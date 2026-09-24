//! The `MeshCompiler` implementation: graph in, meshes out.
//!
//! Evaluation is iterative. The graph forbids non-prior references,
//! so cycles are structurally impossible, but depth is unbounded and
//! recursion would risk a stack overflow on adversarial input.

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, ExecutionOptions, ExecutionTarget, GeomError,
    GeomResult, Operation, ScratchRequirement,
};
use axiolid_core::{PlaneFrame, Point3, Scalar, Tolerance, Transform3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::MeshBoolean;
use axiolid_mesh_compile_contract::{CompileOutcome, MeshCompiler};
use axiolid_model::{GeometryGraph, GeometryNode, NodeId, SolidOperation};

use axiolid_construct::extrude::extrude_profile;
use axiolid_construct::profile::profile_rings;

use crate::channels::{self, Built};

/// Scalar reference compiler.
///
/// Generic over the boolean provider so this crate never depends on a
/// particular one: `axiolid-mesh-boolean-boolmesh` is an adapter, and a different provider
/// swaps in without touching this code.
#[derive(Debug, Clone)]
pub struct ReferenceMeshCompiler<B> {
    boolean: B,
}

impl<B> ReferenceMeshCompiler<B> {
    /// Bind a boolean provider.
    pub const fn new(boolean: B) -> Self {
        Self { boolean }
    }

    /// The bound provider.
    pub const fn boolean(&self) -> &B {
        &self.boolean
    }
}

impl<B: MeshBoolean> Backend for ReferenceMeshCompiler<B> {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(
            BackendId::new("scalar-compile"),
            ExecutionTarget::PortableCpu,
        )
    }
}

/// One entry in the explicit evaluation stack.
///
/// `Enter` schedules dependency discovery; `Exit` runs after every dependency
/// already has a mesh, which is what replaces the recursive call.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct EvalKey {
    id: NodeId,
    linear_bits: u64,
    angular_bits: u64,
}

impl EvalKey {
    fn new(id: NodeId, tolerance: Tolerance) -> Self {
        let bits = |value: Scalar| if value == 0.0 { 0 } else { value.to_bits() };
        Self {
            id,
            linear_bits: bits(tolerance.linear()),
            angular_bits: bits(tolerance.angular()),
        }
    }

    fn tolerance(self) -> Tolerance {
        Tolerance::new(
            Scalar::from_bits(self.linear_bits),
            Scalar::from_bits(self.angular_bits),
        )
        .expect("evaluation keys only contain validated tolerances")
    }
}

#[derive(Debug, Clone, Copy)]
enum Step {
    Enter(EvalKey),
    Exit(EvalKey),
}

/// Cached meshes keyed by node and effective local tolerance.
///
/// A transformed instance changes the local chord budget. Keying only by node
/// would incorrectly reuse a coarse source mesh for a larger instance.
type Cache = std::collections::HashMap<EvalKey, Built>;

impl<B: MeshBoolean> ReferenceMeshCompiler<B> {
    /// Resolve a node handle, blaming the graph rather than panicking.
    /// Resolve a closed 2D boundary curve into rings.
    ///
    /// Only a closed polyline boundary is handled: its points ARE the ring.
    /// An analytic boundary needs the curve evaluator, and an open one does
    /// not bound anything, so both are refused rather than guessed at.
    fn boundary_rings(
        &self,
        graph: &GeometryGraph,
        id: NodeId,
    ) -> GeomResult<axiolid_construct::profile::Rings> {
        let node = self.node(graph, id)?;
        let GeometryNode::Curve2(curve) = node else {
            return Err(GeomError::InvalidInput(format!(
                "half-space boundary {id:?} is not a Curve2 node"
            )));
        };
        match curve {
            axiolid_curve::Curve2::Polyline(p) => {
                let mut pts = p.points.clone();
                // A closed polyline may or may not repeat its first point.
                // Dropping the duplicate keeps the ring's edge count honest.
                if pts.len() >= 2 && pts[0] == pts[pts.len() - 1] {
                    pts.pop();
                }
                if pts.len() < 3 {
                    return Err(GeomError::InvalidInput(
                        "half-space boundary needs at least 3 distinct points".to_owned(),
                    ));
                }
                Ok(axiolid_construct::profile::Rings {
                    outer: pts,
                    holes: Vec::new(),
                })
            }
            _ => Err(GeomError::Unsupported {
                backend: self.descriptor().id,
                operation: Operation::CurveEvaluation,
            }),
        }
    }

    /// Resolve a node that must be a profile, into flattened rings.
    /// Resolve a node that must be a surface.
    ///
    /// Reported by node id rather than by position so a malformed graph
    /// names the offending node instead of the operation that reached it.
    fn surface_of(
        graph: &GeometryGraph,
        id: axiolid_model::NodeId,
    ) -> GeomResult<&axiolid_surface::Surface> {
        match graph.get(id) {
            Some(GeometryNode::Surface(surface)) => Ok(surface),
            Some(_) => Err(GeomError::InvalidInput(format!(
                "reference surface {id:?} is not a Surface node"
            ))),
            None => Err(GeomError::InvalidInput(format!(
                "reference surface {id:?} does not belong to this graph"
            ))),
        }
    }

    fn surface_normals(
        &self,
        graph: &GeometryGraph,
        id: NodeId,
        path: &[Point3],
        options: &ExecutionOptions,
    ) -> GeomResult<Vec<axiolid_core::Vec3>> {
        if let Some(GeometryNode::SurfaceRelation(
            axiolid_model::SurfaceRelation::LinearExtrusion { direction, .. },
        )) = graph.get(id)
        {
            return axiolid_construct::sweep::linear_extrusion_normals(path, *direction);
        }
        let surface = Self::surface_of(graph, id)?;
        path.iter()
            .map(|point| {
                let (u, v) =
                    axiolid_reference::surface::invert(surface, *point, options.tolerance())?;
                axiolid_reference::surface::normal(surface, u, v)
            })
            .collect()
    }

    fn rings_of(
        &self,
        graph: &GeometryGraph,
        id: NodeId,
        options: &ExecutionOptions,
        what: &str,
    ) -> GeomResult<axiolid_construct::profile::Rings> {
        let node = self.node(graph, id)?;
        let GeometryNode::Profile(shape) = node else {
            return Err(GeomError::InvalidInput(format!(
                "{what} {id:?} is not a Profile node"
            )));
        };
        profile_rings(shape, chord_error(options), options.tolerance())
    }

    /// Sample a directrix curve into a polyline.
    ///
    /// Only a polyline directrix is handled here: its points ARE the samples,
    /// so no evaluator is involved. An analytic directrix needs the curve
    /// provider, which the compiler does not yet hold, and is refused with the
    /// named capability rather than silently approximated by its control
    /// points.
    fn directrix_points(
        &self,
        graph: &GeometryGraph,
        id: NodeId,
        range: Option<(Scalar, Scalar)>,
        options: &ExecutionOptions,
    ) -> GeomResult<Vec<Point3>> {
        crate::directrix::points(graph, id, range, options)
    }

    fn node<'g>(&self, graph: &'g GeometryGraph, id: NodeId) -> GeomResult<&'g GeometryNode> {
        graph.get(id).ok_or_else(|| {
            GeomError::InvalidInput(format!("node {id:?} does not belong to this graph"))
        })
    }

    /// Nodes that must have meshes before `id` can be evaluated.
    ///
    /// Only mesh-producing dependencies are listed. A profile referenced by an
    /// extrusion is consumed as 2D data, not as a mesh, so it is deliberately
    /// absent: compiling it standalone would be meaningless.
    fn mesh_dependencies(
        &self,
        graph: &GeometryGraph,
        node: &GeometryNode,
        key: EvalKey,
    ) -> GeomResult<Vec<EvalKey>> {
        let tolerance = key.tolerance();
        let same = |id| EvalKey::new(id, tolerance);
        Ok(match node {
            GeometryNode::Instance(instance) => vec![EvalKey::new(
                instance.source,
                instance_local_tolerance(instance.transform, tolerance)?,
            )],
            GeometryNode::Collection(members) => members.iter().copied().map(same).collect(),
            GeometryNode::SolidOperation(
                operation @ SolidOperation::Boolean { left, right, .. },
            ) => {
                if is_subject_bounded_half_space_boolean(graph, operation) {
                    vec![same(*left)]
                } else {
                    vec![same(*left), same(*right)]
                }
            }
            _ => Vec::new(),
        })
    }

    /// Iterative post-order evaluation of one root.
    fn evaluate(
        &self,
        graph: &GeometryGraph,
        root: NodeId,
        options: &ExecutionOptions,
        cache: &mut Cache,
    ) -> GeomResult<Built> {
        let root_key = EvalKey::new(root, options.tolerance());
        let mut stack = vec![Step::Enter(root_key)];
        while let Some(step) = stack.pop() {
            match step {
                Step::Enter(key) => {
                    if cache.contains_key(&key) {
                        continue;
                    }
                    let node = self.node(graph, key.id)?;
                    let deps = self.mesh_dependencies(graph, node, key)?;
                    // Exit runs after every dependency, so push it first.
                    stack.push(Step::Exit(key));
                    for dep in deps {
                        if !cache.contains_key(&dep) {
                            stack.push(Step::Enter(dep));
                        }
                    }
                }
                Step::Exit(key) => {
                    if cache.contains_key(&key) {
                        continue;
                    }
                    let local_options = options.clone().with_tolerance(key.tolerance());
                    let built = self.build(graph, key, &local_options, cache)?;
                    cache.insert(key, built);
                }
            }
        }
        cache
            .get(&root_key)
            .cloned()
            .ok_or_else(|| GeomError::InvalidInput(format!("root {root:?} produced no mesh")))
    }

    /// Convert an authored polygon mesh to triangles (#160).
    ///
    /// A plain triangle keeps its exact corner order. Any other face -- an
    /// n-gon, a concave face, a face with holes (IFC4
    /// `IfcIndexedPolygonalFaceWithVoids`) -- is triangulated in its own
    /// plane by [`crate::planar::triangulate_polygon`]. Positions are kept
    /// as authored and shared, so the output welds exactly where the input
    /// did; no corner is moved or added.
    ///
    /// A face that is not planar within the linear tolerance, has no area,
    /// or whose rings cross is refused with an error naming its index: a
    /// non-planar n-gon has no unique triangulation, so picking one would
    /// invent geometry.
    fn compile_authored_polygons(
        &self,
        mesh: &axiolid_mesh::PolygonMesh,
        options: &ExecutionOptions,
    ) -> GeomResult<TriMesh> {
        let position_count = mesh.positions.len();
        if let Some(index) = mesh
            .faces
            .iter()
            .flat_map(|face| face_rings(face).flatten().copied())
            .find(|&index| index as usize >= position_count)
        {
            return Err(GeomError::InvalidInput(format!(
                "authored polygon index {index} exceeds position count {position_count}"
            )));
        }

        let mut positions = Vec::new();
        positions
            .try_reserve_exact(position_count)
            .map_err(|_| GeomError::BudgetExceeded {
                resource: "authored polygon positions",
            })?;
        positions.extend_from_slice(&mesh.positions);

        // A polygon with k corners (holes included) gives at most
        // k + 2 * holes - 2 triangles; bound by corner count times three.
        let index_bound = mesh
            .faces
            .iter()
            .try_fold(0_usize, |total, face| {
                let corners = face_rings(face).map(Vec::len).sum::<usize>();
                corners
                    .checked_add(2 * face.holes.len())
                    .and_then(|triangles| triangles.checked_mul(3))
                    .and_then(|indices| total.checked_add(indices))
            })
            .ok_or(GeomError::BudgetExceeded {
                resource: "authored polygon indices",
            })?;
        let mut indices = Vec::new();
        indices
            .try_reserve_exact(index_bound)
            .map_err(|_| GeomError::BudgetExceeded {
                resource: "authored polygon indices",
            })?;

        let linear = options.tolerance().linear();
        for (face_index, face) in mesh.faces.iter().enumerate() {
            if face.outer.len() == 3 && face.holes.is_empty() {
                indices.extend_from_slice(&face.outer);
                continue;
            }
            // A repeated corner (an exporter's closing point, a doubled
            // corner) needs no handling here: earcut drops coincident
            // consecutive points itself (`a_repeated_closing_corner_is_not_a_triangle`).
            let rings: Vec<&Vec<u32>> = face_rings(face).collect();
            let points: Vec<Vec<axiolid_core::Point3>> = rings
                .iter()
                .map(|ring| ring.iter().map(|&i| mesh.positions[i as usize]).collect())
                .collect();
            let views: Vec<&[axiolid_core::Point3]> = points.iter().map(Vec::as_slice).collect();
            let local = crate::planar::triangulate_polygon(&views, linear)
                .map_err(|refusal| crate::planar::face_error(face_index, refusal))?;
            let corners: Vec<u32> = rings.into_iter().flatten().copied().collect();
            indices.extend(local.into_iter().map(|corner| corners[corner]));
        }

        let triangles = TriMesh::new(positions, indices);
        triangles.validate_structure().map_err(|error| {
            GeomError::InvalidInput(format!("invalid authored polygon mesh: {error}"))
        })?;
        Ok(triangles)
    }

    /// Build one node, assuming its mesh dependencies are already cached.
    fn build(
        &self,
        graph: &GeometryGraph,
        key: EvalKey,
        options: &ExecutionOptions,
        cache: &Cache,
    ) -> GeomResult<Built> {
        let id = key.id;
        let node = self.node(graph, id)?;
        match node {
            GeometryNode::TriMesh(mesh) => Ok(authored_mesh(mesh.clone(), options)),
            GeometryNode::PolygonMesh(mesh) => self
                .compile_authored_polygons(mesh, options)
                .map(|mesh| authored_mesh(mesh, options)),
            GeometryNode::Instance(instance) => {
                let source_tolerance =
                    instance_local_tolerance(instance.transform, options.tolerance())?;
                let source = self.cached(cache, instance.source, source_tolerance)?;
                Ok(channels::transform(source, instance.transform))
            }
            GeometryNode::Collection(members) => {
                let members = members
                    .iter()
                    .map(|&member| self.cached(cache, member, options.tolerance()))
                    .collect::<GeomResult<Vec<_>>>()?;
                Ok(channels::merge(&members))
            }
            GeometryNode::SolidOperation(SolidOperation::Boolean {
                left,
                right,
                operator,
            }) => self.build_boolean(graph, *left, *right, *operator, options, cache),
            GeometryNode::SolidOperation(operation) => {
                self.build_solid(graph, operation, options).map(Built::leaf)
            }
            GeometryNode::BRep(brep) => crate::brep::tessellate(brep, graph, options.tolerance())
                .map(|(mesh, closure)| Built::with_closure(mesh, closure)),
            // CSG primitives are analytic solids: no surface evaluation,
            // no trim curves, just a closed mesh at the caller's tolerance.
            GeometryNode::Primitive(primitive) => {
                axiolid_reference::primitive::tessellate_primitive(primitive, options.tolerance())
                    .map(Built::leaf)
            }
            other => Err(GeomError::Unsupported {
                backend: self.descriptor().id,
                operation: unsupported_operation(other),
            }),
        }
    }

    /// Read an already-built dependency.
    fn cached<'c>(
        &self,
        cache: &'c Cache,
        id: NodeId,
        tolerance: Tolerance,
    ) -> GeomResult<&'c Built> {
        cache.get(&EvalKey::new(id, tolerance)).ok_or_else(|| {
            GeomError::InvalidInput(format!(
                "dependency {id:?} was not evaluated first at tolerance {tolerance:?}"
            ))
        })
    }

    /// A boolean, with the operands' channel fates composed onto the
    /// provider's (#115): a channel lost upstream stays reported lost.
    fn build_boolean(
        &self,
        graph: &GeometryGraph,
        left: NodeId,
        right: NodeId,
        operator: axiolid_core::BooleanOperator,
        options: &ExecutionOptions,
        cache: &Cache,
    ) -> GeomResult<Built> {
        let subject = self.cached(cache, left, options.tolerance())?;
        refuse_surface_operand(subject, "subject")?;
        let bounded_tool = match self.node(graph, right)? {
            GeometryNode::HalfSpace(hs) => Some(axiolid_construct::half_space::for_subject(
                &subject.mesh,
                *hs,
                options.tolerance(),
            )?),
            _ => None,
        };
        // A bounded half-space is built here from the subject, so it has no
        // upstream fates of its own.
        let (tool, tool_fates) = match bounded_tool.as_ref() {
            Some(tool) => (tool, None),
            None => {
                let built = self.cached(cache, right, options.tolerance())?;
                refuse_surface_operand(built, "tool")?;
                (&built.mesh, Some(&built.fates))
            }
        };
        let outcome = self
            .boolean
            .boolean(&subject.mesh, tool, operator, options)?;
        Ok(channels::after_boolean(
            outcome.mesh,
            &subject.fates,
            tool_fates,
            outcome.evidence.attribute_fates,
        ))
    }

    /// Every non-boolean solid family; unsupported ones are explicitly
    /// refused. Booleans go through [`Self::build_boolean`], which needs the
    /// operands' cached fates.
    fn build_solid(
        &self,
        graph: &GeometryGraph,
        operation: &SolidOperation,
        options: &ExecutionOptions,
    ) -> GeomResult<TriMesh> {
        match operation {
            SolidOperation::Extrusion {
                profile,
                direction,
                depth,
            } => {
                let node = self.node(graph, *profile)?;
                let GeometryNode::Profile(shape) = node else {
                    return Err(GeomError::InvalidInput(format!(
                        "extrusion profile {profile:?} is not a Profile node"
                    )));
                };
                let rings = profile_rings(shape, chord_error(options), options.tolerance())?;
                extrude_profile(&rings, *direction, *depth, options.tolerance())
            }
            SolidOperation::Revolution {
                profile,
                axis_origin,
                axis_direction,
                angle,
            } => {
                let node = self.node(graph, *profile)?;
                let GeometryNode::Profile(shape) = node else {
                    return Err(GeomError::InvalidInput(format!(
                        "revolution profile {profile:?} is not a Profile node"
                    )));
                };
                let rings = profile_rings(shape, chord_error(options), options.tolerance())?;
                axiolid_construct::revolve::revolve(
                    &rings,
                    *axis_origin,
                    *axis_direction,
                    *angle,
                    options.tolerance(),
                )
            }
            SolidOperation::TaperedExtrusion {
                start_profile,
                end_profile,
                direction,
                depth,
            } => {
                let a = self.rings_of(graph, *start_profile, options, "taper start profile")?;
                let b = self.rings_of(graph, *end_profile, options, "taper end profile")?;
                axiolid_construct::sweep::tapered_extrude(&a, &b, *direction, *depth)
            }
            SolidOperation::TaperedRevolution {
                start_profile,
                end_profile,
                axis_origin,
                axis_direction,
                angle,
            } => {
                let a = self.rings_of(graph, *start_profile, options, "taper start profile")?;
                let b = self.rings_of(graph, *end_profile, options, "taper end profile")?;
                axiolid_construct::sweep::tapered_revolve(
                    &a,
                    &b,
                    *axis_origin,
                    *axis_direction,
                    *angle,
                    options.tolerance(),
                )
            }
            SolidOperation::SweptDisk {
                directrix,
                radius,
                inner_radius,
                parameter_range,
                fillet_radius,
            } => {
                let path = self.directrix_points(graph, *directrix, *parameter_range, options)?;
                axiolid_construct::sweep::swept_disk(
                    &path,
                    *radius,
                    *inner_radius,
                    *fillet_radius,
                    options.tolerance(),
                )
            }
            SolidOperation::FixedReferenceSweep {
                profile,
                directrix,
                reference_direction,
                parameter_range,
            } => {
                let rings = self.rings_of(graph, *profile, options, "sweep profile")?;
                let path = self.directrix_points(graph, *directrix, *parameter_range, options)?;
                axiolid_construct::sweep::fixed_reference_sweep(&rings, &path, *reference_direction)
            }
            SolidOperation::SurfaceCurveSweep {
                profile,
                directrix,
                reference_surface,
                parameter_range,
            } => {
                let rings =
                    self.rings_of(graph, *profile, options, "surface curve sweep profile")?;
                let path = self.directrix_points(graph, *directrix, *parameter_range, options)?;
                let normals = self.surface_normals(graph, *reference_surface, &path, options)?;
                axiolid_construct::sweep::surface_curve_sweep(&rings, &path, &normals)
            }
            SolidOperation::SectionedSpine { spine, sections } => {
                let path = self.directrix_points(graph, *spine, None, options)?;
                if sections.len() != path.len() {
                    return Err(GeomError::InvalidInput(format!(
                        "a sectioned spine needs one section per spine point: {} sections, {} points",
                        sections.len(),
                        path.len()
                    )));
                }
                let mut placed = Vec::with_capacity(sections.len());
                for (section, origin) in sections.iter().zip(&path) {
                    let rings =
                        self.rings_of(graph, section.profile, options, "spine section profile")?;
                    // The section's own placement positions its profile;
                    // the spine point supplies the station origin.
                    let pts = rings
                        .outer
                        .iter()
                        .chain(rings.holes.iter().flatten())
                        .map(|p| {
                            section
                                .placement
                                .transform_point3(Point3::new(p.x, p.y, 0.0))
                                + *origin
                        })
                        .collect();
                    placed.push((rings, pts));
                }
                axiolid_construct::sweep::sectioned_spine(&placed)
            }
            SolidOperation::BoundedHalfSpace {
                half_space,
                boundary,
                placement,
            } => {
                let node = self.node(graph, *half_space)?;
                let GeometryNode::HalfSpace(hs) = node else {
                    return Err(GeomError::InvalidInput(format!(
                        "half-space {half_space:?} is not a HalfSpace node"
                    )));
                };
                let rings = self.boundary_rings(graph, *boundary)?;
                // The declared margin is the contract's own knob for how far
                // an unbounded half-space extends before it can be meshed.
                let margin = axiolid_primitive::ClipMargin::new(2.0)
                    .expect("2.0 is a valid positive clip margin");
                // `placement` is the boundary's OWN frame, independent of the
                // clip plane: the profile is authored in it, so it has to
                // orient the profile before the slab is built. Applying it
                // only to the finished mesh cannot express an in-plane
                // rotation, because by then the boundary has already been
                // framed against the plane normal.
                let frame = PlaneFrame::new(
                    placement.translation,
                    placement.matrix3.x_axis.normalize_or_zero(),
                    placement.matrix3.y_axis.normalize_or_zero(),
                    options.tolerance(),
                )
                .map_err(|error| {
                    GeomError::InvalidInput(format!(
                        "bounded half-space placement is not a usable boundary frame: {error}"
                    ))
                })?;
                let mesh = axiolid_construct::half_space::bounded_half_space_in_frame(
                    &rings,
                    hs.boundary,
                    frame,
                    hs.agreement,
                    margin,
                    options.tolerance(),
                )?;
                Ok(mesh)
            }
            // Routed to `build_boolean` by `build`; unreachable here.
            SolidOperation::Boolean { .. } => Err(GeomError::InvalidInput(
                "boolean must be built through build_boolean".into(),
            )),
            // Naming the capability lets a caller register a provider for it
            // rather than guess. `Unsupported` names only the operation, which
            // collapses revolution, swept disk, fixed-reference sweep and the
            // rest into one indistinguishable answer -- a caller learns that
            // "a sweep" failed but not WHICH provider is missing, so it cannot
            // act. `UnsupportedInput` carries the family, matching what the
            // exact path already reports.
            other => Err(GeomError::UnsupportedInput {
                backend: self.descriptor().id,
                operation: Operation::Sweep,
                input: crate::exact::solid_operation_family(other),
            }),
        }
    }
}

/// Convert a world-space tolerance into conservative instance-local space.
fn instance_local_tolerance(transform: Transform3, tolerance: Tolerance) -> GeomResult<Tolerance> {
    let m = transform.matrix3;
    let sx = m.x_axis.length();
    let sy = m.y_axis.length();
    let sz = m.z_axis.length();
    let max_scale = sx.max(sy).max(sz);
    if !max_scale.is_finite() || max_scale == 0.0 {
        return Err(GeomError::InvalidInput(
            "instance transform has no finite scale".into(),
        ));
    }
    let eps = 32.0 * f64::EPSILON;
    let orthogonal = m.x_axis.dot(m.y_axis).abs() <= eps * sx * sy
        && m.x_axis.dot(m.z_axis).abs() <= eps * sx * sz
        && m.y_axis.dot(m.z_axis).abs() <= eps * sy * sz;
    let stretch = if orthogonal {
        max_scale * (1.0 + 3.0 * eps)
    } else {
        (sx * sx + sy * sy + sz * sz).sqrt()
    };
    if !stretch.is_finite() {
        return Err(GeomError::InvalidInput(
            "instance transform scale overflow".into(),
        ));
    }
    Tolerance::new(tolerance.linear() / stretch, tolerance.angular())
        .map_err(|error| GeomError::InvalidInput(error.to_string()))
}

fn chord_error(options: &ExecutionOptions) -> Scalar {
    options.tolerance().linear()
}

fn is_subject_bounded_half_space_boolean(
    graph: &GeometryGraph,
    operation: &SolidOperation,
) -> bool {
    let SolidOperation::Boolean {
        right, operator, ..
    } = operation
    else {
        return false;
    };
    // These operators stay within the finite left-hand subject, so a prism
    // covering its bounds is an exact finite stand-in for the half-space.
    // Union and XOR are unbounded and must remain unsupported.
    matches!(
        operator,
        axiolid_core::BooleanOperator::Difference | axiolid_core::BooleanOperator::Intersection
    ) && matches!(graph.get(*right), Some(GeometryNode::HalfSpace(_)))
}

/// Which capability a node family would need.
///
/// Reporting the real missing capability lets a caller register a provider
/// that supplies it, instead of guessing from a generic failure.
fn unsupported_operation(node: &GeometryNode) -> Operation {
    match node {
        GeometryNode::Curve2(_) | GeometryNode::Curve3(_) | GeometryNode::OpenProfile(_) => {
            Operation::CurveEvaluation
        }
        GeometryNode::Surface(_) => Operation::SurfaceEvaluation,
        GeometryNode::Profile(_) => Operation::ProfileTriangulation,
        GeometryNode::BRep(_) | GeometryNode::PolygonMesh(_) => Operation::Tessellation,
        _ => Operation::GraphCompilation,
    }
}

impl<B: MeshBoolean> MeshCompiler for ReferenceMeshCompiler<B> {
    /// Bounded by the peak mesh size, which is data-dependent, so the honest
    /// answer is unbounded rather than an invented constant.
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::Unbounded
    }

    fn compile_mesh(
        &self,
        graph: &GeometryGraph,
        root: NodeId,
        options: &ExecutionOptions,
    ) -> GeomResult<TriMesh> {
        self.admit_budget(options)?;
        let mut cache = Cache::new();
        self.evaluate(graph, root, options, &mut cache)
            .map(|built| built.mesh)
    }

    /// Every channel on the result, and every channel an input had that the
    /// result lost, is reported with its fate (#115).
    fn compile_mesh_reported(
        &self,
        graph: &GeometryGraph,
        root: NodeId,
        options: &ExecutionOptions,
    ) -> GeomResult<CompileOutcome> {
        self.admit_budget(options)?;
        let mut cache = Cache::new();
        let built = self.evaluate(graph, root, options, &mut cache)?;
        Ok(CompileOutcome::tracked(built.mesh, built.fates.into_vec()).with_closure(built.closure))
    }

    /// Overriding the `_into` seam gives both call shapes one shared cache,
    /// so a subtree referenced by several roots is compiled once per batch.
    fn compile_mesh_batch_into(
        &self,
        graph: &GeometryGraph,
        roots: &[NodeId],
        options: &ExecutionOptions,
        destination: &mut Vec<TriMesh>,
    ) -> GeomResult<()> {
        self.admit_budget(options)?;
        destination.reserve(roots.len());
        let mut cache = Cache::new();
        for &root in roots {
            destination.push(self.evaluate(graph, root, options, &mut cache)?.mesh);
        }
        Ok(())
    }
}

impl<B: MeshBoolean> ReferenceMeshCompiler<B> {
    /// Refuse a compilation whose memory budget this compiler cannot honour.
    ///
    /// Graph compilation caches every intermediate mesh, and peak size is
    /// data-dependent, so `scratch_requirement` honestly reports
    /// [`ScratchRequirement::Unbounded`]. By contract an unbounded requirement
    /// never fits a DECLARED budget -- admitting it anyway would make the
    /// budget advisory, which is the failure the type exists to prevent.
    ///
    /// Before this check the field was simply ignored here: a caller could set
    /// a budget and watch the compiler allocate past it without a word. A
    /// caller that sets no budget is unaffected.
    fn admit_budget(&self, options: &ExecutionOptions) -> GeomResult<()> {
        if self.scratch_requirement().fits_budget(options, 0) {
            return Ok(());
        }
        Err(GeomError::BudgetExceeded {
            resource: "graph compilation intermediate meshes",
        })
    }

    /// The boolean provider this compiler dispatches to.
    ///
    /// Exposed so an application can apply its own source-format set
    /// operations with the same provider the compiler uses, rather than
    /// constructing a second one that might differ.
    pub const fn boolean_provider(&self) -> &B {
        &self.boolean
    }
}

/// A face's outer ring followed by its holes.
fn face_rings(face: &axiolid_mesh::PolygonFace) -> impl Iterator<Item = &Vec<u32>> {
    std::iter::once(&face.outer).chain(face.holes.iter())
}

/// A boolean needs two volumes; a surface model has none (#161).
///
/// Refused rather than handed to the provider, which would see a closed
/// surface model as a valid solid and return a confident, meaningless
/// result.
fn refuse_surface_operand(built: &Built, role: &'static str) -> GeomResult<()> {
    if built.closure == axiolid_mesh_compile_contract::MeshClosure::Solid {
        return Ok(());
    }
    Err(GeomError::InvalidInput(format!(
        "boolean {role} is a surface model: it encloses no volume to combine"
    )))
}

/// An authored mesh, with its closure read from its structure (#161).
///
/// A mesh node carries no "this is a solid" declaration the way a B-rep
/// does, so the geometry is the only evidence: a closed, consistently
/// wound two-manifold bounds a solid, anything else (an open face set, a
/// single sheet) is a surface.
fn authored_mesh(mesh: TriMesh, options: &ExecutionOptions) -> Built {
    let closure = if axiolid_mesh::audit_mesh(&mesh, options.tolerance()).is_closed_two_manifold() {
        axiolid_mesh_compile_contract::MeshClosure::Solid
    } else {
        axiolid_mesh_compile_contract::MeshClosure::Surface
    };
    Built::with_closure(mesh, closure)
}
