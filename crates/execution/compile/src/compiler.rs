//! The `MeshCompiler` implementation: graph in, meshes out.
//!
//! Evaluation is iterative. The graph forbids non-prior references,
//! so cycles are structurally impossible, but depth is unbounded and
//! recursion would risk a stack overflow on adversarial input.

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, ExecutionOptions, ExecutionTarget, GeomError,
    GeomResult, Operation, ScratchRequirement,
};
use axiolid_core::{Point3, Scalar, Tolerance, Transform3, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::MeshBoolean;
use axiolid_mesh_compile_contract::{
    CompileOutcome, MeshCompiler, PolylineCompiler, PolylineOutcome,
};
use axiolid_model::{GeometryGraph, GeometryNode, NodeId, SolidOperation};

use crate::half_space_boundary::PrismAxis;

use axiolid_construct::extrude::extrude_profile;
use axiolid_construct::profile::profile_rings;

use crate::channels::{self, Built};
use crate::station::{PlacementAccuracy, UNCERTIFIED_STATION_FRAME};

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

/// Why the mesh cut of a boolean was refused.
#[derive(Debug)]
enum CutRefusal {
    /// The result touches itself (#194) or holds two bodies within the
    /// tolerance (ADR 0096, D2): the operands meet where the mesh boolean,
    /// which has no tolerance, cannot decide. The exact referee may.
    Contact(GeomError),
    /// Any other refusal.
    Other(GeomError),
}

impl CutRefusal {
    fn into_error(self) -> GeomError {
        match self {
            Self::Contact(error) | Self::Other(error) => error,
        }
    }
}

/// A boolean's mesh cut, and whether its faces include a pair that nearly
/// share a plane ([`axiolid_inspect::near_coplanar`]): contact the mesh
/// boolean did not decide, whatever it returned (#391).
struct MeshCut {
    built: Built,
    unsettled: bool,
    /// The result holds a skin thinner than the rounding reach (#400, ADR
    /// 0102): the faces a clip left a residue apart, which the snap could
    /// not close.
    skin: Option<PendingSkin>,
}

/// A skin within rounding a cut's result holds: the refusal that names it
/// and the rounding reach it was read at.
#[derive(Debug, Clone)]
pub(crate) struct PendingSkin {
    pub error: GeomError,
    pub reach: Scalar,
}

/// Who decided a boolean's result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Decided {
    /// The mesh boolean.
    Mesh,
    /// The exact compiler, tessellated ([`crate::referee`]).
    Referee,
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
    /// The chord budget in this node's local space (#165). Part of the
    /// key because an instance scales it with its transform, and two
    /// budgets tessellate one node differently.
    chord_bits: u64,
    /// The placed coordinate scale (#312): part of the key because a
    /// boolean's snap reaches the rounding of the coordinates its result is
    /// placed at, so a source placed near and far compiles once per scale.
    placed_bits: u64,
}

/// The evaluation budgets that vary per node: the tolerance, and the chord
/// budget curved geometry is flattened to (#165). Both rescale together
/// through an instance transform. And the placed coordinate scale a
/// boolean's snap reaches (#312, [`crate::snap::placed_scale`]): `0` where
/// nothing places the node, else the world coordinate magnitude its
/// placement adds, in the node's own units.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Budget {
    tolerance: Tolerance,
    chord: Scalar,
    placed: Scalar,
}

impl Budget {
    fn of(options: &ExecutionOptions) -> Self {
        let tolerance = options.tolerance();
        Self {
            tolerance,
            chord: options.chord_error().unwrap_or(tolerance.linear()),
            placed: 0.0,
        }
    }
}

impl EvalKey {
    fn new(id: NodeId, budget: Budget) -> Self {
        let bits = |value: Scalar| if value == 0.0 { 0 } else { value.to_bits() };
        Self {
            id,
            linear_bits: bits(budget.tolerance.linear()),
            angular_bits: bits(budget.tolerance.angular()),
            chord_bits: bits(budget.chord),
            placed_bits: bits(budget.placed),
        }
    }

    fn budget(self) -> Budget {
        Budget {
            tolerance: self.tolerance(),
            chord: Scalar::from_bits(self.chord_bits),
            placed: Scalar::from_bits(self.placed_bits),
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

/// Whether a compilation measures what no construction proves: a
/// boolean's deviation against its exact result (#235) costs an exact
/// compilation and a certified search, so only a deviation report pays it,
/// and only for the booleans whose result is emitted: an operand of
/// another boolean is replaced by that boolean's own measurement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Measure<'m> {
    Skip,
    /// Measure the booleans in this set
    /// ([`crate::deviation::emitted_booleans`]).
    Deviation(&'m std::collections::HashSet<NodeId>),
}

/// Cached meshes keyed by node and effective local tolerance.
///
/// A transformed instance changes the local chord budget. Keying only by node
/// would incorrectly reuse a coarse source mesh for a larger instance.
///
/// It also remembers which nodes are regions ([`Regions`]), so a batch of
/// roots decides each node once.
#[derive(Debug, Default)]
struct Cache {
    built: std::collections::HashMap<EvalKey, Built>,
    regions: crate::regions::RegionMemo,
}

impl Cache {
    fn new() -> Self {
        Self::default()
    }
}

impl std::ops::Deref for Cache {
    type Target = std::collections::HashMap<EvalKey, Built>;

    fn deref(&self) -> &Self::Target {
        &self.built
    }
}

impl std::ops::DerefMut for Cache {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.built
    }
}

impl<B: MeshBoolean> ReferenceMeshCompiler<B> {
    /// Resolve a bounded half-space's boundary into rings.
    ///
    /// A closed polyline's points ARE the ring, taken as wound. A profile
    /// (#277) is flattened as an extrusion's is, every arc chorded within
    /// `chord` by the certified flattener, the rings oriented outer
    /// counter-clockwise; an open, self-crossing or degenerate one is
    /// refused there or by the cap triangulation, by name. Any other curve
    /// needs the curve evaluator, and an open one does not bound anything,
    /// so both are refused rather than guessed at.
    fn boundary_rings(
        &self,
        graph: &GeometryGraph,
        id: NodeId,
        chord: Scalar,
        options: &ExecutionOptions,
    ) -> GeomResult<axiolid_construct::profile::Rings> {
        let node = self.node(graph, id)?;
        let curve = match node {
            GeometryNode::Curve2(curve) => curve,
            GeometryNode::Profile(shape) => {
                crate::half_space_boundary::check(shape, options.tolerance())?;
                return profile_rings(shape, chord, options.tolerance());
            }
            _ => {
                return Err(GeomError::InvalidInput(format!(
                    "half-space boundary {id:?} is neither a Curve2 nor a Profile node"
                )))
            }
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

    /// The slab of a bounded half-space, `depth` deep on its kept side.
    ///
    /// `placement` is the boundary's OWN frame, independent of the clip
    /// plane: the profile is authored in it, so it orients the profile
    /// before the slab is built. Applying it only to the finished mesh
    /// cannot express an in-plane rotation, because by then the boundary has
    /// already been framed against the plane normal. The whole frame is
    /// honoured (#322): along the plane normal its `x` is projected into the
    /// plane and its scale, shear and handedness kept
    /// ([`BoundaryPrism::along_plane_normal`]); along the placement's `z`
    /// (#321) it is read as written ([`BoundaryPrism::along_placement_z`]),
    /// the prism at any angle to the plane. Either way the boundary is
    /// chorded within the budget over the frame's longest axis, so the
    /// chords stay within it once placed.
    ///
    /// [`BoundaryPrism::along_plane_normal`]: axiolid_construct::half_space::BoundaryPrism::along_plane_normal
    /// [`BoundaryPrism::along_placement_z`]: axiolid_construct::half_space::BoundaryPrism::along_placement_z
    fn bounded_half_space(
        &self,
        graph: &GeometryGraph,
        (half_space, boundary): (NodeId, NodeId),
        placement: Transform3,
        prism_axis: PrismAxis,
        depth: axiolid_construct::half_space::SlabDepth,
        options: &ExecutionOptions,
    ) -> GeomResult<TriMesh> {
        use axiolid_construct::half_space::{bounded_half_space_prism, BoundaryPrism};
        let node = self.node(graph, half_space)?;
        let GeometryNode::HalfSpace(hs) = node else {
            return Err(GeomError::InvalidInput(format!(
                "half-space {half_space:?} is not a HalfSpace node"
            )));
        };
        let tolerance = options.tolerance();
        let prism = match prism_axis {
            PrismAxis::PlaneNormal => {
                BoundaryPrism::along_plane_normal(placement, hs.boundary, tolerance)?
            }
            PrismAxis::PlacementZ => BoundaryPrism::along_placement_z(placement, tolerance)?,
        };
        let stretch = prism.x.length().max(prism.y.length());
        let chord = chord_error(options) / stretch.max(1.0);
        let rings = self.boundary_rings(graph, boundary, chord, options)?;
        bounded_half_space_prism(&rings, hs.boundary, prism, hs.agreement, depth, tolerance)
    }

    /// The finite stand-in for a half-space `tool`, placed by `placement`,
    /// over `subject`: the tool is unbounded on its kept side, and a
    /// difference or an intersection with it stays within the subject, so
    /// a slab reaching past the subject is the half-space for them.
    ///
    /// The subject is brought into the tool's own frame. An unbounded
    /// half-space covers it as before (`for_subject`); a bounded one is cut
    /// off past the subject's farthest point on its kept side, by a tenth of
    /// the subject's span plus the tolerance (#325): sized by its boundary
    /// alone, a small boundary over a tall subject ended inside it and cut
    /// it too little, with no error. The slab is placed back, its boundary
    /// chorded within the budget over the placement's stretch.
    fn clip_tool_mesh(
        &self,
        graph: &GeometryGraph,
        tool: ClipTool,
        placement: Transform3,
        subject: &TriMesh,
        options: &ExecutionOptions,
    ) -> GeomResult<TriMesh> {
        let tolerance = options.tolerance();
        let placed = placement != Transform3::IDENTITY;
        let inverse = placement.inverse();
        if placed
            && !(placement.matrix3.determinant() != 0.0
                && inverse.matrix3.is_finite()
                && inverse.translation.is_finite())
        {
            return Err(GeomError::InvalidInput(
                "half-space placed by a singular transform".to_owned(),
            ));
        }
        let local_subject;
        let subject = if placed {
            local_subject = TriMesh::new(
                subject
                    .positions
                    .iter()
                    .map(|&p| inverse.transform_point3(p))
                    .collect(),
                subject.indices.clone(),
            );
            &local_subject
        } else {
            subject
        };
        let local = match tool {
            ClipTool::Unbounded(half_space) => {
                axiolid_construct::half_space::for_subject(subject, half_space, tolerance)?
            }
            ClipTool::Bounded {
                half_space,
                boundary,
                placement: frame,
                prism_axis,
            } => {
                let Some(GeometryNode::HalfSpace(hs)) = graph.get(half_space) else {
                    return Err(GeomError::InvalidInput(format!(
                        "half-space {half_space:?} is not a HalfSpace node"
                    )));
                };
                let normal = hs.boundary.normal.normalize_or_zero();
                if normal == Vec3::ZERO || subject.positions.is_empty() {
                    return Err(GeomError::InvalidInput(
                        "bounded half-space boolean needs a plane and a subject with bounds".into(),
                    ));
                }
                let kept = if hs.agreement { normal } else { -normal };
                let (mut lo, mut hi) = (
                    Vec3::splat(Scalar::INFINITY),
                    Vec3::splat(Scalar::NEG_INFINITY),
                );
                let mut farthest = Scalar::NEG_INFINITY;
                for &point in &subject.positions {
                    if !point.is_finite() {
                        return Err(GeomError::InvalidInput(
                            "bounded half-space subject contains non-finite points".into(),
                        ));
                    }
                    lo = lo.min(point);
                    hi = hi.max(point);
                    farthest = farthest.max((point - hs.boundary.origin).dot(kept));
                }
                let span = (hi - lo).max_element().max(tolerance.linear());
                let reach = farthest.max(0.0) + 0.1 * span + tolerance.linear();
                let local_options;
                let options = if placed {
                    let chord = chord_error(options) / crate::deviation::stretch(placement);
                    local_options = options.clone().with_chord_error(chord).ok_or_else(|| {
                        GeomError::InvalidInput(format!(
                            "chord budget {chord} is not positive and finite"
                        ))
                    })?;
                    &local_options
                } else {
                    options
                };
                self.bounded_half_space(
                    graph,
                    (half_space, boundary),
                    frame,
                    prism_axis,
                    axiolid_construct::half_space::SlabDepth::Reach(reach),
                    options,
                )?
            }
        };
        Ok(if placed {
            channels::transform(&Built::leaf(local), placement, 0.0).mesh
        } else {
            local
        })
    }

    /// The finite stand-in for a region over `subject` (#326): the first
    /// leaf's slab is sized by the subject, as a single clip's is, and each
    /// further leaf's by the stand-in so far, then cut together. `None`
    /// when the leaves share nothing within the subject's reach: an empty
    /// tool. With one leaf it is `clip_tool_mesh`, unchanged.
    fn region_tool_mesh(
        &self,
        graph: &GeometryGraph,
        leaves: &[(ClipTool, Transform3)],
        subject: &TriMesh,
        options: &ExecutionOptions,
        placed: Scalar,
    ) -> GeomResult<Option<Built>> {
        let mut tool: Option<Built> = None;
        for &(leaf, placement) in leaves {
            let around = tool.as_ref().map_or(subject, |built| &built.mesh);
            let side = self.clip_tool_mesh(graph, leaf, placement, around, options)?;
            tool = Some(match tool {
                None => Built::leaf(side),
                Some(so_far) => {
                    let cut = self
                        .boolean_of(
                            &so_far,
                            (&side, None, None),
                            axiolid_core::BooleanOperator::Intersection,
                            options,
                            placed,
                        )
                        .map_err(CutRefusal::into_error)?;
                    if let Some(skin) = cut.skin {
                        return Err(skin.error);
                    }
                    if cut.built.mesh.indices.is_empty() {
                        return Ok(None);
                    }
                    cut.built
                }
            });
        }
        Ok(tool)
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
            // A relation is a model node with no `Surface` value to
            // evaluate; a refusal by name, not a malformed graph (#402).
            Some(GeometryNode::SurfaceRelation(_)) => Err(GeomError::UnsupportedInput {
                backend: crate::BACKEND_ID,
                operation: axiolid_contracts::Operation::SurfaceEvaluation,
                input: "a surface relation as a reference surface",
            }),
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

    /// A profile's rings, chorded to `chord`.
    fn rings_of(
        &self,
        graph: &GeometryGraph,
        id: NodeId,
        chord: Scalar,
        options: &ExecutionOptions,
        what: &str,
    ) -> GeomResult<axiolid_construct::profile::Rings> {
        let node = self.node(graph, id)?;
        let GeometryNode::Profile(shape) = node else {
            return Err(GeomError::InvalidInput(format!(
                "{what} {id:?} is not a Profile node"
            )));
        };
        profile_rings(shape, chord, options.tolerance())
    }

    /// A directrix sampled to the chord budget `chord`, with its end
    /// tangents when it is one smooth conic (#231).
    fn directrix_sampled(
        &self,
        graph: &GeometryGraph,
        id: NodeId,
        range: Option<(Scalar, Scalar)>,
        options: &ExecutionOptions,
        chord: Scalar,
    ) -> GeomResult<axiolid_construct::sweep::SampledPath> {
        let options = options.clone().with_chord_error(chord).ok_or_else(|| {
            GeomError::InvalidInput(format!(
                "directrix chord budget {chord} is not positive and finite"
            ))
        })?;
        crate::directrix::sampled(graph, id, range, &options)
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
        regions: &Regions<'_>,
        node: &GeometryNode,
        key: EvalKey,
    ) -> GeomResult<Vec<EvalKey>> {
        let budget = key.budget();
        let same = |id| EvalKey::new(id, budget);
        // Placed by a transform, or in a station's frame (#264): a rigid
        // motion resolved from the graph. The source is meshed with the
        // initial room for the placement's rounding (#295); a retry with
        // more is scheduled by `evaluate`.
        if let Some((source, transform)) = placed_source(graph, node, key.id)? {
            return Ok(vec![EvalKey::new(
                source,
                instance_local_budget(transform, budget, 0.0)?,
            )]);
        }
        Ok(match node {
            GeometryNode::Collection(members) => members.iter().copied().map(same).collect(),
            GeometryNode::SolidOperation(SolidOperation::Boolean {
                left,
                right,
                operator,
            }) => {
                if let Some(clip) = clip_of(graph, regions, (*left, *right), *operator)? {
                    // The region is built from its subject in `build_boolean`.
                    vec![same(clip.subject)]
                } else if regions.has(*left) || regions.has(*right) {
                    // A region has no mesh of its own: `build_boolean` refuses
                    // the boolean by name (#326).
                    Vec::new()
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
        measure: Measure<'_>,
    ) -> GeomResult<Built> {
        let root_key = EvalKey::new(root, Budget::of(options));
        let mut stack = vec![Step::Enter(root_key)];
        // The room an instance's source was re-meshed with for its
        // placement's rounding (#295), once its first mesh rounded by more
        // than the initial room left.
        let mut rooms: std::collections::HashMap<EvalKey, Scalar> =
            std::collections::HashMap::new();
        // How far each node placed at a station lies from its exact
        // placement (#405), read once: the retry and the build see one.
        let mut accuracies: std::collections::HashMap<EvalKey, PlacementAccuracy> =
            std::collections::HashMap::new();
        while let Some(step) = stack.pop() {
            match step {
                Step::Enter(key) => {
                    if cache.contains_key(&key) {
                        continue;
                    }
                    let node = self.node(graph, key.id)?;
                    let deps = self.mesh_dependencies(
                        graph,
                        &Regions::new(graph, &cache.regions),
                        node,
                        key,
                    )?;
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
                    let budget = key.budget();
                    let node = self.node(graph, key.id)?;
                    if matches!(node, GeometryNode::InstanceAtStation(_)) {
                        if let std::collections::hash_map::Entry::Vacant(slot) =
                            accuracies.entry(key)
                        {
                            slot.insert(crate::station::placement_accuracy(
                                graph,
                                key.id,
                                crate::station::CERTIFICATION_SHARE * budget.chord,
                            )?);
                        }
                    }
                    let accuracy = accuracies.get(&key).copied();
                    if let std::collections::hash_map::Entry::Vacant(slot) = rooms.entry(key) {
                        if let Some((source, transform, room)) =
                            self.placement_retry(graph, node, key, cache, accuracy)?
                        {
                            slot.insert(room);
                            stack.push(Step::Exit(key));
                            stack.push(Step::Enter(EvalKey::new(
                                source,
                                instance_local_budget(transform, budget, room)?,
                            )));
                            continue;
                        }
                    }
                    let room = rooms.get(&key).copied().unwrap_or(0.0);
                    let local_options = options
                        .clone()
                        .with_tolerance(budget.tolerance)
                        .with_chord_error(budget.chord)
                        .ok_or_else(|| {
                            GeomError::InvalidInput(format!(
                                "chord budget {} is not positive and finite",
                                budget.chord
                            ))
                        })?;
                    let built = self.build(
                        graph,
                        (key, room, accuracy),
                        &local_options,
                        (cache, &Regions::new(graph, &cache.regions)),
                        (measure, root),
                    )?;
                    cache.insert(key, built);
                }
            }
        }
        cache
            .get(&root_key)
            .cloned()
            .ok_or_else(|| GeomError::InvalidInput(format!("root {root:?} produced no mesh")))
    }

    /// The room for its placement's rounding an instance's source needs
    /// beyond the initial one (#295), with the source and the transform:
    /// `Some` when the source, meshed with the initial room
    /// ([`placement_room`]), rounds by more once placed and its bound would
    /// carry the instance over the chord budget. The room is an eighth more
    /// than the rounding shown, kept below the budget: a mesh finer by
    /// that little has nearly the same coordinates. `None` for any other
    /// node, an unbounded source, or one that fits. A source whose placed
    /// vertices round by the whole budget is refused
    /// ([`placed_rounding`]).
    fn placement_retry(
        &self,
        graph: &GeometryGraph,
        node: &GeometryNode,
        key: EvalKey,
        cache: &Cache,
        accuracy: Option<PlacementAccuracy>,
    ) -> GeomResult<Option<(NodeId, Transform3, Scalar)>> {
        let Some((source_id, transform)) = placed_source(graph, node, key.id)? else {
            return Ok(None);
        };
        let budget = key.budget();
        let source = self.cached(
            cache,
            source_id,
            instance_local_budget(transform, budget, 0.0)?,
        )?;
        // What the placement adds to every bound: its rounding, and where
        // the station's frame has a certified error (#405), that.
        let rounding = placed_rounding(&source.mesh, transform, budget.chord)?;
        let shown = rounding + frame_term(accuracy, &source.mesh).unwrap_or(0.0);
        if shown <= placement_room(transform, budget.chord, 0.0)? || shown >= budget.chord {
            return Ok(None);
        }
        let Some(bound) = source.deviation.bound() else {
            return Ok(None);
        };
        let factor = crate::deviation::stretch(transform);
        if crate::deviation::grown(bound, factor, shown) <= budget.chord {
            return Ok(None);
        }
        let room = (shown * 1.125).min(shown + (budget.chord - shown) * 0.5);
        Ok(Some((source_id, transform, room)))
    }

    /// Build one node, assuming its mesh dependencies are already cached.
    /// `room` is the room an instance's source was re-meshed with for its
    /// placement's rounding (#295), `0` for the initial one, and `accuracy`
    /// how far a node placed at a station lies from its exact placement
    /// (#405).
    fn build(
        &self,
        graph: &GeometryGraph,
        (key, room, accuracy): (EvalKey, Scalar, Option<PlacementAccuracy>),
        options: &ExecutionOptions,
        (cache, regions): (&Cache, &Regions<'_>),
        (measure, root): (Measure<'_>, NodeId),
    ) -> GeomResult<Built> {
        let id = key.id;
        // The key's budget, placed scale included: what the dependencies
        // were evaluated at (`mesh_dependencies`).
        let budget = key.budget();
        let node = self.node(graph, id)?;
        match node {
            GeometryNode::TriMesh(mesh) => Ok(authored_mesh(mesh.clone(), options)),
            GeometryNode::PolygonMesh(mesh) => {
                compile_authored_polygons(mesh, options).map(|(mesh, warp, steiner_offset)| {
                    let built = authored_mesh(mesh, options);
                    let mut deviation = built.deviation.clone();
                    for (detail, bound) in [
                        (crate::deviation::WARPED_AUTHORED_FACE, warp),
                        (crate::deviation::KEYHOLE_STEINER_POINT, steiner_offset),
                    ] {
                        if let Some(bound) = bound {
                            deviation.add(
                                crate::deviation::DeviationPath::AuthoredMesh,
                                detail,
                                crate::deviation::DeviationBound::Certified(bound),
                            );
                        }
                    }
                    built.with_deviation(deviation)
                })
            }
            GeometryNode::Instance(instance) => {
                let source_budget = instance_local_budget(instance.transform, budget, room)?;
                let source = self.cached(cache, instance.source, source_budget)?;
                let rounding = placed_rounding(&source.mesh, instance.transform, budget.chord)?;
                Ok(channels::transform(source, instance.transform, rounding))
            }
            // An instance in a station's frame (#264): exact on a line;
            // elsewhere its frame's certified error grows every bound
            // (#405), or its placement is named unbounded.
            GeometryNode::InstanceAtStation(placed) => {
                let placement = crate::station::placement(graph, id)?;
                let source_budget = instance_local_budget(placement.transform, budget, room)?;
                let source = self.cached(cache, placed.source, source_budget)?;
                let rounding = placed_rounding(&source.mesh, placement.transform, budget.chord)?;
                let accuracy = accuracy.ok_or_else(|| {
                    GeomError::InvalidInput(format!(
                        "node {id:?} placed at a station was built before its accuracy was read"
                    ))
                })?;
                let unbounded = |reason| {
                    let mut built = channels::transform(source, placement.transform, rounding);
                    built.deviation.add(
                        crate::deviation::DeviationPath::StationPlacement,
                        crate::deviation::NUMERICAL_STATION_FRAME,
                        crate::deviation::DeviationBound::Unbounded(reason),
                    );
                    built
                };
                Ok(match accuracy {
                    PlacementAccuracy::Exact => {
                        channels::transform(source, placement.transform, rounding)
                    }
                    PlacementAccuracy::Certified(_) => {
                        match frame_term(Some(accuracy), &source.mesh) {
                            Some(frame) => channels::transform_within(
                                source,
                                placement.transform,
                                rounding,
                                frame,
                            ),
                            None => unbounded(UNCERTIFIED_STATION_FRAME),
                        }
                    }
                    PlacementAccuracy::Unbounded(reason) => unbounded(reason),
                })
            }
            GeometryNode::Collection(members) => {
                let members = members
                    .iter()
                    .map(|&member| self.cached(cache, member, budget))
                    .collect::<GeomResult<Vec<_>>>()?;
                Ok(channels::merge(&members))
            }
            GeometryNode::SolidOperation(SolidOperation::Boolean {
                left,
                right,
                operator,
            }) => {
                let (built, decided) = self.build_boolean(
                    graph,
                    (id, (*left, *right)),
                    *operator,
                    options,
                    budget,
                    (cache, regions, root),
                )?;
                Ok(match measure {
                    // The exact compiler's own result states its bound
                    // (#391): measuring it against itself would be circular.
                    Measure::Deviation(_) | Measure::Skip if decided == Decided::Referee => built,
                    // Measured against the exact result only when a report
                    // asks for it: it compiles the exact boolean (#235).
                    // Operands snapped within tolerance stay reported (#276).
                    Measure::Deviation(emitted) if emitted.contains(&id) => {
                        let mut deviation =
                            crate::deviation::of_boolean(graph, id, options, &built.mesh);
                        deviation.carry_perturbations(&built.deviation);
                        built.with_deviation(deviation)
                    }
                    Measure::Deviation(_) | Measure::Skip => built,
                })
            }
            GeometryNode::SolidOperation(operation) => {
                let (mesh, deviation) = self.build_solid_within(graph, operation, options)?;
                Ok(Built::leaf(mesh).with_deviation(deviation))
            }
            GeometryNode::BRep(brep) => {
                crate::brep::tessellate(brep, graph, options.tolerance(), chord_error(options)).map(
                    |(mesh, closure, deviation)| {
                        Built::with_closure(mesh, closure).with_deviation(deviation)
                    },
                )
            }
            // A plane trimmed by boundary curves: a surface (#192).
            GeometryNode::SurfaceRelation(axiolid_model::SurfaceRelation::CurveBounded {
                basis,
                boundaries,
                implicit_outer,
            }) => crate::bounded::curve_bounded(
                self.descriptor().id,
                graph,
                *basis,
                boundaries,
                *implicit_outer,
                options,
            )
            .map(|(mesh, merged)| {
                Built::with_closure(mesh, axiolid_mesh_compile_contract::MeshClosure::Surface)
                    .with_deviation(crate::deviation::of_curve_bounded(
                        graph, *basis, boundaries, options, merged,
                    ))
            }),
            // Open sections at stations along a directrix: a sheet (#241).
            GeometryNode::SurfaceRelation(axiolid_model::SurfaceRelation::SectionedSurface {
                directrix,
                sections,
                frame,
            }) => {
                let sections: Vec<_> = sections.iter().cloned().map(Into::into).collect();
                crate::station::sectioned_surface(graph, *directrix, &sections, *frame, options)
                    .map(|(mesh, bound)| {
                        Built::with_closure(
                            mesh,
                            axiolid_mesh_compile_contract::MeshClosure::Surface,
                        )
                        .with_deviation(crate::deviation::of_sectioned_surface(bound))
                    })
            }
            // The general form: oriented sections, joined by tag (#246).
            GeometryNode::SurfaceRelation(
                axiolid_model::SurfaceRelation::OpenSectionsAtStations {
                    directrix,
                    sections,
                    frame,
                },
            ) => crate::station::sectioned_surface(graph, *directrix, sections, *frame, options)
                .map(|(mesh, bound)| {
                    Built::with_closure(mesh, axiolid_mesh_compile_contract::MeshClosure::Surface)
                        .with_deviation(crate::deviation::of_sectioned_surface(bound))
                }),
            // A surface bounded by a parameter rectangle: a sheet (#311).
            GeometryNode::SurfaceRelation(axiolid_model::SurfaceRelation::RectangularTrimmed {
                basis,
                u,
                v,
                u_sense,
                v_sense,
            }) => crate::trimmed::rectangular(
                graph,
                *basis,
                crate::trimmed::Trim {
                    u: *u,
                    v: *v,
                    u_sense: *u_sense,
                    v_sense: *v_sense,
                },
                options,
            )
            .map(|(mesh, family, bound)| {
                Built::with_closure(mesh, axiolid_mesh_compile_contract::MeshClosure::Surface)
                    .with_deviation(crate::deviation::Deviation::one(
                        crate::deviation::DeviationPath::TrimmedSurface,
                        family,
                        bound,
                    ))
            }),
            // A surface ruled between two sides: a sheet (#402).
            GeometryNode::SurfaceRelation(axiolid_model::SurfaceRelation::Ruled {
                first,
                second,
            }) => crate::ruled::ruled(graph, first, second, options).map(|(mesh, bound)| {
                Built::with_closure(mesh, axiolid_mesh_compile_contract::MeshClosure::Surface)
                    .with_deviation(crate::deviation::Deviation::one(
                        crate::deviation::DeviationPath::RuledSurface,
                        "ruled",
                        bound,
                    ))
            }),
            // CSG primitives are analytic solids: no surface evaluation,
            // no trim curves, just a closed mesh at the caller's tolerance.
            GeometryNode::Primitive(primitive) => {
                axiolid_reference::primitive::tessellate_primitive(
                    primitive,
                    chord_tolerance(options)?,
                )
                .map(|mesh| {
                    let deviation = crate::deviation::of_primitive(primitive, chord_error(options));
                    Built::leaf(mesh).with_deviation(deviation)
                })
            }
            // Unbounded on its kept side, a half-space has a finite stand-in
            // only against a subject (#326): it is a clip's tool.
            GeometryNode::HalfSpace(_) => Err(unsupported_input(HALF_SPACE_WITHOUT_SUBJECT)),
            other => Err(GeomError::Unsupported {
                backend: self.descriptor().id,
                operation: unsupported_operation(other),
            }),
        }
    }

    /// Read an already-built dependency.
    fn cached<'c>(&self, cache: &'c Cache, id: NodeId, budget: Budget) -> GeomResult<&'c Built> {
        cache.get(&EvalKey::new(id, budget)).ok_or_else(|| {
            GeomError::InvalidInput(format!(
                "dependency {id:?} was not evaluated first at {budget:?}"
            ))
        })
    }

    /// A boolean, with the operands' channel fates composed onto the
    /// provider's (#115): a channel lost upstream stays reported lost.
    /// `budget` is the boolean's own, the placed scale included (#312).
    ///
    /// Operands placed by the same linear part are cut in the frame they
    /// share and the result is placed (#356,
    /// [`Self::boolean_in_shared_frame`]): placing each operand rounds the
    /// vertices of faces that lie on one plane in exact arithmetic onto
    /// planes a residue apart, which the snap, landing only on axis-aligned
    /// planes, cannot close under a general rotation. That cut is never
    /// worse: where it is refused, for any reason (the unmerged-contact
    /// guard included), the operands are cut in the world frame as before.
    ///
    /// A region (a half-space, or an intersection of half-spaces, read
    /// through instances) is a clip's tool (#326): of a difference or an
    /// intersection, or, an intersection commuting, on the left of one with
    /// a solid. It has no mesh of its own and is built here from its
    /// subject ([`Self::region_tool_mesh`]); any other boolean with one is
    /// refused by name.
    ///
    /// An operand with no triangles (a clip that left nothing, an
    /// intersection of disjoint solids) decides the boolean by an identity
    /// before any frame or cut (#326, [`crate::empty`]): the subject, the
    /// tool or nothing, as built, with no deviation of a cut. A skin the
    /// subject has pending passes through to the root of its chain.
    ///
    /// A boolean the world-frame mesh cut cannot settle is decided by the
    /// exact compiler (#391, #423, [`crate::referee`]): when faces of its operands
    /// nearly share a plane within the linear tolerance, or when the cut is
    /// refused as a pinch or as unmerged contact. The referee is never
    /// worse than the mesh path: where it refuses, a valid mesh result is
    /// returned as before, and a contact refusal keeps its message and names
    /// the referee's reason.
    fn build_boolean(
        &self,
        graph: &GeometryGraph,
        (id, (left, right)): (NodeId, (NodeId, NodeId)),
        operator: axiolid_core::BooleanOperator,
        options: &ExecutionOptions,
        budget: Budget,
        (cache, regions, root): (&Cache, &Regions<'_>, NodeId),
    ) -> GeomResult<(Built, Decided)> {
        // A region is the tool of a clip, or the left operand of an
        // intersection with a solid, whose operands are then read the other
        // way round; any other boolean with one is unbounded (#326).
        let clip = clip_of(graph, regions, (left, right), operator)?;
        if clip.is_none() && (regions.has(left) || regions.has(right)) {
            return Err(unsupported_input(UNCLIPPED_HALF_SPACE));
        }
        let (left, right) = clip
            .as_ref()
            .map_or((left, right), |clip| (clip.subject, clip.tool));
        let subject = self.cached(cache, left, budget)?;
        refuse_surface_operand(subject, "subject")?;
        // Any other tool is a node of its own, built before this one.
        let built_tool = match &clip {
            Some(_) => None,
            None => {
                let built = self.cached(cache, right, budget)?;
                refuse_surface_operand(built, "tool")?;
                Some(built)
            }
        };
        // An empty operand decides the boolean by an identity, before any
        // frame or cut (#326, `crate::empty`).
        let identity = crate::empty::identity(
            operator,
            subject.mesh.indices.is_empty(),
            built_tool.is_some_and(|built| built.mesh.indices.is_empty()),
        );
        if identity.is_none() {
            if let Ok(Some(built)) = self.boolean_in_shared_frame(
                graph,
                (left, right),
                operator,
                options,
                budget,
                (cache, regions),
            ) {
                return Ok((built, Decided::Mesh));
            }
            // A clip is a constructed stand-in with no coordinates of its
            // own: cut in its subject's frame, never worse than the world cut
            // below (#400, ADR 0102).
            if let Ok(Some(built)) = self.clip_in_subject_frame(
                graph,
                (left, right),
                operator,
                options,
                budget,
                (cache, regions),
            ) {
                return Ok((built, Decided::Mesh));
            }
        }
        let region_tool = match &clip {
            Some(clip) if identity.is_none() => {
                self.region_tool_mesh(graph, &clip.leaves, &subject.mesh, options, budget.placed)?
            }
            _ => None,
        };
        // A region that leaves nothing of the subject's reach is an empty
        // tool (#326).
        let identity = match (&clip, &region_tool, identity) {
            (Some(_), None, None) => crate::empty::identity(operator, false, true),
            (.., identity) => identity,
        };
        // A half-space is built here from the subject, so it has no
        // upstream fates of its own; the deviation of a tool cut from
        // several leaves is carried (#326).
        let nothing = TriMesh::default();
        let (tool, tool_fates, tool_deviation) = match (&region_tool, built_tool, &clip) {
            (Some(built), _, Some(clip)) => (
                &built.mesh,
                None,
                (clip.leaves.len() > 1).then_some(&built.deviation),
            ),
            (None, Some(built), _) => (&built.mesh, Some(&built.fates), Some(&built.deviation)),
            _ => (&nothing, None, None),
        };
        let cut = match identity {
            Some(identity) => {
                identity_cut(identity, subject, built_tool).map_err(CutRefusal::Other)
            }
            None => self.boolean_of(
                subject,
                (tool, tool_fates, tool_deviation),
                operator,
                options,
                budget.placed,
            ),
        };
        // A boolean the mesh cut cannot settle is decided by the exact
        // compiler (#391, #423, ADR 0098). It is never worse than the mesh path:
        // where the referee refuses, an unsettled but valid mesh result is
        // returned as before, and a contact refusal keeps its message and
        // names the referee's reason.
        let half_space_tool = clip.is_some();
        let referee = |reason: &mut Option<String>,
                       against: Option<crate::referee::Against<'_>>,
                       rounding: Option<Scalar>| {
            // A half-space is clipped to the subject here and never built as
            // a node of its own, so it has no built operand to hand the
            // referee. The extra-shell trigger of #423 leaves it to the mesh
            // result; a skin within rounding is decided for it too (#400,
            // ADR 0102).
            if half_space_tool && rounding.is_none() {
                *reason = Some("the tool is a half-space clipped to the subject".to_owned());
                return Ok(None);
            }
            match crate::referee::boolean(
                graph,
                id,
                operator,
                options,
                (subject, built_tool),
                against,
                rounding,
            ) {
                Ok(built) => Ok::<_, GeomError>(Some(built)),
                Err(why) => {
                    *reason = Some(why);
                    Ok(None)
                }
            }
        };
        let mut reason = None;
        match cut {
            // A skin within rounding is decided by the exact referee, once,
            // at the root of the chain of clips it sits in: a chain defers
            // it, and one exact compile of the whole chain decides it
            // (#400, ADR 0102).
            Ok(MeshCut {
                mut built,
                skin: Some(skin),
                ..
            }) => {
                if clip_chain_continues(graph, regions, id, root) {
                    built.pending_skin = Some(skin);
                    return Ok((built, Decided::Mesh));
                }
                match referee(&mut reason, None, Some(skin.reach))? {
                    Some(decided) => Ok((decided, Decided::Referee)),
                    None => Err(refused_by_referee(skin.error, reason)),
                }
            }
            Ok(cut) if !cut.unsettled => Ok((cut.built, Decided::Mesh)),
            Ok(cut) => Ok(
                match referee(
                    &mut reason,
                    (operator != axiolid_core::BooleanOperator::Union).then(|| {
                        crate::referee::Against {
                            mesh: &cut.built.mesh,
                            reach: crate::snap::reach(
                                &subject.mesh,
                                tool,
                                options.tolerance(),
                                budget.placed,
                            ),
                        }
                    }),
                    None,
                )? {
                    Some(built) => (built, Decided::Referee),
                    None => (cut.built, Decided::Mesh),
                },
            ),
            // A pinch against the mesh a deferred skin left (#400) is decided
            // by the referee, at the rounding reach.
            Err(CutRefusal::Contact(error))
                if operator != axiolid_core::BooleanOperator::Union
                    && subject.pending_skin.is_some() =>
            {
                let reach =
                    crate::snap::reach(&subject.mesh, tool, options.tolerance(), budget.placed);
                match referee(&mut reason, None, Some(reach))? {
                    Some(decided) => Ok((decided, Decided::Referee)),
                    None => Err(refused_by_referee(error, reason)),
                }
            }
            // A pinched difference or an intersection has no mesh result to
            // hold a void: its refusal stands (ADR 0096, D7).
            Err(CutRefusal::Contact(error)) if operator != axiolid_core::BooleanOperator::Union => {
                Err(error)
            }
            Err(CutRefusal::Contact(error)) => match referee(&mut reason, None, None)? {
                Some(built) => Ok((built, Decided::Referee)),
                None => Err(refused_by_referee(error, reason)),
            },
            Err(CutRefusal::Other(error)) => Err(error),
        }
    }

    /// The operand of a boolean as a frame sees it (#356): an `Instance`'s
    /// transform with its source as cached for it, or the frame a boolean
    /// was itself cut in with its unplaced result. `None` for any other
    /// node, and for an instance whose source is re-meshed for its
    /// placement's rounding ([`Self::placement_retry`]): the world frame
    /// handles that.
    fn frame_of<'c>(
        &self,
        graph: &GeometryGraph,
        id: NodeId,
        budget: Budget,
        cache: &'c Cache,
    ) -> GeomResult<Option<(Transform3, &'c Built)>> {
        let node = self.node(graph, id)?;
        if let GeometryNode::Instance(instance) = node {
            if self
                .placement_retry(graph, node, EvalKey::new(id, budget), cache, None)?
                .is_some()
            {
                return Ok(None);
            }
            let source_budget = instance_local_budget(instance.transform, budget, 0.0)?;
            let source = self.cached(cache, instance.source, source_budget)?;
            return Ok(Some((instance.transform, source)));
        }
        Ok(self
            .cached(cache, id, budget)?
            .local
            .as_deref()
            .map(|local| (local.transform, &local.built)))
    }

    /// The boolean cut in the frame its operands share and placed by the
    /// subject's transform (#356): `None` when they share none, that is
    /// when either is not framed ([`Self::frame_of`]) or their linear parts
    /// differ in a bit ([`crate::frame::shared_linear`]; a union is then
    /// refereed, [`crate::referee`]), or
    /// when the tool is a half-space clip, or when the instance chains
    /// beneath the frames compose to linear parts that differ in a bit
    /// (the sources were themselves placed by differently composed turns):
    /// the world path re-reads it and the referee decides. The tool is moved by its offset
    /// in the subject's frame, which is exactly zero when the translations
    /// are equal.
    ///
    /// The result keeps its unplaced form in [`Built::local`], so an
    /// enclosing boolean in the same frame cuts against it: a chain of
    /// corbels on a column stays in the column's frame.
    ///
    /// Its operands are the caller's, already refused if they are surface
    /// models and not empty (#326).
    fn boolean_in_shared_frame(
        &self,
        graph: &GeometryGraph,
        (left, right): (NodeId, NodeId),
        operator: axiolid_core::BooleanOperator,
        options: &ExecutionOptions,
        budget: Budget,
        (cache, regions): (&Cache, &Regions<'_>),
    ) -> GeomResult<Option<Built>> {
        if clip_of(graph, regions, (left, right), operator)?.is_some() {
            return Ok(None);
        }
        let Some((subject_frame, subject)) = self.frame_of(graph, left, budget, cache)? else {
            return Ok(None);
        };
        let Some((tool_frame, tool)) = self.frame_of(graph, right, budget, cache)? else {
            return Ok(None);
        };
        if !crate::frame::shared_linear(subject_frame, tool_frame) {
            return Ok(None);
        }
        // The outer frames agree, but a source placed beneath them by turns
        // computed differently leaves faces that nearly share a plane in
        // the frame as well: the world path re-reads and referees it (#391).
        if let (Some(a), Some(b)) = (
            chain_linear(graph, left, subject_frame),
            chain_linear(graph, right, tool_frame),
        ) {
            if !crate::frame::shared_linear(a, b) {
                return Ok(None);
            }
        }
        let offset = crate::frame::relative_offset(subject_frame, tool_frame).ok_or_else(|| {
            GeomError::InvalidInput("operand translations have no finite offset".into())
        })?;
        let shifted;
        let tool = if offset == Vec3::ZERO {
            tool
        } else {
            let shift = Transform3::from_translation(offset);
            let rounding = crate::deviation::placement_rounding(&tool.mesh.positions, shift);
            shifted = channels::transform(tool, shift, rounding);
            &shifted
        };
        let local_budget = instance_local_budget(subject_frame, budget, 0.0)?;
        let tool_budget = instance_local_budget(tool_frame, budget, 0.0)?;
        let local_options = local_options(options, local_budget)?;
        let local = self
            .boolean_of(
                subject,
                (&tool.mesh, Some(&tool.fates), Some(&tool.deviation)),
                operator,
                &local_options,
                local_budget.placed.max(tool_budget.placed),
            )
            .map_err(CutRefusal::into_error)?;
        // A skin within rounding is the world cut's to decide (#400).
        if let Some(skin) = local.skin {
            return Err(skin.error);
        }
        let local = local.built;
        let rounding = placed_rounding(&local.mesh, subject_frame, budget.chord)?;
        let mut placed = channels::transform(&local, subject_frame, rounding);
        // Defence in depth: the local tolerance is the world one over the
        // frame's stretch, so a union that still leaves bodies within the
        // world tolerance once placed is refused here, and `build_boolean`
        // cuts it in the world frame instead.
        if operator == axiolid_core::BooleanOperator::Union {
            let eps = options.tolerance().linear();
            if let Some(found) = crate::contact::unmerged(&placed.mesh, eps) {
                return Err(crate::contact::refusal(&found, eps));
            }
        }
        placed.local = Some(Box::new(channels::LocalBoolean {
            transform: subject_frame,
            built: local,
        }));
        Ok(Some(placed))
    }

    /// A difference or intersection with a half-space clip, cut in the frame
    /// its subject is placed by (#400): `None` when the operator does not
    /// clip, the tool is no clip, the subject is not framed
    /// ([`Self::frame_of`]) or not a similarity, or the tool's placement has
    /// no finite form in that frame ([`crate::frame::relative`]).
    ///
    /// Unlike [`Self::boolean_in_shared_frame`] the linear parts need not be
    /// bit-equal (ADR 0096, D3): a clip has no coordinates of its own to
    /// keep. The tool is rebuilt against the subject's local mesh, where the
    /// subject's faces are exactly axis-aligned; where the linear parts
    /// differ, its faces lie a residue off them and the snap lands them. The
    /// result is placed by the frame and keeps its unplaced form in
    /// [`Built::local`]. A region (#326) is placed leaf by leaf, each
    /// relative to the frame; one with a leaf that has no finite form there,
    /// or that leaves nothing of the subject, is the world cut's. The subject
    /// is the caller's, already refused if it is a surface model and not
    /// empty.
    fn clip_in_subject_frame(
        &self,
        graph: &GeometryGraph,
        (left, right): (NodeId, NodeId),
        operator: axiolid_core::BooleanOperator,
        options: &ExecutionOptions,
        budget: Budget,
        (cache, regions): (&Cache, &Regions<'_>),
    ) -> GeomResult<Option<Built>> {
        let Some(clip) = clip_of(graph, regions, (left, right), operator)? else {
            return Ok(None);
        };
        let Some((frame, subject)) = self.frame_of(graph, left, budget, cache)? else {
            return Ok(None);
        };
        // Each leaf of a region (#326) is placed in the subject's frame; one
        // that has no finite form there leaves the world cut to decide.
        let mut leaves = Vec::with_capacity(clip.leaves.len());
        for &(leaf, chain) in &clip.leaves {
            let Some(relative) = crate::frame::relative(frame, chain) else {
                return Ok(None);
            };
            leaves.push((leaf, relative));
        }
        let local_budget = instance_local_budget(frame, budget, 0.0)?;
        let local_options = local_options(options, local_budget)?;
        // A region that leaves nothing of the subject is an empty tool: the
        // world cut decides it by an identity.
        let Some(tool) = self.region_tool_mesh(
            graph,
            &leaves,
            &subject.mesh,
            &local_options,
            local_budget.placed,
        )?
        else {
            return Ok(None);
        };
        let local = self
            .boolean_of(
                subject,
                (
                    &tool.mesh,
                    None,
                    (leaves.len() > 1).then_some(&tool.deviation),
                ),
                operator,
                &local_options,
                local_budget.placed,
            )
            .map_err(CutRefusal::into_error)?;
        if let Some(skin) = local.skin {
            return Err(skin.error);
        }
        let local = local.built;
        let rounding = placed_rounding(&local.mesh, frame, budget.chord)?;
        let mut placed = channels::transform(&local, frame, rounding);
        placed.local = Some(Box::new(channels::LocalBoolean {
            transform: frame,
            built: local,
        }));
        Ok(Some(placed))
    }

    /// The cut of two operands in one frame: settled onto each other
    /// ([`crate::snap`]), cut, and reported. `placed` is the scale the
    /// result is placed at (#312).
    fn boolean_of(
        &self,
        subject: &Built,
        (tool, tool_fates, tool_deviation): (
            &TriMesh,
            Option<&channels::Fates>,
            Option<&crate::deviation::Deviation>,
        ),
        operator: axiolid_core::BooleanOperator,
        options: &ExecutionOptions,
        placed: Scalar,
    ) -> Result<MeshCut, CutRefusal> {
        // Faces of the two operands a rounding residue apart are moved
        // together first, so a skin or sliver that thin is not cut as one
        // (#276, narrowed to rounding in #291; `crate::snap`): the rounding
        // of the operands' coordinates, or of the coordinates the result is
        // placed at, whichever is coarser (#312).
        // A union reads contact at the tolerance (#357), the others at the
        // rounding.
        let reach = if operator == axiolid_core::BooleanOperator::Union {
            crate::snap::contact_reach(options.tolerance())
        } else {
            crate::snap::reach(&subject.mesh, tool, options.tolerance(), placed)
        };
        let settled = crate::snap::settle(
            &subject.mesh,
            tool,
            reach,
            operator != axiolid_core::BooleanOperator::Union,
        );
        let snapped = if settled.moved > 0.0 {
            // A snap never turns a working boolean into a refusal (#291):
            // where the snapped operands are refused, or their result touches
            // itself or holds a skin, they are cut as given below, and no
            // snap is reported.
            self.cut(
                settled.subject.as_ref().unwrap_or(&subject.mesh),
                settled.tool.as_ref().unwrap_or(tool),
                operator,
                options,
                reach,
            )
            .ok()
            .filter(|(_, skin)| skin.is_none())
        } else {
            None
        };
        let ((outcome, skin), moved) = match snapped {
            Some(cut) => (cut, settled.moved),
            None => (
                self.cut(&subject.mesh, tool, operator, options, reach)?,
                0.0,
            ),
        };
        // Faces of the pair actually cut that nearly share a plane, a
        // residue off exact contact, are not decided by a boolean without a
        // tolerance whatever its result looks like (#391): the exact
        // compiler referees them.
        let snapped_pair = moved > 0.0;
        let cut_subject = settled
            .subject
            .as_ref()
            .filter(|_| snapped_pair)
            .unwrap_or(&subject.mesh);
        let cut_tool = settled
            .tool
            .as_ref()
            .filter(|_| snapped_pair)
            .unwrap_or(tool);
        let eps = options.tolerance().linear();
        // Near contact calls the referee. A union is decided by the exact
        // compiler then; a difference or an intersection only where its mesh
        // holds a shell the exact result lacks (#423), as a skin thinner
        // than the tolerance stays as authored (ADR 0096, D7, #291).
        // The cheap conditions come first: the exact compile of the whole
        // subtree is dear, and a mesh of one component holds no void.
        let unsettled = eps > 0.0
            && (operator == axiolid_core::BooleanOperator::Union
                || crate::deviation::geometric_components(&outcome.mesh) > 1)
            && axiolid_inspect::near_coplanar(cut_subject, cut_tool, eps).is_some();
        let mut built = channels::after_boolean(
            outcome.mesh,
            &subject.fates,
            tool_fates,
            outcome.evidence.attribute_fates,
        );
        built.deviation.carry_perturbations(&subject.deviation);
        if let Some(tool) = tool_deviation {
            built.deviation.carry_perturbations(tool);
        }
        if moved > 0.0 {
            built.deviation.add(
                crate::deviation::DeviationPath::Boolean,
                crate::deviation::SNAPPED_OPERANDS,
                crate::deviation::DeviationBound::Certified(moved),
            );
        }
        Ok(MeshCut {
            built,
            unsettled,
            skin,
        })
    }

    /// The mesh boolean of `subject` and `tool`, refused where its result
    /// touches itself (#194), and, for a union, where it holds two bodies
    /// within the linear tolerance of each other (#356, #357): the operands
    /// meet there and the boolean, which has no tolerance, left them apart.
    /// That union is refused by name, [`crate::UNMERGED_UNION_CONTACT`]. A
    /// difference or an intersection is refused where its result holds a
    /// skin thinner than the rounding `reach`, [`crate::SKIN_WITHIN_ROUNDING`]
    /// (#400).
    fn cut(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        operator: axiolid_core::BooleanOperator,
        options: &ExecutionOptions,
        reach: Scalar,
    ) -> Result<
        (
            axiolid_mesh_boolean_contract::BooleanOutcome,
            Option<PendingSkin>,
        ),
        CutRefusal,
    > {
        let outcome = self
            .boolean
            .boolean(subject, tool, operator, options)
            .map_err(CutRefusal::Other)?;
        if let Some(pinch) = crate::pinch::find(&outcome.mesh) {
            return Err(CutRefusal::Contact(GeomError::Degenerate(match pinch {
                crate::pinch::Pinch::Edge { from, to } => format!(
                    "{operator:?} result touches itself along the edge {from:?} to {to:?}: \
                     the operands meet tangentially there (a void tangent to its host's face), \
                     leaving no material between two faces"
                ),
                crate::pinch::Pinch::Vertex { at } => format!(
                    "{operator:?} result touches itself at the point {at:?}: the operands \
                     meet tangentially there, leaving no material between two faces"
                ),
            })));
        }
        if operator == axiolid_core::BooleanOperator::Union {
            let eps = options.tolerance().linear();
            if let Some(found) = crate::contact::unmerged(&outcome.mesh, eps) {
                return Err(CutRefusal::Contact(crate::contact::refusal(&found, eps)));
            }
        }
        // A difference or an intersection leaves no skin its coordinates
        // cannot carry (#400): two faces facing each other within the
        // rounding are residue, and a Solid with them is not decided here.
        let skin = if matches!(
            operator,
            axiolid_core::BooleanOperator::Difference | axiolid_core::BooleanOperator::Intersection
        ) {
            axiolid_inspect::thin_sheet(&outcome.mesh, reach).map(|found| PendingSkin {
                error: crate::skin::refusal(&found, &outcome.mesh, reach, operator),
                reach,
            })
        } else {
            None
        };
        Ok((outcome, skin))
    }

    /// A non-boolean solid and its deviation. A swept disk along a smooth
    /// directrix places its stations by a second-order estimate, so its
    /// certified bound can land a little over the budget: it is rebuilt
    /// with the construction's chord shrunk by the overshoot, and
    /// certified again against the requested budget, at most
    /// [`REFINE_ROUNDS`] times. A bound still over it is reported as it
    /// is, never trimmed (#252).
    fn build_solid_within(
        &self,
        graph: &GeometryGraph,
        operation: &SolidOperation,
        options: &ExecutionOptions,
    ) -> GeomResult<(TriMesh, crate::deviation::Deviation)> {
        // A run of sections at stations is bounded as it is sampled (#247).
        if let Some(run) = crate::deviation::stationed_spine(graph, operation, options) {
            return run;
        }
        let mut mesh = self.build_solid(graph, operation, options)?;
        let mut deviation = crate::deviation::of_solid(graph, operation, options, options, &mesh)?;
        if !matches!(operation, SolidOperation::SweptDisk { .. }) {
            return Ok((mesh, deviation));
        }
        let budget = chord_error(options);
        let mut chord = budget;
        for _ in 0..REFINE_ROUNDS {
            let Some(over) = deviation.certified_over(budget) else {
                break;
            };
            chord *= REFINE_MARGIN * budget / over;
            let Some(finer) = options.clone().with_chord_error(chord) else {
                break;
            };
            let refined = self.build_solid(graph, operation, &finer)?;
            let certified =
                crate::deviation::of_solid(graph, operation, options, &finer, &refined)?;
            match certified.certified_over(budget) {
                // No better: the overshoot is not the stations' to remove.
                Some(next) if next >= over => break,
                _ => (mesh, deviation) = (refined, certified),
            }
        }
        Ok((mesh, deviation))
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
                // Doubly curved (#231): half the budget across the profile,
                // half round the axis, so the two chordings add to at most
                // the whole (see `axiolid_construct::revolve`).
                let half = 0.5 * chord_error(options);
                let rings = profile_rings(shape, half, options.tolerance())?;
                axiolid_construct::revolve::revolve(
                    &rings,
                    *axis_origin,
                    *axis_direction,
                    *angle,
                    half_tolerance(options)?,
                )
            }
            SolidOperation::TaperedExtrusion {
                start_profile,
                end_profile,
                direction,
                depth,
            } => {
                let chord = chord_error(options);
                let a =
                    self.rings_of(graph, *start_profile, chord, options, "taper start profile")?;
                let b = self.rings_of(graph, *end_profile, chord, options, "taper end profile")?;
                axiolid_construct::sweep::tapered_extrude(&a, &b, *direction, *depth)
            }
            SolidOperation::TaperedRevolution {
                start_profile,
                end_profile,
                axis_origin,
                axis_direction,
                angle,
            } => {
                let half = 0.5 * chord_error(options);
                let a =
                    self.rings_of(graph, *start_profile, half, options, "taper start profile")?;
                let b = self.rings_of(graph, *end_profile, half, options, "taper end profile")?;
                axiolid_construct::sweep::tapered_revolve(
                    &a,
                    &b,
                    *axis_origin,
                    *axis_direction,
                    *angle,
                    half_tolerance(options)?,
                )
            }
            // Sweeps along a directrix are doubly curved too: the
            // construction splits the budget between the section and the
            // directrix and refines the directrix until its walls fit.
            // Lines, polylines and composites of lines and arcs are swept
            // piece by piece, corners mitred (#245, beside an arc #248) or
            // filleted (#232); a closed polyline is swept as a loop (#248).
            SolidOperation::SweptDisk {
                directrix,
                radius,
                inner_radius,
                parameter_range,
                fillet_radius,
            } => match crate::directrix::pieces(graph, *directrix, *parameter_range, options)? {
                Some(_) if fillet_radius.is_some() && parameter_range.is_some() => {
                    Err(GeomError::UnsupportedInput {
                        backend: self.descriptor().id,
                        operation: axiolid_contracts::Operation::Sweep,
                        input: "a parameter range on a swept disk with a fillet radius",
                    })
                }
                Some(path) if path.closed => axiolid_construct::pipe::swept_disk_along_loop(
                    &path.pieces,
                    *radius,
                    *inner_radius,
                    *fillet_radius,
                    chord_error(options),
                ),
                Some(path) => axiolid_construct::pipe::swept_disk_along_pieces(
                    &path.pieces,
                    *radius,
                    *inner_radius,
                    *fillet_radius,
                    chord_error(options),
                ),
                None => {
                    // The conic path checks only the radii; a disk wider
                    // than a lone arc's bend would mesh into a tube that
                    // folds through the bend's axis (#122).
                    if let Some(path) = crate::directrix::exact_pieces(
                        graph,
                        *directrix,
                        *parameter_range,
                        options,
                    )? {
                        if let [arc @ axiolid_construct::pipe::PathPiece::Arc { .. }] =
                            path.pieces.as_slice()
                        {
                            axiolid_construct::pipe::check_disk_along(arc, *radius)?;
                        }
                    }
                    axiolid_construct::sweep::swept_disk_within(
                        |chord| {
                            self.directrix_sampled(
                                graph,
                                *directrix,
                                *parameter_range,
                                options,
                                chord,
                            )
                        },
                        *radius,
                        *inner_radius,
                        *fillet_radius,
                        chord_error(options),
                    )
                }
            },
            SolidOperation::FixedReferenceSweep {
                profile,
                directrix,
                reference_direction,
                parameter_range,
            } => {
                let half = 0.5 * chord_error(options);
                let rings = self.rings_of(graph, *profile, half, options, "sweep profile")?;
                // A tangent-continuous chain of segments and arcs is swept
                // piece by piece with a proven bound, a corner refused (#122);
                // an arc oblique to the reference, or any other family, is
                // swept as sampled.
                if let Some(path) = crate::directrix::fixed_reference_pieces(
                    graph,
                    *directrix,
                    *parameter_range,
                    options,
                )? {
                    let plan = axiolid_construct::pipe::fixed_reference_plan(
                        &path.pieces,
                        path.closed,
                        *reference_direction,
                        axiolid_construct::pipe::ring_reach(&rings),
                        chord_error(options),
                        options.tolerance(),
                    )?;
                    if let axiolid_construct::pipe::FixedReferencePlan::Untwisted(chain) = plan {
                        return axiolid_construct::pipe::fixed_reference_along_pieces(
                            &rings, &chain,
                        );
                    }
                }
                axiolid_construct::sweep::fixed_reference_sweep_within(
                    &rings,
                    |chord| {
                        self.directrix_sampled(graph, *directrix, *parameter_range, options, chord)
                    },
                    *reference_direction,
                    half,
                )
            }
            SolidOperation::SurfaceCurveSweep {
                profile,
                directrix,
                reference_surface,
                parameter_range,
            } => {
                let half = 0.5 * chord_error(options);
                let rings = self.rings_of(
                    graph,
                    *profile,
                    half,
                    options,
                    "surface curve sweep profile",
                )?;
                axiolid_construct::sweep::surface_curve_sweep_within(
                    &rings,
                    |chord| {
                        let path = self.directrix_sampled(
                            graph,
                            *directrix,
                            *parameter_range,
                            options,
                            chord,
                        )?;
                        let normals =
                            self.surface_normals(graph, *reference_surface, &path.points, options)?;
                        Ok((path, normals))
                    },
                    half,
                )
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
                    let rings = self.rings_of(
                        graph,
                        section.profile,
                        chord_error(options),
                        options,
                        "spine section profile",
                    )?;
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
            // Sections at stations along the directrix (#241).
            SolidOperation::StationedSpine {
                directrix,
                sections,
                frame,
            } => {
                let sections: Vec<_> = sections.iter().copied().map(Into::into).collect();
                crate::station::stationed_spine(graph, *directrix, &sections, *frame, options)
            }
            // The general form: oriented sections, matched by tag (#246).
            SolidOperation::SectionsAtStations {
                directrix,
                sections,
                frame,
            } => crate::station::stationed_spine(graph, *directrix, sections, *frame, options),
            // Unbounded on its kept side, a bounded half-space has a finite
            // stand-in only against a subject (#325): a slab sized by its
            // boundary alone ends wherever it ends, and a subject reaching
            // further would be cut too little, silently. Its booleans are
            // routed to `build_boolean`, which sizes it; here it has none.
            SolidOperation::BoundedHalfSpace { .. }
            | SolidOperation::BoundedHalfSpaceAlongZ { .. } => Err(GeomError::UnsupportedInput {
                backend: self.descriptor().id,
                operation: Operation::GraphCompilation,
                input: "bounded half-space without a subject to size it by: it is unbounded \
                        on its kept side, and is compiled only as the tool of a difference \
                        or an intersection",
            }),
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

/// The cut of a boolean with an empty operand (#326): the subject or the
/// tool as it is, or nothing, by the identity table. Nothing is cut, so no
/// snap, skin or contact is read; a skin the subject carries stays pending,
/// for the root of its chain to decide.
fn identity_cut(
    identity: crate::empty::Identity,
    subject: &Built,
    tool: Option<&Built>,
) -> GeomResult<MeshCut> {
    use crate::empty::Identity;
    let built = match identity {
        Identity::Subject => {
            let mut built = subject.clone();
            let skin = built.pending_skin.take();
            return Ok(MeshCut {
                built,
                unsettled: false,
                skin,
            });
        }
        Identity::Tool => tool
            .ok_or_else(|| {
                GeomError::InvalidInput("boolean identity names a tool that was not built".into())
            })?
            .clone(),
        Identity::Empty => channels::after_boolean(
            TriMesh::default(),
            &subject.fates,
            tool.map(|built| &built.fates),
            Vec::new(),
        ),
    };
    Ok(MeshCut {
        built,
        unsettled: false,
        skin: None,
    })
}

/// `error` with the reason the exact referee refused to decide it.
fn refused_by_referee(error: GeomError, reason: Option<String>) -> GeomError {
    match error {
        GeomError::Degenerate(text) => GeomError::Degenerate(format!(
            "{text}; the exact referee refused: {}",
            reason.unwrap_or_default()
        )),
        other => other,
    }
}

/// Whether the difference or intersection `id` is only the subject of
/// further differences and intersections, and no root asked for: a chain of
/// clips goes on past it, so a skin within rounding in its result waits for
/// the root of the chain (#400, ADR 0102). A node any other node reads, as a
/// tool, in a union or an instance, ends its chain.
fn clip_chain_continues(
    graph: &GeometryGraph,
    regions: &Regions<'_>,
    id: NodeId,
    root: NodeId,
) -> bool {
    if id == root || graph.roots().contains(&id) {
        return false;
    }
    let mut continues = false;
    for (_, node) in graph.iter() {
        if !node.references().contains(&id) {
            continue;
        }
        match node {
            GeometryNode::SolidOperation(SolidOperation::Boolean {
                left,
                right,
                operator:
                    axiolid_core::BooleanOperator::Difference
                    | axiolid_core::BooleanOperator::Intersection,
            }) if *left == id && *right != id => continues = true,
            // The subject of an intersection with a region on its left (#326).
            GeometryNode::SolidOperation(SolidOperation::Boolean {
                left,
                right,
                operator: axiolid_core::BooleanOperator::Intersection,
            }) if *right == id && regions.has(*left) => {
                continues = true;
            }
            _ => return false,
        }
    }
    continues
}

/// `options` read in the frame a boolean is cut in (#356, #400): its
/// tolerance and chord budget are the frame's `local_budget`.
fn local_options(options: &ExecutionOptions, local_budget: Budget) -> GeomResult<ExecutionOptions> {
    options
        .clone()
        .with_tolerance(local_budget.tolerance)
        .with_chord_error(local_budget.chord)
        .ok_or_else(|| {
            GeomError::InvalidInput(format!(
                "chord budget {} is not positive and finite",
                local_budget.chord
            ))
        })
}

/// Convert world-space budgets into conservative instance-local space.
///
/// The tolerance and the chord budget are both lengths, so both shrink by
/// the transform's stretch bound ([`crate::deviation::stretch`], the same
/// factor the instance's deviation report grows by): a chord within `c` of
/// the curve in local space is within `stretch * c` of it in world space.
/// The chord budget is shrunk until that product, rounded, is within the
/// world budget ([`crate::deviation::shrunk_budget`]), so a source meshed
/// to its local budget reports at most the world one (#258). The chord
/// budget also leaves room for the rounding of placing the source's
/// vertices ([`placement_room`], at least `room`), which the instance's
/// bound adds (#295), and a budget at or below the floor at the
/// translation's scale is refused. The placed coordinate scale grows by
/// the translation, in the source's units ([`crate::snap::placed_scale`],
/// #312).
fn instance_local_budget(
    transform: Transform3,
    budget: Budget,
    room: Scalar,
) -> GeomResult<Budget> {
    let (tolerance, stretch) = instance_local_tolerance(transform, budget.tolerance)?;
    let room = placement_room(transform, budget.chord, room)?;
    Ok(Budget {
        tolerance,
        chord: crate::deviation::shrunk_budget_within(budget.chord, stretch, room),
        placed: crate::snap::placed_scale(transform, budget.placed).ok_or_else(|| {
            GeomError::InvalidInput("instance translation has no finite magnitude".into())
        })?,
    })
}

/// The share of what the floor leaves of a chord budget an instance's
/// source first leaves for the rounding of its own coordinates once placed
/// (#295): `2^-20`, enough for a source within about `10^9` budgets of the
/// origin, so most sources are meshed once.
const ROOM_SHARE: Scalar = 1.0 / 1_048_576.0;

/// The room an instance's source leaves in `chord` for its placement's
/// rounding (#295): the floor every placement by `transform` reaches, the
/// rounding at the origin ([`crate::deviation::placement_rounding`]), plus
/// [`ROOM_SHARE`] of what it leaves, and at least `shown`, the rounding a
/// first mesh of the source showed. A `chord` at or below the floor is
/// refused by name.
pub(crate) fn placement_room(
    transform: Transform3,
    chord: Scalar,
    shown: Scalar,
) -> GeomResult<Scalar> {
    let floor = crate::deviation::placement_rounding(&[Point3::ZERO], transform);
    refuse_below_rounding(chord, floor)?;
    Ok((floor + (chord - floor) * ROOM_SHARE).max(shown))
}

/// How far a source lies from its exact placement at a station whose frame
/// has a certified error (#405): the bound's displacement of the source's
/// farthest vertex, `point_error + axis_error R`, `R` the largest `|x| +
/// |y| + |z|` of its vertices, which covers every point of the mesh (`|q|_1`
/// is convex). `None` for any other accuracy, or a term that is not finite.
fn frame_term(accuracy: Option<PlacementAccuracy>, source: &TriMesh) -> Option<Scalar> {
    let Some(PlacementAccuracy::Certified(bound)) = accuracy else {
        return None;
    };
    let reach = source
        .positions
        .iter()
        .map(|p| p.x.abs() + p.y.abs() + p.z.abs())
        .fold(0.0, Scalar::max)
        * (1.0 + 1e-12);
    let frame = bound.displacement(Point3::new(reach, 0.0, 0.0));
    frame.is_finite().then_some(frame)
}

/// The rounding of placing `source`'s vertices with `transform` (#295),
/// refused by name when it reaches the chord budget: no mesh placed there
/// is within it.
fn placed_rounding(source: &TriMesh, transform: Transform3, chord: Scalar) -> GeomResult<Scalar> {
    let rounding = crate::deviation::placement_rounding(&source.positions, transform);
    refuse_below_rounding(chord, rounding)?;
    Ok(rounding)
}

/// [`crate::deviation::BUDGET_BELOW_PLACEMENT_ROUNDING`] unless `rounding`
/// is below `chord`.
pub(crate) fn refuse_below_rounding(chord: Scalar, rounding: Scalar) -> GeomResult<()> {
    if rounding < chord {
        return Ok(());
    }
    Err(GeomError::UnsupportedInput {
        backend: BackendId::new("scalar-compile"),
        operation: Operation::GraphCompilation,
        input: crate::deviation::BUDGET_BELOW_PLACEMENT_ROUNDING,
    })
}

/// The source and the transform of a node that places one: an `Instance`,
/// or an `InstanceAtStation` with its placement resolved (#264). `None` for
/// any other node.
fn placed_source(
    graph: &GeometryGraph,
    node: &GeometryNode,
    id: NodeId,
) -> GeomResult<Option<(NodeId, Transform3)>> {
    Ok(match node {
        GeometryNode::Instance(instance) => Some((instance.source, instance.transform)),
        GeometryNode::InstanceAtStation(placed) => Some((
            placed.source,
            crate::station::placement(graph, id)?.transform,
        )),
        _ => None,
    })
}

pub(crate) fn instance_local_tolerance(
    transform: Transform3,
    tolerance: Tolerance,
) -> GeomResult<(Tolerance, Scalar)> {
    let m = transform.matrix3;
    let max_scale = m
        .x_axis
        .length()
        .max(m.y_axis.length())
        .max(m.z_axis.length());
    if !max_scale.is_finite() || max_scale == 0.0 {
        return Err(GeomError::InvalidInput(
            "instance transform has no finite scale".into(),
        ));
    }
    let stretch = crate::deviation::stretch(transform);
    if !stretch.is_finite() {
        return Err(GeomError::InvalidInput(
            "instance transform scale overflow".into(),
        ));
    }
    Tolerance::new(tolerance.linear() / stretch, tolerance.angular())
        .map(|local| (local, stretch))
        .map_err(|error| GeomError::InvalidInput(error.to_string()))
}

/// How far a chord may sit from the curve it replaces (#165).
///
/// The caller's `ExecutionOptions::with_chord_error`, else the linear
/// tolerance. Kept separate from the tolerance because the tolerance is a
/// coincidence test: at `Tolerance::MILLIMETRE` it leaves a 5 mm arc four
/// chords per half turn, percent-level area error on a slot or a gutter.
/// Rebuilds a certified sweep may take to bring its bound within budget.
const REFINE_ROUNDS: usize = 3;

/// Each rebuild's chord is the last one shrunk by the overshoot and this
/// margin, so one round usually settles it.
const REFINE_MARGIN: Scalar = 0.9;

pub(crate) fn chord_error(options: &ExecutionOptions) -> Scalar {
    options
        .chord_error()
        .unwrap_or_else(|| options.tolerance().linear())
}

/// The tolerance with the chord budget as its linear part, for providers
/// that take one number for both (the CSG primitive tessellator).
fn chord_tolerance(options: &ExecutionOptions) -> GeomResult<Tolerance> {
    Tolerance::new(chord_error(options), options.tolerance().angular())
        .map_err(|error| GeomError::InvalidInput(error.to_string()))
}

/// The tolerance with HALF the chord budget as its linear part: the share
/// a revolution's steps round the axis get, the other half having chorded
/// its profile (#231).
fn half_tolerance(options: &ExecutionOptions) -> GeomResult<Tolerance> {
    Tolerance::new(0.5 * chord_error(options), options.tolerance().angular())
        .map_err(|error| GeomError::InvalidInput(error.to_string()))
}

/// The message of a half-space built on its own (#326): it has a finite
/// stand-in only against a subject.
const HALF_SPACE_WITHOUT_SUBJECT: &str = "half-space without a subject to size it by: it is \
     unbounded on its kept side, and is compiled only as the tool of a difference or an \
     intersection";

/// The message of a boolean with a region operand that is no clip (#326).
const UNCLIPPED_HALF_SPACE: &str = "boolean of a half-space that is not a clip: a half-space, or \
     an intersection of half-spaces, is unbounded, and is compiled only as the tool of a \
     difference or an intersection, or as one operand of an intersection with a solid";

fn unsupported_input(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: crate::BACKEND_ID,
        operation: Operation::GraphCompilation,
        input,
    }
}

/// A clip: a difference or an intersection with a region as its tool
/// (#326), the operands read in the order that makes the region the tool.
struct Clip {
    subject: NodeId,
    tool: NodeId,
    leaves: Vec<(ClipTool, Transform3)>,
}

/// Which nodes of `graph` are regions (#326): a half-space, an instance of a
/// region, or an intersection of regions, decided on demand and remembered
/// in the compilation's cache ([`crate::regions`]).
struct Regions<'a> {
    graph: &'a GeometryGraph,
    memo: &'a crate::regions::RegionMemo,
}

impl<'a> Regions<'a> {
    fn new(graph: &'a GeometryGraph, memo: &'a crate::regions::RegionMemo) -> Self {
        Self { graph, memo }
    }

    fn has(&self, id: NodeId) -> bool {
        self.memo
            .has(self.graph, id, |leaf| clip_tool(self.graph, leaf).is_some())
    }
}

/// The clip a boolean is, or `None`. These operators stay within the
/// finite subject, so a slab covering its bounds is an exact finite
/// stand-in for the region (#325). Union and XOR are unbounded and remain
/// unsupported. An intersection commutes, so a region on its left clips the
/// solid on its right.
fn clip_of(
    graph: &GeometryGraph,
    regions: &Regions<'_>,
    (left, right): (NodeId, NodeId),
    operator: axiolid_core::BooleanOperator,
) -> GeomResult<Option<Clip>> {
    use axiolid_core::BooleanOperator::{Difference, Intersection};
    if !matches!(operator, Difference | Intersection) {
        return Ok(None);
    }
    let (subject, tool) = match (regions.has(left), regions.has(right)) {
        (false, true) => (left, right),
        (true, false) if operator == Intersection => (right, left),
        _ => return Ok(None),
    };
    let leaves = region(graph, tool)?;
    Ok(Some(Clip {
        subject,
        tool,
        leaves,
    }))
}

/// The half-spaces of the region `id` (#326), each with the placement above
/// it, left operand first. `id` must be a region ([`Regions`]). Iterative,
/// with the leaves bounded by the graph's size: a graph that shares a
/// sub-region is read through many paths.
fn region(graph: &GeometryGraph, id: NodeId) -> GeomResult<Vec<(ClipTool, Transform3)>> {
    let mut leaves = Vec::new();
    let mut work = vec![(id, Transform3::IDENTITY)];
    while let Some((mut id, mut outer)) = work.pop() {
        loop {
            let node = graph.get(id).ok_or_else(|| {
                GeomError::InvalidInput(format!("operand {id:?} does not belong to this graph"))
            })?;
            match node {
                GeometryNode::Instance(instance) => {
                    outer *= instance.transform;
                    id = instance.source;
                }
                GeometryNode::SolidOperation(SolidOperation::Boolean {
                    left,
                    right,
                    operator: axiolid_core::BooleanOperator::Intersection,
                }) => {
                    work.push((*right, outer));
                    work.push((*left, outer));
                    break;
                }
                _ => {
                    let tool = clip_tool(graph, id).ok_or_else(|| {
                        GeomError::InvalidInput(format!("operand {id:?} is not a half-space"))
                    })?;
                    if leaves.len() >= graph.len() {
                        return Err(unsupported_input(
                            "half-space intersection with more half-spaces than its graph has nodes",
                        ));
                    }
                    leaves.push((tool, outer));
                    break;
                }
            }
        }
    }
    Ok(leaves)
}

/// A boolean tool that is a half-space, bounded or not.
#[derive(Debug, Clone, Copy)]
enum ClipTool {
    /// One side of a plane.
    Unbounded(axiolid_primitive::HalfSpace),
    /// One side of a plane within the prism of a boundary.
    Bounded {
        half_space: NodeId,
        boundary: NodeId,
        placement: Transform3,
        prism_axis: PrismAxis,
    },
}

/// The half-space `id` is, bounded or not; `None` when it is none. The
/// instances above one are the region's to read ([`region`]).
fn clip_tool(graph: &GeometryGraph, id: NodeId) -> Option<ClipTool> {
    match graph.get(id)? {
        GeometryNode::HalfSpace(half_space) => Some(ClipTool::Unbounded(*half_space)),
        GeometryNode::SolidOperation(operation) => {
            let (half_space, boundary, frame, prism_axis) =
                crate::half_space_boundary::bounded(operation)?;
            Some(ClipTool::Bounded {
                half_space,
                boundary,
                placement: frame,
                prism_axis,
            })
        }
        _ => None,
    }
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
        self.evaluate(graph, root, options, &mut cache, Measure::Skip)
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
        let built = self.evaluate(graph, root, options, &mut cache, Measure::Skip)?;
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
            destination.push(
                self.evaluate(graph, root, options, &mut cache, Measure::Skip)?
                    .mesh,
            );
        }
        Ok(())
    }
}

/// Curve nodes to polylines (#311): `Curve2` (lifted to `z = 0`),
/// `Curve3` and `OpenProfile`, and their placements, each bisected until
/// every span's certified chord bound fits the budget. See
/// `crate::polyline` for what is bounded how and what is refused.
impl<B: MeshBoolean> PolylineCompiler for ReferenceMeshCompiler<B> {
    fn compile_polyline(
        &self,
        graph: &GeometryGraph,
        root: NodeId,
        options: &ExecutionOptions,
    ) -> GeomResult<PolylineOutcome> {
        crate::polyline::compile(graph, root, options)
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

    /// Compile one root and report how far its exact surface may lie from
    /// the mesh (#232).
    ///
    /// The outcome is [`MeshCompiler::compile_mesh_reported`]'s. The report
    /// carries a certified upper bound on the distance from every point of
    /// the exact surface to the mesh's triangles, in world units, the
    /// paths that contributed, and whether the bound is within the chord
    /// budget asked for; see [`crate::DeviationReport`] and the
    /// [`deviation`](crate::deviation) module for what each path's bound
    /// rests on. A path with no certified bound is named, and makes the
    /// overall bound `None`. A boolean is measured against the result
    /// [`crate::ReferenceExactCompiler`] builds for it, where it builds one
    /// (#235): that costs an exact compilation and a certified search for
    /// each boolean whose result the root emits (not for the inner
    /// booleans of a chain).
    ///
    /// An inherent method rather than a [`MeshCompiler`] one: the bound is
    /// this compiler's claim about its own meshes, and another compiler
    /// would have to prove its own before it could make it.
    ///
    /// # Errors
    ///
    /// As [`MeshCompiler::compile_mesh_reported`].
    pub fn compile_mesh_with_deviation(
        &self,
        graph: &GeometryGraph,
        root: NodeId,
        options: &ExecutionOptions,
    ) -> GeomResult<(CompileOutcome, crate::DeviationReport)> {
        self.admit_budget(options)?;
        let mut cache = Cache::new();
        let emitted = crate::deviation::emitted_booleans(graph, root);
        let built = self.evaluate(
            graph,
            root,
            options,
            &mut cache,
            Measure::Deviation(&emitted),
        )?;
        let report = built.deviation.report(chord_error(options));
        Ok((
            CompileOutcome::tracked(built.mesh, built.fates.into_vec()).with_closure(built.closure),
            report,
        ))
    }
}

/// Convert an authored polygon mesh to triangles (#160).
///
/// A plain triangle keeps its exact corner order. Any other face -- an
/// n-gon, a concave face, a face with holes (IFC4
/// `IfcIndexedPolygonalFaceWithVoids`) -- is triangulated in its own
/// plane by [`crate::planar::triangulate_authored_polygon`]. Positions
/// are kept as authored and shared, so the output welds exactly where the
/// input did; no position is moved. The only positions added are the
/// Steiner points a keyhole face places inside itself so that its seam is
/// no triangle edge (#283), each used by that face alone. A corner that
/// one face's
/// triangulation inserts into a ring edge (a hole touching the outer
/// edge at an interior point) is inserted into every face sharing that
/// edge by index, a triangle face included, so the faces still meet
/// along the same pieces of it ([`crate::weld`], #265).
///
/// A face with no area, or whose rings cross in its plane, is refused
/// with an error naming its index. A face off its plane by more than
/// the linear tolerance has no unique surface; it is triangulated in
/// its fit plane all the same (#254), and the largest warp of any such
/// face (the width of the slab its corners span about that plane, #261)
/// is returned with the mesh, `None` when every face is planar.
/// [`crate::planar`] states what that warp bounds; the caller reports
/// it as a certified deviation, never as exact. So is the largest
/// [`crate::planar::steiner_offset`] of a keyhole face whose Steiner
/// points lie off an axis plane (#283), `None` when there is none.
pub(crate) fn compile_authored_polygons(
    mesh: &axiolid_mesh::PolygonMesh,
    options: &ExecutionOptions,
) -> GeomResult<(TriMesh, Option<Scalar>, Option<Scalar>)> {
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
    let splits = authored_edge_splits(mesh);
    let mut warp: Option<Scalar> = None;
    let mut steiner_offset: Option<Scalar> = None;
    for (face_index, face) in mesh.faces.iter().enumerate() {
        let mut rings: Vec<Vec<u32>> = face_rings(face).cloned().collect();
        // Weld through every corner a face inserts into a shared edge
        // (#265): this face takes each such corner on its own copy of
        // the edge, a triangle face included, which then has four.
        splits.split(&mut rings, |i| (i, mesh.positions[i as usize]), |i, _| i);
        if rings.len() == 1 && rings[0].len() == 3 {
            indices.extend_from_slice(&rings[0]);
            continue;
        }
        // A repeated corner (an exporter's closing point, a doubled
        // corner) needs no handling here: `planar::clip_projected` drops
        // coincident consecutive points itself
        // (`a_repeated_closing_corner_is_not_a_triangle`).
        let points: Vec<Vec<axiolid_core::Point3>> = rings
            .iter()
            .map(|ring| ring.iter().map(|&i| mesh.positions[i as usize]).collect())
            .collect();
        let views: Vec<&[axiolid_core::Point3]> = points.iter().map(Vec::as_slice).collect();
        let local = crate::planar::triangulate_authored_polygon(&views, linear)
            .map_err(|refusal| crate::planar::face_error(face_index, refusal))?;
        if let Some(face_warp) = local.warp {
            warp = Some(warp.map_or(face_warp, |w| w.max(face_warp)));
        }
        if let Some(offset) = local.steiner_offset {
            steiner_offset = Some(steiner_offset.map_or(offset, |w| w.max(offset)));
        }
        let mut corners: Vec<u32> = rings.into_iter().flatten().collect();
        // Steiner points follow the corners, each a new position of this
        // face alone (#283).
        for point in local.steiner {
            let index = u32::try_from(positions.len()).map_err(|_| GeomError::BudgetExceeded {
                resource: "authored polygon positions",
            })?;
            positions.push(point);
            corners.push(index);
        }
        indices.extend(local.indices.into_iter().map(|corner| corners[corner]));
    }

    let triangles = TriMesh::new(positions, indices);
    triangles.validate_structure().map_err(|error| {
        GeomError::InvalidInput(format!("invalid authored polygon mesh: {error}"))
    })?;
    Ok((triangles, warp, steiner_offset))
}

/// The corners the authored polygon faces of `mesh` insert into their own
/// ring edges, recorded against those edges by position index (#265,
/// [`crate::weld`]). Every index is in range: the caller checked.
fn authored_edge_splits(mesh: &axiolid_mesh::PolygonMesh) -> crate::weld::EdgeSplits<u32> {
    let mut splits = crate::weld::EdgeSplits::new();
    for face in &mesh.faces {
        // One triangle has no corner inside its own edge.
        if face.outer.len() == 3 && face.holes.is_empty() {
            continue;
        }
        let points: Vec<Vec<axiolid_core::Point3>> = face_rings(face)
            .map(|ring| ring.iter().map(|&i| mesh.positions[i as usize]).collect())
            .collect();
        let views: Vec<&[axiolid_core::Point3]> = points.iter().map(Vec::as_slice).collect();
        let touches = crate::planar::authored_edge_touches(&views);
        if !touches.is_empty() {
            let corners: Vec<(u32, axiolid_core::Point3)> = face_rings(face)
                .flatten()
                .map(|&i| (i, mesh.positions[i as usize]))
                .collect();
            splits.record(&corners, &touches);
        }
    }
    splits
}

/// A face's outer ring followed by its holes.
fn face_rings(face: &axiolid_mesh::PolygonFace) -> impl Iterator<Item = &Vec<u32>> {
    std::iter::once(&face.outer).chain(face.holes.iter())
}

/// A boolean needs two volumes; a surface model has none (#161), and
/// neither has a declared solid whose mesh does not close (#265).
///
/// Refused rather than handed to the provider, which would see a closed
/// surface model as a valid solid and return a confident, meaningless
/// result.
fn refuse_surface_operand(built: &Built, role: &'static str) -> GeomResult<()> {
    match built.closure {
        axiolid_mesh_compile_contract::MeshClosure::Solid => Ok(()),
        axiolid_mesh_compile_contract::MeshClosure::OpenSolid => {
            Err(GeomError::InvalidInput(format!(
                "boolean {role} is a solid whose mesh is not closed: it encloses no volume to \
                 combine"
            )))
        }
        _ => Err(GeomError::InvalidInput(format!(
            "boolean {role} is a surface model: it encloses no volume to combine"
        ))),
    }
}

/// An authored mesh, with its closure read from its structure (#161).
///
/// A mesh node carries no "this is a solid" declaration the way a B-rep
/// does, so the geometry is the only evidence: a closed, consistently
/// wound two-manifold bounds a solid, anything else (an open face set, a
/// single sheet) is a surface.
///
/// Closure is a property of the authored connectivity, so it is read from
/// the corner indices ([`axiolid_mesh::EdgeAdjacency`]), not from
/// `audit_mesh`. The audit drops every triangle below an area threshold
/// before it counts edges, and real exports carry legitimate small faces
/// (a 1 mm reveal, the end of a thin lining); dropping them opened
/// otherwise closed shells, so 522 authored-closed face sets on real
/// models were reported as surfaces. Only a triangle that repeats a
/// corner index is skipped: it covers no edge pair at all.
fn authored_mesh(mesh: TriMesh, _options: &ExecutionOptions) -> Built {
    let closure = if authored_mesh_is_closed(&mesh) {
        axiolid_mesh_compile_contract::MeshClosure::Solid
    } else {
        axiolid_mesh_compile_contract::MeshClosure::Surface
    };
    // An authored mesh is its own exact surface (#232).
    Built::with_closure(mesh, closure).with_deviation(crate::deviation::Deviation::one(
        crate::deviation::DeviationPath::AuthoredMesh,
        "",
        crate::deviation::DeviationBound::Proven(0.0),
    ))
}

/// Every edge of the index connectivity is shared by exactly two
/// triangles traversing it in opposite directions, the mesh has at least
/// one such triangle, and every position is finite.
fn authored_mesh_is_closed(mesh: &TriMesh) -> bool {
    if mesh.indices.is_empty() || mesh.positions.iter().any(|p| !p.is_finite()) {
        return false;
    }
    let adjacency = axiolid_mesh::EdgeAdjacency::build(mesh);
    adjacency.edge_count() > 0 && adjacency.is_closed_two_manifold()
}

/// The linear part an operand is placed by through its whole chain of
/// instances, `outer` being its own transform: the product down to the
/// first node that is not an instance. `None` when the operand is not an
/// instance (a boolean's frame was checked when it was cut).
fn chain_linear(graph: &GeometryGraph, id: NodeId, outer: Transform3) -> Option<Transform3> {
    let Some(GeometryNode::Instance(first)) = graph.get(id) else {
        return None;
    };
    let mut composed = outer;
    let mut source = first.source;
    while let Some(GeometryNode::Instance(next)) = graph.get(source) {
        composed *= next.transform;
        source = next.source;
    }
    Some(composed)
}
