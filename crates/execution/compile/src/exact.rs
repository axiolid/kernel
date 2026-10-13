//! Reference graph-to-exact-B-rep compiler.
//!
//! Exact and mesh compilation are separate result domains. One exact batch owns
//! a `NodeId -> ExactBRep` memo table; a discrete value cannot enter that cache.
//!
//! It never falls back to mesh compilation. Extrusions (either way along
//! the profile normal, #275), revolutions,
//! booleans of sharp rectangle prisms along +z, differences of placed
//! extrusions ([`boolean`], #228), clips of those by half-spaces ([`clip`],
//! #234), disks swept along one segment or one arc, runs of sections at
//! stations with a closed form (#247: a constant section along a line or a
//! horizontal circle, a linear taper along a line; `crate::station`), and
//! instances of any of these under a rigid transform are compiled
//! exactly; every other family is refused with `GeomError::UnsupportedInput`
//! naming it.
//!
//! An instance places its source's exact B-rep with
//! [`ExactBRep::transformed`] (#223): a rotation, a reflection and a
//! translation are exact, and a scale or shear is refused, never
//! approximated. An instance at a station (#264) is placed the same way
//! when its frame is exact (a line basis) and refused by name otherwise.
//!
//! [`ReferenceExactCompiler::compile_exact_with_report`] also says, per
//! body, whether any boolean on the way read a feature within the
//! tolerance (#236): an exact [`BooleanReport`] means the body is the exact
//! result of its operands as given; otherwise it was built from operands
//! moved by at most the report's magnitudes, within the tolerance.
//!
//! Only for a boolean's deviation report (#266, crate-internal), a
//! compilation may flatten mesh operands: an authored mesh, or a B-rep
//! of planes and warped faces without a surface, becomes the exact
//! polyhedron of the triangles the mesh compiler makes of it
//! (`polyhedral`), and the largest slab width a warped face was flattened
//! by is reported with the body. That is the perturbed operand of
//! ADR 0080, not the exact surface of a warped face, so the public
//! compiler never does it and refuses meshes by name.

use std::collections::{HashMap, HashSet};

use axiolid_brep::{ExactBRep, TransformError};
use axiolid_brep_boolean::BooleanReport;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_construct::swept_disk_exact::{
    swept_disk_along_arc_exact, swept_disk_along_line_exact,
};
use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, ExecutionOptions, ExecutionTarget, GeomError,
    GeomResult, Operation,
};
use axiolid_core::Scalar;
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_model::{GeometryGraph, GeometryNode, NodeId, SolidOperation};

mod boolean;
mod clip;
mod polyhedral;

/// Scalar reference implementation of the exact-compilation capability.
#[derive(Debug, Clone, Copy, Default)]
pub struct ReferenceExactCompiler;

impl ReferenceExactCompiler {
    /// Stable identity of this backend.
    pub const ID: BackendId = BackendId::new("scalar-exact-compile");

    /// Construct the reference exact compiler.
    pub const fn new() -> Self {
        Self
    }

    /// [`ExactCompiler::compile_exact`], with the within-tolerance
    /// decisions the body's booleans took (#236), merged over every boolean
    /// in its graph (each kind once, with its worst magnitude).
    ///
    /// An exact report ([`BooleanReport::is_exact`]) means no boolean read
    /// anything within the tolerance: the body is the exact boolean of its
    /// operands as given, and may be cited as exact. Otherwise it is the
    /// exact result of operands moved by at most [`BooleanReport::linear`]
    /// and turned by at most [`BooleanReport::angular`], and a consumer
    /// must widen distances by that much. Bodies without booleans, and the
    /// integer-exact prism path, always report exact.
    ///
    /// Either way, constructed points closer than
    /// [`BooleanReport::rounding_floor`] count as one point and are not
    /// reported, so a consumer widens distances by that floor too (#244).
    /// The merged floor is the largest floor of any boolean or clip beneath
    /// the body; an instance is rigid and keeps its source's floor, in
    /// model units, unchanged. Zero means no general boolean or clip ran.
    ///
    /// # Errors
    ///
    /// As [`ExactCompiler::compile_exact`].
    pub fn compile_exact_with_report(
        &self,
        graph: &GeometryGraph,
        root: NodeId,
        options: &ExecutionOptions,
    ) -> GeomResult<(ExactBRep, BooleanReport)> {
        let mut compilation = ExactCompilation::new(graph, options);
        let body = compilation.compile(root)?;
        Ok((body, compilation.report_of(root)))
    }

    /// [`Self::compile_exact_with_report`] for a boolean's deviation
    /// report (#266): mesh operands (authored meshes, and B-reps whose
    /// faces are all planes or warped faces without a surface) are taken
    /// as the exact polyhedra of their triangle meshes, and a union or
    /// intersection is taken when an operand is such a mesh. Also returns
    /// the largest slab width any flattened operand beneath `root` was
    /// moved by (zero when none was warped). The public compiler keeps
    /// refusing meshes: their triangles are one reading of a warped face,
    /// not its exact surface.
    pub(crate) fn compile_flattened_with_report(
        &self,
        graph: &GeometryGraph,
        root: NodeId,
        options: &ExecutionOptions,
    ) -> GeomResult<(ExactBRep, BooleanReport, Scalar)> {
        let mut compilation = ExactCompilation::new(graph, options);
        compilation.flatten = true;
        let body = compilation.compile(root)?;
        Ok((
            body,
            compilation.report_of(root),
            compilation.flattened_of(root),
        ))
    }

    /// [`Self::compile_exact_with_report`] for several roots in one
    /// compilation, sharing its cache, in root order.
    ///
    /// # Errors
    ///
    /// As [`ExactCompiler::compile_exact_batch_into`].
    pub fn compile_exact_batch_with_reports(
        &self,
        graph: &GeometryGraph,
        roots: &[NodeId],
        options: &ExecutionOptions,
    ) -> GeomResult<Vec<(ExactBRep, BooleanReport)>> {
        let mut compilation = ExactCompilation::new(graph, options);
        let mut out = Vec::with_capacity(roots.len());
        for &root in roots {
            let body = compilation.compile(root)?;
            out.push((body, compilation.report_of(root)));
        }
        Ok(out)
    }
}

impl Backend for ReferenceExactCompiler {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(Self::ID, ExecutionTarget::PortableCpu)
    }
}

impl ExactCompiler for ReferenceExactCompiler {
    fn compile_exact(
        &self,
        graph: &GeometryGraph,
        root: NodeId,
        options: &ExecutionOptions,
    ) -> GeomResult<ExactBRep> {
        ExactCompilation::new(graph, options).compile(root)
    }

    fn compile_exact_batch_into(
        &self,
        graph: &GeometryGraph,
        roots: &[NodeId],
        options: &ExecutionOptions,
        destination: &mut Vec<ExactBRep>,
    ) -> GeomResult<()> {
        destination.reserve(roots.len());
        let mut compilation = ExactCompilation::new(graph, options);
        for &root in roots {
            destination.push(compilation.compile(root)?);
        }
        Ok(())
    }
}

struct ExactCompilation<'a> {
    graph: &'a GeometryGraph,
    options: &'a ExecutionOptions,
    cache: HashMap<NodeId, ExactBRep>,
    active: HashSet<NodeId>,
    cache_hits: usize,
    evaluated_nodes: usize,
    /// Per compiled node with a boolean beneath it: what its booleans read
    /// within tolerance (#236), and their rounding floor (#244). Absent
    /// means exact, with no floor.
    reports: HashMap<NodeId, BooleanReport>,
    /// The report of the general boolean `compile_boolean` just ran.
    boolean_report: Option<BooleanReport>,
    /// Whether mesh operands are taken as the exact polyhedra of their
    /// triangle meshes (#266,
    /// [`ReferenceExactCompiler::compile_flattened_with_report`]).
    flatten: bool,
    /// Per compiled node with a warped mesh beneath it: the largest slab
    /// width a flattened operand was moved by (#266). Absent means none.
    flattened: HashMap<NodeId, Scalar>,
}

impl<'a> ExactCompilation<'a> {
    fn new(graph: &'a GeometryGraph, options: &'a ExecutionOptions) -> Self {
        Self {
            graph,
            options,
            cache: HashMap::new(),
            active: HashSet::new(),
            cache_hits: 0,
            evaluated_nodes: 0,
            reports: HashMap::new(),
            boolean_report: None,
            flatten: false,
            flattened: HashMap::new(),
        }
    }

    /// The largest slab width of the meshes flattened beneath a node.
    fn flattened_of(&self, node: NodeId) -> Scalar {
        self.flattened.get(&node).copied().unwrap_or(0.0)
    }

    /// Record the largest slab width flattened beneath a node, if any.
    fn set_flattened(&mut self, node: NodeId, width: Scalar) {
        if width > 0.0 {
            self.flattened.insert(node, width);
        }
    }

    /// A mesh operand as the exact polyhedron of its triangle mesh, with
    /// the slab width its warped faces were flattened by (#266).
    fn flattened_mesh(&self, node: &GeometryNode) -> GeomResult<(ExactBRep, Scalar)> {
        match node {
            GeometryNode::PolygonMesh(mesh) => {
                // A keyhole face's Steiner points off an axis plane (#283)
                // move the mesh off the cover by at most their slab width,
                // the same kind of bound as a warp: flattened by either.
                let (triangles, warp, steiner_offset) =
                    crate::compiler::compile_authored_polygons(mesh, self.options)?;
                let width = warp.unwrap_or(0.0).max(steiner_offset.unwrap_or(0.0));
                Ok((polyhedral::from_triangles(&triangles)?, width))
            }
            GeometryNode::TriMesh(mesh) => Ok((polyhedral::from_triangles(mesh)?, 0.0)),
            GeometryNode::BRep(brep) => {
                let (triangles, closure, deviation) = crate::brep::tessellate(
                    brep,
                    self.graph,
                    self.options.tolerance(),
                    crate::compiler::chord_error(self.options),
                )?;
                if closure != axiolid_mesh_compile_contract::MeshClosure::Solid {
                    return Err(unsupported(
                        "exact boolean over a B-rep operand that bounds no solid",
                    ));
                }
                let warp = deviation.polyhedral_warp().ok_or_else(|| {
                    unsupported(
                        "exact boolean over a B-rep operand with a face that is neither a plane \
                         with straight edges nor a warped face without a surface",
                    )
                })?;
                Ok((polyhedral::from_triangles(&triangles)?, warp))
            }
            other => Err(unsupported(exact_input_family(other))),
        }
    }

    /// What the booleans beneath a compiled node read within tolerance.
    fn report_of(&self, node: NodeId) -> BooleanReport {
        self.reports.get(&node).cloned().unwrap_or_default()
    }

    /// Record a node's report, if anything was read or a rounding floor
    /// applied (#244).
    fn set_report(&mut self, node: NodeId, report: BooleanReport) {
        if !report.is_exact() || report.rounding_floor() > 0.0 {
            self.reports.insert(node, report);
        }
    }

    fn compile(&mut self, root: NodeId) -> GeomResult<ExactBRep> {
        if let Some(cached) = self.cache.get(&root) {
            self.cache_hits += 1;
            return Ok(cached.clone());
        }
        if !self.active.insert(root) {
            return Err(GeomError::InvalidInput(format!(
                "exact compilation cycle reached at node {root:?}"
            )));
        }

        self.evaluated_nodes += 1;
        let result = self.compile_uncached(root);
        self.active.remove(&root);
        let exact = result?;
        self.cache.insert(root, exact.clone());
        Ok(exact)
    }

    fn compile_uncached(&mut self, root: NodeId) -> GeomResult<ExactBRep> {
        let node = self.graph.get(root).ok_or_else(|| {
            GeomError::InvalidInput(format!("node {root:?} does not belong to this graph"))
        })?;

        match node {
            GeometryNode::SolidOperation(SolidOperation::Extrusion {
                profile,
                direction,
                depth,
            }) => {
                let profile_node = self.graph.get(*profile).ok_or_else(|| {
                    GeomError::InvalidInput(format!(
                        "extrusion profile {profile:?} does not belong to this graph"
                    ))
                })?;
                let GeometryNode::Profile(profile) = profile_node else {
                    return Err(unsupported("extrusion profile reference"));
                };
                extrude_profile_exact(profile, *direction, *depth, self.options.tolerance())
                    .map_err(remap_construction_error)
            }
            GeometryNode::SolidOperation(SolidOperation::Revolution {
                profile,
                axis_origin,
                axis_direction,
                angle,
            }) => {
                let profile_node = self.graph.get(*profile).ok_or_else(|| {
                    GeomError::InvalidInput(format!(
                        "revolution profile {profile:?} does not belong to this graph"
                    ))
                })?;
                let GeometryNode::Profile(profile) = profile_node else {
                    return Err(unsupported("revolution profile reference"));
                };
                revolve_profile_exact(
                    profile,
                    *axis_origin,
                    *axis_direction,
                    *angle,
                    self.options.tolerance(),
                )
                .map_err(remap_construction_error)
            }
            GeometryNode::Instance(instance) => {
                let instance = *instance;
                let source = self.compile(instance.source)?;
                self.set_report(root, self.report_of(instance.source));
                self.set_flattened(root, self.flattened_of(instance.source));
                source
                    .transformed(&instance.transform)
                    .map_err(|error| match error {
                        TransformError::NotRigid => {
                            unsupported("exact instance under a scaled or sheared transform")
                        }
                        TransformError::Unsupported(what) => unsupported(what),
                        other => GeomError::InvalidInput(format!("instance transform: {other}")),
                    })
            }
            // An instance in a station's frame (#264): exact only where the
            // frame is (a line basis); any other frame is read at a distance
            // resolved by an estimate, which an exact body must not carry.
            GeometryNode::InstanceAtStation(placed) => {
                let source_id = placed.source;
                let placement = crate::station::placement(self.graph, root)?;
                if !placement.exact {
                    return Err(unsupported(
                        "exact instance at a station whose frame is not exact: the basis is \
                         not a line, so its distance is read by an estimate",
                    ));
                }
                let source = self.compile(source_id)?;
                self.set_report(root, self.report_of(source_id));
                self.set_flattened(root, self.flattened_of(source_id));
                source
                    .transformed(&placement.transform)
                    .map_err(|error| match error {
                        TransformError::Unsupported(what) => unsupported(what),
                        other => GeomError::InvalidInput(format!("instance at a station: {other}")),
                    })
            }
            // A single segment or arc has no corners, so a fillet radius
            // has nothing to round there; one with corners is refused.
            GeometryNode::SolidOperation(SolidOperation::SweptDisk {
                directrix,
                radius,
                inner_radius,
                parameter_range,
                fillet_radius: _,
            }) => {
                let tolerance = self.options.tolerance();
                match crate::directrix::exact(
                    self.graph,
                    *directrix,
                    *parameter_range,
                    self.options,
                    unsupported,
                )? {
                    crate::directrix::ExactDirectrix::Segment(start, end) => {
                        swept_disk_along_line_exact(start, end, *radius, *inner_radius, tolerance)
                    }
                    crate::directrix::ExactDirectrix::Arc(circle, span) => {
                        swept_disk_along_arc_exact(&circle, span, *radius, *inner_radius, tolerance)
                    }
                }
                .map_err(remap_construction_error)
            }
            GeometryNode::SolidOperation(SolidOperation::Boolean {
                left,
                right,
                operator,
            }) => {
                let (left, right, operator) = (*left, *right, *operator);
                let body = self.compile_boolean(left, right, operator)?;
                let own = self.boolean_report.take().unwrap_or_default();
                let report = self
                    .report_of(left)
                    .merged(&self.report_of(right))
                    .merged(&own);
                self.set_report(root, report);
                let width = self.flattened_of(left).max(self.flattened_of(right));
                self.set_flattened(root, width);
                Ok(body)
            }
            // Only for a boolean's report (#266): a mesh operand flattened
            // to the exact polyhedron of its triangles.
            GeometryNode::PolygonMesh(_) | GeometryNode::TriMesh(_) | GeometryNode::BRep(_)
                if self.flatten =>
            {
                let (body, width) = self.flattened_mesh(node)?;
                self.set_flattened(root, width);
                Ok(body)
            }
            // A constant section along a line or a horizontal circle, and a
            // linear taper along a line, have closed forms (#247).
            GeometryNode::SolidOperation(SolidOperation::StationedSpine {
                directrix,
                sections,
                frame,
            }) => {
                let sections: Vec<_> = sections.iter().copied().map(Into::into).collect();
                crate::station::stationed_spine_exact(
                    self.graph,
                    *directrix,
                    &sections,
                    *frame,
                    self.options,
                    unsupported,
                )
            }
            GeometryNode::SolidOperation(SolidOperation::SectionsAtStations {
                directrix,
                sections,
                frame,
            }) => crate::station::stationed_spine_exact(
                self.graph,
                *directrix,
                sections,
                *frame,
                self.options,
                unsupported,
            ),
            _ => Err(unsupported(exact_input_family(node))),
        }
    }
}

fn remap_construction_error(error: GeomError) -> GeomError {
    match error {
        GeomError::UnsupportedInput { input, .. } => unsupported(input),
        error => error,
    }
}

/// Resolve a swept disk's directrix exactly as [`ReferenceExactCompiler`]
/// reads it for `SolidOperation::SweptDisk` (#230).
///
/// `directrix` and `parameter_range` are the swept disk's own fields. The
/// result is one straight segment or one circular arc, read with the same
/// trim, sense and range conventions as the mesh compiler's sampled path, so
/// an exact boundary built from it matches what both compilers produce.
/// Accepted: a bounded line, a two-point polyline, a circle or a sub-range of
/// one, a trim of a line or circle (across the seam included), and a
/// one-segment composite. Refused with `GeomError::UnsupportedInput` naming
/// the input: directrices with corners, other curve families and unbounded
/// lines; invalid or empty ranges are `InvalidInput` or `Degenerate`.
pub fn exact_directrix(
    graph: &GeometryGraph,
    directrix: NodeId,
    parameter_range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<crate::directrix::ExactDirectrix> {
    crate::directrix::exact(graph, directrix, parameter_range, options, unsupported)
}

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: ReferenceExactCompiler::ID,
        operation: Operation::GraphCompilation,
        input,
    }
}

fn exact_input_family(node: &GeometryNode) -> &'static str {
    match node {
        GeometryNode::Point2(_) => "2D point",
        GeometryNode::Point3(_) => "3D point",
        GeometryNode::Vector2(_) => "2D vector",
        GeometryNode::Vector3(_) => "3D vector",
        GeometryNode::Frame2(_) => "2D frame",
        GeometryNode::Frame3(_) => "3D frame",
        GeometryNode::Transform(_) => "transform",
        GeometryNode::PointList2(_) => "2D point list",
        GeometryNode::PointList3(_) => "3D point list",
        GeometryNode::Primitive(_) => "primitive",
        GeometryNode::Curve2(_) => "2D curve",
        GeometryNode::Curve3(_) => "3D curve",
        GeometryNode::CurveRelation(_) => "curve relation",
        GeometryNode::PointOnCurve(_) => "point on curve",
        GeometryNode::Surface(_) => "surface",
        GeometryNode::SurfaceRelation(_) => "surface relation",
        GeometryNode::PointOnSurface(_) => "point on surface",
        GeometryNode::Profile(_) => "profile",
        GeometryNode::OpenProfile(_) => "open profile",
        GeometryNode::HalfSpace(_) => "half-space",
        GeometryNode::BRep(_) => "B-rep",
        GeometryNode::PolygonMesh(_) => "polygon mesh",
        GeometryNode::TriMesh(_) => "triangle mesh",
        GeometryNode::BoundingBox(_) => "bounding box",
        GeometryNode::SolidOperation(operation) => solid_operation_family(operation),
        GeometryNode::Instance(_) => "instance",
        GeometryNode::Collection(_) => "collection",
        GeometryNode::CurveStation(_) => "curve station",
        GeometryNode::OrientedCurveStation(_) => "oriented curve station",
        GeometryNode::InstanceAtStation(_) => "instance at a station",
        _ => "unknown geometry node",
    }
}

pub(crate) fn solid_operation_family(operation: &SolidOperation) -> &'static str {
    match operation {
        SolidOperation::Extrusion { .. } => "extrusion",
        SolidOperation::TaperedExtrusion { .. } => "tapered extrusion",
        SolidOperation::Revolution { .. } => "revolution",
        SolidOperation::TaperedRevolution { .. } => "tapered revolution",
        SolidOperation::SweptDisk { .. } => "swept disk",
        SolidOperation::FixedReferenceSweep { .. } => "fixed-reference sweep",
        SolidOperation::SurfaceCurveSweep { .. } => "surface-curve sweep",
        SolidOperation::SectionedSpine { .. } => "sectioned spine",
        SolidOperation::StationedSpine { .. } | SolidOperation::SectionsAtStations { .. } => {
            "station-placed spine"
        }
        SolidOperation::Boolean { .. } => "boolean",
        SolidOperation::BoundedHalfSpace { .. } => "bounded half-space",
        _ => "unknown solid operation",
    }
}

#[cfg(test)]
mod family_name_tests {
    use super::{solid_operation_family, SOLID_FAMILY_NAMES};
    use axiolid_core::Vec3;
    use axiolid_model::{GeometryGraphBuilder, GeometryNode, SolidOperation};
    use axiolid_profile::{Profile, RectangleProfile};

    /// The exported table matches what the naming function actually returns.
    ///
    /// The table is hand-written, so it can drift from the match arm it
    /// mirrors. Walking the real function for a constructed variant turns
    /// drift into a failing test rather than a stale diagnostic, and the
    /// count assertion catches a family added without a name.
    #[test]
    fn the_exported_table_matches_the_naming_function() {
        assert_eq!(
            SOLID_FAMILY_NAMES.len(),
            11,
            "the declared family count changed; update the table"
        );
        for name in SOLID_FAMILY_NAMES {
            assert!(!name.is_empty());
            assert_ne!(*name, "unknown solid operation");
        }

        // NodeId has no public constructor, so build a real one.
        let mut builder = GeometryGraphBuilder::new();
        let profile = builder
            .push(GeometryNode::Profile(Profile::Rectangle(
                RectangleProfile {
                    x: 1.0,
                    y: 1.0,
                    thickness: None,
                    outer_radius: None,
                    inner_radius: None,
                },
            )))
            .expect("a rectangle profile is a valid node");

        let extrusion = SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 1.0,
        };
        assert!(
            SOLID_FAMILY_NAMES.contains(&solid_operation_family(&extrusion)),
            "the function returned a name the table does not carry"
        );
    }
}

#[cfg(test)]
mod tests {
    use axiolid_core::{Tolerance, Vec3};
    use axiolid_model::GeometryGraphBuilder;
    use axiolid_profile::{Profile, RectangleProfile};

    use super::*;

    #[test]
    fn duplicate_roots_hit_the_exact_result_cache() {
        let mut builder = GeometryGraphBuilder::new();
        let profile = builder
            .push(GeometryNode::Profile(Profile::Rectangle(
                RectangleProfile {
                    x: 2.0,
                    y: 1.0,
                    thickness: None,
                    outer_radius: None,
                    inner_radius: None,
                },
            )))
            .unwrap();
        let root = builder
            .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
                profile,
                direction: Vec3::Z,
                depth: 1.0,
            }))
            .unwrap();
        let graph = builder.finish(vec![root]).unwrap();
        let options = ExecutionOptions::new(Tolerance::METRE);
        let mut compilation = ExactCompilation::new(&graph, &options);

        let first = compilation.compile(root).unwrap();
        let second = compilation.compile(root).unwrap();

        assert_eq!(first, second);
        assert_eq!(compilation.evaluated_nodes, 1);
        assert_eq!(compilation.cache_hits, 1);
        assert_eq!(compilation.cache.len(), 1);
    }
}

/// Every diagnostic name `solid_operation_family` can return.
///
/// Kept beside that function so the two cannot drift: a family added there
/// without a name here fails the contract test. Exposed for that test rather
/// than for callers, who receive the name inside a refusal.
pub const SOLID_FAMILY_NAMES: &[&str] = &[
    "extrusion",
    "tapered extrusion",
    "revolution",
    "tapered revolution",
    "swept disk",
    "fixed-reference sweep",
    "surface-curve sweep",
    "sectioned spine",
    "station-placed spine",
    "boolean",
    "bounded half-space",
];
