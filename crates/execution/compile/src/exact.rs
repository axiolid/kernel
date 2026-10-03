//! Reference graph-to-exact-B-rep compiler.
//!
//! Exact and mesh compilation are separate result domains. One exact batch owns
//! a `NodeId -> ExactBRep` memo table; a discrete value cannot enter that cache.
//!
//! It never falls back to mesh compilation. Extrusions, revolutions,
//! booleans of sharp rectangle prisms along +z, differences of placed
//! extrusions ([`boolean`], #228), clips of those by half-spaces ([`clip`],
//! #234), disks swept along one segment or one arc, and instances of any
//! of these under a rigid transform are compiled
//! exactly; every other family is refused with `GeomError::UnsupportedInput`
//! naming it.
//!
//! An instance places its source's exact B-rep with
//! [`ExactBRep::transformed`] (#223): a rotation, a reflection and a
//! translation are exact, and a scale or shear is refused, never
//! approximated.
//!
//! [`ReferenceExactCompiler::compile_exact_with_report`] also says, per
//! body, whether any boolean on the way read a feature within the
//! tolerance (#236): an exact [`BooleanReport`] means the body is the exact
//! result of its operands as given; otherwise it was built from operands
//! moved by at most the report's magnitudes, within the tolerance.

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
                Ok(body)
            }
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
            10,
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
    "boolean",
    "bounded half-space",
];
