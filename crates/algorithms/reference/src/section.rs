//! Portable scalar plane section of a closed oriented triangle mesh.
//!
//! Topology is classified from the exact binary64 plane equation. Geometry is
//! emitted in the caller's plane frame, with source-edge connectivity preserved;
//! no tolerance-based point welding is used.
//!
//! Two entry points share the exact vertex classification and segment
//! building:
//!
//! - `section` is the legacy path. A face lying in the plane is refused as
//!   degenerate, and its messages are part of its contract.
//! - `section_request` (ADR 0091) reports faces lying in the plane as
//!   regions with a side, keeps cut curves to transversal crossings, and
//!   counts touching edges. `request` orchestrates it, `regions` traces the
//!   region boundaries and `chains` assembles open and closed chains.

mod chains;
mod cut_regions;
mod regions;
mod request;
mod ring_nesting;

use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, ExecutionOptions,
    ExecutionTarget, GeomError, GeomResult, ScratchRequirement, Sign,
};
use axiolid_core::{Frame3, Point2, Point3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::{
    MeshPlaneSection, PlaneSection, SectionEvidence, SectionLimits, SectionOutcome, SectionRequest,
};

use crate::orient3::orient3d;
use crate::orientation::orient2d;
use chains::assemble_contours;

/// Portable deterministic mesh plane-section oracle.
#[derive(Debug, Default, Clone, Copy)]
pub struct ScalarSection;

impl ScalarSection {
    /// Stable identity for this reference provider.
    pub const ID: BackendId = BackendId::new("scalar-section");

    /// Construct the provider.
    #[must_use]
    pub const fn new() -> Self {
        Self
    }
}

impl Backend for ScalarSection {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(Self::ID, ExecutionTarget::PortableCpu)
    }
}

impl MeshPlaneSection for ScalarSection {
    fn scratch_requirement(&self) -> ScratchRequirement {
        // Signs/distances, source-edge keys, segment adjacency, and cycle state.
        // A closed triangular manifold has O(triangles) vertices and edges.
        ScratchRequirement::PerElement {
            bytes_per_element: 512,
        }
    }

    fn request_scratch_requirement(&self) -> ScratchRequirement {
        // Per element, where an element is a position or a triangle: the
        // legacy bound (512); for a sheet, the edge-use map (3 entries of
        // about 40 bytes with tree overhead, 120) and the boundary edge and
        // vertex sets (at most 3 and 3 entries of about 24 bytes, 144); for
        // in-plane triangles, the directed-edge map (3 entries of about 40
        // bytes, 120), the triangles, boxes and flat copies (about 112), the
        // sweep order and component roots (16), the boundary set (72) and
        // the cycle and ring state (about 100). About 1.2 KiB in all,
        // rounded up to 1536 bytes.
        ScratchRequirement::PerElement {
            bytes_per_element: 1536,
        }
    }

    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::Incremental
    }

    fn section_request(
        &self,
        mesh: &TriMesh,
        request: &SectionRequest,
        options: &ExecutionOptions,
    ) -> GeomResult<PlaneSection> {
        request::plane_section(mesh, request, options)
    }

    fn section(
        &self,
        mesh: &TriMesh,
        frame: Frame3,
        limits: SectionLimits,
        options: &ExecutionOptions,
    ) -> GeomResult<SectionOutcome> {
        options.check_cancelled()?;
        check_source_limits(mesh, limits)?;

        let plane = ExactSectionPlane::new(frame)?;
        let classifications = classify_vertices(mesh, plane)?;

        let mut segments = BTreeSet::<Segment>::new();
        let mut on_plane_edges = BTreeMap::<EdgeKey, Vec<Sign>>::new();
        for (triangle_index, triangle) in mesh.triangles().enumerate() {
            options.check_cancelled()?;
            let vertices = [
                mesh_index(triangle[0])?,
                mesh_index(triangle[1])?,
                mesh_index(triangle[2])?,
            ];
            let signs = vertices.map(|index| classifications[index].sign);
            let zero_count = signs.iter().filter(|&&sign| sign == Sign::Zero).count();
            match zero_count {
                3 => {
                    return Err(GeomError::Degenerate(format!(
                        "section plane contains source triangle {triangle_index}; a two-dimensional overlap is not a curve"
                    )))
                }
                2 => {
                    let mut zeros = triangle
                        .into_iter()
                        .zip(signs)
                        .filter_map(|(index, sign)| (sign == Sign::Zero).then_some(index));
                    let left = zeros.next().ok_or_else(internal_topology_error)?;
                    let right = zeros.next().ok_or_else(internal_topology_error)?;
                    let third = signs
                        .into_iter()
                        .find(|&sign| sign != Sign::Zero)
                        .ok_or_else(internal_topology_error)?;
                    on_plane_edges
                        .entry(EdgeKey::new(left, right))
                        .or_default()
                        .push(third);
                }
                1 => {
                    let zero_corner = (0..3)
                        .find(|&corner| signs[corner] == Sign::Zero)
                        .ok_or_else(internal_topology_error)?;
                    let (first, second) = match zero_corner {
                        0 => (1, 2),
                        1 => (0, 2),
                        2 => (0, 1),
                        _ => return Err(internal_topology_error()),
                    };
                    if opposite(signs[first], signs[second]) {
                        let segment = Segment::new(
                            NodeKey::Vertex(triangle[zero_corner]),
                            NodeKey::Edge(EdgeKey::new(triangle[first], triangle[second])),
                        )?;
                        insert_segment(&mut segments, segment, limits)?;
                    }
                }
                0 => {
                    let mut crossing = [None, None];
                    let mut crossing_count = 0usize;
                    for (left, right) in [(0, 1), (1, 2), (2, 0)] {
                        if opposite(signs[left], signs[right]) {
                            if crossing_count >= crossing.len() {
                                return Err(internal_topology_error());
                            }
                            crossing[crossing_count] = Some(NodeKey::Edge(EdgeKey::new(
                                triangle[left],
                                triangle[right],
                            )));
                            crossing_count += 1;
                        }
                    }
                    match (crossing[0], crossing[1]) {
                        (Some(first), Some(second)) => insert_segment(
                            &mut segments,
                            Segment::new(first, second)?,
                            limits,
                        )?,
                        (None, None) => {}
                        _ => return Err(internal_topology_error()),
                    }
                }
                _ => return Err(internal_topology_error()),
            }
        }

        for (edge, incident_signs) in on_plane_edges {
            options.check_cancelled()?;
            if incident_signs.len() != 2 {
                return Err(GeomError::NotManifold(format!(
                    "on-plane mesh edge {:?} has {} incident triangles",
                    edge,
                    incident_signs.len()
                )));
            }
            if opposite(incident_signs[0], incident_signs[1]) {
                insert_segment(
                    &mut segments,
                    Segment::new(NodeKey::Vertex(edge.0), NodeKey::Vertex(edge.1))?,
                    limits,
                )?;
            }
        }

        let contours =
            assemble_contours(mesh, frame, &classifications, &segments, limits, options)?;
        let output_vertices = contours.iter().map(|contour| contour.points.len()).sum();
        let evidence =
            SectionEvidence::input_mesh(mesh.triangle_count(), output_vertices, contours.len());
        Ok(SectionOutcome::new(frame, contours, evidence))
    }
}

#[derive(Debug, Clone, Copy)]
struct Classification {
    sign: Sign,
    distance: f64,
}

#[derive(Debug, Clone, Copy)]
struct ExactSectionPlane {
    origin: Point3,
    x_point: Point3,
    y_point: Point3,
    normal: axiolid_core::Vec3,
}

impl ExactSectionPlane {
    fn new(frame: Frame3) -> GeomResult<Self> {
        let x_point = frame.origin + frame.x;
        let y_point = frame.origin + frame.y;
        if !x_point.is_finite()
            || !y_point.is_finite()
            || x_point == frame.origin
            || y_point == frame.origin
            || x_point == y_point
        {
            return Err(GeomError::Degenerate(
                "section frame cannot resolve a finite affine plane at this coordinate magnitude"
                    .into(),
            ));
        }
        let normal = (x_point - frame.origin).cross(y_point - frame.origin);
        let normal_length = normal.length();
        if !normal_length.is_finite() || normal_length == 0.0 {
            return Err(GeomError::Degenerate(
                "section affine plane has no finite normal".into(),
            ));
        }
        Ok(Self {
            origin: frame.origin,
            x_point,
            y_point,
            normal: normal / normal_length,
        })
    }

    fn classify(self, point: Point3) -> GeomResult<Classification> {
        let sign = match orient3d(self.origin, self.x_point, self.y_point, point) {
            axiolid_contracts::Certified::Certain { sign, .. } => sign,
            _ => {
                return Err(GeomError::Degenerate(
                    "certified plane-side predicate returned an uncertain sign".into(),
                ))
            }
        };
        let distance = self.normal.dot(point - self.origin);
        if !distance.is_finite() {
            return Err(GeomError::Degenerate(
                "section signed distance is not finite".into(),
            ));
        }
        if sign != Sign::Zero && distance == 0.0 {
            return Err(GeomError::Degenerate(
                "section interpolation lost a certified nonzero plane offset".into(),
            ));
        }
        Ok(Classification { sign, distance })
    }
}

/// Exact plane side and signed distance of every source position.
fn classify_vertices(mesh: &TriMesh, plane: ExactSectionPlane) -> GeomResult<Vec<Classification>> {
    let mut classifications = Vec::new();
    classifications
        .try_reserve_exact(mesh.positions.len())
        .map_err(|_| GeomError::BudgetExceeded { resource: "memory" })?;
    for &point in &mesh.positions {
        let classification = plane.classify(point)?;
        classifications.push(classification);
    }
    Ok(classifications)
}

fn opposite(left: Sign, right: Sign) -> bool {
    matches!(
        (left, right),
        (Sign::Negative, Sign::Positive) | (Sign::Positive, Sign::Negative)
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct EdgeKey(u32, u32);

impl EdgeKey {
    fn new(left: u32, right: u32) -> Self {
        Self(left.min(right), left.max(right))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum NodeKey {
    Vertex(u32),
    Edge(EdgeKey),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Segment(NodeKey, NodeKey);

impl Segment {
    fn new(left: NodeKey, right: NodeKey) -> GeomResult<Self> {
        if left == right {
            return Err(GeomError::Degenerate(
                "plane section collapsed a segment to one source-topology node".into(),
            ));
        }
        Ok(Self(left.min(right), left.max(right)))
    }
}

fn insert_segment(
    segments: &mut BTreeSet<Segment>,
    segment: Segment,
    limits: SectionLimits,
) -> GeomResult<()> {
    if segments.contains(&segment) {
        return Err(GeomError::NotManifold(
            "two source faces produced the same non-coplanar section segment".into(),
        ));
    }
    if segments.len() >= limits.max_output_vertices {
        return Err(GeomError::BudgetExceeded {
            resource: "section output vertices",
        });
    }
    segments.insert(segment);
    Ok(())
}

fn mesh_index(index: u32) -> GeomResult<usize> {
    usize::try_from(index)
        .map_err(|_| GeomError::InvalidInput("mesh index does not fit usize".into()))
}

fn node_point(
    mesh: &TriMesh,
    classifications: &[Classification],
    node: NodeKey,
) -> GeomResult<Point3> {
    match node {
        NodeKey::Vertex(index) => Ok(mesh.positions[mesh_index(index)?]),
        NodeKey::Edge(EdgeKey(left, right)) => {
            let a = mesh.positions[mesh_index(left)?];
            let b = mesh.positions[mesh_index(right)?];
            let da = classifications[mesh_index(left)?].distance.abs();
            let db = classifications[mesh_index(right)?].distance.abs();
            if !(da > 0.0 && db > 0.0 && da.is_finite() && db.is_finite()) {
                return Err(GeomError::Degenerate(
                    "crossing edge has no representable endpoint distance".into(),
                ));
            }
            let t = if da >= db {
                1.0 / (1.0 + db / da)
            } else {
                let ratio = da / db;
                ratio / (1.0 + ratio)
            };
            let point = a + (b - a) * t;
            if !point.is_finite() {
                return Err(GeomError::Degenerate(
                    "edge-plane intersection exceeded finite arithmetic".into(),
                ));
            }
            Ok(point)
        }
    }
}

fn simplify_collinear(points: &mut Vec<Point2>, tolerance: f64) {
    loop {
        if points.len() <= 3 {
            return;
        }
        let mut removed = false;
        for index in 0..points.len() {
            let previous = points[(index + points.len() - 1) % points.len()];
            let current = points[index];
            let next = points[(index + 1) % points.len()];
            let chord = next - previous;
            let scale = chord.length();
            let cross = (current - previous).perp_dot(chord).abs();
            let within_tolerance = scale > 0.0 && cross <= tolerance * scale;
            if orient2d(previous, current, next).sign() == Some(Sign::Zero) || within_tolerance {
                points.remove(index);
                removed = true;
                break;
            }
        }
        if !removed {
            return;
        }
    }
}

fn signed_area(points: &[Point2]) -> f64 {
    let mut twice = 0.0;
    for index in 0..points.len() {
        let current = points[index];
        let next = points[(index + 1) % points.len()];
        twice += current.x * next.y - current.y * next.x;
    }
    twice * 0.5
}

/// The exact sign of the signed area of the polygon through `points`: the sum
/// of `x_i y_{i+1} - x_{i+1} y_i` in expansion arithmetic, with no rounding.
/// Positive is counter-clockwise.
///
/// Products are exact only inside the binary64 range, so a nonzero coordinate
/// outside `2^-450 ..= 2^450` in magnitude (subnormal included) is refused:
/// its products could underflow or overflow, and a wrong sign is worse than a
/// refusal.
fn exact_area_sign(points: &[Point2]) -> GeomResult<Sign> {
    use axiolid_predicates::arithmetic::{expansion_sign, expansion_sum};
    use axiolid_predicates::two_product;

    let in_range = |value: f64| {
        let encoded = ((value.to_bits() >> 52) & 0x7ff) as i32;
        value == 0.0 || (encoded != 0 && (-450..=450).contains(&(encoded - 1023)))
    };
    if !points.iter().all(|p| in_range(p.x) && in_range(p.y)) {
        return Err(GeomError::Degenerate(
            "section chain coordinates are outside the exact area range".into(),
        ));
    }
    let mut twice = vec![0.0];
    for index in 0..points.len() {
        let (current, next) = (points[index], points[(index + 1) % points.len()]);
        if !(current.is_finite() && next.is_finite()) {
            return Err(GeomError::Degenerate(
                "section chain has a non-finite point".into(),
            ));
        }
        let (forward, forward_error) = two_product(current.x, next.y);
        let (backward, backward_error) = two_product(next.x, current.y);
        twice = expansion_sum(&twice, &[forward_error, forward]);
        twice = expansion_sum(&twice, &[-backward_error, -backward]);
    }
    Ok(expansion_sign(&twice))
}

fn point_order(left: &Point2, right: &Point2) -> std::cmp::Ordering {
    left.x
        .total_cmp(&right.x)
        .then_with(|| left.y.total_cmp(&right.y))
}

/// Lexicographic order of point sequences; a prefix sorts first.
fn sequence_order(left: &[Point2], right: &[Point2]) -> std::cmp::Ordering {
    left.iter()
        .zip(right)
        .map(|(a, b)| point_order(a, b))
        .find(|order| order.is_ne())
        .unwrap_or_else(|| left.len().cmp(&right.len()))
}

fn check_source_limits(mesh: &TriMesh, limits: SectionLimits) -> GeomResult<()> {
    if mesh.positions.len() > limits.max_source_vertices {
        return Err(GeomError::BudgetExceeded {
            resource: "section source vertices",
        });
    }
    if mesh.triangle_count() > limits.max_source_triangles {
        return Err(GeomError::BudgetExceeded {
            resource: "section source triangles",
        });
    }
    Ok(())
}

fn internal_topology_error() -> GeomError {
    GeomError::BackendContractViolation {
        backend: ScalarSection::ID,
        detail: "internal section topology state is inconsistent".into(),
    }
}

/// Plane-local coordinates of a model-space point. Chain points and region
/// vertices both go through this one function, so a chain end that lies on a
/// region boundary is bit-identical to the region vertex.
fn project_local(frame: Frame3, world: Point3) -> GeomResult<Point2> {
    let local = world - frame.origin;
    let point = Point2::new(local.dot(frame.x), local.dot(frame.y));
    if point.is_finite() {
        Ok(point)
    } else {
        Err(GeomError::Degenerate(
            "section projection exceeded finite arithmetic".into(),
        ))
    }
}

/// The exact sign of `orient2d`, refusing an uncertain predicate.
fn orient_sign(a: Point2, b: Point2, c: Point2) -> GeomResult<Sign> {
    orient2d(a, b, c).sign().ok_or_else(|| {
        GeomError::Degenerate("certified orientation predicate returned an uncertain sign".into())
    })
}

/// Product of two signs.
fn sign_times(left: Sign, right: Sign) -> Sign {
    if left == Sign::Zero || right == Sign::Zero {
        Sign::Zero
    } else if left == right {
        Sign::Positive
    } else {
        Sign::Negative
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_area_sign_is_exact_where_the_float_area_is_wrong() {
        // A thin triangle 1e8 from the origin: counter-clockwise, twice its
        // area about 0.41, while the float sum of cross products is -2.
        let points = [
            Point2::new(100000002.228125, 100000001.5671875),
            Point2::new(100000001.9375, 100000001.275),
            Point2::new(100000002.615625, 100000000.5578125),
        ];
        assert!(signed_area(&points) < 0.0);
        assert_eq!(exact_area_sign(&points).expect("sign"), Sign::Positive);
        let reversed = [points[2], points[1], points[0]];
        assert_eq!(exact_area_sign(&reversed).expect("sign"), Sign::Negative);
        let flat = [points[0], points[0], points[1]];
        assert_eq!(exact_area_sign(&flat).expect("sign"), Sign::Zero);
    }

    #[test]
    fn the_area_sign_refuses_coordinates_whose_products_leave_the_range() {
        // Products near 2^-1060 underflow: refused, not a wrong sign.
        let tiny = [
            Point2::new(-8.036185191855754e-160, -8.029313620019926e-160),
            Point2::new(-8.036193993089111e-160, -8.028301829579562e-160),
            Point2::new(-8.019240586512659e-160, -8.033664284419165e-160),
        ];
        assert!(matches!(
            exact_area_sign(&tiny),
            Err(GeomError::Degenerate(_))
        ));
        // Products near 1e400 overflow.
        let huge = [
            Point2::new(1e200, 1e200),
            Point2::new(2e200, 1e200),
            Point2::new(1e200, 3e200),
        ];
        assert!(matches!(
            exact_area_sign(&huge),
            Err(GeomError::Degenerate(_))
        ));
        // Large but inside the range stays exact.
        let large = [
            Point2::new(1e100, 1e100),
            Point2::new(2e100, 1e100),
            Point2::new(1e100, 3e100),
        ];
        assert_eq!(exact_area_sign(&large).expect("sign"), Sign::Positive);
    }

    #[test]
    fn exact_plane_sign_keeps_a_tiny_nonzero_binary64_offset() {
        let frame = Frame3 {
            origin: Point3::ZERO,
            x: Point3::X,
            y: Point3::Y,
            z: Point3::Z,
        };
        let plane = ExactSectionPlane::new(frame).expect("resolvable plane");
        let point = Point3::new(1.0, -1.0, f64::from_bits(1));
        let certified = orient3d(plane.origin, plane.x_point, plane.y_point, point);
        assert!(matches!(
            certified,
            axiolid_contracts::Certified::Certain {
                precision: axiolid_contracts::Precision::Exact,
                ..
            }
        ));
        let positive = plane.classify(point).expect("finite subnormal");
        assert_ne!(positive.sign, Sign::Zero);
    }
}
