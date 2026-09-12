//! Certified bounded intersections between clamped NURBS surfaces.
//!
//! The complete trace certificate implemented here is deliberately narrow: both
//! inputs must be single-span polynomial affine patches with continuous clamped
//! axes. General patch pairs are still bounded and conservatively returned as
//! unresolved candidates rather than passed through heuristic marching.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar};
use axiolid_curve::{BSplineCurve3, KnotSpec};
use axiolid_surface::BSplineSurface;

use crate::{
    certified_bezier::Interval,
    certified_curve_surface_intersection::{
        intersect_curve_surface_certified, CertifiedCurveSurfaceIntersection3,
        CertifiedCurveSurfaceIntersectionOptions, TransverseCurveSurfaceIntersection3,
    },
    certified_projection::ParameterInterval,
    certified_refinement::RefinementBudget,
    certified_surface_bezier::{piecewise_bezier_patches, Patch},
};

const MAX_REFINEMENT_WORK: u32 = 100_000;
const MAX_BOUNDARY_NODES: u32 = 100_000;
const MAX_DEPTH: u16 = 64;
const BOUNDARY_QUERY_COUNT: u8 = 8;

/// Bounded policy for certified surface/surface intersection and affine tracing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CertifiedSurfaceSurfaceIntersectionOptions {
    parameter_tolerance: Scalar,
    max_refinement_work: u32,
    max_boundary_nodes: u32,
    max_depth: u16,
}

impl CertifiedSurfaceSurfaceIntersectionOptions {
    /// Construct a finite policy.
    ///
    /// `max_refinement_work` is shared by both tensor refinements and the initial
    /// patch-pair product. At most eight boundary curve/surface queries are run,
    /// each independently capped by `max_boundary_nodes` and `max_depth`.
    pub fn new(
        parameter_tolerance: Scalar,
        max_refinement_work: u32,
        max_boundary_nodes: u32,
        max_depth: u16,
    ) -> GeomResult<Self> {
        if !parameter_tolerance.is_finite() || parameter_tolerance <= 0.0 {
            return Err(GeomError::InvalidInput(
                "surface/surface parameter tolerance must be finite and positive".to_owned(),
            ));
        }
        if max_refinement_work == 0 || max_refinement_work > MAX_REFINEMENT_WORK {
            return Err(GeomError::InvalidInput(format!(
                "surface/surface max_refinement_work must be in 1..={MAX_REFINEMENT_WORK}"
            )));
        }
        if max_boundary_nodes == 0 || max_boundary_nodes > MAX_BOUNDARY_NODES {
            return Err(GeomError::InvalidInput(format!(
                "surface/surface max_boundary_nodes must be in 1..={MAX_BOUNDARY_NODES}"
            )));
        }
        if max_depth == 0 || max_depth > MAX_DEPTH {
            return Err(GeomError::InvalidInput(format!(
                "surface/surface max_depth must be in 1..={MAX_DEPTH}"
            )));
        }
        Ok(Self {
            parameter_tolerance,
            max_refinement_work,
            max_boundary_nodes,
            max_depth,
        })
    }
}

impl Default for CertifiedSurfaceSurfaceIntersectionOptions {
    fn default() -> Self {
        Self {
            parameter_tolerance: 1.0e-8,
            max_refinement_work: MAX_REFINEMENT_WORK,
            max_boundary_nodes: MAX_BOUNDARY_NODES,
            max_depth: MAX_DEPTH,
        }
    }
}

/// Native four-parameter box for a pair of tensor-product surfaces.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct SurfaceSurfaceParameterBox {
    /// First-surface `u` enclosure.
    pub first_u: ParameterInterval,
    /// First-surface `v` enclosure.
    pub first_v: ParameterInterval,
    /// Second-surface `u` enclosure.
    pub second_u: ParameterInterval,
    /// Second-surface `v` enclosure.
    pub second_v: ParameterInterval,
}

/// Certified endpoint enclosure for an affine surface/surface trace.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct SurfaceSurfaceTraceEndpoint3 {
    /// Native parameter enclosure on both surfaces.
    pub parameters: SurfaceSurfaceParameterBox,
    /// Representative midpoint of the two evaluated images.
    pub point: Point3,
    /// Conservative endpoint residual upper bound.
    pub residual_upper_bound: Scalar,
}

/// Certificate for one complete regular affine intersection segment.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct TransverseSurfaceSurfaceTrace3 {
    /// Lexicographically ordered first endpoint certificate.
    pub start: SurfaceSurfaceTraceEndpoint3,
    /// Lexicographically ordered second endpoint certificate.
    pub end: SurfaceSurfaceTraceEndpoint3,
    /// Positive interval lower bound for `|(S1_u x S1_v) x (S2_u x S2_v)|^2`.
    pub normal_cross_squared_lower_bound: Scalar,
}

/// Certified surface/surface query outcome.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum CertifiedSurfaceSurfaceIntersection3 {
    /// Every patch pair was excluded or represented by a complete certified trace.
    Complete {
        /// Certified regular affine trace segments.
        traces: Vec<TransverseSurfaceSurfaceTrace3>,
        /// Number of initial patch pairs visited by the broad phase.
        visited_patch_pairs: u32,
        /// Number of certified boundary queries executed.
        boundary_queries: u8,
    },
    /// One or more conservative patch-pair candidates remain unresolved.
    Unresolved {
        /// Trace segments proved independently of unresolved candidates.
        traces: Vec<TransverseSurfaceSurfaceTrace3>,
        /// Conservative native boxes that may contain untraced intersections.
        candidate_boxes: Vec<SurfaceSurfaceParameterBox>,
        /// Number of initial patch pairs visited by the broad phase.
        visited_patch_pairs: u32,
        /// Number of certified boundary queries executed.
        boundary_queries: u8,
    },
}

/// Certify bounded intersections of two clamped, internally continuous NURBS surfaces.
///
/// Coordinate-hull exclusion is available for general polynomial and positive-weight
/// rational inputs. Complete traced segments currently require single-span polynomial
/// affine patches. Non-affine, multispan, coincident, tangential, corner-owned, and
/// proof-insufficient candidates remain [`Unresolved`](CertifiedSurfaceSurfaceIntersection3::Unresolved).
pub fn intersect_surface_surface_certified(
    first: &BSplineSurface,
    second: &BSplineSurface,
    options: CertifiedSurfaceSurfaceIntersectionOptions,
) -> GeomResult<CertifiedSurfaceSurfaceIntersection3> {
    let options = CertifiedSurfaceSurfaceIntersectionOptions::new(
        options.parameter_tolerance,
        options.max_refinement_work,
        options.max_boundary_nodes,
        options.max_depth,
    )?;
    let mut budget = RefinementBudget::new(
        options.max_refinement_work,
        "certified surface/surface refinement budget",
    );
    let first_patches = piecewise_bezier_patches(first, &mut budget)?;
    let second_patches = piecewise_bezier_patches(second, &mut budget)?;
    let pair_count = first_patches
        .len()
        .checked_mul(second_patches.len())
        .ok_or(GeomError::BudgetExceeded {
            resource: "certified surface/surface refinement budget",
        })?;
    budget.charge(u128::try_from(pair_count).ok())?;
    let visited_patch_pairs = u32::try_from(pair_count).map_err(|_| GeomError::BudgetExceeded {
        resource: "certified surface/surface refinement budget",
    })?;
    let mut candidates = Vec::new();
    candidates
        .try_reserve_exact(pair_count)
        .map_err(|_| allocation_error("certified surface/surface candidate allocation"))?;
    for first_patch in &first_patches {
        for second_patch in &second_patches {
            if !patches_are_disjoint(first_patch, second_patch)? {
                candidates.push(parameter_box(first_patch, second_patch));
            }
        }
    }
    if candidates.is_empty() {
        return Ok(CertifiedSurfaceSurfaceIntersection3::Complete {
            traces: Vec::new(),
            visited_patch_pairs,
            boundary_queries: 0,
        });
    }
    if candidates.len() != 1
        || first_patches.len() != 1
        || second_patches.len() != 1
        || !is_exact_single_span_affine(first)
        || !is_exact_single_span_affine(second)
    {
        return unresolved(candidates, visited_patch_pairs, 0);
    }
    let normal_lower = normal_cross_squared_lower_bound(&first_patches[0], &second_patches[0])?;
    if normal_lower <= 0.0 {
        return unresolved(candidates, visited_patch_pairs, 0);
    }
    match trace_affine_pair(first, second, options, normal_lower)? {
        AffineTraceOutcome::Complete(traces) => {
            Ok(CertifiedSurfaceSurfaceIntersection3::Complete {
                traces,
                visited_patch_pairs,
                boundary_queries: BOUNDARY_QUERY_COUNT,
            })
        }
        AffineTraceOutcome::Unresolved(traces) => {
            Ok(CertifiedSurfaceSurfaceIntersection3::Unresolved {
                traces,
                candidate_boxes: candidates,
                visited_patch_pairs,
                boundary_queries: BOUNDARY_QUERY_COUNT,
            })
        }
    }
}

fn unresolved(
    candidate_boxes: Vec<SurfaceSurfaceParameterBox>,
    visited_patch_pairs: u32,
    boundary_queries: u8,
) -> GeomResult<CertifiedSurfaceSurfaceIntersection3> {
    Ok(CertifiedSurfaceSurfaceIntersection3::Unresolved {
        traces: Vec::new(),
        candidate_boxes,
        visited_patch_pairs,
        boundary_queries,
    })
}

fn allocation_error(resource: &'static str) -> GeomError {
    GeomError::BudgetExceeded { resource }
}

fn parameter_box(first: &Patch, second: &Patch) -> SurfaceSurfaceParameterBox {
    SurfaceSurfaceParameterBox {
        first_u: ParameterInterval {
            start: first.u_start,
            end: first.u_end,
        },
        first_v: ParameterInterval {
            start: first.v_start,
            end: first.v_end,
        },
        second_u: ParameterInterval {
            start: second.u_start,
            end: second.u_end,
        },
        second_v: ParameterInterval {
            start: second.v_start,
            end: second.v_end,
        },
    }
}

fn patches_are_disjoint(first: &Patch, second: &Patch) -> GeomResult<bool> {
    let first_bounds = first.coordinate_intervals()?;
    let second_bounds = second.coordinate_intervals()?;
    Ok((0..3).any(|axis| {
        first_bounds[axis].upper() < second_bounds[axis].lower()
            || second_bounds[axis].upper() < first_bounds[axis].lower()
    }))
}

fn is_exact_single_span_affine(surface: &BSplineSurface) -> bool {
    if surface.u_degree != 1
        || surface.v_degree != 1
        || surface.control_points.len() != 2
        || surface.control_points.iter().any(|row| row.len() != 2)
        || surface.weights.is_some()
        || surface.u_knots.len() != 2
        || surface.v_knots.len() != 2
    {
        return false;
    }
    let p00 = surface.control_points[0][0];
    let p01 = surface.control_points[0][1];
    let p10 = surface.control_points[1][0];
    let p11 = surface.control_points[1][1];
    [
        [p11.x, -p10.x, -p01.x, p00.x],
        [p11.y, -p10.y, -p01.y, p00.y],
        [p11.z, -p10.z, -p01.z, p00.z],
    ]
    .into_iter()
    .all(exact_sum_is_zero)
}

// Shewchuk-style nonoverlapping expansion addition. Each binary64 input is an
// exact dyadic rational, so an all-zero final expansion proves the affine
// cross-term identity without a tolerance or rounded cancellation decision.
fn exact_sum_is_zero(values: [Scalar; 4]) -> bool {
    let mut expansion = [0.0; 8];
    let mut length = 1usize;
    expansion[0] = values[0];
    for value in values.into_iter().skip(1) {
        let mut next = [0.0; 8];
        let mut next_length = 0usize;
        let mut accumulator = value;
        for component in expansion.iter().take(length).copied() {
            let (sum, error) = two_sum(accumulator, component);
            if error != 0.0 {
                next[next_length] = error;
                next_length += 1;
            }
            accumulator = sum;
        }
        if accumulator != 0.0 || next_length == 0 {
            next[next_length] = accumulator;
            next_length += 1;
        }
        expansion = next;
        length = next_length;
    }
    expansion
        .iter()
        .take(length)
        .all(|component| *component == 0.0)
}

fn two_sum(left: Scalar, right: Scalar) -> (Scalar, Scalar) {
    let sum = left + right;
    let right_virtual = sum - left;
    let left_virtual = sum - right_virtual;
    let right_roundoff = right - right_virtual;
    let left_roundoff = left - left_virtual;
    (sum, left_roundoff + right_roundoff)
}

fn normal_cross_squared_lower_bound(first: &Patch, second: &Patch) -> GeomResult<Scalar> {
    let first_normal = cross_intervals(first.partial_u_intervals()?, first.partial_v_intervals()?)?;
    let second_normal =
        cross_intervals(second.partial_u_intervals()?, second.partial_v_intervals()?)?;
    let cross = cross_intervals(first_normal, second_normal)?;
    let mut squared = Interval::exact(0.0)?;
    for component in cross {
        let lower = Interval::exact(component.absolute_lower_bound())?;
        squared = squared.add(lower.multiply(lower)?)?;
    }
    Ok(squared.lower().max(0.0))
}

fn cross_intervals(left: [Interval; 3], right: [Interval; 3]) -> GeomResult<[Interval; 3]> {
    Ok([
        left[1]
            .multiply(right[2])?
            .subtract(left[2].multiply(right[1])?)?,
        left[2]
            .multiply(right[0])?
            .subtract(left[0].multiply(right[2])?)?,
        left[0]
            .multiply(right[1])?
            .subtract(left[1].multiply(right[0])?)?,
    ])
}

#[derive(Debug, Clone, Copy)]
enum Boundary {
    UStart,
    UEnd,
    VStart,
    VEnd,
}

const BOUNDARIES: [Boundary; 4] = [
    Boundary::UStart,
    Boundary::UEnd,
    Boundary::VStart,
    Boundary::VEnd,
];

enum AffineTraceOutcome {
    Complete(Vec<TransverseSurfaceSurfaceTrace3>),
    Unresolved(Vec<TransverseSurfaceSurfaceTrace3>),
}

fn trace_affine_pair(
    first: &BSplineSurface,
    second: &BSplineSurface,
    options: CertifiedSurfaceSurfaceIntersectionOptions,
    normal_lower: Scalar,
) -> GeomResult<AffineTraceOutcome> {
    let boundary_options = CertifiedCurveSurfaceIntersectionOptions::new(
        options.parameter_tolerance,
        options.max_boundary_nodes,
        options.max_depth,
    )?;
    let mut endpoints = Vec::new();
    endpoints
        .try_reserve_exact(BOUNDARY_QUERY_COUNT.into())
        .map_err(|_| allocation_error("certified surface/surface endpoint allocation"))?;
    let mut any_unresolved = false;
    for boundary in BOUNDARIES {
        if let Err(error) = collect_boundary_roots(
            first,
            second,
            boundary,
            true,
            boundary_options,
            &mut endpoints,
            &mut any_unresolved,
        ) {
            if matches!(error, GeomError::BudgetExceeded { .. }) {
                any_unresolved = true;
                continue;
            }
            return Err(error);
        }
    }
    for boundary in BOUNDARIES {
        if let Err(error) = collect_boundary_roots(
            second,
            first,
            boundary,
            false,
            boundary_options,
            &mut endpoints,
            &mut any_unresolved,
        ) {
            if matches!(error, GeomError::BudgetExceeded { .. }) {
                any_unresolved = true;
                continue;
            }
            return Err(error);
        }
    }
    // A chord ending ON a boundary of BOTH patches is discovered twice:
    // once scanning the first surface's boundaries, once the second's.
    // Both reports denote one point, so fold them before counting.
    // Step 1 made these boundary roots certifiable, which is what
    // exposed the duplication -- previously they were never found.
    dedupe_coincident_endpoints(&mut endpoints);
    if any_unresolved || endpoints.len() == 1 || endpoints.len() > 2 {
        return Ok(AffineTraceOutcome::Unresolved(Vec::new()));
    }
    if endpoints.is_empty() {
        return Ok(AffineTraceOutcome::Complete(Vec::new()));
    }
    endpoints.sort_by(|left, right| lexicographic_point_cmp(left.point, right.point));
    if endpoint_boxes_overlap(&endpoints[0], &endpoints[1]) {
        return Ok(AffineTraceOutcome::Unresolved(Vec::new()));
    }
    let mut traces = Vec::new();
    traces
        .try_reserve_exact(1)
        .map_err(|_| allocation_error("certified surface/surface trace allocation"))?;
    traces.push(TransverseSurfaceSurfaceTrace3 {
        start: endpoints.remove(0),
        end: endpoints.remove(0),
        normal_cross_squared_lower_bound: normal_lower,
    });
    Ok(AffineTraceOutcome::Complete(traces))
}

fn collect_boundary_roots(
    owner: &BSplineSurface,
    other: &BSplineSurface,
    boundary: Boundary,
    owner_is_first: bool,
    options: CertifiedCurveSurfaceIntersectionOptions,
    endpoints: &mut Vec<SurfaceSurfaceTraceEndpoint3>,
    any_unresolved: &mut bool,
) -> GeomResult<()> {
    let curve = boundary_curve(owner, boundary)?;
    match intersect_curve_surface_certified(&curve, other, options)? {
        CertifiedCurveSurfaceIntersection3::Complete { intersections, .. } => {
            for intersection in intersections {
                push_endpoint(
                    endpoints,
                    map_endpoint(owner, boundary, owner_is_first, intersection),
                )?;
            }
        }
        CertifiedCurveSurfaceIntersection3::Unresolved { .. } => *any_unresolved = true,
    }
    Ok(())
}

fn push_endpoint(
    endpoints: &mut Vec<SurfaceSurfaceTraceEndpoint3>,
    endpoint: SurfaceSurfaceTraceEndpoint3,
) -> GeomResult<()> {
    if endpoints.len() == endpoints.capacity() {
        endpoints
            .try_reserve_exact(1)
            .map_err(|_| allocation_error("certified surface/surface endpoint allocation"))?;
    }
    endpoints.push(endpoint);
    Ok(())
}

fn boundary_curve(surface: &BSplineSurface, boundary: Boundary) -> GeomResult<BSplineCurve3> {
    let (degree, knots, multiplicities, controls) = match boundary {
        Boundary::UStart => (
            surface.v_degree,
            &surface.v_knots,
            &surface.v_multiplicities,
            clone_points(surface.control_points[0].iter().copied())?,
        ),
        Boundary::UEnd => (
            surface.v_degree,
            &surface.v_knots,
            &surface.v_multiplicities,
            clone_points(surface.control_points[1].iter().copied())?,
        ),
        Boundary::VStart => (
            surface.u_degree,
            &surface.u_knots,
            &surface.u_multiplicities,
            clone_points(surface.control_points.iter().map(|row| row[0]))?,
        ),
        Boundary::VEnd => (
            surface.u_degree,
            &surface.u_knots,
            &surface.u_multiplicities,
            clone_points(surface.control_points.iter().map(|row| row[1]))?,
        ),
    };
    Ok(BSplineCurve3 {
        degree,
        control_points: controls,
        knots: clone_scalars(knots)?,
        multiplicities: clone_multiplicities(multiplicities)?,
        weights: None,
        knot_spec: KnotSpec::Unspecified,
        closed: false,
        self_intersect: None,
    })
}

fn map_endpoint(
    owner: &BSplineSurface,
    boundary: Boundary,
    owner_is_first: bool,
    root: TransverseCurveSurfaceIntersection3,
) -> SurfaceSurfaceTraceEndpoint3 {
    let fixed_u = |value: Scalar| ParameterInterval {
        start: value,
        end: value,
    };
    let (owner_u, owner_v) = match boundary {
        Boundary::UStart => (fixed_u(owner.u_knots[0]), root.curve_parameter),
        Boundary::UEnd => (
            fixed_u(owner.u_knots[owner.u_knots.len() - 1]),
            root.curve_parameter,
        ),
        Boundary::VStart => (root.curve_parameter, fixed_u(owner.v_knots[0])),
        Boundary::VEnd => (
            root.curve_parameter,
            fixed_u(owner.v_knots[owner.v_knots.len() - 1]),
        ),
    };
    let other_u = root.surface_u_parameter;
    let other_v = root.surface_v_parameter;
    let parameters = if owner_is_first {
        SurfaceSurfaceParameterBox {
            first_u: owner_u,
            first_v: owner_v,
            second_u: other_u,
            second_v: other_v,
        }
    } else {
        SurfaceSurfaceParameterBox {
            first_u: other_u,
            first_v: other_v,
            second_u: owner_u,
            second_v: owner_v,
        }
    };
    SurfaceSurfaceTraceEndpoint3 {
        parameters,
        point: root.point,
        residual_upper_bound: root.residual_upper_bound,
    }
}

/// Fold endpoint reports that denote the same intersection point.
///
/// Two reports are the same point when their certified parameter boxes
/// overlap on ALL FOUR axes -- both surfaces' u and v. That is the same
/// test the trace builder already uses to reject a degenerate chord, so
/// coincidence is decided one way in this file.
///
/// Conservative by construction: boxes that merely touch are folded, and
/// two genuinely distinct roots cannot have overlapping boxes on every
/// axis without violating the certificate that isolated each of them.
fn dedupe_coincident_endpoints(endpoints: &mut Vec<SurfaceSurfaceTraceEndpoint3>) {
    let mut kept: Vec<SurfaceSurfaceTraceEndpoint3> = Vec::new();
    for endpoint in endpoints.drain(..) {
        if kept
            .iter()
            .any(|existing| endpoint_boxes_overlap(existing, &endpoint))
        {
            continue;
        }
        kept.push(endpoint);
    }
    *endpoints = kept;
}

fn endpoint_boxes_overlap(
    first: &SurfaceSurfaceTraceEndpoint3,
    second: &SurfaceSurfaceTraceEndpoint3,
) -> bool {
    let left = first.parameters;
    let right = second.parameters;
    intervals_overlap(left.first_u, right.first_u)
        && intervals_overlap(left.first_v, right.first_v)
        && intervals_overlap(left.second_u, right.second_u)
        && intervals_overlap(left.second_v, right.second_v)
}

fn intervals_overlap(first: ParameterInterval, second: ParameterInterval) -> bool {
    first.start <= second.end && second.start <= first.end
}

fn lexicographic_point_cmp(left: Point3, right: Point3) -> std::cmp::Ordering {
    left.x
        .total_cmp(&right.x)
        .then_with(|| left.y.total_cmp(&right.y))
        .then_with(|| left.z.total_cmp(&right.z))
}

fn clone_points(values: impl ExactSizeIterator<Item = Point3>) -> GeomResult<Vec<Point3>> {
    let count = values.len();
    let mut output = Vec::new();
    output
        .try_reserve_exact(count)
        .map_err(|_| allocation_error("certified surface/surface boundary allocation"))?;
    for value in values {
        output.push(value);
    }
    Ok(output)
}

fn clone_scalars(values: &[Scalar]) -> GeomResult<Vec<Scalar>> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(values.len())
        .map_err(|_| allocation_error("certified surface/surface boundary allocation"))?;
    output.extend(values.iter().copied());
    Ok(output)
}

fn clone_multiplicities(values: &[u32]) -> GeomResult<Vec<u32>> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(values.len())
        .map_err(|_| allocation_error("certified surface/surface boundary allocation"))?;
    output.extend(values.iter().copied());
    Ok(output)
}
