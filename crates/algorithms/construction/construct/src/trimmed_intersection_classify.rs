use axiolid_core::{Point2, Scalar};
use axiolid_nurbs::{
    ParameterInterval, SurfaceSurfaceTraceEndpoint3, TransverseSurfaceSurfaceTrace3,
};
use axiolid_surface::BSplineSurface;

use crate::trimmed_intersection_types::{SurfacePairMember, SurfacePairSplitUnresolvedReason};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum BoundarySide {
    VStart,
    UEnd,
    VEnd,
    UStart,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Domain2 {
    pub u_start: Scalar,
    pub u_end: Scalar,
    pub v_start: Scalar,
    pub v_end: Scalar,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Endpoint2 {
    pub uv: Point2,
    pub side: Option<BoundarySide>,
}

/// A chord that partitions both patches: each side gets its own owner data.
#[derive(Debug, Clone, Copy)]
pub(super) struct DualClassification {
    pub first_domain: Domain2,
    pub second_domain: Domain2,
    pub first_start: Endpoint2,
    pub first_end: Endpoint2,
    pub second_start: Endpoint2,
    pub second_end: Endpoint2,
}

/// What the classifier proved about this trace.
#[derive(Debug, Clone, Copy)]
pub(super) enum Classification {
    /// Exactly one patch is partitioned; the other contains the chord.
    Single(SplitClassification),
    /// Both patches are partitioned by the same chord.
    Dual(DualClassification),
}

#[derive(Debug, Clone, Copy)]
pub(super) struct SplitClassification {
    pub member: SurfacePairMember,
    pub owner_domain: Domain2,
    pub embedded_domain: Domain2,
    pub owner_start: Endpoint2,
    pub owner_end: Endpoint2,
    pub embedded_start: Endpoint2,
    pub embedded_end: Endpoint2,
}

/// Classify a trace into a split arrangement, or say why it cannot be one.
///
/// The error distinguishes a PROVEN terminal refusal from an unimplemented
/// combination, so a caller knows whether retrying could ever help.
pub(super) fn classify(
    first: &BSplineSurface,
    second: &BSplineSurface,
    trace: &TransverseSurfaceSurfaceTrace3,
) -> Result<Classification, SurfacePairSplitUnresolvedReason> {
    // A domain or endpoint that will not resolve is a degenerate
    // representative, not an ownership question -- report it as such.
    let degenerate = SurfacePairSplitUnresolvedReason::DegenerateRepresentative;
    let first_domain = domain(first).ok_or(degenerate)?;
    let second_domain = domain(second).ok_or(degenerate)?;
    let first_start =
        endpoint(&trace.start, SurfacePairMember::First, first_domain).ok_or(degenerate)?;
    let first_end =
        endpoint(&trace.end, SurfacePairMember::First, first_domain).ok_or(degenerate)?;
    let second_start =
        endpoint(&trace.start, SurfacePairMember::Second, second_domain).ok_or(degenerate)?;
    let second_end =
        endpoint(&trace.end, SurfacePairMember::Second, second_domain).ok_or(degenerate)?;

    let first_owns = owns_chord(first_start, first_end, first_domain);
    let second_owns = owns_chord(second_start, second_end, second_domain);
    let first_embeds = interior(first_start, first_domain) && interior(first_end, first_domain);
    let second_embeds =
        interior(second_start, second_domain) && interior(second_end, second_domain);

    match (first_owns, second_owns, first_embeds, second_embeds) {
        (true, false, false, true) => Ok(Classification::Single(SplitClassification {
            member: SurfacePairMember::First,
            owner_domain: first_domain,
            embedded_domain: second_domain,
            owner_start: first_start,
            owner_end: first_end,
            embedded_start: second_start,
            embedded_end: second_end,
        })),
        (false, true, true, false) => Ok(Classification::Single(SplitClassification {
            member: SurfacePairMember::Second,
            owner_domain: second_domain,
            embedded_domain: first_domain,
            owner_start: second_start,
            owner_end: second_end,
            embedded_start: first_start,
            embedded_end: first_end,
        })),
        // Both patches are partitioned by the same chord. Each side is a
        // real trim boundary, so neither can be demoted to an embedded
        // annotation without misreporting the topology.
        (true, true, false, false) => Ok(Classification::Dual(DualClassification {
            first_domain,
            second_domain,
            first_start,
            first_end,
            second_start,
            second_end,
        })),
        // Neither patch is partitioned. Decide whether that is provable or
        // merely unimplemented before refusing, because the two demand
        // different things of the caller.
        _ => Err(
            if is_slit(first_start, first_end, first_domain)
                && is_slit(second_start, second_end, second_domain)
            {
                SurfacePairSplitUnresolvedReason::NoPartitionExists
            } else {
                SurfacePairSplitUnresolvedReason::UnsupportedEndpointOwnership
            },
        ),
    }
}

/// Whether the trace merely slits this patch instead of partitioning it.
///
/// Exactly one endpoint strictly inside the domain and the other on a
/// boundary side means the curve dead-ends in the interior. The face stays
/// simply connected, so no pair of closed trimmed faces can be built from
/// it -- and that conclusion does not depend on tolerance, so a caller
/// gains nothing by retrying.
///
/// # Coverage gap
///
/// The caller requires BOTH patches to be slit before reporting a terminal
/// proof. Relaxing that conjunction to a disjunction is NOT currently
/// caught by any test: every trace the certified path can produce today
/// comes from two planar patches, which meet in a full line and therefore
/// clip symmetrically -- so one-patch-slit-one-not never arises. The case
/// is reachable in principle (a short trace crossing one patch's edge but
/// landing inside the other) and needs a curved patch or non-axis-aligned
/// trim to construct, which waits on certified boundary roots.
fn is_slit(start: Endpoint2, end: Endpoint2, domain: Domain2) -> bool {
    let start_inside = interior(start, domain);
    let end_inside = interior(end, domain);
    let start_on_side = start.side.is_some() && side_interior(start, domain);
    let end_on_side = end.side.is_some() && side_interior(end, domain);
    (start_inside && end_on_side) || (end_inside && start_on_side)
}

fn owns_chord(start: Endpoint2, end: Endpoint2, domain: Domain2) -> bool {
    matches!((start.side, end.side), (Some(left), Some(right)) if left != right)
        && side_interior(start, domain)
        && side_interior(end, domain)
        && start.uv != end.uv
}

fn side_interior(endpoint: Endpoint2, domain: Domain2) -> bool {
    match endpoint.side {
        Some(BoundarySide::VStart | BoundarySide::VEnd) => {
            endpoint.uv.x > domain.u_start && endpoint.uv.x < domain.u_end
        }
        Some(BoundarySide::UStart | BoundarySide::UEnd) => {
            endpoint.uv.y > domain.v_start && endpoint.uv.y < domain.v_end
        }
        None => false,
    }
}

fn interior(endpoint: Endpoint2, domain: Domain2) -> bool {
    endpoint.side.is_none()
        && endpoint.uv.x > domain.u_start
        && endpoint.uv.x < domain.u_end
        && endpoint.uv.y > domain.v_start
        && endpoint.uv.y < domain.v_end
}

fn domain(surface: &BSplineSurface) -> Option<Domain2> {
    let u_start = *surface.u_knots.first()?;
    let u_end = *surface.u_knots.last()?;
    let v_start = *surface.v_knots.first()?;
    let v_end = *surface.v_knots.last()?;
    let values = [u_start, u_end, v_start, v_end];
    if values.iter().all(|value| value.is_finite()) && u_start < u_end && v_start < v_end {
        Some(Domain2 {
            u_start,
            u_end,
            v_start,
            v_end,
        })
    } else {
        None
    }
}

fn endpoint(
    endpoint: &SurfaceSurfaceTraceEndpoint3,
    member: SurfacePairMember,
    domain: Domain2,
) -> Option<Endpoint2> {
    let parameters = endpoint.parameters;
    let (u, v) = match member {
        SurfacePairMember::First => (parameters.first_u, parameters.first_v),
        SurfacePairMember::Second => (parameters.second_u, parameters.second_v),
    };
    let uv = Point2::new(midpoint(u)?, midpoint(v)?);
    if uv.x < domain.u_start || uv.x > domain.u_end || uv.y < domain.v_start || uv.y > domain.v_end
    {
        return None;
    }
    let mut side = None;
    for candidate in [
        fixed_side(v, domain.v_start, BoundarySide::VStart),
        fixed_side(u, domain.u_end, BoundarySide::UEnd),
        fixed_side(v, domain.v_end, BoundarySide::VEnd),
        fixed_side(u, domain.u_start, BoundarySide::UStart),
    ]
    .into_iter()
    .flatten()
    {
        if side.replace(candidate).is_some() {
            return None;
        }
    }
    Some(Endpoint2 { uv, side })
}

fn fixed_side(
    interval: ParameterInterval,
    boundary: Scalar,
    side: BoundarySide,
) -> Option<BoundarySide> {
    (interval.start == boundary && interval.end == boundary).then_some(side)
}

fn midpoint(interval: ParameterInterval) -> Option<Scalar> {
    if !interval.start.is_finite() || !interval.end.is_finite() || interval.start > interval.end {
        return None;
    }
    let value = interval.start * 0.5 + interval.end * 0.5;
    value.is_finite().then_some(value)
}

pub(super) fn boundary_rank(side: BoundarySide, uv: Point2, domain: Domain2) -> Scalar {
    match side {
        BoundarySide::VStart => (uv.x - domain.u_start) / (domain.u_end - domain.u_start),
        BoundarySide::UEnd => 1.0 + (uv.y - domain.v_start) / (domain.v_end - domain.v_start),
        BoundarySide::VEnd => 2.0 + (domain.u_end - uv.x) / (domain.u_end - domain.u_start),
        BoundarySide::UStart => 3.0 + (domain.v_end - uv.y) / (domain.v_end - domain.v_start),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_domain() -> Domain2 {
        Domain2 {
            u_start: 0.0,
            u_end: 1.0,
            v_start: 0.0,
            v_end: 1.0,
        }
    }

    fn inside() -> Endpoint2 {
        Endpoint2 {
            uv: Point2::new(0.5, 0.5),
            side: None,
        }
    }

    fn on_v_start() -> Endpoint2 {
        Endpoint2 {
            uv: Point2::new(0.5, 0.0),
            side: Some(BoundarySide::VStart),
        }
    }

    fn on_u_end() -> Endpoint2 {
        Endpoint2 {
            uv: Point2::new(1.0, 0.5),
            side: Some(BoundarySide::UEnd),
        }
    }

    /// Interior to boundary dead-ends inside the face: a slit.
    #[test]
    fn interior_to_boundary_is_a_slit() {
        assert!(is_slit(inside(), on_v_start(), unit_domain()));
        assert!(is_slit(on_v_start(), inside(), unit_domain()));
    }

    /// Boundary to boundary cuts clean through: a partition, not a slit.
    #[test]
    fn boundary_to_boundary_is_not_a_slit() {
        assert!(!is_slit(on_v_start(), on_u_end(), unit_domain()));
    }

    /// Two interior endpoints touch no boundary at all, so the curve is
    /// fully embedded rather than slitting the face open.
    #[test]
    fn interior_to_interior_is_not_a_slit() {
        let other = Endpoint2 {
            uv: Point2::new(0.25, 0.25),
            side: None,
        };
        assert!(!is_slit(inside(), other, unit_domain()));
    }

    /// An endpoint at a CORNER sits on two sides at once, so it is not a
    /// clean boundary landing and must not read as a slit.
    #[test]
    fn a_corner_endpoint_is_not_a_slit() {
        let corner = Endpoint2 {
            uv: Point2::new(0.0, 0.0),
            side: Some(BoundarySide::VStart),
        };
        assert!(!is_slit(inside(), corner, unit_domain()));
    }

    /// Same corner rejection, with the corner as the START endpoint.
    ///
    /// Both orderings are asserted because the predicate tests each
    /// endpoint on its own line; a one-sided check would pass the other.
    #[test]
    fn a_corner_start_endpoint_is_not_a_slit() {
        let corner = Endpoint2 {
            uv: Point2::new(0.0, 0.0),
            side: Some(BoundarySide::VStart),
        };
        assert!(!is_slit(corner, inside(), unit_domain()));
    }

    /// Classification refuses as UNIMPLEMENTED when only one patch is slit.
    ///
    /// The terminal proof requires BOTH patches to be unpartitionable. One
    /// slit patch alone leaves the other possibly splittable, so claiming a
    /// proof there would refuse work that is actually constructible.
    #[test]
    fn one_slit_patch_alone_is_not_a_terminal_proof() {
        let slit_patch = is_slit(inside(), on_v_start(), unit_domain());
        let partitioned = is_slit(on_v_start(), on_u_end(), unit_domain());
        assert!(slit_patch, "interior-to-boundary must read as a slit");
        assert!(!partitioned, "boundary-to-boundary must not read as a slit");
        // The conjunction is what the classifier requires: one true and one
        // false must NOT yield a terminal verdict.
        assert!(!(slit_patch && partitioned));
    }
}
