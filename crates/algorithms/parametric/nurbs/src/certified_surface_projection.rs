use std::cmp::Ordering;

use axiolid_contracts::{BackendId, GeomError, GeomResult, Operation};
use axiolid_core::{Point3, Scalar};
use axiolid_evaluate::surface::bspline_jet;
use axiolid_surface::BSplineSurface;

use crate::{
    certified_bezier::{
        distance_to_box_lower, distance_to_point_interval_upper, next_up, representative_distance,
    },
    certified_projection::{
        CertifiedSurfaceProjection3, CertifiedSurfaceProjectionOptions, ParameterInterval,
        SurfaceParameterBox, SurfaceProjectionCertificate3, SurfaceProjectionUnresolvedReason,
    },
    certified_refinement::RefinementBudget,
    certified_surface_bezier::{
        piecewise_bezier_patches, piecewise_periodic_bezier_patches, Patch,
    },
    PeriodicBSplineSurface,
};

const DIMENSIONS: usize = 3;

/// How a certified surface search ended.
#[derive(Debug, Clone)]
pub(crate) enum SurfaceSearch {
    /// The search ran to its normal completion or sound stop.
    Finished(CertifiedSurfaceProjection3),
    /// The certified distance lower bound already exceeded the caller's
    /// separation threshold, so the search stopped without refining to the
    /// parameter resolution (#320).
    Separated {
        /// Certified lower bound on the distance to the whole surface.
        distance_lower_bound: Scalar,
    },
}

/// Unwrap a search that was given no separation threshold.
fn finished(search: SurfaceSearch) -> GeomResult<CertifiedSurfaceProjection3> {
    match search {
        SurfaceSearch::Finished(projection) => Ok(projection),
        SurfaceSearch::Separated { .. } => Err(GeomError::Degenerate(
            "surface projection stopped at a separation threshold it was not given".to_owned(),
        )),
    }
}

#[derive(Debug, Clone, Copy)]
struct Pending {
    patch_index: usize,
    domain: SurfaceParameterBox,
    lower: Scalar,
    depth: u16,
    serial: u32,
}

#[derive(Debug, Clone, Copy)]
struct Candidate {
    u: Scalar,
    v: Scalar,
    point: Point3,
    distance: Scalar,
    upper: Scalar,
}

#[derive(Debug, Clone, Copy)]
enum SplitAxis {
    U,
    V,
}

/// Exhaustively bound the global point-to-surface minimum over the closed
/// native domain of a finite clamped NURBS surface.
///
/// This first certified slice rejects closed/periodic axes because closed
/// metadata alone does not define a verified periodic control-net topology.
/// `Complete` proves both the requested global distance gap and the requested
/// native U/V widths for every retained minimizer box. Depth/no-progress
/// termination returns a sound `Unresolved`; shared work or allocation
/// exhaustion returns `GeomError::BudgetExceeded`.
pub fn project_surface_certified(
    surface: &BSplineSurface,
    target: Point3,
    options: CertifiedSurfaceProjectionOptions,
) -> GeomResult<CertifiedSurfaceProjection3> {
    validate_query(surface, target)?;
    finished(project_surface_certified_with_modes(
        surface, target, options, false, false, None,
    )?)
}

/// [`project_surface_certified`] that stops as soon as the certified distance
/// lower bound exceeds `separation` (#320).
pub(crate) fn search_surface_certified(
    surface: &BSplineSurface,
    target: Point3,
    options: CertifiedSurfaceProjectionOptions,
    separation: Scalar,
) -> GeomResult<SurfaceSearch> {
    validate_query(surface, target)?;
    project_surface_certified_with_modes(surface, target, options, false, false, Some(separation))
}

/// Exhaustively bound the global point-to-surface minimum over one canonical
/// period of an explicitly validated cyclic B-spline surface.
///
/// Periodic axes are searched over their complete native period, including the
/// seam-equivalent endpoints. The returned minimizer boxes therefore form a
/// sound cover on the quotient domain; callers can canonicalize witness
/// parameters with [`PeriodicBSplineSurface::wrap_parameters`].
pub fn project_periodic_surface_certified(
    surface: &PeriodicBSplineSurface,
    target: Point3,
    options: CertifiedSurfaceProjectionOptions,
) -> GeomResult<CertifiedSurfaceProjection3> {
    validate_target(target)?;
    finished(project_surface_certified_with_modes(
        surface.as_bspline_surface(),
        target,
        options,
        surface.u_is_periodic(),
        surface.v_is_periodic(),
        None,
    )?)
}

/// [`project_periodic_surface_certified`] that stops as soon as the certified
/// distance lower bound exceeds `separation` (#320).
pub(crate) fn search_periodic_surface_certified(
    surface: &PeriodicBSplineSurface,
    target: Point3,
    options: CertifiedSurfaceProjectionOptions,
    separation: Scalar,
) -> GeomResult<SurfaceSearch> {
    validate_target(target)?;
    project_surface_certified_with_modes(
        surface.as_bspline_surface(),
        target,
        options,
        surface.u_is_periodic(),
        surface.v_is_periodic(),
        Some(separation),
    )
}

fn project_surface_certified_with_modes(
    surface: &BSplineSurface,
    target: Point3,
    options: CertifiedSurfaceProjectionOptions,
    u_periodic: bool,
    v_periodic: bool,
    separation: Option<Scalar>,
) -> GeomResult<SurfaceSearch> {
    let target_array = target.to_array();
    let mut budget = RefinementBudget::new(options.max_work(), "certified surface projection");
    let patches = if u_periodic || v_periodic {
        piecewise_periodic_bezier_patches(surface, u_periodic, v_periodic, &mut budget)?
    } else {
        piecewise_bezier_patches(surface, &mut budget)?
    };
    if patches.is_empty() {
        return Err(GeomError::InvalidInput(
            "certified surface projection requires at least one Bezier patch".to_owned(),
        ));
    }

    let root_work = u128::try_from(patches.len()).map_err(|_| search_overflow())?;
    budget.charge(Some(root_work))?;
    let mut pending = Vec::new();
    pending
        .try_reserve_exact(patches.len())
        .map_err(|_| allocation_error())?;
    let mut candidate = None;
    let mut visited_nodes = 0_u32;
    let mut serial = 0_u32;

    for (patch_index, patch) in patches.iter().enumerate() {
        let domain = patch_domain(patch);
        budget.charge(Some(patch.representative_bound_work()?))?;
        let lower = patch_lower(patch, target_array)?;
        update_candidate(
            &mut candidate,
            sample(surface, patch, domain, target_array)?,
        );
        pending.push(Pending {
            patch_index,
            domain,
            lower,
            depth: 0,
            serial,
        });
        serial = serial.checked_add(1).ok_or_else(search_overflow)?;
        visited_nodes = visited_nodes.checked_add(1).ok_or_else(search_overflow)?;
    }

    let mut candidate = candidate.ok_or_else(|| {
        GeomError::Degenerate("surface projection could not construct an upper witness".to_owned())
    })?;

    loop {
        pending.retain(|record| can_contain_global_minimizer(record.lower, candidate.upper));
        if pending.is_empty() {
            return Err(GeomError::Degenerate(
                "outward surface projection bounds excluded every candidate".to_owned(),
            ));
        }

        let lower = global_lower(&pending)?;
        // Sound: a box is dropped only when its lower bound exceeds an
        // attained upper bound, so every global minimiser lies in a retained
        // box and the minimum over `pending` is at most the true distance.
        if separation.is_some_and(|threshold| lower > threshold) {
            return Ok(SurfaceSearch::Separated {
                distance_lower_bound: lower,
            });
        }
        let parameter_ready = pending.iter().try_fold(true, |ready, record| {
            Ok::<_, GeomError>(
                ready
                    && interval_width(record.domain.u)? <= options.parameter_tolerance()
                    && interval_width(record.domain.v)? <= options.parameter_tolerance(),
            )
        })?;
        let gap = certified_gap(candidate.upper, lower)?;
        if gap <= options.distance_tolerance().linear() && parameter_ready {
            let boxes = take_sorted_boxes(pending)?;
            return Ok(SurfaceSearch::Finished(
                CertifiedSurfaceProjection3::Complete(certificate(
                    candidate,
                    lower,
                    boxes,
                    visited_nodes,
                )),
            ));
        }

        let selected = select_refinable(
            &pending,
            options.max_depth(),
            gap > options.distance_tolerance().linear(),
        )?;
        let Some(selected_index) = selected else {
            let reason = if pending
                .iter()
                .any(|record| record.depth >= options.max_depth())
            {
                SurfaceProjectionUnresolvedReason::DepthLimit
            } else {
                SurfaceProjectionUnresolvedReason::FloatingPointNoProgress
            };
            let boxes = take_sorted_boxes(pending)?;
            return Ok(SurfaceSearch::Finished(
                CertifiedSurfaceProjection3::Unresolved {
                    certificate: certificate(candidate, lower, boxes, visited_nodes),
                    reason,
                },
            ));
        };

        let record = pending.swap_remove(selected_index);
        let (axis, midpoint) = split_choice(record.domain)?.ok_or_else(|| {
            GeomError::Degenerate("selected surface parameter box cannot advance".to_owned())
        })?;
        let child_depth = record.depth.checked_add(1).ok_or_else(search_overflow)?;
        let patch = patches.get(record.patch_index).ok_or_else(|| {
            GeomError::Degenerate("surface projection patch index escaped its catalog".to_owned())
        })?;
        let child_node_work = patch
            .restriction_bound_work()?
            .checked_add(patch.representative_bound_work()?)
            .and_then(|work| work.checked_add(1))
            .ok_or_else(search_overflow)?;
        let child_work = child_node_work.checked_mul(2).ok_or_else(search_overflow)?;
        budget.charge(Some(child_work))?;
        let child_domains = split_domain(record.domain, axis, midpoint);
        let children = [
            make_pending(
                &patches,
                record.patch_index,
                child_domains.0,
                child_depth,
                serial,
                target_array,
            )?,
            make_pending(
                &patches,
                record.patch_index,
                child_domains.1,
                child_depth,
                serial.checked_add(1).ok_or_else(search_overflow)?,
                target_array,
            )?,
        ];
        serial = serial.checked_add(2).ok_or_else(search_overflow)?;
        visited_nodes = visited_nodes.checked_add(2).ok_or_else(search_overflow)?;

        for child in &children {
            update_resolved_candidate(
                &mut candidate,
                sample(
                    surface,
                    &patches[child.patch_index],
                    child.domain,
                    target_array,
                )?,
            );
        }
        for child in children {
            if can_contain_global_minimizer(child.lower, candidate.upper) {
                try_push(&mut pending, child)?;
            }
        }
    }
}

fn validate_query(surface: &BSplineSurface, target: Point3) -> GeomResult<()> {
    validate_target(target)?;
    if surface.u_closed || surface.v_closed {
        return Err(GeomError::Unsupported {
            backend: BackendId::new("axiolid-nurbs"),
            operation: Operation::SpatialQuery,
        });
    }
    Ok(())
}

fn validate_target(target: Point3) -> GeomResult<()> {
    if !target.is_finite() {
        return Err(GeomError::InvalidInput(
            "surface projection target must be finite".to_owned(),
        ));
    }
    Ok(())
}

fn patch_domain(patch: &Patch) -> SurfaceParameterBox {
    SurfaceParameterBox {
        u: ParameterInterval {
            start: patch.u_start,
            end: patch.u_end,
        },
        v: ParameterInterval {
            start: patch.v_start,
            end: patch.v_end,
        },
    }
}

fn patch_lower(patch: &Patch, target: [Scalar; 3]) -> GeomResult<Scalar> {
    let bounds = patch.coordinate_intervals()?;
    distance_to_box_lower(
        target,
        [bounds[0].lower(), bounds[1].lower(), bounds[2].lower()],
        [bounds[0].upper(), bounds[1].upper(), bounds[2].upper()],
        DIMENSIONS,
    )
}

fn restricted_lower(
    patch: &Patch,
    domain: SurfaceParameterBox,
    target: [Scalar; 3],
) -> GeomResult<Scalar> {
    let restricted = patch.restrict(domain.u.start, domain.u.end, domain.v.start, domain.v.end)?;
    patch_lower(&restricted, target)
}

fn sample(
    surface: &BSplineSurface,
    patch: &Patch,
    domain: SurfaceParameterBox,
    target: [Scalar; 3],
) -> GeomResult<Candidate> {
    let u = representative_parameter(domain.u);
    let v = representative_parameter(domain.v);
    let enclosed = patch.point_at(u, v)?.euclidean()?;
    let upper = distance_to_point_interval_upper(target, enclosed, DIMENSIONS)?;
    let point = bspline_jet(surface, u, v)?.point;
    if !point.is_finite() {
        return Err(GeomError::Degenerate(
            "surface projection scalar representative is non-finite".to_owned(),
        ));
    }
    let distance = representative_distance(point.to_array(), target, DIMENSIONS)?;
    Ok(Candidate {
        u,
        v,
        point,
        distance,
        upper,
    })
}

fn representative_parameter(interval: ParameterInterval) -> Scalar {
    let middle = interval.start * 0.5 + interval.end * 0.5;
    if middle > interval.start && middle < interval.end {
        middle
    } else {
        interval.start
    }
}

fn update_candidate(current: &mut Option<Candidate>, next: Candidate) {
    let replace = current.is_none_or(|best| candidate_order(next, best) == Ordering::Less);
    if replace {
        *current = Some(next);
    }
}

fn update_resolved_candidate(current: &mut Candidate, next: Candidate) {
    if candidate_order(next, *current) == Ordering::Less {
        *current = next;
    }
}

fn candidate_order(first: Candidate, second: Candidate) -> Ordering {
    first
        .upper
        .total_cmp(&second.upper)
        .then_with(|| first.u.total_cmp(&second.u))
        .then_with(|| first.v.total_cmp(&second.v))
}

fn make_pending(
    patches: &[Patch],
    patch_index: usize,
    domain: SurfaceParameterBox,
    depth: u16,
    serial: u32,
    target: [Scalar; 3],
) -> GeomResult<Pending> {
    let patch = patches.get(patch_index).ok_or_else(|| {
        GeomError::Degenerate("surface projection patch index escaped its catalog".to_owned())
    })?;
    Ok(Pending {
        patch_index,
        domain,
        lower: restricted_lower(patch, domain, target)?,
        depth,
        serial,
    })
}

fn select_refinable(
    pending: &[Pending],
    max_depth: u16,
    improve_gap: bool,
) -> GeomResult<Option<usize>> {
    let mut selected = None;
    for (index, record) in pending.iter().enumerate() {
        if record.depth >= max_depth || split_choice(record.domain)?.is_none() {
            continue;
        }
        selected = match selected {
            None => Some(index),
            Some(before) => {
                let ordering = if improve_gap {
                    lower_order(record, &pending[before])?
                } else {
                    widest_order(record, &pending[before])?
                };
                Some(if ordering == Ordering::Less {
                    index
                } else {
                    before
                })
            }
        };
    }
    Ok(selected)
}

fn lower_order(first: &Pending, second: &Pending) -> GeomResult<Ordering> {
    Ok(first
        .lower
        .total_cmp(&second.lower)
        .then_with(|| first.depth.cmp(&second.depth))
        .then_with(|| first.patch_index.cmp(&second.patch_index))
        .then_with(|| first.domain.u.start.total_cmp(&second.domain.u.start))
        .then_with(|| first.domain.v.start.total_cmp(&second.domain.v.start))
        .then_with(|| first.serial.cmp(&second.serial)))
}

fn widest_order(first: &Pending, second: &Pending) -> GeomResult<Ordering> {
    let first_width = interval_width(first.domain.u)?.max(interval_width(first.domain.v)?);
    let second_width = interval_width(second.domain.u)?.max(interval_width(second.domain.v)?);
    Ok(second_width
        .total_cmp(&first_width)
        .then_with(|| first.patch_index.cmp(&second.patch_index))
        .then_with(|| first.domain.u.start.total_cmp(&second.domain.u.start))
        .then_with(|| first.domain.v.start.total_cmp(&second.domain.v.start))
        .then_with(|| first.serial.cmp(&second.serial)))
}

fn split_choice(domain: SurfaceParameterBox) -> GeomResult<Option<(SplitAxis, Scalar)>> {
    let u_width = interval_width(domain.u)?;
    let v_width = interval_width(domain.v)?;
    let u_midpoint = advancing_midpoint(domain.u);
    let v_midpoint = advancing_midpoint(domain.v);
    Ok(match (u_midpoint, v_midpoint) {
        (Some(u), Some(_)) if u_width >= v_width => Some((SplitAxis::U, u)),
        (Some(_), Some(v)) => Some((SplitAxis::V, v)),
        (Some(u), None) => Some((SplitAxis::U, u)),
        (None, Some(v)) => Some((SplitAxis::V, v)),
        (None, None) => None,
    })
}

fn advancing_midpoint(interval: ParameterInterval) -> Option<Scalar> {
    let middle = interval.start * 0.5 + interval.end * 0.5;
    (middle > interval.start && middle < interval.end).then_some(middle)
}

fn split_domain(
    domain: SurfaceParameterBox,
    axis: SplitAxis,
    midpoint: Scalar,
) -> (SurfaceParameterBox, SurfaceParameterBox) {
    match axis {
        SplitAxis::U => (
            SurfaceParameterBox {
                u: ParameterInterval {
                    start: domain.u.start,
                    end: midpoint,
                },
                v: domain.v,
            },
            SurfaceParameterBox {
                u: ParameterInterval {
                    start: midpoint,
                    end: domain.u.end,
                },
                v: domain.v,
            },
        ),
        SplitAxis::V => (
            SurfaceParameterBox {
                u: domain.u,
                v: ParameterInterval {
                    start: domain.v.start,
                    end: midpoint,
                },
            },
            SurfaceParameterBox {
                u: domain.u,
                v: ParameterInterval {
                    start: midpoint,
                    end: domain.v.end,
                },
            },
        ),
    }
}

fn interval_width(interval: ParameterInterval) -> GeomResult<Scalar> {
    let width = interval.end - interval.start;
    if width.is_finite() && width >= 0.0 {
        Ok(width)
    } else {
        Err(GeomError::Degenerate(
            "surface projection native parameter width is non-finite".to_owned(),
        ))
    }
}

fn global_lower(pending: &[Pending]) -> GeomResult<Scalar> {
    let lower = pending
        .iter()
        .map(|record| record.lower)
        .fold(Scalar::INFINITY, Scalar::min);
    if lower.is_finite() {
        Ok(lower)
    } else {
        Err(GeomError::Degenerate(
            "surface projection global lower bound is non-finite".to_owned(),
        ))
    }
}

fn can_contain_global_minimizer(lower: Scalar, attained_upper: Scalar) -> bool {
    lower <= attained_upper
}

fn certified_gap(upper: Scalar, lower: Scalar) -> GeomResult<Scalar> {
    if lower > upper {
        return Err(GeomError::Degenerate(
            "surface projection lower bound exceeds its attained upper bound".to_owned(),
        ));
    }
    let gap = next_up((upper - lower).max(0.0));
    if gap.is_finite() {
        Ok(gap)
    } else {
        Err(GeomError::Degenerate(
            "surface projection distance gap is non-finite".to_owned(),
        ))
    }
}

fn take_sorted_boxes(pending: Vec<Pending>) -> GeomResult<Vec<SurfaceParameterBox>> {
    let mut boxes = Vec::new();
    boxes
        .try_reserve_exact(pending.len())
        .map_err(|_| allocation_error())?;
    for record in pending {
        boxes.push(record.domain);
    }
    boxes.sort_by(|first, second| {
        first
            .u
            .start
            .total_cmp(&second.u.start)
            .then_with(|| first.v.start.total_cmp(&second.v.start))
            .then_with(|| first.u.end.total_cmp(&second.u.end))
            .then_with(|| first.v.end.total_cmp(&second.v.end))
    });
    Ok(boxes)
}

fn certificate(
    candidate: Candidate,
    lower: Scalar,
    boxes: Vec<SurfaceParameterBox>,
    visited_nodes: u32,
) -> SurfaceProjectionCertificate3 {
    SurfaceProjectionCertificate3 {
        u: candidate.u,
        v: candidate.v,
        point: candidate.point,
        distance: candidate.distance,
        distance_lower_bound: lower,
        distance_upper_bound: candidate.upper,
        possible_minimizer_boxes: boxes,
        visited_nodes,
    }
}

fn try_push<T>(values: &mut Vec<T>, value: T) -> GeomResult<()> {
    values.try_reserve(1).map_err(|_| allocation_error())?;
    values.push(value);
    Ok(())
}

fn allocation_error() -> GeomError {
    GeomError::BudgetExceeded {
        resource: "certified surface projection allocation",
    }
}

fn search_overflow() -> GeomError {
    GeomError::BudgetExceeded {
        resource: "certified surface projection search nodes",
    }
}

#[cfg(test)]
mod tests {
    use super::can_contain_global_minimizer;

    #[test]
    fn equality_is_retained() {
        assert!(can_contain_global_minimizer(1.0, 1.0));
    }
}
