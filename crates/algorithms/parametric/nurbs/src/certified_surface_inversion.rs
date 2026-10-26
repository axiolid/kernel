//! Globally certified closest-point inversion for B-spline surfaces (#7).
//!
//! # Why inversion is not projection
//!
//! `project_surface_certified` answers "how far is this point from the
//! surface". Inversion answers "which parameters name this point", and that
//! is a strictly stronger claim: it requires the answer to be UNIQUE.
//!
//! A projection certificate retains EVERY box that could hold a global
//! minimizer. On a sphere-like patch a pole point is equidistant from a whole
//! parameter circle, so many boxes survive and no single (u, v) names the
//! point. Returning one of them would be an arbitrary choice dressed up as a
//! certified answer. This module therefore refuses unless the surviving
//! cover is a single connected box: ambiguity is reported, never resolved by
//! picking a representative.
//!
//! Inversion also asserts the point lies ON the surface. That is a separate
//! obligation from uniqueness, checked against the certified distance LOWER
//! bound so a miss cannot be absorbed by a loose representative evaluation.
//! A point is refused off-surface as soon as that separation is proven,
//! without refining to the parameter resolution (#320).

use axiolid_contracts::GeomResult;
use axiolid_core::{Point3, Scalar};
use axiolid_surface::BSplineSurface;

use crate::certified_projection::{
    CertifiedSurfaceProjection3, CertifiedSurfaceProjectionOptions, ParameterInterval,
    SurfaceParameterBox, SurfaceProjectionCertificate3, SurfaceProjectionUnresolvedReason,
};
use crate::certified_surface_projection::{
    search_periodic_surface_certified, search_surface_certified, SurfaceSearch,
};
use crate::periodic_surface::PeriodicBSplineSurface;

/// Why a sound inversion query did not yield unique parameters.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceInversionRefusal {
    /// The point is off the surface by more than the linear tolerance.
    ///
    /// Reported with the certified LOWER bound, so this is a proof of
    /// separation rather than one representative evaluation missing.
    OffSurface {
        /// Certified lower bound on the distance to the whole surface.
        distance_lower_bound: Scalar,
    },
    /// Several disjoint parameter regions attain the same minimum distance.
    ///
    /// A pole, a seam, or a self-touching patch. The retained cover is
    /// returned so callers can inspect the ambiguity instead of guessing.
    Ambiguous {
        /// Every closed box that may contain a global minimizer.
        candidates: Vec<SurfaceParameterBox>,
    },
    /// The search stayed sound but could not resolve to the requested
    /// accuracy within the configured budget.
    Unresolved {
        /// Exact reason the underlying projection stopped.
        reason: SurfaceProjectionUnresolvedReason,
    },
}

/// Unique native parameters proven to name a point on the surface.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceInversionCertificate3 {
    /// Native U parameter of the unique closest point.
    pub u: Scalar,
    /// Native V parameter of the unique closest point.
    pub v: Scalar,
    /// The single retained box proven to contain every global minimizer.
    pub enclosure: SurfaceParameterBox,
    /// Certified upper bound on how far the point is from the surface.
    ///
    /// Retained rather than discarded: an on-surface claim is only as good
    /// as the residual it was accepted with.
    pub residual_upper_bound: Scalar,
    /// Full underlying global projection certificate.
    pub projection: SurfaceProjectionCertificate3,
}

/// Globally invert a point against a clamped B-spline surface.
///
/// Returns `Ok(Ok(..))` only when the point is proven on the surface AND the
/// minimizer cover proves the parameters are unique.
///
/// A structurally sound refusal is `Ok(Err(..))`, not `Err`: an ambiguous or
/// off-surface point is a real answer about the geometry. `Err` is reserved
/// for invalid input and exhausted budgets.
pub fn invert_surface_certified(
    surface: &BSplineSurface,
    point: Point3,
    options: CertifiedSurfaceProjectionOptions,
) -> GeomResult<Result<SurfaceInversionCertificate3, SurfaceInversionRefusal>> {
    decide(
        search_surface_certified(
            surface,
            point,
            options,
            options.distance_tolerance().linear(),
        )?,
        options,
    )
}

/// Globally invert a point against one canonical period of a cyclic surface.
///
/// The seam is searched as part of the quotient domain, so a point exactly ON
/// the seam yields one enclosure rather than two rival endpoint boxes.
pub fn invert_periodic_surface_certified(
    surface: &PeriodicBSplineSurface,
    point: Point3,
    options: CertifiedSurfaceProjectionOptions,
) -> GeomResult<Result<SurfaceInversionCertificate3, SurfaceInversionRefusal>> {
    decide(
        search_periodic_surface_certified(
            surface,
            point,
            options,
            options.distance_tolerance().linear(),
        )?,
        options,
    )
}

/// Turn a certified surface search into an inversion verdict.
///
/// Order matters. The on-surface check comes BEFORE uniqueness: a plainly
/// off-surface point must be named as such, not mislabelled ambiguous
/// because the search stopped at its distance gap with a wide cover.
fn decide(
    search: SurfaceSearch,
    options: CertifiedSurfaceProjectionOptions,
) -> GeomResult<Result<SurfaceInversionCertificate3, SurfaceInversionRefusal>> {
    let outcome = match search {
        SurfaceSearch::Finished(outcome) => outcome,
        SurfaceSearch::Separated {
            distance_lower_bound,
        } => {
            return Ok(Err(SurfaceInversionRefusal::OffSurface {
                distance_lower_bound,
            }))
        }
    };
    let certificate = match outcome {
        CertifiedSurfaceProjection3::Complete(certificate) => certificate,
        // An unresolved projection has sound bounds but has NOT proven the
        // cover is final. Uniqueness cannot be concluded from it -- but a
        // multi-component cover already PROVES ambiguity, and refining
        // further can only split boxes, never reconnect them. A pole must
        // therefore be named as ambiguous, not excused as slow.
        CertifiedSurfaceProjection3::Unresolved {
            certificate,
            reason,
        } => {
            // A wide or disconnected cover already PROVES the answer is not
            // unique. Further refinement can only split boxes, never shrink
            // the span of the whole region, so this verdict is final.
            if unique_enclosure(&certificate.possible_minimizer_boxes, options).is_none() {
                return Ok(Err(SurfaceInversionRefusal::Ambiguous {
                    candidates: certificate.possible_minimizer_boxes,
                }));
            }
            return Ok(Err(SurfaceInversionRefusal::Unresolved { reason }));
        }
    };

    // On-surface FIRST. The search legitimately stops once the distance
    // gap is met, which for a far-away point leaves a wide parameter box;
    // judging uniqueness first would then mislabel a plainly off-surface
    // point as 'ambiguous' and hide the real reason.
    //
    // The verdict uses the certified LOWER bound: the representative
    // distance is one lucky evaluation and could accept a point the
    // surface provably never reaches.
    let linear = options.distance_tolerance().linear();
    if certificate.distance_lower_bound > linear {
        return Ok(Err(SurfaceInversionRefusal::OffSurface {
            distance_lower_bound: certificate.distance_lower_bound,
        }));
    }

    // Uniqueness needs BOTH connectedness and localization, so it is one
    // shared predicate: subdivision splits a genuine minimizer across
    // adjacent cells (so counting boxes is wrong), while a pole yields a
    // single CONNECTED strip spanning the whole domain (so connectivity
    // alone is not enough either).
    let Some(enclosure) = unique_enclosure(&certificate.possible_minimizer_boxes, options) else {
        return Ok(Err(SurfaceInversionRefusal::Ambiguous {
            candidates: certificate.possible_minimizer_boxes,
        }));
    };

    Ok(Ok(SurfaceInversionCertificate3 {
        u: certificate.u,
        v: certificate.v,
        enclosure,
        residual_upper_bound: certificate.distance_upper_bound,
        projection: certificate,
    }))
}

/// The single connected, localized region owning every global minimizer.
///
/// Returns `None` when the cover has two separated components (rival
/// minimizers) OR when the single component is wider than the parameter
/// tolerance (a pole, whose whole family ties).
///
/// Connectivity is computed by PAIRWISE touching, then transitively closed.
/// Growing a running hull instead would be unsound in the dangerous
/// direction: a hull spans the gap between two disjoint regions, so a third
/// box touching only that empty span would fuse rival minimizers and let an
/// ambiguous point be reported as uniquely invertible.
///
/// A leftover box therefore means two separated minimizer regions really
/// exist, and the point has no unique inverse.
fn unique_enclosure(
    boxes: &[SurfaceParameterBox],
    options: CertifiedSurfaceProjectionOptions,
) -> Option<SurfaceParameterBox> {
    let enclosure = single_component(boxes)?;
    let parameter = options.parameter_tolerance();
    // Each retained cell is already within the parameter tolerance, and a
    // point can lie on at most ONE cell boundary per axis, so a genuinely
    // unique minimizer spans at most two cells per axis. A pole spans the
    // whole family and is far wider, so this bound separates the two.
    let span = 2.0 * parameter;
    if width(enclosure.u) > span || width(enclosure.v) > span {
        return None;
    }
    Some(enclosure)
}

fn single_component(boxes: &[SurfaceParameterBox]) -> Option<SurfaceParameterBox> {
    let (first, rest) = boxes.split_first()?;
    let mut component = vec![*first];
    let mut remaining: Vec<SurfaceParameterBox> = rest.to_vec();
    let mut cursor = 0;
    while cursor < component.len() {
        let current = component[cursor];
        remaining.retain(|candidate| {
            if boxes_touch(&current, candidate) {
                component.push(*candidate);
                return false;
            }
            true
        });
        cursor += 1;
    }
    // Anything unreachable from the first box is a separate minimizer.
    if !remaining.is_empty() {
        return None;
    }
    component
        .into_iter()
        .reduce(|left, right| SurfaceParameterBox {
            u: hull(left.u, right.u),
            v: hull(left.v, right.v),
        })
}

/// Two closed boxes share at least a boundary point.
fn boxes_touch(left: &SurfaceParameterBox, right: &SurfaceParameterBox) -> bool {
    intervals_touch(left.u, right.u) && intervals_touch(left.v, right.v)
}

fn intervals_touch(left: ParameterInterval, right: ParameterInterval) -> bool {
    left.start <= right.end && right.start <= left.end
}

fn hull(left: ParameterInterval, right: ParameterInterval) -> ParameterInterval {
    ParameterInterval {
        start: left.start.min(right.start),
        end: left.end.max(right.end),
    }
}

fn width(interval: ParameterInterval) -> Scalar {
    (interval.end - interval.start).max(0.0)
}
