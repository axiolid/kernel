//! Certified positions of elevated and banked curves over an intrinsic
//! plan (#393).
//!
//! An [`Elevated3`](axiolid_curve::Elevated3) is `(p(d), z(d))`: the plan
//! `p` read by its own arc length and a closed-form height `z`. A
//! [`Banked3`](axiolid_curve::Banked3) raises that point by the closed-form
//! pivot `e(d)`. Only `p` is an integral, so the point's error is the plan's:
//! [`ElevatedPositionTable`] reads `p` from an
//! [`IntrinsicPositionTable`] and adds the height (and the pivot) in closed
//! form. The closed-form height and pivot are not charged their own
//! rounding, as for every closed-form point in the tube certificate of
//! `axiolid-mesh-compile`.
//!
//! The evaluator's own point (`crate::arc_length::intrinsic_point`) is an
//! estimate with no certified error, off by 1.4 cm at 96.5 m on a fast
//! harmonic: [`chord_bound3`](super::chord_bound3) and the swept-disk
//! certificate read this table instead of trusting it.
//!
//! # Evaluated ends (#393)
//!
//! [`chord_bound3`](super::chord_bound3) bounds the distance from the arc to
//! the chord of the *evaluated* ends. For an elevated or banked curve over
//! an intrinsic plan the exact-chord bound of `crate::elevated` (`h^2 / 8
//! sup |k|` composed with the profile's) is not that: the evaluated plan
//! point is off the exact one, and the chord of two moved ends moves by
//! at most the larger of the two moves. `evaluated_chord_bound3` adds that
//! larger gap, each the evaluator's point against a table point plus the
//! table point's own error, to the exact-chord bound. Over a chain plan
//! (#438) and under an intrinsic profile (#481) the exact-chord bound is
//! returned as it stands, and UNDER-reports the evaluated ends (a chain of a
//! line and `k = 0.01 sin 50 s` gives 1.95e-5 over a span whose chord lies
//! 2.8e-3 off): `certifies_flattening3` claims neither, so no consumer takes
//! it as a certificate.

use axiolid_core::{Point3, Scalar};
use axiolid_curve::{Banked3, Curve2, Curve3, Elevated3, ElevationLaw};

use super::interval::Iv;
use super::{frame_stretch2, inflate, ordered, IntrinsicPositionTable};

/// Certified positions of an elevated or banked curve over an intrinsic
/// plan, in plan distance, read in any order: the plan from an
/// [`IntrinsicPositionTable`], the height and the pivot in closed form.
#[derive(Debug, Clone)]
pub struct ElevatedPositionTable {
    plan: IntrinsicPositionTable,
    curve: Curve3,
    end: Scalar,
}

impl ElevatedPositionTable {
    /// The table of `curve` over `[0, end]`, its plan points within about
    /// `error` ([`IntrinsicPositionTable::new2`]). `Some` only for an
    /// elevated or banked curve whose plan is an intrinsic curve, whose
    /// profile is closed form (no intrinsic profile), and for a banked
    /// curve whose cant and pivot laws are well formed with no held-rail
    /// piece in the cant; `None` otherwise, for an error or an end that is
    /// not finite and positive, or an end past the plan's length.
    #[must_use]
    pub fn new(curve: &Curve3, error: Scalar, end: Scalar) -> Option<Self> {
        let base: &Elevated3 = match curve {
            Curve3::Elevated(e) => e,
            Curve3::Banked(b) => {
                if !(b.cant.is_well_formed() && b.pivot.is_well_formed())
                    || b.cant.has_rail_pieces()
                {
                    return None;
                }
                &b.base
            }
            _ => return None,
        };
        let Curve2::Intrinsic(plan) = base.plan.as_ref() else {
            return None;
        };
        if !crate::elevated::certifies_elevated(base) {
            return None;
        }
        Some(Self {
            plan: IntrinsicPositionTable::new2(plan, error, end)?,
            curve: curve.clone(),
            end,
        })
    }

    /// The point at plan distance `d` and a certified bound on its
    /// distance from the exact point: the plan table's error at `d`, with
    /// the closed-form height (and pivot) exact. `None` outside `[0, end]`
    /// or where the profile or the pivot has no value.
    #[must_use]
    pub fn point3(&self, d: Scalar) -> Option<(Point3, Scalar)> {
        if !(d.is_finite() && (0.0..=self.end).contains(&d)) {
            return None;
        }
        let (plan, error) = self.plan.point2(d)?;
        let (law, pivot) = match &self.curve {
            Curve3::Elevated(e) => (&e.elevation, 0.0),
            Curve3::Banked(b) => (&b.base.elevation, pivot_at(b, d)?),
            _ => return None,
        };
        let height = crate::elevation::elevation_height(law, d).ok()?;
        let z = height + pivot;
        let out = Point3::new(plan.x, plan.y, z);
        out.is_finite().then_some((out, error))
    }

    /// A certified bound on the error of every point the table returns:
    /// [`IntrinsicPositionTable::max_error`] of its plan.
    #[must_use]
    pub fn max_error(&self) -> Scalar {
        self.plan.max_error()
    }
}

fn pivot_at(curve: &Banked3, d: Scalar) -> Option<Scalar> {
    curve.pivot_at(d).ok().map(|(e, _)| e)
}

/// Whether any piece of `law` is an intrinsic (vertical clothoid) profile.
fn has_intrinsic_piece(law: &ElevationLaw) -> bool {
    match law {
        ElevationLaw::Intrinsic { .. } => true,
        ElevationLaw::Piecewise { laws, .. } => laws.iter().any(has_intrinsic_piece),
        _ => false,
    }
}

/// Certified chord bound of an elevated or banked curve over plan
/// distances `[a, b]`, for the *evaluated* ends (see the module
/// documentation). `None` where the exact-chord bound is `None`, or where
/// the table, the evaluator or the pivot refuses an end.
pub(super) fn evaluated_chord_bound3(curve: &Curve3, a: Scalar, b: Scalar) -> Option<Scalar> {
    let (law, base) = match curve {
        Curve3::Elevated(e) => (crate::elevated::elevated_chord_bound(e, a, b)?, e),
        Curve3::Banked(banked) => (
            crate::elevated::banked_chord_bound(banked, a, b)?,
            &banked.base,
        ),
        _ => return None,
    };
    let Curve2::Intrinsic(plan) = base.plan.as_ref() else {
        // TODO(#438): the chain plan's evaluated pieces carry no certified
        // gap yet, so this under-reports the evaluated ends. Not a
        // certificate: `certifies_flattening3` refuses a chain plan.
        return Some(law);
    };
    if has_intrinsic_piece(&base.elevation) {
        // TODO(#481): a vertical clothoid's heights carry no certified
        // error; `certifies_flattening3` does not claim it either.
        return Some(law);
    }
    let (lo, hi) = ordered(a, b);
    let stretch = frame_stretch2(&plan.start);
    let origin = plan.start.origin.abs().max_element();
    let floor = 64.0 * Scalar::EPSILON * (1.0 + origin + stretch * hi);
    let table = ElevatedPositionTable::new(curve, (law / 8.0).max(floor), hi)?;
    let gap = |s: Scalar| {
        let (exact, error) = table.point3(s)?;
        let read = crate::curve::evaluate3(curve, s).ok()?;
        let d = [read.x - exact.x, read.y - exact.y, read.z - exact.z]
            .map(|c| Iv::point(c).sq())
            .into_iter()
            .fold(Iv::point(0.0), Iv::add)
            .sqrt();
        Some(Iv::point(d.hi).add(Iv::point(error)).hi)
    };
    let g = gap(lo)?.max(gap(hi)?);
    let out = inflate(Iv::point(law).add(Iv::point(g)).hi);
    out.is_finite().then_some(out)
}
