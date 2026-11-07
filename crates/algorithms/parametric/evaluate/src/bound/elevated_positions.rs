//! Certified positions of elevated and banked curves over an intrinsic
//! plan (#393) or a chain plan (#438).
//!
//! An [`Elevated3`](axiolid_curve::Elevated3) is `(p(d), z(d))`: the plan
//! `p` read by its own arc length and a closed-form height `z`. A
//! [`Banked3`](axiolid_curve::Banked3) raises that point by the closed-form
//! pivot `e(d)`. Only `p` is an integral, so the point's error is the plan's:
//! [`ElevatedPositionTable`] reads `p` from an
//! [`IntrinsicPositionTable`] (or, over a chain plan, a
//! [`ChainPositionTable`]) and adds the height (and the pivot) in closed
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
//! (#438) the ends are read as a chain's chord bound reads them
//! ([`crate::chain::chain_chord_bound`]): the lower one as the piece that
//! starts there, the upper one as the piece that ends there. Under an
//! intrinsic profile (#481) the exact-chord bound is returned as it
//! stands, and UNDER-reports the evaluated ends: `certifies_flattening3`
//! does not claim it, so no consumer takes it as a certificate.
//!
//! # Jets over a chain plan (#438)
//!
//! A swept tube reads the tangent and the second derivative of its centre
//! line too. Over an intrinsic plan those are closed forms of the heading,
//! exact but for rounding; over a chain plan a parametric piece's are read at
//! the evaluator's parameter. [`ElevatedPositionTable::jet3_on`] reads them
//! from the chain table instead (the plan's `p'` and `p''`, the profile's
//! `z'` and `z''` in closed form), and [`ElevatedPositionTable::jet_errors`]
//! bounds how far each is from the exact one.

use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_curve::{Curve2, Curve3, Elevated3, ElevationLaw};

use super::interval::Iv;
use super::{
    frame_stretch2, inflate, ordered, ChainPositionTable, IntrinsicPositionTable, SeamSide,
};

/// The certified positions of the plan of an [`ElevatedPositionTable`].
#[derive(Debug, Clone)]
enum PlanPositions {
    Intrinsic(IntrinsicPositionTable),
    Chain(ChainPositionTable),
}

/// Certified positions of an elevated or banked curve over an intrinsic or
/// a chain plan, in plan distance, read in any order: the plan from an
/// [`IntrinsicPositionTable`] or a [`ChainPositionTable`], the height and
/// the pivot in closed form.
#[derive(Debug, Clone)]
pub struct ElevatedPositionTable {
    plan: PlanPositions,
    curve: Curve3,
    end: Scalar,
}

impl ElevatedPositionTable {
    /// The table of `curve` over `[0, end]`, its plan points within about
    /// `error` ([`IntrinsicPositionTable::new2`],
    /// [`ChainPositionTable::new`]). `Some` only for an elevated or banked
    /// curve whose plan is an intrinsic curve or a chain a
    /// [`ChainPositionTable`] covers, whose profile is closed form (no
    /// intrinsic profile), and for a banked curve whose cant and pivot laws
    /// are well formed with no held-rail piece in the cant; `None`
    /// otherwise, for an error or an end that is not finite and positive,
    /// or an end past the plan's length.
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
        if !crate::elevated::certifies_elevated(base) {
            return None;
        }
        let plan = match base.plan.as_ref() {
            Curve2::Intrinsic(plan) => {
                PlanPositions::Intrinsic(IntrinsicPositionTable::new2(plan, error, end)?)
            }
            Curve2::Chain(chain) => {
                PlanPositions::Chain(ChainPositionTable::new(chain, error, end)?)
            }
            _ => return None,
        };
        Some(Self {
            plan,
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
        self.point3_on(d, SeamSide::Outgoing)
    }

    /// [`Self::point3`] on `side` of a join of a chain plan (the side makes
    /// no difference over an intrinsic plan, or off a join).
    #[must_use]
    pub fn point3_on(&self, d: Scalar, side: SeamSide) -> Option<(Point3, Scalar)> {
        if !(d.is_finite() && (0.0..=self.end).contains(&d)) {
            return None;
        }
        let (plan, error) = match &self.plan {
            PlanPositions::Intrinsic(table) => table.point2(d)?,
            PlanPositions::Chain(table) => table.point2_on(d, side)?,
        };
        let (z, _, _) = self.vertical(d)?;
        let out = Point3::new(plan.x, plan.y, z);
        out.is_finite().then_some((out, error))
    }

    /// The height, grade and second derivative of the profile (and the
    /// pivot) at `d`, in closed form, as the evaluator reads them.
    fn vertical(&self, d: Scalar) -> Option<(Scalar, Scalar, Scalar)> {
        let (law, pivot) = match &self.curve {
            Curve3::Elevated(e) => (&e.elevation, None),
            Curve3::Banked(b) => (&b.base.elevation, Some(b)),
            _ => return None,
        };
        let height = crate::elevation::elevation_height(law, d).ok()?;
        let grade = crate::elevation::elevation_grade(law, d).ok()?;
        let bend = crate::elevation::elevation_second(law, d).ok()?;
        let (e, rate, extra) = match pivot {
            Some(b) => {
                let (e, rate) = b.pivot_at(d).ok()?;
                (e, rate, crate::banked::pivot_jet(b, d, false)?[2])
            }
            None => (0.0, 0.0, 0.0),
        };
        let out = (height + e, grade + rate, bend + extra);
        (out.0.is_finite() && out.1.is_finite() && out.2.is_finite()).then_some(out)
    }

    /// The point, the derivative and the second derivative at plan distance
    /// `d` on `side` of a join, `(c, c', c'')`: over a chain plan the plan
    /// parts (`p'`, `p''`) from the chain table and the profile's and the
    /// pivot's in closed form; over an intrinsic plan the curve's own
    /// closed-form derivatives, bit for bit. [`Self::jet_errors`] bounds
    /// their distances from the exact ones.
    #[must_use]
    pub fn jet3_on(&self, d: Scalar, side: SeamSide) -> Option<(Point3, Vec3, Vec3)> {
        let (point, _) = self.point3_on(d, side)?;
        let (first, second) = match &self.plan {
            PlanPositions::Intrinsic(_) => (
                crate::curve::derivative3(&self.curve, d).ok()?,
                crate::curve::second_derivative3(&self.curve, d).ok()?,
            ),
            PlanPositions::Chain(table) => {
                let (_, grade, bend) = self.vertical(d)?;
                let (tangent, _) = table.tangent2_on(d, side)?;
                let (curvature, _) = table.curvature2_on(d, side)?;
                (
                    Vec3::new(tangent.x, tangent.y, grade),
                    Vec3::new(curvature.x, curvature.y, bend),
                )
            }
        };
        (first.is_finite() && second.is_finite()).then_some((point, first, second))
    }

    /// Certified bounds, uniform over the table, on how far a point, a
    /// derivative and a second derivative it returns are from the exact
    /// ones: `[e0, e1, e2]`. The plan's `max_error`, `max_tangent_error` and
    /// `max_curvature_error` over a chain plan (the closed-form profile adds
    /// none); `[max_error, 0, 0]` over an intrinsic plan, whose derivatives
    /// are exact closed forms of the heading.
    #[must_use]
    pub fn jet_errors(&self) -> [Scalar; 3] {
        match &self.plan {
            PlanPositions::Intrinsic(table) => [table.max_error(), 0.0, 0.0],
            PlanPositions::Chain(table) => [
                table.max_error(),
                table.max_tangent_error(),
                table.max_curvature_error(),
            ],
        }
    }

    /// A certified bound on the error of every point the table returns:
    /// the `max_error` of its plan's table.
    #[must_use]
    pub fn max_error(&self) -> Scalar {
        self.jet_errors()[0]
    }
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
    let (start, length) = match base.plan.as_ref() {
        Curve2::Intrinsic(plan) => (plan.start, plan.length),
        Curve2::Chain(chain) => (chain.start, chain.length()?),
        _ => return Some(law),
    };
    if has_intrinsic_piece(&base.elevation) {
        // TODO(#481): a vertical clothoid's heights carry no certified
        // error; `certifies_flattening3` does not claim it either.
        return Some(law);
    }
    let (lo, hi) = ordered(a, b);
    let stretch = frame_stretch2(&start);
    let origin = start.origin.abs().max_element();
    let floor = 64.0 * Scalar::EPSILON * (1.0 + origin + stretch * hi);
    let table = ElevatedPositionTable::new(curve, (law / 8.0).max(floor), hi.min(length))?;
    let gap = |s: Scalar, side: SeamSide| {
        let (exact, error) = table.point3_on(s.min(length), side)?;
        let read = crate::curve::evaluate3(curve, s).ok()?;
        let d = [read.x - exact.x, read.y - exact.y, read.z - exact.z]
            .map(|c| Iv::point(c).sq())
            .into_iter()
            .fold(Iv::point(0.0), Iv::add)
            .sqrt();
        Some(Iv::point(d.hi).add(Iv::point(error)).hi)
    };
    let g = gap(lo, SeamSide::Outgoing)?.max(gap(hi, SeamSide::Incoming)?);
    let out = inflate(Iv::point(law).add(Iv::point(g)).hi);
    out.is_finite().then_some(out)
}
