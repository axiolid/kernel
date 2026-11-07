//! Frames read along intrinsic curves, certified (#361, ADR 0082
//! amendment).
//!
//! An intrinsic curve's heading `theta(s)` is the closed-form integral of
//! its curvature law, so every axis of a frame on it is closed form. Only
//! its POINT is a quadrature (`crate::arc_length::intrinsic_point`, an
//! estimate with no error bound), or on a space curve a Magnus integration
//! (`crate::frenet`). As for the other curves read by quadrature
//! ([`super::reading`]), the frame law these bounds are about is the EXACT
//! one, and each frame the kernel reads is compared with it:
//!
//! - the exact point comes from a certified
//!   [`IntrinsicPositionTable`](crate::bound::IntrinsicPositionTable), and
//!   its error is part of the charge;
//! - the exact tangent is the closed-form heading through the start frame
//!   (on an elevated curve lifted by the profile's grade), and the exact
//!   lateral and up axes follow from it as the kernel frames it
//!   ([`framed`]).
//!
//! The measure is exact: arc length (plan distance on an elevated curve)
//! is the curve's parameter. The families: a 2D intrinsic curve, a level
//! planar intrinsic space curve, an elevated curve over an intrinsic
//! plan with a profile of closed-form pieces, and a banked curve over one
//! (#396), each with an orthonormal start frame. A banked curve's point is
//! the plan's table point raised by the profile's height and the pivot, its
//! tangent the plan heading lifted by the grade and the pivot's rate, and
//! its axes the rolled section the kernel frames
//! ([`crate::banked::rolled_axes`]), the roll from the cant law; a held-rail
//! cant is not read. A chain, a stretched or skewed start frame, a tilted or
//! twisted space curve and a placement that tilts `+Z` are not read here.

use axiolid_contracts::GeomResult;
use axiolid_core::{Frame2, Scalar, Vec2, Vec3};
use axiolid_curve::{Banked3, Curve2, Curve3, Elevated3, Intrinsic2, Intrinsic3, SeamSide};

use super::reading::{framed, AXIS_ROUNDING};
use super::POINT_ROUNDING;
use crate::bound::{intrinsic_derivative_bounds2, CurveDerivativeBounds, IntrinsicPositionTable};
use crate::station::composite::{StationCurve, StationPiece};
use crate::station::{SectionFrame, KEEPS_UP_TOLERANCE};

/// How far an orthonormal start frame's axes may stray from unit length
/// and from each other.
const NEAR: Scalar = 1e-12;

/// What a piece's frames are read along.
#[derive(Clone, Copy)]
enum Base<'c> {
    /// A 2D intrinsic curve.
    Plan(&'c Intrinsic2),
    /// A level planar space curve, with its start frame's `x`, `y` plane
    /// read as a 2D one.
    Level(&'c Intrinsic3),
    /// An elevated curve over an intrinsic plan.
    Elevated(&'c Elevated3, &'c Intrinsic2),
    /// A banked curve over an intrinsic plan.
    Banked(&'c Banked3, &'c Intrinsic2),
}

fn orthonormal(frame: &Frame2) -> bool {
    (frame.x.length() - 1.0).abs() <= NEAR
        && (frame.y.length() - 1.0).abs() <= NEAR
        && frame.x.dot(frame.y).abs() <= NEAR
}

impl<'c> Base<'c> {
    /// The family of `curve`, with an orthonormal start frame; `None` for
    /// everything else.
    fn of(curve: StationCurve<'c>) -> Option<Self> {
        match curve {
            StationCurve::Two(Curve2::Intrinsic(plan)) => {
                orthonormal(&plan.start).then_some(Self::Plan(plan))
            }
            StationCurve::Three(Curve3::Elevated(elevated)) => {
                let Curve2::Intrinsic(plan) = elevated.plan.as_ref() else {
                    return None;
                };
                (orthonormal(&plan.start) && crate::elevated::certifies_elevated(elevated))
                    .then_some(Self::Elevated(elevated, plan))
            }
            StationCurve::Three(Curve3::Banked(banked)) => {
                let Curve2::Intrinsic(plan) = banked.base.plan.as_ref() else {
                    return None;
                };
                (orthonormal(&plan.start)
                    && crate::elevated::certifies_elevated(&banked.base)
                    && banked.cant.is_well_formed()
                    && banked.pivot.is_well_formed()
                    && !banked.cant.has_rail_pieces())
                .then_some(Self::Banked(banked, plan))
            }
            StationCurve::Three(Curve3::Intrinsic(space)) => {
                let start = &space.start;
                let level = start.x.z == 0.0 && start.y.z == 0.0;
                (space.is_planar() && level && orthonormal(&plane_frame(space)))
                    .then_some(Self::Level(space))
            }
            _ => None,
        }
    }
}

/// A level space curve's start frame in its plane.
fn plane_frame(space: &Intrinsic3) -> Frame2 {
    Frame2 {
        origin: axiolid_core::Point2::new(space.start.origin.x, space.start.origin.y),
        x: Vec2::new(space.start.x.x, space.start.x.y),
        y: Vec2::new(space.start.y.x, space.start.y.y),
    }
}

/// Whether the frames of `curve` are read through an
/// [`IntrinsicReader`]: an intrinsic family of the module documentation.
pub(crate) fn reads_intrinsic(curve: StationCurve<'_>) -> bool {
    Base::of(curve).is_some()
}

/// A piece of an intrinsic curve whose frames are bounded against the
/// exact frame law, with the certified positions of the curve.
pub(crate) struct IntrinsicReader<'c> {
    base: Base<'c>,
    /// The plan curve in the start frame's plane.
    local: Intrinsic2,
    table: IntrinsicPositionTable,
    end: Scalar,
}

impl<'c> IntrinsicReader<'c> {
    /// The reader of `piece`, its positions certified to `tolerance`;
    /// `None` for a piece that is not an intrinsic family of the
    /// [module documentation](self), is placed in a frame that tilts `+Z`,
    /// or whose positions the table cannot certify.
    pub(crate) fn new(piece: &StationPiece<'c>, tolerance: Scalar) -> GeomResult<Option<Self>> {
        let Some(base) = Base::of(piece.curve) else {
            return Ok(None);
        };
        if let Some(rigid) = piece.placement {
            if (rigid.transform_vector3(Vec3::Z) - Vec3::Z).length() > KEEPS_UP_TOLERANCE {
                return Ok(None);
            }
        }
        let (local, table) = match base {
            Base::Plan(plan) | Base::Elevated(_, plan) | Base::Banked(_, plan) => (
                plan.clone(),
                IntrinsicPositionTable::new2(plan, tolerance, piece.end),
            ),
            Base::Level(space) => (
                Intrinsic2::new(plane_frame(space), space.curvature.clone(), space.length),
                IntrinsicPositionTable::new3(space, tolerance, piece.end),
            ),
        };
        Ok(table.map(|table| Self {
            base,
            local,
            table,
            end: piece.end,
        }))
    }

    /// Whether the curve is an elevated or a banked one, measured in plan
    /// distance.
    pub(crate) fn is_plan_measured(&self) -> bool {
        matches!(self.base, Base::Elevated(..) | Base::Banked(..))
    }

    /// The banked curve, if it is one.
    pub(crate) fn banked(&self) -> Option<&'c Banked3> {
        match self.base {
            Base::Banked(banked, _) => Some(banked),
            _ => None,
        }
    }

    /// Certified suprema of the curve's derivatives in its own measure over
    /// `[lo, hi]` (`third` infinite across a seam of its law), and its
    /// derivative at the middle: the unit tangent of a plan or level curve,
    /// an elevated or banked curve's `(p', z')` in plan distance.
    pub(crate) fn derivatives(
        &self,
        lo: Scalar,
        hi: Scalar,
    ) -> Option<(CurveDerivativeBounds, Vec3)> {
        let mid = 0.5 * (lo + hi);
        match self.base {
            Base::Elevated(elevated, _) => {
                let bounds = crate::elevated::elevated_derivative_bounds(elevated, lo, hi)?;
                let middle = crate::elevated::elevated_derivative(elevated, mid).ok()?;
                Some((bounds, middle))
            }
            Base::Banked(banked, _) => {
                let bounds = crate::elevated::banked_derivative_bounds(banked, lo, hi)?;
                let middle = crate::banked::banked_derivative(banked, mid).ok()?;
                Some((bounds, middle))
            }
            Base::Plan(_) | Base::Level(_) => {
                let bounds = intrinsic_derivative_bounds2(&self.local, lo, hi)?;
                let (sin, cos) = self.local.heading_at(mid)?.sin_cos();
                let t = self.local.start.x * cos + self.local.start.y * sin;
                Some((bounds, Vec3::new(t.x, t.y, 0.0)))
            }
        }
    }

    /// The curve's closed-form jet `(O', O'')` at its measure `m`, per unit
    /// of that measure (#361): the unit tangent of a plan or level curve
    /// and `k` times its normal, or an elevated curve's `(p', z')` and
    /// `(k n, z'')` in plan distance, a banked curve's with the pivot's rate
    /// and bend added to the profile's. `None` where the law has no value at
    /// `m`.
    pub(crate) fn jet(&self, m: Scalar) -> Option<(Vec3, Vec3)> {
        match self.base {
            Base::Elevated(elevated, _) => Some((
                crate::elevated::elevated_derivative(elevated, m).ok()?,
                crate::elevated::elevated_second_derivative(elevated, m).ok()?,
            )),
            Base::Banked(banked, _) => Some((
                crate::banked::banked_derivative(banked, m).ok()?,
                crate::banked::banked_second_derivative(banked, m).ok()?,
            )),
            Base::Plan(_) | Base::Level(_) => {
                let (sin, cos) = self.local.heading_at(m)?.sin_cos();
                let start = &self.local.start;
                let tangent = start.x * cos + start.y * sin;
                let normal = start.y * cos - start.x * sin;
                let kappa = crate::frenet::law_value(&self.local.curvature, m)?;
                Some((
                    Vec3::new(tangent.x, tangent.y, 0.0),
                    kappa * Vec3::new(normal.x, normal.y, 0.0),
                ))
            }
        }
    }

    /// The exact frame at the curve's measure `m` read from `side`, the
    /// certified error of its point and the roll of its axes (zero off a
    /// banked curve).
    fn exact(
        &self,
        m: Scalar,
        side: SeamSide,
    ) -> GeomResult<Option<(SectionFrame, Scalar, Scalar)>> {
        let slack = 1e-9 * self.end.max(1.0);
        if !m.is_finite() || m < -slack || m > self.end + slack {
            return Ok(None);
        }
        let m = m.clamp(0.0, self.end);
        let Some(theta) = self.local.heading_at(m) else {
            return Ok(None);
        };
        let (sin, cos) = theta.sin_cos();
        let plane = self.local.start.x * cos + self.local.start.y * sin;
        let plane = plane / plane.length();
        let (point, error, tangent) = match self.base {
            Base::Banked(banked, _) => {
                let Some((p, e)) = self.table.point2(m) else {
                    return Ok(None);
                };
                // The laws the kernel reads at `m` from `side`: the pieces
                // ending there on the incoming side.
                let truncated;
                let read = if side == SeamSide::Incoming {
                    truncated = crate::station::seam::banked_ending_at(banked, m);
                    &truncated
                } else {
                    banked
                };
                let Ok((pivot, rate)) = read.pivot_at(m) else {
                    return Ok(None);
                };
                let height = crate::elevation::elevation_height(&read.base.elevation, m)?;
                let grade = crate::elevation::elevation_grade(&read.base.elevation, m)?;
                let derivative = Vec3::new(plane.x, plane.y, grade + rate);
                let axes = crate::banked::rolled_axes(read, m, derivative)?;
                let frame = SectionFrame {
                    point: axiolid_core::Point3::new(p.x, p.y, height + pivot),
                    tangent: axes.tangent,
                    lateral: axes.lateral,
                    up: axes.up,
                };
                return Ok(Some((frame, e, axes.roll)));
            }
            Base::Plan(_) => {
                let Some((p, e)) = self.table.point2(m) else {
                    return Ok(None);
                };
                (
                    axiolid_core::Point3::new(p.x, p.y, 0.0),
                    e,
                    Vec3::new(plane.x, plane.y, 0.0),
                )
            }
            Base::Level(_) => {
                let Some((p, e)) = self.table.point3(m) else {
                    return Ok(None);
                };
                (p, e, Vec3::new(plane.x, plane.y, 0.0))
            }
            Base::Elevated(elevated, _) => {
                let Some((p, e)) = self.table.point2(m) else {
                    return Ok(None);
                };
                let height = crate::elevation::elevation_height(&elevated.elevation, m)?;
                // A kink of the profile has two grades; the read side's.
                let grade = if side == SeamSide::Incoming {
                    crate::elevation::elevation_grade_before(&elevated.elevation, m)?
                } else {
                    crate::elevation::elevation_grade(&elevated.elevation, m)?
                };
                let size = (1.0 + grade * grade).sqrt();
                (
                    axiolid_core::Point3::new(p.x, p.y, height),
                    e,
                    Vec3::new(plane.x / size, plane.y / size, grade / size),
                )
            }
        };
        Ok(Some((framed(point, tangent)?, error, 0.0)))
    }

    /// The charge of the read `frame` at the distance `u` into `piece`
    /// from `side`: its point's distance from the exact point plus the
    /// table's error and the points' rounding, and the largest distance of
    /// an axis from the exact axis plus the heading's rounding (condition
    /// one for an orthonormal start frame, and the lateral and up axes of
    /// a leaning tangent amplified by the horizontal part it keeps).
    fn compare(
        &self,
        piece: &StationPiece<'c>,
        frame: &SectionFrame,
        u: Scalar,
        side: SeamSide,
    ) -> GeomResult<Option<[Scalar; 2]>> {
        let (m, curve_side) = piece.read_at(u, side);
        let Some((exact, table, roll)) = self.exact(m, curve_side)? else {
            return Ok(None);
        };
        let scale = exact.point.abs().max_element().max(1.0);
        let point_error = (frame.point - exact.point).length() + table + POINT_ROUNDING * scale;
        let r_t = self.table.heading_rounding() * (1.0 + 4.0 * NEAR) + AXIS_ROUNDING;
        let (r_l, r_u) = if self.is_plan_measured() {
            let flat = exact.tangent.x.hypot(exact.tangent.y) - r_t;
            if flat.is_nan() || flat <= 0.0 {
                return Ok(None);
            }
            let r_n = 2.0 * r_t / flat;
            let r_u = r_n + r_t;
            if self.banked().is_some() {
                // The horizontal normal and the square up carry the heading's
                // error; the roll, from the cant law at the same distance,
                // only its rounding, through `asin`. Each rolled axis mixes
                // the two.
                let cosine = roll.cos();
                if cosine.is_nan() || cosine <= 0.0 {
                    return Ok(None);
                }
                let r_roll = 16.0 * Scalar::EPSILON / cosine + r_n + r_u;
                (r_roll, r_roll)
            } else {
                (r_n, r_u)
            }
        } else {
            (r_t, 0.0)
        };
        let axis_error = ((frame.tangent - exact.tangent).length() + r_t)
            .max((frame.lateral - exact.lateral).length() + r_l)
            .max((frame.up - exact.up).length() + r_u)
            + AXIS_ROUNDING;
        Ok(Some([point_error, axis_error]))
    }

    /// How far the frames `piece` reads at `from` (from the outgoing side)
    /// and at `to` (the incoming) lie from the exact frames: `[point_error,
    /// axis_error]`, the larger of both ends; `None` where the exact
    /// frame is not known well enough.
    ///
    /// # Errors
    ///
    /// The evaluators' refusals.
    pub(crate) fn errors(
        &self,
        piece: &StationPiece<'c>,
        from: Scalar,
        to: Scalar,
    ) -> GeomResult<Option<[Scalar; 2]>> {
        let mut out: [Scalar; 2] = [0.0; 2];
        for (u, side) in [(from, SeamSide::Outgoing), (to, SeamSide::Incoming)] {
            let (_, frame) = piece.raw_section(u, side, None)?;
            let Some(charge) = self.compare(piece, &frame, u, side)? else {
                return Ok(None);
            };
            out = [out[0].max(charge[0]), out[1].max(charge[1])];
        }
        Ok(out.iter().all(|v| v.is_finite()).then_some(out))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::{Point2, Vec2};
    use axiolid_curve::CurvatureLaw;

    fn clothoid() -> Curve2 {
        Curve2::Intrinsic(Intrinsic2::new(
            Frame2 {
                origin: Point2::new(2.0, 1.0),
                x: Vec2::new(0.6, 0.8),
                y: Vec2::new(-0.8, 0.6),
            },
            CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
            60.0,
        ))
    }

    /// A read off the exact frame by a known amount is charged at least
    /// that much: its point, and each of its axes.
    #[test]
    fn a_read_off_the_exact_frame_is_charged_its_distance() {
        let curve = clothoid();
        let piece = StationPiece::between(StationCurve::Two(&curve), 0.0, 60.0).unwrap();
        let reader = IntrinsicReader::new(&piece, 1e-7).unwrap().unwrap();
        let u = 41.0;
        let side = SeamSide::Outgoing;
        let (_, frame) = piece.raw_section(u, side, None).unwrap();
        let off = 1e-6;
        let charge =
            |frame: SectionFrame| reader.compare(&piece, &frame, u, side).unwrap().unwrap();
        let [point, axis] = charge(frame);
        assert!(point <= 1e-7 && axis <= 1e-12, "{point} {axis}");
        let moved = SectionFrame {
            point: frame.point + off * Vec3::Y,
            ..frame
        };
        assert!(charge(moved)[0] >= off);
        let turned = |axis: fn(&mut SectionFrame) -> &mut Vec3| {
            let mut out = frame;
            *axis(&mut out) += off * Vec3::Z;
            charge(out)[1]
        };
        assert!(turned(|f| &mut f.tangent) >= off);
        assert!(turned(|f| &mut f.lateral) >= off);
        assert!(turned(|f| &mut f.up) >= off);
    }

    /// A banked read off the exact frame by a known amount is charged at
    /// least that much, its rolled axes included, and a read on the exact
    /// frame only its rounding.
    #[test]
    fn a_banked_read_off_the_exact_frame_is_charged_its_distance() {
        use axiolid_curve::{BankConvention, CantLaw, CantPiece, ElevationLaw};
        let curve = Curve3::Banked(Banked3::new(
            Elevated3::new(clothoid(), ElevationLaw::constant_grade(5.0, 0.03)),
            CantLaw::new(vec![CantPiece::linear(60.0, 0.0, 0.12)]),
            CantLaw::new(vec![CantPiece::linear(60.0, 0.0, 0.06)]),
            1.5,
            BankConvention::TangentRotation,
        ));
        let piece = StationPiece::between(StationCurve::Three(&curve), 0.0, 60.0).unwrap();
        let reader = IntrinsicReader::new(&piece, 1e-7).unwrap().unwrap();
        assert!(reader.banked().is_some() && reader.is_plan_measured());
        let u = 41.0;
        let side = SeamSide::Outgoing;
        let (_, frame) = piece.raw_section(u, side, None).unwrap();
        let off = 1e-6;
        let charge =
            |frame: SectionFrame| reader.compare(&piece, &frame, u, side).unwrap().unwrap();
        let [point, axis] = charge(frame);
        assert!(point <= 1e-7 && axis <= 1e-12, "{point} {axis}");
        let moved = SectionFrame {
            point: frame.point + off * Vec3::Z,
            ..frame
        };
        assert!(charge(moved)[0] >= off);
        let turned = |axis: fn(&mut SectionFrame) -> &mut Vec3| {
            let mut out = frame;
            *axis(&mut out) += off * Vec3::Z;
            charge(out)[1]
        };
        assert!(turned(|f| &mut f.tangent) >= off);
        assert!(turned(|f| &mut f.lateral) >= off);
        assert!(turned(|f| &mut f.up) >= off);
        // The exact point is raised by the height and the pivot, as read.
        let (exact, _, roll) = reader.exact(u, side).unwrap().unwrap();
        assert!((exact.point.z - frame.point.z).abs() <= 1e-9 && roll > 0.0);
    }
}
