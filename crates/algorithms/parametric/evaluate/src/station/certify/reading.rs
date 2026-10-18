//! Frames read by quadrature along a piece, certified (#316, ADR 0082
//! amendment).
//!
//! On an ellipse or a B-spline, and on an offset beside a curve in closed
//! form (a line, a circle, an ellipse or a B-spline), the kernel reads a
//! frame at a distance through estimates: the arc-length inverse, an
//! offset's own length table, and an offset's tangent by a difference
//! quotient of its displacement. The frame law these bounds are about is
//! the EXACT one: the curve the length certificates of this module measure
//! (for an offset by distances, its law runs linearly in the base's true
//! arc length between the points the kernel reads as the piece's ends, the
//! `beside` module's curve), framed by its own tangent, at its true arc
//! length from the piece's start.
//!
//! A frame the kernel reads at a distance `u` into the piece stands at a
//! native parameter `t`. There, the exact curve's point `P(t)` and
//! velocity `P'(t)` are known in closed form to a certified error (the
//! `panels::Measured` reading the length certificate rests on), so the
//! read frame's distance from the exact frame at `t` -- its point, and
//! each axis from the exact axes framed as the kernel frames them -- is
//! computed and bounded, whatever estimate produced it. The exact frame at
//! `t` stands at the exact curve's true arc length from the piece's start
//! to `t`, which a table of chord panels over the piece bounds
//! (`beside::BaseTable`); its distance from `u` is the read's measure.
//!
//! Over a span between two reads, at `t_a` and `t_b`, the exact frame law
//! is bounded in its own arc length `s` from the suprema `M2 >= |P''|` and
//! `M3 >= |P'''|` over `[t_a, t_b]` and the floor `sigma` of the speed
//! there (the midpoint's, less its slip and `h M2 / 2`). The unit tangent
//! `u = P' / |P'|` has `|u_t| <= x = M2 / sigma` and `|u_tt| <= M3 /
//! sigma + 3 x^2` (`panels::unit_jet`), and `dt/ds = 1 / |P'|` has a
//! derivative of at most `M2 / sigma^3`; so `|O'| = 1`, `|T'| = |O''| <=
//! M2 / sigma^2` and `|T''| <= |u_tt| / sigma^2 + M2^2 / sigma^4`. The axes
//! follow from the tangent as on every other curve
//! (`super::super::frame_bound`).
//!
//! At a B-spline's corner knot the kernel reads the side it names in
//! closed form at the knot's parameter; there the read axes are the exact
//! ones when the read point is the knot's, rounding aside.

use axiolid_contracts::GeomResult;
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_curve::{OffsetLaw, SeamSide};

use super::beside::{BaseTable, Beside};
use super::panels::{Measured, PANEL_SHARE};
use super::{Native, PlanModel, POINT_ROUNDING};
use crate::station::composite::{StationCurve, StationPiece};
use crate::station::offset::{OffsetTable, StationOffset};
use crate::station::{reference_up_frame, SectionFrame, KEEPS_UP_TOLERANCE};

/// The rounding allowed a computed unit axis.
const AXIS_ROUNDING: Scalar = 1e-14;

/// What a piece's frames are read along.
enum Source<'c> {
    /// An ellipse or a B-spline at its native parameter.
    Curve(Box<Native<'c>>),
    /// An offset of a closed-form base.
    Offset(Box<(StationOffset<'c>, Beside<'c>)>),
}

/// The certified frame law of a piece over a span, as the
/// [module documentation](self) derives it: `[|O'|, |O''|]`, `[|T'|,
/// |T''|]` in the exact curve's arc length, the unit tangent at the
/// middle (in the curve's own direction), whether every tangent is
/// horizontal, how far a read's true distance lies from the one asked
/// for, and how far a read point and a read axis lie from the exact ones.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ReadLaw {
    pub(crate) point: [Scalar; 2],
    pub(crate) tangent: [Scalar; 2],
    pub(crate) middle: Vec3,
    pub(crate) planar: bool,
    pub(crate) measure: Scalar,
    pub(crate) point_error: Scalar,
    pub(crate) axis_error: Scalar,
}

/// A piece whose frames are read by quadrature, with its own length
/// measured once by chord panels.
pub(crate) struct PieceReader<'c> {
    source: Source<'c>,
    table: BaseTable,
    /// Whether the distance into the piece grows with the native parameter.
    increasing: bool,
    /// Whether every tangent of the exact curve is horizontal.
    planar: bool,
    /// Native parameters where the curve's tangent may jump.
    corners: Vec<Scalar>,
}

/// The exact frame at `point` with unit tangent `tangent`, as the kernel
/// frames it: the planar frame where the tangent is horizontal, the
/// reference-up frame elsewhere.
fn framed(point: Point3, tangent: Vec3) -> GeomResult<SectionFrame> {
    if tangent.z == 0.0 {
        return Ok(SectionFrame {
            point,
            tangent,
            lateral: Vec3::new(-tangent.y, tangent.x, 0.0),
            up: Vec3::Z,
        });
    }
    reference_up_frame(point, tangent)
}

impl<'c> PieceReader<'c> {
    /// The reader of `piece` (an offset's through `table`), its length
    /// measured to `tolerance`; `None` for a piece whose frames are not
    /// read by quadrature here: one read from stored data, an offset of a
    /// base parameterised by its plan length or of a banked curve, a
    /// placed offset, and a curve placed in a frame that tilts `+Z`.
    ///
    /// # Errors
    ///
    /// The certificate's own refusals.
    pub(crate) fn new(
        piece: &StationPiece<'c>,
        table: Option<&OffsetTable>,
        tolerance: Scalar,
    ) -> GeomResult<Option<Self>> {
        let target = PANEL_SHARE * tolerance;
        match piece.curve {
            StationCurve::Offset(offset) => {
                if piece.placement.is_some() || !matches!(PlanModel::of(&offset), Ok(None)) {
                    return Ok(None);
                }
                let Some(table) = table else {
                    return Ok(None);
                };
                let a = offset.base_measure(piece.start, table)?;
                let b = offset.base_measure(piece.end, table)?;
                let Some(beside) = Beside::new(&offset, a.min(b), a.max(b), tolerance)? else {
                    return Ok(None);
                };
                let (t0, t1) = (beside.t0, beside.t1);
                let lengths = beside
                    .measured(|curve| BaseTable::build(curve, t0, t1, &beside.breaks, target))?;
                let base = offset.base();
                let level = match offset.law() {
                    OffsetLaw::Planar { .. } => true,
                    OffsetLaw::Linear { start, end, .. } => start.vertical == end.vertical,
                    _ => false,
                };
                Ok(Some(Self {
                    planar: matches!(base.curve, StationCurve::Two(_)) && level,
                    increasing: base.reversed == piece.reversed,
                    source: Source::Offset(Box::new((offset, beside))),
                    table: lengths,
                    corners: Vec::new(),
                }))
            }
            curve if curve.measure_is_exact() => Ok(None),
            curve => {
                if let Some(rigid) = piece.placement {
                    if (rigid.transform_vector3(Vec3::Z) - Vec3::Z).length() > KEEPS_UP_TOLERANCE {
                        return Ok(None);
                    }
                }
                let native = Native::of(curve, piece.start, piece.end)?;
                let (ta, tb) = (native.parameter(piece.start)?, native.parameter(piece.end)?);
                let (t0, t1) = (ta.min(tb), ta.max(tb));
                let corners = native.breaks(1);
                let mut breaks = native.breaks(2);
                breaks.extend(&corners);
                let lengths =
                    native.measured(|curve| BaseTable::build(curve, t0, t1, &breaks, target))?;
                Ok(Some(Self {
                    planar: matches!(curve, StationCurve::Two(_)),
                    increasing: !piece.reversed,
                    source: Source::Curve(Box::new(native)),
                    table: lengths,
                    corners,
                }))
            }
        }
    }

    fn measured<R>(&self, run: impl FnOnce(&Measured<'_>) -> GeomResult<R>) -> GeomResult<R> {
        match &self.source {
            Source::Curve(native) => native.measured(run),
            Source::Offset(offset) => offset.1.measured(run),
        }
    }

    /// The direction the piece's measure runs along the native parameter.
    fn travel(&self) -> Scalar {
        match &self.source {
            Source::Curve(_) => 1.0,
            Source::Offset(offset) => offset.1.travel(),
        }
    }

    /// The native parameter the curve's measure `m` is read at.
    fn parameter(&self, m: Scalar, table: Option<&OffsetTable>) -> GeomResult<Option<Scalar>> {
        match &self.source {
            Source::Curve(native) => native.parameter(m).map(Some),
            Source::Offset(boxed) => {
                let (offset, beside) = boxed.as_ref();
                let Some(table) = table else {
                    return Ok(None);
                };
                let v = offset.base_measure(m, table)?;
                let base = offset.base();
                let mb = if base.reversed {
                    base.end - v
                } else {
                    base.start + v
                };
                beside.native.parameter(mb).map(Some)
            }
        }
    }

    /// The certified frame law of `piece` between the distances `from <
    /// to` into it, read from inside the span (an offset through
    /// `table`); `None` where it is not bounded.
    ///
    /// # Errors
    ///
    /// The evaluators' refusals.
    pub(crate) fn read(
        &self,
        piece: &StationPiece<'c>,
        from: Scalar,
        to: Scalar,
        table: Option<&OffsetTable>,
    ) -> GeomResult<Option<ReadLaw>> {
        let travel = self.travel();
        let mut reads = Vec::with_capacity(2);
        for (u, side) in [(from, SeamSide::Outgoing), (to, SeamSide::Incoming)] {
            let (m, frame) = piece.raw_section(u, side, table)?;
            let Some(t) = self.parameter(m, table)? else {
                return Ok(None);
            };
            reads.push((u, frame, t));
        }
        let (ta, tb) = (reads[0].2, reads[1].2);
        self.measured(|curve| {
            let mut out = ReadLaw {
                point: [0.0; 2],
                tangent: [0.0; 2],
                middle: Vec3::ZERO,
                planar: self.planar,
                measure: 0.0,
                point_error: 0.0,
                axis_error: 0.0,
            };
            for (u, read, t) in &reads {
                let Some((measure, point_error, axis_error)) = self.compare(curve, read, *u, *t)?
                else {
                    return Ok(None);
                };
                out.measure = out.measure.max(measure);
                out.point_error = out.point_error.max(point_error);
                out.axis_error = out.axis_error.max(axis_error);
            }
            let (lo, hi) = (ta.min(tb), ta.max(tb));
            let (second, third) = (curve.bounds)(lo, hi)?;
            let Some(third) = third.filter(|m| m.is_finite()) else {
                return Ok(None);
            };
            let (velocity, slip) = (curve.velocity)(0.5 * (lo + hi))?;
            let sigma = velocity.length() - slip - 0.5 * (hi - lo) * second;
            if !(sigma.is_finite() && sigma > 0.0) {
                return Ok(None);
            }
            let x = second / sigma;
            let unit_second = third / sigma + 3.0 * x * x;
            let k1 = second / (sigma * sigma);
            let k2 = unit_second / (sigma * sigma) + x * x / (sigma * sigma);
            out.point = [1.0, k1];
            out.tangent = [k1, k2];
            out.middle = (travel * velocity) / velocity.length();
            let finite = [k1, k2, out.measure, out.point_error, out.axis_error]
                .iter()
                .all(|v| v.is_finite());
            Ok(finite.then_some(out))
        })
    }

    /// The read `frame` at the distance `u` into the piece, standing at the
    /// native parameter `t`, against the exact frame there: the true
    /// distance's distance from `u`, the point's error and the largest
    /// axis's; `None` where the exact tangent is not known well enough.
    fn compare(
        &self,
        curve: &Measured<'_>,
        frame: &SectionFrame,
        u: Scalar,
        t: Scalar,
    ) -> GeomResult<Option<(Scalar, Scalar, Scalar)>> {
        let (point, error) = (curve.point)(t)?;
        let (velocity, slip) = (curve.velocity)(t)?;
        let speed = velocity.length();
        if !(speed.is_finite() && speed > 2.0 * slip) {
            return Ok(None);
        }
        let tangent = (self.travel() * velocity) / speed;
        let exact = framed(point, tangent)?;
        let scale = curve.scale.max(point.abs().max_element()).max(1.0);
        let point_error = (frame.point - point).length() + error + POINT_ROUNDING * scale;
        let corner =
            self.corners.contains(&t) && (frame.point - point).length() <= POINT_ROUNDING * scale;
        let axis_error = if corner {
            AXIS_ROUNDING
        } else {
            let r_t = 2.0 * slip / speed + AXIS_ROUNDING;
            let (r_l, r_u) = if self.planar {
                (r_t, 0.0)
            } else {
                let flat = tangent.x.hypot(tangent.y) - r_t;
                if flat.is_nan() || flat <= 0.0 {
                    return Ok(None);
                }
                let r_l = 2.0 * r_t / flat;
                (r_l, r_l + r_t)
            };
            ((frame.tangent - exact.tangent).length() + r_t)
                .max((frame.lateral - exact.lateral).length() + r_l)
                .max((frame.up - exact.up).length() + r_u)
                + AXIS_ROUNDING
        };
        let (lo, hi) = self.table.at(curve, t)?;
        let (into_lo, into_hi) = if self.increasing {
            (lo, hi)
        } else {
            let (whole_lo, whole_hi) = self.table.total();
            (whole_lo - hi, whole_hi - lo)
        };
        let measure = (u - into_lo).max(into_hi - u).max(0.0);
        Ok(Some((measure, point_error, axis_error)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::{Frame2, Point2, Vec2};
    use axiolid_curve::{Curve2, Ellipse2};

    /// A read off the exact frame by a known amount is charged at least
    /// that much: its point, and each of its axes.
    #[test]
    fn a_read_off_the_exact_frame_is_charged_its_distance() {
        let curve = Curve2::Ellipse(Ellipse2 {
            frame: Frame2 {
                origin: Point2::new(30.0, 0.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            semi_axis_x: 10.0,
            semi_axis_y: 6.0,
        });
        let piece = StationPiece::between(StationCurve::Two(&curve), 0.0, 20.0).unwrap();
        let reader = PieceReader::new(&piece, None, 1e-7).unwrap().unwrap();
        let u = 7.0;
        let (m, frame) = piece.raw_section(u, SeamSide::Outgoing, None).unwrap();
        let t = reader.parameter(m, None).unwrap().unwrap();
        let off = 1e-6;
        let charge = |frame: SectionFrame| {
            reader
                .measured(|curve| reader.compare(curve, &frame, u, t))
                .unwrap()
                .unwrap()
        };
        let (_, point, axis) = charge(frame);
        assert!(point <= 1e-12 && axis <= 1e-8, "{point} {axis}");
        let moved = SectionFrame {
            point: frame.point + off * Vec3::Z,
            ..frame
        };
        assert!(charge(moved).1 >= off);
        let turned = |axis: fn(&mut SectionFrame) -> &mut Vec3| {
            let mut out = frame;
            *axis(&mut out) += off * Vec3::Z;
            charge(out).2
        };
        assert!(turned(|f| &mut f.tangent) >= off);
        assert!(turned(|f| &mut f.lateral) >= off);
        assert!(turned(|f| &mut f.up) >= off);
    }
}
