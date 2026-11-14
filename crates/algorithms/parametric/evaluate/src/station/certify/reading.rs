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

use std::cell::RefCell;
use std::collections::HashMap;

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
pub(crate) const AXIS_ROUNDING: Scalar = 1e-14;

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

/// The exact frame law's jet at the native parameter a read stands at, in
/// the piece's own measure and direction (#406): the unit tangent `O'` and
/// the curvature vector `O'' = (P'' - (P''.V) V) / |P'|^2`, `V = P' / |P'|`,
/// each with a certified error.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ReadJet {
    pub(crate) tangent: Vec3,
    pub(crate) curvature: Vec3,
    pub(crate) tangent_error: Scalar,
    pub(crate) curvature_error: Scalar,
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
    /// Native parameters where its curvature may jump: `P''` is continuous
    /// between them (#406).
    bend_breaks: Vec<Scalar>,
    /// What a read at the distance `u` into the piece from a side stood at
    /// and was charged, by `u`'s bits and whether the side is the outgoing
    /// one: a span and its halves share their ends, and each read costs an
    /// arc-length inverse and a frame (#406).
    ends: RefCell<HashMap<(u64, bool), Option<End>>>,
    /// The native parameter a measure is read at, by the measure's bits.
    parameters: RefCell<HashMap<u64, Option<Scalar>>>,
    /// The native parameters the piece spans, ascending.
    range: (Scalar, Scalar),
}

/// A read's native parameter, its measure and the errors of its point and
/// axes.
type End = (Scalar, Scalar, Scalar, Scalar);

/// The exact frame at `point` with unit tangent `tangent`, as the kernel
/// frames it: the planar frame where the tangent is horizontal, the
/// reference-up frame elsewhere.
pub(super) fn framed(point: Point3, tangent: Vec3) -> GeomResult<SectionFrame> {
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
                    bend_breaks: beside.breaks.clone(),
                    ends: RefCell::default(),
                    parameters: RefCell::default(),
                    range: (t0, t1),
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
                    bend_breaks: breaks,
                    ends: RefCell::default(),
                    parameters: RefCell::default(),
                    range: (t0, t1),
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

    /// The native parameter the curve's measure `m` is read at: an
    /// arc-length inverse, kept for the span and its halves that share it.
    fn parameter(&self, m: Scalar, table: Option<&OffsetTable>) -> GeomResult<Option<Scalar>> {
        let known = self.parameters.borrow().get(&m.to_bits()).copied();
        if let Some(t) = known {
            return Ok(t);
        }
        let t = self.uncached_parameter(m, table)?;
        self.parameters.borrow_mut().insert(m.to_bits(), t);
        Ok(t)
    }

    fn uncached_parameter(
        &self,
        m: Scalar,
        table: Option<&OffsetTable>,
    ) -> GeomResult<Option<Scalar>> {
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

    /// Native parameters where the curve's curvature may jump: the exact
    /// law's curvature is continuous between them (#406).
    pub(crate) fn bend_breaks(&self) -> &[Scalar] {
        &self.bend_breaks
    }

    /// The native parameters the reads at the distances `from` (read from
    /// the outgoing side) and `to` (from the incoming one) into the piece
    /// stand at; `None` where the offset's length table is missing (#406).
    ///
    /// # Errors
    ///
    /// The evaluators' refusals.
    pub(crate) fn span_parameters(
        &self,
        piece: &StationPiece<'c>,
        from: Scalar,
        to: Scalar,
        table: Option<&OffsetTable>,
    ) -> GeomResult<Option<(Scalar, Scalar)>> {
        let mut out = [0.0; 2];
        for (slot, (u, side)) in out
            .iter_mut()
            .zip([(from, SeamSide::Outgoing), (to, SeamSide::Incoming)])
        {
            let (m, _) = piece.read_at(u, side);
            let Some(t) = self.parameter(m, table)? else {
                return Ok(None);
            };
            *slot = t;
        }
        Ok(Some((out[0], out[1])))
    }

    /// The exact frame law's jet at the distance `u` into the piece, read
    /// from the outgoing side as a frame is, at the native parameter that
    /// read stands at (#406); `None` where the speed is not certified
    /// positive.
    ///
    /// # Errors
    ///
    /// The evaluators' refusals.
    pub(crate) fn jet(
        &self,
        piece: &StationPiece<'c>,
        u: Scalar,
        table: Option<&OffsetTable>,
    ) -> GeomResult<Option<ReadJet>> {
        let (m, _) = piece.read_at(u, SeamSide::Outgoing);
        let Some(t) = self.parameter(m, table)? else {
            return Ok(None);
        };
        self.measured(|curve| {
            let (velocity, slip) = (curve.velocity)(t)?;
            let (acceleration, accel_error) = (curve.acceleration)(t)?;
            let speed = velocity.length();
            if !(speed.is_finite() && speed > 2.0 * slip) {
                return Ok(None);
            }
            let low = speed - slip;
            let unit = velocity / speed;
            let size = acceleration.length();
            let tangent_error = 2.0 * slip / low;
            let rounding = 64.0 * Scalar::EPSILON;
            let curvature_error =
                ((accel_error + 2.0 * size * tangent_error + 2.0 * accel_error * tangent_error)
                    / (low * low)
                    + 3.0 * size * slip / (low * low * low))
                    * (1.0 + rounding)
                    + rounding * size / (low * low);
            let out = ReadJet {
                tangent: self.travel() * unit,
                curvature: (acceleration - unit * acceleration.dot(unit)) / (speed * speed),
                tangent_error,
                curvature_error,
            };
            let finite = out.tangent.is_finite()
                && out.curvature.is_finite()
                && out.tangent_error.is_finite()
                && out.curvature_error.is_finite();
            Ok(finite.then_some(out))
        })
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
        let mut reads = Vec::with_capacity(2);
        for (u, side) in [(from, SeamSide::Outgoing), (to, SeamSide::Incoming)] {
            let key = (u.to_bits(), side == SeamSide::Outgoing);
            let known = self.ends.borrow().get(&key).copied();
            let end = match known {
                Some(end) => end,
                None => {
                    let (m, frame) = piece.raw_section(u, side, table)?;
                    let end = match self.parameter(m, table)? {
                        Some(t) => self
                            .measured(|curve| self.compare(curve, &frame, u, t))?
                            .map(|(measure, point, axis)| (t, measure, point, axis)),
                        None => None,
                    };
                    self.ends.borrow_mut().insert(key, end);
                    end
                }
            };
            let Some(end) = end else {
                return Ok(None);
            };
            reads.push(end);
        }
        let (ta, tb) = (reads[0].0, reads[1].0);
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
            for &(_, measure, point_error, axis_error) in &reads {
                out.measure = out.measure.max(measure);
                out.point_error = out.point_error.max(point_error);
                out.axis_error = out.axis_error.max(axis_error);
            }
            let Some(law) = self.constants(curve, ta.min(tb), ta.max(tb))? else {
                return Ok(None);
            };
            out.point = law.point;
            out.tangent = law.tangent;
            out.middle = law.middle;
            let finite = [out.measure, out.point_error, out.axis_error]
                .iter()
                .chain(&out.point)
                .chain(&out.tangent)
                .all(|v| v.is_finite());
            Ok(finite.then_some(out))
        })
    }

    /// The frame law's constants over the native parameters `[lo, hi]`, with
    /// no read's errors: `[|O'|, |O''|]`, `[|T'|, |T''|]` in the exact
    /// curve's arc length and the unit tangent at the middle.
    fn constants(
        &self,
        curve: &Measured<'_>,
        lo: Scalar,
        hi: Scalar,
    ) -> GeomResult<Option<ReadLaw>> {
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
        Ok(Some(ReadLaw {
            point: [1.0, k1],
            tangent: [k1, k2],
            middle: (self.travel() * velocity) / velocity.length(),
            planar: self.planar,
            measure: 0.0,
            point_error: 0.0,
            axis_error: 0.0,
        }))
    }

    /// The frame law's constants over every true distance of a span read
    /// between `from` and `to` into the piece, each read standing within
    /// `measure` of its distance, without reading a frame at their ends
    /// (#406): the native parameters of the two reads, padded by `2 measure`
    /// of exact arc length, which at a speed of at least half the floor over
    /// the span is `4 measure / sigma` of parameter. `None` where the speed
    /// floor does not hold over the padded range or the reads have no
    /// parameters or the padded range strictly contains a point where its
    /// curvature may jump.
    ///
    /// # Errors
    ///
    /// The evaluators' refusals.
    pub(crate) fn read_padded(
        &self,
        piece: &StationPiece<'c>,
        (from, to): (Scalar, Scalar),
        table: Option<&OffsetTable>,
        measure: Scalar,
    ) -> GeomResult<Option<ReadLaw>> {
        let Some((ta, tb)) = self.span_parameters(piece, from, to, table)? else {
            return Ok(None);
        };
        let (lo, hi) = (ta.min(tb), ta.max(tb));
        self.measured(|curve| {
            // The speed floor over a range, as `constants` takes it.
            let floor = |lo: Scalar, hi: Scalar| -> GeomResult<Scalar> {
                let (velocity, slip) = (curve.velocity)(0.5 * (lo + hi))?;
                let (second, _) = (curve.bounds)(lo, hi)?;
                Ok(velocity.length() - slip - 0.5 * (hi - lo) * second)
            };
            let sigma = floor(lo, hi)?;
            if !(sigma.is_finite() && sigma > 0.0) {
                return Ok(None);
            }
            let pad = 4.0 * measure / sigma;
            let (start, end) = self.range;
            let (lo, hi) = ((lo - pad).max(start), (hi + pad).min(end));
            // The padding assumed half the floor of the unpadded range.
            if floor(lo, hi)? < 0.5 * sigma {
                return Ok(None);
            }
            // Its Lipschitz constants hold where `P''` is continuous only.
            if self.bend_breaks.iter().any(|&k| lo < k && k < hi) {
                return Ok(None);
            }
            self.constants(curve, lo, hi)
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
