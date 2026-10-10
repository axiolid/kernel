//! Offset pieces of a station basis: the curve beside a span of an atomic
//! curve, measured in its own length (#289, ADR 0082 amendment).
//!
//! A [`StationOffset`] is the curve `C(v) = B(v) + D(v)` beside a base
//! piece `B` (a [`StationPiece`] of an atomic curve, reversed or placed),
//! `v` the base's station measure into the piece and `D` the displacement
//! its [`OffsetLaw`] gives in the base's frame at `v`:
//!
//! - [`OffsetLaw::Planar`]: `d` along the base's left lateral, on a 2D
//!   base (an offset curve 2D: positive to the left, the orthogonal
//!   complement of the tangent);
//! - [`OffsetLaw::Directed`]: `d` along `normalise(V x T)`, `T` the base's
//!   unit tangent (an offset curve 3D); a tangent parallel to `V` is
//!   refused by name;
//! - [`OffsetLaw::Linear`]: lateral, vertical and longitudinal offsets
//!   interpolated linearly in `v` between the piece's ends, in the base's
//!   section or plan frame (one interval of an offset by distances).
//!
//! # Distance convention
//!
//! A station on an offset is measured in the OFFSET's own length from its
//! start (the base piece's start), not in the base's: its arc length, or,
//! where the base is plan-measured (an elevated or banked curve), its own
//! plan length (the arc length of its plan projection), so an offset of an
//! alignment keeps the alignment's convention and joins plan-measured
//! pieces. The offset's frame is its own section frame: the planar frame
//! of its own tangent where it runs horizontally in a horizontal plane,
//! else the reference-up frame against `+Z` of its own point and tangent;
//! an offset of a banked curve keeps the base's rolled lateral, made
//! perpendicular to its own tangent.
//!
//! # Exactness
//!
//! A constant offset of a line is a line beside it, and a constant offset
//! (no longitudinal part) of a circle whose offset stays in planes normal
//! to its axis (a 2D or horizontal circle under a planar, section or plan
//! law, any circle under a reference direction along its axis) is a
//! circle of radius `r - a`, `a` the displacement towards the centre: its
//! length is the base's scaled by `(r - a) / r`, read in closed form.
//! Every other offset is read numerically: its speed is the base's
//! tangent plus the derivative of the displacement by a five-point
//! difference, integrated by adaptive Gauss-Legendre to
//! [`OFFSET_TOLERANCE`] and inverted by a safeguarded Newton step within a
//! panel. That is an estimate, not a proof: a frame on such an offset is
//! never exact, and a mesh placed in it is `Unbounded` by name.
//!
//! # Refusals
//!
//! By name: a base that is itself an offset, a base span with a seam of
//! its curve inside it (split it there, as [`offset_pieces`] does), a 2D
//! law on a 3D base, a circle whose offset radius collapses to zero or
//! below, an offset whose speed in the base's direction falls to
//! [`CUSP_TOLERANCE`] of the base's or below (a cusp, or a reversal), and
//! an offset lying in a horizontal plane that crosses itself within the
//! piece.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Transform3, Vec3};
use axiolid_curve::{Curve2, Curve3, OffsetFrame, OffsetLaw, SeamSide};
use axiolid_curve_evaluate_contract::DistanceConvention;

use super::composite::{StationCurve, StationPiece};
use super::{invalid, reference_up_frame, SectionFrame};
use crate::arc_parameter::{NODES, WEIGHTS};

/// Relative tolerance of a numerical offset's own length and its inverse:
/// `OFFSET_TOLERANCE * max(1, L)` in the length unit. An estimate, as the
/// quadrature's agreement test is.
pub const OFFSET_TOLERANCE: Scalar = 1e-9;

/// Smallest speed of an offset in its base's direction, relative to the
/// base's, that is not a cusp; also the smallest scale `(r - a) / r` of an
/// offset circle that has not collapsed.
pub const CUSP_TOLERANCE: Scalar = 1e-6;

/// Most quadrature panels an offset's length may use.
const MAX_OFFSET_PANELS: usize = 1 << 12;

/// Root-find iterations within a panel.
const MAX_ITERATIONS: usize = 100;

fn degenerate(detail: String) -> GeomError {
    GeomError::Degenerate(detail)
}

/// The atomic curve under an offset's base.
#[derive(Debug, Clone, Copy, PartialEq)]
enum BaseCurve<'c> {
    Two(&'c Curve2),
    Three(&'c Curve3),
}

/// An offset of one base piece, measured in its own length; see the
/// [module documentation](self).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StationOffset<'c> {
    curve: BaseCurve<'c>,
    start: Scalar,
    end: Scalar,
    reversed: bool,
    placement: Option<Transform3>,
    placement_exact: bool,
    law: OffsetLaw,
}

/// How an offset's own measure maps to its base's: closed form, or panels.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct OffsetTable {
    /// The speed ratio of a closed form (a line or a circle), else `None`.
    closed: Option<Scalar>,
    /// Panel bounds in the base's measure into the piece, ascending.
    knots: Vec<Scalar>,
    /// The offset's own measure at each knot.
    sums: Vec<Scalar>,
}

impl OffsetTable {
    /// The offset's own length.
    pub(crate) fn length(&self, base: Scalar) -> Scalar {
        match self.closed {
            Some(scale) => scale * base,
            None => self.sums[self.sums.len() - 1],
        }
    }
}

impl<'c> StationOffset<'c> {
    /// The offset of `base` by `law`.
    ///
    /// # Errors
    ///
    /// A base that is itself an offset, a base span with a seam of its
    /// curve strictly inside it, a planar law on a 3D base, a non-finite
    /// law or a zero reference direction, each by name.
    pub fn new(base: StationPiece<'c>, law: OffsetLaw) -> GeomResult<Self> {
        let curve =
            match base.curve {
                StationCurve::Two(curve) => BaseCurve::Two(curve),
                StationCurve::Three(curve) => BaseCurve::Three(curve),
                _ => return Err(invalid(
                    "station: an offset of an offset as a station basis: offset the base once, \
                     by the summed law"
                        .into(),
                )),
            };
        match law {
            OffsetLaw::Planar { distance } => {
                if !distance.is_finite() {
                    return Err(invalid("station: an offset distance is not finite".into()));
                }
                if matches!(curve, BaseCurve::Three(_)) {
                    return Err(invalid(
                        "station: a planar (2D) offset of a 3D curve: a 3D offset needs a \
                         reference direction"
                            .into(),
                    ));
                }
            }
            OffsetLaw::Directed {
                distance,
                reference_direction,
            } => {
                if !distance.is_finite()
                    || !reference_direction.is_finite()
                    || reference_direction.length() <= 1e-12
                {
                    return Err(invalid(
                        "station: a 3D offset's distance or reference direction is not finite, \
                         or the direction is zero"
                            .into(),
                    ));
                }
            }
            OffsetLaw::Linear { start, end, .. } => {
                let finite = |o: axiolid_curve::PathOffsets| {
                    o.lateral.is_finite() && o.vertical.is_finite() && o.longitudinal.is_finite()
                };
                if !(finite(start) && finite(end)) {
                    return Err(invalid(
                        "station: an offset by distances has an offset that is not finite".into(),
                    ));
                }
            }
            _ => {
                return Err(invalid(
                    "station: an offset law this evaluator does not know".into(),
                ))
            }
        }
        if !base.seams()?.is_empty() {
            return Err(invalid(
                "station: an offset piece whose base crosses a seam of its curve: split the base \
                 there, one offset piece per span between seams"
                    .into(),
            ));
        }
        Ok(Self {
            curve,
            start: base.start,
            end: base.end,
            reversed: base.reversed,
            placement: base.placement,
            placement_exact: base.placement_exact,
            law,
        })
    }

    /// The base piece.
    #[must_use]
    pub fn base(&self) -> StationPiece<'c> {
        let curve = match self.curve {
            BaseCurve::Two(curve) => StationCurve::Two(curve),
            BaseCurve::Three(curve) => StationCurve::Three(curve),
        };
        StationPiece {
            curve,
            start: self.start,
            end: self.end,
            reversed: self.reversed,
            placement: self.placement,
            placement_exact: self.placement_exact,
        }
    }

    /// The law.
    #[must_use]
    pub fn law(&self) -> OffsetLaw {
        self.law
    }

    /// The convention the offset is measured in: its base's (plan distance
    /// beside an elevated or banked curve, arc length beside any other).
    #[must_use]
    pub fn convention(&self) -> DistanceConvention {
        self.base().convention()
    }

    /// The offset's own length in its convention, from its start.
    ///
    /// # Errors
    ///
    /// The refusals of the [module documentation](self#refusals), and the
    /// base's evaluators'.
    pub fn length(&self) -> GeomResult<Scalar> {
        Ok(self.table()?.length(self.base_length()))
    }

    fn base_length(&self) -> Scalar {
        self.end - self.start
    }

    fn plan_measured(&self) -> bool {
        self.convention() == DistanceConvention::PlanDistance
    }

    fn banked(&self) -> bool {
        matches!(self.curve, BaseCurve::Three(Curve3::Banked(_)))
    }

    /// Whether the offset is a line beside a line, its frame read exactly.
    pub(crate) fn frame_is_exact(&self) -> bool {
        self.law.is_constant()
            && matches!(
                self.curve,
                BaseCurve::Two(Curve2::Line(_)) | BaseCurve::Three(Curve3::Line(_))
            )
            && self.placement_exact
    }

    /// Whether its measure is read in closed form (a line or a circle).
    pub(crate) fn measure_is_exact(&self) -> bool {
        matches!(self.closed_form(), Ok(Some(_)))
    }

    /// The base's section at `v` into the base piece, read from inside.
    fn base_section(&self, v: Scalar) -> GeomResult<SectionFrame> {
        self.base().section_on(v, SeamSide::Outgoing)
    }

    /// The displacement at `v` from the base's section there.
    fn displacement(&self, section: &SectionFrame, v: Scalar) -> GeomResult<Vec3> {
        match self.law {
            OffsetLaw::Planar { distance } => Ok(distance * section.lateral),
            OffsetLaw::Directed {
                distance,
                reference_direction,
            } => {
                let reference = reference_direction.normalize();
                let normal = reference.cross(section.tangent);
                let size = normal.length();
                if !(size.is_finite() && size > 1e-9) {
                    return Err(invalid(
                        "station: a 3D offset whose base tangent is parallel to its reference \
                         direction: the offset direction is undefined there"
                            .into(),
                    ));
                }
                Ok(distance / size * normal)
            }
            OffsetLaw::Linear { start, end, frame } => {
                let length = self.base_length();
                let at = start.lerp(end, (v / length).clamp(0.0, 1.0));
                let frame = match frame {
                    OffsetFrame::Section => *section,
                    OffsetFrame::Plan => section.plan()?,
                    _ => {
                        return Err(invalid(
                            "station: an offset frame this evaluator does not know".into(),
                        ))
                    }
                };
                Ok(at.lateral * frame.lateral
                    + at.vertical * frame.up
                    + at.longitudinal * frame.tangent)
            }
            _ => Err(invalid(
                "station: an offset law this evaluator does not know".into(),
            )),
        }
    }

    /// The displacement at `v`, its base section read afresh.
    fn displacement_at(&self, v: Scalar) -> GeomResult<Vec3> {
        let section = self.base_section(v)?;
        self.displacement(&section, v)
    }

    /// The step of the displacement's difference quotient.
    fn step(&self) -> Scalar {
        (self.base_length() / 8.0).min(1e-2)
    }

    /// The offset's velocity `dC/dv` at `v`, and the base's section there.
    fn velocity(&self, v: Scalar) -> GeomResult<(Vec3, SectionFrame)> {
        let section = self.base_section(v)?;
        // The base's own rate: unit in arc length, `1 / |T_xy|` in plan
        // distance.
        let rate = if self.plan_measured() {
            let horizontal = Vec3::new(section.tangent.x, section.tangent.y, 0.0).length();
            if horizontal <= 1e-12 {
                return Err(invalid(
                    "station: a plan-measured offset base runs vertically, so it has no plan \
                     rate"
                        .into(),
                ));
            }
            1.0 / horizontal
        } else {
            1.0
        };
        let length = self.base_length();
        let h = self.step();
        let d = |x: Scalar| self.displacement_at(x);
        // Five-point differences, one-sided at the piece's ends so that
        // nothing is read beyond the base piece.
        let derivative = if v - 2.0 * h < 0.0 {
            (-25.0 * self.displacement(&section, v)? + 48.0 * d(v + h)? - 36.0 * d(v + 2.0 * h)?
                + 16.0 * d(v + 3.0 * h)?
                - 3.0 * d(v + 4.0 * h)?)
                / (12.0 * h)
        } else if v + 2.0 * h > length {
            (25.0 * self.displacement(&section, v)? - 48.0 * d(v - h)? + 36.0 * d(v - 2.0 * h)?
                - 16.0 * d(v - 3.0 * h)?
                + 3.0 * d(v - 4.0 * h)?)
                / (12.0 * h)
        } else {
            (-d(v + 2.0 * h)? + 8.0 * d(v + h)? - 8.0 * d(v - h)? + d(v - 2.0 * h)?) / (12.0 * h)
        };
        Ok((rate * section.tangent + derivative, section))
    }

    /// The speed of the offset's own measure along `v`, refusing a cusp.
    fn speed(&self, v: Scalar) -> GeomResult<Scalar> {
        let (velocity, section) = self.velocity(v)?;
        let (along, forward) = if self.plan_measured() {
            let flat = Vec3::new(velocity.x, velocity.y, 0.0);
            let heading = Vec3::new(section.tangent.x, section.tangent.y, 0.0).normalize();
            (flat.length(), flat.dot(heading))
        } else {
            (velocity.length(), velocity.dot(section.tangent))
        };
        if !(forward.is_finite() && forward > CUSP_TOLERANCE) {
            return Err(degenerate(format!(
                "station: an offset with a cusp or turning back at {v} along its base: its speed \
                 in the base's direction is {forward}, the offset reaches the base's radius of \
                 curvature there"
            )));
        }
        Ok(along)
    }

    /// The speed ratio of a closed form: `1` beside a line, `(r - a) / r`
    /// beside a circle; `None` where the offset is read numerically.
    fn closed_form(&self) -> GeomResult<Option<Scalar>> {
        if !self.law.is_constant() {
            return Ok(None);
        }
        let longitudinal = matches!(
            self.law,
            OffsetLaw::Linear { start, .. } if start.longitudinal != 0.0
        );
        let place_point = |p: Point3| self.placement.map_or(p, |r| r.transform_point3(p));
        let place_vector = |v: Vec3| self.placement.map_or(v, |r| r.transform_vector3(v));
        let (centre, axis, radius) = match self.curve {
            BaseCurve::Two(Curve2::Line(_)) | BaseCurve::Three(Curve3::Line(_)) => {
                return Ok(Some(1.0))
            }
            BaseCurve::Two(Curve2::Circle(circle)) => (
                Point3::new(circle.frame.origin.x, circle.frame.origin.y, 0.0),
                Vec3::Z,
                circle.radius,
            ),
            BaseCurve::Three(Curve3::Circle(circle)) => {
                (circle.frame.origin, circle.frame.z, circle.radius)
            }
            _ => return Ok(None),
        };
        if longitudinal {
            return Ok(None);
        }
        let (centre, axis) = (place_point(centre), place_vector(axis).normalize());
        let upright = axis.cross(Vec3::Z).length() <= 1e-12;
        let invariant = match self.law {
            OffsetLaw::Directed {
                reference_direction,
                ..
            } => reference_direction.normalize().cross(axis).length() <= 1e-12,
            _ => upright,
        };
        if !invariant {
            return Ok(None);
        }
        let section = self.base_section(0.0)?;
        let towards = centre - section.point;
        let inward = towards - towards.dot(axis) * axis;
        let inward = inward / inward.length();
        let scale = (radius - self.displacement(&section, 0.0)?.dot(inward)) / radius;
        if !(scale.is_finite() && scale > CUSP_TOLERANCE) {
            return Err(degenerate(format!(
                "station: an offset of a circle of radius {radius} whose radius collapses to {} \
                 (at or below zero): the offset passes through or beyond the centre",
                scale * radius
            )));
        }
        Ok(Some(scale))
    }

    /// The panels of the offset's own measure; see the
    /// [module documentation](self#exactness).
    pub(crate) fn table(&self) -> GeomResult<OffsetTable> {
        if let Some(scale) = self.closed_form()? {
            return Ok(OffsetTable {
                closed: Some(scale),
                knots: Vec::new(),
                sums: Vec::new(),
            });
        }
        let length = self.base_length();
        let rule = |a: Scalar, b: Scalar| -> GeomResult<Scalar> {
            let (half, mid) = (0.5 * (b - a), 0.5 * (a + b));
            let mut sum = 0.0;
            for (node, weight) in NODES.iter().zip(WEIGHTS.iter()) {
                sum += weight * self.speed(mid + half * node)?;
            }
            Ok(sum * half)
        };
        let coarse = rule(0.0, length)?;
        let tolerance = OFFSET_TOLERANCE * coarse.abs().max(1.0);
        let floor = 4.0 * self.step() * 1e-3;
        let mut knots = vec![0.0];
        let mut sums = vec![0.0];
        let mut stack = vec![(0.0, length, coarse)];
        while let Some((lo, hi, whole)) = stack.pop() {
            if knots.len() > MAX_OFFSET_PANELS {
                return Err(GeomError::BudgetExceeded {
                    resource: "offset length quadrature panels",
                });
            }
            let mid = 0.5 * (lo + hi);
            let (left, right) = (rule(lo, mid)?, rule(mid, hi)?);
            let share = tolerance * (hi - lo) / length;
            if (left + right - whole).abs() <= share || hi - lo <= floor || !(mid > lo && mid < hi)
            {
                let run = sums[sums.len() - 1];
                knots.extend([mid, hi]);
                sums.extend([run + left, run + left + right]);
            } else {
                stack.push((mid, hi, right));
                stack.push((lo, mid, left));
            }
        }
        let table = OffsetTable {
            closed: None,
            knots,
            sums,
        };
        self.refuse_self_crossing(&table)?;
        Ok(table)
    }

    /// Refuse an offset lying in one horizontal plane that crosses itself:
    /// its points at the panel bounds and quarter points, as a polyline.
    fn refuse_self_crossing(&self, table: &OffsetTable) -> GeomResult<()> {
        let mut points = Vec::new();
        for pair in table.knots.windows(2) {
            for k in 0..4 {
                let v = pair[0] + (pair[1] - pair[0]) * Scalar::from(k) / 4.0;
                points.push(self.point_at_base(v)?);
            }
        }
        points.push(self.point_at_base(self.base_length())?);
        let z = points[0].z;
        if points.iter().any(|p| p.z != z) {
            return Ok(());
        }
        let cross =
            |a: Point3, b: Point3, c: Point3| (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        let n = points.len();
        for i in 0..n.saturating_sub(1) {
            for j in i + 2..n - 1 {
                let (a, b, c, d) = (points[i], points[i + 1], points[j], points[j + 1]);
                let (d1, d2) = (cross(a, b, c), cross(a, b, d));
                let (d3, d4) = (cross(c, d, a), cross(c, d, b));
                if d1 * d2 < 0.0 && d3 * d4 < 0.0 {
                    return Err(degenerate(
                        "station: an offset that crosses itself within its piece: the offset \
                         distance exceeds the room the base's curvature leaves"
                            .into(),
                    ));
                }
            }
        }
        Ok(())
    }

    /// The base's measure `v` into the piece at the offset's own measure
    /// `m`.
    fn base_measure(&self, m: Scalar, table: &OffsetTable) -> GeomResult<Scalar> {
        let length = self.base_length();
        if let Some(scale) = table.closed {
            return Ok((m / scale).clamp(0.0, length));
        }
        let last = table.sums.len() - 1;
        let m = m.clamp(0.0, table.sums[last]);
        let i = table.sums[1..last]
            .partition_point(|sum| *sum <= m)
            .min(last - 1);
        let (lo, hi) = (table.knots[i], table.knots[i + 1]);
        let target = m - table.sums[i];
        let span = table.sums[i + 1] - table.sums[i];
        let tolerance = 1e-3 * OFFSET_TOLERANCE * m.max(1.0);
        let travelled = |v: Scalar| -> GeomResult<Scalar> {
            let (half, mid) = (0.5 * (v - lo), 0.5 * (lo + v));
            let mut sum = 0.0;
            for (node, weight) in NODES.iter().zip(WEIGHTS.iter()) {
                sum += weight * self.speed(mid + half * node)?;
            }
            Ok(sum * half)
        };
        let (mut a, mut b) = (lo, hi);
        let mut v = if span > 0.0 {
            lo + (hi - lo) * (target / span)
        } else {
            lo
        };
        for _ in 0..MAX_ITERATIONS {
            let residual = travelled(v)? - target;
            if residual.abs() <= tolerance {
                return Ok(v);
            }
            if residual < 0.0 {
                a = v;
            } else {
                b = v;
            }
            let newton = v - residual / self.speed(v)?;
            v = if newton > a && newton < b {
                newton
            } else {
                0.5 * (a + b)
            };
            if !(v > a && v < b) {
                return Ok(v);
            }
        }
        Ok(v)
    }

    /// The offset's point at the base's measure `v`.
    fn point_at_base(&self, v: Scalar) -> GeomResult<Point3> {
        let section = self.base_section(v)?;
        Ok(section.point + self.displacement(&section, v)?)
    }

    /// The point at the offset's own measure `m`, through `table` (built
    /// afresh when `None`).
    pub(crate) fn point_at(&self, m: Scalar, table: Option<&OffsetTable>) -> GeomResult<Point3> {
        let built;
        let table = match table {
            Some(table) => table,
            None => {
                built = self.table()?;
                &built
            }
        };
        self.point_at_base(self.base_measure(m, table)?)
    }

    /// The section frame at the offset's own measure `m`, through `table`
    /// (built afresh when `None`); see the
    /// [module documentation](self#distance-convention).
    pub(crate) fn section_at(
        &self,
        m: Scalar,
        table: Option<&OffsetTable>,
    ) -> GeomResult<SectionFrame> {
        let built;
        let table = match table {
            Some(table) => table,
            None => {
                built = self.table()?;
                &built
            }
        };
        let v = self.base_measure(m, table)?;
        let (point, tangent, base) = if table.closed.is_some() {
            // A line or a circle beside its base: the base's tangent.
            let section = self.base_section(v)?;
            let point = section.point + self.displacement(&section, v)?;
            (point, section.tangent, section)
        } else {
            let (velocity, section) = self.velocity(v)?;
            let point = section.point + self.displacement(&section, v)?;
            (point, velocity, section)
        };
        own_frame(point, tangent, &base, self.banked())
    }
}

/// The offset's own section frame at `point` with tangent `direction`:
/// the base's rolled lateral made perpendicular to it beside a banked
/// curve; else the planar frame where it runs horizontally at `z` exactly
/// level, the reference-up frame elsewhere.
fn own_frame(
    point: Point3,
    direction: Vec3,
    base: &SectionFrame,
    banked: bool,
) -> GeomResult<SectionFrame> {
    let size = direction.length();
    if !(size.is_finite() && size > 0.0) {
        return Err(invalid(
            "station: the offset has no tangent direction there".into(),
        ));
    }
    let tangent = direction / size;
    if banked {
        let lateral = base.lateral - base.lateral.dot(tangent) * tangent;
        let magnitude = lateral.length();
        if !(magnitude.is_finite() && magnitude > 1e-12) {
            return Err(invalid(
                "station: the offset's tangent lies along its banked base's lateral".into(),
            ));
        }
        let lateral = lateral / magnitude;
        return Ok(SectionFrame {
            point,
            tangent,
            lateral,
            up: tangent.cross(lateral),
        });
    }
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

/// The offset of pieces laid end to end (a composite's, in its measure)
/// by `law`: one [`StationOffset`] piece per span of a base piece between
/// the seams of its curve, so every seam and joint of the base is a joint
/// of the offset. A [`OffsetLaw::Linear`] law runs from its `start` at the
/// first piece's start to its `end` at the last piece's end, linearly in
/// the base's measure.
///
/// The pieces are not joined here: [`super::CompositeBasis::new`] refuses
/// two that do not meet (an offset across a corner of its base, or a jump
/// of a banked base's roll) by name.
///
/// # Errors
///
/// No pieces, and the refusals of [`StationOffset::new`] and of the
/// offset's length, by name.
pub fn offset_pieces<'c>(
    base: &[StationPiece<'c>],
    law: OffsetLaw,
) -> GeomResult<Vec<StationPiece<'c>>> {
    if base.is_empty() {
        return Err(invalid("station: an offset of no pieces".into()));
    }
    let total: Scalar = base.iter().map(StationPiece::length).sum();
    let mut out = Vec::new();
    let mut run = 0.0;
    for piece in base {
        let mut cuts = vec![0.0];
        cuts.extend(piece.seams()?.into_iter().map(|seam| seam.distance));
        cuts.push(piece.length());
        for pair in cuts.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            let mut span = *piece;
            if piece.reversed {
                span.start = piece.end - b;
                span.end = piece.end - a;
            } else {
                span.start = piece.start + a;
                span.end = piece.start + b;
            }
            let law = match law {
                OffsetLaw::Linear { start, end, frame } => OffsetLaw::Linear {
                    start: start.lerp(end, (run + a) / total),
                    end: start.lerp(end, (run + b) / total),
                    frame,
                },
                other => other,
            };
            out.push(StationPiece::offset(StationOffset::new(span, law)?)?);
        }
        run += piece.length();
    }
    Ok(out)
}
