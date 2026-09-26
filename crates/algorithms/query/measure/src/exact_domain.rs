//! A face's domain in its surface parameters: the pcurves of its loops,
//! joined across seams and poles.
//!
//! Shared by the mass-property integral ([`crate::exact_face`]) and the
//! certified distance ([`crate::exact_distance`]). See `exact_face` for why
//! a periodic surface's loops need not close in the parameter plane, and how
//! a net winding reaches a pole.

use axiolid_brep::ExactBRep;
use axiolid_core::{Point2, Point3, Scalar, Vec2};
use axiolid_curve::Curve2;
use axiolid_evaluate::surface::evaluate;
use axiolid_evaluate::{derivative2, evaluate2};
use axiolid_surface::Surface;
use axiolid_topology::{Face, Orientation};
use core::f64::consts::{FRAC_PI_2, TAU};

use crate::exact::ExactMeasureError;

/// How a face's surface is parameterised: its periods and its poles.
pub(crate) struct Chart {
    pub(crate) u_period: Option<Scalar>,
    pub(crate) v_period: Option<Scalar>,
    /// `v` values where the surface collapses to one point for every `u`.
    pub(crate) poles: Vec<Scalar>,
}

pub(crate) fn chart(surface: &Surface) -> Result<Chart, ExactMeasureError> {
    let chart = match surface {
        Surface::Plane(_) | Surface::BSpline(_) => Chart {
            u_period: None,
            v_period: None,
            poles: Vec::new(),
        },
        Surface::Cylinder(_) | Surface::EllipticalCylinder(_) => Chart {
            u_period: Some(TAU),
            v_period: None,
            poles: Vec::new(),
        },
        Surface::Cone(cone) => {
            let slope = cone.semi_angle.tan();
            let poles = if slope.is_finite() && slope != 0.0 {
                vec![-cone.radius / slope]
            } else {
                Vec::new()
            };
            Chart {
                u_period: Some(TAU),
                v_period: None,
                poles,
            }
        }
        Surface::Sphere(_) => Chart {
            u_period: Some(TAU),
            v_period: None,
            poles: vec![-FRAC_PI_2, FRAC_PI_2],
        },
        Surface::Torus(_) => Chart {
            u_period: Some(TAU),
            v_period: Some(TAU),
            poles: Vec::new(),
        },
        _ => {
            return Err(ExactMeasureError::NonPlanarFace(crate::exact::family(
                surface,
            )))
        }
    };
    Ok(chart)
}

/// One stretch of a loop in the parameter plane.
pub(crate) enum Piece<'a> {
    /// A pcurve over `[start, end]`, shifted by whole periods.
    Curve {
        curve: &'a Curve2,
        start: Scalar,
        end: Scalar,
        offset: Vec2,
    },
    /// A straight closing stretch: across a pole, or a gap within tolerance.
    Segment { from: Point2, to: Point2 },
}

impl Piece<'_> {
    pub(crate) fn span(&self) -> (Scalar, Scalar) {
        match self {
            Self::Curve { start, end, .. } => (*start, *end),
            Self::Segment { .. } => (0.0, 1.0),
        }
    }

    pub(crate) fn at(&self, t: Scalar) -> Result<(Point2, Vec2), ExactMeasureError> {
        match self {
            Self::Curve { curve, offset, .. } => {
                let point = evaluate2(curve, t).map_err(|_| ExactMeasureError::Evaluation)?;
                let tangent = derivative2(curve, t).map_err(|_| ExactMeasureError::Evaluation)?;
                Ok((point + *offset, tangent))
            }
            Self::Segment { from, to } => Ok((*from + (*to - *from) * t, *to - *from)),
        }
    }

    pub(crate) fn end_point(&self) -> Result<Point2, ExactMeasureError> {
        Ok(self.at(self.span().1)?.0)
    }
}

/// A face's boundary, assembled in the parameter plane.
pub(crate) struct Boundary<'a> {
    pub(crate) pieces: Vec<Piece<'a>>,
    /// Net whole periods each loop winds, summed over the face.
    pub(crate) winding: [i64; 2],
    /// Whether any single loop winds in `u` / in `v`.
    pub(crate) wraps: [bool; 2],
    /// A point on the boundary, for choosing the reference line.
    pub(crate) anchor: Point2,
}

/// The pole the face's domain reaches, used as the reference line so the
/// pole itself -- the boundary no loop states -- contributes zero.
pub(crate) fn pole_on_domain_side(
    chart: &Chart,
    boundary: &Boundary<'_>,
) -> Result<Scalar, ExactMeasureError> {
    // A loop running +u keeps its domain on its left, towards +v.
    let upward = match boundary.winding[0] {
        1 => true,
        -1 => false,
        _ => {
            return Err(ExactMeasureError::ParameterDomain(
                "face boundary winds around the surface more than once",
            ))
        }
    };
    chart
        .poles
        .iter()
        .copied()
        .filter(|pole| (*pole > boundary.anchor.y) == upward)
        .min_by(|a, b| {
            (a - boundary.anchor.y)
                .abs()
                .total_cmp(&(b - boundary.anchor.y).abs())
        })
        .ok_or(ExactMeasureError::ParameterDomain(
            "face boundary winds around a surface with no pole on the domain side",
        ))
}

/// Assemble every bound's pcurves into one boundary, joining each use to the
/// next through whole periods, across a pole, or over a gap within tolerance.
pub(crate) fn assemble<'a>(
    brep: &'a ExactBRep,
    face: &Face<axiolid_brep::SurfaceId>,
    surface: &Surface,
    chart: &Chart,
    linear: Scalar,
) -> Result<Boundary<'a>, ExactMeasureError> {
    let topology = brep.topology();
    let mut pieces = Vec::new();
    let mut winding = [0_i64; 2];
    let mut wraps = [false; 2];
    let mut anchor = None;

    for bound in &face.bounds {
        let wire = topology
            .loops()
            .get(bound.loop_id.index())
            .ok_or(ExactMeasureError::DanglingReference)?;
        let reversed = bound.orientation == Orientation::Reversed;

        // Each use's pcurve interval already runs in the use's traversal
        // (ADR 0024); a reversed bound walks the loop backwards.
        let mut uses = Vec::with_capacity(wire.edges.len());
        for (index, use_) in wire.edges.iter().enumerate() {
            let pcurve = use_.pcurve.ok_or(ExactMeasureError::DanglingReference)?;
            let curve = brep
                .curves2()
                .get(pcurve.index())
                .ok_or(ExactMeasureError::DanglingReference)?;
            let interval = brep
                .pcurve_interval(bound.loop_id, index)
                .ok_or(ExactMeasureError::DanglingReference)?;
            let (start, end) = if reversed {
                (interval.end, interval.start)
            } else {
                (interval.start, interval.end)
            };
            uses.push((curve, start, end));
        }
        if reversed {
            uses.reverse();
        }
        let Some(&(first_curve, first_start, _)) = uses.first() else {
            continue;
        };

        let loop_start =
            evaluate2(first_curve, first_start).map_err(|_| ExactMeasureError::Evaluation)?;
        anchor.get_or_insert(loop_start);
        let mut offset = Vec2::ZERO;
        let mut cursor = loop_start;
        for (curve, start, end) in uses {
            let raw = evaluate2(curve, start).map_err(|_| ExactMeasureError::Evaluation)?;
            let (shift, bridge) = join(surface, chart, cursor, raw + offset, linear)?;
            offset += shift;
            if let Some(segment) = bridge {
                pieces.push(segment);
            }
            let piece = Piece::Curve {
                curve,
                start,
                end,
                offset,
            };
            cursor = piece.end_point()?;
            pieces.push(piece);
        }

        // Close the loop back onto its own start.
        let (shift, bridge) = join(surface, chart, cursor, loop_start, linear)?;
        if let Some(segment) = bridge {
            pieces.push(segment);
        }
        // `shift` moved the start onto the end: the loop ended that many
        // periods on from where it began.
        let turns = [
            periods(shift.x, chart.u_period),
            periods(shift.y, chart.v_period),
        ];
        for axis in 0..2 {
            winding[axis] += turns[axis];
            wraps[axis] |= turns[axis] != 0;
        }
    }

    let anchor = anchor.ok_or(ExactMeasureError::Degenerate)?;
    Ok(Boundary {
        pieces,
        winding,
        wraps,
        anchor,
    })
}

/// Whole periods in `shift`, which is already a multiple of `period`.
fn periods(shift: Scalar, period: Option<Scalar>) -> i64 {
    period.map_or(0, |period| (shift / period).round() as i64)
}

/// Join a loop that has reached `from` to the next use starting at `to`.
///
/// Returns the whole-period shift to apply to `to` and everything after it,
/// and the closing segment, if one is needed. The two must be the same point
/// on the surface; the parameter gap between them may be whole periods (a
/// seam), a stretch along a pole, or a residue within tolerance.
fn join<'a>(
    surface: &Surface,
    chart: &Chart,
    from: Point2,
    to: Point2,
    linear: Scalar,
) -> Result<(Vec2, Option<Piece<'a>>), ExactMeasureError> {
    let there = point(surface, from)?;
    let here = point(surface, to)?;
    if (there - here).length() > linear {
        return Err(ExactMeasureError::ParameterDomain(
            "consecutive pcurves of a face loop do not meet on the surface",
        ));
    }

    // Along a pole `u` is free: keep the stretch as it is, so its `H du`
    // is integrated, rather than folding it into periods.
    let on_pole = chart.poles.iter().any(|pole| {
        (from.y - pole).abs() <= parameter_slack(*pole)
            && (to.y - pole).abs() <= parameter_slack(*pole)
    });
    let shift = if on_pole {
        Vec2::ZERO
    } else {
        let wrap = |gap: Scalar, period: Option<Scalar>| {
            period.map_or(0.0, |period| (gap / period).round() * period)
        };
        Vec2::new(
            wrap(from.x - to.x, chart.u_period),
            wrap(from.y - to.y, chart.v_period),
        )
    };
    let landed = to + shift;
    let bridge = (landed != from).then_some(Piece::Segment { from, to: landed });
    Ok((shift, bridge))
}

fn parameter_slack(value: Scalar) -> Scalar {
    1e-9 * value.abs().max(1.0)
}

pub(crate) fn point(surface: &Surface, at: Point2) -> Result<Point3, ExactMeasureError> {
    evaluate(surface, at.x, at.y).map_err(|_| ExactMeasureError::Evaluation)
}

/// One stretch of a boundary piece over which both parameters are monotone.
///
/// Its box is the box of its two ends, exactly, and a line `u = c` crosses
/// it at most once -- which is what lets a crossing be decided from the
/// ends alone.
struct MonoArc {
    piece: usize,
    t0: Scalar,
    t1: Scalar,
    a: Point2,
    b: Point2,
}

/// A face's parameter domain, able to say which points lie in it.
///
/// Decisions are certified or refused: a point too close to the boundary
/// for the crossing count to be sure is `None`, never a guess.
pub(crate) struct Domain<'a> {
    boundary: Boundary<'a>,
    arcs: Vec<MonoArc>,
    /// Loops wind round the tube (`v`): coordinates are swapped so the ray
    /// is always cast in the second one.
    transpose: bool,
    /// Period of the first (possibly swapped) coordinate.
    period: Option<Scalar>,
    /// The pole the domain reaches and the net winding that reaches it.
    pole: Option<(Scalar, i64)>,
    /// Box of the domain, in unswapped `(u, v)`.
    pub(crate) min: Point2,
    pub(crate) max: Point2,
}

/// Parameters in `(lo, hi)` where one coordinate of a pcurve turns.
fn turning_points(curve: &Curve2, lo: Scalar, hi: Scalar) -> Option<Vec<Scalar>> {
    let mut out = Vec::new();
    let mut add_trig = |a: Scalar, b: Scalar| {
        // `a cos t + b sin t` turns where its derivative vanishes.
        if a == 0.0 && b == 0.0 {
            return;
        }
        let base = b.atan2(a);
        let first = ((lo - base) / core::f64::consts::PI).floor() as i64;
        let last = ((hi - base) / core::f64::consts::PI).ceil() as i64;
        for k in first..=last {
            let t = base + k as Scalar * core::f64::consts::PI;
            if t > lo && t < hi {
                out.push(t);
            }
        }
    };
    match curve {
        Curve2::Line(_) => {}
        Curve2::Circle(c) => {
            add_trig(c.radius * c.frame.x.x, c.radius * c.frame.y.x);
            add_trig(c.radius * c.frame.x.y, c.radius * c.frame.y.y);
        }
        Curve2::Ellipse(e) => {
            add_trig(e.semi_axis_x * e.frame.x.x, e.semi_axis_y * e.frame.y.x);
            add_trig(e.semi_axis_x * e.frame.x.y, e.semi_axis_y * e.frame.y.y);
        }
        Curve2::Sinusoid(w) => add_trig(w.cosine, w.sine),
        Curve2::Polyline(_) => {
            let mut k = lo.floor() + 1.0;
            while k < hi {
                out.push(k);
                k += 1.0;
            }
        }
        _ => return None,
    }
    out.sort_by(Scalar::total_cmp);
    Some(out)
}

fn slack(value: Scalar) -> Scalar {
    1e-9 * (1.0 + value.abs())
}

impl<'a> Domain<'a> {
    /// The face's domain, or `None` when a pcurve family has no monotone
    /// split here (a B-spline or intrinsic trim): such a face is still
    /// bounded, just never classified.
    pub(crate) fn new(
        brep: &'a ExactBRep,
        face: &Face<axiolid_brep::SurfaceId>,
        surface: &Surface,
        linear: Scalar,
    ) -> Result<Option<Self>, ExactMeasureError> {
        let chart = chart(surface)?;
        let boundary = assemble(brep, face, surface, &chart, linear)?;
        let transpose = boundary.wraps[1];
        if transpose && (boundary.wraps[0] || boundary.winding[1] != 0) {
            return Err(ExactMeasureError::ParameterDomain(
                "face boundary winds around the surface in both directions",
            ));
        }
        let pole = if !transpose && boundary.winding[0] != 0 {
            Some((pole_on_domain_side(&chart, &boundary)?, boundary.winding[0]))
        } else {
            None
        };
        let swap = |p: Point2| if transpose { Point2::new(p.y, p.x) } else { p };

        let mut arcs = Vec::new();
        for (index, piece) in boundary.pieces.iter().enumerate() {
            let (start, end) = piece.span();
            let (lo, hi) = (start.min(end), start.max(end));
            let mut cuts = match piece {
                Piece::Curve { curve, .. } => match turning_points(curve, lo, hi) {
                    Some(cuts) => cuts,
                    None => return Ok(None),
                },
                Piece::Segment { .. } => Vec::new(),
            };
            if start > end {
                cuts.reverse();
            }
            let mut ts = vec![start];
            ts.extend(cuts);
            ts.push(end);
            for pair in ts.windows(2) {
                let a = swap(piece.at(pair[0])?.0);
                let b = swap(piece.at(pair[1])?.0);
                arcs.push(MonoArc {
                    piece: index,
                    t0: pair[0],
                    t1: pair[1],
                    a,
                    b,
                });
            }
        }
        if arcs.is_empty() {
            return Ok(None);
        }

        let mut min = Point2::splat(Scalar::INFINITY);
        let mut max = Point2::splat(Scalar::NEG_INFINITY);
        for arc in &arcs {
            min = min.min(arc.a.min(arc.b));
            max = max.max(arc.a.max(arc.b));
        }
        if let Some((pole, _)) = pole {
            min.y = min.y.min(pole);
            max.y = max.y.max(pole);
        }
        let period = if transpose {
            chart.v_period
        } else {
            chart.u_period
        };
        let (min, max) = if transpose {
            (Point2::new(min.y, min.x), Point2::new(max.y, max.x))
        } else {
            (min, max)
        };
        Ok(Some(Self {
            boundary,
            arcs,
            transpose,
            period,
            pole,
            min,
            max,
        }))
    }

    fn swap(&self, p: Point2) -> Point2 {
        if self.transpose {
            Point2::new(p.y, p.x)
        } else {
            p
        }
    }

    /// Whole-period shifts `s` of the first coordinate for which the span
    /// `[lo, hi]`, moved to `[lo - s, hi - s]`, can meet the boundary's own
    /// range `[min, max]`: `s` in `[lo - max, hi - min]`.
    fn shifts(&self, lo: Scalar, hi: Scalar) -> Vec<Scalar> {
        let Some(period) = self.period else {
            return vec![0.0];
        };
        let (min, max) = self.arcs.iter().fold(
            (Scalar::INFINITY, Scalar::NEG_INFINITY),
            |(min, max), arc| (min.min(arc.a.x.min(arc.b.x)), max.max(arc.a.x.max(arc.b.x))),
        );
        let first = ((lo - max) / period).floor() as i64;
        let last = ((hi - min) / period).ceil() as i64;
        (first..=last).map(|k| k as Scalar * period).collect()
    }

    /// Whether the boundary may pass through the rectangle.
    ///
    /// `false` is certain. A monotone arc's box is exact at its ends but
    /// covers a whole triangle beside a diagonal, so an arc whose box meets
    /// the rectangle is bisected until its halves clear the rectangle or
    /// one lies inside it.
    pub(crate) fn touches(&self, lo: Point2, hi: Point2) -> Result<bool, ExactMeasureError> {
        let (lo, hi) = {
            let (a, b) = (self.swap(lo), self.swap(hi));
            (a.min(b), a.max(b))
        };
        for shift in self.shifts(lo.x, hi.x) {
            let (lo, hi) = (
                Point2::new(lo.x - shift, lo.y),
                Point2::new(hi.x - shift, hi.y),
            );
            for arc in &self.arcs {
                if self.arc_touches(arc, lo, hi, 0)? {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }

    fn arc_touches(
        &self,
        arc: &MonoArc,
        lo: Point2,
        hi: Point2,
        depth: u32,
    ) -> Result<bool, ExactMeasureError> {
        let (a_lo, a_hi) = (arc.a.min(arc.b), arc.a.max(arc.b));
        let (sx, sy) = (
            slack(a_hi.x.abs().max(a_lo.x.abs())),
            slack(a_hi.y.abs().max(a_lo.y.abs())),
        );
        if a_lo.x - sx > hi.x || a_hi.x + sx < lo.x || a_lo.y - sy > hi.y || a_hi.y + sy < lo.y {
            return Ok(false);
        }
        // A straight stretch lying along one of the rectangle's own sides
        // (a rim or seam of a wall, cut exactly where the patch was cut)
        // does not enter it: the rectangle is still wholly in or out.
        let along = |a: Scalar, b: Scalar, edge: Scalar, s: Scalar| {
            (a - b).abs() <= s && (a - edge).abs() <= s && (b - edge).abs() <= s
        };
        if along(arc.a.x, arc.b.x, lo.x, sx)
            || along(arc.a.x, arc.b.x, hi.x, sx)
            || along(arc.a.y, arc.b.y, lo.y, sy)
            || along(arc.a.y, arc.b.y, hi.y, sy)
        {
            return Ok(false);
        }
        let inside = |p: Point2| p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y;
        if inside(arc.a) || inside(arc.b) || depth >= 40 {
            return Ok(true);
        }
        let tm = 0.5 * (arc.t0 + arc.t1);
        let m = self.swap(self.boundary.pieces[arc.piece].at(tm)?.0);
        let first = MonoArc {
            piece: arc.piece,
            t0: arc.t0,
            t1: tm,
            a: arc.a,
            b: m,
        };
        let second = MonoArc {
            piece: arc.piece,
            t0: tm,
            t1: arc.t1,
            a: m,
            b: arc.b,
        };
        Ok(self.arc_touches(&first, lo, hi, depth + 1)?
            || self.arc_touches(&second, lo, hi, depth + 1)?)
    }

    /// Whether `p` lies in the domain: `Some` when certain, `None` when it
    /// is too close to the boundary to say.
    pub(crate) fn contains(&self, p: Point2) -> Result<Option<bool>, ExactMeasureError> {
        // Outside the domain's box in a coordinate that does not wrap is
        // outside, certainly -- and it is where a ray cast along the box's
        // own edge could not decide.
        let beyond = |value: Scalar, lo: Scalar, hi: Scalar| {
            value < lo - slack(lo) || value > hi + slack(hi)
        };
        let (u_wraps, v_wraps) = if self.transpose {
            (false, self.period.is_some())
        } else {
            (self.period.is_some(), false)
        };
        if (!u_wraps && beyond(p.x, self.min.x, self.max.x))
            || (!v_wraps && beyond(p.y, self.min.y, self.max.y))
        {
            return Ok(Some(false));
        }
        let q = self.swap(p);
        let mut total: i64 = 0;
        for shift in self.shifts(q.x, q.x) {
            for arc in &self.arcs {
                match self.crossing(arc, q.x - shift, q.y, 0)? {
                    Some(weight) => total += weight,
                    None => return Ok(None),
                }
            }
        }
        if let Some((pole, winding)) = self.pole {
            if (pole - q.y).abs() <= slack(pole) {
                return Ok(None);
            }
            if pole > q.y {
                total += winding;
            }
        }
        // Swapping the coordinates mirrors the plane: an anticlockwise loop
        // counts clockwise.
        if self.transpose {
            total = -total;
        }
        Ok(match total {
            0 => Some(false),
            1 => Some(true),
            _ => None,
        })
    }

    /// Signed crossing of the ray `x = c, y > y0` with one monotone arc:
    /// `+1` where the boundary runs towards `-x`, `-1` towards `+x`.
    fn crossing(
        &self,
        arc: &MonoArc,
        c: Scalar,
        y0: Scalar,
        depth: u32,
    ) -> Result<Option<i64>, ExactMeasureError> {
        let (a, b) = (arc.a, arc.b);
        let dx = slack(c);
        if c < a.x.min(b.x) - dx || c > a.x.max(b.x) + dx {
            return Ok(Some(0));
        }
        if (c - a.x).abs() <= dx || (c - b.x).abs() <= dx {
            // Through an end, or along a piece parallel to the ray.
            return Ok(None);
        }
        let weight = if b.x < a.x { 1 } else { -1 };
        let dy = slack(y0);
        if a.y.min(b.y) > y0 + dy {
            return Ok(Some(weight));
        }
        if a.y.max(b.y) < y0 - dy {
            return Ok(Some(0));
        }
        if depth >= 48 {
            return Ok(None);
        }
        // Monotone: the half that spans `c` holds the crossing.
        let tm = 0.5 * (arc.t0 + arc.t1);
        let m = self.swap(self.boundary.pieces[arc.piece].at(tm)?.0);
        let half = if (c - a.x) * (c - m.x) < 0.0 {
            MonoArc {
                piece: arc.piece,
                t0: arc.t0,
                t1: tm,
                a,
                b: m,
            }
        } else {
            MonoArc {
                piece: arc.piece,
                t0: tm,
                t1: arc.t1,
                a: m,
                b,
            }
        };
        self.crossing(&half, c, y0, depth + 1)
    }
}

/// Which parameters of one exact face lie inside it, decided with a
/// certificate or not at all.
///
/// The public face of the crate-private `Domain`: the pcurves of the face's loops, joined
/// across seams and poles, split into pieces monotone in both parameters
/// so that a ray crossing is decided from the pieces' ends. `contains`
/// answers `Some(true)` or `Some(false)` only when certain, and `None` for
/// a point too close to the boundary to say.
pub struct FaceDomain<'a> {
    domain: Domain<'a>,
}

impl<'a> FaceDomain<'a> {
    /// The domain of `face`, or `None` when a pcurve family has no monotone
    /// split here (a B-spline or intrinsic trim).
    ///
    /// # Errors
    ///
    /// A missing surface, a dangling handle, a boundary that encloses no
    /// domain in the surface's parameters, or an evaluation failure.
    pub fn new(
        brep: &'a ExactBRep,
        face: axiolid_topology::FaceId,
        tolerance: axiolid_core::Tolerance,
    ) -> Result<Option<Self>, ExactMeasureError> {
        let topology = brep.topology();
        let record = topology
            .faces()
            .get(face.index())
            .ok_or(ExactMeasureError::DanglingReference)?;
        let surface = record
            .surface
            .and_then(|id| brep.surfaces().get(id.index()))
            .ok_or(ExactMeasureError::MissingSurface)?;
        let linear = tolerance.linear().max(1e-12);
        Ok(Domain::new(brep, record, surface, linear)?.map(|domain| Self { domain }))
    }

    /// Whether the parameters `at` lie in the face: `None` when too close to
    /// its boundary to decide.
    ///
    /// # Errors
    ///
    /// A pcurve that cannot be evaluated where the decision needs it.
    pub fn contains(&self, at: Point2) -> Result<Option<bool>, ExactMeasureError> {
        self.domain.contains(at)
    }

    /// A box in `(u, v)` holding the whole domain.
    #[must_use]
    pub fn bounds(&self) -> (Point2, Point2) {
        (self.domain.min, self.domain.max)
    }
}
