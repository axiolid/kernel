//! Sections through a point where the two surfaces touch, at a smooth edge
//! (#249).
//!
//! A round web hole touching an I-beam's flange also touches each root
//! fillet: the fillet is a quarter cylinder running tangent into the
//! flange, so at the fillet/flange edge the hole and the fillet cylinder
//! (axes crossing at right angles, skew by `r_hole - r_fillet`) share
//! their tangent plane. Their section is a quartic with a double point
//! there: two loops round the fillet cylinder, `z = z0 +- g(phi)` in the
//! fillet's angle `phi`, with `g` vanishing like `sqrt(r_hole r_fillet)
//! |phi - phi0|` at the edge, so two branches cross there transversally.
//! Three things are needed to build it, each placed here:
//!
//! - **The trace** ([`trace`]). The section is traced as an implicit curve
//!   (ADR 0077) on one cylinder, in its face's parameter box. The fillet
//!   face ends at the double point, so that box's edge runs through it,
//!   where the trace cannot certify the pieces about the crossing and
//!   refuses (`Undecided`). Then the trace is taken again over the whole
//!   turn of the carrier: the double point is inside, the trace finds the
//!   saddle of the field and ends both loops at it (a vertex of the
//!   section, ADR 0077 "Singular points"). Pieces off the face are dropped
//!   as for any section. On the carrier's own face each piece keeps the
//!   traced curve itself as its pcurve (`crate::split`), so no second trace
//!   over the face's box is needed.
//! - **Cuts at a smooth edge** ([`smooth_edge_cuts`]). A section on the
//!   fillet is cut where it crosses the fillet's edges. Along a smooth edge
//!   (the fillet runs tangent into the flange or web there) the section
//!   touches the adjacent surface instead of crossing it: the exact
//!   crossing is a double root, which rounding splits or loses. The edge
//!   itself lies on the fillet, so the section crosses it exactly where the
//!   edge crosses the section's other surface (the hole): a line against a
//!   cylinder, transversal at the web edge and a double root exactly at the
//!   double point (closed form, decided exactly). Whether an edge is smooth
//!   only chooses which of these two exact formulations is evaluated; it
//!   decides no geometry.
//! - **One contact** (#243). Under a general placement, or a fraction of
//!   the tolerance off touching, the flange plane is read as touching the
//!   hole (`PlaneTouchesCylinder`), and every crossing of a curve in the
//!   flange plane with the hole is placed on the contact ruling. The
//!   fillet/flange edge is such a curve, so the double point is placed
//!   there too, and the section must pass through it within tolerance.
//!   Where it does (the residue is rounding: the trace sees one crossing),
//!   the cut is that point. Where it does not, the hole and fillet meet in
//!   two separate arcs that the moved flange would have to join, which no
//!   single move of one operand does (moving the hole onto the flange
//!   makes it touch the fillet too; moving the flange alone leaves the
//!   fillet crossing it): that is refused by name
//!   ([`BooleanError::UnsupportedContact`]). A traced section on a cylinder
//!   that some plane is read as touching meets the contact ruling where it
//!   crosses the plane through the ruling and the axis
//!   ([`traced_on_ruling`]).
//!
//! An exactly tangent pair (exact placements, dyadic sizes) needs no
//! reading: the double point is the exact double root of the edge against
//! the hole, and the trace's crossing lies on it to rounding. (A fillet
//! lowered from a profile carries its radius from the arc's bulge, an ulp
//! or so off; the trace's crossing absorbs that residue, as it does any
//! below its field's rounding, ADR 0077.)

use axiolid_core::{Frame3, Interval, Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::{Curve3, ImplicitSection3, Line3};
use axiolid_evaluate::surface::{locate, normal};
use axiolid_evaluate::{curve::locate3, evaluate3};
use axiolid_nurbs::extrema::{minimum_distance, Piece as DistancePiece};
use axiolid_nurbs::{
    exact_curve_surface_intersection, implicit_surface_intersection, ExactCurveIntersection,
    ExactIntersectionRefusal,
};
use axiolid_surface::{Plane, Surface};
use core::f64::consts::TAU;

use crate::contact::Contacts;
use crate::report::{self, ToleranceDecisionKind};
use crate::support::{periods, window};
use crate::BooleanError;

/// How far apart two unit normals may be for an edge to count as smooth:
/// a fixed slack far below any tolerance. It only chooses which exact
/// formulation of a crossing is evaluated.
const SMOOTH: Scalar = 1e-9;

/// The traced section of `carrier` by `other` over the face's parameter box
/// `[lo, hi]`; where that cannot be decided (a crossing of the section on
/// the box's edge), over a whole turn of a periodic carrier.
pub(crate) fn trace(
    carrier: &Surface,
    other: &Surface,
    lo: Point2,
    hi: Point2,
) -> Result<Vec<ImplicitSection3>, ExactIntersectionRefusal> {
    match implicit_surface_intersection(carrier, other, Some(window(carrier, lo, hi))) {
        Err(ExactIntersectionRefusal::Undecided) => match whole_turn(carrier, lo, hi) {
            Some(turn) => implicit_surface_intersection(carrier, other, Some(turn)),
            None => Err(ExactIntersectionRefusal::Undecided),
        },
        other => other,
    }
}

/// The trace window of a whole turn of `surface`'s angle about the face's
/// parameter box `[lo, hi]`, when the angle is periodic and the box is
/// less than a turn wide.
fn whole_turn(surface: &Surface, lo: Point2, hi: Point2) -> Option<(Point2, Point2)> {
    if !periods(surface).0 || hi.x - lo.x >= TAU - 1e-9 {
        return None;
    }
    let mid = 0.5 * (lo.x + hi.x);
    Some(window(
        surface,
        Point2::new(mid - 0.5 * TAU, lo.y),
        Point2::new(mid + 0.5 * TAU, hi.y),
    ))
}

/// Where `curve`, a section lying on `own` and `meets`, crosses an edge
/// (`edge` over `span`) of `own`'s face whose adjacent surface `cutter`
/// runs tangent into `own` along it: where the edge crosses `meets`, as
/// parameters on `curve` (see the module docs). `None` when the edge is
/// not smooth, or its crossing with `meets` has no closed form here; the
/// caller then cuts against `cutter`.
///
/// # Errors
///
/// [`BooleanError::UnsupportedContact`] when the edge lies in a plane read
/// as touching `meets` and the section misses the contact point by more
/// than the tolerance while passing within the contact's reach of it.
#[allow(clippy::too_many_arguments)]
pub(crate) fn smooth_edge_cuts(
    own: &Surface,
    cutter: &Surface,
    meets: &Surface,
    edge: &Curve3,
    span: Interval,
    curve: &Curve3,
    contacts: &Contacts,
    tolerance: Tolerance,
) -> Result<Option<Vec<Scalar>>, BooleanError> {
    if !smooth(own, cutter, edge, span, tolerance)? {
        return Ok(None);
    }
    let (points, placed) = match contacts.crossing([own, cutter], meets, edge, tolerance)? {
        Some(params) => (
            params
                .into_iter()
                .map(|t| evaluate3(edge, t).map_err(|_| BooleanError::Evaluation))
                .collect::<Result<Vec<_>, _>>()?,
            true,
        ),
        None => match exact_curve_surface_intersection(edge, meets) {
            Ok(ExactCurveIntersection::Points(hits)) => {
                (hits.into_iter().map(|h| h.point).collect(), false)
            }
            _ => return Ok(None),
        },
    };
    let mut out = Vec::new();
    for point in points {
        // A point the edge's curve cannot locate is off it: the far root
        // of an edge running nearly parallel to `meets`.
        if !crate::section::on_span(edge, span, point, tolerance).unwrap_or(false) {
            continue;
        }
        if let Ok(s) = locate3(curve, point, report::floored(tolerance)) {
            let on = evaluate3(curve, s).map_err(|_| BooleanError::Evaluation)?;
            if report::near(
                ToleranceDecisionKind::IncidentPoint,
                (on - point).length(),
                tolerance,
            ) {
                out.push(s);
                continue;
            }
        }
        // The contact puts the crossing where the section does not pass,
        // though the section passes within the contact's reach of it: the
        // two readings disagree, and no single move of one operand makes
        // them agree. A distance that cannot be bounded is refused too.
        if placed && distance(curve, point).is_none_or(|d| d <= reach(meets, tolerance)) {
            return Err(BooleanError::UnsupportedContact);
        }
    }
    Ok(Some(out))
}

/// A lower bound on how far `point` lies from `curve` (a line, a conic or
/// a traced section), or `None` where it cannot be bounded here.
fn distance(curve: &Curve3, point: Point3) -> Option<Scalar> {
    let span = match curve {
        Curve3::Line(line) => {
            let d = line.direction.normalize_or_zero();
            let off = point - line.origin;
            return (d != Vec3::ZERO).then(|| (off - d * off.dot(d)).length());
        }
        Curve3::Circle(_) | Curve3::Ellipse(_) => Interval::new(0.0, TAU),
        Curve3::ImplicitSection(s) => Interval::new(0.0, s.curve.end()),
        _ => return None,
    };
    minimum_distance(
        &DistancePiece::Curve { curve, span },
        &DistancePiece::Point(point),
        1e-9,
    )
    .ok()
    .map(|e| e.lower)
}

/// How far from a contact point a section crossing the plane read as
/// touching a cylinder can lie: the exact crossings of a plane up to the
/// tolerance into the cylinder are within `2 sqrt(2 r eps)` of the contact
/// ruling, and a curve through them within twice that of its point.
fn reach(cylinder: &Surface, tolerance: Tolerance) -> Scalar {
    let r = match cylinder {
        Surface::Cylinder(c) => c.radius,
        _ => 0.0,
    };
    let eps = tolerance.linear();
    4.0 * (2.0 * r * eps).sqrt() + eps
}

/// Whether `own` and `cutter` share their tangent plane along the edge, to
/// [`SMOOTH`]: sampled at the edge's middle.
fn smooth(
    own: &Surface,
    cutter: &Surface,
    edge: &Curve3,
    span: Interval,
    tolerance: Tolerance,
) -> Result<bool, BooleanError> {
    let mid =
        evaluate3(edge, 0.5 * (span.start + span.end)).map_err(|_| BooleanError::Evaluation)?;
    let at = |s: &Surface| -> Result<Vec3, BooleanError> {
        let (u, v) =
            locate(s, mid, report::floored(tolerance)).map_err(|_| BooleanError::Evaluation)?;
        Ok(normal(s, u, v)
            .map_err(|_| BooleanError::Evaluation)?
            .normalize_or_zero())
    };
    let (a, b) = (at(own)?, at(cutter)?);
    Ok(a != Vec3::ZERO && b != Vec3::ZERO && a.cross(b).length() <= SMOOTH)
}

/// Where a traced section on `cylinder` meets the contact `ruling` of a
/// plane read as touching it (`crate::contact`): where it crosses the plane
/// through the ruling and the axis, kept where that crossing lies on the
/// ruling within tolerance (that plane holds the opposite ruling too).
///
/// # Errors
///
/// [`BooleanError::UnsupportedContact`] when that crossing cannot be found
/// (the curve runs in the plane, or the exact intersection refuses).
pub(crate) fn traced_on_ruling(
    cylinder: &Surface,
    ruling: &Line3,
    curve: &Curve3,
    tolerance: Tolerance,
) -> Result<Vec<Scalar>, BooleanError> {
    let Surface::Cylinder(c) = cylinder else {
        return Err(BooleanError::UnsupportedContact);
    };
    let axis = c.frame.z.normalize_or_zero();
    let off = ruling.origin - c.frame.origin;
    let radial = (off - axis * off.dot(axis)).normalize_or_zero();
    let across = axis.cross(radial);
    if across == Vec3::ZERO {
        return Err(BooleanError::UnsupportedContact);
    }
    let plane = Surface::Plane(Plane {
        frame: Frame3 {
            origin: ruling.origin,
            x: axis,
            y: radial,
            z: across,
        },
    });
    let Ok(ExactCurveIntersection::Points(hits)) = exact_curve_surface_intersection(curve, &plane)
    else {
        return Err(BooleanError::UnsupportedContact);
    };
    let d = ruling.direction.normalize_or_zero();
    let mut out = Vec::new();
    for hit in hits {
        let rel = hit.point - ruling.origin;
        let gap = (rel - d * rel.dot(d)).length();
        if report::near(ToleranceDecisionKind::TangentCrossing, gap, tolerance) {
            out.push(hit.parameter.approx());
        }
    }
    Ok(out)
}
