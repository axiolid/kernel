//! The boundary of a curved face whose edges carry no p-curve (#310).
//!
//! Exchange formats usually give an edge only its 3D curve: IFC's
//! `IfcEdgeCurve` without an `IfcPcurve`, for one. The face still needs its
//! boundary in its surface's parameters, so the edge's 3D samples are
//! inverted onto the surface instead of a p-curve being fitted. The samples
//! are the ones every face sharing the edge uses (`state.edge_samples`): the
//! edge pass ([`super::edges`], #336) samples its 3D curve as finely as the
//! strictest face using it asks, by [`sample_edge`] run for each face, and
//! every face inverts the same points, so shared vertices stay shared.
//! The edge is sampled in pieces between the derivative breaks of its 3D
//! curve (#342, [`EdgeTrace::breaks`]): each break is a sample, so no span
//! crosses a corner, and inside a piece the spans are uniform, a power of
//! two of them, refined per piece. The pieces' counts together stay within
//! the segment cap.
//!
//! - Analytic surfaces invert in closed form
//!   ([`axiolid_reference::surface::invert`]), which checks the residual
//!   against the linear tolerance. B-spline surfaces invert by the certified
//!   global search ([`axiolid_nurbs::invert_surface_certified`], the
//!   periodic form for a surface whose closed axes validate as periodic),
//!   which refuses an off-surface point as soon as it proves the separation
//!   (#320) and parameters that are not unique. A face over a surface
//!   relation inverts by the certified search along its swept parameter
//!   ([`RelationInverter`], #454, ADR 0110), with the same refusals; every
//!   step below applies to it as to an atomic face.
//! - An edge may end at a pole: a cone's apex, a sphere's pole, a
//!   revolution's profile on its axis or a ruled apex, where every `u`
//!   names the point (#463, ADR 0112). That end is not inverted:
//!   [`SurfaceInverter::pole`] reads its `v`, and it takes the `u` of its
//!   neighbouring sample on the edge, so the edge's last span runs along a
//!   meridian (a ruling) in parameters. An edge with one span between two
//!   poles is sampled again until it has two. Any other sample at, or
//!   within the linear tolerance of, a pole is refused by name (the edge
//!   runs through the pole; the producer splits it at the pole's vertex),
//!   and so is one on a collapsed ruling. A span into a pole is halved
//!   only while its middle stays off the pole.
//! - Periodic parameters are unwrapped along the ring by continuity. That
//!   reading is unambiguous only while consecutive samples turn less than a
//!   quarter period, so a new edge is sampled until they do; an edge whose
//!   samples still jump further at the segment cap, or a shared edge whose
//!   samples jump further, runs into a singularity and is refused.
//! - The reported deviation adds, per span, the edge's certified 3D chord
//!   bound `d` ([`axiolid_reference::bound::chord_bound3`]) and a bound `e`
//!   on how far the surface image of the span's parameter chord `D` strays
//!   from the 3D chord: `q(D) / 8` from the surface's interpolation
//!   coefficients over the chord's box (the one-dimensional case of the
//!   triangle bound in the parent module's notes), or `|S_u| |D_u| + |S_v|
//!   |D_v|` where a full-multiplicity knot line crosses it. Every point of
//!   the exact edge is within `d` of the 3D chord, a triangle edge, and the
//!   image of the parameter chord within `e` of it. The face's exact domain
//!   differs from the triangulated one by the lens between the edge's
//!   preimage and the parameter chord. A line across the lens, normal to
//!   the chord, meets its boundary at two points each on the edge or on the
//!   chord, so each lens point lies on a segment whose ends map within
//!   `max(d, e)` of the 3D chord, and its image is within `q / 8` of that
//!   segment's chord. The segment is at most twice the lens's width `W`, and
//!   `W <= (d + e) / sigma` where `sigma` bounds the singular values of the
//!   surface's derivative from below over a box the preimage provably stays
//!   in (`FaceBound::lens`). The span's term is `d + e + M W^2 / 2`, and the
//!   report is unbounded by name where `sigma` cannot be shown positive.
//!   Each vertex's distance from the surface at its inverted parameters is
//!   added by the triangle bound as for any vertex. A new edge is sampled
//!   until `d + e` fits the budget, where that is cheap.
//! - Near a pole the lens cannot be had: `|S_u|` vanishes there, so the
//!   scaled box reaches the pole and `sigma` is not positive. A span into a
//!   pole, and a span near one whose lens fails, are bounded along the
//!   meridians instead (#463, ADR 0112). The span's parameter chord is a
//!   segment; its preimage runs from one end to the other (into the pole
//!   row's point `(u*, v_pole)` for a span into a pole, where the pole row
//!   closes the region). Take a point of the region between them and the
//!   meridian (`u` fixed) through it: within the region it is a segment
//!   whose two ends lie on the preimage, on the parameter chord or, for a
//!   span into a pole, on the pole row. Their surface points are edge
//!   points, within `d` of the mesh chord; points of the chord's image,
//!   within `e` of it; or the exact pole, within `o`, the pole vertex's
//!   distance from `S(u_neighbour, v_pole)`, of the chord's end. Distance
//!   to a segment is convex, so the 3D chord between the two ends is within
//!   the larger of those of the mesh chord, and the meridian arc between
//!   them is within `|S_vv| dv^2 / 8` of its chord, `dv` the segment's
//!   length in `v`. On a support turning about its pole, `S_u` is
//!   orthogonal to `S_v` and `|S_v|` does not depend on `u`, so the edge's
//!   `v` moves at most `|C'| / s` per unit of its parameter, `s` a lower
//!   bound on `|S_v|` over a box of every `u` holding the preimage: the
//!   preimage's `v` stays within `L / s` of its start, `L` the span's
//!   length, and `dv <= |v_b - v_a| + L / s`. The box is grown as the
//!   lens's is, clamped to the support's range of `v` and, for a span into
//!   a pole, to its side of the pole. The term is `d + o + dvv dv^2 / 8`
//!   into a pole and `d + e + dvv dv^2 / 8` near one; on a ruled surface's
//!   straight rulings `S_vv = 0` and it is `d + o` or `d + e`. No quarter
//!   turn is needed: the meridians bound the region however far the
//!   preimage swings in `u` as it reaches the pole. A span near a pole is
//!   bounded this way only where its ends show it off the poles (their
//!   distances from them sum to more than its length), so no span is
//!   bounded across a pole it passes; where `s` cannot be shown positive
//!   or the box is not smooth (a tangent break of the profile, #498), the
//!   face is unbounded by name ("the surface is not certifiably regular
//!   along an edge into a pole").
//! - The certified B-spline search proves uniqueness to a coarse parameter
//!   resolution; Gauss-Newton then sharpens the parameters inside the
//!   proven enclosure, and the residual left reaches the report as above.
//!
//! An edge runs from its start vertex to its end vertex in its curve's own
//! sense, and its vertices bound it (#319, #333): a curve whose whole domain
//! does not already run between them is trimmed at their parameters
//! ([`super::trim`]). Everything downstream (the shared samples, their span
//! bounds and the report) covers that span only.
//!
//! An edge with no 3D curve is straight (the topology's contract). An edge
//! shared with a face that sampled it through a p-curve is inverted at those
//! samples; its 3D chord bound is then known only for a straight edge, and a
//! curved one is reported unbounded by name.

use super::support::{PoleKind, Support};
use super::view::Pcurve;
use crate::surface_relation::RelationInverter;
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Scalar, Tolerance, Vec3};
use axiolid_nurbs::{
    CertifiedSurfaceProjectionOptions, PeriodicBSplineSurface, SurfaceInversionCertificate3,
    SurfaceInversionRefusal, MAX_CERTIFIED_SURFACE_PROJECTION_DEPTH,
    MAX_CERTIFIED_SURFACE_PROJECTION_WORK,
};
use axiolid_surface::{BSplineSurface, Surface};
use axiolid_topology::{EdgeUse, Orientation};

use super::trim::CurveSpan;
use super::{
    edge_samples, localized_parameter, surface_periods, CurvedMeshState, FaceBound, FaceContext,
    MAX_CURVED_EDGE_SEGMENTS,
};

/// A certified B-spline inversion proves its parameters unique to this
/// fraction of the surface's parameter span; Gauss-Newton then sharpens
/// them inside that enclosure. The parameters' remaining error reaches the
/// report as the vertex's distance from the surface there, so this only
/// trades search work against polishing.
const BSPLINE_PARAMETER_RESOLUTION: Scalar = 1e-4;

/// Gauss-Newton steps sharpening certified B-spline parameters.
const POLISH_STEPS: usize = 8;

/// How a face's support names a point on it.
pub(super) struct SurfaceInverter<'s> {
    target: Target<'s>,
    /// Periods along which parameters are unwrapped by continuity.
    periods: (Option<Scalar>, Option<Scalar>),
}

/// What inverts a sample: an atomic surface's own inversion, or a surface
/// relation's certified one ([`RelationInverter`], ADR 0110).
enum Target<'s> {
    Atomic(AtomicInverter<'s>),
    Relation(RelationInverter<'s>),
}

impl<'s> SurfaceInverter<'s> {
    pub(super) fn new(support: Support<'s>, tolerance: Tolerance) -> GeomResult<Self> {
        match support {
            Support::Surface(surface) => {
                let atomic = AtomicInverter::new(surface, tolerance)?;
                Ok(Self {
                    periods: atomic.periods,
                    target: Target::Atomic(atomic),
                })
            }
            Support::Relation { surface, .. } => Ok(Self {
                periods: surface.periods(),
                target: Target::Relation(RelationInverter::new(surface, tolerance)?),
            }),
        }
    }

    /// The parameters naming `point`, or a refusal by name.
    fn invert(&self, point: Vec3) -> GeomResult<Point2> {
        match &self.target {
            Target::Atomic(inverter) => inverter.invert(point),
            Target::Relation(inverter) => inverter.invert(point),
        }
    }

    /// The `v` of the pole within the linear tolerance of `point` (#463):
    /// `None` exactly where [`Self::invert`] raises no pole refusal (R1,
    /// R2, "cone apex", "sphere pole"); where it raises one, the pole's `v`,
    /// or that refusal when `point` is near the axis but not near a pole.
    pub(super) fn pole(&self, point: Vec3) -> GeomResult<Option<Scalar>> {
        match &self.target {
            Target::Atomic(inverter) => inverter.pole(point),
            Target::Relation(inverter) => inverter.pole(point),
        }
    }

    /// See [`AtomicInverter::check_spline_parameters`]. A relation
    /// evaluates exactly periodically along its periodic axes and nothing
    /// is unwrapped along the others; an edge crossing the seam of a closed
    /// swept curve without a period is refused instead
    /// ([`RelationInverter::check_spans`]).
    fn check_spline_parameters(&self, params: &[Point2], points: &[Vec3]) -> GeomResult<()> {
        match &self.target {
            Target::Atomic(inverter) => inverter.check_spline_parameters(params, points),
            Target::Relation(inverter) => inverter.check_spans(params),
        }
    }
}

/// How an atomic surface names a point on it.
struct AtomicInverter<'s> {
    surface: &'s Surface,
    tolerance: Tolerance,
    kind: InverterKind,
    /// Periods along which parameters are unwrapped by continuity.
    periods: (Option<Scalar>, Option<Scalar>),
}

enum InverterKind {
    Analytic,
    BSpline(SplineSearch),
    Periodic(Box<PeriodicBSplineSurface>, SplineSearch),
}

/// The two policies a B-spline sample is inverted under.
///
/// The certified search stops once its distance gap meets the distance
/// tolerance, and then judges uniqueness on the boxes it kept. At the
/// model's linear tolerance it stops while neighbouring boxes still tie, so
/// a point on a knot line reads as ambiguous. The `strict` pass therefore
/// runs at a thousandth of it, refining until those ties separate; a point
/// the strict pass finds off the surface but within the linear tolerance is
/// searched again at that tolerance, `model`, and accepted only when that
/// pass proves it unique. The certified search refuses a point as soon as
/// it proves it off the surface, so a point far off never refines to the
/// parameter resolution (#320). A point on the surface can still overrun the
/// strict pass's budget (where four patches meet); `model` is then tried,
/// and a strict overrun stands unless `model` accepts the point or proves
/// it off the surface.
#[derive(Clone, Copy)]
struct SplineSearch {
    strict: CertifiedSurfaceProjectionOptions,
    model: CertifiedSurfaceProjectionOptions,
}

impl<'s> AtomicInverter<'s> {
    fn new(surface: &'s Surface, tolerance: Tolerance) -> GeomResult<Self> {
        let Surface::BSpline(spline) = surface else {
            return Ok(Self {
                surface,
                tolerance,
                kind: InverterKind::Analytic,
                periods: surface_periods(surface),
            });
        };
        let options = SplineSearch {
            strict: bspline_options(spline, strict_tolerance(tolerance)?)?,
            model: bspline_options(spline, tolerance)?,
        };
        // A surface whose closed axes are not the expanded periodic form is
        // a clamped surface that meets itself: its seam points are refused
        // as ambiguous by the plain search, never resolved by a guess.
        if spline.u_closed || spline.v_closed {
            if let Ok(periodic) = PeriodicBSplineSurface::new(spline.clone()) {
                let span = |periodic: bool, (a, b): (Scalar, Scalar)| periodic.then_some(b - a);
                let periods = (
                    span(periodic.u_is_periodic(), periodic.u_domain()),
                    span(periodic.v_is_periodic(), periodic.v_domain()),
                );
                return Ok(Self {
                    surface,
                    tolerance,
                    kind: InverterKind::Periodic(Box::new(periodic), options),
                    periods,
                });
            }
        }
        Ok(Self {
            surface,
            tolerance,
            kind: InverterKind::BSpline(options),
            periods: (None, None),
        })
    }

    /// The parameters naming `point`, or a refusal by name.
    fn invert(&self, point: Vec3) -> GeomResult<Point2> {
        self.check_singularity(point)?;
        type Outcome = Result<SurfaceInversionCertificate3, SurfaceInversionRefusal>;
        let run = |options: CertifiedSurfaceProjectionOptions| -> GeomResult<Outcome> {
            match (&self.kind, self.surface) {
                (InverterKind::BSpline(_), Surface::BSpline(spline)) => {
                    axiolid_nurbs::invert_surface_certified(spline, point, options)
                }
                (InverterKind::Periodic(periodic, _), _) => {
                    axiolid_nurbs::invert_periodic_surface_certified(periodic, point, options)
                }
                _ => Err(GeomError::InvalidInput(
                    "B-spline inverter built for another surface".to_owned(),
                )),
            }
        };
        let search = match &self.kind {
            InverterKind::Analytic => {
                let (u, v) =
                    axiolid_reference::surface::invert(self.surface, point, self.tolerance)?;
                return Ok(Point2::new(u, v));
            }
            InverterKind::BSpline(search) | InverterKind::Periodic(_, search) => *search,
        };
        let outcome = match run(search.strict) {
            Ok(Err(SurfaceInversionRefusal::OffSurface {
                distance_lower_bound,
            })) if distance_lower_bound <= self.tolerance.linear() => run(search.model)?,
            Ok(outcome) => outcome,
            Err(exhausted @ GeomError::BudgetExceeded { .. }) => match run(search.model) {
                Ok(Ok(certificate)) => Ok(certificate),
                Ok(Err(off @ SurfaceInversionRefusal::OffSurface { .. })) => Err(off),
                _ => return Err(exhausted),
            },
            Err(error) => return Err(error),
        };
        let certificate = certified_parameters(outcome)?;
        self.polish(point, certificate)
    }

    /// Refine the certified parameters inside their enclosure by
    /// Gauss-Newton on `|S - p|^2`. The search proved the enclosure holds
    /// the unique minimizer; this only sharpens the representative, and
    /// the better of the two is kept. Whatever residual remains reaches the
    /// report as the vertex's distance from the surface.
    fn polish(&self, point: Vec3, certificate: SurfaceInversionCertificate3) -> GeomResult<Point2> {
        let evaluate = |uv: Point2| -> GeomResult<Scalar> {
            Ok((axiolid_reference::surface::evaluate(self.surface, uv.x, uv.y)? - point).length())
        };
        let enclosure = certificate.enclosure;
        let grow = |start: Scalar, end: Scalar| {
            let width = end - start;
            (start - width, end + width)
        };
        // Kept inside the domain on an open axis: the evaluator clamps
        // there, so a parameter past it names the same point, but the
        // derivative bounds over a box reaching past it are not smooth.
        let (u_domain, v_domain) = self.open_domain();
        let (u_lo, u_hi) = grow(enclosure.u.start, enclosure.u.end);
        let (v_lo, v_hi) = grow(enclosure.v.start, enclosure.v.end);
        let (u_lo, u_hi) = (u_lo.max(u_domain.0), u_hi.min(u_domain.1));
        let (v_lo, v_hi) = (v_lo.max(v_domain.0), v_hi.min(v_domain.1));
        if !(u_lo <= u_hi && v_lo <= v_hi) {
            return Err(GeomError::Degenerate(
                "a certified B-spline enclosure lies outside the surface's domain".to_owned(),
            ));
        }
        let mut best = Point2::new(
            certificate.u.clamp(u_lo, u_hi),
            certificate.v.clamp(v_lo, v_hi),
        );
        let mut best_residual = evaluate(best)?;
        let mut current = best;
        for _ in 0..POLISH_STEPS {
            if best_residual == 0.0 {
                break;
            }
            let on = axiolid_reference::surface::evaluate(self.surface, current.x, current.y)?;
            let (su, sv) =
                axiolid_reference::surface::partials(self.surface, current.x, current.y)?;
            let r = on - point;
            let (a, b, c) = (su.dot(su), su.dot(sv), sv.dot(sv));
            let det = a * c - b * b;
            if !(det.is_finite() && det > 0.0) {
                break;
            }
            let (gu, gv) = (su.dot(r), sv.dot(r));
            let step = Point2::new((c * gu - b * gv) / det, (a * gv - b * gu) / det);
            current = Point2::new(
                (current.x - step.x).clamp(u_lo, u_hi),
                (current.y - step.y).clamp(v_lo, v_hi),
            );
            let residual = evaluate(current)?;
            if residual < best_residual {
                best = current;
                best_residual = residual;
            } else {
                break;
            }
        }
        Ok(best)
    }

    /// The native domain of each open B-spline axis; unbounded on a
    /// periodic axis.
    fn open_domain(&self) -> ((Scalar, Scalar), (Scalar, Scalar)) {
        let unbounded = (Scalar::NEG_INFINITY, Scalar::INFINITY);
        let Surface::BSpline(spline) = self.surface else {
            return (unbounded, unbounded);
        };
        let axis = |knots: &[Scalar], period: Option<Scalar>| match (knots.first(), knots.last()) {
            (Some(&a), Some(&b)) if period.is_none() => (a, b),
            _ => unbounded,
        };
        (
            axis(&spline.u_knots, self.periods.0),
            axis(&spline.v_knots, self.periods.1),
        )
    }

    /// See [`SurfaceInverter::pole`]: a sphere's pole is at `v = +-pi/2` by
    /// the side of the centre the point lies on, a cone's apex where its
    /// radius `r + v tan(a)` vanishes (moved by an ulp or two to where it
    /// is not negative, which its evaluator refuses).
    fn pole(&self, point: Vec3) -> GeomResult<Option<Scalar>> {
        let Err(refusal) = self.check_singularity(point) else {
            return Ok(None);
        };
        let v = match self.surface {
            Surface::Sphere(sphere) => {
                let axis = sphere.frame.z.normalize_or_zero();
                if (point - sphere.frame.origin).dot(axis) >= 0.0 {
                    core::f64::consts::FRAC_PI_2
                } else {
                    -core::f64::consts::FRAC_PI_2
                }
            }
            Surface::Cone(cone) => {
                let slope = cone.semi_angle.tan();
                if !(slope.is_finite() && slope != 0.0) {
                    return Err(refusal);
                }
                let mut v = -cone.radius / slope;
                for _ in 0..4 {
                    if cone.radius + v * slope >= 0.0 {
                        break;
                    }
                    v = if slope > 0.0 {
                        v.next_up()
                    } else {
                        v.next_down()
                    };
                }
                v
            }
            _ => return Err(refusal),
        };
        let Ok(at) = axiolid_reference::surface::evaluate(self.surface, 0.0, v) else {
            return Err(refusal);
        };
        let off = (at - point).length();
        if off.is_nan() || off > self.tolerance.linear() {
            return Err(refusal);
        }
        Ok(Some(v))
    }

    /// Refuse a point at, or within the linear tolerance of, a cone's apex
    /// or a sphere's pole: every `u` names it, so the boundary's parameters
    /// there are not determined by the point.
    fn check_singularity(&self, point: Vec3) -> GeomResult<()> {
        let (frame, what) = match self.surface {
            Surface::Cone(cone) => (&cone.frame, "cone apex"),
            Surface::Sphere(sphere) => (&sphere.frame, "sphere pole"),
            _ => return Ok(()),
        };
        let offset = point - frame.origin;
        let axis = frame.z.normalize_or_zero();
        let radial = (offset - axis * offset.dot(axis)).length();
        if radial.is_nan() || radial <= self.tolerance.linear() {
            return Err(GeomError::Degenerate(format!(
                "an edge without a p-curve runs into the {what} of its face, where every u \
                 names the point: its parameters there are not determined"
            )));
        }
        Ok(())
    }
}

/// A thousandth of the linear tolerance, for [`SplineSearch::strict`].
fn strict_tolerance(tolerance: Tolerance) -> GeomResult<Tolerance> {
    Tolerance::new(tolerance.linear() * 1e-3, tolerance.angular())
        .map_err(|error| GeomError::InvalidInput(format!("strict inversion tolerance: {error:?}")))
}

fn bspline_options(
    spline: &BSplineSurface,
    tolerance: Tolerance,
) -> GeomResult<CertifiedSurfaceProjectionOptions> {
    let span = |knots: &[Scalar]| match (knots.first(), knots.last()) {
        (Some(a), Some(b)) => (b - a).abs(),
        _ => 0.0,
    };
    let extent = span(&spline.u_knots).max(span(&spline.v_knots));
    let resolution = if extent.is_finite() && extent > 0.0 {
        extent * BSPLINE_PARAMETER_RESOLUTION
    } else {
        BSPLINE_PARAMETER_RESOLUTION
    };
    CertifiedSurfaceProjectionOptions::new(
        tolerance,
        resolution,
        MAX_CERTIFIED_SURFACE_PROJECTION_WORK,
        MAX_CERTIFIED_SURFACE_PROJECTION_DEPTH,
    )
}

fn certified_parameters(
    outcome: Result<SurfaceInversionCertificate3, SurfaceInversionRefusal>,
) -> GeomResult<SurfaceInversionCertificate3> {
    match outcome {
        Ok(certificate) => Ok(certificate),
        Err(SurfaceInversionRefusal::OffSurface {
            distance_lower_bound,
        }) => Err(GeomError::Degenerate(format!(
            "an edge without a p-curve leaves its B-spline face's surface by at least \
             {distance_lower_bound:e}: its samples must lie on the surface"
        ))),
        Err(SurfaceInversionRefusal::Ambiguous { .. }) => Err(GeomError::Degenerate(
            "an edge sample without a p-curve has no unique parameters on its B-spline \
             face's surface (a pole, a seam or a patch touching itself)"
                .to_owned(),
        )),
        Err(SurfaceInversionRefusal::Unresolved { .. }) => Err(GeomError::BudgetExceeded {
            resource: "certified B-spline inversion of an edge sample",
        }),
        Err(_) => Err(GeomError::Degenerate(
            "an edge sample without a p-curve could not be inverted onto its B-spline face"
                .to_owned(),
        )),
    }
}

/// How an edge is traced in 3D, in its own start-to-end direction.
pub(super) enum EdgeTrace<'g> {
    /// No curve: the segment between its vertices.
    Straight(Vec3, Vec3),
    /// Its curve over the span its vertices bound ([`super::trim`]).
    Curve(&'g axiolid_curve::Curve3, CurveSpan),
}

impl EdgeTrace<'_> {
    /// The point at fraction `s` of the edge's parameter span.
    pub(super) fn at(&self, s: Scalar) -> GeomResult<Vec3> {
        match self {
            Self::Straight(a, b) => Ok(*a + (*b - *a) * s),
            Self::Curve(curve, span) => span.evaluate(curve, s),
        }
    }

    /// The fractions of the edge's span where its curve has a derivative
    /// break (#342): strictly inside, ascending. An edge with more pieces
    /// than the segment cap allows is not cut: it is sampled uniformly, its
    /// span bounds still split at the breaks.
    pub(super) fn breaks(&self) -> Vec<Scalar> {
        match self {
            Self::Straight(..) => Vec::new(),
            Self::Curve(_, span) => {
                let breaks = span.break_fractions();
                if breaks.len() + 1 > MAX_CURVED_EDGE_SEGMENTS {
                    Vec::new()
                } else {
                    breaks
                }
            }
        }
    }

    /// The certified bounds of the span between fractions `a` and `b`,
    /// whose ends are the points `pa` and `pb`.
    fn span(&self, a: Scalar, b: Scalar, pa: Vec3, pb: Vec3) -> Option<SpanBound> {
        let Self::Curve(curve, span) = self else {
            return Some(SpanBound::straight(pa, pb));
        };
        span.bound(curve, a, b, pa, pb)
    }
}

/// Certified bounds of one span of an edge on its 3D curve.
#[derive(Debug, Clone, Copy)]
pub(super) struct SpanBound {
    /// How far the span strays from its chord.
    pub(super) chord: Scalar,
    /// How long the span is.
    pub(super) arc: Scalar,
}

impl SpanBound {
    pub(super) fn straight(a: Vec3, b: Vec3) -> Self {
        Self {
            chord: 0.0,
            arc: (b - a).length(),
        }
    }
}

impl FaceBound<'_> {
    /// How far the surface image of the parameter chord `a`-`b` may stray
    /// from the 3D chord between `S(a)` and `S(b)`: `q(D) / 8` with the
    /// interpolation coefficients over the chord's box (either set; the
    /// smaller is kept), or the first-order bound where a knot line of full
    /// multiplicity crosses it. `None` without a derivative bound.
    fn chord_image(&self, a: Point2, b: Point2) -> Option<Scalar> {
        let (lo, hi) = (a.min(b), a.max(b));
        let bounds = self.oracle.as_ref()?.bounds((lo.x, hi.x), (lo.y, hi.y))?;
        let d = b - a;
        if !bounds.smooth {
            return Some(bounds.du * d.x.abs() + bounds.dv * d.y.abs());
        }
        let q = |[x, y, z]: [Scalar; 3]| {
            0.125 * (x * d.x * d.x + 2.0 * y * (d.x * d.y).abs() + z * d.y * d.y)
        };
        Some(q(bounds.interpolation).min(q([bounds.duu, bounds.duv, bounds.dvv])))
    }

    /// The trim terms of one edge use without a pcurve: per span, the
    /// edge's chord bound plus [`Self::chord_image`] of its parameter chord
    /// plus the lens between the two (the module notes give the argument);
    /// a span into a pole (`at_pole`, at the first and the last sample) is
    /// bounded by [`Self::pole_span`] instead.
    pub(super) fn inverted_edge(
        &mut self,
        spans: &[SpanBound],
        params: &[Point2],
        points: &[Vec3],
        at_pole: [bool; 2],
    ) -> GeomResult<()> {
        if self.unbounded.is_some() {
            return Ok(());
        }
        for (i, (span, ends)) in spans.iter().zip(params.windows(2)).enumerate() {
            if let Some((at, pole_v, pole)) = pole_side(i, params, points, at_pole) {
                let offset = (self.surface.evaluate(at.x, pole_v)? - pole).length();
                let Some(term) = self.pole_span(at, pole_v, span, offset) else {
                    self.unbounded(
                        "the surface is not certifiably regular along an edge into a pole",
                    );
                    return Ok(());
                };
                self.trim = self.trim.max(term);
                continue;
            }
            let Some(image) = self.chord_image(ends[0], ends[1]) else {
                self.unbounded("surface has no derivative bound over an edge chord");
                return Ok(());
            };
            let near = span.chord + image;
            let Some(term) = self
                .lens(ends[0], ends[1], span.arc, near)
                .map(|lens| near + lens)
                .or_else(|| {
                    self.meridian_lens(
                        (ends[0], ends[1]),
                        (points[i], points[i + 1]),
                        span.arc,
                        near,
                    )
                })
            else {
                self.unbounded(
                    "the surface is not certifiably invertible along an edge without a pcurve",
                );
                return Ok(());
            };
            self.trim = self.trim.max(term);
        }
        Ok(())
    }

    /// The lens term of one span, or `None` when the surface cannot be
    /// shown invertible along it. `arc` bounds the span's 3D length and
    /// `near` how far each edge point is from the surface image of the
    /// chord `a`-`b`.
    ///
    /// Distances are taken in parameters scaled by `alpha = |S_u(m)|` and
    /// `beta = |S_v(m)|`, `m` the chord's midpoint, where a unit step in
    /// either is about as long on the surface. Over a box `B` every
    /// singular value of the scaled derivative is at least `sigma(m) -
    /// |DS - DS(m)|` (the last by the second derivatives over `B`), and a
    /// `sigma` below that is assumed. While inside `B` the edge's preimage
    /// moves at most `|C'| / sigma` per unit of its parameter, so it never
    /// leaves `B` when `B` holds every point within `arc / sigma` of the
    /// chord's start (`B` is clamped to the support's open domain, a
    /// B-spline's or the ruled unit square, which the preimage never
    /// leaves). Inside `B`, `|S(x) - S(y)| >= sigma |x - y|`,
    /// so the preimage is within `W = near / sigma` of the chord, a segment
    /// across the lens is at most `2 W` long, and the surface strays from
    /// that segment's chord by `q(2 W) / 8` at most.
    fn lens(&self, a: Point2, b: Point2, arc: Scalar, near: Scalar) -> Option<Scalar> {
        let oracle = self.oracle.as_ref()?;
        let m = (a + b) * 0.5;
        let (su, sv) = self.surface.partials(m.x, m.y).ok()?;
        let (alpha, beta) = (su.length(), sv.length());
        if !(alpha > 0.0 && beta > 0.0 && alpha.is_finite() && beta.is_finite()) {
            return None;
        }
        // The scaled Gram matrix is [[1, c], [c, 1]]: eigenvalues 1 +- c.
        let c = (su.dot(sv) / (alpha * beta)).abs();
        let sigma_m = (1.0 - c).max(0.0).sqrt();
        let ((u0, u1), (v0, v1)) = self.surface.open_domain();
        let mut sigma = 0.5 * sigma_m;
        for _ in 0..LENS_ATTEMPTS {
            if sigma.is_nan() || sigma <= 0.0 {
                return None;
            }
            let reach = arc / sigma;
            let (gu, gv) = (reach / alpha, reach / beta);
            let lo = Point2::new((a.x.min(b.x) - gu).max(u0), (a.y.min(b.y) - gv).max(v0));
            let hi = Point2::new((a.x.max(b.x) + gu).min(u1), (a.y.max(b.y) + gv).min(v1));
            let bounds = oracle.bounds((lo.x, hi.x), (lo.y, hi.y))?;
            // TODO(#498): a box across a tangent break is never shown
            // invertible, so a face touching one without pcurves is
            // unbounded by name; bound each side of the break instead.
            if !bounds.smooth {
                return None;
            }
            let (ru, rv) = ((m.x - lo.x).max(hi.x - m.x), (m.y - lo.y).max(hi.y - m.y));
            let du = (bounds.duu * ru + bounds.duv * rv) / alpha;
            let dv = (bounds.duv * ru + bounds.dvv * rv) / beta;
            let on_box = sigma_m - (du * du + dv * dv).sqrt();
            if on_box >= sigma {
                let w = near / sigma;
                let (wu, wv) = (2.0 * w / alpha, 2.0 * w / beta);
                let q = |[x, y, z]: [Scalar; 3]| x * wu * wu + 2.0 * y * wu * wv + z * wv * wv;
                let spread = q(bounds.interpolation).min(q([bounds.duu, bounds.duv, bounds.dvv]));
                return Some(0.125 * spread);
            }
            // The box was grown for a larger `sigma` than it has: grow it
            // for the one it has and check again.
            sigma = on_box;
        }
        None
    }
}

impl FaceBound<'_> {
    /// The trim term of a span from the sample at `at` into a pole at `v =
    /// pole_v`, the end's `u` taken from `at` (#463, ADR 0112):
    /// `span.chord + offset + dvv dv^2 / 8`, `offset` how far the pole
    /// vertex lies from `S(at.u, pole_v)`. `None` where the surface cannot
    /// be shown regular along the span. The module notes give the argument.
    pub(super) fn pole_span(
        &self,
        at: Point2,
        pole_v: Scalar,
        span: &SpanBound,
        offset: Scalar,
    ) -> Option<Scalar> {
        let near = span.chord + offset;
        match self.surface.pole_kind()? {
            // Rulings are straight: `S_vv = 0`.
            PoleKind::Apex => Some(near),
            PoleKind::Rotational => {
                Some(near + self.meridian_spread(at.x, (at.y, pole_v), true, span.arc)?)
            }
        }
    }

    /// The term of a span of a support with poles by its meridians, for a
    /// span whose [`Self::lens`] cannot be had (near a pole, where `|S_u|`
    /// vanishes): `near`, the span's chord bound plus the image of its
    /// parameter chord, plus `dvv dv^2 / 8` on a rotational support and
    /// nothing along a ruled surface's straight rulings. A meridian through
    /// a lens point leaves the lens at two points on the edge's preimage or
    /// on the parameter chord, each within `near` of the mesh chord, and the
    /// module notes give the rest of the argument. The span's ends `pa` and
    /// `pb` must show the edge off the poles over it (a span through a pole
    /// is refused where a sample lands at it, never bounded across it).
    /// `None` off a support with poles or where it cannot be shown regular.
    fn meridian_lens(
        &self,
        (a, b): (Point2, Point2),
        (pa, pb): (Vec3, Vec3),
        arc: Scalar,
        near: Scalar,
    ) -> Option<Scalar> {
        let clear = self.surface.pole_distance(pa)? + self.surface.pole_distance(pb)?;
        if clear.is_nan() || clear <= arc {
            return None;
        }
        match self.surface.pole_kind()? {
            PoleKind::Apex => Some(near),
            PoleKind::Rotational => {
                Some(near + self.meridian_spread(0.5 * (a.x + b.x), (a.y, b.y), false, arc)?)
            }
        }
    }

    /// `dvv dv^2 / 8` for an edge span of 3D length at most `arc` from `v =
    /// from` to `v = to` of a rotational support, `u` near `at_u`:
    /// `dv = |to - from| + arc / s`, `s` a lower bound on `|S_v|` over a
    /// box of every `u` that holds the span's preimage, `dvv` the bound on
    /// `|S_vv|` there. The box stays in the support's [`Support::v_range`], and on
    /// `from`'s side of `to` when `to` is a pole, which the span reaches at
    /// its end only. `None` where `s` cannot be shown positive or the box
    /// has no smooth bound.
    fn meridian_spread(
        &self,
        at_u: Scalar,
        (from, to): (Scalar, Scalar),
        pole: bool,
        arc: Scalar,
    ) -> Option<Scalar> {
        let oracle = self.oracle.as_ref()?;
        let period = self.surface.periods().0?;
        let (v_lo, v_hi) = (from.min(to), from.max(to));
        let middle = 0.5 * (v_lo + v_hi);
        // `|S_v|` does not depend on `u` on a rotational support.
        let speed = self.surface.partials(at_u, middle).ok()?.1.length();
        if !(speed > 0.0 && speed.is_finite()) {
            return None;
        }
        let (mut v0, mut v1) = self.surface.v_range();
        if pole {
            if from <= to {
                v1 = v1.min(to);
            } else {
                v0 = v0.max(to);
            }
        }
        let u = (at_u - 0.5 * period, at_u + 0.5 * period);
        let mut sigma = 0.5 * speed;
        for _ in 0..LENS_ATTEMPTS {
            if sigma.is_nan() || sigma <= 0.0 {
                return None;
            }
            let reach = arc / sigma;
            let (lo, hi) = ((v_lo - reach).max(v0), (v_hi + reach).min(v1));
            let bounds = oracle.bounds(u, (lo, hi))?;
            if !bounds.smooth {
                return None;
            }
            let on_box = speed - bounds.dvv * (middle - lo).max(hi - middle);
            if on_box >= sigma {
                let dv = (v_hi - v_lo) + reach;
                let spread = 0.125 * bounds.dvv * dv * dv;
                return spread.is_finite().then_some(spread);
            }
            sigma = on_box;
        }
        None
    }
}

/// Tries at growing the lens box until the preimage provably stays in it.
const LENS_ATTEMPTS: usize = 4;

/// The native domain of a B-spline's open axes, which its preimages stay
/// in; unbounded otherwise.
pub(super) fn parameter_domain(surface: &Surface) -> ((Scalar, Scalar), (Scalar, Scalar)) {
    let unbounded = (Scalar::NEG_INFINITY, Scalar::INFINITY);
    let Surface::BSpline(spline) = surface else {
        return (unbounded, unbounded);
    };
    let axis = |knots: &[Scalar], closed: bool| match (knots.first(), knots.last()) {
        (Some(&a), Some(&b)) if !closed => (a, b),
        _ => unbounded,
    };
    (
        axis(&spline.u_knots, spline.u_closed),
        axis(&spline.v_knots, spline.v_closed),
    )
}

/// The parameters and 3D points of one edge use with no p-curve, in the
/// use's direction, both endpoints included, the mesh indices of the
/// points, shared with every face that uses the edge, and whether its first
/// and its last sample are a pole of the surface (#463).
pub(super) fn inverted_edge_use(
    state: &mut CurvedMeshState<'_>,
    ctx: &FaceContext<'_>,
    use_: &EdgeUse<Pcurve<'_>>,
    inverter: &SurfaceInverter<'_>,
    bound: &mut FaceBound<'_>,
    previous_end: Option<Point2>,
) -> GeomResult<(Vec<Point2>, Vec<u32>, [bool; 2])> {
    let edge = ctx
        .brep
        .edges()
        .get(use_.edge.index())
        .ok_or_else(|| GeomError::InvalidInput("edge missing".to_string()))?;
    let (start_vertex, end_vertex) = if use_.orientation == Orientation::Forward {
        (edge.start, edge.end)
    } else {
        (edge.end, edge.start)
    };
    let trace = edge_trace(ctx, edge)?;
    let (points, mut params, spans, at_pole) = match state.edge_samples.get(&use_.edge) {
        Some(cached) => {
            let (indices, spans) = cached.in_direction(use_.orientation);
            let points: Vec<Vec3> = indices
                .iter()
                .map(|&index| state.mesh.positions[index as usize])
                .collect();
            let last = points.len().checked_sub(1).ok_or_else(|| {
                GeomError::Degenerate("curved edge needs both endpoint samples".to_owned())
            })?;
            // An end at a pole is not inverted: only its `v` is read (#463).
            let poles = [inverter.pole(points[0])?, inverter.pole(points[last])?];
            let at_pole = poles.map(|pole| pole.is_some());
            if last == 1 && at_pole == [true, true] {
                return Err(GeomError::Degenerate(
                    "an edge without a p-curve runs from a pole to a pole in one span: no end \
                     has a neighbouring sample to take its u from"
                        .to_owned(),
                ));
            }
            let mut params = points
                .iter()
                .enumerate()
                .map(|(k, &point)| match (k, poles) {
                    (0, [Some(v), _]) => Ok(Point2::new(0.0, v)),
                    (k, [_, Some(v)]) if k == last => Ok(Point2::new(0.0, v)),
                    _ => inverter.invert(point),
                })
                .collect::<GeomResult<Vec<_>>>()?;
            pole_ends(&mut params, at_pole);
            unwrap(inverter, &mut params);
            if !turns_within_quarter(inverter, &params) {
                return Err(singular_between_samples());
            }
            let straight = match trace {
                EdgeTrace::Straight(..) => true,
                EdgeTrace::Curve(curve, _) => straight_curve(curve),
            };
            let spans = match spans {
                Ok(spans) => Ok(spans),
                Err(_) if straight => Ok(points
                    .windows(2)
                    .map(|ends| SpanBound::straight(ends[0], ends[1]))
                    .collect()),
                Err(reason) => Err(reason),
            };
            (points, params, spans, at_pole)
        }
        None => {
            let (mut points, mut params, spans, _) =
                sample_edge(&trace, inverter, bound, ctx, &trace.breaks())?;
            let mut spans = spans.ok_or(NO_CHORD_BOUND);
            if use_.orientation == Orientation::Reversed {
                points.reverse();
                params.reverse();
                if let Ok(spans) = spans.as_mut() {
                    spans.reverse();
                }
            }
            let ends = [points[0], points[points.len() - 1]];
            let at_pole = [
                inverter.pole(ends[0])?.is_some(),
                inverter.pole(ends[1])?.is_some(),
            ];
            (points, params, spans, at_pole)
        }
    };
    let n = points.len().checked_sub(1).ok_or_else(|| {
        GeomError::Degenerate("curved edge needs both endpoint samples".to_owned())
    })?;
    state.reserve_face_vertices(n)?;
    align(inverter, &mut params, previous_end);
    inverter.check_spline_parameters(&params, &points)?;
    match &spans {
        Ok(spans) => bound.inverted_edge(spans, &params, &points, at_pole)?,
        Err(reason) => bound.unbounded(reason),
    }
    let shared = edge_samples(
        state,
        use_.edge,
        use_.orientation,
        start_vertex,
        end_vertex,
        &points,
        spans,
    )?;
    Ok((params, shared, at_pole))
}

/// Why a curved edge sampled through a pcurve has no span bounds on its 3D
/// curve.
pub(super) const THROUGH_PCURVE: &str =
    "curved edge without a pcurve, sampled through a neighbour's pcurve, has no chord bound \
     on its 3D curve";

/// Why an edge's curve family has no span bounds.
pub(super) const NO_CHORD_BOUND: &str = "edge curve family with no chord bound";

pub(super) fn edge_trace<'g>(
    ctx: &FaceContext<'g>,
    edge: &axiolid_topology::Edge<&'g axiolid_curve::Curve3>,
) -> GeomResult<EdgeTrace<'g>> {
    let position = |vertex: axiolid_topology::VertexId| -> GeomResult<Vec3> {
        Ok(ctx
            .brep
            .vertices()
            .get(vertex.index())
            .ok_or_else(|| GeomError::InvalidInput("vertex missing".to_string()))?
            .position)
    };
    let (start, end) = (position(edge.start)?, position(edge.end)?);
    let Some(curve) = edge.curve else {
        return Ok(EdgeTrace::Straight(start, end));
    };
    let span =
        super::trim::vertex_span(curve, (edge.start, start), (edge.end, end), ctx.tolerance)?;
    Ok(EdgeTrace::Curve(curve, span))
}

pub(super) fn straight_curve(curve: &axiolid_curve::Curve3) -> bool {
    match curve {
        axiolid_curve::Curve3::Line(_) => true,
        axiolid_curve::Curve3::Polyline(p) => p.points.len() == 2 && !p.closed,
        _ => false,
    }
}

/// Sample a new edge in its own direction, cut at its derivative breaks
/// (#342): uniform spans inside each piece between breaks, a power of two
/// of them, each piece doubled until every one of its spans' 3D midpoints is
/// within the budget of its chord and its parameters turn less than a
/// quarter period, then while the spans' certified terms miss the budget.
/// A break is a sample, so a straight piece (a polyline segment) is one
/// span with a chord bound of zero. The pieces' span counts together stay
/// within the segment cap. Returns the points, their parameters, the span
/// bounds and the span count of each piece.
///
/// TODO(#497): the slowest stretch of a piece sets the spacing of all of it;
/// grade the samples by the surface's local speed (ADR 0110).
#[allow(clippy::type_complexity)]
pub(super) fn sample_edge(
    trace: &EdgeTrace<'_>,
    inverter: &SurfaceInverter<'_>,
    bound: &FaceBound<'_>,
    ctx: &FaceContext<'_>,
    breaks: &[Scalar],
) -> GeomResult<(Vec<Vec3>, Vec<Point2>, Option<Vec<SpanBound>>, Vec<usize>)> {
    let budget = ctx.chord_error;
    if !(budget.is_finite() && budget > 0.0) {
        return Err(GeomError::InvalidInput(
            "curved-edge chord tolerance must be positive and finite".to_owned(),
        ));
    }
    if breaks.len() >= MAX_CURVED_EDGE_SEGMENTS {
        return Err(GeomError::BudgetExceeded {
            resource: "curved edge samples",
        });
    }
    let ends: Vec<Scalar> = core::iter::once(0.0)
        .chain(breaks.iter().copied())
        .chain(core::iter::once(1.0))
        .collect();
    let corners = ends
        .iter()
        .map(|&s| trace.at(s))
        .collect::<GeomResult<Vec<_>>>()?;
    // The edge's ends may lie at a pole, where only `v` is read (#463).
    let last = corners.len() - 1;
    let poles = [inverter.pole(corners[0])?, inverter.pole(corners[last])?];
    let corner_params = corners
        .iter()
        .enumerate()
        .map(|(k, &point)| match (k, poles) {
            (0, [Some(v), _]) => Ok(Point2::new(0.0, v)),
            (k, [_, Some(v)]) if k == last => Ok(Point2::new(0.0, v)),
            _ => inverter.invert(point),
        })
        .collect::<GeomResult<Vec<_>>>()?;
    let at_pole = poles.map(|pole| pole.is_some());
    // Each piece's points and their raw parameters, both ends included.
    let mut pieces: Vec<Piece> = (0..ends.len() - 1)
        .map(|k| Piece {
            points: vec![corners[k], corners[k + 1]],
            raw: vec![corner_params[k], corner_params[k + 1]],
        })
        .collect();
    // Per piece, the coarsest sampling meeting the midpoint and quarter-turn
    // tests, kept as a stride into the finer samples should the certified
    // terms never fit.
    let mut basic: Vec<Option<usize>> = vec![None; pieces.len()];
    loop {
        let counts: Vec<usize> = pieces.iter().map(Piece::spans).collect();
        let fractions = piece_fractions(breaks, &counts);
        let (points, mut params) = joined(&pieces);
        pole_ends(&mut params, at_pole);
        unwrap(inverter, &mut params);
        let spans = spans_of(trace, &fractions, &points);
        // The new midpoints, at the odd fractions of the doubled sampling, so
        // the samples and the bounds name the same fractions.
        let doubled: Vec<usize> = counts.iter().map(|n| 2 * n).collect();
        let fine = piece_fractions(breaks, &doubled);
        let mut mids = Vec::with_capacity(points.len() - 1);
        let mut unfinished: Vec<usize> = Vec::new();
        let mut midpoint_failed = vec![false; pieces.len()];
        // The pieces whose current sampling meets the midpoint and
        // quarter-turn tests.
        let mut passed = vec![false; pieces.len()];
        let mut offset = 0;
        for (k, &n) in counts.iter().enumerate() {
            let mut measured: Scalar = 0.0;
            for i in offset..offset + n {
                let mid = trace.at(fine[2 * i + 1])?;
                measured = measured.max((mid - 0.5 * (points[i] + points[i + 1])).length());
                mids.push(mid);
            }
            midpoint_failed[k] = measured > budget;
            let quarter = turns_within_quarter(inverter, &params[offset..=offset + n]);
            // One span between two poles has no neighbour to take `u` from.
            let lone = n == 1 && counts.len() == 1 && at_pole == [true, true];
            let mut fits = false;
            if !midpoint_failed[k] && quarter && !lone {
                passed[k] = true;
                basic[k].get_or_insert(n);
                let certified = spans.as_ref().and_then(|spans| {
                    (offset..offset + n)
                        .map(|i| {
                            if let Some((at, pole_v, pole)) =
                                pole_side(i, &params, &points, at_pole)
                            {
                                let offset =
                                    (bound.surface.evaluate(at.x, pole_v).ok()? - pole).length();
                                // Not shown regular yet: split as the lens is.
                                return Some(
                                    bound
                                        .pole_span(at, pole_v, &spans[i], offset)
                                        .unwrap_or(Scalar::INFINITY),
                                );
                            }
                            let (a, b) = (params[i], params[i + 1]);
                            let near = spans[i].chord + bound.chord_image(a, b)?;
                            // A span too long to show the surface invertible
                            // along it is split like one that misses the
                            // budget.
                            Some(
                                bound
                                    .lens(a, b, spans[i].arc, near)
                                    .map(|lens| near + lens)
                                    .or_else(|| {
                                        bound.meridian_lens(
                                            (a, b),
                                            (points[i], points[i + 1]),
                                            spans[i].arc,
                                            near,
                                        )
                                    })
                                    .unwrap_or(Scalar::INFINITY),
                            )
                        })
                        .try_fold(0.0, |worst: Scalar, term| term.map(|t| worst.max(t)))
                });
                // No certified term to steer by: the report names it.
                fits = certified.is_none_or(|worst| worst <= budget);
            }
            if !fits {
                unfinished.push(k);
            }
            offset += n;
        }
        if unfinished.is_empty() {
            return Ok((points, params, spans, counts));
        }
        let added: usize = unfinished.iter().map(|&k| counts[k]).sum();
        let mut start = 0;
        let mut starts = Vec::with_capacity(counts.len());
        for &n in &counts {
            starts.push(start);
            start += n;
        }
        // A span into a pole is halved until its middle would lie within
        // the tolerance of the pole, where it names no `u` (#463): the edge
        // is then as fine as it can be sampled.
        let last_piece = counts.len() - 1;
        let mut crowded = false;
        for &k in &unfinished {
            let ends = [
                (k == 0 && at_pole[0]).then(|| mids[starts[k]]),
                (k == last_piece && at_pole[1]).then(|| mids[starts[k] + counts[k] - 1]),
            ];
            for mid in ends.into_iter().flatten() {
                crowded |= !matches!(inverter.pole(mid), Ok(None));
            }
        }
        if crowded || points.len() - 1 + added > MAX_CURVED_EDGE_SEGMENTS {
            // Out of segments: each unfinished piece falls back to its
            // coarsest passing sampling, or the edge is refused. Crowded
            // against a pole, a piece whose sampling passes keeps it: it is
            // as fine as the pole lets it be, and its terms the smallest.
            for &k in &unfinished {
                if crowded && passed[k] {
                    continue;
                }
                let Some(coarse) = basic[k] else {
                    if midpoint_failed[k] {
                        return Err(GeomError::BudgetExceeded {
                            resource: "curved edge samples",
                        });
                    }
                    return Err(singular_between_samples());
                };
                pieces[k].coarsen(counts[k] / coarse);
            }
            let counts: Vec<usize> = pieces.iter().map(Piece::spans).collect();
            let fractions = piece_fractions(breaks, &counts);
            let (points, mut params) = joined(&pieces);
            pole_ends(&mut params, at_pole);
            unwrap(inverter, &mut params);
            let spans = spans_of(trace, &fractions, &points);
            return Ok((points, params, spans, counts));
        }
        // Halve every span of the unfinished pieces, inverting only the new
        // midpoints.
        for &k in &unfinished {
            pieces[k].halve(&mids[starts[k]..starts[k] + counts[k]], inverter)?;
        }
    }
}

/// Give an edge's end at a pole the `u` of its neighbouring sample (#463):
/// every `u` names the pole, and the neighbour's runs along the meridian
/// (or ruling) the edge's last span is bounded along. Run after the
/// samples are joined and before they are unwrapped.
fn pole_ends(params: &mut [Point2], at_pole: [bool; 2]) {
    let n = params.len() - 1;
    if at_pole[0] {
        params[0].x = params[1].x;
    }
    if at_pole[1] {
        params[n].x = params[n - 1].x;
    }
}

/// The neighbouring sample's parameters, the pole's `v` and the pole's
/// point, when span `i` of an edge ends at a pole.
fn pole_side(
    i: usize,
    params: &[Point2],
    points: &[Vec3],
    at_pole: [bool; 2],
) -> Option<(Point2, Scalar, Vec3)> {
    let n = params.len() - 1;
    if i == 0 && at_pole[0] {
        Some((params[1], params[0].y, points[0]))
    } else if i + 1 == n && at_pole[1] {
        Some((params[n - 1], params[n].y, points[n]))
    } else {
        None
    }
}

/// The samples of one piece of an edge between breaks.
struct Piece {
    /// Both ends included.
    points: Vec<Vec3>,
    /// The points' parameters on the face's surface, before unwrapping.
    raw: Vec<Point2>,
}

impl Piece {
    fn spans(&self) -> usize {
        self.points.len() - 1
    }

    /// Keep every `stride`-th sample.
    fn coarsen(&mut self, stride: usize) {
        let keep = |i: &usize| i.is_multiple_of(stride);
        let kept = |len: usize| (0..len).filter(keep).collect::<Vec<_>>();
        self.points = kept(self.points.len())
            .iter()
            .map(|&i| self.points[i])
            .collect();
        self.raw = kept(self.raw.len()).iter().map(|&i| self.raw[i]).collect();
    }

    /// Insert the `mids` of its spans.
    fn halve(&mut self, mids: &[Vec3], inverter: &SurfaceInverter<'_>) -> GeomResult<()> {
        let n = self.spans();
        let mut points = Vec::with_capacity(2 * n + 1);
        let mut raw = Vec::with_capacity(2 * n + 1);
        for (i, &mid) in mids.iter().enumerate() {
            points.push(self.points[i]);
            raw.push(self.raw[i]);
            points.push(mid);
            raw.push(inverter.invert(mid)?);
        }
        points.push(self.points[n]);
        raw.push(self.raw[n]);
        self.points = points;
        self.raw = raw;
        Ok(())
    }
}

/// The points and raw parameters of `pieces`, a shared end once.
fn joined(pieces: &[Piece]) -> (Vec<Vec3>, Vec<Point2>) {
    let (mut points, mut raw) = (Vec::new(), Vec::new());
    for (k, piece) in pieces.iter().enumerate() {
        let skip = usize::from(k > 0);
        points.extend_from_slice(&piece.points[skip..]);
        raw.extend_from_slice(&piece.raw[skip..]);
    }
    (points, raw)
}

/// The fractions of an edge's samples when it is cut at `breaks` into pieces
/// of `counts` uniform spans: every break exactly, and the ends 0 and 1.
pub(super) fn piece_fractions(breaks: &[Scalar], counts: &[usize]) -> Vec<Scalar> {
    debug_assert_eq!(counts.len(), breaks.len() + 1);
    let end = |k: usize| match k {
        0 => 0.0,
        _ => breaks.get(k - 1).copied().unwrap_or(1.0),
    };
    let mut fractions = Vec::with_capacity(counts.iter().sum::<usize>() + 1);
    for (k, &n) in counts.iter().enumerate() {
        let (from, to) = (end(k), end(k + 1));
        fractions.extend((0..n).map(|i| from + (to - from) * (i as Scalar / n as Scalar)));
    }
    fractions.push(1.0);
    fractions
}

/// The bounds of the spans between `points`, sampled at `fractions` of the
/// edge's parameter; `None` where a family has none.
pub(super) fn spans_of(
    trace: &EdgeTrace<'_>,
    fractions: &[Scalar],
    points: &[Vec3],
) -> Option<Vec<SpanBound>> {
    debug_assert_eq!(fractions.len(), points.len());
    (0..points.len().checked_sub(1)?)
        .map(|i| trace.span(fractions[i], fractions[i + 1], points[i], points[i + 1]))
        .collect()
}

/// Unwrap periodic parameters along one edge by continuity.
fn unwrap(inverter: &SurfaceInverter<'_>, params: &mut [Point2]) {
    let (pu, pv) = inverter.periods;
    for i in 1..params.len() {
        params[i] = Point2::new(
            localized_parameter(params[i - 1].x, params[i].x, pu),
            localized_parameter(params[i - 1].y, params[i].y, pv),
        );
    }
}

/// Whether every span of an unwrapped edge turns less than a quarter
/// period: only then is continuity a reading of the edge rather than a
/// guess across a singularity.
fn turns_within_quarter(inverter: &SurfaceInverter<'_>, params: &[Point2]) -> bool {
    let (pu, pv) = inverter.periods;
    let within = |d: Scalar, period: Option<Scalar>| period.is_none_or(|p| d.abs() < 0.25 * p);
    params
        .windows(2)
        .all(|w| within(w[1].x - w[0].x, pu) && within(w[1].y - w[0].y, pv))
}

fn singular_between_samples() -> GeomError {
    GeomError::Degenerate(
        "an edge without a p-curve turns a quarter period or more between consecutive \
         samples: it runs into a singularity of its face's surface"
            .to_owned(),
    )
}

/// Move the edge's parameters by whole periods so it starts where the
/// previous edge of the ring ended.
fn align(inverter: &SurfaceInverter<'_>, params: &mut [Point2], previous_end: Option<Point2>) {
    let (Some(end), Some(&first)) = (previous_end, params.first()) else {
        return;
    };
    let (pu, pv) = inverter.periods;
    let shift = Point2::new(
        localized_parameter(end.x, first.x, pu),
        localized_parameter(end.y, first.y, pv),
    ) - first;
    for p in params.iter_mut() {
        *p += shift;
    }
}

impl AtomicInverter<'_> {
    /// A B-spline's evaluator does not wrap a parameter moved by whole
    /// periods back into its domain, so each unwrapped sample must still
    /// name its point there.
    fn check_spline_parameters(&self, params: &[Point2], points: &[Vec3]) -> GeomResult<()> {
        if matches!(self.kind, InverterKind::Analytic) {
            return Ok(());
        }
        for (uv, point) in params.iter().zip(points) {
            let on = axiolid_reference::surface::evaluate(self.surface, uv.x, uv.y)?;
            let off = (on - *point).length();
            if off.is_nan() || off > self.tolerance.linear() {
                return Err(GeomError::Degenerate(
                    "an edge without a p-curve crosses the seam of its periodic B-spline \
                     face: its unwrapped parameters leave the surface's domain"
                        .to_owned(),
                ));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pieces of power-of-two counts between breaks put a sample on every
    /// break, bit for bit, and on the ends.
    #[test]
    fn piece_fractions_land_on_the_breaks() {
        let breaks = [0.1, 0.35];
        let fractions = piece_fractions(&breaks, &[1, 4, 2]);
        assert_eq!(fractions.len(), 8);
        assert_eq!(fractions[0], 0.0);
        assert_eq!(fractions[1], 0.1);
        assert_eq!(fractions[5], 0.35);
        assert_eq!(fractions[7], 1.0);
        assert!(fractions.windows(2).all(|w| w[0] < w[1]));
        // No breaks is the plain uniform sampling.
        assert_eq!(piece_fractions(&[], &[4]), vec![0.0, 0.25, 0.5, 0.75, 1.0]);
    }

    /// Continuity reads an edge only while each span turns less than a
    /// quarter period.
    #[test]
    fn the_quarter_turn_guard_is_a_quarter_period() {
        let surface = Surface::Cylinder(axiolid_surface::Cylinder {
            frame: axiolid_core::Frame3 {
                origin: Vec3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius: 0.3,
        });
        let inverter =
            SurfaceInverter::new(Support::Surface(&surface), Tolerance::MILLIMETRE).unwrap();
        let quarter = core::f64::consts::FRAC_PI_2;
        let turn = |d: Scalar| [Point2::new(0.0, 0.0), Point2::new(d, 0.1)];
        assert!(turns_within_quarter(&inverter, &turn(0.99 * quarter)));
        assert!(!turns_within_quarter(&inverter, &turn(1.01 * quarter)));
        assert!(!turns_within_quarter(&inverter, &turn(-1.01 * quarter)));
    }

    /// `pole` is `None` exactly where `invert` raises no "cone apex" or
    /// "sphere pole" refusal, and the pole's `v` within the tolerance of it.
    #[test]
    fn a_pole_is_named_exactly_where_a_sample_is_refused_there() {
        let frame = axiolid_core::Frame3 {
            origin: Vec3::new(0.2, -0.1, 0.3),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        };
        let tau = Tolerance::MILLIMETRE.linear();
        let check =
            |name: &str, inverter: &SurfaceInverter<'_>, point: Vec3, want: Option<Scalar>| {
                let pole = inverter.pole(point);
                let at_pole = matches!(inverter.invert(point), Err(GeomError::Degenerate(text))
                if text.contains("cone apex") || text.contains("sphere pole"));
                assert_eq!(!matches!(pole, Ok(None)), at_pole, "{name}: {pole:?}");
                if let Some(want) = want {
                    assert_eq!(pole.unwrap(), Some(want), "{name}");
                }
            };
        let off = |k: Scalar| Vec3::new(0.6, 0.8, 0.0) * (k * tau);
        let sphere = Surface::Sphere(axiolid_surface::Sphere { frame, radius: 0.5 });
        let inverter =
            SurfaceInverter::new(Support::Surface(&sphere), Tolerance::MILLIMETRE).unwrap();
        for (z, v) in [
            (0.5, core::f64::consts::FRAC_PI_2),
            (-0.5, -core::f64::consts::FRAC_PI_2),
        ] {
            let pole = frame.origin + Vec3::Z * z;
            check("sphere pole", &inverter, pole, Some(v));
            check("near the pole", &inverter, pole + off(0.6), Some(v));
            check("beyond it", &inverter, pole + off(1.5), None);
        }
        check("the centre", &inverter, frame.origin, None);
        assert!(inverter.pole(frame.origin).is_err());

        let (radius, semi_angle) = (0.4, -0.4_f64);
        let cone = Surface::Cone(axiolid_surface::Cone {
            frame,
            radius,
            semi_angle,
        });
        let inverter =
            SurfaceInverter::new(Support::Surface(&cone), Tolerance::MILLIMETRE).unwrap();
        let apex_v = -radius / semi_angle.tan();
        let apex = frame.origin + Vec3::Z * apex_v;
        check("cone apex", &inverter, apex, Some(apex_v));
        check("near the apex", &inverter, apex + off(0.6), None);
        check("beyond it", &inverter, apex + off(1.5), None);
        check("on the axis, off", &inverter, frame.origin, None);
        assert!(inverter.pole(frame.origin).is_err());
        let v = inverter.pole(apex).unwrap().unwrap();
        assert!(radius + v * semi_angle.tan() >= 0.0);
    }

    /// I6 against dense closed-form sampling: every surface point between
    /// a small circle's last span into a sphere's north pole and the
    /// meridian the span is meshed along (its end takes its neighbour's
    /// `u`), closed by the pole row, lies within the span's term of the
    /// mesh chord. The circles are tilted from a meridian by `alpha`, so
    /// their preimage swings in `u` as they reach the pole. The term needs
    /// its chord bound: the rest of it falls short of a measured point. (On
    /// a sphere no lens point lies farther from the chord than the edge
    /// does, as every curve on it bends at least as much as its meridians;
    /// the meridian term is what the argument needs on any other support.)
    #[test]
    fn a_span_into_a_pole_bounds_its_lens_along_the_meridians() {
        let radius = 1.0;
        let sphere = Surface::Sphere(axiolid_surface::Sphere {
            frame: axiolid_core::Frame3 {
                origin: Vec3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius,
        });
        let bound = FaceBound::new(Support::Surface(&sphere));
        let north = Vec3::Z * radius;
        let segment = |p: Vec3, a: Vec3, b: Vec3| {
            let d = b - a;
            let t = ((p - a).dot(d) / d.dot(d)).clamp(0.0, 1.0);
            (p - (a + d * t)).length()
        };
        let invert = |p: Vec3| Point2::new(p.y.atan2(p.x), (p.z / radius).clamp(-1.0, 1.0).asin());
        let pole_v = core::f64::consts::FRAC_PI_2;
        let mut beyond_chord = false;
        for alpha in [0.1_f64, 0.4, 0.8, 1.2] {
            // The circle through the pole in the plane of normal `n`.
            let n = Vec3::new(alpha.sin(), 0.0, alpha.cos());
            let centre = n * n.dot(north);
            let rho = (radius * radius - centre.length_squared()).sqrt();
            let e1 = (north - centre).normalize();
            let e2 = n.cross(e1);
            let circle = |theta: Scalar| centre + (e1 * theta.cos() + e2 * theta.sin()) * rho;
            for h in [0.2, 0.05] {
                let arc = rho * h;
                let start = circle(h);
                let at = invert(start);
                let span = SpanBound {
                    chord: rho * (1.0 - (0.5 * h).cos()),
                    arc,
                };
                let term = bound.pole_span(at, pole_v, &span, 0.0).unwrap();
                let mut measured: Scalar = 0.0;
                for i in 0..=200 {
                    let mut uv = invert(circle(h * (1.0 - i as Scalar / 200.0)));
                    if i == 200 {
                        uv = Point2::new(at.x, pole_v);
                    }
                    let u = localized_parameter(at.x, uv.x, Some(core::f64::consts::TAU));
                    for j in 0..=40 {
                        let t = j as Scalar / 40.0;
                        let q = axiolid_reference::surface::evaluate(
                            &sphere,
                            at.x + (u - at.x) * t,
                            uv.y,
                        )
                        .unwrap();
                        measured = measured.max(segment(q, start, north));
                    }
                }
                assert!(
                    measured <= term,
                    "alpha {alpha}, h {h}: {measured:e} above {term:e}"
                );
                beyond_chord |= measured > term - span.chord;
            }
        }
        assert!(beyond_chord, "the chord bound is never needed");
    }

    /// The bound along the meridians against dense closed-form sampling on
    /// a curved profile, the canopy's cubic from the axis to the axis
    /// turned about it: a span near the top pole oblique to the meridians
    /// (its preimage bowed away from its parameter chord), and a span into
    /// the pole whose preimage swings in `u`. Every surface point between
    /// the preimage and the chord (and the pole row) lies within the term
    /// of the mesh chord. The meridian term is not what holds it here
    /// either: the farthest points measured are the edge's own (6.7e-3
    /// against `d + e` = 1.8e-2 near the pole, 2.21e-3 against `d` =
    /// 2.23e-3 into it, for the widest spans), so leaving the term out
    /// shows in no sample; the argument needs it.
    #[test]
    fn a_span_near_a_curved_profile_s_pole_is_bounded_along_the_meridians() {
        use axiolid_model::{GeometryGraphBuilder, GeometryNode, SurfaceRelation};
        let profile = axiolid_curve::Curve3::BSpline(axiolid_curve::BSplineCurve3 {
            degree: 3,
            control_points: vec![
                Vec3::ZERO,
                Vec3::new(0.6, 0.0, 0.1),
                Vec3::new(0.5, 0.0, 0.9),
                Vec3::new(0.0, 0.0, 1.0),
            ],
            knots: vec![0.0, 1.0],
            multiplicities: vec![4, 4],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: axiolid_curve::KnotSpec::Unspecified,
        });
        let mut builder = GeometryGraphBuilder::new();
        let curve = builder.push(GeometryNode::Curve3(profile)).unwrap();
        let graph = builder.finish(vec![curve]).unwrap();
        let surface = crate::surface_relation::RelationSurface::revolution(
            &graph,
            curve,
            Vec3::ZERO,
            Vec3::Z,
        )
        .unwrap();
        let node = SurfaceRelation::LinearExtrusion {
            swept_curve: curve,
            direction: Vec3::Z,
        };
        let support = Support::Relation {
            node: &node,
            surface,
        };
        let bound = FaceBound::new(support);
        let at = |p: Point2| surface.evaluate(p.x, p.y).unwrap();
        let segment = |p: Vec3, a: Vec3, b: Vec3| {
            let d = b - a;
            let t = ((p - a).dot(d) / d.dot(d)).clamp(0.0, 1.0);
            (p - (a + d * t)).length()
        };
        // The edge's chord bound and length, read densely and widened.
        let measure = |path: &dyn Fn(Scalar) -> Point2| {
            let (pa, pb) = (at(path(0.0)), at(path(1.0)));
            let (mut chord, mut arc): (Scalar, Scalar) = (0.0, 0.0);
            let mut previous = pa;
            for i in 1..=2000 {
                let p = at(path(i as Scalar / 2000.0));
                chord = chord.max(segment(p, pa, pb));
                arc += (p - previous).length();
                previous = p;
            }
            SpanBound {
                chord: chord * 1.01,
                arc: arc * 1.01,
            }
        };
        for (du, bow) in [(0.4, 0.01), (0.2, 0.005), (0.1, 0.0025)] {
            // Near the pole, oblique to the meridians.
            let (a, b) = (Point2::new(0.2, 0.9), Point2::new(0.2 + du, 0.93));
            let chord = move |s: Scalar| a + (b - a) * s;
            let path = move |s: Scalar| {
                chord(s) + Point2::new(0.0, bow * (core::f64::consts::PI * s).sin())
            };
            let span = measure(&path);
            let near = span.chord + bound.chord_image(a, b).unwrap();
            let (pa, pb) = (at(a), at(b));
            let term = bound
                .meridian_lens((a, b), (pa, pb), span.arc, near)
                .unwrap();
            let mut measured: Scalar = 0.0;
            for i in 0..=200 {
                let s = i as Scalar / 200.0;
                for j in 0..=40 {
                    let t = j as Scalar / 40.0;
                    let q = at(chord(s) + (path(s) - chord(s)) * t);
                    measured = measured.max(segment(q, pa, pb));
                }
            }
            assert!(measured <= term, "near, {du}: {measured:e} above {term:e}");
            // Into the pole at `v = 1`, from `v = 0.95`, sweeping `du` in
            // `u` as it goes; meshed along the meridian at its start.
            let start = Point2::new(0.3, 1.0 - 0.05 * (du / 0.8).min(1.0));
            let path =
                move |s: Scalar| Point2::new(start.x + du * s * s, start.y + (1.0 - start.y) * s);
            let span = measure(&path);
            let pole = at(Point2::new(start.x, 1.0));
            let term = bound.pole_span(start, 1.0, &span, 0.0).unwrap();
            let x1 = at(start);
            let mut measured: Scalar = 0.0;
            for i in 0..=200 {
                let p = path(i as Scalar / 200.0);
                for j in 0..=40 {
                    let t = j as Scalar / 40.0;
                    let q = at(Point2::new(start.x + (p.x - start.x) * t, p.y));
                    measured = measured.max(segment(q, x1, pole));
                }
            }
            assert!(measured <= term, "into, {du}: {measured:e} above {term:e}");
        }
    }

    /// The certified search proves the B-spline parameters unique to a
    /// coarse resolution; the polish brings them to rounding.
    #[test]
    fn certified_bspline_parameters_are_polished() {
        let net: Vec<Vec<Vec3>> = (0..4)
            .map(|i| {
                (0..4)
                    .map(|j| {
                        let (x, y) = (i as Scalar / 3.0, j as Scalar / 3.0);
                        Vec3::new(x, y, 0.2 * (x - 0.5) * (y - 0.3) + 0.1 * x * x)
                    })
                    .collect()
            })
            .collect();
        let surface = Surface::BSpline(BSplineSurface {
            u_degree: 3,
            v_degree: 3,
            control_points: net,
            u_knots: vec![0.0, 1.0],
            u_multiplicities: vec![4, 4],
            v_knots: vec![0.0, 1.0],
            v_multiplicities: vec![4, 4],
            weights: None,
            u_closed: false,
            v_closed: false,
            knot_spec: axiolid_curve::KnotSpec::Unspecified,
            self_intersect: None,
        });
        let inverter =
            SurfaceInverter::new(Support::Surface(&surface), Tolerance::MILLIMETRE).unwrap();
        for (u, v) in [(0.31, 0.77), (0.5, 0.0), (0.9, 0.42)] {
            let point = axiolid_reference::surface::evaluate(&surface, u, v).unwrap();
            let uv = inverter.invert(point).unwrap();
            let on = axiolid_reference::surface::evaluate(&surface, uv.x, uv.y).unwrap();
            assert!((on - point).length() < 1e-12, "{u} {v}: {uv:?}");
        }
    }

    /// Where four Bezier patches meet the strict search can exhaust its
    /// budget although the point is on the surface; the model-tolerance
    /// search proves it unique, so the sample is accepted (#320).
    #[test]
    fn a_strict_budget_overrun_is_retried_at_the_model_tolerance() {
        let net: Vec<Vec<Vec3>> = (0..5)
            .map(|i| {
                (0..5)
                    .map(|j| {
                        let (x, y) = (i as Scalar * 0.25, j as Scalar * 0.25);
                        let lift = if (1..4).contains(&i) && (1..4).contains(&j) {
                            0.12
                        } else {
                            0.0
                        };
                        Vec3::new(x, y, lift + 0.05 * x * y)
                    })
                    .collect()
            })
            .collect();
        let bump = Surface::BSpline(BSplineSurface {
            u_degree: 3,
            v_degree: 3,
            control_points: net,
            u_knots: vec![0.0, 0.5, 1.0],
            u_multiplicities: vec![4, 1, 4],
            v_knots: vec![0.0, 0.5, 1.0],
            v_multiplicities: vec![4, 1, 4],
            weights: None,
            u_closed: false,
            v_closed: false,
            knot_spec: axiolid_curve::KnotSpec::Unspecified,
            self_intersect: None,
        });
        let inverter =
            SurfaceInverter::new(Support::Surface(&bump), Tolerance::new(1e-4, 1e-9).unwrap())
                .unwrap();
        let point = axiolid_reference::surface::evaluate(&bump, 0.5, 0.5).unwrap();
        let uv = inverter
            .invert(point)
            .expect("the model pass proves the four-patch corner unique");
        let on = axiolid_reference::surface::evaluate(&bump, uv.x, uv.y).unwrap();
        assert!((on - point).length() < 1e-9, "{uv:?}");
    }
}
