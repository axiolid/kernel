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
//!
//! - Analytic surfaces invert in closed form
//!   ([`axiolid_reference::surface::invert`]), which checks the residual
//!   against the linear tolerance. B-spline surfaces invert by the certified
//!   global search ([`axiolid_nurbs::invert_surface_certified`], the
//!   periodic form for a surface whose closed axes validate as periodic),
//!   which refuses an off-surface point and parameters that are not unique.
//! - A sample at a cone's apex or a sphere's pole, where every `u` names
//!   the point, is refused by name, as is a sample within the linear
//!   tolerance of it.
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

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Scalar, Tolerance, Vec3};
use axiolid_model::NodeId;
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

/// How a face's surface names a point on it.
pub(super) struct SurfaceInverter<'s> {
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
/// pass proves it unique.
#[derive(Clone, Copy)]
struct SplineSearch {
    strict: CertifiedSurfaceProjectionOptions,
    model: CertifiedSurfaceProjectionOptions,
}

impl<'s> SurfaceInverter<'s> {
    pub(super) fn new(surface: &'s Surface, tolerance: Tolerance) -> GeomResult<Self> {
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
        let outcome = match run(search.strict)? {
            Err(SurfaceInversionRefusal::OffSurface {
                distance_lower_bound,
            }) if distance_lower_bound <= self.tolerance.linear() => run(search.model)?,
            outcome => outcome,
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
    /// plus the lens between the two (the module notes give the argument).
    pub(super) fn inverted_edge(&mut self, spans: &[SpanBound], params: &[Point2]) {
        if self.unbounded.is_some() {
            return;
        }
        for (span, ends) in spans.iter().zip(params.windows(2)) {
            let Some(image) = self.chord_image(ends[0], ends[1]) else {
                self.unbounded("surface has no derivative bound over an edge chord");
                return;
            };
            let near = span.chord + image;
            let Some(lens) = self.lens(ends[0], ends[1], span.arc, near) else {
                self.unbounded(
                    "the surface is not certifiably invertible along an edge without a pcurve",
                );
                return;
            };
            self.trim = self.trim.max(near + lens);
        }
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
    /// chord's start (`B` is clamped to a B-spline's domain, which the
    /// preimage never leaves). Inside `B`, `|S(x) - S(y)| >= sigma |x - y|`,
    /// so the preimage is within `W = near / sigma` of the chord, a segment
    /// across the lens is at most `2 W` long, and the surface strays from
    /// that segment's chord by `q(2 W) / 8` at most.
    fn lens(&self, a: Point2, b: Point2, arc: Scalar, near: Scalar) -> Option<Scalar> {
        let oracle = self.oracle.as_ref()?;
        let m = (a + b) * 0.5;
        let (su, sv) = axiolid_reference::surface::partials(self.surface, m.x, m.y).ok()?;
        let (alpha, beta) = (su.length(), sv.length());
        if !(alpha > 0.0 && beta > 0.0 && alpha.is_finite() && beta.is_finite()) {
            return None;
        }
        // The scaled Gram matrix is [[1, c], [c, 1]]: eigenvalues 1 +- c.
        let c = (su.dot(sv) / (alpha * beta)).abs();
        let sigma_m = (1.0 - c).max(0.0).sqrt();
        let ((u0, u1), (v0, v1)) = parameter_domain(self.surface);
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

/// Tries at growing the lens box until the preimage provably stays in it.
const LENS_ATTEMPTS: usize = 4;

/// The native domain of a B-spline's open axes, which its preimages stay
/// in; unbounded otherwise.
fn parameter_domain(surface: &Surface) -> ((Scalar, Scalar), (Scalar, Scalar)) {
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
/// use's direction, both endpoints included, and the mesh indices of the
/// points, shared with every face that uses the edge.
pub(super) fn inverted_edge_use(
    state: &mut CurvedMeshState<'_>,
    ctx: &FaceContext<'_>,
    use_: &EdgeUse<NodeId>,
    inverter: &SurfaceInverter<'_>,
    bound: &mut FaceBound<'_>,
    previous_end: Option<Point2>,
) -> GeomResult<(Vec<Point2>, Vec<u32>)> {
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
    let (points, mut params, spans) = match state.edge_samples.get(&use_.edge) {
        Some(cached) => {
            let (indices, spans) = cached.in_direction(use_.orientation);
            let points: Vec<Vec3> = indices
                .iter()
                .map(|&index| state.mesh.positions[index as usize])
                .collect();
            let mut params = points
                .iter()
                .map(|&point| inverter.invert(point))
                .collect::<GeomResult<Vec<_>>>()?;
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
            (points, params, spans)
        }
        None => {
            let (mut points, mut params, spans) = sample_edge(&trace, inverter, bound, ctx)?;
            let mut spans = spans.ok_or(NO_CHORD_BOUND);
            if use_.orientation == Orientation::Reversed {
                points.reverse();
                params.reverse();
                if let Ok(spans) = spans.as_mut() {
                    spans.reverse();
                }
            }
            (points, params, spans)
        }
    };
    let n = points.len().checked_sub(1).ok_or_else(|| {
        GeomError::Degenerate("curved edge needs both endpoint samples".to_owned())
    })?;
    state.reserve_face_vertices(n)?;
    align(inverter, &mut params, previous_end);
    inverter.check_spline_parameters(&params, &points)?;
    match &spans {
        Ok(spans) => bound.inverted_edge(spans, &params),
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
    Ok((params, shared))
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
    edge: &axiolid_topology::Edge<NodeId>,
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
    let Some(id) = edge.curve else {
        return Ok(EdgeTrace::Straight(start, end));
    };
    let Some(axiolid_model::GeometryNode::Curve3(curve)) = ctx.graph.get(id) else {
        return Err(GeomError::InvalidInput(
            "edge curve must reference a Curve3 node".to_string(),
        ));
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

/// Sample a new edge in its own direction: uniform spans of its parameter,
/// doubled until every span's 3D midpoint is within the budget of its chord
/// and its parameters turn less than a quarter period, then while the
/// span's certified terms miss the budget, up to the segment cap.
#[allow(clippy::type_complexity)]
pub(super) fn sample_edge(
    trace: &EdgeTrace<'_>,
    inverter: &SurfaceInverter<'_>,
    bound: &FaceBound<'_>,
    ctx: &FaceContext<'_>,
) -> GeomResult<(Vec<Vec3>, Vec<Point2>, Option<Vec<SpanBound>>)> {
    let budget = ctx.chord_error;
    if !(budget.is_finite() && budget > 0.0) {
        return Err(GeomError::InvalidInput(
            "curved-edge chord tolerance must be positive and finite".to_owned(),
        ));
    }
    let mut points = vec![trace.at(0.0)?, trace.at(1.0)?];
    let mut raw = points
        .iter()
        .map(|&point| inverter.invert(point))
        .collect::<GeomResult<Vec<_>>>()?;
    // The coarsest sampling meeting the midpoint and quarter-turn tests,
    // kept as a stride into the finer samples should the certified terms
    // never fit.
    let mut basic: Option<usize> = None;
    let mut n = 1usize;
    loop {
        let mut params = raw.clone();
        unwrap(inverter, &mut params);
        let fraction = |i: usize| i as Scalar / n as Scalar;
        let mut mids = Vec::with_capacity(n);
        let mut measured: Scalar = 0.0;
        for i in 0..n {
            let mid = trace.at(0.5 * (fraction(i) + fraction(i + 1)))?;
            measured = measured.max((mid - 0.5 * (points[i] + points[i + 1])).length());
            mids.push(mid);
        }
        let quarter = turns_within_quarter(inverter, &params);
        let spans = spans_of(trace, &points);
        if measured <= budget && quarter {
            basic.get_or_insert(n);
            let certified = spans.as_ref().and_then(|spans| {
                (0..n)
                    .map(|i| {
                        let (a, b) = (params[i], params[i + 1]);
                        let near = spans[i].chord + bound.chord_image(a, b)?;
                        // A span too long to show the surface invertible
                        // along it is split like one that misses the budget.
                        let lens = bound
                            .lens(a, b, spans[i].arc, near)
                            .unwrap_or(Scalar::INFINITY);
                        Some(near + lens)
                    })
                    .try_fold(0.0, |worst: Scalar, term| term.map(|t| worst.max(t)))
            });
            // No certified term to steer by: the report names it.
            if certified.is_none_or(|worst| worst <= budget) {
                return Ok((points, params, spans));
            }
        }
        if n >= MAX_CURVED_EDGE_SEGMENTS {
            let Some(coarse) = basic else {
                if measured > budget {
                    return Err(GeomError::BudgetExceeded {
                        resource: "curved edge samples",
                    });
                }
                return Err(singular_between_samples());
            };
            let stride = n / coarse;
            let keep = |i: &usize| i % stride == 0;
            let points: Vec<Vec3> = (0..=n).filter(keep).map(|i| points[i]).collect();
            let params: Vec<Point2> = (0..=n).filter(keep).map(|i| params[i]).collect();
            let spans = spans_of(trace, &points);
            return Ok((points, params, spans));
        }
        // Halve every span, inverting only the new midpoints.
        let mut next_points = Vec::with_capacity(2 * n + 1);
        let mut next_raw = Vec::with_capacity(2 * n + 1);
        for i in 0..n {
            next_points.push(points[i]);
            next_raw.push(raw[i]);
            next_points.push(mids[i]);
            next_raw.push(inverter.invert(mids[i])?);
        }
        next_points.push(points[n]);
        next_raw.push(raw[n]);
        points = next_points;
        raw = next_raw;
        n *= 2;
    }
}

/// The bounds of the uniform spans between `points`, sampled at equal
/// fractions of the edge's parameter; `None` where a family has none.
pub(super) fn spans_of(trace: &EdgeTrace<'_>, points: &[Vec3]) -> Option<Vec<SpanBound>> {
    let n = points.len().checked_sub(1)?;
    let fraction = |i: usize| i as Scalar / n as Scalar;
    (0..n)
        .map(|i| trace.span(fraction(i), fraction(i + 1), points[i], points[i + 1]))
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

impl SurfaceInverter<'_> {
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
        let inverter = SurfaceInverter::new(&surface, Tolerance::MILLIMETRE).unwrap();
        let quarter = core::f64::consts::FRAC_PI_2;
        let turn = |d: Scalar| [Point2::new(0.0, 0.0), Point2::new(d, 0.1)];
        assert!(turns_within_quarter(&inverter, &turn(0.99 * quarter)));
        assert!(!turns_within_quarter(&inverter, &turn(1.01 * quarter)));
        assert!(!turns_within_quarter(&inverter, &turn(-1.01 * quarter)));
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
        let inverter = SurfaceInverter::new(&surface, Tolerance::MILLIMETRE).unwrap();
        for (u, v) in [(0.31, 0.77), (0.5, 0.0), (0.9, 0.42)] {
            let point = axiolid_reference::surface::evaluate(&surface, u, v).unwrap();
            let uv = inverter.invert(point).unwrap();
            let on = axiolid_reference::surface::evaluate(&surface, uv.x, uv.y).unwrap();
            assert!((on - point).length() < 1e-12, "{u} {v}: {uv:?}");
        }
    }
}
