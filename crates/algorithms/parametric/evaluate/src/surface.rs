//! Scalar reference implementation of surface evaluation (ADR 0012).
//!
//! # What this closes
//!
//! `axiolid-surface` declared six surface families and a `SurfaceEvaluator`
//! trait. Nothing implemented it, so a B-rep face on any curved surface could
//! not be tessellated, which is most faces in a real building model. This
//! reference that makes the declaration executable.
//!
//! # Parameterisation
//!
//! Each family uses the conventional parameterisation, chosen so `u` is the
//! angular direction wherever one exists (matching the curve module, where a
//! full turn is `[0, tau]`):
//!
//! | family   | `u`                  | `v`                     |
//! |----------|----------------------|-------------------------|
//! | Plane    | local x offset       | local y offset          |
//! | Cylinder | angle about z        | height along z          |
//! | Cone     | angle about z        | height along z          |
//! | Sphere   | azimuth about z      | polar, `-pi/2 .. pi/2`  |
//! | Torus    | angle about z        | angle around the tube   |
//! | BSpline  | first knot axis      | second knot axis        |
//!
//! Normals point outward for closed families (away from the axis for a
//! cylinder, away from the centre for a sphere, away from the tube centre for
//! a torus). A caller that needs the opposite convention negates; the kernel
//! does not guess.
//!
//! # What it does not do
//!
//! No surface-surface intersection, no trimming, no blending. Those are the
//! parts of a NURBS kernel this crate deliberately does not attempt: for a
//! tessellate-and-check pipeline the useful operation is evaluation, and the
//! boolean stack works on meshes.

use axiolid_contracts::BackendId;
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame3, Point3, Scalar, SpaceFrame, Tolerance, Vec3};
use axiolid_surface::{
    BSplineSurface, Cone, Cylinder, EllipticalCylinder, Plane, Sphere, Surface, Torus,
};

use crate::curve::{de_boor_recurrence, eval_homogeneous, span_in};
use crate::nurbs::SplineAxis;

/// A finite parameter rectangle for tessellation.
///
/// Elementary surfaces are infinite (a plane, a cylinder) or only periodic in
/// one direction, so a caller must say which patch it wants. Returning a
/// default would invent geometry the source never declared.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Patch {
    /// Start of the `u` interval.
    pub u_start: Scalar,
    /// End of the `u` interval.
    pub u_end: Scalar,
    /// Start of the `v` interval.
    pub v_start: Scalar,
    /// End of the `v` interval.
    pub v_end: Scalar,
}

/// Position and all first/second partial derivatives of a surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceJet {
    /// Position at `(u, v)`.
    pub point: Point3,
    /// First partial with respect to `u`.
    pub du: Vec3,
    /// First partial with respect to `v`.
    pub dv: Vec3,
    /// Second partial with respect to `u` twice.
    pub duu: Vec3,
    /// Mixed second partial.
    pub duv: Vec3,
    /// Second partial with respect to `v` twice.
    pub dvv: Vec3,
}

impl Patch {
    /// Construct a patch, rejecting an empty or non-finite rectangle.
    pub fn new(u_start: Scalar, u_end: Scalar, v_start: Scalar, v_end: Scalar) -> GeomResult<Self> {
        let all = [u_start, u_end, v_start, v_end];
        if !all.iter().all(|value| value.is_finite()) {
            return Err(GeomError::InvalidInput(format!(
                "patch bounds must be finite, got {all:?}"
            )));
        }
        if !(u_end > u_start && v_end > v_start) {
            return Err(GeomError::Degenerate(format!(
                "patch must have positive extent, got u {u_start}..{u_end}, v {v_start}..{v_end}"
            )));
        }
        Ok(Self {
            u_start,
            u_end,
            v_start,
            v_end,
        })
    }

    /// The full closed patch for a family that is periodic in `u`.
    pub fn full_turn(v_start: Scalar, v_end: Scalar) -> GeomResult<Self> {
        Self::new(0.0, core::f64::consts::TAU, v_start, v_end)
    }
}

/// Map a local-frame point into world coordinates.
fn place(frame: &Frame3, local: Vec3) -> Point3 {
    frame.origin + frame.x * local.x + frame.y * local.y + frame.z * local.z
}

/// Map a local-frame direction into world coordinates (no translation).
fn direct(frame: &Frame3, local: Vec3) -> Vec3 {
    frame.x * local.x + frame.y * local.y + frame.z * local.z
}

fn finite(value: Scalar, what: &str) -> GeomResult<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(GeomError::InvalidInput(format!(
            "{what} must be finite, got {value}"
        )))
    }
}

fn positive(value: Scalar, what: &str) -> GeomResult<()> {
    finite(value, what)?;
    if value > 0.0 {
        Ok(())
    } else {
        Err(GeomError::InvalidInput(format!(
            "{what} must be positive, got {value}"
        )))
    }
}

fn finite_surface_frame(surface: &Surface) -> GeomResult<()> {
    let frame = match surface {
        Surface::Plane(value) => Some(&value.frame),
        Surface::Cylinder(value) => Some(&value.frame),
        Surface::Cone(value) => Some(&value.frame),
        Surface::Sphere(value) => Some(&value.frame),
        Surface::Torus(value) => Some(&value.frame),
        Surface::EllipticalCylinder(value) => Some(&value.frame),
        _ => None,
    };
    if frame.is_none_or(|frame| {
        frame.origin.is_finite()
            && frame.x.is_finite()
            && frame.y.is_finite()
            && frame.z.is_finite()
    }) {
        Ok(())
    } else {
        Err(GeomError::InvalidInput(
            "surface frame must be finite".to_owned(),
        ))
    }
}

/// Position on a surface at `(u, v)`.
pub fn evaluate(surface: &Surface, u: Scalar, v: Scalar) -> GeomResult<Point3> {
    finite(u, "surface parameter u")?;
    finite(v, "surface parameter v")?;
    finite_surface_frame(surface)?;
    let point = match surface {
        Surface::Plane(p) => Ok(plane_point(p, u, v)),
        Surface::Cylinder(c) => cylinder_point(c, u, v),
        Surface::EllipticalCylinder(c) => elliptical_cylinder_point(c, u, v),
        Surface::Cone(c) => cone_point(c, u, v),
        Surface::Sphere(s) => sphere_point(s, u, v),
        Surface::Torus(t) => torus_point(t, u, v),
        Surface::BSpline(b) => bspline_point(b, u, v),
        _ => Err(GeomError::Unsupported {
            backend: ScalarSurface::ID,
            operation: axiolid_contracts::Operation::SurfaceEvaluation,
        }),
    }?;
    if point.is_finite() {
        Ok(point)
    } else {
        Err(GeomError::Degenerate(
            "surface point is non-finite".to_owned(),
        ))
    }
}

/// Analytic first partial derivatives `(∂S/∂u, ∂S/∂v)` at `(u, v)`.
///
/// Rational B-spline derivatives are evaluated in homogeneous space and
/// projected with the quotient rule. This remains stable when valid imported
/// knot domains have large offsets that make finite-difference steps vanish.
pub fn partials(surface: &Surface, u: Scalar, v: Scalar) -> GeomResult<(Vec3, Vec3)> {
    finite(u, "surface parameter u")?;
    finite(v, "surface parameter v")?;
    finite_surface_frame(surface)?;
    let value = match surface {
        Surface::Plane(p) => Ok((p.frame.x, p.frame.y)),
        Surface::Cylinder(c) => {
            positive(c.radius, "cylinder radius")?;
            let (s, co) = u.sin_cos();
            Ok((
                direct(&c.frame, Vec3::new(-c.radius * s, c.radius * co, 0.0)),
                c.frame.z,
            ))
        }
        Surface::Cone(c) => {
            finite(c.radius, "cone radius")?;
            finite(c.semi_angle, "cone semi-angle")?;
            let slope = c.semi_angle.tan();
            let radius = c.radius + v * slope;
            if radius < 0.0 {
                return Err(GeomError::Degenerate(format!(
                    "cone radius is negative at v = {v}: the patch crosses the apex"
                )));
            }
            let (s, co) = u.sin_cos();
            Ok((
                direct(&c.frame, Vec3::new(-radius * s, radius * co, 0.0)),
                direct(&c.frame, Vec3::new(slope * co, slope * s, 1.0)),
            ))
        }
        Surface::Sphere(sphere) => {
            positive(sphere.radius, "sphere radius")?;
            let (su, cu) = u.sin_cos();
            let (sv, cv) = v.sin_cos();
            Ok((
                direct(
                    &sphere.frame,
                    Vec3::new(-sphere.radius * cv * su, sphere.radius * cv * cu, 0.0),
                ),
                direct(
                    &sphere.frame,
                    Vec3::new(
                        -sphere.radius * sv * cu,
                        -sphere.radius * sv * su,
                        sphere.radius * cv,
                    ),
                ),
            ))
        }
        Surface::Torus(torus) => {
            positive(torus.major_radius, "torus major radius")?;
            positive(torus.minor_radius, "torus minor radius")?;
            let (su, cu) = u.sin_cos();
            let (sv, cv) = v.sin_cos();
            let ring = torus.major_radius + torus.minor_radius * cv;
            Ok((
                direct(&torus.frame, Vec3::new(-ring * su, ring * cu, 0.0)),
                direct(
                    &torus.frame,
                    Vec3::new(
                        -torus.minor_radius * sv * cu,
                        -torus.minor_radius * sv * su,
                        torus.minor_radius * cv,
                    ),
                ),
            ))
        }
        Surface::EllipticalCylinder(c) => {
            positive(c.semi_axis_x, "elliptical cylinder semi-axis x")?;
            positive(c.semi_axis_y, "elliptical cylinder semi-axis y")?;
            let (su, cu) = u.sin_cos();
            // The u-partial is NOT radial: its components carry different
            // semi-axes, which is exactly why the normal of an elliptical
            // cylinder is not its radial direction.
            Ok((
                direct(
                    &c.frame,
                    Vec3::new(-c.semi_axis_x * su, c.semi_axis_y * cu, 0.0),
                ),
                direct(&c.frame, Vec3::Z),
            ))
        }
        Surface::BSpline(b) => bspline_partials(b, u, v),
        _ => Err(GeomError::Unsupported {
            backend: ScalarSurface::ID,
            operation: axiolid_contracts::Operation::SurfaceEvaluation,
        }),
    }?;
    if value.0.is_finite() && value.1.is_finite() {
        Ok(value)
    } else {
        Err(GeomError::Degenerate(
            "surface partial is non-finite".to_owned(),
        ))
    }
}

/// Second-order differential jet at `(u, v)`.
///
/// All values are analytic in the surface's native parameterisation. Rational
/// B-splines are differentiated in homogeneous space before projection.
pub fn jet(surface: &Surface, u: Scalar, v: Scalar) -> GeomResult<SurfaceJet> {
    finite(u, "surface parameter u")?;
    finite(v, "surface parameter v")?;
    finite_surface_frame(surface)?;
    let value = match surface {
        Surface::Plane(p) => SurfaceJet {
            point: plane_point(p, u, v),
            du: p.frame.x,
            dv: p.frame.y,
            duu: Vec3::ZERO,
            duv: Vec3::ZERO,
            dvv: Vec3::ZERO,
        },
        Surface::Cylinder(c) => {
            positive(c.radius, "cylinder radius")?;
            let (s, co) = u.sin_cos();
            SurfaceJet {
                point: cylinder_point(c, u, v)?,
                du: direct(&c.frame, Vec3::new(-c.radius * s, c.radius * co, 0.0)),
                dv: c.frame.z,
                duu: direct(&c.frame, Vec3::new(-c.radius * co, -c.radius * s, 0.0)),
                duv: Vec3::ZERO,
                dvv: Vec3::ZERO,
            }
        }
        Surface::Cone(c) => {
            finite(c.radius, "cone radius")?;
            finite(c.semi_angle, "cone semi-angle")?;
            let slope = c.semi_angle.tan();
            let radius = c.radius + v * slope;
            if radius < 0.0 {
                return Err(GeomError::Degenerate(format!(
                    "cone radius is negative at v = {v}: the patch crosses the apex"
                )));
            }
            let (s, co) = u.sin_cos();
            SurfaceJet {
                point: cone_point(c, u, v)?,
                du: direct(&c.frame, Vec3::new(-radius * s, radius * co, 0.0)),
                dv: direct(&c.frame, Vec3::new(slope * co, slope * s, 1.0)),
                duu: direct(&c.frame, Vec3::new(-radius * co, -radius * s, 0.0)),
                duv: direct(&c.frame, Vec3::new(-slope * s, slope * co, 0.0)),
                dvv: Vec3::ZERO,
            }
        }
        Surface::Sphere(sphere) => {
            positive(sphere.radius, "sphere radius")?;
            let r = sphere.radius;
            let (su, cu) = u.sin_cos();
            let (sv, cv) = v.sin_cos();
            SurfaceJet {
                point: sphere_point(sphere, u, v)?,
                du: direct(&sphere.frame, Vec3::new(-r * cv * su, r * cv * cu, 0.0)),
                dv: direct(&sphere.frame, Vec3::new(-r * sv * cu, -r * sv * su, r * cv)),
                duu: direct(&sphere.frame, Vec3::new(-r * cv * cu, -r * cv * su, 0.0)),
                duv: direct(&sphere.frame, Vec3::new(r * sv * su, -r * sv * cu, 0.0)),
                dvv: direct(
                    &sphere.frame,
                    Vec3::new(-r * cv * cu, -r * cv * su, -r * sv),
                ),
            }
        }
        Surface::Torus(torus) => {
            positive(torus.major_radius, "torus major radius")?;
            positive(torus.minor_radius, "torus minor radius")?;
            let r = torus.minor_radius;
            let (su, cu) = u.sin_cos();
            let (sv, cv) = v.sin_cos();
            let ring = torus.major_radius + r * cv;
            SurfaceJet {
                point: torus_point(torus, u, v)?,
                du: direct(&torus.frame, Vec3::new(-ring * su, ring * cu, 0.0)),
                dv: direct(&torus.frame, Vec3::new(-r * sv * cu, -r * sv * su, r * cv)),
                duu: direct(&torus.frame, Vec3::new(-ring * cu, -ring * su, 0.0)),
                duv: direct(&torus.frame, Vec3::new(r * sv * su, -r * sv * cu, 0.0)),
                dvv: direct(&torus.frame, Vec3::new(-r * cv * cu, -r * cv * su, -r * sv)),
            }
        }
        Surface::BSpline(b) => bspline_jet(b, u, v)?,
        _ => {
            return Err(GeomError::Unsupported {
                backend: ScalarSurface::ID,
                operation: axiolid_contracts::Operation::SurfaceEvaluation,
            })
        }
    };
    if [
        value.point,
        value.du,
        value.dv,
        value.duu,
        value.duv,
        value.dvv,
    ]
    .iter()
    .all(|vector| vector.is_finite())
    {
        Ok(value)
    } else {
        Err(GeomError::Degenerate(
            "surface differential jet is non-finite".to_owned(),
        ))
    }
}

/// Unit normal at `(u, v)`.
///
/// Computed from the exact analytic partial derivatives rather than by
/// differencing evaluated points: a finite difference loses precision exactly
/// where it matters most, at high curvature.
pub fn normal(surface: &Surface, u: Scalar, v: Scalar) -> GeomResult<Vec3> {
    finite(u, "surface parameter u")?;
    finite(v, "surface parameter v")?;
    finite_surface_frame(surface)?;
    let n = match surface {
        Surface::Plane(p) => p.frame.z,
        Surface::Cylinder(c) => {
            positive(c.radius, "cylinder radius")?;
            let (s, co) = u.sin_cos();
            direct(&c.frame, Vec3::new(co, s, 0.0))
        }
        Surface::Cone(c) => cone_normal(c, u)?,
        Surface::EllipticalCylinder(c) => {
            positive(c.semi_axis_x, "elliptical cylinder semi-axis x")?;
            positive(c.semi_axis_y, "elliptical cylinder semi-axis y")?;
            // NOT the radial direction. For a circular cylinder the two
            // coincide, but for a 3:1 ellipse they differ by up to 53
            // degrees, agreeing only at the four axis points. Taking the
            // cross product of the partials is the definition and is right
            // everywhere.
            let (su, cu) = u.sin_cos();
            let along_u = direct(
                &c.frame,
                Vec3::new(-c.semi_axis_x * su, c.semi_axis_y * cu, 0.0),
            );
            along_u.cross(c.frame.z)
        }
        Surface::Sphere(s) => {
            positive(s.radius, "sphere radius")?;
            let (su, cu) = u.sin_cos();
            let (sv, cv) = v.sin_cos();
            direct(&s.frame, Vec3::new(cv * cu, cv * su, sv))
        }
        Surface::Torus(t) => {
            positive(t.minor_radius, "torus minor radius")?;
            let (su, cu) = u.sin_cos();
            let (sv, cv) = v.sin_cos();
            direct(&t.frame, Vec3::new(cv * cu, cv * su, sv))
        }
        Surface::BSpline(b) => bspline_normal(b, u, v)?,
        _ => {
            return Err(GeomError::Unsupported {
                backend: ScalarSurface::ID,
                operation: axiolid_contracts::Operation::SurfaceEvaluation,
            })
        }
    };
    let length = n.length();
    if !(length > 0.0 && n.is_finite()) {
        return Err(GeomError::Degenerate(format!(
            "surface normal is not orientable at ({u}, {v})"
        )));
    }
    Ok(n / length)
}

fn plane_point(p: &Plane, u: Scalar, v: Scalar) -> Point3 {
    place(&p.frame, Vec3::new(u, v, 0.0))
}

fn cylinder_point(c: &Cylinder, u: Scalar, v: Scalar) -> GeomResult<Point3> {
    positive(c.radius, "cylinder radius")?;
    let (s, co) = u.sin_cos();
    Ok(place(&c.frame, Vec3::new(c.radius * co, c.radius * s, v)))
}

fn elliptical_cylinder_point(c: &EllipticalCylinder, u: Scalar, v: Scalar) -> GeomResult<Point3> {
    positive(c.semi_axis_x, "elliptical cylinder semi-axis x")?;
    positive(c.semi_axis_y, "elliptical cylinder semi-axis y")?;
    let (s, co) = u.sin_cos();
    Ok(place(
        &c.frame,
        Vec3::new(c.semi_axis_x * co, c.semi_axis_y * s, v),
    ))
}

fn cone_point(c: &Cone, u: Scalar, v: Scalar) -> GeomResult<Point3> {
    finite(c.radius, "cone radius")?;
    finite(c.semi_angle, "cone semi-angle")?;
    // Radius shrinks with height at the semi-angle; a negative radius means
    // the surface has passed through the apex, which is not a valid patch.
    let r = c.radius + v * c.semi_angle.tan();
    if r < 0.0 {
        return Err(GeomError::Degenerate(format!(
            "cone radius is negative at v = {v}: the patch crosses the apex"
        )));
    }
    let (s, co) = u.sin_cos();
    Ok(place(&c.frame, Vec3::new(r * co, r * s, v)))
}

fn cone_normal(c: &Cone, u: Scalar) -> GeomResult<Vec3> {
    finite(c.semi_angle, "cone semi-angle")?;
    let (s, co) = u.sin_cos();
    // Outward radial component, tilted by the semi-angle: the normal leans
    // toward the axis as the cone narrows.
    let (sa, ca) = c.semi_angle.sin_cos();
    Ok(direct(&c.frame, Vec3::new(ca * co, ca * s, -sa)))
}

fn sphere_point(s: &Sphere, u: Scalar, v: Scalar) -> GeomResult<Point3> {
    positive(s.radius, "sphere radius")?;
    let (su, cu) = u.sin_cos();
    let (sv, cv) = v.sin_cos();
    Ok(place(
        &s.frame,
        Vec3::new(s.radius * cv * cu, s.radius * cv * su, s.radius * sv),
    ))
}

fn torus_point(t: &Torus, u: Scalar, v: Scalar) -> GeomResult<Point3> {
    positive(t.major_radius, "torus major radius")?;
    positive(t.minor_radius, "torus minor radius")?;
    let (su, cu) = u.sin_cos();
    let (sv, cv) = v.sin_cos();
    let ring = t.major_radius + t.minor_radius * cv;
    Ok(place(
        &t.frame,
        Vec3::new(ring * cu, ring * su, t.minor_radius * sv),
    ))
}

// --- tensor-product B-spline ------------------------------------------------

type Axis = SplineAxis;

/// Validate the control net and both axes together.
fn bspline_axes(b: &BSplineSurface) -> GeomResult<(Axis, Axis)> {
    let rows = b.control_points.len();
    if rows == 0 {
        return Err(GeomError::InvalidInput(
            "B-spline surface has no control points".to_owned(),
        ));
    }
    let cols = b.control_points[0].len();
    if cols == 0 {
        return Err(GeomError::InvalidInput(
            "B-spline surface control net has an empty row".to_owned(),
        ));
    }
    // A ragged net is a data error, not something to paper over: evaluating it
    // would silently read a different surface than the source declared.
    if b.control_points.iter().any(|row| row.len() != cols) {
        return Err(GeomError::InvalidInput(
            "B-spline surface control net is ragged".to_owned(),
        ));
    }
    if b.control_points
        .iter()
        .flatten()
        .any(|point| !point.is_finite())
    {
        return Err(GeomError::InvalidInput(
            "B-spline surface control points must be finite".to_owned(),
        ));
    }
    if let Some(w) = &b.weights {
        if w.len() != rows || w.iter().any(|row| row.len() != cols) {
            return Err(GeomError::InvalidInput(
                "B-spline surface weight net does not match the control net".to_owned(),
            ));
        }
        if w.iter()
            .flatten()
            .any(|weight| !weight.is_finite() || *weight <= 0.0)
        {
            return Err(GeomError::InvalidInput(
                "B-spline surface weights must be finite and strictly positive".to_owned(),
            ));
        }
    }
    let u = Axis::new(&b.u_knots, &b.u_multiplicities, b.u_degree, rows, "u")?;
    let v = Axis::new(&b.v_knots, &b.v_multiplicities, b.v_degree, cols, "v")?;
    Ok((u, v))
}

/// Tensor-product de Boor: evaluate along `v` per influencing row, then along
/// `u` through those results.
///
/// Rational surfaces interpolate in homogeneous space throughout; projecting
/// per row and averaging afterwards is the classic wrong answer.
fn bspline_point(b: &BSplineSurface, u: Scalar, v: Scalar) -> GeomResult<Point3> {
    let (ua, va) = bspline_axes(b)?;
    let (uc, vc) = (ua.clamp(u), va.clamp(v));
    let uspan = span_in(&ua.knots, ua.count, ua.degree, uc);
    let vspan = span_in(&va.knots, va.count, va.degree, vc);

    // Stage one: collapse each influencing row along v, staying homogeneous.
    let mut row_points: Vec<[Scalar; 3]> = Vec::with_capacity(ua.degree + 1);
    let mut row_weights: Vec<Scalar> = Vec::with_capacity(ua.degree + 1);
    for i in 0..=ua.degree {
        let row = uspan - ua.degree + i;
        let mut pts: Vec<[Scalar; 3]> = Vec::with_capacity(va.degree + 1);
        let mut wts: Vec<Scalar> = Vec::with_capacity(va.degree + 1);
        for j in 0..=va.degree {
            let col = vspan - va.degree + j;
            let w = b.weights.as_ref().map_or(1.0, |ws| ws[row][col]);
            let p = b.control_points[row][col];
            let homogeneous = [p.x * w, p.y * w, p.z * w];
            if homogeneous.iter().any(|value| !value.is_finite()) {
                return Err(GeomError::Degenerate(
                    "B-spline surface homogeneous control point overflowed".to_owned(),
                ));
            }
            pts.push(homogeneous);
            wts.push(w);
        }
        de_boor_recurrence(&va.knots, vspan, va.degree, vc, &mut pts, &mut wts);
        row_points.push(pts[va.degree]);
        row_weights.push(wts[va.degree]);
    }

    // Stage two: collapse the row results along u.
    de_boor_recurrence(
        &ua.knots,
        uspan,
        ua.degree,
        uc,
        &mut row_points,
        &mut row_weights,
    );

    let w = row_weights[ua.degree];
    if !w.is_finite() || w == 0.0 {
        return Err(GeomError::Degenerate(
            "B-spline surface weight collapsed to zero".to_owned(),
        ));
    }
    let p = row_points[ua.degree];
    Ok(Point3::new(p[0] / w, p[1] / w, p[2] / w))
}

/// Borrowed knot axes for one homogeneous tensor-product evaluation.
#[derive(Clone, Copy)]
struct HomogeneousAxes<'a> {
    u_knots: &'a [Scalar],
    u_degree: usize,
    v_knots: &'a [Scalar],
    v_degree: usize,
}

/// Evaluate one homogeneous tensor-product control net without projecting.
fn eval_tensor_homogeneous(
    axes: HomogeneousAxes<'_>,
    points: &[Vec<[Scalar; 3]>],
    weights: &[Vec<Scalar>],
    u: Scalar,
    v: Scalar,
) -> ([Scalar; 3], Scalar) {
    let mut row_points = Vec::with_capacity(points.len());
    let mut row_weights = Vec::with_capacity(points.len());
    for (row_points_h, row_weights_h) in points.iter().zip(weights) {
        let (point, weight) =
            eval_homogeneous(axes.v_knots, axes.v_degree, row_points_h, row_weights_h, v);
        row_points.push(point);
        row_weights.push(weight);
    }
    eval_homogeneous(axes.u_knots, axes.u_degree, &row_points, &row_weights, u)
}

type HomogeneousPointNet = Vec<Vec<[Scalar; 3]>>;
type HomogeneousWeightNet = Vec<Vec<Scalar>>;

/// Homogeneous point and weight control nets for a validated surface.
fn homogeneous_control_net(
    b: &BSplineSurface,
) -> GeomResult<(HomogeneousPointNet, HomogeneousWeightNet)> {
    let mut points = Vec::with_capacity(b.control_points.len());
    let mut weights = Vec::with_capacity(b.control_points.len());
    for (i, row) in b.control_points.iter().enumerate() {
        let mut point_row = Vec::with_capacity(row.len());
        let mut weight_row = Vec::with_capacity(row.len());
        for (j, point) in row.iter().enumerate() {
            let weight = b.weights.as_ref().map_or(1.0, |net| net[i][j]);
            let homogeneous = [point.x * weight, point.y * weight, point.z * weight];
            if homogeneous.iter().any(|value| !value.is_finite()) {
                return Err(GeomError::Degenerate(
                    "B-spline surface homogeneous control point overflowed".to_owned(),
                ));
            }
            point_row.push(homogeneous);
            weight_row.push(weight);
        }
        points.push(point_row);
        weights.push(weight_row);
    }
    Ok((points, weights))
}

/// Differentiate a homogeneous control net along `u`.
fn derivative_net_u(
    points: &[Vec<[Scalar; 3]>],
    weights: &[Vec<Scalar>],
    knots: &[Scalar],
    degree: usize,
) -> (Vec<Vec<[Scalar; 3]>>, Vec<Vec<Scalar>>) {
    let rows = points.len() - 1;
    let cols = points[0].len();
    let mut derivative_points = Vec::with_capacity(rows);
    let mut derivative_weights = Vec::with_capacity(rows);
    for i in 0..rows {
        let denominator = knots[i + degree + 1] - knots[i + 1];
        let factor = if denominator.abs() > 0.0 {
            degree as Scalar / denominator
        } else {
            0.0
        };
        let mut point_row = Vec::with_capacity(cols);
        let mut weight_row = Vec::with_capacity(cols);
        for j in 0..cols {
            point_row.push(core::array::from_fn(|k| {
                factor * (points[i + 1][j][k] - points[i][j][k])
            }));
            weight_row.push(factor * (weights[i + 1][j] - weights[i][j]));
        }
        derivative_points.push(point_row);
        derivative_weights.push(weight_row);
    }
    (derivative_points, derivative_weights)
}

/// Differentiate a homogeneous control net along `v`.
fn derivative_net_v(
    points: &[Vec<[Scalar; 3]>],
    weights: &[Vec<Scalar>],
    knots: &[Scalar],
    degree: usize,
) -> (Vec<Vec<[Scalar; 3]>>, Vec<Vec<Scalar>>) {
    let rows = points.len();
    let cols = points[0].len() - 1;
    let mut derivative_points = Vec::with_capacity(rows);
    let mut derivative_weights = Vec::with_capacity(rows);
    for i in 0..rows {
        let mut point_row = Vec::with_capacity(cols);
        let mut weight_row = Vec::with_capacity(cols);
        for j in 0..cols {
            let denominator = knots[j + degree + 1] - knots[j + 1];
            let factor = if denominator.abs() > 0.0 {
                degree as Scalar / denominator
            } else {
                0.0
            };
            point_row.push(core::array::from_fn(|k| {
                factor * (points[i][j + 1][k] - points[i][j][k])
            }));
            weight_row.push(factor * (weights[i][j + 1] - weights[i][j]));
        }
        derivative_points.push(point_row);
        derivative_weights.push(weight_row);
    }
    (derivative_points, derivative_weights)
}

/// Project one homogeneous derivative with the rational quotient rule.
fn project_derivative(
    point: [Scalar; 3],
    weight: Scalar,
    derivative: [Scalar; 3],
    derivative_weight: Scalar,
    axis: &str,
) -> GeomResult<Vec3> {
    if !weight.is_finite() || weight == 0.0 {
        return Err(GeomError::Degenerate(
            "B-spline surface weight collapsed to a non-finite or zero value".to_owned(),
        ));
    }
    let value = Vec3::new(
        (derivative[0] - point[0] * derivative_weight / weight) / weight,
        (derivative[1] - point[1] * derivative_weight / weight) / weight,
        (derivative[2] - point[2] * derivative_weight / weight) / weight,
    );
    if !value.is_finite() {
        return Err(GeomError::Degenerate(format!(
            "B-spline surface {axis} derivative is non-finite"
        )));
    }
    Ok(value)
}

/// Exact first partials of a rational tensor-product B-spline.
fn bspline_partials(b: &BSplineSurface, u: Scalar, v: Scalar) -> GeomResult<(Vec3, Vec3)> {
    let (ua, va) = bspline_axes(b)?;
    let (uc, vc) = (ua.clamp(u), va.clamp(v));
    let (points, weights) = homogeneous_control_net(b)?;
    let (point, weight) = eval_tensor_homogeneous(
        HomogeneousAxes {
            u_knots: &ua.knots,
            u_degree: ua.degree,
            v_knots: &va.knots,
            v_degree: va.degree,
        },
        &points,
        &weights,
        uc,
        vc,
    );

    let (u_points, u_weights) = derivative_net_u(&points, &weights, &ua.knots, ua.degree);
    let (du, du_weight) = eval_tensor_homogeneous(
        HomogeneousAxes {
            u_knots: &ua.knots[1..ua.knots.len() - 1],
            u_degree: ua.degree - 1,
            v_knots: &va.knots,
            v_degree: va.degree,
        },
        &u_points,
        &u_weights,
        uc,
        vc,
    );

    let (v_points, v_weights) = derivative_net_v(&points, &weights, &va.knots, va.degree);
    let (dv, dv_weight) = eval_tensor_homogeneous(
        HomogeneousAxes {
            u_knots: &ua.knots,
            u_degree: ua.degree,
            v_knots: &va.knots[1..va.knots.len() - 1],
            v_degree: va.degree - 1,
        },
        &v_points,
        &v_weights,
        uc,
        vc,
    );

    Ok((
        project_derivative(point, weight, du, du_weight, "u")?,
        project_derivative(point, weight, dv, dv_weight, "v")?,
    ))
}

/// Full second-order jet of a rational tensor-product B-spline.
/// Second-order differential jet of a B-spline surface without enum wrapping.
pub fn bspline_jet(b: &BSplineSurface, u: Scalar, v: Scalar) -> GeomResult<SurfaceJet> {
    let (ua, va) = bspline_axes(b)?;
    let (uc, vc) = (ua.clamp(u), va.clamp(v));
    let (points, weights) = homogeneous_control_net(b)?;
    let base_axes = HomogeneousAxes {
        u_knots: &ua.knots,
        u_degree: ua.degree,
        v_knots: &va.knots,
        v_degree: va.degree,
    };
    let (point, weight) = eval_tensor_homogeneous(base_axes, &points, &weights, uc, vc);
    if !weight.is_finite() || weight == 0.0 {
        return Err(GeomError::Degenerate(
            "B-spline surface weight collapsed to a non-finite or zero value".to_owned(),
        ));
    }
    let position = Point3::new(point[0] / weight, point[1] / weight, point[2] / weight);

    let (u_points, u_weights) = derivative_net_u(&points, &weights, &ua.knots, ua.degree);
    let u_knots = &ua.knots[1..ua.knots.len() - 1];
    let (du_h, du_weight) = eval_tensor_homogeneous(
        HomogeneousAxes {
            u_knots,
            u_degree: ua.degree - 1,
            v_knots: &va.knots,
            v_degree: va.degree,
        },
        &u_points,
        &u_weights,
        uc,
        vc,
    );
    let du = project_derivative(point, weight, du_h, du_weight, "u")?;

    let (v_points, v_weights) = derivative_net_v(&points, &weights, &va.knots, va.degree);
    let v_knots = &va.knots[1..va.knots.len() - 1];
    let (dv_h, dv_weight) = eval_tensor_homogeneous(
        HomogeneousAxes {
            u_knots: &ua.knots,
            u_degree: ua.degree,
            v_knots,
            v_degree: va.degree - 1,
        },
        &v_points,
        &v_weights,
        uc,
        vc,
    );
    let dv = project_derivative(point, weight, dv_h, dv_weight, "v")?;

    let (duu_h, duu_weight) = if ua.degree >= 2 {
        let (net, net_weights) = derivative_net_u(&u_points, &u_weights, u_knots, ua.degree - 1);
        eval_tensor_homogeneous(
            HomogeneousAxes {
                u_knots: &u_knots[1..u_knots.len() - 1],
                u_degree: ua.degree - 2,
                v_knots: &va.knots,
                v_degree: va.degree,
            },
            &net,
            &net_weights,
            uc,
            vc,
        )
    } else {
        ([0.0; 3], 0.0)
    };
    let duu = project_second(point, weight, du, duu_h, du_weight, duu_weight, "uu")?;

    let (dvv_h, dvv_weight) = if va.degree >= 2 {
        let (net, net_weights) = derivative_net_v(&v_points, &v_weights, v_knots, va.degree - 1);
        eval_tensor_homogeneous(
            HomogeneousAxes {
                u_knots: &ua.knots,
                u_degree: ua.degree,
                v_knots: &v_knots[1..v_knots.len() - 1],
                v_degree: va.degree - 2,
            },
            &net,
            &net_weights,
            uc,
            vc,
        )
    } else {
        ([0.0; 3], 0.0)
    };
    let dvv = project_second(point, weight, dv, dvv_h, dv_weight, dvv_weight, "vv")?;

    let (uv_points, uv_weights) = derivative_net_v(&u_points, &u_weights, &va.knots, va.degree);
    let (duv_h, duv_weight) = eval_tensor_homogeneous(
        HomogeneousAxes {
            u_knots,
            u_degree: ua.degree - 1,
            v_knots,
            v_degree: va.degree - 1,
        },
        &uv_points,
        &uv_weights,
        uc,
        vc,
    );
    let duv = project_mixed(
        point, weight, du, du_weight, dv, dv_weight, duv_h, duv_weight,
    )?;

    Ok(SurfaceJet {
        point: position,
        du,
        dv,
        duu,
        duv,
        dvv,
    })
}

fn project_second(
    point: [Scalar; 3],
    weight: Scalar,
    first: Vec3,
    second: [Scalar; 3],
    first_weight: Scalar,
    second_weight: Scalar,
    axis: &str,
) -> GeomResult<Vec3> {
    let position = Vec3::new(point[0], point[1], point[2]) / weight;
    let value = (Vec3::new(second[0], second[1], second[2])
        - 2.0 * first_weight * first
        - second_weight * position)
        / weight;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(GeomError::Degenerate(format!(
            "B-spline surface {axis} second derivative is non-finite"
        )))
    }
}

#[allow(clippy::too_many_arguments)]
fn project_mixed(
    point: [Scalar; 3],
    weight: Scalar,
    du: Vec3,
    du_weight: Scalar,
    dv: Vec3,
    dv_weight: Scalar,
    mixed: [Scalar; 3],
    mixed_weight: Scalar,
) -> GeomResult<Vec3> {
    let position = Vec3::new(point[0], point[1], point[2]) / weight;
    let value = (Vec3::new(mixed[0], mixed[1], mixed[2])
        - du_weight * dv
        - dv_weight * du
        - mixed_weight * position)
        / weight;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(GeomError::Degenerate(
            "B-spline surface uv mixed derivative is non-finite".to_owned(),
        ))
    }
}

/// Normal from analytic rational tensor-product partial derivatives.
fn bspline_normal(b: &BSplineSurface, u: Scalar, v: Scalar) -> GeomResult<Vec3> {
    let (du, dv) = bspline_partials(b, u, v)?;
    Ok(du.cross(dv))
}

/// The [`axiolid_surface::SurfaceEvaluator`] implementation, so a caller can dispatch through
/// the trait rather than the free functions.
#[derive(Debug, Default, Clone, Copy)]
pub struct ScalarSurface;

impl ScalarSurface {
    /// Identity reported in structured errors, matching `ScalarBoolean`.
    pub const ID: BackendId = BackendId::new("scalar-reference");
}

impl axiolid_surface::SurfaceEvaluator<Surface> for ScalarSurface {
    type Error = GeomError;

    fn evaluate(
        &self,
        surface: &Surface,
        u: Scalar,
        v: Scalar,
        _tolerance: axiolid_core::Tolerance,
    ) -> Result<Point3, Self::Error> {
        evaluate(surface, u, v)
    }

    fn normal(
        &self,
        surface: &Surface,
        u: Scalar,
        v: Scalar,
        _tolerance: axiolid_core::Tolerance,
    ) -> Result<Vec3, Self::Error> {
        normal(surface, u, v)
    }
}

/// Surface parameters `(u, v)` whose evaluation reproduces `point`.
///
/// This is the exact inverse of [`evaluate`] for the analytic surfaces,
/// derived from each parameterisation rather than found by iteration, so
/// it neither needs a seed nor converges to a nearby-but-wrong branch.
///
/// The point must already lie ON the surface: this answers "which
/// parameters name this point", not "which point is nearest". A sweep
/// directrix that has drifted off its reference surface is a modelling
/// error, and silently projecting it would tilt every section frame by an
/// amount nothing downstream can detect. The residual is therefore checked
/// against `tolerance` and a miss is reported rather than absorbed.
///
/// Parameters that no unique answer exists for are refused, not guessed:
/// at a cone apex or a sphere pole the whole `u` circle maps to one point,
/// so any choice would be arbitrary and would rotate the swept section.
pub fn invert(
    surface: &Surface,
    point: Point3,
    tolerance: axiolid_core::Tolerance,
) -> GeomResult<(Scalar, Scalar)> {
    let (u, v) = match surface {
        Surface::Plane(p) => {
            let local = to_local(&p.frame, point, tolerance)?;
            (local.x, local.y)
        }
        Surface::Cylinder(c) => {
            positive(c.radius, "cylinder radius")?;
            let local = to_local(&c.frame, point, tolerance)?;
            (angle_about_axis(local, "cylinder")?, local.z)
        }
        Surface::Cone(c) => {
            finite(c.radius, "cone radius")?;
            finite(c.semi_angle, "cone semi-angle")?;
            let local = to_local(&c.frame, point, tolerance)?;
            // At the apex the radius vanishes and every u names the same
            // point, so the angle is unrecoverable rather than merely
            // imprecise.
            (angle_about_axis(local, "cone")?, local.z)
        }
        Surface::Sphere(s) => {
            positive(s.radius, "sphere radius")?;
            let local = to_local(&s.frame, point, tolerance)?;
            // Latitude first: it is well defined even at the poles, which
            // the angle lookup then rejects.
            let sin_v = (local.z / s.radius).clamp(-1.0, 1.0);
            (angle_about_axis(local, "sphere")?, sin_v.asin())
        }
        Surface::Torus(t) => {
            positive(t.major_radius, "torus major radius")?;
            positive(t.minor_radius, "torus minor radius")?;
            let local = to_local(&t.frame, point, tolerance)?;
            let ring = (local.x * local.x + local.y * local.y).sqrt();
            (
                angle_about_axis(local, "torus")?,
                (local.z).atan2(ring - t.major_radius),
            )
        }
        // A B-spline has no closed-form inverse; recovering parameters
        // needs iterative closest-point with its own seeding and
        // convergence contract. `Surface` is non-exhaustive, so any
        // future variant lands here too and is refused by name rather
        // than silently taking an analytic branch that does not fit it.
        _ => {
            return Err(GeomError::Unsupported {
                backend: ScalarSurface::ID,
                operation: axiolid_contracts::Operation::SurfaceEvaluation,
            });
        }
    };
    // The parameters are only meaningful if they reproduce the point.
    // This is what turns a silent mis-parameterisation into an error.
    let round_trip = evaluate(surface, u, v)?;
    let residual = (round_trip - point).length();
    if residual > tolerance.linear() {
        return Err(GeomError::Degenerate(format!(
            "point is {residual} from the surface, beyond the {} tolerance: \
             inversion names a point ON the surface and does not project",
            tolerance.linear()
        )));
    }
    Ok((u, v))
}

/// Parameters of the closest point on a surface to an arbitrary point.
///
/// This is the counterpart to [`invert`], and the distinction matters.
/// `invert` names a point that is ALREADY on the surface and refuses one
/// that is not. Projection accepts a point anywhere and answers where the
/// surface is nearest to it.
///
/// Refinement needs projection, not inversion: the midpoint of a chord
/// across a faceted cylinder lies strictly inside the cylinder, so
/// inversion correctly refuses it while projection is exactly the question
/// being asked.
///
/// Every arm here is a closed form, so the result is exact rather than the
/// stopping point of an iteration. Where the closest point is genuinely
/// ambiguous -- a point on a cylinder's axis is equidistant from every
/// point of the surface -- this refuses by name instead of returning one
/// arbitrary member of the tie.
///
/// # Errors
///
/// Refuses a non-orthonormal frame, a non-finite point, a degenerate
/// radius, an ambiguous (equidistant) configuration, and any surface with
/// no closed-form projection.
pub fn project(
    surface: &Surface,
    point: Point3,
    tolerance: axiolid_core::Tolerance,
) -> GeomResult<(Scalar, Scalar)> {
    let ambiguous = |what: &str| {
        GeomError::Degenerate(format!(
            "{what}: the closest point is not unique, so no projection names it"
        ))
    };
    let (u, v) = match surface {
        // A plane's closest point is the orthogonal foot, which the local
        // frame already gives directly.
        Surface::Plane(p) => {
            let local = to_local(&p.frame, point, tolerance)?;
            (local.x, local.y)
        }
        // Radial projection: slide along the axis-perpendicular direction
        // to the radius. Undefined exactly on the axis.
        Surface::Cylinder(c) => {
            positive(c.radius, "cylinder radius")?;
            let local = to_local(&c.frame, point, tolerance)?;
            let ring = (local.x * local.x + local.y * local.y).sqrt();
            if ring <= tolerance.linear() {
                return Err(ambiguous("point lies on the cylinder axis"));
            }
            (local.y.atan2(local.x), local.z)
        }
        // Radial projection from the centre. Undefined exactly at it.
        Surface::Sphere(s) => {
            positive(s.radius, "sphere radius")?;
            let local = to_local(&s.frame, point, tolerance)?;
            let distance = (local.x * local.x + local.y * local.y + local.z * local.z).sqrt();
            if distance <= tolerance.linear() {
                return Err(ambiguous("point lies at the sphere centre"));
            }
            let ring = (local.x * local.x + local.y * local.y).sqrt();
            if ring <= tolerance.linear() {
                return Err(ambiguous("point lies on the sphere's polar axis"));
            }
            let sin_v = (local.z / distance).clamp(-1.0, 1.0);
            (local.y.atan2(local.x), sin_v.asin())
        }
        // The nearest point on a cone is along the SLANT, not the radius:
        // the generator is a line in the (rho, z) half-plane, so project
        // onto that line rather than onto a circle of constant z.
        Surface::Cone(c) => {
            finite(c.radius, "cone radius")?;
            finite(c.semi_angle, "cone semi-angle")?;
            let local = to_local(&c.frame, point, tolerance)?;
            let ring = (local.x * local.x + local.y * local.y).sqrt();
            if ring <= tolerance.linear() {
                return Err(ambiguous("point lies on the cone axis"));
            }
            let slope = c.semi_angle.tan();
            if !slope.is_finite() {
                return Err(GeomError::Degenerate(
                    "cone semi-angle is a right angle: the surface degenerates to a plane".into(),
                ));
            }
            // Generator through (radius, 0) with direction (slope, 1),
            // normalised so the dot product is a true arc position.
            let length = (slope * slope + 1.0).sqrt();
            let (dr, dz) = (slope / length, 1.0 / length);
            let step = (ring - c.radius) * dr + local.z * dz;
            let foot_radius = c.radius + step * dr;
            // Past the apex the foot crosses onto the mirrored nappe,
            // which is a different sheet of the surface.
            if foot_radius < 0.0 {
                return Err(GeomError::Degenerate(
                    "closest point on the cone lies beyond the apex, on the opposite nappe".into(),
                ));
            }
            (local.y.atan2(local.x), step * dz)
        }
        // Reduce to the tube's cross-section circle: project onto the
        // major circle first, then onto the tube around it. Both steps are
        // radial, so the composition is closed form.
        Surface::Torus(t) => {
            positive(t.major_radius, "torus major radius")?;
            positive(t.minor_radius, "torus minor radius")?;
            let local = to_local(&t.frame, point, tolerance)?;
            let ring = (local.x * local.x + local.y * local.y).sqrt();
            if ring <= tolerance.linear() {
                return Err(ambiguous("point lies on the torus axis"));
            }
            let planar = ring - t.major_radius;
            // On the tube's centre circle every cross-section angle is
            // equidistant, the same tie the axis case has one dimension up.
            if planar.abs() <= tolerance.linear() && local.z.abs() <= tolerance.linear() {
                return Err(ambiguous("point lies on the torus tube centre circle"));
            }
            (local.y.atan2(local.x), local.z.atan2(planar))
        }
        // A B-spline needs iterative closest-point with its own seeding and
        // convergence contract, which `axiolid-nurbs` provides as a
        // certified projection. `Surface` is non-exhaustive, so a future
        // variant is refused here by name rather than silently taking an
        // analytic branch that does not describe it.
        _ => {
            return Err(GeomError::Unsupported {
                backend: ScalarSurface::ID,
                operation: axiolid_contracts::Operation::SurfaceEvaluation,
            });
        }
    };
    // Projection promises a point ON the surface, so the parameters must
    // evaluate to one. `invert` is the independent check: it refuses a
    // point further than tolerance from the surface, so a wrong arm here
    // is caught rather than inherited by the caller as bad geometry.
    let landed = evaluate(surface, u, v)?;
    invert(surface, landed, tolerance).map_err(|_| {
        GeomError::Degenerate(
            "projection produced parameters that do not name a point on the surface".into(),
        )
    })?;
    Ok((u, v))
}

/// Express a world point in a frame's local coordinates.
///
/// `place` maps local to world by scaling the frame axes, so the inverse
/// is a projection onto those axes -- but ONLY when they are orthonormal.
/// `Frame3` stores three free vectors and the core documents that
/// algorithms validate orthonormality explicitly, so this checks rather
/// than assumes: on a skewed or scaled frame a dot-product projection is
/// silently wrong, and every parameter derived from it would be wrong by
/// an amount that still round-trips through the same bad frame.
fn to_local(frame: &Frame3, point: Point3, tolerance: Tolerance) -> GeomResult<Vec3> {
    // Validation lives in the core SpaceFrame so surface evaluation, sampled
    // fields, and sectioning cannot drift apart on what a valid frame is.
    // Handedness matters here and was previously unchecked: a mirrored basis
    // passes every unit-length and perpendicularity test while reflecting the
    // surface it parameterises.
    let validated = SpaceFrame::new(frame.origin, frame.x, frame.y, frame.z, tolerance)
        .map_err(|error| GeomError::Degenerate(format!("surface frame is invalid: {error}")))?;
    Ok(validated.to_local(point))
}

/// Angle of a local point about the frame's z axis.
///
/// Refuses points ON the axis. There the whole `u` circle collapses to a
/// single location -- a cone apex, a sphere pole -- so no angle is more
/// correct than any other. Returning zero would look successful and would
/// rotate a swept section arbitrarily about its own path.
fn angle_about_axis(local: Vec3, surface: &str) -> GeomResult<Scalar> {
    let radial = (local.x * local.x + local.y * local.y).sqrt();
    if radial <= 1e-12 {
        return Err(GeomError::Degenerate(format!(
            "{surface} point lies on the axis, where every u names it: \
             the angular parameter is not recoverable"
        )));
    }
    Ok(local.y.atan2(local.x))
}
