//! Fields over a surface's parameters (ADR 0077): building them from two
//! analytic surfaces, and bounding them over parameter boxes.
//!
//! Every analytic surface's point is a [`Field2`] in each world coordinate:
//! powers of a linear parameter, harmonics of an angle. Every analytic
//! surface also has an implicit equation, polynomial in the world
//! coordinates. Substituting the first into the second gives the second
//! surface's equation read in the first one's parameters -- the field whose
//! zero set is the section -- with the product-to-sum rules keeping it a
//! finite sum of the same terms.
//!
//! Bounds over parameter boxes live with the field itself
//! (`axiolid_curve::implicit`), where measure uses them too.

use axiolid_core::{Frame3, Scalar, Vec3};
use axiolid_curve::implicit::grow;
use axiolid_curve::{Basis, Carrier, Field2, RuledCarrier, SeriesField2, TorusCarrier};
use axiolid_surface::Surface;

// --- Algebra ---------------------------------------------------------------

/// Products of two basis terms along one parameter, as terms.
fn product(basis: Basis, a: usize, b: usize) -> Vec<(usize, Scalar)> {
    match basis {
        Basis::Power => vec![(a + b, 1.0)],
        Basis::Fourier => {
            if a == 0 {
                return vec![(b, 1.0)];
            }
            if b == 0 {
                return vec![(a, 1.0)];
            }
            let (wa, ca) = (a.div_ceil(2), a % 2 == 1);
            let (wb, cb) = (b.div_ceil(2), b % 2 == 1);
            // Term index of cos(w x) / sin(w x) for a signed frequency w.
            let cos = |w: i64| -> (usize, Scalar) {
                if w == 0 {
                    (0, 1.0)
                } else {
                    (2 * w.unsigned_abs() as usize - 1, 1.0)
                }
            };
            let sin = |w: i64| -> (usize, Scalar) {
                if w == 0 {
                    (0, 0.0)
                } else {
                    (2 * w.unsigned_abs() as usize, w.signum() as Scalar)
                }
            };
            let (p, q) = (wa as i64, wb as i64);
            let half = |(k, s): (usize, Scalar), f: Scalar| (k, 0.5 * s * f);
            match (ca, cb) {
                // cos p cos q = (cos(p - q) + cos(p + q)) / 2
                (true, true) => vec![half(cos(p - q), 1.0), half(cos(p + q), 1.0)],
                // cos p sin q = (sin(p + q) - sin(p - q)) / 2
                (true, false) => vec![half(sin(p + q), 1.0), half(sin(p - q), -1.0)],
                // sin p cos q = (sin(p + q) + sin(p - q)) / 2
                (false, true) => vec![half(sin(p + q), 1.0), half(sin(p - q), 1.0)],
                // sin p sin q = (cos(p - q) - cos(p + q)) / 2
                (false, false) => vec![half(cos(p - q), 1.0), half(cos(p + q), -1.0)],
            }
        }
    }
}

pub(crate) fn constant(u: Basis, v: Basis, c: Scalar) -> SeriesField2 {
    SeriesField2 {
        u,
        v,
        coefficients: vec![vec![c]],
    }
}

pub(crate) fn add(a: &SeriesField2, b: &SeriesField2, scale_b: Scalar) -> SeriesField2 {
    let mut out = a.coefficients.clone();
    for (i, row) in b.coefficients.iter().enumerate() {
        for (j, &c) in row.iter().enumerate() {
            if c != 0.0 {
                grow(&mut out, i, j);
                out[i][j] += scale_b * c;
            }
        }
    }
    SeriesField2 {
        u: a.u,
        v: a.v,
        coefficients: out,
    }
}

pub(crate) fn mul(a: &SeriesField2, b: &SeriesField2) -> SeriesField2 {
    let mut out: Vec<Vec<Scalar>> = vec![vec![0.0]];
    for (i1, r1) in a.coefficients.iter().enumerate() {
        for (j1, &c1) in r1.iter().enumerate() {
            if c1 == 0.0 {
                continue;
            }
            for (i2, r2) in b.coefficients.iter().enumerate() {
                for (j2, &c2) in r2.iter().enumerate() {
                    if c2 == 0.0 {
                        continue;
                    }
                    for (iu, fu) in product(a.u, i1, i2) {
                        if fu == 0.0 {
                            continue;
                        }
                        for (jv, fv) in product(a.v, j1, j2) {
                            if fv == 0.0 {
                                continue;
                            }
                            grow(&mut out, iu, jv);
                            out[iu][jv] += c1 * c2 * fu * fv;
                        }
                    }
                }
            }
        }
    }
    SeriesField2 {
        u: a.u,
        v: a.v,
        coefficients: out,
    }
}

// --- Surfaces --------------------------------------------------------------

/// The carrier form of an analytic surface, or `None` for a B-spline.
pub(crate) fn carrier_of(surface: &Surface) -> Option<Carrier> {
    Some(match surface {
        Surface::Plane(p) => Carrier::Plane(p.frame),
        Surface::Cylinder(c) => Carrier::Ruled(RuledCarrier {
            frame: c.frame,
            x_radius: c.radius,
            y_radius: c.radius,
            slope: 0.0,
        }),
        Surface::EllipticalCylinder(c) => Carrier::Ruled(RuledCarrier {
            frame: c.frame,
            x_radius: c.semi_axis_x,
            y_radius: c.semi_axis_y,
            slope: 0.0,
        }),
        Surface::Cone(c) => Carrier::Ruled(RuledCarrier {
            frame: c.frame,
            x_radius: c.radius,
            y_radius: c.radius,
            slope: c.semi_angle.tan(),
        }),
        Surface::Sphere(s) => Carrier::Sphere {
            frame: s.frame,
            radius: s.radius,
        },
        Surface::Torus(t) => Carrier::Torus(TorusCarrier {
            frame: t.frame,
            major_radius: t.major_radius,
            minor_radius: t.minor_radius,
        }),
        Surface::BSpline(b) => Carrier::Spline(Box::new(b.clone())),
        _ => return None,
    })
}

/// The bases of a carrier's parameters.
pub(crate) fn bases(carrier: &Carrier) -> (Basis, Basis) {
    match carrier {
        Carrier::Plane(_) | Carrier::Spline(_) => (Basis::Power, Basis::Power),
        Carrier::Ruled(_) => (Basis::Fourier, Basis::Power),
        Carrier::Sphere { .. } | Carrier::Torus(_) => (Basis::Fourier, Basis::Fourier),
    }
}

/// The carrier's point as three fields, one per world coordinate.
fn world(carrier: &Carrier) -> [SeriesField2; 3] {
    let (bu, bv) = bases(carrier);
    // Terms: (u term, v term, vector coefficient).
    let (origin, terms): (Vec3, Vec<(usize, usize, Vec3)>) = match carrier {
        Carrier::Plane(f) => (f.origin, vec![(1, 0, f.x), (0, 1, f.y)]),
        Carrier::Ruled(k) => {
            let f = &k.frame;
            (
                f.origin,
                vec![
                    (0, 1, f.z),
                    (1, 0, f.x * k.x_radius),
                    (1, 1, f.x * k.slope),
                    (2, 0, f.y * k.y_radius),
                    (2, 1, f.y * k.slope),
                ],
            )
        }
        Carrier::Sphere {
            frame: f,
            radius: r,
        } => (
            f.origin,
            vec![(1, 1, f.x * *r), (2, 1, f.y * *r), (0, 2, f.z * *r)],
        ),
        Carrier::Torus(t) => {
            let f = &t.frame;
            let (big, small) = (t.major_radius, t.minor_radius);
            (
                f.origin,
                vec![
                    (1, 0, f.x * big),
                    (1, 1, f.x * small),
                    (2, 0, f.y * big),
                    (2, 1, f.y * small),
                    (0, 2, f.z * small),
                ],
            )
        }
        // Not a series: `section_field` takes a spline carrier elsewhere.
        Carrier::Spline(_) => (Vec3::ZERO, Vec::new()),
    };
    let one = |axis: usize| {
        let mut c = vec![vec![0.0; 3]; 3];
        c[0][0] = origin[axis];
        for (i, j, v) in &terms {
            c[*i][*j] += v[axis];
        }
        SeriesField2 {
            u: bu,
            v: bv,
            coefficients: c,
        }
    };
    [one(0), one(1), one(2)]
}

/// A world point's coordinate along `axis` from `origin`, as a field.
fn local(point: &[SeriesField2; 3], origin: Vec3, axis: Vec3) -> SeriesField2 {
    let axis = axis.normalize();
    let mut out = constant(point[0].u, point[0].v, -origin.dot(axis));
    for (k, field) in point.iter().enumerate() {
        out = add(&out, field, axis[k]);
    }
    out
}

fn frame_locals(point: &[SeriesField2; 3], f: &Frame3) -> [SeriesField2; 3] {
    [
        local(point, f.origin, f.x),
        local(point, f.origin, f.y),
        local(point, f.origin, f.z),
    ]
}

/// `other`'s implicit equation read in `carrier`'s parameters: zero exactly
/// where the carrier's point lies on `other` (for a cone, on either nappe).
/// `None` for a B-spline `other`.
pub(crate) fn section_field(carrier: &Carrier, other: &Surface) -> Option<Field2> {
    if let Carrier::Spline(b) = carrier {
        return crate::spline_field::spline_section_field(b, other);
    }
    let p = world(carrier);
    let (bu, bv) = bases(carrier);
    let square = |f: &SeriesField2| mul(f, f);
    Some(Field2::Series(match other {
        Surface::Plane(q) => local(&p, q.frame.origin, q.frame.z),
        Surface::Cylinder(c) => {
            let [x, y, _] = frame_locals(&p, &c.frame);
            add(
                &add(&square(&x), &square(&y), 1.0),
                &constant(bu, bv, c.radius * c.radius),
                -1.0,
            )
        }
        Surface::EllipticalCylinder(c) => {
            let [x, y, _] = frame_locals(&p, &c.frame);
            let (a2, b2) = (c.semi_axis_x.powi(2), c.semi_axis_y.powi(2));
            let mut f = add(&constant(bu, bv, 0.0), &square(&x), b2);
            f = add(&f, &square(&y), a2);
            add(&f, &constant(bu, bv, a2 * b2), -1.0)
        }
        Surface::Cone(c) => {
            let [x, y, z] = frame_locals(&p, &c.frame);
            let slope = c.semi_angle.tan();
            let radius = add(&constant(bu, bv, c.radius), &z, slope);
            let f = add(&square(&x), &square(&y), 1.0);
            add(&f, &square(&radius), -1.0)
        }
        Surface::Sphere(s) => {
            let [x, y, z] = frame_locals(&p, &s.frame);
            let f = add(&add(&square(&x), &square(&y), 1.0), &square(&z), 1.0);
            add(&f, &constant(bu, bv, s.radius * s.radius), -1.0)
        }
        Surface::Torus(t) => {
            let [x, y, z] = frame_locals(&p, &t.frame);
            let (big, small) = (t.major_radius, t.minor_radius);
            let planar = add(&square(&x), &square(&y), 1.0);
            let all = add(&planar, &square(&z), 1.0);
            let inner = add(&all, &constant(bu, bv, big * big - small * small), 1.0);
            add(&square(&inner), &planar, -4.0 * big * big)
        }
        _ => return None,
    }))
}
