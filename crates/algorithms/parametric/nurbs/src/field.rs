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
//! Bounds are interval arithmetic over the terms (a power's or a harmonic's
//! exact range over an interval), tightened by the mean-value form, and
//! widened by a margin that covers the rounding of the sums.

use axiolid_core::{Frame3, Point2, Scalar, Vec3};
use axiolid_curve::{Basis, Carrier, Field2, RuledCarrier, TorusCarrier};
use axiolid_surface::Surface;
use core::f64::consts::{PI, TAU};

/// A closed interval of reals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Range {
    pub(crate) lo: Scalar,
    pub(crate) hi: Scalar,
}

impl Range {
    pub(crate) fn point(x: Scalar) -> Self {
        Self { lo: x, hi: x }
    }

    pub(crate) fn new(a: Scalar, b: Scalar) -> Self {
        Self {
            lo: a.min(b),
            hi: a.max(b),
        }
    }

    fn add(self, o: Self) -> Self {
        Self {
            lo: self.lo + o.lo,
            hi: self.hi + o.hi,
        }
    }

    fn mul(self, o: Self) -> Self {
        let p = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        Self {
            lo: p.iter().copied().fold(Scalar::INFINITY, Scalar::min),
            hi: p.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),
        }
    }

    fn scale(self, c: Scalar) -> Self {
        Self::new(self.lo * c, self.hi * c)
    }

    fn intersect(self, o: Self) -> Self {
        Self {
            lo: self.lo.max(o.lo),
            hi: self.hi.min(o.hi),
        }
    }

    fn widen(self, by: Scalar) -> Self {
        Self {
            lo: self.lo - by,
            hi: self.hi + by,
        }
    }

    /// Whether zero lies in the interval.
    pub(crate) fn straddles_zero(self) -> bool {
        self.lo <= 0.0 && self.hi >= 0.0
    }
}

/// The range of `cos(x)` over `[a, b]`.
fn cos_range(a: Scalar, b: Scalar) -> Range {
    if b - a >= TAU {
        return Range { lo: -1.0, hi: 1.0 };
    }
    let (ca, cb) = (a.cos(), b.cos());
    let mut r = Range::new(ca, cb);
    // A maximum at 2 pi m, a minimum at pi + 2 pi m.
    let m = (a / TAU).ceil();
    if m * TAU <= b {
        r.hi = 1.0;
    }
    let m = ((a - PI) / TAU).ceil();
    if PI + m * TAU <= b {
        r.lo = -1.0;
    }
    r
}

/// The range of basis term `k` over `[a, b]`.
fn term_range(basis: Basis, k: usize, a: Scalar, b: Scalar) -> Range {
    match basis {
        Basis::Power => {
            if k == 0 {
                return Range::point(1.0);
            }
            let (pa, pb) = (a.powi(k as i32), b.powi(k as i32));
            // Monotone unless an even power straddles zero.
            if k % 2 == 1 || a >= 0.0 || b <= 0.0 {
                Range::new(pa, pb)
            } else {
                Range {
                    lo: 0.0,
                    hi: pa.max(pb),
                }
            }
        }
        Basis::Fourier => {
            if k == 0 {
                return Range::point(1.0);
            }
            let w = k.div_ceil(2) as Scalar;
            if k % 2 == 1 {
                cos_range(w * a, w * b)
            } else {
                // sin(y) = cos(y - pi/2).
                cos_range(w * a - 0.5 * PI, w * b - 0.5 * PI)
            }
        }
    }
}

/// A box in the parameter plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Cell {
    pub(crate) lo: Point2,
    pub(crate) hi: Point2,
}

impl Cell {
    pub(crate) fn centre(&self) -> Point2 {
        (self.lo + self.hi) * 0.5
    }
}

/// Naive interval bound of the field over the box.
fn naive(field: &Field2, cell: &Cell) -> Range {
    let mut total = Range::point(0.0);
    let (n, m) = size(field);
    let ru: Vec<Range> = (0..n)
        .map(|k| term_range(field.u, k, cell.lo.x, cell.hi.x))
        .collect();
    let rv: Vec<Range> = (0..m)
        .map(|k| term_range(field.v, k, cell.lo.y, cell.hi.y))
        .collect();
    for (i, row) in field.coefficients.iter().enumerate() {
        for (j, &c) in row.iter().enumerate() {
            if c != 0.0 {
                total = total.add(ru[i].mul(rv[j]).scale(c));
            }
        }
    }
    total
}

/// The margin covering rounding in a sum over the field's terms at the
/// given parameter magnitudes.
fn margin(field: &Field2, cell: &Cell) -> Scalar {
    let (n, m) = size(field);
    // Harmonics never exceed one; powers grow with the parameter.
    let reach = |basis: Basis, x: Scalar, terms: usize| match basis {
        Basis::Fourier => 1.0,
        Basis::Power => (1.0 + x.abs()).powi(terms.saturating_sub(1) as i32),
    };
    let scale = field.magnitude()
        * reach(field.u, cell.lo.x.abs().max(cell.hi.x.abs()), n)
        * reach(field.v, cell.lo.y.abs().max(cell.hi.y.abs()), m);
    64.0 * Scalar::EPSILON * scale
}

/// A bound certain to hold the field's values over the box (up to the
/// rounding margin included in it).
pub(crate) fn bound(field: &Field2, du: &Field2, dv: &Field2, cell: &Cell) -> Range {
    let direct = naive(field, cell);
    let c = cell.centre();
    let (hu, hv) = (0.5 * (cell.hi.x - cell.lo.x), 0.5 * (cell.hi.y - cell.lo.y));
    let mean = Range::point(field.value(c))
        .add(naive(du, cell).mul(Range { lo: -hu, hi: hu }))
        .add(naive(dv, cell).mul(Range { lo: -hv, hi: hv }));
    direct.intersect(mean).widen(margin(field, cell))
}

/// A bound of the field over the box without the mean-value tightening.
pub(crate) fn bound_simple(field: &Field2, cell: &Cell) -> Range {
    naive(field, cell).widen(margin(field, cell))
}

fn size(field: &Field2) -> (usize, usize) {
    (
        field.coefficients.len(),
        field.coefficients.iter().map(Vec::len).max().unwrap_or(0),
    )
}

/// The partial derivative of a field along `u` (`along_u`) or `v`.
pub(crate) fn partial(field: &Field2, along_u: bool) -> Field2 {
    let (n, m) = size(field);
    let mut out = vec![vec![0.0; m]; n];
    let basis = if along_u { field.u } else { field.v };
    let map = |k: usize| -> Option<(usize, Scalar)> {
        match basis {
            Basis::Power => (k > 0).then(|| (k - 1, k as Scalar)),
            Basis::Fourier => {
                if k == 0 {
                    None
                } else {
                    let w = k.div_ceil(2) as Scalar;
                    if k % 2 == 1 {
                        // cos(w x) -> -w sin(w x)
                        Some((k + 1, -w))
                    } else {
                        // sin(w x) -> w cos(w x)
                        Some((k - 1, w))
                    }
                }
            }
        }
    };
    for (i, row) in field.coefficients.iter().enumerate() {
        for (j, &c) in row.iter().enumerate() {
            if c == 0.0 {
                continue;
            }
            if along_u {
                if let Some((k, f)) = map(i) {
                    grow(&mut out, k, j);
                    out[k][j] += c * f;
                }
            } else if let Some((k, f)) = map(j) {
                grow(&mut out, i, k);
                out[i][k] += c * f;
            }
        }
    }
    Field2 {
        u: field.u,
        v: field.v,
        coefficients: out,
    }
}

fn grow(c: &mut Vec<Vec<Scalar>>, i: usize, j: usize) {
    if c.len() <= i {
        let m = c.first().map_or(0, Vec::len);
        c.resize(i + 1, vec![0.0; m]);
    }
    if c[0].len() <= j {
        for row in c.iter_mut() {
            row.resize(j + 1, 0.0);
        }
    }
    for row in c.iter_mut() {
        if row.len() <= j {
            row.resize(j + 1, 0.0);
        }
    }
}

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

pub(crate) fn constant(u: Basis, v: Basis, c: Scalar) -> Field2 {
    Field2 {
        u,
        v,
        coefficients: vec![vec![c]],
    }
}

pub(crate) fn add(a: &Field2, b: &Field2, scale_b: Scalar) -> Field2 {
    let mut out = a.coefficients.clone();
    for (i, row) in b.coefficients.iter().enumerate() {
        for (j, &c) in row.iter().enumerate() {
            if c != 0.0 {
                grow(&mut out, i, j);
                out[i][j] += scale_b * c;
            }
        }
    }
    Field2 {
        u: a.u,
        v: a.v,
        coefficients: out,
    }
}

pub(crate) fn mul(a: &Field2, b: &Field2) -> Field2 {
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
    Field2 {
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
        _ => return None,
    })
}

/// The bases of a carrier's parameters.
pub(crate) fn bases(carrier: &Carrier) -> (Basis, Basis) {
    match carrier {
        Carrier::Plane(_) => (Basis::Power, Basis::Power),
        Carrier::Ruled(_) => (Basis::Fourier, Basis::Power),
        Carrier::Sphere { .. } | Carrier::Torus(_) => (Basis::Fourier, Basis::Fourier),
    }
}

/// The carrier's point as three fields, one per world coordinate.
fn world(carrier: &Carrier) -> [Field2; 3] {
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
    };
    let one = |axis: usize| {
        let mut c = vec![vec![0.0; 3]; 3];
        c[0][0] = origin[axis];
        for (i, j, v) in &terms {
            c[*i][*j] += v[axis];
        }
        Field2 {
            u: bu,
            v: bv,
            coefficients: c,
        }
    };
    [one(0), one(1), one(2)]
}

/// A world point's coordinate along `axis` from `origin`, as a field.
fn local(point: &[Field2; 3], origin: Vec3, axis: Vec3) -> Field2 {
    let axis = axis.normalize();
    let mut out = constant(point[0].u, point[0].v, -origin.dot(axis));
    for (k, field) in point.iter().enumerate() {
        out = add(&out, field, axis[k]);
    }
    out
}

fn frame_locals(point: &[Field2; 3], f: &Frame3) -> [Field2; 3] {
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
    let p = world(carrier);
    let (bu, bv) = bases(carrier);
    let square = |f: &Field2| mul(f, f);
    Some(match other {
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
    })
}
