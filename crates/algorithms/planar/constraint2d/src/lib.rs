#![forbid(unsafe_code)]

//! Apollonius and tangent-circle constructions (#159, ledger row B14).
//!
//! [`tangent_circles`] finds every circle tangent to three given points,
//! lines or circles, in any combination -- the classical Apollonius problem
//! (three circles), its ten degenerate variants (mixing in points and
//! lines), and everything in between.
//!
//! # Method
//!
//! A tangency to a candidate circle `(c, r)` is one algebraic equation per
//! constraint, all sharing the same shape after multiplying through:
//!
//! `a * (|c|^2 - r^2) + b . c + d * r + f = 0`
//!
//! where `a = 1` for a point or circle constraint (`a = 0` for a line), and
//! a per-constraint sign choice `s` in `{-1, 1}` selects which of the two
//! tangency senses `d` encodes (internal/external for a circle, which side
//! for a line; a point has only one sense, so its sign does not change the
//! equation). Three such equations, one per constraint, are eliminated in
//! closed form: subtracting the quadratic-type equations from each other
//! cancels their shared `|c|^2 - r^2` term, leaving two equations linear in
//! `(cx, cy, r)`; solving those for `cx` and `cy` as affine functions of `r`
//! and substituting into the remaining quadratic equation gives a quadratic
//! in `r` alone (linear when every constraint is a line). Every sign choice
//! (up to eight) is tried; every real root with a positive radius is a
//! solution, deduplicated and returned in a deterministic order.
//!
//! This is the standard "linearization" solution to Apollonius' problem
//! (see e.g. Coakley 1860, or Gisch and Ribando's 2004 survey for the
//! algebraic form used here), extended uniformly to points and lines rather
//! than special-cased.
//!
//! # Scope
//!
//! Not covered: two circle- or point-type constraints centred at the same
//! point. Their equations then differ only in the radius term, with no
//! `cx`/`cy` coefficient left to eliminate against, which is singular under
//! every sign choice in the elimination this crate uses (a solution family
//! can still exist; finding it needs a different elimination order). This
//! is refused as [`TangencyError::Degenerate`], not answered wrong.

use axiolid_core::{Frame2, Point2, Scalar, Tolerance, Vec2};
use axiolid_curve::Circle2;
use axiolid_linear::Line2;

/// A point, line or circle a solution circle must be tangent to (a point
/// counts as a zero-radius circle, tangent when the solution circle passes
/// through it).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Tangent {
    /// A point the solution circle must pass through.
    Point(Point2),
    /// A line the solution circle must touch.
    Line(Line2),
    /// A circle the solution circle must touch, internally or externally.
    Circle(Circle2),
}

/// Why [`tangent_circles`] could not be solved.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum TangencyError {
    /// A coordinate, direction or radius is not finite.
    NonFinite,
    /// Constraint `index` is a circle whose radius is not positive.
    InvalidRadius {
        /// The constraint's index (0, 1 or 2).
        index: usize,
        /// The radius given.
        radius: Scalar,
    },
    /// Constraint `index` is a line whose direction is shorter than the
    /// linear tolerance.
    ZeroDirection {
        /// The constraint's index (0, 1 or 2).
        index: usize,
    },
    /// Every sign choice's elimination system was singular: the three
    /// constraints do not pin down a family of tangent circles (for
    /// example, three coincident points, three identical circles, or three
    /// parallel lines).
    Degenerate,
}

impl core::fmt::Display for TangencyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFinite => f.write_str("a coordinate, direction or radius is not finite"),
            Self::InvalidRadius { index, radius } => {
                write!(f, "constraint {index} has non-positive radius {radius}")
            }
            Self::ZeroDirection { index } => {
                write!(f, "constraint {index} is a line with no direction")
            }
            Self::Degenerate => f.write_str("no sign choice pins down a family of tangent circles"),
        }
    }
}

impl core::error::Error for TangencyError {}

/// One tangency equation after a sign choice:
/// `a * (cx^2 + cy^2 - r^2) + b * cx + c * cy + d * r + f = 0`.
#[derive(Debug, Clone, Copy)]
struct Equation {
    a: Scalar,
    b: Scalar,
    c: Scalar,
    d: Scalar,
    f: Scalar,
}

fn equation(constraint: &Tangent, sign: Scalar) -> Equation {
    match *constraint {
        Tangent::Point(p) => Equation {
            a: 1.0,
            b: -2.0 * p.x,
            c: -2.0 * p.y,
            d: 0.0,
            f: p.x * p.x + p.y * p.y,
        },
        Tangent::Circle(circle) => {
            let o = circle.frame.origin;
            let rho = circle.radius;
            Equation {
                a: 1.0,
                b: -2.0 * o.x,
                c: -2.0 * o.y,
                d: -2.0 * sign * rho,
                f: o.x * o.x + o.y * o.y - rho * rho,
            }
        }
        Tangent::Line(line) => {
            let n = line.direction.perp().normalize();
            Equation {
                a: 0.0,
                b: n.x,
                c: n.y,
                d: -sign,
                f: -(n.dot(line.origin)),
            }
        }
    }
}

fn sub(lhs: &Equation, rhs: &Equation) -> Equation {
    Equation {
        a: lhs.a - rhs.a,
        b: lhs.b - rhs.b,
        c: lhs.c - rhs.c,
        d: lhs.d - rhs.d,
        f: lhs.f - rhs.f,
    }
}

/// Every circle tangent to `constraints[0]`, `constraints[1]` and
/// `constraints[2]`, in any combination of point, line and circle.
///
/// Solutions are returned in ascending order of centre `x`, then `y`, then
/// radius, and deduplicated: two candidates within the linear tolerance of
/// each other in both centre and radius are the same solution. A well-posed
/// configuration with no tangent circle (for example, three mutually
/// separate circles too far apart for any sign choice to reach) answers
/// with an empty, not an error: only a configuration that pins down no
/// family of circles at all is refused.
///
/// # Errors
///
/// [`TangencyError::NonFinite`], [`TangencyError::InvalidRadius`] for a
/// circle constraint with a non-positive radius,
/// [`TangencyError::ZeroDirection`] for a line constraint shorter than the
/// linear tolerance, and [`TangencyError::Degenerate`] when no sign choice
/// yields a solvable elimination system.
pub fn tangent_circles(
    constraints: [Tangent; 3],
    tolerance: Tolerance,
) -> Result<Vec<Circle2>, TangencyError> {
    validate(&constraints, tolerance)?;
    let linear = tolerance.linear();

    let mut candidates: Vec<(Point2, Scalar)> = Vec::new();
    let mut any_solvable = false;
    for bits in 0u8..8 {
        let signs = [
            if bits & 1 == 0 { -1.0 } else { 1.0 },
            if bits & 2 == 0 { -1.0 } else { 1.0 },
            if bits & 4 == 0 { -1.0 } else { 1.0 },
        ];
        let eqs = [
            equation(&constraints[0], signs[0]),
            equation(&constraints[1], signs[1]),
            equation(&constraints[2], signs[2]),
        ];
        if let Some(roots) = solve_system(&eqs) {
            any_solvable = true;
            for (cx, cy, r) in roots {
                if r.is_finite() && cx.is_finite() && cy.is_finite() && r > linear {
                    candidates.push((Point2::new(cx, cy), r));
                }
            }
        }
    }
    if !any_solvable {
        return Err(TangencyError::Degenerate);
    }

    candidates.sort_by(|a, b| {
        a.0.x
            .total_cmp(&b.0.x)
            .then_with(|| a.0.y.total_cmp(&b.0.y))
            .then_with(|| a.1.total_cmp(&b.1))
    });
    let mut deduped: Vec<(Point2, Scalar)> = Vec::new();
    for (centre, radius) in candidates {
        let same = deduped
            .iter()
            .any(|&(c, r)| (c - centre).length() <= linear && (r - radius).abs() <= linear);
        if !same {
            deduped.push((centre, radius));
        }
    }

    Ok(deduped
        .into_iter()
        .map(|(centre, radius)| Circle2 {
            frame: Frame2 {
                origin: centre,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius,
        })
        .collect())
}

fn validate(constraints: &[Tangent; 3], tolerance: Tolerance) -> Result<(), TangencyError> {
    let linear = tolerance.linear();
    for (index, constraint) in constraints.iter().enumerate() {
        match *constraint {
            Tangent::Point(p) => {
                if !p.is_finite() {
                    return Err(TangencyError::NonFinite);
                }
            }
            Tangent::Line(line) => {
                if !line.origin.is_finite() || !line.direction.is_finite() {
                    return Err(TangencyError::NonFinite);
                }
                if line.direction.length() <= linear {
                    return Err(TangencyError::ZeroDirection { index });
                }
            }
            Tangent::Circle(circle) => {
                if !circle.frame.origin.is_finite() || !circle.radius.is_finite() {
                    return Err(TangencyError::NonFinite);
                }
                if circle.radius <= 0.0 {
                    return Err(TangencyError::InvalidRadius {
                        index,
                        radius: circle.radius,
                    });
                }
            }
        }
    }
    Ok(())
}

/// Solve the 3-equation elimination system for one sign choice: `None` when
/// it is singular, `Some` with zero, one or two `(cx, cy, r)` roots
/// otherwise.
fn solve_system(eqs: &[Equation; 3]) -> Option<Vec<(Scalar, Scalar, Scalar)>> {
    let quad: Vec<usize> = (0..3).filter(|&i| eqs[i].a != 0.0).collect();
    if quad.is_empty() {
        return solve_linear3(eqs).map(|root| vec![root]);
    }
    let base = quad[0];
    let mut lines: Vec<Equation> = Vec::with_capacity(2);
    for i in 0..3 {
        if i == base {
            continue;
        }
        lines.push(if eqs[i].a != 0.0 {
            sub(&eqs[base], &eqs[i])
        } else {
            eqs[i]
        });
    }
    let (m1, n1, m2, n2) = solve_cxy_in_r(&lines[0], &lines[1])?;

    let e = &eqs[base];
    let a_coef = m1 * m1 + m2 * m2 - 1.0;
    let b_coef = 2.0 * (m1 * n1 + m2 * n2) + e.b * m1 + e.c * m2 + e.d;
    let c_coef = n1 * n1 + n2 * n2 + e.b * n1 + e.c * n2 + e.f;
    Some(
        solve_quadratic(a_coef, b_coef, c_coef)
            .into_iter()
            .map(|r| (m1 * r + n1, m2 * r + n2, r))
            .collect(),
    )
}

/// Solve two linear equations `b*cx + c*cy + d*r + f = 0` for `cx` and `cy`
/// as affine functions of `r`: `cx = m1*r + n1`, `cy = m2*r + n2`. `None`
/// when the two equations' `(cx, cy)` coefficients are (near-)parallel.
fn solve_cxy_in_r(l0: &Equation, l1: &Equation) -> Option<(Scalar, Scalar, Scalar, Scalar)> {
    let det = l0.b * l1.c - l1.b * l0.c;
    let scale = (l0.b.hypot(l0.c) * l1.b.hypot(l1.c)).max(1.0);
    if det.abs() <= Scalar::EPSILON * 1e8 * scale {
        return None;
    }
    let m1 = (l1.d * l0.c - l0.d * l1.c) / det;
    let n1 = (l1.f * l0.c - l0.f * l1.c) / det;
    let m2 = (l1.b * l0.d - l0.b * l1.d) / det;
    let n2 = (l1.b * l0.f - l0.b * l1.f) / det;
    Some((m1, n1, m2, n2))
}

/// Solve the 3x3 linear system `b_i*cx + c_i*cy + d_i*r + f_i = 0` for
/// `i = 0, 1, 2` directly (the all-lines case, which has no quadratic
/// equation to substitute into).
fn solve_linear3(eqs: &[Equation; 3]) -> Option<(Scalar, Scalar, Scalar)> {
    let rows = [
        [eqs[0].b, eqs[0].c, eqs[0].d],
        [eqs[1].b, eqs[1].c, eqs[1].d],
        [eqs[2].b, eqs[2].c, eqs[2].d],
    ];
    let rhs = [-eqs[0].f, -eqs[1].f, -eqs[2].f];
    let det = det3(rows);
    let scale = rows
        .iter()
        .map(|row| row[0].hypot(row[1]).hypot(row[2]))
        .fold(1.0, Scalar::max);
    if det.abs() <= Scalar::EPSILON * 1e8 * scale * scale * scale {
        return None;
    }
    let cx = det3(replace_col(rows, 0, rhs)) / det;
    let cy = det3(replace_col(rows, 1, rhs)) / det;
    let r = det3(replace_col(rows, 2, rhs)) / det;
    Some((cx, cy, r))
}

fn det3(m: [[Scalar; 3]; 3]) -> Scalar {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn replace_col(mut m: [[Scalar; 3]; 3], col: usize, values: [Scalar; 3]) -> [[Scalar; 3]; 3] {
    for row in 0..3 {
        m[row][col] = values[row];
    }
    m
}

/// Real roots of `a*x^2 + b*x + c = 0`, falling back to the linear and
/// constant cases as `a`, then `b`, vanish.
fn solve_quadratic(a: Scalar, b: Scalar, c: Scalar) -> Vec<Scalar> {
    let scale = a.abs().max(b.abs()).max(c.abs()).max(1.0);
    if a.abs() <= Scalar::EPSILON * 1e8 * scale {
        if b.abs() <= Scalar::EPSILON * 1e8 * scale {
            return Vec::new();
        }
        return vec![-c / b];
    }
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        return Vec::new();
    }
    if discriminant <= Scalar::EPSILON * 1e8 * scale * scale {
        return vec![-b / (2.0 * a)];
    }
    let sqrt_d = discriminant.sqrt();
    // The numerically stable form (Press et al.): avoids cancellation when
    // `b` and `sqrt_d` are close in magnitude and the same sign.
    let q = if b >= 0.0 {
        -0.5 * (b + sqrt_d)
    } else {
        -0.5 * (b - sqrt_d)
    };
    if q.abs() <= Scalar::EPSILON * 1e8 * scale {
        return vec![-b / (2.0 * a)];
    }
    let mut roots = vec![q / a, c / q];
    roots.sort_by(Scalar::total_cmp);
    roots
}
