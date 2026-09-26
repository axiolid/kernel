//! The general boolean with B-spline faces (#167, ADR 0075 stage 3, ADR
//! 0077): a box with a bi-quadratic Bezier roof against pipes, boxes and a
//! sphere.
//!
//! The roof is `z = sum h_ij B_i(x/2) B_j(y/2)` over `[0, 2]^2`, a
//! polynomial, so every volume below is an integral of a polynomial over a
//! disc or a rectangle, evaluated here by Gauss quadrature exactly.

mod common;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::boolean;
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_core::{BooleanOperator, Point2, Tolerance};
use axiolid_measure::exact_properties;
use axiolid_overlay::ArcRing;
use common::{spline_box, spline_box_at, spline_box_grid, spline_box_sized};

const PI: f64 = std::f64::consts::PI;
const H: [[f64; 3]; 3] = [[1.0, 1.3, 1.1], [1.2, 1.8, 1.4], [0.9, 1.5, 1.2]];

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, tol())
        .expect("measurable")
        .signed_volume
}

fn audited(brep: &ExactBRep, what: &str) {
    let health = geometric_audit(brep, tol());
    assert!(health.is_consistent(), "{what}: {:?}", health.defects());
    let topology = axiolid_topology::audit_brep(brep.topology());
    assert!(topology.is_closed_manifold(), "{what}: {topology:?}");
}

fn run(a: &ExactBRep, b: &ExactBRep, op: BooleanOperator) -> ExactBRep {
    let result = boolean(a, b, op, tol()).unwrap_or_else(|e| panic!("{op:?}: {e}"));
    audited(&result, &format!("{op:?}"));
    result
}

fn close(what: &str, got: f64, want: f64, rel: f64) {
    assert!(
        (got - want).abs() <= rel * want.abs().max(1.0),
        "{what}: expected {want}, got {got}"
    );
}

fn check(a: &ExactBRep, b: &ExactBRep, want: f64, rel: f64) {
    let u = volume(&run(a, b, BooleanOperator::Union));
    let i = volume(&run(a, b, BooleanOperator::Intersection));
    let d = volume(&run(a, b, BooleanOperator::Difference));
    close("intersection", i, want, rel);
    close("identity", u + i, volume(a) + volume(b), rel);
    close("difference", d, volume(a) - i, rel);
}

/// The roof height over `(x, y)`.
fn roof(x: f64, y: f64) -> f64 {
    roof_of(&H, x, y)
}

/// A roof with control heights `h` over `(x, y)`, in its own box's frame.
#[allow(clippy::needless_range_loop)]
fn roof_of(h: &[[f64; 3]; 3], x: f64, y: f64) -> f64 {
    let b = |t: f64| [(1.0 - t) * (1.0 - t), 2.0 * t * (1.0 - t), t * t];
    let (bu, bv) = (b(x / 2.0), b(y / 2.0));
    let mut z = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            z += h[i][j] * bu[i] * bv[j];
        }
    }
    z
}

/// Gauss-Legendre with 12 points: exact for the degree-8 polynomials here.
fn gauss(a: f64, b: f64, f: impl Fn(f64) -> f64) -> f64 {
    const X: [f64; 6] = [
        0.125_233_408_511_468_9,
        0.367_831_498_998_180_2,
        0.587_317_954_286_617_4,
        0.769_902_674_194_304_7,
        0.904_117_256_370_474_9,
        0.981_560_634_246_719_3,
    ];
    const W: [f64; 6] = [
        0.249_147_045_813_402_8,
        0.233_492_536_538_354_8,
        0.203_167_426_723_065_9,
        0.160_078_328_543_346_2,
        0.106_939_325_995_318_4,
        0.047_175_336_386_511_8,
    ];
    let (m, h) = (0.5 * (a + b), 0.5 * (b - a));
    let mut s = 0.0;
    for k in 0..6 {
        s += W[k] * (f(m + h * X[k]) + f(m - h * X[k]));
    }
    s * h
}

fn pipe(x: f64, y: f64, r: f64, bottom: f64, top: f64) -> ExactBRep {
    let p = ArcPrism {
        section: ArcRing::circle(Point2::new(x, y), r),
        bottom,
        top,
    };
    boolean_arc_prisms_exact(&p, &p, BooleanOperator::Intersection, tol()).expect("a pipe")
}

fn block(x: (f64, f64), y: (f64, f64), z: (f64, f64)) -> ExactBRep {
    let p = ArcPrism {
        section: ArcRing::from_points(&[
            Point2::new(x.0, y.0),
            Point2::new(x.1, y.0),
            Point2::new(x.1, y.1),
            Point2::new(x.0, y.1),
        ]),
        bottom: z.0,
        top: z.1,
    };
    boolean_arc_prisms_exact(&p, &p, BooleanOperator::Intersection, tol()).expect("a box")
}

#[test]
fn the_spline_roofed_box_measures_and_audits() {
    let solid = spline_box(H);
    audited(&solid, "spline box");
    let want = gauss(0.0, 2.0, |x| gauss(0.0, 2.0, |y| roof(x, y)));
    close("volume", volume(&solid), want, 1e-12);
}

#[test]
fn a_pipe_through_the_spline_roof() {
    let solid = spline_box(H);
    let (cx, cy, r) = (1.1, 0.9, 0.5);
    let bore = pipe(cx, cy, r, -1.0, 3.0);
    // The box over the pipe's disc: the roof height integrated over it.
    let want = gauss(0.0, 2.0 * PI, |t| {
        gauss(0.0, r, |s| roof(cx + s * t.cos(), cy + s * t.sin()) * s)
    });
    check(&solid, &bore, want, 1e-9);
}

#[test]
fn a_block_cutting_the_spline_roof() {
    let solid = spline_box(H);
    // x in [0.5, 1.5], y in [0.4, 3], z in [0.2, 1.5]: the roof dips below
    // z = 1.5 in places, so the block meets it along a traced curve.
    let cut = block((0.5, 1.5), (0.4, 3.0), (0.2, 1.5));
    // Over each x the roof is a quadratic in y: split the y range where it
    // crosses z = 1.5, so every piece is a polynomial and Gauss is exact on
    // it; the outer integral in x is split finely around the kinks.
    let inner = |x: f64| {
        let q = |y: f64| roof(x, y) - 1.5;
        // q(y) = a y^2 + b y + c from three values.
        let (q0, q1, q2) = (q(0.0), q(1.0), q(2.0));
        let a = 0.5 * (q2 - 2.0 * q1 + q0);
        let b = q1 - q0 - a;
        let c = q0;
        let mut cuts = vec![0.4, 2.0];
        let d = b * b - 4.0 * a * c;
        if a.abs() > 1e-14 && d > 0.0 {
            for r in [(-b + d.sqrt()) / (2.0 * a), (-b - d.sqrt()) / (2.0 * a)] {
                if r > 0.4 && r < 2.0 {
                    cuts.push(r);
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.windows(2)
            .map(|w| {
                let mid = 0.5 * (w[0] + w[1]);
                if roof(x, mid) < 1.5 {
                    gauss(w[0], w[1], |y| roof(x, y) - 0.2)
                } else {
                    1.3 * (w[1] - w[0])
                }
            })
            .sum::<f64>()
    };
    let panels = 400;
    let want: f64 = (0..panels)
        .map(|k| {
            let (a, b) = (
                0.5 + k as f64 / panels as f64,
                0.5 + (k + 1) as f64 / panels as f64,
            );
            gauss(a, b, inner)
        })
        .sum();
    check(&solid, &cut, want, 1e-8);
}

#[test]
fn a_ball_dipping_into_the_spline_roof() {
    // The sphere/roof section is traced on the roof; on the sphere it is
    // the same space curve read in the sphere's parameters. No closed form:
    // the identities, and the intersection below the sphere's own volume.
    let solid = spline_box(H);
    let ball = common::sphere(axiolid_core::Point3::new(1.0, 1.1, 1.9), 0.6);
    let u = volume(&run(&solid, &ball, BooleanOperator::Union));
    let i = volume(&run(&solid, &ball, BooleanOperator::Intersection));
    let d = volume(&run(&solid, &ball, BooleanOperator::Difference));
    close("identity", u + i, volume(&solid) + volume(&ball), 1e-9);
    close("difference", d, volume(&solid) - i, 1e-9);
    assert!(i > 0.0 && i < volume(&ball));
    // Directly: the ball below the roof, by quadrature over its shadow
    // (the roof never dips below the ball's lowest point here).
    let (cx, cy, cz, r) = (1.0, 1.1, 1.9, 0.6);
    let want = gauss(0.0, 2.0 * PI, |t| {
        // Split the radius finely: the integrand has a kink where the roof
        // leaves the ball.
        let panels = 40;
        (0..panels)
            .map(|k| {
                let (a, b) = (
                    r * k as f64 / panels as f64,
                    r * (k + 1) as f64 / panels as f64,
                );
                gauss(a, b, |s| {
                    let (x, y) = (cx + s * t.cos(), cy + s * t.sin());
                    let half = (r * r - s * s).max(0.0).sqrt();
                    let top = roof(x, y).min(cz + half);
                    (top - (cz - half)).max(0.0) * s
                })
            })
            .sum::<f64>()
    });
    close("intersection", i, want, 1e-5);
}

#[test]
fn two_spline_roofed_boxes_meet_roof_to_roof() {
    // The second box is moved by (0.5, 0.3, -0.2) and its roof crosses the
    // first one's: the roofs meet along a traced pair section (ADR 0077).
    const G: [[f64; 3]; 3] = [[1.6, 1.0, 1.5], [1.1, 0.9, 1.3], [1.7, 1.2, 1.0]];
    let (dx, dy, dz) = (0.5, 0.3, -0.2);
    let first = spline_box(H);
    let second = spline_box_at(axiolid_core::Vec3::new(dx, dy, dz), G);
    let other = |x: f64, y: f64| roof_of(&G, x - dx, y - dy) + dz;
    // Over each x both roofs are quadratics in y, so their difference is
    // one too: split at its roots, and every piece is a polynomial.
    let inner = |x: f64| {
        let q = |y: f64| roof(x, y) - other(x, y);
        let (q0, q1, q2) = (q(0.0), q(1.0), q(2.0));
        let a = 0.5 * (q2 - 2.0 * q1 + q0);
        let b = q1 - q0 - a;
        let c = q0;
        let mut cuts = vec![dy, 2.0];
        let d = b * b - 4.0 * a * c;
        if a.abs() > 1e-14 && d > 0.0 {
            for r in [(-b + d.sqrt()) / (2.0 * a), (-b - d.sqrt()) / (2.0 * a)] {
                if r > dy && r < 2.0 {
                    cuts.push(r);
                }
            }
        } else if a.abs() <= 1e-14 && b != 0.0 {
            let r = -c / b;
            if r > dy && r < 2.0 {
                cuts.push(r);
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.windows(2)
            .map(|w| gauss(w[0], w[1], |y| roof(x, y).min(other(x, y))))
            .sum::<f64>()
    };
    let panels = 600;
    let want: f64 = (0..panels)
        .map(|k| {
            let (a, b) = (
                dx + (2.0 - dx) * k as f64 / panels as f64,
                dx + (2.0 - dx) * (k + 1) as f64 / panels as f64,
            );
            gauss(a, b, inner)
        })
        .sum();
    check(&first, &second, want, 1e-7);
}

#[test]
fn a_spline_bump_pierces_the_spline_roof_in_a_closed_loop() {
    // A smaller box under the roof, its own roof a bump poking through the
    // first one's inside the face: the roofs meet in one closed pair section
    // touching no edge, which leaves a hole in each roof.
    const G: [[f64; 3]; 3] = [[0.8, 0.8, 0.8], [0.8, 4.0, 0.8], [0.8, 0.8, 0.8]];
    let (x0, y0, z0, step) = (0.5, 0.5, 0.1, 0.5);
    let first = spline_box(H);
    let bump = spline_box_sized(axiolid_core::Vec3::new(x0, y0, z0), step, G);
    let other = |x: f64, y: f64| roof_of(&G, (x - x0) / step, (y - y0) / step) + z0;
    let (x1, y1) = (x0 + 2.0 * step, y0 + 2.0 * step);
    let inner = |x: f64| {
        let q = |y: f64| roof(x, y) - other(x, y);
        let (q0, q1, q2) = (q(0.0), q(1.0), q(2.0));
        let a = 0.5 * (q2 - 2.0 * q1 + q0);
        let b = q1 - q0 - a;
        let c = q0;
        let mut cuts = vec![y0, y1];
        let d = b * b - 4.0 * a * c;
        if a.abs() > 1e-14 && d > 0.0 {
            for r in [(-b + d.sqrt()) / (2.0 * a), (-b - d.sqrt()) / (2.0 * a)] {
                if r > y0 && r < y1 {
                    cuts.push(r);
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.windows(2)
            .map(|w| gauss(w[0], w[1], |y| roof(x, y).min(other(x, y)) - z0))
            .sum::<f64>()
    };
    let panels = 600;
    let want: f64 = (0..panels)
        .map(|k| {
            let (a, b) = (
                x0 + (x1 - x0) * k as f64 / panels as f64,
                x0 + (x1 - x0) * (k + 1) as f64 / panels as f64,
            );
            gauss(a, b, inner)
        })
        .sum();
    check(&first, &bump, want, 1e-7);
}

#[test]
fn a_roof_tangent_to_a_block_top_along_a_line_and_crossing_it() {
    // The roof z = 1 + (x - 1)^3 / 2 over [0, 2]^2 (cubic in x; the
    // Bernstein coefficients of s^3 over [-1, 1] are -1, 1, -1, 1) is
    // tangent to the plane z = 1 along the line x = 1 and crosses it there.
    // The block [0.5, 1.5]^2 x [0.2, 1] has that plane as its top: the
    // section is the line of contact, traced as a second derivative's
    // zeros (ADR 0077).
    let c = [-1.0, 1.0, -1.0, 1.0];
    let h: Vec<Vec<f64>> = c.iter().map(|ci| vec![1.0 + 0.5 * ci; 2]).collect();
    let roofed = spline_box_grid(axiolid_core::Vec3::ZERO, 2.0, &h);
    let top = block((0.5, 1.5), (0.5, 1.5), (0.2, 1.0));
    // Under the roof where it is lower (x < 1), under the top elsewhere.
    let want = 0.4 + 0.5 * (-(0.5f64.powi(4)) / 4.0) + 0.4;
    check(&roofed, &top, want, 1e-9);
}
