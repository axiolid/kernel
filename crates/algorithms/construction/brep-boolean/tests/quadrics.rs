//! The general boolean between spheres and cylinders at any angle (#167
//! stage 2): pipe tees, Steinmetz solids, bored and overlapping spheres.
//!
//! Every volume has an oracle independent of the kernel: a closed form, or
//! a one- or two-dimensional integral evaluated here by Gauss-Legendre
//! quadrature on a smooth integrand.

mod common;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::boolean;
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_core::{BooleanOperator, Point2, Point3, Tolerance, Vec3};
use axiolid_measure::exact_properties;
use axiolid_overlay::ArcRing;
use common::{moved, sphere, Motion};

const PI: f64 = std::f64::consts::PI;

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

/// All three operators; the intersection against `want`, the rest against
/// the identities.
fn check(a: &ExactBRep, b: &ExactBRep, want: f64, rel: f64) {
    let u = volume(&run(a, b, BooleanOperator::Union));
    let i = volume(&run(a, b, BooleanOperator::Intersection));
    let d = volume(&run(a, b, BooleanOperator::Difference));
    close("intersection", i, want, rel);
    let (va, vb) = (volume(a), volume(b));
    close("identity", u + i, va + vb, rel);
    close("difference", d, va - i, rel);
}

/// A vertical pipe (world z) of radius `r` about `(x, y)`.
fn pipe(x: f64, y: f64, r: f64, bottom: f64, top: f64) -> ExactBRep {
    let p = ArcPrism {
        section: ArcRing::circle(Point2::new(x, y), r),
        bottom,
        top,
    };
    boolean_arc_prisms_exact(&p, &p, BooleanOperator::Intersection, tol()).expect("a pipe")
}

/// A pipe along world x, through the origin at height `z`.
fn cross_pipe(r: f64, half_length: f64, z: f64) -> ExactBRep {
    let upright = pipe(0.0, 0.0, r, -half_length, half_length);
    // Turn z onto x: a quarter turn about y, exactly.
    moved(
        &upright,
        Motion {
            axis: Vec3::Y,
            angle: Motion::QUARTER_Y,
            shift: Vec3::new(0.0, 0.0, z),
        },
    )
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

/// Gauss-Legendre nodes and weights on `[-1, 1]` by Newton on `P_n`.
fn gauss(n: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-16 {
                break;
            }
        }
        let (mut p0, mut p1) = (1.0, x);
        for k in 2..=n {
            let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
            p0 = p1;
            p1 = p2;
        }
        let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
        out.push((x, 2.0 / ((1.0 - x * x) * dp * dp)));
    }
    out
}

fn integrate(a: f64, b: f64, f: impl Fn(f64) -> f64) -> f64 {
    // Composite: 16 panels of 40 points.
    let rule = gauss(40);
    let mut total = 0.0;
    let panels = 16;
    for k in 0..panels {
        let (lo, hi) = (
            a + (b - a) * k as f64 / panels as f64,
            a + (b - a) * (k + 1) as f64 / panels as f64,
        );
        let (m, h) = (0.5 * (lo + hi), 0.5 * (hi - lo));
        total += rule.iter().map(|(x, w)| w * f(m + h * x)).sum::<f64>() * h;
    }
    total
}

#[test]
fn a_sphere_fixture_measures_and_audits() {
    let r = 1.3;
    let ball = sphere(Point3::new(0.2, -0.1, 0.4), r);
    audited(&ball, "sphere");
    close("sphere", volume(&ball), 4.0 / 3.0 * PI * r * r * r, 1e-12);
}

#[test]
fn half_a_sphere() {
    let r = 1.0;
    let ball = sphere(Point3::ZERO, r);
    let upper = block((-2.0, 2.0), (-2.0, 2.0), (0.0, 2.0));
    check(&ball, &upper, 2.0 / 3.0 * PI, 1e-9);
}

#[test]
fn a_sphere_bored_along_its_axis_leaves_a_napkin_ring() {
    let (big, r) = (1.0, 0.5);
    let ball = sphere(Point3::ZERO, big);
    let bore = pipe(0.0, 0.0, r, -2.0, 2.0);
    let h = 2.0 * (big * big - r * r).sqrt();
    let ring = PI * h * h * h / 6.0;
    check(&ball, &bore, 4.0 / 3.0 * PI - ring, 1e-9);
}

#[test]
fn a_sphere_bored_off_its_axis() {
    // The pipe (radius 0.4 about (0.3, 0.2)) stays inside the sphere's
    // shadow, so the sphere over the pipe's disc is
    // int 2 sqrt(R^2 - x^2 - y^2) over that disc, in polar coordinates about
    // the pipe's centre.
    let (big, r, cx, cy) = (1.0, 0.4, 0.3, 0.2);
    let ball = sphere(Point3::ZERO, big);
    let bore = pipe(cx, cy, r, -2.0, 2.0);
    let want = integrate(0.0, 2.0 * PI, |t| {
        integrate(0.0, r, |s| {
            let (x, y) = (cx + s * t.cos(), cy + s * t.sin());
            2.0 * (big * big - x * x - y * y).sqrt() * s
        })
    });
    check(&ball, &bore, want, 1e-8);
}

#[test]
fn two_overlapping_spheres_share_a_lens() {
    let (big, r, d) = (1.0, 0.8, 1.1);
    let a = sphere(Point3::ZERO, big);
    let b = sphere(Point3::new(d, 0.0, 0.0), r);
    let lens = PI
        * (big + r - d).powi(2)
        * (d * d + 2.0 * d * r - 3.0 * r * r + 2.0 * d * big + 6.0 * r * big - 3.0 * big * big)
        / (12.0 * d);
    check(&a, &b, lens, 1e-9);
}

#[test]
fn a_steinmetz_solid() {
    // Two unit pipes at right angles: their intersection is 16/3.
    let upright = pipe(0.0, 0.0, 1.0, -3.0, 3.0);
    let across = cross_pipe(1.0, 3.0, 0.0);
    check(&upright, &across, 16.0 / 3.0, 1e-9);
}

#[test]
fn a_pipe_tee_with_unequal_radii() {
    // Radius 1 upright, 0.6 across: the intersection is
    // 8 int_0^{pi/2} r^2 cos^2 t sqrt(R^2 - r^2 sin^2 t) dt.
    let (big, r) = (1.0, 0.6);
    let upright = pipe(0.0, 0.0, big, -3.0, 3.0);
    let across = cross_pipe(r, 3.0, 0.4);
    let want = 8.0
        * integrate(0.0, 0.5 * PI, |t| {
            r * r * t.cos().powi(2) * (big * big - r * r * t.sin().powi(2)).sqrt()
        });
    check(&upright, &across, want, 1e-8);
}

#[test]
fn a_pipe_through_a_sphere_at_an_angle() {
    // No closed form: the identities, and the intersection between the
    // bore through the centre and the whole pipe.
    let ball = sphere(Point3::new(0.1, 0.0, 0.2), 1.0);
    let upright = pipe(0.0, 0.0, 0.45, -3.0, 3.0);
    let tilted = moved(
        &upright,
        Motion {
            axis: Vec3::new(1.0, 0.4, 0.0),
            angle: 0.6,
            shift: Vec3::new(0.2, -0.1, 0.0),
        },
    );
    let u = volume(&run(&ball, &tilted, BooleanOperator::Union));
    let i = volume(&run(&ball, &tilted, BooleanOperator::Intersection));
    let d = volume(&run(&ball, &tilted, BooleanOperator::Difference));
    close("identity", u + i, volume(&ball) + volume(&tilted), 1e-9);
    close("difference", d, volume(&ball) - i, 1e-9);
    assert!(i > 0.0 && i < PI * 0.45 * 0.45 * 2.0);
}
