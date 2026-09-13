//! Variable-radius fillet surface (ADR 0051).
//!
//! The blend is checked by EVALUATING the rational surface and comparing
//! against the closed-form fillet: at height v the section is a circle of
//! radius r(v) centred on the bisector at r(v)/sin(theta/2). A surface that
//! merely looks plausible fails this.

use axiolid_construct::feature::BlendCorner;
use axiolid_construct::fillet_variable::tapered_blend_surface;
use axiolid_core::{Point2, Point3};

/// de Boor evaluation of a rational Bezier patch (u quadratic, v linear).
///
/// Written out rather than calling a library evaluator: the point is to
/// check the control net independently, and reusing the same evaluation
/// code that built it would not do that.
fn evaluate(net: &[Vec<Point3>], weights: &[Vec<f64>], u: f64, v: f64) -> Point3 {
    let bu = [(1.0 - u) * (1.0 - u), 2.0 * u * (1.0 - u), u * u];
    let bv = [1.0 - v, v];
    let mut num = [0.0f64; 3];
    let mut den = 0.0f64;
    for (i, bu_i) in bu.iter().enumerate() {
        for (j, bv_j) in bv.iter().enumerate() {
            let w = bu_i * bv_j * weights[i][j];
            num[0] += w * net[i][j].x;
            num[1] += w * net[i][j].y;
            num[2] += w * net[i][j].z;
            den += w;
        }
    }
    Point3::new(num[0] / den, num[1] / den, num[2] / den)
}

#[test]
fn blend_matches_closed_form() {
    let theta: f64 = std::f64::consts::FRAC_PI_2;
    let half = theta / 2.0;
    let bottom_radius = 0.30_f64;
    let top_radius = 0.70_f64;
    let height = 2.0_f64;

    // Wedge with its corner at the origin, bisector along +x.
    let corner = Point2::new(0.0, 0.0);
    let centre_at = |r: f64| Point2::new(r / half.sin(), 0.0);

    let bottom = BlendCorner {
        centre: centre_at(bottom_radius),
        start: Point2::new(
            bottom_radius / half.tan() * half.cos(),
            bottom_radius / half.tan() * half.sin(),
        ),
        end: Point2::new(
            bottom_radius / half.tan() * half.cos(),
            -bottom_radius / half.tan() * half.sin(),
        ),
        sweep: std::f64::consts::PI - theta,
    };
    let top = BlendCorner {
        centre: centre_at(top_radius),
        start: Point2::new(
            top_radius / half.tan() * half.cos(),
            top_radius / half.tan() * half.sin(),
        ),
        end: Point2::new(
            top_radius / half.tan() * half.cos(),
            -top_radius / half.tan() * half.sin(),
        ),
        sweep: std::f64::consts::PI - theta,
    };
    let _ = corner;

    let surface = tapered_blend_surface(&bottom, &top, height).expect("tapered blend");
    let weights = surface.weights.clone().expect("rational");

    // At each height the section must be a circle of radius r(v) centred on
    // the bisector -- the closed form, computed here independently.
    let mut worst = 0.0_f64;
    for vi in 0..=8 {
        let v = f64::from(vi) / 8.0;
        let radius = bottom_radius + (top_radius - bottom_radius) * v;
        let centre = centre_at(radius);
        for ui in 0..=8 {
            let u = f64::from(ui) / 8.0;
            let point = evaluate(&surface.control_points, &weights, u, v);
            let dx = point.x - centre.x;
            let dy = point.y - centre.y;
            worst = worst.max((dx.hypot(dy) - radius).abs());
            worst = worst.max((point.z - height * v).abs());
        }
    }
    assert!(
        worst < 1e-12,
        "the tapered blend must equal the closed-form fillet, worst error {worst:e}"
    );
}

/// Build a wedge blend corner at the given radius and height.
fn wedge(theta: f64, radius: f64) -> BlendCorner {
    let half = theta / 2.0;
    let setback = radius / half.tan();
    BlendCorner {
        centre: Point2::new(radius / half.sin(), 0.0),
        start: Point2::new(setback * half.cos(), setback * half.sin()),
        end: Point2::new(setback * half.cos(), -setback * half.sin()),
        // Signed sweep from start to end about the centre. Derived rather
        // than assumed: the sign convention is what the first draft of this
        // test got wrong, and hardcoding it hid a real disagreement.
        sweep: std::f64::consts::PI - theta,
    }
}

#[test]
fn blend_is_tangent_to_both_walls() {
    // Tangency is the property that makes a fillet a fillet: at the seam the
    // blend must touch each wall plane, not merely come close.
    for degrees in [60.0_f64, 90.0, 120.0, 150.0] {
        let theta = degrees.to_radians();
        let half = theta / 2.0;
        let bottom_radius = 0.25_f64;
        let top_radius = 0.55_f64;
        let surface =
            tapered_blend_surface(&wedge(theta, bottom_radius), &wedge(theta, top_radius), 1.5)
                .expect("tapered blend");
        let weights = surface.weights.clone().expect("rational");

        // Wall planes through the origin, normals perpendicular to each face.
        let upper = (-half.sin(), half.cos());
        let lower = (-half.sin(), -half.cos());

        for vi in 0..=6 {
            let v = f64::from(vi) / 6.0;
            let radius = bottom_radius + (top_radius - bottom_radius) * v;
            // u = 0 is the start seam, u = 1 the end seam.
            let start = evaluate(&surface.control_points, &weights, 0.0, v);
            let end = evaluate(&surface.control_points, &weights, 1.0, v);
            let d_start = (start.x * upper.0 + start.y * upper.1).abs();
            let d_end = (end.x * lower.0 + end.y * lower.1).abs();
            assert!(
                d_start < 1e-12 && d_end < 1e-12,
                "seams must lie on the wall planes at theta={degrees}, v={v}: {d_start:e} {d_end:e}"
            );
            let centre = (radius / half.sin(), 0.0);
            let rs = (start.x - centre.0).hypot(start.y - centre.1);
            assert!(
                (rs - radius).abs() < 1e-12,
                "seam must sit at radius {radius} from the centre, got {rs}"
            );
        }
    }
}

#[test]
fn a_reversed_taper_is_accepted_symmetrically() {
    // Shrinking upward is as valid as growing upward; refusing one direction
    // would be an arbitrary restriction, so both must build.
    let theta = std::f64::consts::FRAC_PI_2;
    let growing = tapered_blend_surface(&wedge(theta, 0.2), &wedge(theta, 0.6), 1.0);
    let shrinking = tapered_blend_surface(&wedge(theta, 0.6), &wedge(theta, 0.2), 1.0);
    assert!(growing.is_ok() && shrinking.is_ok());
}

#[test]
fn a_constant_radius_taper_is_a_cylinder_section() {
    // The degenerate taper must still be exact: equal radii means every
    // height has the same circle, which is the constant-radius fillet.
    let theta = std::f64::consts::FRAC_PI_2;
    let radius = 0.4_f64;
    let surface =
        tapered_blend_surface(&wedge(theta, radius), &wedge(theta, radius), 1.0).expect("blend");
    let weights = surface.weights.clone().expect("rational");
    let centre = (radius / (theta / 2.0).sin(), 0.0);
    for vi in 0..=4 {
        let v = f64::from(vi) / 4.0;
        for ui in 0..=4 {
            let u = f64::from(ui) / 4.0;
            let point = evaluate(&surface.control_points, &weights, u, v);
            let r = (point.x - centre.0).hypot(point.y - centre.1);
            assert!(
                (r - radius).abs() < 1e-12,
                "constant taper must stay at radius {radius}, got {r}"
            );
        }
    }
}

#[test]
fn the_end_seam_is_the_start_rotated_by_the_sweep() {
    // An asymmetric corner: with a symmetric wedge a dropped rotation can
    // still land on the circle, so symmetry alone hides the bug.
    let centre = Point2::new(1.3, -0.4);
    let start = Point2::new(1.3 + 0.5, -0.4);
    let sweep = 1.1_f64;
    let bottom = BlendCorner {
        centre,
        start,
        end: Point2::new(0.0, 0.0),
        sweep,
    };
    let top = BlendCorner {
        centre: Point2::new(1.3, -0.4),
        start: Point2::new(1.3 + 0.9, -0.4),
        end: Point2::new(0.0, 0.0),
        sweep,
    };
    let surface = tapered_blend_surface(&bottom, &top, 1.0).expect("blend");
    let weights = surface.weights.clone().expect("rational");

    // Closed-form end point: start rotated about the centre by the sweep.
    let (sin, cos) = sweep.sin_cos();
    let radial = (start.x - centre.x, start.y - centre.y);
    let expected = (
        centre.x + radial.0 * cos - radial.1 * sin,
        centre.y + radial.0 * sin + radial.1 * cos,
    );
    let got = evaluate(&surface.control_points, &weights, 1.0, 0.0);
    let error = (got.x - expected.0).hypot(got.y - expected.1);
    assert!(
        error < 1e-12,
        "u=1 must be the start rotated by the sweep, off by {error:e}"
    );
}
