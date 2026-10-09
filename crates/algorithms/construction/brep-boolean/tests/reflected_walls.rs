//! Booleans cut by a reflected curved wall (#288).
//!
//! An extrusion along `-z` is the reflection `z -> -z` of the prism along
//! `+z` (#275, `ExactBRep::transformed`). Its cylinder walls must cut a
//! box exactly as the same walls built directly along `+z` do: same
//! closed-form volume, clean audit, closed two-manifold.

use std::f64::consts::{FRAC_PI_2, PI};

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::boolean;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{BooleanOperator, Frame2, Interval, Point2, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::exact_properties;
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile,
};

fn tol() -> Tolerance {
    Tolerance::METRE
}

/// Wall length, thickness and height, and where the cut starts.
const L: f64 = 6.0;
const T: f64 = 0.3;
const H: f64 = 3.0;
const ZC: f64 = 2.0;
/// The arc radius.
const R: f64 = 1.2;

fn line(a: Point2, b: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: a,
            direction: b - a,
        }),
        domain: Interval::new(0.0, 1.0),
        same_sense: true,
    }
}

fn arc(centre: Point2, radius: f64, (from, to): (f64, f64)) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: centre,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius,
        }),
        domain: Interval::new(from.min(to), from.max(to)),
        same_sense: from < to,
    }
}

fn contour(segments: Vec<ProfileSegment>) -> Profile {
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

/// The rectangle `[0, 4] x [cy - 1, cy + 2]` less the disk of radius `R`
/// about `(0, cy)`, the arc split where it crosses `y = cy`, clockwise.
fn bitten(cy: f64) -> Profile {
    let c = Point2::new(0.0, cy);
    let foot = (R * R - 1.0).sqrt();
    let below = -(1.0 / R).asin();
    contour(vec![
        line(Point2::new(foot, cy - 1.0), Point2::new(4.0, cy - 1.0)),
        line(Point2::new(4.0, cy - 1.0), Point2::new(4.0, cy + 2.0)),
        line(Point2::new(4.0, cy + 2.0), Point2::new(0.0, cy + 2.0)),
        line(Point2::new(0.0, cy + 2.0), Point2::new(0.0, cy + R)),
        arc(c, R, (FRAC_PI_2, 0.0)),
        arc(c, R, (0.0, below)),
    ])
}

/// The rectangle `[-4, 0] x [cy - 2, cy + 2]` with the right half of the
/// disk added, the arc split at `y = cy`, counter-clockwise.
fn bulging(cy: f64) -> Profile {
    let c = Point2::new(0.0, cy);
    contour(vec![
        line(Point2::new(-4.0, cy - 2.0), Point2::new(0.0, cy - 2.0)),
        line(Point2::new(0.0, cy - 2.0), Point2::new(0.0, cy - R)),
        arc(c, R, (-FRAC_PI_2, 0.0)),
        arc(c, R, (0.0, FRAC_PI_2)),
        line(Point2::new(0.0, cy + R), Point2::new(0.0, cy + 2.0)),
        line(Point2::new(0.0, cy + 2.0), Point2::new(-4.0, cy + 2.0)),
        line(Point2::new(-4.0, cy + 2.0), Point2::new(-4.0, cy - 2.0)),
    ])
}

/// `integral of sqrt(R^2 - u^2) du` over the strip `|y| <= T/2`, `u = y - cy`.
fn disk_strip(cy: f64) -> f64 {
    let f = |u: f64| 0.5 * u * (R * R - u * u).sqrt() + 0.5 * R * R * (u / R).asin();
    let t = T / 2.0;
    f(t - cy) - f(-t - cy)
}

/// The wall's plan area each boundary covers.
fn covered(profile: &str, cy: f64) -> f64 {
    match profile {
        "bitten" => L / 2.0 * T - disk_strip(cy),
        "bulging" => L / 2.0 * T + disk_strip(cy),
        _ => unreachable!(),
    }
}

fn boundary(profile: &str, cy: f64) -> Profile {
    match profile {
        "bitten" => bitten(cy),
        "bulging" => bulging(cy),
        _ => unreachable!(),
    }
}

fn rect(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// The wall `[-L/2, L/2] x [-T/2, T/2] x [0, H]`.
fn wall() -> ExactBRep {
    extrude_profile_exact(&rect(L, T), Vec3::Z, H, tol()).expect("a wall")
}

/// `profile` swept over `z in [ZC, top]`: down from `top` (a reflected
/// prism) or up from `ZC` (built directly), then moved by `placement`.
fn tool(profile: &Profile, top: f64, down: bool, placement: Transform3) -> ExactBRep {
    let depth = top - ZC;
    let (direction, base) = if down { (-Vec3::Z, top) } else { (Vec3::Z, ZC) };
    extrude_profile_exact(profile, direction, depth, tol())
        .expect("a prism")
        .transformed(&(placement * Transform3::from_translation(Vec3::new(0.0, 0.0, base))))
        .expect("rigid")
}

fn run(a: &ExactBRep, b: &ExactBRep, op: BooleanOperator, what: &str) -> ExactBRep {
    let result = boolean(a, b, op, tol()).unwrap_or_else(|e| panic!("{what} {op:?}: {e}"));
    let health = geometric_audit(&result, tol());
    assert!(
        health.is_consistent(),
        "{what} {op:?}: {:?}",
        health.defects()
    );
    let topology = axiolid_topology::audit_brep(result.topology());
    assert!(topology.is_closed_manifold(), "{what} {op:?}: {topology:?}");
    result
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, tol())
        .expect("measurable")
        .signed_volume
}

fn close(what: &str, got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-9 * want.abs().max(1.0),
        "{what}: expected {want}, got {got}"
    );
}

fn placements() -> [(&'static str, Transform3); 2] {
    [
        ("identity", Transform3::IDENTITY),
        (
            "general",
            Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
                * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7),
        ),
    ]
}

#[test]
fn a_box_minus_a_downward_arc_prism_matches_the_upward_one() {
    for profile in ["bitten", "bulging"] {
        for cy in [0.0, 0.1, 0.3, -0.2] {
            for (name, p) in placements() {
                let box_ = wall().transformed(&p).expect("rigid");
                let shape = boundary(profile, cy);
                let want = L * T * H - (H - ZC) * covered(profile, cy);
                for down in [false, true] {
                    let what = format!("{profile} cy {cy} {name} down {down}");
                    let cut = run(
                        &box_,
                        &tool(&shape, 5.0, down, p),
                        BooleanOperator::Difference,
                        &what,
                    );
                    close(&what, volume(&cut), want);
                }
            }
        }
    }
}

#[test]
fn a_box_minus_a_downward_circle_matches_the_upward_one() {
    let radius = 0.1;
    for (cx, cy) in [(0.0, 0.0), (1.0, 0.1), (-1.5, 0.15)] {
        for (name, p) in placements() {
            let box_ = wall().transformed(&p).expect("rigid");
            let shape = Profile::Derived {
                basis: Box::new(Profile::Circle(CircleProfile {
                    radius,
                    thickness: None,
                })),
                transform: axiolid_core::Transform2::from_translation(Vec2::new(cx, cy)),
            };
            // The disk's part inside the strip `|y| <= T/2`.
            let f = |u: f64| {
                let u = u.clamp(-radius, radius);
                u * (radius * radius - u * u).sqrt() + radius * radius * (u / radius).asin()
            };
            let inside = f(T / 2.0 - cy) - f(-T / 2.0 - cy);
            let want = L * T * H - (H - ZC) * inside;
            for down in [false, true] {
                let what = format!("circle ({cx}, {cy}) {name} down {down}");
                let cut = run(
                    &box_,
                    &tool(&shape, 5.0, down, p),
                    BooleanOperator::Difference,
                    &what,
                );
                close(&what, volume(&cut), want);
            }
        }
    }
}

#[test]
fn a_reflected_cylinder_walled_prism_against_a_box_for_every_operator() {
    // A full disk standing through the box's top: the box is
    // `[-1, 1]^2 x [0, 1]`, the tool the disk of radius 0.5 about the
    // origin over `z in [0.5, 1.5]`, swept down from 1.5.
    let block = extrude_profile_exact(&rect(2.0, 2.0), Vec3::Z, 1.0, tol()).expect("a box");
    let r = 0.5;
    let disk = Profile::Circle(CircleProfile {
        radius: r,
        thickness: None,
    });
    let bite = boundary("bitten", 0.1);
    let cases = [
        ("disk", disk, PI * r * r),
        // The bitten profile clipped to the box `[-1, 1]^2`: `x in [0, 1]`,
        // `y in [-0.9, 1]`, outside the disk of radius `R` about `(0, 0.1)`.
        ("bitten", bite, bite_area_in_unit_box(0.1)),
    ];
    for (shape_name, shape, area_in_box) in cases {
        let block_area = 4.0;
        let shape_area = exact_properties(
            &extrude_profile_exact(&shape, Vec3::Z, 1.0, tol()).expect("a prism"),
            tol(),
        )
        .expect("measurable")
        .signed_volume;
        for (name, p) in placements() {
            let block = block.transformed(&p).expect("rigid");
            for down in [false, true] {
                let what = format!("{shape_name} {name} down {down}");
                let depth = 1.0;
                let base = if down { 1.5 } else { 0.5 };
                let direction = if down { -Vec3::Z } else { Vec3::Z };
                let prism = extrude_profile_exact(&shape, direction, depth, tol())
                    .expect("a prism")
                    .transformed(&(p * Transform3::from_translation(Vec3::new(0.0, 0.0, base))))
                    .expect("rigid");
                let common = 0.5 * area_in_box;
                for (op, want) in [
                    (BooleanOperator::Difference, block_area - common),
                    (BooleanOperator::Intersection, common),
                    (BooleanOperator::Union, block_area + shape_area - common),
                ] {
                    let result = run(&block, &prism, op, &what);
                    close(&format!("{what} {op:?}"), volume(&result), want);
                }
                // The tool minus the box too: the reflected solid kept.
                let result = run(&prism, &block, BooleanOperator::Difference, &what);
                close(
                    &format!("{what} reversed"),
                    volume(&result),
                    shape_area - common,
                );
            }
        }
    }
}

/// The area of `bitten(cy)` inside `[-1, 1]^2`: `[0, 1] x [cy - 1, 1]`
/// less the disk of radius `R` about `(0, cy)`.
fn bite_area_in_unit_box(cy: f64) -> f64 {
    // Over `x in [0, 1]`, the disk covers `y in [cy - s, cy + s]` with
    // `s = sqrt(R^2 - x^2)`, clipped to `[cy - 1, 1]`.
    let (lo, hi) = (cy - 1.0, 1.0);
    let rect = (hi - lo) * 1.0;
    // `s(x) >= 1` on `x <= sqrt(R^2 - 1)`: the disk spans the bottom line.
    let g = |x: f64| 0.5 * x * (R * R - x * x).sqrt() + 0.5 * R * R * (x / R).asin();
    let foot = (R * R - 1.0).sqrt();
    let upper_cap = hi - cy; // `0.9`: the disk reaches past `y = 1` while `s > 0.9`.
    let top = (R * R - upper_cap * upper_cap).sqrt();
    // Covered length at `x`: `min(cy + s, hi) - max(cy - s, lo)`.
    // Below: `max(cy - s, lo)` is `lo` for `x <= foot`, else `cy - s`.
    // Above: `min(cy + s, hi)` is `hi` for `x <= top`, else `cy + s`.
    let below = foot * 1.0 + (g(1.0) - g(foot));
    let above = top * upper_cap + (g(1.0) - g(top));
    rect - (below + above)
}

#[test]
#[ignore = "TODO(#287): a tilted elliptical-cylinder wall is not cut yet, either way up"]
fn a_box_minus_an_oblique_downward_arc_prism_has_the_closed_form_volume() {
    // Swept down from `z = 5` along `(a, 0, -1)`: at height `z` the profile
    // is moved `s = a (5 - z)` along `x`, so the arc's cut over the wall's
    // strip moves with it, `T s` of covered area per unit height, `s` at
    // its mean `a (5 - (ZC + H) / 2)` over the clipped height.
    let a = 0.2;
    let top = 5.0;
    let mean = a * (top - 0.5 * (ZC + H));
    for profile in ["bitten", "bulging"] {
        for cy in [0.0, 0.1, 0.3] {
            let shape = boundary(profile, cy);
            let moved = match profile {
                "bitten" => -T * mean,
                _ => T * mean,
            };
            let want = L * T * H - (H - ZC) * (covered(profile, cy) + moved);
            let direction = Vec3::new(a, 0.0, -1.0);
            let depth = (top - ZC) * direction.length();
            let what = format!("oblique {profile} cy {cy}");
            let prism = match extrude_profile_exact(&shape, direction, depth, tol()) {
                Ok(prism) => prism,
                Err(e) => panic!("{what}: {e}"),
            }
            .transformed(&Transform3::from_translation(Vec3::new(0.0, 0.0, top)))
            .expect("rigid");
            let cut = run(&wall(), &prism, BooleanOperator::Difference, &what);
            close(&what, volume(&cut), want);
        }
    }
}

// --- every curved family under a reflection --------------------------------

mod common;

/// The identity checks of all three operators, with the intersection's
/// closed form.
fn check(a: &ExactBRep, b: &ExactBRep, intersection: f64, what: &str) {
    let i = volume(&run(a, b, BooleanOperator::Intersection, what));
    let u = volume(&run(a, b, BooleanOperator::Union, what));
    let d = volume(&run(a, b, BooleanOperator::Difference, what));
    close(&format!("{what} intersection"), i, intersection);
    close(&format!("{what} identity"), u + i, volume(a) + volume(b));
    close(&format!("{what} difference"), d, volume(a) - i);
}

/// A reflection in the plane `x = 0`, then a turn about `axis`.
fn mirror(axis: Vec3, angle: f64) -> Transform3 {
    Transform3::from_axis_angle(axis.normalize(), angle)
        * Transform3::from_mat3(axiolid_core::Mat3::from_diagonal(Vec3::new(-1.0, 1.0, 1.0)))
}

fn block(x: (f64, f64), y: (f64, f64), z: (f64, f64)) -> ExactBRep {
    extrude_profile_exact(&rect(x.1 - x.0, y.1 - y.0), Vec3::Z, z.1 - z.0, tol())
        .expect("a box")
        .transformed(&Transform3::from_translation(Vec3::new(
            0.5 * (x.0 + x.1),
            0.5 * (y.0 + y.1),
            z.0,
        )))
        .expect("rigid")
}

fn revolve(segments: Vec<ProfileSegment>) -> ExactBRep {
    let profile = contour(segments);
    axiolid_construct::revolve_exact::revolve_profile_exact(
        &profile,
        axiolid_core::Point3::ZERO,
        Vec3::Y,
        std::f64::consts::TAU,
        tol(),
    )
    .expect("revolves")
}

fn reflected(solid: &ExactBRep, m: Transform3, what: &str) -> ExactBRep {
    let out = solid
        .transformed(&m)
        .unwrap_or_else(|e| panic!("{what}: {e}"));
    let health = geometric_audit(&out, tol());
    assert!(health.is_consistent(), "{what}: {:?}", health.defects());
    close(&format!("{what} volume"), volume(&out), volume(solid));
    out
}

#[test]
fn a_reflected_sphere_cut_by_a_box() {
    // Centred on the origin, every reflection keeps the ball: the box
    // beyond `x = d` holds the cap of height `r - d`.
    let r = 1.5;
    let ball = common::sphere(axiolid_core::Point3::ZERO, r);
    for m in [mirror(Vec3::Z, 0.0), mirror(Vec3::new(1.0, 2.0, 3.0), 0.7)] {
        let ball = reflected(&ball, m, "sphere");
        for d in [0.0, 0.4] {
            let h = r - d;
            let cap = PI * h * h * (3.0 * r - h) / 3.0;
            check(
                &ball,
                &block((d, 4.0), (-4.0, 4.0), (-4.0, 4.0)),
                cap,
                "sphere",
            );
        }
    }
}

#[test]
fn a_reflected_cone_and_torus_cut_by_a_box() {
    // Solids of revolution about `y`: a reflection in `x = 0` and a turn
    // about `y` keep them as sets.
    let ms = [mirror(Vec3::Y, 0.0), mirror(Vec3::Y, 0.3)];
    let p = Point2::new;
    let cone = revolve(vec![
        line(p(2.0, 0.0), p(4.0, 0.0)),
        line(p(4.0, 0.0), p(3.0, 2.0)),
        line(p(3.0, 2.0), p(2.0, 2.0)),
        line(p(2.0, 2.0), p(2.0, 0.0)),
    ]);
    // Frustum above y = 0.5: outer radius 3.75 -> 3, inner 2, height 1.5.
    let frustum = PI * 1.5 / 3.0 * (3.75 * 3.75 + 3.75 * 3.0 + 3.0 * 3.0) - PI * 4.0 * 1.5;
    let (big, r) = (4.0, 1.0);
    let c = p(big, 0.0);
    let q = FRAC_PI_2;
    let quarter = |from: f64| ProfileSegment {
        domain: Interval::new(from, from + q),
        same_sense: true,
        ..arc(c, r, (0.0, 1.0))
    };
    let torus = revolve(vec![
        quarter(0.0),
        quarter(q),
        quarter(PI),
        quarter(3.0 * q),
    ]);
    // Above y = 0.3: the tube's segment beyond that chord, by Pappus.
    let segment = r * r * (0.3_f64 / r).acos() - 0.3 * (r * r - 0.09).sqrt();
    for m in ms {
        let cone = reflected(&cone, m, "cone");
        let half = 0.5 * volume(&cone);
        check(
            &cone,
            &block((0.0, 6.0), (-6.0, 6.0), (-6.0, 6.0)),
            half,
            "cone half",
        );
        check(
            &cone,
            &block((-6.0, 6.0), (0.5, 6.0), (-6.0, 6.0)),
            frustum,
            "cone slab",
        );
        let torus = reflected(&torus, m, "torus");
        let half = PI * PI * big * r * r;
        check(
            &torus,
            &block((0.0, 6.0), (-6.0, 6.0), (-6.0, 6.0)),
            half,
            "torus half",
        );
        let slab = std::f64::consts::TAU * big * segment;
        check(
            &torus,
            &block((-6.0, 6.0), (0.3, 6.0), (-6.0, 6.0)),
            slab,
            "torus slab",
        );
    }
}

/// A disk of radius `r` swept down and sideways over `z in [-2, 0]`: an
/// oblique circular cylinder, the reflection of the one swept up.
fn oblique_disk(r: f64) -> (ExactBRep, Vec3) {
    let disk = Profile::Circle(CircleProfile {
        radius: r,
        thickness: None,
    });
    let direction = Vec3::new(0.3, -0.2, -1.0);
    let prism = extrude_profile_exact(&disk, direction, 2.0 * direction.length(), tol())
        .unwrap_or_else(|e| panic!("oblique disk: {e}"));
    assert!(prism
        .surfaces()
        .iter()
        .any(|s| matches!(s, axiolid_surface::Surface::EllipticalCylinder(_))));
    (prism, direction)
}

#[test]
#[ignore = "TODO(#287): the boolean takes no elliptical-cylinder face yet, either way up"]
fn a_reflected_elliptical_cylinder_cut_by_a_box() {
    let r = 0.5;
    let (prism, direction) = oblique_disk(r);
    close("prism", volume(&prism), PI * r * r * 2.0);
    // A box inside the prism about its axis at `z = -1`: no face is cut,
    // and every region of the reflected wall is classified.
    let h = 0.1;
    let c = direction;
    let inner = block((c.x - h, c.x + h), (c.y - h, c.y + h), (-1.0 - h, -1.0 + h));
    check(
        &prism,
        &inner,
        volume(&inner),
        "elliptical cylinder around a box",
    );
    // The slab `z in [-1.5, -0.5]` crosses the prism whole, so it holds
    // the disk's area times its thickness.
    let slab = block((-3.0, 3.0), (-3.0, 3.0), (-1.5, -0.5));
    check(
        &prism,
        &slab,
        PI * r * r,
        "elliptical cylinder cut by a slab",
    );
}
