//! Exact rounded and hollow rectangles (kernel#111, `IfcRoundedRectangleProfileDef`
//! and `IfcRectangleHollowProfileDef`).
//!
//! Areas are checked against the closed form `x*y - (4 - pi)*r^2`, computed
//! independently of the kernel. A corner sampled into chords, dropped, or
//! rounded the wrong way changes that area; so does a hollow core that lost
//! its hole.

use axiolid_brep_audit::geometric_audit;
use axiolid_construct::contour_lower::{
    arc_ring_signed_area, contour_to_arc_ring, orient_arc_ring,
};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_construct::section_lower::rectangle_contour;
use axiolid_contracts::GeomError;
use axiolid_core::{Point3, Tolerance, Vec3};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_surface::Surface;

const PI: f64 = core::f64::consts::PI;
const TAU: f64 = core::f64::consts::TAU;

fn rect(
    x: f64,
    y: f64,
    thickness: Option<f64>,
    outer: Option<f64>,
    inner: Option<f64>,
) -> RectangleProfile {
    RectangleProfile {
        x,
        y,
        thickness,
        outer_radius: outer,
        inner_radius: inner,
    }
}

/// Net section area from the kernel's own exact rings: outer counter-clockwise,
/// each hole re-oriented clockwise, exactly as the extruder does it.
fn net_area(rectangle: &RectangleProfile) -> f64 {
    let contour = rectangle_contour(rectangle).expect("rectangle lowers");
    let tol = Tolerance::METRE;
    let outer = orient_arc_ring(&contour_to_arc_ring(&contour.outer, tol).unwrap(), true).unwrap();
    let mut area = arc_ring_signed_area(&outer);
    for hole in &contour.holes {
        let ring = orient_arc_ring(&contour_to_arc_ring(hole, tol).unwrap(), false).unwrap();
        area += arc_ring_signed_area(&ring);
    }
    area
}

fn rounded_area(x: f64, y: f64, r: f64) -> f64 {
    x * y - (4.0 - PI) * r * r
}

fn cylinder_radii(solid: &axiolid_brep::ExactBRep) -> Vec<f64> {
    solid
        .surfaces()
        .iter()
        .filter_map(|s| match s {
            Surface::Cylinder(c) => Some(c.radius),
            _ => None,
        })
        .collect()
}

fn extrude_clean(rectangle: RectangleProfile) -> axiolid_brep::ExactBRep {
    let solid = extrude_profile_exact(
        &Profile::Rectangle(rectangle),
        Vec3::Z,
        2.0,
        Tolerance::METRE,
    )
    .expect("a valid rectangle extrudes exactly");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "solid must audit clean: {:?}",
        health.defects()
    );
    solid
}

#[test]
fn a_rounded_rectangle_extrudes_with_four_cylinder_corners() {
    let (x, y, r) = (3.0, 2.0, 0.4);
    let solid = extrude_clean(rect(x, y, None, Some(r), None));
    let radii = cylinder_radii(&solid);
    assert_eq!(
        radii.len(),
        4,
        "one cylinder wall per corner, got {radii:?}"
    );
    for value in &radii {
        assert!((value - r).abs() < 1e-12, "corner radius {value}, want {r}");
    }
    let got = net_area(&rect(x, y, None, Some(r), None));
    let want = rounded_area(x, y, r);
    assert!(
        (got - want).abs() < 1e-12,
        "section area {got}, want {want}"
    );
}

#[test]
fn a_hollow_rectangle_with_both_radii_keeps_its_rounded_core() {
    // IfcRectangleHollowProfileDef: 0.4 x 0.3 tube, 0.02 wall, outer radius
    // 0.03, inner radius 0.01. A dropped hole or a chorded corner on either
    // ring shows up as an area error.
    let (x, y, t, ro, ri) = (0.4, 0.3, 0.02, 0.03, 0.01);
    let profile = rect(x, y, Some(t), Some(ro), Some(ri));
    let solid = extrude_clean(profile);
    let mut radii = cylinder_radii(&solid);
    radii.sort_by(f64::total_cmp);
    assert_eq!(
        radii.len(),
        8,
        "four outer and four inner corners, got {radii:?}"
    );
    assert!(
        radii[..4].iter().all(|v| (v - ri).abs() < 1e-12),
        "inner: {radii:?}"
    );
    assert!(
        radii[4..].iter().all(|v| (v - ro).abs() < 1e-12),
        "outer: {radii:?}"
    );

    let want = rounded_area(x, y, ro) - rounded_area(x - 2.0 * t, y - 2.0 * t, ri);
    let got = net_area(&profile);
    assert!((got - want).abs() < 1e-12, "net section {got}, want {want}");
}

#[test]
fn a_hollow_rectangle_with_sharp_corners_still_takes_the_dedicated_path() {
    // No radii: the existing planar-only builder stays in charge, so this
    // change cannot have altered the shape it already produced.
    let solid = extrude_clean(rect(2.0, 1.0, Some(0.1), None, None));
    assert!(
        cylinder_radii(&solid).is_empty(),
        "a sharp tube has no cylinder walls"
    );
}

#[test]
fn a_radius_filling_the_short_side_makes_a_stadium() {
    // r = y/2 consumes both short edges entirely: two semicircular ends and
    // two straight sides, the extreme the router must accept.
    let (x, y) = (4.0, 2.0);
    let r = y / 2.0;
    let got = net_area(&rect(x, y, None, Some(r), None));
    let want = (x - y) * y + PI * r * r;
    assert!(
        (got - want).abs() < 1e-12,
        "stadium area {got}, want {want}"
    );
    extrude_clean(rect(x, y, None, Some(r), None));
}

#[test]
fn a_zero_radius_is_the_sharp_rectangle() {
    let got = net_area(&rect(3.0, 2.0, None, Some(0.0), None));
    assert!((got - 6.0).abs() < 1e-12, "got {got}");
}

#[test]
fn invalid_radii_are_refused_rather_than_clamped() {
    let cases = [
        ("negative", rect(3.0, 2.0, None, Some(-0.1), None)),
        (
            "wider than the half extent",
            rect(3.0, 2.0, None, Some(1.01), None),
        ),
        ("non-finite", rect(3.0, 2.0, None, Some(f64::NAN), None)),
        (
            "inner radius on a filled rectangle",
            rect(3.0, 2.0, None, None, Some(0.1)),
        ),
        // Outer 0.5 on a 0.02 wall with a sharp core: the corner wall vanishes.
        (
            "corner wall vanishes",
            rect(2.0, 2.0, Some(0.02), Some(0.5), None),
        ),
    ];
    for (what, profile) in cases {
        let result =
            extrude_profile_exact(&Profile::Rectangle(profile), Vec3::Z, 1.0, Tolerance::METRE);
        assert!(
            matches!(
                result,
                Err(GeomError::InvalidInput(_) | GeomError::Degenerate(_))
            ),
            "{what}: expected a refusal, got {result:?}"
        );
    }
}

#[test]
fn a_rounded_rectangle_revolves_into_tori_at_the_arc_centres() {
    // Section centred 5 from the axis. Its rounded corners sit on the inner
    // and outer walls, so their arc centres are at 5 -+ (x/2 - r) from the
    // axis: each corner sweeps a torus with THAT major radius and minor
    // radius r. `exact_properties` refuses toroidal faces, so the check is
    // on the surfaces themselves rather than a volume that cannot be taken.
    let (x, y, r, centre) = (2.0, 3.0, 0.5, 5.0);
    let profile = Profile::Rectangle(rect(x, y, None, Some(r), None));
    let solid = revolve_profile_exact(
        &profile,
        Point3::new(-centre, 0.0, 0.0),
        Vec3::Y,
        TAU,
        Tolerance::METRE,
    )
    .expect("a rounded rectangle clear of the axis revolves");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "revolved solid must audit clean: {:?}",
        health.defects()
    );

    let mut tori: Vec<(f64, f64)> = solid
        .surfaces()
        .iter()
        .filter_map(|s| match s {
            Surface::Torus(t) => Some((t.major_radius, t.minor_radius)),
            _ => None,
        })
        .collect();
    tori.sort_by(|a, b| a.0.total_cmp(&b.0));
    assert_eq!(
        tori.len(),
        4,
        "four rounded corners sweep four tori, got {tori:?}"
    );
    let inner = centre - (x / 2.0 - r);
    let outer = centre + (x / 2.0 - r);
    for (index, (major, minor)) in tori.iter().enumerate() {
        let want = if index < 2 { inner } else { outer };
        assert!(
            (major - want).abs() < 1e-12 && (minor - r).abs() < 1e-12,
            "torus {index}: ({major}, {minor}), want ({want}, {r}); all {tori:?}"
        );
    }
}

#[test]
fn a_hollow_rectangle_revolution_still_refuses_by_name() {
    let profile = Profile::Rectangle(rect(2.0, 3.0, Some(0.2), Some(0.3), None));
    let result = revolve_profile_exact(
        &profile,
        Point3::new(-5.0, 0.0, 0.0),
        Vec3::Y,
        TAU,
        Tolerance::METRE,
    );
    assert!(
        matches!(&result, Err(GeomError::UnsupportedInput { input, .. }) if *input == "hollow rectangle exact revolution"),
        "got {result:?}"
    );
}
