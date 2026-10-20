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
use axiolid_curve::Curve2;
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

/// buildingSMART IFC4 Annex E `air-terminal-element.ifc` (#339):
/// `IFCRECTANGLEHOLLOWPROFILEDEF(.AREA.,$,$,24.,24.,2.,10.,0.)` in inches,
/// lowered to metres. The 10 in inner fillet is exactly half the 20 in
/// opening, a fully rounded (circular) core, but in metres it lands one ulp
/// above `x / 2 - thickness`.
#[test]
fn the_annex_e_tube_builds_with_a_fully_rounded_core() {
    let inch = 0.0254;
    let (x, t, ri) = (24.0 * inch, 2.0 * inch, 10.0 * inch);
    let profile = rect(x, x, Some(t), Some(0.0), Some(ri));
    let contour = rectangle_contour(&profile).expect("the Annex E profile lowers");
    // Fully rounded: the hole is four quarter arcs and no straight run.
    let hole = &contour.holes[0];
    assert_eq!(hole.segments.len(), 4, "{:?}", hole.segments);
    let core = x / 2.0 - t;
    for segment in &hole.segments {
        let Curve2::Circle(arc) = segment.curve else {
            panic!("a straight run in a fully rounded core: {segment:?}");
        };
        assert!((arc.radius - core).abs() < 1e-15, "radius {}", arc.radius);
        let centre = arc.frame.origin;
        assert!(centre.x.hypot(centre.y) < 1e-15, "centre {centre:?}");
    }
    let got = net_area(&profile);
    let want = x * x - PI * core * core;
    assert!((got - want).abs() < 1e-12, "net section {got}, want {want}");
    extrude_clean(profile);
}

/// The inner limit's rounding is that of the half-extent it is cut from,
/// not of the opening: 161 x 161 in with an 80 in wall leaves a 1 in core,
/// and its 0.5 in fillet lands about 100 epsilons of the core above it in
/// metres, 0.63 of the half-extent (#339).
#[test]
fn a_narrow_core_rounds_at_the_scale_of_its_half_extent() {
    let inch = 0.0254;
    let (x, t, ri) = (161.0 * inch, 80.0 * inch, 0.5 * inch);
    assert!(ri > x / 2.0 - t, "the case must round past its limit");
    let contour = rectangle_contour(&rect(x, x, Some(t), Some(0.0), Some(ri)))
        .expect("a fully rounded narrow core lowers");
    assert_eq!(contour.holes[0].segments.len(), 4);
}

/// The outer radius takes the same reach: a radius one ulp above the
/// half-extent is the stadium it states (#339).
#[test]
fn an_outer_radius_one_ulp_past_the_half_extent_is_the_stadium() {
    let r: f64 = 0.254;
    let half = f64::from_bits(r.to_bits() - 1);
    let (x, y) = (4.0 * half, 2.0 * half);
    assert!(r > y / 2.0);
    let got = net_area(&rect(x, y, None, Some(r), None));
    let want = (x - y) * y + PI * half * half;
    assert!((got - want).abs() < 1e-12, "stadium {got}, want {want}");
    extrude_clean(rect(x, y, None, Some(r), None));
}

/// The reach is rounding only: a radius a micrometre past its limit is a
/// modelling error and is still refused, inner and outer alike (#339).
#[test]
fn a_radius_a_micrometre_past_its_limit_is_refused() {
    let inch = 0.0254;
    let (x, t) = (24.0 * inch, 2.0 * inch);
    let core = x / 2.0 - t;
    let cases = [
        ("inner", rect(x, x, Some(t), Some(0.0), Some(core + 1e-6))),
        ("outer", rect(x, x, None, Some(x / 2.0 + 1e-6), None)),
    ];
    for (which, profile) in cases {
        let result = rectangle_contour(&profile);
        assert!(
            matches!(&result, Err(GeomError::InvalidInput(message)) if message.contains(which)),
            "{which}: expected a refusal by name, got {result:?}"
        );
    }
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
    // radius r. The surfaces are checked here; `curved_measure.rs` checks
    // the volume against Pappus.
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
fn a_hollow_rounded_rectangle_revolves_with_its_opening_as_a_void() {
    // The opening sweeps a closed ring-shaped cavity (#111). Pappus: the
    // section is symmetric about its centre, 5 from the axis.
    let (x, y, t, r) = (2.0, 3.0, 0.2, 0.3);
    let profile = Profile::Rectangle(rect(x, y, Some(t), Some(r), None));
    let solid = revolve_profile_exact(
        &profile,
        Point3::new(-5.0, 0.0, 0.0),
        Vec3::Y,
        TAU,
        Tolerance::METRE,
    )
    .expect("a hollow rectangle revolves");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    assert_eq!(solid.topology().solids()[0].voids.len(), 1, "one cavity");
    let section = x * y - (4.0 - core::f64::consts::PI) * r * r - (x - 2.0 * t) * (y - 2.0 * t);
    let expected = TAU * 5.0 * section;
    let volume = axiolid_measure::exact_properties(&solid, Tolerance::METRE)
        .expect("measurable")
        .signed_volume;
    assert!(
        (volume - expected).abs() < 1e-10 * expected,
        "expected {expected}, got {volume}"
    );
}
