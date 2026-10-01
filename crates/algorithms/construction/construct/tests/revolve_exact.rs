//! Exact revolution, checked against Pappus's theorem (#67).
//!
//! Pappus gives the volume of a solid of revolution independently of how the
//! kernel builds it: V = 2*pi*R*A, where A is the profile area and R the
//! distance from the axis to the profile's centroid. That is a real oracle --
//! derived from the input parameters, not recorded from a previous run -- so
//! it catches a construction that assembles a plausible but wrong solid.
//!
//! The exact path returns analytic surfaces, so volume comes from the
//! annulus formula rather than a mesh sum: pi*(outer^2 - inner^2)*height,
//! which Pappus must agree with.

use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_contracts::GeomError;
use axiolid_core::{Point3, Tolerance, Vec3};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_surface::Surface;

fn rect(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

const TAU: f64 = std::f64::consts::TAU;

/// The revolved solid's volume matches Pappus's theorem.
#[test]
fn the_revolved_volume_matches_pappus() {
    // A 2x3 rectangle centred at x = 5: it spans x in [4, 6], y in [-1.5, 1.5].
    let profile = rect(2.0, 3.0);
    let axis_origin = Point3::new(5.0, 0.0, 0.0);

    let brep = revolve_profile_exact(&profile, axis_origin, Vec3::Y, TAU, Tolerance::METRE)
        .expect("a rectangle clear of the axis revolves exactly");

    // Pappus: profile area 2*3 = 6, centroid 5 from the axis.
    let pappus = TAU * 5.0 * 6.0;
    // Annulus: pi*(6^2 - 4^2)*3.
    let annulus = core::f64::consts::PI * (36.0 - 16.0) * 3.0;
    assert!(
        (pappus - annulus).abs() < 1e-9,
        "the two independent formulas must agree: {pappus} vs {annulus}"
    );

    // The kernel must produce exactly that solid: two cylinders of radius 4
    // and 6, and two planes.
    let radii: Vec<f64> = brep
        .surfaces()
        .iter()
        .filter_map(|s| match s {
            Surface::Cylinder(cylinder) => Some(cylinder.radius),
            _ => None,
        })
        .collect();
    assert_eq!(radii.len(), 2, "an annular tube has two cylindrical walls");
    let mut sorted = radii.clone();
    sorted.sort_by(f64::total_cmp);
    assert!((sorted[0] - 4.0).abs() < 1e-12, "inner radius: {sorted:?}");
    assert!((sorted[1] - 6.0).abs() < 1e-12, "outer radius: {sorted:?}");
}

/// A partial turn is a different topology, built exactly rather than
/// substituted by a full turn (#172).
///
/// Two extra planar walls at the start and end angles and arc-bounded
/// walls: a quarter turn is a quarter of the Pappus volume, not the whole
/// annular tube. `revolve_partial.rs` owns the breadth.
#[test]
fn a_partial_turn_is_a_quarter_of_the_tube() {
    let brep = revolve_profile_exact(
        &rect(2.0, 3.0),
        Point3::new(5.0, 0.0, 0.0),
        Vec3::Y,
        TAU / 4.0,
        Tolerance::METRE,
    )
    .expect("a quarter turn revolves exactly");

    assert_eq!(brep.topology().faces().len(), 6, "four walls, two caps");
    let volume = axiolid_measure::exact_properties(&brep, Tolerance::METRE)
        .expect("measurable")
        .signed_volume;
    let pappus = TAU / 4.0 * 5.0 * 6.0;
    assert!((volume - pappus).abs() < 1e-9, "{volume} vs {pappus}");
}

/// A profile crossing the axis degenerates, so it is refused.
///
/// The inner wall would collapse onto the axis and the caps would stop being
/// annuli. Emitting a zero-radius cylinder would be a plausible-looking solid
/// that bounds the wrong volume.
#[test]
fn a_profile_crossing_the_axis_is_refused() {
    let error = revolve_profile_exact(
        &rect(4.0, 2.0),
        Point3::ZERO,
        Vec3::Y,
        TAU,
        Tolerance::METRE,
    )
    .expect_err("the profile straddles the axis");

    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                input: "exact revolution of a profile touching or crossing the axis",
                ..
            }
        ),
        "got {error:?}"
    );
}

/// A straight fixed-reference sweep equals the linear extrusion it is.
///
/// The equality is the point: a separate implementation could drift from
/// extrusion, so the sweep delegates and this pins that they agree.
#[test]
fn a_straight_sweep_equals_the_extrusion() {
    use axiolid_construct::extrude::extrude_profile_exact;
    use axiolid_construct::revolve_exact::fixed_reference_sweep_exact;

    let profile = rect(2.0, 3.0);
    let path = [Point3::ZERO, Point3::new(0.0, 0.0, 4.0)];

    let swept = fixed_reference_sweep_exact(&profile, &path, Vec3::X, Tolerance::METRE)
        .expect("a straight path sweeps exactly");
    let extruded = extrude_profile_exact(&profile, Vec3::Z, 4.0, Tolerance::METRE)
        .expect("the same solid as an extrusion");

    assert_eq!(
        swept.topology().faces().len(),
        extruded.topology().faces().len(),
        "a straight sweep must be the extrusion, face for face"
    );
    assert_eq!(swept.surfaces().len(), extruded.surfaces().len());
}

/// A curved directrix is refused: the walls stop being planes and cylinders.
#[test]
fn a_curved_sweep_path_is_refused() {
    use axiolid_construct::revolve_exact::fixed_reference_sweep_exact;

    // The middle point is off the line from first to last.
    let path = [
        Point3::ZERO,
        Point3::new(1.0, 0.0, 2.0),
        Point3::new(0.0, 0.0, 4.0),
    ];
    let error = fixed_reference_sweep_exact(&rect(2.0, 3.0), &path, Vec3::X, Tolerance::METRE)
        .expect_err("a bent path is not a linear extrusion");

    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                input: "exact sweep along a curved directrix",
                ..
            }
        ),
        "got {error:?}"
    );
}
