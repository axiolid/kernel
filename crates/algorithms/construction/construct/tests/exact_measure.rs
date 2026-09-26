//! Exact B-rep mass properties, without tessellating (#72).
//!
//! The gap this closes: v0.6's revolution and chamfer tests hand-rolled
//! divergence sums over `ExactBRep` vertices because no provider could
//! measure an exact solid. `MeshMeasure` needs triangles; an exact B-rep has
//! none.
//!
//! # Why this test lives in `axiolid-construct`
//!
//! It needs an exact B-rep to measure, and the only way to build one is
//! `extrude_profile_exact`. `axiolid-construct` already depends on
//! `axiolid-measure`, so a dev-dependency the other way would close a
//! cycle. The consumer side is the correct home for the test.

use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{Tolerance, Vec3};
use axiolid_measure::{exact_properties, ExactMeasureError};
use axiolid_profile::{Profile, RectangleProfile};

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("tolerance")
}

fn rectangle(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// An exact prism measures at machine precision, not tessellation fidelity.
#[test]
fn an_exact_prism_matches_its_closed_form() {
    let (x, y, depth) = (3.0, 5.0, 7.0);
    let brep = extrude_profile_exact(&rectangle(x, y), Vec3::Z, depth, tol())
        .expect("a rectangle extrudes exactly");

    let props = exact_properties(&brep, tol()).expect("a closed prism is measurable");

    let expected_volume = x * y * depth;
    // SIGNED, not absolute: an outward-oriented solid must measure POSITIVE.
    // Comparing magnitudes would let a provider that ignores face
    // orientation pass, and orientation is exactly what distinguishes a
    // solid from its inside-out twin.
    assert!(
        (props.signed_volume - expected_volume).abs() < 1e-12 * expected_volume,
        "expected signed volume {expected_volume}, got {}",
        props.signed_volume
    );

    let expected_area = 2.0 * (x * y) + 2.0 * depth * (x + y);
    assert!(
        (props.area - expected_area).abs() < 1e-12 * expected_area,
        "expected area {expected_area}, got {}",
        props.area
    );
}

/// The centroid of a symmetric prism sits on its axis, at mid-height.
#[test]
fn the_centroid_is_where_symmetry_demands() {
    let depth = 4.0;
    let brep = extrude_profile_exact(&rectangle(2.0, 2.0), Vec3::Z, depth, tol())
        .expect("a rectangle extrudes exactly");
    let props = exact_properties(&brep, tol()).expect("measurable");

    assert!(
        props.centroid.x.abs() < 1e-12 && props.centroid.y.abs() < 1e-12,
        "a profile centred on the origin has a centroid on the axis: {:?}",
        props.centroid
    );
    assert!(
        (props.centroid.z - depth / 2.0).abs() < 1e-12,
        "expected mid-height {}, got {}",
        depth / 2.0,
        props.centroid.z
    );
}

/// A curved face is integrated over its own parameters, not refused and not
/// sampled (#125): a cylinder measures `pi r^2 h` to machine precision rather
/// than to the chord error of a tessellation. `curved_measure.rs` holds the
/// per-family oracles.
#[test]
fn a_curved_face_is_integrated_not_refused() {
    use axiolid_profile::CircleProfile;

    let circle = Profile::Circle(CircleProfile {
        radius: 1.0,
        thickness: None,
    });
    let brep =
        extrude_profile_exact(&circle, Vec3::Z, 2.0, tol()).expect("a circle extrudes exactly");

    let props = exact_properties(&brep, tol()).expect("a cylinder is measurable");
    let expected = core::f64::consts::PI * 2.0;
    assert!(
        (props.signed_volume - expected).abs() < 1e-12 * expected,
        "expected {expected}, got {}",
        props.signed_volume
    );
}

/// What cannot be integrated is still refused by name, pointing at the
/// approximate path.
#[test]
fn a_refusal_names_the_approximate_path() {
    let text = ExactMeasureError::NotConverged.to_string();
    assert!(
        text.contains("MeshMeasure"),
        "refusal should point at the approximate path, got: {text}"
    );
}

/// The exact and mesh paths agree on the same solid.
///
/// A differential check across two implementations that share no code: the
/// exact path integrates over B-rep boundary polygons, the mesh path over a
/// tessellation. For a planar-faced prism the tessellation is exact, so they
/// must agree to near machine precision rather than to mesh fidelity.
#[test]
fn the_exact_and_mesh_paths_agree() {
    use axiolid_construct::extrude::extrude_profile;
    use axiolid_construct::profile::profile_rings;
    use axiolid_measure::{Measure, MeshMeasure};

    let (x, y, depth) = (3.0, 5.0, 7.0);
    let profile = rectangle(x, y);

    let brep = extrude_profile_exact(&profile, Vec3::Z, depth, tol()).expect("exact extrusion");
    let exact = exact_properties(&brep, tol()).expect("measurable");

    let rings = profile_rings(&profile, 1e-5, tol()).expect("profile rings");
    let mesh = extrude_profile(&rings, Vec3::Z, depth, tol()).expect("mesh extrusion");
    let meshed = MeshMeasure.measure(&mesh, tol()).expect("measurable");

    assert!(
        (exact.signed_volume - meshed.signed_volume).abs() < 1e-9,
        "exact {} vs mesh {}",
        exact.signed_volume,
        meshed.signed_volume
    );
    assert!(
        (exact.area - meshed.area).abs() < 1e-9,
        "exact area {} vs mesh area {}",
        exact.area,
        meshed.area
    );
    assert!(
        (exact.centroid - meshed.centroid).length() < 1e-9,
        "exact centroid {:?} vs mesh {:?}",
        exact.centroid,
        meshed.centroid
    );
}

/// A solid that does not touch `z = 0` measures its true volume.
///
/// Every other prism here starts at `z = 0`, where the bottom cap adds
/// nothing to `int z n_z dA` whichever way it is wound -- so a measure that
/// ignored `Face::orientation` still got them right. This one sits at
/// `2 <= z <= 3`: its bottom cap contributes `-2 * area` only when the
/// `Reversed` face is honoured, and `+2 * area` (volume 7/3 for a unit
/// cube, not 1) when it is not.
#[test]
fn a_raised_solid_measures_the_same_as_one_on_the_ground() {
    use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
    use axiolid_core::{BooleanOperator, Point2};
    use axiolid_overlay::{ArcRing, ArcVertex};

    let square = |half: f64| {
        ArcRing::new(
            [(-half, -half), (half, -half), (half, half), (-half, half)]
                .into_iter()
                .map(|(x, y)| ArcVertex::straight(Point2::new(x, y)))
                .collect(),
        )
    };
    let block = |bottom: f64, top: f64| ArcPrism {
        section: square(0.5),
        bottom,
        top,
    };
    let cover = |bottom: f64, top: f64| ArcPrism {
        section: square(1.0),
        bottom,
        top,
    };
    let at = |bottom: f64| {
        let solid = boolean_arc_prisms_exact(
            &block(bottom, bottom + 1.0),
            &cover(bottom, bottom + 1.0),
            BooleanOperator::Intersection,
            tol(),
        )
        .expect("a unit cube");
        exact_properties(&solid, tol()).expect("all-planar")
    };
    let ground = at(0.0);
    let raised = at(2.0);
    assert!(
        (ground.signed_volume - 1.0).abs() < 1e-12,
        "{}",
        ground.signed_volume
    );
    assert!(
        (raised.signed_volume - 1.0).abs() < 1e-12,
        "{}",
        raised.signed_volume
    );
    // Translating a solid moves its centroid by the same amount.
    assert!(
        (raised.centroid.z - 2.5).abs() < 1e-12,
        "{:?}",
        raised.centroid
    );
    assert!(
        (ground.centroid.z - 0.5).abs() < 1e-12,
        "{:?}",
        ground.centroid
    );
}

/// A straight-edged planar face with a hole: summed as vectors, the hole's
/// opposite winding subtracts its area. Summing triangle magnitudes instead
/// added it (a 4 x 4 plate with a 2 x 2 hole read 20 per cap, not 12).
#[test]
fn a_planar_hole_subtracts_its_area() {
    use axiolid_core::{Interval, Point2};
    use axiolid_curve::{Curve2, Line2};
    use axiolid_profile::{Contour, ContourProfile, ProfileSegment};

    let square = |half: f64| {
        let corners = [
            Point2::new(-half, -half),
            Point2::new(half, -half),
            Point2::new(half, half),
            Point2::new(-half, half),
        ];
        Contour::new(
            (0..4)
                .map(|i| ProfileSegment {
                    curve: Curve2::Line(Line2 {
                        origin: corners[i],
                        direction: corners[(i + 1) % 4] - corners[i],
                    }),
                    domain: Interval::UNIT,
                    same_sense: true,
                })
                .collect(),
        )
    };
    let depth = 3.0;
    let profile = Profile::Contour(ContourProfile {
        outer: square(2.0),
        holes: vec![square(1.0)],
    });
    let brep = extrude_profile_exact(&profile, Vec3::Z, depth, tol()).expect("extrudes");
    let props = exact_properties(&brep, tol()).expect("measurable");
    let section = 16.0 - 4.0;
    assert!(
        (props.signed_volume - section * depth).abs() < 1e-12,
        "{}",
        props.signed_volume
    );
    let area = 2.0 * section + (16.0 + 8.0) * depth;
    assert!(
        (props.area - area).abs() < 1e-12,
        "expected {area}, got {}",
        props.area
    );
}
