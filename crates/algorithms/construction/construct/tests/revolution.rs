//! Revolution against Pappus's theorem.
//!
//! A full revolution of a profile with area `A` whose centroid is `R` from
//! the axis has volume `2*pi*R*A`. That constant comes from outside this
//! crate, so it is a real oracle rather than the code agreeing with itself.

use axiolid_construct::profile::profile_rings;
use axiolid_construct::profile::Rings;
use axiolid_construct::revolve::revolve;
use axiolid_contracts::GeomError;
use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_profile::{Profile, RectangleProfile};

fn tol_for(chord: Scalar) -> Tolerance {
    Tolerance::new(chord, 1e-9).expect("tolerance")
}

fn volume(mesh: &TriMesh, tol: Tolerance) -> Scalar {
    volume_properties(mesh, tol)
        .expect("a revolved solid must be closed and two-manifold")
        .signed_volume
}

/// A rectangle offset from the axis: the textbook Pappus case, a torus of
/// rectangular section.
fn rect(x: Scalar, y: Scalar) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// Move every ring point away from the axis along +x.
///
/// The profile is authored in its own z = 0 plane. Revolving about the Y
/// axis sweeps the x offset into the major radius, so the offset is what
/// makes this a torus rather than a degenerate self-intersecting sweep.
fn shift(
    rings: &axiolid_construct::profile::Rings,
    dx: Scalar,
) -> axiolid_construct::profile::Rings {
    let f = |p: &axiolid_core::Point2| axiolid_core::Point2::new(p.x + dx, p.y);
    axiolid_construct::profile::Rings {
        outer: rings.outer.iter().map(f).collect(),
        holes: rings
            .holes
            .iter()
            .map(|h| h.iter().map(f).collect())
            .collect(),
    }
}

#[test]
fn a_full_revolution_matches_pappus() {
    // Rectangle 1 x 2 centred at x = 4, revolved about the Y axis through
    // the origin. Pappus: V = 2*pi*R*A = 2*pi*4*2.
    let chord = 1e-5;
    let rings = profile_rings(&rect(1.0, 2.0), chord, tol_for(chord)).expect("rings");
    let shifted = shift(&rings, 4.0);
    let mesh = revolve(
        &shifted,
        Point3::ZERO,
        Vec3::Y,
        core::f64::consts::TAU,
        tol_for(chord),
    )
    .expect("revolve");
    let want = core::f64::consts::TAU * 4.0 * 2.0;
    let got = volume(&mesh, tol_for(chord));
    // A revolved polygon is inscribed in the true torus, so it
    // under-estimates. Asserting the direction catches an outward-bulging
    // sweep that a symmetric bound would accept.
    assert!(
        got < want,
        "inscribed sweep must under-estimate: {got} vs {want}"
    );
    assert!((want - got) / want < 1e-4, "pappus volume {got} vs {want}");
}

#[test]
fn a_full_revolution_is_closed_and_seamless() {
    let chord = 1e-4;
    let rings = profile_rings(&rect(1.0, 2.0), chord, tol_for(chord)).expect("rings");
    let shifted = shift(&rings, 4.0);
    let mesh = revolve(
        &shifted,
        Point3::ZERO,
        Vec3::Y,
        core::f64::consts::TAU,
        tol_for(chord),
    )
    .expect("revolve");
    // volume_properties refuses anything not closed and two-manifold, so a
    // split seam cannot reach the assertion above.
    assert!(volume(&mesh, tol_for(chord)) > 0.0, "outward wound");
}

#[test]
fn a_half_revolution_is_half_the_volume() {
    // A partial turn is capped at both ends. Half the angle must give half
    // the solid, which also proves the caps are not double counted.
    let chord = 1e-5;
    let rings = profile_rings(&rect(1.0, 2.0), chord, tol_for(chord)).expect("rings");
    let shifted = shift(&rings, 4.0);
    let half = revolve(
        &shifted,
        Point3::ZERO,
        Vec3::Y,
        core::f64::consts::PI,
        tol_for(chord),
    )
    .expect("revolve");
    let want = core::f64::consts::TAU * 4.0 * 2.0 / 2.0;
    let got = volume(&half, tol_for(chord));
    assert!(
        (want - got) / want < 1e-4,
        "half revolution {got} vs {want}"
    );
}

#[test]
fn a_zero_angle_is_refused() {
    let chord = 1e-4;
    let rings = profile_rings(&rect(1.0, 2.0), chord, tol_for(chord)).expect("rings");
    let shifted = shift(&rings, 4.0);
    assert!(revolve(&shifted, Point3::ZERO, Vec3::Y, 0.0, tol_for(chord)).is_err());
    assert!(revolve(&shifted, Point3::ZERO, Vec3::ZERO, 1.0, tol_for(chord)).is_err());
}

// --- Which side of the axis the rings lie on (#351) ---

fn square(x0: Scalar, x1: Scalar, y0: Scalar, y1: Scalar) -> Vec<Point2> {
    vec![
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ]
}

fn rings(outer: Vec<Point2>, holes: Vec<Vec<Point2>>) -> Rings {
    Rings { outer, holes }
}

fn turn(r: &Rings, angle: Scalar) -> Result<TriMesh, GeomError> {
    revolve(r, Point3::ZERO, Vec3::Y, angle, tol_for(1e-3))
}

const FULL: Scalar = core::f64::consts::TAU;

#[test]
fn a_profile_crossing_the_axis_is_invalid_input() {
    let r = rings(square(-1.0, 2.0, 0.0, 1.0), vec![]);
    for angle in [FULL, 1.0, -1.0] {
        assert!(matches!(turn(&r, angle), Err(GeomError::InvalidInput(_))));
    }
}

#[test]
fn a_hole_crossing_the_axis_is_invalid_input() {
    // The outer ring is one-sided; only the hole reaches across. The hole
    // pokes out of the outer ring, which is a malformed input; the check
    // judges every ring's vertices and refuses it as crossing regardless.
    let r = rings(
        square(1.0, 6.0, 0.0, 4.0),
        vec![square(-0.5, 2.0, 1.0, 2.0)],
    );
    assert!(matches!(turn(&r, FULL), Err(GeomError::InvalidInput(_))));
}

#[test]
fn a_vertex_on_the_axis_is_unsupported() {
    let r = rings(
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 1.0),
        ],
        vec![],
    );
    assert!(matches!(
        turn(&r, FULL),
        Err(GeomError::UnsupportedInput { .. })
    ));
}

#[test]
fn an_edge_on_the_axis_is_unsupported() {
    let r = rings(square(0.0, 2.0, 0.0, 1.0), vec![]);
    for angle in [FULL, 1.0] {
        assert!(matches!(
            turn(&r, angle),
            Err(GeomError::UnsupportedInput { .. })
        ));
    }
}

#[test]
fn a_hole_touching_the_axis_is_unsupported() {
    let clear = rings(
        square(1.0, 6.0, 0.0, 4.0),
        vec![vec![
            Point2::new(1.5, 1.0),
            Point2::new(2.0, 1.0),
            Point2::new(2.0, 2.0),
        ]],
    );
    assert!(turn(&clear, FULL).is_ok());
    let touching = rings(
        square(1.0, 6.0, 0.0, 4.0),
        vec![vec![
            Point2::new(0.0, 1.0),
            Point2::new(2.0, 1.0),
            Point2::new(2.0, 2.0),
        ]],
    );
    assert!(matches!(
        turn(&touching, FULL),
        Err(GeomError::UnsupportedInput { .. })
    ));
}

#[test]
fn a_barely_crossing_vertex_is_decided_exactly() {
    // 1e-300 past the axis is still across it: no tolerance band.
    let r = rings(
        vec![
            Point2::new(-1e-300, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 1.0),
        ],
        vec![],
    );
    assert!(matches!(turn(&r, FULL), Err(GeomError::InvalidInput(_))));
}

#[test]
fn a_one_sided_profile_revolves_on_either_side_and_on_a_tilted_axis() {
    let right = rings(square(1.0, 2.0, 0.0, 1.0), vec![]);
    let left = rings(square(-2.0, -1.0, 0.0, 1.0), vec![]);
    for r in [&right, &left] {
        for angle in [FULL, 1.0, -1.0] {
            let mesh = turn(r, angle).expect("one-sided profile revolves");
            assert!(!mesh.indices.is_empty());
        }
    }
    // The axis y = x leaves the square (x in 2..3, y in 0..1) on one side.
    let tilted = revolve(
        &rings(square(2.0, 3.0, 0.0, 1.0), vec![]),
        Point3::ZERO,
        Vec3::new(1.0, 1.0, 0.0),
        FULL,
        tol_for(1e-3),
    );
    assert!(tilted.is_ok());
}

#[test]
fn a_tilted_axis_is_judged_by_its_own_line() {
    // Right of x = 0 throughout, yet the axis y = x runs through the square.
    let tilt = |r: &Rings| {
        revolve(
            r,
            Point3::ZERO,
            Vec3::new(1.0, 1.0, 0.0),
            FULL,
            tol_for(1e-3),
        )
    };
    let across = rings(square(2.0, 3.0, 0.0, 3.0), vec![]);
    assert!(matches!(tilt(&across), Err(GeomError::InvalidInput(_))));
    let touching = rings(
        vec![
            Point2::new(1.0, 1.0),
            Point2::new(3.0, 0.0),
            Point2::new(3.0, 1.0),
        ],
        vec![],
    );
    assert!(matches!(
        tilt(&touching),
        Err(GeomError::UnsupportedInput { .. })
    ));
}

#[test]
fn a_tilted_axis_with_a_large_offset_is_decided_without_rounding() {
    // `origin + direction` rounds in f64 here and tilts the line enough to
    // put a one-sided triangle across it. Checked with rationals: all three
    // vertices lie strictly right of the axis through (1000.1, 0) along
    // (0.6, 0.8).
    let r = rings(
        vec![
            Point2::new(1006.1, 8.0),
            Point2::new(1010.0, 0.0),
            Point2::new(1010.0, 8.0),
        ],
        vec![],
    );
    let mesh = revolve(
        &r,
        Point3::new(1000.1, 0.0, 0.0),
        Vec3::new(0.6, 0.8, 0.0),
        FULL,
        tol_for(1e-3),
    );
    assert!(mesh.is_ok(), "{mesh:?}");
}

#[test]
fn an_axis_parallel_to_the_plane_and_offset_is_checked_for_crossing_only() {
    let at =
        |z: Scalar, r: &Rings| revolve(r, Point3::new(0.0, 0.0, z), Vec3::Y, FULL, tol_for(1e-3));
    // Crossing is still refused above the plane.
    let across = rings(square(-1.0, 2.0, 0.0, 1.0), vec![]);
    assert!(matches!(at(0.5, &across), Err(GeomError::InvalidInput(_))));
    // A vertex on the projected axis is 0.5 from the real one: no touch.
    let on_projection = rings(square(0.0, 2.0, 0.0, 1.0), vec![]);
    assert!(at(0.5, &on_projection).is_ok());
    // Placement noise does not switch the check off or make it offset.
    assert!(matches!(
        at(1e-17, &across),
        Err(GeomError::InvalidInput(_))
    ));
    assert!(matches!(
        at(1e-17, &on_projection),
        Err(GeomError::UnsupportedInput { .. })
    ));
}

#[test]
fn an_axis_leaning_out_of_the_plane_is_not_judged_by_its_projection() {
    // The ring's edge x = 0 lies on the axis's projection into the plane,
    // but the axis itself passes 1 m or more clear of it. Whether such a
    // sweep crosses itself is not decided here (TODO(#362)).
    let r = rings(square(0.0, 2.0, 1.0, 2.0), vec![]);
    let mesh = revolve(
        &r,
        Point3::ZERO,
        Vec3::new(0.0, 1.0, 0.5),
        FULL,
        tol_for(1e-3),
    );
    assert!(mesh.is_ok(), "{mesh:?}");
}
