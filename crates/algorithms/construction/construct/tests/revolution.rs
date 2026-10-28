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

// --- Profiles touching the axis (#363) ---

use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::audit_mesh;

const AXIS_CHORD: Scalar = 1e-3;

/// Assert the structural facts of a touching revolution: closed, consistently
/// wound, no degenerate triangle, no duplicate vertex and none unreferenced.
fn assert_sound(mesh: &TriMesh) {
    let health = audit_mesh(mesh, tol_for(1e-9));
    assert!(health.is_closed_two_manifold(), "{health:?}");
    let mut seen = std::collections::HashSet::new();
    for p in &mesh.positions {
        assert!(
            seen.insert([p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]),
            "duplicate vertex {p:?}"
        );
    }
    let mut used = vec![false; mesh.positions.len()];
    for i in &mesh.indices {
        used[*i as usize] = true;
    }
    assert!(used.iter().all(|u| *u), "an unreferenced vertex");
}

/// A closed full-turn wall of `S` stations over a four-point profile has
/// `4 S` vertices and `8 S` triangles: nothing was welded or collapsed.
fn assert_unwelded(mesh: &TriMesh) {
    assert_eq!(mesh.positions.len() * 2, mesh.indices.len() / 3);
}

/// Number of vertices on the axis x = z = 0.
fn on_axis(mesh: &TriMesh) -> usize {
    mesh.positions
        .iter()
        .filter(|p| p.x == 0.0 && p.z == 0.0)
        .count()
}

/// Pappus volume of the full turn about the Y axis of a polygon on one side:
/// `2*pi*|integral of x dA|`, with `integral x dA = sum (x0 + x1) cross / 6`.
fn pappus(ring: &[Point2]) -> Scalar {
    let mut moment = 0.0;
    for k in 0..ring.len() {
        let (p, q) = (ring[k], ring[(k + 1) % ring.len()]);
        moment += (p.x + q.x) * (p.x * q.y - q.x * p.y) / 6.0;
    }
    FULL * moment.abs()
}

fn revolve_touching(ring: &[Point2], angle: Scalar) -> TriMesh {
    turn_with(&rings(ring.to_vec(), vec![]), angle, AXIS_CHORD)
        .expect("a touching profile revolves")
}

fn turn_with(r: &Rings, angle: Scalar, chord: Scalar) -> Result<TriMesh, GeomError> {
    revolve(r, Point3::ZERO, Vec3::Y, angle, tol_for(chord))
}

/// The inscribed polygon of the right half of the unit disk, closed along
/// the axis.
fn half_disk(k: usize) -> Vec<Point2> {
    (0..=k)
        .map(|i| {
            let t = core::f64::consts::PI * i as Scalar / k as Scalar;
            Point2::new(t.sin(), -t.cos())
        })
        .map(|p| Point2::new(if p.x < 1e-12 { 0.0 } else { p.x }, p.y))
        .collect()
}

/// A right triangle with one leg on the axis.
fn cone() -> Vec<Point2> {
    vec![
        Point2::new(0.0, 0.0),
        Point2::new(2.0, 0.0),
        Point2::new(0.0, 1.0),
    ]
}

fn pie_cylinder() -> Vec<Point2> {
    square(0.0, 2.0, 0.0, 1.0)
}

/// Largest distance from the surface the rings sweep, sampled off the mesh's
/// own stations, to the mesh. Edges lying on the axis sweep nothing.
fn surface_deviation(ring: &[Point2], angle: Scalar, mesh: &TriMesh) -> Scalar {
    let triangles: Vec<[Point3; 3]> = mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| t.map(|i| mesh.positions[i as usize]))
        .collect();
    let mut worst: Scalar = 0.0;
    for k in 0..ring.len() {
        let (p, q) = (ring[k], ring[(k + 1) % ring.len()]);
        if p.x == 0.0 && q.x == 0.0 {
            continue;
        }
        for i in 0..4 {
            let x = (i as Scalar + 0.37) / 4.0;
            let (r, y) = (p.x + (q.x - p.x) * x, p.y + (q.y - p.y) * x);
            for j in 0..48 {
                let t = angle * (j as Scalar + 0.41) / 48.0;
                let point = Point3::new(r * t.cos(), y, -r * t.sin());
                let near = triangles
                    .iter()
                    .map(|tri| {
                        point.distance(closest_point_on_triangle(point, *tri).expect("closest"))
                    })
                    .fold(Scalar::INFINITY, Scalar::min);
                worst = worst.max(near);
            }
        }
    }
    worst
}

#[test]
fn a_half_disk_revolves_to_a_closed_sphere_like_mesh() {
    let ring = half_disk(16);
    let mesh = revolve_touching(&ring, FULL);
    assert_sound(&mesh);
    // The two poles are one vertex each, not one per station.
    assert_eq!(on_axis(&mesh), 2);
    let (want, got) = (pappus(&ring), volume(&mesh, tol_for(AXIS_CHORD)));
    assert!(got < want, "inscribed sweep: {got} vs {want}");
    assert!((want - got) / want < 3.0 * AXIS_CHORD, "{got} vs {want}");
    assert!(surface_deviation(&ring, FULL, &mesh) <= AXIS_CHORD);
}

#[test]
fn a_rectangle_on_the_axis_revolves_to_a_closed_cylinder() {
    let ring = pie_cylinder();
    let mesh = revolve_touching(&ring, FULL);
    assert_sound(&mesh);
    // The axis edge sweeps nothing: both ends are single vertices and no
    // wall triangle spans them.
    assert_eq!(on_axis(&mesh), 2);
    let (want, got) = (pappus(&ring), volume(&mesh, tol_for(AXIS_CHORD)));
    assert!(got < want && (want - got) / want < 3.0 * AXIS_CHORD);
    assert!(surface_deviation(&ring, FULL, &mesh) <= AXIS_CHORD);
}

#[test]
fn a_right_triangle_with_a_leg_on_the_axis_revolves_to_a_cone() {
    let ring = cone();
    let mesh = revolve_touching(&ring, FULL);
    assert_sound(&mesh);
    assert_eq!(on_axis(&mesh), 2);
    let want = core::f64::consts::PI * 4.0 / 3.0;
    assert!((pappus(&ring) - want).abs() < 1e-12);
    let got = volume(&mesh, tol_for(AXIS_CHORD));
    assert!(got < pappus(&ring) && (pappus(&ring) - got) / want < 3.0 * AXIS_CHORD);
    assert!(surface_deviation(&ring, FULL, &mesh) <= AXIS_CHORD);
}

#[test]
fn a_profile_touching_at_one_vertex_is_closed_and_meets_the_axis_once() {
    // The apex sits between the two rim points: a cylinder of radius 1 and
    // height 2 with a double cone taken out. Edge-manifold: the apex is one
    // shared vertex of a fan on each nappe.
    let ring = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, -1.0),
        Point2::new(1.0, 1.0),
    ];
    let mesh = revolve_touching(&ring, FULL);
    assert_sound(&mesh);
    assert_eq!(on_axis(&mesh), 1);
    let (want, got) = (pappus(&ring), volume(&mesh, tol_for(AXIS_CHORD)));
    assert!(got < want && (want - got) / want < 3.0 * AXIS_CHORD);
    assert!(surface_deviation(&ring, FULL, &mesh) <= AXIS_CHORD);
}

#[test]
fn partial_turns_of_touching_profiles_are_closed_and_share_the_axis_vertices() {
    for (ring, axis) in [(pie_cylinder(), 2), (half_disk(16), 2), (cone(), 2)] {
        for angle in [core::f64::consts::FRAC_PI_2, 1.0, 4.0] {
            let mesh = revolve_touching(&ring, angle);
            assert_sound(&mesh);
            // Start cap, wall and end cap meet at one vertex per axis point.
            assert_eq!(on_axis(&mesh), axis);
            let want = pappus(&ring) * angle / FULL;
            let got = volume(&mesh, tol_for(AXIS_CHORD));
            // The caps and the axis make the partial solid exact up to the
            // chord term of the arc.
            assert!(
                (want - got).abs() / want < 3.0 * AXIS_CHORD,
                "{angle}: {got} vs {want}"
            );
            assert!(surface_deviation(&ring, angle, &mesh) <= AXIS_CHORD);
        }
    }
}

#[test]
fn a_negative_angle_winds_a_touching_solid_outward_too() {
    for ring in [pie_cylinder(), half_disk(16), cone()] {
        for angle in [-FULL, -1.0] {
            let mesh = revolve_touching(&ring, angle);
            assert_sound(&mesh);
            let want = pappus(&ring) * angle.abs() / FULL;
            assert!((want - volume(&mesh, tol_for(AXIS_CHORD))).abs() / want < 3.0 * AXIS_CHORD);
        }
    }
}

#[test]
fn a_profile_with_a_hole_clear_of_the_axis_is_unaffected() {
    let r = rings(
        square(1.0, 6.0, 0.0, 4.0),
        vec![vec![
            Point2::new(1.5, 1.0),
            Point2::new(2.0, 1.0),
            Point2::new(2.0, 2.0),
        ]],
    );
    let mesh = turn(&r, FULL).expect("revolves");
    assert_sound(&mesh);
    assert_eq!(on_axis(&mesh), 0);
}

#[test]
fn a_hole_vertex_on_the_axis_follows_the_same_exact_rule() {
    // The hole pokes out of the outer ring, a malformed input; the exact
    // rule still welds its on-axis vertex to one mesh vertex.
    let r = rings(
        square(1.0, 6.0, 0.0, 4.0),
        vec![vec![
            Point2::new(0.0, 1.0),
            Point2::new(2.0, 1.0),
            Point2::new(2.0, 2.0),
        ]],
    );
    if let Ok(mesh) = turn(&r, FULL) {
        assert_eq!(on_axis(&mesh), 1);
    }
}

#[test]
fn touching_profiles_keep_the_budget_of_the_widest_radius() {
    // A point on the axis traces no circle, so the step count is that of
    // the same profile with the axis edge moved off by nothing: widest
    // radius governs. Halving the budget grows the stations, not the poles.
    let ring = half_disk(16);
    let coarse = turn_with(&rings(ring.clone(), vec![]), FULL, 1e-2).expect("coarse");
    let fine = turn_with(&rings(ring, vec![]), FULL, 1e-3).expect("fine");
    assert!(fine.positions.len() > coarse.positions.len());
    assert_eq!(on_axis(&coarse), 2);
    assert_eq!(on_axis(&fine), 2);
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
    // The vertex (1, 1) lies on the axis: meshed once, not refused (#363).
    let mesh = tilt(&touching).expect("a touching profile revolves");
    assert_sound(&mesh);
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
    let clear = at(0.5, &on_projection).expect("no touch");
    assert_sound(&clear);
    assert_unwelded(&clear);
    // Placement noise does not switch the check off or make it offset.
    assert!(matches!(
        at(1e-17, &across),
        Err(GeomError::InvalidInput(_))
    ));
    let touching = at(1e-17, &on_projection).expect("an on-axis edge is meshed");
    assert_sound(&touching);
    assert_eq!(on_axis(&touching), 2);
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
    let mesh = mesh.expect("revolves");
    assert_sound(&mesh);
    assert_unwelded(&mesh);
}
