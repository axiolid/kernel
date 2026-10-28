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
    assert_vertex_manifold(mesh);
}

/// The triangles around every vertex form a single fan: their link edges
/// make one cycle. The edge audit cannot see a pinch, where two fans meet.
fn assert_vertex_manifold(mesh: &TriMesh) {
    let mut link: Vec<Vec<(u32, u32)>> = vec![Vec::new(); mesh.positions.len()];
    for t in mesh.indices.as_chunks::<3>().0 {
        for k in 0..3 {
            link[t[k] as usize].push((t[(k + 1) % 3], t[(k + 2) % 3]));
        }
    }
    for (v, edges) in link.iter().enumerate() {
        let mut next = std::collections::HashMap::new();
        for (a, b) in edges {
            assert!(next.insert(*a, *b).is_none(), "vertex {v}: not a fan");
        }
        let mut at = edges[0].0;
        let mut walked = 0;
        loop {
            at = next[&at];
            walked += 1;
            if at == edges[0].0 {
                break;
            }
            assert!(walked <= edges.len(), "vertex {v}: not a cycle");
        }
        assert_eq!(walked, edges.len(), "vertex {v}: pinched into several fans");
    }
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
fn a_profile_touching_at_one_vertex_is_refused_as_a_pinch_in_a_full_turn() {
    // The apex sits between the two rim points: the two nappes of a double
    // cone meet there, so the boundary is not a two-manifold.
    let r = rings(
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, -1.0),
            Point2::new(1.0, 1.0),
        ],
        vec![],
    );
    for angle in [FULL, -FULL] {
        assert!(matches!(
            turn(&r, angle),
            Err(GeomError::UnsupportedInput { .. })
        ));
    }
}

#[test]
fn a_profile_touching_at_one_vertex_is_a_sector_point_in_a_partial_turn() {
    let ring = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, -1.0),
        Point2::new(1.0, 1.0),
    ];
    for angle in [1.0, -1.0, 4.0] {
        let mesh = revolve_touching(&ring, angle);
        assert_sound(&mesh);
        assert_eq!(on_axis(&mesh), 1);
        let want = pappus(&ring) * angle.abs() / FULL;
        let got = volume(&mesh, tol_for(AXIS_CHORD));
        assert!(
            (want - got).abs() / want < 3.0 * AXIS_CHORD,
            "{got} vs {want}"
        );
        assert!(surface_deviation(&ring, angle, &mesh) <= AXIS_CHORD);
    }
}

#[test]
fn an_isolated_touch_on_a_later_vertex_or_a_hole_is_refused_in_a_full_turn() {
    let outer = rings(
        vec![
            Point2::new(1.0, -1.0),
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 0.0),
        ],
        vec![],
    );
    assert!(matches!(
        turn(&outer, FULL),
        Err(GeomError::UnsupportedInput { .. })
    ));
    let hole = rings(
        square(1.0, 6.0, 0.0, 4.0),
        vec![vec![
            Point2::new(0.0, 1.0),
            Point2::new(2.0, 1.0),
            Point2::new(2.0, 2.0),
        ]],
    );
    assert!(matches!(
        turn(&hole, FULL),
        Err(GeomError::UnsupportedInput { .. })
    ));
}

#[test]
fn consecutive_axis_vertices_leave_no_unreferenced_vertex() {
    // Three consecutive vertices on the axis: the middle one has no wall
    // quad and, in a full turn, no cap.
    let ring = vec![
        Point2::new(0.0, 0.0),
        Point2::new(2.0, 0.0),
        Point2::new(2.0, 2.0),
        Point2::new(0.0, 2.0),
        Point2::new(0.0, 1.0),
    ];
    for angle in [FULL, -FULL] {
        let mesh = revolve_touching(&ring, angle);
        assert_sound(&mesh);
        assert_eq!(on_axis(&mesh), 2);
        let want = pappus(&ring);
        assert!((want - volume(&mesh, tol_for(AXIS_CHORD))).abs() / want < 3.0 * AXIS_CHORD);
    }
}

#[test]
fn a_c_profile_with_two_axis_segments_gives_a_shell_around_a_void() {
    let ring = vec![
        Point2::new(0.0, 0.0),
        Point2::new(2.0, 0.0),
        Point2::new(2.0, 4.0),
        Point2::new(0.0, 4.0),
        Point2::new(0.0, 3.0),
        Point2::new(1.0, 3.0),
        Point2::new(1.0, 1.0),
        Point2::new(0.0, 1.0),
    ];
    let mesh = revolve_touching(&ring, FULL);
    assert_sound(&mesh);
    assert_eq!(on_axis(&mesh), 4);
    let want = pappus(&ring);
    let got = volume(&mesh, tol_for(AXIS_CHORD));
    assert!(
        (want - got).abs() / want < 3.0 * AXIS_CHORD,
        "{got} vs {want}"
    );
}

#[test]
fn z_noise_in_the_axis_refuses_a_touching_profile() {
    // Within the tolerance the axis counts as in the plane, but a vertex on
    // its projection is then off the real axis: welding it would collapse a
    // real orbit, so the touch is refused instead (#363).
    let r = rings(
        vec![
            Point2::new(0.0, 100.0),
            Point2::new(1.0, 100.0),
            Point2::new(0.0, 101.0),
        ],
        vec![],
    );
    let noise = 0.999e-3;
    let tilted = revolve(
        &r,
        Point3::ZERO,
        Vec3::new(0.0, 1.0, noise),
        FULL,
        tol_for(1e-3),
    );
    let raised = revolve(
        &r,
        Point3::new(0.0, 0.0, noise),
        Vec3::Y,
        FULL,
        tol_for(1e-3),
    );
    for result in [tilted, raised] {
        assert!(matches!(result, Err(GeomError::UnsupportedInput { .. })));
    }
    // A crossing is still named as one under the same noise.
    let across = rings(square(-1.0, 2.0, 100.0, 101.0), vec![]);
    assert!(matches!(
        revolve(
            &across,
            Point3::new(0.0, 0.0, noise),
            Vec3::Y,
            FULL,
            tol_for(1e-3)
        ),
        Err(GeomError::InvalidInput(_))
    ));
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
            Point2::new(2.0, 2.0),
            Point2::new(3.0, 1.0),
        ],
        vec![],
    );
    // The edge from (1, 1) to (2, 2) lies on the axis: meshed, not refused
    // (#363).
    let mesh = tilt(&touching).expect("a touching profile revolves");
    assert_sound(&mesh);
    // The shared vertices are the profile points themselves, not rotated
    // copies of them.
    for p in [Point3::new(1.0, 1.0, 0.0), Point3::new(2.0, 2.0, 0.0)] {
        assert!(mesh.positions.contains(&p), "{p:?}");
    }
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
    // ... but z noise makes the on-axis edge a real thin tube, which is not
    // welded: the touch is refused (#363).
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
    let mesh = mesh.expect("revolves");
    assert_sound(&mesh);
    assert_unwelded(&mesh);
}

// An arc ending on the axis (#365): its end is the joint with the next
// segment, not a trigonometric sample, so it lies on the axis to the bit.

fn arc_segment(radius: Scalar, from: Scalar, to: Scalar) -> axiolid_profile::ProfileSegment {
    use axiolid_core::{Frame2, Interval, Vec2};
    axiolid_profile::ProfileSegment {
        curve: axiolid_curve::Curve2::Circle(axiolid_curve::Circle2 {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius,
        }),
        domain: Interval::new(from, to),
        same_sense: true,
    }
}

fn line_segment(from: Point2, to: Point2) -> axiolid_profile::ProfileSegment {
    axiolid_profile::ProfileSegment {
        curve: axiolid_curve::Curve2::Line(axiolid_curve::Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: axiolid_core::Interval::UNIT,
        same_sense: true,
    }
}

fn contour_of(segments: Vec<axiolid_profile::ProfileSegment>) -> Profile {
    Profile::Contour(axiolid_profile::ContourProfile {
        outer: axiolid_profile::Contour::new(segments),
        holes: Vec::new(),
    })
}

/// The half-disc of radius 0.5: an arc from (0,-0.5) over (0.5,0) to
/// (0,0.5), and a line back along x = 0.
fn arc_half_disc() -> Profile {
    let half = core::f64::consts::FRAC_PI_2;
    contour_of(vec![
        arc_segment(0.5, -half, half),
        line_segment(Point2::new(0.0, 0.5), Point2::new(0.0, -0.5)),
    ])
}

fn assert_no_degenerate_triangle(mesh: &TriMesh) {
    for t in mesh.indices.as_chunks::<3>().0 {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        assert!(
            (b - a).cross(c - a).length() > 1e-12,
            "degenerate triangle {a:?} {b:?} {c:?}"
        );
    }
}

#[test]
fn an_arc_ending_on_the_axis_ends_exactly_on_it() {
    let rings = profile_rings(&arc_half_disc(), AXIS_CHORD, tol_for(AXIS_CHORD)).expect("rings");
    let ring = &rings.outer;
    let on_axis: Vec<_> = ring.iter().filter(|p| p.x == 0.0).collect();
    assert_eq!(on_axis.len(), 2, "{ring:?}");
    assert!(on_axis.contains(&&Point2::new(0.0, 0.5)));
    assert!(on_axis.contains(&&Point2::new(0.0, -0.5)));
    assert!(ring.iter().all(|p| p.x >= 0.0));
}

#[test]
fn an_arc_half_disc_revolves_to_a_closed_sphere() {
    let rings = profile_rings(&arc_half_disc(), AXIS_CHORD, tol_for(AXIS_CHORD)).expect("rings");
    let mesh = turn_with(&rings, FULL, AXIS_CHORD).expect("a full turn");
    assert_sound(&mesh);
    assert_no_degenerate_triangle(&mesh);
    assert_eq!(on_axis(&mesh), 2);
    let want = 4.0 / 3.0 * core::f64::consts::PI * 0.5_f64.powi(3);
    let got = volume(&mesh, tol_for(AXIS_CHORD));
    assert!(
        got < want && (want - got) / want < 3.0 * AXIS_CHORD / 0.5,
        "{got} vs {want}"
    );
}

#[test]
fn half_and_quarter_turns_of_an_arc_half_disc_are_sound() {
    let rings = profile_rings(&arc_half_disc(), AXIS_CHORD, tol_for(AXIS_CHORD)).expect("rings");
    let sphere = 4.0 / 3.0 * core::f64::consts::PI * 0.5_f64.powi(3);
    for (angle, share) in [(FULL / 2.0, 0.5), (FULL / 4.0, 0.25)] {
        let mesh = turn_with(&rings, angle, AXIS_CHORD).expect("a partial turn");
        assert_sound(&mesh);
        assert_no_degenerate_triangle(&mesh);
        assert_eq!(on_axis(&mesh), 2);
        let (want, got) = (sphere * share, volume(&mesh, tol_for(AXIS_CHORD)));
        assert!(
            (want - got).abs() / want < 3.0 * AXIS_CHORD / 0.5,
            "{angle}: {got} vs {want}"
        );
    }
}

#[test]
fn a_quarter_disc_revolves_to_a_hemisphere() {
    let half = core::f64::consts::FRAC_PI_2;
    let profile = contour_of(vec![
        arc_segment(0.5, 0.0, half),
        line_segment(Point2::new(0.0, 0.5), Point2::ZERO),
        line_segment(Point2::ZERO, Point2::new(0.5, 0.0)),
    ]);
    let rings = profile_rings(&profile, AXIS_CHORD, tol_for(AXIS_CHORD)).expect("rings");
    let mesh = turn_with(&rings, FULL, AXIS_CHORD).expect("a full turn");
    assert_sound(&mesh);
    assert_no_degenerate_triangle(&mesh);
    let want = 2.0 / 3.0 * core::f64::consts::PI * 0.5_f64.powi(3);
    let got = volume(&mesh, tol_for(AXIS_CHORD));
    assert!(
        got < want && (want - got) / want < 3.0 * AXIS_CHORD / 0.5,
        "{got} vs {want}"
    );
}

// Joints rank by how exactly each side states the point (#365): a line's
// origin, a polyline vertex, a line's far end, a curved end.

const AXIS_X: Scalar = 0.3;

fn axis_x_vertices(mesh: &TriMesh) -> usize {
    mesh.positions
        .iter()
        .filter(|p| p.x == AXIS_X && p.z == 0.0)
        .count()
}

/// The profile revolved a full turn about the Y axis through `AXIS_X`.
fn turn_about_offset_axis(profile: &Profile) -> TriMesh {
    let rings = profile_rings(profile, AXIS_CHORD, tol_for(AXIS_CHORD)).expect("rings");
    revolve(
        &rings,
        Point3::new(AXIS_X, 0.0, 0.0),
        Vec3::Y,
        FULL,
        tol_for(AXIS_CHORD),
    )
    .expect("a profile touching the axis at its ends revolves")
}

/// An arc through (0.3, -0.5) and (1.1, 0), centre (0.45, 0.15), swept
/// counter-clockwise: the angles are rounded, so its ends are inexact.
fn offset_arc() -> axiolid_profile::ProfileSegment {
    use axiolid_core::{Frame2, Interval, Vec2};
    let centre = Point2::new(0.45, 0.15);
    let from = Point2::new(0.3, -0.5) - centre;
    let to = Point2::new(1.1, 0.0) - centre;
    axiolid_profile::ProfileSegment {
        curve: axiolid_curve::Curve2::Circle(axiolid_curve::Circle2 {
            frame: Frame2 {
                origin: centre,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: from.length(),
        }),
        domain: Interval::new(from.y.atan2(from.x), to.y.atan2(to.x)),
        same_sense: true,
    }
}

fn pole_points() -> (Point2, Point2, Point2) {
    (
        Point2::new(AXIS_X, 0.5),
        Point2::new(AXIS_X, -0.5),
        Point2::new(1.1, 0.0),
    )
}

#[test]
fn a_line_joint_takes_the_next_lines_origin_not_the_rounded_far_end() {
    let (top, bottom, right) = pole_points();
    for middle in [line_segment(bottom, right), offset_arc()] {
        let profile = contour_of(vec![
            line_segment(top, bottom),
            middle,
            line_segment(right, top),
        ]);
        let mesh = turn_about_offset_axis(&profile);
        assert_sound(&mesh);
        assert_no_degenerate_triangle(&mesh);
        assert_eq!(axis_x_vertices(&mesh), 2);
    }
}

#[test]
fn a_reversed_line_states_its_ends_by_the_oriented_parameters() {
    let (top, bottom, right) = pole_points();
    // The closing line stored from the top to the right and read backwards.
    let mut closing = line_segment(top, right);
    closing.same_sense = false;
    let profile = contour_of(vec![
        line_segment(top, bottom),
        line_segment(bottom, right),
        closing,
    ]);
    let mesh = turn_about_offset_axis(&profile);
    assert_sound(&mesh);
    assert_eq!(axis_x_vertices(&mesh), 2);
}

#[test]
fn a_polyline_vertex_beats_a_far_line_end_and_a_curved_end() {
    let (top, bottom, right) = pole_points();
    let corner = axiolid_profile::ProfileSegment {
        curve: axiolid_curve::Curve2::Polyline(axiolid_curve::Polyline2 {
            points: vec![top, bottom, right],
            closed: false,
        }),
        domain: axiolid_core::Interval::new(0.0, 2.0),
        same_sense: true,
    };
    let profile = contour_of(vec![corner, line_segment(right, top)]);
    let mesh = turn_about_offset_axis(&profile);
    assert_sound(&mesh);
    assert_eq!(axis_x_vertices(&mesh), 2);
}

#[test]
fn two_arcs_meeting_off_the_axis_still_close_a_ring() {
    let ring = profile_rings(
        &contour_of(vec![
            arc_segment(1.0, 0.0, core::f64::consts::PI),
            arc_segment(1.0, core::f64::consts::PI, FULL),
        ]),
        AXIS_CHORD,
        tol_for(AXIS_CHORD),
    )
    .expect("rings")
    .outer;
    for k in 0..ring.len() {
        assert_ne!(ring[k], ring[(k + 1) % ring.len()]);
    }
    let area: Scalar = (0..ring.len())
        .map(|k| {
            let (p, q) = (ring[k], ring[(k + 1) % ring.len()]);
            (p.x * q.y - q.x * p.y) / 2.0
        })
        .sum();
    assert!((area - core::f64::consts::PI).abs() < 1e-2, "{area}");
}

#[test]
fn a_hole_of_two_arcs_is_joined_like_an_outer_ring() {
    let square = |s: Scalar| {
        let c = [(-s, -s), (s, -s), (s, s), (-s, s)].map(|(x, y)| Point2::new(x, y));
        (0..4).map(|k| line_segment(c[k], c[(k + 1) % 4])).collect()
    };
    let profile = Profile::Contour(axiolid_profile::ContourProfile {
        outer: axiolid_profile::Contour::new(square(2.0)),
        holes: vec![axiolid_profile::Contour::new(vec![
            arc_segment(0.5, 0.0, core::f64::consts::PI),
            arc_segment(0.5, core::f64::consts::PI, FULL),
        ])],
    });
    let rings = profile_rings(&profile, AXIS_CHORD, tol_for(AXIS_CHORD)).expect("rings");
    assert_eq!(rings.holes.len(), 1);
    let hole = &rings.holes[0];
    for k in 0..hole.len() {
        assert_ne!(hole[k], hole[(k + 1) % hole.len()]);
    }
}

#[test]
fn the_snap_distance_joins_the_reported_deviation() {
    use axiolid_construct::profile::{profile_deviation, ProfileDeviation};
    // The second line starts 5e-4 above where the first ends; within the
    // 1e-3 linear tolerance they are one joint, and the ring moves by that.
    let profile = contour_of(vec![
        line_segment(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)),
        line_segment(Point2::new(1.0, 0.0005), Point2::new(1.0, 1.0)),
        line_segment(Point2::new(1.0, 1.0), Point2::new(0.0, 0.0)),
    ]);
    let tolerance = Tolerance::new(1e-3, 1e-9).expect("tolerance");
    match profile_deviation(&profile, 1e-3, tolerance).expect("deviation") {
        ProfileDeviation::Bounded(d) => assert!(d >= 5e-4, "{d}"),
        other => panic!("{other:?}"),
    }
}
