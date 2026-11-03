//! The caps of `extrude_profile` avoid sliver triangles wherever a
//! triangulation avoids them, from every start vertex (#415).
//!
//! A sliver has a doubled area within `tolerance^2`, which `audit_mesh`
//! reads as degenerate; one that is also flat to rounding is placed onto a
//! line or inverted by a georeferenced placement. Each ring is extruded,
//! placed at a national-grid translation and turned, and audited closed with
//! no degenerate triangle. Ear clipping alone leaves a sliver as the last
//! triangle of a polygon cut down to a quadrilateral; the flip pass removes
//! it.

use axiolid_construct::extrude::extrude_profile;
use axiolid_construct::profile::Rings;
use axiolid_core::{Point2, Point3, Tolerance, Vec3};
use axiolid_mesh::{audit_mesh, TriMesh};

const X_AXIS: [f64; 2] = [0.999215271103513, 0.0396086100934191];
const SHIFT: [f64; 3] = [643692.128712773, 5649156.75663248, 213.66];
const TOL: Tolerance = Tolerance::MILLIMETRE;

fn rings(ring: &[[f64; 2]]) -> Rings {
    Rings {
        outer: ring.iter().map(|q| Point2::new(q[0], q[1])).collect(),
        holes: Vec::new(),
    }
}

fn placed(mesh: &TriMesh) -> TriMesh {
    let positions = mesh
        .positions
        .iter()
        .map(|p| {
            Point3::new(
                X_AXIS[0] * p.x - X_AXIS[1] * p.y + SHIFT[0],
                X_AXIS[1] * p.x + X_AXIS[0] * p.y + SHIFT[1],
                p.z + SHIFT[2],
            )
        })
        .collect();
    TriMesh::new(positions, mesh.indices.clone())
}

/// The ring extruded, placed and audited: closed, no degenerate triangle.
fn assert_closed(ring: &[[f64; 2]], what: &str) {
    let mesh = extrude_profile(&rings(ring), Vec3::Z, 3.6, TOL)
        .unwrap_or_else(|e| panic!("{what}: refused: {e}"));
    let health = audit_mesh(&placed(&mesh), TOL);
    assert!(
        health.is_closed_two_manifold() && health.degenerate_triangles == 0,
        "{what}: {health:?}"
    );
}

/// `ring` started at each of its vertices.
fn rotations(ring: &[[f64; 2]]) -> impl Iterator<Item = (usize, Vec<[f64; 2]>)> + '_ {
    (0..ring.len()).map(|k| {
        let mut r = ring.to_vec();
        r.rotate_left(k);
        (k, r)
    })
}

fn assert_every_start_closed(ring: &[[f64; 2]], what: &str) {
    for (k, r) in rotations(ring) {
        assert_closed(&r, &format!("{what}, starting at vertex {k}"));
    }
}

#[test]
fn a_triangular_room_with_a_base_vertex_off_the_side_closes() {
    // Cutting the corner at the start leaves the sliver as the last
    // triangle; the diagonal from the base vertex to the apex avoids it.
    assert_every_start_closed(
        &[[2.0, 3.0], [0.0, 0.0], [1.7, -2e-13], [4.0, 0.0]],
        "triangular room",
    );
}

#[test]
fn a_quad_with_a_sliver_wider_than_the_flat_band_closes() {
    // 5e-8 off the side: a sliver by the audit, not flat to rounding, so it
    // was meshed (open) and not refused.
    assert_every_start_closed(
        &[[1.0, 1.0], [0.0, 0.0], [1.0, -5e-8], [2.0, 0.0]],
        "quad, h = 5e-8",
    );
}

#[test]
fn an_l_shape_with_near_collinear_vertices_closes_from_every_start() {
    let l = [
        [0.0, 0.0],
        [3.0, -1e-13],
        [6.0, 0.0],
        [6.0, 2.0],
        [4.0, 2.0 + 1e-13],
        [2.5, 2.0],
        [2.5, 5.0],
        [2.5 - 1e-13, 6.5],
        [0.0, 5.0],
    ];
    assert_every_start_closed(&l, "L shape");
}

#[test]
fn the_issue_ring_closes_from_every_start() {
    let ring = [
        [132.8, -8.0],
        [132.4, -8.0],
        [132.4, -7.60000000000037],
        [132.8, -7.60000000000037],
        [132.8, -3.96250000000036],
        [127.4, -3.96250000000034],
        [127.3375, -3.96250000000034],
        [127.2622304, -3.96250000000034],
        [125.549671244606, -3.96250000000034],
        [125.549671244606, -5.36699040000169],
        [125.549671244606, -5.52000000000169],
        [125.549671244606, -9.98250000000169],
        [125.549671244606, -10.1355096000017],
        [125.549671244606, -12.2419904000017],
        [125.549671244606, -12.3950000000017],
        [125.549671244606, -14.2200000000018],
        [125.549671244606, -14.9450000000018],
        [125.549671244606, -15.6000000000008],
        [127.4, -15.6000000000008],
        [127.4, -14.4530096000007],
        [127.4, -14.2375000000007],
        [127.5377696, -14.2375000000007],
        [132.8, -14.2375000000008],
    ];
    assert_every_start_closed(&ring, "issue ring");
}

#[test]
fn a_long_side_with_a_vertex_50_nanometres_off_it_still_extrudes_closed() {
    // Relatively flat (5e-8 over 100 m) but a doubled area of 5e-6, over the
    // audit's floor: the old triangulation meshed it closed, so it is not
    // refused.
    assert_every_start_closed(
        &[
            [50.0, -5e-8],
            [100.0, 0.0],
            [100.0, 10.0],
            [0.0, 10.0],
            [0.0, 0.0],
        ],
        "100 m side",
    );
}

#[test]
fn a_cap_that_must_be_flat_is_refused_by_name() {
    let error = extrude_profile(
        &rings(&[[0.0, 0.0], [10.0, 0.0], [5.0, 1e-13]]),
        Vec3::Z,
        1.0,
        TOL,
    )
    .expect_err("a flat sliver cap is refused");
    assert!(error.to_string().contains("cap triangle"), "{error}");
}

#[test]
fn a_sliver_that_is_not_flat_and_cannot_be_avoided_is_still_meshed() {
    // The whole profile is the sliver: a fine chord's cap is legitimately
    // this thin, so it is meshed below the audit's floor, not refused.
    let mesh = extrude_profile(
        &rings(&[[0.0, 0.0], [10.0, 0.0], [5.0, 5e-8]]),
        Vec3::Z,
        1.0,
        TOL,
    )
    .expect("meshed");
    assert!(audit_mesh(&mesh, TOL).degenerate_triangles > 0);
}
