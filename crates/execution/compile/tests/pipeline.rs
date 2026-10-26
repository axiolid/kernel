//! The composability gate: meshes this crate produces must be accepted by the
//! boolean provider, unmodified.
//!
//! Edge-manifoldness is a local proxy; this is the real thing. If extrusion
//! output cannot enter `axiolid-mesh-boolean-boolmesh`, the pipeline does not exist no matter
//! how good each half looks alone.

use axiolid_construct::extrude::extrude;
use axiolid_construct::profile::{profile_rings, triangulate, Rings};
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Tolerance, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
use axiolid_profile::{Profile, RectangleProfile};

fn volume(m: &TriMesh) -> f64 {
    m.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let (a, b, c) = (
                m.positions[t[0] as usize],
                m.positions[t[1] as usize],
                m.positions[t[2] as usize],
            );
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0
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

#[allow(clippy::single_range_in_vec_init)]
fn loops(rings: &Rings) -> Vec<core::ops::Range<usize>> {
    let mut out = vec![0..rings.outer.len()];
    let mut start = rings.outer.len();
    for h in &rings.holes {
        out.push(start..start + h.len());
        start += h.len();
    }
    out
}

/// Extrude a box centred on (cx, cy) from z0, of the given size.
fn box_solid(cx: f64, cy: f64, z0: f64, sx: f64, sy: f64, sz: f64) -> TriMesh {
    let rings = profile_rings(&rect(sx, sy), 1e-4, Tolerance::METRE).expect("rings");
    let (pts, tris) = triangulate(&rings).expect("triangulate");
    let mut mesh = extrude(&pts, &tris, &loops(&rings), Vec3::Z, sz).expect("extrude");
    for p in &mut mesh.positions {
        *p = Point3::new(p.x + cx, p.y + cy, p.z + z0);
    }
    mesh
}

/// The IFC-dominant case, built entirely from this crate's output.
///
/// Wall 4 x 0.2 x 3 minus one 1.0 x 0.4 x 1.2 opening.
#[test]
fn an_extruded_wall_can_be_cut_by_the_boolean_provider() {
    let wall = box_solid(2.0, 0.1, 0.0, 4.0, 0.2, 3.0);
    let opening = box_solid(1.5, 0.1, 0.3, 1.0, 0.4, 1.2);

    assert!(volume(&wall) > 0.0, "wall must be outward-oriented");
    assert!((volume(&wall) - 2.4).abs() < 1e-9);

    let provider = BoolmeshBoolean::new();
    let options = ExecutionOptions::new(Tolerance::METRE);
    let cut = provider
        .boolean(&wall, &opening, BooleanOperator::Difference, &options)
        .expect("the provider must accept compiler output")
        .mesh;

    // The opening spans the full 0.2 m thickness, so the removed volume is
    // 1.0 x 0.2 x 1.2 = 0.24, not the opening's own 0.48.
    let expected = 2.4 - 0.24;
    let actual = volume(&cut);
    assert!(
        (actual - expected).abs() < 1e-9,
        "expected {expected}, got {actual}"
    );
}

/// Conservation across the seam: the two halves must reconstitute the wall.
#[test]
fn difference_and_intersection_partition_an_extruded_wall() {
    let wall = box_solid(2.0, 0.1, 0.0, 4.0, 0.2, 3.0);
    let tool = box_solid(1.5, 0.1, 0.3, 1.0, 0.4, 1.2);
    let provider = BoolmeshBoolean::new();
    let options = ExecutionOptions::new(Tolerance::METRE);

    let diff = provider
        .boolean(&wall, &tool, BooleanOperator::Difference, &options)
        .expect("difference")
        .mesh;
    let isect = provider
        .boolean(&wall, &tool, BooleanOperator::Intersection, &options)
        .expect("intersection")
        .mesh;

    let sum = volume(&diff) + volume(&isect);
    assert!(
        (sum - volume(&wall)).abs() < 1e-9,
        "a\\b + a^b = {sum} must equal a = {}",
        volume(&wall)
    );
}
