//! Reconstruction identity: splitting a solid and re-uniting the pieces
//! must return one solid, not two shells.
//!
//! Ignored pending ADR 0048. This is a KNOWN upstream limitation shared
//! with crates.io boolmesh 0.1.9, recorded as axiolid/kernel#100. The test
//! is committed failing-but-ignored so the defect stays visible and the
//! fix has a gate to turn green, rather than living only in prose.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Tolerance};
use axiolid_mesh::{component_count, EdgeAdjacency, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;

/// Oriented box: centre, half-extents, rotation about z.
fn obb(centre: [f64; 3], half: [f64; 3], angle: f64) -> TriMesh {
    let (s, c) = angle.sin_cos();
    TriMesh::new(
        (0..8)
            .map(|i| {
                let sx = if i & 1 == 0 { -half[0] } else { half[0] };
                let sy = if i & 2 == 0 { -half[1] } else { half[1] };
                let sz = if i & 4 == 0 { -half[2] } else { half[2] };
                Point3::new(
                    centre[0] + sx * c - sy * s,
                    centre[1] + sx * s + sy * c,
                    centre[2] + sz,
                )
            })
            .collect(),
        vec![
            0, 2, 1, 1, 2, 3, 4, 5, 6, 5, 7, 6, 0, 1, 4, 1, 5, 4, 2, 6, 3, 3, 6, 7, 0, 4, 2, 2, 4,
            6, 1, 3, 5, 3, 7, 5,
        ],
    )
}

fn topology(mesh: &TriMesh) -> (i64, usize) {
    let adjacency = EdgeAdjacency::build(mesh);
    (adjacency.euler_characteristic(), component_count(mesh))
}

/// `(A-B) u (A^B)` must reconstruct `A` in topology, not only in volume.
///
/// The operands are the exactness-suite pair: an axis-aligned slab and a
/// box rotated 30 degrees that passes clean through it. The rotation
/// matters -- the axis-aligned equivalent already reconstructs exactly,
/// because its coordinates are exactly representable in binary floating
/// point and the two sides retriangulate identically.
#[test]
#[ignore = "known upstream limitation, see ADR 0048 and axiolid/kernel#100"]
fn reconstruction_preserves_topology() {
    let a = obb([2.0, 0.2, 1.5], [2.0, 0.2, 1.5], 0.0);
    let b = obb(
        [1.0, 0.1, 1.5],
        [0.6, 0.5, 0.7],
        std::f64::consts::FRAC_PI_6,
    );
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let provider = BoolmeshBoolean::new();
    let op = |s: &TriMesh, t: &TriMesh, o: BooleanOperator| {
        provider
            .boolean(s, t, o, &options)
            .expect("operands are valid closed solids")
            .mesh
    };

    let difference = op(&a, &b, BooleanOperator::Difference);
    let intersection = op(&a, &b, BooleanOperator::Intersection);
    let rebuilt = op(&difference, &intersection, BooleanOperator::Union);

    let (want_chi, want_components) = topology(&a);
    let (got_chi, got_components) = topology(&rebuilt);

    assert_eq!(
        (got_chi, got_components),
        (want_chi, want_components),
        "reconstruction changed the shape: expected chi={want_chi} comps={want_components}, \
         got chi={got_chi} comps={got_components}"
    );
}

/// The same reconstruction axis-aligned, which must pass today.
///
/// Guards the claim that the defect is specific to inexactly-representable
/// coordinates. If this ever fails, the diagnosis in ADR 0048 is wrong.
#[test]
fn axis_aligned_reconstruction_preserves_topology() {
    let a = obb([2.0, 0.2, 1.5], [2.0, 0.2, 1.5], 0.0);
    let b = obb([1.0, 0.2, 1.5], [0.6, 0.5, 0.7], 0.0);
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let provider = BoolmeshBoolean::new();
    let op = |s: &TriMesh, t: &TriMesh, o: BooleanOperator| {
        provider
            .boolean(s, t, o, &options)
            .expect("operands are valid closed solids")
            .mesh
    };

    let difference = op(&a, &b, BooleanOperator::Difference);
    let intersection = op(&a, &b, BooleanOperator::Intersection);
    let rebuilt = op(&difference, &intersection, BooleanOperator::Union);

    assert_eq!(
        topology(&rebuilt),
        topology(&a),
        "axis-aligned reconstruction must be exact"
    );
}
