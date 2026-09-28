//! Generated solids must be accepted by a boolean provider, unmodified.
//!
//! `axiolid-construct` sits at L2 and cannot depend on a provider, so this
//! seam is verified here, at L3, where both sides are already in scope. The
//! claim is not "extrusion looks right" but "extrusion output satisfies the
//! preconditions a real boolean implementation enforces".

use axiolid_construct::extrude::extrude_profile;
use axiolid_construct::profile::profile_rings;
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Scalar, Tolerance, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

/// Volume measured through the audited oracle.
///
/// `volume_properties` refuses a mesh that is not closed and two-manifold, so
/// calling it is itself an assertion about watertightness.
fn volume_at(mesh: &TriMesh, tolerance: Tolerance) -> Scalar {
    axiolid_measure::volume_properties(mesh, tolerance)
        .expect("generated output must be closed and two-manifold")
        .signed_volume
}

/// A tolerance proportional to a chord budget, floored at f64 sanity.
///
/// The mesh audit calls a triangle degenerate when its doubled area is at
/// most `tolerance.linear()²`. A cylinder flattened at chord `c` has side
/// quads about `sqrt(8 r c)` wide and cap slivers far smaller, so a fixed
/// `Tolerance::MILLIMETRE` rejects correct geometry as soon as the chord
/// budget goes sub-millimetre. The measuring tolerance has to follow the
/// budget that drove flattening.
fn tolerance_for(chord: Scalar) -> Tolerance {
    Tolerance::new((chord * 1e-3).max(1e-12), 1e-9).expect("valid tolerance")
}

/// Every interior edge must be shared by exactly two oppositely-wound faces.
fn assert_closed_manifold(mesh: &TriMesh, what: &str) {
    use std::collections::HashMap;
    let mut edges: HashMap<(u32, u32), i32> = HashMap::new();
    for c in mesh.indices.chunks_exact(3) {
        for (a, b) in [(c[0], c[1]), (c[1], c[2]), (c[2], c[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            *edges.entry(key).or_default() += if a < b { 1 } else { -1 };
        }
    }
    for (edge, balance) in edges {
        assert_eq!(
            balance, 0,
            "{what}: edge {edge:?} is not shared by two oppositely-wound faces"
        );
    }
}

#[test]
fn an_extruded_solid_is_accepted_by_the_conformance_gated_boolean() {
    let disc = Profile::Circle(CircleProfile {
        radius: 2.0,
        thickness: None,
    });
    let rings = profile_rings(&disc, 1e-4, Tolerance::MILLIMETRE).unwrap();
    let cylinder = extrude_profile(&rings, Vec3::Z, 4.0, Tolerance::MILLIMETRE).unwrap();

    let bar = Profile::Rectangle(RectangleProfile {
        x: 1.0,
        y: 8.0,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    });
    let bar_rings = profile_rings(&bar, 1e-4, Tolerance::MILLIMETRE).unwrap();
    let block = extrude_profile(&bar_rings, Vec3::Z, 2.0, Tolerance::MILLIMETRE).unwrap();

    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let outcome = BoolmeshBoolean::new()
        .boolean(&cylinder, &block, BooleanOperator::Difference, &options)
        .expect("generated solids must satisfy the boolean preconditions");

    // The cut removed material but did not annihilate the cylinder.
    let cut = volume_at(&outcome.mesh, tolerance_for(1e-4));
    let whole = volume_at(&cylinder, tolerance_for(1e-4));
    assert!(
        cut > 0.0 && cut < whole,
        "difference volume {cut} must be between 0 and {whole}"
    );
    assert_closed_manifold(&outcome.mesh, "boolean result");
}
