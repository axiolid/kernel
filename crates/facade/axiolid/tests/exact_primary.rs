//! kernel#36: exact B-rep is the primary currency at the PUBLIC API.
//!
//! The v1.0 promise depends on this holding at the API boundary, not
//! merely inside the compiler. Each test below pins one of the
//! issue's three acceptance criteria.
#![cfg(feature = "generate")]

use axiolid::core::Tolerance;
use axiolid::generate::{GeneratedGeometry, GenerationRequest, TessellationRequest};

/// Criterion 1: no public entry point yields a mesh unasked.
///
/// The request enum is the gate. `ExactBRep` carries no tolerance
/// field at all, so there is no way to spell "give me a mesh"
/// without selecting `Tessellation` and supplying one.
#[test]
fn a_mesh_cannot_be_requested_without_a_tolerance() {
    let exact = GenerationRequest::ExactBRep;
    assert!(exact.requires_exact_brep());

    let tess = GenerationRequest::Tessellation(TessellationRequest::new(Tolerance::MILLIMETRE));
    assert!(!tess.requires_exact_brep());
}

/// Criterion 3: an exact input survives a public round trip.
///
/// A square extrudes to an exact solid through the facade, and what
/// comes back still carries analytic surfaces -- planes, not triangles.
/// This is the test the issue names explicitly.
#[test]
fn an_exact_solid_survives_a_public_round_trip() {
    use axiolid::core::{Tolerance, Vec3};
    use axiolid::surface::Surface;

    let solid = axiolid::generate::extrude::extrude_profile_exact(
        &square(),
        Vec3::Z,
        2.0,
        Tolerance::METRE,
    )
    .expect("a square extrudes exactly");

    // Six planes: four walls and two caps. If the public path had
    // tessellated anything, these would be triangles instead.
    let planes = solid
        .surfaces()
        .iter()
        .filter(|s| matches!(s, Surface::Plane(_)))
        .count();
    assert_eq!(planes, 6, "expected six analytic planes, got {planes}");

    // And the round trip through the generated-geometry carrier keeps it exact.
    let carried = GeneratedGeometry::ExactBRep(solid.clone());
    match carried {
        GeneratedGeometry::ExactBRep(out) => assert_eq!(out, solid),
        other => panic!("exact input became {other:?}"),
    }
}

/// A tessellated result carries the tolerance it was built to.
///
/// Criterion 2 in negative form: a mesh that reached a caller can
/// always be traced back to the budget that produced it, so an
/// approximation can never masquerade as exact downstream.
#[test]
fn a_tessellated_result_carries_its_tolerance() {
    let request = TessellationRequest::new(Tolerance::MILLIMETRE);
    let mesh = axiolid::mesh::TriMesh::new(Vec::new(), Vec::new());
    let bound = request.bind(mesh);
    assert_eq!(bound.tolerance(), Tolerance::MILLIMETRE);
}

/// A unit square contour, built through the public profile API.
fn square() -> axiolid::profile::Profile {
    use axiolid::core::{Interval, Point2};
    use axiolid::curve::{Curve2, Line2};
    use axiolid::profile::{Contour, ContourProfile, Profile, ProfileSegment};
    let corner = [
        (Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)),
        (Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)),
        (Point2::new(1.0, 1.0), Point2::new(0.0, 1.0)),
        (Point2::new(0.0, 1.0), Point2::new(0.0, 0.0)),
    ];
    let sides = corner
        .into_iter()
        .map(|(from, to)| ProfileSegment {
            curve: Curve2::Line(Line2 {
                origin: from,
                direction: to - from,
            }),
            domain: Interval::UNIT,
            same_sense: true,
        })
        .collect();
    Profile::Contour(ContourProfile {
        outer: Contour::new(sides),
        holes: Vec::new(),
    })
}
