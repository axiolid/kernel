//! A skin in a difference's or intersection's result, thinner than rounding
//! (#400).
//!
//! The snap ([`crate::snap`]) closes a skin or sliver a rounding residue
//! thick only where the faces are axis-aligned; a clip turned by a placement
//! lands its faces on planes no `f64` coordinate can hold exactly, and the
//! mesh boolean, which has no tolerance, leaves a lid between them. The
//! result is closed and consistently wound, so the closure check labels it a
//! Solid, yet two of its faces face each other within the rounding of the
//! coordinates they were placed at: residue, not authored geometry (ADR 0096,
//! D7). [`axiolid_inspect::thin_sheet`] finds that pair on the result, and the
//! cut is refused by name instead of returned.
//!
//! A skin above the rounding but below the tolerance stays authored (ADR 0096,
//! D7): the guard reads the rounding reach of [`crate::snap::reach`], never
//! the tolerance. A union is not guarded here (#434).

use axiolid_contracts::GeomError;
use axiolid_core::{BooleanOperator, Scalar};
use axiolid_inspect::ThinSheet;
use axiolid_mesh::TriMesh;

/// The start of the `GeomError::Degenerate` message a difference or
/// intersection is refused with when its result holds a skin thinner than
/// its coordinates' rounding (#400): two faces face each other on one plane,
/// closer than the snap's reach.
pub const SKIN_WITHIN_ROUNDING: &str = "boolean result holds a skin thinner than its coordinates' rounding: two of its faces face each other on one plane";

/// The refusal for `found` in `mesh`, the result of `operator`, read at the
/// rounding `reach`: the constant, then the two triangles, their offset, the
/// reach and where the first one lies.
pub(crate) fn refusal(
    found: &ThinSheet,
    mesh: &TriMesh,
    reach: Scalar,
    operator: BooleanOperator,
) -> GeomError {
    let centroid = mesh
        .indices
        .get(found.first * 3..found.first * 3 + 3)
        .and_then(|corner| {
            let at = |i: u32| mesh.positions.get(i as usize).copied();
            Some((at(corner[0])? + at(corner[1])? + at(corner[2])?) / 3.0)
        });
    GeomError::Degenerate(format!(
        "{SKIN_WITHIN_ROUNDING}: {operator:?} result, triangles {} and {} {:e} apart, \
         within the rounding reach {reach:e}, near {centroid:?}",
        found.first, found.second, found.offset
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::Point3;

    #[test]
    fn the_refusal_starts_with_the_constant_and_names_the_pair() {
        // A fin folded back onto its face about a shared edge.
        let mesh = TriMesh::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 1e-12),
            ],
            vec![0, 1, 2, 1, 0, 3],
        );
        let found = axiolid_inspect::thin_sheet(&mesh, 1e-9).expect("a fin on its face");
        let GeomError::Degenerate(text) = refusal(&found, &mesh, 1e-9, BooleanOperator::Difference)
        else {
            panic!("a degenerate refusal");
        };
        assert!(text.starts_with(SKIN_WITHIN_ROUNDING), "{text}");
        assert!(text.contains("triangles 0 and 1"), "{text}");
        assert!(text.contains("Difference"), "{text}");
    }
}
