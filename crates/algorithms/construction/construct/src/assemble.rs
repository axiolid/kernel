//! Put several exact solids into one `ExactBRep`.
//!
//! A composite profile whose members do not touch is several bodies acting
//! as one section (#111). The graph compiles one node to one `ExactBRep`, so
//! the bodies travel together: one `Solid` per piece, each with its own
//! outer shell and voids. Nothing is unioned or glued here; the pieces must
//! already be disjoint.

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_topology::Solid;

/// One `ExactBRep` holding every solid of every piece, in order.
pub(crate) fn merge_solids(mut pieces: Vec<ExactBRep>) -> GeomResult<ExactBRep> {
    if pieces.len() == 1 {
        return Ok(pieces.remove(0));
    }
    let mut builder = ExactBRepBuilder::default();
    for piece in &pieces {
        let shells = builder.append(piece, false);
        for solid in piece.topology().solids() {
            let outer = shells[solid.outer.index()];
            let voids = solid
                .voids
                .iter()
                .map(|void| shells[void.index()])
                .collect();
            builder.topology_mut().add_solid(Solid { outer, voids });
        }
    }
    builder.finish().map_err(|error| {
        GeomError::Degenerate(format!("disjoint solids did not assemble: {error}"))
    })
}
