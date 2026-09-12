//! Recovering which input face each boolean output face came from.
//!
//! An exact 2D overlay never invents an edge: every edge of the result
//! lies on an edge of one input or the other. So a boolean output face is
//! always a FRAGMENT of an input face, and can say which -- no bookkeeping
//! through the overlay required, because the geometry itself carries the
//! answer.
//!
//! This module recovers that mapping after the fact, with exact
//! predicates. A side wall of the result is a fragment of the input wall
//! whose supporting line it lies on; the caps are fragments of the input
//! caps. Where no input edge supports an output edge -- which should not
//! happen for an exact overlay, but is not proven here -- the face is
//! reported anonymous rather than guessed at.

use axiolid_brep::{FaceName, Operand, SweptFace};
use axiolid_core::Point2;
use axiolid_guarantees::{Certified, Sign};
use axiolid_predicates::orient2d;

/// One ring of an operand cross-section, with the operand it belongs to.
pub(crate) struct OperandRings<'a> {
    /// Which operand these rings came from.
    pub operand: Operand,
    /// Cross-section rings: outer first, then holes.
    pub rings: &'a [Vec<Point2>],
}

/// Whether `point` lies exactly on the line through `a` and `b`.
///
/// Uses the filtered-then-exact `orient2d` cascade, so a point that is
/// nearly-but-not-quite collinear is correctly rejected rather than
/// absorbed by a tolerance band. That distinction is the whole reason this
/// recovery can be trusted: a wrong answer here would attach a fragment to
/// the wrong input face and carry the wrong material through.
fn on_supporting_line(a: Point2, b: Point2, point: Point2) -> bool {
    matches!(
        orient2d(a, b, point),
        Certified::Certain {
            sign: Sign::Zero,
            ..
        }
    )
}

/// Name the input wall that supports the output edge `start -> end`.
///
/// Both endpoints must lie on the input edge's supporting line. Testing
/// both is what rejects an edge that merely crosses the line at a point.
///
/// Returns `None` when no input edge supports it, which the caller must
/// treat as unnamed rather than substituting a default.
pub(crate) fn name_side_fragment(
    start: Point2,
    end: Point2,
    operands: &[OperandRings<'_>],
) -> Option<FaceName> {
    for operand in operands {
        // Ring 0 is the outer boundary and the rest are holes, matching the
        // order the sweep names its walls in, so the wall ordinal is
        // counted across rings in exactly that sequence.
        let mut ordinal = 0u32;
        for ring in operand.rings {
            for index in 0..ring.len() {
                let a = ring[index];
                let b = ring[(index + 1) % ring.len()];
                if on_supporting_line(a, b, start) && on_supporting_line(a, b, end) {
                    return Some(
                        FaceName::swept(SweptFace::Side(ordinal)).fragment(operand.operand),
                    );
                }
                ordinal = ordinal.checked_add(1)?;
            }
        }
    }
    None
}
