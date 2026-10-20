//! A memory budget for the straight-edge boolean (#338).
//!
//! [`crate::overlay_within`] and [`crate::union_soup_within`] read the
//! caller's [`axiolid_contracts::ExecutionOptions::memory_budget_bytes`].
//! The subdivision charges what it is about to hold before it allocates
//! it, and refuses with [`OverlayRefusal::BudgetExceeded`] as soon as the
//! charges pass the budget:
//!
//! - the operands' rings, by vertex, before they are snapped and reduced;
//! - the edges, before their boxes, trees and winding parts are built;
//! - every point a pair of edges meets at, before it is recorded (the
//!   pairs come from the box tree, so this stops a soup whose edges cross
//!   far more often than its output has corners before those crossings
//!   are held);
//! - the pieces, counted from the recorded points, before any is built;
//! - each piece's rings and carriers, as they are found.
//!
//! The charges are upper bounds on the heap each holds, measured with a
//! counting allocator by `tests/overlay_budget.rs` (which also checks that
//! a call never holds more than its budget, refused or not), so a call
//! admitted holds at most its budget in the subdivision, and a refused
//! call stops before it. Without a budget nothing is refused: the charges
//! are counted and never compared. The rounding and settling of the
//! output hold about the output, which no budget refuses.

use crate::OverlayError;

/// Bytes charged per input vertex: the copies the snap and the chain
/// reduction take, and the rings handed to the subdivision.
pub(crate) const PER_VERTEX: usize = 320;
/// Bytes charged per edge of the subdivision: the edge, its boxes in the
/// edge tree, its winding part and tree, and its first two stops.
pub(crate) const PER_EDGE: usize = 2_000;
/// Bytes charged per point recorded on an edge: the meeting, with room
/// for its vector to double.
pub(crate) const PER_MEETING: usize = 400;
/// Bytes charged per crossing computed once for repeated segments.
pub(crate) const PER_SHARED_CROSSING: usize = 200;
/// Bytes charged per piece: the piece, its raw edge, its stop, its
/// interned vertex and its share of the linking.
pub(crate) const PER_PIECE: usize = 600;
/// Bytes charged per ring containing or carrying a piece.
pub(crate) const PER_MEMBER: usize = 40;

/// Why a budgeted overlay gave no result.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum OverlayRefusal {
    /// An operand failed the overlay's checks, as [`crate::overlay_with`]
    /// and [`crate::union_soup_with`] refuse it.
    Invalid(OverlayError),
    /// The subdivision would hold more than the memory budget: refused
    /// before allocating past it. `needed` is what had been charged when
    /// the budget was passed, a lower bound on what the operation needs;
    /// `budget` is the caller's.
    BudgetExceeded {
        /// Bytes charged when the budget was passed.
        needed: usize,
        /// The caller's memory budget in bytes.
        budget: usize,
    },
}

impl From<OverlayError> for OverlayRefusal {
    fn from(error: OverlayError) -> Self {
        Self::Invalid(error)
    }
}

impl std::fmt::Display for OverlayRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Invalid(error) => write!(f, "invalid overlay operand: {error:?}"),
            Self::BudgetExceeded { needed, budget } => write!(
                f,
                "overlay exceeded its memory budget: {needed} bytes charged, {budget} allowed"
            ),
        }
    }
}

impl std::error::Error for OverlayRefusal {}

/// The running charge against an optional budget.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Charge {
    budget: Option<usize>,
    used: usize,
}

impl Charge {
    /// No budget: every charge is admitted.
    pub(crate) const fn unlimited() -> Self {
        Self {
            budget: None,
            used: 0,
        }
    }

    /// A charge against `budget` bytes, if any.
    pub(crate) const fn within(budget: Option<usize>) -> Self {
        Self { budget, used: 0 }
    }

    /// Charge `count` items of `each` bytes before allocating them.
    pub(crate) fn add(&mut self, count: usize, each: usize) -> Result<(), OverlayRefusal> {
        self.used = self.used.saturating_add(count.saturating_mul(each));
        match self.budget {
            Some(budget) if self.used > budget => Err(OverlayRefusal::BudgetExceeded {
                needed: self.used,
                budget,
            }),
            _ => Ok(()),
        }
    }
}

/// The error of an unbudgeted call: an unlimited charge never refuses.
pub(crate) fn unbudgeted(refusal: OverlayRefusal) -> OverlayError {
    match refusal {
        OverlayRefusal::Invalid(error) => error,
        OverlayRefusal::BudgetExceeded { .. } => {
            unreachable!("an unlimited charge admits every charge")
        }
    }
}
