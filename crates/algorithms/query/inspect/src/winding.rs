//! Certificate that a closed oriented cycle's winding number is 0 or 1
//! (#420, ADR 0099 D6).
//!
//! `axiolid_measure::cycle_volume_properties` returns the integral of the
//! surface's integer winding number `w`, which is the enclosed volume only
//! where `w` is 0 or 1: duplicated or overlapping shells count twice and an
//! inverted shell subtracts. This module bounds how much volume breaks that.
//!
//! # Moments of the winding number
//!
//! The prism decomposition of [`crate::overlap`] needs no solid, only a
//! closed oriented 2-cycle. Above a plane `z = z0` below every vertex, the
//! winding number at a point is the signed count of faces above it,
//!
//! ```text
//! w = sum over faces f of  s_f * 1{ z0 < z < h_f(x, y) }
//! ```
//!
//! with `s_f` the sign of the face's shadow `f'` (`+1` facing up). Then
//!
//! ```text
//! integral of w   = sum over f of     s_f * integral over f' of (h_f - z0)
//! integral of w^2 = sum over f, g of  s_f s_g * integral over f' and g' of (min(h_f, h_g) - z0)
//! ```
//!
//! The first is the signed volume and the second the unclamped pair sum
//! behind `intersection_volume`, here of the mesh against itself over the
//! faces as given (no flip of an inward-wound mesh, which would change the
//! sign of `w`).
//!
//! # The excess
//!
//! `w` is an integer, so `w^2 - w = w (w - 1)` is 0 where `w` is 0 or 1 and
//! at least 2 elsewhere. Hence
//!
//! ```text
//! vol{ w not in {0, 1} }  <=  E := (integral of w^2 - integral of w) / 2
//!                              =   integral of w (w - 1) / 2
//! ```
//!
//! `E` weighs a point by `w (w - 1) / 2`, which is 0 at `w` = 0 or 1, 1 at
//! `w` = -1 or 2 and 3 at `w` = 3: it
//! equals the volume only where `w` is in {-1..2}, and a triple-welded unit
//! cube has `E` = 3 over a bad volume of 1. `E` is 0 exactly when `w` is 0
//! or 1 almost everywhere. Both
//! integrals are outward-rounded intervals, so the excess is an interval
//! clamped below at 0; its upper end bounds the volume where the winding
//! number leaves {0, 1}. The returned [`VolumeInterval`] brackets `E`, not
//! that volume, so its "true volume" wording reads as `E` here. Edge-touching shells and a cavity (an inverted
//! inner shell, `w` = 0 inside it) keep the excess at 0; a duplicated,
//! overlapping or inverted shell makes it positive.

use crate::overlap::{clamp, faces_as_given, lowest, pair_sum, signed_volume, Iv, VolumeInterval};
use axiolid_core::Tolerance;
use axiolid_mesh::{audit_mesh_cycle, MeshCycleHealth, TriangleMeshView};

/// Why the winding number was not certified.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum WindingError {
    /// A coordinate is not finite.
    NonFinite,
    /// The mesh is not a closed oriented cycle; the audit says why.
    NotClosedOrientedCycle(Box<MeshCycleHealth>),
    /// The winding excess `E = integral of w(w-1)/2` (an upper bound on the
    /// volume with `w` outside {0, 1}, equal to it where `w` is in
    /// {-1..2}) is certified greater than the budget.
    Exceeds {
        /// Bracket of `E`, its lower end above the budget.
        excess: VolumeInterval,
    },
    /// The bracket straddles the budget: neither certified within nor over.
    Undecided {
        /// Bracket of `E`.
        excess: VolumeInterval,
    },
    /// The budget is not a finite positive number.
    InvalidBudget,
}

impl core::fmt::Display for WindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFinite => f.write_str("a coordinate is not finite"),
            Self::NotClosedOrientedCycle(_) => {
                f.write_str("the mesh is not a closed oriented cycle")
            }
            Self::Exceeds { excess } => write!(
                f,
                "the winding excess (an upper bound on the volume with winding number outside 0 and 1) is at least {}, over the budget",
                excess.lower
            ),
            Self::Undecided { excess } => write!(
                f,
                "the winding excess (an upper bound on the volume with winding number outside 0 and 1) lies in [{}, {}], straddling the budget",
                excess.lower, excess.upper
            ),
            Self::InvalidBudget => f.write_str("the budget is not a finite positive number"),
        }
    }
}

impl core::error::Error for WindingError {}

/// Certifies that the winding excess `E = integral of w(w-1)/2` of a
/// closed oriented cycle, an upper bound on the volume with `w` not in
/// {0, 1} that equals it where `w` is in {-1..2}, is at most `max_excess`,
/// and returns the bracket of `E` (`upper <= max_excess`; the interval's
/// "volume" is `E`, not the bad volume).
///
/// `max_excess` must be a small positive volume scaled to the model: the
/// bracket is outward-rounded, so a budget of 0 can never be met.
///
/// # Errors
///
/// [`WindingError::InvalidBudget`] for a zero, negative or non-finite
/// `max_excess`, [`WindingError::NonFinite`] for a non-finite position,
/// [`WindingError::NotClosedOrientedCycle`] with the audit when the mesh is
/// not a closed oriented cycle, [`WindingError::Exceeds`] when the excess is
/// certified above the budget and [`WindingError::Undecided`] when its
/// bracket straddles it.
pub fn certify_unit_winding<M: TriangleMeshView + ?Sized>(
    mesh: &M,
    tolerance: Tolerance,
    max_excess: f64,
) -> Result<VolumeInterval, WindingError> {
    if !max_excess.is_finite() || max_excess <= 0.0 {
        return Err(WindingError::InvalidBudget);
    }
    if (0..mesh.position_count()).any(|i| !mesh.position(i).is_finite()) {
        return Err(WindingError::NonFinite);
    }
    let audit = audit_mesh_cycle(mesh, tolerance);
    if !audit.is_closed_oriented_cycle() {
        return Err(WindingError::NotClosedOrientedCycle(Box::new(audit)));
    }
    let base = lowest(mesh);
    let faces = faces_as_given(mesh);
    let second = pair_sum(&faces, &faces, base);
    let first = signed_volume(&faces, base);
    let excess = clamp(second.sub(first).mul(Iv::point(0.5)), 0.0, f64::INFINITY);
    if excess.upper <= max_excess {
        Ok(excess)
    } else if excess.lower > max_excess {
        Err(WindingError::Exceeds { excess })
    } else {
        Err(WindingError::Undecided { excess })
    }
}
