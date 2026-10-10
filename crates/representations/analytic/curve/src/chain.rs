//! Arc-length chains: plane curves parameterised by cumulative arc length,
//! built from pieces placed rigidly end to end.
//!
//! A horizontal alignment is authored as a run of segments, each written in
//! its own frame and each starting where the previous one ends, in the
//! direction it ends in. Most segment kinds have a natural equation (line,
//! arc, clothoid) and fit a [`CurvatureLaw`] piece. Some do not: a cubic
//! parabola `y = x^3 / (6 R L)` is an exact polynomial graph, but its
//! curvature is not elementary in arc length, and its length along the
//! curve inverts an elliptic integral. A [`Chain2`] holds both kinds at once:
//! [`ChainPiece2::Intrinsic`] carries a curvature law, and
//! [`ChainPiece2::Parametric`] carries any [`Curve2`] in its own local frame,
//! read from a start parameter over a stated arc length.
//!
//! # The parameter is arc length
//!
//! The chain's parameter is distance along it, from `0` at its start to the
//! sum of the piece lengths at its end. That is what makes a chain usable as
//! the plan of an [`Elevated3`](crate::Elevated3), whose elevation law is
//! written against plan distance (ADR 0060): a parametric piece is read by
//! arc length, never by its own parameter.
//!
//! # Rigid placement
//!
//! Each piece starts at the end point of the previous one with the previous
//! end tangent as its local `+x` axis, so a chain is tangent-continuous at
//! every join by construction. A parametric piece's curve must itself start
//! at its local origin with tangent `+x` at `start`; an evaluator refuses one
//! that does not rather than move it there. The first piece is placed by
//! [`Chain2::start`], read like [`Intrinsic2::start`](crate::Intrinsic2).
//!
//! # Exact data, evaluated to tolerance
//!
//! Like [`Intrinsic2`](crate::Intrinsic2), a chain is exact stored data and
//! this module computes no points. Where a piece ends depends on the arc
//! length of its curve, which is a quadrature, and a parametric piece is
//! read at a distance by inverting that quadrature; both belong to an
//! evaluator that can state its tolerance (`axiolid-evaluate`). Whether a
//! parametric piece's curve is long enough for its `length` is therefore an
//! evaluator's refusal, not a structural one.

use axiolid_core::{Frame2, Scalar};

use crate::{CurvatureLaw, Curve2};

/// One piece of a [`Chain2`], written in its own local frame: origin at the
/// piece start, `+x` along the start tangent.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub enum ChainPiece2 {
    /// A natural-equation piece: curvature as a function of arc length from
    /// the piece start, over `length`.
    Intrinsic {
        /// Curvature law, in the piece's own arc length.
        curvature: CurvatureLaw,
        /// Arc length of the piece.
        length: Scalar,
    },
    /// A parametric curve read by arc length.
    ///
    /// `curve` is written in the piece's local frame: it passes through the
    /// origin at parameter `start` with its tangent along `+x` there, and the
    /// piece is the stretch of it from `start`, in increasing parameter,
    /// whose arc length is `length`.
    Parametric {
        /// The curve, in the piece's local frame.
        curve: Curve2,
        /// Curve parameter where the piece begins.
        start: Scalar,
        /// Arc length of the piece along `curve`.
        length: Scalar,
    },
}

impl ChainPiece2 {
    /// Arc length of the piece, as stored.
    #[must_use]
    pub fn length(&self) -> Scalar {
        match self {
            Self::Intrinsic { length, .. } | Self::Parametric { length, .. } => *length,
        }
    }
}

/// A plane curve parameterised by cumulative arc length, its pieces placed
/// rigidly end to end. See the [module documentation](self).
///
/// Dirty imported data stays representable, as everywhere else in this
/// crate: [`Self::is_well_formed`] names what an evaluator refuses
/// structurally.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct Chain2 {
    /// Start frame: origin at the chain start, `x` along the start tangent.
    pub start: Frame2,
    /// Pieces in order of travel.
    pub pieces: Vec<ChainPiece2>,
}

impl Chain2 {
    /// A chain of `pieces` starting at `start`.
    #[must_use]
    pub const fn new(start: Frame2, pieces: Vec<ChainPiece2>) -> Self {
        Self { start, pieces }
    }

    /// Whether the chain is structurally sound: at least one piece, a finite
    /// start frame with independent axes, every length finite and positive,
    /// every curvature law well formed, every parametric start finite.
    ///
    /// Whether a parametric piece's curve has `length` of arc after `start`
    /// needs a quadrature, so it is an evaluator's check, not this one.
    #[must_use]
    pub fn is_well_formed(&self) -> bool {
        let frame = &self.start;
        !self.pieces.is_empty()
            && frame.origin.is_finite()
            && frame.x.is_finite()
            && frame.y.is_finite()
            && frame.x.perp_dot(frame.y) != 0.0
            && self.pieces.iter().all(|piece| {
                let length = piece.length();
                length.is_finite()
                    && length > 0.0
                    && match piece {
                        ChainPiece2::Intrinsic { curvature, .. } => curvature.is_well_formed(),
                        ChainPiece2::Parametric { start, .. } => start.is_finite(),
                    }
            })
    }

    /// Total arc length: the sum of the piece lengths, `None` when one is
    /// not finite and positive.
    #[must_use]
    pub fn length(&self) -> Option<Scalar> {
        self.pieces.iter().try_fold(0.0, |total, piece| {
            let length = piece.length();
            (length.is_finite() && length > 0.0).then_some(total + length)
        })
    }

    /// Arc length at each interior join, ascending: where piece `i + 1`
    /// begins. Empty for one piece; `None` when a length is not finite and
    /// positive.
    #[must_use]
    pub fn joins(&self) -> Option<Vec<Scalar>> {
        let mut total = 0.0;
        let mut out = Vec::with_capacity(self.pieces.len().saturating_sub(1));
        for (index, piece) in self.pieces.iter().enumerate() {
            let length = piece.length();
            if !(length.is_finite() && length > 0.0) {
                return None;
            }
            total += length;
            if index + 1 < self.pieces.len() {
                out.push(total);
            }
        }
        Some(out)
    }

    /// The piece covering arc length `s`, its index, and `s` rebased to the
    /// piece start.
    ///
    /// Pieces are half-open, so a join belongs to the piece that starts
    /// there; the last piece is closed at the chain end. `None` outside
    /// `[0, length]`, for a non-finite `s` and for a malformed length.
    #[must_use]
    pub fn piece_at(&self, s: Scalar) -> Option<(usize, &ChainPiece2, Scalar)> {
        let total = self.length()?;
        if !s.is_finite() || s < 0.0 || s > total {
            return None;
        }
        let mut begin = 0.0;
        let last = self.pieces.len() - 1;
        for (index, piece) in self.pieces.iter().enumerate() {
            let end = begin + piece.length();
            if s < end || index == last {
                return Some((index, piece, s - begin));
            }
            begin = end;
        }
        None
    }
}
