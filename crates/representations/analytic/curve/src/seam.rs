//! Which side of a seam a position on it is read from (#263).
//!
//! A composite curve -- a polyline, an arc-length chain, an elevated curve
//! whose profile is piecewise, a banked curve whose cant law is -- is made
//! of pieces laid end to end. Where two pieces meet, the position is shared
//! but the tangent, or the roll of a banked section, may jump: the curve
//! has two frames there. Every evaluator in this workspace that picks a
//! piece by a distance or a parameter takes the piece that STARTS at the
//! seam ([`ElevationLaw::piece_at`](crate::ElevationLaw::piece_at),
//! [`Chain2::piece_at`](crate::Chain2::piece_at),
//! [`CantLaw`](crate::CantLaw), a polyline's spans): that is
//! [`SeamSide::Outgoing`], the default. [`SeamSide::Incoming`] reads the
//! piece that ENDS there instead, as some authoring rules require (the
//! previous segment's tangent governs at a seam).
//!
//! At the curve's own start and end only one piece exists, and both sides
//! read it.

/// The piece a position exactly on a seam is read from.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub enum SeamSide {
    /// The piece that starts at the seam: the frame the curve is about to
    /// take. What every evaluator reads by default.
    #[default]
    Outgoing,
    /// The piece that ends at the seam: the frame the curve arrived with.
    Incoming,
}
