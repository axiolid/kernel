//! Centre-line profiles: an open path plus a constant width.

use axiolid_core::Scalar;

use crate::contour::Contour;

/// An open path swept by a constant width, centred on the path.
///
/// This is a *representation*, not a contour: the closed area it denotes is
/// the set of points within `half_width` of the path. Resolving it into a
/// boundary is an offsetting algorithm, which belongs above this tier.
///
/// Keeping it exact matters because offsetting is lossy in a way that cannot
/// be undone. Once a centre line is flattened into an outline, the fact that
/// it *was* a constant-width path is gone, and with it the ability to change
/// the width, re-offset at a different tolerance, or round-trip the source
/// unchanged. Sheet metal, curtain-wall mullions and pipe centre lines are all
/// authored this way.
///
/// The path is held as a [`Contour`] because a centre line is a sequence of
/// bounded curve segments exactly like a boundary is; the difference is that
/// this one is not closed and does not enclose anything on its own.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct CenterLineProfile {
    /// Open path the width is measured from.
    pub path: Contour,
    /// Half the constant width, measured perpendicular to the path.
    ///
    /// Stored as a half-width rather than a full thickness so the two offset
    /// sides are symmetric by construction: a consumer offsets by `+half` and
    /// `-half` without having to remember to halve anything. Halving at the
    /// boundary is the kind of factor-of-two error that produces geometry
    /// which looks plausible and is wrong.
    pub half_width: Scalar,
}

impl CenterLineProfile {
    /// Construct from an open path and the full width across it.
    pub fn from_width(path: Contour, width: Scalar) -> Self {
        Self {
            path,
            half_width: width / 2.0,
        }
    }

    /// The full width across the path.
    pub fn width(&self) -> Scalar {
        self.half_width * 2.0
    }
}
