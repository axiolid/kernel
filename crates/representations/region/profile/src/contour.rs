//! Arbitrary closed contours.

use axiolid_core::Interval;
use axiolid_curve::Curve2;

/// One bounded oriented curve segment.
#[derive(Debug, Clone, PartialEq)]
pub struct ProfileSegment {
    /// Exact supporting curve.
    pub curve: Curve2,
    /// Parameter interval on `curve`.
    ///
    /// Always a finite span: a contour edge is bounded, so an unbounded
    /// support such as a line contributes only this interval, never itself.
    pub domain: Interval,
    /// Whether parameter direction follows contour orientation.
    pub same_sense: bool,
}

/// One oriented closed contour assembled from bounded exact curve segments.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Contour {
    /// Consecutive bounded segments; closure is validated separately.
    pub segments: Vec<ProfileSegment>,
}

impl Contour {
    /// Construct from oriented segments without silently repairing gaps.
    pub fn new(segments: Vec<ProfileSegment>) -> Self {
        Self { segments }
    }

    /// Number of exact segments.
    pub fn len(&self) -> usize {
        self.segments.len()
    }

    /// Whether this contour has no segments.
    pub fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }
}

/// Closed outer contour and zero or more closed holes.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ContourProfile {
    /// Outer boundary.
    pub outer: Contour,
    /// Inner boundaries.
    pub holes: Vec<Contour>,
}
