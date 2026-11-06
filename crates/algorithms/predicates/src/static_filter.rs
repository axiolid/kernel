//! Static filters: bounds computed once from a coordinate range.
//!
//! The dynamic filter in each predicate computes a *permanent* -- a sum of
//! absolute products -- on every call, which costs roughly as much as the
//! determinant itself. When a caller can state an upper bound on coordinate
//! magnitude up front (a model's bounding box, a quantisation grid), the error
//! bound can be precomputed and the per-call work drops to one comparison.
//!
//! The trade is coverage, not correctness: a static bound is necessarily
//! looser than a per-input one, so it defers more often. Its bound covers
//! every rounding error of the f64 evaluation, so a sign it certifies is the
//! exact sign. It does not defer wherever the dynamic filter does: the dynamic
//! `orient2d` filter declines differences outside [1e-90, 1e90], which the
//! static one may still settle. A caller that exceeds the declared range gets
//! `None` rather than a silently invalid answer, and a declared range so small
//! that a bound would sit below the normal range, where products lose their
//! relative accuracy, is refused by [`StaticFilter::new`] (#429).

use axiolid_core::{Point2, Point3};
use axiolid_guarantees::Sign;

/// Machine epsilon for binary64.
const EPSILON: f64 = f64::EPSILON / 2.0;

/// Precomputed error bounds for a declared coordinate range.
///
/// Construct once per model, reuse across every predicate call.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StaticFilter {
    /// Maximum absolute coordinate value this filter is valid for.
    bound: f64,
    /// Absolute error bound for the `orient2d` determinant.
    orient2d: f64,
    /// Absolute error bound for the `orient3d` determinant.
    orient3d: f64,
}

impl StaticFilter {
    /// Build bounds valid for coordinates with `|x| <= bound`.
    ///
    /// Returns `None` for a non-finite or non-positive bound, and for a bound
    /// so large that the derived error bound is not finite: in both cases no
    /// sign could be certified, and returning a filter that always defers
    /// would hide the configuration error. A bound so small that a derived
    /// error bound falls below the normal range is refused as well: the
    /// determinant's products then underflow, the bound no longer covers the
    /// rounding error, and a bound of zero would certify every nonzero sign.
    #[must_use]
    pub fn new(bound: f64) -> Option<Self> {
        if !bound.is_finite() || bound <= 0.0 {
            return None;
        }
        // A coordinate difference is at most 2*bound, so a 2x2 determinant
        // term is at most (2*bound)^2 and a 3x3 term (2*bound)^3.
        let span = 2.0 * bound;
        let orient2d = (3.0 + 16.0 * EPSILON) * EPSILON * (2.0 * span * span);
        let orient3d = (7.0 + 56.0 * EPSILON) * EPSILON * (6.0 * span * span * span);
        let usable = |bound: f64| bound.is_finite() && bound >= f64::MIN_POSITIVE;
        if !usable(orient2d) || !usable(orient3d) {
            return None;
        }
        Some(Self {
            bound,
            orient2d,
            orient3d,
        })
    }

    /// The coordinate range this filter was built for.
    #[must_use]
    pub const fn bound(self) -> f64 {
        self.bound
    }

    /// Whether every coordinate of a 2D point is inside the declared range.
    #[must_use]
    fn covers2(self, p: Point2) -> bool {
        p.x.abs() <= self.bound && p.y.abs() <= self.bound
    }

    /// Whether every coordinate of a 3D point is inside the declared range.
    #[must_use]
    fn covers3(self, p: Point3) -> bool {
        p.x.abs() <= self.bound && p.y.abs() <= self.bound && p.z.abs() <= self.bound
    }
}

impl StaticFilter {
    /// Try to settle `orient2d` with the precomputed bound.
    ///
    /// `None` means "not settled": either a point lies outside the declared
    /// range, or the determinant is too small for the static bound. The caller
    /// must fall back to the full predicate, which is always correct.
    #[must_use]
    pub fn orient2d(self, a: Point2, b: Point2, c: Point2) -> Option<Sign> {
        if !(self.covers2(a) && self.covers2(b) && self.covers2(c)) {
            return None;
        }
        let determinant = (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x);
        decide(determinant, self.orient2d)
    }

    /// Try to settle `orient3d` with the precomputed bound.
    #[must_use]
    pub fn orient3d(self, a: Point3, b: Point3, c: Point3, d: Point3) -> Option<Sign> {
        if !(self.covers3(a) && self.covers3(b) && self.covers3(c) && self.covers3(d)) {
            return None;
        }
        let (adx, ady, adz) = (a.x - d.x, a.y - d.y, a.z - d.z);
        let (bdx, bdy, bdz) = (b.x - d.x, b.y - d.y, b.z - d.z);
        let (cdx, cdy, cdz) = (c.x - d.x, c.y - d.y, c.z - d.z);
        let determinant = adz * (bdx * cdy - cdx * bdy)
            + bdz * (cdx * ady - adx * cdy)
            + cdz * (adx * bdy - bdx * ady);
        decide(determinant, self.orient3d)
    }
}

/// Certify a sign only when the magnitude strictly exceeds the bound.
///
/// Strict, not `>=`: a determinant exactly equal to its own error bound could
/// have true value zero, so the sign is not proven.
#[inline]
#[must_use]
fn decide(determinant: f64, bound: f64) -> Option<Sign> {
    if determinant > bound {
        Some(Sign::Positive)
    } else if determinant < -bound {
        Some(Sign::Negative)
    } else {
        None
    }
}
