//! Certified positions of intrinsic curves (#317).
//!
//! An intrinsic curve's point is the integral of its unit tangent `f =
//! (cos theta, sin theta)` in the start frame, `theta` the law's closed-form
//! heading. [`IntrinsicPositions`] integrates it by the composite midpoint
//! rule and bounds the error of every point it returns, so a polyline
//! through such points can carry a bound rather than the evaluator's
//! quadrature estimate (`crate::arc_length`).
//!
//! # The remainder
//!
//! On a panel of width `w` with midpoint `m`, `|int f - w f(m)| <= w^3 / 24
//! sup |f''|` (the vector Taylor remainder). `f'' = theta'' f^perp -
//! theta'^2 f`, so `|f''| <= |k'| + k^2`, bounded per seam-free piece of
//! the law from its closed-form suprema (`crate::chain::curvature_bound`),
//! as the certified station lengths bound their speeds (#294). Panels never
//! straddle a seam, where `k` may jump. The panel width is chosen so the
//! remainders over the whole curve sum to half the requested error; each
//! point reports the remainders actually summed up to it, plus:
//!
//! - the heading's rounding, at most `64 eps (1 + K L) (1 + seams)` radians
//!   (`K` bounding `|k|`, `L` the length), which turns each tangent by as
//!   much, so a point at `s` moves by at most `s` times it;
//! - the sums' rounding, `16 eps s` (compensated sums);
//! - the start frame's stretch `sigma`, applied to all of these, and
//!   `8 eps (|origin| + sigma s)` for placing the point.
//!
//! A planar space curve (torsion identically zero) is the same integral in
//! its start frame's `x`, `y` plane: its tangent is `cos theta x + sin
//! theta y`.
//!
//! Points are marched forward from the last one read, so reading them in
//! increasing arc length costs one pass; reading backwards starts again.

use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Vec2};
use axiolid_curve::{Intrinsic2, Intrinsic3};

use super::{frame_stretch2, frame_stretch3, inflate};
use crate::chain::curvature_bound;

/// Most midpoint panels one curve's positions may use.
pub const MAX_POSITION_PANELS: usize = 1 << 22;

/// Where a marched point is placed: by a 2D curve's start frame, or in the
/// plane of a planar space curve's start frame.
#[derive(Debug, Clone, Copy)]
enum Placement {
    Plane(Frame2),
    Space(Frame3),
}

/// Certified positions of an intrinsic curve: the unit tangent of its
/// closed-form heading integrated by the composite midpoint rule, each
/// point with a certified bound on its error (the remainder `w^3 / 24
/// (|k'| + k^2)` per panel, panels split at seams, and the heading's, the
/// sums' and the placement's rounding).
#[derive(Debug, Clone)]
pub struct IntrinsicPositions {
    /// The curve in its own start frame.
    local: Intrinsic2,
    placement: Placement,
    stretch: Scalar,
    origin: Scalar,
    /// Seam-free pieces `(lo, hi, sup |f''|)`.
    pieces: Vec<(Scalar, Scalar, Scalar)>,
    step: Scalar,
    heading_rounding: Scalar,
    at: Scalar,
    sum: [(Scalar, Scalar); 2],
    remainder: Scalar,
}

/// Neumaier's compensated addition into `(sum, carry)`.
fn add(acc: &mut (Scalar, Scalar), value: Scalar) {
    let total = acc.0 + value;
    acc.1 += if acc.0.abs() >= value.abs() {
        (acc.0 - total) + value
    } else {
        (value - total) + acc.0
    };
    acc.0 = total;
}

fn unit_frame() -> Frame2 {
    Frame2 {
        origin: Point2::ZERO,
        x: Vec2::X,
        y: Vec2::Y,
    }
}

impl IntrinsicPositions {
    /// Positions of a 2D intrinsic curve within about `error`, half of it
    /// spent on the midpoint remainders; `None` for a non-positive or
    /// non-finite error or length, a law with no closed-form bound, or more
    /// than [`MAX_POSITION_PANELS`] panels.
    #[must_use]
    pub fn new2(curve: &Intrinsic2, error: Scalar) -> Option<Self> {
        let local = Intrinsic2::new(unit_frame(), curve.curvature.clone(), curve.length);
        Self::build(
            local,
            Placement::Plane(curve.start),
            frame_stretch2(&curve.start),
            curve.start.origin.abs().max_element(),
            error,
        )
    }

    /// Positions of a planar space curve (its torsion law identically
    /// zero, [`Intrinsic3::is_planar`]); `None` as for [`Self::new2`], or
    /// for a curve with torsion.
    #[must_use]
    pub fn new3(curve: &Intrinsic3, error: Scalar) -> Option<Self> {
        if !curve.is_planar() {
            return None;
        }
        let local = Intrinsic2::new(unit_frame(), curve.curvature.clone(), curve.length);
        Self::build(
            local,
            Placement::Space(curve.start),
            frame_stretch3(&curve.start),
            curve.start.origin.abs().max_element(),
            error,
        )
    }

    fn build(
        local: Intrinsic2,
        placement: Placement,
        stretch: Scalar,
        origin: Scalar,
        error: Scalar,
    ) -> Option<Self> {
        let length = local.length;
        if !(error.is_finite() && error > 0.0 && length.is_finite() && length > 0.0) {
            return None;
        }
        let law = &local.curvature;
        let rate = law.derivative();
        let seams = law.seams_within(length);
        let mut ends = vec![0.0];
        ends.extend(seams.iter().copied());
        ends.push(length);
        let mut pieces = Vec::with_capacity(ends.len() - 1);
        let mut bend: Scalar = 0.0;
        for pair in ends.windows(2) {
            let (lo, hi) = (pair[0], pair[1]);
            let k = curvature_bound(law, lo, hi)?;
            let m = inflate(curvature_bound(&rate, lo, hi)? + k * k);
            bend = bend.max(m);
            pieces.push((lo, hi, m));
        }
        let budget = 0.5 * error / stretch;
        let step = if bend > 0.0 {
            (24.0 * budget / (length * bend)).sqrt().min(length)
        } else {
            length
        };
        if !(step.is_finite() && step > 0.0) || length / step > MAX_POSITION_PANELS as Scalar {
            return None;
        }
        let turning = curvature_bound(law, 0.0, length)? * length;
        let heading_rounding =
            64.0 * Scalar::EPSILON * (1.0 + turning) * (1.0 + seams.len() as Scalar);
        Some(Self {
            local,
            placement,
            stretch,
            origin,
            pieces,
            step,
            heading_rounding,
            at: 0.0,
            sum: [(0.0, 0.0); 2],
            remainder: 0.0,
        })
    }

    /// March the local integral to `s`.
    fn march(&mut self, s: Scalar) -> Option<()> {
        if !(s.is_finite() && (0.0..=self.local.length).contains(&s)) {
            return None;
        }
        if s < self.at {
            self.at = 0.0;
            self.sum = [(0.0, 0.0); 2];
            self.remainder = 0.0;
        }
        for &(lo, hi, bend) in &self.pieces {
            let (a, b) = (lo.max(self.at), hi.min(s));
            if b <= a {
                continue;
            }
            let n = ((b - a) / self.step).ceil().max(1.0);
            if n > MAX_POSITION_PANELS as Scalar {
                return None;
            }
            let w = (b - a) / n;
            for i in 0..n as usize {
                let theta = self.local.heading_at(a + (i as Scalar + 0.5) * w)?;
                add(&mut self.sum[0], w * theta.cos());
                add(&mut self.sum[1], w * theta.sin());
            }
            self.remainder += n * w * w * w / 24.0 * bend;
        }
        self.at = s;
        Some(())
    }

    /// The marched point in the start frame's coordinates, and its error
    /// there.
    fn local_point(&self) -> (Point2, Scalar) {
        let p = Point2::new(self.sum[0].0 + self.sum[0].1, self.sum[1].0 + self.sum[1].1);
        let error = self.remainder + self.at * (self.heading_rounding + 16.0 * Scalar::EPSILON);
        (p, error)
    }

    /// A local error through the start frame, with the placement's rounding.
    fn placed_error(&self, local: Scalar) -> Scalar {
        inflate(
            self.stretch * local + 8.0 * Scalar::EPSILON * (self.origin + self.stretch * self.at),
        )
    }

    /// The point at arc length `s` of a 2D curve and a certified bound on
    /// its distance from the exact point; `None` outside `[0, length]`, for
    /// a law that does not integrate there, or for a space curve.
    pub fn point2(&mut self, s: Scalar) -> Option<(Point2, Scalar)> {
        let Placement::Plane(frame) = self.placement else {
            return None;
        };
        self.march(s)?;
        let (p, local) = self.local_point();
        let out = Point2::new(
            frame.origin.x + frame.x.x * p.x + frame.y.x * p.y,
            frame.origin.y + frame.x.y * p.x + frame.y.y * p.y,
        );
        Some((out, self.placed_error(local)))
    }

    /// [`Self::point2`] for a planar space curve.
    pub fn point3(&mut self, s: Scalar) -> Option<(Point3, Scalar)> {
        let Placement::Space(frame) = self.placement else {
            return None;
        };
        self.march(s)?;
        let (p, local) = self.local_point();
        Some((
            frame.origin + frame.x * p.x + frame.y * p.y,
            self.placed_error(local),
        ))
    }

    /// The point at arc length `s` of the 2D curve's offset by `distance`
    /// to its left (the curve of [`super::offset_point2`]) and a certified
    /// bound on its distance from the exact one: the point's, plus `|d|`
    /// times the turn of the left normal by the heading's rounding through
    /// the start frame (its stretch squared over its determinant bounds its
    /// condition) and the normalisation's rounding.
    pub fn offset_point2(&mut self, s: Scalar, distance: Scalar) -> Option<(Point2, Scalar)> {
        let Placement::Plane(frame) = self.placement else {
            return None;
        };
        let (p, error) = self.point2(s)?;
        let theta = self.local.heading_at(s)?;
        let (c, n) = (theta.cos(), theta.sin());
        let t = Vec2::new(frame.x.x * c + frame.y.x * n, frame.x.y * c + frame.y.y * n);
        let size = t.length();
        let det = (frame.x.x * frame.y.y - frame.x.y * frame.y.x).abs();
        if !(size.is_finite() && size > 0.0 && det > 0.0 && distance.is_finite()) {
            return None;
        }
        let left = Vec2::new(-t.y, t.x) / size;
        let condition = self.stretch * self.stretch / det;
        let turn = condition * (self.heading_rounding + 16.0 * Scalar::EPSILON);
        Some((
            p + left * distance,
            inflate(error + distance.abs() * (turn + 8.0 * Scalar::EPSILON)),
        ))
    }
}
