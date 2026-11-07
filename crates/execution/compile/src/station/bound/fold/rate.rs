//! A lower bound of the rate along the turned section plane over a span
//! (#406).
//!
//! The section point `p = q + x y^ + y e` (components in the BASE frame
//! `E = [T, L, U]`, `q` the offsets) advances along the plane's normal `n` at
//!
//! `F = sigma n_T + rise n_U + q'.n + w.(p x n) + x c_x + y c_y`
//!
//! (see [`turned`](super::turned) and the [module above](super)), `w =
//! (tau, -beta_U, beta_L)` the frame's Darboux vector. Over a sub-span every
//! factor is an interval:
//!
//! - the rates from `FrameLaw::rates`: a certified range, or for the
//!   witness the middle value `+-` its error;
//! - the turned frame in closed form at the span's middle, widened by half
//!   the span times the supremum of its derivative (`Turn`);
//! - `q`, `x` and `y` between their values at the span's ends (they are
//!   linear), and `p` from those where the orientation is constant, else
//!   from the interval sum.
//!
//! The lower bound of `F`, less the profile's flattening bound times the
//! gradient of `F` in `(x, y)` and the rounding of the arithmetic, is a
//! certified floor of the vertex's rate over the span.

use axiolid_core::{Scalar, Vec3};
use axiolid_model::{StationOffsets, StationOrientation};
use axiolid_reference::station::FrameRates;

use super::turned::{turned_frame, TurnedFrame};
use super::{lerp_offsets, turn, SectionLaw, ROUNDING};
use crate::station::lerp;

/// A closed interval of reals.
#[derive(Debug, Clone, Copy)]
struct Range {
    lo: Scalar,
    hi: Scalar,
}

impl Range {
    fn point(x: Scalar) -> Self {
        Self { lo: x, hi: x }
    }

    fn around(middle: Scalar, radius: Scalar) -> Self {
        Self {
            lo: middle - radius,
            hi: middle + radius,
        }
    }

    fn hull(p: Scalar, q: Scalar) -> Self {
        Self {
            lo: p.min(q),
            hi: p.max(q),
        }
    }

    fn from(pair: [Scalar; 2]) -> Self {
        Self {
            lo: pair[0],
            hi: pair[1],
        }
    }

    fn add(self, other: Self) -> Self {
        Self {
            lo: self.lo + other.lo,
            hi: self.hi + other.hi,
        }
    }

    fn mul(self, other: Self) -> Self {
        let products = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        Self {
            lo: products.iter().copied().fold(Scalar::INFINITY, Scalar::min),
            hi: products
                .iter()
                .copied()
                .fold(Scalar::NEG_INFINITY, Scalar::max),
        }
    }

    fn magnitude(self) -> Scalar {
        self.lo.abs().max(self.hi.abs())
    }
}

/// The five rates of the frame law a span's floors are built from.
#[derive(Debug, Clone, Copy)]
pub(super) struct Coefficients {
    forward: Range,
    rise: Range,
    lateral: Range,
    up: Range,
    twist: Range,
}

impl Coefficients {
    /// The certified ranges over the span: what clears it.
    pub(super) fn ranges(rates: &FrameRates) -> Self {
        let range = |r: &axiolid_reference::station::RateBounds| Range::from(r.range);
        Self {
            forward: range(&rates.forward),
            rise: range(&rates.rise),
            lateral: range(&rates.lateral),
            up: range(&rates.up),
            twist: range(&rates.twist),
        }
    }

    /// The middle values `+-` their errors: what refuses.
    pub(super) fn middles(rates: &FrameRates) -> Self {
        let around = |r: &axiolid_reference::station::RateBounds| Range::around(r.middle, r.error);
        Self {
            forward: around(&rates.forward),
            rise: around(&rates.rise),
            lateral: around(&rates.lateral),
            up: around(&rates.up),
            twist: around(&rates.twist),
        }
    }
}

/// A vertex's lower bound of the rate along the turned plane's normal.
#[derive(Debug, Clone, Copy)]
pub(super) struct Floor {
    /// The certified lower bound.
    pub(super) floor: Scalar,
    /// The middle of the interval the rate lies in, before the margin and the
    /// rounding are taken off.
    pub(super) centre: Scalar,
    /// The vertex `(T, L, U)` in the base frame, at the span's middle.
    pub(super) point: [Scalar; 3],
}

/// The terms of `F` as `(sign, factors)`, the lower bound of each product
/// summed with its sign.
struct Sum {
    lo: Scalar,
    hi: Scalar,
    magnitude: Scalar,
}

impl Sum {
    fn new() -> Self {
        Self {
            lo: 0.0,
            hi: 0.0,
            magnitude: 0.0,
        }
    }

    fn term(&mut self, sign: Scalar, factors: &[Range]) {
        let product = factors[1..]
            .iter()
            .fold(factors[0], |acc, factor| acc.mul(*factor));
        self.magnitude += factors.iter().map(|f| f.magnitude()).product::<Scalar>();
        if sign > 0.0 {
            self.lo += product.lo;
            self.hi += product.hi;
        } else {
            self.lo -= product.hi;
            self.hi -= product.lo;
        }
    }
}

impl SectionLaw<'_, '_> {
    /// Each vertex's floor of the rate along the turned plane's normal over
    /// the distances `[lo, hi]` of station interval `i`, from `c`, less
    /// `margin` times the rate's gradient in the section's plane; `None`
    /// where the orientation's turn cannot be bounded.
    pub(super) fn floors(
        &self,
        i: usize,
        (lo, hi): (Scalar, Scalar),
        c: &Coefficients,
        margin: Scalar,
    ) -> Option<Vec<Floor>> {
        let (d0, d1) = (self.distances[i], self.distances[i + 1]);
        let delta = d1 - d0;
        let fraction = |s: Scalar| ((s - d0) / delta).clamp(0.0, 1.0);
        let (u0, u1) = (fraction(lo), fraction(hi));
        let unit = |o: &StationOrientation| o.unit_axes().ok();
        let (ta, tb) = (
            unit(&self.orientations[i])?,
            unit(&self.orientations[i + 1])?,
        );
        let turn = turn(ta, tb, u0, u1, delta)?;
        let both_base = self.orientations[i].is_base() && self.orientations[i + 1].is_base();
        let frame = if both_base {
            TurnedFrame::BASE
        } else {
            turned_frame(ta, tb, 0.5 * (u0 + u1), delta)?
        };
        // The orientation does not change over the interval.
        let constant = ta == tb;
        let half = 0.5 * (hi - lo);
        let spread = |v: Vec3, rate: Scalar| {
            let radius = half * rate + ROUNDING;
            [
                Range::around(v.x, radius),
                Range::around(v.y, radius),
                Range::around(v.z, radius),
            ]
        };
        let rate = |bound: Scalar| if both_base { 0.0 } else { bound };
        let n = spread(frame.normal, rate(turn.tangent[0]));
        let l = spread(frame.lateral, rate(turn.lateral[0]));
        let e = spread(frame.up, rate(turn.up[0]));
        let c_x = Range::around(
            frame.lateral_rate,
            half * (turn.lateral[1] + turn.lateral[0] * turn.tangent[0]) + frame.rate_error[0],
        );
        let c_y = Range::around(
            frame.up_rate,
            half * (turn.up[1] + turn.up[0] * turn.tangent[0]) + frame.rate_error[1],
        );
        let (a, b) = (&self.offsets[i], &self.offsets[i + 1]);
        let (at_lo, at_hi) = (lerp_offsets(a, b, u0), lerp_offsets(a, b, u1));
        // The offsets `q = (T, L, U)` and their rates.
        let components = |o: &StationOffsets| [o.longitudinal, o.lateral, o.vertical];
        let q = [0, 1, 2].map(|k| Range::hull(components(&at_lo)[k], components(&at_hi)[k]));
        let rates = [0, 1, 2].map(|k| (components(b)[k] - components(a)[k]) / delta);
        // The gradient of `F` in `(x, y)`: `G = -w.e + c_x`, `H = w.y^ + c_y`.
        let w = [c.twist, Range::from([-c.up.hi, -c.up.lo]), c.lateral];
        let dot = |u: &[Range; 3]| {
            w.iter()
                .zip(u)
                .fold(Range::point(0.0), |acc, (w, u)| acc.add(w.mul(*u)))
        };
        let (w_e, w_l) = (dot(&e), dot(&l));
        let g = Range::from([-w_e.hi, -w_e.lo]).add(c_x);
        let h = w_l.add(c_y);
        let reach = margin * g.magnitude().hypot(h.magnitude());

        let mut floors = Vec::with_capacity(self.locals[i].len());
        for (p, r) in self.locals[i].iter().zip(&self.locals[i + 1]) {
            let (x_lo, x_hi) = (lerp(p.x, r.x, u0), lerp(p.x, r.x, u1));
            let (y_lo, y_hi) = (lerp(p.y, r.y, u0), lerp(p.y, r.y, u1));
            let (x, y) = (Range::hull(x_lo, x_hi), Range::hull(y_lo, y_hi));
            // The vertex in the base frame.
            let point: [Range; 3] = if constant {
                // Linear in the distance: from its ends.
                let at = |o: &StationOffsets, x: Scalar, y: Scalar| {
                    let q = components(o);
                    [0, 1, 2].map(|k| q[k] + x * frame.lateral[k] + y * frame.up[k])
                };
                let (start, end) = (at(&at_lo, x_lo, y_lo), at(&at_hi, x_hi, y_hi));
                [0, 1, 2].map(|k| Range::hull(start[k], end[k]))
            } else {
                [0, 1, 2].map(|k| q[k].add(x.mul(l[k])).add(y.mul(e[k])))
            };
            let mut sum = Sum::new();
            sum.term(1.0, &[c.forward, n[0]]);
            sum.term(1.0, &[c.rise, n[2]]);
            for k in 0..3 {
                sum.term(1.0, &[Range::point(rates[k]), n[k]]);
            }
            // `w.(p x n)`, `w = (tau, -beta_U, beta_L)`.
            sum.term(1.0, &[c.twist, point[1], n[2]]);
            sum.term(-1.0, &[c.twist, point[2], n[1]]);
            sum.term(-1.0, &[c.up, point[2], n[0]]);
            sum.term(1.0, &[c.up, point[0], n[2]]);
            sum.term(1.0, &[c.lateral, point[0], n[1]]);
            sum.term(-1.0, &[c.lateral, point[1], n[0]]);
            sum.term(1.0, &[x, c_x]);
            sum.term(1.0, &[y, c_y]);
            let floor = sum.lo - reach - ROUNDING * (sum.magnitude + reach);
            floors.push(Floor {
                floor,
                centre: 0.5 * (sum.lo + sum.hi),
                point: point.map(|range| 0.5 * (range.lo + range.hi)),
            });
        }
        Some(floors)
    }
}
