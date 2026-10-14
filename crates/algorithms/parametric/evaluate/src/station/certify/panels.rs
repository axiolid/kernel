//! Panels of a certified length, refined where they are widest (#294,
//! #307).
//!
//! A certified length is a sum of panels, each bounded from both sides.
//! [`adapt`] keeps every panel in a queue by its width and halves the
//! widest until their widths sum to the panel share of the tolerance
//! ([`PANEL_SHARE`]), so the panels are fine only where the bound is loose
//! and the budget of [`MAX_CERTIFIED_PANELS`] goes where it is needed.
//!
//! # A chord's deficit
//!
//! A curve `P` measured by chords over `[a, b]`, `h = b - a`, midpoint
//! `m`, falls short of its arc length `L` by `D = L - c = int |P'| (1 -
//! cos theta)`, `theta` the angle of `P'` to the chord's direction `u`.
//! With `g = P' - (P'.u) u`, the velocity across the chord,
//! `|P'| (1 - cos theta) = |g|^2 / (|P'| (1 + cos theta))`, and `g` has
//! zero mean over the panel (the chord is the integral of `P'`). So
//!
//! `int |g|^2 / (2 S) <= D <= int |g|^2 / (s (1 + cos 2 beta))`
//!
//! with `s <= |P'| <= S` on the panel and `beta` the turning of the
//! tangent from its midpoint's (the chord lies within `beta` of it, so
//! every tangent within `2 beta` of the chord). The integral of `|g|^2` is
//! bounded three ways:
//!
//! - by `M2^2 h^3 / 12`, `M2 >= |P''| >= |g'|`: for a zero-mean `g`,
//!   `int |g|^2 = 1 / (2 h) int int |g(x) - g(y)|^2`, and `|g(x) - g(y)| <=
//!   M2 |x - y|`;
//! - from above by `int |g - g(m)|^2 <= int (alpha |t| + M3 t^2 / 2)^2 =
//!   alpha^2 h^3 / 12 + alpha M3 h^4 / 32 + M3^2 h^5 / 320`, `M3 >=
//!   |P'''|` and `alpha >= |g'(m)|` (a zero-mean `g` is closest to zero
//!   among its shifts);
//! - from below by Cauchy-Schwarz against `t`: `int |g|^2 >= |int t g|^2 /
//!   (h^3 / 12)`, with `|int t g| >= |g'(m)| h^3 / 12 - M3 h^4 / 64`.
//!
//! `g'(m) = P''(m)` across the chord is read from `G = (P'(b) - P'(a)) /
//! h`, within `M3 h / 4` of `P''(m)`, so no further derivative is
//! evaluated; the computed chord direction's own error and the velocities'
//! rounding widen it. The bounds from `M3` are a full order tighter: a
//! panel's width falls as `h^4`, not `h^3`, which is what lets an ordinary
//! curve certify to `1e-8` and below within the budget. The older bound
//! `c / cos(beta)` is kept as a ceiling, so no panel is wider than it was.
//!
//! `M3` needs `P''` continuous on the panel, so a caller cuts its span
//! where it may not be (`continuity_breaks*(.., 2)`, and `3` beside an
//! offset, whose `P''` holds the base's third derivative).

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use axiolid_contracts::GeomResult;
use axiolid_core::{Point3, Scalar, Vec3};

use super::{budget, uncertified, DistanceInterval, MAX_CERTIFIED_PANELS, POINT_ROUNDING};

/// The relative rounding allowed a computed velocity: generous, since it
/// widens a panel only through `h^2` times it.
pub(super) const VELOCITY_ROUNDING: Scalar = 1e-9;

/// The share of a tolerance the panels' own widths are refined to: the
/// rest is left to the points' rounding, and the widths come out narrower
/// than the uniform refinement before #307 left them.
pub(super) const PANEL_SHARE: Scalar = 0.125;

/// Suprema of a function and its derivatives over a panel, `[|f|, |f'|,
/// ...]`.
pub(super) type Jet<const N: usize> = [Scalar; N];

/// `a b`, zero when either is: a derivative known to vanish does not
/// meet an unbounded one.
fn times(a: Scalar, b: Scalar) -> Scalar {
    if a == 0.0 || b == 0.0 {
        0.0
    } else {
        a * b
    }
}

/// Leibniz: suprema of `(f g)^(k)`.
pub(super) fn product<const N: usize>(f: Jet<N>, g: Jet<N>) -> Jet<N> {
    let binomial = |n: usize, k: usize| -> Scalar {
        (0..k).fold(1.0, |acc, i| acc * (n - i) as Scalar / (i + 1) as Scalar)
    };
    core::array::from_fn(|k| {
        (0..=k)
            .map(|j| binomial(k, j) * times(f[j], g[k - j]))
            .sum()
    })
}

/// `[1, |u'|, |u''|, |u'''|]` of a unit vector `u = w / |w|`, from a lower
/// bound `r` of `|w|` and suprema `w1`, `w2`, `w3` of `|w'|`, `|w''|`,
/// `|w'''|`; unbounded where `r` is not positive.
///
/// With `x = w1 / r`, `y = w2 / r`, `z = w3 / r`: `u' = P w' / |w|`, `P =
/// I - u u^T`, so `|u'| <= x`. Differentiating `u'' = (P w'' + P' w' - u'
/// s') / s` (`s = |w|`, `|P'| <= 2 |u'|`, `|s'| <= |w'|`) gives `|u''| <=
/// y + 3 x^2`. Once more, `u''' = (P w''' + 2 P' w'' + P'' w' - 2 u'' s' -
/// u' s'') / s` with `|P''| <= 2 |u''| + 2 |u'|^2` and `|s''| <= |w'|^2 /
/// s + |w''|` gives `|u'''| <= z + 9 x y + 15 x^3`.
pub(super) fn unit_jet(r: Scalar, w1: Scalar, w2: Scalar, w3: Scalar) -> Jet<4> {
    if r.is_nan() || r <= 0.0 {
        return [1.0, Scalar::INFINITY, Scalar::INFINITY, Scalar::INFINITY];
    }
    let (x, y, z) = (w1 / r, w2 / r, w3 / r);
    [1.0, x, y + 3.0 * x * x, z + 9.0 * x * y + 15.0 * x * x * x]
}

/// A unit vector `w / |w|` and its derivative, from `w` and `w'`.
pub(super) fn unit_rate(w: Vec3, dw: Vec3) -> GeomResult<(Vec3, Vec3)> {
    let r = w.length();
    if r.is_nan() || r <= 0.0 {
        return Err(uncertified());
    }
    let u = w / r;
    Ok((u, (dw - u * u.dot(dw)) / r))
}

// --- sums ---------------------------------------------------------------------

/// A sum of many positive terms, compensated (Neumaier): its error is a
/// few units in the last place of the sum, whatever the number of terms.
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct Compensated {
    sum: Scalar,
    carry: Scalar,
}

impl Compensated {
    pub(super) fn add(&mut self, value: Scalar) {
        let total = self.sum + value;
        self.carry += if self.sum.abs() >= value.abs() {
            (self.sum - total) + value
        } else {
            (value - total) + self.sum
        };
        self.sum = total;
    }

    pub(super) fn value(self) -> Scalar {
        self.sum + self.carry
    }
}

/// A compensated sum of panel bounds as an interval, widened for its own
/// rounding; refused when wider than `tolerance`.
pub(super) fn finish(
    lower: Scalar,
    upper: Scalar,
    tolerance: Scalar,
) -> GeomResult<DistanceInterval> {
    let rounding = 4.0 * Scalar::EPSILON;
    let out = DistanceInterval::new(
        (lower * (1.0 - rounding)).max(0.0),
        upper * (1.0 + rounding),
    );
    if !(out.lower.is_finite() && out.upper.is_finite()) || out.width() > tolerance {
        return Err(uncertified());
    }
    Ok(out)
}

// --- adaptive refinement ----------------------------------------------------------

/// A queued panel: its width, then its slot.
#[derive(Debug, Clone, Copy)]
struct Queued(Scalar, usize);

impl PartialEq for Queued {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Queued {}

impl PartialOrd for Queued {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Queued {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0).then(other.1.cmp(&self.1))
    }
}

/// Panels refined by [`adapt`], in no particular order, and whether the
/// budget ran out before their widths reached the target.
#[derive(Debug, Clone)]
pub(super) struct Refined<P> {
    pub(super) panels: Vec<P>,
    pub(super) exhausted: bool,
}

/// Halve the widest of `initial` panels (by `width`, an unbounded or NaN
/// width first) through `split` until the widths sum to `target`, or the
/// panel budget is spent.
///
/// # Errors
///
/// More initial panels than the budget, and `split`'s own.
pub(super) fn adapt<P>(
    initial: Vec<P>,
    target: Scalar,
    width: impl Fn(&P) -> Scalar,
    mut split: impl FnMut(&P) -> GeomResult<(P, P)>,
) -> GeomResult<Refined<P>> {
    let mut count = initial.len();
    if count > MAX_CERTIFIED_PANELS {
        return Err(budget());
    }
    let measure = |p: &P| {
        let w = width(p);
        if w.is_nan() {
            Scalar::INFINITY
        } else {
            w
        }
    };
    let mut slots: Vec<Option<P>> = Vec::with_capacity(2 * count + 16);
    let mut queue = BinaryHeap::with_capacity(2 * count + 16);
    let (mut total, mut unbounded) = (0.0, 0_usize);
    let enqueue = |p: P,
                   slots: &mut Vec<Option<P>>,
                   queue: &mut BinaryHeap<Queued>,
                   total: &mut Scalar,
                   unbounded: &mut usize| {
        let w = measure(&p);
        if w.is_finite() {
            *total += w;
        } else {
            *unbounded += 1;
        }
        queue.push(Queued(w, slots.len()));
        slots.push(Some(p));
    };
    for p in initial {
        enqueue(p, &mut slots, &mut queue, &mut total, &mut unbounded);
    }
    let mut exhausted = false;
    loop {
        if unbounded == 0 && total <= target {
            // The running total drifts with its subtractions: recount it.
            let mut exact = Compensated::default();
            for p in slots.iter().flatten() {
                exact.add(measure(p));
            }
            total = exact.value();
            if total <= target {
                break;
            }
        }
        if count + 2 > MAX_CERTIFIED_PANELS {
            exhausted = true;
            break;
        }
        let Some(Queued(w, slot)) = queue.pop() else {
            break;
        };
        let Some(parent) = slots[slot].take() else {
            continue;
        };
        if w.is_finite() {
            total -= w;
        } else {
            unbounded -= 1;
        }
        let (left, right) = split(&parent)?;
        count += 2;
        enqueue(left, &mut slots, &mut queue, &mut total, &mut unbounded);
        enqueue(right, &mut slots, &mut queue, &mut total, &mut unbounded);
    }
    Ok(Refined {
        panels: slots.into_iter().flatten().collect(),
        exhausted,
    })
}

// --- chords -------------------------------------------------------------------------

/// A curve measured by chords: its point and its velocity at a native
/// parameter, each with an absolute error beyond rounding, and suprema of
/// `|P''|` and, where bounded, `|P'''|` over a parameter interval; and the
/// coordinate magnitude a point's rounding is relative to.
pub(super) struct Measured<'a> {
    pub(super) point: &'a dyn Fn(Scalar) -> GeomResult<(Point3, Scalar)>,
    pub(super) velocity: &'a dyn Fn(Scalar) -> GeomResult<(Vec3, Scalar)>,
    pub(super) bounds: &'a dyn Fn(Scalar, Scalar) -> GeomResult<(Scalar, Option<Scalar>)>,
    pub(super) scale: Scalar,
}

/// A point of a measured curve with its velocity, each with its error
/// beyond rounding.
#[derive(Debug, Clone, Copy)]
pub(super) struct Node {
    pub(super) t: Scalar,
    pub(super) point: Point3,
    pub(super) error: Scalar,
    pub(super) velocity: Vec3,
    pub(super) slip: Scalar,
}

/// A chord panel: its ends, the computed chord, and the arc's excess over
/// the true chord, between `low` and `high`.
#[derive(Debug, Clone, Copy)]
pub(super) struct Panel {
    pub(super) a: Node,
    pub(super) b: Node,
    pub(super) chord: Scalar,
    pub(super) low: Scalar,
    pub(super) high: Scalar,
}

impl Panel {
    pub(super) fn width(&self) -> Scalar {
        self.high - self.low
    }
}

impl Measured<'_> {
    pub(super) fn node(&self, t: Scalar) -> GeomResult<Node> {
        let (point, error) = (self.point)(t)?;
        let (velocity, slip) = (self.velocity)(t)?;
        Ok(Node {
            t,
            point,
            error,
            velocity,
            slip,
        })
    }

    /// The bound of a point's error: its rounding and its own error.
    pub(super) fn point_error(&self, node: &Node) -> Scalar {
        POINT_ROUNDING * self.scale.max(node.point.abs().max_element()).max(1.0) + node.error
    }

    /// The panel from `a` to `b`; see the [module documentation](self).
    pub(super) fn panel(&self, a: Node, b: Node) -> GeomResult<Panel> {
        let h = b.t - a.t;
        let mid = 0.5 * (a.t + b.t);
        let chord = (b.point - a.point).length();
        let (second, third) = (self.bounds)(a.t, b.t)?;
        let (velocity, slip) = (self.velocity)(mid)?;
        let speed = velocity.length();
        // Every tangent on the panel lies within `turn` of the midpoint's,
        // so the chord is at least `cos(turn)` of the arc.
        let floor = speed - slip - 0.5 * h * second;
        let turn = if floor > 0.0 {
            0.5 * h * second / floor + 1e-12
        } else {
            Scalar::INFINITY
        };
        // `c / cos(turn) - c`, without cancellation.
        let mut high = if turn < 1.0 {
            let half = (0.5 * turn).sin();
            chord * (1.0 + 4.0 * Scalar::EPSILON) * 2.0 * half * half / turn.cos()
        } else {
            Scalar::INFINITY
        };
        let mut low = 0.0;
        if 2.0 * turn < 1.0 {
            let h3 = h * h * h / 12.0;
            let mut upper = second * second * h3;
            let mut lower: Scalar = 0.0;
            let e = self.point_error(&a).max(self.point_error(&b));
            if let Some(third) = third.filter(|m| m.is_finite()) {
                if chord > 4.0 * e {
                    let u = (b.point - a.point) / chord;
                    let drift = 4.0 * e / (chord - 2.0 * e);
                    let g = (b.velocity - a.velocity) / h;
                    let error = 0.25 * third * h + (a.slip + b.slip) / h + 2.0 * g.length() * drift;
                    let across = (g - u * g.dot(u)).length();
                    let alpha = across + error;
                    let h4 = h3 * 12.0 * h;
                    upper = upper.min(
                        alpha * alpha * h3
                            + alpha * third * h4 / 32.0
                            + third * third * h4 * h / 320.0,
                    );
                    let moment = (across - error) * h3 - third * h4 / 64.0;
                    if moment > 0.0 {
                        lower = moment * moment / h3;
                    }
                }
            }
            let ceiling = speed + slip + 0.5 * h * second;
            let deficit = upper / (floor * (1.0 + (2.0 * turn).cos())) * (1.0 + 1e-12);
            if deficit.is_finite() {
                high = high.min(deficit);
            }
            low = (lower / (2.0 * ceiling) * (1.0 - 1e-12)).max(0.0);
            if low > high {
                low = 0.0;
            }
        }
        Ok(Panel {
            a,
            b,
            chord,
            low,
            high,
        })
    }

    /// Chord panels over `[t0, t1]` cut at `breaks`, refined until their
    /// widths sum to `target`, in order along the curve.
    ///
    /// # Errors
    ///
    /// A panel that cannot be halved further, and the evaluators' own.
    pub(super) fn chords(
        &self,
        t0: Scalar,
        t1: Scalar,
        breaks: &[Scalar],
        target: Scalar,
    ) -> GeomResult<Refined<Panel>> {
        let span = t1 - t0;
        if span.is_nan() || span <= 0.0 {
            return Ok(Refined {
                panels: Vec::new(),
                exhausted: false,
            });
        }
        let mut cuts = vec![t0];
        let mut inside: Vec<Scalar> = breaks
            .iter()
            .copied()
            .filter(|&t| t > t0 && t < t1)
            .collect();
        inside.sort_by(Scalar::total_cmp);
        inside.dedup();
        cuts.extend(inside);
        cuts.push(t1);
        let nodes = cuts
            .iter()
            .map(|&t| self.node(t))
            .collect::<GeomResult<Vec<_>>>()?;
        let initial = nodes
            .windows(2)
            .map(|pair| self.panel(pair[0], pair[1]))
            .collect::<GeomResult<Vec<_>>>()?;
        let mut refined = adapt(initial, target, Panel::width, |p| {
            let mid = 0.5 * (p.a.t + p.b.t);
            if !(mid > p.a.t && mid < p.b.t) {
                return Err(uncertified());
            }
            let m = self.node(mid)?;
            Ok((self.panel(p.a, m)?, self.panel(m, p.b)?))
        })?;
        refined.panels.sort_by(|p, q| p.a.t.total_cmp(&q.a.t));
        Ok(refined)
    }

    /// The certified length of the curve over `[t0, t1]` by chords (see
    /// the [module documentation](self)), split at `breaks` and refined
    /// adaptively within `tolerance`.
    ///
    /// # Errors
    ///
    /// An interval wider than `tolerance` (by
    /// [`super::UNCERTIFIED_LENGTH`], or the budget's name where the
    /// panels ran out first), and the evaluators' own.
    pub(super) fn chord_turning(
        &self,
        t0: Scalar,
        t1: Scalar,
        breaks: &[Scalar],
        tolerance: Scalar,
    ) -> GeomResult<DistanceInterval> {
        let refined = self.chords(t0, t1, breaks, PANEL_SHARE * tolerance)?;
        if refined.panels.is_empty() {
            return Ok(DistanceInterval::point(0.0));
        }
        let (lower, upper) = self.chain(&refined.panels);
        finish(lower, upper, tolerance).map_err(
            |error| {
                if refined.exhausted {
                    budget()
                } else {
                    error
                }
            },
        )
    }

    /// The bounds of the length of a chain of chord panels, in order:
    /// their chords and excesses, widened for the points' rounding.
    pub(super) fn chain(&self, panels: &[Panel]) -> (Scalar, Scalar) {
        let mut rounding = Rounding::new(self, panels);
        let (mut chords, mut lows, mut highs) = (
            Compensated::default(),
            Compensated::default(),
            Compensated::default(),
        );
        for panel in panels {
            chords.add(panel.chord);
            lows.add(panel.low);
            highs.add(panel.high);
            rounding.panel(panel);
        }
        let error = rounding.bound();
        let sum = chords.value();
        (sum + lows.value() - error, sum + highs.value() + error)
    }
}

/// The rounding of a chain of chords between computed points (#294).
///
/// Each point is off its curve by at most `e`: its rounding,
/// [`POINT_ROUNDING`] relative to the largest coordinate magnitude, plus
/// its own error (#307: an offset by distances read through its base's
/// certified length). Writing a chord `|D + d|`, `D` the true one and `d`
/// the difference of its ends' errors (`|d| <= 2e`): `|D + d| - |D| = d.u +
/// q`, `u = D / |D|`, with `0 <= q <= |d|^2 / (2 |D|)`. The parts `d.u`
/// telescope along the chain: their sum is at most `e (2 + sum |u_i -
/// u_(i-1)|)`, the directions' total variation, which the chords
/// themselves bound. So the rounding of the whole chain stays a few `e`,
/// however many chords it has; each chord's own norm adds `4 eps |c|`.
#[derive(Debug, Clone, Copy)]
struct Rounding {
    e: Scalar,
    previous: Option<Vec3>,
    variation: Scalar,
    quadratic: Scalar,
    norms: Scalar,
}

impl Rounding {
    fn new(curve: &Measured<'_>, panels: &[Panel]) -> Self {
        let e = panels
            .iter()
            .flat_map(|p| [curve.point_error(&p.a), curve.point_error(&p.b)])
            .fold(POINT_ROUNDING * curve.scale.max(1.0), Scalar::max);
        Self {
            e,
            previous: None,
            variation: 0.0,
            quadratic: 0.0,
            norms: 0.0,
        }
    }

    fn panel(&mut self, panel: &Panel) {
        let (e, chord) = (self.e, panel.chord);
        if chord > 4.0 * e {
            let direction = (panel.b.point - panel.a.point) / chord;
            if let Some(previous) = self.previous {
                self.variation += (direction - previous).length();
            }
            self.previous = Some(direction);
            // The chord's second-order part, and the error of the computed
            // direction against the true one (at most `4 e / |D|`, twice).
            self.quadratic += 10.0 * e * e / (chord - 2.0 * e);
        } else {
            // Too short to have a direction: count its whole error.
            self.variation += 4.0;
            self.previous = None;
        }
        self.norms += 4.0 * Scalar::EPSILON * chord;
    }

    /// A bound on the difference of the computed and the true chord sums.
    fn bound(&self) -> Scalar {
        self.e * (2.0 + self.variation) + self.quadratic + self.norms
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A unit vector's derivative bounds hold for a curve whose `|w|`
    /// dips: the third by differences of the exact second.
    #[test]
    fn the_unit_vector_bounds_hold() {
        let w = |t: Scalar| Vec3::new(1.0 + 0.5 * t.cos(), 0.8 * t.sin(), 0.3 * t);
        let dw = |t: Scalar| Vec3::new(-0.5 * t.sin(), 0.8 * t.cos(), 0.3);
        let unit = |t: Scalar| {
            let v = w(t);
            v / v.length()
        };
        let r = 0.45;
        let jet = unit_jet(r, 1.0, 1.0, 1.0);
        let h = 1e-3;
        for i in 0..2000 {
            let t = Scalar::from(i) * 3e-3;
            let (_, d) = unit_rate(w(t), dw(t)).unwrap();
            assert!(d.length() <= jet[1]);
            let (_, up) = unit_rate(w(t + h), dw(t + h)).unwrap();
            let (_, down) = unit_rate(w(t - h), dw(t - h)).unwrap();
            assert!(((up - down) / (2.0 * h)).length() <= jet[2]);
            let third = (unit(t + 2.0 * h) - 2.0 * unit(t + h) + 2.0 * unit(t - h)
                - unit(t - 2.0 * h))
                / (2.0 * h * h * h);
            assert!(
                third.length() <= jet[3],
                "{t}: {} > {}",
                third.length(),
                jet[3]
            );
        }
    }

    #[test]
    fn a_compensated_sum_of_many_terms_stays_exact() {
        let mut sum = Compensated::default();
        for _ in 0..1_000_000 {
            sum.add(0.1);
        }
        assert!((sum.value() - 100_000.0).abs() <= 1e-9);
    }

    /// The widest panel is halved first, and the refinement stops once the
    /// widths sum to the target.
    #[test]
    fn the_widest_panel_is_halved_first() {
        // A panel `[a, b]` whose width is `(b - a)^2` on the left half of
        // `[0, 1]` and `(b - a)^3` on the right.
        let width = |p: &(Scalar, Scalar)| {
            let h = p.1 - p.0;
            if p.1 <= 0.5 {
                h * h
            } else {
                h * h * h
            }
        };
        let refined = adapt(vec![(0.0, 0.5), (0.5, 1.0)], 1e-3, width, |p| {
            let mid = 0.5 * (p.0 + p.1);
            Ok(((p.0, mid), (mid, p.1)))
        })
        .unwrap();
        assert!(!refined.exhausted);
        let total: Scalar = refined.panels.iter().map(width).sum();
        assert!(total <= 1e-3);
        let left = refined.panels.iter().filter(|p| p.1 <= 0.5).count();
        let right = refined.panels.len() - left;
        assert!(left > 4 * right, "{left} {right}");
    }
}
