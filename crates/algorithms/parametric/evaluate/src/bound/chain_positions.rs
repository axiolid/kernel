//! Certified positions of an arc-length chain (#438).
//!
//! The evaluator reads a chain ([`Chain2`]) piece by piece: an intrinsic
//! piece by quadrature, a parametric piece at the parameter of an
//! arc-length inverse whose `1e-12` tolerance is an estimate, and every
//! piece is placed at the end of the one before it, so each error is carried
//! on. [`ChainPositionTable`] reads the *exact* chain with a certified bound
//! on the distance of every point, tangent and curvature vector it returns,
//! so a chord bound can charge the evaluated ends against it
//! ([`chain_chord_bound`](crate::chain::chain_chord_bound), #414's rule for
//! chains).
//!
//! # The exact chain
//!
//! It is the representation's chain, nothing is re-fitted:
//! `placement_0` is the identity and `placement_{k+1} = placement_k o
//! (q_k(L_k), unit q_k'(L_k))`, with `q_k` the exact local curve of piece
//! `k`: the integral of `e^(i theta)` for an intrinsic piece, `c(u*(m))`
//! with `lambda(u*) = m` measured from `start` for a parametric one. A
//! parametric curve may start up to `PIECE_FRAME_TOLERANCE max(1, L)` off its
//! local origin, so the exact chain may jump by that at a join: a read
//! carries a [`SeamSide`], and `Incoming` at a join reads the piece before
//! it at its end.
//!
//! # A local read
//!
//! A piece at local distance `m` gives a point `p~` within `e` of the exact
//! `p*`, a unit tangent `t~` within `phi` and a curvature vector `w~` (`p''`
//! in arc length) within `omega`:
//!
//! - **Intrinsic.** [`IntrinsicPositionTable`] of the piece in its identity
//!   frame, built at the piece's share of the budget; the tangent is the
//!   closed-form heading, off by its `heading_rounding`; `w~ = k(m)
//!   perp(t~)` with `k` the law's value, off by `K (Phi_k + rho)`, `K` the
//!   law's supremum over the piece and `Phi_k` the turn error of the
//!   placement.
//! - **Line.** `u = start + m / |d|`: a distance off by `delta = |(u -
//!   start) |d| - m|` plus a few roundings, which moves a point of a unit
//!   speed curve by at most `delta`. Straight: no curvature.
//! - **Circle, ellipse, B-spline.** The parameter is *enclosed*, not
//!   trusted: [`CurveLengths`](crate::station::CurveLengths) bounds
//!   `lambda(u)` between `lo` and `hi` by chord panels, the evaluator's
//!   inverse is only a guess, and a guess whose enclosure misses `m` is
//!   bisected until it holds `m` (a piece longer than its curve is `None`).
//!   `delta = max(m - lo, hi - m, 0)` bounds `|lambda(u) - m|`, which is the
//!   arc between `u` and the exact `u*`, which bounds `|c(u) - c(u*)|` and
//!   turns the unit tangent by at most `K delta`, `K` the curvature supremum
//!   over the panels the piece reaches; the curvature vector moves by `K'
//!   delta`, `K'` bounding `|p'''|`. A read whose `delta` exceeds the
//!   piece's share is `None`, so `max_error` is a true supremum by
//!   refusal.
//!
//! # Carried errors
//!
//! With `E_k` the position error and `Phi_k` the turn error of piece `k`'s
//! placement (`E_0 = Phi_0 = 0`): a point `P = O_k + X_k p~` is off by `E_k +
//! e + Phi_k (|p~| + e)` plus the rounding of the sum, because `|R(X) p~ -
//! R(X*) p*| <= |p~ - p*| + |X - X*| |p*|`; and the next placement carries
//! `E_{k+1} = E_k + e_end + Phi_k (|p~_end| + e_end) + 8 eps (|O_k| +
//! |p~_end|)`, `Phi_{k+1} = Phi_k + phi_end + 8 eps`, the angles of
//! successive turns adding. The start frame maps all of it, stretching by at
//! most [`frame_stretch2`]. Each piece's share of the budget is `error / (2
//! sigma n)`, `n` the pieces up to the end, so the sum of the pieces'
//! own errors is half the request; what the carried turns add on top is
//! rounding.
//!
//! `max_error`, `max_tangent_error` and `max_curvature_error` are suprema
//! over every read: a read never reports more (it refuses first).

use axiolid_core::{Frame2, Point2, Scalar, Vec2};
use axiolid_curve::{Chain2, ChainPiece2, Curve2, Intrinsic2};

use super::interval::Iv;
use super::{frame_stretch2, inflate, IntrinsicPositionTable};
use crate::arc_length::intrinsic_tangent;
use crate::arc_parameter::{parameter_at_arc_length2, ARC_LENGTH_TOLERANCE};
use crate::chain::{
    check_local_start, check_well_formed, curvature_bound, local_frame, parametric_bend,
    parametric_plan_jet, unit, Placement,
};
use crate::curve::{derivative2, domain2, evaluate2};
use crate::station::{CurveLengths, POINT_ROUNDING};

/// The relative rounding of a computed derivative against the curve's
/// derivative bound, as [`POINT_ROUNDING`] is against coordinates.
const DERIVATIVE_ROUNDING: Scalar = 1e-14;

/// Most halvings the search for the enclosed parameter makes.
const BISECTIONS: usize = 128;

/// The rounding of one placement step, in units of epsilon.
const PLACE: Scalar = 8.0 * Scalar::EPSILON;

/// Which piece a read at a join belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeamSide {
    /// The piece that starts at the join (what the evaluator reads).
    Outgoing,
    /// The piece that ends at the join: its end, before any jump a
    /// parametric piece's start makes.
    Incoming,
}

/// One local read of a piece.
#[derive(Debug, Clone, Copy)]
struct Local {
    p: Vec2,
    e: Scalar,
    t: Vec2,
    phi: Scalar,
    w: Vec2,
    /// The curvature vector's error before the placement's turn is added.
    omega: Scalar,
}

/// How a piece is read.
#[derive(Debug, Clone)]
enum Kind {
    Intrinsic {
        local: Intrinsic2,
        table: IntrinsicPositionTable,
    },
    Parametric {
        curve: Curve2,
        start: Scalar,
        /// The parameter the evaluator reads at the piece end.
        end: Scalar,
        /// A line's `|d|`; `None` for a curve measured by its panels.
        rate: Option<Scalar>,
        lengths: Option<CurveLengths>,
        high: Scalar,
        /// `K'`, bounding `|p'''|`.
        jerk: Scalar,
        /// A parameter magnitude the line's rounding is relative to.
        reach: Scalar,
    },
}

/// A piece with its placement and the errors carried into it.
#[derive(Debug, Clone)]
struct Placed {
    length: Scalar,
    placement: Placement,
    position: Scalar,
    turn: Scalar,
    share: Scalar,
    /// `K`, bounding `|p''|` (the curvature) over the piece.
    curvature: Scalar,
    /// The tangent's rounding.
    rho: Scalar,
    /// The largest error of a local point.
    emax: Scalar,
    kind: Kind,
}

/// Certified positions, unit tangents and curvature vectors of an
/// arc-length chain over `[0, end]`, read in any order; see the module
/// documentation.
#[derive(Debug, Clone)]
pub struct ChainPositionTable {
    start: Frame2,
    stretch: Scalar,
    total: Scalar,
    end: Scalar,
    /// Where each built piece starts, summed as `Chain2::piece_at` sums.
    begins: Vec<Scalar>,
    pieces: Vec<Placed>,
    max_point: Scalar,
    max_tangent: Scalar,
    max_curvature: Scalar,
}

/// The world read of a piece at a local distance.
struct Read {
    point: Point2,
    point_error: Scalar,
    tangent: Vec2,
    tangent_error: Scalar,
    curvature: Vec2,
    curvature_error: Scalar,
}

impl ChainPositionTable {
    /// The table over `[0, end]`, its points within about `error`; `None`
    /// for a malformed chain (or a parametric piece off its local frame), an
    /// error or end that is not finite and positive, an end past the chain,
    /// an intrinsic piece whose law has no closed-form bound of `k` to
    /// `k'''`, a parametric piece of a family other than a line, a circle,
    /// an ellipse or a B-spline, a piece longer than its curve, a curve
    /// without a positive speed floor on a panel, or a budget exhausted.
    #[must_use]
    pub fn new(chain: &Chain2, error: Scalar, end: Scalar) -> Option<Self> {
        if !(error.is_finite() && error > 0.0 && end.is_finite() && end > 0.0) {
            return None;
        }
        check_well_formed(chain).ok()?;
        let total = chain.length()?;
        if end > total + ARC_LENGTH_TOLERANCE * total.max(1.0) {
            return None;
        }
        let end = end.min(total);
        let stretch = frame_stretch2(&chain.start);
        let mut begins = vec![0.0];
        for piece in &chain.pieces {
            let last = *begins.last()?;
            begins.push(last + piece.length());
        }
        let built = begins[..chain.pieces.len()]
            .iter()
            .take_while(|&&begin| begin <= end)
            .count()
            .max(1);
        begins.truncate(built);
        let budget = 0.5 * error / (stretch * built as Scalar);
        let mut pieces: Vec<Placed> = Vec::with_capacity(built);
        let mut placement = Placement::IDENTITY;
        let (mut position, mut turn) = (0.0, 0.0);
        for (index, piece) in chain.pieces[..built].iter().enumerate() {
            let placed = Placed::new(index, piece, budget, placement, position, turn)?;
            // The end of the piece: where the next one is placed.
            let local = placed.read(placed.length)?;
            let end_turn = turn + local.phi + PLACE;
            position = position
                + local.e
                + turn * (local.p.length() + local.e)
                + PLACE * (placement.origin.abs().max_element() + local.p.length());
            turn = end_turn;
            placement = Placement {
                origin: placement.place(local.p),
                x: unit(placement.turn(local.t), "chain join").ok()?,
            };
            pieces.push(placed);
        }
        let mut table = Self {
            start: chain.start,
            stretch,
            total,
            end,
            begins,
            pieces,
            max_point: 0.0,
            max_tangent: 0.0,
            max_curvature: 0.0,
        };
        table.bound_suprema();
        let finite = table.max_point.is_finite()
            && table.max_tangent.is_finite()
            && table.max_curvature.is_finite();
        finite.then_some(table)
    }

    /// The suprema of the three errors over every read.
    fn bound_suprema(&mut self) {
        let origin = self.start.origin.abs().max_element();
        let mut worst = (0.0, 0.0, 0.0);
        for piece in &self.pieces {
            let size = piece.length * (1.0 + PIECE_SLACK) + PIECE_SLACK + piece.emax;
            let corner = piece.placement.origin.abs().max_element();
            let local = piece.position
                + piece.emax
                + piece.turn * (size + piece.emax)
                + PLACE * (corner + size);
            let placed = (piece.placement.origin.length() + size).max(1.0);
            let point = self.stretch * local + PLACE * (origin + self.stretch * placed);
            let tangent = self.stretch * (piece.turn + piece.phi_max() + PLACE);
            let omega = piece.omega_max() + piece.curvature * piece.turn;
            let curvature = self.stretch * (omega + PLACE * piece.curvature);
            worst.0 = Scalar::max(worst.0, point);
            worst.1 = Scalar::max(worst.1, tangent);
            worst.2 = Scalar::max(worst.2, curvature);
        }
        self.max_point = inflate(worst.0);
        self.max_tangent = inflate(worst.1);
        self.max_curvature = inflate(worst.2);
    }

    /// The piece holding `s` on `side`, and the local distance in it.
    fn locate(&self, s: Scalar, side: SeamSide) -> Option<(&Placed, Scalar)> {
        let slack = ARC_LENGTH_TOLERANCE * self.total.max(1.0);
        // A distance a rounding step past the end is the end, as the
        // evaluator reads it.
        let s = if s > self.total && s <= self.total + slack {
            self.total
        } else {
            s
        };
        if !(s.is_finite() && s >= 0.0 && s <= self.end) {
            return None;
        }
        let index = self
            .begins
            .partition_point(|&begin| begin <= s)
            .saturating_sub(1)
            .min(self.pieces.len() - 1);
        if side == SeamSide::Incoming && index > 0 && s == self.begins[index] {
            let piece = &self.pieces[index - 1];
            return Some((piece, piece.length));
        }
        let piece = &self.pieces[index];
        Some((piece, (s - self.begins[index]).clamp(0.0, piece.length)))
    }

    fn read(&self, s: Scalar, side: SeamSide) -> Option<Read> {
        let (piece, m) = self.locate(s, side)?;
        let local = piece.read(m)?;
        let origin = self.start.origin.abs().max_element();
        let turned = piece.placement.turn(local.p);
        let p = piece.placement.origin + turned;
        let size = local.p.length();
        let at = piece.position
            + local.e
            + piece.turn * (size + local.e)
            + PLACE * (piece.placement.origin.abs().max_element() + size);
        let world = |v: Vec2| {
            let x = self.start.x * v.x + self.start.y * v.y;
            Vec2::new(x.x, x.y)
        };
        let w = world(p);
        let point_error =
            inflate(self.stretch * at + PLACE * (origin + self.stretch * w.length().max(1.0)));
        let omega = local.omega + piece.curvature * piece.turn;
        let out = Read {
            point: Point2::new(self.start.origin.x + w.x, self.start.origin.y + w.y),
            point_error,
            tangent: world(piece.placement.turn(local.t)),
            tangent_error: inflate(self.stretch * (piece.turn + local.phi + PLACE)),
            curvature: world(piece.placement.turn(local.w)),
            curvature_error: inflate(self.stretch * (omega + PLACE * piece.curvature)),
        };
        (point_error <= self.max_point
            && out.tangent_error <= self.max_tangent
            && out.curvature_error <= self.max_curvature)
            .then_some(out)
    }

    /// The point at arc length `s` as the evaluator reads it (a join belongs
    /// to the piece that starts there), and a certified bound on its
    /// distance from the exact point. `None` outside `[0, end]`, or where a
    /// piece refuses the read.
    #[must_use]
    pub fn point2(&self, s: Scalar) -> Option<(Point2, Scalar)> {
        self.point2_on(s, SeamSide::Outgoing)
    }

    /// [`Self::point2`] on `side` of a join.
    #[must_use]
    pub fn point2_on(&self, s: Scalar, side: SeamSide) -> Option<(Point2, Scalar)> {
        self.read(s, side).map(|r| (r.point, r.point_error))
    }

    /// The tangent at `s` on `side`, and a certified bound on its distance
    /// from the exact tangent.
    #[must_use]
    pub fn tangent2_on(&self, s: Scalar, side: SeamSide) -> Option<(Vec2, Scalar)> {
        self.read(s, side).map(|r| (r.tangent, r.tangent_error))
    }

    /// `p''(s)` in arc length (`k n`) on `side`, and a certified bound on its
    /// distance from the exact one. Where a piece's curvature law has a seam
    /// inside, the value is the left one, as the evaluator reads it.
    #[must_use]
    pub fn curvature2_on(&self, s: Scalar, side: SeamSide) -> Option<(Vec2, Scalar)> {
        self.read(s, side).map(|r| (r.curvature, r.curvature_error))
    }

    /// A certified bound on the error of every point the table returns.
    #[must_use]
    pub fn max_error(&self) -> Scalar {
        self.max_point
    }

    /// A certified bound on the error of every tangent the table returns.
    #[must_use]
    pub fn max_tangent_error(&self) -> Scalar {
        self.max_tangent
    }

    /// A certified bound on the error of every curvature vector the table
    /// returns.
    #[must_use]
    pub fn max_curvature_error(&self) -> Scalar {
        self.max_curvature
    }
}

/// Certified bound on the distance from a chain over `[a, b]` (inside one
/// piece) to the chord of its evaluated ends; see
/// [`chain_chord_bound`](crate::chain::chain_chord_bound).
///
/// The exact-chord bound plus the larger distance of an evaluated end from
/// the table's point (and that point's error): the lower end read as the
/// piece that starts there, the upper end as the piece that ends there. The
/// table is built at an eighth of the exact bound, or a few ulps of the
/// placement if more.
pub(crate) fn evaluated_chain_chord_bound(chain: &Chain2, a: Scalar, b: Scalar) -> Option<Scalar> {
    let (lo, hi) = super::ordered(a, b);
    let law = crate::chain::exact_chain_chord_bound(chain, lo, hi)?;
    let total = chain.length()?;
    let stretch = frame_stretch2(&chain.start);
    let origin = chain.start.origin.abs().max_element();
    let floor = 64.0 * Scalar::EPSILON * (1.0 + origin + stretch * hi);
    let table = ChainPositionTable::new(chain, (law / 8.0).max(floor), hi.min(total))?;
    let gap = |s: Scalar, side: SeamSide| {
        let (exact, error) = table.point2_on(s.min(total), side)?;
        let read = crate::chain::chain_point(chain, s).ok()?;
        let d = [read.x - exact.x, read.y - exact.y]
            .map(|c| Iv::point(c).sq())
            .into_iter()
            .fold(Iv::point(0.0), Iv::add)
            .sqrt();
        Some(Iv::point(d.hi).add(Iv::point(error)).hi)
    };
    let g = gap(lo, SeamSide::Outgoing)?.max(gap(hi, SeamSide::Incoming)?);
    let out = inflate(Iv::point(law).add(Iv::point(g)).hi);
    out.is_finite().then_some(out)
}

/// How far a parametric piece's curve may sit from its local origin, per
/// `max(1, length)`: [`crate::chain::PIECE_FRAME_TOLERANCE`].
const PIECE_SLACK: Scalar = crate::chain::PIECE_FRAME_TOLERANCE;

enum Side {
    Below,
    Inside(Scalar, Scalar),
    Above,
}

impl Placed {
    fn new(
        index: usize,
        piece: &ChainPiece2,
        budget: Scalar,
        placement: Placement,
        position: Scalar,
        turn: Scalar,
    ) -> Option<Self> {
        let length = piece.length();
        match piece {
            ChainPiece2::Intrinsic { curvature, .. } => {
                let local = Intrinsic2::new(local_frame(), curvature.clone(), length);
                let table = IntrinsicPositionTable::new2(&local, budget, length)?;
                let k = curvature_bound(curvature, 0.0, length)?;
                Some(Self {
                    length,
                    placement,
                    position,
                    turn,
                    share: budget,
                    curvature: k,
                    rho: table.heading_rounding(),
                    emax: table.max_error(),
                    kind: Kind::Intrinsic { local, table },
                })
            }
            ChainPiece2::Parametric {
                curve,
                start,
                length,
            } => {
                check_local_start(index, curve, *start, *length).ok()?;
                // At least what the evaluator's own inverse is good for.
                let share = budget.max(64.0 * ARC_LENGTH_TOLERANCE * length.max(1.0));
                let end = parameter_at_arc_length2(curve, *start, *length).ok()?;
                match curve {
                    Curve2::Line(line) => {
                        let rate = line.direction.length();
                        if !(rate > 0.0 && rate.is_finite()) {
                            return None;
                        }
                        let o = line.origin.abs().max_element();
                        let d = line.direction.abs().max_element();
                        let reach = (o + d * (start.abs() + end.abs())).max(1.0);
                        Some(Self {
                            length: *length,
                            placement,
                            position,
                            turn,
                            share,
                            curvature: 0.0,
                            rho: 2.0 * DERIVATIVE_ROUNDING,
                            emax: share + POINT_ROUNDING * reach,
                            kind: Kind::Parametric {
                                curve: curve.clone(),
                                start: *start,
                                end,
                                rate: Some(rate),
                                lengths: None,
                                high: end,
                                jerk: 0.0,
                                reach,
                            },
                        })
                    }
                    Curve2::Circle(_) | Curve2::Ellipse(_) | Curve2::BSpline(_) => {
                        let high = match curve {
                            Curve2::BSpline(_) => {
                                let domain = domain2(curve);
                                domain.start.max(domain.end)
                            }
                            _ => start + core::f64::consts::TAU,
                        };
                        if high.partial_cmp(start) != Some(core::cmp::Ordering::Greater) {
                            return None;
                        }
                        let lengths = CurveLengths::new(curve, *start, high, 0.25 * share).ok()?;
                        // Even its least length falls short of the piece:
                        // the piece is longer than its curve.
                        if lengths.total().1 < *length {
                            return None;
                        }
                        let (mut k, mut jerk, mut speed, mut first): (
                            Scalar,
                            Scalar,
                            Scalar,
                            Scalar,
                        ) = (0.0, 0.0, Scalar::INFINITY, 0.0);
                        for (a, b) in lengths.panels() {
                            if lengths.at(a).ok()?.0 > *length + share {
                                break;
                            }
                            let (p2, p3, sigma) = parametric_plan_jet(curve, a, b)?;
                            let bounds = crate::bound::curve_derivative_bounds2(curve, a, b)?;
                            k = k.max(p2);
                            jerk = jerk.max(p3);
                            speed = speed.min(sigma);
                            first = first.max(bounds.first);
                        }
                        if !(k.is_finite() && jerk.is_finite() && speed > 0.0 && first > 0.0) {
                            return None;
                        }
                        let rho = 2.0 * DERIVATIVE_ROUNDING * first / speed;
                        let extent = lengths.extent();
                        Some(Self {
                            length: *length,
                            placement,
                            position,
                            turn,
                            share,
                            curvature: k,
                            rho,
                            emax: share + POINT_ROUNDING * extent,
                            kind: Kind::Parametric {
                                curve: curve.clone(),
                                start: *start,
                                end,
                                rate: None,
                                lengths: Some(lengths),
                                high,
                                jerk,
                                reach: extent,
                            },
                        })
                    }
                    _ => None,
                }
            }
            _ => None,
        }
    }

    /// The tangent's largest error in a local read.
    fn phi_max(&self) -> Scalar {
        match &self.kind {
            Kind::Intrinsic { .. } => self.rho,
            Kind::Parametric { rate: Some(_), .. } => self.rho,
            Kind::Parametric { .. } => self.curvature * self.share + self.rho,
        }
    }

    /// The curvature vector's largest error in a local read, before the
    /// placement's turn.
    fn omega_max(&self) -> Scalar {
        match &self.kind {
            Kind::Intrinsic { .. } => self.curvature * self.rho,
            Kind::Parametric { rate: Some(_), .. } => 0.0,
            Kind::Parametric { jerk, .. } => {
                self.curvature * self.rho
                    + jerk * self.share
                    + 4.0 * DERIVATIVE_ROUNDING * self.curvature
            }
        }
    }

    /// The local read at distance `m` of the piece.
    fn read(&self, m: Scalar) -> Option<Local> {
        if !(m.is_finite() && (0.0..=self.length).contains(&m)) {
            return None;
        }
        match &self.kind {
            Kind::Intrinsic { local, table } => {
                let (p, e) = table.point2(m)?;
                let t = intrinsic_tangent(local, m).ok()?;
                let kappa = crate::frenet::law_value(&local.curvature, m)?;
                let out = Local {
                    p: Vec2::new(p.x, p.y),
                    e,
                    t,
                    phi: self.rho,
                    w: Vec2::new(-t.y, t.x) * kappa,
                    omega: self.curvature * self.rho,
                };
                (e <= self.emax && out.w.is_finite()).then_some(out)
            }
            Kind::Parametric {
                curve,
                start,
                end,
                rate,
                lengths,
                high,
                reach,
                ..
            } => {
                let (u, delta) = match (rate, lengths) {
                    (Some(rate), _) => {
                        let u = if m == self.length {
                            *end
                        } else {
                            start + m / rate
                        };
                        let delta = ((u - start) * rate - m).abs()
                            + 4.0 * Scalar::EPSILON * (rate * (u.abs() + start.abs()) + m);
                        (u, delta)
                    }
                    (None, Some(lengths)) => {
                        let guess = if m == self.length {
                            *end
                        } else {
                            parameter_at_arc_length2(curve, *start, m).ok()?
                        };
                        solve(lengths, *start, *high, m, guess)?
                    }
                    (None, None) => return None,
                };
                if !(delta.is_finite() && delta <= self.share && u.is_finite()) {
                    return None;
                }
                let p = evaluate2(curve, u).ok()?;
                let t = unit(derivative2(curve, u).ok()?, "chain piece curve").ok()?;
                let rounding = POINT_ROUNDING * reach;
                let e = delta + rounding;
                let (w, omega) = if rate.is_some() {
                    (Vec2::ZERO, 0.0)
                } else {
                    (
                        parametric_bend(curve, u).ok()?,
                        self.curvature * self.rho
                            + self.jerk() * delta
                            + 4.0 * DERIVATIVE_ROUNDING * self.curvature,
                    )
                };
                let out = Local {
                    p: Vec2::new(p.x, p.y),
                    e,
                    t,
                    phi: if rate.is_some() {
                        self.rho
                    } else {
                        self.curvature * delta + self.rho
                    },
                    w,
                    omega,
                };
                (e <= self.emax && out.p.is_finite() && out.w.is_finite()).then_some(out)
            }
        }
    }

    fn jerk(&self) -> Scalar {
        match &self.kind {
            Kind::Parametric { jerk, .. } => *jerk,
            Kind::Intrinsic { .. } => 0.0,
        }
    }
}

/// The parameter `u` with `lambda(u)` enclosing `m`, and the width `delta`
/// of the enclosure about it; `None` where the piece is longer than the
/// curve or no enclosure is found.
///
/// `guess` is the evaluator's parameter, tried first and never trusted.
fn solve(
    lengths: &CurveLengths,
    start: Scalar,
    high: Scalar,
    m: Scalar,
    guess: Scalar,
) -> Option<(Scalar, Scalar)> {
    let classify = |u: Scalar| -> Option<Side> {
        let (lo, hi) = lengths.at(u).ok()?;
        Some(if hi < m {
            Side::Below
        } else if lo > m {
            Side::Above
        } else {
            Side::Inside(lo, hi)
        })
    };
    let inside = |lo: Scalar, hi: Scalar| (m - lo).max(hi - m).max(0.0);
    let guess = guess.clamp(start, high);
    let (mut below, mut above) = match classify(guess)? {
        Side::Inside(lo, hi) => return Some((guess, inside(lo, hi))),
        Side::Below => match classify(high)? {
            Side::Below => return None,
            Side::Inside(lo, hi) => return Some((high, inside(lo, hi))),
            Side::Above => (guess, high),
        },
        Side::Above => match classify(start)? {
            Side::Above => return None,
            Side::Inside(lo, hi) => return Some((start, inside(lo, hi))),
            Side::Below => (start, guess),
        },
    };
    for _ in 0..BISECTIONS {
        let mid = 0.5 * (below + above);
        if !(mid > below && mid < above) {
            return None;
        }
        match classify(mid)? {
            Side::Inside(lo, hi) => return Some((mid, inside(lo, hi))),
            Side::Below => below = mid,
            Side::Above => above = mid,
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::Point2;
    use axiolid_curve::{BSplineCurve2, KnotSpec};

    /// The cubic `y = x^3 / (6 R X)` of radius 300 over 60 m, as a Bezier.
    fn cubic() -> Curve2 {
        let (reach, a) = (60.0, 1.0 / (6.0 * 300.0 * 60.0));
        Curve2::BSpline(BSplineCurve2 {
            degree: 3,
            control_points: vec![
                Point2::new(0.0, 0.0),
                Point2::new(reach / 3.0, 0.0),
                Point2::new(2.0 * reach / 3.0, 0.0),
                Point2::new(reach, a * reach * reach * reach),
            ],
            knots: vec![0.0, 1.0],
            multiplicities: vec![4, 4],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::PiecewiseBezier,
        })
    }

    #[test]
    fn a_parametric_read_off_its_distance_is_charged() {
        // A coarse enclosure holds both the guess's distance and `m`, so the
        // read is charged at least the gap between them.
        let curve = cubic();
        let lengths = CurveLengths::new(&curve, 0.0, 1.0, 1e-3).unwrap();
        let m = 30.0;
        let off = parameter_at_arc_length2(&curve, 0.0, m + 1e-6).unwrap();
        let (u, delta) = solve(&lengths, 0.0, 1.0, m, off).unwrap();
        assert_eq!(u, off);
        assert!(delta >= 1e-6 - 1e-9, "{delta}");
        let (lo, hi) = lengths.at(u).unwrap();
        assert!(lo <= m && m <= hi);
    }

    #[test]
    fn a_bad_guess_is_bisected() {
        let curve = cubic();
        let lengths = CurveLengths::new(&curve, 0.0, 1.0, 1e-9).unwrap();
        let m = 30.0;
        for far in [m + 1e-3, m + 10.0, m - 10.0, 59.0, 0.5] {
            let guess = parameter_at_arc_length2(&curve, 0.0, far).unwrap();
            let (u, delta) = solve(&lengths, 0.0, 1.0, m, guess).unwrap();
            let (lo, hi) = lengths.at(u).unwrap();
            assert!(lo <= m && m <= hi, "guess for {far}");
            assert!(delta <= 1e-9, "{delta}");
        }
        // Past the curve: no parameter has that length.
        assert!(solve(&lengths, 0.0, 1.0, 61.0, 1.0).is_none());
    }
}
