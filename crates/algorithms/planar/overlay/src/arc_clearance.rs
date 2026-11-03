//! Opt-in clearance between the edges of an arc polygon (#388).
//!
//! [`check_arc_polygon`](crate::check_arc_polygon) decides contact exactly
//! and with no tolerance: rings that miss each other by `2^-40` are valid. A
//! consumer that treats anything closer than a distance as touching (the
//! straight overlay's `validate_ring` does, at the linear tolerance) asks
//! for [`ArcCheckOptions::with_clearance`] and gets an
//! [`ArcPolygonDefect::Clearance`](crate::ArcPolygonDefect::Clearance) for
//! every pair of edges that comes closer than it, by the contract below.
//!
//! # What is reported
//!
//! A pair of edges is reported when
//!
//! - their distance has a local minimum strictly between zero and the
//!   clearance, away from the points where they meet, or
//! - they form an exact cusp: they meet at a point and leave it in the same
//!   direction (a zero interior angle, decided exactly, at a vertex or at a
//!   tangent touch inside both edges), so they stay within the clearance of
//!   each other for a stretch however small it is.
//!
//! Not reported: the neighbourhood of a meeting point at a nonzero angle. A
//! narrow wedge, two edges closing in on a shared vertex at a small angle,
//! has all its points near each other but no local minimum of the distance
//! other than the vertex itself, so it is found only once an end of one edge
//! comes within the clearance of the other. Edges that leave a joint in
//! opposite directions (a smooth joint, a rounded corner) are no cusp.
//!
//! # Which pairs
//!
//! Every pair of edges over all usable rings, in ring then edge order,
//! including consecutive edges and pairs that meet. A meeting point is not a
//! near miss, but it does not hide one: the distance tests below look for
//! critical points at a distance strictly between zero and the clearance, so
//! the joint two consecutive edges share and a vertex touching an edge are
//! left out as the zero-distance points they are (they are reported as
//! contacts and self-contacts by the exact check), while a part of the same
//! two edges that misses by less than the clearance elsewhere, or two edges
//! of a thin sliver, is found. The cusp test below adds the meeting points
//! whose edges leave along one direction.
//!
//! A ring that fails `validate_arc_ring` has no usable edges and is left
//! out. An edge shorter than the clearance, or two edges that close in on
//! each other within it, put their ends within it, and that is reported.
//!
//! # How a pair is decided
//!
//! 1. **Boxes.** The broad phase queries a box tree with every edge's box
//!    grown by the clearance. A pair whose boxes are `clearance` or more
//!    apart cannot come closer (`Bounds::gap` is a certified lower bound,
//!    shrunk past its own rounding) and is dropped.
//! 2. **Exact distance.** The shortest segment between two edges, left out
//!    of where they meet, ends at an end point of one of them, or runs
//!    between interior points of both along a common normal. So some pair
//!    of points is at a distance in `(0, t)` exactly when one of these
//!    holds, each with a distance that is not zero:
//!    - an end point of one edge is within `t` of the other edge: for a
//!      segment, the foot of the perpendicular if it lies on the segment,
//!      else the nearer end; for an arc, `| |X - c| - r |` if the ray from
//!      the centre through `X` crosses the arc, else the nearer end;
//!    - a segment and an arc: the foot of the perpendicular from the arc's
//!      centre lies on the segment, and the arc point on that perpendicular
//!      (`c +- r n`) lies on the arc and within `t` of the line;
//!    - two arcs with distinct centres: the arc points on the line of
//!      centres (`c1 +- r1 w`, `c2 +- r2 w`) lie on their arcs and are
//!      within `t` of each other. Concentric arcs are covered by the end
//!      point tests: overlapping angular ranges put an end of one in the
//!      other.
//!
//!    Every test is an exact sign. An arc's centre is `K / k` with `k = 4b`,
//!    its radius `|D| (1 + b^2) / |k|`; coordinates are scaled by `|k|` so no
//!    division is needed, and each distance comparison is the sign of a sum
//!    of at most three square roots of dyadics, decided by `axiolid-exact`'s
//!    tower, which also decides the exact zero that marks a meeting point.
//!    No `f64` rounding enters a verdict.
//! 3. **Cusps.** An end of one edge on the other edge, with the two edges
//!    leaving it in the same direction (`cross = 0`, `dot > 0` of exact
//!    leaving vectors; an arc's is the perpendicular of its centre-relative
//!    point), or a zero common-normal distance strictly inside both edges (a
//!    line tangent to an arc, two tangent arcs), is a cusp.
//! 4. **Fail closed.** A test the exact tier cannot decide (zero or not
//!    included) makes the pair a finding with `undecided` set; a pair is
//!    dropped only when every test is decided and none is below.
//!
//! A finding names the two edges and carries a certified lower bound on
//! their distance (the gap of their boxes), so a consumer can tell a pair
//! that is certainly closer than the clearance (`undecided` false) from one
//! reported only because it could not be cleared.

use axiolid_core::{Point2, Scalar};
use axiolid_exact::{Arith, Dyadic, Tower};
use axiolid_guarantees::Sign;

use crate::arc::ArcRing;
use crate::exact_arc::boxes::BoxTree;
use crate::exact_arc::edges_of;

/// Options of [`check_arc_polygon_with`](crate::check_arc_polygon_with).
///
/// The default is the behaviour of [`check_arc_polygon`](crate::check_arc_polygon):
/// contact is decided exactly, with no clearance.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[non_exhaustive]
pub struct ArcCheckOptions {
    /// Report edges that miss each other by less than this distance, or
    /// `None` (the default) to leave nearness unchecked.
    pub clearance: Option<Scalar>,
}

impl ArcCheckOptions {
    /// The default options: no clearance check.
    #[must_use]
    pub const fn new() -> Self {
        Self { clearance: None }
    }

    /// These options with a clearance check at `distance`.
    ///
    /// A distance that is not finite or is negative is reported by the check
    /// as [`ArcPolygonDefect::InvalidClearance`](crate::ArcPolygonDefect::InvalidClearance),
    /// never replaced; zero finds nothing, since only touching edges are at
    /// distance zero.
    #[must_use]
    pub const fn with_clearance(mut self, distance: Scalar) -> Self {
        self.clearance = Some(distance);
        self
    }
}

/// One pair of edges closer than the clearance, by flat ring index
/// (`0` the outer ring, `k + 1` hole `k`) and edge index.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct NearPair {
    pub(crate) edges: [(usize, usize); 2],
    pub(crate) lower_bound: f64,
    pub(crate) undecided: bool,
}

/// The pairs of edges of `rings` closer than `clearance` (finite, positive).
/// `usable[i]` says whether ring `i` may be tested.
pub(crate) fn near_pairs(rings: &[&ArcRing], usable: &[bool], clearance: f64) -> Vec<NearPair> {
    near_pairs_by(rings, usable, clearance, |a, b| a.closer_than(b, clearance))
}

/// [`near_pairs`] with the exact distance test given, so a test can make it
/// undecided.
fn near_pairs_by(
    rings: &[&ArcRing],
    usable: &[bool],
    clearance: f64,
    decide: impl Fn(&Shape, &Shape) -> Verdict,
) -> Vec<NearPair> {
    let mut edges = Vec::new();
    let mut shapes = Vec::new();
    let mut owner: Vec<(usize, usize)> = Vec::new();
    for (index, ring) in rings.iter().enumerate() {
        if !usable[index] {
            continue;
        }
        let n = ring.vertices.len();
        for (e, edge) in edges_of(ring).into_iter().enumerate() {
            let from = ring.vertices[e];
            let to = ring.vertices[(e + 1) % n];
            shapes.push(Shape::new(from.point, to.point, from.bulge));
            edges.push(edge);
            owner.push((index, e));
        }
    }
    let tree = BoxTree::new(edges.iter().map(|e| e.bounds).collect());
    let (mut hits, mut stack) = (Vec::new(), Vec::new());
    let mut found = Vec::new();
    for i in 0..edges.len() {
        let grown = edges[i].bounds.grown(clearance);
        tree.query_into(|b| b.overlaps(&grown), &mut hits, &mut stack);
        for &j in hits.iter().filter(|&&j| j > i) {
            let ((ra, ea), (rb, eb)) = (owner[i], owner[j]);
            let lower_bound = edges[i].bounds.gap(&edges[j].bounds);
            if lower_bound >= clearance {
                continue;
            }
            let undecided = match decide(&shapes[i], &shapes[j]) {
                Verdict::Clear => continue,
                Verdict::Below => false,
                Verdict::Undecided => true,
            };
            found.push(NearPair {
                edges: [(ra, ea), (rb, eb)],
                lower_bound,
                undecided,
            });
        }
    }
    found
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Verdict {
    /// The edges are certainly closer than the distance.
    Below,
    /// The edges are certainly not.
    Clear,
    /// An exact test could not decide.
    Undecided,
}

type V = (Dyadic, Dyadic);

fn dy(value: f64) -> Dyadic {
    Dyadic::from_f64(value)
}

fn point(p: Point2) -> V {
    (dy(p.x), dy(p.y))
}

fn sub(a: &V, b: &V) -> V {
    (a.0.sub(&b.0), a.1.sub(&b.1))
}

fn dot(a: &V, b: &V) -> Dyadic {
    a.0.mul(&b.0).add(&a.1.mul(&b.1))
}

fn cross(a: &V, b: &V) -> Dyadic {
    a.0.mul(&b.1).sub(&a.1.mul(&b.0))
}

fn scale(a: &V, by: &Dyadic) -> V {
    (a.0.mul(by), a.1.mul(by))
}

fn is_zero(a: &V) -> bool {
    sign(&a.0) == Sign::Zero && sign(&a.1) == Sign::Zero
}

/// Sign of an exact dyadic. A dyadic always decides.
fn sign(value: &Dyadic) -> Sign {
    value.sign().expect("a dyadic has a sign")
}

/// Whether two nonzero vectors point the same way.
fn same_direction(a: &V, b: &V) -> bool {
    sign(&cross(a, b)) == Sign::Zero && sign(&dot(a, b)) == Sign::Positive
}

fn at_least_zero(value: &Dyadic) -> bool {
    sign(value) != Sign::Negative
}

/// One edge as exact data.
#[derive(Debug, Clone)]
struct Shape {
    from: V,
    to: V,
    kind: Kind,
}

#[derive(Debug, Clone)]
enum Kind {
    Segment,
    Arc(Box<ArcShape>),
}

/// An arc in coordinates scaled by `ka = |k|`, `k = 4b`, and centred on its
/// own centre `c = K / k`: a point `Y` becomes `ka (Y - c) = ka Y - ks K`.
/// There the circle has radius `rho` with `rho^2 = r2`, so no division is
/// needed anywhere.
#[derive(Debug, Clone)]
struct ArcShape {
    ka: Dyadic,
    ks: Dyadic,
    k_centre: V,
    /// `rho^2 = |D|^2 (1 + b^2)^2`.
    r2: Dyadic,
    /// The ends in scaled coordinates, ordered counter-clockwise.
    start: V,
    end: V,
    major: bool,
}

impl Shape {
    fn new(from: Point2, to: Point2, bulge: f64) -> Self {
        let (f, t) = (point(from), point(to));
        if bulge == 0.0 {
            return Self {
                from: f,
                to: t,
                kind: Kind::Segment,
            };
        }
        let b = dy(bulge);
        let k = b.mul(&dy(4.0));
        let one_minus = dy(1.0).sub(&b.square());
        let d = sub(&t, &f);
        let half = dy(0.5);
        let mid = (f.0.add(&t.0).mul(&half), f.1.add(&t.1).mul(&half));
        let k_centre = (
            k.mul(&mid.0).sub(&one_minus.mul(&d.1)),
            k.mul(&mid.1).add(&one_minus.mul(&d.0)),
        );
        let one_plus = dy(1.0).add(&b.square());
        let r2 = dot(&d, &d).mul(&one_plus.square());
        let positive = sign(&k) == Sign::Positive;
        let (ka, ks) = if positive {
            (k, dy(1.0))
        } else {
            (k.neg(), dy(-1.0))
        };
        let mut arc = ArcShape {
            ka,
            ks,
            k_centre,
            r2,
            start: (dy(0.0), dy(0.0)),
            end: (dy(0.0), dy(0.0)),
            major: sign(&b.square().sub(&dy(1.0))) == Sign::Positive,
        };
        let (a, e) = (arc.scaled(&f), arc.scaled(&t));
        // A clockwise arc is the counter-clockwise one between the same ends,
        // read the other way.
        (arc.start, arc.end) = if positive { (a, e) } else { (e, a) };
        Self {
            from: f,
            to: t,
            kind: Kind::Arc(Box::new(arc)),
        }
    }

    /// Whether points of the edges are at a distance in `(0, t)`.
    fn closer_than(&self, other: &Self, t: f64) -> Verdict {
        let t = dy(t);
        let mut undecided = false;
        let mut hit = |result: Option<bool>| match result {
            Some(true) => true,
            Some(false) => false,
            None => {
                undecided = true;
                false
            }
        };
        // End point of one against the other edge.
        for end in [&self.from, &self.to] {
            if hit(other.near_point(end, &t)) {
                return Verdict::Below;
            }
        }
        for end in [&other.from, &other.to] {
            if hit(self.near_point(end, &t)) {
                return Verdict::Below;
            }
        }
        // Interior against interior, along a common normal.
        let interior = match (&self.kind, &other.kind) {
            (Kind::Segment, Kind::Segment) => Some(Dist::Far),
            (Kind::Segment, Kind::Arc(arc)) => segment_arc(self, arc, &t),
            (Kind::Arc(arc), Kind::Segment) => segment_arc(other, arc, &t),
            (Kind::Arc(a), Kind::Arc(b)) => arc_arc(a, b, &t),
        };
        // A zero here is a tangent touch of the two curves: a cusp.
        if hit(interior.map(|dist| dist != Dist::Far)) {
            return Verdict::Below;
        }
        // A meeting point at which the edges leave along one direction.
        if self.cusp_with(other) || other.cusp_with(self) {
            return Verdict::Below;
        }
        if undecided {
            Verdict::Undecided
        } else {
            Verdict::Clear
        }
    }

    /// Whether an end of this edge lies on `other` and both edges leave it in
    /// the same direction (zero interior angle, decided exactly).
    fn cusp_with(&self, other: &Self) -> bool {
        [&self.from, &self.to].into_iter().any(|end| {
            other.holds_point(end)
                && self.leaving(end).iter().any(|mine| {
                    other
                        .leaving(end)
                        .iter()
                        .any(|theirs| same_direction(mine, theirs))
                })
        })
    }

    /// Whether the exact point `x` lies on this edge, ends included.
    fn holds_point(&self, x: &V) -> bool {
        match &self.kind {
            Kind::Segment => {
                let d = sub(&self.to, &self.from);
                let w = sub(x, &self.from);
                let along = dot(&w, &d);
                sign(&cross(&d, &w)) == Sign::Zero
                    && at_least_zero(&along)
                    && at_least_zero(&dot(&d, &d).sub(&along))
            }
            Kind::Arc(arc) => {
                let u = arc.scaled(x);
                dot(&u, &u) == arc.r2 && arc.holds(&u)
            }
        }
    }

    /// The directions in which this edge leaves `x`, a point on it: one at
    /// an end, both ways inside.
    fn leaving(&self, x: &V) -> Vec<V> {
        let travel = match &self.kind {
            Kind::Segment => sub(&self.to, &self.from),
            Kind::Arc(arc) => {
                let u = arc.scaled(x);
                if arc.ccw() {
                    (u.1.neg(), u.0)
                } else {
                    (u.1, u.0.neg())
                }
            }
        };
        let back = (travel.0.neg(), travel.1.neg());
        match (*x == self.from, *x == self.to) {
            (true, _) => vec![travel],
            (_, true) => vec![back],
            _ => vec![travel, back],
        }
    }

    /// Whether `x` is within `t` of this edge. `None` when undecided.
    fn near_point(&self, x: &V, t: &Dyadic) -> Option<bool> {
        let t2 = t.square();
        let ends_near = |p: &V| {
            let w = sub(x, p);
            // Not the end itself: that is a meeting point.
            let d2 = dot(&w, &w);
            sign(&d2) != Sign::Zero && sign(&t2.sub(&d2)) == Sign::Positive
        };
        match &self.kind {
            Kind::Segment => {
                let d = sub(&self.to, &self.from);
                let w = sub(x, &self.from);
                let along = dot(&w, &d);
                let length2 = dot(&d, &d);
                if !at_least_zero(&along) {
                    return Some(ends_near(&self.from));
                }
                if sign(&along.sub(&length2)) != Sign::Negative {
                    return Some(ends_near(&self.to));
                }
                // Perpendicular distance: 0 < cross^2 < t^2 |D|^2; on the
                // segment itself it is a meeting point.
                let c = cross(&d, &w).square();
                Some(sign(&c) != Sign::Zero && sign(&t2.mul(&length2).sub(&c)) == Sign::Positive)
            }
            Kind::Arc(arc) => {
                let u = arc.scaled(x);
                let tk = t.mul(&arc.ka);
                if is_zero(&u) {
                    // At the centre: every arc point is `r` away.
                    return Some(sign(&tk.square().sub(&arc.r2)) == Sign::Positive);
                }
                if arc.holds(&u) {
                    // | |u| - rho | < tk
                    // On the arc itself it is a meeting point, not a near miss.
                    return abs_below(
                        &dy(0.0),
                        &[(dy(1.0), dot(&u, &u)), (dy(-1.0), arc.r2.clone())],
                        &(tk, None),
                    )
                    .map(|dist| dist == Dist::Near);
                }
                Some(ends_near(&self.from) || ends_near(&self.to))
            }
        }
    }
}

impl ArcShape {
    /// Whether the arc runs counter-clockwise from its first end.
    fn ccw(&self) -> bool {
        sign(&self.ks) == Sign::Positive
    }

    /// Whether the direction `v` from the centre meets the arc strictly
    /// between its ends.
    fn holds_open(&self, v: &V) -> bool {
        self.holds(v) && !same_direction(v, &self.start) && !same_direction(v, &self.end)
    }

    /// `ka Y - ks K`.
    fn scaled(&self, y: &V) -> V {
        sub(&scale(y, &self.ka), &scale(&self.k_centre, &self.ks))
    }

    /// Whether the direction `v` (not zero) from the centre meets the arc,
    /// ends included.
    fn holds(&self, v: &V) -> bool {
        let (a, e) = (&self.start, &self.end);
        if self.major {
            // Not strictly inside the complementary arc from `end` to `start`.
            !(sign(&cross(e, v)) == Sign::Positive && sign(&cross(v, a)) == Sign::Positive)
        } else {
            at_least_zero(&cross(a, v)) && at_least_zero(&cross(v, e))
        }
    }
}

/// Where a distance stands against the clearance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Dist {
    /// At least the clearance.
    Far,
    /// Strictly between zero and the clearance.
    Near,
    /// Exactly zero: a meeting point.
    Zero,
}

/// Where `| p + sum q_i sqrt(c_i) |` stands against `tau_c * sqrt(tau_e)`
/// (`tau_e` one when absent), over exact dyadics. `None` when the tower
/// cannot decide.
fn abs_below(
    p: &Dyadic,
    roots: &[(Dyadic, Dyadic)],
    tau: &(Dyadic, Option<Dyadic>),
) -> Option<Dist> {
    let mut tower = Tower::<Dyadic>::new();
    let mut sum = tower.value(p.clone());
    for (q, c) in roots {
        let radicand = tower.value(c.clone());
        let root = tower.sqrt(&radicand).ok()?;
        sum = tower.add(&sum, &tower.mul(&tower.value(q.clone()), &root));
    }
    let mut bound = tower.value(tau.0.clone());
    if let Some(e) = &tau.1 {
        let radicand = tower.value(e.clone());
        let root = tower.sqrt(&radicand).ok()?;
        bound = tower.mul(&bound, &root);
    }
    let below = tower.sign(&tower.sub(&sum, &bound))?;
    let above = tower.sign(&tower.add(&sum, &bound))?;
    if tower.sign(&sum)? == Sign::Zero {
        return Some(Dist::Zero);
    }
    let near = below == Sign::Negative && above == Sign::Positive;
    Some(if near { Dist::Near } else { Dist::Far })
}

/// Interior points of a segment and an arc along a common normal.
fn segment_arc(segment: &Shape, arc: &ArcShape, t: &Dyadic) -> Option<Dist> {
    let d = sub(&segment.to, &segment.from);
    let length2 = dot(&d, &d);
    if sign(&length2) == Sign::Zero {
        return Some(Dist::Far);
    }
    // The centre against the segment's line, scaled by `ka`.
    let c = sub(
        &scale(&arc.k_centre, &arc.ks),
        &scale(&segment.from, &arc.ka),
    );
    let along = dot(&c, &d);
    if !at_least_zero(&along) || !at_least_zero(&arc.ka.mul(&length2).sub(&along)) {
        return Some(Dist::Far);
    }
    let height = cross(&d, &c);
    let normal = (d.1.neg(), d.0.clone());
    let tk = t.mul(&arc.ka);
    let inside =
        sign(&along) == Sign::Positive && sign(&arc.ka.mul(&length2).sub(&along)) == Sign::Positive;
    let mut best = Dist::Far;
    let mut undecided = false;
    for sigma in [dy(1.0), dy(-1.0)] {
        let direction = scale(&normal, &sigma);
        if !arc.holds(&direction) {
            continue;
        }
        // A touch at an end of either edge is judged by its directions
        // (`Shape::cusp_with`); only one inside both is a tangency here.
        let interior = inside && arc.holds_open(&direction);
        // | height + sigma rho sqrt(L) | against tk sqrt(L), with rho^2 = r2;
        // zero is a tangent touch of the line and the circle.
        let roots = [(sigma, arc.r2.mul(&length2))];
        match abs_below(&height, &roots, &(tk.clone(), Some(length2.clone()))) {
            Some(Dist::Zero) if !interior => {}
            Some(dist) => best = best.max(dist),
            None => undecided = true,
        }
    }
    decided_unless(best, undecided)
}

/// `best`, or `None` (undecided) when nothing below or at zero was found and
/// some test could not decide: a decided hit stands whatever else is open.
fn decided_unless(best: Dist, undecided: bool) -> Option<Dist> {
    if best == Dist::Far && undecided {
        None
    } else {
        Some(best)
    }
}

/// Interior points of two arcs on the line through their centres.
fn arc_arc(a: &ArcShape, b: &ArcShape, t: &Dyadic) -> Option<Dist> {
    // `ka_a ka_b (c_b - c_a)`.
    let w = sub(
        &scale(&b.k_centre, &a.ka.mul(&b.ks)),
        &scale(&a.k_centre, &b.ka.mul(&a.ks)),
    );
    if is_zero(&w) {
        return Some(Dist::Far);
    }
    let w2 = dot(&w, &w);
    let tk = t.mul(&a.ka).mul(&b.ka);
    let mut best = Dist::Far;
    let mut undecided = false;
    for sa in [dy(1.0), dy(-1.0)] {
        let da = scale(&w, &sa);
        if !a.holds(&da) {
            continue;
        }
        for sb in [dy(1.0), dy(-1.0)] {
            let db = scale(&w, &sb);
            if !b.holds(&db) {
                continue;
            }
            let interior = a.holds_open(&da) && b.holds_open(&db);
            // | |w| + sb r_b - sa r_a | < t, times ka_a ka_b.
            let roots = [
                (dy(1.0), w2.clone()),
                (sb, a.ka.square().mul(&b.r2)),
                (sa.neg(), b.ka.square().mul(&a.r2)),
            ];
            match abs_below(&dy(0.0), &roots, &(tk.clone(), None)) {
                Some(Dist::Zero) if !interior => {}
                Some(dist) => best = best.max(dist),
                None => undecided = true,
            }
        }
    }
    decided_unless(best, undecided)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64) -> Point2 {
        Point2::new(x, y)
    }

    #[test]
    fn a_scaled_arc_end_lies_on_the_scaled_circle() {
        for bulge in [0.5, -0.5, 1.0, 2.0, -3.0, 0.1] {
            let shape = Shape::new(p(0.0, 0.0), p(2.0, 1.0), bulge);
            let Kind::Arc(arc) = &shape.kind else {
                panic!("an arc");
            };
            assert_eq!(dot(&arc.start, &arc.start), arc.r2, "bulge {bulge}");
            assert_eq!(dot(&arc.end, &arc.end), arc.r2, "bulge {bulge}");
        }
    }

    #[test]
    fn a_negative_radicand_is_undecided_not_a_guess() {
        let result = abs_below(&dy(0.0), &[(dy(1.0), dy(-1.0))], &(dy(1.0), None));
        assert_eq!(result, None);
    }

    fn square_and_hole(gap: f64) -> Vec<ArcRing> {
        let ring = |x0: f64, y0: f64, x1: f64, y1: f64| {
            ArcRing::from_points(&[p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)])
        };
        vec![ring(0.0, 0.0, 10.0, 10.0), ring(gap, 1.0, 5.0, 4.0)]
    }

    #[test]
    fn an_undecided_pair_is_reported_as_undecided_and_a_clear_one_is_not() {
        let rings = square_and_hole(0.0005);
        let rings: Vec<&ArcRing> = rings.iter().collect();
        let run = |verdict: Verdict| near_pairs_by(&rings, &[true, true], 0.001, |_, _| verdict);
        // Every pair whose boxes are within the clearance is a candidate.
        let undecided = run(Verdict::Undecided);
        assert!(!undecided.is_empty());
        assert!(undecided.iter().all(|pair| pair.undecided));
        let below = run(Verdict::Below);
        assert_eq!(below.len(), undecided.len());
        assert!(below.iter().all(|pair| !pair.undecided));
        assert!(run(Verdict::Clear).is_empty());
    }

    #[test]
    fn a_ring_that_is_not_usable_is_left_out() {
        let rings = square_and_hole(0.0005);
        let rings: Vec<&ArcRing> = rings.iter().collect();
        let found = near_pairs_by(&rings, &[true, false], 0.001, |_, _| Verdict::Below);
        assert!(!found.is_empty());
        assert!(found
            .iter()
            .all(|pair| pair.edges.iter().all(|edge| edge.0 == 0)));
    }

    #[test]
    fn parallel_segments_are_below_only_within_their_distance() {
        let a = Shape::new(p(0.0, 0.0), p(1.0, 0.0), 0.0);
        let b = Shape::new(p(0.0, 1.0), p(1.0, 1.0), 0.0);
        assert_eq!(a.closer_than(&b, 2.0), Verdict::Below);
        assert_eq!(a.closer_than(&b, 0.5), Verdict::Clear);
    }
}
