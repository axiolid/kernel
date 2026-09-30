//! Every intersection among many segments: a Bentley–Ottmann sweep (#146).
//!
//! [`segment_intersections`] reports each point where two or more input
//! segments meet, with every segment through it, and each collinear piece
//! that two or more segments share, in `O((n + k) log n)` time for `n`
//! segments and `k` reported incidences (segment–point pairs; `k` is twice
//! the number of crossings when segments are in general position). The
//! pairwise alternative, calling a segment/segment test on every pair, is
//! `O(n^2)`.
//!
//! # Exactness
//!
//! Every decision is an exact sign: the order of two event points, which
//! side of a segment a point lies on, whether two segments cross, and the
//! order of two directions. Each is asked of an outward-rounded interval
//! first and of exact dyadic arithmetic (`axiolid-exact`) only when the
//! interval cannot decide, so no tolerance enters and the answer does not
//! depend on drawing units. A crossing point is a rational number, kept
//! exactly as homogeneous dyadic coordinates ([`ExactPoint2`]) and rounded
//! once, correctly (nearest double, ties to even), when a caller asks for
//! [`ExactPoint2::rounded`]. A rounded value never feeds a decision.
//!
//! # Sweep order and degeneracies
//!
//! Points are swept in lexicographic order, `x` first and `y` breaking ties,
//! which is a sweep line turned by an infinitesimal angle. A vertical
//! segment is then an ordinary segment that the line crosses from its lower
//! end to its upper end, and needs no special case. Along the line,
//! segments through the current point are ordered by direction, and
//! segments with the same direction through the same point lie on one line
//! and form a bundle; a bundle of two or more is an overlap.
//!
//! - Shared endpoints, T-junctions and many segments through one point are
//!   found by locating the event point in the sweep status: all active
//!   segments through it form one contiguous run, whatever their number.
//! - Only proper crossings (interior to both segments) are scheduled as new
//!   events; every endpoint is an event from the start, so a crossing at an
//!   endpoint is already one.
//! - A zero-length segment is a point. It is reported where it lies on
//!   another segment or on another zero-length segment, with
//!   [`SegmentLocation::Degenerate`], and it takes no part in overlaps.
//! - Non-finite coordinates are refused with the segment's index.
//!
//! # What is reported
//!
//! A point is reported when at least two segments contain it and it is an
//! endpoint of one of them or a proper crossing of two. Equivalently: every
//! point on two or more segments, except those inside an overlap where all
//! segments through the point are collinear and pass through it, which say
//! nothing the overlap does not.
//! Overlaps are maximal collinear pieces covered by a constant set of at
//! least two segments; a new set starts only where a segment of the line
//! starts or ends, so both ends of an overlap are input endpoints and are
//! also reported as points.
//!
//! The model is CGAL's `Surface_sweep_2` (`compute_intersection_points`,
//! `compute_subcurves`), restricted to segments and kept exact.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::sync::{Arc, OnceLock};

use axiolid_core::Point2;
use axiolid_exact::{certify, Arith, Dyadic, Interval, SignExpr};
use axiolid_guarantees::Sign;

use crate::exact_arc::{orient_doubles, round_ratio};

/// Why the sweep refused its input.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SegmentSweepError {
    /// A coordinate of the segment at this index was NaN or infinite.
    NonFinite {
        /// Index of the segment in the input slice.
        segment: usize,
    },
}

impl core::fmt::Display for SegmentSweepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFinite { segment } => {
                write!(f, "segment {segment} has a non-finite coordinate")
            }
        }
    }
}

impl std::error::Error for SegmentSweepError {}

/// A planar point with exact rational coordinates.
///
/// Either an input endpoint, which is a pair of doubles, or a crossing of
/// two segments, `(x / w, y / w)` with dyadic `x`, `y` and `w > 0`.
/// Equality and order are exact; the order is lexicographic (`x`, then
/// `y`), which is the sweep order.
#[derive(Debug, Clone)]
pub struct ExactPoint2 {
    repr: Repr,
}

#[derive(Debug, Clone)]
enum Repr {
    Input(Point2),
    Rational(Arc<Rational>),
}

#[derive(Debug)]
struct Rational {
    x: Dyadic,
    y: Dyadic,
    w: Dyadic,
    /// Sound enclosures of `x / w` and `y / w`: the fast tier.
    bx: Interval,
    by: Interval,
    rounded: OnceLock<Point2>,
}

impl ExactPoint2 {
    /// The point `(x / w, y / w)`, or `None` when `w` is zero.
    #[must_use]
    pub fn from_homogeneous(x: Dyadic, y: Dyadic, w: Dyadic) -> Option<Self> {
        let flip = match w.sign()? {
            Sign::Positive => false,
            Sign::Negative => true,
            _ => return None,
        };
        let (x, y, w) = if flip {
            (x.neg(), y.neg(), w.neg())
        } else {
            (x, y, w)
        };
        let we = w.enclosure();
        Some(Self {
            repr: Repr::Rational(Arc::new(Rational {
                bx: x.enclosure().quotient(we),
                by: y.enclosure().quotient(we),
                x,
                y,
                w,
                rounded: OnceLock::new(),
            })),
        })
    }

    fn input(point: Point2) -> Self {
        Self {
            repr: Repr::Input(point),
        }
    }

    /// The input endpoint this point is, when it is one.
    ///
    /// Every point the sweep reports that is not a proper crossing is an
    /// input endpoint, and then these are exactly its doubles.
    #[must_use]
    pub fn as_input(&self) -> Option<Point2> {
        match &self.repr {
            Repr::Input(p) => Some(*p),
            Repr::Rational(_) => None,
        }
    }

    /// Exact homogeneous coordinates `(x, y, w)` with `w > 0`.
    #[must_use]
    pub fn homogeneous(&self) -> (Dyadic, Dyadic, Dyadic) {
        match &self.repr {
            Repr::Input(p) => (dy(p.x), dy(p.y), dy(1.0)),
            Repr::Rational(r) => (r.x.clone(), r.y.clone(), r.w.clone()),
        }
    }

    /// Each coordinate correctly rounded to the nearest double, ties to
    /// even. An input endpoint comes back bit for bit.
    ///
    /// For output only: rounding moves the point, so a rounded crossing
    /// may lie off both segments it is on.
    #[must_use]
    pub fn rounded(&self) -> Point2 {
        match &self.repr {
            Repr::Input(p) => *p,
            Repr::Rational(r) => *r.rounded.get_or_init(|| {
                let guess = |n: &Dyadic| {
                    let (mn, en) = n.approx_parts();
                    let (mw, ew) = r.w.approx_parts();
                    scale2(mn / mw, en - ew)
                };
                Point2::new(
                    round_ratio(&r.x, &r.w, guess(&r.x)),
                    round_ratio(&r.y, &r.w, guess(&r.y)),
                )
            }),
        }
    }

    fn bound(&self, axis: Axis) -> Interval {
        match (&self.repr, axis) {
            (Repr::Input(p), Axis::X) => Interval::point(p.x),
            (Repr::Input(p), Axis::Y) => Interval::point(p.y),
            (Repr::Rational(r), Axis::X) => r.bx,
            (Repr::Rational(r), Axis::Y) => r.by,
        }
    }

    fn exact(&self, axis: Axis) -> (Dyadic, Dyadic) {
        match (&self.repr, axis) {
            (Repr::Input(p), Axis::X) => (dy(p.x), dy(1.0)),
            (Repr::Input(p), Axis::Y) => (dy(p.y), dy(1.0)),
            (Repr::Rational(r), Axis::X) => (r.x.clone(), r.w.clone()),
            (Repr::Rational(r), Axis::Y) => (r.y.clone(), r.w.clone()),
        }
    }
}

impl PartialEq for ExactPoint2 {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for ExactPoint2 {}

impl PartialOrd for ExactPoint2 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ExactPoint2 {
    fn cmp(&self, other: &Self) -> Ordering {
        if let (Repr::Input(a), Repr::Input(b)) = (&self.repr, &other.repr) {
            // Finite doubles: `total_cmp` would order -0 before +0, which
            // are the same point, so compare as reals.
            return cmp_f64(a.x, b.x).then_with(|| cmp_f64(a.y, b.y));
        }
        cmp_axis(self, other, Axis::X).then_with(|| cmp_axis(self, other, Axis::Y))
    }
}

/// Where on a segment an incidence lies, in the segment's own direction
/// (`Start` is the segment's first point as given).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum SegmentLocation {
    /// The segment's first point.
    Start,
    /// The segment's second point.
    End,
    /// Strictly between its endpoints.
    Interior,
    /// The segment has zero length and this is its only point.
    Degenerate,
}

/// One segment through a reported point.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Incidence {
    /// Index of the segment in the input slice.
    pub segment: usize,
    /// Where on that segment the point lies.
    pub location: SegmentLocation,
}

/// A point where two or more segments meet.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct IntersectionPoint {
    /// The exact point.
    pub point: ExactPoint2,
    /// [`ExactPoint2::rounded`]: correctly rounded, for output.
    pub rounded: Point2,
    /// Every segment containing the point, by increasing index.
    pub incidences: Vec<Incidence>,
}

/// A collinear piece shared by two or more segments.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct SegmentOverlap {
    /// The lexicographically smaller end. Always an input endpoint, so it
    /// is exact.
    pub from: Point2,
    /// The lexicographically larger end, also an input endpoint.
    pub to: Point2,
    /// The segments covering the whole piece, by increasing index. The set
    /// changes at both ends (or the piece ends there).
    pub segments: Vec<usize>,
}

/// What the sweep did, for measurement.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SweepEvidence {
    /// Input segments.
    pub segments: usize,
    /// Distinct event points swept: endpoints and proper crossings.
    pub events: usize,
    /// Proper crossings found between neighbouring segments.
    pub crossings: usize,
}

/// The result of [`segment_intersections`].
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct SegmentIntersections {
    /// Every reported point, in lexicographic order (`x`, then `y`).
    pub points: Vec<IntersectionPoint>,
    /// Every overlap, ordered by its lexicographically smaller end, then
    /// by where the sweep first met it.
    pub overlaps: Vec<SegmentOverlap>,
    /// Counters.
    pub evidence: SweepEvidence,
}

/// Every intersection among `segments`, exactly, by a Bentley–Ottmann
/// sweep in `O((n + k) log n)` time for `k` reported incidences.
///
/// A point is reported when at least two segments contain it and it is an
/// endpoint of one of them or a proper crossing of two, with every segment
/// through it. Each maximal collinear piece covered by a constant set of
/// two or more segments is reported as a [`SegmentOverlap`]; its ends are
/// input endpoints and are reported as points too. A zero-length segment is
/// a point ([`SegmentLocation::Degenerate`]) and takes no part in overlaps.
///
/// Every decision is an exact sign, so no tolerance enters. Crossings are
/// [`ExactPoint2`] rationals; [`IntersectionPoint::rounded`] is their
/// correctly rounded value, and an input endpoint comes back bit for bit.
/// Output order is deterministic: points in lexicographic order (`x`, then
/// `y`), incidences by segment index.
///
/// # Errors
///
/// [`SegmentSweepError::NonFinite`] when a coordinate is NaN or infinite.
pub fn segment_intersections(
    segments: &[[Point2; 2]],
) -> Result<SegmentIntersections, SegmentSweepError> {
    let segs = prepare(segments)?;
    Ok(Sweep::new(segs).run())
}

// ------------------------------------------------------------ arithmetic

#[derive(Debug, Clone, Copy)]
enum Axis {
    X,
    Y,
}

fn dy(value: f64) -> Dyadic {
    Dyadic::from_f64(value)
}

fn to_i8(sign: Sign) -> i8 {
    match sign {
        Sign::Positive => 1,
        Sign::Negative => -1,
        _ => 0,
    }
}

fn exact_sign(value: &Dyadic) -> i8 {
    to_i8(value.sign().expect("dyadic signs are always decided"))
}

fn cmp_f64(a: f64, b: f64) -> Ordering {
    a.partial_cmp(&b)
        .expect("coordinates were validated finite")
}

/// `x * 2^n` without intermediate overflow; for a rounding guess only.
fn scale2(mut x: f64, mut n: i64) -> f64 {
    while n > 1000 && x.is_finite() && x != 0.0 {
        x *= 2f64.powi(1000);
        n -= 1000;
    }
    while n < -1000 && x != 0.0 {
        x *= 2f64.powi(-1000);
        n += 1000;
    }
    x * 2f64.powi(n as i32)
}

/// Exact order of one coordinate: disjoint enclosures first, then
/// `a / wa - b / wb` with the denominators cleared.
fn cmp_axis(p: &ExactPoint2, q: &ExactPoint2, axis: Axis) -> Ordering {
    let (pb, qb) = (p.bound(axis), q.bound(axis));
    if pb.hi() < qb.lo() {
        return Ordering::Less;
    }
    if qb.hi() < pb.lo() {
        return Ordering::Greater;
    }
    let (pn, pw) = p.exact(axis);
    let (qn, qw) = q.exact(axis);
    exact_sign(&pn.mul(&qw).sub(&qn.mul(&pw))).cmp(&0)
}

/// Which side of the line `a -> b` the point lies on: positive for left.
fn orient_point(a: Point2, b: Point2, p: &ExactPoint2) -> i8 {
    match &p.repr {
        Repr::Input(c) => {
            if *c == a || *c == b {
                return 0;
            }
            to_i8(orient_doubles(a, b, *c))
        }
        Repr::Rational(r) => {
            let (ax, ay) = (Interval::point(a.x), Interval::point(a.y));
            let bax = Interval::point(b.x).sub(&ax);
            let bay = Interval::point(b.y).sub(&ay);
            let value = bax.mul(&r.by.sub(&ay)).sub(&bay.mul(&r.bx.sub(&ax)));
            if let Some(sign) = value.sign() {
                return to_i8(sign);
            }
            // (b - a) x (p - a) * w, and w > 0.
            let (ax, ay) = (dy(a.x), dy(a.y));
            let bax = dy(b.x).sub(&ax);
            let bay = dy(b.y).sub(&ay);
            let value = bax
                .mul(&r.y.sub(&ay.mul(&r.w)))
                .sub(&bay.mul(&r.x.sub(&ax.mul(&r.w))));
            exact_sign(&value)
        }
    }
}

/// `cross(u1 - u0, v1 - v0)`: positive when `v` turns left of `u`.
struct CrossDirections([Point2; 4]);

impl SignExpr for CrossDirections {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let [u0, u1, v0, v1] = self.0;
        let ux = T::from_f64(u1.x).sub(&T::from_f64(u0.x));
        let uy = T::from_f64(u1.y).sub(&T::from_f64(u0.y));
        let vx = T::from_f64(v1.x).sub(&T::from_f64(v0.x));
        let vy = T::from_f64(v1.y).sub(&T::from_f64(v0.y));
        ux.mul(&vy).sub(&uy.mul(&vx)).sign()
    }
}

// ------------------------------------------------------------- segments

/// An input segment with its ends in sweep order.
#[derive(Debug, Clone, Copy)]
struct Seg {
    lo: Point2,
    hi: Point2,
    /// `lo` is the input's second point.
    reversed: bool,
    degenerate: bool,
}

impl Seg {
    fn lo_location(&self) -> SegmentLocation {
        if self.degenerate {
            SegmentLocation::Degenerate
        } else if self.reversed {
            SegmentLocation::End
        } else {
            SegmentLocation::Start
        }
    }

    fn hi_location(&self) -> SegmentLocation {
        if self.reversed {
            SegmentLocation::Start
        } else {
            SegmentLocation::End
        }
    }
}

fn prepare(segments: &[[Point2; 2]]) -> Result<Vec<Seg>, SegmentSweepError> {
    segments
        .iter()
        .enumerate()
        .map(|(index, [a, b])| {
            if ![a.x, a.y, b.x, b.y].iter().all(|v| v.is_finite()) {
                return Err(SegmentSweepError::NonFinite { segment: index });
            }
            let reversed = cmp_f64(a.x, b.x).then_with(|| cmp_f64(a.y, b.y)) == Ordering::Greater;
            let (lo, hi) = if reversed { (*b, *a) } else { (*a, *b) };
            Ok(Seg {
                lo,
                hi,
                reversed,
                degenerate: lo == hi,
            })
        })
        .collect()
}

/// The order of two directions through one point, along the sweep line:
/// positive when `t` lies above `s` just after the point.
fn cross_sign(s: &Seg, t: &Seg) -> i8 {
    let expr = CrossDirections([s.lo, s.hi, t.lo, t.hi]);
    to_i8(certify(&expr).expect("inputs were validated finite"))
}

fn boxes_meet(s: &Seg, t: &Seg) -> bool {
    let (sy0, sy1) = (s.lo.y.min(s.hi.y), s.lo.y.max(s.hi.y));
    let (ty0, ty1) = (t.lo.y.min(t.hi.y), t.lo.y.max(t.hi.y));
    s.lo.x <= t.hi.x && t.lo.x <= s.hi.x && sy0 <= ty1 && ty0 <= sy1
}

/// Whether the segments cross at one point interior to both.
fn proper_crossing(s: &Seg, t: &Seg) -> bool {
    if !boxes_meet(s, t) {
        return false;
    }
    let o = |a: Point2, b: Point2, c: Point2| to_i8(orient_doubles(a, b, c));
    let (o1, o2) = (o(s.lo, s.hi, t.lo), o(s.lo, s.hi, t.hi));
    if o1 == 0 || o2 == 0 || o1 == o2 {
        return false;
    }
    let (o3, o4) = (o(t.lo, t.hi, s.lo), o(t.lo, t.hi, s.hi));
    o3 != 0 && o4 != 0 && o3 != o4
}

/// The crossing of the lines through `s` and `t`, which are not parallel:
/// `s.lo + (num / den) * (s.hi - s.lo)` with the denominator cleared.
fn crossing_point(s: &Seg, t: &Seg) -> ExactPoint2 {
    let (ax, ay) = (dy(s.lo.x), dy(s.lo.y));
    let (ux, uy) = (dy(s.hi.x).sub(&ax), dy(s.hi.y).sub(&ay));
    let (cx, cy) = (dy(t.lo.x), dy(t.lo.y));
    let (vx, vy) = (dy(t.hi.x).sub(&cx), dy(t.hi.y).sub(&cy));
    let den = ux.mul(&vy).sub(&uy.mul(&vx));
    let (acx, acy) = (cx.sub(&ax), cy.sub(&ay));
    let num = acx.mul(&vy).sub(&acy.mul(&vx));
    let x = ax.mul(&den).add(&num.mul(&ux));
    let y = ay.mul(&den).add(&num.mul(&uy));
    ExactPoint2::from_homogeneous(x, y, den).expect("proper crossings are not parallel")
}

// ---------------------------------------------------------------- status

const NIL: u32 = u32::MAX;

#[derive(Debug, Clone, Copy)]
struct Node {
    left: u32,
    right: u32,
    parent: u32,
    priority: u64,
}

/// The sweep status: a treap over segment indices, in order along the
/// sweep line from below. A node's index is its segment's, so a segment is
/// removed by handle and never needs comparing to be found.
///
/// Priorities come from a fixed-seed generator: expected depth
/// `O(log n)`, and a reproducible shape.
struct Status {
    nodes: Vec<Node>,
    root: u32,
    state: u64,
}

impl Status {
    fn new(count: usize) -> Self {
        let empty = Node {
            left: NIL,
            right: NIL,
            parent: NIL,
            priority: 0,
        };
        Self {
            nodes: vec![empty; count],
            root: NIL,
            state: 0x9E37_79B9_7F4A_7C15,
        }
    }

    fn random(&mut self) -> u64 {
        // xorshift64*
        self.state ^= self.state >> 12;
        self.state ^= self.state << 25;
        self.state ^= self.state >> 27;
        self.state.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn node(&self, id: u32) -> &Node {
        &self.nodes[id as usize]
    }

    fn node_mut(&mut self, id: u32) -> &mut Node {
        &mut self.nodes[id as usize]
    }

    /// The first node, in order, for which `pred` holds; `pred` must be
    /// false on a prefix and true on the rest.
    fn first_where(&self, mut pred: impl FnMut(u32) -> bool) -> u32 {
        let (mut current, mut found) = (self.root, NIL);
        while current != NIL {
            if pred(current) {
                found = current;
                current = self.node(current).left;
            } else {
                current = self.node(current).right;
            }
        }
        found
    }

    fn min_from(&self, mut id: u32) -> u32 {
        while self.node(id).left != NIL {
            id = self.node(id).left;
        }
        id
    }

    fn max_from(&self, mut id: u32) -> u32 {
        while self.node(id).right != NIL {
            id = self.node(id).right;
        }
        id
    }

    fn last(&self) -> u32 {
        if self.root == NIL {
            NIL
        } else {
            self.max_from(self.root)
        }
    }

    fn next(&self, mut id: u32) -> u32 {
        if self.node(id).right != NIL {
            return self.min_from(self.node(id).right);
        }
        loop {
            let parent = self.node(id).parent;
            if parent == NIL || self.node(parent).left == id {
                return parent;
            }
            id = parent;
        }
    }

    fn prev(&self, mut id: u32) -> u32 {
        if self.node(id).left != NIL {
            return self.max_from(self.node(id).left);
        }
        loop {
            let parent = self.node(id).parent;
            if parent == NIL || self.node(parent).right == id {
                return parent;
            }
            id = parent;
        }
    }

    /// Replace `old` by `new` as a child of `parent` (or as the root).
    fn relink(&mut self, parent: u32, old: u32, new: u32) {
        if parent == NIL {
            self.root = new;
        } else if self.node(parent).left == old {
            self.node_mut(parent).left = new;
        } else {
            self.node_mut(parent).right = new;
        }
        if new != NIL {
            self.node_mut(new).parent = parent;
        }
    }

    /// Rotate `id` above its parent, keeping the in-order sequence.
    fn rotate_up(&mut self, id: u32) {
        let parent = self.node(id).parent;
        let grand = self.node(parent).parent;
        if self.node(parent).left == id {
            let middle = self.node(id).right;
            self.node_mut(parent).left = middle;
            if middle != NIL {
                self.node_mut(middle).parent = parent;
            }
            self.node_mut(id).right = parent;
        } else {
            let middle = self.node(id).left;
            self.node_mut(parent).right = middle;
            if middle != NIL {
                self.node_mut(middle).parent = parent;
            }
            self.node_mut(id).left = parent;
        }
        self.node_mut(parent).parent = id;
        self.relink(grand, parent, id);
    }

    /// Insert `id` immediately before `before`, or last when `before` is
    /// [`NIL`].
    fn insert_before(&mut self, id: u32, before: u32) {
        let priority = self.random();
        *self.node_mut(id) = Node {
            left: NIL,
            right: NIL,
            parent: NIL,
            priority,
        };
        if self.root == NIL {
            self.root = id;
            return;
        }
        if before == NIL {
            let last = self.max_from(self.root);
            self.node_mut(last).right = id;
            self.node_mut(id).parent = last;
        } else if self.node(before).left == NIL {
            self.node_mut(before).left = id;
            self.node_mut(id).parent = before;
        } else {
            let last = self.max_from(self.node(before).left);
            self.node_mut(last).right = id;
            self.node_mut(id).parent = last;
        }
        while self.node(id).parent != NIL
            && self.node(self.node(id).parent).priority < self.node(id).priority
        {
            self.rotate_up(id);
        }
    }

    fn remove(&mut self, id: u32) {
        loop {
            let Node { left, right, .. } = *self.node(id);
            if left == NIL || right == NIL {
                let child = if left == NIL { right } else { left };
                let parent = self.node(id).parent;
                self.relink(parent, id, child);
                return;
            }
            let up = if self.node(left).priority > self.node(right).priority {
                left
            } else {
                right
            };
            self.rotate_up(up);
        }
    }
}

// ----------------------------------------------------------------- sweep

struct Piece {
    from: ExactPoint2,
    to: Option<ExactPoint2>,
    segments: Vec<usize>,
}

struct Sweep {
    segs: Vec<Seg>,
    /// Event points, each with the segments that start there.
    queue: BTreeMap<ExactPoint2, Vec<u32>>,
    status: Status,
    /// The open overlap each segment is part of.
    open: Vec<Option<usize>>,
    pieces: Vec<Piece>,
    /// Pieces that ended at the current event, which may continue.
    closed_here: Vec<usize>,
    points: Vec<IntersectionPoint>,
    evidence: SweepEvidence,
}

impl Sweep {
    // The only interior mutability in a point is its cached rounding,
    // which order and equality never read.
    #[allow(clippy::mutable_key_type)]
    fn new(segs: Vec<Seg>) -> Self {
        let mut queue: BTreeMap<ExactPoint2, Vec<u32>> = BTreeMap::new();
        for (index, seg) in segs.iter().enumerate() {
            queue
                .entry(ExactPoint2::input(seg.lo))
                .or_default()
                .push(index as u32);
            if !seg.degenerate {
                queue.entry(ExactPoint2::input(seg.hi)).or_default();
            }
        }
        let count = segs.len();
        Self {
            status: Status::new(count),
            open: vec![None; count],
            pieces: Vec::new(),
            closed_here: Vec::new(),
            points: Vec::new(),
            evidence: SweepEvidence {
                segments: count,
                ..SweepEvidence::default()
            },
            queue,
            segs,
        }
    }

    fn run(mut self) -> SegmentIntersections {
        while let Some((p, starts)) = self.queue.pop_first() {
            self.evidence.events += 1;
            self.event(&p, &starts);
        }
        let mut overlaps: Vec<SegmentOverlap> = self
            .pieces
            .into_iter()
            .map(|piece| {
                let to = piece.to.expect("every overlap ends at an event");
                SegmentOverlap {
                    from: endpoint(&piece.from),
                    to: endpoint(&to),
                    segments: piece.segments,
                }
            })
            .collect();
        overlaps
            .sort_by(|a, b| cmp_f64(a.from.x, b.from.x).then_with(|| cmp_f64(a.from.y, b.from.y)));
        SegmentIntersections {
            points: self.points,
            overlaps,
            evidence: self.evidence,
        }
    }

    fn event(&mut self, p: &ExactPoint2, starts: &[u32]) {
        // The active segments through `p` are one contiguous run: the
        // status is ordered along the sweep line just before `p`.
        let segs = &self.segs;
        let side = |id: u32| {
            let seg = &segs[id as usize];
            orient_point(seg.lo, seg.hi, p)
        };
        let first = self.status.first_where(|id| side(id) <= 0);
        let mut block = Vec::new();
        let mut above = first;
        while above != NIL && side(above) == 0 {
            block.push(above);
            above = self.status.next(above);
        }
        let below = if first == NIL {
            self.status.last()
        } else {
            self.status.prev(first)
        };

        // A crossing event is never an input endpoint (endpoints are
        // queued from the start and deduplicate it), so only an input
        // event can end a segment.
        let at = p.as_input();
        let ended: Vec<bool> = block
            .iter()
            .map(|&id| at.is_some_and(|q| segs[id as usize].hi == q))
            .collect();

        self.report(p, starts, &block, &ended);
        self.close_overlaps(p, &block);

        for &id in &block {
            self.status.remove(id);
        }
        let mut through: Vec<u32> = starts
            .iter()
            .copied()
            .filter(|&id| !self.segs[id as usize].degenerate)
            .chain(
                block
                    .iter()
                    .zip(&ended)
                    .filter(|(_, &ended)| !ended)
                    .map(|(&id, _)| id),
            )
            .collect();
        let segs = &self.segs;
        through.sort_by(
            |&s, &t| match cross_sign(&segs[s as usize], &segs[t as usize]) {
                1 => Ordering::Less,
                -1 => Ordering::Greater,
                _ => s.cmp(&t),
            },
        );
        for &id in &through {
            self.status.insert_before(id, above);
        }
        self.open_overlaps(p, &through);

        match (through.first(), through.last()) {
            (Some(&low), Some(&high)) => {
                self.check(below, low, p);
                self.check(high, above, p);
            }
            _ => self.check(below, above, p),
        }
    }

    fn report(&mut self, p: &ExactPoint2, starts: &[u32], block: &[u32], ended: &[bool]) {
        if starts.len() + block.len() < 2 {
            return;
        }
        let mut incidences = Vec::with_capacity(starts.len() + block.len());
        for &id in starts {
            incidences.push(Incidence {
                segment: id as usize,
                location: self.segs[id as usize].lo_location(),
            });
        }
        for (&id, &ends) in block.iter().zip(ended) {
            let location = if ends {
                self.segs[id as usize].hi_location()
            } else {
                SegmentLocation::Interior
            };
            incidences.push(Incidence {
                segment: id as usize,
                location,
            });
        }
        // Every event is an endpoint (a segment starts or ends here) or a
        // proper crossing of two segments that are not collinear, so the
        // interior of an overlap is never an event and needs no filter.
        incidences.sort_unstable();
        self.points.push(IntersectionPoint {
            rounded: p.rounded(),
            point: p.clone(),
            incidences,
        });
    }

    fn close_overlaps(&mut self, p: &ExactPoint2, block: &[u32]) {
        self.closed_here.clear();
        for &id in block {
            if let Some(index) = self.open[id as usize].take() {
                let piece = &mut self.pieces[index];
                if piece.to.is_none() {
                    piece.to = Some(p.clone());
                    self.closed_here.push(index);
                }
            }
        }
    }

    /// Runs of equal direction in `through` share a line from `p` on; each
    /// run of two or more is an overlap. A run with the same segments as a
    /// piece that just ended here continues that piece.
    fn open_overlaps(&mut self, p: &ExactPoint2, through: &[u32]) {
        let mut start = 0;
        while start < through.len() {
            let head = &self.segs[through[start] as usize];
            let mut end = start + 1;
            while end < through.len() && cross_sign(head, &self.segs[through[end] as usize]) == 0 {
                end += 1;
            }
            if end - start >= 2 {
                let mut members: Vec<usize> =
                    through[start..end].iter().map(|&id| id as usize).collect();
                members.sort_unstable();
                let resumed = self
                    .closed_here
                    .iter()
                    .copied()
                    .find(|&index| self.pieces[index].segments == members);
                let index = if let Some(index) = resumed {
                    self.pieces[index].to = None;
                    index
                } else {
                    self.pieces.push(Piece {
                        from: p.clone(),
                        to: None,
                        segments: members.clone(),
                    });
                    self.pieces.len() - 1
                };
                for &member in &members {
                    self.open[member] = Some(index);
                }
            }
            start = end;
        }
    }

    /// Schedule the crossing of two neighbours if it lies ahead.
    fn check(&mut self, low: u32, high: u32, p: &ExactPoint2) {
        if low == NIL || high == NIL {
            return;
        }
        let (s, t) = (&self.segs[low as usize], &self.segs[high as usize]);
        if !proper_crossing(s, t) {
            return;
        }
        let q = crossing_point(s, t);
        if q > *p {
            if !self.queue.contains_key(&q) {
                self.evidence.crossings += 1;
            }
            self.queue.entry(q).or_default();
        }
    }
}

/// An overlap end, which is always an input endpoint (see the module
/// documentation).
fn endpoint(point: &ExactPoint2) -> Point2 {
    debug_assert!(point.as_input().is_some(), "overlap ends are input points");
    point.as_input().unwrap_or_else(|| point.rounded())
}
