//! Exact points of the plane in homogeneous coordinates.
//!
//! A point is `(x, y, w)` with `w > 0` and stands for `(x / w, y / w)`. The
//! coordinates are [`Dyadic`] numbers, so a point can be the crossing of two
//! segments or of an edge with a plane without being a binary64 position.
//!
//! The predicates are polynomials in these values. `orient` is the 3 x 3
//! determinant and `compare` is `c_a w_b - c_b w_a` for a coordinate `c`. When
//! every point has `w = 1` and binary64 coordinates (a "flat" point) they are
//! the exact `orient2d` and a plain comparison. Otherwise each is evaluated
//! twice at most: in `Interval`s of the stored enclosures, which decides
//! nearly every call, and, when that straddles zero, in `Dyadic`s of the
//! stored exact values. `Dyadic` has no range limit, so products of tiny
//! coordinates are no reason to refuse.
//!
//! The boxes in [`ExactPoints::bounds`] are sound enclosures of the points.
//! They only prune candidates and never decide.

use std::cell::Cell;
use std::cmp::Ordering;

use axiolid_contracts::{GeomError, GeomResult, Sign};
use axiolid_core::Point2;
use axiolid_exact::{Arith, Dyadic, Interval};
use axiolid_predicates::orient2d;

use crate::ratio::round_ratio;

/// Homogeneous coordinates `(x, y, w)`.
type Homogeneous<T> = [T; 3];

/// A point of an [`ExactPoints`] store.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PointId(pub usize);

#[derive(Debug, Clone)]
struct ExactPoint {
    exact: Homogeneous<Dyadic>,
    enclosure: Homogeneous<Interval>,
    /// `[x_lo, y_lo, x_hi, y_hi]` of the point.
    bounds: [f64; 4],
    /// The coordinates when `w = 1` and they are binary64.
    flat: Option<Point2>,
}

/// Where a point lies relative to a ring.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Location {
    /// In the interior.
    Inside,
    /// Outside the closed ring.
    Outside,
    /// On the ring.
    Boundary,
}

/// The extent of a ring, for pruning candidate pairs. Never decides.
#[derive(Debug, Clone, Copy)]
pub struct RingBox {
    /// `[min_x, min_y, max_x, max_y]` of the exact points, with every vertex
    /// box at its outer limit: contains the ring.
    bounds: [f64; 4],
    /// The same with every vertex box at its inner limit: the extremes of the
    /// ring reach at least this far.
    core: [f64; 4],
}

impl RingBox {
    /// Whether this ring can lie in `other` by their boxes: its smallest
    /// coordinate cannot be below the box of `other`, nor its largest above.
    pub fn possibly_within(&self, other: &Self) -> bool {
        self.core[0] >= other.bounds[0]
            && self.core[1] >= other.bounds[1]
            && self.core[2] <= other.bounds[2]
            && self.core[3] <= other.bounds[3]
    }
}

/// A store of exact points and the exact predicates on them.
///
/// Orientation is counter-clockwise in `(x, y)` times the frame sign given to
/// [`ExactPoints::with_frame_sign`]; a caller whose `(x, y)` are the
/// coordinates of a mirrored frame passes `Sign::Negative` and reads
/// "counter-clockwise" in its own frame.
#[derive(Debug, Clone)]
pub struct ExactPoints {
    frame_sign: Sign,
    points: Vec<ExactPoint>,
    escalations: Cell<usize>,
}

impl Default for ExactPoints {
    fn default() -> Self {
        Self::new()
    }
}

fn not_finite() -> GeomError {
    GeomError::Degenerate("an exact point is not finite".into())
}

pub(crate) fn undecided() -> GeomError {
    GeomError::Degenerate("an exact sign could not be decided".into())
}

pub(crate) fn inconsistent() -> GeomError {
    GeomError::Degenerate("the exact arrangement state is inconsistent".into())
}

fn exact(value: f64) -> GeomResult<Dyadic> {
    Dyadic::try_from_f64(value).ok_or_else(not_finite)
}

/// `p_i q_j - q_i p_j`.
fn minor<T: Arith>(p: &Homogeneous<T>, q: &Homogeneous<T>, i: usize, j: usize) -> T {
    p[i].mul(&q[j]).sub(&q[i].mul(&p[j]))
}

/// `p x q`: the line through two points, or the point where two lines meet.
fn cross<T: Arith>(p: &Homogeneous<T>, q: &Homogeneous<T>) -> Homogeneous<T> {
    [minor(p, q, 1, 2), minor(p, q, 2, 0), minor(p, q, 0, 1)]
}

/// `det[a; b; c]`: positive when the points turn counter-clockwise (their
/// weights are positive).
fn determinant<T: Arith>(a: &Homogeneous<T>, b: &Homogeneous<T>, c: &Homogeneous<T>) -> T {
    a[0].mul(&minor(b, c, 1, 2))
        .sub(&a[1].mul(&minor(b, c, 0, 2)))
        .add(&a[2].mul(&minor(b, c, 0, 1)))
}

/// `c_a w_b - c_b w_a` for coordinate `axis` (0 = x, 1 = y).
fn difference<T: Arith>(a: &Homogeneous<T>, b: &Homogeneous<T>, axis: usize) -> T {
    a[axis].mul(&b[2]).sub(&b[axis].mul(&a[2]))
}

pub(crate) fn ordering(sign: Sign) -> Ordering {
    match sign {
        Sign::Positive => Ordering::Greater,
        Sign::Negative => Ordering::Less,
        _ => Ordering::Equal,
    }
}

fn coordinate(point: Point2, axis: usize) -> f64 {
    if axis == 0 {
        point.x
    } else {
        point.y
    }
}

/// Product of two signs.
fn sign_times(left: Sign, right: Sign) -> Sign {
    if left == Sign::Zero || right == Sign::Zero {
        Sign::Zero
    } else if left == right {
        Sign::Positive
    } else {
        Sign::Negative
    }
}

/// `numerator / weight` as a binary64 number, when that is exact. `hint` is a
/// candidate, such as an end of an edge that does not move in this coordinate.
fn exact_quotient(numerator: &Dyadic, weight: &Dyadic, hint: Option<f64>) -> Option<f64> {
    let candidate = match hint {
        Some(hint) => hint,
        // Both convert exactly up to 53 bits, so the division is then
        // correctly rounded and equals the quotient whenever that is
        // representable.
        None if numerator.bits() <= 53 && weight.bits() <= 53 => {
            numerator.to_f64() / weight.to_f64()
        }
        None => return None,
    };
    let exact_candidate = Dyadic::try_from_f64(candidate)?;
    (exact_candidate.mul(weight) == *numerator).then_some(candidate)
}

impl ExactPoints {
    /// An empty store whose orientation is counter-clockwise in `(x, y)`.
    pub fn new() -> Self {
        Self::with_frame_sign(Sign::Positive)
    }

    /// An empty store whose orientation is `frame_sign` times the sense of
    /// `(x, y)`. `frame_sign` must not be `Sign::Zero`.
    pub fn with_frame_sign(frame_sign: Sign) -> Self {
        Self {
            frame_sign,
            points: Vec::new(),
            escalations: Cell::new(0),
        }
    }

    /// The frame sign the store was made with.
    pub fn frame_sign(&self) -> Sign {
        self.frame_sign
    }

    /// The number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether the store has no point.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// A binary64 position.
    pub fn push_flat(&mut self, flat: Point2) -> GeomResult<PointId> {
        let point = [exact(flat.x)?, exact(flat.y)?, exact(1.0)?];
        let enclosure = [
            Interval::point(flat.x),
            Interval::point(flat.y),
            Interval::point(1.0),
        ];
        self.points.push(ExactPoint {
            exact: point,
            enclosure,
            bounds: [flat.x, flat.y, flat.x, flat.y],
            flat: Some(flat),
        });
        Ok(PointId(self.points.len() - 1))
    }

    /// The point `(x / w, y / w)` of homogeneous coordinates with `w != 0`.
    /// A negative `w` is normalised. A point that is a binary64 position is
    /// stored as one, so it takes the fast path of the predicates.
    pub fn push(&mut self, homogeneous: [Dyadic; 3]) -> GeomResult<PointId> {
        self.push_hinted(homogeneous, [None, None])
    }

    /// Like [`ExactPoints::push`], with a candidate value for each
    /// coordinate to try first, such as the coordinate of an edge that does
    /// not move in it.
    pub fn push_hinted(
        &mut self,
        homogeneous: [Dyadic; 3],
        hints: [Option<f64>; 2],
    ) -> GeomResult<PointId> {
        let mut point = homogeneous;
        match point[2].sign() {
            Some(Sign::Positive) => {}
            Some(Sign::Negative) => point = point.map(|value| value.neg()),
            _ => return Err(inconsistent()),
        }
        if let (Some(x), Some(y)) = (
            exact_quotient(&point[0], &point[2], hints[0]),
            exact_quotient(&point[1], &point[2], hints[1]),
        ) {
            return self.push_flat(Point2::new(x, y));
        }
        let enclosure = [
            point[0].enclosure(),
            point[1].enclosure(),
            point[2].enclosure(),
        ];
        let (x, y) = (
            enclosure[0].quotient(enclosure[2]),
            enclosure[1].quotient(enclosure[2]),
        );
        self.points.push(ExactPoint {
            exact: point,
            enclosure,
            bounds: [x.lo(), y.lo(), x.hi(), y.hi()],
            flat: None,
        });
        Ok(PointId(self.points.len() - 1))
    }

    /// The coordinates of a point that is a binary64 position (`w = 1`).
    pub fn flat(&self, point: PointId) -> Option<Point2> {
        self.points[point.0].flat
    }

    /// The exact homogeneous coordinates of a point, `w > 0`.
    pub fn exact(&self, point: PointId) -> &[Dyadic; 3] {
        &self.points[point.0].exact
    }

    /// The point where the segments `ends[0]`-`ends[1]` and
    /// `ends[2]`-`ends[3]` cross properly.
    ///
    /// The lines `L = a x b` and the meeting point `X = L1 x L2` have degrees
    /// 2 and (3, 3, 2) in the coordinates; `w` is made positive. When the
    /// point is a binary64 position it takes the fast path. Parallel segments
    /// have no crossing.
    pub fn segment_crossing(&mut self, ends: [PointId; 4]) -> GeomResult<PointId> {
        let [a, b, c, d] = ends.map(|end| &self.points[end.0]);
        let point = cross(&cross(&a.exact, &b.exact), &cross(&c.exact, &d.exact));
        // An edge that does not move in a coordinate crosses at that value.
        let fixed = |axis: usize| -> Option<f64> {
            [(a, b), (c, d)].into_iter().find_map(|(from, to)| {
                let (from, to) = (coordinate(from.flat?, axis), coordinate(to.flat?, axis));
                (from == to).then_some(from)
            })
        };
        let hints = [fixed(0), fixed(1)];
        self.push_hinted(point, hints)
    }

    /// The parameter `t` along the first of two crossing segments (`a` to `b`
    /// and `c` to `d`, as in [`ExactPoints::segment_crossing`]) at which they
    /// cross: `a + t (b - a)`. For output only.
    pub fn segment_parameter(&self, ends: [PointId; 4]) -> GeomResult<f64> {
        let [a, b, c, d] = ends.map(|end| &self.points[end.0].exact);
        let difference = |to: &Homogeneous<Dyadic>, from: &Homogeneous<Dyadic>| {
            [to[0].sub(&from[0]), to[1].sub(&from[1])]
        };
        let turn = |p: [Dyadic; 2], q: [Dyadic; 2]| p[0].mul(&q[1]).sub(&p[1].mul(&q[0]));
        let along = difference(d, c);
        let numerator = turn(difference(c, a), along.clone());
        let denominator = turn(difference(b, a), along);
        if denominator.sign() == Some(Sign::Zero) {
            return Err(inconsistent());
        }
        let (m_n, e_n) = numerator.approx_parts();
        let (m_d, e_d) = denominator.approx_parts();
        let shift = (e_n - e_d).clamp(-1100, 1100) as i32;
        let t = m_n / m_d * 2f64.powi(shift / 2) * 2f64.powi(shift - shift / 2);
        if t.is_finite() {
            Ok(t)
        } else {
            Err(not_finite())
        }
    }

    /// Exact orientation of `a`, `b`, `c`: positive when they turn
    /// counter-clockwise in the store's sense (see [`ExactPoints`]).
    pub fn orient(&self, a: PointId, b: PointId, c: PointId) -> GeomResult<Sign> {
        let (a, b, c) = (&self.points[a.0], &self.points[b.0], &self.points[c.0]);
        if let (Some(fa), Some(fb), Some(fc)) = (a.flat, b.flat, c.flat) {
            let sign = orient2d(fa, fb, fc).sign().ok_or_else(undecided)?;
            return Ok(sign_times(sign, self.frame_sign));
        }
        let sign = match determinant(&a.enclosure, &b.enclosure, &c.enclosure).sign() {
            Some(sign) => sign,
            None => {
                self.escalations.set(self.escalations.get() + 1);
                determinant(&a.exact, &b.exact, &c.exact)
                    .sign()
                    .ok_or_else(undecided)?
            }
        };
        Ok(sign_times(sign, self.frame_sign))
    }

    /// Order of the coordinates of `a` and `b` along `axis` (0 = x, 1 = y).
    /// The stored axes, not the frame's.
    pub fn compare(&self, a: PointId, b: PointId, axis: usize) -> GeomResult<Ordering> {
        let (a, b) = (&self.points[a.0], &self.points[b.0]);
        if let (Some(fa), Some(fb)) = (a.flat, b.flat) {
            return coordinate(fa, axis)
                .partial_cmp(&coordinate(fb, axis))
                .ok_or_else(not_finite);
        }
        let sign = match difference(&a.enclosure, &b.enclosure, axis).sign() {
            Some(sign) => sign,
            None => {
                self.escalations.set(self.escalations.get() + 1);
                difference(&a.exact, &b.exact, axis)
                    .sign()
                    .ok_or_else(undecided)?
            }
        };
        Ok(ordering(sign))
    }

    /// Whether `a` and `b` are the same point.
    pub fn same(&self, a: PointId, b: PointId) -> GeomResult<bool> {
        Ok(a == b
            || (self.compare(a, b, 0)? == Ordering::Equal
                && self.compare(a, b, 1)? == Ordering::Equal))
    }

    /// Order by `x`, then `y`.
    pub fn lexicographic(&self, a: PointId, b: PointId) -> GeomResult<Ordering> {
        match self.compare(a, b, 0)? {
            Ordering::Equal => self.compare(a, b, 1),
            order => Ok(order),
        }
    }

    /// `[x_lo, y_lo, x_hi, y_hi]`, a sound box of the point.
    pub fn bounds(&self, point: PointId) -> [f64; 4] {
        self.points[point.0].bounds
    }

    /// The axis along which two distinct points of a line differ.
    pub fn dominant(&self, a: PointId, b: PointId) -> GeomResult<usize> {
        Ok(if self.compare(a, b, 0)? == Ordering::Equal {
            1
        } else {
            0
        })
    }

    /// `(low, high)` of two points along `axis`.
    pub fn ordered(&self, a: PointId, b: PointId, axis: usize) -> GeomResult<(PointId, PointId)> {
        Ok(if self.compare(a, b, axis)? == Ordering::Greater {
            (b, a)
        } else {
            (a, b)
        })
    }

    /// Whether collinear `at` lies strictly between `before` and `after`.
    pub fn strictly_between(
        &self,
        before: PointId,
        at: PointId,
        after: PointId,
    ) -> GeomResult<bool> {
        let axis = self.dominant(before, after)?;
        let (low, high) = self.ordered(before, after, axis)?;
        Ok(self.compare(low, at, axis)? == Ordering::Less
            && self.compare(at, high, axis)? == Ordering::Less)
    }

    /// Whether collinear `point` lies on the closed segment `start`-`end`.
    pub fn within_segment(&self, point: PointId, start: PointId, end: PointId) -> GeomResult<bool> {
        let axis = self.dominant(start, end)?;
        let (low, high) = self.ordered(start, end, axis)?;
        Ok(self.compare(low, point, axis)? != Ordering::Greater
            && self.compare(point, high, axis)? != Ordering::Greater)
    }

    /// Where `point` lies relative to the closed `ring`: a crossing-number
    /// walk with exact orientation tests.
    ///
    /// `orient` is in the frame's sense and the `y` comparisons are in the
    /// stored axes. When the frame mirrors them the walk counts the crossings
    /// of the opposite ray, which is the winding number negated; only zero or
    /// nonzero is used.
    pub fn locate(&self, point: PointId, ring: &[PointId]) -> GeomResult<Location> {
        let mut winding = 0i32;
        for index in 0..ring.len() {
            let (a, b) = (ring[index], ring[(index + 1) % ring.len()]);
            let side = self.orient(a, b, point)?;
            if side == Sign::Zero && self.within_segment(point, a, b)? {
                return Ok(Location::Boundary);
            }
            if self.compare(a, point, 1)? != Ordering::Greater {
                if self.compare(b, point, 1)? == Ordering::Greater && side == Sign::Positive {
                    winding += 1;
                }
            } else if self.compare(b, point, 1)? != Ordering::Greater && side == Sign::Negative {
                winding -= 1;
            }
        }
        Ok(if winding == 0 {
            Location::Outside
        } else {
            Location::Inside
        })
    }

    /// The boxes of a ring of exact points.
    pub fn ring_box(&self, ring: &[PointId]) -> RingBox {
        let mut bounds = [
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ];
        let mut core = bounds;
        for &node in ring {
            let [x_lo, y_lo, x_hi, y_hi] = self.bounds(node);
            bounds = [
                bounds[0].min(x_lo),
                bounds[1].min(y_lo),
                bounds[2].max(x_hi),
                bounds[3].max(y_hi),
            ];
            core = [
                core[0].min(x_hi),
                core[1].min(y_hi),
                core[2].max(x_lo),
                core[3].max(y_lo),
            ];
        }
        RingBox { bounds, core }
    }

    /// Each coordinate `x / w` and `y / w` rounded once to the nearest
    /// binary64 number, ties to even. `Degenerate` when one is beyond the
    /// binary64 range.
    pub fn round(&self, point: PointId) -> GeomResult<Point2> {
        let point = &self.points[point.0];
        if let Some(flat) = point.flat {
            return Ok(flat);
        }
        let [x, y, w] = &point.exact;
        match (round_ratio(x, w), round_ratio(y, w)) {
            (Some(x), Some(y)) => Ok(Point2::new(x, y)),
            _ => Err(not_finite()),
        }
    }

    /// The number of predicate calls that left the interval filter for exact
    /// arithmetic. Not counted by the flat fast path.
    pub fn escalations(&self) -> usize {
        self.escalations.get()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flats(points: &mut ExactPoints, corners: &[(f64, f64)]) -> Vec<PointId> {
        corners
            .iter()
            .map(|&(x, y)| points.push_flat(Point2::new(x, y)).expect("point"))
            .collect()
    }

    #[test]
    fn the_diagonals_of_the_unit_square_meet_at_its_centre_on_the_fast_path() {
        let mut points = ExactPoints::new();
        let at = flats(
            &mut points,
            &[(0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0), (0.5, 0.5)],
        );
        let crossing = points
            .segment_crossing([at[0], at[1], at[2], at[3]])
            .expect("crossing");
        assert_eq!(points.flat(crossing), Some(Point2::new(0.5, 0.5)));
        assert!(points.same(crossing, at[4]).expect("same"));
    }

    #[test]
    fn a_segment_crossing_does_not_depend_on_the_direction_of_its_segments() {
        let mut points = ExactPoints::new();
        let at = flats(
            &mut points,
            &[(0.0, 0.0), (0.5, 1.0), (1.0, 0.0), (0.0, 1.0)],
        );
        let forward = points
            .segment_crossing([at[0], at[1], at[2], at[3]])
            .expect("crossing");
        let swapped = points
            .segment_crossing([at[2], at[3], at[0], at[1]])
            .expect("crossing");
        let reversed = points
            .segment_crossing([at[1], at[0], at[3], at[2]])
            .expect("crossing");
        assert_eq!(points.flat(forward), None);
        let exact = points.exact(forward);
        // (1/3, 2/3): the numerators differ by a factor of two.
        assert_eq!(exact[0].add(&exact[0]), exact[1]);
        assert_eq!(exact[2].sign(), Some(Sign::Positive));
        for other in [swapped, reversed] {
            assert!(points.same(forward, other).expect("same"));
            assert_eq!(points.exact(forward), points.exact(other));
        }
    }

    #[test]
    fn crossings_on_one_segment_are_collinear_in_either_frame_sense() {
        for (frame_sign, turn) in [
            (Sign::Positive, Sign::Positive),
            (Sign::Negative, Sign::Negative),
        ] {
            let mut points = ExactPoints::with_frame_sign(frame_sign);
            // The segment 0-1 from (0.1, 0.1) to (3.7, 2.9), crossed by three
            // segments whose crossing points are not binary64.
            let mut corners = vec![(0.1, 0.1), (3.7, 2.9), (2.0, 0.0)];
            for x in [0.7, 1.3, 2.9] {
                corners.extend([(x, -1.0), (x + 0.3, 4.0)]);
            }
            let at = flats(&mut points, &corners);
            let crossing: Vec<PointId> = (0..3)
                .map(|k| {
                    points
                        .segment_crossing([at[0], at[1], at[3 + 2 * k], at[4 + 2 * k]])
                        .expect("crossing")
                })
                .collect();
            assert!(crossing.iter().all(|&p| points.flat(p).is_none()));
            assert_eq!(
                points
                    .orient(crossing[0], crossing[1], crossing[2])
                    .expect("orient"),
                Sign::Zero
            );
            // (2, 0) lies below the segment, to its right going from the
            // first crossing to the second.
            let side = points
                .orient(crossing[0], crossing[1], at[2])
                .expect("orient");
            assert_eq!(side, sign_times(Sign::Negative, turn));
        }
    }

    #[test]
    fn a_mirrored_frame_reverses_the_orientation() {
        let corners = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
        let mut counter = ExactPoints::new();
        let mut mirrored = ExactPoints::with_frame_sign(Sign::Negative);
        let (a, b) = (
            flats(&mut counter, &corners),
            flats(&mut mirrored, &corners),
        );
        assert_eq!(counter.orient(a[0], a[1], a[2]).unwrap(), Sign::Positive);
        assert_eq!(mirrored.orient(b[0], b[1], b[2]).unwrap(), Sign::Negative);
    }

    #[test]
    fn the_parameter_of_a_segment_crossing_is_its_position_along_the_first_segment() {
        let mut points = ExactPoints::new();
        let at = flats(
            &mut points,
            &[(0.0, 0.0), (0.5, 1.0), (1.0, 0.0), (0.0, 1.0)],
        );
        let t = points
            .segment_parameter([at[0], at[1], at[2], at[3]])
            .expect("parameter");
        assert!((t - 2.0 / 3.0).abs() < 1e-15, "{t}");
        let back = points
            .segment_parameter([at[2], at[3], at[0], at[1]])
            .expect("parameter");
        assert!((back - 2.0 / 3.0).abs() < 1e-15, "{back}");
        assert!(points
            .segment_parameter([at[0], at[1], at[0], at[1]])
            .is_err());
    }

    #[test]
    fn parallel_segments_have_no_crossing() {
        let mut points = ExactPoints::new();
        let at = flats(
            &mut points,
            &[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        );
        assert!(points
            .segment_crossing([at[0], at[1], at[2], at[3]])
            .is_err());
    }

    #[test]
    fn a_homogeneous_point_with_a_binary64_position_takes_the_fast_path() {
        let mut points = ExactPoints::new();
        let three = Dyadic::try_from_f64(3.0).unwrap();
        let (x, y) = (
            Dyadic::try_from_f64(1.5).unwrap(),
            Dyadic::try_from_f64(0.75).unwrap(),
        );
        // (3 * 0.5, 3 * 0.25, 3) is (0.5, 0.25); a negative weight is the same.
        let half = points
            .push([
                x.mul(&Dyadic::try_from_f64(1.0).unwrap()),
                y.mul(&Dyadic::try_from_f64(1.0).unwrap()),
                three.clone(),
            ])
            .expect("point");
        let negated = points.push([x.neg(), y.neg(), three.neg()]).expect("point");
        assert_eq!(points.flat(half), Some(Point2::new(0.5, 0.25)));
        assert!(points.same(half, negated).expect("same"));
        assert!(points.push([x, y, Dyadic::zero()]).is_err());
    }

    #[test]
    fn a_crossing_rounds_once_to_the_nearest_binary64() {
        let mut points = ExactPoints::new();
        let at = flats(
            &mut points,
            &[(0.0, 0.0), (0.5, 1.0), (1.0, 0.0), (0.0, 1.0)],
        );
        let crossing = points
            .segment_crossing([at[0], at[1], at[2], at[3]])
            .expect("crossing");
        let rounded = points.round(crossing).expect("rounded");
        assert_eq!(rounded.x.to_bits(), (1.0f64 / 3.0).to_bits());
        assert_eq!(rounded.y.to_bits(), (2.0f64 / 3.0).to_bits());
        assert_eq!(points.round(at[1]).unwrap(), Point2::new(0.5, 1.0));
    }

    #[test]
    fn locate_tells_inside_outside_and_boundary() {
        let mut points = ExactPoints::new();
        let at = flats(
            &mut points,
            &[
                (0.0, 0.0),
                (2.0, 0.0),
                (2.0, 2.0),
                (0.0, 2.0),
                (1.0, 1.0),
                (3.0, 1.0),
                (1.0, 0.0),
            ],
        );
        let ring = &at[..4];
        assert_eq!(points.locate(at[4], ring).unwrap(), Location::Inside);
        assert_eq!(points.locate(at[5], ring).unwrap(), Location::Outside);
        assert_eq!(points.locate(at[6], ring).unwrap(), Location::Boundary);
        assert!(points
            .ring_box(ring)
            .possibly_within(&points.ring_box(ring)));
    }

    #[test]
    fn escalations_count_the_calls_the_interval_filter_could_not_decide() {
        let mut points = ExactPoints::new();
        let at = flats(
            &mut points,
            &[(0.0, 0.0), (0.5, 1.0), (1.0, 0.0), (0.0, 1.0)],
        );
        let crossing = points
            .segment_crossing([at[0], at[1], at[2], at[3]])
            .expect("crossing");
        assert_eq!(points.escalations(), 0);
        // The crossing lies exactly on its own segment: the enclosure
        // straddles zero and the exact value decides.
        assert_eq!(
            points.orient(at[0], crossing, at[1]).expect("orient"),
            Sign::Zero
        );
        assert!(points.escalations() >= 1);
    }
}
