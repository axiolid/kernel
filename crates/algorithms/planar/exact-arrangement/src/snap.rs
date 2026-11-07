//! Snap rounding of exact edges to the binary64 grid (ADR 0100).
//!
//! A *hot cell* is the set of points that [`ExactPoints::round`] rounds to
//! the rounded position of an end of an input edge. The rounding is to
//! nearest with ties to even, so the cells are half-open boxes of the
//! binary64 tensor grid that partition the plane: the cell of an even centre
//! owns both of its ties.
//!
//! [`snap_round`] routes every edge through the centres of the hot cells that
//! some point of the closed exact edge rounds to, in order along the edge,
//! and arranges the routed fragments exactly. Where nothing but the edge's
//! own ends rounds onto an edge, the fragment is the edge with its ends
//! rounded once. The result is certified, not assumed: a crossing left among
//! the fragments is refused.
//!
//! The displacement bound is the larger half-extent of the cells a route
//! visits (zero for a binary64 end), so a point farther from the exact
//! boundary than the bound has the same winding number before and after.
//!
//! [`Snapped::class_boundary`] nets the routes of a class's edges per
//! fragment: opposite traversals cancel, so a part that became zero-width
//! vanishes and a zero-width gap merges its two sides.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use axiolid_contracts::{DisplacementBound, ExecutionOptions, GeomError, GeomResult, Sign};
use axiolid_core::Point2;
use axiolid_exact::{Arith, Dyadic};

use crate::arrangement::{
    arrange, charge, sort_exact, Arrangement, PairBudget, PieceId, Segment, VertexId,
};
use crate::points::{inconsistent, ExactPoints, PointId};
use crate::ring::point_order;
use crate::trace::HalfEdge;

/// Bytes charged for one hot cell: its exact bounds, the point and its
/// index entries.
pub const HOT_POINT_SCRATCH_BYTES: usize = 256;

/// Bytes charged for one routed fragment: the segment and its path entry.
pub const FRAGMENT_SCRATCH_BYTES: usize = 64;

fn exact(value: f64) -> GeomResult<Dyadic> {
    Dyadic::try_from_f64(value).ok_or_else(|| GeomError::Degenerate("a cell is not finite".into()))
}

/// The cell of one centre: per axis its exact bounds, whether it owns them,
/// and the extents used for the displacement bound and for pruning.
#[derive(Debug, Clone)]
struct Cell {
    center: Point2,
    /// `[x, y]` lower bounds: the midpoints to the previous binary64 number.
    lo: [Dyadic; 2],
    /// `[x, y]` upper bounds.
    hi: [Dyadic; 2],
    /// Per axis, whether the cell owns its bounds (an even centre).
    closed: [bool; 2],
    /// The largest distance from the centre to a point of the cell, per
    /// coordinate, over both axes.
    half_width: f64,
    /// `[x_lo, y_lo, x_hi, y_hi]` from the previous to the next binary64
    /// number; contains the cell. Only prunes.
    reach: [f64; 4],
}

/// The gaps from a finite `value` down and up to its neighbours; at the
/// largest magnitude the missing neighbour mirrors the finite gap.
fn gaps(value: f64) -> (f64, f64) {
    let (down, up) = (value.next_down(), value.next_up());
    match (down.is_finite(), up.is_finite()) {
        (true, true) => (value - down, up - value),
        (false, _) => (up - value, up - value),
        (_, false) => (value - down, value - down),
    }
}

fn cell(center: Point2) -> GeomResult<Cell> {
    let center = Point2::new(center.x + 0.0, center.y + 0.0);
    let half = exact(0.5)?;
    let mut lo = [Dyadic::zero(), Dyadic::zero()];
    let mut hi = [Dyadic::zero(), Dyadic::zero()];
    let mut closed = [false; 2];
    let mut reach = [0.0; 4];
    let mut half_width = 0.0f64;
    for axis in 0..2 {
        let c = if axis == 0 { center.x } else { center.y };
        let (down, up) = gaps(c);
        let c_exact = exact(c)?;
        lo[axis] = c_exact.sub(&exact(down)?.mul(&half));
        hi[axis] = c_exact.add(&exact(up)?.mul(&half));
        closed[axis] = c.abs().to_bits() & 1 == 0;
        reach[axis] = c - down;
        reach[axis + 2] = c + up;
        let gap = up.max(down);
        let halved = gap * 0.5;
        half_width = half_width.max(if halved > 0.0 { halved } else { gap });
    }
    Ok(Cell {
        center,
        lo,
        hi,
        closed,
        half_width,
        reach,
    })
}

/// A parameter `num / den` along an edge (`den > 0`), open when the edge
/// meets the cell only after it.
#[derive(Debug, Clone)]
struct Param {
    num: Dyadic,
    den: Dyadic,
    open: bool,
}

fn compare_values(a: &Param, b: &Param) -> GeomResult<Ordering> {
    let sign = a
        .num
        .mul(&b.den)
        .sub(&b.num.mul(&a.den))
        .sign()
        .ok_or_else(inconsistent)?;
    Ok(match sign {
        Sign::Positive => Ordering::Greater,
        Sign::Negative => Ordering::Less,
        _ => Ordering::Equal,
    })
}

/// Order by value, a closed parameter before an open one of equal value.
fn compare_params(a: &Param, b: &Param) -> GeomResult<Ordering> {
    Ok(compare_values(a, b)?.then(a.open.cmp(&b.open)))
}

/// Where the closed exact segment `from`-`to` (homogeneous, `w > 0`) first
/// meets `cell`: the lower end of the parameter interval of the meeting, or
/// `None` when it does not meet it.
fn meets(from: &[Dyadic; 3], to: &[Dyadic; 3], cell: &Cell) -> GeomResult<Option<Param>> {
    let mut lower = Param {
        num: Dyadic::zero(),
        den: exact(1.0)?,
        open: false,
    };
    let mut upper = Param {
        num: exact(1.0)?,
        den: exact(1.0)?,
        open: false,
    };
    for axis in 0..2 {
        let open = !cell.closed[axis];
        // `sigma * (P_axis - bound * P_w)` is at least zero on the cell's
        // side of the bound; it is linear along the segment.
        for (bound, negate) in [(&cell.lo[axis], false), (&cell.hi[axis], true)] {
            let at = |p: &[Dyadic; 3]| {
                let value = p[axis].sub(&bound.mul(&p[2]));
                if negate {
                    value.neg()
                } else {
                    value
                }
            };
            let (h0, h1) = (at(from), at(to));
            let satisfied = |h: &Dyadic| -> GeomResult<bool> {
                Ok(match h.sign().ok_or_else(inconsistent)? {
                    Sign::Positive => true,
                    Sign::Zero => !open,
                    _ => false,
                })
            };
            match (satisfied(&h0)?, satisfied(&h1)?) {
                (true, true) => {}
                (false, false) => return Ok(None),
                (true, false) => {
                    let candidate = Param {
                        num: h0.clone(),
                        den: h0.sub(&h1),
                        open,
                    };
                    // An open end excludes its value: it wins a tie, also
                    // against the initial closed end at the edge's far end.
                    match compare_values(&candidate, &upper)? {
                        Ordering::Less => upper = candidate,
                        Ordering::Equal if candidate.open => upper = candidate,
                        _ => {}
                    }
                }
                (false, true) => {
                    let candidate = Param {
                        num: h0.neg(),
                        den: h1.sub(&h0),
                        open,
                    };
                    match compare_values(&candidate, &lower)? {
                        Ordering::Greater => lower = candidate,
                        Ordering::Equal if candidate.open => lower = candidate,
                        _ => {}
                    }
                }
            }
        }
    }
    // `lower <= s <= upper`, each end excluded when open.
    Ok(match compare_values(&lower, &upper)? {
        Ordering::Less => Some(lower),
        Ordering::Equal if !lower.open && !upper.open => Some(lower),
        _ => None,
    })
}

/// The edges of a snap-rounded set, as fragments of an exact arrangement.
#[derive(Debug, Clone)]
pub struct Snapped {
    /// The hot centres, all flat, with the frame sign of the input points
    /// (the arrangement may add points only if it was refused).
    pub points: ExactPoints,
    /// The arrangement of the fragments; no crossings.
    pub arrangement: Arrangement,
    /// Per input edge, its pieces from `ends[0]` to `ends[1]`, each with
    /// whether that direction runs from `ends[0]` to `ends[1]` of the piece.
    /// Empty when both ends round to the same cell.
    pub routes: Vec<Vec<(PieceId, bool)>>,
    /// How far a boundary point moved at most, in the maximum norm.
    pub displacement: DisplacementBound,
    /// The edges routed through a cell other than their ends' cells.
    pub rerouted: usize,
}

/// The net boundary of a class after snapping.
#[derive(Debug, Clone)]
pub struct ClassBoundary {
    /// The half-edges with the class on their left, ascending by piece.
    pub half_edges: Vec<HalfEdge>,
    /// Per element of the chain, whether some piece of its route carries a
    /// surviving half-edge.
    pub survives: Vec<bool>,
}

fn refused(message: &str) -> GeomError {
    GeomError::Degenerate(message.into())
}

/// Snap-round `edges`, straight segments between distinct points of
/// `points`.
///
/// The edges must be pairwise interior-disjoint with no end inside another
/// (a boundary subgraph of an arrangement). Every examined edge and cell
/// pair counts against `budget`; cells and fragments are charged against
/// the memory budget on top of `declared_bytes`, as in [`arrange`].
pub fn snap_round(
    points: &ExactPoints,
    edges: &[Segment],
    budget: &mut PairBudget,
    declared_bytes: usize,
    options: &ExecutionOptions,
) -> GeomResult<Snapped> {
    // The hot cells: the distinct rounded ends, by centre.
    let rounded = |id: PointId| -> GeomResult<Point2> {
        let at = points.round(id)?;
        Ok(Point2::new(at.x + 0.0, at.y + 0.0))
    };
    let mut centers: Vec<Point2> = Vec::new();
    for edge in edges {
        for end in edge.ends {
            centers.push(rounded(end)?);
        }
    }
    centers.sort_by(point_order);
    centers.dedup_by(|a, b| point_order(a, b).is_eq());
    let mut charged = 0usize;
    let mut store = ExactPoints::with_frame_sign(points.frame_sign());
    let mut cells = Vec::with_capacity(centers.len());
    let mut ids = Vec::with_capacity(centers.len());
    for &center in &centers {
        charge(
            options,
            declared_bytes,
            &mut charged,
            HOT_POINT_SCRATCH_BYTES,
        )?;
        let hot = cell(center)?;
        ids.push(store.push_flat(hot.center)?);
        cells.push(hot);
    }
    let find = |at: Point2| {
        centers
            .binary_search_by(|probe| point_order(probe, &at))
            .map_err(|_| inconsistent())
    };

    // Per edge the route: the start cell, the cells met in order, the end.
    let mut routes: Vec<Vec<usize>> = Vec::with_capacity(edges.len());
    let mut displacement = 0.0f64;
    let mut rerouted = 0usize;
    for edge in edges {
        options.check_cancelled()?;
        let [from, to] = edge.ends;
        let start = find(rounded(from)?)?;
        let end = find(rounded(to)?)?;
        if start == end {
            // Both ends in one convex cell: the whole edge is inside it.
            for (id, at) in [(from, start), (to, end)] {
                if points.flat(id).is_none() {
                    displacement = displacement.max(cells[at].half_width);
                }
            }
            routes.push(Vec::new());
            continue;
        }
        let (a, b) = (points.bounds(from), points.bounds(to));
        let reach = [
            a[0].min(b[0]),
            a[1].min(b[1]),
            a[2].max(b[2]),
            a[3].max(b[3]),
        ];
        let first = cells.partition_point(|c| c.reach[2] < reach[0]);
        let mut found: Vec<(usize, Param)> = Vec::new();
        for (index, candidate) in cells.iter().enumerate().skip(first) {
            if candidate.reach[0] > reach[2] {
                break;
            }
            if index == start
                || index == end
                || candidate.reach[3] < reach[1]
                || candidate.reach[1] > reach[3]
            {
                continue;
            }
            budget.count()?;
            if let Some(param) = meets(points.exact(from), points.exact(to), candidate)? {
                found.push((index, param));
            }
        }
        let mut order: Vec<usize> = (0..found.len()).collect();
        sort_exact(&mut order, |l, r| compare_params(&found[l].1, &found[r].1))?;
        for pair in order.windows(2) {
            if compare_params(&found[pair[0]].1, &found[pair[1]].1)?.is_eq() {
                return Err(inconsistent());
            }
        }
        for (id, at) in [(from, start), (to, end)] {
            if points.flat(id).is_none() {
                displacement = displacement.max(cells[at].half_width);
            }
        }
        for &(index, _) in &found {
            displacement = displacement.max(cells[index].half_width);
        }
        rerouted += usize::from(!found.is_empty());
        let mut route = Vec::with_capacity(found.len() + 2);
        route.push(start);
        route.extend(order.iter().map(|&i| found[i].0));
        route.push(end);
        routes.push(route);
    }

    // The fragments between consecutive centres, arranged exactly.
    let mut segments = Vec::new();
    let mut spans = Vec::with_capacity(routes.len());
    for route in &routes {
        let first = segments.len();
        for pair in route.windows(2) {
            charge(
                options,
                declared_bytes,
                &mut charged,
                FRAGMENT_SCRATCH_BYTES,
            )?;
            segments.push(Segment {
                ends: [ids[pair[0]], ids[pair[1]]],
            });
        }
        spans.push(first..segments.len());
    }
    let arrangement = arrange(
        &mut store,
        &segments,
        budget,
        declared_bytes.saturating_add(charged),
        options,
    )?;
    if arrangement.crossings != 0 {
        return Err(refused("snap rounding left two fragments crossing"));
    }
    let routes = spans
        .into_iter()
        .map(|span| {
            span.flat_map(|segment| arrangement.paths[segment].iter().copied())
                .collect()
        })
        .collect();
    Ok(Snapped {
        points: store,
        arrangement,
        routes,
        displacement: DisplacementBound::new(displacement).ok_or_else(inconsistent)?,
        rerouted,
    })
}

impl Snapped {
    /// The binary64 position of a vertex of the arrangement.
    pub fn point(&self, vertex: VertexId) -> GeomResult<Point2> {
        self.points
            .flat(self.arrangement.vertices[vertex.0].point)
            .ok_or_else(inconsistent)
    }

    /// The net boundary of a class whose boundary is `chain`: elements
    /// `(edge, forward)` of the input edges, `forward` when the class runs
    /// along the edge from `ends[0]` to `ends[1]` with the class on its
    /// left. Opposite traversals of a piece cancel. A piece traversed twice
    /// the same way is refused.
    pub fn class_boundary(&self, chain: &[(usize, bool)]) -> GeomResult<ClassBoundary> {
        let route_of = |&(edge, forward): &(usize, bool)| -> Vec<(PieceId, bool)> {
            let route = &self.routes[edge];
            if forward {
                route.clone()
            } else {
                route.iter().rev().map(|&(p, f)| (p, !f)).collect()
            }
        };
        let mut net: BTreeMap<PieceId, i32> = BTreeMap::new();
        for element in chain {
            for (piece, forward) in route_of(element) {
                *net.entry(piece).or_default() += if forward { 1 } else { -1 };
            }
        }
        let mut half_edges = Vec::new();
        for (&piece, &count) in &net {
            match count {
                0 => {}
                1 => half_edges.push(HalfEdge {
                    piece,
                    forward: true,
                }),
                -1 => half_edges.push(HalfEdge {
                    piece,
                    forward: false,
                }),
                _ => return Err(refused("snap rounding folded a boundary onto itself")),
            }
        }
        let survives = chain
            .iter()
            .map(|element| {
                route_of(element)
                    .iter()
                    .any(|(piece, _)| net.get(piece).is_some_and(|&count| count != 0))
            })
            .collect();
        Ok(ClassBoundary {
            half_edges,
            survives,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat(x: f64, y: f64) -> Point2 {
        Point2::new(x, y)
    }

    #[test]
    fn the_cell_at_one_is_asymmetric_by_a_binade() {
        let c = cell(flat(1.0, 1.0)).unwrap();
        let one = exact(1.0).unwrap();
        // Below 1.0 the spacing is 2^-53, above it 2^-52.
        assert_eq!(
            one.sub(&c.lo[0]),
            exact(2f64.powi(-54)).unwrap(),
            "lower half-width"
        );
        assert_eq!(
            c.hi[0].sub(&one),
            exact(2f64.powi(-53)).unwrap(),
            "upper half-width"
        );
        assert_eq!(c.half_width, 2f64.powi(-53));
        assert!(c.closed[0], "1.0 has an even significand");
    }

    #[test]
    fn cells_at_zero_the_smallest_subnormal_and_the_largest_number() {
        let tiny = f64::from_bits(1);
        let zero = cell(flat(0.0, -0.0)).unwrap();
        assert_eq!(zero.center.y.to_bits(), 0, "-0.0 is +0.0");
        assert_eq!(zero.hi[0], exact(tiny).unwrap().mul(&exact(0.5).unwrap()));
        assert_eq!(zero.lo[0], zero.hi[0].neg());
        assert!(zero.closed[0]);
        assert_eq!(zero.half_width, tiny, "halving underflows");

        let small = cell(flat(tiny, tiny)).unwrap();
        assert!(!small.closed[0], "the smallest subnormal is odd");
        assert_eq!(small.lo[0], zero.hi[0]);

        let big = cell(flat(f64::MAX, -f64::MAX)).unwrap();
        let gap = f64::MAX - f64::MAX.next_down();
        assert_eq!(
            big.hi[0].sub(&exact(f64::MAX).unwrap()),
            exact(gap).unwrap().mul(&exact(0.5).unwrap()),
            "the missing neighbour mirrors the gap"
        );
        assert_eq!(big.reach[2], f64::MAX + gap);
    }

    fn segment_through(points: &ExactPoints, a: (f64, f64), b: (f64, f64), cell: &Cell) -> bool {
        let mut store = points.clone();
        let from = store.push_flat(flat(a.0, a.1)).unwrap();
        let to = store.push_flat(flat(b.0, b.1)).unwrap();
        meets(store.exact(from), store.exact(to), cell)
            .unwrap()
            .is_some()
    }

    #[test]
    fn a_segment_on_a_midpoint_meets_only_the_even_neighbour() {
        let odd = 1.0 + f64::EPSILON; // 0x3ff0...01: odd significand
        assert_eq!(1.0f64.to_bits() & 1, 0);
        assert_eq!(odd.to_bits() & 1, 1);
        // The vertical line x = tie runs between the cells of 1.0 and odd.
        let tie = exact(1.0)
            .unwrap()
            .add(&exact(odd).unwrap())
            .mul(&exact(0.5).unwrap());
        let (even, odd_cell) = (cell(flat(1.0, 1.0)).unwrap(), cell(flat(odd, 1.0)).unwrap());
        // Homogeneous points on x = tie: (tie * 2, y * 2, 2) and a second.
        let two = exact(2.0).unwrap();
        let line = |y: f64| [tie.mul(&two), exact(y).unwrap().mul(&two), two.clone()];
        let (p, q) = (line(0.9), line(1.1));
        assert!(meets(&p, &q, &even).unwrap().is_some());
        assert!(meets(&p, &q, &odd_cell).unwrap().is_none());
    }

    #[test]
    fn a_point_meets_its_own_cell_and_no_neighbour() {
        let mut points = ExactPoints::new();
        let at = flat(0.3, -2.5);
        let id = points.push_flat(at).unwrap();
        let own = cell(at).unwrap();
        assert!(meets(points.exact(id), points.exact(id), &own)
            .unwrap()
            .is_some());
        for neighbour in [
            flat(at.x.next_up(), at.y),
            flat(at.x.next_down(), at.y),
            flat(at.x, at.y.next_up()),
            flat(at.x, at.y.next_down()),
        ] {
            let other = cell(neighbour).unwrap();
            assert!(meets(points.exact(id), points.exact(id), &other)
                .unwrap()
                .is_none());
        }
    }

    #[test]
    fn a_segment_through_a_cell_is_found_and_one_beside_it_is_not() {
        let points = ExactPoints::new();
        let c = cell(flat(1.0, 1.0)).unwrap();
        assert!(segment_through(&points, (0.0, 0.0), (2.0, 2.0), &c));
        assert!(!segment_through(&points, (0.0, 0.5), (2.0, 2.5), &c));
    }

    /// A small seeded generator; the tests need spread, not quality.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 11
        }

        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// The point at parameter `num / den` of the segment, rounded once.
    fn rounded_at(from: &[Dyadic; 3], to: &[Dyadic; 3], num: &Dyadic, den: &Dyadic) -> Point2 {
        let rest = den.sub(num);
        let mix = |i: usize| rest.mul(&from[i]).add(&num.mul(&to[i]));
        let mut store = ExactPoints::new();
        let id = store.push([mix(0), mix(1), mix(2)]).unwrap();
        let at = store.round(id).unwrap();
        Point2::new(at.x + 0.0, at.y + 0.0)
    }

    /// `meets` against rounding the points of the segment one by one: the
    /// lower end it returns is where the segment enters the cell, nothing
    /// sampled before it rounds into the cell, and every sample that does
    /// lies at or after it.
    #[test]
    fn meets_agrees_with_rounding_the_points_of_the_segment() {
        let mut rng = Lcg(0x5eed);
        let (mut hits, mut misses, mut ties) = (0, 0, 0);
        let tiny = exact(2f64.powi(-200)).unwrap();
        for case in 0..5000 {
            let scale = [1.0, 0.1, 1e-3, 3.7e5, 2.0][rng.below(5) as usize];
            let center = Point2::new(
                (rng.next() as f64 / (1u64 << 53) as f64 * 8.0 - 4.0) * scale,
                (rng.next() as f64 / (1u64 << 53) as f64 * 8.0 - 4.0) * scale,
            );
            let center = Point2::new(center.x + 0.0, center.y + 0.0);
            let cell = cell(center).unwrap();
            let ulp = |c: f64| exact(c.next_up() - c).unwrap();
            let end = |rng: &mut Lcg| {
                let w = [1u64, 3, 5, 7][rng.below(4) as usize];
                let weight = exact(w as f64).unwrap();
                let quarter = exact(0.25).unwrap();
                let offset = |rng: &mut Lcg, c: f64| {
                    let m = rng.below(41) as f64 - 20.0;
                    exact(c)
                        .unwrap()
                        .mul(&weight)
                        .add(&ulp(c).mul(&exact(m).unwrap()).mul(&quarter))
                };
                [offset(rng, center.x), offset(rng, center.y), weight.clone()]
            };
            // A third of the ends sit exactly on a corner of the cell, where
            // ties of both axes meet, whichever way they go.
            let corner = |rng: &mut Lcg| {
                let pick = |rng: &mut Lcg, axis: usize| {
                    if rng.below(2) == 0 {
                        cell.lo[axis].clone()
                    } else {
                        cell.hi[axis].clone()
                    }
                };
                [pick(rng, 0), pick(rng, 1), exact(1.0).unwrap()]
            };
            let either = |rng: &mut Lcg| {
                if rng.below(3) == 0 {
                    corner(rng)
                } else {
                    end(rng)
                }
            };
            let (from, to) = (either(&mut rng), either(&mut rng));
            let found = meets(&from, &to, &cell).unwrap();
            let inside =
                |num: &Dyadic, den: &Dyadic| rounded_at(&from, &to, num, den) == cell.center;
            match &found {
                Some(lower) => {
                    hits += 1;
                    if lower.open {
                        assert!(!inside(&lower.num, &lower.den), "case {case}: open");
                        assert!(
                            inside(&lower.num.add(&lower.den.mul(&tiny)), &lower.den),
                            "case {case}"
                        );
                    } else {
                        assert!(inside(&lower.num, &lower.den), "case {case}: closed");
                    }
                    let before = lower.num.sub(&lower.den.mul(&tiny));
                    if before.sign() == Some(Sign::Positive) {
                        assert!(!inside(&before, &lower.den), "case {case}: before");
                    }
                    ties += usize::from(lower.open);
                }
                None => misses += 1,
            }
            for k in 0..=16 {
                let num = exact(f64::from(k)).unwrap();
                let den = exact(16.0).unwrap();
                if inside(&num, &den) {
                    let lower = found
                        .as_ref()
                        .unwrap_or_else(|| panic!("case {case}: missed s={k}/16"));
                    let at = Param {
                        num,
                        den,
                        open: false,
                    };
                    assert!(
                        compare_values(lower, &at).unwrap().is_le(),
                        "case {case}: s={k}/16 before the lower end"
                    );
                }
            }
        }
        assert!(hits > 500 && misses > 500, "{hits} hits, {misses} misses");
        assert!(ties > 0, "no open lower end was exercised");
    }
}
