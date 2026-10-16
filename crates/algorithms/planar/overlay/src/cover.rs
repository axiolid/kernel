//! Fewer and smaller edge hulls for a Minkowski sum (#323).
//!
//! [`crate::minkowski`] builds `S + P` for a set `S` and a convex piece `P`
//! from one hull `e + P` per boundary edge `e` of `S`. Along a wall of
//! many short edges -- overlay output along tessellated or flush bodies,
//! with sub-millimetre teeth or exactly collinear corners -- each hull
//! overlaps hundreds of its neighbours. Their boundaries cross about that
//! many times each, every piece of the subdivision lies in about that many
//! of them, and time and memory grew with the square of the corners
//! (1,202 corners: 5 s and 1.6 GB). Two exact steps remove most of that.
//!
//! # Collinear edges merged
//!
//! A vertex lying exactly on the line through its neighbours, strictly
//! between them (an exact orientation and coordinate comparisons, no
//! tolerance), only splits one straight edge in two. [`merged`] drops it,
//! so the edge `a -> c` gives one hull: `(ab + P) union (bc + P) = ac + P`.
//! The hull's vertices are vertex sums `a + k` and `c + k`, rounded once
//! as before.
//!
//! # Chords cover the hulls
//!
//! A segment `ab` between two vertices of one ring that lies in the closed
//! set `S` gives a term `ab + P` contained in `S + P`. Across a toothed
//! wall it runs along the valleys, and its hull covers nearly all of every
//! tooth's hull. Each such edge hull then keeps only the vertices outside
//! the chord's hull, plus the first and last of each run of covered ones:
//! the cut-off cap is the hull of covered vertices, so it lies in the
//! chord's hull, and the union of the terms is the same set. Only a sliver
//! about the depth of the teeth is left of each edge hull, and those
//! overlap a few neighbours instead of hundreds.
//!
//! [`chords`] proposes chords greedily: a stack over the ring's vertices
//! pops a vertex bending away from `S` while the chord stays within
//! `depth` of every vertex it spans (an `f64` bound, a choice only). Every
//! chord used is then proved to lie in the closed set, exactly:
//!
//! 1. every vertex it spans lies on the far side of its line or on it,
//!    and one on the chord turns away from `S` (`S` holds both directions
//!    along the chord there);
//! 2. at each end, the chord leaves inside the closed sector of `S`
//!    between the two boundary edges meeting there;
//! 3. no other boundary edge, of any ring of `S`, crosses it properly, and
//!    no other boundary vertex lies on it.
//!
//! Between the points where the chord meets the boundary, which (3) limits
//! to its own vertices, each open stretch lies wholly in `S` or wholly
//! outside it, and (1) and (2) put the start of each in `S`.
//!
//! The trimmed hulls and the chord hulls are rounded exactly as the edge
//! hulls were: each vertex a vertex sum, rounded once, and containment in
//! the chord's hull is decided by exact orientations of those rounded
//! points. A trimmed hull that tidying would change keeps its full hull.

use axiolid_core::Point2;
use axiolid_exact::{certify, Arith, SignExpr};
use axiolid_guarantees::Sign;

use crate::exact_arc::boxes::BoxTree;
use crate::exact_arc::{orient_doubles as orient, Bounds};

/// A ring with its exactly collinear vertices merged.
#[derive(Debug, Clone)]
pub(crate) struct Boundary {
    /// The vertices kept, in ring order.
    pub(crate) points: Vec<Point2>,
    /// Each kept vertex's index in the ring as given.
    origin: Vec<usize>,
}

/// Whether `b`, exactly collinear with `a` and `c`, lies strictly between
/// them.
fn between(a: Point2, b: Point2, c: Point2) -> bool {
    let strictly = |p: f64, q: f64, r: f64| (p < q && q < r) || (p > q && q > r);
    strictly(a.x, b.x, c.x) || (a.x == b.x && b.x == c.x && strictly(a.y, b.y, c.y))
}

/// Whether `p` lies on the closed segment `a b`, exactly.
fn on_segment(a: Point2, b: Point2, p: Point2) -> bool {
    orient(a, b, p) == Sign::Zero
        && p.x >= a.x.min(b.x)
        && p.x <= a.x.max(b.x)
        && p.y >= a.y.min(b.y)
        && p.y <= a.y.max(b.y)
}

/// The ring with every vertex dropped that lies exactly on the segment
/// between its neighbours. Removing them all at once is sound: a run of
/// such vertices lies on one line, in order, between the two kept
/// vertices around it.
pub(crate) fn merged(ring: &[Point2]) -> Boundary {
    let n = ring.len();
    let keep = |i: usize| {
        let (a, b, c) = (ring[(i + n - 1) % n], ring[i], ring[(i + 1) % n]);
        !(orient(a, b, c) == Sign::Zero && between(a, b, c))
    };
    let origin: Vec<usize> = (0..n).filter(|&i| keep(i)).collect();
    if origin.len() < 3 {
        return Boundary {
            points: ring.to_vec(),
            origin: (0..n).collect(),
        };
    }
    Boundary {
        points: origin.iter().map(|&i| ring[i]).collect(),
        origin,
    }
}

/// Every boundary edge of the set, boxed, for the chords' exact checks.
pub(crate) struct Index<'a> {
    rings: &'a [Vec<Point2>],
    /// `(ring, edge)`: the edge leaving `rings[ring][edge]`.
    edges: Vec<(usize, usize)>,
    tree: BoxTree,
}

impl<'a> Index<'a> {
    pub(crate) fn new(rings: &'a [Vec<Point2>]) -> Self {
        let mut edges = Vec::new();
        let mut boxes = Vec::new();
        for (r, ring) in rings.iter().enumerate() {
            for e in 0..ring.len() {
                edges.push((r, e));
                boxes.push(Bounds::around(&[ring[e], ring[(e + 1) % ring.len()]], 0.0));
            }
        }
        Self {
            rings,
            edges,
            tree: BoxTree::new(boxes),
        }
    }
}

/// A chord from kept vertex `from` to kept vertex `to` of a [`Boundary`],
/// spanning the `span` edges after `from`, proved to lie in the set.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Chord {
    pub(crate) from: usize,
    pub(crate) span: usize,
}

/// The distance from `p` to the line through `a` and `b`; infinite when
/// they coincide.
fn off_line(p: Point2, a: Point2, b: Point2) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len = dx.hypot(dy);
    if len == 0.0 {
        return f64::INFINITY;
    }
    ((p.x - a.x) * dy - (p.y - a.y) * dx).abs() / len
}

/// The chords of ring `ring` (its [`Boundary`] `boundary`), each spanning
/// two edges or more, proved to lie in the set. `side` is the side of the
/// ring's travel the set lies on: `Positive` left, `Negative` right.
pub(crate) fn chords(
    ring: usize,
    boundary: &Boundary,
    side: Sign,
    depth: f64,
    index: &Index<'_>,
) -> Vec<Chord> {
    let p = &boundary.points;
    let n = p.len();
    if side == Sign::Zero || n < 4 || depth.is_nan() || depth <= 0.0 {
        return Vec::new();
    }
    let away = side.flip();
    let turn = |i: usize| orient(p[(i + n - 1) % n], p[i], p[(i + 1) % n]);
    // Start at a vertex that bends away from the set: no chord passes it.
    let start = (0..n).find(|&i| turn(i) == away).unwrap_or(0);
    let at = |t: usize| p[(start + t) % n];
    // Positions along the ring from `start`, with a bound on how far the
    // chord ending at each lies from the vertices it spans.
    let mut stack: Vec<(usize, f64)> = vec![(0, 0.0)];
    for t in 1..=n {
        let v = at(t);
        let mut bound = 0.0_f64;
        while stack.len() >= 2 {
            let (u, _) = stack[stack.len() - 2];
            let (w, below) = stack[stack.len() - 1];
            if orient(at(u), at(w), v) == away {
                break;
            }
            let wider = below.max(bound) + off_line(at(w), at(u), v);
            if wider.is_nan() || wider > depth {
                break;
            }
            stack.pop();
            bound = wider;
        }
        stack.push((t, bound));
    }
    stack
        .windows(2)
        .filter(|w| w[1].0 - w[0].0 >= 2 && w[1].0 - w[0].0 < n)
        .map(|w| Chord {
            from: (start + w[0].0) % n,
            span: w[1].0 - w[0].0,
        })
        .filter(|&chord| proved(ring, boundary, chord, side, index))
        .collect()
}

/// Whether, at `at`, the direction to `other` lies in the closed sector
/// swept from the ray to `next` towards `side` until the ray to `prev`.
/// `next` is a vertex the chord spans, so check 1 has put `other` towards
/// `side` of the ray to it, or ahead on it: the sector holds the chord
/// unless it closes before.
fn leaves_inside(at: Point2, next: Point2, prev: Point2, other: Point2, side: Sign) -> bool {
    match orient(at, next, prev) {
        s if s == side => orient(at, prev, other) != side,
        // The boundary doubles back along `next`: merging has removed every
        // vertex it runs straight through.
        Sign::Zero => false,
        _ => true,
    }
}

/// The exact proof that `chord` lies in the closed set (module docs).
fn proved(ring: usize, boundary: &Boundary, chord: Chord, side: Sign, index: &Index<'_>) -> bool {
    let p = &boundary.points;
    let n = p.len();
    let to = (chord.from + chord.span) % n;
    let (a, b) = (p[chord.from], p[to]);
    if a == b {
        return false;
    }
    // 1. The spanned vertices lie on the far side or on the chord, and one
    //    on the chord turns away from the set.
    for k in 1..chord.span {
        let i = (chord.from + k) % n;
        match orient(a, b, p[i]) {
            Sign::Zero => {
                if !on_segment(a, b, p[i])
                    || orient(p[(i + n - 1) % n], p[i], p[(i + 1) % n]) != side.flip()
                {
                    return false;
                }
            }
            s if s == side => return false,
            _ => {}
        }
    }
    // 2. Both ends leave inside the set's sector.
    let next = |i: usize| p[(i + 1) % n];
    let prev = |i: usize| p[(i + n - 1) % n];
    if !leaves_inside(a, next(chord.from), prev(chord.from), b, side)
        || !leaves_inside(b, prev(to), next(to), a, side.flip())
    {
        return false;
    }
    // 3. No other boundary crosses it or has a vertex on it.
    let length = index.rings[ring].len();
    let (first, last) = (boundary.origin[chord.from], boundary.origin[to]);
    let run = (last + length - first) % length;
    let query = Bounds::around(&[a, b], 0.0);
    let (mut found, mut scratch) = (Vec::new(), Vec::new());
    index
        .tree
        .query_into(|q| q.overlaps(&query), &mut found, &mut scratch);
    found.into_iter().all(|item| {
        let (r, e) = index.edges[item];
        let points = &index.rings[r];
        let f = (e + 1) % points.len();
        if r == ring && (e + length - first) % length < run {
            return true;
        }
        let end = |v: usize| r == ring && (v == first || v == last);
        let (c, d) = (points[e], points[f]);
        if (!end(e) && on_segment(a, b, c)) || (!end(f) && on_segment(a, b, d)) {
            return false;
        }
        let opposite = |s: Sign, t: Sign| s != Sign::Zero && t != Sign::Zero && s != t;
        !(opposite(orient(a, b, c), orient(a, b, d)) && opposite(orient(c, d, a), orient(c, d, b)))
    })
}

/// The sign of the cross product `(b - a) x (d - c)`, exactly: zero when
/// the two directions are parallel.
struct Cross {
    a: Point2,
    b: Point2,
    c: Point2,
    d: Point2,
}

impl SignExpr for Cross {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let f = T::from_f64;
        let (ux, uy) = (f(self.b.x).sub(&f(self.a.x)), f(self.b.y).sub(&f(self.a.y)));
        let (vx, vy) = (f(self.d.x).sub(&f(self.c.x)), f(self.d.y).sub(&f(self.c.y)));
        ux.mul(&vy).sub(&uy.mul(&vx)).sign()
    }
}

fn parallel(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    certify(&Cross { a, b, c, d }) == Ok(Sign::Zero)
}

/// Whether `q` lies in the closed convex polygon `polygon`
/// (counter-clockwise).
fn holds(polygon: &[Point2], q: Point2) -> bool {
    (0..polygon.len())
        .all(|i| orient(polygon[i], polygon[(i + 1) % polygon.len()], q) != Sign::Negative)
}

/// The convex `hull` (counter-clockwise, as tidied) less its part inside
/// the convex `cover`, the hull of the chord `a -> b`: `None` when the
/// cover holds it all. Each run of vertices inside the cover keeps its
/// first and last; the cap between them is the hull of covered vertices,
/// so the cover holds it.
///
/// A cut parallel to the chord keeps one more covered vertex. The teeth of
/// a wall are translates of one another along it, so their cuts parallel
/// to it would lie on one line and overlap there, a degenerate pile-up for
/// the subdivision; tilted, they are parallel but apart. A trimmed ring
/// that `tidy` would change is returned whole instead.
pub(crate) fn trim(
    hull: Vec<Point2>,
    cover: &[Point2],
    (a, b): (Point2, Point2),
    tidy: impl Fn(Vec<Point2>) -> Option<Vec<Point2>>,
) -> Option<Vec<Point2>> {
    let n = hull.len();
    let covered: Vec<bool> = hull.iter().map(|&q| holds(cover, q)).collect();
    if covered.iter().all(|&c| c) {
        return None;
    }
    let mut keep: Vec<bool> = (0..n)
        .map(|i| !(covered[i] && covered[(i + n - 1) % n] && covered[(i + 1) % n]))
        .collect();
    let first = (0..n).find(|&i| !covered[i]).unwrap_or(0);
    let mut i = first;
    loop {
        // From a kept vertex, the next kept one.
        let mut j = (i + 1) % n;
        while !keep[j] {
            j = (j + 1) % n;
        }
        if (j + n - i) % n > 2 && parallel(hull[i], hull[j], a, b) {
            keep[(i + 1) % n] = true;
        }
        i = j;
        if i == first {
            break;
        }
    }
    let kept: Vec<Point2> = (0..n).filter(|&i| keep[i]).map(|i| hull[i]).collect();
    if kept.len() == n {
        return Some(hull);
    }
    match tidy(kept.clone()) {
        Some(same) if same == kept => Some(kept),
        _ => Some(hull),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn points(p: &[(f64, f64)]) -> Vec<Point2> {
        p.iter().map(|&(x, y)| Point2::new(x, y)).collect()
    }

    /// A room's bottom wall with two teeth into the room, from `(0, 0)`
    /// to `(4, 0)`, then the given corners back.
    fn room(back: &[(f64, f64)]) -> Vec<(f64, f64)> {
        let mut ring = vec![(0.0, 0.0), (1.0, 0.01), (2.0, 0.0), (3.0, 0.01), (4.0, 0.0)];
        ring.extend_from_slice(back);
        ring
    }

    /// Whether the chord from vertex 0 over `span` edges of the first ring
    /// is proved, the set lying right of its travel (outside a
    /// counter-clockwise room).
    fn proves(rings: &[Vec<(f64, f64)>], span: usize) -> bool {
        let rings: Vec<Vec<Point2>> = rings.iter().map(|r| points(r)).collect();
        let index = Index::new(&rings);
        let boundary = merged(&rings[0]);
        assert_eq!(boundary.points.len(), rings[0].len());
        proved(
            0,
            &boundary,
            Chord { from: 0, span },
            Sign::Negative,
            &index,
        )
    }

    #[test]
    fn merging_drops_only_vertices_exactly_between_their_neighbours() {
        let square = merged(&points(&[
            (0.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (0.0, 1.0),
        ]));
        assert_eq!(
            square.points,
            points(&[(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])
        );
        assert_eq!(square.origin, vec![0, 2, 3, 4]);
        // Off the line by far less than any tolerance, and so kept.
        let near = points(&[(0.0, 0.0), (1.0, 1e-17), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)]);
        assert_eq!(merged(&near).points, near);
        // On the line, but the boundary turns back there.
        let back = points(&[(0.0, 0.0), (2.0, 0.0), (1.0, 0.0), (0.5, 1.0)]);
        assert_eq!(merged(&back).points, back);
    }

    #[test]
    fn a_chord_along_the_valleys_of_a_toothed_wall_is_found_and_proved() {
        let ring = points(&room(&[(4.0, 3.0), (0.0, 3.0)]));
        let rings = vec![ring];
        let index = Index::new(&rings);
        let boundary = merged(&rings[0]);
        let found = chords(0, &boundary, Sign::Negative, 0.05, &index);
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!((found[0].from, found[0].span), (0, 4));
        // Too deep a pocket for the chord to be taken at all.
        assert!(chords(0, &boundary, Sign::Negative, 0.001, &index).is_empty());
        // The set on the other side, inside: a chord across the peaks.
        let inside = chords(0, &boundary, Sign::Positive, 0.05, &index);
        assert_eq!(inside.len(), 1, "{inside:?}");
        assert_eq!((inside[0].from, inside[0].span), (1, 2));
    }

    #[test]
    fn a_chord_is_refused_when_a_vertex_lies_on_the_set_side() {
        assert!(proves(&[room(&[(4.0, 3.0), (0.0, 3.0)])], 4));
        let mut dip = room(&[(4.0, 3.0), (0.0, 3.0)]);
        dip[2].1 = -0.01;
        assert!(!proves(&[dip], 4));
    }

    #[test]
    fn a_chord_is_refused_where_the_boundary_turns_towards_the_set_on_it() {
        // (2, 0) lies on the chord, and the boundary bends towards the set
        // there: the region's side holds both directions along the chord.
        let ring = vec![
            (0.0, 0.0),
            (1.5, 1.0),
            (2.0, 0.0),
            (1.0, 1.0),
            (3.0, 1.0),
            (4.0, 0.0),
            (4.0, 3.0),
            (0.0, 3.0),
        ];
        assert!(!proves(&[ring], 5));
    }

    #[test]
    fn a_chord_is_refused_when_it_runs_on_past_its_end() {
        // (5, 0) lies on the chord's line beyond its end.
        let ring = vec![
            (0.0, 0.0),
            (1.0, 1.0),
            (5.0, 0.0),
            (4.5, 1.0),
            (4.0, 0.0),
            (4.0, -1.0),
            (0.0, -1.0),
        ];
        assert!(!proves(&[ring], 4));
    }

    #[test]
    fn a_chord_is_refused_when_an_end_leaves_outside_the_set() {
        // The edge into the start runs just above the chord.
        assert!(!proves(&[room(&[(4.0, 3.0), (10.0, 0.05)])], 4));
        // The edge out of the end does.
        assert!(!proves(&[room(&[(-6.0, 0.05), (0.0, 3.0)])], 4));
    }

    #[test]
    fn a_chord_is_refused_at_a_spike() {
        // The edge into the start doubles back along the first spanned one.
        assert!(!proves(&[room(&[(4.0, 3.0), (2.0, 0.02)])], 4));
    }

    #[test]
    fn a_chord_is_refused_when_another_boundary_crosses_or_touches_it() {
        let wall = room(&[(4.0, 3.0), (0.0, 3.0)]);
        let crossing = vec![(1.5, -0.5), (2.5, -0.5), (2.0, 0.005)];
        assert!(!proves(&[wall.clone(), crossing], 4));
        let touching = vec![(1.5, -0.5), (2.5, -0.5), (1.5, 0.0)];
        assert!(!proves(&[wall.clone(), touching], 4));
        let apart = vec![(1.5, -0.5), (2.5, -0.5), (1.5, -0.1)];
        assert!(proves(&[wall, apart], 4));
        // Its own ring, past the run, counts too.
        assert!(!proves(
            &[room(&[(4.0, 3.0), (2.5, 3.0), (2.5, -1.0), (0.0, -1.0)])],
            4
        ));
    }

    #[test]
    fn trimming_keeps_what_the_cover_does_not_hold() {
        let hull = points(&[
            (0.0, 0.0),
            (1.0, -0.5),
            (2.0, 0.0),
            (2.0, 2.0),
            (1.0, 2.5),
            (0.0, 2.0),
        ]);
        let cover = points(&[(-1.0, -1.0), (3.0, -1.0), (3.0, 1.0), (-1.0, 1.0)]);
        let same = |p: Vec<Point2>| Some(p);
        let tilted = (Point2::new(0.0, 0.0), Point2::new(1.0, 1.0));
        let kept = trim(hull.clone(), &cover, tilted, same).unwrap();
        assert_eq!(
            kept,
            points(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (1.0, 2.5), (0.0, 2.0)])
        );
        // Every point of the hull stays in the trimmed ring or the cover:
        // each edge's midpoint, and the cut-off vertex.
        let held = |p: &[Point2], q: Point2| {
            (0..p.len()).all(|i| orient(p[i], p[(i + 1) % p.len()], q) != Sign::Negative)
        };
        for i in 0..hull.len() {
            let (a, b) = (hull[i], hull[(i + 1) % hull.len()]);
            let mid = Point2::new(0.5 * a.x + 0.5 * b.x, 0.5 * a.y + 0.5 * b.y);
            assert!(held(&kept, mid) || held(&cover, mid), "{mid:?}");
        }
        assert!(trim(
            hull.clone(),
            &points(&[(-1.0, -1.0), (3.0, -1.0), (3.0, 3.0), (-1.0, 3.0)]),
            tilted,
            same
        )
        .is_none());
        assert_eq!(
            trim(
                hull.clone(),
                &points(&[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0)]),
                tilted,
                same
            ),
            Some(hull.clone())
        );
        // A ring tidying would change is kept whole.
        assert_eq!(
            trim(hull.clone(), &cover, tilted, |_| None),
            Some(hull.clone())
        );
        assert_eq!(
            trim(hull.clone(), &cover, tilted, |p| Some(p[1..].to_vec())),
            Some(hull)
        );
    }

    #[test]
    fn a_cut_parallel_to_the_chord_is_tilted() {
        // Cut from (0, 0) to (2, 0) along a horizontal chord: one more
        // covered vertex is kept, so that the cuts of a wall's teeth do not
        // share a line.
        let hull = points(&[
            (0.0, 0.0),
            (0.7, -0.4),
            (1.3, -0.4),
            (2.0, 0.0),
            (2.0, 2.0),
            (0.0, 2.0),
        ]);
        let cover = points(&[(-1.0, -1.0), (3.0, -1.0), (3.0, 1.0), (-1.0, 1.0)]);
        let chord = (Point2::new(-5.0, -1.0), Point2::new(5.0, -1.0));
        let kept = trim(hull, &cover, chord, Some).unwrap();
        assert_eq!(
            kept,
            points(&[(0.0, 0.0), (0.7, -0.4), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)])
        );
    }
}
