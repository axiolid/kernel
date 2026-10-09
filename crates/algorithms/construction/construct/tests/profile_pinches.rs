//! Rings that touch at single points triangulate on 2D and surface paths
//! and are refused for solids (#262).
//!
//! A plan footprint that is the union of two shadows touching at a corner,
//! a hole whose corner sits on the outer ring or on another hole: each
//! bounds a valid region. `triangulate_with(.., PinchPolicy::Accept)` must
//! tile it exactly -- strictly counter-clockwise triangles, exact area, no
//! vertex inside a triangle or on one of its edges, every triangle edge
//! either twinned or running along a ring edge in its direction and
//! covering it exactly -- while `triangulate` and `extrude_profile`, which
//! build solids, refuse the same rings by name: the extrusion would share
//! one wall edge between four faces.
//!
//! Coordinates are dyadic, so areas are exact in `f64` and compared with
//! `==`.

use std::collections::HashMap;

use axiolid_construct::extrude::extrude_profile;
use axiolid_construct::profile::{ring_touches, triangulate, triangulate_with, PinchPolicy, Rings};
use axiolid_contracts::GeomError;
use axiolid_core::{Point2, Tolerance, Vec3};
use axiolid_guarantees::Sign;
use axiolid_predicates::orient2d;

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Point2> {
    vec![p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)]
}

fn reversed(mut ring: Vec<Point2>) -> Vec<Point2> {
    ring.reverse();
    ring
}

fn twice_area(a: Point2, b: Point2, c: Point2) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)
}

fn cross(a: Point2, b: Point2, c: Point2) -> f64 {
    twice_area(a, b, c)
}

/// Whether `c` lies on the closed segment `a b` (exact for dyadic input).
fn on_segment(a: Point2, b: Point2, c: Point2) -> bool {
    cross(a, b, c) == 0.0
        && a.x.min(b.x) <= c.x
        && c.x <= a.x.max(b.x)
        && a.y.min(b.y) <= c.y
        && c.y <= a.y.max(b.y)
}

/// Every check #262 asks of an accepted triangulation; returns twice the
/// covered area. `exact` compares coordinates exactly; off-grid input
/// skips the parts that need exact products.
fn assert_tiles(rings: &Rings, what: &str, exact: bool) -> (f64, usize) {
    let (points, triangles) = triangulate_with(rings, PinchPolicy::Accept)
        .unwrap_or_else(|e| panic!("{what}: refused: {e:?}"));
    let count: usize = rings.outer.len() + rings.holes.iter().map(Vec::len).sum::<usize>();
    assert_eq!(points.len(), count, "{what}: vertex layout changed");
    let first = |v: u32| {
        points
            .iter()
            .position(|&q| q == points[v as usize])
            .unwrap() as u32
    };

    let mut covered = 0.0;
    for t in &triangles {
        for &v in t {
            assert_eq!(
                first(v),
                v,
                "{what}: {t:?} does not use a point's first index"
            );
        }
        let (a, b, c) = (
            points[t[0] as usize],
            points[t[1] as usize],
            points[t[2] as usize],
        );
        let doubled = twice_area(a, b, c);
        assert_eq!(
            orient2d(a, b, c).sign(),
            Some(Sign::Positive),
            "{what}: {t:?} is not strictly counter-clockwise"
        );
        covered += doubled;
        if exact {
            for (i, &q) in points.iter().enumerate() {
                if t.contains(&(i as u32)) || q == a || q == b || q == c {
                    continue;
                }
                let inside =
                    cross(a, b, q) >= 0.0 && cross(b, c, q) >= 0.0 && cross(c, a, q) >= 0.0;
                assert!(!inside, "{what}: vertex {i} lies in or on triangle {t:?}");
            }
        }
    }

    // Edge parity in coordinates: an edge not twinned runs along a ring
    // edge in its own direction, and those pieces cover each ring edge once.
    let key = |q: Point2| ((q.x + 0.0).to_bits(), (q.y + 0.0).to_bits());
    let mut directed: HashMap<(Key, Key), (Point2, Point2, u32)> = HashMap::new();
    for t in &triangles {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let (pa, pb) = (points[a as usize], points[b as usize]);
            directed.entry((key(pa), key(pb))).or_insert((pa, pb, 0)).2 += 1;
        }
    }
    let ring_edges: Vec<(Point2, Point2)> = core::iter::once(&rings.outer)
        .chain(&rings.holes)
        .flat_map(|ring| (0..ring.len()).map(move |k| (ring[k], ring[(k + 1) % ring.len()])))
        .collect();
    let mut along: Vec<f64> = vec![0.0; ring_edges.len()];
    for (&(ka, kb), &(a, b, uses)) in &directed {
        assert_eq!(uses, 1, "{what}: edge {a:?}->{b:?} used {uses} times");
        if directed.contains_key(&(kb, ka)) {
            continue;
        }
        if !exact {
            continue;
        }
        // A boundary edge: on some ring edge, same direction (either
        // ring orientation is allowed, so either direction of the ring
        // edge, but all pieces of one ring edge must agree).
        let host = ring_edges
            .iter()
            .position(|&(u, v)| on_segment(u, v, a) && on_segment(u, v, b));
        let Some(host) = host else {
            panic!("{what}: edge {a:?}->{b:?} is neither twinned nor on a ring edge");
        };
        let (u, v) = ring_edges[host];
        let length = (b.x - a.x).abs() + (b.y - a.y).abs();
        let same = (b.x - a.x) * (v.x - u.x) + (b.y - a.y) * (v.y - u.y) > 0.0;
        along[host] += if same { length } else { -length };
    }
    if exact {
        for (k, &(u, v)) in ring_edges.iter().enumerate() {
            let length = (v.x - u.x).abs() + (v.y - u.y).abs();
            assert_eq!(
                along[k].abs(),
                length,
                "{what}: ring edge {u:?}->{v:?} is not covered exactly once from one side"
            );
        }
    }
    (covered, triangles.len())
}

/// A point's exact coordinates, as a hash key.
type Key = (u64, u64);

fn region(outer: Vec<Point2>, holes: Vec<Vec<Point2>>) -> Rings {
    Rings { outer, holes }
}

/// `rings` with the `k`-th zero coordinate (in ring order, `x` before `y`)
/// written as `-0.0` where bit `k % 64` of `pattern` is set: the same
/// points, as projecting onto plane axes can write them (#269).
fn signed_zeros(rings: &Rings, pattern: u64) -> Rings {
    let mut k = 0;
    let mut flip = |v: f64| {
        if v != 0.0 {
            return v;
        }
        let negate = pattern >> (k % 64) & 1 == 1;
        k += 1;
        if negate {
            -0.0
        } else {
            0.0
        }
    };
    let mut ring = |r: &Vec<Point2>| -> Vec<Point2> {
        r.iter()
            .map(|q| {
                let x = flip(q.x);
                p(x, flip(q.y))
            })
            .collect()
    };
    let outer = ring(&rings.outer);
    let holes = rings.holes.iter().map(&mut ring).collect();
    region(outer, holes)
}

/// The zero-sign patterns every fixture is checked under: as written, every
/// zero negative, and alternating in runs of one, two and four zeros, so
/// any two of a fixture's first eight zeros differ in sign under one of
/// them.
const SIGN_PATTERNS: [u64; 8] = [
    0,
    u64::MAX,
    0x5555_5555_5555_5555,
    0xAAAA_AAAA_AAAA_AAAA,
    0x3333_3333_3333_3333,
    0xCCCC_CCCC_CCCC_CCCC,
    0x0F0F_0F0F_0F0F_0F0F,
    0xF0F0_F0F0_F0F0_F0F0,
];

/// Accepted with `triangles` triangles of twice-area `twice`, from every
/// start vertex of the outer ring, either way round and under every
/// [`SIGN_PATTERNS`] of its zeros; refused by `triangulate` and the
/// extrusion with a message containing `refusal`.
fn check(rings: Rings, triangles: usize, twice: f64, refusal: &str, what: &str) {
    for pattern in SIGN_PATTERNS {
        let signed = signed_zeros(&rings, pattern);
        check_signed(
            signed,
            triangles,
            twice,
            refusal,
            &format!("{what}, zeros {pattern:x}"),
        );
    }
}

fn check_signed(rings: Rings, triangles: usize, twice: f64, refusal: &str, what: &str) {
    let n = rings.outer.len();
    for start in 0..n {
        for flip in [false, true] {
            let mut outer: Vec<Point2> = (0..n).map(|k| rings.outer[(start + k) % n]).collect();
            if flip {
                outer.reverse();
            }
            let holes = rings
                .holes
                .iter()
                .map(|h| if flip { reversed(h.clone()) } else { h.clone() })
                .collect();
            let turned = region(outer, holes);
            let what = format!("{what} (start {start}, flipped {flip})");
            let (covered, count) = assert_tiles(&turned, &what, true);
            assert_eq!(covered, twice, "{what}: area");
            assert_eq!(count, triangles, "{what}: triangle count");
        }
    }
    match triangulate(&rings) {
        Err(GeomError::InvalidInput(message)) => {
            assert!(message.contains(refusal), "{what}: {message}");
        }
        other => panic!("{what}: a solid cap must refuse touching rings, got {other:?}"),
    }
    match extrude_profile(&rings, Vec3::Z, 1.0, Tolerance::METRE) {
        Err(GeomError::InvalidInput(message)) => {
            assert!(message.contains(refusal), "{what}: {message}");
        }
        other => panic!("{what}: the extrusion must be refused, got {other:?}"),
    }
}

#[test]
fn two_footprints_touching_at_a_corner_triangulate() {
    // The union of two unit squares diagonal to each other, one ring
    // through the shared corner twice: n = 8, two parts, 8 - 4 = 4.
    let ring = vec![
        p(0.0, 0.0),
        p(1.0, 0.0),
        p(1.0, 1.0),
        p(2.0, 1.0),
        p(2.0, 2.0),
        p(1.0, 2.0),
        p(1.0, 1.0),
        p(0.0, 1.0),
    ];
    check(
        region(ring, vec![]),
        4,
        4.0,
        "outer ring intersects itself",
        "diagonal squares",
    );
    // Round the origin, the second visit written as -0.0: the same point.
    let ring = vec![
        p(-1.0, -1.0),
        p(0.0, -1.0),
        p(0.0, 0.0),
        p(1.0, 0.0),
        p(1.0, 1.0),
        p(0.0, 1.0),
        p(-0.0, 0.0),
        p(-1.0, 0.0),
    ];
    check(
        region(ring, vec![]),
        4,
        4.0,
        "outer ring intersects itself",
        "diagonal squares at a signed zero",
    );
}

#[test]
fn an_l_shaped_room_and_a_column_touching_at_a_corner_triangulate() {
    // An L (3 units of area) and a unit square meeting at the L's corner
    // (1, 2): n = 10, two parts, 10 - 4 = 6.
    let ring = vec![
        p(0.0, 0.0),
        p(2.0, 0.0),
        p(2.0, 1.0),
        p(1.0, 1.0),
        p(1.0, 2.0),
        p(2.0, 2.0),
        p(2.0, 3.0),
        p(1.0, 3.0),
        p(1.0, 2.0),
        p(0.0, 2.0),
    ];
    check(
        region(ring, vec![]),
        6,
        8.0,
        "outer ring intersects itself",
        "L and column",
    );
    // The same with a diamond column whose corner sits at the middle of
    // the L's top edge: one ring through (0.5, 2) twice.
    let ring = vec![
        p(0.0, 0.0),
        p(2.0, 0.0),
        p(2.0, 1.0),
        p(1.0, 1.0),
        p(1.0, 2.0),
        p(0.5, 2.0),
        p(1.0, 2.5),
        p(0.5, 3.0),
        p(0.0, 2.5),
        p(0.5, 2.0),
        p(0.0, 2.0),
    ];
    check(
        region(ring, vec![]),
        7,
        7.0,
        "outer ring intersects itself",
        "L and diamond",
    );
}

#[test]
fn a_hole_touching_the_outer_ring_at_a_corner_triangulates() {
    // n = 7 in one part with no hole left: 5 triangles.
    let hole = vec![p(0.0, 0.0), p(2.0, 1.0), p(1.0, 2.0)];
    check(
        region(rect(0.0, 0.0, 4.0, 4.0), vec![hole]),
        5,
        29.0,
        "hole 0 touches or crosses the outer ring",
        "hole at the outer corner",
    );
}

#[test]
fn a_hole_touching_the_outer_ring_inside_an_edge_triangulates() {
    // The apex (2, 0) splits the bottom edge: n = 8, one part, 6.
    let hole = vec![p(2.0, 0.0), p(3.0, 1.0), p(1.0, 1.0)];
    check(
        region(rect(0.0, 0.0, 4.0, 4.0), vec![hole]),
        6,
        30.0,
        "hole 0 touches or crosses the outer ring",
        "hole on the outer edge",
    );
    // A diamond touching the bottom and the top edge splits the region in
    // two: n = 10, two parts, 6.
    let hole = vec![p(2.0, 0.0), p(3.0, 2.0), p(2.0, 4.0), p(1.0, 2.0)];
    check(
        region(rect(0.0, 0.0, 4.0, 4.0), vec![hole]),
        6,
        24.0,
        "hole 0 touches or crosses the outer ring",
        "hole across the region",
    );
}

#[test]
fn two_holes_touching_at_a_corner_triangulate() {
    // n = 12; the holes join into one hole cycle: 12 + 2 - 2.
    check(
        region(
            rect(0.0, 0.0, 4.0, 4.0),
            vec![rect(1.0, 1.0, 2.0, 2.0), rect(2.0, 2.0, 3.0, 3.0)],
        ),
        12,
        28.0,
        "holes 0 and 1 overlap or touch",
        "holes at a corner",
    );
    // A triangle's apex on the middle of the other hole's side.
    check(
        region(
            rect(0.0, 0.0, 4.0, 4.0),
            vec![
                rect(1.0, 1.0, 2.0, 3.0),
                vec![p(2.0, 2.0), p(3.0, 1.5), p(3.0, 2.5)],
            ],
        ),
        12,
        27.0,
        "holes 0 and 1 overlap or touch",
        "hole on a hole's side",
    );
}

#[test]
fn a_hole_with_every_corner_on_the_outer_ring_triangulates() {
    // Each corner of the hole splits a wall, and every vertex of the hole
    // is then a vertex of the outer ring too: its side is read from the
    // sector it leaves the shared corner into. Three corner parts and the
    // top: 10 loop vertices, 3 parts, 10 - 6 = 4 triangles.
    let hole = vec![p(2.0, 0.0), p(4.0, 2.0), p(0.0, 2.0)];
    check(
        region(rect(0.0, 0.0, 4.0, 4.0), vec![hole]),
        4,
        24.0,
        "hole 0 touches or crosses the outer ring",
        "hole with every corner on a wall",
    );
}

#[test]
fn two_holes_sharing_their_extreme_corner_triangulate() {
    // Both holes leave their leftmost (then their rightmost) corner: the
    // joined hole cycle visits it twice, once round the wedge between
    // them and once round the outside, so it must not read as an outer
    // cycle there, and its bridge must leave from the outside visit.
    let left = vec![
        vec![p(1.0, 2.0), p(2.0, 1.0), p(3.0, 1.0)],
        vec![p(1.0, 2.0), p(3.0, 3.0), p(2.0, 3.0)],
    ];
    let right: Vec<Vec<Point2>> = left
        .iter()
        .map(|h| reversed(h.iter().map(|q| p(4.0 - q.x, q.y)).collect()))
        .collect();
    for (holes, what) in [(left, "leftmost"), (right, "rightmost")] {
        check(
            region(rect(0.0, 0.0, 4.0, 4.0), holes),
            10,
            30.0,
            "holes 0 and 1 overlap or touch",
            what,
        );
    }
}

#[test]
fn holes_touching_in_a_ring_close_off_a_part_with_its_own_hole() {
    // Four holes round [2, 4]^2, each touching the next at a corner of
    // it, make that square a part of its own, bounded by hole edges; the
    // fifth hole inside belongs to it, not to the outer ring.
    let holes = vec![
        rect(2.0, 1.0, 4.0, 2.0),
        rect(4.0, 2.0, 5.0, 4.0),
        rect(2.0, 4.0, 4.0, 5.0),
        rect(1.0, 2.0, 2.0, 4.0),
        rect(2.5, 2.5, 3.5, 3.5),
    ];
    // 24 loop vertices, 2 outer cycles, 2 hole cycles: 24 triangles over
    // 36 - 8 - 1 = 27 units.
    check(
        region(rect(0.0, 0.0, 6.0, 6.0), holes),
        24,
        54.0,
        "overlap or touch",
        "ring of holes",
    );
}

#[test]
fn an_outer_ring_pinched_around_a_hole_triangulates() {
    // The hole drawn as a loop of the outer ring itself, through (0, 0).
    let ring = vec![
        p(0.0, 0.0),
        p(1.0, 2.0),
        p(2.0, 1.0),
        p(0.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 4.0),
        p(0.0, 4.0),
    ];
    check(
        region(ring, vec![]),
        5,
        29.0,
        "outer ring intersects itself",
        "pinched loop",
    );
}

// The consumer's call sites (axioval, #262): plan regions re-fed from
// overlay output -- unallocated space, cap coverage, free floor, span
// sections, envelope coverage, room minus door zones, walkable plan
// join / subtract / intersect -- are where pinches come from.

#[test]
fn a_union_of_two_rooms_meeting_at_one_corner_triangulates() {
    // Two 2 x 1 rooms, the second up and to the right, as the union
    // reports them: one ring through the shared corner (2, 1) twice.
    let ring = vec![
        p(0.0, 0.0),
        p(2.0, 0.0),
        p(2.0, 1.0),
        p(4.0, 1.0),
        p(4.0, 2.0),
        p(2.0, 2.0),
        p(2.0, 1.0),
        p(0.0, 1.0),
    ];
    check(
        region(ring, vec![]),
        4,
        8.0,
        "outer ring intersects itself",
        "rooms at a corner",
    );
    // The same with a straight corner on the left wall: under the zero
    // sign patterns `(-0.0, 0.5)` comes before the corner `(0, 0)` in
    // `total_cmp` order, and the left room's cycle read at it turns
    // neither way (#269). n = 9, two parts: 5 triangles.
    let ring = vec![
        p(0.0, 0.0),
        p(2.0, 0.0),
        p(2.0, 1.0),
        p(4.0, 1.0),
        p(4.0, 2.0),
        p(2.0, 2.0),
        p(2.0, 1.0),
        p(0.0, 1.0),
        p(0.0, 0.5),
    ];
    check(
        region(ring, vec![]),
        5,
        8.0,
        "outer ring intersects itself",
        "rooms at a corner, straight left wall",
    );
}

#[test]
fn an_l_shaped_room_wrapped_round_a_column_corner_triangulates() {
    // A room [0, 3] x [0, 3] minus a column [2, 3] x [2, 3] at its corner
    // and minus [1, 2] x [1, 2], the column's diagonal neighbour cut away
    // too: the ring touches itself at the column's corner (2, 2).
    let ring = vec![
        p(0.0, 0.0),
        p(3.0, 0.0),
        p(3.0, 2.0),
        p(2.0, 2.0),
        p(2.0, 1.0),
        p(1.0, 1.0),
        p(1.0, 2.0),
        p(2.0, 2.0),
        p(2.0, 3.0),
        p(0.0, 3.0),
    ];
    // n = 10 in one part with the pinched-off loop as a hole cycle:
    // the loop [1, 2]^2 bounds a hole touching the rest at (2, 2), so
    // 10 + 0 - 2 = 8 triangles over 9 - 1 - 1 = 7 units.
    check(
        region(ring, vec![]),
        8,
        14.0,
        "outer ring intersects itself",
        "L round a column",
    );
}

#[test]
fn a_corridor_eroded_to_a_point_triangulates() {
    // A 6 x 2 corridor narrowed to nothing at (3, 1): two lobes of five
    // units each, joined at one vertex.
    let ring = vec![
        p(0.0, 0.0),
        p(2.0, 0.0),
        p(3.0, 1.0),
        p(4.0, 0.0),
        p(6.0, 0.0),
        p(6.0, 2.0),
        p(4.0, 2.0),
        p(3.0, 1.0),
        p(2.0, 2.0),
        p(0.0, 2.0),
    ];
    check(
        region(ring, vec![]),
        6,
        20.0,
        "outer ring intersects itself",
        "eroded corridor",
    );
}

#[test]
fn an_even_odd_pinch_triangulates() {
    // The same corridor as an even-odd fill reports it: the ring runs
    // straight through (3, 1) twice, crossing itself there, the right lobe
    // clockwise. Inside an odd number of times is the same two lobes.
    let ring = vec![
        p(0.0, 0.0),
        p(2.0, 0.0),
        p(3.0, 1.0),
        p(4.0, 2.0),
        p(6.0, 2.0),
        p(6.0, 0.0),
        p(4.0, 0.0),
        p(3.0, 1.0),
        p(2.0, 2.0),
        p(0.0, 2.0),
    ];
    check(
        region(ring, vec![]),
        6,
        20.0,
        "outer ring intersects itself",
        "even-odd corridor",
    );
    // Two rooms at a corner the same way: a figure eight through (2, 2).
    let eight = vec![
        p(0.0, 0.0),
        p(2.0, 2.0),
        p(4.0, 4.0),
        p(4.0, 0.0),
        p(2.0, 2.0),
        p(0.0, 4.0),
    ];
    check(
        region(eight, vec![]),
        2,
        16.0,
        "outer ring intersects itself",
        "even-odd rooms",
    );
}

#[test]
fn a_room_minus_door_zones_touching_its_walls_triangulates() {
    // A 4 x 3 room minus two door swing zones: one a triangle whose corner
    // is the room's corner, one whose corners sit inside the right and the
    // top wall, touching the first at (2, 1).
    let room = rect(0.0, 0.0, 4.0, 3.0);
    let zones = vec![
        vec![p(0.0, 0.0), p(2.0, 1.0), p(1.0, 2.0)],
        vec![p(2.0, 1.0), p(4.0, 1.0), p(3.0, 3.0)],
    ];
    // 12 loop vertices (the room's walls split at (4, 1) and (3, 3)); the
    // zones cut the room into three parts (below them, above them, and the
    // corner beyond the second): 12 - 6 = 6 triangles over
    // 12 - 1.5 - 2 = 8.5 units.
    check(
        region(room, zones),
        6,
        17.0,
        "touch",
        "room minus door zones",
    );
}

#[test]
fn rings_crossing_each_other_at_shared_vertices_are_refused_by_name() {
    // Two holes crossing only through two shared corners: the second
    // leaves the first's square at (1, 1), runs round outside it and
    // comes back in at (3, 3). The square starts at (1, 3), outside the
    // second, so neither reads as inside the other and the crossing is
    // found where the spokes round (1, 1) do not alternate.
    let holes = vec![
        vec![p(1.0, 3.0), p(1.0, 1.0), p(3.0, 1.0), p(3.0, 3.0)],
        vec![
            p(1.0, 1.0),
            p(3.5, 0.5),
            p(3.5, 3.5),
            p(3.0, 3.0),
            p(2.0, 1.5),
        ],
    ];
    match triangulate_with(
        &region(rect(0.0, 0.0, 4.0, 4.0), holes),
        PinchPolicy::Accept,
    ) {
        Err(GeomError::InvalidInput(message)) => {
            assert!(
                message.contains("cross at their shared vertex"),
                "{message}"
            );
        }
        other => panic!("expected a refusal, got {other:?}"),
    }
    // Two holes whose corners meet at (2, 2) with their sectors
    // interleaved: each runs through the other's corner.
    let holes = vec![
        vec![p(1.0, 1.0), p(3.0, 3.0), p(1.0, 2.0)],
        vec![p(2.0, 1.0), p(2.0, 3.0), p(3.0, 2.0)],
    ];
    match triangulate_with(
        &region(rect(0.0, 0.0, 4.0, 4.0), holes),
        PinchPolicy::Accept,
    ) {
        Err(GeomError::InvalidInput(message)) => {
            assert!(message.contains("cross"), "{message}");
        }
        other => panic!("expected a refusal, got {other:?}"),
    }
}

#[test]
fn rings_that_overlap_or_cross_stay_refused_when_pinches_are_accepted() {
    let refused = |rings: Rings| match triangulate_with(&rings, PinchPolicy::Accept) {
        Err(GeomError::InvalidInput(message)) => message,
        other => panic!("expected a refusal, got {other:?}"),
    };
    let square = || rect(0.0, 0.0, 4.0, 4.0);
    // Sharing a stretch of edge.
    let message = refused(region(
        square(),
        vec![rect(1.0, 1.0, 2.0, 2.0), rect(2.0, 1.5, 3.0, 2.5)],
    ));
    assert!(message.contains("holes 0 and 1 overlap"), "{message}");
    // Sharing a stretch of edge from a common corner.
    let message = refused(region(
        square(),
        vec![
            rect(1.0, 1.0, 2.0, 2.0),
            vec![p(1.0, 1.0), p(2.0, 0.5), p(3.0, 1.0)],
        ],
    ));
    assert!(message.contains("holes 0 and 1 overlap"), "{message}");
    // A hole along the outer ring's edge.
    let message = refused(region(square(), vec![rect(0.0, 1.0, 1.0, 2.0)]));
    assert!(
        message.contains("hole 0 overlaps the outer ring"),
        "{message}"
    );
    // A proper crossing.
    let message = refused(region(square(), vec![rect(3.0, 1.0, 5.0, 2.0)]));
    assert!(
        message.contains("hole 0 crosses the outer ring"),
        "{message}"
    );
    let message = refused(region(
        square(),
        vec![rect(1.0, 1.0, 2.5, 2.5), rect(2.0, 2.0, 3.0, 3.0)],
    ));
    assert!(message.contains("holes 0 and 1 cross"), "{message}");
    // A hole outside, touching the outer ring at a corner.
    let message = refused(region(square(), vec![rect(4.0, 4.0, 5.0, 5.0)]));
    assert!(
        message.contains("hole 0 lies outside the outer ring"),
        "{message}"
    );
    // A hole inside another, touching it at a corner.
    let message = refused(region(
        square(),
        vec![rect(1.0, 1.0, 3.0, 3.0), rect(1.0, 1.0, 2.0, 2.0)],
    ));
    assert!(message.contains("overlap"), "{message}");
    let message = refused(region(
        square(),
        vec![
            rect(1.0, 1.0, 3.0, 3.0),
            vec![p(1.0, 1.0), p(2.0, 1.5), p(1.5, 2.0)],
        ],
    ));
    assert!(message.contains("hole 1 lies inside hole 0"), "{message}");
    // A self-crossing outer ring away from any vertex.
    let bowtie = vec![p(0.0, 0.0), p(2.0, 2.0), p(2.0, 0.0), p(0.0, 2.0)];
    let message = refused(region(bowtie, vec![]));
    assert!(
        message.contains("outer ring intersects itself"),
        "{message}"
    );
}

/// A deterministic LCG, as in `profile_holes.rs`.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// A random region full of pinches on an `n x n` grid: square holes on
/// a checkerboard (touching each other at corners), diamonds in the
/// bottom row (touching each other, the outer ring inside its bottom edge,
/// and squares above them inside their bottom edges), rings started
/// anywhere and wound either way.
fn pinched_layout(rng: &mut Lcg, scale: f64, shift: Point2) -> (Rings, f64) {
    let n = 4 + rng.below(6) as i64;
    let at = |x: f64, y: f64| p(x * scale + shift.x, y * scale + shift.y);
    let mut holes: Vec<Vec<Point2>> = Vec::new();
    let mut twice = 2.0 * (n * n) as f64;
    for i in 1..n - 1 {
        for j in 1..n - 1 {
            if (i + j) % 2 == 0 && rng.below(3) != 0 {
                let (x, y) = (i as f64, j as f64);
                holes.push(vec![
                    at(x, y),
                    at(x + 1.0, y),
                    at(x + 1.0, y + 1.0),
                    at(x, y + 1.0),
                ]);
                twice -= 2.0;
            }
        }
    }
    for i in 0..n {
        if rng.below(2) == 0 {
            let x = i as f64;
            holes.push(vec![
                at(x + 0.5, 0.0),
                at(x + 1.0, 0.5),
                at(x + 0.5, 1.0),
                at(x, 0.5),
            ]);
            twice -= 1.0;
        }
    }
    let size = n as f64;
    let mut outer = vec![at(0.0, 0.0), at(size, 0.0), at(size, size), at(0.0, size)];
    if rng.below(2) == 0 {
        outer.reverse();
    }
    for hole in &mut holes {
        let k = rng.below(hole.len() as u64) as usize;
        hole.rotate_left(k);
        if rng.below(2) == 0 {
            hole.reverse();
        }
    }
    // Shuffle the hole order.
    for k in (1..holes.len()).rev() {
        let j = rng.below(k as u64 + 1) as usize;
        holes.swap(k, j);
    }
    (region(outer, holes), twice * scale * scale)
}

#[test]
fn random_pinched_layouts_tile_exactly() {
    let mut rng = Lcg(0x262);
    let mut pinched = 0;
    for round in 0..300 {
        let (rings, twice) = pinched_layout(&mut rng, 1.0, p(0.0, 0.0));
        let what = format!("round {round}");
        let (covered, _) = assert_tiles(&rings, &what, true);
        assert_eq!(covered, twice, "{what}: area");
        // A layout where nothing happens to touch is a solid cap too, and
        // then triangulates identically; any touch is refused for a solid.
        match triangulate(&rings) {
            Ok(solid) => assert_eq!(
                Ok(solid),
                triangulate_with(&rings, PinchPolicy::Accept),
                "{what}"
            ),
            Err(GeomError::InvalidInput(_)) => pinched += 1,
            Err(other) => panic!("{what}: a solid cap refused with {other:?}"),
        }
    }
    assert!(
        pinched >= 250,
        "only {pinched} of 300 layouts touch anywhere"
    );
}

#[test]
fn random_pinched_layouts_with_signed_zeros_tile_exactly() {
    // As above, with zero coordinates written as -0.0 at random (#269);
    // shifted by one cell too, so the first column of square holes and
    // the diamonds' corners sit on the axes.
    let mut rng = Lcg(0x269);
    for round in 0..300 {
        let shift = if round % 2 == 0 {
            p(0.0, 0.0)
        } else {
            p(-1.0, -1.0)
        };
        let (rings, twice) = pinched_layout(&mut rng, 1.0, shift);
        let pattern = rng.next() << 31 ^ rng.next();
        let rings = signed_zeros(&rings, pattern);
        let what = format!("signed-zero round {round}");
        let (covered, _) = assert_tiles(&rings, &what, true);
        assert_eq!(covered, twice, "{what}: area");
        if let Ok(solid) = triangulate(&rings) {
            assert_eq!(
                Ok(solid),
                triangulate_with(&rings, PinchPolicy::Accept),
                "{what}"
            );
        }
    }
}

#[test]
fn a_straight_corner_on_negative_zero_does_not_flip_the_ring() {
    // #269: a 1 x 3 rectangle, counter-clockwise, with a straight corner on
    // its left side, which is -0.0 as projecting onto plane axes writes it.
    // `total_cmp` read the orientation at (-0.0, 2), the straight corner,
    // and reversed the ring.
    let rings = region(
        vec![
            p(0.0, 0.0),
            p(1.0, 0.0),
            p(1.0, 3.0),
            p(-0.0, 3.0),
            p(-0.0, 2.0),
        ],
        vec![],
    );
    for policy in [PinchPolicy::Refuse, PinchPolicy::Accept] {
        let (_, triangles) = triangulate_with(&rings, policy).expect("a rectangle triangulates");
        assert_eq!(triangles.len(), 3);
    }
    let (covered, count) = assert_tiles(&rings, "rectangle at -0.0", true);
    assert_eq!((covered, count), (6.0, 3));
    // Every start vertex, either way round, under every zero sign.
    check_accepted(&rings, 3, 6.0, "rectangle at -0.0");
}

#[test]
fn a_hole_whose_lowest_corner_has_a_negative_zero_twin_triangulates() {
    // The same ring as a hole of a 4 x 6 square: a misread hole is
    // bridged the wrong way round and finds no ear either.
    let hole = vec![
        p(0.0, 0.0),
        p(1.0, 0.0),
        p(1.0, 3.0),
        p(-0.0, 3.0),
        p(-0.0, 2.0),
    ];
    for hole in [hole.clone(), reversed(hole)] {
        let rings = region(rect(-2.0, -2.0, 2.0, 4.0), vec![hole]);
        // n = 9 with one hole: 9 triangles over 24 - 3 units.
        for policy in [PinchPolicy::Refuse, PinchPolicy::Accept] {
            let (_, triangles) = triangulate_with(&rings, policy).expect("triangulates");
            assert_eq!(triangles.len(), 9);
        }
        check_accepted(&rings, 9, 42.0, "hole at -0.0");
    }
}

/// [`check`]'s acceptance half alone, for rings that touch nowhere.
fn check_accepted(rings: &Rings, triangles: usize, twice: f64, what: &str) {
    for pattern in SIGN_PATTERNS {
        let signed = signed_zeros(rings, pattern);
        let n = signed.outer.len();
        for start in 0..n {
            for flip in [false, true] {
                let mut outer: Vec<Point2> =
                    (0..n).map(|k| signed.outer[(start + k) % n]).collect();
                let mut holes = signed.holes.clone();
                if flip {
                    outer.reverse();
                    holes = holes.into_iter().map(reversed).collect();
                }
                let turned = region(outer, holes);
                let what = format!("{what} (zeros {pattern:x}, start {start}, flipped {flip})");
                let (covered, count) = assert_tiles(&turned, &what, true);
                assert_eq!((covered, count), (twice, triangles), "{what}");
                let solid = triangulate(&turned).unwrap_or_else(|e| panic!("{what}: {e:?}"));
                assert_eq!(solid.1.len(), triangles, "{what}");
            }
        }
    }
}

#[test]
fn random_pinched_layouts_off_the_grid_tile() {
    // Scaled and shifted by non-dyadic amounts: every contact survives
    // (it lies on an axis-parallel line or is a shared vertex), but the
    // coordinates are no longer exact sums.
    let mut rng = Lcg(0x2620);
    for round in 0..200 {
        let (rings, twice) = pinched_layout(&mut rng, 0.3, p(0.1, -7.7));
        let what = format!("off-grid round {round}");
        let (covered, _) = assert_tiles(&rings, &what, false);
        assert!(
            (covered - twice).abs() <= 1e-12 * twice,
            "{what}: covered {covered}, expected {twice}"
        );
    }
}

#[test]
fn rings_touching_nowhere_triangulate_as_for_a_solid() {
    let mut rng = Lcg(0x2621);
    for round in 0..100 {
        let (mut rings, _) = pinched_layout(&mut rng, 1.0, p(0.0, 0.0));
        // Shrink every hole about its centre so nothing touches.
        for hole in &mut rings.holes {
            let c = hole
                .iter()
                .fold(p(0.0, 0.0), |s, q| p(s.x + q.x, s.y + q.y));
            let c = p(c.x / hole.len() as f64, c.y / hole.len() as f64);
            for q in hole.iter_mut() {
                *q = p(c.x + (q.x - c.x) * 0.5, c.y + (q.y - c.y) * 0.5);
            }
        }
        assert_eq!(
            triangulate_with(&rings, PinchPolicy::Accept).expect("accepted"),
            triangulate(&rings).expect("a solid cap too"),
            "round {round}"
        );
    }
}

/// `ring_touches` names exactly the ring edges the accepted triangulation
/// splits (#265): a ring edge is a triangle edge unless a vertex lies
/// inside it, and then it is not. A caller welding a shared edge through
/// that vertex relies on both directions.
#[test]
fn ring_touches_are_the_edges_the_triangulation_splits() {
    let hole = vec![p(2.0, 0.0), p(3.0, 1.0), p(1.0, 1.0)];
    let rings = region(rect(0.0, 0.0, 4.0, 4.0), vec![hole]);
    let touches = ring_touches(&rings).expect("valid rings");
    assert_eq!(touches.len(), 1, "{touches:?}");
    let touch = touches[0];
    assert_eq!((touch.ring, touch.edge, touch.vertex), (0, 0, 4));
    assert!(ring_touches(&region(rect(0.0, 0.0, 4.0, 4.0), vec![]))
        .expect("valid")
        .is_empty());
    let crossing = vec![p(-1.0, 1.0), p(1.0, 1.0), p(1.0, 2.0)];
    assert!(ring_touches(&region(rect(0.0, 0.0, 4.0, 4.0), vec![crossing])).is_err());

    let key = |q: Point2| ((q.x + 0.0).to_bits(), (q.y + 0.0).to_bits());
    let mut rng = Lcg(0x2650);
    for round in 0..200 {
        let (rings, _) = pinched_layout(&mut rng, 1.0, p(0.0, 0.0));
        let (points, triangles) = triangulate_with(&rings, PinchPolicy::Accept).expect("accepted");
        let edges: std::collections::HashSet<(Key, Key)> = triangles
            .iter()
            .flat_map(|t| (0..3).map(move |i| (t[i], t[(i + 1) % 3])))
            .map(|(a, b)| {
                let (a, b) = (key(points[a as usize]), key(points[b as usize]));
                (a.min(b), a.max(b))
            })
            .collect();
        let split: std::collections::HashSet<(usize, usize)> = ring_touches(&rings)
            .expect("accepted rings are valid")
            .iter()
            .map(|touch| (touch.ring, touch.edge))
            .collect();
        for (r, ring) in std::iter::once(&rings.outer)
            .chain(&rings.holes)
            .enumerate()
        {
            for k in 0..ring.len() {
                let (a, b) = (key(ring[k]), key(ring[(k + 1) % ring.len()]));
                assert_eq!(
                    edges.contains(&(a.min(b), a.max(b))),
                    !split.contains(&(r, k)),
                    "round {round}: ring {r} edge {k}"
                );
            }
        }
    }
}
