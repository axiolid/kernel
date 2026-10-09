//! Profiles with several holes triangulate to closed, consistently wound
//! extrusions with exact area and volume (#253).
//!
//! Two holes side by side in one horizontal band used to leave the cap
//! with a T-junction: a triangle edge ran along the band's bottom line
//! straight past the inner corners of both holes, so the extrusion had
//! eight boundary edges. Every case below checks the triangulation the
//! strict way -- each triangle strictly counter-clockwise, every ring edge
//! used exactly once in its own direction, every other edge exactly once in
//! each direction -- and then measures the extrusion, which
//! `volume_properties` refuses unless it is a closed two-manifold.
//!
//! Coordinates are dyadic (multiples of a power of two), so the area sums
//! below are exact in `f64` and compared with `==`.

use std::collections::HashMap;

use axiolid_construct::extrude::{extrude_profile, outward_orientation};
use axiolid_construct::profile::{profile_rings, triangulate, Rings};
use axiolid_contracts::GeomError;
use axiolid_core::{Frame2, Interval, Point2, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::volume_properties;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};

const DEPTH: f64 = 2.0;

/// An axis-aligned rectangle, counter-clockwise.
fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Point2> {
    vec![
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ]
}

/// The same ring, clockwise: the hole convention.
fn cw(mut ring: Vec<Point2>) -> Vec<Point2> {
    ring.reverse();
    ring
}

/// Twice the signed area of a ring, by the shoelace formula.
fn twice_ring_area(ring: &[Point2]) -> f64 {
    (0..ring.len())
        .map(|i| {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            a.x * b.y - b.x * a.y
        })
        .sum()
}

/// Twice the signed area of one triangle.
fn twice_area(a: Point2, b: Point2, c: Point2) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)
}

/// Every check #253 asks of a profile with holes; returns twice the area
/// the triangles cover.
fn assert_sound(rings: &Rings, what: &str) -> f64 {
    let (points, triangles) =
        triangulate(rings).unwrap_or_else(|e| panic!("{what}: triangulation refused: {e:?}"));
    let vertex_count: usize = rings.outer.len() + rings.holes.iter().map(Vec::len).sum::<usize>();
    assert_eq!(points.len(), vertex_count, "{what}: vertex layout changed");
    assert_eq!(
        triangles.len(),
        vertex_count + 2 * rings.holes.len() - 2,
        "{what}: a polygon with n vertices and h holes has n + 2h - 2 triangles"
    );

    // Strictly counter-clockwise triangles, and their exact area.
    let mut covered = 0.0;
    for t in &triangles {
        let (a, b, c) = (
            points[t[0] as usize],
            points[t[1] as usize],
            points[t[2] as usize],
        );
        let doubled = twice_area(a, b, c);
        assert!(
            doubled > 0.0,
            "{what}: triangle {t:?} is not counter-clockwise"
        );
        covered += doubled;
    }

    // Edge parity: every ring edge once in its own direction, never
    // reversed; every other edge once in each direction.
    let mut directed: HashMap<(u32, u32), u32> = HashMap::new();
    for t in &triangles {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *directed.entry((a, b)).or_default() += 1;
        }
    }
    let mut ring_edges = Vec::new();
    let mut start = 0u32;
    for ring in core::iter::once(&rings.outer).chain(&rings.holes) {
        let n = ring.len() as u32;
        for k in 0..n {
            ring_edges.push((start + k, start + (k + 1) % n));
        }
        start += n;
    }
    for &(a, b) in &ring_edges {
        assert_eq!(
            directed.get(&(a, b)).copied().unwrap_or(0),
            1,
            "{what}: ring edge {a}->{b} is not covered exactly once"
        );
        assert_eq!(
            directed.get(&(b, a)).copied().unwrap_or(0),
            0,
            "{what}: ring edge {a}->{b} is covered from outside"
        );
    }
    for (&(a, b), &count) in &directed {
        assert_eq!(count, 1, "{what}: edge {a}->{b} is used {count} times");
        if !ring_edges.contains(&(a, b)) {
            assert_eq!(
                directed.get(&(b, a)).copied().unwrap_or(0),
                1,
                "{what}: interior edge {a}->{b} has no twin (a T-junction or gap)"
            );
        }
    }

    // The extrusion is closed, outward and of volume area * depth.
    let mesh = extrude_profile(rings, Vec3::Z, DEPTH, Tolerance::METRE)
        .unwrap_or_else(|e| panic!("{what}: extrusion refused: {e:?}"));
    let volume = volume_properties(&mesh, Tolerance::METRE)
        .unwrap_or_else(|e| panic!("{what}: extrusion is not a closed two-manifold: {e:?}"))
        .signed_volume;
    assert_eq!(outward_orientation(&mesh), Some(true), "{what}: inside out");
    let expected = covered / 2.0 * DEPTH;
    assert!(
        (volume - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{what}: volume {volume}, expected {expected}"
    );
    covered
}

/// The polygonal area outer - sum(holes), doubled, from the rings alone.
fn twice_expected_area(rings: &Rings) -> f64 {
    twice_ring_area(&rings.outer).abs()
        - rings
            .holes
            .iter()
            .map(|h| twice_ring_area(h).abs())
            .sum::<f64>()
}

fn check(outer: Vec<Point2>, holes: Vec<Vec<Point2>>, what: &str) {
    let rings = Rings {
        outer,
        holes: holes.into_iter().map(cw).collect(),
    };
    let covered = assert_sound(&rings, what);
    assert_eq!(covered, twice_expected_area(&rings), "{what}: area");
}

#[test]
fn two_holes_side_by_side_close_the_extrusion() {
    // The exact profile of the issue.
    check(
        rect(-2.0, -2.0, 2.0, 2.0),
        vec![rect(-1.5, -0.5, -0.5, 0.5), rect(0.5, -0.5, 1.5, 0.5)],
        "side by side",
    );
}

#[test]
fn two_holes_stacked_close_the_extrusion() {
    check(
        rect(-2.0, -2.0, 2.0, 2.0),
        vec![rect(-0.5, -1.5, 0.5, -0.5), rect(-0.5, 0.5, 0.5, 1.5)],
        "stacked",
    );
}

#[test]
fn three_holes_on_one_ray_close_the_extrusion() {
    check(
        rect(-2.0, -2.0, 2.0, 2.0),
        vec![
            rect(-1.75, -0.25, -1.25, 0.25),
            rect(-0.25, -0.25, 0.25, 0.25),
            rect(1.25, -0.25, 1.75, 0.25),
        ],
        "three in a row",
    );
}

#[test]
fn a_hole_vertex_on_another_holes_ray_closes_the_extrusion() {
    // The square's right corners lie at y = +-0.5. The diamond's left
    // vertex lies exactly on the ray from the upper one, and the
    // triangle's apex exactly on the ray from the lower one.
    let diamond = vec![
        Point2::new(0.5, 0.5),
        Point2::new(1.0, 0.0),
        Point2::new(1.5, 0.5),
        Point2::new(1.0, 1.0),
    ];
    let triangle = vec![
        Point2::new(0.25, -1.5),
        Point2::new(1.5, -1.5),
        Point2::new(1.0, -0.5),
    ];
    check(
        rect(-2.0, -2.0, 2.0, 2.0),
        vec![rect(-1.5, -0.5, -0.5, 0.5), diamond, triangle],
        "vertex on a ray",
    );
}

#[test]
fn holes_close_to_the_outer_boundary_close_the_extrusion() {
    let gap = 1.0 / 1024.0;
    check(
        rect(-2.0, -2.0, 2.0, 2.0),
        vec![
            rect(-2.0 + gap, -2.0 + gap, -1.0, -1.0),
            rect(1.0, -2.0 + gap, 2.0 - gap, -1.0),
            rect(1.0, 1.0, 2.0 - gap, 2.0 - gap),
            rect(-2.0 + gap, 1.0, -1.0, 2.0 - gap),
            rect(-0.5, -2.0 + gap, 0.5, 2.0 - gap),
        ],
        "near the boundary",
    );
}

#[test]
fn holes_given_counter_clockwise_and_a_clockwise_outer_still_triangulate() {
    // `triangulate` takes rings either way round; only the triangles' own
    // winding is fixed.
    let rings = Rings {
        outer: cw(rect(-2.0, -2.0, 2.0, 2.0)),
        holes: vec![rect(-1.5, -0.5, -0.5, 0.5), rect(0.5, -0.5, 1.5, 0.5)],
    };
    let (points, triangles) = triangulate(&rings).expect("triangulates");
    let covered: f64 = triangles
        .iter()
        .map(|t| {
            let doubled = twice_area(
                points[t[0] as usize],
                points[t[1] as usize],
                points[t[2] as usize],
            );
            assert!(doubled > 0.0, "triangle {t:?} is not counter-clockwise");
            doubled
        })
        .sum();
    assert_eq!(covered, 28.0);
    assert_eq!(triangles.len(), 12 + 4 - 2);
}

#[test]
fn a_vertex_in_line_with_an_edge_past_its_end_does_not_touch_it() {
    // The triangle's top-left vertex is in line with the square's right
    // side, half a unit past its top, and the triangle's edge from there
    // overlaps that side's bounding box. Then the same turned a quarter.
    let square = rect(-1.0, -1.0, -0.5, 0.0);
    let triangle = vec![
        Point2::new(-0.5, 0.5),
        Point2::new(0.0, -0.5),
        Point2::new(0.5, 0.5),
    ];
    let turn =
        |ring: &[Point2]| -> Vec<Point2> { ring.iter().map(|p| Point2::new(-p.y, p.x)).collect() };
    check(
        rect(-2.0, -2.0, 2.0, 2.0),
        vec![square.clone(), triangle.clone()],
        "in line with a vertical side",
    );
    check(
        rect(-2.0, -2.0, 2.0, 2.0),
        vec![turn(&square), turn(&triangle)],
        "in line with a horizontal side",
    );
}

/// Every cyclic rotation of `ring`, so each vertex takes a turn as the
/// first one: the walk's start, and the vertex orientation is not read at.
fn rotations(ring: &[Point2]) -> Vec<Vec<Point2>> {
    (0..ring.len())
        .map(|r| ring[r..].iter().chain(&ring[..r]).copied().collect())
        .collect()
}

#[test]
fn a_straight_vertex_blocks_the_ear_across_it_from_every_start() {
    // A triangle with a vertex in the middle of its base: the ear at the
    // apex would cut straight along the base past it. And an L whose
    // reflex corner comes first in one rotation.
    let triangle = vec![
        Point2::new(-1.0, -1.0),
        Point2::new(0.0, -1.0),
        Point2::new(1.0, -1.0),
        Point2::new(0.0, 1.0),
    ];
    let ell = vec![
        Point2::new(-1.0, -1.0),
        Point2::new(1.0, -1.0),
        Point2::new(1.0, 0.0),
        Point2::new(0.0, 0.0),
        Point2::new(0.0, 1.0),
        Point2::new(-1.0, 1.0),
    ];
    for ring in [&triangle, &ell] {
        for (r, outer) in rotations(ring).into_iter().enumerate() {
            check(outer, Vec::new(), &format!("outer rotation {r}"));
        }
        for (r, hole) in rotations(ring).into_iter().enumerate() {
            let hole = hole
                .iter()
                .map(|p| Point2::new(p.x / 2.0, p.y / 2.0))
                .collect();
            check(
                rect(-2.0, -2.0, 2.0, 2.0),
                vec![hole],
                &format!("hole rotation {r}"),
            );
        }
    }
}

fn line(from: Point2, to: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}

fn polygon(points: &[Point2]) -> Contour {
    Contour::new(
        (0..points.len())
            .map(|i| line(points[i], points[(i + 1) % points.len()]))
            .collect(),
    )
}

/// A full circle as four quarter arcs.
fn circle(centre: Point2, radius: f64) -> Contour {
    let frame = Frame2 {
        origin: centre,
        x: Vec2::X,
        y: Vec2::Y,
    };
    let quarter = core::f64::consts::FRAC_PI_2;
    Contour::new(
        (0..4)
            .map(|index| ProfileSegment {
                curve: Curve2::Circle(Circle2 { frame, radius }),
                domain: Interval::new(quarter * index as f64, quarter * (index + 1) as f64),
                same_sense: true,
            })
            .collect(),
    )
}

#[test]
fn round_holes_side_by_side_close_the_extrusion() {
    // Flattened arcs: two round holes in one band and a third on the same
    // horizontal line, the case that left the cap open for polygons.
    let profile = Profile::Contour(ContourProfile {
        outer: polygon(&rect(-3.0, -1.5, 3.0, 1.5)),
        holes: vec![
            circle(Point2::new(-1.75, 0.0), 0.75),
            circle(Point2::new(0.0, 0.0), 0.75),
            circle(Point2::new(1.75, 0.0), 0.75),
        ],
    });
    for chord in [1e-2, 1e-3] {
        let rings = profile_rings(&profile, chord, Tolerance::METRE).expect("rings");
        let covered = assert_sound(&rings, "round holes");
        let expected = twice_expected_area(&rings);
        assert!(
            (covered - expected).abs() <= 1e-12 * expected,
            "chord {chord}: covered {covered}, rings {expected}"
        );
    }
}

/// A small deterministic generator: no dependency, same layouts every run.
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

/// One hole shape, in grid units, placed with its lower-left box corner at
/// the origin and fitting a `w x h` box. Every shape has its box's corners
/// or edge points as vertices, so neighbouring holes share rays, columns
/// and edge lines constantly.
fn shape(kind: u64, w: f64, h: f64) -> Vec<Point2> {
    match kind {
        0 => rect(0.0, 0.0, w, h),
        // A triangle with its apex on the box's top edge.
        1 => vec![
            Point2::new(0.0, 0.0),
            Point2::new(w, 0.0),
            Point2::new(w / 2.0, h),
        ],
        // An L, reflex at its inner corner.
        2 => vec![
            Point2::new(0.0, 0.0),
            Point2::new(w, 0.0),
            Point2::new(w, h / 2.0),
            Point2::new(w / 2.0, h / 2.0),
            Point2::new(w / 2.0, h),
            Point2::new(0.0, h),
        ],
        // A diamond: vertices at the box's edge midpoints.
        3 => vec![
            Point2::new(w / 2.0, 0.0),
            Point2::new(w, h / 2.0),
            Point2::new(w / 2.0, h),
            Point2::new(0.0, h / 2.0),
        ],
        // A rectangle with a collinear midpoint on every side.
        _ => vec![
            Point2::new(0.0, 0.0),
            Point2::new(w / 2.0, 0.0),
            Point2::new(w, 0.0),
            Point2::new(w, h / 2.0),
            Point2::new(w, h),
            Point2::new(w / 2.0, h),
            Point2::new(0.0, h),
            Point2::new(0.0, h / 2.0),
        ],
    }
}

/// Random hole layouts on a coarse grid inside a 16 x 16 outer square:
/// boxes never overlap, but they share rows, columns and lines all the
/// time, and some come within a quarter of the outer boundary.
#[test]
fn random_hole_layouts_triangulate_soundly() {
    let mut rng = Lcg(0x253);
    let mut checked = 0;
    for round in 0..400 {
        let (outer, holes) = hole_layout(&mut rng, round);
        check(outer, holes, &format!("round {round}"));
        checked += 1;
    }
    assert_eq!(checked, 400);
}

/// The layouts above moved so that one hole's lowest corner is the origin,
/// with zero coordinates written as -0.0 at random, as projecting onto
/// plane axes writes them. Orientation is read at a ring's lowest corner,
/// which a `-0.0` twin of its `x` used to displace (#269).
#[test]
fn random_hole_layouts_with_signed_zeros_triangulate_soundly() {
    let mut rng = Lcg(0x269);
    let mut moved = 0;
    for round in 0..400 {
        let (mut outer, mut holes) = hole_layout(&mut rng, round);
        if holes.is_empty() {
            continue;
        }
        let pick = rng.below(holes.len() as u64) as usize;
        let origin = holes[pick]
            .iter()
            .copied()
            .min_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)))
            .expect("a hole has corners");
        let mut sign = rng.next() << 31 ^ rng.next();
        let mut place = |q: Point2| {
            let mut coordinate = |v: f64| {
                if v != 0.0 {
                    return v;
                }
                sign = sign.rotate_left(1);
                if sign & 1 == 1 {
                    -0.0
                } else {
                    0.0
                }
            };
            let x = coordinate(q.x - origin.x);
            Point2::new(x, coordinate(q.y - origin.y))
        };
        outer = outer.into_iter().map(&mut place).collect();
        for hole in &mut holes {
            *hole = hole.iter().map(|&q| place(q)).collect();
        }
        moved += usize::from(
            holes
                .iter()
                .flatten()
                .any(|q| q.x.to_bits() == (-0.0f64).to_bits()),
        );
        let what = format!("signed-zero round {round}");
        // Wound the other way too: a straight corner taken for the lowest
        // one turns neither way, which reads as clockwise, so only rings
        // given counter-clockwise show the misreading.
        let flipped = Rings {
            outer: cw(outer.clone()),
            holes: holes.clone(),
        };
        let (points, triangles) =
            triangulate(&flipped).unwrap_or_else(|e| panic!("{what}, flipped: {e:?}"));
        let n = flipped.outer.len() + flipped.holes.iter().map(Vec::len).sum::<usize>();
        assert_eq!(triangles.len(), n + 2 * holes.len() - 2, "{what}, flipped");
        let covered: f64 = triangles
            .iter()
            .map(|t| {
                twice_area(
                    points[t[0] as usize],
                    points[t[1] as usize],
                    points[t[2] as usize],
                )
            })
            .sum();
        assert_eq!(
            covered,
            twice_expected_area(&flipped),
            "{what}, flipped: area"
        );
        check(outer, holes, &what);
    }
    assert!(moved >= 150, "only {moved} layouts have a -0.0 corner");
}

/// One random layout of [`random_hole_layouts_triangulate_soundly`]: the
/// outer ring and its holes, counter-clockwise.
fn hole_layout(rng: &mut Lcg, round: usize) -> (Vec<Point2>, Vec<Vec<Point2>>) {
    // Quarter-unit grid, so every coordinate and every product below is
    // exact in f64.
    let q = 0.25;
    let mut boxes: Vec<(f64, f64, f64, f64)> = Vec::new();
    let target = 1 + rng.below(9) as usize;
    for _ in 0..target * 8 {
        if boxes.len() == target {
            break;
        }
        let w = (2 + rng.below(16)) as f64 * q;
        let h = (2 + rng.below(16)) as f64 * q;
        let x = (1 + rng.below(63)) as f64 * q - 8.0;
        let y = (1 + rng.below(63)) as f64 * q - 8.0;
        if x + w > 8.0 - q || y + h > 8.0 - q {
            continue;
        }
        // Boxes keep at least a quarter apart, so even shapes that fill
        // their box (rectangles) never touch.
        let clear = boxes.iter().all(|&(bx, by, bw, bh)| {
            x >= bx + bw + q || bx >= x + w + q || y >= by + bh + q || by >= y + h + q
        });
        if clear {
            boxes.push((x, y, w, h));
        }
    }
    let holes: Vec<Vec<Point2>> = boxes
        .iter()
        .map(|&(x, y, w, h)| {
            let kind = rng.below(5);
            shape(kind, w, h)
                .into_iter()
                .map(|p| Point2::new(p.x + x, p.y + y))
                .collect()
        })
        .collect();
    // The outer ring sometimes carries collinear points too.
    let outer = if round % 3 == 0 {
        shape(4, 16.0, 16.0)
            .into_iter()
            .map(|p| Point2::new(p.x - 8.0, p.y - 8.0))
            .collect()
    } else {
        rect(-8.0, -8.0, 8.0, 8.0)
    };
    (outer, holes)
}

/// Random regular polygons at arbitrary (non-dyadic) positions and turns
/// inside a 96-gon: the predicates see no grid, rings of up to 40
/// vertices, and holes as close as a hundredth of their size.
#[test]
fn random_off_grid_hole_layouts_triangulate_soundly() {
    let mut rng = Lcg(0x2530);
    let mut unit = move || rng.next() as f64 / (1u64 << 31) as f64;
    let regular = |centre: Point2, radius: f64, sides: usize, turn: f64| -> Vec<Point2> {
        (0..sides)
            .map(|k| {
                let t = turn + k as f64 * core::f64::consts::TAU / sides as f64;
                Point2::new(centre.x + radius * t.cos(), centre.y + radius * t.sin())
            })
            .collect()
    };
    for round in 0..200 {
        let outer = regular(Point2::new(0.0, 0.0), 10.0, 96, unit());
        let mut discs: Vec<(Point2, f64)> = Vec::new();
        let target = 1 + (unit() * 12.0) as usize;
        for _ in 0..target * 10 {
            if discs.len() == target {
                break;
            }
            let radius = 0.2 + 1.8 * unit();
            let centre = Point2::new(16.0 * unit() - 8.0, 16.0 * unit() - 8.0);
            // Stay inside the 96-gon's inscribed circle, and clear of other
            // discs by a hundredth of the radius.
            let inside = centre.x.hypot(centre.y) + radius * 1.01 < 10.0 * 0.999;
            let clear = discs
                .iter()
                .all(|&(c, r)| (c.x - centre.x).hypot(c.y - centre.y) > (r + radius) * 1.01);
            if inside && clear {
                discs.push((centre, radius));
            }
        }
        let holes: Vec<Vec<Point2>> = discs
            .iter()
            .map(|&(c, r)| cw(regular(c, r, 3 + (unit() * 38.0) as usize, unit())))
            .collect();
        let rings = Rings { outer, holes };
        let what = format!("off-grid round {round}");
        let covered = assert_sound(&rings, &what);
        let expected = twice_expected_area(&rings);
        assert!(
            (covered - expected).abs() <= 1e-12 * expected,
            "{what}: covered {covered}, rings {expected}"
        );
    }
}

fn refused(outer: Vec<Point2>, holes: Vec<Vec<Point2>>) -> String {
    let rings = Rings {
        outer,
        holes: holes.into_iter().map(cw).collect(),
    };
    match triangulate(&rings) {
        Err(GeomError::InvalidInput(message)) => message,
        other => panic!("expected an invalid-input refusal, got {other:?}"),
    }
}

#[test]
fn overlapping_holes_are_refused_by_name() {
    let message = refused(
        rect(-2.0, -2.0, 2.0, 2.0),
        vec![rect(-1.0, -0.5, 0.5, 0.5), rect(0.0, -0.5, 1.5, 0.5)],
    );
    assert!(message.contains("holes 0 and 1 overlap"), "{message}");
}

#[test]
fn holes_touching_at_a_vertex_are_refused_by_name() {
    let message = refused(
        rect(-2.0, -2.0, 2.0, 2.0),
        vec![rect(-1.0, -0.5, 0.0, 0.5), rect(0.0, 0.5, 1.0, 1.5)],
    );
    assert!(
        message.contains("holes 0 and 1 overlap or touch"),
        "{message}"
    );
}

#[test]
fn a_hole_inside_another_hole_is_refused_by_name() {
    let message = refused(
        rect(-2.0, -2.0, 2.0, 2.0),
        vec![rect(-1.5, -1.5, 1.5, 1.5), rect(-0.5, -0.5, 0.5, 0.5)],
    );
    assert!(message.contains("hole 1 lies inside hole 0"), "{message}");
}

#[test]
fn a_hole_outside_the_outer_ring_is_refused_by_name() {
    let message = refused(
        rect(-2.0, -2.0, 2.0, 2.0),
        vec![rect(-1.5, -0.5, -0.5, 0.5), rect(3.0, -0.5, 4.0, 0.5)],
    );
    assert!(
        message.contains("hole 1 lies outside the outer ring"),
        "{message}"
    );
}

#[test]
fn a_hole_crossing_the_outer_ring_is_refused_by_name() {
    let message = refused(rect(-2.0, -2.0, 2.0, 2.0), vec![rect(1.5, -0.5, 2.5, 0.5)]);
    assert!(
        message.contains("hole 0 touches or crosses the outer ring"),
        "{message}"
    );
}

#[test]
fn a_hole_touching_the_outer_ring_at_one_vertex_is_refused_by_name() {
    // A triangle's apex on the outer ring's left side, then on its right:
    // once the outer edge is swept first, once the hole's.
    for (apex, base_x) in [(-2.0, -1.0), (2.0, 1.0)] {
        let triangle = vec![
            Point2::new(apex, 0.0),
            Point2::new(base_x, -0.5),
            Point2::new(base_x, 0.5),
        ];
        // Counter-clockwise either side, so `refused` turns it clockwise.
        let triangle = if apex < 0.0 { triangle } else { cw(triangle) };
        let message = refused(rect(-2.0, -2.0, 2.0, 2.0), vec![triangle]);
        assert!(
            message.contains("hole 0 touches or crosses the outer ring"),
            "apex at x = {apex}: {message}"
        );
    }
}

#[test]
fn a_self_intersecting_outer_ring_is_refused_by_name() {
    let bowtie = vec![
        Point2::new(-2.0, -2.0),
        Point2::new(2.0, 2.0),
        Point2::new(2.0, -2.0),
        Point2::new(-2.0, 2.0),
    ];
    let message = refused(bowtie, Vec::new());
    assert!(
        message.contains("outer ring intersects itself"),
        "{message}"
    );
}

#[test]
fn a_ring_folding_back_on_itself_is_refused_by_name() {
    // The second vertex runs past the first along the same line.
    let spike = vec![
        Point2::new(-2.0, -2.0),
        Point2::new(2.0, -2.0),
        Point2::new(0.0, -2.0),
        Point2::new(0.0, 2.0),
    ];
    let message = refused(spike, Vec::new());
    assert!(
        message.contains("outer ring folds back on itself"),
        "{message}"
    );
}

#[test]
fn a_repeated_vertex_is_refused_by_name() {
    let mut outer = rect(-2.0, -2.0, 2.0, 2.0);
    outer.insert(1, outer[1]);
    let message = refused(outer, Vec::new());
    assert!(message.contains("outer ring repeats vertex 1"), "{message}");
}

#[test]
fn a_non_finite_vertex_is_refused_by_name() {
    let mut outer = rect(-2.0, -2.0, 2.0, 2.0);
    outer[2].x = f64::NAN;
    let message = refused(outer, Vec::new());
    assert!(
        message.contains("outer ring has a non-finite vertex"),
        "{message}"
    );
}
