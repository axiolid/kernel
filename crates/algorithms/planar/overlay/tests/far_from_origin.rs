//! Small rings far from the origin keep their area and orientation (#274).
//!
//! A ring's orientation and its `ZeroArea` check used to come from a
//! shoelace sum over the coordinates as given. At georeferenced plan
//! positions (about 6e5 and 5.6e6 m) each term is the size of the
//! coordinates squared, and its rounding swamps the area of a small ring:
//! a valid triangle was refused as `ZeroArea`, and a clockwise sliver read
//! as counter-clockwise made a union fail with `SelfIntersection`. Both
//! are now exact signs of the area taken relative to a ring vertex, so a
//! ring decides the same far from the origin as near it.
//!
//! Every ring here is built near the origin from dyadic coordinates and
//! translated by offsets that keep it exact, so the far ring is the same
//! ring, and its result must be the near result translated.

use axiolid_core::{Frame2, Point2, Tolerance, Vec2};
use axiolid_overlay::{
    arc_overlay, arc_ring_area, overlay, polygon_area, reverse_arc_ring, union_soup,
    validate_arc_ring, ArcArrangement, ArcRing, ArcVertex, FillRule, OverlayError, OverlayInput,
    OverlayOperation, Polygon, Region, Ring,
};

fn tol() -> Tolerance {
    Tolerance::new(1e-9, 1e-9).unwrap()
}

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(0.0, 0.0),
        x: Vec2::new(1.0, 0.0),
        y: Vec2::new(0.0, 1.0),
    }
}

fn poly(points: &[(f64, f64)]) -> Polygon {
    Polygon {
        outer: Ring {
            points: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
        },
        holes: Vec::new(),
    }
}

fn union(polygons: Vec<Polygon>, tolerance: Tolerance) -> Result<Vec<Polygon>, OverlayError> {
    let input = OverlayInput {
        frame: frame(),
        polygons,
    };
    overlay(
        &input,
        &input,
        OverlayOperation::Union,
        FillRule::NonZero,
        tolerance,
    )
    .map(|result| result.polygons)
}

/// Offsets, as `(x, y)`, at which every ring below translates exactly.
const OFFSETS: [(f64, f64); 5] = [
    (0.0, 0.0),
    (1e5, 1e5),
    (5e6, 5e6),
    (1e7, 1e7),
    (600_000.0, 5_600_000.0),
];

fn shifted(points: &[(f64, f64)], (dx, dy): (f64, f64)) -> Vec<(f64, f64)> {
    points
        .iter()
        .map(|&(x, y)| {
            let (sx, sy) = (x + dx, y + dy);
            assert_eq!((sx - dx, sy - dy), (x, y), "translation must be exact");
            (sx, sy)
        })
        .collect()
}

/// Every start vertex, in both windings.
fn variants(points: &[(f64, f64)]) -> Vec<Vec<(f64, f64)>> {
    let mut out = Vec::new();
    for reversed in [false, true] {
        let mut ring = points.to_vec();
        if reversed {
            ring.reverse();
        }
        for k in 0..ring.len() {
            let mut rotated = ring.clone();
            rotated.rotate_left(k);
            out.push(rotated);
        }
    }
    out
}

fn back(polygons: &[Polygon], (dx, dy): (f64, f64)) -> Vec<Polygon> {
    let ring = |r: &Ring| Ring {
        points: r
            .points
            .iter()
            .map(|p| Point2::new(p.x - dx, p.y - dy))
            .collect(),
    };
    polygons
        .iter()
        .map(|p| Polygon {
            outer: ring(&p.outer),
            holes: p.holes.iter().map(ring).collect(),
        })
        .collect()
}

/// The ring's union with itself and its validated region at every offset,
/// translated back, equal its own near the origin, for every start vertex
/// and both windings. Except for a ring on one line, which has no
/// orientation, the results do not depend on the start vertex or the
/// winding either. Returns the near union and the near region.
fn same_everywhere(
    points: &[(f64, f64)],
    tolerance: Tolerance,
) -> (
    Result<Vec<Polygon>, OverlayError>,
    Result<Vec<Polygon>, OverlayError>,
) {
    let region_of = |points: &[(f64, f64)]| {
        Region::new(vec![poly(points)], tolerance).map(|r| r.polygons().to_vec())
    };
    let near = union(vec![poly(points)], tolerance);
    let near_region = region_of(points);
    let degenerate = matches!(&near, Ok(polygons) if polygons.is_empty());
    for ring in variants(points) {
        let near_union = union(vec![poly(&ring)], tolerance);
        let near_ring = region_of(&ring);
        if !degenerate {
            assert_eq!(near_union, near, "union of {ring:?}");
            assert_eq!(near_ring, near_region, "region of {ring:?}");
        }
        for offset in OFFSETS {
            let far = union(vec![poly(&shifted(&ring, offset))], tolerance);
            let far_region = region_of(&shifted(&ring, offset));
            assert_eq!(
                far.map(|f| back(&f, offset)),
                near_union,
                "union of {ring:?} at {offset:?}"
            );
            assert_eq!(
                far_region.map(|f| back(&f, offset)),
                near_ring,
                "region of {ring:?} at {offset:?}"
            );
        }
    }
    (near, near_region)
}

const P: f64 = 1.0 / 1024.0;

#[test]
fn the_issue_triangle_is_accepted_far_from_the_origin() {
    let near = poly(&[(0.0, 0.0), (0.02, 0.0), (0.0, 0.01)]);
    let far = poly(&[
        (600_000.0, 5_600_000.0),
        (600_000.02, 5_600_000.0),
        (600_000.0, 5_600_000.01),
    ]);
    for t in [tol(), Tolerance::new(1e-6, 1e-6).unwrap()] {
        let near = union(vec![near.clone()], t).unwrap();
        let far = union(vec![far.clone()], t).expect("a 1e-4 m^2 triangle is not ZeroArea");
        assert_eq!((near.len(), far.len()), (1, 1));
        let (a, b) = (polygon_area(&near[0]), polygon_area(&far[0]));
        assert!(
            (a - 1e-4).abs() < 1e-12 && (b - 1e-4).abs() < 1e-9,
            "{a} {b}"
        );
    }
}

#[test]
fn the_issue_fan_is_united_far_from_the_origin() {
    let a = (600_026.029859079, 5_599_990.725323705);
    let b = (600_026.1426650629, 5_599_990.741066861);
    let c1 = (600_026.086448961, 5_599_990.731856141);
    let c2 = (600_026.2534703149, 5_599_990.767434134);
    let offset = (600_026.0, 5_599_990.0);
    let local = |(x, y): (f64, f64)| (x - offset.0, y - offset.1);
    let far = union(vec![poly(&[a, b, c1]), poly(&[a, b, c2])], tol())
        .expect("two triangles sharing an edge unite");
    let near = union(
        vec![
            poly(&[local(a), local(b), local(c1)]),
            poly(&[local(a), local(b), local(c2)]),
        ],
        tol(),
    )
    .unwrap();
    assert_eq!(far.len(), 1);
    assert_eq!(back(&far, offset), near);
    // The clockwise sliver alone keeps its orientation too.
    let sliver = union(vec![poly(&[a, b, c1])], tol()).unwrap();
    let sliver_near = union(vec![poly(&[local(a), local(b), local(c1)])], tol()).unwrap();
    assert_eq!(back(&sliver, offset), sliver_near);
    let soup = union_soup(&[poly(&[a, b, c1]).outer, poly(&[a, b, c2]).outer], tol()).unwrap();
    assert_eq!(soup, far);
}

#[test]
fn thin_triangles_decide_the_same_at_every_offset() {
    // Slivers down to 2^-24 m wide and 2^-31 m^2 in area.
    for (base, apex_x, height) in [
        (40.0 * P, 13.0 * P, P / 1024.0),
        (64.0 * P, 70.0 * P, P / 16384.0),
        (8.0 * P, 3.0 * P, P / 2048.0),
        (512.0 * P, -5.0 * P, P / 4096.0),
    ] {
        let triangle = [(0.0, 0.0), (base, 0.0), (apex_x, height)];
        let (near, region) = same_everywhere(&triangle, tol());
        assert_eq!(near.expect("a thin triangle is valid").len(), 1);
        assert_eq!(region.expect("a thin triangle is valid").len(), 1);
    }
}

#[test]
fn small_rings_decide_the_same_at_every_offset() {
    let rings: [&[(f64, f64)]; 3] = [
        // A small L.
        &[
            (0.0, 0.0),
            (3.0 * P, 0.0),
            (3.0 * P, P),
            (P, P),
            (P, 2.0 * P),
            (0.0, 2.0 * P),
        ],
        // A thin quadrilateral with a reflex vertex.
        &[
            (0.0, 0.0),
            (32.0 * P, 0.0),
            (16.0 * P, P / 512.0),
            (0.0, P / 64.0),
        ],
        // A small square.
        &[(0.0, 0.0), (P, 0.0), (P, P), (0.0, P)],
    ];
    for ring in rings {
        let (near, region) = same_everywhere(ring, tol());
        assert_eq!(near.expect("a small ring is valid").len(), 1);
        assert_eq!(region.expect("a small ring is valid").len(), 1);
    }
}

#[test]
fn rings_on_one_line_decide_the_same_at_every_offset() {
    // A ring on one line encloses nothing and a boolean leaves it out
    // (#219), wherever it lies; before, its rounded area decided whether
    // it was refused as ZeroArea, so the same ring passed at some offsets
    // and not at others.
    for ring in [
        [(0.0, 0.0), (P, P / 2.0), (2.0 * P, P)],
        [(0.0, 0.0), (2.0 * P, 0.0), (P, 0.0)],
        [(0.0, 0.0), (0.0, 3.0 * P), (0.0, P)],
    ] {
        let (near, region) = same_everywhere(&ring, tol());
        assert_eq!(near, Ok(Vec::new()), "{ring:?}");
        assert!(region.is_ok(), "{ring:?}");
    }
    // Off the line by one step in the last place at 1e7, the triangle has
    // an area of 2^-40 m^2 and is a ring like any other: ZeroArea
    // only at a tolerance whose square exceeds it.
    let thin = [(0.0, 0.0), (P, P / 2.0), (2.0 * P, P + P / 524_288.0)];
    assert_eq!(same_everywhere(&thin, tol()).1.map(|p| p.len()), Ok(1));
    let coarse = Tolerance::new(P / 64.0, P / 64.0).unwrap();
    assert_eq!(
        same_everywhere(&thin, coarse).1,
        Err(OverlayError::ZeroArea)
    );
}

#[test]
fn the_zero_area_threshold_is_exact_at_every_offset() {
    // Tolerance 2^-10: a ring of area at most 2^-20 is ZeroArea. This
    // triangle's area is exactly 2^-20; raised by 2^-28 in height it is
    // 2^-33 above the threshold. (The overlay then snaps the apex, which
    // lies within the tolerance of the base, onto it: #222.)
    let t = Tolerance::new(P, P).unwrap();
    let at = [(0.0, 0.0), (64.0 * P, 0.0), (32.0 * P, P / 32.0)];
    assert_eq!(same_everywhere(&at, t).1, Err(OverlayError::ZeroArea));
    let above = [
        (0.0, 0.0),
        (64.0 * P, 0.0),
        (32.0 * P, P / 32.0 + P / 262_144.0),
    ];
    assert_eq!(same_everywhere(&above, t).1.map(|p| p.len()), Ok(1));
}

/// A counter-clockwise triangle whose doubled area is exactly `+1` while
/// the `f64` fan from its first vertex gives `0`: `2^27 * 2^27` is exact,
/// `(2^27 + 1) * (2^27 - 1) = 2^54 - 1` rounds to `2^54`. Only an exact
/// sign tells it from a degenerate ring.
fn cancelling() -> [(f64, f64); 3] {
    let k = 134_217_728.0; // 2^27
    [(0.0, 0.0), (k, k - 1.0), (k + 1.0, k)]
}

#[test]
fn orientation_and_area_are_exact_where_the_rounded_fan_cancels() {
    let t = tol();
    let (near, region) = same_everywhere(&cancelling(), t);
    let expected = vec![poly(&cancelling())];
    assert_eq!(near.as_ref(), Ok(&expected), "kept counter-clockwise");
    assert_eq!(region, Ok(expected.clone()));
    // Given clockwise, it comes back counter-clockwise all the same.
    let mut clockwise = cancelling();
    clockwise.reverse();
    assert_eq!(union(vec![poly(&clockwise)], t), Ok(expected.clone()));
    assert_eq!(
        union_soup(&[poly(&clockwise).outer], t),
        Ok(expected.clone())
    );
    assert_eq!(union_soup(&[poly(&cancelling()).outer], t), Ok(expected));
}

fn arc_ring(points: &[(f64, f64)]) -> ArcRing {
    ArcRing::from_points(
        &points
            .iter()
            .map(|&(x, y)| Point2::new(x, y))
            .collect::<Vec<_>>(),
    )
}

#[test]
fn straight_arc_rings_take_the_exact_orientation() {
    let t = tol();
    let mut clockwise = cancelling();
    clockwise.reverse();
    for given in [cancelling(), clockwise] {
        let ring = arc_ring(&given);
        validate_arc_ring(&ring, t).expect("valid");
        let result = arc_overlay(&ring, &ring, OverlayOperation::Union, t).unwrap();
        assert_eq!(result.regions.len(), 1);
        let mut got: Vec<(f64, f64)> = result.regions[0]
            .outer
            .vertices
            .iter()
            .map(|v| (v.point.x, v.point.y))
            .collect();
        let k = (0..got.len())
            .min_by(|&a, &b| got[a].partial_cmp(&got[b]).unwrap())
            .unwrap();
        got.rotate_left(k);
        assert_eq!(got, cancelling().to_vec(), "counter-clockwise");

        let arrangement = ArcArrangement::new(std::slice::from_ref(&ring), t).unwrap();
        let regions = arrangement.regions(|inside| inside[0]).unwrap();
        assert_eq!(regions.len(), 1);
        assert!(regions[0].holes.is_empty());
    }
}

#[test]
fn small_arc_rings_keep_their_area_far_from_the_origin() {
    // A quarter-disc of radius 2^-7 m: two radii, then the quarter arc
    // back (bulge tan(pi / 8)).
    let r = 8.0 * P;
    let near_area = std::f64::consts::PI * r * r / 4.0;
    let bulge = (std::f64::consts::PI / 8.0).tan();
    for (dx, dy) in OFFSETS {
        let quarter = ArcRing::new(vec![
            ArcVertex::straight(Point2::new(dx, dy)),
            ArcVertex::bulged(Point2::new(dx + r, dy), bulge),
            ArcVertex::straight(Point2::new(dx, dy + r)),
        ]);
        validate_arc_ring(&quarter, tol()).expect("a small quarter-disc is valid");
        let area = arc_ring_area(&quarter);
        assert!((area - near_area).abs() <= 1e-12, "{area} at {dx}");
        let reversed = reverse_arc_ring(&quarter);
        assert!((arc_ring_area(&reversed) + near_area).abs() <= 1e-12);
        let result = arc_overlay(&reversed, &quarter, OverlayOperation::Union, tol()).unwrap();
        assert_eq!(result.regions.len(), 1, "at {dx}");
        assert!(result.regions[0].holes.is_empty());
        assert!(arc_ring_area(&result.regions[0].outer) > 0.0);
    }
}
