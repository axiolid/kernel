//! Every region an operation returns is a region `Region::new` accepts,
//! and can be fed into the next operation (#188 follow-up).
//!
//! Unions of dilated and visible regions share collinear edges, touch at
//! vertices and pinch holes off against their outer rings; the backends
//! also leave edges shorter than the tolerance. Outputs are settled so
//! that none of that reaches the caller.

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{Polygon, Region, Ring};

fn t() -> Tolerance {
    Tolerance::MILLIMETRE
}

fn ring(points: &[(f64, f64)]) -> Ring {
    Ring {
        points: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
    }
}

fn rect(x: f64, y: f64, w: f64, h: f64) -> Region {
    Region::new(
        vec![Polygon {
            outer: ring(&[(x, y), (x + w, y), (x + w, y + h), (x, y + h)]),
            holes: Vec::new(),
        }],
        t(),
    )
    .unwrap()
}

/// The region is accepted as it is, with the same area.
fn round_trips(what: &str, r: &Region) {
    let back = Region::new(r.polygons().to_vec(), t())
        .unwrap_or_else(|e| panic!("{what}: Region::new refused an operation's own output: {e:?}"));
    assert!(
        (back.area() - r.area()).abs() <= 1e-9 * (1.0 + r.area()),
        "{what}"
    );
}

#[test]
fn touching_collinear_and_pinched_shapes_round_trip() {
    let cases = [
        (
            "corner to corner",
            rect(0.0, 0.0, 1.0, 1.0).union(&rect(1.0, 1.0, 1.0, 1.0), t()),
        ),
        (
            "sharing an edge",
            rect(0.0, 0.0, 1.0, 1.0).union(&rect(1.0, 0.0, 1.0, 1.0), t()),
        ),
        (
            "sharing part of an edge",
            rect(0.0, 0.0, 2.0, 1.0).union(&rect(1.0, 1.0, 2.0, 1.0), t()),
        ),
        ("a frame round a hole", {
            let a = rect(0.0, 0.0, 3.0, 1.0)
                .union(&rect(0.0, 2.0, 3.0, 1.0), t())
                .unwrap();
            a.union(&rect(0.0, 0.0, 1.0, 3.0), t())
                .and_then(|r| r.union(&rect(2.0, 0.0, 1.0, 3.0), t()))
        }),
        ("a hole touching the outer ring", {
            let notch = Region::new(
                vec![Polygon {
                    outer: ring(&[(0.0, 2.0), (2.0, 1.0), (2.0, 3.0)]),
                    holes: Vec::new(),
                }],
                t(),
            )
            .unwrap();
            rect(0.0, 0.0, 4.0, 4.0).difference(&notch, t())
        }),
        ("two holes touching at a corner", {
            rect(0.0, 0.0, 6.0, 6.0)
                .difference(&rect(1.0, 1.0, 2.0, 2.0), t())
                .and_then(|r| r.difference(&rect(3.0, 3.0, 2.0, 2.0), t()))
        }),
        ("a ring closed on itself at a vertex", {
            // A U whose arms meet at one corner: the enclosed square is a
            // hole touching the outer boundary.
            rect(0.0, 0.0, 3.0, 1.0)
                .union(&rect(0.0, 0.0, 1.0, 3.0), t())
                .and_then(|r| r.union(&rect(2.0, 0.0, 1.0, 2.0), t()))
                .and_then(|r| r.union(&rect(1.0, 2.0, 1.0, 1.0), t()))
        }),
    ];
    for (what, r) in cases {
        round_trips(what, &r.unwrap());
    }
}

#[test]
fn chains_of_operations_accept_their_own_outputs() {
    let mut s = 0x1234_5678_9abc_def1u64;
    let mut next = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    for case in 0..120 {
        let mut acc = Region::empty();
        let mut cells = std::collections::BTreeSet::new();
        for _ in 0..2 + (next() * 4.0) as usize {
            let (x, y) = ((next() * 8.0).round() * 0.5, (next() * 8.0).round() * 0.5);
            let (w, h) = (
                0.5 + (next() * 4.0).round() * 0.5,
                0.5 + (next() * 4.0).round() * 0.5,
            );
            let mut piece = rect(x, y, w, h);
            match case % 3 {
                0 => {
                    for i in 0..(w * 2.0) as i64 {
                        for j in 0..(h * 2.0) as i64 {
                            cells.insert(((x * 2.0) as i64 + i, (y * 2.0) as i64 + j));
                        }
                    }
                }
                1 => {
                    piece = piece
                        .dilate(0.25 + (next() * 4.0).round() * 0.25, t())
                        .unwrap()
                }
                _ => piece = piece.dilate_outer(0.5, t()).unwrap(),
            }
            round_trips("dilation", &piece);
            let both = acc.intersection(&piece, t()).unwrap();
            let union = acc.union(&piece, t()).unwrap();
            round_trips("intersection", &both);
            round_trips("union", &union);
            // Inclusion-exclusion, up to what settling may move: a
            // tolerance along the boundary.
            let slack = 1e-3 * 200.0;
            assert!((union.area() + both.area() - acc.area() - piece.area()).abs() < slack);
            acc = union;
        }
        if case % 3 == 0 {
            // Grid-aligned rectangles: the union's area is exact.
            assert!(
                (acc.area() - 0.25 * cells.len() as f64).abs() < 1e-9,
                "case {case}"
            );
        }
        let hole = rect(1.0, 1.0, 2.0, 2.0);
        round_trips("difference", &acc.difference(&hole, t()).unwrap());
        round_trips("erosion", &acc.erode(0.3, t()).unwrap());
        if case % 10 == 0 {
            round_trips("inner erosion", &acc.erode_inner(0.3, t()).unwrap());
        }
        if let Ok(v) = acc.visibility_polygon(Point2::new(next() * 10.0, next() * 10.0), t()) {
            round_trips("visibility", &v);
            round_trips(
                "visibility grown",
                &acc.union(&v.dilate(0.3, t()).unwrap(), t()).unwrap(),
            );
        }
    }
}

#[test]
fn a_gap_under_the_tolerance_closes_round_a_hole() {
    // Three arms and a top bar stopping 0.5 mm short of the right arm:
    // closer than the tolerance, so the ring is closed there and the
    // square inside is a hole -- not filled, not a second polygon.
    let r = rect(0.0, 0.0, 3.0, 1.0)
        .union(&rect(0.0, 0.0, 1.0, 3.0), t())
        .and_then(|r| r.union(&rect(2.0, 0.0, 1.0, 2.0), t()))
        .and_then(|r| r.union(&rect(1.0, 2.0, 0.9995, 1.0), t()))
        .unwrap();
    round_trips("closed C", &r);
    assert!((r.area() - 6.9995).abs() < 1e-3, "{}", r.area());
}
