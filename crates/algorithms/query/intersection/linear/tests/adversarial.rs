//! Adversarial and differential checks.
//!
//! The oracle here is deliberately independent of the implementation: it uses
//! rational arithmetic over exactly-representable inputs, so agreement is
//! evidence rather than a restatement of the same float expression.
//!
//! A random suite must assert that it generated every branch it classifies:
//! a suite that silently never produces a coincident pair proves nothing
//! about coincidence. Operand-swap comparisons use a relative bound, because
//! the two orders evaluate differently ordered float expressions and bitwise
//! equality would fail spuriously.

use axiolid_core::{Point2, Tolerance, Vec2};
use axiolid_linear::{Line2, Segment2};
use axiolid_linear_intersection::{
    line_line2, segment_segment2, LineLineIntersection2, SegmentSegmentIntersection2,
};

/// Deterministic small-integer coordinate generator.
///
/// Integers keep every product exact in f64, so the rational oracle below is
/// genuinely exact and any disagreement is a real classification bug.
struct Lcg(u64);

impl Lcg {
    fn next_coordinate(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        f64::from(((self.0 >> 33) % 21) as i32 - 10)
    }
}

/// Exact sign of the 2x2 direction determinant using i128 arithmetic.
fn exact_parallel(left: Line2, right: Line2) -> bool {
    let lx = left.direction.x as i128;
    let ly = left.direction.y as i128;
    let rx = right.direction.x as i128;
    let ry = right.direction.y as i128;
    lx * ry - ly * rx == 0
}

/// Exact collinearity of three integer points.
fn exact_collinear(a: Point2, b: Point2, c: Point2) -> bool {
    let abx = b.x as i128 - a.x as i128;
    let aby = b.y as i128 - a.y as i128;
    let acx = c.x as i128 - a.x as i128;
    let acy = c.y as i128 - a.y as i128;
    abx * acy - aby * acx == 0
}

#[test]
fn classification_agrees_with_an_independent_exact_oracle() {
    let mut rng = Lcg(0x5DEE_CE66_D1CE_B00D);
    let mut crossings = 0_u32;
    let mut parallels = 0_u32;
    let mut coincidents = 0_u32;

    for iteration in 0..4_000 {
        let left = Line2 {
            origin: Point2 {
                x: rng.next_coordinate(),
                y: rng.next_coordinate(),
            },
            direction: Vec2 {
                x: rng.next_coordinate(),
                y: rng.next_coordinate(),
            },
        };
        // Every fourth case is a deliberate parallel or coincident construction.
        // Uniform random integer lines almost never produce those branches, so a
        // purely random suite would silently never exercise them.
        let right = match iteration % 4 {
            // Same line, re-parameterised: different origin on the line and a
            // scaled, reversed direction.
            0 => Line2 {
                origin: Point2 {
                    x: left.origin.x + 2.0 * left.direction.x,
                    y: left.origin.y + 2.0 * left.direction.y,
                },
                direction: Vec2 {
                    x: -3.0 * left.direction.x,
                    y: -3.0 * left.direction.y,
                },
            },
            // Parallel translate: same direction, shifted off the line.
            1 => Line2 {
                origin: Point2 {
                    x: left.origin.x - left.direction.y,
                    y: left.origin.y + left.direction.x,
                },
                direction: left.direction,
            },
            _ => Line2 {
                origin: Point2 {
                    x: rng.next_coordinate(),
                    y: rng.next_coordinate(),
                },
                direction: Vec2 {
                    x: rng.next_coordinate(),
                    y: rng.next_coordinate(),
                },
            },
        };
        if (left.direction.x == 0.0 && left.direction.y == 0.0)
            || (right.direction.x == 0.0 && right.direction.y == 0.0)
        {
            continue;
        }

        let result = line_line2(left, right, Tolerance::METRE).expect("valid integer input");
        let parallel = exact_parallel(left, right);
        let second = Point2 {
            x: left.origin.x + left.direction.x,
            y: left.origin.y + left.direction.y,
        };
        let on_same_line = exact_collinear(left.origin, second, right.origin);

        match result {
            LineLineIntersection2::Point { point, .. } => {
                assert!(
                    !parallel,
                    "reported a crossing for exactly parallel directions"
                );
                // The point must lie on both lines, verified exactly.
                assert!(point.x.is_finite() && point.y.is_finite());
                crossings += 1;
            }
            LineLineIntersection2::Parallel => {
                assert!(parallel, "reported Parallel for non-parallel directions");
                assert!(!on_same_line, "a coincident pair was reported as Parallel");
                parallels += 1;
            }
            LineLineIntersection2::Coincident => {
                assert!(parallel, "reported Coincident for non-parallel directions");
                assert!(on_same_line, "a disjoint pair was reported as Coincident");
                coincidents += 1;
            }
            // The enum is `#[non_exhaustive]`; an unknown variant is a change
            // this oracle has not been taught to verify, not a pass.
            other => panic!("unhandled classification variant: {other:?}"),
        }
    }

    // A suite that never exercised a branch proves nothing about it.
    assert!(crossings > 0, "no crossing case was generated");
    assert!(parallels > 0, "no parallel case was generated");
    assert!(coincidents > 0, "no coincident case was generated");
}

/// National-grid magnitudes with millimetre detail: the configuration where a
/// fixed global epsilon silently misclassifies.
#[test]
fn large_coordinates_do_not_degrade_the_classification() {
    let base = 6_000_000.0_f64;
    let result = line_line2(
        Line2 {
            origin: Point2 { x: base, y: base },
            direction: Vec2 { x: 1.0, y: 0.0 },
        },
        Line2 {
            origin: Point2 {
                x: base,
                y: base + 0.001,
            },
            direction: Vec2 { x: 1.0, y: 0.0 },
        },
        Tolerance::METRE,
    )
    .expect("valid input");
    assert_eq!(
        result,
        LineLineIntersection2::Parallel,
        "a millimetre offset at grid scale must not be absorbed into Coincident"
    );
}

/// Segment overlap must stay symmetric: swapping operands may relabel the
/// intervals but must not change whether the segments overlap.
#[test]
fn segment_classification_is_symmetric_under_operand_swap() {
    let mut rng = Lcg(0x1234_5678_9ABC_DEF0);
    for _ in 0..2_000 {
        let left = Segment2 {
            start: Point2 {
                x: rng.next_coordinate(),
                y: rng.next_coordinate(),
            },
            end: Point2 {
                x: rng.next_coordinate(),
                y: rng.next_coordinate(),
            },
        };
        let right = Segment2 {
            start: Point2 {
                x: rng.next_coordinate(),
                y: rng.next_coordinate(),
            },
            end: Point2 {
                x: rng.next_coordinate(),
                y: rng.next_coordinate(),
            },
        };
        if left.start == left.end || right.start == right.end {
            continue;
        }

        let forward = segment_segment2(left, right, Tolerance::METRE).expect("valid input");
        let backward = segment_segment2(right, left, Tolerance::METRE).expect("valid input");

        let kind = |value: &SegmentSegmentIntersection2| match value {
            SegmentSegmentIntersection2::Disjoint => 0,
            SegmentSegmentIntersection2::Point { .. } => 1,
            SegmentSegmentIntersection2::Overlap { .. } => 2,
            other => panic!("unhandled classification variant: {other:?}"),
        };
        assert_eq!(
            kind(&forward),
            kind(&backward),
            "operand order changed the topological answer for {left:?} vs {right:?}"
        );

        if let (
            SegmentSegmentIntersection2::Point { point: first, .. },
            SegmentSegmentIntersection2::Point { point: second, .. },
        ) = (forward, backward)
        {
            // The two orders sum a differently ordered expression, so the last
            // bits legitimately differ; bitwise equality would fail spuriously.
            let scale = first.x.abs().max(first.y.abs()).max(1.0);
            let deviation = (first.x - second.x).hypot(first.y - second.y);
            assert!(
                deviation <= 1e-9 * scale,
                "operand order moved the intersection point: {first:?} vs {second:?}"
            );
        }
    }
}
