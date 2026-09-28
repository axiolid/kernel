//! Exact straight-edge boolean (#173).
//!
//! The polygon path used to hand its rings to an integer backend that maps
//! the operands' bounding box onto an `i32` grid and back, so every output
//! coordinate, even an input vertex the operation does not move, came back
//! snapped to a step of about 1.5e-8 of the extent. It now runs on the
//! exact subdivision the arc path uses (`exact_arc::arrangement`; a
//! straight ring is an arc ring without arcs):
//!
//! 1. Cut the plane by every ring of both operands at once. Where
//!    boundaries cross, which pieces coincide and which rings contain the
//!    region on each side of a piece are exact signs.
//! 2. Count, on each side of every piece, each operand's winding number:
//!    the rings containing that side, each counted `+1` when given
//!    counter-clockwise and `-1` when clockwise. The fill rule reads the
//!    count, the operation combines the two operands.
//! 3. Keep the pieces with the result on exactly one side, facing it, and
//!    link them into regions. A vertex where the boundary goes straight on
//!    (it only marked where another ring touched) is dropped, exactly.
//! 4. Round once: an input vertex is a double already and comes back
//!    bit-identical, and a crossing is the double nearest to it.
//!
//! The tolerance does not enter any decision. It validates operands, and
//! the result is settled like every other output ([`crate::settle`]), which
//! merges edges that rounding left shorter than the tolerance.

use crate::arc::{ArcRing, ArcVertex};
use crate::exact_arc::arrangement;
use crate::{signed, FillRule, OverlayError, OverlayOperation, Polygon, Ring};

/// Which operand a ring belongs to, and how it counts.
struct Counted {
    clip: bool,
    /// `+1` for a counter-clockwise ring as given, `-1` for clockwise.
    weight: i64,
}

fn filled(fill: FillRule, winding: i64) -> bool {
    match fill {
        FillRule::EvenOdd => winding % 2 != 0,
        FillRule::NonZero => winding != 0,
        FillRule::Positive => winding > 0,
        FillRule::Negative => winding < 0,
    }
}

fn combined(operation: OverlayOperation, subject: bool, clip: bool) -> bool {
    match operation {
        OverlayOperation::Intersection => subject && clip,
        OverlayOperation::Union => subject || clip,
        OverlayOperation::Difference => subject && !clip,
        OverlayOperation::Xor => subject != clip,
    }
}

/// The boolean of two sets of polygons, each ring counted by its own
/// orientation. Rings must be simple (validated); they may overlap each
/// other, within an operand too.
///
/// Returns the result's rings as `(outer, holes)`, before canonical
/// orientation and settling.
pub(crate) fn boolean(
    subject: &[Polygon],
    clip: &[Polygon],
    operation: OverlayOperation,
    fill: FillRule,
) -> Result<Vec<(Ring, Vec<Ring>)>, OverlayError> {
    let mut rings: Vec<ArcRing> = Vec::new();
    let mut counted: Vec<Counted> = Vec::new();
    for (is_clip, polygons) in [(false, subject), (true, clip)] {
        for ring in polygons
            .iter()
            .flat_map(|p| std::iter::once(&p.outer).chain(&p.holes))
        {
            let positive = signed(ring) > 0.0;
            let mut points = ring.points.clone();
            if !positive {
                points.reverse();
            }
            rings.push(ArcRing {
                vertices: points.into_iter().map(ArcVertex::straight).collect(),
            });
            counted.push(Counted {
                clip: is_clip,
                weight: if positive { 1 } else { -1 },
            });
        }
    }
    if rings.is_empty() {
        return Ok(Vec::new());
    }

    let raw = arrangement::build(&rings);
    // Whether the result holds the region on one side of a piece: the
    // rings containing that side are the ones containing the piece, and
    // each carrier whose inside faces that side.
    let result = |edge: &arrangement::RawEdge, left: bool| {
        let (mut subject, mut clip) = (0, 0);
        let rings = edge.inside.iter().copied().chain(
            edge.sources
                .iter()
                .filter(|source| source.2 == left)
                .map(|source| source.0),
        );
        for ring in rings {
            let count = &counted[ring];
            if count.clip {
                clip += count.weight;
            } else {
                subject += count.weight;
            }
        }
        combined(operation, filled(fill, subject), filled(fill, clip))
    };
    let keep: Vec<Option<bool>> = raw
        .edges
        .iter()
        .map(|edge| {
            let (left, right) = (result(edge, true), result(edge, false));
            // Kept pieces face the result: reversed when it lies right.
            (left != right).then_some(right)
        })
        .collect();

    let positions = raw.vertex_positions();
    let to_ring = |uses: &arrangement::RingUses| {
        let n = uses.len();
        let points = (0..n)
            .filter(|&i| !raw.straight_on(uses[(i + n - 1) % n], uses[i]))
            .map(|i| {
                let (piece, reversed) = uses[i];
                let edge = &raw.edges[piece];
                positions[if reversed { edge.to } else { edge.from }]
            })
            .collect();
        Ring { points }
    };
    Ok(raw
        .regions(&keep)?
        .iter()
        .map(|(outer, holes)| (to_ring(outer), holes.iter().map(to_ring).collect()))
        .collect())
}

#[cfg(test)]
mod tests {
    //! The exact boolean against the integer backend it replaced: the same
    //! regions, up to that backend's grid.
    use super::*;
    use axiolid_core::Point2;
    use i_overlay::core::{fill_rule::FillRule as Fill, overlay_rule::OverlayRule};
    use i_overlay::float::single::SingleFloatOverlay;

    /// A small deterministic generator, so failures reproduce.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// A star-shaped ring around `(cx, cy)`: simple by construction.
    fn star(rng: &mut Lcg, cx: f64, cy: f64, n: usize, clockwise: bool) -> Ring {
        let mut points: Vec<Point2> = (0..n)
            .map(|i| {
                let angle = std::f64::consts::TAU * (i as f64 + 0.8 * rng.next()) / n as f64;
                let radius = 0.4 + rng.next();
                Point2::new(cx + radius * angle.cos(), cy + radius * angle.sin())
            })
            .collect();
        if clockwise {
            points.reverse();
        }
        Ring { points }
    }

    /// Even-odd membership over all rings of a result.
    fn inside(polygons: &[(Ring, Vec<Ring>)], p: Point2) -> bool {
        let crosses = |ring: &Ring| {
            let n = ring.points.len();
            (0..n)
                .filter(|&i| {
                    let (a, b) = (ring.points[i], ring.points[(i + 1) % n]);
                    (a.y > p.y) != (b.y > p.y)
                        && p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x
                })
                .count()
                % 2
                == 1
        };
        polygons
            .iter()
            .flat_map(|(outer, holes)| std::iter::once(outer).chain(holes))
            .filter(|ring| crosses(ring))
            .count()
            % 2
            == 1
    }

    fn backend(polygons: &[Polygon]) -> Vec<Vec<Vec<[f64; 2]>>> {
        polygons
            .iter()
            .map(|p| {
                std::iter::once(&p.outer)
                    .chain(&p.holes)
                    .map(|r| r.points.iter().map(|q| [q.x, q.y]).collect())
                    .collect()
            })
            .collect()
    }

    fn distance_to_boundary(polygons: &[Polygon], p: Point2) -> f64 {
        let mut best = f64::INFINITY;
        for ring in polygons
            .iter()
            .flat_map(|q| std::iter::once(&q.outer).chain(&q.holes))
        {
            let n = ring.points.len();
            for i in 0..n {
                let (a, b) = (ring.points[i], ring.points[(i + 1) % n]);
                let t = ((p - a).dot(b - a) / (b - a).length_squared()).clamp(0.0, 1.0);
                best = best.min((a + (b - a) * t - p).length());
            }
        }
        best
    }

    #[test]
    fn agrees_with_the_integer_backend() {
        let mut rng = Lcg(173);
        let operations = [
            (OverlayOperation::Intersection, OverlayRule::Intersect),
            (OverlayOperation::Union, OverlayRule::Union),
            (OverlayOperation::Difference, OverlayRule::Difference),
            (OverlayOperation::Xor, OverlayRule::Xor),
        ];
        let fills = [
            (FillRule::EvenOdd, Fill::EvenOdd),
            (FillRule::NonZero, Fill::NonZero),
            (FillRule::Positive, Fill::Positive),
            (FillRule::Negative, Fill::Negative),
        ];
        let mut compared = 0;
        for case in 0..256 {
            let mut operand = |count: usize| -> Vec<Polygon> {
                (0..count)
                    .map(|_| {
                        let (cx, cy) = (2.0 * rng.next(), 2.0 * rng.next());
                        // Every eighth case has long rings, which the
                        // subdivision indexes rather than scans.
                        let n = if case % 8 == 7 {
                            40 + (rng.next() * 40.0) as usize
                        } else {
                            3 + (rng.next() * 9.0) as usize
                        };
                        let clockwise = rng.next() < 0.3;
                        Polygon {
                            outer: star(&mut rng, cx, cy, n, clockwise),
                            holes: Vec::new(),
                        }
                    })
                    .collect()
            };
            let subject = operand(1 + case % 3);
            let clip = operand(1 + case % 2);
            let (operation, rule) = operations[case % 4];
            let (fill, backend_fill) = fills[(case / 4) % 4];
            let exact = boolean(&subject, &clip, operation, fill).expect("simple rings link");
            let grid: Vec<(Ring, Vec<Ring>)> = backend(&subject)
                .overlay(&backend(&clip), rule, backend_fill)
                .into_iter()
                .map(|shape| {
                    let mut rings = shape.into_iter().map(|r| Ring {
                        points: r.into_iter().map(|p| Point2::new(p[0], p[1])).collect(),
                    });
                    let outer = rings.next().expect("a shape has an outer ring");
                    (outer, rings.collect())
                })
                .collect();
            let all: Vec<Polygon> = subject.iter().chain(&clip).cloned().collect();
            for _ in 0..200 {
                let p = Point2::new(4.0 * rng.next() - 1.0, 4.0 * rng.next() - 1.0);
                // Away from every input boundary, the grid cannot matter.
                if distance_to_boundary(&all, p) < 1e-6 {
                    continue;
                }
                compared += 1;
                assert_eq!(
                    inside(&exact, p),
                    inside(&grid, p),
                    "case {case}: {operation:?} {fill:?} at {p:?}"
                );
            }
        }
        assert!(compared > 40_000, "{compared}");
    }
}
