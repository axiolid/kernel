//! Wall-clock cost of the chord bound of a long arc-length chain (#438).
//!
//! `chord_bound2` of a chain charges the evaluated ends of its span against
//! a certified `ChainPositionTable` built over the chain's prefix up to the
//! span, so a call costs in proportion to the pieces before it, and a
//! flatten pays that at every span it tests. A 2 km chain of 30 pieces
//! (alternating tangent lines, clothoid transitions and arcs, with
//! parametric lines among them) makes that measurable:
//!
//! - `bound`: `chord_bound2` over a 1 m span inside the last piece, the
//!   worst case for the prefix;
//! - `bound/first`: the same inside the first piece, the prefix of one;
//! - `flatten`: `flatten2` of the whole chain at a 1 cm tolerance.
//!
//! Each case is validated before it is timed: the bound exists, and the
//! flattened polyline has more than the chain's pieces as points.
//!
//! ```bash
//! cargo bench -p axiolid-benchmark --bench chain_chord_bound
//! ```

use std::hint::black_box;

use axiolid::evaluate::bound::chord_bound2;
use axiolid::evaluate::flatten2;
use axiolid_core::{Frame2, Interval, Point2, Scalar, Vec2};
use axiolid_curve::{Chain2, ChainPiece2, CurvatureLaw, Curve2, Line2};
use criterion::{criterion_group, criterion_main, Criterion};

/// Pieces of the chain and the length of each, in metres.
const PIECES: usize = 30;
const PIECE: Scalar = 2000.0 / PIECES as Scalar;

/// A chain of 30 pieces: a tangent line (intrinsic or parametric), a
/// clothoid into a bend, an arc, a clothoid out, repeating.
fn chain() -> Chain2 {
    let k = 1.0 / 400.0;
    let pieces = (0..PIECES)
        .map(|i| match i % 5 {
            0 => ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::circular(0.0),
                length: PIECE,
            },
            1 => ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::clothoid(0.0, k, PIECE),
                length: PIECE,
            },
            2 => ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::circular(k),
                length: PIECE,
            },
            3 => ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::clothoid(k, 0.0, PIECE),
                length: PIECE,
            },
            _ => ChainPiece2::Parametric {
                curve: Curve2::Line(Line2 {
                    origin: Vec2::ZERO,
                    direction: Vec2::X,
                }),
                start: 0.0,
                length: PIECE,
            },
        })
        .collect();
    Chain2::new(
        Frame2 {
            origin: Point2::new(100.0, -50.0),
            x: Vec2::new(0.8, 0.6),
            y: Vec2::new(-0.6, 0.8),
        },
        pieces,
    )
}

fn bench(c: &mut Criterion) {
    let chain = chain();
    let curve = Curve2::Chain(chain.clone());
    let last = (PIECES as Scalar - 1.0) * PIECE;
    let spans = [(last + 10.0, last + 11.0), (10.0, 11.0)];
    for (lo, hi) in spans {
        let bound = chord_bound2(&curve, lo, hi).expect("a bound inside one piece");
        assert!(bound.is_finite() && bound >= 0.0);
        eprintln!("bound [{lo}, {hi}]: {bound:e}");
    }
    let total = PIECES as Scalar * PIECE;
    let points = flatten2(&curve, Interval::new(0.0, total), 1e-2, 30).expect("a polyline");
    assert!(points.len() > PIECES, "{} points", points.len());
    eprintln!("flatten: {} points", points.len());

    let mut group = c.benchmark_group("chain/chord_bound");
    group.sample_size(10);
    group.bench_function("bound last piece", |b| {
        b.iter(|| black_box(chord_bound2(black_box(&curve), spans[0].0, spans[0].1)))
    });
    group.bench_function("bound first piece", |b| {
        b.iter(|| black_box(chord_bound2(black_box(&curve), spans[1].0, spans[1].1)))
    });
    group.bench_function("flatten 2 km at 1 cm", |b| {
        b.iter(|| {
            black_box(flatten2(
                black_box(&curve),
                Interval::new(0.0, total),
                1e-2,
                30,
            ))
        })
    });
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
