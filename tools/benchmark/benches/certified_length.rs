//! Wall-clock cost of certified joint distances along an offset by
//! distances beside an ellipse (#307).
//!
//! The deck of openbimrs/ifc#443: an ellipse of semi-axes 10 and 6 about
//! `(30, 0)`, offset by `0.5` to station 2, by distances widening to `0.75`
//! at station 20, by `0.75` on to the end. A windowed station certifies the
//! joints before it at `1e-3` of its window, so a `1e-5` window asks for
//! `1e-8`; before #307 the varying piece was refused, and the issue reports
//! uniform panels running out of budget at such tolerances. Two cases,
//! each at `1e-6`, `1e-8` and `1e-9`:
//!
//! - `joints`: both interior joints, the varying piece's length included,
//!   by `CompositeBasis::joint_intervals`;
//! - `station 30`: the distance a station 30 along with a window of a
//!   thousand times the tolerance is read at, which certifies the joints
//!   before it (`CompositeBasis::snapped_distance`).
//!
//! Measured on a shared 4-job machine (medians, noisy): `joints` 6 ms at
//! `1e-6`, 17 ms at `1e-8`, 34 ms at `1e-9`; `station 30` the same, its
//! certification dominating. A thousand times tighter costs under six
//! times as much: a panel's width falls as the fourth power of its size.
//!
//! Each case is validated before it is timed: every interval is no wider
//! than its tolerance and its middle lies within the tolerance plus the
//! offset quadrature's own tolerance (`1e-9` relative) of the kernel's
//! nominal joint, and the station is read where it stands (no seam lies
//! within its window).
//!
//! ```bash
//! cargo bench -p axiolid-benchmark --bench certified_length
//! ```

use std::hint::black_box;

use axiolid::evaluate::station::{offset_pieces, CompositeBasis, StationCurve, StationPiece};
use axiolid_core::{Frame2, Point2, Scalar, Vec2};
use axiolid_curve::{Curve2, Ellipse2, OffsetFrame, OffsetLaw, PathOffsets};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

/// The tolerances measured.
const TOLERANCES: [Scalar; 3] = [1e-6, 1e-8, 1e-9];

fn ellipse() -> Curve2 {
    Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::new(30.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: 10.0,
        semi_axis_y: 6.0,
    })
}

fn lateral(start: Scalar, end: Scalar) -> OffsetLaw {
    OffsetLaw::Linear {
        start: PathOffsets::new(start, 0.0, 0.0),
        end: PathOffsets::new(end, 0.0, 0.0),
        frame: OffsetFrame::Section,
    }
}

/// The deck's offset as a composite basis over `curve`.
fn deck(curve: &Curve2) -> CompositeBasis<'_> {
    let length = StationPiece::whole(StationCurve::Two(curve))
        .expect("an ellipse is a station basis")
        .length();
    let span = |a, b| StationPiece::between(StationCurve::Two(curve), a, b).expect("a span");
    let mut pieces = Vec::new();
    for (piece, law) in [
        (span(0.0, 2.0), lateral(0.5, 0.5)),
        (span(2.0, 20.0), lateral(0.5, 0.75)),
        (span(20.0, length), lateral(0.75, 0.75)),
    ] {
        pieces.extend(offset_pieces(&[piece], law).expect("an offset of the ellipse"));
    }
    CompositeBasis::new(pieces).expect("the offset's pieces meet")
}

/// The checks every case passes before it is timed.
fn validate(basis: &CompositeBasis<'_>, tolerance: Scalar) {
    let joints = basis
        .joint_intervals(tolerance)
        .expect("the deck's joints certify");
    let mut nominal = 0.0;
    for (k, interval) in joints.iter().enumerate() {
        nominal += basis.pieces()[k].length();
        let middle = 0.5 * (interval.lower + interval.upper);
        assert!(
            interval.width() <= tolerance,
            "{interval:?} at {tolerance:e}"
        );
        assert!(
            (middle - nominal).abs() <= tolerance + 1e-9 * nominal,
            "joint {k}: {interval:?} against {nominal}"
        );
    }
    let window = tolerance * 1e3;
    assert_eq!(
        basis
            .snapped_distance(30.0, window)
            .expect("the station certifies"),
        30.0
    );
}

fn certified_length(c: &mut Criterion) {
    let curve = ellipse();
    let basis = deck(&curve);
    let mut group = c.benchmark_group("certified_length");
    group.sample_size(10);
    for tolerance in TOLERANCES {
        validate(&basis, tolerance);
        group.bench_with_input(
            BenchmarkId::new("joints", format!("{tolerance:e}")),
            &tolerance,
            |b, &tolerance| b.iter(|| black_box(basis.joint_intervals(tolerance).unwrap())),
        );
        group.bench_with_input(
            BenchmarkId::new("station 30", format!("{tolerance:e}")),
            &tolerance,
            |b, &tolerance| {
                b.iter(|| black_box(basis.snapped_distance(30.0, tolerance * 1e3).unwrap()))
            },
        );
    }
    group.finish();
}

criterion_group!(benches, certified_length);
criterion_main!(benches);
