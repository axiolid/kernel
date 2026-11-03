//! Wall-clock cost of certifying a run of sections along a clothoid plan
//! (#361).
//!
//! A 5 km alignment (line, clothoid, arc, clothoid, line) under a sag and
//! a crest, the deck of openbimrs/ifc alignments through
//! `rheonda/cad#139`. Its frames are certified by comparing every frame
//! read with a certified position table; a run's refinement reads them in
//! bisecting order, backward as often as forward. Two cases:
//!
//! - `positions`: the certified positions of the alignment at `1e-6`
//!   (a thousandth of a 1 mm budget), 512 points read in a bisecting
//!   order. `midpoint` is `IntrinsicPositions` (#317), whose forward march
//!   restarts on every backward read; `gauss2` is the checkpointed
//!   two-point Gauss-Legendre `IntrinsicPositionTable`, built once;
//! - `run`: a 12 m deck widening to 12 m along the gradient curve,
//!   compiled with its deviation report at a 1 mm chord error, certified;
//! - `placement`: 1000 blocks placed along the gradient curve, one every
//!   5 m (#405), compiled with their deviation report at a 1 mm chord
//!   error, every placement certified.
//!
//! Measured on a shared 4-job machine (medians, noisy): `positions`
//! 198 ms for the midpoint rule against 3.9 ms for the Gauss-Legendre
//! table, build included, about fifty times; `run` 367 ms with a certified
//! bound of 0.49 mm against the 1 mm request. `placement` (#405), on a
//! machine under a load of 70: 97 ms for the 1000 blocks (median of 10, 93
//! to 106 ms), every placement certified with a frame error of 0.32 um
//! against the 1 mm request; a measurement, not a performance claim.
//!
//! Each case is validated before it is timed: both tables agree to their
//! errors, and the run and the placements are `Certified` within their
//! request.
//!
//! ```bash
//! cargo bench -p axiolid-benchmark --bench station_clothoid
//! ```

use std::hint::black_box;

use axiolid::evaluate::bound::{IntrinsicPositionTable, IntrinsicPositions};
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame2, Point2, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Intrinsic2};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, ReferenceMeshCompiler};
use axiolid_model::{
    CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode, InstanceAtStation, NodeId,
    SectionAtStation, SolidOperation, Station, StationFrame, StationOffsets,
};
use axiolid_profile::{Profile, RectangleProfile};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

/// The chord error requested of the run, in metres.
const REQUEST: Scalar = 1e-3;

/// The error asked of the positions: a thousandth of the request.
const POSITION_ERROR: Scalar = 1e-6;

/// Line 1000, clothoid 0 to 1/1500 over 250, arc 2500, clothoid back 250,
/// line 1000: 5 km.
fn alignment() -> Intrinsic2 {
    let law = CurvatureLaw::piecewise(
        vec![1000.0, 1250.0, 3750.0, 4000.0],
        vec![
            CurvatureLaw::Constant { curvature: 0.0 },
            CurvatureLaw::clothoid(0.0, 1.0 / 1500.0, 250.0),
            CurvatureLaw::Constant {
                curvature: 1.0 / 1500.0,
            },
            CurvatureLaw::clothoid(1.0 / 1500.0, 0.0, 250.0),
            CurvatureLaw::Constant { curvature: 0.0 },
        ],
    );
    Intrinsic2::new(
        Frame2 {
            origin: Point2::new(3.0, -1.0),
            x: Vec2::new(0.6, 0.8),
            y: Vec2::new(-0.8, 0.6),
        },
        law,
        5000.0,
    )
}

/// `count` distances in the order a bisection reads them: the ends, then
/// the middle, then the quarters, and so on.
fn bisecting(length: Scalar, count: usize) -> Vec<Scalar> {
    let mut out = vec![0.0, length];
    let mut level = 1usize;
    'fill: loop {
        for i in 0..level {
            if out.len() >= count {
                break 'fill;
            }
            out.push(length * (2 * i + 1) as Scalar / (2 * level) as Scalar);
        }
        level *= 2;
    }
    out
}

fn positions(c: &mut Criterion) {
    let curve = alignment();
    let reads = bisecting(curve.length, 512);
    // Both agree to their errors before either is timed.
    let mut midpoint = IntrinsicPositions::new2(&curve, POSITION_ERROR).expect("midpoint");
    let table =
        IntrinsicPositionTable::new2(&curve, POSITION_ERROR, curve.length).expect("a table");
    eprintln!("gauss2 panels: {}", table.panels());
    for &s in reads.iter().take(64) {
        let (p, e) = midpoint.point2(s).unwrap();
        let (q, f) = table.point2(s).unwrap();
        assert!((p - q).length() <= e + f, "{s}");
    }
    let mut group = c.benchmark_group("positions");
    group.sample_size(10);
    group.bench_function(BenchmarkId::new("midpoint", "512 reads"), |b| {
        b.iter(|| {
            let mut positions = IntrinsicPositions::new2(&curve, POSITION_ERROR).unwrap();
            for &s in &reads {
                black_box(positions.point2(black_box(s)).unwrap());
            }
        })
    });
    group.bench_function(BenchmarkId::new("gauss2", "512 reads"), |b| {
        b.iter(|| {
            let table = IntrinsicPositionTable::new2(&curve, POSITION_ERROR, curve.length).unwrap();
            for &s in &reads {
                black_box(table.point2(black_box(s)).unwrap());
            }
        })
    });
    group.bench_function(BenchmarkId::new("gauss2 reads only", "512 reads"), |b| {
        b.iter(|| {
            for &s in &reads {
                black_box(table.point2(black_box(s)).unwrap());
            }
        })
    });
    group.finish();
}

fn rectangle(b: &mut GeometryGraphBuilder, x: Scalar) -> NodeId {
    b.push_value(Profile::Rectangle(RectangleProfile {
        x,
        y: 0.4,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    }))
    .unwrap()
}

/// A deck 10 m wide widening to 12 m along the gradient curve over the
/// alignment, at the stations 100, 1600 and 4700.
fn deck() -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let curve = Curve3::Elevated(Elevated3::new(
        Curve2::Intrinsic(alignment()),
        ElevationLaw::parabolic(10.0, 0.01, -0.01, 5000.0),
    ));
    let basis = b.push_value(curve).unwrap();
    let (narrow, wide) = (rectangle(&mut b, 10.0), rectangle(&mut b, 12.0));
    let sections = vec![
        SectionAtStation::new(narrow, Station::new(100.0, StationOffsets::ZERO)),
        SectionAtStation::new(
            wide,
            Station::new(1600.0, StationOffsets::new(1.5, 0.2, 0.0)),
        ),
        SectionAtStation::new(narrow, Station::new(4700.0, StationOffsets::ZERO)),
    ];
    let root = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix: basis,
                sections,
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    (b.finish(vec![root]).unwrap(), root)
}

fn run(c: &mut Criterion) {
    let (graph, root) = deck();
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(REQUEST)
        .expect("a chord error");
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    let (_, report) = compiler
        .compile_mesh_with_deviation(&graph, root, &options)
        .expect("the run compiles");
    let run = report
        .contributions
        .iter()
        .find(|c| c.path == DeviationPath::StationedSpine)
        .expect("the run reports its path");
    let DeviationBound::Certified(bound) = run.bound else {
        panic!("a certified bound expected, got {report:?}");
    };
    assert!(report.meets_requested(), "{bound:e} past {REQUEST:e}");
    eprintln!("run: bound {bound:e}");
    let mut group = c.benchmark_group("station/clothoid");
    group.sample_size(10);
    group.bench_function("run 5 km at 1 mm", |b| {
        b.iter(|| {
            black_box(
                compiler
                    .compile_mesh_with_deviation(black_box(&graph), root, &options)
                    .expect("a mesh and its report"),
            )
        })
    });
    group.finish();
}

/// 1000 blocks of 4 x 2 x 3 m along the gradient curve over the alignment,
/// one every 5 m from 2.5 m, collected.
fn blocks() -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let curve = Curve3::Elevated(Elevated3::new(
        Curve2::Intrinsic(alignment()),
        ElevationLaw::parabolic(10.0, 0.01, -0.01, 5000.0),
    ));
    let basis = b.push_value(curve).unwrap();
    let profile = b
        .push_value(Profile::Rectangle(RectangleProfile {
            x: 4.0,
            y: 2.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }))
        .unwrap();
    let block = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 3.0,
        }))
        .unwrap();
    let members: Vec<NodeId> = (0..1000)
        .map(|i| {
            let station = CurveStation::new(basis, Station::at(2.5 + 5.0 * Scalar::from(i)));
            b.push_value(InstanceAtStation::new(block, station.into()))
                .unwrap()
        })
        .collect();
    let root = b.push(GeometryNode::Collection(members)).unwrap();
    (b.finish(vec![root]).unwrap(), root)
}

fn placement(c: &mut Criterion) {
    let (graph, root) = blocks();
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(REQUEST)
        .expect("a chord error");
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    let (_, report) = compiler
        .compile_mesh_with_deviation(&graph, root, &options)
        .expect("the blocks compile");
    let frame = report
        .contributions
        .iter()
        .find(|c| c.path == DeviationPath::StationPlacement)
        .expect("the placements report their path");
    let DeviationBound::Certified(bound) = frame.bound else {
        panic!("a certified frame expected, got {report:?}");
    };
    assert!(report.meets_requested(), "{report:?}");
    eprintln!("placement: frame error {bound:e}");
    let mut group = c.benchmark_group("station/placement");
    group.sample_size(10);
    group.bench_function("1000 blocks at 1 mm", |b| {
        b.iter(|| {
            black_box(
                compiler
                    .compile_mesh_with_deviation(black_box(&graph), root, &options)
                    .expect("a mesh and its report"),
            )
        })
    });
    group.finish();
}

criterion_group!(benches, positions, run, placement);
criterion_main!(benches);
