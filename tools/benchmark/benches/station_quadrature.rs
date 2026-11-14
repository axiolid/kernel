//! Wall-clock cost of compiling a run of sections along a frame read by
//! quadrature (#406).
//!
//! The fold check of a stationed run decides every span of a station
//! interval from the frame law's certified rates, and a frame read by
//! quadrature (an ellipse, a B-spline, an offset beside them) pays for a
//! certified second derivative and an arc-length inverse per read. Three
//! cases, each a 1 m by 1 m rectangle compiled with its deviation report:
//!
//! - `ellipse`: stations 0 and 10 along the 20 by 10 ellipse, 1.5 m to the
//!   side, at a 10 mm chord error: certified;
//! - `offset across a knot, near`: stations 1 and 22 along the offset by
//!   1 m of a cubic B-spline (one interior knot, where the offset's
//!   curvature jumps), 0.3 m to the side, at 50 mm: certified, the spans
//!   across the knot left to the conservative bound;
//! - `offset across a knot, far`: the same 1 m to the side: reported
//!   unbounded by name, every span across the knot bisected to the
//!   fallback depth.
//!
//! Each case is validated before it is timed: the outcome is the one named.
//!
//! Measured on a shared machine under a load of 18 (best of three wall-clock
//! runs of the same compile, a measurement, not a performance claim): the
//! ellipse 10 ms (18 ms before the fold check kept the reads at the ends of
//! its spans); near the knot 20.6 s against 46.4 s, far from it 11.7 s
//! against 52.0 s. In the near case the fold decision alone is 6.5 s against
//! 47 s (83 spans with rates at about 50 ms each); the rest is the mesh's
//! own cells.
//!
//! ```bash
//! cargo bench -p axiolid-benchmark --bench station_quadrature
//! ```

use std::hint::black_box;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame2, Point2, Scalar, Tolerance, Vec2};
use axiolid_curve::{BSplineCurve2, Curve2, Ellipse2, KnotSpec};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, ReferenceMeshCompiler};
use axiolid_model::{
    CurveRelation, GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, SectionAtStation,
    SolidOperation, Station, StationFrame, StationOffsets,
};
use axiolid_profile::{Profile, RectangleProfile};
use criterion::{criterion_group, criterion_main, Criterion};

/// What a case must report.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Expect {
    Certified,
    Unbounded,
}

fn ellipse() -> Curve2 {
    Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: 20.0,
        semi_axis_y: 10.0,
    })
}

/// A cubic plan B-spline with one interior knot, about 25 long.
fn cubic() -> Curve2 {
    Curve2::BSpline(BSplineCurve2 {
        degree: 3,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(6.0, 1.0),
            Point2::new(12.0, -2.0),
            Point2::new(18.0, 3.0),
            Point2::new(24.0, 2.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![4, 1, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// A rectangle run along `curve` (or its offset by 1 m) between `at`,
/// `lateral` to the side.
fn run(curve: Curve2, offset: bool, lateral: Scalar, at: [Scalar; 2]) -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let mut basis = b.push_value(curve).unwrap();
    if offset {
        basis = b
            .push(GeometryNode::CurveRelation(CurveRelation::Offset {
                basis,
                distance: 1.0,
                reference_direction: None,
            }))
            .unwrap();
    }
    let profile = b
        .push_value(Profile::Rectangle(RectangleProfile {
            x: 1.0,
            y: 1.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }))
        .unwrap();
    let offsets = StationOffsets::new(lateral, 0.0, 0.0);
    let root = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix: basis,
                sections: at
                    .iter()
                    .map(|&d| SectionAtStation::new(profile, Station::new(d, offsets)))
                    .collect(),
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    (b.finish(vec![root]).unwrap(), root)
}

fn compile(
    compiler: &ReferenceMeshCompiler<BoolmeshBoolean>,
    graph: &GeometryGraph,
    root: NodeId,
    options: &ExecutionOptions,
) -> DeviationBound {
    let (_, report) = compiler
        .compile_mesh_with_deviation(graph, root, options)
        .expect("the run compiles");
    report
        .contributions
        .iter()
        .find(|c| c.path == DeviationPath::StationedSpine)
        .expect("the run reports its path")
        .bound
}

fn quadrature(c: &mut Criterion) {
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    let cases = [
        (
            "ellipse",
            run(ellipse(), false, 1.5, [0.0, 10.0]),
            0.01,
            Expect::Certified,
        ),
        (
            "offset across a knot, near",
            run(cubic(), true, 0.3, [1.0, 22.0]),
            0.05,
            Expect::Certified,
        ),
        (
            "offset across a knot, far",
            run(cubic(), true, 1.0, [1.0, 22.0]),
            0.05,
            Expect::Unbounded,
        ),
    ];
    let mut group = c.benchmark_group("station/quadrature");
    group.sample_size(10);
    for (name, (graph, root), chord, expect) in cases {
        let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
            .with_chord_error(chord)
            .expect("a chord error");
        let bound = compile(&compiler, &graph, root, &options);
        let ok = match bound {
            DeviationBound::Certified(_) => expect == Expect::Certified,
            DeviationBound::Unbounded(_) => expect == Expect::Unbounded,
            _ => false,
        };
        assert!(ok, "{name}: unexpected {bound:?}");
        group.bench_function(name, |b| {
            b.iter(|| black_box(compile(&compiler, black_box(&graph), root, &options)))
        });
    }
    group.finish();
}

criterion_group!(benches, quadrature);
criterion_main!(benches);
