//! Wall-clock cost of a boolean's measured deviation report (#314).
//!
//! Two walls at a 1 mm request, each compiled with and without its
//! deviation report, so the difference is the cost of measuring the mesh
//! against the exact result:
//!
//! - **planar**: the wall `3.65 x 0.25 x 3.67` less two door openings
//!   `1.01 x 2.26` (#312, #314). Every face is a plane, so the search runs
//!   to the measured distance, near rounding.
//! - **curved**: a wall `6 x 0.3 x 3` less a round window of radius `0.4`
//!   under a building placement (#235). Its cylinder faces keep the rule
//!   that stops once the bound is within the request.
//!
//! Each case is validated before it is timed: the report must be one
//! certified boolean bound within the request, and the planar one within
//! `1e-11`.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Tolerance, Transform2, Transform3, Vec2, Vec3};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::hint::black_box;

const REQUEST: f64 = 1e-3;

fn rectangle(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

struct Graph(GeometryGraphBuilder);

impl Graph {
    fn push(&mut self, node: GeometryNode) -> NodeId {
        self.0.push(node).expect("a valid node")
    }

    fn extrusion(&mut self, profile: Profile, depth: f64) -> NodeId {
        let profile = self.push(GeometryNode::Profile(profile));
        self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
    }

    fn place(&mut self, source: NodeId, transform: Transform3) -> NodeId {
        self.push(GeometryNode::Instance(Instance { source, transform }))
    }

    fn difference(&mut self, left: NodeId, right: NodeId) -> NodeId {
        self.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator: BooleanOperator::Difference,
        }))
    }

    fn finish(self, root: NodeId) -> (GeometryGraph, NodeId) {
        (self.0.finish(vec![root]).expect("a valid graph"), root)
    }
}

/// The wall less two doors, as #312 builds it.
fn planar() -> (GeometryGraph, NodeId) {
    let mut g = Graph(GeometryGraphBuilder::new());
    let derived = |x: f64, y: f64, transform: Transform2| Profile::Derived {
        basis: Box::new(rectangle(x, y)),
        transform,
    };
    let wall = derived(
        3.65,
        0.25,
        Transform2::from_cols(Vec2::NEG_X, Vec2::NEG_Y, Vec2::new(1.825, 0.0)),
    );
    let mut net = g.extrusion(wall, 3.67);
    for x in [3.325, 1.35] {
        let door = derived(
            1.01,
            2.26,
            Transform2::from_cols(Vec2::Y, Vec2::NEG_X, Vec2::new(1.13, 0.505)),
        );
        let door = g.extrusion(door, 0.25);
        let door = g.place(
            door,
            Transform3::from_cols(Vec3::Z, Vec3::NEG_X, Vec3::NEG_Y, Vec3::new(x, 0.125, 0.2)),
        );
        net = g.difference(net, door);
    }
    g.finish(net)
}

/// The wall with a round window, under a building placement.
fn curved() -> (GeometryGraph, NodeId) {
    let building =
        Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2)) * Transform3::from_rotation_z(0.6);
    let mut g = Graph(GeometryGraphBuilder::new());
    let wall = g.extrusion(rectangle(6.0, 0.3), 3.0);
    let wall = g.place(wall, building);
    let window = g.extrusion(
        Profile::Circle(CircleProfile {
            radius: 0.4,
            thickness: None,
        }),
        1.0,
    );
    let across = Transform3::from_translation(Vec3::new(0.8, 0.5, 1.5))
        * Transform3::from_rotation_x(std::f64::consts::FRAC_PI_2);
    let window = g.place(window, building * across);
    let cut = g.difference(wall, window);
    g.finish(cut)
}

fn deviation_report(c: &mut Criterion) {
    let mut group = c.benchmark_group("deviation/boolean");
    group.sample_size(10);
    let options = ExecutionOptions::new(Tolerance::new(REQUEST, 1e-9).expect("a tolerance"))
        .with_chord_error(REQUEST)
        .expect("a chord error");
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    for (name, (graph, root), ceiling) in
        [("planar", planar(), 1e-11), ("curved", curved(), REQUEST)]
    {
        let (_, report) = compiler
            .compile_mesh_with_deviation(&graph, root, &options)
            .expect("the boolean compiles");
        let [contribution] = report.contributions.as_slice() else {
            panic!("{name}: one contribution expected, got {report:?}");
        };
        assert_eq!(contribution.path, DeviationPath::Boolean, "{name}");
        let DeviationBound::Certified(bound) = contribution.bound else {
            panic!("{name}: a certified bound expected, got {report:?}");
        };
        assert!(
            bound <= ceiling,
            "{name}: bound {bound:e} above {ceiling:e}"
        );
        eprintln!("{name}: bound {bound:e}");
        group.bench_with_input(BenchmarkId::new("mesh", name), &root, |b, &root| {
            b.iter(|| {
                black_box(
                    compiler
                        .compile_mesh_reported(black_box(&graph), root, &options)
                        .expect("a mesh"),
                )
            })
        });
        group.bench_with_input(
            BenchmarkId::new("mesh_and_report", name),
            &root,
            |b, &root| {
                b.iter(|| {
                    black_box(
                        compiler
                            .compile_mesh_with_deviation(black_box(&graph), root, &options)
                            .expect("a mesh and its report"),
                    )
                })
            },
        );
    }
    group.finish();
}

criterion_group!(benches, deviation_report);
criterion_main!(benches);
