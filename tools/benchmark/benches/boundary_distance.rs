//! Wall-clock cost of the exact boundary distance between building
//! elements that touch or cross (#273).
//!
//! Clash checks measure many pairs of elements whose boundaries meet: a
//! wall standing on a slab, two walls crossing. Their distance is zero, and
//! the search used to refine patch pairs round the contact until they were
//! within the accuracy (or its step budget ran out). A pair of elements 5 cm
//! apart is the control: the early end must leave it as it was.
//!
//! The elements are extrusions under one general rigid placement, compiled
//! exactly by `ReferenceExactCompiler`, so no coordinate is round. Each
//! case is validated before it is timed: a meeting pair must measure
//! `lower == 0` with `upper` within the accuracy, the apart pair an
//! interval holding its 5 cm gap.

use axiolid_brep::ExactBRep;
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Tolerance, Transform3, Vec3};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::boundary_distance;
use axiolid_mesh_compile::ReferenceExactCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, Instance, SolidOperation};
use axiolid_profile::{Profile, RectangleProfile};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::hint::black_box;

/// The accuracy a clash check asks for.
const ACCURACY: f64 = 1e-5;
/// The gap of the control pair.
const GAP: f64 = 0.05;

/// A `x` by `y` rectangle extruded by `depth` along `z`, placed by
/// `transform`, compiled exactly.
fn element(x: f64, y: f64, depth: f64, transform: Transform3) -> ExactBRep {
    let mut builder = GeometryGraphBuilder::new();
    let profile = builder
        .push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x,
                y,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )))
        .expect("a profile");
    let extrusion = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
        .expect("an extrusion");
    let root = builder
        .push(GeometryNode::Instance(Instance {
            source: extrusion,
            transform,
        }))
        .expect("an instance");
    let graph = builder.finish(vec![root]).expect("a valid graph");
    ReferenceExactCompiler::new()
        .compile_exact(&graph, root, &ExecutionOptions::new(Tolerance::METRE))
        .expect("an element compiles exactly")
}

/// The pairs: name, the two elements, and the true distance.
fn pairs() -> Vec<(&'static str, ExactBRep, ExactBRep, f64)> {
    // A storey placed off the axes, as a building's elements are.
    let storey = Transform3::from_translation(Vec3::new(12.0, -7.0, 3.0))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7);
    let at = |x: f64, y: f64, z: f64| storey * Transform3::from_translation(Vec3::new(x, y, z));
    // A 6 x 4 slab, 0.25 thick, its top at z = 0.
    let slab = || element(6.0, 4.0, 0.25, at(0.0, 0.0, -0.25));
    // A 4 m wall, 0.3 thick and 3 high, off the slab's centre.
    let wall = |z: f64| element(4.0, 0.3, 3.0, at(0.5, 0.7, z));
    let across = element(
        0.3,
        4.0,
        3.0,
        at(1.1, 0.4, 0.0) * Transform3::from_rotation_z(0.2),
    );
    vec![
        ("wall_on_slab", slab(), wall(0.0), 0.0),
        ("crossing_walls", wall(0.0), across, 0.0),
        ("wall_above_slab", slab(), wall(GAP), GAP),
    ]
}

fn touching(c: &mut Criterion) {
    let mut group = c.benchmark_group("exact/boundary_distance");
    group.sample_size(10);
    let tolerance = Tolerance::METRE;
    for (name, a, b, expected) in pairs() {
        let bounds = boundary_distance(&a, &b, ACCURACY, tolerance).expect("bounded");
        let width = (bounds.point_a - bounds.point_b).length();
        assert!(
            (width - bounds.upper).abs() <= 1e-12,
            "{name}: witnesses {width} apart, upper {}",
            bounds.upper
        );
        if expected == 0.0 {
            assert_eq!(bounds.lower, 0.0, "{name}: {bounds:?}");
            assert!(bounds.upper <= ACCURACY, "{name}: {bounds:?}");
        } else {
            assert!(
                bounds.lower <= expected + 1e-9 && expected <= bounds.upper + 1e-9,
                "{name}: {bounds:?} must hold {expected}"
            );
            assert!(
                bounds.upper - bounds.lower <= ACCURACY,
                "{name}: {bounds:?}"
            );
        }
        group.bench_with_input(BenchmarkId::new("pair", name), &name, |bench, _| {
            bench.iter(|| {
                black_box(
                    boundary_distance(black_box(&a), black_box(&b), ACCURACY, tolerance)
                        .expect("bounded"),
                )
            })
        });
    }
    group.finish();
}

criterion_group!(benches, touching);
criterion_main!(benches);
