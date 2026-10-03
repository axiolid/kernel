//! Wall-clock cost of cutting placed openings from a placed wall exactly
//! (#228).
//!
//! A 6 x 0.3 x 3 wall under a general rigid placement loses `n` rectangular
//! windows, each an extrusion placed across the wall, as nested differences
//! compiled by `ReferenceExactCompiler`. Every opening adds faces the next
//! boolean must consider, so this shows how the cost per opening grows with
//! the openings already cut.
//!
//! Each case is validated before it is timed: the exact volume must equal
//! the wall's minus the windows'.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Tolerance, Transform3, Vec3};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::exact_properties;
use axiolid_mesh_compile::ReferenceExactCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{Profile, RectangleProfile};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::hint::black_box;

const LENGTH: f64 = 6.0;
const THICKNESS: f64 = 0.3;
const HEIGHT: f64 = 3.0;
const WINDOW: (f64, f64) = (0.3, 0.6);

fn rectangle(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

fn placed_extrusion(
    builder: &mut GeometryGraphBuilder,
    profile: Profile,
    depth: f64,
    transform: Transform3,
) -> NodeId {
    let profile = builder
        .push(GeometryNode::Profile(profile))
        .expect("a profile");
    let extrusion = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
        .expect("an extrusion");
    builder
        .push(GeometryNode::Instance(Instance {
            source: extrusion,
            transform,
        }))
        .expect("an instance")
}

/// The wall with `count` windows along it, and the root of the last cut.
fn wall_with_windows(count: usize) -> (GeometryGraph, NodeId) {
    let placement = Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7);
    let mut builder = GeometryGraphBuilder::new();
    let mut body = placed_extrusion(
        &mut builder,
        rectangle(LENGTH, THICKNESS),
        HEIGHT,
        placement,
    );
    for index in 0..count {
        let x = -LENGTH / 2.0 + 0.3 + 0.55 * index as f64;
        let across = placement
            * Transform3::from_translation(Vec3::new(x, THICKNESS / 2.0 + 0.1, 1.5))
            * Transform3::from_rotation_x(std::f64::consts::FRAC_PI_2);
        let window = placed_extrusion(
            &mut builder,
            rectangle(WINDOW.0, WINDOW.1),
            THICKNESS + 0.2,
            across,
        );
        body = builder
            .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                left: body,
                right: window,
                operator: BooleanOperator::Difference,
            }))
            .expect("a difference");
    }
    let graph = builder.finish(vec![body]).expect("a valid graph");
    (graph, body)
}

fn placed_openings(c: &mut Criterion) {
    let mut group = c.benchmark_group("exact/placed_openings");
    group.sample_size(10);
    let options = ExecutionOptions::new(Tolerance::METRE);
    let compiler = ReferenceExactCompiler::new();
    for count in [1, 3, 10] {
        let (graph, root) = wall_with_windows(count);
        let solid = compiler
            .compile_exact(&graph, root, &options)
            .expect("placed openings compile exactly");
        let volume = exact_properties(&solid, Tolerance::METRE)
            .expect("measurable")
            .signed_volume;
        let expected = THICKNESS * (LENGTH * HEIGHT - count as f64 * WINDOW.0 * WINDOW.1);
        assert!(
            (volume - expected).abs() <= 1e-9 * expected,
            "{count} windows: volume {volume}, expected {expected}"
        );
        group.bench_with_input(BenchmarkId::new("windows", count), &count, |b, _| {
            b.iter(|| {
                black_box(
                    compiler
                        .compile_exact(black_box(&graph), root, &options)
                        .expect("exact"),
                )
            })
        });
    }
    group.finish();
}

criterion_group!(benches, placed_openings);
criterion_main!(benches);
