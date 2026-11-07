//! Wall-clock cost of the exact boundary distance between building
//! elements that touch or cross (#273).
//!
//! Clash checks measure many pairs of elements whose boundaries meet: a
//! wall standing on a slab, two walls crossing. Their distance is zero, and
//! the search used to refine patch pairs round the contact until they were
//! within the accuracy (or its step budget ran out). A pair of elements 5 cm
//! apart is the control: the early end must leave it as it was.
//!
//! Every pair is also measured at georeferenced coordinates (`_far`: the
//! storey at (450012, 5599993, 3)), where the rounding of the margins is
//! about `1e-9` and the old `1e-12` margins left an interval `1e-5` wide
//! (#318), and two regular `n`-gon prisms 1.1 m apart (`ngon_prisms_n`, at
//! the origin and `_far`), built as contours of line segments in world
//! coordinates, whose nearest points are two facing corners.
//!
//! The elements are extrusions under one general rigid placement, compiled
//! exactly by `ReferenceExactCompiler`, so no coordinate is round. Each
//! case is validated before it is timed: a meeting pair must measure
//! `lower == 0` with `upper` within the accuracy, the apart pair an
//! interval holding its 5 cm gap, and far from the origin up to the rounding
//! of the witnesses' coordinates.

use axiolid_brep::ExactBRep;
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Interval, Point2, Tolerance, Transform3, Vec3};
use axiolid_curve::{Curve2, Line2};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::boundary_distance;
use axiolid_mesh_compile::ReferenceExactCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, Instance, SolidOperation};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::hint::black_box;

/// The accuracy a clash check asks for.
const ACCURACY: f64 = 1e-5;
/// The gap of the control pair.
const GAP: f64 = 0.05;

/// The accuracy of the `n`-gon prisms.
const NGON_ACCURACY: f64 = 1e-6;

/// A georeferenced storey origin (#318).
const FAR_STOREY: Vec3 = Vec3::new(450_012.0, 5_599_993.0, 3.0);

/// One case: name, the two elements, the true distance, the accuracy asked.
type Case = (String, ExactBRep, ExactBRep, f64, f64);

/// A `x` by `y` rectangle extruded by `depth` along `z`, placed by
/// `transform`, compiled exactly.
fn element(x: f64, y: f64, depth: f64, transform: Transform3) -> ExactBRep {
    compile(
        Profile::Rectangle(RectangleProfile {
            x,
            y,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }),
        depth,
        transform,
    )
}

/// `profile` extruded by `depth` along `z`, placed by `transform`, compiled
/// exactly.
fn compile(profile: Profile, depth: f64, transform: Transform3) -> ExactBRep {
    let mut builder = GeometryGraphBuilder::new();
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

/// The pairs of a storey placed at `origin` off the axes, as a building's
/// elements are, with `suffix` after their names.
fn storey_pairs(origin: Vec3, suffix: &str) -> Vec<Case> {
    let storey = Transform3::from_translation(origin)
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
        (
            format!("wall_on_slab{suffix}"),
            slab(),
            wall(0.0),
            0.0,
            ACCURACY,
        ),
        (
            format!("crossing_walls{suffix}"),
            wall(0.0),
            across,
            0.0,
            ACCURACY,
        ),
        (
            format!("wall_above_slab{suffix}"),
            slab(),
            wall(GAP),
            GAP,
            ACCURACY,
        ),
    ]
}

/// A regular `n`-gon of circumradius 0.2 about `(cx, cy)`, counter-clockwise
/// from angle `start`; the first corner is `first` when given.
fn ngon(n: usize, cx: f64, cy: f64, start: f64, first: Option<Point2>) -> Vec<Point2> {
    (0..n)
        .map(|k| {
            if let (0, Some(p)) = (k, first) {
                return p;
            }
            let theta = start + core::f64::consts::TAU * k as f64 / n as f64;
            Point2::new(cx + 0.2 * theta.cos(), cy + 0.2 * theta.sin())
        })
        .collect()
}

/// A 3 m prism over `corners`, in world coordinates, as a contour of line
/// segments under the identity placement.
fn prism(corners: &[Point2]) -> ExactBRep {
    let segments = (0..corners.len())
        .map(|k| {
            let (from, to) = (corners[k], corners[(k + 1) % corners.len()]);
            ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: from,
                    direction: to - from,
                }),
                domain: Interval {
                    start: 0.0,
                    end: 1.0,
                },
                same_sense: true,
            }
        })
        .collect();
    let profile = Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    });
    compile(profile, 3.0, Transform3::IDENTITY)
}

/// Two `n`-gon prisms 1.1 m apart about `(x, y)`: b's first corner stands
/// exactly 0.7 from a's, which are the nearest points.
fn ngon_pair(n: usize, x: f64, y: f64, suffix: &str) -> Case {
    let a = ngon(n, x, y, 0.0, None);
    let b = ngon(
        n,
        x + 1.1,
        y,
        core::f64::consts::PI,
        Some(Point2::new(x + 1.1 - 0.2, y)),
    );
    let truth = b[0].x - a[0].x;
    (
        format!("ngon_prisms_{n}{suffix}"),
        prism(&a),
        prism(&b),
        truth,
        NGON_ACCURACY,
    )
}

/// The cases: the near and far storeys, and the `n`-gon prisms.
fn pairs() -> Vec<Case> {
    let mut cases = storey_pairs(Vec3::new(12.0, -7.0, 3.0), "");
    cases.extend(storey_pairs(FAR_STOREY, "_far"));
    for n in [8, 128] {
        cases.push(ngon_pair(n, 0.0, 0.0, ""));
        cases.push(ngon_pair(n, 450_000.0, 5_600_000.0, "_far"));
    }
    cases
}

fn touching(c: &mut Criterion) {
    let mut group = c.benchmark_group("exact/boundary_distance");
    group.sample_size(10);
    let tolerance = Tolerance::METRE;
    for (name, a, b, expected, accuracy) in pairs() {
        let bounds = boundary_distance(&a, &b, accuracy, tolerance).expect("bounded");
        let width = (bounds.point_a - bounds.point_b).length();
        assert!(
            (width - bounds.upper).abs() <= 1e-12,
            "{name}: witnesses {width} apart, upper {}",
            bounds.upper
        );
        // The witnesses carry the rounding of their coordinates.
        let slack =
            1e-9 + 32.0 * f64::EPSILON * (bounds.point_a.length() + bounds.point_b.length());
        if expected == 0.0 {
            assert_eq!(bounds.lower, 0.0, "{name}: {bounds:?}");
            assert!(bounds.upper <= accuracy, "{name}: {bounds:?}");
        } else {
            assert!(
                bounds.lower <= expected + slack && expected <= bounds.upper + slack,
                "{name}: {bounds:?} must hold {expected}"
            );
            assert!(
                bounds.upper - bounds.lower <= accuracy,
                "{name}: {bounds:?}"
            );
        }
        group.bench_with_input(BenchmarkId::new("pair", &name), &name, |bench, _| {
            bench.iter(|| {
                black_box(
                    boundary_distance(black_box(&a), black_box(&b), accuracy, tolerance)
                        .expect("bounded"),
                )
            })
        });
    }
    group.finish();
}

criterion_group!(benches, touching);
criterion_main!(benches);
