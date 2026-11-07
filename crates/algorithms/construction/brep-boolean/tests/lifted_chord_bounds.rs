//! The trims of the general boolean with B-spline faces carry certified
//! chord bounds (#414, ADR 0077): every lifted pcurve of a result, cut at
//! its continuity breaks and split into eight spans each, is bounded by
//! `chord_bound2`, and no sampled point lies further from the chord of the
//! evaluated ends than the bound.
//!
//! The solids are those of `splines.rs`: a box with a bi-quadratic Bezier
//! roof against a pipe, a ball, another roofed box, and a bump through its
//! roof.

mod common;

use axiolid_brep::ExactBRep;
use axiolid_brep_boolean::boolean;
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_core::{BooleanOperator, Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::Curve2;
use axiolid_evaluate::bound::{chord_bound2, continuity_breaks2};
use axiolid_evaluate::curve::evaluate2;
use axiolid_overlay::ArcRing;
use common::{spline_box, spline_box_at, spline_box_sized};

const H: [[f64; 3]; 3] = [[1.0, 1.3, 1.1], [1.2, 1.8, 1.4], [0.9, 1.5, 1.2]];
const SAMPLES: usize = 400;
const SPANS: usize = 8;

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn pipe(x: f64, y: f64, r: f64, bottom: f64, top: f64) -> ExactBRep {
    let p = ArcPrism {
        section: ArcRing::circle(Point2::new(x, y), r),
        bottom,
        top,
    };
    boolean_arc_prisms_exact(&p, &p, BooleanOperator::Intersection, tol()).expect("a pipe")
}

/// The first solid, the second and a name, for each case.
fn cases() -> Vec<(&'static str, ExactBRep, ExactBRep)> {
    const G: [[f64; 3]; 3] = [[1.6, 1.0, 1.5], [1.1, 0.9, 1.3], [1.7, 1.2, 1.0]];
    const BUMP: [[f64; 3]; 3] = [[0.8, 0.8, 0.8], [0.8, 4.0, 0.8], [0.8, 0.8, 0.8]];
    vec![
        (
            "a pipe through the roof",
            spline_box(H),
            pipe(1.1, 0.9, 0.5, -1.0, 3.0),
        ),
        (
            "a ball dipping into the roof",
            spline_box(H),
            common::sphere(Point3::new(1.0, 1.1, 1.9), 0.6),
        ),
        (
            "roof to roof",
            spline_box(H),
            spline_box_at(Vec3::new(0.5, 0.3, -0.2), G),
        ),
        (
            "a bump through the roof",
            spline_box(H),
            spline_box_sized(Vec3::new(0.5, 0.5, 0.1), 0.5, BUMP),
        ),
    ]
}

/// The largest distance of the sampled points from the chord of the ends.
fn sampled(curve: &Curve2, lo: Scalar, hi: Scalar) -> Scalar {
    let point = |t: Scalar| evaluate2(curve, t).expect("evaluates");
    let (p, q) = (point(lo), point(hi));
    let chord = q - p;
    (1..SAMPLES)
        .map(|i| {
            let r = point(lo + (hi - lo) * i as Scalar / SAMPLES as Scalar) - p;
            let s = if chord.length_squared() == 0.0 {
                0.0
            } else {
                (r.dot(chord) / chord.length_squared()).clamp(0.0, 1.0)
            };
            (r - chord * s).length()
        })
        .fold(0.0, Scalar::max)
}

#[test]
fn every_lifted_trim_of_a_spline_boolean_has_a_certified_chord_bound() {
    let mut lifted = 0;
    for (name, a, b) in cases() {
        // The sections are the same for every operator; a union keeps the
        // outer parts of both operands, a difference the first's outside the
        // second and the second's inside the first.
        for op in [BooleanOperator::Union, BooleanOperator::Difference] {
            let result = boolean(&a, &b, op, tol()).unwrap_or_else(|e| panic!("{name}: {e}"));
            let topology = result.topology();
            for face in topology.faces() {
                for bound in &face.bounds {
                    let wire = &topology.loops()[bound.loop_id.index()];
                    for (k, edge_use) in wire.edges.iter().enumerate() {
                        let Some(handle) = edge_use.pcurve else {
                            continue;
                        };
                        let curve = &result.curves2()[handle.index()];
                        if !matches!(curve, Curve2::Lifted(_)) {
                            continue;
                        }
                        lifted += 1;
                        let span = result.pcurve_interval(bound.loop_id, k).expect("a span");
                        let (lo, hi) = (span.start.min(span.end), span.start.max(span.end));
                        let mut cuts = vec![lo];
                        cuts.extend(
                            continuity_breaks2(curve, 1)
                                .into_iter()
                                .filter(|&t| t > lo && t < hi),
                        );
                        cuts.push(hi);
                        for piece in cuts.windows(2) {
                            for i in 0..SPANS {
                                let at = |j: usize| {
                                    piece[0] + (piece[1] - piece[0]) * j as Scalar / SPANS as Scalar
                                };
                                let (a, b) = (at(i), at(i + 1));
                                let what = format!("{name} {op:?} [{a}, {b}]");
                                let bound = chord_bound2(curve, a, b)
                                    .unwrap_or_else(|| panic!("{what}: no bound"));
                                let sample = sampled(curve, a, b);
                                assert!(sample <= bound, "{what}: {sample} above {bound}");
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(lifted > 0, "no lifted pcurve in any result");
}
