//! Geometric consistency auditing for exact boundary representations.
//!
//! # Why this exists next to `audit_brep`
//!
//! [`audit_brep`](axiolid_topology::audit_brep) is pure topology: handles and
//! adjacency, no coordinates and no tolerance. That makes it exact and
//! reproducible, and it is deliberately kept that way.
//!
//! It also means a whole class of defect passes it silently. A pcurve that
//! has nothing to do with the 3D edge it trims still resolves, still closes
//! its loop, and still balances its edge uses. Every variant of
//! `ExactBRepError` is likewise about PRESENCE -- is there a pcurve, does the
//! handle resolve -- never about AGREEMENT.
//!
//! That gap is not hypothetical. An earlier arc-aware overlay built cap loops
//! whose pcurves were straight chords across arc edges. The solid closed,
//! validated, audited clean, and reported a plausible area, while the cap
//! boundary disagreed with the wall boundary along the same edge.
//!
//! This module closes that gap by EVALUATING. Each check maps parameters to
//! points through the real evaluators and compares positions, so it is
//! necessarily tolerance-dependent -- which is exactly why it is separate
//! from the exact topological audit rather than folded into it.

#![forbid(unsafe_code)]

mod report;

pub use report::{GeometricDefect, GeometricHealth};

use axiolid_brep::ExactBRep;
use axiolid_core::{Point3, Scalar, Tolerance};
use axiolid_evaluate::{curve, surface};
use axiolid_topology::Orientation;

/// Audit the geometric consistency of an exact B-rep.
///
/// Complements the topological audit: this one evaluates curves and surfaces
/// and compares positions, so it needs a tolerance and can only report
/// agreement to within it.
///
/// Checks performed:
///
/// - every edge's start and end vertex lies on the edge's own 3D curve;
/// - every pcurve, mapped through its face's surface, lands on the 3D curve
///   of the edge it trims.
#[must_use]
pub fn geometric_audit(brep: &ExactBRep, tolerance: Tolerance) -> GeometricHealth {
    let mut health = GeometricHealth::default();
    check_vertices_on_curves(brep, tolerance, &mut health);
    check_pcurves_against_curves(brep, tolerance, &mut health);
    health
}

/// Sample parameters used to compare a pcurve against its 3D curve.
///
/// Endpoints catch a pcurve that trims the wrong portion; the interior
/// samples catch one that shares endpoints but takes a different path --
/// a straight chord standing in for an arc, for instance, which agrees
/// exactly at both ends.
const SAMPLES: [Scalar; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];

/// The worst deviation of a pcurve from its edge when the two share no
/// parameterisation: each lifted sample is projected onto the edge (by
/// inversion), and must land within the edge's span, at the use's start and
/// end for the end samples, and in order along the use between them. A
/// sample that fails any of these counts as infinitely far.
fn follows_in_order(
    pcurve: &axiolid_curve::Curve2,
    pcurve_interval: axiolid_core::Interval,
    support: &axiolid_surface::Surface,
    curve3: &axiolid_curve::Curve3,
    edge_interval: axiolid_core::Interval,
    orientation: Orientation,
    tolerance: Tolerance,
) -> Scalar {
    const DENSE: usize = 16;
    let (lo, hi) = (
        edge_interval.start.min(edge_interval.end),
        edge_interval.start.max(edge_interval.end),
    );
    let periodic = matches!(
        curve3,
        axiolid_curve::Curve3::Circle(_) | axiolid_curve::Curve3::Ellipse(_)
    );
    // An edge that closes on itself (one vertex) reads its end as its start.
    let closed = match (curve::evaluate3(curve3, lo), curve::evaluate3(curve3, hi)) {
        (Ok(a), Ok(b)) => distance(a, b) <= tolerance.linear(),
        _ => false,
    };
    let slack = 1e-9 * (1.0 + lo.abs().max(hi.abs()));
    let mut worst: Scalar = 0.0;
    let mut previous: Option<Scalar> = None;
    for i in 0..=DENSE {
        let sample = i as Scalar / DENSE as Scalar;
        let Ok(uv) = curve::evaluate2(pcurve, lerp(pcurve_interval, sample)) else {
            return Scalar::INFINITY;
        };
        let Ok(lifted) = surface::evaluate(support, uv.x, uv.y) else {
            return Scalar::INFINITY;
        };
        let Ok(t) = curve::invert3(curve3, lifted, tolerance) else {
            return Scalar::INFINITY;
        };
        // Into the edge's span, by whole turns for a closed conic.
        let candidates: &[Scalar] = if periodic {
            &[0.0, 1.0, -1.0, 2.0]
        } else {
            &[0.0]
        };
        let Some(t) = candidates
            .iter()
            .map(|k| t + k * core::f64::consts::TAU)
            .find(|t| *t >= lo - slack && *t <= hi + slack)
        else {
            return Scalar::INFINITY;
        };
        // Where along the use, from 0 at its start to 1 at its end.
        let along = (t - edge_interval.start) / (edge_interval.end - edge_interval.start);
        let along = match orientation {
            Orientation::Forward => along,
            Orientation::Reversed => 1.0 - along,
        };
        // The end samples at the use's ends (a closed edge's end may read
        // as its start).
        let at_end = |target: Scalar| {
            (along - target).abs() <= 1e-6 || (closed && (along - (1.0 - target)).abs() <= 1e-6)
        };
        if (i == 0 && !at_end(0.0)) || (i == DENSE && !at_end(1.0)) {
            return Scalar::INFINITY;
        }
        if i > 0 && i < DENSE {
            if let Some(last) = previous {
                if along < last - 1e-9 {
                    return Scalar::INFINITY;
                }
            }
            previous = Some(along);
        } else if i == 0 {
            previous = Some(0.0);
        }
        let Ok(on) = curve::evaluate3(curve3, t) else {
            return Scalar::INFINITY;
        };
        worst = worst.max(distance(lifted, on));
    }
    worst
}

fn lerp(interval: axiolid_core::Interval, t: Scalar) -> Scalar {
    interval.start + (interval.end - interval.start) * t
}

fn distance(a: Point3, b: Point3) -> Scalar {
    (a - b).length()
}

/// Every edge's endpoints must lie on the curve that edge claims to follow.
fn check_vertices_on_curves(brep: &ExactBRep, tolerance: Tolerance, health: &mut GeometricHealth) {
    let topology = brep.topology();
    for (index, edge) in topology.edges().iter().enumerate() {
        let Some(curve_id) = edge.curve else { continue };
        let Some(curve3) = brep.curves3().get(curve_id.index()) else {
            continue;
        };
        let Some(edge_id) = topology.edge_id_at(index) else {
            continue;
        };
        let Some(interval) = brep.edge_interval(edge_id) else {
            continue;
        };
        let vertices = topology.vertices();
        let (Some(start), Some(end)) = (
            vertices.get(edge.start.index()),
            vertices.get(edge.end.index()),
        ) else {
            continue;
        };
        for (parameter, vertex) in [(interval.start, start), (interval.end, end)] {
            let Ok(point) = curve::evaluate3(curve3, parameter) else {
                health.push(GeometricDefect::UnevaluableCurve { edge: index });
                continue;
            };
            let error = distance(point, vertex.position);
            if error > tolerance.linear() {
                health.push(GeometricDefect::VertexOffCurve { edge: index, error });
            }
        }
    }
}

/// Every pcurve, lifted through its face surface, must follow the 3D edge.
fn check_pcurves_against_curves(
    brep: &ExactBRep,
    tolerance: Tolerance,
    health: &mut GeometricHealth,
) {
    let topology = brep.topology();
    for face in topology.faces() {
        let Some(surface_id) = face.surface else {
            continue;
        };
        let Some(support) = brep.surfaces().get(surface_id.index()) else {
            continue;
        };
        for bound in &face.bounds {
            let Some(lp) = topology.loops().get(bound.loop_id.index()) else {
                continue;
            };
            for (use_index, use_) in lp.edges.iter().enumerate() {
                let Some(pcurve_id) = use_.pcurve else {
                    continue;
                };
                let Some(pcurve) = brep.curves2().get(pcurve_id.index()) else {
                    continue;
                };
                let Some(pcurve_interval) = brep.pcurve_interval(bound.loop_id, use_index) else {
                    continue;
                };
                let Some(edge) = topology.edges().get(use_.edge.index()) else {
                    continue;
                };
                let Some(curve_id) = edge.curve else { continue };
                let Some(curve3) = brep.curves3().get(curve_id.index()) else {
                    continue;
                };
                let Some(edge_interval) = brep.edge_interval(use_.edge) else {
                    continue;
                };

                let mut worst: Scalar = 0.0;
                for sample in SAMPLES {
                    // The pcurve runs in loop-traversal order; the 3D curve
                    // runs from its own start vertex. A reversed use walks
                    // the edge backwards, so the sample must be mirrored or
                    // every reversed edge would look like a mismatch.
                    let along = match use_.orientation {
                        Orientation::Forward => sample,
                        Orientation::Reversed => 1.0 - sample,
                    };
                    let parametric = lerp(pcurve_interval, sample);
                    let spatial = lerp(edge_interval, along);

                    let Ok(uv) = curve::evaluate2(pcurve, parametric) else {
                        health.push(GeometricDefect::UnevaluablePcurve {
                            loop_id: bound.loop_id.index(),
                            use_index,
                        });
                        break;
                    };
                    let Ok(lifted) = surface::evaluate(support, uv.x, uv.y) else {
                        health.push(GeometricDefect::UnevaluableSurface {
                            loop_id: bound.loop_id.index(),
                            use_index,
                        });
                        break;
                    };
                    let Ok(expected) = curve::evaluate3(curve3, spatial) else {
                        health.push(GeometricDefect::UnevaluableCurve {
                            edge: use_.edge.index(),
                        });
                        break;
                    };
                    // Track the WORST sample rather than reporting the first
                    // failure: the first sample over tolerance is rarely the
                    // largest deviation, and an under-reported error invites
                    // someone to widen the tolerance just past it.
                    let error = distance(lifted, expected);
                    if error > worst {
                        worst = error;
                    }
                }
                // An implicit pcurve (ADR 0077) is parameterised by its own
                // cells, not in proportion to the edge: it follows the edge
                // when every lifted sample lies on the edge's span, in the
                // order the use runs, from the use's start to its end.
                if worst > tolerance.linear()
                    && matches!(pcurve, axiolid_curve::Curve2::Implicit(_))
                {
                    worst = follows_in_order(
                        pcurve,
                        pcurve_interval,
                        support,
                        curve3,
                        edge_interval,
                        use_.orientation,
                        tolerance,
                    );
                }
                if worst > tolerance.linear() {
                    health.push(GeometricDefect::PcurveOffCurve {
                        loop_id: bound.loop_id.index(),
                        use_index,
                        error: worst,
                    });
                }
            }
        }
    }
}
