// SPDX-License-Identifier: MPL-2.0

//! Behavioural tests for constrained Delaunay triangulation.
//!
//! These assert the properties the crate promises, not a fixed triangle list:
//! a different but equally valid triangulation must not fail the suite.

use axiolid_core::Point2;
use axiolid_triangulate::{triangulate, triangulate_refined, Constraint, Quality, RefineOutcome};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

/// Every triangle must be counter-clockwise and non-degenerate.
fn assert_all_ccw(tri: &axiolid_triangulate::Triangulation) {
    for t in 0..tri.triangle_count() {
        let a = tri.points()[tri.triangles()[3 * t] as usize];
        let b = tri.points()[tri.triangles()[3 * t + 1] as usize];
        let c = tri.points()[tri.triangles()[3 * t + 2] as usize];
        let cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        assert!(
            cross > 0.0,
            "triangle {t} is not counter-clockwise: {cross}"
        );
    }
}

/// Smallest interior angle across the whole mesh, in degrees.
fn worst_angle(tri: &axiolid_triangulate::Triangulation) -> f64 {
    let mut worst: f64 = 180.0;
    for t in 0..tri.triangle_count() {
        let a = tri.points()[tri.triangles()[3 * t] as usize];
        let b = tri.points()[tri.triangles()[3 * t + 1] as usize];
        let c = tri.points()[tri.triangles()[3 * t + 2] as usize];
        for (u, v, w) in [(a, b, c), (b, c, a), (c, a, b)] {
            let ux = v.x - u.x;
            let uy = v.y - u.y;
            let vx = w.x - u.x;
            let vy = w.y - u.y;
            let dot = ux * vx + uy * vy;
            let nu = (ux * ux + uy * uy).sqrt();
            let nv = (vx * vx + vy * vy).sqrt();
            if nu == 0.0 || nv == 0.0 {
                continue;
            }
            let angle = (dot / (nu * nv)).clamp(-1.0, 1.0).acos().to_degrees();
            worst = worst.min(angle);
        }
    }
    worst
}

#[test]
fn a_square_triangulates_into_two_triangles() {
    let pts = [p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)];
    let tri = triangulate(&pts, &[]).expect("square is triangulable");
    assert_eq!(tri.triangle_count(), 2);
    assert_all_ccw(&tri);
}

#[test]
fn fewer_than_three_points_is_an_error() {
    let pts = [p(0.0, 0.0), p(1.0, 0.0)];
    assert!(triangulate(&pts, &[]).is_err());
}

#[test]
fn collinear_input_is_reported_not_silently_empty() {
    let pts = [p(0.0, 0.0), p(1.0, 0.0), p(2.0, 0.0), p(3.0, 0.0)];
    let err = triangulate(&pts, &[]).expect_err("collinear input has no triangulation");
    assert_eq!(
        err,
        axiolid_triangulate::TriangulationError::AllPointsCollinear
    );
}

#[test]
fn an_out_of_range_constraint_is_rejected() {
    let pts = [p(0.0, 0.0), p(1.0, 0.0), p(0.0, 1.0)];
    let err = triangulate(&pts, &[Constraint::new(0, 99)]).expect_err("index 99 does not exist");
    assert!(matches!(
        err,
        axiolid_triangulate::TriangulationError::ConstraintOutOfRange { index: 99 }
    ));
}

#[test]
fn a_constraint_edge_survives_triangulation() {
    // Five points where the unconstrained Delaunay triangulation would not
    // choose the 0-3 diagonal; the constraint forces it.
    let pts = [
        p(0.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 4.0),
        p(0.0, 4.0),
        p(2.0, 2.05),
    ];
    let c = Constraint::new(0, 2);
    let tri = triangulate(&pts, &[c]).expect("triangulable");
    let present = (0..tri.triangle_count()).any(|t| {
        let v = &tri.triangles()[3 * t..3 * t + 3];
        v.contains(&0) && v.contains(&2)
    });
    assert!(present, "constraint edge (0, 2) was lost");
    assert_all_ccw(&tri);
}

#[test]
fn constraint_direction_does_not_matter() {
    // (a, b) and (b, a) must be the same constraint; a directional compare
    // is a classic source of "constraint lost" bugs.
    assert_eq!(Constraint::new(7, 2), Constraint::new(2, 7));
}

/// Refinement improves mesh quality where it can, and reports honestly where
/// it cannot.
///
/// The original version of this test asserted that the worst angle always
/// improves. That assertion was wrong, and the implementation proved it: in
/// this fixture the worst triangle is `(9.9, 1.4) - (10.0, 3.0) - (9.9, 1.6)`,
/// whose three corners are all INPUT vertices. It sits in the corridor
/// between the opening and the boundary, and interior Steiner insertion
/// cannot open it -- only moving an input point or splitting the boundary
/// would, and both are forbidden. Ruppert's termination guarantee explicitly
/// excludes small input angles for exactly this reason.
///
/// So the property worth asserting is not "always improves" but "never
/// claims more than it delivered".
#[test]
fn refinement_reports_honestly_on_an_input_pinned_sliver() {
    let pts = [
        p(0.0, 0.0),
        p(10.0, 0.0),
        p(10.0, 3.0),
        p(0.0, 3.0),
        // A thin opening pushed against the right edge.
        p(9.7, 1.4),
        p(9.9, 1.4),
        p(9.9, 1.6),
        p(9.7, 1.6),
    ];
    let quality = Quality {
        min_angle_degrees: 20.0,
        max_steiner_points: 64,
    };
    let (after, outcome) = triangulate_refined(&pts, &[], quality).expect("triangulable");
    let achieved = worst_angle(&after);

    assert_all_ccw(&after);
    match outcome {
        RefineOutcome::Achieved { .. } => {
            assert!(
                achieved >= 20.0 - 1e-9,
                "reported Achieved but the worst angle is {achieved}"
            );
        }
        RefineOutcome::Capped {
            worst_angle_degrees,
            inserted,
        } => {
            assert!(
                worst_angle_degrees < 20.0,
                "reported Capped but the bound was met"
            );
            // The reported figure must be the real one, not a stale estimate.
            assert!(
                (worst_angle_degrees - achieved).abs() < 1e-9,
                "reported {worst_angle_degrees} but the mesh's worst angle is {achieved}"
            );
            assert!(inserted > 0, "capped without attempting any insertion");
        }
        other => panic!("unexpected refine outcome: {other:?}"),
    }
}

/// Refinement must improve a sliver that interior insertion CAN fix.
///
/// The fixture matters more than it looks. An earlier version used a point
/// at `(5.0, 0.15)` -- close to the bottom edge -- which produces a triangle
/// whose small angle is pinned between two INPUT vertices and the boundary,
/// exactly the unfixable case above. A genuinely fixable sliver needs its
/// bad angle to be interior, so the circumcentre lands inside the hull and
/// the insertion is legal.
#[test]
fn refinement_improves_a_fixable_interior_sliver() {
    let pts = [
        p(0.0, 0.0),
        p(12.0, 0.0),
        p(12.0, 12.0),
        p(0.0, 12.0),
        // Two interior points close together, far from every boundary: the
        // thin triangle they form has room around it for Steiner insertion.
        p(5.8, 6.0),
        p(6.2, 6.05),
    ];
    let before = triangulate(&pts, &[]).expect("triangulable");
    let rough = worst_angle(&before);

    let quality = Quality {
        min_angle_degrees: 20.0,
        max_steiner_points: 256,
    };
    let (after, outcome) = triangulate_refined(&pts, &[], quality).expect("triangulable");
    let improved = worst_angle(&after);

    assert_all_ccw(&after);
    assert!(
        improved > rough,
        "interior refinement did not improve a fixable sliver: {rough} -> {improved}"
    );
    // And whatever it reports must match the mesh it produced.
    if let RefineOutcome::Capped {
        worst_angle_degrees,
        ..
    } = outcome
    {
        assert!((worst_angle_degrees - improved).abs() < 1e-9);
    }
}

/// A capped run must say so rather than pretending it met the bound.
#[test]
fn an_impossible_bound_reports_capped_not_achieved() {
    let pts = [p(0.0, 0.0), p(10.0, 0.0), p(10.0, 0.02), p(0.0, 0.05)];
    // 59 degrees is unreachable for a sliver quad, and a tiny budget makes
    // sure the call cannot brute-force its way there.
    let quality = Quality {
        min_angle_degrees: 59.0,
        max_steiner_points: 8,
    };
    let (_, outcome) = triangulate_refined(&pts, &[], quality).expect("triangulable");
    assert!(
        !outcome.achieved(),
        "an unreachable bound must not report Achieved"
    );
    match outcome {
        RefineOutcome::Capped {
            worst_angle_degrees,
            ..
        } => assert!(worst_angle_degrees < 59.0),
        RefineOutcome::Achieved { .. } => panic!("unreachable bound reported as achieved"),
        other => panic!("unexpected refine outcome: {other:?}"),
    }
}

/// Refinement is monotone: it never returns a worse mesh than it received.
///
/// Regression test for a real defect found while building this crate.
/// Inserting the circumcentre of an input-pinned sliver produced an even
/// thinner triangle beside it -- the worst angle went 0.51 -> 0.22 degrees.
/// A "refine" that can degrade its input forces every caller to compare
/// before and after, which defeats the point of the call.
#[test]
fn refinement_never_makes_the_mesh_worse() {
    let fixtures: [&[Point2]; 3] = [
        &[
            p(0.0, 0.0),
            p(10.0, 0.0),
            p(10.0, 3.0),
            p(0.0, 3.0),
            p(9.7, 1.4),
            p(9.9, 1.4),
            p(9.9, 1.6),
            p(9.7, 1.6),
        ],
        &[p(0.0, 0.0), p(10.0, 0.0), p(10.0, 0.02), p(0.0, 0.05)],
        &[
            p(0.0, 0.0),
            p(4.0, 0.0),
            p(4.0, 4.0),
            p(0.0, 4.0),
            p(2.0, 2.05),
        ],
    ];
    let quality = Quality {
        min_angle_degrees: 25.0,
        max_steiner_points: 64,
    };
    for (i, pts) in fixtures.iter().enumerate() {
        let before = triangulate(pts, &[]).expect("triangulable");
        let rough = worst_angle(&before);
        let (after, _) = triangulate_refined(pts, &[], quality).expect("triangulable");
        let refined = worst_angle(&after);
        assert!(
            refined >= rough - 1e-9,
            "fixture {i}: refinement degraded the mesh, {rough} -> {refined}"
        );
    }
}

#[test]
fn refinement_preserves_constraints() {
    let pts = [
        p(0.0, 0.0),
        p(8.0, 0.0),
        p(8.0, 5.0),
        p(0.0, 5.0),
        p(4.0, 0.1),
    ];
    let c = Constraint::new(0, 2);
    let quality = Quality {
        min_angle_degrees: 22.0,
        max_steiner_points: 256,
    };
    let (tri, _) = triangulate_refined(&pts, &[c], quality).expect("triangulable");
    let present = (0..tri.triangle_count()).any(|t| {
        let v = &tri.triangles()[3 * t..3 * t + 3];
        v.contains(&0) && v.contains(&2)
    });
    assert!(present, "refinement destroyed a constraint edge");
}
