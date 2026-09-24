//! Booleans whose result falls apart into several solids (#120).
//!
//! `boolean_prisms_exact` and `boolean_arc_prisms_exact` return one
//! `ExactBRep`, which is one solid, so they refuse a result with several
//! separate pieces rather than silently drop material. The `_solids`
//! variants return every piece. The oracle is hand-derived per piece:
//! volume for all-planar pieces, extent and wall radius for curved ones.

use axiolid_brep::ExactBRep;
use axiolid_construct::boolean_exact::{
    boolean_arc_prisms_exact, boolean_arc_prisms_exact_solids, boolean_prisms_exact,
    boolean_prisms_exact_solids, ArcPrism, Prism,
};
use axiolid_contracts::GeomError;
use axiolid_core::{BooleanOperator, Point2, Tolerance};
use axiolid_overlay::{ArcRing, ArcVertex};
use axiolid_surface::Surface;

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Point2> {
    vec![
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ]
}

fn prism(ring: Vec<Point2>, bottom: f64, top: f64) -> Prism {
    Prism {
        rings: vec![ring],
        bottom,
        top,
    }
}

fn arc_rect(x0: f64, y0: f64, x1: f64, y1: f64) -> ArcRing {
    ArcRing {
        vertices: rect(x0, y0, x1, y1)
            .into_iter()
            .map(ArcVertex::straight)
            .collect(),
    }
}

fn volume(solid: &ExactBRep) -> f64 {
    axiolid_measure::exact_properties(solid, Tolerance::METRE)
        .expect("all-planar piece is measurable")
        .signed_volume
}

/// Axis-aligned bounds of a solid's vertices: (min x, max x, min z, max z).
fn extent(solid: &ExactBRep) -> (f64, f64, f64, f64) {
    let mut e = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for v in solid.topology().vertices() {
        e.0 = e.0.min(v.position.x);
        e.1 = e.1.max(v.position.x);
        e.2 = e.2.min(v.position.z);
        e.3 = e.3.max(v.position.z);
    }
    e
}

fn assert_sound(solid: &ExactBRep) {
    let health = axiolid_brep_audit::geometric_audit(solid, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
}

fn is_disconnected_refusal(error: &GeomError) -> bool {
    matches!(error, GeomError::UnsupportedInput { input, .. } if input.contains("disconnected"))
}

#[test]
fn a_wall_cut_through_by_a_full_height_slot_is_two_walls() {
    // Wall 10 x 1 x 3, slot 2 wide through its whole depth and height at
    // x in [4, 6]: two walls of 4 x 1 x 3 = 12 each, left one first.
    let wall = prism(rect(0.0, 0.0, 10.0, 1.0), 0.0, 3.0);
    let slot = prism(rect(4.0, -1.0, 6.0, 2.0), 0.0, 3.0);

    let error = boolean_prisms_exact(&wall, &slot, BooleanOperator::Difference, Tolerance::METRE)
        .expect_err("one ExactBRep cannot hold two solids");
    assert!(is_disconnected_refusal(&error), "{error:?}");

    let pieces =
        boolean_prisms_exact_solids(&wall, &slot, BooleanOperator::Difference, Tolerance::METRE)
            .expect("two pieces are representable");
    assert_eq!(pieces.len(), 2);
    for piece in &pieces {
        assert_sound(piece);
        assert!((volume(piece) - 12.0).abs() < 1e-12, "{}", volume(piece));
    }
    assert_eq!(extent(&pieces[0]).0, 0.0, "left piece first");
    assert_eq!(extent(&pieces[1]).1, 10.0, "right piece second");
}

#[test]
fn three_separate_islands_come_back_in_a_stable_order() {
    // A comb: three teeth joined by a spine, minus the spine, leaves three
    // separate teeth. Ordered by lowest vertex, x first.
    let comb = Prism {
        rings: vec![vec![
            Point2::new(0.0, 0.0),
            Point2::new(5.0, 0.0),
            Point2::new(5.0, 3.0),
            Point2::new(4.0, 3.0),
            Point2::new(4.0, 1.0),
            Point2::new(3.0, 1.0),
            Point2::new(3.0, 3.0),
            Point2::new(2.0, 3.0),
            Point2::new(2.0, 1.0),
            Point2::new(1.0, 1.0),
            Point2::new(1.0, 3.0),
            Point2::new(0.0, 3.0),
        ]],
        bottom: 0.0,
        top: 2.0,
    };
    let spine = prism(rect(-1.0, -1.0, 6.0, 1.0), 0.0, 2.0);
    let pieces =
        boolean_prisms_exact_solids(&comb, &spine, BooleanOperator::Difference, Tolerance::METRE)
            .expect("three teeth");
    assert_eq!(pieces.len(), 3);
    let starts: Vec<f64> = pieces.iter().map(|p| extent(p).0).collect();
    assert_eq!(starts, vec![0.0, 2.0, 4.0]);
    for piece in &pieces {
        assert_sound(piece);
        assert!((volume(piece) - 4.0).abs() < 1e-12, "1 x 2 x 2 tooth");
    }
}

#[test]
fn a_single_piece_matches_the_single_solid_entry_point() {
    let wall = prism(rect(0.0, 0.0, 10.0, 1.0), 0.0, 3.0);
    let notch = prism(rect(4.0, 0.5, 6.0, 2.0), 0.0, 3.0);
    let one = boolean_prisms_exact(&wall, &notch, BooleanOperator::Difference, Tolerance::METRE)
        .expect("a notch leaves one solid");
    let all =
        boolean_prisms_exact_solids(&wall, &notch, BooleanOperator::Difference, Tolerance::METRE)
            .expect("same");
    assert_eq!(all, vec![one]);
}

#[test]
fn an_empty_result_is_an_empty_list_not_an_error() {
    let a = prism(rect(0.0, 0.0, 1.0, 1.0), 0.0, 1.0);
    let far = prism(rect(5.0, 5.0, 6.0, 6.0), 0.0, 1.0);
    let above = prism(rect(0.0, 0.0, 1.0, 1.0), 2.0, 3.0);
    for (tool, what) in [(&far, "disjoint sections"), (&above, "disjoint heights")] {
        let pieces =
            boolean_prisms_exact_solids(&a, tool, BooleanOperator::Intersection, Tolerance::METRE)
                .unwrap_or_else(|e| panic!("{what}: {e:?}"));
        assert!(pieces.is_empty(), "{what}");
    }
    // The single-solid entry point keeps refusing, as before.
    assert!(matches!(
        boolean_prisms_exact(&a, &far, BooleanOperator::Intersection, Tolerance::METRE),
        Err(GeomError::Degenerate(_))
    ));
}

#[test]
fn refusals_other_than_emptiness_are_kept() {
    // An enclosed cavity is still refused, even by the multi-solid path.
    let a = prism(rect(0.0, 0.0, 1.0, 1.0), 0.0, 1.0);
    let buried = prism(rect(0.25, 0.25, 0.75, 0.75), 0.25, 0.75);
    let error =
        boolean_prisms_exact_solids(&a, &buried, BooleanOperator::Difference, Tolerance::METRE)
            .expect_err("cavity");
    assert!(
        matches!(error, GeomError::UnsupportedInput { input, .. } if input.contains("cavity")),
        "{error:?}"
    );
    // Malformed input is still malformed.
    let bad = prism(rect(0.0, 0.0, 1.0, 1.0), 1.0, 0.0);
    assert!(matches!(
        boolean_prisms_exact_solids(&a, &bad, BooleanOperator::Union, Tolerance::METRE),
        Err(GeomError::InvalidInput(_))
    ));
}

#[test]
fn a_disc_cut_by_a_strip_is_two_curved_pieces() {
    // Unit-radius column minus a strip |x| <= 0.25: two D-shaped pieces,
    // each with a cylindrical wall of the column's radius.
    let column = ArcPrism {
        section: ArcRing::circle(Point2::new(0.0, 0.0), 1.0),
        bottom: 0.0,
        top: 2.0,
    };
    let strip = ArcPrism {
        section: arc_rect(-0.25, -2.0, 0.25, 2.0),
        bottom: 0.0,
        top: 2.0,
    };
    let error = boolean_arc_prisms_exact(
        &column,
        &strip,
        BooleanOperator::Difference,
        Tolerance::METRE,
    )
    .expect_err("one ExactBRep cannot hold two solids");
    assert!(is_disconnected_refusal(&error), "{error:?}");

    let pieces = boolean_arc_prisms_exact_solids(
        &column,
        &strip,
        BooleanOperator::Difference,
        Tolerance::METRE,
    )
    .expect("two pieces");
    assert_eq!(pieces.len(), 2);
    let spans: Vec<(f64, f64)> = pieces
        .iter()
        .map(|p| {
            let e = extent(p);
            (e.0, e.1)
        })
        .collect();
    assert!((spans[0].0 + 1.0).abs() < 1e-9 && (spans[0].1 + 0.25).abs() < 1e-9);
    assert!((spans[1].0 - 0.25).abs() < 1e-9 && (spans[1].1 - 1.0).abs() < 1e-9);
    for piece in &pieces {
        assert_sound(piece);
        let e = extent(piece);
        assert_eq!((e.2, e.3), (0.0, 2.0));
        let radii: Vec<f64> = piece
            .surfaces()
            .iter()
            .filter_map(|s| match s {
                Surface::Cylinder(c) => Some(c.radius),
                _ => None,
            })
            .collect();
        assert!(!radii.is_empty(), "the curved wall stays a cylinder");
        assert!(radii.iter().all(|r| (r - 1.0).abs() < 1e-12), "{radii:?}");
    }
}

#[test]
fn raised_curved_pieces_keep_their_height() {
    // A column standing on [1, 2] minus a full-height strip: two halves,
    // both still spanning [1, 2] rather than being rebuilt from z = 0.
    let raised = ArcPrism {
        section: ArcRing::circle(Point2::new(0.0, 0.0), 1.0),
        bottom: 1.0,
        top: 2.0,
    };
    let strip = ArcPrism {
        section: arc_rect(-0.25, -2.0, 0.25, 2.0),
        bottom: 0.0,
        top: 5.0,
    };
    let pieces = boolean_arc_prisms_exact_solids(
        &raised,
        &strip,
        BooleanOperator::Difference,
        Tolerance::METRE,
    )
    .expect("the tool spans the subject's height");
    assert_eq!(pieces.len(), 2);
    for piece in &pieces {
        assert_sound(piece);
        let e = extent(piece);
        assert_eq!((e.2, e.3), (1.0, 2.0));
    }
}

#[test]
fn pieces_sharing_their_lowest_x_are_ordered_by_y() {
    // A tall bar cut by two horizontal slots leaves three pieces, all with
    // their lowest vertex at x = 0: only the y tie-break orders them, and
    // it must give bottom, middle, top whatever order the overlay emits.
    let bar = prism(rect(0.0, 0.0, 1.0, 9.0), 0.0, 1.0);
    let slots = Prism {
        rings: vec![rect(-1.0, 2.0, 2.0, 3.0)],
        bottom: 0.0,
        top: 1.0,
    };
    let first =
        boolean_prisms_exact_solids(&bar, &slots, BooleanOperator::Difference, Tolerance::METRE)
            .expect("two pieces");
    assert_eq!(first.len(), 2);
    let slots_two = Prism {
        rings: vec![vec![
            Point2::new(-1.0, 2.0),
            Point2::new(2.0, 2.0),
            Point2::new(2.0, 6.0),
            Point2::new(-1.0, 6.0),
            Point2::new(-1.0, 5.0),
            Point2::new(1.5, 5.0),
            Point2::new(1.5, 3.0),
            Point2::new(-1.0, 3.0),
        ]],
        bottom: 0.0,
        top: 1.0,
    };
    let pieces = boolean_prisms_exact_solids(
        &bar,
        &slots_two,
        BooleanOperator::Difference,
        Tolerance::METRE,
    )
    .expect("three pieces");
    assert_eq!(pieces.len(), 3);
    let lows: Vec<f64> = pieces
        .iter()
        .map(|piece| {
            piece
                .topology()
                .vertices()
                .iter()
                .map(|v| v.position.y)
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    assert_eq!(lows, vec![0.0, 3.0, 6.0]);
    // Bottom 1x2x1, middle 1x2x1 (y in [3, 5]), top 1x3x1 (y in [6, 9]).
    let volumes: Vec<f64> = pieces.iter().map(volume).collect();
    for (got, want) in volumes.iter().zip([2.0, 2.0, 3.0]) {
        assert!((got - want).abs() < 1e-12, "{volumes:?}");
    }
}

#[test]
fn arc_pieces_are_ordered_left_to_right_whatever_the_operand_order() {
    // Two discs from different operands, the subject's on the right. The
    // overlay emits regions in its own linking order; the result must
    // still come back lowest x first, in both argument orders.
    let right = ArcPrism {
        section: ArcRing::circle(Point2::new(5.0, 0.0), 1.0),
        bottom: 0.0,
        top: 1.0,
    };
    let left = ArcPrism {
        section: ArcRing::circle(Point2::new(0.0, 0.0), 1.0),
        bottom: 0.0,
        top: 1.0,
    };
    for (subject, tool) in [(&right, &left), (&left, &right)] {
        let pieces = boolean_arc_prisms_exact_solids(
            subject,
            tool,
            BooleanOperator::Union,
            Tolerance::METRE,
        )
        .expect("two separate cylinders");
        let lows: Vec<f64> = pieces
            .iter()
            .map(|piece| {
                piece
                    .topology()
                    .vertices()
                    .iter()
                    .map(|v| v.position.x)
                    .fold(f64::INFINITY, f64::min)
            })
            .collect();
        assert_eq!(lows.len(), 2);
        assert!(lows[0] < 0.0 && lows[1] > 3.0, "left first: {lows:?}");
    }
}

#[test]
fn a_slot_cut_through_the_middle_heights_leaves_two_slabs() {
    // A 1 x 1 x 3 block minus a wider plate at z 1..2: two separate 1 x 1 x 1
    // slabs, one above the other. Stepped AND disconnected.
    let block = prism(rect(0.0, 0.0, 1.0, 1.0), 0.0, 3.0);
    let plate = prism(rect(-1.0, -1.0, 2.0, 2.0), 1.0, 2.0);
    let solids = boolean_prisms_exact_solids(
        &block,
        &plate,
        BooleanOperator::Difference,
        Tolerance::METRE,
    )
    .expect("two slabs");
    assert_eq!(solids.len(), 2);
    for solid in &solids {
        assert_sound(solid);
        assert!((volume(solid) - 1.0).abs() < 1e-12, "{}", volume(solid));
    }
    let heights: Vec<_> = solids.iter().map(|s| (extent(s).2, extent(s).3)).collect();
    assert!(
        heights.contains(&(0.0, 1.0)) && heights.contains(&(2.0, 3.0)),
        "{heights:?}"
    );
    // The single-solid entry point refuses rather than dropping a slab.
    let error = boolean_prisms_exact(
        &block,
        &plate,
        BooleanOperator::Difference,
        Tolerance::METRE,
    )
    .expect_err("two pieces");
    assert!(is_disconnected_refusal(&error), "{error:?}");
}
