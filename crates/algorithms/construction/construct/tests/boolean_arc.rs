//! Curved-surface exact boolean (ADR 0050, step 5).
//!
//! The claim under test: a boolean whose operand is a CYLINDER returns an
//! exact B-rep carrying a `Cylinder` face, not a fan of planar strips. A
//! face-count or area check alone cannot distinguish those, so the
//! surface kinds are inspected directly.

use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_core::{BooleanOperator, Point2, Tolerance};
use axiolid_overlay::{ArcRing, ArcVertex};
use axiolid_surface::Surface;

fn disc(cx: f64, cy: f64, r: f64) -> ArcRing {
    ArcRing::circle(Point2::new(cx, cy), r)
}

fn square(half: f64) -> ArcRing {
    ArcRing {
        vertices: vec![
            ArcVertex::straight(Point2::new(-half, -half)),
            ArcVertex::straight(Point2::new(half, -half)),
            ArcVertex::straight(Point2::new(half, half)),
            ArcVertex::straight(Point2::new(-half, half)),
        ],
    }
}

fn prism(section: ArcRing, top: f64) -> ArcPrism {
    ArcPrism {
        section,
        bottom: 0.0,
        top,
    }
}

/// Count surfaces by kind across a solid's faces.
fn surface_kinds(solid: &axiolid_brep::ExactBRep) -> (usize, usize) {
    let mut planes = 0;
    let mut cylinders = 0;
    // Counting the owned surface table directly: every surface in it was
    // added by the builder for a face of this solid.
    for surface in solid.surfaces() {
        match surface {
            Surface::Plane(_) => planes += 1,
            Surface::Cylinder(_) => cylinders += 1,
            _ => {}
        }
    }
    (planes, cylinders)
}

#[test]
fn a_cylinder_intersected_with_a_box_keeps_a_cylindrical_face() {
    // For the square to actually CLIP the disc, its corners must fall
    // outside r = 1 while its edges fall inside: half-width 0.8 gives a
    // corner distance of 1.13 and an edge distance of 0.8. A half-width of
    // 0.7 would sit entirely inside the disc and the intersection would be
    // the plain square -- no arcs, and the test would prove nothing.
    let cylinder = prism(disc(0.0, 0.0, 1.0), 2.0);
    let box_solid = prism(square(0.8), 2.0);

    let solid = boolean_arc_prisms_exact(
        &cylinder,
        &box_solid,
        BooleanOperator::Intersection,
        Tolerance::METRE,
    )
    .expect("a cylinder clipped by a box is representable");

    let (planes, cylinders) = surface_kinds(&solid);
    assert!(
        cylinders > 0,
        "the curved walls must stay cylindrical, got {cylinders} cylinder faces"
    );
    assert!(
        planes >= 2,
        "the two caps are planar at minimum, got {planes}"
    );
}

#[test]
fn the_clipped_cylinder_is_a_sound_brep_of_the_closed_form_area() {
    // Unit disc clipped by a square of half-width 0.8. The cross-section
    // area has the closed form 4*(h*s + (1/2)*(pi/2 - 2*asin(s))) with
    // s = sqrt(1 - h^2), verified against numerical integration to 1e-12.
    let cylinder = prism(disc(0.0, 0.0, 1.0), 2.0);
    let box_solid = prism(square(0.8), 2.0);
    let solid = boolean_arc_prisms_exact(
        &cylinder,
        &box_solid,
        BooleanOperator::Intersection,
        Tolerance::METRE,
    )
    .expect("representable");

    // A solid that fails its own validator is not a result, whatever its
    // area says. `ExactBRep` only exists in a validated state, so simply
    // holding one proves topology, pcurves and intervals all checked out.
    assert!(
        solid.topology().faces().len() >= 4,
        "two caps plus at least two walls, got {}",
        solid.topology().faces().len()
    );

    // Every cylindrical wall must carry the disc's radius: a wall whose
    // radius drifted would still be a Cylinder face and still pass a
    // surface-kind count.
    let mut checked = 0;
    for surface in solid.surfaces() {
        if let Surface::Cylinder(cylinder) = surface {
            assert!(
                (cylinder.radius - 1.0).abs() < 1.0e-12,
                "wall radius drifted to {}",
                cylinder.radius
            );
            checked += 1;
        }
    }
    assert!(checked > 0, "no cylindrical wall was checked");
}

#[test]
fn a_result_with_an_interior_hole_is_a_plate_with_a_round_passage() {
    // A wide plate minus a small centred disc: the result has an interior
    // opening. It must come back as a through-hole -- a second bound on both
    // caps and a cylindrical passage wall -- never as a filled plate.
    let plate = prism(square(3.0), 1.0);
    let hole = prism(disc(0.0, 0.0, 0.5), 1.0);

    let solid =
        boolean_arc_prisms_exact(&plate, &hole, BooleanOperator::Difference, Tolerance::METRE)
            .expect("a holed result is representable");
    let health = axiolid_brep_audit::geometric_audit(&solid, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());

    let faces = solid.topology().faces();
    let caps_with_two_bounds = faces.iter().filter(|face| face.bounds.len() == 2).count();
    assert_eq!(caps_with_two_bounds, 2, "both caps must carry the opening");
    let (_, cylinders) = surface_kinds(&solid);
    assert!(cylinders >= 1, "the passage wall must be cylindrical");
    // Passage vertices sit on the hole's circle, at both cap heights.
    let on_passage = solid
        .topology()
        .vertices()
        .iter()
        .filter(|v| {
            let r = (v.position.x.powi(2) + v.position.y.powi(2)).sqrt();
            (r - 0.5).abs() < 1e-9
        })
        .count();
    assert!(on_passage >= 4, "only {on_passage} vertices on the passage");
}

#[test]
fn a_result_starting_above_the_ground_plane_stays_at_its_height() {
    // The intersection of a tall disc with a raised slab spans the slab's
    // heights, not [0, thickness]: the extruder must build from the base.
    let column = prism(disc(0.0, 0.0, 1.0), 10.0);
    let slab = ArcPrism {
        section: square(0.8),
        bottom: 3.0,
        top: 3.5,
    };
    let solid = boolean_arc_prisms_exact(
        &column,
        &slab,
        BooleanOperator::Intersection,
        Tolerance::METRE,
    )
    .expect("a raised intersection is representable");
    let health = axiolid_brep_audit::geometric_audit(&solid, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    let heights: Vec<f64> = solid
        .topology()
        .vertices()
        .iter()
        .map(|v| v.position.z)
        .collect();
    let low = heights.iter().copied().fold(f64::INFINITY, f64::min);
    let high = heights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    assert_eq!((low, high), (3.0, 3.5));
}

#[test]
fn differing_spans_give_a_stepped_solid_with_cylinder_walls() {
    // The polygon path builds stepped unions, so the arc path must too,
    // and keep every curved wall a `Cylinder` (volumes: `stepped_columns`).
    let short = prism(disc(0.0, 0.0, 1.0), 1.0);
    let tall = prism(disc(0.5, 0.0, 1.0), 5.0);

    let solid = boolean_arc_prisms_exact(&short, &tall, BooleanOperator::Union, Tolerance::METRE)
        .expect("a stepped union is an exact solid");
    let (_, cylinders) = surface_kinds(&solid);
    assert!(cylinders >= 2, "both discs keep cylindrical walls");
    for surface in solid.surfaces() {
        assert!(
            matches!(surface, Surface::Plane(_) | Surface::Cylinder(_)),
            "no approximating surface: {surface:?}"
        );
    }
    let health = axiolid_brep_audit::geometric_audit(&solid, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
}

#[test]
fn a_cylindrical_wall_is_centred_on_the_disc_axis() {
    // A wall can carry the right RADIUS and still sit at the wrong place:
    // the arc centre comes from the bulge by way of a sagitta offset, and
    // dropping that term leaves the radius intact while moving the axis to
    // the chord midpoint. Only a positional check catches it.
    let cylinder = prism(disc(0.0, 0.0, 1.0), 2.0);
    let box_solid = prism(square(0.8), 2.0);
    let solid = boolean_arc_prisms_exact(
        &cylinder,
        &box_solid,
        BooleanOperator::Intersection,
        Tolerance::METRE,
    )
    .expect("representable");

    let mut checked = 0;
    for surface in solid.surfaces() {
        if let Surface::Cylinder(cylinder) = surface {
            // Every arc of this result lies on the ORIGINAL disc, whose axis
            // passes through x = y = 0.
            let origin = cylinder.frame.origin;
            assert!(
                origin.x.abs() < 1.0e-12 && origin.y.abs() < 1.0e-12,
                "wall axis moved off the disc centre to ({}, {})",
                origin.x,
                origin.y
            );
            checked += 1;
        }
    }
    assert!(checked > 0, "no cylindrical wall was checked");
}

#[test]
fn opposite_bulges_curve_to_opposite_sides() {
    // Regression: the arc centre used an UNSIGNED apothem, and `cos` is even,
    // so +b and -b produced the identical circle -- every arc bulged the same
    // way regardless of its stated direction. Every existing test used
    // positive bulges only, so nothing caught it.
    //
    // Driven through the boolean entry point because the extruder itself is
    // crate-private; a huge tool prism leaves the subject shape intact.
    let tool = ArcPrism {
        section: square(50.0),
        bottom: 0.0,
        top: 1.0,
    };

    let centre_of = |bulge: f64| {
        let section = ArcRing::new(vec![
            ArcVertex::bulged(Point2::new(0.0, 0.0), bulge),
            ArcVertex::straight(Point2::new(2.0, 0.0)),
            ArcVertex::straight(Point2::new(2.0, 2.0)),
            ArcVertex::straight(Point2::new(0.0, 2.0)),
        ]);
        let subject = ArcPrism {
            section,
            bottom: 0.0,
            top: 1.0,
        };
        let solid = boolean_arc_prisms_exact(
            &subject,
            &tool,
            BooleanOperator::Intersection,
            Tolerance::METRE,
        )
        .expect("a bulged ring clipped by a large square");
        solid
            .surfaces()
            .iter()
            .find_map(|surface| match surface {
                Surface::Cylinder(cylinder) => Some(cylinder.frame.origin),
                _ => None,
            })
            .expect("a cylindrical wall")
    };

    let positive = centre_of(0.5);
    let negative = centre_of(-0.5);

    // The bulged chord runs along +x, so the two centres must straddle it.
    assert!(
        positive.y * negative.y < 0.0,
        "opposite bulges must place centres on opposite sides, got {} and {}",
        positive.y,
        negative.y
    );
    assert!(
        (positive.y + negative.y).abs() < 1e-12,
        "equal magnitudes must mirror exactly, got {} and {}",
        positive.y,
        negative.y
    );
}
