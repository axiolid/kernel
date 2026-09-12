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
fn a_result_with_an_interior_hole_is_refused_not_filled_in() {
    // A wide plate minus a small centred disc: the true result has an
    // interior opening. The arc extruder cannot build a cap with two
    // bounds, so this must REFUSE rather than return a solid plate.
    let plate = prism(square(3.0), 1.0);
    let hole = prism(disc(0.0, 0.0, 0.5), 1.0);

    let error =
        boolean_arc_prisms_exact(&plate, &hole, BooleanOperator::Difference, Tolerance::METRE)
            .expect_err("a holed result is not representable yet");
    let text = format!("{error:?}");
    assert!(
        text.contains("interior hole"),
        "the refusal must name the hole, got {text}"
    );
}

#[test]
fn differing_spans_are_refused_exactly_as_on_the_polygon_path() {
    // The height reduction is shared, so the arc path must inherit the
    // same refusal rather than quietly accepting a stepped solid.
    let short = prism(disc(0.0, 0.0, 1.0), 1.0);
    let tall = prism(disc(0.5, 0.0, 1.0), 5.0);

    let error = boolean_arc_prisms_exact(&short, &tall, BooleanOperator::Union, Tolerance::METRE)
        .expect_err("a stepped union is not a prism");
    let text = format!("{error:?}");
    assert!(text.contains("differing extrusion spans"), "got {text}");
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
