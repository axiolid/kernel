use std::f64::consts::TAU;

use axiolid_construct::contour_lower::contour_to_arc_ring;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_contracts::{GeomError, Operation};
use axiolid_core::{Frame2, Interval, Point2, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2, Polyline2};
use axiolid_evaluate::{evaluate2, evaluate3};
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile,
};
use axiolid_surface::Surface;
use axiolid_topology::{audit_brep, Orientation};

fn rectangle(thickness: Option<f64>) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x: 4.0,
        y: 2.0,
        thickness,
        outer_radius: None,
        inner_radius: None,
    })
}

fn assert_complete_and_closed(value: &axiolid_construct::ExactBRep) {
    let topology = value.topology();
    assert!(audit_brep(topology).is_closed_manifold());
    assert_eq!(topology.shells().len(), 1);
    assert!(topology.shells()[0].closed);
    assert_eq!(topology.solids().len(), 1);

    for index in 0..topology.edges().len() {
        let edge = topology.edge_id_at(index).expect("dense edge");
        assert!(value.edge_interval(edge).is_some(), "edge {index} interval");
    }
    for (loop_index, wire) in topology.loops().iter().enumerate() {
        let loop_id = topology.loop_id_at(loop_index).expect("dense loop");
        for use_index in 0..wire.edges.len() {
            assert!(
                value.pcurve_interval(loop_id, use_index).is_some(),
                "loop {loop_index} use {use_index} pcurve interval"
            );
        }
    }
}

fn assert_pcurve_surface_agreement(brep: &axiolid_construct::ExactBRep) {
    let topology = brep.topology();
    for face in topology.faces() {
        let surface_id = face.surface.expect("exact face support");
        let surface = &brep.surfaces()[surface_id.index()];
        for bound in &face.bounds {
            let wire = &topology.loops()[bound.loop_id.index()];
            for (use_index, edge_use) in wire.edges.iter().enumerate() {
                let pcurve_id = edge_use.pcurve.expect("exact pcurve");
                let pcurve = &brep.curves2()[pcurve_id.index()];
                let pspan = brep
                    .pcurve_interval(bound.loop_id, use_index)
                    .expect("pcurve span");
                let edge = &topology.edges()[edge_use.edge.index()];
                let curve_id = edge.curve.expect("exact edge support");
                let curve = &brep.curves3()[curve_id.index()];
                let espan = brep.edge_interval(edge_use.edge).expect("edge span");

                for alpha in [0.0, 0.5, 1.0] {
                    let p_t = pspan.start + (pspan.end - pspan.start) * alpha;
                    let e_t = match edge_use.orientation {
                        Orientation::Forward => espan.start + (espan.end - espan.start) * alpha,
                        Orientation::Reversed => espan.end + (espan.start - espan.end) * alpha,
                    };
                    let uv = evaluate2(pcurve, p_t).expect("pcurve evaluation");
                    let from_surface = axiolid_evaluate::surface::evaluate(surface, uv.x, uv.y)
                        .expect("surface evaluation");
                    let from_curve = evaluate3(curve, e_t).expect("edge evaluation");
                    assert!(
                        from_surface.distance(from_curve) <= 1.0e-10,
                        "pcurve/surface mismatch: {from_surface:?} vs {from_curve:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn filled_rectangle_extrusion_is_an_exact_planar_solid() {
    let value = extrude_profile_exact(
        &rectangle(None),
        Vec3::new(0.25, 0.5, 1.0),
        3.0,
        Tolerance::METRE,
    )
    .expect("exact rectangle extrusion");

    assert_complete_and_closed(&value);
    assert_pcurve_surface_agreement(&value);
    assert_eq!(value.topology().vertices().len(), 8);
    assert_eq!(value.topology().edges().len(), 12);
    assert_eq!(value.topology().faces().len(), 6);
    assert_eq!(value.surfaces().len(), 6);
    assert!(value
        .surfaces()
        .iter()
        .all(|surface| matches!(surface, Surface::Plane(_))));
}

#[test]
fn hollow_rectangle_uses_cap_inner_bounds_and_through_passage_faces() {
    let value = extrude_profile_exact(&rectangle(Some(0.25)), Vec3::Z, 2.0, Tolerance::METRE)
        .expect("exact hollow rectangle extrusion");

    assert_complete_and_closed(&value);
    assert_pcurve_surface_agreement(&value);
    assert_eq!(value.topology().vertices().len(), 16);
    assert_eq!(value.topology().edges().len(), 24);
    assert_eq!(value.topology().faces().len(), 10);
    assert_eq!(value.topology().faces()[0].bounds.len(), 2);
    assert_eq!(value.topology().faces()[1].bounds.len(), 2);
    assert_eq!(value.topology().solids()[0].voids.len(), 0);
}

#[test]
fn axial_circle_extrusion_preserves_cylinder_and_distinct_seam_charts() {
    let value = extrude_profile_exact(
        &Profile::Circle(CircleProfile {
            radius: 2.0,
            thickness: None,
        }),
        Vec3::Z,
        5.0,
        Tolerance::METRE,
    )
    .expect("exact cylinder");

    assert_complete_and_closed(&value);
    assert_pcurve_surface_agreement(&value);
    assert_eq!(value.topology().vertices().len(), 2);
    assert_eq!(value.topology().edges().len(), 3);
    assert_eq!(value.topology().faces().len(), 3);
    assert_eq!(
        value
            .surfaces()
            .iter()
            .filter(|surface| matches!(surface, Surface::Cylinder(_)))
            .count(),
        1
    );
    assert_eq!(
        value
            .surfaces()
            .iter()
            .filter(|surface| matches!(surface, Surface::Plane(_)))
            .count(),
        2
    );

    let side = &value.topology().loops()[2];
    let seam_u = side
        .edges
        .iter()
        .filter_map(|edge_use| edge_use.pcurve)
        .filter_map(|id| match &value.curves2()[id.index()] {
            Curve2::Line(line) if line.direction.x == 0.0 => Some(line.origin.x),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(seam_u.len(), 2);
    assert!(seam_u.contains(&0.0));
    assert!(seam_u.contains(&TAU));
}

#[test]
fn exact_extrusion_refuses_families_whose_supports_are_not_populated() {
    // A rounded rectangle used to sit in this list. It now extrudes exactly
    // through the contour path (kernel#111); see tests/rounded_rectangle.rs.
    let annulus = Profile::Circle(CircleProfile {
        radius: 1.0,
        thickness: Some(0.2),
    });
    // A list rather than one case: every family that still refuses belongs
    // here, and a list keeps adding the next one a one-line change.
    let cases = vec![(annulus, "annular circle extrusion")];
    for (profile, family) in cases {
        let error = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
            .expect_err("unsupported exact family must refuse");
        assert!(matches!(
            error,
            GeomError::UnsupportedInput {
                operation: Operation::Sweep,
                input,
                ..
            } if input == family
        ));
    }
}

/// An oblique circle sweeps an oblique circular cylinder (#280), never a
/// right cylinder mislabelled as its support: the wall is the elliptical
/// cylinder with its axis along the direction.
#[test]
fn oblique_circle_extrusion_builds_an_oblique_cylinder_not_a_right_one() {
    let profile = Profile::Circle(CircleProfile {
        radius: 1.0,
        thickness: None,
    });
    let direction = Vec3::new(1.0, 0.0, 1.0);
    let solid = extrude_profile_exact(&profile, direction, 1.0, Tolerance::METRE)
        .expect("an oblique circle extrudes exactly");
    assert_complete_and_closed(&solid);
    assert_pcurve_surface_agreement(&solid);
    let walls: Vec<_> = solid
        .surfaces()
        .iter()
        .filter(|surface| !matches!(surface, Surface::Plane(_)))
        .collect();
    assert_eq!(walls.len(), 1);
    let Surface::EllipticalCylinder(wall) = walls[0] else {
        panic!(
            "an oblique circle's wall is not a right cylinder: {:?}",
            walls[0]
        );
    };
    let unit = direction.normalize();
    assert!((wall.frame.z - unit).length() < 1e-15);
    assert!((wall.semi_axis_x - 1.0).abs() < 1e-15);
    assert!((wall.semi_axis_y - unit.z).abs() < 1e-15);
}

/// A lean below tolerance is still a lean: the exact constructor builds the
/// oblique cylinder it denotes rather than snapping it to the normal.
#[test]
fn sub_tolerance_oblique_circle_is_not_snapped_to_the_normal() {
    let profile = Profile::Circle(CircleProfile {
        radius: 1.0,
        thickness: None,
    });
    let direction = Vec3::new(Tolerance::METRE.linear() * 0.5, 0.0, 1.0);
    let solid = extrude_profile_exact(&profile, direction, 1.0, Tolerance::METRE)
        .expect("a slightly leaning circle extrudes exactly");
    assert_complete_and_closed(&solid);
    assert_pcurve_surface_agreement(&solid);
    assert!(solid
        .surfaces()
        .iter()
        .any(|surface| matches!(surface, Surface::EllipticalCylinder(_))));
}

#[test]
fn nonzero_direction_magnitude_is_not_a_model_length() {
    let direction = Vec3::Z * (Tolerance::METRE.linear() * 0.5);
    let value = extrude_profile_exact(
        &Profile::Circle(CircleProfile {
            radius: 1.0,
            thickness: None,
        }),
        direction,
        1.0,
        Tolerance::METRE,
    )
    .expect("direction magnitude is unitless and normalization preserves the axis");

    assert_complete_and_closed(&value);
    assert_pcurve_surface_agreement(&value);
}

#[test]
fn finite_nonzero_extreme_direction_magnitudes_normalize_safely() {
    let profile = Profile::Circle(CircleProfile {
        radius: 1.0,
        thickness: None,
    });

    for magnitude in [f64::MIN_POSITIVE, f64::MAX] {
        let value = extrude_profile_exact(
            &profile,
            Vec3::new(0.0, 0.0, magnitude),
            1.0,
            Tolerance::METRE,
        )
        .expect("finite nonzero unitless directions must normalize without range loss");

        assert_complete_and_closed(&value);
        assert_pcurve_surface_agreement(&value);
    }
}

#[test]
fn mixed_extreme_direction_refuses_nonzero_component_loss() {
    let error = extrude_profile_exact(
        &rectangle(None),
        Vec3::new(f64::MIN_POSITIVE, 0.0, f64::MAX),
        1.0,
        Tolerance::METRE,
    )
    .expect_err("normalization must not silently erase a nonzero direction component");

    assert!(matches!(error, GeomError::Degenerate(_)));
}

#[test]
fn invalid_exact_extrusion_inputs_are_not_reclassified_as_unsupported() {
    let error = extrude_profile_exact(
        &rectangle(None),
        Vec2::ZERO.extend(0.0),
        1.0,
        Tolerance::METRE,
    )
    .expect_err("zero direction");
    assert!(matches!(error, GeomError::InvalidInput(_)));
}

// ----------------------------------------------------------- polylines (#389)

fn pt(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

/// A rectangle with coordinates that are not dyadic, so a ring vertex that
/// was recomputed rather than copied would differ in the last bit.
fn rect_points() -> Vec<Point2> {
    vec![pt(0.1, 0.2), pt(2.3, 0.2), pt(2.3, 1.7), pt(0.1, 1.7)]
}

fn polyline_segment(
    points: &[Point2],
    closed: bool,
    domain: (f64, f64),
    same_sense: bool,
) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Polyline(Polyline2 {
            points: points.to_vec(),
            closed,
        }),
        domain: Interval::new(domain.0, domain.1),
        same_sense,
    }
}

fn line_segment(from: Point2, to: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}

fn contour_profile(segments: Vec<ProfileSegment>) -> Profile {
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

fn shoelace(points: &[Point2]) -> f64 {
    let n = points.len();
    (0..n)
        .map(|i| points[i].perp_dot(points[(i + 1) % n]))
        .sum::<f64>()
        / 2.0
}

/// Extrude exactly, require a closed B-rep, and return its volume.
fn extruded_volume(profile: &Profile, depth: f64) -> f64 {
    let value = extrude_profile_exact(profile, Vec3::Z, depth, Tolerance::METRE)
        .expect("the polyline contour extrudes exactly");
    assert_complete_and_closed(&value);
    axiolid_measure::exact_properties(&value, Tolerance::METRE)
        .expect("measurable")
        .signed_volume
}

fn ring_points(segments: Vec<ProfileSegment>) -> Vec<(Point2, f64)> {
    contour_to_arc_ring(&Contour::new(segments), Tolerance::METRE)
        .expect("lowers")
        .vertices
        .iter()
        .map(|v| (v.point, v.bulge))
        .collect()
}

fn assert_bit_equal(found: Point2, expected: Point2) {
    assert_eq!(found.x.to_bits(), expected.x.to_bits(), "{found:?}");
    assert_eq!(found.y.to_bits(), expected.y.to_bits(), "{found:?}");
}

#[test]
fn an_open_polyline_and_a_closing_line_extrude_to_the_shoelace_volume() {
    let points = rect_points();
    let depth = 1.5;
    let expected = shoelace(&points) * depth;
    let volume = extruded_volume(
        &contour_profile(vec![
            polyline_segment(&points, false, (0.0, 3.0), true),
            line_segment(points[3], points[0]),
        ]),
        depth,
    );
    assert!((volume - expected).abs() < 1e-12, "{volume} vs {expected}");
}

#[test]
fn one_closed_polyline_is_a_whole_contour_of_straight_edges() {
    let points = rect_points();
    let segment = polyline_segment(&points, true, (0.0, 4.0), true);
    let ring = ring_points(vec![segment.clone()]);
    assert_eq!(ring.len(), 4);
    for ((point, bulge), expected) in ring.iter().zip(&points) {
        assert_bit_equal(*point, *expected);
        assert_eq!(*bulge, 0.0);
    }
    let depth = 2.0;
    let volume = extruded_volume(&contour_profile(vec![segment]), depth);
    assert!(
        (volume - shoelace(&points) * depth).abs() < 1e-12,
        "{volume}"
    );
}

#[test]
fn a_reversed_polyline_is_the_same_solid() {
    let points = rect_points();
    let forward = extruded_volume(
        &contour_profile(vec![polyline_segment(&points, true, (0.0, 4.0), true)]),
        1.0,
    );
    let backward_segment = polyline_segment(&points, true, (0.0, 4.0), false);
    let backward = extruded_volume(&contour_profile(vec![backward_segment.clone()]), 1.0);
    assert!(
        (forward - backward).abs() < 1e-12,
        "{forward} vs {backward}"
    );

    // Traversed from the end: the points in reverse order, closing on the
    // first.
    let ring = ring_points(vec![backward_segment]);
    let order = [0, 3, 2, 1];
    assert_eq!(ring.len(), 4);
    for ((point, _), index) in ring.iter().zip(order) {
        assert_bit_equal(*point, points[index]);
    }
}

#[test]
fn a_partial_domain_starts_and_ends_between_vertices() {
    let points = vec![pt(0.0, 0.0), pt(2.0, 0.0), pt(2.0, 2.0), pt(0.0, 2.0)];
    let ring = ring_points(vec![
        polyline_segment(&points, true, (0.5, 2.5), true),
        line_segment(pt(1.0, 2.0), pt(1.0, 0.0)),
    ]);
    let expected = [pt(1.0, 0.0), pt(2.0, 0.0), pt(2.0, 2.0), pt(1.0, 2.0)];
    assert_eq!(ring.len(), expected.len());
    for ((point, bulge), expected) in ring.iter().zip(expected) {
        assert_bit_equal(*point, expected);
        assert_eq!(*bulge, 0.0);
    }
}

#[test]
fn a_polyline_and_an_arc_make_a_mixed_contour() {
    // The upper half of the unit circle, then the other three sides of a
    // rectangle below it.
    let arch = ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: pt(0.0, 0.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 1.0,
        }),
        domain: Interval::new(0.0, core::f64::consts::PI),
        same_sense: true,
    };
    let base = [pt(-1.0, 0.0), pt(-1.0, -1.0), pt(1.0, -1.0), pt(1.0, 0.0)];
    let depth = 1.0;
    let volume = extruded_volume(
        &contour_profile(vec![arch, polyline_segment(&base, false, (0.0, 3.0), true)]),
        depth,
    );
    let expected = (2.0 + core::f64::consts::FRAC_PI_2) * depth;
    assert!((volume - expected).abs() < 1e-9, "{volume} vs {expected}");
}

#[test]
fn a_polyline_domain_ends_at_its_last_point_and_no_further() {
    let points = rect_points();
    let too_far = contour_profile(vec![polyline_segment(&points, true, (0.0, 5.0), true)]);
    let error = extrude_profile_exact(&too_far, Vec3::Z, 1.0, Tolerance::METRE)
        .expect_err("past the closing point");
    assert!(matches!(error, GeomError::Degenerate(_)), "{error:?}");
}

fn square_points() -> Vec<Point2> {
    vec![pt(0.0, 0.0), pt(2.0, 0.0), pt(2.0, 2.0), pt(0.0, 2.0)]
}

#[test]
fn a_closed_polyline_repeating_its_first_point_is_the_same_square_in_both_senses() {
    let mut points = square_points();
    points.push(points[0]);
    for same_sense in [true, false] {
        let segment = polyline_segment(&points, true, (0.0, 5.0), same_sense);
        assert_eq!(ring_points(vec![segment.clone()]).len(), 4, "{same_sense}");
        let volume = extruded_volume(&contour_profile(vec![segment]), 1.5);
        assert!((volume - 6.0).abs() < 1e-12, "{same_sense}: {volume}");
    }
}

#[test]
fn a_backward_partial_domain_starts_below_its_first_vertex() {
    // Opposite sense on 0.5..3.5, so from 3.5 down to 0.5: the vertices 3, 2 and 1 lie between.
    let points = square_points();
    let segment = polyline_segment(&points, true, (0.5, 3.5), false);
    let closing = line_segment(pt(1.0, 0.0), pt(0.0, 1.0));
    let ring = ring_points(vec![segment.clone(), closing.clone()]);
    let expected = [
        pt(0.0, 1.0),
        pt(0.0, 2.0),
        pt(2.0, 2.0),
        pt(2.0, 0.0),
        pt(1.0, 0.0),
    ];
    assert_eq!(ring.len(), expected.len());
    for ((point, _), expected) in ring.iter().zip(expected) {
        assert_bit_equal(*point, expected);
    }
    let volume = extruded_volume(&contour_profile(vec![segment, closing]), 1.0);
    assert!((volume - 3.5).abs() < 1e-12, "{volume}");
}

#[test]
fn a_domain_crossing_the_seam_ends_on_the_first_point() {
    let points = square_points();
    let segment = polyline_segment(&points, true, (2.5, 4.0), true);
    let volume = extruded_volume(
        &contour_profile(vec![segment, line_segment(pt(0.0, 0.0), pt(1.0, 2.0))]),
        1.0,
    );
    assert!((volume - 1.0).abs() < 1e-12, "{volume}");
}

#[test]
fn a_single_full_circle_contour_extrudes() {
    let circle = ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: pt(0.0, 0.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 1.0,
        }),
        domain: Interval::new(0.0, core::f64::consts::TAU),
        same_sense: true,
    };
    let volume = extruded_volume(&contour_profile(vec![circle]), 2.0);
    assert!(
        (volume - 2.0 * core::f64::consts::PI).abs() < 1e-9,
        "{volume}"
    );
}
