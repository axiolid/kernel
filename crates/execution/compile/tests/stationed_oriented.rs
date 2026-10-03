//! Oriented and tagged stations through the compiler (#246): explicit
//! axis and reference direction against hand-computed frames on a line, an
//! arc, an elevated and a banked curve; spines matched by tag against
//! their index-matched equivalents; refusals by name.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Frame2, Interval, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CantPiece, Circle2, Curve2, Curve3, Elevated3, ElevationLaw,
    Line2, Polyline2,
};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::station::resolve;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::{MeshClosure, MeshCompiler};
use axiolid_model::{
    CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode, GraphError, NodeId,
    OpenProfile, OrientedCurveStation, SectionAtStation, SolidOperation, Station, StationFrame,
    StationOffsets, StationOrientation, StationedSection, SurfaceRelation,
};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile};

const EPS: Scalar = 1e-9;

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

fn compile(graph: &GeometryGraph, root: NodeId) -> GeomResult<TriMesh> {
    compiler().compile_mesh(graph, root, &options())
}

fn close3(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    assert!(
        (actual - expected).abs().max_element() <= eps,
        "{what}: {actual:?} != {expected:?}"
    );
}

fn close(actual: Scalar, expected: Scalar, relative: Scalar, what: &str) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{what}: {actual} != {expected}"
    );
}

fn volume(mesh: &TriMesh) -> Scalar {
    volume_properties(mesh, Tolerance::MILLIMETRE)
        .expect("a stationed spine is closed")
        .signed_volume
}

fn bounds(mesh: &TriMesh) -> (Point3, Point3) {
    mesh.positions.iter().fold(
        (
            Point3::splat(Scalar::INFINITY),
            Point3::splat(-Scalar::INFINITY),
        ),
        |(lo, hi), p| (lo.min(*p), hi.max(*p)),
    )
}

/// A plan line along `+y` from `(2, 1)`: tangent `+y`, lateral `-x`, up
/// `+Z`.
fn plan_line(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve2::Line(Line2 {
        origin: Point2::new(2.0, 1.0),
        direction: Vec2::new(0.0, 2.0),
    }))
    .unwrap()
}

fn orientation(axis: Option<Vec3>, reference: Option<Vec3>) -> StationOrientation {
    StationOrientation::new(axis, reference)
}

/// Resolve an oriented station node on `basis`.
fn resolve_oriented(
    b: GeometryGraphBuilder,
    basis: NodeId,
    station: Station,
    frame: StationFrame,
    turn: StationOrientation,
) -> axiolid_mesh_compile::station::ResolvedStation {
    let mut b = b;
    let id = b
        .push_value(OrientedCurveStation::new(
            CurveStation {
                basis,
                station,
                frame,
            },
            turn,
        ))
        .unwrap();
    let graph = b.finish(vec![id]).unwrap();
    resolve(&graph, id).unwrap()
}

/// The presented frame of tangent `t`, lateral `l` and up `u`: `x = t`,
/// `y = u`, `z = -l`.
fn assert_frame(
    resolved: &axiolid_mesh_compile::station::ResolvedStation,
    [t, l, u]: [Vec3; 3],
    what: &str,
) {
    close3(resolved.frame.x, t, EPS, &format!("{what}: x tangent"));
    close3(resolved.frame.y, u, EPS, &format!("{what}: y up"));
    close3(resolved.frame.z, -l, EPS, &format!("{what}: z right"));
    close3(resolved.frame.origin, resolved.point, 0.0, what);
}

// --- resolved frames ----------------------------------------------------------

#[test]
fn an_axis_rolls_a_line_station_about_its_tangent_and_leaves_its_point() {
    // Base: t = +y, l = -x, u = +Z. Axis (0, -sin, cos) in (t, l, u) is
    // up' = sin * x + cos * z; the default reference keeps t' = t, and
    // l' = up' x t' = (-cos, 0, sin).
    let theta: Scalar = 0.3;
    let (s, c) = theta.sin_cos();
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let resolved = resolve_oriented(
        b,
        basis,
        Station::new(4.0, StationOffsets::new(1.5, 0.25, 0.5)),
        StationFrame::Section,
        orientation(Some(Vec3::new(0.0, -s, c)), None),
    );
    // The offsets stay in the base frame: the point does not move.
    close3(resolved.point, Point3::new(0.5, 5.5, 0.25), EPS, "point");
    assert_frame(
        &resolved,
        [Vec3::Y, Vec3::new(-c, 0.0, s), Vec3::new(s, 0.0, c)],
        "rolled",
    );
    // The base section is reported unturned.
    close3(resolved.section.up, Vec3::Z, EPS, "base up");
}

#[test]
fn a_reference_direction_off_the_axis_is_made_perpendicular_to_it() {
    // Axis (0, 0, 2) is up, not unit; reference (1, 1, 1) loses its up
    // component, so t' = (t + l) / sqrt 2 = (-1, 1, 0) / sqrt 2 and
    // l' = Z x t' = (-1, -1, 0) / sqrt 2.
    let r = core::f64::consts::FRAC_1_SQRT_2;
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let resolved = resolve_oriented(
        b,
        basis,
        Station::at(4.0),
        StationFrame::Section,
        orientation(Some(Vec3::new(0.0, 0.0, 2.0)), Some(Vec3::ONE)),
    );
    close3(resolved.point, Point3::new(2.0, 5.0, 0.0), EPS, "point");
    assert_frame(
        &resolved,
        [Vec3::new(-r, r, 0.0), Vec3::new(-r, -r, 0.0), Vec3::Z],
        "Gram-Schmidt",
    );
}

#[test]
fn an_arc_station_can_stand_its_section_along_the_tangent() {
    // A quarter of the way round a counter-clockwise R = 30 circle: point
    // (0, 30), t = -x, l = -y (towards the centre), u = +Z. Axis (1, 0, 0)
    // and reference (0, 1, 0) make up' = t = -x, t' = l = -y and
    // l' = up' x t' = +Z.
    let mut b = GeometryGraphBuilder::new();
    let radius = 30.0;
    let basis = b
        .push_value(Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius,
        }))
        .unwrap();
    let resolved = resolve_oriented(
        b,
        basis,
        Station::new(
            0.5 * core::f64::consts::PI * radius,
            StationOffsets::new(2.0, 0.0, 0.0),
        ),
        StationFrame::Section,
        orientation(Some(Vec3::X), Some(Vec3::Y)),
    );
    close3(resolved.point, Point3::new(0.0, 28.0, 0.0), 1e-6, "point");
    let tight = |v: Vec3, e: Vec3, what| close3(v, e, 1e-6, what);
    tight(resolved.frame.x, -Vec3::Y, "x");
    tight(resolved.frame.y, -Vec3::X, "y");
    tight(resolved.frame.z, -Vec3::Z, "z");
}

#[test]
fn an_elevated_station_turns_within_its_leaning_section_frame() {
    // 2% grade along +x at plan distance 50: t = (1, 0, 0.02) / k,
    // l = +y, u = (-0.02, 0, 1) / k. Axis (0, 1, 1) is up' = (l + u) / sqrt 2;
    // the default reference is already perpendicular, so t' = t and
    // l' = (l - u) / sqrt 2.
    let k = 1.0004_f64.sqrt();
    let r = core::f64::consts::FRAC_1_SQRT_2;
    let t = Vec3::new(1.0, 0.0, 0.02) / k;
    let l = Vec3::Y;
    let u = Vec3::new(-0.02, 0.0, 1.0) / k;
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push_value(Curve3::Elevated(Elevated3::new(
            Curve2::Line(Line2 {
                origin: Point2::ZERO,
                direction: Vec2::X,
            }),
            ElevationLaw::constant_grade(100.0, 0.02),
        )))
        .unwrap();
    let resolved = resolve_oriented(
        b,
        basis,
        Station::new(50.0, StationOffsets::new(0.0, 1.0, 0.0)),
        StationFrame::Section,
        orientation(Some(Vec3::new(0.0, 1.0, 1.0)), None),
    );
    close3(
        resolved.point,
        Point3::new(50.0, 0.0, 101.0) + u,
        EPS,
        "offset along the base up",
    );
    assert_frame(&resolved, [t, r * (l - u), r * (l + u)], "elevated");
}

#[test]
fn an_axis_against_the_bank_stands_a_banked_section_upright() {
    // Constant 150 mm cant over a 1.5 m gauge on a level straight along
    // +x: sin psi = 0.1, l = (0, cos, sin), u = (0, -sin, cos). The axis
    // (0, sin, cos) in (t, l, u) is world +Z, so the turned frame is
    // (x, y, z) upright while the offset still rolls with the bank.
    let (s, c) = (0.1_f64, (1.0_f64 - 0.01).sqrt());
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push_value(Curve3::Banked(Banked3::new(
            Elevated3::new(
                Curve2::Line(Line2 {
                    origin: Point2::ZERO,
                    direction: Vec2::X,
                }),
                ElevationLaw::level(0.0),
            ),
            CantLaw::new(vec![CantPiece::constant(100.0, 0.15)]),
            CantLaw::zero(100.0),
            1.5,
            BankConvention::TangentRotation,
        )))
        .unwrap();
    let resolved = resolve_oriented(
        b,
        basis,
        Station::new(30.0, StationOffsets::new(0.75, 0.0, 0.0)),
        StationFrame::Section,
        orientation(Some(Vec3::new(0.0, s, c)), None),
    );
    close3(
        resolved.point,
        Point3::new(30.0, 0.75 * c, 0.075),
        EPS,
        "left rail head",
    );
    assert_frame(&resolved, [Vec3::X, Vec3::Y, Vec3::Z], "upright");
}

#[test]
fn parallel_or_degenerate_orientations_are_refused_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let push = |b: &mut GeometryGraphBuilder, axis, reference| {
        b.push_value(OrientedCurveStation::new(
            CurveStation::new(basis, Station::at(1.0)),
            orientation(axis, reference),
        ))
    };
    let refused = |result: Result<NodeId, GraphError>, needle: &str| match result {
        Err(GraphError::InvalidStation { detail }) => {
            assert!(detail.contains(needle), "{detail}");
        }
        other => panic!("expected InvalidStation({needle}), got {other:?}"),
    };
    refused(
        push(&mut b, Some(Vec3::Z), Some(Vec3::new(0.0, 0.0, -3.0))),
        "parallel",
    );
    // Only a reference along the default axis.
    refused(push(&mut b, None, Some(Vec3::Z)), "parallel");
    // Only an axis along the default reference.
    refused(push(&mut b, Some(Vec3::X), None), "parallel");
    refused(push(&mut b, Some(Vec3::ZERO), None), "axis is zero");
    refused(
        push(&mut b, None, Some(Vec3::new(Scalar::NAN, 0.0, 0.0))),
        "reference direction is zero or not finite",
    );
    // A section's orientation is checked the same way.
    let profile = rectangle(&mut b, 1.0, 1.0);
    let spine = SolidOperation::SectionsAtStations {
        directrix: basis,
        sections: vec![
            SectionAtStation::new(profile, Station::at(0.0)),
            SectionAtStation::new(profile, Station::at(1.0))
                .with_orientation(orientation(Some(Vec3::Y), Some(Vec3::Y))),
        ],
        frame: StationFrame::Section,
    };
    refused(b.push(GeometryNode::SolidOperation(spine)), "parallel");
}

// --- oriented spines --------------------------------------------------------

fn rectangle(b: &mut GeometryGraphBuilder, x: Scalar, y: Scalar) -> NodeId {
    b.push_value(Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    }))
    .unwrap()
}

fn sections_spine(
    b: &mut GeometryGraphBuilder,
    directrix: NodeId,
    sections: Vec<SectionAtStation>,
) -> NodeId {
    b.push(GeometryNode::SolidOperation(
        SolidOperation::SectionsAtStations {
            directrix,
            sections,
            frame: StationFrame::Section,
        },
    ))
    .unwrap()
}

#[test]
fn a_leaning_section_sweeps_an_oblique_prism() {
    // Axis (1, 0, 1): up' = (t + u) / sqrt 2, t' = (t - u) / sqrt 2, l' = l.
    // The 1 x 2 section leans 45 degrees towards the tangent, so the prism
    // over 10 m holds A * L * cos 45 = 20 / sqrt 2.
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let profile = rectangle(&mut b, 1.0, 2.0);
    let lean = orientation(Some(Vec3::new(1.0, 0.0, 1.0)), None);
    // Offsets 0.5 left (-x) and 1 up stay in the base frame.
    let offsets = StationOffsets::new(0.5, 1.0, 0.0);
    let at = |d| SectionAtStation::new(profile, Station::new(d, offsets)).with_orientation(lean);
    let root = sections_spine(&mut b, basis, vec![at(0.0), at(10.0)]);
    let graph = b.finish(vec![root]).unwrap();
    let mesh = compile(&graph, root).unwrap();
    close(
        volume(&mesh),
        20.0 * core::f64::consts::FRAC_1_SQRT_2,
        1e-12,
        "oblique prism",
    );
    // The top edge (profile y = 1) leans forward by 1 / sqrt 2 in y.
    let (lo, hi) = bounds(&mesh);
    let r = core::f64::consts::FRAC_1_SQRT_2;
    close3(lo, Point3::new(1.0, 1.0 - r, 1.0 - r), EPS, "min");
    close3(hi, Point3::new(2.0, 11.0 + r, 1.0 + r), EPS, "max");
}

#[test]
fn an_unturned_section_at_stations_matches_the_stationed_spine() {
    // The general form without tags or orientation is the #241 spine.
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let small = rectangle(&mut b, 1.0, 1.0);
    let tall = rectangle(&mut b, 1.0, 3.0);
    let old = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::StationedSpine {
                directrix: basis,
                sections: vec![
                    StationedSection {
                        profile: small,
                        station: Station::at(0.0),
                    },
                    StationedSection {
                        profile: tall,
                        station: Station::new(10.0, StationOffsets::new(0.0, 1.0, 0.0)),
                    },
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let new = sections_spine(
        &mut b,
        basis,
        vec![
            SectionAtStation::new(small, Station::at(0.0)),
            SectionAtStation::new(tall, Station::new(10.0, StationOffsets::new(0.0, 1.0, 0.0))),
        ],
    );
    let graph = b.finish(vec![old, new]).unwrap();
    let (old, new) = (compile(&graph, old).unwrap(), compile(&graph, new).unwrap());
    assert_eq!(old.indices, new.indices);
    for (p, q) in old.positions.iter().zip(&new.positions) {
        close3(*p, *q, 1e-12, "same vertices");
    }
}

#[test]
fn an_orientation_that_flips_between_sections_is_refused_where_it_degenerates() {
    // Axis up at 0 and down at 10 interpolate through zero at 5.
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let profile = rectangle(&mut b, 1.0, 1.0);
    let root = sections_spine(
        &mut b,
        basis,
        vec![
            SectionAtStation::new(profile, Station::at(0.0)),
            SectionAtStation::new(profile, Station::at(10.0))
                .with_orientation(orientation(Some(-Vec3::Z), None)),
        ],
    );
    let graph = b.finish(vec![root]).unwrap();
    let error = compile(&graph, root).unwrap_err();
    assert!(error.to_string().contains("zero or not finite"), "{error}");
}

// --- spines matched by tag --------------------------------------------------

fn line(a: Point2, b: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: a,
            direction: b - a,
        }),
        domain: Interval::new(0.0, 1.0),
        same_sense: true,
    }
}

fn ring(points: &[(Scalar, Scalar)]) -> Contour {
    let points: Vec<Point2> = points.iter().map(|&(x, y)| Point2::new(x, y)).collect();
    Contour::new(
        (0..points.len())
            .map(|k| line(points[k], points[(k + 1) % points.len()]))
            .collect(),
    )
}

fn contour(
    b: &mut GeometryGraphBuilder,
    outer: &[(Scalar, Scalar)],
    holes: &[&[(Scalar, Scalar)]],
) -> NodeId {
    b.push_value(Profile::Contour(ContourProfile {
        outer: ring(outer),
        holes: holes.iter().map(|hole| ring(hole)).collect(),
    }))
    .unwrap()
}

/// A contour through `points` whose segments run against it: each line is
/// parameterised from the next point back, `same_sense` false.
fn contour_against(b: &mut GeometryGraphBuilder, points: &[(Scalar, Scalar)]) -> NodeId {
    let points: Vec<Point2> = points.iter().map(|&(x, y)| Point2::new(x, y)).collect();
    let segments = (0..points.len())
        .map(|k| ProfileSegment {
            same_sense: false,
            ..line(points[(k + 1) % points.len()], points[k])
        })
        .collect();
    b.push_value(Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    }))
    .unwrap()
}

fn tagged(profile: NodeId, station: Station, tags: &[&str]) -> SectionAtStation {
    SectionAtStation::new(profile, station).with_tags(tags.iter().copied())
}

#[test]
fn a_tag_matched_spine_with_reordered_vertices_matches_its_index_matched_twin() {
    // 1 x 1 at 0 to 1 x 3 (lifted 1) at 10, as the index-matched #241
    // spine of two rectangles: volume 20. The tall section is authored
    // clockwise from its north-east corner, its tags following its
    // corners, so only tag matching lines the walls up. Its segments run
    // against the contour, which must not change its vertex order.
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let small = contour(
        &mut b,
        &[(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)],
        &[],
    );
    let tall = contour_against(
        &mut b,
        &[(0.5, 1.5), (0.5, -1.5), (-0.5, -1.5), (-0.5, 1.5)],
    );
    let end = Station::new(10.0, StationOffsets::new(0.0, 1.0, 0.0));
    let by_tag = sections_spine(
        &mut b,
        basis,
        vec![
            tagged(small, Station::at(0.0), &["sw", "se", "ne", "nw"]),
            tagged(tall, end, &["ne", "se", "sw", "nw"]),
        ],
    );
    let rect_small = rectangle(&mut b, 1.0, 1.0);
    let rect_tall = rectangle(&mut b, 1.0, 3.0);
    let by_index = sections_spine(
        &mut b,
        basis,
        vec![
            SectionAtStation::new(rect_small, Station::at(0.0)),
            SectionAtStation::new(rect_tall, end),
        ],
    );
    // The same reordered profiles by index twist the walls.
    let twisted = sections_spine(
        &mut b,
        basis,
        vec![
            SectionAtStation::new(small, Station::at(0.0)),
            SectionAtStation::new(tall, end),
        ],
    );
    let graph = b.finish(vec![by_tag, by_index, twisted]).unwrap();
    let tag_volume = volume(&compile(&graph, by_tag).unwrap());
    close(tag_volume, 20.0, 1e-12, "tag-matched");
    close(
        tag_volume,
        volume(&compile(&graph, by_index).unwrap()),
        1e-12,
        "index-matched twin",
    );
    let (lo, hi) = bounds(&compile(&graph, by_tag).unwrap());
    close3(lo, Point3::new(1.5, 1.0, -0.5), EPS, "start bottom");
    close3(hi, Point3::new(2.5, 11.0, 2.5), EPS, "end top");
    let twisted = compile(&graph, twisted).map(|mesh| volume(&mesh));
    assert!(
        twisted.map_or(true, |v| (v - 20.0).abs() > 1e-6),
        "index matching the reordered profiles must not give the tagged result"
    );
}

#[test]
fn tags_carry_holes_listed_in_another_order() {
    // A 4 x 4 slab with two 1 x 1 voids: the second section lists the
    // voids swapped, starts its outer ring elsewhere and winds a void the
    // other way. Volume (16 - 2) * 10.
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let outer = [(-2.0, -2.0), (2.0, -2.0), (2.0, 2.0), (-2.0, 2.0)];
    // Diagonal voids: the profile triangulator mishandles two voids side by
    // side in one horizontal band (an extrusion of that slab is open too),
    // which is not what this test is about.
    let left = [(-1.5, -1.5), (-1.5, -0.5), (-0.5, -0.5), (-0.5, -1.5)];
    let right = [(0.5, 0.5), (0.5, 1.5), (1.5, 1.5), (1.5, 0.5)];
    let first = contour(&mut b, &outer, &[&left, &right]);
    let outer_rotated = [(2.0, 2.0), (-2.0, 2.0), (-2.0, -2.0), (2.0, -2.0)];
    let right_reversed = [(1.5, 0.5), (1.5, 1.5), (0.5, 1.5), (0.5, 0.5)];
    let second = contour(&mut b, &outer_rotated, &[&right_reversed, &left]);
    let first_tags = [
        "o1", "o2", "o3", "o4", "l1", "l2", "l3", "l4", "r1", "r2", "r3", "r4",
    ];
    let second_tags = [
        "o3", "o4", "o1", "o2", "r4", "r3", "r2", "r1", "l1", "l2", "l3", "l4",
    ];
    let root = sections_spine(
        &mut b,
        basis,
        vec![
            tagged(first, Station::at(0.0), &first_tags),
            tagged(second, Station::at(10.0), &second_tags),
        ],
    );
    let graph = b.finish(vec![root]).unwrap();
    let mesh = compile(&graph, root).unwrap();
    close(volume(&mesh), 140.0, 1e-12, "slab with voids");
}

#[test]
fn a_mirrored_section_is_matched_through_its_transform() {
    // The unit square right of the directrix (profile x in [0, 1]) at both
    // ends, the second authored as the mirror of the square left of it:
    // its corners, wound clockwise before the mirror, keep their tags.
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let right = contour(
        &mut b,
        &[(0.0, -0.5), (1.0, -0.5), (1.0, 0.5), (0.0, 0.5)],
        &[],
    );
    let left = Profile::Contour(ContourProfile {
        outer: ring(&[(0.0, -0.5), (-1.0, -0.5), (-1.0, 0.5), (0.0, 0.5)]),
        holes: Vec::new(),
    });
    let mirrored = b
        .push_value(Profile::Derived {
            basis: Box::new(left),
            transform: axiolid_core::Transform2::from_scale(Vec2::new(-1.0, 1.0)),
        })
        .unwrap();
    let names = ["a", "b", "c", "d"];
    let root = sections_spine(
        &mut b,
        basis,
        vec![
            tagged(right, Station::at(0.0), &names),
            tagged(mirrored, Station::at(10.0), &names),
        ],
    );
    let graph = b.finish(vec![root]).unwrap();
    let mesh = compile(&graph, root).unwrap();
    close(volume(&mesh), 10.0, 1e-12, "unit prism");
    // Profile x runs along the lateral axis, -x in the world: [1, 2].
    let (lo, hi) = bounds(&mesh);
    close3(lo, Point3::new(1.0, 1.0, -0.5), EPS, "min");
    close3(hi, Point3::new(2.0, 11.0, 0.5), EPS, "max");
}

#[test]
fn mismatched_tags_are_refused_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let square = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)];
    let profile = contour(&mut b, &square, &[]);
    let at = |d, tags: &[&str]| tagged(profile, Station::at(d), tags);
    let spine = |sections| {
        GeometryNode::SolidOperation(SolidOperation::SectionsAtStations {
            directrix: basis,
            sections,
            frame: StationFrame::Section,
        })
    };
    let refused = |result: Result<NodeId, GraphError>, needle: &str| match result {
        Err(GraphError::InvalidStation { detail }) => {
            assert!(detail.contains(needle), "{detail}");
        }
        other => panic!("expected InvalidStation({needle}), got {other:?}"),
    };
    // Refused when the graph is built: tag sets and mixed tagging.
    refused(
        b.push(spine(vec![
            at(0.0, &["a", "b", "c", "d"]),
            at(5.0, &["a", "b", "c", "e"]),
        ])),
        "same set of tags",
    );
    refused(
        b.push(spine(vec![
            at(0.0, &["a", "b", "c", "d"]),
            SectionAtStation::new(profile, Station::at(5.0)),
        ])),
        "tagged and others not",
    );
    refused(
        b.push(spine(vec![
            at(0.0, &["a", "b", "c", "c"]),
            at(5.0, &["a", "b", "c", "c"]),
        ])),
        "repeats a tag",
    );
    // Refused when compiled: an order that is no rotation, a tag count
    // that is not the vertex count, a curved contour.
    let crossed = b
        .push(spine(vec![
            at(0.0, &["a", "b", "c", "d"]),
            at(5.0, &["a", "c", "b", "d"]),
        ]))
        .unwrap();
    let short = b
        .push(spine(vec![
            at(0.0, &["a", "b", "c"]),
            at(5.0, &["a", "b", "c"]),
        ]))
        .unwrap();
    let disk = rectangle(&mut b, 1.0, 1.0);
    let parametric = b
        .push(spine(vec![
            tagged(disk, Station::at(0.0), &["a", "b", "c", "d"]),
            tagged(disk, Station::at(5.0), &["a", "b", "c", "d"]),
        ]))
        .unwrap();
    let round = b
        .push_value(Profile::Contour(ContourProfile {
            outer: Contour::new(vec![ProfileSegment {
                curve: Curve2::Circle(Circle2 {
                    frame: Frame2 {
                        origin: Point2::ZERO,
                        x: Vec2::X,
                        y: Vec2::Y,
                    },
                    radius: 1.0,
                }),
                domain: Interval::new(0.0, core::f64::consts::TAU),
                same_sense: true,
            }]),
            holes: Vec::new(),
        }))
        .unwrap();
    let curved = b
        .push(spine(vec![
            tagged(round, Station::at(0.0), &["a"]),
            tagged(round, Station::at(5.0), &["a"]),
        ]))
        .unwrap();
    // A slab with a void whose second section names the void's corners as
    // the outer ring's and back.
    let slab = contour(
        &mut b,
        &[(-2.0, -2.0), (2.0, -2.0), (2.0, 2.0), (-2.0, 2.0)],
        &[&[(-0.5, -0.5), (-0.5, 0.5), (0.5, 0.5), (0.5, -0.5)]],
    );
    let names = ["o1", "o2", "o3", "o4", "h1", "h2", "h3", "h4"];
    let swapped = ["h1", "h2", "h3", "h4", "o1", "o2", "o3", "o4"];
    let inside_out = b
        .push(spine(vec![
            tagged(slab, Station::at(0.0), &names),
            tagged(slab, Station::at(5.0), &swapped),
        ]))
        .unwrap();
    let graph = b
        .finish(vec![crossed, short, parametric, curved, inside_out])
        .unwrap();
    let error = compile(&graph, crossed).unwrap_err();
    assert!(error.to_string().contains("same cyclic order"), "{error}");
    let error = compile(&graph, inside_out).unwrap_err();
    assert!(
        error.to_string().contains("the outer onto the outer"),
        "{error}"
    );
    assert!(matches!(
        compile(&graph, curved).unwrap_err(),
        GeomError::UnsupportedInput { input, .. } if input.contains("curved contour segment")
    ));
    let error = compile(&graph, short).unwrap_err();
    assert!(
        error.to_string().contains("3 tags for 4 contour vertices"),
        "{error}"
    );
    assert!(matches!(
        compile(&graph, parametric).unwrap_err(),
        GeomError::UnsupportedInput { input, .. } if input.contains("not a contour")
    ));
}

// --- open sections ------------------------------------------------------------

fn open_section(b: &mut GeometryGraphBuilder, points: Vec<Point2>) -> NodeId {
    let path = b
        .push_value(Curve2::Polyline(Polyline2 {
            points,
            closed: false,
        }))
        .unwrap();
    b.push_value(OpenProfile::new(path)).unwrap()
}

fn area(mesh: &TriMesh) -> Scalar {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|k| mesh.positions[t[k] as usize]);
            0.5 * (b - a).cross(c - a).length()
        })
        .sum()
}

fn open_sections(
    b: &mut GeometryGraphBuilder,
    directrix: NodeId,
    sections: Vec<SectionAtStation>,
) -> Result<NodeId, GraphError> {
    b.push(GeometryNode::SurfaceRelation(
        SurfaceRelation::OpenSectionsAtStations {
            directrix,
            sections,
            frame: StationFrame::Plan,
        },
    ))
}

#[test]
fn a_reversed_open_section_is_joined_reversed() {
    // The crowned deck of the #241 test, its wide section authored right
    // to left: joined by tag it is the same sheet.
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let narrow = open_section(
        &mut b,
        vec![
            Point2::new(-2.0, -0.04),
            Point2::new(0.0, 0.0),
            Point2::new(2.0, -0.04),
        ],
    );
    let wide_backwards = open_section(
        &mut b,
        vec![
            Point2::new(3.0, -0.06),
            Point2::new(0.0, 0.0),
            Point2::new(-3.0, -0.06),
        ],
    );
    let root = open_sections(
        &mut b,
        basis,
        vec![
            tagged(narrow, Station::at(0.0), &["left", "crown", "right"]),
            tagged(
                wide_backwards,
                Station::at(10.0),
                &["right", "crown", "left"],
            ),
        ],
    )
    .unwrap();
    // Neither the first's order nor its reverse: refused by name.
    let crossed = open_sections(
        &mut b,
        basis,
        vec![
            tagged(narrow, Station::at(0.0), &["left", "crown", "right"]),
            tagged(
                wide_backwards,
                Station::at(10.0),
                &["crown", "right", "left"],
            ),
        ],
    );
    assert!(
        matches!(&crossed, Err(GraphError::InvalidStation { detail }) if detail.contains("in reverse")),
        "{crossed:?}"
    );
    let mixed = open_sections(
        &mut b,
        basis,
        vec![
            tagged(narrow, Station::at(0.0), &["left", "crown", "right"]),
            SectionAtStation::new(wide_backwards, Station::at(10.0)),
        ],
    );
    assert!(
        matches!(&mixed, Err(GraphError::InvalidStation { detail }) if detail.contains("tagged and others not")),
        "{mixed:?}"
    );
    let graph = b.finish(vec![root]).unwrap();
    let outcome = compiler()
        .compile_mesh_with_deviation(&graph, root, &options())
        .unwrap()
        .0;
    assert_eq!(outcome.closure, MeshClosure::Surface);
    let expected = 2.0 * 0.5 * (2.0 + 3.0) * 10.0 * (1.0_f64 + 0.0004).sqrt();
    close(area(&outcome.mesh), expected, 1e-12, "deck area");
}

#[test]
fn an_oriented_open_section_stands_in_its_turned_plane() {
    // Axis (0, 1, 0) is up' = l = -x; t' = t = +y; l' = up' x t' = -Z. A
    // flat section from profile x = -1 to 1 hangs from z = 1 to z = -1 at
    // x = 2, then 2 m wide over 10 m of plan line.
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let flat = open_section(&mut b, vec![Point2::new(-1.0, 0.0), Point2::new(1.0, 0.0)]);
    let turn = orientation(Some(Vec3::Y), None);
    let at = |d| SectionAtStation::new(flat, Station::at(d)).with_orientation(turn);
    let root = open_sections(&mut b, basis, vec![at(0.0), at(10.0)]).unwrap();
    let graph = b.finish(vec![root]).unwrap();
    let mesh = compile(&graph, root).unwrap();
    close(area(&mesh), 20.0, 1e-12, "vertical sheet");
    let (lo, hi) = bounds(&mesh);
    close3(lo, Point3::new(2.0, 1.0, -1.0), EPS, "min");
    close3(hi, Point3::new(2.0, 11.0, 1.0), EPS, "max");
}
