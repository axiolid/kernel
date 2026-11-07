//! Booleans far from the origin (#418, #419).
//!
//! The #356 column with a corbel on each of its four side faces, turned
//! about `z` and placed at a site `6e5, 5.6e6, 200` from the origin. At the
//! origin every union sews into one solid. Far out, the operands'
//! coordinates carry a nanometre ulp, and before the exact recentring
//! (`recentre.rs`) most of these were refused: "kept faces do not sew"
//! (#418) or "a point too close to a face boundary to classify" (#419).
//!
//! Operands whose planar faces are charted in world coordinates (#227: a
//! contour with absolute coordinates, `boolean_arc_prisms_exact` sections)
//! keep parameters as far out as before the move; the recentred cut of
//! those is refused and the operands are cut as given, as without it.

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::{boolean_with_report, BooleanError, BooleanReport, ROUNDING_FACTOR};
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{BooleanOperator, Interval, Point2, Tolerance, Transform3, Vec3};
use axiolid_curve::{Curve2, Line2};
use axiolid_measure::exact_properties;
use axiolid_overlay::ArcRing;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile};

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn rect(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// The corbel's centre in the column's frame and its profile size.
const CORBEL_FACES: [([f64; 2], [f64; 2]); 4] = [
    ([-0.325, 0.0], [0.25, 0.2]),
    ([0.325, 0.0], [0.25, 0.2]),
    ([0.0, -0.375], [0.2, 0.25]),
    ([0.0, 0.375], [0.2, 0.25]),
];

const FAR: Vec3 = Vec3::new(6e5, 5.6e6, 200.0);

/// Sites the column and corbel are placed at: far in every axis, far in
/// `x` and `y` with a height near zero (`z` in `[0.3, 6.3]` and `[-0.3,
/// 5.7]`), and off the grid of the ulp.
const SITES: [Vec3; 4] = [
    FAR,
    Vec3::new(6e5, 5.6e6, 0.3),
    Vec3::new(6e5, 5.6e6, -0.3),
    Vec3::new(6e5 + 0.1, 5.6e6 + 0.37, 3.15),
];

/// The union of the column and the corbel on `face`, the corbel's turn
/// composed from two half turns.
fn union(theta: f64, face: usize, site: Vec3) -> Result<(ExactBRep, BooleanReport), BooleanError> {
    union_at(theta, face, site, tol())
}

/// [`union`] at `tolerance`.
fn union_at(
    theta: f64,
    face: usize,
    site: Vec3,
    tolerance: Tolerance,
) -> Result<(ExactBRep, BooleanReport), BooleanError> {
    let ([cx, cy], [sx, sy]) = CORBEL_FACES[face];
    let column = extrude_profile_exact(&rect(0.4, 0.5), Vec3::Z, 6.0, tol()).expect("a column");
    let corbel = extrude_profile_exact(&rect(sx, sy), Vec3::Z, 0.5, tol()).expect("a corbel");
    let half = Transform3::from_rotation_z(theta / 2.0);
    let column = column
        .transformed(&(Transform3::from_translation(site) * Transform3::from_rotation_z(theta)))
        .expect("placed");
    let corbel = corbel
        .transformed(
            &(Transform3::from_translation(site)
                * half
                * half
                * Transform3::from_translation(Vec3::new(cx, cy, 3.5))),
        )
        .expect("placed");
    boolean_with_report(&column, &corbel, BooleanOperator::Union, tolerance)
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, tol())
        .expect("measurable")
        .signed_volume
}

#[test]
fn a_union_far_from_the_origin_sews_into_one_solid() {
    let mut refused: Vec<String> = Vec::new();
    for site in SITES {
        for face in 0..CORBEL_FACES.len() {
            for theta in [0.001, 0.5, 2.0] {
                let what = format!("{site:?} {face} {theta}");
                match union(theta, face, site) {
                    Ok((brep, _)) => {
                        let health = geometric_audit(&brep, Tolerance::METRE);
                        assert!(health.is_consistent(), "{what}: {:?}", health.defects());
                        assert_eq!(brep.topology().shells().len(), 1, "{what}");
                        // The result is where the operands are, not recentred.
                        assert!(
                            brep.topology()
                                .vertices()
                                .iter()
                                .all(|v| (v.position - site).length() < 10.0),
                            "{what}: the result left the site"
                        );
                        let measured = volume(&brep);
                        assert!(
                            (measured - 1.225).abs() <= 1e-6,
                            "{what}: volume {measured}"
                        );
                    }
                    Err(error) => refused.push(format!("{what}: {error}")),
                }
            }
        }
    }
    assert!(refused.is_empty(), "{refused:#?}");
}

/// A square block with a contour of absolute coordinates under an identity
/// placement (#227), raised to `z0`.
fn contour_block(x0: f64, y0: f64, x1: f64, y1: f64, z0: f64, height: f64) -> ExactBRep {
    let points = [
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ];
    let segments = (0..4)
        .map(|i| {
            let (from, to) = (points[i], points[(i + 1) % 4]);
            ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: from,
                    direction: to - from,
                }),
                domain: Interval::new(0.0, 1.0),
                same_sense: true,
            }
        })
        .collect();
    let profile = Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    });
    extrude_profile_exact(&profile, Vec3::Z, height, tol())
        .expect("a block")
        .transformed(&Transform3::from_translation(Vec3::new(0.0, 0.0, z0)))
        .expect("placed")
}

const OPERATORS: [BooleanOperator; 3] = [
    BooleanOperator::Union,
    BooleanOperator::Intersection,
    BooleanOperator::Difference,
];

/// The volume of `result`, or `None` for a result with nothing left.
fn outcome(result: Result<(ExactBRep, BooleanReport), BooleanError>, what: &str) -> Option<f64> {
    match result {
        Ok((brep, _)) => {
            let health = geometric_audit(&brep, tol());
            assert!(health.is_consistent(), "{what}: {:?}", health.defects());
            Some(volume(&brep))
        }
        Err(BooleanError::EmptyResult) => None,
        Err(error) => panic!("{what}: refused: {error}"),
    }
}

/// Two 2 x 2 x 1 blocks offset by (1, 1), flush or stacked, charted in
/// world coordinates: recentred, they are refused (the classifier cannot
/// decide parameters of 6e5), and are cut as given, as at the origin.
#[test]
fn world_charted_contours_far_from_the_origin_are_cut_as_given() {
    // flush, stacked: union, intersection, difference.
    let expected = [
        [Some(7.0), Some(1.0), Some(3.0)],
        [Some(8.0), None, Some(4.0)],
    ];
    for (x, y, z) in [
        (0.0, 0.0, 0.0),
        (6e5, 5.6e6, 200.0),
        (6e5 + 0.1, 5.6e6 + 0.3, 200.0),
    ] {
        let flush = (
            contour_block(x, y, x + 2.0, y + 2.0, z, 1.0),
            contour_block(x + 1.0, y + 1.0, x + 3.0, y + 3.0, z, 1.0),
        );
        let stacked = (
            contour_block(x, y, x + 2.0, y + 2.0, z, 1.0),
            contour_block(x + 1.0, y + 1.0, x + 3.0, y + 3.0, z + 1.0, 1.0),
        );
        for (row, (a, b)) in [flush, stacked].iter().enumerate() {
            for (column, operator) in OPERATORS.into_iter().enumerate() {
                let what = format!("({x}, {y}, {z}) case {row} {operator:?}");
                let got = outcome(boolean_with_report(a, b, operator, tol()), &what);
                match (got, expected[row][column]) {
                    (Some(got), Some(want)) => {
                        assert!((got - want).abs() <= 1e-6, "{what}: volume {got}");
                    }
                    (None, None) => {}
                    (got, want) => panic!("{what}: {got:?}, expected {want:?}"),
                }
            }
        }
    }
}

fn prism(section: ArcRing, bottom: f64, top: f64) -> ExactBRep {
    let prism = ArcPrism {
        section,
        bottom,
        top,
    };
    boolean_arc_prisms_exact(&prism, &prism, BooleanOperator::Intersection, tol()).expect("a solid")
}

fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> ArcRing {
    ArcRing::from_points(&[
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ])
}

/// Operands in contact, built with their sections in world coordinates at
/// the far site, cut as at the origin: the same volume, or the same
/// refusal.
#[test]
fn contact_far_from_the_origin_is_decided_as_at_the_origin() {
    type Build = fn(f64, f64, f64) -> (ExactBRep, ExactBRep);
    let cases: [(&str, Build); 5] = [
        ("stacked offset", |x, y, z| {
            (
                prism(square(x, y, x + 2.0, y + 2.0), z, z + 1.0),
                prism(square(x + 1.0, y + 1.0, x + 3.0, y + 3.0), z + 1.0, z + 2.0),
            )
        }),
        ("coplanar floors", |x, y, z| {
            (
                prism(square(x, y, x + 2.0, y + 2.0), z, z + 1.0),
                prism(square(x + 1.0, y + 1.0, x + 3.0, y + 3.0), z, z + 1.0),
            )
        }),
        ("touching face", |x, y, z| {
            (
                prism(square(x, y, x + 2.0, y + 2.0), z, z + 1.0),
                prism(square(x + 2.0, y, x + 3.0, y + 2.0), z, z + 1.0),
            )
        }),
        ("edge touch", |x, y, z| {
            (
                prism(square(x, y, x + 1.0, y + 1.0), z, z + 1.0),
                prism(square(x + 1.0, y + 1.0, x + 2.0, y + 2.0), z, z + 1.0),
            )
        }),
        ("gap 5e-7", |x, y, z| {
            (
                prism(square(x, y, x + 2.0, y + 2.0), z, z + 1.0),
                prism(square(x + 2.0 + 5e-7, y, x + 3.0, y + 2.0), z, z + 1.0),
            )
        }),
    ];
    for (name, build) in cases {
        for operator in OPERATORS {
            let (a, b) = build(0.0, 0.0, 0.0);
            let what = format!("{name} {operator:?}");
            let near = boolean_with_report(&a, &b, operator, tol());
            let near_kind = near.as_ref().err().cloned();
            let near = near.ok().map(|(brep, _)| volume(&brep));
            for [x, y, z] in [FAR.to_array(), [6e5 + 0.1, 5.6e6 + 0.3, 200.0]] {
                let (a, b) = build(x, y, z);
                match boolean_with_report(&a, &b, operator, tol()) {
                    Ok((brep, _)) => {
                        let got = volume(&brep);
                        let want =
                            near.unwrap_or_else(|| panic!("{what}: far {got}, near refused"));
                        assert!((got - want).abs() <= 1e-6, "{what} ({x}): {got} vs {want}");
                    }
                    Err(error) => assert!(
                        near.is_none() && near_kind.as_ref() == Some(&error),
                        "{what} ({x}): refused as {error}, near: {near:?}"
                    ),
                }
            }
        }
    }
}

/// The column and the corbel on its +x face of #445, placed in world
/// coordinates at `site`: the corbel's bottom corners sit on the column's
/// face x = site.x + 0.2.
fn world_column_and_corbel(site: Vec3) -> (ExactBRep, ExactBRep) {
    let (x, y, z) = (site.x, site.y, site.z);
    (
        contour_block(x - 0.2, y - 0.25, x + 0.2, y + 0.25, z, 6.0),
        contour_block(x + 0.2, y - 0.1, x + 0.45, y + 0.1, z + 3.5, 0.5),
    )
}

/// The largest distance of an edge end from its curve.
fn vertex_residue(brep: &ExactBRep) -> f64 {
    let topology = brep.topology();
    let mut worst: f64 = 0.0;
    for (index, edge) in topology.edges().iter().enumerate() {
        let curve = &brep.curves3()[edge.curve.expect("a curve").index()];
        let span = brep
            .edge_interval(topology.edge_id_at(index).expect("an edge"))
            .expect("an interval");
        for (parameter, vertex) in [(span.start, edge.start), (span.end, edge.end)] {
            let point = axiolid_evaluate::evaluate3(curve, parameter).expect("evaluates");
            let at = topology.vertices()[vertex.index()].position;
            worst = worst.max((point - at).length());
        }
    }
    worst
}

/// #445: the corbel's bottom corner is the corbel's own vertex. A window
/// crossing 1.2e-6 past the column's face used to replace the feature cut
/// there, so the difference returned a vertex off its edge's curve with an
/// exact report and the union was refused as not sewing. The corners now
/// come out bit for bit, at every tolerance.
#[test]
fn a_world_charted_corbel_far_from_the_origin_is_cut_exactly() {
    let sites = [
        Vec3::ZERO,
        FAR,
        Vec3::new(6e5 + 0.1, 5.6e6 + 0.3, 200.0),
        Vec3::new(6e5, 5.6e6, 0.3),
    ];
    let tolerances = [
        Tolerance::METRE,
        Tolerance::new(1e-7, 1e-9).expect("a tolerance"),
        Tolerance::ZERO,
    ];
    for site in sites {
        let (column, corbel) = world_column_and_corbel(site);
        for tolerance in tolerances {
            let cases = [
                (
                    "union",
                    &column,
                    &corbel,
                    BooleanOperator::Union,
                    Some(1.225),
                ),
                (
                    "union swapped",
                    &corbel,
                    &column,
                    BooleanOperator::Union,
                    Some(1.225),
                ),
                (
                    "intersection",
                    &column,
                    &corbel,
                    BooleanOperator::Intersection,
                    None,
                ),
                (
                    "intersection swapped",
                    &corbel,
                    &column,
                    BooleanOperator::Intersection,
                    None,
                ),
                (
                    "column - corbel",
                    &column,
                    &corbel,
                    BooleanOperator::Difference,
                    Some(1.2),
                ),
                (
                    "corbel - column",
                    &corbel,
                    &column,
                    BooleanOperator::Difference,
                    Some(0.025),
                ),
            ];
            for (name, a, b, operator, want) in cases {
                let what = format!("{site:?} {tolerance:?} {name}");
                let result = boolean_with_report(a, b, operator, tolerance);
                let Some(want) = want else {
                    assert!(
                        matches!(result, Err(BooleanError::EmptyResult)),
                        "{what}: {:?}",
                        result.map(|r| r.1)
                    );
                    continue;
                };
                let (brep, report) = result.unwrap_or_else(|e| panic!("{what}: refused: {e}"));
                assert!(report.is_exact(), "{what}: {report:?}");
                assert_eq!(brep.topology().shells().len(), 1, "{what}");
                if tolerance.linear() > 0.0 {
                    let health = geometric_audit(&brep, Tolerance::METRE);
                    assert!(health.is_consistent(), "{what}: {:?}", health.defects());
                }
                let measured = exact_properties(&brep, tolerance)
                    .expect("measurable")
                    .signed_volume;
                assert!((measured - want).abs() <= 1e-6, "{what}: volume {measured}");
                let residue = vertex_residue(&brep);
                assert!(
                    residue <= ROUNDING_FACTOR * report.extent(),
                    "{what}: vertex {residue:e} off its curve"
                );
                if name == "column - corbel" {
                    // The imprint's corners are the corbel's own vertices.
                    let kept: Vec<_> = brep
                        .topology()
                        .vertices()
                        .iter()
                        .map(|v| v.position)
                        .collect();
                    for v in corbel.topology().vertices() {
                        if v.position.x == site.x + 0.2 {
                            assert!(kept.contains(&v.position), "{what}: {:?}", v.position);
                        }
                    }
                }
            }
        }
    }
}

/// A tolerance below what the coordinates resolve (one ulp at 5.6e6 is
/// 9.3e-10) is refused by name, not answered with a vertex off its curve;
/// at the origin the same cut is fine.
#[test]
fn a_tolerance_finer_than_far_coordinates_resolve_is_refused_by_name() {
    let fine = Tolerance::new(1e-10, 1e-9).expect("a tolerance");
    assert!(matches!(
        union_at(0.001, 0, FAR, fine),
        Err(BooleanError::VertexOffCurve)
    ));
    let (brep, _) = union_at(0.001, 0, Vec3::ZERO, fine).expect("fine at the origin");
    let health = geometric_audit(&brep, fine);
    assert!(health.is_consistent(), "{:?}", health.defects());
}
