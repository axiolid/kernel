//! Exact extrusions along oblique directions of profiles with arcs (#280).
//!
//! An arc swept along a direction `d` that leans off the profile normal
//! sweeps an oblique circular cylinder: a cylinder whose rulings run along
//! `d` and whose cross-section perpendicular to `d` is an ellipse. Every
//! case here is checked against closed forms, never against a second run
//! of the same builder: a boundary that closes and audits clean, a signed
//! volume of `area * depth * d.z / |d|` (a shear keeps volume), top
//! vertices that are bottom vertices moved by the offset, and curved walls
//! whose points, sheared back to the profile plane, land on the arc's own
//! circle.

use std::f64::consts::{FRAC_PI_2, PI};

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_contracts::{GeomError, Operation};
use axiolid_core::{Frame2, Interval, Point2, Point3, Tolerance, Transform2, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_evaluate::surface::evaluate as evaluate_surface;
use axiolid_measure::exact_properties;
use axiolid_profile::{
    CenterLineProfile, CircleProfile, Contour, ContourProfile, EllipseProfile, Profile,
    ProfileSegment, RectangleProfile, SectionProfile,
};
use axiolid_surface::Surface;
use axiolid_topology::audit_brep;

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("tolerance")
}

fn line(from: Point2, to: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}

fn arc(centre: Point2, radius: f64, from: f64, to: f64) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: centre,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius,
        }),
        domain: Interval::new(from, to),
        same_sense: true,
    }
}

/// A full circle as four quarter arcs (ADR 0053 refuses half turns).
fn circle_contour(centre: Point2, radius: f64) -> Contour {
    Contour::new(
        (0..4)
            .map(|index| {
                arc(
                    centre,
                    radius,
                    FRAC_PI_2 * index as f64,
                    FRAC_PI_2 * (index + 1) as f64,
                )
            })
            .collect(),
    )
}

fn square_contour(half: f64) -> Contour {
    Contour::new(vec![
        line(Point2::new(-half, -half), Point2::new(half, -half)),
        line(Point2::new(half, -half), Point2::new(half, half)),
        line(Point2::new(half, half), Point2::new(-half, half)),
        line(Point2::new(-half, half), Point2::new(-half, -half)),
    ])
}

fn rounded_rectangle() -> Profile {
    Profile::Rectangle(RectangleProfile {
        x: 4.0,
        y: 2.0,
        thickness: None,
        outer_radius: Some(0.5),
        inner_radius: None,
    })
}

fn rounded_rectangle_area() -> f64 {
    8.0 - (4.0 - PI) * 0.25
}

/// The corner circles of a rounded rectangle centred at `at`, of radius
/// `radius` about centres `(+-1.5, +-0.5)` from it.
fn corners(at: Point2, radius: f64) -> Vec<Round> {
    [(1.5, 0.5), (-1.5, 0.5), (-1.5, -0.5), (1.5, -0.5)]
        .into_iter()
        .map(|(x, y)| Round {
            centre: at + Vec2::new(x, y),
            radius,
        })
        .collect()
}

/// A curved wall the test can check point by point: points of a surface of
/// this kind, sheared back to the profile plane, lie on this circle.
#[derive(Clone, Copy)]
struct Round {
    centre: Point2,
    radius: f64,
}

/// A profile family with arcs, its area, and the circles of all its arcs.
struct Family {
    name: &'static str,
    profile: Profile,
    area: f64,
    rounds: Vec<Round>,
}

fn families() -> Vec<Family> {
    let quarter_path = Profile::CenterLine(CenterLineProfile::from_width(
        Contour::new(vec![arc(Point2::new(0.0, 0.0), 2.0, 0.0, FRAC_PI_2)]),
        0.4,
    ));
    // A D: a chord closed by a major arc of 270 degrees.
    let d_shape = Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            arc(Point2::new(0.0, 0.0), 1.0, 0.0, 1.5 * PI),
            line(Point2::new(0.0, -1.0), Point2::new(1.0, 0.0)),
        ]),
        holes: Vec::new(),
    });
    vec![
        Family {
            name: "rounded rectangle",
            profile: rounded_rectangle(),
            area: rounded_rectangle_area(),
            rounds: corners(Point2::ZERO, 0.5),
        },
        Family {
            name: "hollow rounded rectangle",
            profile: Profile::Rectangle(RectangleProfile {
                x: 4.0,
                y: 2.0,
                thickness: Some(0.4),
                outer_radius: Some(0.5),
                inner_radius: Some(0.1),
            }),
            area: rounded_rectangle_area() - (3.2 * 1.2 - (4.0 - PI) * 0.01),
            rounds: [corners(Point2::ZERO, 0.5), corners(Point2::ZERO, 0.1)].concat(),
        },
        Family {
            name: "contour with a round hole",
            profile: Profile::Contour(ContourProfile {
                outer: square_contour(2.0),
                holes: vec![circle_contour(Point2::new(-0.75, 0.25), 0.5)],
            }),
            area: 16.0 - PI * 0.25,
            rounds: vec![Round {
                centre: Point2::new(-0.75, 0.25),
                radius: 0.5,
            }],
        },
        Family {
            name: "round contour",
            profile: Profile::Contour(ContourProfile {
                outer: circle_contour(Point2::new(0.5, -0.25), 1.25),
                holes: Vec::new(),
            }),
            area: PI * 1.5625,
            rounds: vec![Round {
                centre: Point2::new(0.5, -0.25),
                radius: 1.25,
            }],
        },
        Family {
            name: "major arc",
            profile: d_shape,
            area: 0.75 * PI + 0.5,
            rounds: vec![Round {
                centre: Point2::new(0.0, 0.0),
                radius: 1.0,
            }],
        },
        Family {
            name: "curved centre line",
            profile: quarter_path,
            // Width times the centre line's length.
            area: 0.4 * 2.0 * FRAC_PI_2,
            rounds: [1.8, 2.2]
                .into_iter()
                .map(|radius| Round {
                    centre: Point2::ZERO,
                    radius,
                })
                .collect(),
        },
        Family {
            name: "circle",
            profile: Profile::Circle(CircleProfile {
                radius: 1.5,
                thickness: None,
            }),
            area: PI * 2.25,
            rounds: vec![Round {
                centre: Point2::new(0.0, 0.0),
                radius: 1.5,
            }],
        },
        Family {
            name: "derived rounded rectangle",
            profile: Profile::Derived {
                basis: Box::new(rounded_rectangle()),
                transform: Transform2::from_translation(Vec2::new(3.0, -1.0)),
            },
            area: rounded_rectangle_area(),
            rounds: corners(Point2::new(3.0, -1.0), 0.5),
        },
        Family {
            name: "composite of a circle and a rounded rectangle",
            profile: Profile::Composite(vec![
                Profile::Circle(CircleProfile {
                    radius: 0.5,
                    thickness: None,
                }),
                Profile::Derived {
                    basis: Box::new(rounded_rectangle()),
                    transform: Transform2::from_translation(Vec2::new(5.0, 0.0)),
                },
            ]),
            area: PI * 0.25 + rounded_rectangle_area(),
            rounds: [
                vec![Round {
                    centre: Point2::ZERO,
                    radius: 0.5,
                }],
                corners(Point2::new(5.0, 0.0), 0.5),
            ]
            .concat(),
        },
    ]
}

/// Oblique directions, leaning a little, a lot, and nearly flat; each is
/// also taken downward.
fn oblique_directions() -> Vec<Vec3> {
    let up = [
        Vec3::new(0.3, -0.2, 1.0),
        Vec3::new(-1.0, 2.0, 1.5),
        Vec3::new(2.0, 0.5, 0.3),
        Vec3::new(0.0, 1e-3, 1.0),
    ];
    up.iter()
        .flat_map(|d| [*d, Vec3::new(d.x, d.y, -d.z)])
        .collect()
}

fn assert_sound(solid: &ExactBRep, what: &str) {
    assert!(
        audit_brep(solid.topology()).is_closed_manifold(),
        "{what}: not a closed manifold"
    );
    let health = geometric_audit(solid, tol());
    assert!(
        health.is_consistent(),
        "{what}: audit found {:?}",
        health.defects()
    );
}

fn signed_volume(solid: &ExactBRep, what: &str) -> f64 {
    exact_properties(solid, tol())
        .unwrap_or_else(|error| panic!("{what}: unmeasurable: {error:?}"))
        .signed_volume
}

/// Every vertex is on the profile plane or on the far cap, and every far
/// vertex is a profile-plane vertex moved by the offset.
fn assert_sheared_vertices(solid: &ExactBRep, offset: Vec3, what: &str) {
    let points: Vec<Point3> = solid
        .topology()
        .vertices()
        .iter()
        .map(|vertex| vertex.position)
        .collect();
    for p in &points {
        if p.z == 0.0 {
            continue;
        }
        assert!(
            (p.z - offset.z).abs() <= 1e-12,
            "{what}: vertex {p:?} is on neither cap"
        );
        let base = *p - offset;
        assert!(
            points
                .iter()
                .any(|q| q.z == 0.0 && (q.truncate() - base.truncate()).length() <= 1e-12),
            "{what}: far vertex {p:?} is no base vertex moved by {offset:?}"
        );
    }
}

/// Points of every curved wall, sheared back along the offset onto the
/// profile plane, lie on one of the family's known circles.
fn assert_walls_on_rounds(solid: &ExactBRep, offset: Vec3, rounds: &[Round], what: &str) {
    for surface in solid.surfaces() {
        if matches!(surface, Surface::Plane(_)) {
            continue;
        }
        for i in 0..8 {
            for j in 0..3 {
                let (u, v) = (i as f64 * 0.8, j as f64 * 0.3);
                let p = evaluate_surface(surface, u, v).expect("wall point");
                let back = p - offset * (p.z / offset.z);
                let on = rounds.iter().any(|round| {
                    ((back.truncate() - round.centre).length() - round.radius).abs() <= 1e-9
                });
                assert!(
                    on,
                    "{what}: wall point {p:?} shears back to {back:?}, on no profile circle"
                );
            }
        }
    }
}

/// The issue's repro: a rounded rectangle along `(0.3, -0.2, 1)` at depth
/// 0.75. Its arc walls used to be right cylinders standing on the profile
/// plane, so the geometric audit found vertex and pcurve errors of about
/// the horizontal shear.
#[test]
fn the_issue_rounded_rectangle_extrudes_obliquely_and_audits_clean() {
    let direction = Vec3::new(0.3, -0.2, 1.0);
    let depth = 0.75;
    let solid = extrude_profile_exact(&rounded_rectangle(), direction, depth, tol())
        .expect("an oblique rounded rectangle extrudes exactly");
    assert_sound(&solid, "issue repro");
    let offset = direction.normalize() * depth;
    assert_sheared_vertices(&solid, offset, "issue repro");
    let expected = rounded_rectangle_area() * offset.z;
    let volume = signed_volume(&solid, "issue repro");
    assert!(
        (volume - expected).abs() <= 1e-9 * expected,
        "signed volume {volume}, closed form {expected}"
    );
    // The four corner walls are oblique cylinders: elliptical cross-section
    // perpendicular to the direction, axis along it.
    let walls: Vec<_> = solid
        .surfaces()
        .iter()
        .filter_map(|surface| match surface {
            Surface::EllipticalCylinder(wall) => Some(*wall),
            _ => None,
        })
        .collect();
    assert_eq!(walls.len(), 4, "one oblique wall per rounded corner");
    let unit = direction.normalize();
    for wall in walls {
        assert!(
            (wall.frame.z - unit).length() <= 1e-15,
            "axis {:?}",
            wall.frame.z
        );
        assert!((wall.semi_axis_x - 0.5).abs() <= 1e-15);
        assert!((wall.semi_axis_y - 0.5 * unit.z).abs() <= 1e-15);
    }
}

/// Every family with arcs over oblique up and down directions either
/// builds a sound solid of the closed-form volume, or is refused by name.
#[test]
fn every_family_with_arcs_extrudes_obliquely_or_refuses_by_name() {
    let depth = 0.75;
    for Family {
        name,
        profile,
        area,
        rounds,
    } in families()
    {
        for direction in oblique_directions() {
            let what = format!("{name} along {direction:?}");
            let solid = extrude_profile_exact(&profile, direction, depth, tol())
                .unwrap_or_else(|error| panic!("{what}: {error:?}"));
            assert_sound(&solid, &what);
            let offset = direction.normalize() * depth;
            assert_sheared_vertices(&solid, offset, &what);
            assert_walls_on_rounds(&solid, offset, &rounds, &what);
            let expected = area * offset.z.abs();
            let volume = signed_volume(&solid, &what);
            assert!(
                (volume - expected).abs() <= 1e-9 * expected,
                "{what}: signed volume {volume}, closed form {expected}"
            );
        }
    }
}

/// A rolled I section with root fillets: its area has no short closed
/// form, so the oblique prism is held to the straight one of the same
/// height, built along the normal.
#[test]
fn a_filleted_i_section_extrudes_obliquely_with_the_straight_volume() {
    let profile = Profile::Section(SectionProfile::I {
        depth: 0.4,
        width: 0.3,
        web_thickness: 0.011,
        flange_thickness: 0.019,
        fillet_radius: Some(0.021),
        flange_edge_radius: None,
        flange_slope: None,
    });
    for direction in oblique_directions() {
        let what = format!("I section along {direction:?}");
        let solid = extrude_profile_exact(&profile, direction, 2.0, tol())
            .unwrap_or_else(|error| panic!("{what}: {error:?}"));
        assert_sound(&solid, &what);
        let offset = direction.normalize() * 2.0;
        assert_sheared_vertices(&solid, offset, &what);
        let straight = extrude_profile_exact(&profile, Vec3::Z, offset.z.abs(), tol())
            .expect("the straight I section");
        let expected = signed_volume(&straight, "straight I section");
        let volume = signed_volume(&solid, &what);
        assert!(
            (volume - expected).abs() <= 1e-9 * expected,
            "{what}: signed volume {volume}, straight {expected}"
        );
    }
}

/// An ellipse swept obliquely is an elliptical cylinder whose principal
/// axes are not the profile's; that is not built, and is refused by name.
#[test]
fn an_oblique_ellipse_is_refused_by_name() {
    for direction in oblique_directions() {
        let error = extrude_profile_exact(
            &Profile::Ellipse(EllipseProfile {
                semi_axis_x: 2.0,
                semi_axis_y: 1.0,
            }),
            direction,
            0.75,
            tol(),
        )
        .expect_err("an oblique ellipse is refused");
        assert!(
            matches!(
                error,
                GeomError::UnsupportedInput {
                    operation: Operation::Sweep,
                    input: "oblique ellipse extrusion",
                    ..
                }
            ),
            "{direction:?}: {error:?}"
        );
    }
}

/// Along the normal the arc walls stay right circular cylinders: the
/// oblique support is only for a direction that leans.
#[test]
fn along_the_normal_arc_walls_stay_right_cylinders() {
    for profile in [
        rounded_rectangle(),
        Profile::Circle(CircleProfile {
            radius: 1.5,
            thickness: None,
        }),
    ] {
        for direction in [Vec3::Z, -Vec3::Z] {
            let solid = extrude_profile_exact(&profile, direction, 0.75, tol()).expect("straight");
            assert!(solid
                .surfaces()
                .iter()
                .all(|surface| matches!(surface, Surface::Plane(_) | Surface::Cylinder(_))));
        }
    }
}
