//! Exact extrusions whose direction points against the profile normal (#275).
//!
//! The solid along `o` with `o.z < 0` is the mirror image, in the profile
//! plane, of the solid along `(o.x, o.y, -o.z)`. Each test checks that
//! against closed forms and against the forward build, never against a
//! second call of the downward path: a signed volume that is positive (so
//! the faces point outward), a boundary that closes and audits clean, and
//! vertices that are the forward vertices with `z` negated.

use axiolid_brep_audit::geometric_audit;
use axiolid_construct::extrude::{extrude_profile, extrude_profile_exact};
use axiolid_construct::profile::Rings;
use axiolid_contracts::{GeomError, Operation};
use axiolid_core::{Frame2, Interval, Point2, Point3, Tolerance, Transform2, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::exact_properties;
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, EllipseProfile, Profile, ProfileSegment,
    RectangleProfile,
};
use axiolid_topology::audit_brep;

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("tolerance")
}

fn rectangle(x: f64, y: f64, thickness: Option<f64>) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness,
        outer_radius: None,
        inner_radius: None,
    })
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

/// A full circle as four quarter arcs (ADR 0053 refuses half turns).
fn circle_contour(centre: Point2, radius: f64) -> Contour {
    let frame = Frame2 {
        origin: centre,
        x: Vec2::X,
        y: Vec2::Y,
    };
    let quarter = core::f64::consts::FRAC_PI_2;
    Contour::new(
        (0..4)
            .map(|index| ProfileSegment {
                curve: Curve2::Circle(Circle2 { frame, radius }),
                domain: Interval::new(quarter * index as f64, quarter * (index + 1) as f64),
                same_sense: true,
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

/// A 4 x 4 square with a round hole of radius 1 and a square hole of side 1.
fn holed_contour() -> Profile {
    Profile::Contour(ContourProfile {
        outer: square_contour(2.0),
        holes: vec![
            circle_contour(Point2::new(-0.75, 0.0), 0.5),
            Contour::new(vec![
                line(Point2::new(0.5, -0.5), Point2::new(0.5, 0.5)),
                line(Point2::new(0.5, 0.5), Point2::new(1.5, 0.5)),
                line(Point2::new(1.5, 0.5), Point2::new(1.5, -0.5)),
                line(Point2::new(1.5, -0.5), Point2::new(0.5, -0.5)),
            ]),
        ],
    })
}

fn holed_contour_area() -> f64 {
    16.0 - core::f64::consts::PI * 0.25 - 1.0
}

/// A 4 x 4 square with two unit square holes: straight edges only.
fn polygon_with_holes() -> Profile {
    let hole = |x: f64| {
        square_ring(&[
            Point2::new(x, -0.5),
            Point2::new(x, 0.5),
            Point2::new(x + 1.0, 0.5),
            Point2::new(x + 1.0, -0.5),
        ])
    };
    Profile::Contour(ContourProfile {
        outer: square_contour(2.0),
        holes: vec![hole(-1.5), hole(0.5)],
    })
}

/// A profile family, its closed-form area, and whether it extrudes along
/// a direction leaning off the normal.
struct Family {
    name: &'static str,
    profile: Profile,
    area: f64,
    leans: bool,
}

fn curved(name: &'static str, profile: Profile, area: f64) -> Family {
    Family {
        name,
        profile,
        area,
        leans: true,
    }
}

fn straight(name: &'static str, profile: Profile, area: f64) -> Family {
    Family {
        name,
        profile,
        area,
        leans: true,
    }
}

/// The profiles every builder family takes.
fn families() -> Vec<Family> {
    vec![
        straight("rectangle", rectangle(4.0, 2.0, None), 8.0),
        straight(
            "hollow rectangle",
            rectangle(4.0, 2.0, Some(0.25)),
            8.0 - 3.5 * 1.5,
        ),
        straight("polygon with holes", polygon_with_holes(), 14.0),
        curved(
            "rounded rectangle",
            Profile::Rectangle(RectangleProfile {
                x: 4.0,
                y: 2.0,
                thickness: None,
                outer_radius: Some(0.5),
                inner_radius: None,
            }),
            8.0 - (4.0 - core::f64::consts::PI) * 0.25,
        ),
        curved(
            "circle",
            Profile::Circle(CircleProfile {
                radius: 1.5,
                thickness: None,
            }),
            core::f64::consts::PI * 2.25,
        ),
        // An ellipse refuses a leaning direction by name.
        Family {
            leans: false,
            ..curved(
                "ellipse",
                Profile::Ellipse(EllipseProfile {
                    semi_axis_x: 2.0,
                    semi_axis_y: 1.0,
                }),
                core::f64::consts::PI * 2.0,
            )
        },
        curved("contour with holes", holed_contour(), holed_contour_area()),
        straight(
            "derived",
            Profile::Derived {
                basis: Box::new(rectangle(2.0, 1.0, None)),
                transform: Transform2::from_translation(Vec2::new(3.0, -1.0)),
            },
            2.0,
        ),
        straight(
            "composite",
            Profile::Composite(vec![
                rectangle(1.0, 1.0, None),
                Profile::Derived {
                    basis: Box::new(rectangle(1.0, 2.0, None)),
                    transform: Transform2::from_translation(Vec2::new(4.0, 0.0)),
                },
            ]),
            3.0,
        ),
    ]
}

fn assert_sound(solid: &axiolid_brep::ExactBRep, what: &str) {
    let topology = solid.topology();
    assert!(
        audit_brep(topology).is_closed_manifold(),
        "{what}: not a closed manifold"
    );
    let health = geometric_audit(solid, tol());
    assert!(
        health.is_consistent(),
        "{what}: audit found {:?}",
        health.defects()
    );
}

fn signed_volume(solid: &axiolid_brep::ExactBRep, what: &str) -> f64 {
    exact_properties(solid, tol())
        .unwrap_or_else(|error| panic!("{what}: unmeasurable: {error:?}"))
        .signed_volume
}

fn sorted_positions(solid: &axiolid_brep::ExactBRep, mirror: bool) -> Vec<Point3> {
    let mut points: Vec<Point3> = solid
        .topology()
        .vertices()
        .iter()
        .map(|vertex| {
            let p = vertex.position;
            if mirror {
                Point3::new(p.x, p.y, -p.z)
            } else {
                p
            }
        })
        .collect();
    points.sort_by(|a, b| {
        a.x.total_cmp(&b.x)
            .then(a.y.total_cmp(&b.y))
            .then(a.z.total_cmp(&b.z))
    });
    points
}

/// The issue's repro, at construct level: a rectangle extruded along `-z`
/// is the prism below the profile plane, outward, of volume `x y depth`.
#[test]
fn a_rectangle_extrudes_down_the_normal() {
    let solid = extrude_profile_exact(&rectangle(2.0, 1.0, None), -Vec3::Z, 0.25, tol())
        .expect("a downward rectangle extrudes exactly");
    assert_sound(&solid, "downward rectangle");
    let volume = signed_volume(&solid, "downward rectangle");
    assert!((volume - 0.5).abs() < 1e-12, "signed volume {volume}");
    // Every vertex is on the profile plane or one depth below it.
    for vertex in solid.topology().vertices() {
        let z = vertex.position.z;
        assert!(z == 0.0 || z == -0.25, "vertex at z = {z}");
    }
}

/// Every builder family along the normal, and every one but the ellipse
/// along oblique directions too: the downward solid is the forward solid
/// mirrored in the profile plane, closed, consistent, and of the same
/// positive closed-form volume. Oblique arc walls are oblique cylinders
/// since #280; an ellipse refuses an oblique direction.
#[test]
fn every_family_extrudes_down_as_the_mirror_of_its_forward_prism() {
    let depth = 0.75;
    for Family {
        name,
        profile,
        area,
        leans,
    } in families()
    {
        let mut directions = vec![Vec3::Z];
        if leans {
            directions.extend([Vec3::new(0.3, -0.2, 1.0), Vec3::new(-1.0, 2.0, 1.5)]);
        }
        for up in directions {
            let down = Vec3::new(up.x, up.y, -up.z);
            let what = format!("{name} along {down:?}");
            let forward = extrude_profile_exact(&profile, up, depth, tol())
                .unwrap_or_else(|error| panic!("{name} along {up:?}: {error:?}"));
            assert_sound(&forward, &format!("{name} along {up:?}"));
            let solid = extrude_profile_exact(&profile, down, depth, tol())
                .unwrap_or_else(|error| panic!("{what}: {error:?}"));
            assert_sound(&solid, &what);

            // Closed form: the prism's height is depth times the normal
            // component of the unit direction.
            let height = depth * up.z / up.length();
            let expected = area * height;
            let up_volume = signed_volume(&forward, &format!("{name} along {up:?}"));
            let down_volume = signed_volume(&solid, &what);
            assert!(
                down_volume > 0.0,
                "{what}: signed volume {down_volume} is not outward"
            );
            assert!(
                (down_volume - expected).abs() <= 1e-9 * expected,
                "{what}: signed volume {down_volume}, closed form {expected}"
            );
            assert!(
                (down_volume - up_volume).abs() <= 1e-9 * expected,
                "{what}: mirror volumes differ, {down_volume} vs {up_volume}"
            );

            // Mirror image, bit for bit: negating z is exact in f64.
            assert_eq!(
                sorted_positions(&solid, true),
                sorted_positions(&forward, false),
                "{what}: vertices are not the forward vertices mirrored"
            );
            assert_eq!(
                solid.topology().faces().len(),
                forward.topology().faces().len(),
                "{what}: face count"
            );
        }
    }
}

/// The mesh path already flips its winding below the plane; the exact path
/// must build the same solid, so the two measure the same volume.
#[test]
fn the_exact_and_mesh_paths_agree_below_the_plane() {
    let rings = Rings {
        outer: vec![
            Point2::new(-2.0, -1.0),
            Point2::new(2.0, -1.0),
            Point2::new(2.0, 1.0),
            Point2::new(-2.0, 1.0),
        ],
        holes: vec![vec![
            Point2::new(-0.5, -0.5),
            Point2::new(-0.5, 0.5),
            Point2::new(0.5, 0.5),
            Point2::new(0.5, -0.5),
        ]],
    };
    let profile = Profile::Contour(ContourProfile {
        outer: square_ring(&rings.outer),
        holes: vec![square_ring(&rings.holes[0])],
    });
    let direction = Vec3::new(0.25, 0.5, -1.0);
    let mesh = extrude_profile(&rings, direction, 3.0, tol()).expect("mesh extrusion");
    let exact = extrude_profile_exact(&profile, direction, 3.0, tol()).expect("exact extrusion");
    let mesh_volume = mesh_signed_volume(&mesh);
    let exact_volume = signed_volume(&exact, "exact");
    let expected = 7.0 * 3.0 / direction.length();
    assert!((mesh_volume - expected).abs() < 1e-9, "mesh {mesh_volume}");
    assert!(
        (exact_volume - expected).abs() < 1e-9,
        "exact {exact_volume}"
    );
    // Same solid: every exact vertex is a mesh vertex.
    for vertex in exact.topology().vertices() {
        assert!(
            mesh.positions
                .iter()
                .any(|p| p.distance(vertex.position) < 1e-12),
            "exact vertex {:?} is not on the mesh",
            vertex.position
        );
    }
}

fn square_ring(points: &[Point2]) -> Contour {
    Contour::new(
        (0..points.len())
            .map(|i| line(points[i], points[(i + 1) % points.len()]))
            .collect(),
    )
}

fn mesh_signed_volume(mesh: &axiolid_mesh::TriMesh) -> f64 {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|k| mesh.positions[t[k] as usize]);
            a.dot(b.cross(c)) / 6.0
        })
        .sum()
}

/// A direction within tolerance of the profile plane bounds no volume on
/// either side, and is refused by name.
#[test]
fn a_direction_in_the_profile_plane_is_refused_by_name() {
    let tolerance = Tolerance::METRE;
    let within = tolerance.linear() * 0.5;
    for direction in [
        Vec3::X,
        Vec3::new(1.0, 1.0, 0.0),
        Vec3::new(1.0, 0.0, within),
        Vec3::new(1.0, 0.0, -within),
        Vec3::new(0.0, 1.0, -tolerance.linear()),
    ] {
        for profile in [rectangle(2.0, 1.0, None), holed_contour()] {
            let error = extrude_profile_exact(&profile, direction, 1.0, tolerance)
                .expect_err("an in-plane direction must refuse");
            assert!(
                matches!(
                    error,
                    GeomError::UnsupportedInput {
                        operation: Operation::Sweep,
                        input: "extrusion direction in the profile plane",
                        ..
                    }
                ),
                "{direction:?}: {error:?}"
            );
        }
    }
    // Just past the tolerance, on either side, it builds.
    for z in [2.0 * tolerance.linear(), -2.0 * tolerance.linear()] {
        let solid = extrude_profile_exact(
            &rectangle(2.0, 1.0, None),
            Vec3::new(1.0, 0.0, z),
            1.0,
            tolerance,
        )
        .expect("a direction past tolerance extrudes");
        assert!(signed_volume(&solid, "near plane") > 0.0);
    }
}

/// An ellipse along `-z` still refuses an oblique direction rather than
/// mislabeling an elliptical cylinder.
#[test]
fn an_oblique_downward_ellipse_still_refuses() {
    let error = extrude_profile_exact(
        &Profile::Ellipse(EllipseProfile {
            semi_axis_x: 2.0,
            semi_axis_y: 1.0,
        }),
        Vec3::new(1.0, 0.0, -1.0),
        1.0,
        tol(),
    )
    .expect_err("oblique ellipse");
    assert!(matches!(
        error,
        GeomError::UnsupportedInput {
            input: "oblique ellipse extrusion",
            ..
        }
    ));
}

/// The downward cap at the profile plane keeps the profile's coordinates
/// bit for bit, so a coplanar cut from a slab's top meets it exactly.
#[test]
fn the_profile_plane_cap_keeps_the_profile_coordinates() {
    let solid = extrude_profile_exact(&holed_contour(), Vec3::new(0.1, 0.2, -1.0), 0.5, tol())
        .expect("downward holed contour");
    let on_plane: Vec<Point3> = solid
        .topology()
        .vertices()
        .iter()
        .map(|vertex| vertex.position)
        .filter(|p| p.z == 0.0)
        .collect();
    assert!(!on_plane.is_empty());
    for corner in [
        Point3::new(-2.0, -2.0, 0.0),
        Point3::new(2.0, 2.0, 0.0),
        Point3::new(0.5, -0.5, 0.0),
        Point3::new(1.5, 0.5, 0.0),
    ] {
        assert!(on_plane.contains(&corner), "missing {corner:?}");
    }
}
