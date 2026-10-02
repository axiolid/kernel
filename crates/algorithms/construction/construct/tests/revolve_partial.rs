//! Exact partial-turn revolution (#172), checked against Pappus.
//!
//! Pappus's theorems hold for any swept angle `theta`, not only a full turn:
//! the volume is `theta * R * A` (`A` the section area, `R` its centroid's
//! distance from the axis) and the swept lateral area is
//! `theta * integral(r ds)` round the section boundary. A partial turn adds
//! two planar end caps, each of area `A`. Every expected value below is
//! derived from the profile's dimensions, never recorded from a run.

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_contracts::GeomError;
use axiolid_core::{Frame2, Interval, Point2, Point3, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::MassProperties;
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, EllipseProfile, Profile, ProfileSegment,
    RectangleProfile,
};
use axiolid_surface::Surface;

const PI: f64 = std::f64::consts::PI;
const TAU: f64 = std::f64::consts::TAU;

/// Quarter, half, three-quarter, and two small turns.
const ANGLES: [f64; 5] = [PI / 2.0, PI, 1.5 * PI, PI / 180.0, 1e-3];

fn rect(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
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

fn contour(segments: Vec<ProfileSegment>) -> Profile {
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

fn try_revolve(
    profile: &Profile,
    axis_x: f64,
    direction: Vec3,
    angle: f64,
) -> Result<ExactBRep, GeomError> {
    revolve_profile_exact(
        profile,
        Point3::new(axis_x, 0.0, 0.0),
        direction,
        angle,
        Tolerance::METRE,
    )
}

/// Revolve, then require a clean geometric audit and a closed manifold.
fn revolve(profile: &Profile, axis_x: f64, angle: f64) -> ExactBRep {
    revolve_about(profile, axis_x, Vec3::Y, angle)
}

fn revolve_about(profile: &Profile, axis_x: f64, direction: Vec3, angle: f64) -> ExactBRep {
    let solid = try_revolve(profile, axis_x, direction, angle)
        .unwrap_or_else(|error| panic!("revolves through {angle}: {error:?}"));
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "partial revolution through {angle} must audit clean: {:?}",
        health.defects()
    );
    let topology = axiolid_topology::audit_brep(solid.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
    solid
}

fn measure(solid: &ExactBRep) -> MassProperties {
    axiolid_measure::exact_properties(solid, Tolerance::METRE).expect("measurable")
}

fn close(what: &str, got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-9 * want.abs().max(1e-12),
        "{what}: expected {want}, got {got}"
    );
}

/// Volume `theta R A` and area `theta L_r + 2 A`, where `L_r` is the
/// boundary's `integral(r ds)`.
fn check_pappus(solid: &ExactBRep, angle: f64, area: f64, centroid_radius: f64, moment: f64) {
    let props = measure(solid);
    close(
        "volume",
        props.signed_volume,
        angle * centroid_radius * area,
    );
    close("area", props.area, angle * moment + 2.0 * area);
}

#[test]
fn a_rectangle_matches_pappus_at_every_angle() {
    // 2 x 3 rectangle, axis 5 to its left: r in [4, 6], h = 3.
    let (a, b, h) = (4.0, 6.0, 3.0);
    for angle in ANGLES {
        let solid = revolve(&rect(2.0, 3.0), -5.0, angle);
        // r ds: the two walls a*h and b*h, the two annular sectors
        // (b^2 - a^2) / 2 each.
        let moment = a * h + b * h + (b * b - a * a);
        check_pappus(&solid, angle, (b - a) * h, 5.0, moment);
        // Two cylinder walls, two planar sectors, two planar caps.
        let count =
            |want: fn(&Surface) -> bool| solid.surfaces().iter().filter(|s| want(s)).count();
        assert_eq!(count(|s| matches!(s, Surface::Cylinder(_))), 2);
        assert_eq!(count(|s| matches!(s, Surface::Plane(_))), 4);
        assert_eq!(solid.topology().faces().len(), 6);
    }
}

#[test]
fn a_circle_matches_pappus_at_every_angle() {
    let (rho, major) = (1.25, 4.0);
    let profile = Profile::Circle(CircleProfile {
        radius: rho,
        thickness: None,
    });
    for angle in ANGLES {
        let solid = revolve(&profile, -major, angle);
        let tori = solid
            .surfaces()
            .iter()
            .filter(|s| matches!(s, Surface::Torus(_)))
            .count();
        assert_eq!(tori, 4, "four quarter arcs, four torus walls");
        // r ds round a circle about its centre: 2 pi rho * major.
        check_pappus(&solid, angle, PI * rho * rho, major, TAU * rho * major);
    }
}

#[test]
fn a_half_disk_with_an_arc_matches_pappus() {
    // Flat side on x = R, arc bulging away from the axis: two quarter arcs,
    // since one profile segment stays below half a turn.
    let (rho, big) = (1.0, 3.0);
    let quarter = |from: f64| ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(big, 0.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: rho,
        }),
        domain: Interval::new(from, from + PI / 2.0),
        same_sense: true,
    };
    let profile = contour(vec![
        quarter(-PI / 2.0),
        quarter(0.0),
        line(Point2::new(big, rho), Point2::new(big, -rho)),
    ]);
    let area = PI * rho * rho / 2.0;
    let centroid = big + 4.0 * rho / (3.0 * PI);
    // Flat side 2 rho at r = R; the arc integrates rho (pi R + 2 rho).
    let moment = 2.0 * rho * big + rho * (PI * big + 2.0 * rho);
    for angle in ANGLES {
        let solid = revolve(&profile, 0.0, angle);
        check_pappus(&solid, angle, area, centroid, moment);
    }
}

#[test]
fn a_hollow_rectangle_keeps_its_opening_as_a_tunnel() {
    // A 4 x 2 rectangle with wall thickness 0.5, axis 5 to its left: the
    // opening sweeps a tunnel from cap to cap, not a void.
    let profile = Profile::Rectangle(RectangleProfile {
        x: 4.0,
        y: 2.0,
        thickness: Some(0.5),
        outer_radius: None,
        inner_radius: None,
    });
    let (outer, inner) = (4.0 * 2.0, 3.0 * 1.0);
    // Outer r in [3, 7], h = 2; opening r in [3.5, 6.5], h = 1.
    let ring = |a: f64, b: f64, h: f64| a * h + b * h + (b * b - a * a);
    let moment = ring(3.0, 7.0, 2.0) + ring(3.5, 6.5, 1.0);
    for angle in ANGLES {
        let solid = revolve(&profile, -5.0, angle);
        assert!(solid.topology().solids()[0].voids.is_empty());
        check_pappus(&solid, angle, outer - inner, 5.0, moment);
    }
}

#[test]
fn a_hollow_circle_keeps_its_bore() {
    let (rho, t, major) = (1.25, 0.5, 4.0);
    let profile = Profile::Circle(CircleProfile {
        radius: rho,
        thickness: Some(t),
    });
    let inner = rho - t;
    for angle in [PI / 2.0, 1.5 * PI] {
        let solid = revolve(&profile, -major, angle);
        check_pappus(
            &solid,
            angle,
            PI * (rho * rho - inner * inner),
            major,
            TAU * (rho + inner) * major,
        );
    }
}

#[test]
fn a_rectangle_on_the_axis_sweeps_a_wedge() {
    // r in [0, 2], h = 3: the inner side lies ON the axis. It sweeps no
    // wall, and the two caps share it as one edge.
    let (b, h) = (2.0, 3.0);
    for angle in ANGLES {
        let solid = revolve(&rect(b, h), -b / 2.0, angle);
        // r ds: outer wall b*h, two sectors b^2/2 each, axis side nothing.
        check_pappus(&solid, angle, b * h, b / 2.0, b * h + b * b);
        assert_eq!(solid.topology().faces().len(), 5, "no wall on the axis");
        let cylinders = solid
            .surfaces()
            .iter()
            .filter(|s| matches!(s, Surface::Cylinder(_)))
            .count();
        assert_eq!(cylinders, 1, "only the outer side sweeps a cylinder");
    }
}

#[test]
fn a_triangle_with_its_apex_on_the_axis_sweeps_a_cone_to_the_apex() {
    // Apex (0, 0) on the axis, base along h = 0 out to r = b, outer side
    // r = b up to h. The hypotenuse sweeps a cone closing at the apex.
    let (b, h) = (2.0, 3.0);
    let profile = contour(vec![
        line(Point2::new(0.0, 0.0), Point2::new(b, 0.0)),
        line(Point2::new(b, 0.0), Point2::new(b, h)),
        line(Point2::new(b, h), Point2::new(0.0, 0.0)),
    ]);
    let slant = (b * b + h * h).sqrt();
    // r ds: base b^2/2, outer side b h, hypotenuse b/2 * slant.
    let moment = b * b / 2.0 + b * h + b * slant / 2.0;
    for angle in ANGLES {
        let solid = revolve(&profile, 0.0, angle);
        check_pappus(&solid, angle, b * h / 2.0, 2.0 * b / 3.0, moment);
        let cones = solid
            .surfaces()
            .iter()
            .filter(|s| matches!(s, Surface::Cone(_)))
            .count();
        assert_eq!(cones, 1);
    }
}

#[test]
fn a_triangle_with_an_edge_on_the_axis_sweeps_a_solid_cone_sector() {
    // (0,0) -> (b,0) -> (0,h): the side x = 0 lies on the axis, so this is
    // a sector of a solid cone with its apex at (0, h).
    let (b, h) = (2.0, 3.0);
    let profile = contour(vec![
        line(Point2::new(0.0, 0.0), Point2::new(b, 0.0)),
        line(Point2::new(b, 0.0), Point2::new(0.0, h)),
        line(Point2::new(0.0, h), Point2::new(0.0, 0.0)),
    ]);
    let slant = (b * b + h * h).sqrt();
    for angle in ANGLES {
        let solid = revolve(&profile, 0.0, angle);
        check_pappus(
            &solid,
            angle,
            b * h / 2.0,
            b / 3.0,
            b * b / 2.0 + b * slant / 2.0,
        );
        // Base sector, cone, two caps: the axis side sweeps nothing.
        assert_eq!(solid.topology().faces().len(), 4);
    }
}

#[test]
fn the_far_side_of_the_axis_and_both_senses_sweep_the_same_amount() {
    // The same rectangle on the other side of the axis, swept the other
    // way, or about the reversed axis: all are the same volume, measured
    // positive, so none of them is built inside out.
    let angle = PI / 2.0;
    let want = angle * 5.0 * 6.0;
    for (axis_x, direction, angle) in [
        (5.0, Vec3::Y, angle),
        (-5.0, Vec3::Y, -angle),
        (-5.0, Vec3::NEG_Y, angle),
        (5.0, Vec3::NEG_Y, -angle),
    ] {
        let solid = revolve_about(&rect(2.0, 3.0), axis_x, direction, angle);
        let props = measure(&solid);
        close("volume", props.signed_volume, want);
        // The solid stays on the profile's own side of the axis: a sweep
        // started from the wrong side has the same volume, mirrored.
        let offset = props.centroid.x - axis_x;
        assert!(
            offset * (0.0 - axis_x) > 0.0,
            "axis at {axis_x}: centroid {:?} left the profile's side",
            props.centroid
        );
    }
    // Mirroring onto the far side flips every arc's turning direction too.
    let circle = Profile::Circle(CircleProfile {
        radius: 1.25,
        thickness: None,
    });
    for axis_x in [-4.0, 4.0] {
        let solid = revolve(&circle, axis_x, angle);
        check_pappus(&solid, angle, PI * 1.25 * 1.25, 4.0, TAU * 1.25 * 4.0);
    }
}

#[test]
fn a_vertex_within_tolerance_of_the_axis_is_on_it() {
    // The inner side sits 1e-8 off the axis, well inside the micrometre
    // tolerance: it must sweep nothing rather than a hair-thin cylinder.
    let (b, h) = (2.0, 3.0);
    let solid = revolve(&rect(b, h), -b / 2.0 - 1e-8, PI / 2.0);
    assert_eq!(solid.topology().faces().len(), 5, "no wall on the axis");
    let volume = measure(&solid).signed_volume;
    let want = PI / 2.0 * (b / 2.0) * b * h;
    assert!((volume - want).abs() < 1e-6, "{volume} vs {want}");
}

#[test]
fn the_sweep_turns_by_the_right_hand_rule_like_the_mesh_path() {
    use axiolid_construct::profile::profile_rings;
    use axiolid_construct::revolve::revolve as revolve_mesh;

    // A quarter turn about +y carries +x towards -z. The mesh path uses
    // Rodrigues' formula; both must put the solid in the same place.
    let chord = 1e-4;
    let tolerance = Tolerance::new(chord, 1e-9).expect("tolerance");
    let profile = rect(2.0, 3.0);
    let rings = profile_rings(&profile, chord, tolerance).expect("rings");
    for (direction, angle) in [
        (Vec3::Y, PI / 2.0),
        (Vec3::Y, -PI / 3.0),
        (Vec3::NEG_Y, 0.75 * PI),
        (Vec3::Y, -PI / 2.0),
        (Vec3::NEG_Y, -0.75 * PI),
    ] {
        let exact = measure(&revolve_about(&profile, -5.0, direction, angle));
        let mesh = revolve_mesh(
            &rings,
            Point3::new(-5.0, 0.0, 0.0),
            direction,
            angle,
            tolerance,
        )
        .expect("mesh revolves");
        let mesh = axiolid_measure::volume_properties(&mesh, tolerance).expect("closed mesh");
        // Signed: the mesh path must orient its faces outward exactly like
        // the exact path does, for either sign of angle and either axis
        // direction (#221). A mismatched sign here means the mesh is
        // inside out even though its magnitude happens to agree.
        assert!(
            (exact.signed_volume - mesh.signed_volume).abs() < 1e-3 * exact.signed_volume,
            "volume: exact {} mesh {}",
            exact.signed_volume,
            mesh.signed_volume
        );
        assert!(
            (exact.centroid - mesh.centroid).length() < 1e-3,
            "centroid: exact {:?} mesh {:?}",
            exact.centroid,
            mesh.centroid
        );
    }
    // And the direction itself, from geometry alone: the centroid of a
    // quarter turn about +y from +x lies at negative z.
    let quarter = measure(&revolve(&profile, -5.0, PI / 2.0));
    assert!(quarter.centroid.z < -1.0, "{:?}", quarter.centroid);
    assert!(quarter.centroid.x > -5.0 + 1.0, "{:?}", quarter.centroid);
}

#[test]
fn the_mesh_path_orients_outward_for_a_full_turn_either_sense() {
    use axiolid_construct::profile::profile_rings;
    use axiolid_construct::revolve::revolve as revolve_mesh;

    // The annular-tube fixture shared with the exact-path tests: axis 5 to
    // the rectangle's left, so r in [4, 6]. Every wall triangle's radial
    // sense is unambiguous regardless of sign or axis direction (#221).
    let (axis_x, b, h) = (-5.0, 2.0, 3.0);
    let (inner, outer) = (4.0, 6.0);
    let chord = 1e-3;
    let tolerance = Tolerance::new(chord, 1e-9).expect("tolerance");
    let profile = rect(b, h);
    let rings = profile_rings(&profile, chord, tolerance).expect("rings");
    let axis_origin = Point3::new(axis_x, 0.0, 0.0);

    for (direction, angle) in [
        (Vec3::Y, TAU),
        (Vec3::Y, -TAU),
        (Vec3::NEG_Y, TAU),
        (Vec3::NEG_Y, -TAU),
    ] {
        let mesh =
            revolve_mesh(&rings, axis_origin, direction, angle, tolerance).expect("mesh revolves");
        let props = axiolid_measure::volume_properties(&mesh, tolerance).expect("closed mesh");
        let want = TAU * 5.0 * (outer - inner) * h;
        assert!(
            (props.signed_volume - want).abs() < 1e-3 * want,
            "direction {direction:?} angle {angle}: volume {} vs {want}",
            props.signed_volume
        );
        assert_wall_normals_radial(&mesh, axis_origin, direction, inner, outer, chord);
    }
}

#[test]
fn the_mesh_path_orients_outward_for_a_partial_turn_either_sense() {
    use axiolid_construct::profile::profile_rings;
    use axiolid_construct::revolve::revolve as revolve_mesh;

    let (axis_x, b, h) = (-5.0, 2.0, 3.0);
    let (inner, outer) = (4.0, 6.0);
    let chord = 1e-3;
    let tolerance = Tolerance::new(chord, 1e-9).expect("tolerance");
    let profile = rect(b, h);
    let rings = profile_rings(&profile, chord, tolerance).expect("rings");
    let axis_origin = Point3::new(axis_x, 0.0, 0.0);

    for (direction, angle) in [
        (Vec3::Y, PI / 2.0),
        (Vec3::Y, -PI / 2.0),
        (Vec3::NEG_Y, PI / 2.0),
        (Vec3::NEG_Y, -PI / 2.0),
    ] {
        let mesh =
            revolve_mesh(&rings, axis_origin, direction, angle, tolerance).expect("mesh revolves");
        let props = axiolid_measure::volume_properties(&mesh, tolerance).expect("closed mesh");
        let want = angle.abs() * 5.0 * (outer - inner) * h;
        assert!(
            props.signed_volume > 0.0,
            "direction {direction:?} angle {angle}: signed volume {} must be positive",
            props.signed_volume
        );
        assert!(
            (props.signed_volume - want).abs() < 1e-3 * want,
            "direction {direction:?} angle {angle}: volume {} vs {want}",
            props.signed_volume
        );
        assert_wall_normals_radial(&mesh, axis_origin, direction, inner, outer, chord);
    }
}

/// Every wall triangle (all three vertices at the same radius) must face
/// away from the axis at the outer radius and towards it at the inner
/// radius. Cap and mixed triangles are skipped: their outward sense is not
/// a pure radial statement.
fn assert_wall_normals_radial(
    mesh: &axiolid_mesh::TriMesh,
    axis_origin: Point3,
    axis_direction: Vec3,
    inner: f64,
    outer: f64,
    tol: f64,
) {
    let dir = axis_direction / axis_direction.length();
    let radial = |p: Point3| -> (f64, Vec3) {
        let v = p - axis_origin;
        let r_vec = v - dir * dir.dot(v);
        (r_vec.length(), r_vec)
    };
    let mut outer_checked = 0;
    let mut inner_checked = 0;
    for i in 0..mesh.triangle_count() {
        let tri = [
            mesh.indices[i * 3] as usize,
            mesh.indices[i * 3 + 1] as usize,
            mesh.indices[i * 3 + 2] as usize,
        ];
        let [a, b, c] = tri.map(|j| mesh.positions[j]);
        let (ra, rvec_a) = radial(a);
        let (rb, _) = radial(b);
        let (rc, _) = radial(c);
        let normal = (b - a).cross(c - a);
        if normal.length() < 1e-12 {
            continue;
        }
        let is_outer =
            (ra - outer).abs() < tol && (rb - outer).abs() < tol && (rc - outer).abs() < tol;
        let is_inner =
            (ra - inner).abs() < tol && (rb - inner).abs() < tol && (rc - inner).abs() < tol;
        if is_outer {
            assert!(
                normal.dot(rvec_a) > 0.0,
                "triangle {i} on the outer wall has an inward normal {normal:?}"
            );
            outer_checked += 1;
        } else if is_inner {
            assert!(
                normal.dot(rvec_a) < 0.0,
                "triangle {i} on the inner wall has an outward normal {normal:?}"
            );
            inner_checked += 1;
        }
    }
    assert!(outer_checked > 0, "no outer-wall triangles found to check");
    assert!(inner_checked > 0, "no inner-wall triangles found to check");
}

#[test]
fn a_full_turn_still_takes_the_closed_path() {
    // Both senses of a full turn are the annular tube: no end caps, no
    // planar walls at the start and end angles.
    for angle in [TAU, -TAU] {
        let solid = revolve(&rect(2.0, 3.0), -5.0, angle);
        assert_eq!(solid.topology().faces().len(), 4);
        close("volume", measure(&solid).signed_volume, TAU * 5.0 * 6.0);
    }
}

#[test]
fn angles_that_state_no_partial_turn_are_refused() {
    for angle in [0.0, -0.0, f64::NAN, f64::INFINITY] {
        let error = try_revolve(&rect(2.0, 3.0), -5.0, Vec3::Y, angle)
            .expect_err("no solid without a finite non-zero angle");
        assert!(matches!(error, GeomError::InvalidInput(_)), "{error:?}");
    }
    for angle in [1.5 * TAU, -3.0 * PI] {
        let error = try_revolve(&rect(2.0, 3.0), -5.0, Vec3::Y, angle)
            .expect_err("beyond a full turn sweeps through itself");
        assert!(
            matches!(
                error,
                GeomError::UnsupportedInput {
                    input: "exact revolution beyond a full turn",
                    ..
                }
            ),
            "{error:?}"
        );
    }
    // A turn so small the outermost point moves less than the tolerance
    // would collapse the start and end caps onto each other.
    let tiny = Tolerance::METRE.linear() / 100.0;
    let error = try_revolve(&rect(2.0, 3.0), -5.0, Vec3::Y, tiny)
        .expect_err("a sweep below the tolerance is no solid");
    assert!(matches!(error, GeomError::Degenerate(_)), "{error:?}");
}

#[test]
fn a_section_crossing_the_axis_is_refused() {
    let error = try_revolve(&rect(4.0, 2.0), 0.0, Vec3::Y, PI / 2.0)
        .expect_err("the section straddles the axis");
    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                input: "exact revolution of a section crossing the axis",
                ..
            }
        ),
        "{error:?}"
    );
}

#[test]
fn an_arc_reaching_the_axis_is_refused() {
    // A circle tangent to the axis sweeps a horn torus.
    let error = try_revolve(
        &Profile::Circle(CircleProfile {
            radius: 1.0,
            thickness: None,
        }),
        -1.0,
        Vec3::Y,
        PI / 2.0,
    )
    .expect_err("the tube reaches the axis");
    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                input: "revolved arc whose tube reaches the axis (spindle torus)",
                ..
            }
        ),
        "{error:?}"
    );
}

#[test]
fn an_ellipse_partial_turn_is_refused() {
    let error = try_revolve(
        &Profile::Ellipse(EllipseProfile {
            semi_axis_x: 1.0,
            semi_axis_y: 0.5,
        }),
        -4.0,
        Vec3::Y,
        PI / 2.0,
    )
    .expect_err("an ellipse sweeps no elementary surface");
    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                input: "ellipse exact revolution",
                ..
            }
        ),
        "{error:?}"
    );
}
