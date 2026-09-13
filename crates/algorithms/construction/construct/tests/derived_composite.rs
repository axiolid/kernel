//! Derived and composite profile extrusion (ADR 0054).
//!
//! Geometry is checked against independently-computed positions and areas,
//! not against another call of the same lowering code.

use axiolid_brep_audit::geometric_audit;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{Point2, Tolerance, Transform2, Vec2, Vec3};
use axiolid_measure::exact_properties;
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};
use axiolid_surface::Surface;

fn rect(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

fn derived(basis: Profile, transform: Transform2) -> Profile {
    Profile::Derived {
        basis: Box::new(basis),
        transform,
    }
}

fn base_points(solid: &axiolid_brep::ExactBRep) -> Vec<Point2> {
    let mut points: Vec<Point2> = solid
        .topology()
        .vertices()
        .iter()
        .filter(|vertex| vertex.position.z.abs() < 1e-12)
        .map(|vertex| Point2::new(vertex.position.x, vertex.position.y))
        .collect();
    points.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));
    points
}

#[test]
fn a_translated_rectangle_lands_where_the_transform_puts_it() {
    let shift = Vec2::new(5.0, -3.0);
    let profile = derived(rect(2.0, 4.0), Transform2::from_translation(shift));
    let solid = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
        .expect("a translated rectangle extrudes");

    // Closed form: the four corners of a 2x4 rectangle, shifted.
    let mut expected = [
        Point2::new(-1.0, -2.0) + shift,
        Point2::new(1.0, -2.0) + shift,
        Point2::new(1.0, 2.0) + shift,
        Point2::new(-1.0, 2.0) + shift,
    ];
    expected.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));

    let found = base_points(&solid);
    assert_eq!(found.len(), 4, "a rectangle has four base corners");
    for (got, want) in found.iter().zip(expected.iter()) {
        assert!(
            (*got - *want).length() < 1e-12,
            "corner {got:?} should be {want:?}"
        );
    }
    assert!(geometric_audit(&solid, Tolerance::METRE).is_consistent());
}

#[test]
fn a_rotated_rectangle_keeps_its_area_and_moves_its_corners() {
    let angle = 0.7_f64;
    let profile = derived(rect(2.0, 4.0), Transform2::from_angle(angle));
    let solid = extrude_profile_exact(&profile, Vec3::Z, 3.0, Tolerance::METRE)
        .expect("a rotated rectangle extrudes");

    // Rotation preserves area, so the volume is unchanged: 2*4*3.
    let properties = exact_properties(&solid, Tolerance::METRE).expect("volume");
    assert!(
        (properties.signed_volume - 24.0).abs() < 1e-9,
        "rotation must preserve volume, got {}",
        properties.signed_volume
    );

    // ...but the corners genuinely moved.
    let (sin, cos) = angle.sin_cos();
    // Rotate the corner (-1, -2) by `angle`, written out as the closed form.
    let (cx, cy) = (-1.0_f64, -2.0_f64);
    let want = Point2::new(cx * cos - cy * sin, cx * sin + cy * cos);
    let found = base_points(&solid);
    assert!(
        found.iter().any(|p| (*p - want).length() < 1e-12),
        "expected a corner at {want:?}, got {found:?}"
    );
    assert!(geometric_audit(&solid, Tolerance::METRE).is_consistent());
}

#[test]
fn a_uniformly_scaled_circle_stays_an_exact_circle() {
    let profile = derived(
        Profile::Circle(CircleProfile {
            radius: 1.5,
            thickness: None,
        }),
        Transform2::from_scale(Vec2::splat(2.0)),
    );
    let solid = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
        .expect("a scaled circle extrudes");

    let radii: Vec<f64> = solid
        .surfaces()
        .iter()
        .filter_map(|surface| match surface {
            Surface::Cylinder(cylinder) => Some(cylinder.radius),
            _ => None,
        })
        .collect();
    assert_eq!(radii.len(), 1, "one cylindrical wall");
    assert!(
        (radii[0] - 3.0).abs() < 1e-12,
        "radius 1.5 scaled by 2 is 3.0, got {}",
        radii[0]
    );
}

#[test]
fn a_sheared_circle_is_refused_because_it_is_an_ellipse() {
    // A shear has determinant 1, so a determinant check would let it pass
    // while the circle is genuinely distorted into an ellipse.
    // Columns (1,0) and (0.5,1): a unit-determinant shear.
    let mut shear = Transform2::IDENTITY;
    shear.matrix2.y_axis = Vec2::new(0.5, 1.0);
    assert!(
        (shear.matrix2.determinant() - 1.0).abs() < 1e-12,
        "the shear must have unit determinant for this test to mean anything"
    );
    let profile = derived(
        Profile::Circle(CircleProfile {
            radius: 1.0,
            thickness: None,
        }),
        shear,
    );
    let error = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
        .expect_err("a sheared circle is an ellipse");
    assert!(
        format!("{error:?}").contains("non-conformal"),
        "got {error:?}"
    );
}

fn placed(x: f64, y: f64, at: Vec2) -> Profile {
    derived(rect(x, y), Transform2::from_translation(at))
}

#[test]
fn overlapping_members_are_unioned_not_double_counted() {
    // Two 2x2 squares overlapping in a 1x2 strip. Concatenating them as
    // rings would count the shared strip twice and emit crossing walls;
    // the union area is 4 + 4 - 2 = 6.
    let profile = Profile::Composite(vec![
        placed(2.0, 2.0, Vec2::new(0.0, 0.0)),
        placed(2.0, 2.0, Vec2::new(1.0, 0.0)),
    ]);
    let depth = 3.0;
    let solid = extrude_profile_exact(&profile, Vec3::Z, depth, Tolerance::METRE)
        .expect("overlapping members union");

    let properties = exact_properties(&solid, Tolerance::METRE).expect("volume");
    assert!(
        (properties.signed_volume - 6.0 * depth).abs() < 1e-9,
        "union area 6 times depth {depth}, got {}",
        properties.signed_volume
    );
    assert!(geometric_audit(&solid, Tolerance::METRE).is_consistent());
}

#[test]
fn members_forming_a_frame_produce_a_genuine_hole() {
    // Four bars arranged as a picture frame: the union has one outer
    // boundary and one hole, and the hole must survive into the solid.
    let profile = Profile::Composite(vec![
        placed(6.0, 2.0, Vec2::new(0.0, 2.0)),
        placed(6.0, 2.0, Vec2::new(0.0, -2.0)),
        placed(2.0, 6.0, Vec2::new(-2.0, 0.0)),
        placed(2.0, 6.0, Vec2::new(2.0, 0.0)),
    ]);
    let depth = 2.0;
    let solid = extrude_profile_exact(&profile, Vec3::Z, depth, Tolerance::METRE)
        .expect("a picture frame extrudes");

    // Outer 6x6 minus a 2x2 hole: 36 - 4 = 32.
    let properties = exact_properties(&solid, Tolerance::METRE).expect("volume");
    assert!(
        (properties.signed_volume - 32.0 * depth).abs() < 1e-9,
        "a dropped hole would read 36*depth; got {}",
        properties.signed_volume
    );
    assert!(geometric_audit(&solid, Tolerance::METRE).is_consistent());
}

#[test]
fn disjoint_members_are_refused_not_silently_reduced_to_one() {
    // Two separated squares are two bodies. `Solid` holds one outer shell
    // plus voids, so returning either square alone would silently discard
    // the other.
    let profile = Profile::Composite(vec![
        placed(2.0, 2.0, Vec2::new(0.0, 0.0)),
        placed(2.0, 2.0, Vec2::new(10.0, 0.0)),
    ]);
    let error = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
        .expect_err("disjoint members are two solids");
    assert!(
        format!("{error:?}").contains("disjoint"),
        "the refusal must name the disjointness, got {error:?}"
    );
}

#[test]
fn an_empty_composite_is_refused() {
    let error = extrude_profile_exact(
        &Profile::Composite(Vec::new()),
        Vec3::Z,
        1.0,
        Tolerance::METRE,
    )
    .expect_err("an empty composite has no section");
    assert!(format!("{error:?}").contains("at least one member"));
}

#[test]
fn a_doubly_derived_profile_composes_its_transforms() {
    // Nesting must compose rather than apply only the outer transform.
    let inner = derived(
        rect(2.0, 2.0),
        Transform2::from_translation(Vec2::new(3.0, 0.0)),
    );
    let outer = derived(inner, Transform2::from_translation(Vec2::new(0.0, 4.0)));
    let solid = extrude_profile_exact(&outer, Vec3::Z, 1.0, Tolerance::METRE)
        .expect("a nested derived profile extrudes");

    // Both shifts must land: centre at (3, 4), so a corner at (2, 3).
    let want = Point2::new(2.0, 3.0);
    let found = base_points(&solid);
    assert!(
        found.iter().any(|p| (*p - want).length() < 1e-12),
        "expected a corner at {want:?} from both shifts, got {found:?}"
    );
}
