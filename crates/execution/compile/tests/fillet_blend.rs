//! The fillet builds a real cylindrical blend (#79).
//!
//! Lives in `axiolid-mesh-compile` because it checks the tessellated
//! result, which is what a mesh consumer of the fillet receives.

use axiolid_construct::feature::{fillet_extruded_profile, EdgeSelector, FeatureSize};
use axiolid_core::{Point2, Tolerance, Vec3};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_surface::Surface;
use axiolid_topology::audit_brep;

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("tolerance")
}

fn rectangle(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// The blend surface is genuinely a cylinder, not a fan of planes.
///
/// This is the claim the v0.6 contract refused to make: a many-segment
/// chamfer is indistinguishable to a caller reading the mesh, so the
/// distinction has to be checked on the B-rep itself.
#[test]
fn the_blend_face_is_a_cylinder() {
    let filleted = fillet_extruded_profile(
        &rectangle(4.0, 3.0),
        Vec3::Z,
        2.0,
        EdgeSelector::NearestCorner(Point2::new(2.0, 1.5)),
        FeatureSize::ConstantRadius(0.5),
        tol(),
    )
    .expect("filletable");

    let cylinders = filleted
        .surfaces()
        .iter()
        .filter(|s| matches!(s, Surface::Cylinder(_)))
        .count();
    assert_eq!(
        cylinders, 1,
        "the fillet must produce exactly one cylindrical blend face"
    );
}

/// The filleted shell is closed and manifold with the blend stitched in.
///
/// `finish_closed` audits this during construction, so reaching here at all
/// means the cylinder shares its edges and vertices with both neighbouring
/// walls -- which is what tangency means topologically. Asserting it again
/// states the invariant rather than trusting the constructor silently.
#[test]
fn the_filleted_shell_is_closed_and_manifold() {
    let filleted = fillet_extruded_profile(
        &rectangle(4.0, 3.0),
        Vec3::Z,
        2.0,
        EdgeSelector::NearestCorner(Point2::new(2.0, 1.5)),
        FeatureSize::ConstantRadius(0.5),
        tol(),
    )
    .expect("filletable");

    let health = audit_brep(filleted.topology());
    assert!(
        health.is_closed_manifold(),
        "filleted shell is not closed manifold: {health:?}"
    );
}

/// Every corner of the rectangle can be filleted, not just one.
#[test]
fn any_corner_can_be_filleted() {
    for target in [
        Point2::new(2.0, 1.5),
        Point2::new(-2.0, 1.5),
        Point2::new(2.0, -1.5),
        Point2::new(-2.0, -1.5),
    ] {
        let filleted = fillet_extruded_profile(
            &rectangle(4.0, 3.0),
            Vec3::Z,
            2.0,
            EdgeSelector::NearestCorner(target),
            FeatureSize::ConstantRadius(0.4),
            tol(),
        )
        .unwrap_or_else(|error| panic!("corner {target:?} must be filletable: {error}"));
        assert!(audit_brep(filleted.topology()).is_closed_manifold());
    }
}

/// The blend is tangent to both neighbouring walls.
///
/// Tangency is meant to follow from construction: the centre sits on the
/// internal bisector, so its perpendicular distance to each adjacent wall is
/// exactly the radius. That is measurable, and measuring it is what stops
/// the claim being decorative -- a centre pushed off the bisector still
/// yields a closed manifold shell with one cylinder in it, so the structural
/// tests alone cannot tell the difference.
#[test]
fn the_blend_is_tangent_to_both_neighbouring_walls() {
    let (x, y, radius) = (4.0, 3.0, 0.5);
    let filleted = fillet_extruded_profile(
        &rectangle(x, y),
        Vec3::Z,
        2.0,
        EdgeSelector::NearestCorner(Point2::new(x / 2.0, y / 2.0)),
        FeatureSize::ConstantRadius(radius),
        tol(),
    )
    .expect("filletable");

    let cylinder = filleted
        .surfaces()
        .iter()
        .find_map(|s| match s {
            Surface::Cylinder(c) => Some(*c),
            _ => None,
        })
        .expect("the fillet emits a cylinder");

    assert!(
        (cylinder.radius - radius).abs() < 1e-12,
        "blend radius: expected {radius}, got {}",
        cylinder.radius
    );

    // The two walls at the +x/+y corner are the planes x = 2 and y = 1.5.
    // Tangency means the axis stands exactly `radius` from each.
    let axis = cylinder.frame.origin;
    let to_x_wall = (x / 2.0) - axis.x;
    let to_y_wall = (y / 2.0) - axis.y;
    assert!(
        (to_x_wall - radius).abs() < 1e-12,
        "blend axis is {to_x_wall} from the x wall, expected {radius}"
    );
    assert!(
        (to_y_wall - radius).abs() < 1e-12,
        "blend axis is {to_y_wall} from the y wall, expected {radius}"
    );
}
