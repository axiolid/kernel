//! A bounded half-space against its closed-form volume.
//!
//! The construction is a prism: the boundary profile swept along the plane
//! normal by a depth derived from the boundary's own extent. Its volume is
//! therefore area x depth, a constant this crate does not get to choose.

use axiolid_construct::half_space::{bounded_half_space, bounded_half_space_in_frame};
use axiolid_construct::profile::Rings;
use axiolid_core::{Plane3, PlaneFrame, Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_primitive::ClipMargin;

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("tolerance")
}

fn volume(mesh: &TriMesh) -> Scalar {
    volume_properties(mesh, tol())
        .expect("a bounded half-space must be closed and two-manifold")
        .signed_volume
}

/// A square boundary of half-width `h`, wound counter-clockwise.
fn square(h: Scalar) -> Rings {
    Rings {
        outer: vec![
            Point2::new(-h, -h),
            Point2::new(h, -h),
            Point2::new(h, h),
            Point2::new(-h, h),
        ],
        holes: Vec::new(),
    }
}

fn plane_z() -> Plane3 {
    Plane3 {
        origin: Point3::ZERO,
        normal: Vec3::Z,
    }
}

#[test]
fn a_bounded_half_space_is_area_times_depth() {
    // Square of half-width 5 => area 100, extent 5, margin 2 => depth 10.
    let margin = ClipMargin::new(2.0).expect("margin");
    let mesh = bounded_half_space(&square(5.0), plane_z(), true, margin, tol()).expect("clip");
    let want = 100.0 * 10.0;
    let got = volume(&mesh);
    assert!(
        (got - want).abs() / want < 1e-9,
        "bounded half-space {got} vs {want}"
    );
}

#[test]
fn the_margin_scales_the_depth_proportionally() {
    // Doubling the margin must double the volume and nothing else: this is
    // what proves the depth is derived from the margin rather than from a
    // constant that happens to fit the first test.
    let a = bounded_half_space(
        &square(5.0),
        plane_z(),
        true,
        ClipMargin::new(2.0).expect("margin"),
        tol(),
    )
    .expect("clip");
    let b = bounded_half_space(
        &square(5.0),
        plane_z(),
        true,
        ClipMargin::new(4.0).expect("margin"),
        tol(),
    )
    .expect("clip");
    let ratio = volume(&b) / volume(&a);
    assert!((ratio - 2.0).abs() < 1e-9, "margin ratio {ratio} vs 2");
}

#[test]
fn the_construction_is_unit_independent() {
    // The same boundary scaled by 1000 must give exactly 1000^3 the volume.
    // Sizing the slab from the boundary's own extent is what guarantees
    // this; an absolute constant would break it.
    let margin = ClipMargin::new(2.0).expect("margin");
    let small = bounded_half_space(&square(5.0), plane_z(), true, margin, tol()).expect("clip");
    let large = bounded_half_space(&square(5000.0), plane_z(), true, margin, tol()).expect("clip");
    let ratio = volume(&large) / volume(&small);
    assert!(
        (ratio - 1e9).abs() / 1e9 < 1e-9,
        "unit scaling {ratio} vs 1e9"
    );
}

#[test]
fn agreement_selects_the_opposite_side() {
    // Both sides must be genuine solids of equal volume, and they must sit
    // on opposite sides of the plane. Equal volume alone would also hold if
    // agreement were ignored, so the z extent is what makes this a real test.
    let margin = ClipMargin::new(2.0).expect("margin");
    let up = bounded_half_space(&square(5.0), plane_z(), true, margin, tol()).expect("clip");
    let down = bounded_half_space(&square(5.0), plane_z(), false, margin, tol()).expect("clip");
    assert!((volume(&up) - volume(&down)).abs() < 1e-9, "equal volume");

    let max_z = |m: &TriMesh| m.positions.iter().fold(Scalar::MIN, |a, p| a.max(p.z));
    let min_z = |m: &TriMesh| m.positions.iter().fold(Scalar::MAX, |a, p| a.min(p.z));
    assert!(
        max_z(&up) > 0.0 && min_z(&up) >= -1e-12,
        "normal side is +z"
    );
    assert!(
        min_z(&down) < 0.0 && max_z(&down) <= 1e-12,
        "opposite side is -z"
    );
}

#[test]
fn a_degenerate_boundary_is_refused() {
    let margin = ClipMargin::new(2.0).expect("margin");
    // Zero extent cannot size a slab. The ring has three distinct points
    // so it passes the vertex-count check and reaches the extent guard,
    // which is the condition actually under test here.
    let flat = Rings {
        outer: vec![
            Point2::new(0.0, 0.0),
            Point2::new(0.0, 0.0),
            Point2::new(0.0, 0.0),
        ],
        holes: Vec::new(),
    };
    assert!(bounded_half_space(&flat, plane_z(), true, margin, tol()).is_err());
    assert!(bounded_half_space(&square(0.0), plane_z(), true, margin, tol()).is_err());
    // A zero normal has no side to select.
    let bad_plane = Plane3 {
        origin: Point3::ZERO,
        normal: Vec3::ZERO,
    };
    assert!(bounded_half_space(&square(5.0), bad_plane, true, margin, tol()).is_err());
    // Two points do not bound a region.
    let sliver = Rings {
        outer: vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)],
        holes: Vec::new(),
    };
    assert!(bounded_half_space(&sliver, plane_z(), true, margin, tol()).is_err());
}

/// An asymmetric L-shaped boundary, wound counter-clockwise.
///
/// A symmetric profile cannot detect a mirrored footprint: reflecting a
/// square across its own axis is the identity. The L is deliberately
/// chiral so a mirror is observable.
fn l_shape() -> Rings {
    Rings {
        outer: vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(4.0, 1.0),
            Point2::new(1.0, 1.0),
            Point2::new(1.0, 3.0),
            Point2::new(0.0, 3.0),
        ],
        holes: Vec::new(),
    }
}

/// Unique in-plane footprint, rounded so exact vertex order does not matter.
fn footprint_xy(mesh: &TriMesh) -> Vec<(i64, i64)> {
    let mut seen: Vec<(i64, i64)> = mesh
        .positions
        .iter()
        .map(|p| ((p.x * 1e6).round() as i64, (p.y * 1e6).round() as i64))
        .collect();
    seen.sort_unstable();
    seen.dedup();
    seen
}

#[test]
fn agreement_moves_the_solid_without_reshaping_the_boundary() {
    // `agreement` selects which side of the plane is material. It is not a
    // transform of the boundary, so the in-plane footprint must be identical
    // for both values; only the swept z extent may differ.
    let margin = ClipMargin::new(2.0).expect("margin");
    let up = bounded_half_space(&l_shape(), plane_z(), true, margin, tol()).expect("clip");
    let down = bounded_half_space(&l_shape(), plane_z(), false, margin, tol()).expect("clip");

    assert_eq!(
        footprint_xy(&up),
        footprint_xy(&down),
        "agreement must not mirror or otherwise reshape the boundary footprint"
    );
}

#[test]
fn both_agreement_values_stay_outward_wound() {
    // The mirrored-basis workaround kept winding correct by reflecting the
    // profile. Any replacement must keep the solid outward-facing on its
    // own, so assert positive volume for both sides.
    let margin = ClipMargin::new(2.0).expect("margin");
    for agreement in [true, false] {
        let mesh =
            bounded_half_space(&l_shape(), plane_z(), agreement, margin, tol()).expect("clip");
        assert!(
            volume(&mesh) > 0.0,
            "agreement={agreement} produced an inside-out solid"
        );
    }
}

/// An authored in-plane frame must orient the boundary footprint.
///
/// The clip plane here is the z = 0 ground plane, so the internal heuristic
/// picks Vec3::X as its reference and lands the boundary axis-aligned. A
/// consumer whose source format authors its own in-plane rotation (IFC
/// IfcPolygonalBoundedHalfSpace.Position) must be able to override that.
///
/// Volume cannot detect this: rotating a profile about the sweep axis
/// preserves area and depth. The footprint geometry is what changes, so the
/// test pins a vertex position rather than a scalar measure.
#[test]
fn an_authored_in_plane_frame_orients_the_boundary() {
    let ground = PlaneFrame::new(Point3::ZERO, Vec3::X, Vec3::Y, tol())
        .expect("the ground plane is a valid frame");
    let angle = std::f64::consts::FRAC_PI_4;
    let rotated = PlaneFrame::new(
        Point3::ZERO,
        Vec3::new(angle.cos(), angle.sin(), 0.0),
        Vec3::new(-angle.sin(), angle.cos(), 0.0),
        tol(),
    )
    .expect("a rotation about the plane normal is a valid frame");

    let plane = Plane3 {
        origin: Point3::ZERO,
        normal: Vec3::Z,
    };
    let boundary = square(1.0);
    let axis_aligned = bounded_half_space_in_frame(
        &boundary,
        plane,
        ground,
        true,
        ClipMargin::new(1.0).unwrap(),
        tol(),
    )
    .expect("an axis-aligned boundary is valid");
    let turned = bounded_half_space_in_frame(
        &boundary,
        plane,
        rotated,
        true,
        ClipMargin::new(1.0).unwrap(),
        tol(),
    )
    .expect("a rotated boundary is valid");

    // The same profile in a frame turned 45 degrees must not produce the same
    // footprint. If the frame were ignored, these meshes would be identical.
    let differs = axis_aligned
        .positions
        .iter()
        .zip(turned.positions.iter())
        .any(|(a, b)| (*a - *b).length() > 1e-9);
    assert!(
        differs,
        "the authored in-plane frame was ignored: both footprints are identical"
    );

    // Stronger than "differs": the footprint must land exactly where the
    // authored frame puts it. Corner (1,1) of the profile maps to
    // origin + x*1 + y*1, which for a 45-degree frame is (0, sqrt(2)).
    let want = rotated.lift(Point2::new(1.0, 1.0));
    let found = turned.positions.iter().any(|p| (*p - want).length() < 1e-9);
    assert!(
        found,
        "no vertex at the authored corner {want:?}; frame was not applied as written"
    );
}

/// The authored frame's origin anchors the footprint; only its offset along
/// the clip normal is dropped (#164).
///
/// The frame sits at (3, -2, 5) while the clip plane passes through the
/// origin with a +Z normal. The in-plane part (3, -2) must move every
/// footprint vertex; the 5 along the normal must not lift the slab off the
/// clip plane. Before the fix the boundary was anchored at the plane origin,
/// so the footprint stayed centred on (0, 0) with no error.
#[test]
fn an_authored_frame_origin_places_the_boundary_in_the_plane() {
    let offset = PlaneFrame::new(Point3::new(3.0, -2.0, 5.0), Vec3::X, Vec3::Y, tol())
        .expect("a translated ground frame is valid");
    let mesh = bounded_half_space_in_frame(
        &square(1.0),
        plane_z(),
        offset,
        true,
        ClipMargin::new(1.0).unwrap(),
        tol(),
    )
    .expect("a translated boundary is valid");

    let near = |want: Point3| mesh.positions.iter().any(|p| (*p - want).length() < 1e-9);
    // Profile corner (1, 1) at the in-plane anchor (3, -2), on the plane z = 0.
    assert!(
        near(Point3::new(4.0, -1.0, 0.0)),
        "the footprint must follow the frame's in-plane origin: {:?}",
        mesh.positions
    );
    let (min_z, max_z) = mesh
        .positions
        .iter()
        .fold((Scalar::INFINITY, Scalar::NEG_INFINITY), |(lo, hi), p| {
            (lo.min(p.z), hi.max(p.z))
        });
    assert_eq!(
        min_z, 0.0,
        "the sweep must start on the clip plane, not at z = 5"
    );
    assert!(max_z > 0.0, "agreement = true keeps the normal side");
    // Translation is rigid, so the solid is still area x depth.
    let unmoved = bounded_half_space(
        &square(1.0),
        plane_z(),
        true,
        ClipMargin::new(1.0).unwrap(),
        tol(),
    )
    .expect("clip");
    assert!((volume(&mesh) - volume(&unmoved)).abs() < 1e-9);
}
