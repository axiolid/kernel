//! A bounded half-space against its closed-form volume.
//!
//! The construction is a prism: the boundary profile swept along the plane
//! normal by a depth derived from the boundary's own extent. Its volume is
//! therefore area x depth, a constant this crate does not get to choose.

use axiolid_construct::half_space::{
    bounded_half_space, bounded_half_space_in_frame, bounded_half_space_prism, BoundaryPrism,
    SlabDepth,
};
use axiolid_construct::profile::Rings;
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Plane3, PlaneFrame, Point2, Point3, Scalar, Tolerance, Transform3, Vec3};
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

// ---------------------------------------------------------------------------
// A prism along the placement's z (#321)
// ---------------------------------------------------------------------------

/// The unit square `[0, 1]^2`, counter-clockwise.
fn unit_square() -> Rings {
    Rings {
        outer: vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ],
        holes: Vec::new(),
    }
}

/// The frame with columns `x`, `y`, `z` at `origin`.
fn placement(x: Vec3, y: Vec3, z: Vec3, origin: Point3) -> Transform3 {
    Transform3::from_cols(x, y, z, origin)
}

fn prism_along_z(
    boundary: &Rings,
    plane: Plane3,
    frame: Transform3,
    agreement: bool,
    reach: Scalar,
) -> GeomResult<TriMesh> {
    let prism = BoundaryPrism::along_placement_z(frame, tol())?;
    bounded_half_space_prism(
        boundary,
        plane,
        prism,
        agreement,
        SlabDepth::Reach(reach),
        tol(),
    )
}

fn has_vertex(mesh: &TriMesh, want: Point3) -> bool {
    mesh.positions.iter().any(|p| (*p - want).length() < 1e-9)
}

/// A prism leaning off the plane normal is an oblique prism: its volume is
/// the area its section perpendicular to the normal has times the reach,
/// and its far cap is the near one moved along the axis.
#[test]
fn an_oblique_prism_is_its_section_times_its_reach() {
    let frame = placement(Vec3::X, Vec3::Y, Vec3::new(0.5, 0.0, 1.0), Point3::ZERO);
    let mesh = prism_along_z(&unit_square(), plane_z(), frame, true, 2.0).expect("clip");
    assert!((volume(&mesh) - 2.0).abs() < 1e-12, "{}", volume(&mesh));
    // Corner (1, 0) on the plane, and moved by (1, 0, 2) to the far cap.
    assert!(has_vertex(&mesh, Point3::new(1.0, 0.0, 0.0)));
    assert!(has_vertex(&mesh, Point3::new(2.0, 0.0, 2.0)));
    // The axis's length and sign do not matter: the prism is infinite both
    // ways, and the slab runs to the kept side whichever way it points.
    for z in [Vec3::new(-1.0, 0.0, -2.0), Vec3::new(2.5, 0.0, 5.0)] {
        let other = prism_along_z(
            &unit_square(),
            plane_z(),
            placement(Vec3::X, Vec3::Y, z, Point3::ZERO),
            true,
            2.0,
        )
        .expect("clip");
        assert_eq!(footprint(&other), footprint(&mesh), "z = {z:?}");
    }
}

/// Every vertex, rounded, sorted: a mesh's shape independent of its order.
fn footprint(mesh: &TriMesh) -> Vec<[i64; 3]> {
    let mut seen: Vec<[i64; 3]> = mesh
        .positions
        .iter()
        .map(|p| [p.x, p.y, p.z].map(|c| (c * 1e9).round() as i64))
        .collect();
    seen.sort_unstable();
    seen.dedup();
    seen
}

/// The #321 shape: a vertical prism under a sloped plane. The slab between
/// the plane and its parallel at the reach is cut from the vertical prism,
/// so its volume is the boundary's area times the reach over the cosine of
/// the slope; the near cap is the boundary lifted vertically onto the
/// plane, not laid into it.
#[test]
fn a_vertical_prism_under_a_sloped_plane_is_cut_vertically() {
    let sloped = Plane3 {
        origin: Point3::ZERO,
        normal: Vec3::new(0.6, 0.0, 0.8),
    };
    let frame = placement(Vec3::X, Vec3::Y, Vec3::Z, Point3::new(0.0, 0.0, 7.0));
    for agreement in [true, false] {
        let mesh = prism_along_z(&unit_square(), sloped, frame, agreement, 1.0).expect("clip");
        assert!((volume(&mesh) - 1.25).abs() < 1e-12, "{}", volume(&mesh));
        // Corner (1, 1) at z = -0.75 on the plane; the frame's own height
        // does not move the prism, which is infinite along z.
        assert!(has_vertex(&mesh, Point3::new(1.0, 1.0, -0.75)));
        let side = if agreement { 1.0 } else { -1.0 };
        assert!(has_vertex(
            &mesh,
            Point3::new(1.0, 1.0, -0.75 + side * 1.25)
        ));
    }
}

/// A mirrored frame, or a clockwise ring, builds the same region outward
/// wound: a set, not an orientation, is what the boundary states.
#[test]
fn a_mirrored_frame_or_a_clockwise_ring_stays_outward_wound() {
    let clockwise = Rings {
        outer: unit_square().outer.into_iter().rev().collect(),
        holes: Vec::new(),
    };
    let upright = placement(Vec3::X, Vec3::Y, Vec3::new(0.5, 0.0, 1.0), Point3::ZERO);
    let mirrored = placement(Vec3::X, -Vec3::Y, Vec3::new(0.5, 0.0, 1.0), Point3::ZERO);
    for (frame, below) in [(upright, false), (mirrored, true)] {
        for ring in [unit_square(), clockwise.clone()] {
            for agreement in [true, false] {
                let mesh = prism_along_z(&ring, plane_z(), frame, agreement, 2.0).expect("clip");
                assert!(
                    (volume(&mesh) - 2.0).abs() < 1e-12,
                    "mirrored {below}, agreement {agreement}: {}",
                    volume(&mesh)
                );
                // The mirrored frame puts the square at y in [-1, 0].
                let y = if below { -1.0 } else { 1.0 };
                assert!(has_vertex(&mesh, Point3::new(0.0, y, 0.0)));
            }
        }
    }
}

/// A prism along an axis in the plane never crosses it, and a frame whose
/// axes do not span space has no prism: both refused by name.
#[test]
fn a_prism_parallel_to_the_plane_or_a_degenerate_frame_is_refused() {
    // The boundary stands upright in x and z; its prism runs along the plane.
    for z in [Vec3::Y, Vec3::new(1.0, 2.0, 1e-12)] {
        let error = prism_along_z(
            &unit_square(),
            plane_z(),
            placement(Vec3::X, Vec3::Z, z, Point3::ZERO),
            true,
            1.0,
        )
        .expect_err("parallel to the plane");
        assert!(
            matches!(&error, GeomError::InvalidInput(m) if m.contains("parallel to the clip plane")),
            "{error:?}"
        );
    }
    for (x, y, z) in [
        (Vec3::X, Vec3::Y, Vec3::new(1.0, 1.0, 0.0)),
        (Vec3::X, Vec3::X * 2.0, Vec3::Z),
        (Vec3::ZERO, Vec3::Y, Vec3::Z),
    ] {
        let error = BoundaryPrism::along_placement_z(placement(x, y, z, Point3::ZERO), tol())
            .expect_err("degenerate");
        assert!(
            matches!(&error, GeomError::InvalidInput(m) if m.contains("degenerate")),
            "{error:?}"
        );
    }
    let frame = placement(Vec3::X, Vec3::Y, Vec3::Z, Point3::ZERO);
    for reach in [0.0, -1.0, Scalar::INFINITY] {
        let error =
            prism_along_z(&unit_square(), plane_z(), frame, true, reach).expect_err("no depth");
        assert!(
            matches!(&error, GeomError::InvalidInput(m) if m.contains("slab depth")),
            "{error:?}"
        );
    }
}
