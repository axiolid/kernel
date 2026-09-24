//! Gates for profile flattening and extrusion.
//!
//! Signed volume is the single check that catches wrong winding, missing caps,
//! and inverted sides simultaneously: it is positive exactly when the solid is
//! closed and outward-oriented, and its magnitude is `area * depth`.

use axiolid_construct::extrude::{extrude, extrude_profile, outward_orientation};
use axiolid_construct::profile::{profile_rings, triangulate, Rings};
use axiolid_core::{Tolerance, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

/// Six times the signed volume, via the divergence theorem.
fn six_volume(m: &TriMesh) -> f64 {
    m.indices
        .chunks_exact(3)
        .map(|t| {
            let (a, b, c) = (
                m.positions[t[0] as usize],
                m.positions[t[1] as usize],
                m.positions[t[2] as usize],
            );
            a.dot(b.cross(c))
        })
        .sum()
}

fn volume(m: &TriMesh) -> f64 {
    six_volume(m) / 6.0
}

/// Loop ranges for a `Rings` laid out as outer ++ holes.
#[allow(clippy::single_range_in_vec_init)]
fn loops(rings: &Rings) -> Vec<core::ops::Range<usize>> {
    let mut out = vec![0..rings.outer.len()];
    let mut start = rings.outer.len();
    for h in &rings.holes {
        out.push(start..start + h.len());
        start += h.len();
    }
    out
}

fn rect(x: f64, y: f64, thickness: Option<f64>) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness,
        outer_radius: None,
        inner_radius: None,
    })
}

/// A solid box: volume must equal x * y * depth, and be POSITIVE.
///
/// A negative result means the solid is inside-out, which `axiolid-mesh-boolean-boolmesh`
/// refuses -- so this assertion is what keeps the two crates composable.
#[test]
fn a_box_extrudes_to_its_exact_volume_with_outward_winding() {
    let rings = profile_rings(&rect(4.0, 0.2, None), 1e-3, Tolerance::METRE).expect("rings");
    let (pts, tris) = triangulate(&rings).expect("triangulate");
    let mesh = extrude(&pts, &tris, &loops(&rings), Vec3::Z, 3.0).expect("extrude");

    let v = volume(&mesh);
    assert!(v > 0.0, "extruded solid must be outward-oriented, got {v}");
    assert!((v - 4.0 * 0.2 * 3.0).abs() < 1e-9, "expected 2.4, got {v}");
    assert!(mesh.validate_structure().is_ok());
}

/// An extrusion direction with a NEGATIVE z component is still outward.
///
/// `extrude` builds its caps and side walls assuming the offset leaves the
/// profile plane towards +z. Real building models routinely author
/// `ExtrudedDirection = (0, 0, -1)`, e.g. an opening body hung from a lintel
/// and extruded downward. Without compensation every triangle is wound
/// inside-out, the signed volume is `-area * depth`, and any boolean that
/// uses the solid refuses it (or, worse, silently inverts). The magnitude
/// must be unchanged and the sign positive whichever way the offset points.
#[test]
fn a_downward_extrusion_is_still_outward() {
    let rings = profile_rings(&rect(4.0, 0.2, None), 1e-3, Tolerance::METRE).expect("rings");
    let (pts, tris) = triangulate(&rings).expect("triangulate");
    for direction in [
        Vec3::NEG_Z,
        Vec3::new(0.0, 0.3, -1.0),
        Vec3::new(0.5, -0.2, -1.0),
    ] {
        let mesh = extrude(&pts, &tris, &loops(&rings), direction, 3.0).expect("extrude");
        let v = volume(&mesh);
        let expected = 4.0 * 0.2 * 3.0 * direction.normalize().z.abs();
        assert!(
            v > 0.0,
            "direction {direction:?}: solid must be outward-oriented, got {v}"
        );
        assert!(
            (v - expected).abs() < 1e-9,
            "direction {direction:?}: expected {expected}, got {v}"
        );
        assert!(mesh.validate_structure().is_ok());
    }
}

/// A hollow section extruded downward keeps the hole subtracted, not added.
///
/// The inner wall faces the opposite way from the outer one, so a fix that
/// flips only the caps or only the outer loop would still get this wrong.
#[test]
fn a_downward_hollow_section_loses_exactly_the_hole_volume() {
    let rings = profile_rings(&rect(4.0, 2.0, Some(0.25)), 1e-3, Tolerance::METRE).expect("rings");
    let (pts, tris) = triangulate(&rings).expect("triangulate");
    let mesh = extrude(&pts, &tris, &loops(&rings), Vec3::NEG_Z, 3.0).expect("extrude");
    let v = volume(&mesh);
    let expected = (4.0 * 2.0 - 3.5 * 1.5) * 3.0;
    assert!(v > 0.0, "hollow section must be outward-oriented, got {v}");
    assert!((v - expected).abs() < 1e-9, "expected {expected}, got {v}");
}

/// A hollow section must lose exactly the hole's volume. This is the case the
/// hand-rolled ear clipper could not do, and the reason earcut was adopted.
#[test]
fn a_hollow_section_loses_exactly_the_hole_volume() {
    let rings = profile_rings(&rect(10.0, 6.0, Some(1.0)), 1e-3, Tolerance::METRE).expect("rings");
    assert_eq!(rings.holes.len(), 1, "hollow rectangle must produce a hole");

    let (pts, tris) = triangulate(&rings).expect("triangulate");
    let mesh = extrude(&pts, &tris, &loops(&rings), Vec3::Z, 2.0).expect("extrude");

    // Outer 10x6, inner 8x4 -> ring area 60 - 32 = 28.
    let expected = 28.0 * 2.0;
    let v = volume(&mesh);
    assert!(v > 0.0, "hollow solid must be outward-oriented, got {v}");
    assert!((v - expected).abs() < 1e-9, "expected {expected}, got {v}");
}

/// A disk's volume converges on pi r^2 h as the chord budget tightens, and the
/// approximation is always an INSCRIBED polygon, so it never overshoots.
#[test]
fn a_disk_converges_from_below_as_the_chord_budget_tightens() {
    let exact = core::f64::consts::PI * 4.0 * 1.5;
    let mut previous = 0.0;
    for chord in [1e-1, 1e-2, 1e-3, 1e-4] {
        let profile = Profile::Circle(CircleProfile {
            radius: 2.0,
            thickness: None,
        });
        let rings = profile_rings(&profile, chord, Tolerance::METRE).expect("rings");
        let (pts, tris) = triangulate(&rings).expect("triangulate");
        let mesh = extrude(&pts, &tris, &loops(&rings), Vec3::Z, 1.5).expect("extrude");
        let v = volume(&mesh);

        assert!(
            v > 0.0 && v < exact,
            "inscribed volume {v} must be under {exact}"
        );
        assert!(v > previous, "tightening the budget must not lose volume");
        previous = v;
    }
    // The segment count is clamped at 512, so the tightest achievable relative
    // volume error is set by that clamp, not by the requested chord budget:
    // an inscribed n-gon loses O(1/n^2) of the area. 512 segments gives ~2.5e-5.
    // Asserting a tighter bound would be asserting against the clamp we chose.
    assert!(
        (exact - previous).abs() / exact < 1e-4,
        "tightest budget should be within 1e-4 relative, got {previous} vs {exact}"
    );
}

/// Every directed edge must appear exactly once.
///
/// This is the gate that volume CANNOT provide: a cap lying in the z = 0 plane
/// contributes nothing to the divergence integral, so a flipped bottom cap is
/// invisible to a volume check. Directed-edge parity is sensitive to winding
/// everywhere, and it is precisely what `boolmesh::Manifold::new` enforces on
/// input -- so this test is the local proxy for "will the boolean accept it".
fn assert_edge_manifold(mesh: &TriMesh, what: &str) {
    use std::collections::HashMap;
    let mut directed: HashMap<(u32, u32), i32> = HashMap::new();
    for t in mesh.indices.chunks_exact(3) {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *directed.entry((a, b)).or_default() += 1;
        }
    }
    for (&(a, b), &count) in &directed {
        assert_eq!(
            count, 1,
            "{what}: directed edge {a}->{b} appears {count} times"
        );
        let opposite = directed.get(&(b, a)).copied().unwrap_or(0);
        assert_eq!(
            opposite, 1,
            "{what}: edge {a}->{b} has {opposite} opposing half-edges, so the \
             surface is not closed and consistently wound"
        );
    }
}

#[test]
fn a_box_is_edge_manifold() {
    let rings = profile_rings(&rect(4.0, 0.2, None), 1e-3, Tolerance::METRE).expect("rings");
    let (pts, tris) = triangulate(&rings).expect("triangulate");
    let mesh = extrude(&pts, &tris, &loops(&rings), Vec3::Z, 3.0).expect("extrude");
    assert_edge_manifold(&mesh, "solid box");
}

#[test]
fn a_hollow_section_is_edge_manifold() {
    let rings = profile_rings(&rect(10.0, 6.0, Some(1.0)), 1e-3, Tolerance::METRE).expect("rings");
    let (pts, tris) = triangulate(&rings).expect("triangulate");
    let mesh = extrude(&pts, &tris, &loops(&rings), Vec3::Z, 2.0).expect("extrude");
    assert_edge_manifold(&mesh, "hollow section");
}

/// Contour profiles with holes, and mirrored placements.
///
/// The committed IFC corpus contains neither, so these paths would otherwise
/// ship unverified: mutation probes on `contour_points` and the mirror
/// re-orientation both survived the corpus gate.
mod contour_and_mirror {
    use super::{assert_edge_manifold, volume};

    use axiolid_construct::extrude::extrude_profile;
    use axiolid_construct::profile::profile_rings;
    use axiolid_core::{Point2, Scalar, Tolerance, Transform2, Vec2, Vec3};
    use axiolid_curve::{Curve2, Polyline2};
    use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile};

    /// One closed polyline segment covering a whole ring.
    fn ring(points: Vec<Point2>) -> Contour {
        // A polyline's parameter is one unit per segment, so the domain of a
        // closed n-point ring is (0, n). The fixture previously declared
        // (0, 1) and the old flattener ignored the domain entirely; the
        // scalar evaluator honours it, so the declared domain must be right.
        let segments = points.len() as Scalar;
        Contour::new(vec![ProfileSegment {
            curve: Curve2::Polyline(Polyline2 {
                points,
                closed: true,
            }),
            domain: axiolid_core::Interval::new(0.0, segments),
            same_sense: true,
        }])
    }

    fn square(cx: Scalar, cy: Scalar, half: Scalar) -> Vec<Point2> {
        vec![
            Point2::new(cx - half, cy - half),
            Point2::new(cx + half, cy - half),
            Point2::new(cx + half, cy + half),
            Point2::new(cx - half, cy + half),
        ]
    }

    /// A ring whose first point is repeated as its last must still close.
    ///
    /// Authoring tools write both conventions; a duplicated closing point is a
    /// zero-length edge that earcut cannot triangulate.
    #[test]
    fn a_ring_that_repeats_its_first_point_is_closed_once() {
        let mut pts = square(0.0, 0.0, 1.0);
        pts.push(pts[0]);
        let profile = Profile::Contour(ContourProfile {
            outer: ring(pts),
            holes: Vec::new(),
        });
        let rings = profile_rings(&profile, 1e-4, Tolerance::MILLIMETRE).expect("rings");
        assert_eq!(rings.outer.len(), 4, "the duplicate closing point must go");
        let mesh = extrude_profile(&rings, Vec3::Z, 2.0, Tolerance::MILLIMETRE).expect("extrude");
        assert!((volume(&mesh) - 8.0).abs() < 1e-9);
        assert_edge_manifold(&mesh, "closed ring solid");
    }

    /// A contour with a hole must lose exactly the hole's volume.
    #[test]
    fn a_contour_hole_is_subtracted() {
        let profile = Profile::Contour(ContourProfile {
            outer: ring(square(0.0, 0.0, 2.0)),
            holes: vec![ring(square(0.0, 0.0, 0.5))],
        });
        let rings = profile_rings(&profile, 1e-4, Tolerance::MILLIMETRE).expect("rings");
        assert_eq!(rings.holes.len(), 1, "the hole must survive flattening");
        let mesh = extrude_profile(&rings, Vec3::Z, 3.0, Tolerance::MILLIMETRE).expect("extrude");
        // outer 4x4 minus hole 1x1, times depth 3.
        assert!((volume(&mesh) - (16.0 - 1.0) * 3.0).abs() < 1e-9);
        assert_edge_manifold(&mesh, "contour solid");
    }

    /// A mirroring 2D placement must still yield an outward-facing solid.
    #[test]
    fn a_mirrored_derived_profile_stays_outward() {
        let basis = Profile::Rectangle(RectangleProfile {
            x: 2.0,
            y: 1.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        });
        let mirror = Transform2::from_scale_angle_translation(
            Vec2::new(-1.0, 1.0),
            0.0,
            Vec2::new(5.0, 0.0),
        );
        let profile = Profile::Derived {
            basis: Box::new(basis),
            transform: mirror,
        };
        let rings = profile_rings(&profile, 1e-4, Tolerance::MILLIMETRE).expect("rings");
        let mesh = extrude_profile(&rings, Vec3::Z, 2.0, Tolerance::MILLIMETRE).expect("extrude");
        assert!(
            volume(&mesh) > 0.0,
            "a mirrored profile produced an inside-out solid"
        );
        assert!((volume(&mesh) - 4.0).abs() < 1e-9);
        assert_edge_manifold(&mesh, "mirrored solid");
    }
}

/// The certified orientation check must agree with the volume sign on solids
/// this crate produces, and must decline the cases volume cannot judge.
///
/// This is ADR 0016's follow-up: `orient3d` becomes load-bearing here rather
/// than staying a proven-but-unused capability.
#[test]
fn certified_orientation_agrees_with_volume_on_real_solids() {
    for (label, profile, depth) in [
        ("solid box", rect(4.0, 0.2, None), 3.0),
        ("hollow section", rect(10.0, 6.0, Some(1.0)), 2.0),
        ("thin plate", rect(50.0, 50.0, None), 0.001),
    ] {
        let rings = profile_rings(&profile, 1e-4, Tolerance::METRE).expect("rings");
        let mesh = extrude_profile(&rings, Vec3::Z, depth, Tolerance::METRE).expect("extrude");

        let certified = outward_orientation(&mesh)
            .unwrap_or_else(|| panic!("{label}: a valid solid must be judgeable"));
        assert!(certified, "{label}: extrusion output must be outward");
        assert!(
            volume(&mesh) > 0.0,
            "{label}: the volume sign must agree with the certified one"
        );

        // Reversing every triangle must flip the verdict, not merely fail.
        let mut inverted = mesh.clone();
        for corner in inverted.indices.chunks_exact_mut(3) {
            corner.swap(1, 2);
        }
        assert_eq!(
            outward_orientation(&inverted),
            Some(false),
            "{label}: an inverted solid must be reported inward"
        );
    }
}

/// Degenerate input must be declined rather than assigned an orientation.
#[test]
fn an_unjudgeable_mesh_is_declined() {
    // Too few faces to bound a volume.
    let sliver = TriMesh::new(
        vec![
            axiolid_core::Point3::new(0.0, 0.0, 0.0),
            axiolid_core::Point3::new(1.0, 0.0, 0.0),
            axiolid_core::Point3::new(0.0, 1.0, 0.0),
        ],
        vec![0, 1, 2],
    );
    assert_eq!(outward_orientation(&sliver), None);

    // Four coincident-plane faces enclose nothing: every tetrahedron about the
    // centroid is degenerate, so no side can be certified.
    let flat = TriMesh::new(
        vec![
            axiolid_core::Point3::new(0.0, 0.0, 0.0),
            axiolid_core::Point3::new(1.0, 0.0, 0.0),
            axiolid_core::Point3::new(1.0, 1.0, 0.0),
            axiolid_core::Point3::new(0.0, 1.0, 0.0),
        ],
        vec![0, 1, 2, 0, 2, 3, 2, 1, 0, 3, 2, 0],
    );
    assert_eq!(outward_orientation(&flat), None);
}

/// A thin solid far from the origin must still be judged correctly.
///
/// This is the case exact summation exists for. Building models are routinely
/// placed on survey coordinates in the 1e6 range, and a 0.1 mm sheet there
/// produces per-tetrahedron terms around 1e12 whose true sum is around 1e-1.
/// Summing those in f64 loses roughly sixteen digits, so the sign of the
/// result is noise. A mutation replacing the exact accumulation with an f64
/// sum survives every other test in this file; it fails here.
#[test]
fn a_thin_plate_on_survey_coordinates_is_still_judged_correctly() {
    let rings = profile_rings(&rect(2.0, 2.0, None), 1e-4, Tolerance::METRE).expect("rings");
    let mut mesh = extrude_profile(&rings, Vec3::Z, 1e-4, Tolerance::METRE).expect("extrude");

    // Translate onto national-grid style coordinates.
    let offset = axiolid_core::Point3::new(4.5e6, 3.2e6, 1.0e6);
    for p in &mut mesh.positions {
        *p = axiolid_core::Point3::new(p.x + offset.x, p.y + offset.y, p.z + offset.z);
    }

    // Precondition: the naive f64 sum has genuinely lost the sign, so this
    // test is exercising the exact path rather than restating the easy case.
    let naive: f64 = mesh
        .indices
        .chunks_exact(3)
        .map(|t| {
            let (a, b, c) = (
                mesh.positions[t[0] as usize],
                mesh.positions[t[1] as usize],
                mesh.positions[t[2] as usize],
            );
            a.dot(b.cross(c))
        })
        .sum();
    let true_volume = 2.0 * 2.0 * 1e-4;
    assert!(
        (naive / 6.0 - true_volume).abs() > true_volume * 0.5,
        "precondition: the naive sum should be badly wrong here, got {}",
        naive / 6.0
    );

    assert_eq!(
        outward_orientation(&mesh),
        Some(true),
        "a thin plate on survey coordinates must still be outward"
    );

    let mut inverted = mesh.clone();
    for corner in inverted.indices.chunks_exact_mut(3) {
        corner.swap(1, 2);
    }
    assert_eq!(outward_orientation(&inverted), Some(false));
}
