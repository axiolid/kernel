//! Extrusion checked against closed-form volume (C3).
//!
//! # Why this file exists
//!
//! `extrude` was already implemented and watertight, but nothing asserted the
//! one identity that makes it *correct* rather than merely closed:
//!
//! ```text
//! volume(extrude(profile, depth)) == area(profile) * depth
//! ```
//!
//! Both sides are computed independently here: the left by summing signed
//! tetrahedra over the produced mesh, the right by the shoelace formula over
//! the flattened rings. A winding error, a dropped cap, or a mis-stitched side
//! wall all break the identity.
//!
//! These tests also cover the profile families that only became reachable once
//! curve evaluation moved into `axiolid-reference`: ellipse and B-spline contours
//! were a hard `Unsupported` before.

use axiolid_construct::extrude::extrude_profile;
use axiolid_construct::profile::profile_rings;
use axiolid_core::{Frame2, Interval, Point2, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{BSplineCurve2, Curve2, Ellipse2, KnotSpec, Polyline2};
use axiolid_mesh::TriMesh;
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, EllipseProfile, Profile, ProfileSegment,
    RectangleProfile,
};
use axiolid_reference::signed_area2;

const TAU: Scalar = core::f64::consts::TAU;

/// Signed volume via the measurement oracle (`axiolid-measure`).
///
/// Deliberately not a local divergence sum. `volume_properties` refuses a mesh
/// that is not closed and two-manifold *before* it integrates, so every call
/// here asserts watertightness as well as volume. A hand-rolled integral
/// cannot: it returns a plausible number for a mesh with a missing face.
fn volume(mesh: &TriMesh) -> Scalar {
    volume_at(mesh, Tolerance::MILLIMETRE)
}

/// Measure with an explicit tolerance.
///
/// `audit_mesh` calls a triangle degenerate when its doubled area is at most
/// `tolerance.linear()`, so the tolerance must match the scale of the mesh
/// being measured. A cylinder flattened at chord budget `h` has side triangles
/// about `h` wide: judging those with a fixed MILLIMETRE tolerance would
/// classify perfectly correct geometry as degenerate as soon as the caller
/// asks for sub-millimetre accuracy. Tolerance is a property of the
/// measurement, not a constant.
fn volume_at(mesh: &TriMesh, tolerance: Tolerance) -> Scalar {
    axiolid_measure::volume_properties(mesh, tolerance)
        .expect("extrusion output must be closed and two-manifold")
        .signed_volume
}

/// A tolerance proportional to a chord budget, floored at f64 sanity.
fn tolerance_for(chord: Scalar) -> Tolerance {
    Tolerance::new((chord * 1e-3).max(1e-12), 1e-9).expect("valid tolerance")
}

/// Cross-section area from the flattened rings: outer minus holes.
///
/// `signed_area2` returns *twice* the signed area (the raw shoelace sum), so
/// the halving is not a fudge factor -- it is the definition.
fn ring_area(rings: &axiolid_construct::profile::Rings) -> Scalar {
    let mut a = signed_area2(&rings.outer).abs();
    for hole in &rings.holes {
        a -= signed_area2(hole).abs();
    }
    a * 0.5
}

/// Every interior edge must be shared by exactly two triangles.
fn assert_closed_manifold(mesh: &TriMesh, what: &str) {
    use std::collections::HashMap;
    let mut edges: HashMap<(u32, u32), i32> = HashMap::new();
    for c in mesh.indices.as_chunks::<3>().0 {
        for (a, b) in [(c[0], c[1]), (c[1], c[2]), (c[2], c[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            *edges.entry(key).or_default() += if a < b { 1 } else { -1 };
        }
    }
    for (edge, balance) in edges {
        assert_eq!(
            balance, 0,
            "{what}: edge {edge:?} is not shared by two oppositely-wound faces"
        );
    }
}

/// The identity this whole file exists to prove.
fn assert_volume_is_area_times_depth(profile: &Profile, depth: Scalar, what: &str) {
    let chord = 1e-5;
    let rings = profile_rings(profile, chord, Tolerance::MILLIMETRE).expect("rings");
    let mesh = extrude_profile(&rings, Vec3::Z, depth, Tolerance::MILLIMETRE).expect("extrude");

    let want = ring_area(&rings) * depth;
    let got = volume_at(&mesh, tolerance_for(chord));
    // Relative tolerance: a flattened curve's area is itself an approximation,
    // bounded by the chord budget, so an absolute epsilon would be wrong.
    let tol = (want.abs() * 1e-6).max(1e-9);
    assert!(
        (got - want).abs() < tol,
        "{what}: volume {got} != area*depth {want}"
    );
    assert!(got > 0.0, "{what}: extrusion must be outward-oriented");
    assert_closed_manifold(&mesh, what);
}

// --- parameterized profiles -------------------------------------------------

#[test]
fn a_rectangle_extrudes_to_its_analytic_volume() {
    let p = Profile::Rectangle(RectangleProfile {
        x: 3.0,
        y: 2.0,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    });
    // Independent of ring_area: 3 * 2 * 4 = 24 exactly.
    let rings = profile_rings(&p, 1e-5, Tolerance::MILLIMETRE).expect("rings");
    let mesh = extrude_profile(&rings, Vec3::Z, 4.0, Tolerance::MILLIMETRE).expect("extrude");
    assert!(
        (volume(&mesh) - 24.0).abs() < 1e-12,
        "got {}",
        volume(&mesh)
    );
    assert_volume_is_area_times_depth(&p, 4.0, "rectangle");
}

#[test]
fn a_circle_extrudes_to_a_cylinder_of_pi_r_squared_h() {
    let r = 2.0;
    let h = 5.0;
    let p = Profile::Circle(CircleProfile {
        radius: r,
        thickness: None,
    });
    let rings = profile_rings(&p, 1e-6, Tolerance::MILLIMETRE).expect("rings");
    let mesh = extrude_profile(&rings, Vec3::Z, h, Tolerance::MILLIMETRE).expect("extrude");
    let want = core::f64::consts::PI * r * r * h;
    let got = volume_at(&mesh, tolerance_for(1e-6));
    // A flattened circle is an inscribed polygon, so it under-estimates.
    // The gap must be small and one-signed -- that is the tolerance working.
    assert!(got < want, "inscribed polygon must under-estimate");
    assert!(
        (want - got) / want < 1e-6,
        "cylinder volume {got} vs {want} (relative {})",
        (want - got) / want
    );
    assert_volume_is_area_times_depth(&p, h, "circle");
}

#[test]
fn an_annulus_loses_exactly_its_hole() {
    let p = Profile::Circle(CircleProfile {
        radius: 3.0,
        thickness: Some(1.0),
    });
    assert_volume_is_area_times_depth(&p, 2.0, "annulus");

    // And the hole is really absent: compare against the filled disk.
    let filled = Profile::Circle(CircleProfile {
        radius: 3.0,
        thickness: None,
    });
    let hollow_rings = profile_rings(&p, 1e-6, Tolerance::MILLIMETRE).unwrap();
    let filled_rings = profile_rings(&filled, 1e-6, Tolerance::MILLIMETRE).unwrap();
    let hollow = volume_at(
        &extrude_profile(&hollow_rings, Vec3::Z, 2.0, Tolerance::MILLIMETRE).unwrap(),
        tolerance_for(1e-6),
    );
    let solid = volume_at(
        &extrude_profile(&filled_rings, Vec3::Z, 2.0, Tolerance::MILLIMETRE).unwrap(),
        tolerance_for(1e-6),
    );
    let inner_volume = core::f64::consts::PI * 2.0 * 2.0 * 2.0;
    assert!(
        ((solid - hollow) - inner_volume).abs() / inner_volume < 1e-5,
        "removed {}, expected {inner_volume}",
        solid - hollow
    );
}

// --- families unlocked by C1 ------------------------------------------------

#[test]
fn an_ellipse_profile_now_extrudes_at_all() {
    // Before curve evaluation moved into axiolid-reference this was a hard
    // `Unsupported`, so this test is the capability claim.
    let p = Profile::Ellipse(EllipseProfile {
        semi_axis_x: 3.0,
        semi_axis_y: 1.5,
    });
    let rings = profile_rings(&p, 1e-6, Tolerance::MILLIMETRE).expect("ellipse must be supported");
    let mesh = extrude_profile(&rings, Vec3::Z, 2.0, Tolerance::MILLIMETRE).expect("extrude");
    let want = core::f64::consts::PI * 3.0 * 1.5 * 2.0;
    let got = volume_at(&mesh, tolerance_for(1e-6));
    assert!(got < want, "inscribed polygon under-estimates");
    assert!(
        (want - got) / want < 1e-6,
        "elliptic cylinder {got} vs {want}"
    );
    assert_volume_is_area_times_depth(&p, 2.0, "ellipse");
}

#[test]
fn an_elliptical_arc_contour_extrudes() {
    // A contour built from an exact Ellipse2 curve segment, closed by a line.
    let e = Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::new(1.0, 0.0),
            y: Vec2::new(0.0, 1.0),
        },
        semi_axis_x: 2.0,
        semi_axis_y: 1.0,
    });
    let contour = Contour::new(vec![ProfileSegment {
        curve: e,
        domain: Interval::new(0.0, TAU),
        same_sense: true,
    }]);
    let p = Profile::Contour(ContourProfile {
        outer: contour,
        holes: Vec::new(),
    });
    assert_volume_is_area_times_depth(&p, 1.5, "elliptical contour");
}

#[test]
fn a_bspline_contour_extrudes() {
    // A closed rational quarter-arc plus straight returns: exercises the
    // de Boor path inside a real profile.
    let w = 1.0 / Scalar::sqrt(2.0);
    let arc = Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![3, 3],
        weights: Some(vec![1.0, w, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::PiecewiseBezier,
    });
    let back = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 1.0),
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
        ],
        closed: false,
    });
    let p = Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            ProfileSegment {
                curve: arc,
                domain: Interval::new(0.0, 1.0),
                same_sense: true,
            },
            ProfileSegment {
                curve: back,
                domain: Interval::new(0.0, 2.0),
                same_sense: true,
            },
        ]),
        holes: Vec::new(),
    });
    let rings = profile_rings(&p, 1e-6, Tolerance::MILLIMETRE).expect("bspline contour");
    // Quarter disk of radius 1: area pi/4.
    let area = ring_area(&rings);
    assert!(
        (area - core::f64::consts::FRAC_PI_4).abs() < 1e-5,
        "quarter-disc area {area}"
    );
    assert_volume_is_area_times_depth(&p, 3.0, "bspline contour");
}

// --- tolerance is honoured end to end ---------------------------------------

#[test]
fn a_tighter_chord_budget_converges_on_the_true_cylinder() {
    // The whole point of adaptive flattening: error must fall with tolerance,
    // and stay proportional to the budget the caller asked for.
    //
    // Measured behaviour (unit disc, depth 1):
    //
    // ```text
    //   chord     vertices   volume error   error/chord
    //   1e-2            32       2.015e-2         2.015
    //   1e-3           128       1.261e-3         1.261
    //   1e-4           256       3.154e-4         3.154
    //   1e-5          1024       1.971e-5         1.971
    //   1e-6          4096       1.232e-6         1.232
    //   1e-7          8192       3.080e-7         3.080
    // ```
    //
    // Error is O(chord) with a small constant, which is what an inscribed
    // polygon gives. Asserting a fixed final magnitude would be asserting a
    // vertex count; asserting the ratio is asserting the contract.
    let p = Profile::Circle(CircleProfile {
        radius: 1.0,
        thickness: None,
    });
    let exact = core::f64::consts::PI;
    let mut previous = Scalar::INFINITY;
    for chord in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6] {
        let rings = profile_rings(&p, chord, Tolerance::MILLIMETRE).unwrap();
        let mesh = extrude_profile(&rings, Vec3::Z, 1.0, Tolerance::MILLIMETRE).unwrap();
        let error = (exact - volume_at(&mesh, tolerance_for(chord))).abs();
        assert!(
            error < previous,
            "error must shrink monotonically: {error} !< {previous} at chord {chord}"
        );
        // The requested budget must actually bound the outcome, within the
        // small constant an inscribed polygon costs.
        assert!(
            error < chord * 5.0,
            "chord {chord} did not bound the volume error {error}"
        );
        previous = error;
    }
}

// --- refusals ---------------------------------------------------------------

#[test]
fn a_zero_radius_circle_is_refused() {
    let p = Profile::Circle(CircleProfile {
        radius: 0.0,
        thickness: None,
    });
    assert!(profile_rings(&p, 1e-5, Tolerance::MILLIMETRE).is_err());
}

#[test]
fn a_full_circle_ring_does_not_repeat_its_closing_vertex() {
    // A duplicated closing point is a zero-length edge; the extruder would
    // emit a degenerate side quad from it.
    let p = Profile::Circle(CircleProfile {
        radius: 1.0,
        thickness: None,
    });
    let rings = profile_rings(&p, 1e-3, Tolerance::MILLIMETRE).unwrap();
    let first = rings.outer[0];
    let last = *rings.outer.last().unwrap();
    assert!(
        (first - last).length() > 1e-12,
        "ring repeats its first vertex"
    );
}

/// The centroid of a prism sits at exactly half its depth.
///
/// Volume is blind to defects the centroid is not: reflecting a cap through
/// the mid-plane, or emitting the two caps at swapped heights, both preserve
/// the divergence integral while moving the centre of mass.
#[test]
fn an_extruded_prism_centres_at_half_its_depth() {
    let depth = 3.0;
    let p = Profile::Rectangle(RectangleProfile {
        x: 4.0,
        y: 2.0,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    });
    let rings = profile_rings(&p, 1e-4, Tolerance::MILLIMETRE).expect("rings");
    let mesh = extrude_profile(&rings, Vec3::Z, depth, Tolerance::MILLIMETRE).expect("extrude");
    let props = axiolid_measure::volume_properties(&mesh, tolerance_for(1e-4))
        .expect("prism must be closed and two-manifold");
    assert!(
        (props.centroid.z - depth * 0.5).abs() < 1e-12,
        "prism centroid z {} != half depth {}",
        props.centroid.z,
        depth * 0.5
    );
    // A rectangle centred on the origin keeps the centroid on the axis.
    assert!(props.centroid.x.abs() < 1e-12 && props.centroid.y.abs() < 1e-12);
}

/// A cylinder's total surface area converges on 2*pi*r^2 + 2*pi*r*h.
///
/// Volume exercises the caps and the interior; area exercises the *sides*.
/// An extruder emitting one triangle per boundary edge instead of two would
/// still be closed and could still integrate to a believable volume, but its
/// lateral area would be visibly short.
#[test]
fn a_cylinder_surface_area_converges_on_the_analytic_value() {
    let (r, h) = (1.0, 2.0);
    let p = Profile::Circle(CircleProfile {
        radius: r,
        thickness: None,
    });
    let rings = profile_rings(&p, 1e-6, Tolerance::MILLIMETRE).expect("rings");
    let mesh = extrude_profile(&rings, Vec3::Z, h, Tolerance::MILLIMETRE).expect("extrude");
    let props = axiolid_measure::surface_properties(&mesh, tolerance_for(1e-6))
        .expect("cylinder must be surface-usable");
    let analytic = 2.0 * core::f64::consts::PI * r * r + core::f64::consts::TAU * r * h;
    let rel = (props.area - analytic).abs() / analytic;
    assert!(
        rel < 1e-5,
        "cylinder area {} vs analytic {analytic} (rel {rel:.3e})",
        props.area
    );
}

/// The measurement oracle refuses an open mesh; a raw integral would not.
///
/// This is the whole reason these tests call `axiolid-measure` instead of
/// summing tetrahedra locally. Deleting one triangle leaves a mesh whose
/// divergence sum is still a finite, plausible number -- and which is not a
/// solid. `volume_properties` audits closed-two-manifoldness first, so the
/// defect surfaces as a refusal instead of a slightly-wrong volume.
#[test]
fn the_oracle_refuses_an_open_mesh_a_raw_integral_would_accept() {
    let p = Profile::Rectangle(RectangleProfile {
        x: 2.0,
        y: 2.0,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    });
    let rings = profile_rings(&p, 1e-4, Tolerance::MILLIMETRE).expect("rings");
    let mesh = extrude_profile(&rings, Vec3::Z, 1.0, Tolerance::MILLIMETRE).expect("extrude");

    // Sound to begin with.
    let whole = axiolid_measure::volume_properties(&mesh, tolerance_for(1e-4))
        .expect("intact extrusion is a solid")
        .signed_volume;
    assert!((whole - 4.0).abs() < 1e-9, "2x2x1 box, got {whole}");

    // Remove one triangle: still finite under a raw integral, no longer a solid.
    let mut holed = mesh.clone();
    holed.indices.truncate(holed.indices.len() - 3);

    let raw: Scalar = holed
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|c| {
            let a = holed.positions[c[0] as usize];
            let b = holed.positions[c[1] as usize];
            let d = holed.positions[c[2] as usize];
            a.dot(b.cross(d)) / 6.0
        })
        .sum();
    assert!(
        raw.is_finite(),
        "a raw integral reports a plausible number for a holed mesh: {raw}"
    );

    let refused = axiolid_measure::volume_properties(&holed, tolerance_for(1e-4));
    assert!(
        refused.is_err(),
        "the oracle must refuse an open mesh, got {refused:?}"
    );
}
