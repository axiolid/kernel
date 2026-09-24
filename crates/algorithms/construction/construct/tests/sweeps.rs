//! The six remaining sweep families against closed-form volume.
//!
//! Each family is checked against a constant from outside this crate rather
//! than against its own output: a tapered extrusion is a prismatoid, a
//! swept disk on a straight path is a cylinder, and a sectioned spine with
//! equal sections is a plain extrusion.

use axiolid_construct::profile::profile_rings;
use axiolid_construct::sweep;
use axiolid_core::{Point3, Scalar, Tolerance, Vec3};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_profile::{Profile, RectangleProfile};

fn tol_for(chord: Scalar) -> Tolerance {
    Tolerance::new(chord, 1e-9).expect("tolerance")
}

fn volume(mesh: &TriMesh, tol: Tolerance) -> Scalar {
    volume_properties(mesh, tol)
        .expect("a swept solid must be closed and two-manifold")
        .signed_volume
}

fn rect(x: Scalar, y: Scalar) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

#[test]
fn a_tapered_extrusion_matches_the_prismatoid_formula() {
    // Prismatoid: V = h(A1 + 4Am + A2)/6. For a linear taper between two
    // rectangles the mid-section is the average, so this is a real
    // constraint and not a restatement of the taper.
    let chord = 1e-6;
    let tol = tol_for(chord);
    let a = profile_rings(&rect(2.0, 2.0), chord, tol).expect("start");
    let b = profile_rings(&rect(4.0, 4.0), chord, tol).expect("end");
    let mesh = sweep::tapered_extrude(&a, &b, Vec3::Z, 3.0).expect("taper");
    // A1 = 4, A2 = 16, mid rectangle is 3x3 = 9.
    let want = 3.0 * (4.0 + 4.0 * 9.0 + 16.0) / 6.0;
    let got = volume(&mesh, tol);
    assert!(
        (got - want).abs() / want < 1e-9,
        "prismatoid volume {got} vs {want}"
    );
}

#[test]
fn a_tapered_revolution_averages_its_two_profiles() {
    // A linear taper between equal profiles must reproduce the plain
    // revolution exactly: same Pappus volume, no taper contribution.
    let chord = 1e-5;
    let tol = tol_for(chord);
    let mut a = profile_rings(&rect(1.0, 2.0), chord, tol).expect("start");
    for p in &mut a.outer {
        p.x += 4.0;
    }
    let b = a.clone();
    let mesh = sweep::tapered_revolve(&a, &b, Point3::ZERO, Vec3::Y, core::f64::consts::PI, tol)
        .expect("taper");
    // Pappus for the half turn: 2*pi*R*A / 2 with R = 4, A = 2.
    let want = core::f64::consts::TAU * 4.0 * 2.0 / 2.0;
    let got = volume(&mesh, tol);
    assert!(
        (want - got) / want < 1e-4 && got > 0.0,
        "tapered revolution {got} vs {want}"
    );
}

#[test]
fn a_swept_disk_on_a_straight_path_is_a_cylinder() {
    let chord = 1e-6;
    let tol = tol_for(chord);
    let path = [Point3::ZERO, Point3::new(0.0, 0.0, 5.0)];
    let mesh = sweep::swept_disk(&path, 2.0, None, None, tol).expect("disk");
    let want = core::f64::consts::PI * 4.0 * 5.0;
    let got = volume(&mesh, tol);
    // Inscribed, so it under-estimates; assert the direction too.
    assert!(
        (want - got) / want < 1e-4 && got < want,
        "swept disk {got} vs {want}"
    );
}

/// A path of straight legs joined by finely sampled circular bends.
///
/// `legs` are the unit directions of successive straight runs, each
/// `leg_length` long; consecutive legs are joined by a circular bend of
/// radius `bend` sampled every `step` radians. Returns the points and the
/// exact polyline length, which is what a sweep's volume scales with.
fn bent_path(
    legs: &[Vec3],
    leg_length: Scalar,
    bend: Scalar,
    step: Scalar,
) -> (Vec<Point3>, Scalar) {
    let mut points = vec![Point3::ZERO];
    let mut here = Point3::ZERO;
    for (i, dir) in legs.iter().enumerate() {
        here += *dir * leg_length;
        points.push(here);
        let Some(next) = legs.get(i + 1) else { break };
        // Bend in the plane of `dir` and `next`, centred on the inside.
        let angle = dir.angle_between(*next);
        let inward = (*next - *dir * dir.dot(*next)).normalize();
        let centre = here + inward * bend;
        let count = (angle / step).ceil().max(1.0) as usize;
        for k in 1..=count {
            let a = angle * k as Scalar / count as Scalar;
            here = centre - inward * (bend * a.cos()) + *dir * (bend * a.sin());
            points.push(here);
        }
    }
    let length = points.windows(2).map(|w| (w[1] - w[0]).length()).sum();
    (points, length)
}

/// Closed-form volume of the swept disk the compiler builds: its inscribed
/// n-gon area (n from the same chord budget) times the path length.
fn disk_volume(radius: Scalar, inner: Option<Scalar>, length: Scalar, chord: Scalar) -> Scalar {
    let ngon = |r: Scalar| {
        let per = 2.0
            * (1.0 - (chord / r).min(1.0))
                .clamp(-1.0, 1.0)
                .acos()
                .max(1e-9);
        let n = (core::f64::consts::TAU / per).ceil().clamp(2.0, 4096.0);
        0.5 * n * r * r * (core::f64::consts::TAU / n).sin()
    };
    (ngon(radius) - inner.map_or(0.0, ngon)) * length
}

fn assert_bent_disk(legs: &[Vec3], inner: Option<Scalar>, what: &str) {
    let chord = 1e-4;
    let tol = tol_for(chord);
    let (path, length) = bent_path(legs, 3.0, 1.0, 0.02);
    let mesh = sweep::swept_disk(&path, 0.25, inner, None, tol)
        .unwrap_or_else(|e| panic!("{what}: a bent disk must sweep: {e}"));
    let want = disk_volume(0.25, inner, length, chord);
    let got = volume(&mesh, tol);
    assert!(
        (got - want).abs() / want < 1e-3,
        "{what}: swept {got} vs n-gon area x path length {want}"
    );
}

/// A pipe that starts along +X and turns onto +Y.
///
/// A path starting along X seeds the reference on Y. A single fixed
/// reference is then tangent to the second leg and the sweep was refused,
/// although a circle needs no particular orientation (axiolid/kernel#169).
#[test]
fn a_swept_disk_turning_onto_its_seed_axis_sweeps() {
    assert_bent_disk(&[Vec3::X, Vec3::Y], None, "L, solid");
    assert_bent_disk(&[Vec3::X, Vec3::Y], Some(0.1), "L, hollow");
}

/// A U: out along X, across along Y, back along -X.
#[test]
fn a_u_shaped_swept_disk_sweeps() {
    assert_bent_disk(&[Vec3::X, Vec3::Y, -Vec3::X], None, "U, solid");
    assert_bent_disk(&[Vec3::X, Vec3::Y, -Vec3::X], Some(0.1), "U, hollow");
}

/// A leg ALMOST along the seed axis is the silent case.
///
/// Exactly parallel is refused; off by 1e-9 it is not, and projecting the
/// fixed reference onto that leg's normal plane leaves a residue of 1e-9
/// whose direction is arbitrary. The ring rotates about the path between
/// neighbouring stations, the loft connects vertex k to vertex k across the
/// rotation, and the volume collapses with no error. A bent rebar in a real
/// Revit model lost up to 64 % of its volume this way.
#[test]
fn a_leg_nearly_along_the_seed_axis_does_not_twist() {
    let nearly_y = Vec3::new(1e-9, 1.0, 0.0).normalize();
    assert_bent_disk(&[Vec3::X, nearly_y], None, "nearly-parallel leg");
    assert_bent_disk(
        &[Vec3::Z, Vec3::X, nearly_y],
        Some(0.1),
        "3D, nearly-parallel leg",
    );
}

/// A helix has torsion; a rotation-minimising frame must still close the
/// volume to the analytic value (no twist-induced shrinkage).
#[test]
fn a_helical_swept_disk_keeps_its_volume() {
    let chord = 1e-4;
    let tol = tol_for(chord);
    let path: Vec<Point3> = (0..=400)
        .map(|k| {
            let a = k as Scalar * 0.02;
            Point3::new(2.0 * a.cos(), 2.0 * a.sin(), 0.3 * a)
        })
        .collect();
    let length: Scalar = path.windows(2).map(|w| (w[1] - w[0]).length()).sum();
    let mesh = sweep::swept_disk(&path, 0.25, None, None, tol).expect("helix sweeps");
    let want = disk_volume(0.25, None, length, chord);
    let got = volume(&mesh, tol);
    assert!(
        (got - want).abs() / want < 1e-3,
        "helix: swept {got} vs {want}"
    );
}

#[test]
fn a_hollow_swept_disk_subtracts_its_bore() {
    let chord = 1e-6;
    let tol = tol_for(chord);
    let path = [Point3::ZERO, Point3::new(0.0, 0.0, 5.0)];
    let mesh = sweep::swept_disk(&path, 2.0, Some(1.0), None, tol).expect("pipe");
    // An annulus of radii 2 and 1: pi(4 - 1)*5.
    let want = core::f64::consts::PI * 3.0 * 5.0;
    let got = volume(&mesh, tol);
    assert!(
        (want - got).abs() / want < 1e-3 && got > 0.0,
        "hollow swept disk {got} vs {want}"
    );
}

#[test]
fn a_fillet_radius_is_refused_not_silently_sharpened() {
    // The model's own docs require this: a consumer that cannot round
    // corners must refuse rather than drop the request, because a silently
    // sharpened pipe run builds, renders, and is wrong.
    let tol = tol_for(1e-6);
    let path = [Point3::ZERO, Point3::new(0.0, 0.0, 5.0)];
    assert!(sweep::swept_disk(&path, 2.0, None, Some(0.5), tol).is_err());
}

#[test]
fn a_fixed_reference_sweep_on_a_straight_path_is_an_extrusion() {
    let chord = 1e-6;
    let tol = tol_for(chord);
    let rings = profile_rings(&rect(2.0, 3.0), chord, tol).expect("rings");
    let path = [Point3::ZERO, Point3::new(0.0, 0.0, 4.0)];
    let mesh = sweep::fixed_reference_sweep(&rings, &path, Vec3::X).expect("sweep");
    let want = 2.0 * 3.0 * 4.0;
    let got = volume(&mesh, tol);
    assert!(
        (got - want).abs() / want < 1e-9,
        "fixed reference sweep {got} vs {want}"
    );
}

#[test]
fn a_parallel_reference_is_refused() {
    // A reference parallel to the tangent cannot orient the profile.
    // Substituting a fallback axis would rotate the section by an
    // arbitrary angle, so it must be refused.
    let chord = 1e-6;
    let tol = tol_for(chord);
    let rings = profile_rings(&rect(2.0, 3.0), chord, tol).expect("rings");
    let path = [Point3::ZERO, Point3::new(0.0, 0.0, 4.0)];
    assert!(sweep::fixed_reference_sweep(&rings, &path, Vec3::Z).is_err());
}

#[test]
fn a_surface_curve_sweep_takes_its_up_from_the_surface() {
    let chord = 1e-6;
    let tol = tol_for(chord);
    let rings = profile_rings(&rect(2.0, 3.0), chord, tol).expect("rings");
    let path = [Point3::ZERO, Point3::new(0.0, 0.0, 4.0)];
    let normals = [Vec3::X, Vec3::X];
    let mesh = sweep::surface_curve_sweep(&rings, &path, &normals).expect("sweep");
    let want = 2.0 * 3.0 * 4.0;
    let got = volume(&mesh, tol);
    assert!(
        (got - want).abs() / want < 1e-9,
        "surface curve sweep {got} vs {want}"
    );
    // One normal per directrix point, or there is no correspondence.
    assert!(sweep::surface_curve_sweep(&rings, &path, &[Vec3::X]).is_err());
}

#[test]
fn a_sectioned_spine_with_equal_sections_is_an_extrusion() {
    let chord = 1e-6;
    let tol = tol_for(chord);
    let rings = profile_rings(&rect(2.0, 3.0), chord, tol).expect("rings");
    let place = |z: Scalar| -> Vec<Point3> {
        rings
            .outer
            .iter()
            .map(|p| Point3::new(p.x, p.y, z))
            .collect()
    };
    let sections = vec![
        (rings.clone(), place(0.0)),
        (rings.clone(), place(2.0)),
        (rings.clone(), place(4.0)),
    ];
    let mesh = sweep::sectioned_spine(&sections).expect("spine");
    // Three equal sections spanning 4 must give the same solid as one
    // extrusion of depth 4: intermediate stations add vertices, not volume.
    let want = 2.0 * 3.0 * 4.0;
    let got = volume(&mesh, tol);
    assert!(
        (got - want).abs() / want < 1e-9,
        "sectioned spine {got} vs {want}"
    );
    assert!(sweep::sectioned_spine(&sections[..1]).is_err());
}
