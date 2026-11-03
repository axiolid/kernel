//! The six remaining sweep families against closed-form volume.
//!
//! Each family is checked against a constant from outside this crate rather
//! than against its own output: a tapered extrusion is a prismatoid, a
//! swept disk on a straight path is a cylinder, and a sectioned spine with
//! equal sections is a plain extrusion.

use axiolid_construct::profile::{profile_rings, Rings};
use axiolid_construct::sweep;
use axiolid_contracts::GeomError;
use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec3};
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
fn a_tapered_revolution_orients_outward_for_every_angle_sign_and_axis_direction() {
    // Same fixture as `a_tapered_revolution_averages_its_two_profiles`
    // (equal start and end profiles, so Pappus applies exactly), but now
    // covering every combination of angle sign and axis direction (#221):
    // `tapered_revolve` shares `revolve::revolve`'s unconditional station
    // reversal, which only winds the walls outward for a positive angle.
    let chord = 1e-5;
    let tol = tol_for(chord);
    let mut a = profile_rings(&rect(1.0, 2.0), chord, tol).expect("start");
    for p in &mut a.outer {
        p.x += 4.0;
    }
    let b = a.clone();
    let (inner, outer) = (3.5, 4.5);
    let axis_origin = Point3::ZERO;

    for (direction, angle) in [
        (Vec3::Y, core::f64::consts::PI / 2.0),
        (Vec3::Y, -core::f64::consts::PI / 2.0),
        (Vec3::NEG_Y, core::f64::consts::PI / 2.0),
        (Vec3::NEG_Y, -core::f64::consts::PI / 2.0),
    ] {
        let mesh =
            sweep::tapered_revolve(&a, &b, axis_origin, direction, angle, tol).expect("taper");
        let want = angle.abs() * 4.0 * 2.0;
        let got = volume(&mesh, tol);
        assert!(
            got > 0.0,
            "direction {direction:?} angle {angle}: signed volume {got} must be positive"
        );
        assert!(
            (got - want).abs() / want < 1e-4,
            "direction {direction:?} angle {angle}: volume {got} vs {want}"
        );
        assert_wall_normals_radial(&mesh, axis_origin, direction, inner, outer, chord);
    }
}

/// Every wall triangle (all three vertices at the same radius from the
/// axis) must face away from the axis at the outer radius and towards it
/// at the inner radius. Cap and mixed triangles are skipped: their outward
/// sense is not a pure radial statement.
fn assert_wall_normals_radial(
    mesh: &TriMesh,
    axis_origin: Point3,
    axis_direction: Vec3,
    inner: Scalar,
    outer: Scalar,
    tol: Scalar,
) {
    let dir = axis_direction / axis_direction.length();
    let radial = |p: Point3| -> (Scalar, Vec3) {
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
        let [p, q, r] = tri.map(|j| mesh.positions[j]);
        let (rp, rvec_p) = radial(p);
        let (rq, _) = radial(q);
        let (rr, _) = radial(r);
        let normal = (q - p).cross(r - p);
        if normal.length() < 1e-12 {
            continue;
        }
        let is_outer =
            (rp - outer).abs() < tol && (rq - outer).abs() < tol && (rr - outer).abs() < tol;
        let is_inner =
            (rp - inner).abs() < tol && (rq - inner).abs() < tol && (rr - inner).abs() < tol;
        if is_outer {
            assert!(
                normal.dot(rvec_p) > 0.0,
                "triangle {i} on the outer wall has an inward normal {normal:?}"
            );
            outer_checked += 1;
        } else if is_inner {
            assert!(
                normal.dot(rvec_p) < 0.0,
                "triangle {i} on the inner wall has an outward normal {normal:?}"
            );
            inner_checked += 1;
        }
    }
    assert!(outer_checked > 0, "no outer-wall triangles found to check");
    assert!(inner_checked > 0, "no inner-wall triangles found to check");
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

#[test]
fn a_tapered_revolution_crossing_or_touching_the_axis_is_refused() {
    // Both ends are judged, and together: a section that is clear of the
    // axis at each end but on opposite sides crosses it in between (#351).
    let tol = tol_for(1e-3);
    let at = |x: Scalar| {
        let mut rings = profile_rings(&rect(1.0, 2.0), 1e-3, tol).expect("rings");
        for p in &mut rings.outer {
            p.x += x;
        }
        rings
    };
    let turn = |a: &_, b: &_| {
        sweep::tapered_revolve(a, b, Point3::ZERO, Vec3::Y, core::f64::consts::PI, tol)
    };
    assert!(turn(&at(4.0), &at(5.0)).is_ok());
    assert!(matches!(
        turn(&at(0.0), &at(4.0)),
        Err(GeomError::InvalidInput(_))
    ));
    assert!(matches!(
        turn(&at(4.0), &at(0.0)),
        Err(GeomError::InvalidInput(_))
    ));
    assert!(matches!(
        turn(&at(4.0), &at(-4.0)),
        Err(GeomError::InvalidInput(_))
    ));
    // The end ring's left edge lies on the axis.
    assert!(matches!(
        turn(&at(4.0), &at(0.5)),
        Err(GeomError::UnsupportedInput { .. })
    ));
}

fn polygon(points: &[(Scalar, Scalar)]) -> Rings {
    Rings {
        outer: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
        holes: Vec::new(),
    }
}

/// A tapered turn whose start and end rings touch the axis (#364) is meshed
/// closed and holds the volume the swept section gives: a meridian plane at
/// angle `a` cuts the solid in the blended section, so the volume is
/// `angle * integral of M(t) dt` with `M` its first moment about the axis.
fn assert_tapered_revolution(
    start: &Rings,
    end: &Rings,
    moment_integral: Scalar,
    chord: Scalar,
    angle: Scalar,
) {
    let tol = tol_for(chord);
    let mesh = sweep::tapered_revolve(start, end, Point3::ZERO, Vec3::Y, angle, tol)
        .expect("a tapered revolution touching the axis at both ends");
    let health = axiolid_mesh::audit_mesh(&mesh, Tolerance::new(1e-9, 1e-9).expect("tolerance"));
    assert!(health.is_closed_two_manifold(), "{health:?}");
    // A vertex on the axis is one mesh vertex: no two coincide, and no
    // triangle collapsed onto a line survives.
    let mut seen = std::collections::HashSet::new();
    for p in &mesh.positions {
        assert!(
            seen.insert([p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]),
            "{p:?} twice"
        );
    }
    for t in mesh.indices.chunks(3) {
        let [a, b, c] = [t[0], t[1], t[2]].map(|i| mesh.positions[i as usize]);
        assert!(
            (b - a).cross(c - a).length() > 0.0,
            "degenerate triangle {t:?}"
        );
    }
    let (got, want) = (volume(&mesh, tol), angle * moment_integral);
    assert!(
        (got - want).abs() / want < 2.0 * chord,
        "tapered revolution {got} vs {want}"
    );
}

#[test]
fn a_tapered_revolution_of_a_cone_is_closed_and_has_the_swept_volume() {
    // The right triangle (0,0) (r,0) (0,2) grows from r = 1 to 2 with its
    // edge on the axis: first moment r^2 / 3, and the integral of
    // (1 + t)^2 / 3 over [0, 1] is 7 / 9.
    let start = polygon(&[(0.0, 0.0), (1.0, 0.0), (0.0, 2.0)]);
    let end = polygon(&[(0.0, 0.0), (2.0, 0.0), (0.0, 2.0)]);
    assert_tapered_revolution(&start, &end, 7.0 / 9.0, 1e-5, core::f64::consts::PI);
}

#[test]
fn a_tapered_revolution_of_a_double_cone_is_closed_and_has_the_swept_volume() {
    // (0,-1) (k,0) (0,1) with k from 1 to 2: both apexes stay on the axis.
    // First moment k^2 / 3, as for the cone; the integral is 7 / 9.
    let start = polygon(&[(0.0, -1.0), (1.0, 0.0), (0.0, 1.0)]);
    let end = polygon(&[(0.0, -1.0), (2.0, 0.0), (0.0, 1.0)]);
    assert_tapered_revolution(&start, &end, 7.0 / 9.0, 1e-5, core::f64::consts::PI);
}

#[test]
fn a_tapered_revolution_whose_axis_vertex_slides_along_the_axis_is_refused_by_name() {
    // The apexes are on the axis at both ends but move along it: their
    // stations are distinct points, and the caps' axis edges differ in
    // extent. Welding them would sweep the wrong solid.
    let start = polygon(&[(0.0, -1.0), (1.0, 0.0), (0.0, 1.0)]);
    let end = polygon(&[(0.0, -2.0), (2.0, 0.0), (0.0, 2.0)]);
    let error = sweep::tapered_revolve(&start, &end, Point3::ZERO, Vec3::Y, 1.0, tol_for(1e-3))
        .expect_err("a sliding vertex is not welded");
    assert!(
        matches!(error, GeomError::UnsupportedInput { .. }),
        "{error:?}"
    );
    assert!(format!("{error:?}").contains("slides along"), "{error:?}");
}

#[test]
fn a_tapered_revolution_touching_the_axis_at_one_end_only_is_refused_by_name() {
    let start = polygon(&[(0.0, 0.0), (1.0, 0.0), (1.0, 2.0), (0.0, 2.0), (0.0, 1.0)]);
    let end = polygon(&[(0.5, 0.0), (1.0, 0.0), (1.0, 2.0), (0.0, 2.0), (0.0, 1.0)]);
    let tol = tol_for(1e-3);
    for (a, b) in [(&start, &end), (&end, &start)] {
        let error = sweep::tapered_revolve(a, b, Point3::ZERO, Vec3::Y, 1.0, tol)
            .expect_err("a vertex leaves the axis through the blend");
        assert!(
            matches!(error, GeomError::UnsupportedInput { .. }),
            "{error:?}"
        );
        assert!(format!("{error:?}").contains("one end only"), "{error:?}");
    }
}

#[test]
fn a_partial_tapered_revolution_closes_an_isolated_touch_into_a_sector_point() {
    // The apex (0, 0) is the only vertex on the axis; a partial turn welds
    // it as `revolve` does, a full one pinches there. The triangle
    // (0,0) (r,-1) (r,1) has first moment 2 r^2 / 3, and r runs from 1 to 2.
    let start = polygon(&[(0.0, 0.0), (1.0, -1.0), (1.0, 1.0)]);
    let end = polygon(&[(0.0, 0.0), (2.0, -1.0), (2.0, 1.0)]);
    assert_tapered_revolution(&start, &end, 14.0 / 9.0, 1e-4, 1.0);
    let error = sweep::tapered_revolve(
        &start,
        &end,
        Point3::ZERO,
        Vec3::Y,
        core::f64::consts::TAU,
        tol_for(1e-3),
    )
    .expect_err("a full turn pinches at the apex");
    assert!(
        matches!(error, GeomError::UnsupportedInput { .. }),
        "{error:?}"
    );
}

#[test]
fn a_tapered_revolution_about_a_leaning_axis_is_judged_over_the_whole_blend() {
    // The axis leans out of the plane and meets it at the origin; its
    // projection N is the y axis.
    let lean = Vec3::new(0.0, 1.0, 0.5);
    let turn =
        |a: &Rings, b: &Rings| sweep::tapered_revolve(a, b, Point3::ZERO, lean, 1.0, tol_for(1e-3));
    let column = |lo: Scalar, hi: Scalar| polygon(&[(0.0, lo), (1.0, lo), (1.0, hi), (0.0, hi)]);
    // Each end is clear of the origin, but the edge on N slides through it.
    let error = turn(&column(1.0, 2.0), &column(-2.0, -1.0)).expect_err("slides through X");
    assert!(
        matches!(error, GeomError::UnsupportedInput { .. }),
        "{error:?}"
    );
    // Sliding along N clear of the origin is fine.
    turn(&column(1.0, 2.0), &column(3.0, 4.0)).expect("clear of X throughout");
    // A vertex on N at one end only touches it there: that end is judged.
    let off = polygon(&[(0.5, 3.0), (1.0, 3.0), (1.0, 4.0), (0.5, 4.0)]);
    turn(&column(1.0, 2.0), &off).expect("leaves N clear of X");
    assert!(matches!(
        turn(&column(-1.0, 1.0), &off),
        Err(GeomError::UnsupportedInput { .. })
    ));
    assert!(matches!(
        turn(&off, &column(-1.0, 1.0)),
        Err(GeomError::UnsupportedInput { .. })
    ));
    // A lone vertex on N slides through the origin between the ends.
    let rising = polygon(&[(0.0, 1.0), (1.0, 0.0), (1.0, 2.0)]);
    let falling = polygon(&[(0.0, -1.0), (1.0, -2.0), (1.0, 0.0)]);
    assert!(matches!(
        turn(&rising, &falling),
        Err(GeomError::UnsupportedInput { .. })
    ));
    // A lone vertex on N at the origin, at either end, is judged too.
    let apex = polygon(&[(0.0, 0.0), (1.0, -1.0), (1.0, 1.0)]);
    let far = polygon(&[(0.5, 3.0), (1.5, 2.0), (1.5, 4.0)]);
    for (a, b) in [(&apex, &far), (&far, &apex)] {
        assert!(matches!(
            turn(a, b),
            Err(GeomError::UnsupportedInput { .. })
        ));
    }
    // A vertex crossing N between the ends puts the rings on both sides.
    let right = polygon(&[(0.5, 1.0), (1.5, 1.0), (1.5, 2.0), (0.5, 2.0)]);
    let left = polygon(&[(-1.5, 1.0), (-0.5, 1.0), (-0.5, 2.0), (-1.5, 2.0)]);
    let error = turn(&right, &left).expect_err("crosses N");
    assert!(matches!(error, GeomError::InvalidInput(_)), "{error:?}");
}

#[test]
fn a_tapered_revolution_welds_apexes_within_the_tolerance_of_the_axis() {
    // The double cone of #364 with both apexes 1e-5 off the axis, inside
    // the 1e-3 tolerance (#367): they are projected onto it and welded.
    let start = polygon(&[(1e-5, -1.0), (1.0, 0.0), (1e-5, 1.0)]);
    let end = polygon(&[(1e-5, -1.0), (2.0, 0.0), (1e-5, 1.0)]);
    let tol = tol_for(1e-3);
    let mesh = sweep::tapered_revolve(&start, &end, Point3::ZERO, Vec3::Y, 1.0, tol)
        .expect("near-axis apexes are welded");
    let health = axiolid_mesh::audit_mesh(&mesh, Tolerance::new(1e-9, 1e-9).expect("tolerance"));
    assert!(health.is_closed_two_manifold(), "{health:?}");
    let near: Vec<_> = mesh
        .positions
        .iter()
        .filter(|p| p.x.hypot(p.z) < 1e-4)
        .collect();
    assert_eq!(near.len(), 2, "{near:?}");
    assert!(near.iter().all(|p| p.x == 0.0 && p.z == 0.0), "{near:?}");
    // The projection is charged to the chord budget: apexes 0.9 tol off the
    // axis leave a tenth of it for the turn, so it takes more steps.
    let apart = |x: Scalar| {
        let start = polygon(&[(x, -1.0), (1.0, 0.0), (x, 1.0)]);
        let end = polygon(&[(x, -1.0), (2.0, 0.0), (x, 1.0)]);
        sweep::tapered_revolve(&start, &end, Point3::ZERO, Vec3::Y, 1.0, tol)
            .expect("welded")
            .positions
            .len()
    };
    assert!(apart(0.9e-3) > 2 * apart(0.0));
}

#[test]
fn a_tapered_revolution_about_a_tilted_axis_welds_apexes_at_one_place_along_it() {
    // The axis through (0.1, 0.3) along (0.6, 0.8): the apexes lie 1e-9 and
    // 2e-9 off it at the same place along it, so their projections are one
    // point though not bitwise equal.
    let at = |s: Scalar, o: Scalar| Point2::new(0.1 + 0.6 * s - 0.8 * o, 0.3 + 0.8 * s + 0.6 * o);
    let ring = |o: Scalar, far: Scalar| Rings {
        outer: vec![at(1.0, o), at(2.0, o), at(1.5, far)],
        holes: Vec::new(),
    };
    let (start, end) = (ring(1e-9, 0.5), ring(2e-9, 0.9));
    let mesh = sweep::tapered_revolve(
        &start,
        &end,
        Point3::new(0.1, 0.3, 0.0),
        Vec3::new(0.6, 0.8, 0.0),
        1.0,
        tol_for(1e-3),
    )
    .expect("apexes at one place along the axis are welded");
    let health = axiolid_mesh::audit_mesh(&mesh, Tolerance::new(1e-9, 1e-9).expect("tolerance"));
    assert!(health.is_closed_two_manifold(), "{health:?}");
}

#[test]
fn a_tapered_revolution_with_one_end_in_the_band_is_refused() {
    let t = 1e-3;
    let start = polygon(&[(0.9 * t, -1.0), (1.0, 0.0), (0.9 * t, 1.0)]);
    let end = polygon(&[(1.1 * t, -1.0), (2.0, 0.0), (1.1 * t, 1.0)]);
    let error = sweep::tapered_revolve(&start, &end, Point3::ZERO, Vec3::Y, 1.0, tol_for(t))
        .expect_err("touching at one end only");
    assert!(format!("{error:?}").contains("one end only"), "{error:?}");
}

#[test]
fn a_tapered_ring_inside_the_band_is_refused() {
    let t = 1e-3;
    let flat = polygon(&[(0.2 * t, 0.0), (0.6 * t, 0.0), (0.5 * t, 1.0)]);
    let error = sweep::tapered_revolve(&flat, &flat, Point3::ZERO, Vec3::Y, 1.0, tol_for(t))
        .expect_err("a flattened ring");
    assert!(
        matches!(error, GeomError::UnsupportedInput { .. }),
        "{error:?}"
    );
}

#[test]
fn a_slide_within_the_tolerance_is_charged_to_the_budget() {
    // The apexes move 0.9 tol along the axis between the ends: they are
    // welded, and the slide takes that much of the chord budget.
    let t = 1e-3;
    let turn = |slide: Scalar| {
        let start = polygon(&[(0.0, -1.0), (1.0, 0.0), (0.0, 1.0)]);
        let end = polygon(&[(0.0, -1.0 + slide), (2.0, 0.0), (0.0, 1.0 - slide)]);
        sweep::tapered_revolve(&start, &end, Point3::ZERO, Vec3::Y, 1.0, tol_for(t))
            .expect("welded")
            .positions
            .len()
    };
    assert!(turn(0.9 * t) > 2 * turn(0.0));
}
