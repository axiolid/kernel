//! Booleans cut by tilted elliptical-cylinder walls (#287).
//!
//! An exact extrusion of a profile with arcs along a leaning direction
//! sweeps each arc into an oblique circular cylinder, the
//! `EllipticalCylinder` whose axis is the sweep direction (#280). A plane
//! cuts it in an ellipse (a circle where it is exactly one), or in rulings
//! where it is parallel to the axis; two such walls swept along one
//! direction meet in rulings. Each is built in closed form, so a slab less
//! an oblique round shaft or rounded opening is exact: a horizontal section
//! of an oblique prism is its profile, moved, so the slab loses the
//! profile's area times its depth.

use std::f64::consts::PI;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::{boolean, boolean_with_report, section_edges, split_face, PieceSource};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{BooleanOperator, Tolerance, Transform3, Vec3};
use axiolid_curve::{Curve2, Curve3};
use axiolid_measure::exact_properties;
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};
use axiolid_surface::Surface;

fn tol() -> Tolerance {
    Tolerance::METRE
}

/// The slab `[-2, 2] x [-1.5, 1.5] x [0, D]`.
const D: f64 = 0.25;

fn rect(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

fn slab() -> ExactBRep {
    extrude_profile_exact(&rect(4.0, 3.0), Vec3::Z, D, tol()).expect("a slab")
}

fn circle(radius: f64) -> Profile {
    Profile::Circle(CircleProfile {
        radius,
        thickness: None,
    })
}

/// `0.8 x 0.5` with corners of radius `0.15`.
fn rounded() -> Profile {
    Profile::Rectangle(RectangleProfile {
        x: 0.8,
        y: 0.5,
        thickness: None,
        outer_radius: Some(0.15),
        inner_radius: None,
    })
}

const ROUNDED_AREA: f64 = 0.4 - (4.0 - PI) * 0.0225;

/// `profile` swept along `lean` (upward, or downward: a reflected prism)
/// far enough to pass through the slab whole, its profile plane at
/// `(x, y, z)`.
fn oblique(profile: &Profile, lean: Vec3, at: Vec3) -> ExactBRep {
    let depth = 1.0 * lean.length() / lean.z.abs();
    let prism = extrude_profile_exact(profile, lean, depth, tol())
        .unwrap_or_else(|e| panic!("an oblique prism: {e}"));
    assert!(prism
        .surfaces()
        .iter()
        .any(|s| matches!(s, Surface::EllipticalCylinder(_))));
    prism
        .transformed(&Transform3::from_translation(at))
        .expect("rigid")
}

fn run(a: &ExactBRep, b: &ExactBRep, op: BooleanOperator, what: &str) -> ExactBRep {
    let result = boolean(a, b, op, tol()).unwrap_or_else(|e| panic!("{what} {op:?}: {e}"));
    let health = geometric_audit(&result, tol());
    assert!(
        health.is_consistent(),
        "{what} {op:?}: {:?}",
        health.defects()
    );
    let topology = axiolid_topology::audit_brep(result.topology());
    assert!(topology.is_closed_manifold(), "{what} {op:?}: {topology:?}");
    result
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, tol())
        .expect("measurable")
        .signed_volume
}

fn close(what: &str, got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-9 * want.abs().max(1.0),
        "{what}: expected {want}, got {got}"
    );
}

/// Every section edge of the pair is a line or a conic in closed form, and
/// at least one lies on an elliptical-cylinder wall as `kind`.
fn closed_form_sections(a: &ExactBRep, b: &ExactBRep, kind: fn(&Curve3) -> bool, what: &str) {
    let edges = section_edges(a, b, tol()).unwrap_or_else(|e| panic!("{what}: {e}"));
    assert!(!edges.is_empty(), "{what}: no section");
    for edge in &edges {
        assert!(
            matches!(
                edge.curve,
                Curve3::Line(_) | Curve3::Circle(_) | Curve3::Ellipse(_)
            ),
            "{what}: a section is not closed form: {:?}",
            edge.curve
        );
    }
    let walled = |s: &Surface| matches!(s, Surface::EllipticalCylinder(_));
    assert!(
        edges
            .iter()
            .any(|e| (walled(&e.other_a) || walled(&e.other_b)) && kind(&e.curve)),
        "{what}: no expected section on the wall"
    );
}

fn conic(curve: &Curve3) -> bool {
    matches!(curve, Curve3::Circle(_) | Curve3::Ellipse(_))
}

fn line(curve: &Curve3) -> bool {
    matches!(curve, Curve3::Line(_))
}

fn leans() -> [Vec3; 4] {
    [
        Vec3::new(0.3, -0.2, 1.0),
        Vec3::new(-0.4, 0.25, 1.0),
        // Downward: the reflection of the prism swept up.
        Vec3::new(0.3, -0.2, -1.0),
        Vec3::new(-0.4, 0.25, -1.0),
    ]
}

/// Where a prism along `lean` starts so that it passes through the slab.
fn start(lean: Vec3, x: f64, y: f64) -> Vec3 {
    let z = if lean.z > 0.0 { -0.3 } else { D + 0.3 };
    Vec3::new(x, y, z)
}

#[test]
fn a_slab_minus_an_oblique_round_shaft_loses_the_disk_times_its_depth() {
    let r = 0.3;
    for lean in leans() {
        let what = format!("shaft along {lean:?}");
        let shaft = oblique(&circle(r), lean, start(lean, -1.0, 0.4));
        closed_form_sections(&slab(), &shaft, conic, &what);
        let (cut, report) =
            boolean_with_report(&slab(), &shaft, BooleanOperator::Difference, tol())
                .unwrap_or_else(|e| panic!("{what}: {e}"));
        assert!(report.is_exact(), "{what}: {report:?}");
        let health = geometric_audit(&cut, tol());
        assert!(health.is_consistent(), "{what}: {:?}", health.defects());
        close(&what, volume(&cut), 4.0 * 3.0 * D - PI * r * r * D);
        // Exact decisions only: the same at zero tolerance.
        let zero = boolean(
            &slab(),
            &shaft,
            BooleanOperator::Difference,
            Tolerance::ZERO,
        )
        .unwrap_or_else(|e| panic!("{what} at zero tolerance: {e}"));
        close(&format!("{what} zero"), volume(&zero), volume(&cut));
    }
}

#[test]
fn a_slab_minus_an_oblique_rounded_opening_loses_its_area_times_its_depth() {
    for lean in leans() {
        let what = format!("rounded opening along {lean:?}");
        let opening = oblique(&rounded(), lean, start(lean, 0.9, -0.5));
        closed_form_sections(&slab(), &opening, conic, &what);
        let (cut, report) =
            boolean_with_report(&slab(), &opening, BooleanOperator::Difference, tol())
                .unwrap_or_else(|e| panic!("{what}: {e}"));
        assert!(report.is_exact(), "{what}: {report:?}");
        close(&what, volume(&cut), 4.0 * 3.0 * D - ROUNDED_AREA * D);
    }
}

#[test]
fn a_slab_with_both_oblique_openings_under_a_general_placement() {
    let p = Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7);
    let place = |b: ExactBRep| b.transformed(&p).expect("rigid");
    let r = 0.3;
    let slab = place(slab());
    let up = Vec3::new(0.3, -0.2, 1.0);
    let down = Vec3::new(-0.4, 0.25, -1.0);
    let shaft = place(oblique(&circle(r), up, start(up, -1.0, 0.4)));
    let opening = place(oblique(&rounded(), down, start(down, 0.9, -0.5)));
    let once = run(&slab, &shaft, BooleanOperator::Difference, "placed shaft");
    let twice = run(
        &once,
        &opening,
        BooleanOperator::Difference,
        "placed opening",
    );
    close(
        "placed slab",
        volume(&twice),
        4.0 * 3.0 * D - (PI * r * r + ROUNDED_AREA) * D,
    );
}

#[test]
fn a_blind_oblique_hole_ends_on_its_cap() {
    // Swept from `z = 0.1` inside the slab up through its top: the cap
    // circle is a hole's floor, and the top face cuts the wall.
    let r = 0.3;
    for lean in [Vec3::new(0.3, -0.2, 1.0), Vec3::new(-0.2, 0.5, 1.0)] {
        let what = format!("blind hole along {lean:?}");
        let hole = oblique(&circle(r), lean, Vec3::new(0.2, -0.1, 0.1));
        let cut = run(&slab(), &hole, BooleanOperator::Difference, &what);
        close(&what, volume(&cut), 4.0 * 3.0 * D - PI * r * r * (D - 0.1));
        let common = run(&slab(), &hole, BooleanOperator::Intersection, &what);
        close(
            &format!("{what} common"),
            volume(&common),
            PI * r * r * (D - 0.1),
        );
        let union = run(&slab(), &hole, BooleanOperator::Union, &what);
        close(
            &format!("{what} union"),
            volume(&union) + volume(&common),
            volume(&slab()) + volume(&hole),
        );
    }
}

#[test]
fn a_plane_along_the_lean_cuts_the_wall_in_rulings() {
    // Leaning along `x`, the wall's axis lies in every plane `y = c`: a
    // block whose face is `y = c` cuts it in two rulings, and each
    // horizontal section keeps the disk's segment beyond `y = c`.
    let r: f64 = 0.3;
    let c: f64 = 0.1;
    let block = extrude_profile_exact(&rect(4.0, 1.5 - c), Vec3::Z, D, tol())
        .expect("a block")
        .transformed(&Transform3::from_translation(Vec3::new(
            0.0,
            0.5 * (1.5 + c),
            0.0,
        )))
        .expect("rigid");
    let segment = r * r * (c / r).acos() - c * (r * r - c * c).sqrt();
    for lean in [Vec3::new(0.5, 0.0, 1.0), Vec3::new(-0.5, 0.0, -1.0)] {
        let what = format!("rulings along {lean:?}");
        let shaft = oblique(&circle(r), lean, start(lean, 0.0, 0.0));
        closed_form_sections(&block, &shaft, line, &what);
        let (common, report) =
            boolean_with_report(&block, &shaft, BooleanOperator::Intersection, tol())
                .unwrap_or_else(|e| panic!("{what}: {e}"));
        assert!(report.is_exact(), "{what}: {report:?}");
        close(&what, volume(&common), segment * D);
        let cut = run(&block, &shaft, BooleanOperator::Difference, &what);
        close(
            &format!("{what} difference"),
            volume(&cut),
            volume(&block) - segment * D,
        );
    }
}

#[test]
fn an_upright_face_across_the_lean_cuts_the_wall_in_ellipses() {
    // A block `x >= x0` over the slab's height keeps, at each height, the
    // disk's segment beyond `x0`. The disk's centre moves along `x` at
    // `k = lean.x / lean.z` per unit height, so the segment's depth `d`
    // is linear in `z`, and the volume is `[F(d)] / (dd/dz)` with
    // `F' = r^2 acos(d / r) - d sqrt(r^2 - d^2)`. The face's normal has a
    // part along the wall's frame `x`, which a level face's has not.
    let r: f64 = 0.3;
    let x0 = 0.1;
    let block = extrude_profile_exact(&rect(2.0 - x0, 3.0), Vec3::Z, D, tol())
        .expect("a block")
        .transformed(&Transform3::from_translation(Vec3::new(
            0.5 * (2.0 + x0),
            0.0,
            0.0,
        )))
        .expect("rigid");
    let f = |d: f64| {
        r * r * (d * (d / r).acos() - (r * r - d * d).sqrt()) + (r * r - d * d).powf(1.5) / 3.0
    };
    for lean in [Vec3::new(0.3, -0.2, 1.0), Vec3::new(0.3, -0.2, -1.0)] {
        let what = format!("upright face across {lean:?}");
        let at = start(lean, 0.0, 0.0);
        let shaft = oblique(&circle(r), lean, at);
        closed_form_sections(&block, &shaft, conic, &what);
        let k = lean.x / lean.z;
        let depth = |z: f64| x0 - k * (z - at.z);
        let want = (f(depth(D)) - f(depth(0.0))) / -k;
        let (common, report) =
            boolean_with_report(&block, &shaft, BooleanOperator::Intersection, tol())
                .unwrap_or_else(|e| panic!("{what}: {e}"));
        assert!(report.is_exact(), "{what}: {report:?}");
        close(&what, volume(&common), want);
        let cut = run(&block, &shaft, BooleanOperator::Difference, &what);
        close(
            &format!("{what} difference"),
            volume(&cut),
            volume(&block) - want,
        );
    }
}

#[test]
fn two_shafts_along_one_lean_meet_in_rulings() {
    // Disks of radius `r` with centres `s` apart, swept along one
    // direction over one height: every horizontal section is their lens.
    let r: f64 = 0.3;
    let s: f64 = 0.2;
    let lens = 2.0 * r * r * (s / (2.0 * r)).acos() - 0.5 * s * (4.0 * r * r - s * s).sqrt();
    for lean in [Vec3::new(0.4, 0.0, 1.0), Vec3::new(0.4, 0.0, -1.0)] {
        let what = format!("parallel shafts along {lean:?}");
        let height = 1.0;
        let a = oblique(&circle(r), lean, Vec3::ZERO);
        let b = oblique(&circle(r), lean, Vec3::new(0.0, s, 0.0));
        closed_form_sections(&a, &b, line, &what);
        let common = run(&a, &b, BooleanOperator::Intersection, &what);
        close(&what, volume(&common), lens * height);
        let union = run(&a, &b, BooleanOperator::Union, &what);
        close(
            &format!("{what} union"),
            volume(&union),
            2.0 * PI * r * r * height - lens * height,
        );
        let cut = run(&a, &b, BooleanOperator::Difference, &what);
        close(
            &format!("{what} difference"),
            volume(&cut),
            PI * r * r * height - lens * height,
        );
    }
}

/// The section pieces split into the shaft's wall, with their pcurves.
fn wall_pieces(slab: &ExactBRep, shaft: &ExactBRep) -> Vec<(Curve3, Curve2)> {
    let edges = section_edges(slab, shaft, tol()).expect("sections");
    let cuts: Vec<_> = edges.iter().flat_map(|e| [e.start, e.end]).collect();
    let mut out = Vec::new();
    for index in 0..shaft.topology().faces().len() {
        let face = shaft.topology().face_id_at(index).expect("a face");
        let record = &shaft.topology().faces()[index];
        let surface = &shaft.surfaces()[record.surface.expect("a surface").index()];
        if !matches!(surface, Surface::EllipticalCylinder(_)) {
            continue;
        }
        let mine: Vec<_> = edges.iter().filter(|e| e.face_b == face).cloned().collect();
        let regions = split_face(shaft, face, &mine, false, &cuts, tol()).expect("split");
        for region in regions {
            for piece in region.outer.iter().chain(region.holes.iter().flatten()) {
                if matches!(piece.source, PieceSource::Section(_)) {
                    out.push((piece.curve.clone(), piece.pcurve.clone()));
                }
            }
        }
    }
    out
}

#[test]
fn a_cut_wall_carries_closed_form_pcurves() {
    // A slab face's ellipse is a sinusoid in the wall's angle, a ruling a
    // vertical line: never a traced implicit pcurve.
    let r = 0.3;
    for lean in leans() {
        let shaft = oblique(&circle(r), lean, start(lean, -1.0, 0.4));
        let pieces = wall_pieces(&slab(), &shaft);
        assert!(!pieces.is_empty(), "{lean:?}");
        for (curve, pcurve) in &pieces {
            assert!(conic(curve), "{lean:?}: {curve:?}");
            assert!(
                matches!(pcurve, Curve2::Sinusoid(_)),
                "{lean:?}: {pcurve:?}"
            );
        }
    }
    let c = 0.1;
    let block = extrude_profile_exact(&rect(4.0, 1.5 - c), Vec3::Z, D, tol())
        .expect("a block")
        .transformed(&Transform3::from_translation(Vec3::new(
            0.0,
            0.5 * (1.5 + c),
            0.0,
        )))
        .expect("rigid");
    for lean in [Vec3::new(0.5, 0.0, 1.0), Vec3::new(-0.5, 0.0, -1.0)] {
        let shaft = oblique(&circle(r), lean, start(lean, 0.0, 0.0));
        let pieces = wall_pieces(&block, &shaft);
        assert!(
            pieces
                .iter()
                .any(|(curve, pcurve)| line(curve) && matches!(pcurve, Curve2::Line(_))),
            "{lean:?}: {pieces:?}"
        );
        assert!(
            pieces
                .iter()
                .all(|(_, pcurve)| matches!(pcurve, Curve2::Line(_) | Curve2::Sinusoid(_))),
            "{lean:?}: {pieces:?}"
        );
    }
}

#[test]
fn a_flush_oblique_shaft_runs_along_the_slab_faces() {
    // Swept from the slab's bottom face exactly to its top: each cap lies
    // in a slab face, and the face meets the wall along the cap's rim.
    let r = 0.3;
    for lean in [Vec3::new(0.3, -0.2, 1.0), Vec3::new(0.5, 0.0, -1.0)] {
        let what = format!("flush shaft along {lean:?}");
        let base = if lean.z > 0.0 { 0.0 } else { D };
        let depth = D * lean.length() / lean.z.abs();
        let shaft = extrude_profile_exact(&circle(r), lean, depth, tol())
            .expect("a prism")
            .transformed(&Transform3::from_translation(Vec3::new(-0.5, 0.2, base)))
            .expect("rigid");
        let cut = run(&slab(), &shaft, BooleanOperator::Difference, &what);
        close(&what, volume(&cut), 4.0 * 3.0 * D - PI * r * r * D);
    }
}
