//! A difference whose tool stops a rounding error short of its host's
//! face, or reaches that far past it (#276).
//!
//! The wall is `3.65 x 0.25 x 3.67`, centred on `x` and `y`, standing on
//! `z = 0`; two door openings `1.01 x 2.26` stand on its floor face and run
//! along `+y` from `y0 = -0.125 + gap` to well past its `+y` face. With a
//! positive `gap` a skin that thin is left on the `-y` side; with a
//! negative one the doors reach past the `-y` face. The doors are placed by
//! an exact axis matrix, so `gap` is the distance the numbers state (up to
//! the rounding of `-0.125 + gap`), and the flush case `gap = 0` is exact.
//!
//! Subtracted in the wall's own frame and then placed at georeferenced
//! coordinates, as `ifc-geometry` lowers openings: a skin kept as two faces
//! `4.5e-15` apart crosses itself once its coordinates round to about
//! `1e-9`, and the volume kernel refuses it as self-intersecting.

use axiolid_brep::ExactBRep;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Scalar, Tolerance, Transform3, Vec3};
use axiolid_inspect::{enclosed_volume, VolumeInterval};
use axiolid_measure::exact_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::deviation::SNAPPED_OPERANDS;
use axiolid_mesh_compile::{
    BooleanReport, DeviationBound, DeviationReport, ReferenceExactCompiler, ReferenceMeshCompiler,
};
use axiolid_mesh_compile_contract::MeshClosure;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{Profile, RectangleProfile};

const LENGTH: Scalar = 3.65;
const THICKNESS: Scalar = 0.25;
const HEIGHT: Scalar = 3.67;
const DOOR_WIDTH: Scalar = 1.01;
const DOOR_HEIGHT: Scalar = 2.26;
const DOORS: [Scalar; 2] = [-0.9, 0.8];

/// The wall less both doors cut through it.
fn net() -> Scalar {
    LENGTH * THICKNESS * HEIGHT - DOORS.len() as Scalar * DOOR_WIDTH * DOOR_HEIGHT * THICKNESS
}

/// The skin both doors leave when they stop `gap` short.
fn skin(gap: Scalar) -> Scalar {
    DOORS.len() as Scalar * DOOR_WIDTH * DOOR_HEIGHT * gap
}

/// More than the wall's surface area, which bounds how far moving every
/// point of it by `d` changes its volume (`area * d`).
const AREA: Scalar = 40.0;

fn rect(x: Scalar, y: Scalar) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// The wall less its doors, both `gap` short, in its own frame, under
/// `placement` when one is given.
fn wall(gap: Scalar, placement: Option<Transform3>) -> (GeometryGraph, NodeId) {
    wall_with([gap, gap], placement)
}

/// The wall less its doors, each its own `gap` short.
fn wall_with(gaps: [Scalar; 2], placement: Option<Transform3>) -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let mut push = |node| b.push(node).expect("a valid node");
    let profile = push(GeometryNode::Profile(rect(LENGTH, THICKNESS)));
    let mut body = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: HEIGHT,
    }));
    for (cx, gap) in DOORS.into_iter().zip(gaps) {
        let y0 = -THICKNESS / 2.0 + gap;
        let profile = push(GeometryNode::Profile(rect(DOOR_WIDTH, DOOR_HEIGHT)));
        let door = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: THICKNESS + 0.2,
        }));
        // Profile x along -x, profile y up, extruded along +y: entries 0
        // and +-1 only, so every corner lands exactly where it is stated.
        let across = Transform3::from_cols(
            Vec3::NEG_X,
            Vec3::Z,
            Vec3::Y,
            Vec3::new(cx, y0, DOOR_HEIGHT / 2.0),
        );
        let door = push(GeometryNode::Instance(Instance {
            source: door,
            transform: across,
        }));
        body = push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: body,
            right: door,
            operator: BooleanOperator::Difference,
        }));
    }
    if let Some(transform) = placement {
        body = push(GeometryNode::Instance(Instance {
            source: body,
            transform,
        }));
    }
    (b.finish(vec![body]).expect("a valid graph"), body)
}

/// Georeferenced: the site of the consumer model, turned by 2.3 degrees.
fn site() -> Transform3 {
    Transform3::from_translation(Vec3::new(6.0e5, 5.6e6, 0.0))
        * Transform3::from_rotation_z(2.3_f64.to_radians())
}

/// Far out, and only moved.
fn far() -> Transform3 {
    Transform3::from_translation(Vec3::new(4.0e5, 4.0e5, 0.0))
}

/// The caller's tolerance, with a chord budget of its own so a zero
/// tolerance still flattens.
fn options(tolerance: Tolerance) -> ExecutionOptions {
    ExecutionOptions::new(tolerance)
        .with_chord_error(1e-3)
        .expect("a valid budget")
}

fn mesh(
    gap: Scalar,
    placement: Option<Transform3>,
    tolerance: Tolerance,
) -> (TriMesh, DeviationReport) {
    let (graph, root) = wall(gap, placement);
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options(tolerance))
        .expect("the wall compiles");
    assert_eq!(outcome.closure, MeshClosure::Solid);
    (outcome.mesh, report)
}

fn exact(
    gap: Scalar,
    placement: Option<Transform3>,
    tolerance: Tolerance,
) -> Result<(ExactBRep, BooleanReport), GeomError> {
    let (graph, root) = wall(gap, placement);
    ReferenceExactCompiler::new().compile_exact_with_report(&graph, root, &options(tolerance))
}

/// The largest snap the report names, `None` when nothing was snapped.
fn snapped(report: &DeviationReport) -> Option<Scalar> {
    let snaps: Vec<Scalar> = report
        .contributions
        .iter()
        .filter(|c| c.detail == SNAPPED_OPERANDS)
        .map(|c| match c.bound {
            DeviationBound::Certified(moved) => moved,
            other => panic!("a snap is certified: {other:?}"),
        })
        .collect();
    assert!(snaps.len() <= 1, "{report:?}");
    snaps.first().copied()
}

/// `expected` within `slack` of the certified interval.
fn assert_volume(volume: VolumeInterval, expected: Scalar, slack: Scalar, what: &str) {
    assert!(
        volume.lower - slack <= expected && expected <= volume.upper + slack,
        "{what}: [{}, {}] does not hold {expected} within {slack:e}",
        volume.lower,
        volume.upper
    );
}

/// The rounding of a coordinate placed `size` from the origin.
fn ulp(size: Scalar) -> Scalar {
    size * Scalar::EPSILON
}

fn triangles(mesh: &TriMesh) -> usize {
    mesh.indices.len() / 3
}

fn exact_volume(brep: &ExactBRep) -> Scalar {
    exact_properties(brep, Tolerance::METRE)
        .expect("measurable")
        .signed_volume
}

fn closed(brep: &ExactBRep) {
    let topology = axiolid_topology::audit_brep(brep.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
}

/// Short of the face (a skin) and past it (an overlap), from the exported
/// residue to a thousandth of the tolerance.
const THIN: [Scalar; 6] = [4.5e-15, -4.5e-15, 1e-12, -1e-12, 1e-9, -1e-9];

#[test]
fn a_skin_thinner_than_the_tolerance_is_snapped_away_and_reported() {
    let tolerance = Tolerance::METRE;
    let (flush, flush_report) = mesh(0.0, None, tolerance);
    assert_eq!(snapped(&flush_report), None, "flush needs no snap");
    let flush_volume = enclosed_volume(&flush).expect("a valid solid");
    assert!(flush_volume.contains(net()));
    for gap in THIN {
        let (mesh, report) = mesh(gap, None, tolerance);
        let moved = snapped(&report).unwrap_or_else(|| panic!("gap {gap}: {report:?}"));
        // The doors' ends moved onto the wall's face: by the gap, never
        // beyond the tolerance.
        assert!(
            moved >= 0.5 * gap.abs() && moved <= tolerance.linear(),
            "{gap}: {moved}"
        );
        assert_eq!(triangles(&mesh), triangles(&flush), "gap {gap}");
        let volume = enclosed_volume(&mesh).expect("a valid solid");
        assert_volume(volume, net(), AREA * moved, &format!("gap {gap}"));
        // The report's bound covers the snap.
        assert!(report.bound.is_some_and(|b| b >= moved), "{report:?}");
    }
}

#[test]
fn an_inner_booleans_snap_stays_reported_under_the_outer_one() {
    // Only the first door is short; the second, cut after it, is flush and
    // needs no snap. The outer difference is the one measured against its
    // exact result, and the first door's snap must survive that.
    let (graph, root) = wall_with([4.5e-15, 0.0], Some(site()));
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options(Tolerance::METRE))
        .expect("the wall compiles");
    let moved = snapped(&report).unwrap_or_else(|| panic!("{report:?}"));
    assert!((2e-15..=1e-14).contains(&moved), "{moved}");
    let volume = enclosed_volume(&outcome.mesh).expect("a valid placed solid");
    assert_volume(volume, net(), AREA * (moved + 4.0 * ulp(5.7e6)), "placed");
}

#[test]
fn a_snapped_wall_placed_at_georeferenced_coordinates_keeps_a_certified_volume() {
    let tolerance = Tolerance::METRE;
    for gap in THIN {
        for (name, placement) in [("site", site()), ("far", far())] {
            let (mesh, report) = mesh(gap, Some(placement), tolerance);
            let moved = snapped(&report).expect("snapped");
            let volume =
                enclosed_volume(&mesh).unwrap_or_else(|e| panic!("gap {gap}, {name}: {e:?}"));
            // Each placed coordinate rounds by about one unit in the last
            // place of the placement's size.
            let slack = AREA * (moved + 4.0 * ulp(5.7e6));
            assert_volume(volume, net(), slack, &format!("gap {gap}, {name}"));
        }
    }
}

#[test]
fn a_thin_skin_has_an_exact_boundary_whose_report_names_the_reading() {
    let tolerance = Tolerance::METRE;
    let (flush, flush_report) = exact(0.0, None, tolerance).expect("flush");
    assert!(flush_report.is_exact(), "{flush_report:?}");
    closed(&flush);
    for gap in THIN {
        let (brep, report) = exact(gap, None, tolerance).unwrap_or_else(|e| panic!("{gap}: {e}"));
        closed(&brep);
        // Read within tolerance, and reported: never silently.
        assert!(!report.is_exact(), "gap {gap}: {report:?}");
        let linear = report
            .decisions()
            .iter()
            .map(|d| d.linear)
            .fold(0.0, Scalar::max);
        assert!(
            linear >= 0.5 * gap.abs() && linear <= tolerance.linear(),
            "{report:?}"
        );
        let measured = exact_volume(&brep);
        assert!(
            (measured - net()).abs() <= AREA * linear + 1e-12,
            "gap {gap}: {measured} for {}",
            net()
        );
        for placement in [site(), far()] {
            let (placed, placed_report) =
                exact(gap, Some(placement), tolerance).unwrap_or_else(|e| panic!("{gap}: {e}"));
            closed(&placed);
            assert_eq!(placed_report.decisions(), report.decisions());
        }
    }
}

#[test]
fn a_skin_thicker_than_the_tolerance_is_kept() {
    // Ten tolerances, and a micrometre at a tenth of one.
    for (gap, tolerance) in [
        (1e-5, Tolerance::METRE),
        (1e-6, Tolerance::new(1e-7, 1e-9).expect("valid")),
    ] {
        let (kept, report) = mesh(gap, None, tolerance);
        assert_eq!(snapped(&report), None, "gap {gap}: {report:?}");
        let volume = enclosed_volume(&kept).expect("a valid solid");
        assert_volume(volume, net() + skin(gap), 0.0, &format!("gap {gap}"));
        assert!(!volume.contains(net()), "gap {gap}: the skin is kept");
        let (placed, _) = mesh(gap, Some(site()), tolerance);
        let volume = enclosed_volume(&placed).expect("a valid placed solid");
        assert_volume(
            volume,
            net() + skin(gap),
            AREA * 4.0 * ulp(5.7e6),
            &format!("gap {gap}, placed"),
        );
        let (brep, _) = exact(gap, None, tolerance).expect("exact");
        closed(&brep);
        let measured = exact_volume(&brep);
        assert!(
            (measured - net() - skin(gap)).abs() <= 1e-12,
            "gap {gap}: {measured}"
        );
    }
}

#[test]
fn at_zero_tolerance_nothing_is_snapped() {
    for gap in [4.5e-15, -4.5e-15, 1e-9, 1e-6] {
        let (kept, report) = mesh(gap, None, Tolerance::ZERO);
        assert_eq!(snapped(&report), None, "gap {gap}: {report:?}");
        if gap >= 1e-9 {
            let volume = enclosed_volume(&kept).expect("a valid solid");
            assert_volume(volume, net() + skin(gap), 0.0, "zero tolerance");
            assert!(!volume.contains(net()), "gap {gap}: the skin is kept");
        }
        // The exact boolean reads nothing within tolerance: its result is
        // the exact difference of the operands as given, skin and all, or
        // a refusal by name where a piece is too thin to classify exactly.
        match exact(gap, None, Tolerance::ZERO) {
            Ok((brep, report)) => {
                assert!(report.is_exact(), "gap {gap}: {report:?}");
                closed(&brep);
                let measured = exact_volume(&brep);
                assert!(
                    (measured - net() - skin(gap.max(0.0))).abs() <= 1e-12,
                    "gap {gap}: {measured}"
                );
            }
            Err(GeomError::UnsupportedInput { .. }) => assert!(gap < 1e-6, "gap {gap}"),
            Err(other) => panic!("gap {gap}: {other}"),
        }
    }
}
