//! A trimmed conic directrix runs the way its sense says, across the seam.
//!
//! ISO 10303-42 `trimmed_curve`: the curve runs from `trim_1` to `trim_2`,
//! in the basis's own direction when `sense_agreement` is true and against
//! it when false. On a periodic basis that may cross the seam: a circle
//! trimmed `315 deg -> 45 deg` with the sense is the 90 degree arc through
//! 0, not the 270 degree arc through 180. Sorting the two trims picks the
//! wrong one. Standalone that is silently wrong geometry; inside a
//! composite the ends no longer meet and the directrix is refused as
//! gappy (axiolid/kernel#168; 1,494 bent rebars in a real Revit model).

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame3, Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::{Circle3, Curve3, Ellipse3, Line3};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    CurveRelation, CurveSegment, GeometryGraphBuilder, GeometryNode, NodeId, SolidOperation,
    Transition, TrimSelector, TrimmingPreference,
};

const RADIUS: Scalar = 1.0;
const TUBE: Scalar = 0.05;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

fn frame() -> Frame3 {
    Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

fn circle() -> Curve3 {
    Curve3::Circle(Circle3 {
        frame: frame(),
        radius: RADIUS,
    })
}

fn trimmed(
    b: &mut GeometryGraphBuilder,
    basis: NodeId,
    t1: Scalar,
    t2: Scalar,
    sense: bool,
) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
        basis,
        start: vec![TrimSelector::Parameter(t1)],
        end: vec![TrimSelector::Parameter(t2)],
        sense_agreement: sense,
        preference: TrimmingPreference::Parameter,
    }))
    .unwrap()
}

fn disk(
    mut b: GeometryGraphBuilder,
    directrix: NodeId,
    range: Option<(Scalar, Scalar)>,
) -> axiolid_contracts::GeomResult<TriMesh> {
    let sweep = b
        .push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
            directrix,
            radius: TUBE,
            inner_radius: None,
            parameter_range: range,
            fillet_radius: None,
        }))
        .unwrap();
    let graph = b.finish(vec![sweep]).unwrap();
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(&graph, sweep, &options())
}

/// Area of the inscribed n-gon the compiler meshes the disk as (the same
/// chord rule as `axiolid-construct`), so volumes compare path length only.
fn section_area() -> Scalar {
    let chord = Tolerance::MILLIMETRE.linear();
    let per = 2.0 * (1.0 - (chord / TUBE).min(1.0)).clamp(-1.0, 1.0).acos();
    let n = (core::f64::consts::TAU / per).ceil();
    0.5 * n * TUBE * TUBE * (core::f64::consts::TAU / n).sin()
}

fn volume(mesh: &TriMesh) -> Scalar {
    volume_properties(mesh, Tolerance::MILLIMETRE)
        .expect("a swept disk is closed")
        .signed_volume
}

fn assert_length(mesh: &TriMesh, length: Scalar, what: &str) {
    let want = section_area() * length;
    let got = volume(mesh);
    assert!(
        (got - want).abs() / want < 0.01,
        "{what}: volume {got} is not section x {length} = {want}"
    );
}

/// Direction of the mean vertex: the mid-angle of an arc about the origin.
fn mean_angle_deg(mesh: &TriMesh) -> Scalar {
    let n = mesh.positions.len() as Scalar;
    let (x, y) = mesh
        .positions
        .iter()
        .fold((0.0, 0.0), |(x, y), p| (x + p.x / n, y + p.y / n));
    y.atan2(x).to_degrees().rem_euclid(360.0)
}

fn angle_gap(a: Scalar, b: Scalar) -> Scalar {
    let d = (a - b).rem_euclid(360.0);
    d.min(360.0 - d)
}

#[test]
fn a_trimmed_circle_runs_from_trim1_the_way_its_sense_says() {
    // (trim1, trim2, sense, swept degrees, mid-angle)
    let cases: &[(Scalar, Scalar, bool, Scalar, Scalar)] = &[
        (315.0, 45.0, true, 90.0, 0.0),
        (45.0, 315.0, false, 90.0, 0.0),
        (270.0, 45.0, true, 135.0, 337.5),
        (45.0, 270.0, false, 135.0, 337.5),
        // As Revit writes a quarter bend: just past the seam.
        (270.0, 360.00000000000034, true, 90.0, 315.0),
        (270.0, 15.0, true, 105.0, 322.5),
        // The long way round is what these trims SAY.
        (45.0, 315.0, true, 270.0, 180.0),
        (315.0, 45.0, false, 270.0, 180.0),
        // Not crossing the seam: unchanged.
        (90.0, 180.0, true, 90.0, 135.0),
        (180.0, 90.0, false, 90.0, 135.0),
    ];
    for &(t1, t2, sense, span, mid) in cases {
        let what = format!("({t1}, {t2}, sense {sense})");
        let mut b = GeometryGraphBuilder::new();
        let basis = b.push(GeometryNode::Curve3(circle())).unwrap();
        let arc = trimmed(&mut b, basis, t1.to_radians(), t2.to_radians(), sense);
        let mesh = disk(b, arc, None).unwrap_or_else(|e| panic!("{what}: {e}"));
        assert_length(&mesh, RADIUS * span.to_radians(), &what);
        let got = mean_angle_deg(&mesh);
        assert!(
            angle_gap(got, mid) < 3.0,
            "{what}: arc centred at {got} deg, expected {mid} deg"
        );
    }
}

/// A full turn written with rounding past the seam stays a full turn.
///
/// `0 -> 360.0000000000003` degrees travels just over one period. Reduced
/// modulo the period it would be a sliver of 3e-13 degrees and the ring
/// would vanish; within slack it is the whole circle.
#[test]
fn a_full_turn_rounded_past_the_seam_stays_a_full_turn() {
    let turns: [(Scalar, Scalar); 2] = [(0.0, 360.00000000000034), (90.0, 450.00000000000006)];
    for (t1, t2) in turns {
        let what = format!("({t1}, {t2})");
        let mut b = GeometryGraphBuilder::new();
        let basis = b.push(GeometryNode::Curve3(circle())).unwrap();
        let ring = trimmed(&mut b, basis, t1.to_radians(), t2.to_radians(), true);
        let mesh = disk(b, ring, None).unwrap_or_else(|e| panic!("{what}: {e}"));
        assert_length(&mesh, RADIUS * core::f64::consts::TAU, &what);
    }
}

/// The same rule on an ellipse, the other periodic conic.
#[test]
fn a_trimmed_ellipse_crosses_its_seam_too() {
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push(GeometryNode::Curve3(Curve3::Ellipse(Ellipse3 {
            frame: frame(),
            semi_axis_x: 2.0,
            semi_axis_y: 1.0,
        })))
        .unwrap();
    let arc = trimmed(&mut b, basis, 315f64.to_radians(), 45f64.to_radians(), true);
    let mesh = disk(b, arc, None).expect("the ellipse arc through 0 sweeps");
    // The arc through 0 keeps x >= 2 cos 45; the wrong one reaches x = -2.
    let min_x = mesh
        .positions
        .iter()
        .map(|p| p.x)
        .fold(Scalar::MAX, Scalar::min);
    assert!(
        min_x > 2.0 * 45f64.to_radians().cos() - TUBE - 1e-6,
        "the ellipse arc went round the far side: min x {min_x}"
    );
}

/// A rebar leg: a bend across the seam, then a straight leg off its end.
///
/// This is the Revit shape that was refused as a gap: the bend ran round
/// the wrong side, so its end was nowhere near the leg's start.
fn bend_then_leg(b: &mut GeometryGraphBuilder, reversed: bool, against: bool) -> NodeId {
    let basis = b.push(GeometryNode::Curve3(circle())).unwrap();
    // The same bend written two ways: 270 -> 45 with the basis, or its
    // reverse, 45 -> 270 against the basis, used with `same_sense = false`
    // so the composite still walks it 270 -> 45 into the leg. Getting the
    // trimmed curve's own direction wrong leaves a gap at the leg.
    let bend = if against {
        trimmed(b, basis, 45f64.to_radians(), 270f64.to_radians(), false)
    } else {
        trimmed(b, basis, 270f64.to_radians(), 45f64.to_radians(), true)
    };
    let bend_sense = !against;
    let end = 45f64.to_radians();
    let leg_line = b
        .push(GeometryNode::Curve3(Curve3::Line(Line3 {
            origin: Point3::new(end.cos(), end.sin(), 0.0),
            direction: Vec3::new(-end.sin(), end.cos(), 0.0),
        })))
        .unwrap();
    let leg = trimmed(b, leg_line, 0.0, 1.0, true);
    let segment = |curve, same_sense| CurveSegment {
        curve,
        same_sense,
        transition: Transition::Continuous,
    };
    let segments = if reversed {
        // The same path walked backwards: leg tip -> bend start.
        vec![segment(leg, false), segment(bend, !bend_sense)]
    } else {
        vec![segment(bend, bend_sense), segment(leg, true)]
    };
    b.push(GeometryNode::CurveRelation(CurveRelation::Composite {
        segments,
    }))
    .unwrap()
}

#[test]
fn a_composite_bend_across_the_seam_meets_its_leg() {
    for (reversed, against) in [(false, false), (true, false), (false, true), (true, true)] {
        let what = format!("reversed={reversed} against={against}");
        let mut b = GeometryGraphBuilder::new();
        let path = bend_then_leg(&mut b, reversed, against);
        let mesh =
            disk(b, path, None).unwrap_or_else(|e| panic!("{what}: the bent bar must sweep: {e}"));
        assert_length(&mesh, RADIUS * 135f64.to_radians() + 1.0, &what);
    }
}

/// A sweep range on a trimmed conic is read in the trimmed curve's own
/// unwrapped interval, so a sub-range through the seam is expressible.
#[test]
fn a_range_on_a_seam_crossing_trim_selects_within_it() {
    let arc_of = |range: (Scalar, Scalar)| {
        let mut b = GeometryGraphBuilder::new();
        let basis = b.push(GeometryNode::Curve3(circle())).unwrap();
        let arc = trimmed(&mut b, basis, 315f64.to_radians(), 45f64.to_radians(), true);
        disk(b, arc, Some((range.0.to_radians(), range.1.to_radians())))
    };
    // 330 -> 390 and 330 -> 30 name the same 60 degree sub-arc through 0.
    for range in [(330.0, 390.0), (330.0, 30.0), (-30.0, 30.0)] {
        let mesh = arc_of(range).unwrap_or_else(|e| panic!("{range:?}: {e}"));
        assert_length(&mesh, RADIUS * 60f64.to_radians(), &format!("{range:?}"));
        let mid = mean_angle_deg(&mesh);
        assert!(angle_gap(mid, 0.0) < 3.0, "{range:?}: centred at {mid} deg");
    }
    // Outside the trimmed arc, whichever way it is written.
    assert!(
        arc_of((100.0, 200.0)).is_err(),
        "a range off the arc is refused"
    );
}
