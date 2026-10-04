//! Round holes tangent to a planar face (#243).
//!
//! A round hole touching a plane meets it along a ruling of its cylinder.
//! Placed with exact axis matrices and dyadic sizes, the contact is exact
//! and the exact predicates decide it: the difference succeeds at
//! `Tolerance::ZERO` with an empty report. Under a general placement, or
//! with the hole a fraction of the tolerance into or off the plane (#234),
//! the plane is read as touching (`PlaneTouchesCylinder`), and every face
//! pair the two meet in agrees with that one reading.
//!
//! Cases:
//!
//! - an I-beam (no root fillets) `[-B/2, B/2] x [-H/2, H/2] x [0, LEN]`
//!   whose round web hole, its axis along `x`, touches the top flange's
//!   inner face `y = H/2 - TF`: past the flange tips, past the web, or
//!   flush with the web's faces;
//! - a wall whose round hole touches its top face from inside, or whose
//!   round tool touches it from outside;
//! - a round column clipped by a plane a fraction of the tolerance into it
//!   (the #234 roof case);
//! - the I-beam with root fillets (#249): the hole touches each top fillet
//!   where it meets the flange, a cylinder/cylinder section with a double
//!   point; its volume is checked against a slice-by-slice reference.
//!
//! Each case checks the exact volume against its closed form, the
//! certified distance to a probe in the hole, a clean audit, a closed
//! two-manifold, and the report.

use std::f64::consts::{FRAC_PI_2, PI};

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Mat3, Plane3, Point3, Tolerance, Transform3, Vec3};
use axiolid_measure::{boundary_distance, exact_properties};
use axiolid_mesh_compile::{BooleanReport, ReferenceExactCompiler, ToleranceDecisionKind};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{CircleProfile, Profile, RectangleProfile, SectionProfile};

const EPS: f64 = 1e-6;

/// An I section and its hole: dyadic, so every sum is exact.
#[derive(Clone, Copy)]
struct Beam {
    h: f64,
    b: f64,
    tw: f64,
    tf: f64,
    len: f64,
    r: f64,
}

const DYADIC: Beam = Beam {
    h: 0.5,
    b: 0.25,
    tw: 0.125,
    tf: 0.0625,
    len: 2.0,
    r: 0.125,
};

/// IPE 300 sizes without root fillets, a hole of radius 50 mm: nothing is
/// dyadic, so even unplaced the tangency is only within rounding.
const IPE: Beam = Beam {
    h: 0.3,
    b: 0.15,
    tw: 0.0071,
    tf: 0.0107,
    len: 2.0,
    r: 0.05,
};

impl Beam {
    fn section(self) -> Profile {
        Profile::Section(SectionProfile::I {
            depth: self.h,
            width: self.b,
            web_thickness: self.tw,
            flange_thickness: self.tf,
            fillet_radius: None,
            flange_edge_radius: None,
            flange_slope: None,
        })
    }

    /// The hole's axis height: touching the top flange's inner face.
    fn axis(self) -> f64 {
        self.h / 2.0 - self.tf - self.r
    }

    fn volume(self) -> f64 {
        let area = 2.0 * self.b * self.tf + (self.h - 2.0 * self.tf) * self.tw;
        area * self.len - PI * self.r * self.r * self.tw
    }
}

/// How far the hole runs along `x`: past the flange tips, past the web
/// only, or flush with the web's faces.
#[derive(Clone, Copy, Debug)]
enum Run {
    PastFlanges,
    PastWeb,
    Flush,
}

impl Run {
    const ALL: [Run; 3] = [Run::PastFlanges, Run::PastWeb, Run::Flush];

    /// Where it starts and how deep it runs.
    fn span(self, beam: Beam) -> (f64, f64) {
        match self {
            Run::PastFlanges => (-0.75 * beam.b, 1.5 * beam.b),
            Run::PastWeb => (-beam.tw, 2.0 * beam.tw),
            Run::Flush => (-beam.tw / 2.0, beam.tw),
        }
    }
}

struct Graph(GeometryGraphBuilder);

impl Graph {
    fn new() -> Self {
        Self(GeometryGraphBuilder::new())
    }

    fn push(&mut self, node: GeometryNode) -> NodeId {
        self.0.push(node).expect("a valid node")
    }

    fn extrusion(&mut self, profile: Profile, depth: f64) -> NodeId {
        let profile = self.push(GeometryNode::Profile(profile));
        self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
    }

    fn place(&mut self, source: NodeId, transform: Transform3) -> NodeId {
        self.push(GeometryNode::Instance(Instance { source, transform }))
    }

    fn minus(&mut self, left: NodeId, right: NodeId) -> NodeId {
        self.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator: BooleanOperator::Difference,
        }))
    }

    /// A box `[min, max]` under `placement`.
    fn probe(&mut self, min: Vec3, max: Vec3, placement: Transform3) -> NodeId {
        let size = max - min;
        let block = self.extrusion(rect(size.x, size.y), size.z);
        let centre = Vec3::new((min.x + max.x) / 2.0, (min.y + max.y) / 2.0, min.z);
        self.place(block, placement * Transform3::from_translation(centre))
    }

    fn finish(self, roots: Vec<NodeId>) -> GeometryGraph {
        self.0.finish(roots).expect("a valid graph")
    }
}

fn rect(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

fn circle(radius: f64) -> Profile {
    Profile::Circle(CircleProfile {
        radius,
        thickness: None,
    })
}

/// The extrusion's `+z` turned onto `+x` with entries `0` and `+-1`.
fn exact_onto_x() -> Transform3 {
    Transform3::from_mat3(Mat3::from_cols(-Vec3::Z, Vec3::Y, Vec3::X))
}

/// The extrusion's `+z` turned onto `-y` with entries `0` and `+-1`.
fn exact_onto_minus_y() -> Transform3 {
    Transform3::from_mat3(Mat3::from_cols(Vec3::X, Vec3::Z, -Vec3::Y))
}

/// Placements with exact entries: none, and a quarter or half turn with a
/// dyadic offset.
fn exact_placements() -> [Transform3; 3] {
    let quarter = Mat3::from_cols(Vec3::Y, -Vec3::X, Vec3::Z);
    let half = Mat3::from_cols(-Vec3::X, -Vec3::Y, Vec3::Z);
    [
        Transform3::IDENTITY,
        Transform3::from_mat3_translation(quarter, Vec3::new(12.5, -4.0, 3.25)),
        Transform3::from_mat3_translation(half, Vec3::new(-7.75, 2.5, 0.0)),
    ]
}

/// A building placement: turned about `z` and moved.
fn building() -> Transform3 {
    Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2)) * Transform3::from_rotation_z(0.6)
}

/// A rigid placement with no axis left aligned.
fn general() -> Transform3 {
    Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7)
}

fn audited(brep: &ExactBRep) {
    let health = geometric_audit(brep, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topology = axiolid_topology::audit_brep(brep.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, Tolerance::METRE)
        .expect("measurable")
        .signed_volume
}

/// A compiled body with its report and its probe, checked for audit,
/// volume (within `slack`) and the probe's certified distance.
fn checked(
    graph: &GeometryGraph,
    (body, probe): (NodeId, NodeId),
    tolerance: Tolerance,
    (expected, slack): (f64, f64),
    distance: f64,
) -> Result<BooleanReport, GeomError> {
    let mut solids = ReferenceExactCompiler::new().compile_exact_batch_with_reports(
        graph,
        &[body, probe],
        &ExecutionOptions::new(tolerance),
    )?;
    let (probe, _) = solids.pop().expect("probe");
    let (body, report) = solids.pop().expect("body");
    audited(&body);
    let measured = volume(&body);
    assert!(
        (measured - expected).abs() <= slack + 1e-12 * expected,
        "volume {measured}, expected {expected} (slack {slack})"
    );
    let bounds = boundary_distance(&body, &probe, 1e-9, Tolerance::METRE).expect("bounded");
    assert!(
        bounds.lower <= distance + 1e-12 + slack && distance <= bounds.upper + 1e-12 + slack,
        "[{}, {}] must contain {distance}",
        bounds.lower,
        bounds.upper
    );
    Ok(report)
}

/// The beam under `placement` minus its web hole, `gap` above touching
/// the flange (positive: into it), the hole turned onto `x` by `turn`,
/// with a square probe of half side `r / 2` along the hole's axis through
/// the web: the certified distance from the cut beam to it is
/// `r - r / sqrt 2` exactly when the hole is where it should be.
fn beam_case(
    beam: Beam,
    run: Run,
    gap: f64,
    turn: Transform3,
    placement: Transform3,
) -> (GeometryGraph, (NodeId, NodeId), f64) {
    let mut g = Graph::new();
    let body = g.extrusion(beam.section(), beam.len);
    let body = g.place(body, placement);
    let (start, depth) = run.span(beam);
    let hole = g.extrusion(circle(beam.r), depth);
    let local =
        Transform3::from_translation(Vec3::new(start, beam.axis() + gap, beam.len / 2.0)) * turn;
    let hole = g.place(hole, placement * local);
    let body = g.minus(body, hole);
    let a = beam.r / 2.0;
    let (y, z) = (beam.axis() + gap, beam.len / 2.0);
    let probe = g.probe(
        Vec3::new(-beam.b, y - a, z - a),
        Vec3::new(beam.b, y + a, z + a),
        placement,
    );
    let graph = g.finish(vec![body, probe]);
    (graph, (body, probe), beam.r - a * 2.0_f64.sqrt())
}

fn decision(report: &BooleanReport, kind: ToleranceDecisionKind) -> Option<f64> {
    report
        .decisions()
        .iter()
        .find(|d| d.kind == kind)
        .map(|d| d.linear)
}

#[test]
fn a_web_hole_touching_the_flange_is_exact_at_zero_tolerance() {
    // Exactly tangent: the exact predicates decide the contact, nothing is
    // read within tolerance, and the volume is the closed form's.
    for placement in exact_placements() {
        for run in Run::ALL {
            let (graph, roots, distance) = beam_case(DYADIC, run, 0.0, exact_onto_x(), placement);
            let report = checked(
                &graph,
                roots,
                Tolerance::ZERO,
                (DYADIC.volume(), 0.0),
                distance,
            )
            .unwrap_or_else(|e| panic!("{run:?}: {e}"));
            assert!(report.is_exact(), "{run:?}: {report:?}");
            // A positive tolerance changes nothing for an exact contact.
            let report = checked(
                &graph,
                roots,
                Tolerance::METRE,
                (DYADIC.volume(), 0.0),
                distance,
            )
            .unwrap_or_else(|e| panic!("{run:?}: {e}"));
            assert!(report.is_exact(), "{run:?}: {report:?}");
        }
    }
}

#[test]
fn a_web_hole_touching_the_flange_under_a_general_placement_reads_the_contact() {
    // Placed separately under a general rotation, or turned by sines and
    // cosines, the hole touches the flange only up to rounding: read as
    // touching, reported, within the tolerance.
    let rounded = Transform3::from_rotation_y(FRAC_PI_2);
    for beam in [DYADIC, IPE] {
        for placement in [Transform3::IDENTITY, building(), general()] {
            for turn in [exact_onto_x(), rounded] {
                for run in Run::ALL {
                    let (graph, roots, distance) = beam_case(beam, run, 0.0, turn, placement);
                    let report = checked(
                        &graph,
                        roots,
                        Tolerance::METRE,
                        (beam.volume(), 1e-12),
                        distance,
                    )
                    .unwrap_or_else(|e| panic!("{run:?}: {e}"));
                    assert!(report.linear() <= EPS, "{run:?}: {report:?}");
                    assert!(report.angular() <= Tolerance::METRE.angular());
                    // Placed and turned by exact matrices, the pair may be
                    // exactly tangent on the numbers given (it is for both
                    // beams here) and need no reading; otherwise it is read.
                    let exact = placement == Transform3::IDENTITY && turn == exact_onto_x();
                    if !exact {
                        assert!(
                            report.contains(ToleranceDecisionKind::PlaneTouchesCylinder),
                            "{run:?}: {report:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn a_rounded_web_hole_touching_the_flange_is_refused_at_zero_tolerance() {
    // Nothing may be read within a zero tolerance: a contact that holds
    // only up to rounding is refused by name, never guessed.
    let (graph, (body, _), _) = beam_case(DYADIC, Run::PastWeb, 0.0, exact_onto_x(), general());
    let error = ReferenceExactCompiler::new()
        .compile_exact_with_report(&graph, body, &ExecutionOptions::new(Tolerance::ZERO))
        .expect_err("refused");
    assert!(
        matches!(&error, GeomError::UnsupportedInput { input, .. } if input.contains("exact boolean")),
        "{error:?}"
    );
}

/// The beam's root fillets (#249): radius `1/32` for the dyadic beam,
/// IPE 300's `15 mm`.
fn fillet(beam: Beam) -> f64 {
    if beam.h == DYADIC.h {
        0.03125
    } else {
        0.015
    }
}

/// The area of the hole's cross-section above the line `d` above its axis.
fn segment(r: f64, d: f64) -> f64 {
    if d >= r {
        0.0
    } else if d <= -r {
        PI * r * r
    } else {
        r * r * (d / r).acos() - d * (r * r - d * d).sqrt()
    }
}

/// The volume of the filleted beam less a hole `gap` above touching the
/// flange whose run covers both top fillets and no more than the flanges.
///
/// Sliced across the hole's axis: through the web the whole disc; through
/// a fillet of radius `rf`, at `x = x_c - rf sin t`, the material lies
/// above `y_c + rf cos t`, so the slice is a circular segment, summed by
/// Simpson's rule in `t` (smooth there); under the flange the material
/// lies above the flange's inner face.
fn filleted_volume(beam: Beam, run: Run, gap: f64) -> f64 {
    let rf = fillet(beam);
    let (top, x_c) = (beam.h / 2.0 - beam.tf, beam.tw / 2.0 + rf);
    let y = beam.axis() + gap;
    let (start, depth) = run.span(beam);
    assert!(start <= -x_c && start + depth >= x_c);
    let n = 200_000;
    let step = 0.5 * PI / n as f64;
    let mut slices = 0.0;
    for i in 0..=n {
        let t = i as f64 * step;
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        slices += w * segment(beam.r, top - rf + rf * t.cos() - y) * rf * t.cos();
    }
    let flange = (2.0 * (start + depth).min(beam.b / 2.0) - 2.0 * x_c).max(0.0);
    let removed = PI * beam.r * beam.r * beam.tw
        + 2.0 * slices * step / 3.0
        + flange * segment(beam.r, top - y);
    let area = 2.0 * beam.b * beam.tf
        + (beam.h - 2.0 * beam.tf) * beam.tw
        + 4.0 * rf * rf * (1.0 - PI / 4.0);
    area * beam.len - removed
}

/// [`beam_case`] for the beam with root fillets.
fn filleted_case(
    beam: Beam,
    run: Run,
    gap: f64,
    turn: Transform3,
    placement: Transform3,
) -> (GeometryGraph, (NodeId, NodeId), f64) {
    let mut g = Graph::new();
    let section = Profile::Section(SectionProfile::I {
        depth: beam.h,
        width: beam.b,
        web_thickness: beam.tw,
        flange_thickness: beam.tf,
        fillet_radius: Some(fillet(beam)),
        flange_edge_radius: None,
        flange_slope: None,
    });
    let body = g.extrusion(section, beam.len);
    let body = g.place(body, placement);
    let (start, depth) = run.span(beam);
    let hole = g.extrusion(circle(beam.r), depth);
    let local =
        Transform3::from_translation(Vec3::new(start, beam.axis() + gap, beam.len / 2.0)) * turn;
    let hole = g.place(hole, placement * local);
    let body = g.minus(body, hole);
    // A probe along the hole's axis inside the web, square of half side
    // `r / 2` across it: the nearest material is the hole's wall in the
    // web, `r - r / sqrt 2` from its long edges. It stays inside the web
    // (half the web's thickness long), so the distance has a short
    // plateau to certify, not one running past the fillets.
    let a = beam.r / 2.0;
    let (y, z) = (beam.axis() + gap, beam.len / 2.0);
    let probe = g.probe(
        Vec3::new(-beam.tw / 4.0, y - a, z - a),
        Vec3::new(beam.tw / 4.0, y + a, z + a),
        placement,
    );
    let graph = g.finish(vec![body, probe]);
    (graph, (body, probe), beam.r - a * 2.0_f64.sqrt())
}

#[test]
fn a_web_hole_touching_the_flange_of_a_filleted_beam_is_exact_at_zero_tolerance() {
    // With root fillets the hole's top ruling, tangent to the flange, is
    // also tangent to each fillet where the fillet meets the flange: the
    // hole and the fillet cylinder (axes crossing at right angles) meet in
    // a quartic with a double point there. The trace ends both of its
    // loops at that point, which is the exact double root of the
    // fillet/flange edge against the hole: nothing is read within
    // tolerance, and the volume is the reference's.
    for placement in exact_placements() {
        for run in [Run::PastFlanges, Run::PastWeb] {
            let (graph, roots, distance) =
                filleted_case(DYADIC, run, 0.0, exact_onto_x(), placement);
            let expected = filleted_volume(DYADIC, run, 0.0);
            for tolerance in [Tolerance::ZERO, Tolerance::METRE] {
                let report = checked(&graph, roots, tolerance, (expected, 0.0), distance)
                    .unwrap_or_else(|e| panic!("{run:?}: {e}"));
                assert!(
                    report.linear() <= report.rounding_floor(),
                    "{run:?}: {report:?}"
                );
                if tolerance == Tolerance::ZERO {
                    assert!(report.is_exact(), "{run:?}: {report:?}");
                }
            }
        }
    }
}

#[test]
fn a_web_hole_touching_the_flange_of_a_filleted_beam_under_a_general_placement() {
    // The flange is read as touching the hole (#243), and the fillets'
    // double points are placed on that contact: the dyadic beam and an
    // IPE 300 with its 15 mm root fillets.
    let rounded = Transform3::from_rotation_y(FRAC_PI_2);
    for beam in [DYADIC, IPE] {
        for placement in [building(), general()] {
            for turn in [exact_onto_x(), rounded] {
                let run = Run::PastFlanges;
                let (graph, roots, distance) = filleted_case(beam, run, 0.0, turn, placement);
                let report = checked(
                    &graph,
                    roots,
                    Tolerance::METRE,
                    (filleted_volume(beam, run, 0.0), 1e-12),
                    distance,
                )
                .unwrap_or_else(|e| panic!("{run:?} {placement:?} {turn:?}: {e}"));
                assert!(report.linear() <= EPS, "{report:?}");
                assert!(
                    report.contains(ToleranceDecisionKind::PlaneTouchesCylinder),
                    "{report:?}"
                );
            }
        }
    }
}

#[test]
fn a_filleted_beam_whose_hole_is_a_fraction_of_the_tolerance_off_the_flange_is_refused() {
    // Read as touching the flange, the hole would touch each fillet too,
    // where it meets it in two arcs about `sqrt(2 r gap)` from the
    // contact: no single move within the tolerance explains both, so it
    // is refused by name, never sewn inconsistently.
    for fraction in [0.5, -0.5] {
        let (graph, (body, _), _) = filleted_case(
            DYADIC,
            Run::PastFlanges,
            fraction * EPS,
            exact_onto_x(),
            Transform3::IDENTITY,
        );
        let error = ReferenceExactCompiler::new()
            .compile_exact_with_report(&graph, body, &ExecutionOptions::new(Tolerance::METRE))
            .expect_err("refused");
        assert!(
            matches!(&error, GeomError::UnsupportedInput { input, .. } if input.contains("contact")),
            "{fraction}: {error:?}"
        );
    }
}

#[test]
fn a_web_hole_a_fraction_of_the_tolerance_into_or_off_the_flange_reads_as_touching() {
    // The #234 open item on an I-beam: the hole reaches a fraction of the
    // tolerance into the flange (or stops short of it). Read as touching,
    // the flange plane moves by that fraction, and every face pair agrees:
    // the web circle, the cap chords and the flange's edges all meet the
    // one contact ruling.
    for placement in [Transform3::IDENTITY, general()] {
        for run in Run::ALL {
            for fraction in [0.1, 0.5, 0.9, -0.1, -0.5, -0.9] {
                let gap = fraction * EPS;
                let (graph, roots, distance) =
                    beam_case(DYADIC, run, gap, exact_onto_x(), placement);
                // The flange moved by `gap` changes the volume by at most
                // `gap` times the beam's surface; so does reading the
                // vertices it moves as one.
                let surface = 2.0 * (DYADIC.b + DYADIC.h + DYADIC.b) * DYADIC.len;
                let report = checked(
                    &graph,
                    roots,
                    Tolerance::METRE,
                    (DYADIC.volume(), gap.abs() * surface),
                    distance,
                )
                .unwrap_or_else(|e| panic!("{run:?} {fraction}: {e}"));
                let moved = decision(&report, ToleranceDecisionKind::PlaneTouchesCylinder)
                    .unwrap_or_else(|| panic!("{run:?} {fraction}: {report:?}"));
                assert!(
                    (moved - gap.abs()).abs() <= 1e-9 * EPS + 1e-15,
                    "{run:?} {fraction}: moved {moved}"
                );
                assert!(report.linear() <= EPS, "{run:?} {fraction}: {report:?}");
            }
        }
    }
}

#[test]
fn a_web_hole_ten_tolerances_into_the_flange_cuts_its_groove() {
    // Further than the tolerance the exact predicates decide: the hole
    // cuts a groove `10 eps` deep into the flange's underside wherever it
    // runs past the web, and nothing is read as touching.
    let beam = DYADIC;
    let depth = 10.0 * EPS;
    let (r, d) = (beam.r, beam.r - depth);
    let segment = r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    for (run, groove) in [
        (Run::PastFlanges, beam.b - beam.tw),
        (Run::PastWeb, beam.tw),
    ] {
        for placement in [Transform3::IDENTITY, general()] {
            let (graph, roots, distance) = beam_case(beam, run, depth, exact_onto_x(), placement);
            let expected = beam.volume() - segment * groove;
            let report = checked(
                &graph,
                roots,
                Tolerance::METRE,
                (expected, 1e-3 * segment * groove),
                distance,
            )
            .unwrap_or_else(|e| panic!("{run:?}: {e}"));
            assert!(
                !report.contains(ToleranceDecisionKind::PlaneTouchesCylinder),
                "{run:?}: {report:?}"
            );
        }
    }
}

/// Wall length, thickness and height, and a round hole's radius.
const L: f64 = 6.0;
const T: f64 = 0.25;
const H: f64 = 3.0;
const R: f64 = 0.5;

/// A wall `[-L/2, L/2] x [-T/2, T/2] x [0, H]` under `placement` minus a
/// round hole across it (its axis along `y`, at `x = 0.5`), touching the
/// top face from inside (`inside`) or from outside, `gap` above that; with
/// a square probe of half side `R / 2` along the hole's axis.
fn wall_case(
    inside: bool,
    gap: f64,
    placement: Transform3,
) -> (GeometryGraph, (NodeId, NodeId), f64) {
    let mut g = Graph::new();
    let wall = g.extrusion(rect(L, T), H);
    let wall = g.place(wall, placement);
    let hole = g.extrusion(circle(R), 1.0);
    let z = if inside { H - R } else { H + R } + gap;
    let local = Transform3::from_translation(Vec3::new(0.5, 0.5, z)) * exact_onto_minus_y();
    let hole = g.place(hole, placement * local);
    let body = g.minus(wall, hole);
    let a = R / 2.0;
    let probe = g.probe(
        Vec3::new(0.5 - a, -1.0, z - a),
        Vec3::new(0.5 + a, 1.0, z + a),
        placement,
    );
    let graph = g.finish(vec![body, probe]);
    // Inside, the probe's corners are nearest the hole's surface; outside
    // the wall, its bottom is nearest the wall's top face.
    let distance = if inside {
        R - a * 2.0_f64.sqrt()
    } else {
        R - a
    };
    (graph, (body, probe), distance)
}

#[test]
fn a_hole_touching_a_wall_face_from_inside_or_outside() {
    let through = L * T * H - PI * R * R * T;
    for (inside, expected) in [(true, through), (false, L * T * H)] {
        // Exactly tangent at zero tolerance: exact, nothing read.
        for placement in exact_placements() {
            let (graph, roots, distance) = wall_case(inside, 0.0, placement);
            let report = checked(&graph, roots, Tolerance::ZERO, (expected, 0.0), distance)
                .unwrap_or_else(|e| panic!("inside {inside}: {e}"));
            assert!(report.is_exact(), "inside {inside}: {report:?}");
        }
        // Under general placements, and a fraction of the tolerance either
        // way: read as touching, and reported.
        for placement in [building(), general()] {
            for gap in [0.0, 0.5 * EPS, -0.5 * EPS] {
                let (graph, roots, distance) = wall_case(inside, gap, placement);
                let report = checked(
                    &graph,
                    roots,
                    Tolerance::METRE,
                    (expected, gap.abs() * 2.0 * (L * T + L * H + T * H)),
                    distance,
                )
                .unwrap_or_else(|e| panic!("inside {inside} gap {gap}: {e}"));
                assert!(
                    report.contains(ToleranceDecisionKind::PlaneTouchesCylinder),
                    "inside {inside} gap {gap}: {report:?}"
                );
                assert!(report.linear() <= EPS);
            }
        }
    }
}

/// A round column of radius `r`, `H` tall, under a general placement,
/// clipped by the vertical plane `x = r - inset` from outside (#234).
fn column_clip(r: f64, inset: f64) -> Result<(ExactBRep, BooleanReport), GeomError> {
    let placement = general();
    let mut g = Graph::new();
    let column = g.extrusion(circle(r), H);
    let column = g.place(column, placement);
    let cut = g.push(GeometryNode::HalfSpace(HalfSpace {
        boundary: Plane3 {
            origin: placement.transform_point3(Point3::new(r - inset, 0.0, 0.0)),
            normal: placement.transform_vector3(Vec3::X),
        },
        agreement: true,
    }));
    let body = g.minus(column, cut);
    let graph = g.finish(vec![body]);
    ReferenceExactCompiler::new().compile_exact_with_report(
        &graph,
        body,
        &ExecutionOptions::new(Tolerance::METRE),
    )
}

#[test]
fn a_roof_plane_a_fraction_of_the_tolerance_into_a_column_reads_as_touching() {
    // The #234 open item: the plane reaches a fraction of the tolerance
    // into the column. Read as touching, the plane moves out by that much
    // and keeps the whole column, its caps' chords and its side agreeing.
    let r = 0.3;
    let whole = PI * r * r * H;
    for fraction in [0.1, 0.5, 0.9] {
        let inset = fraction * EPS;
        let (solid, report) = column_clip(r, inset).unwrap_or_else(|e| panic!("{fraction}: {e}"));
        audited(&solid);
        let measured = volume(&solid);
        assert!(
            (measured - whole).abs() <= 2.0 * r * H * EPS,
            "{fraction}: volume {measured}"
        );
        let moved = decision(&report, ToleranceDecisionKind::PlaneTouchesCylinder)
            .unwrap_or_else(|| panic!("{fraction}: {report:?}"));
        assert!((moved - inset).abs() <= 1e-9 * EPS, "{fraction}: {moved}");
        assert!(report.linear() <= EPS, "{fraction}: {report:?}");
    }
}
