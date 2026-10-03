//! A boolean's mesh is measured against its exact result (#235).
//!
//! Each case compiles a boolean with `compile_mesh_with_deviation` (the
//! mesh from the boolmesh provider), compiles the same node exactly with
//! `ReferenceExactCompiler`, samples that exact result densely (a lattice
//! of each face's parameter box kept where `axiolid_measure::FaceDomain`
//! certifies the point inside, plus every trim curve, where the cut lies),
//! and asserts that no sample is further from the mesh than the report's
//! bound: the report's one-sided quantity, exact surface to mesh, relative
//! to the exact compiler's result. It also asserts the bound is useful:
//! within the requested budget, or not absurdly above the measured maximum
//! (the search stops once it is within the budget).
//!
//! Cases: a placed wall with a round window, an I-beam with root fillets
//! and round holes through its web, a slab with a round shaft, and a gable
//! wall clipped by two roof half-spaces (#234) with a round window. A union
//! of placed operands, which the exact compiler refuses, stays unbounded
//! with the refusal's name.

use std::collections::HashMap;
use std::f64::consts::FRAC_PI_2;
use std::time::Instant;

use axiolid_brep::ExactBRep;
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Plane3, Point2, Point3, Scalar, Tolerance, Transform3, Vec3};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_measure::FaceDomain;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{
    DeviationBound, DeviationPath, DeviationReport, ReferenceExactCompiler, ReferenceMeshCompiler,
};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{CircleProfile, Profile, RectangleProfile, SectionProfile};

fn options(budget: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(budget.min(1e-3), 1e-9).unwrap())
        .with_chord_error(budget)
        .unwrap()
}

struct Graph(GeometryGraphBuilder);

impl Graph {
    fn push(&mut self, node: GeometryNode) -> NodeId {
        self.0.push(node).expect("a valid node")
    }

    fn extrusion(&mut self, profile: Profile, depth: Scalar) -> NodeId {
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

    fn boolean(&mut self, left: NodeId, right: NodeId, operator: BooleanOperator) -> NodeId {
        self.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator,
        }))
    }
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

fn circle(radius: Scalar) -> Profile {
    Profile::Circle(CircleProfile {
        radius,
        thickness: None,
    })
}

/// A wall placement as a building model has one: turned about `z` and moved.
fn building() -> Transform3 {
    Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2)) * Transform3::from_rotation_z(0.6)
}

/// Compile `root` to a mesh with its report, and exactly.
fn compile(
    build: impl FnOnce(&mut Graph) -> NodeId,
    options: &ExecutionOptions,
) -> (GeometryGraph, NodeId, TriMesh, DeviationReport, Scalar) {
    let mut g = Graph(GeometryGraphBuilder::new());
    let root = build(&mut g);
    let graph = g.0.finish(vec![root]).expect("a valid graph");
    let start = Instant::now();
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, options)
        .expect("the body compiles");
    let seconds = start.elapsed().as_secs_f64();
    (graph, root, outcome.mesh, report, seconds)
}

/// Triangles binned on a uniform grid, each in every cell its box touches
/// once grown by `reach`.
struct TriangleGrid {
    triangles: Vec<[Point3; 3]>,
    cells: HashMap<[i64; 3], Vec<usize>>,
    cell: Scalar,
    reach: Scalar,
}

impl TriangleGrid {
    fn new(mesh: &TriMesh, reach: Scalar) -> Self {
        let triangles: Vec<[Point3; 3]> = mesh
            .indices
            .chunks_exact(3)
            .map(|t| [0, 1, 2].map(|i| mesh.positions[t[i] as usize]))
            .collect();
        let (lo, hi) = mesh.positions.iter().fold(
            (
                Vec3::splat(Scalar::INFINITY),
                Vec3::splat(Scalar::NEG_INFINITY),
            ),
            |(lo, hi), p| (lo.min(*p), hi.max(*p)),
        );
        let cell = ((hi - lo).max_element() / 64.0).max(4.0 * reach);
        let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        for (i, t) in triangles.iter().enumerate() {
            let a = t[0].min(t[1]).min(t[2]) - Vec3::splat(reach);
            let b = t[0].max(t[1]).max(t[2]) + Vec3::splat(reach);
            let (ka, kb) = (key(a, cell), key(b, cell));
            for x in ka[0]..=kb[0] {
                for y in ka[1]..=kb[1] {
                    for z in ka[2]..=kb[2] {
                        cells.entry([x, y, z]).or_default().push(i);
                    }
                }
            }
        }
        Self {
            triangles,
            cells,
            cell,
            reach,
        }
    }

    fn distance(&self, p: Point3) -> Scalar {
        let mut best = self.reach;
        for &i in self.cells.get(&key(p, self.cell)).into_iter().flatten() {
            if let Ok(q) = closest_point_on_triangle(p, self.triangles[i]) {
                best = best.min(p.distance(q));
            }
        }
        best
    }
}

fn key(p: Point3, cell: Scalar) -> [i64; 3] {
    [p.x, p.y, p.z].map(|v| (v / cell).floor() as i64)
}

/// Points of the exact result: a `n x n` lattice of each face's parameter
/// box where the face certainly contains it, and `n` points along every
/// pcurve span of every loop.
fn exact_samples(exact: &ExactBRep, n: usize) -> Vec<Point3> {
    let topology = exact.topology();
    let mut out = Vec::new();
    for (index, face) in topology.faces().iter().enumerate() {
        let surface = &exact.surfaces()[face.surface.expect("a surface").index()];
        let domain = FaceDomain::new(exact, topology.face_id_at(index).unwrap(), Tolerance::METRE)
            .expect("a domain")
            .expect("a monotone split");
        let (lo, hi) = domain.bounds();
        for i in 0..n {
            for j in 0..n {
                let at = Point2::new(
                    lo.x + (hi.x - lo.x) * (i as Scalar + 0.5) / n as Scalar,
                    lo.y + (hi.y - lo.y) * (j as Scalar + 0.5) / n as Scalar,
                );
                if domain.contains(at).expect("decidable") == Some(true) {
                    out.push(axiolid_reference::surface::evaluate(surface, at.x, at.y).unwrap());
                }
            }
        }
        for bound in &face.bounds {
            let wire = &topology.loops()[bound.loop_id.index()];
            for (k, edge_use) in wire.edges.iter().enumerate() {
                let curve = &exact.curves2()[edge_use.pcurve.expect("a pcurve").index()];
                let span = exact.pcurve_interval(bound.loop_id, k).expect("a span");
                for s in 0..=n {
                    let t = span.start + (span.end - span.start) * s as Scalar / n as Scalar;
                    let uv = axiolid_reference::evaluate2(curve, t).unwrap();
                    out.push(axiolid_reference::surface::evaluate(surface, uv.x, uv.y).unwrap());
                }
            }
        }
    }
    out
}

/// Assert the report is certified, measured against the exact result, and
/// that no exact sample lies further from the mesh than it says. Returns
/// `(reported, measured)`.
fn assert_measured(
    name: &str,
    case: (GeometryGraph, NodeId, TriMesh, DeviationReport, Scalar),
    options: &ExecutionOptions,
    ceiling: Scalar,
) -> (Scalar, Scalar) {
    let (graph, root, mesh, report, seconds) = case;
    assert_eq!(report.contributions.len(), 1, "{report:?}");
    let contribution = report.contributions[0];
    assert_eq!(contribution.path, DeviationPath::Boolean);
    // Within tolerance, or the exact boolean of the given operands (#236).
    assert!(
        contribution
            .detail
            .starts_with("measured against the exact compiler's result, "),
        "{report:?}"
    );
    let DeviationBound::Certified(bound) = contribution.bound else {
        panic!("{name}: expected a certified bound, got {report:?}");
    };
    assert_eq!(report.bound, Some(bound));
    let exact = ReferenceExactCompiler::new()
        .compile_exact(&graph, root, options)
        .expect("the exact result");
    let samples = exact_samples(&exact, 400);
    let grid = TriangleGrid::new(&mesh, 2.0 * bound);
    let worst = samples
        .iter()
        .map(|p| grid.distance(*p))
        .fold(0.0, Scalar::max);
    let ratio = bound / worst;
    eprintln!(
        "{name}: reported {bound:e}, measured {worst:e}, ratio {ratio:.2}, requested {:e}, \
         {} triangles, {} samples, {seconds:.3} s",
        report.requested,
        mesh.indices.len() / 3,
        samples.len()
    );
    assert!(
        worst <= bound,
        "{name}: an exact point lies {worst:e} from the mesh, above the reported {bound:e}"
    );
    assert!(
        ratio <= ceiling || bound <= report.requested,
        "{name}: bound {bound:e} is {ratio:.1}x the measured {worst:e}, above the budget"
    );
    (bound, worst)
}

/// A wall `6 x 0.3 x 3` with a round window of radius `0.4` cut across it,
/// both under a building placement.
fn wall_with_round_window(g: &mut Graph) -> NodeId {
    let wall = g.extrusion(rect(6.0, 0.3), 3.0);
    let wall = g.place(wall, building());
    let window = g.extrusion(circle(0.4), 1.0);
    let across = Transform3::from_translation(Vec3::new(0.8, 0.5, 1.5))
        * Transform3::from_rotation_x(FRAC_PI_2);
    let window = g.place(window, building() * across);
    g.boolean(wall, window, BooleanOperator::Difference)
}

#[test]
fn a_wall_with_a_round_window_is_certified_against_its_exact_result() {
    for budget in [1e-2, 1e-3] {
        let options = options(budget);
        let case = compile(wall_with_round_window, &options);
        let (bound, _) = assert_measured("wall, round window", case, &options, 2.0);
        assert!(bound <= 2.0 * budget, "{bound} for a budget of {budget}");
    }
}

/// An IPE 300-like beam, 4 m long, with root fillets, and three round holes
/// of radius 0.05 through its web.
fn beam_with_round_holes(g: &mut Graph) -> NodeId {
    let section = Profile::Section(SectionProfile::I {
        depth: 0.3,
        width: 0.15,
        web_thickness: 0.0071,
        flange_thickness: 0.0107,
        fillet_radius: Some(0.015),
        flange_edge_radius: None,
        flange_slope: None,
    });
    let mut beam = g.extrusion(section, 4.0);
    let hole = g.extrusion(circle(0.05), 0.2);
    for z in [0.8, 2.0, 3.2] {
        let across = Transform3::from_translation(Vec3::new(-0.1, 0.02, z))
            * Transform3::from_rotation_y(FRAC_PI_2);
        let placed = g.place(hole, across);
        beam = g.boolean(beam, placed, BooleanOperator::Difference);
    }
    beam
}

#[test]
fn an_i_beam_with_round_holes_is_certified_against_its_exact_result() {
    let budget = 1e-3;
    let options = options(budget);
    let case = compile(beam_with_round_holes, &options);
    let (bound, _) = assert_measured("I-beam, round holes", case, &options, 2.0);
    assert!(bound <= 2.0 * budget, "{bound}");
}

/// An IPE 300-like beam without root fillets, 4 m long, under a building
/// placement, whose round web hole of radius 0.05, flush with the web's
/// faces, touches the top flange's inner face (#243).
fn beam_with_a_hole_touching_its_flange(g: &mut Graph) -> NodeId {
    let (depth, web, flange, radius) = (0.3, 0.0071, 0.0107, 0.05);
    let section = Profile::Section(SectionProfile::I {
        depth,
        width: 0.15,
        web_thickness: web,
        flange_thickness: flange,
        fillet_radius: None,
        flange_edge_radius: None,
        flange_slope: None,
    });
    let beam = g.extrusion(section, 4.0);
    let beam = g.place(beam, building());
    let hole = g.extrusion(circle(radius), web);
    let across =
        Transform3::from_translation(Vec3::new(-web / 2.0, depth / 2.0 - flange - radius, 2.0))
            * Transform3::from_rotation_y(FRAC_PI_2);
    let hole = g.place(hole, building() * across);
    g.boolean(beam, hole, BooleanOperator::Difference)
}

#[test]
fn an_i_beam_with_a_hole_touching_its_flange_is_certified_against_its_exact_result() {
    // Its exact result reads the hole as touching the flange (#243), so the
    // mesh is measured against it rather than left unbounded.
    let budget = 1e-3;
    let options = options(budget);
    let case = compile(beam_with_a_hole_touching_its_flange, &options);
    let (bound, _) = assert_measured("I-beam, hole touching the flange", case, &options, 2.0);
    assert!(bound <= 2.0 * budget, "{bound}");
}

#[test]
fn a_slab_with_a_round_shaft_is_certified_against_its_exact_result() {
    let budget = 1e-3;
    let options = options(budget);
    let case = compile(
        |g| {
            let slab = g.extrusion(rect(4.0, 3.0), 0.25);
            let slab = g.place(slab, building());
            let shaft = g.extrusion(circle(0.3), 0.75);
            let local = Transform3::from_translation(Vec3::new(-1.0, 0.4, -0.25));
            let shaft = g.place(shaft, building() * local);
            g.boolean(slab, shaft, BooleanOperator::Difference)
        },
        &options,
    );
    assert_measured("slab, round shaft", case, &options, 2.0);
}

/// A gable wall under two roof planes `z = 2.4 -+ 0.3 x` (world
/// half-spaces, as the mesh compiler reads them) with a round window.
#[test]
fn a_roof_clipped_wall_is_certified_against_its_exact_result() {
    let budget = 1e-3;
    let options = options(budget);
    let case = compile(
        |g| {
            let wall = g.extrusion(rect(6.0, 0.3), 3.0);
            let wall = g.place(wall, building());
            let mut body = wall;
            for s in [0.3, -0.3] {
                let roof = g.push(GeometryNode::HalfSpace(HalfSpace {
                    boundary: Plane3 {
                        origin: building().transform_point3(Point3::new(0.0, 0.0, 2.4)),
                        normal: building().transform_vector3(Vec3::new(-s, 0.0, 1.0)),
                    },
                    agreement: true,
                }));
                body = g.boolean(body, roof, BooleanOperator::Difference);
            }
            let window = g.extrusion(circle(0.3), 1.0);
            let across = Transform3::from_translation(Vec3::new(-1.2, 0.5, 1.2))
                * Transform3::from_rotation_x(FRAC_PI_2);
            let window = g.place(window, building() * across);
            g.boolean(body, window, BooleanOperator::Difference)
        },
        &options,
    );
    let (bound, _) = assert_measured("gable wall, round window", case, &options, 2.0);
    assert!(bound <= budget, "{bound}");
}

#[test]
fn an_instanced_cut_wall_scales_the_measured_bound() {
    let budget = 1e-3;
    let options = options(budget);
    let (_, _, _, report, _) = compile(
        |g| {
            let wall = wall_with_round_window(g);
            g.place(wall, Transform3::from_translation(Vec3::new(1.0, 2.0, 3.0)))
        },
        &options,
    );
    let bound = report.bound.expect("certified");
    assert!(matches!(
        report.contributions[0].bound,
        DeviationBound::Certified(_)
    ));
    assert!(bound <= 2.0 * budget, "{bound}");
}

#[test]
fn a_union_the_exact_compiler_refuses_stays_unbounded_by_name() {
    let options = options(1e-3);
    let (_, _, _, report, _) = compile(
        |g| {
            let a = g.extrusion(rect(1.0, 1.0), 1.0);
            let b = g.extrusion(circle(0.4), 2.0);
            let b = g.place(b, Transform3::from_translation(Vec3::new(0.5, 0.0, -0.5)));
            g.boolean(a, b, BooleanOperator::Union)
        },
        &options,
    );
    assert_eq!(report.bound, None);
    let contribution = report.contributions[0];
    assert_eq!(contribution.path, DeviationPath::Boolean);
    assert_eq!(contribution.detail, "no exact result");
    assert_eq!(
        contribution.bound,
        DeviationBound::Unbounded("exact union or intersection of placed operands")
    );
}

/// A wall `6 x 0.25 x 3` minus a window `1.25 x 1.5` flush with both its
/// faces, both under `placement`, the window turned across the wall by
/// `across`.
fn wall_with_flush_window(g: &mut Graph, placement: Transform3, across: Transform3) -> NodeId {
    let wall = g.extrusion(rect(6.0, 0.25), 3.0);
    let wall = g.place(wall, placement);
    let window = g.extrusion(rect(1.25, 1.5), 0.25);
    let local = Transform3::from_translation(Vec3::new(0.75, 0.125, 1.25)) * across;
    let window = g.place(window, placement * local);
    g.boolean(wall, window, BooleanOperator::Difference)
}

/// The boolean contribution's detail for a flush window.
fn flush_window_detail(
    placement: Transform3,
    across: Transform3,
    tolerance: Tolerance,
) -> &'static str {
    let options = ExecutionOptions::new(tolerance)
        .with_chord_error(1e-3)
        .unwrap();
    let (_, _, _, report, _) = compile(|g| wall_with_flush_window(g, placement, across), &options);
    let contribution = report.contributions[0];
    assert_eq!(contribution.path, DeviationPath::Boolean);
    assert!(
        matches!(contribution.bound, DeviationBound::Certified(_)),
        "{report:?}"
    );
    contribution.detail
}

#[test]
fn the_detail_says_when_the_bound_holds_for_the_given_operands() {
    // Placed by exact axis matrices, the window's caps lie exactly in the
    // wall's faces: no decision within tolerance, at zero tolerance or a
    // positive one, so the bound is to the exact boolean of the operands
    // as given (#236).
    let exact_across =
        Transform3::from_mat3(axiolid_core::Mat3::from_cols(Vec3::X, Vec3::Z, -Vec3::Y));
    let grid = Transform3::from_mat3_translation(
        axiolid_core::Mat3::from_cols(Vec3::Y, -Vec3::X, Vec3::Z),
        Vec3::new(12.5, -4.0, 3.25),
    );
    for tolerance in [Tolerance::ZERO, Tolerance::METRE] {
        assert_eq!(
            flush_window_detail(grid, exact_across, tolerance),
            "measured against the exact compiler's result, the exact boolean of the given operands"
        );
    }
    // Turned from sines and cosines under a building placement, the caps
    // meet the faces only up to rounding: read within tolerance.
    assert_eq!(
        flush_window_detail(
            building(),
            Transform3::from_rotation_x(FRAC_PI_2),
            Tolerance::METRE
        ),
        "measured against the exact compiler's result, operands within tolerance"
    );
}
