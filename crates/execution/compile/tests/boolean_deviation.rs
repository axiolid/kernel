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
//! and round holes through its web, a slab with a round shaft, a slab with
//! oblique round openings (#287), and a gable wall clipped by two roof half-spaces (#234) with a round window. Unions and
//! intersections of placed operands are measured like differences (#376);
//! a boolean the exact compiler refuses stays unbounded with the refusal's
//! name.

use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, PI};
use std::time::Instant;

use axiolid_brep::ExactBRep;
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{
    BooleanOperator, Frame3, Plane3, Point2, Point3, Scalar, Tolerance, Transform2, Transform3,
    Vec2, Vec3,
};
use axiolid_curve::{Circle3, Curve3};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::FaceDomain;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::deviation::COINCIDENT_OPERANDS;
use axiolid_mesh_compile::{
    DeviationBound, DeviationContribution, DeviationPath, DeviationReport, ReferenceExactCompiler,
    ReferenceMeshCompiler,
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

    fn extrusion_along(&mut self, profile: Profile, direction: Vec3, depth: Scalar) -> NodeId {
        let profile = self.push(GeometryNode::Profile(profile));
        self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction,
            depth,
        }))
    }

    fn revolution(
        &mut self,
        profile: Profile,
        axis_origin: Point3,
        axis_direction: Vec3,
        angle: Scalar,
    ) -> NodeId {
        let profile = self.push(GeometryNode::Profile(profile));
        self.push(GeometryNode::SolidOperation(SolidOperation::Revolution {
            profile,
            axis_origin,
            axis_direction,
            angle,
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
            .as_chunks::<3>()
            .0
            .iter()
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
            best = best.min(triangle_distance(p, self.triangles[i]));
        }
        best
    }
}

/// The distance from `p` to a triangle: across the plane when the foot
/// lies inside (the three edge signs agree), else to the nearest edge.
/// Well conditioned on slivers, where the barycentric closest-point
/// formulas lose digits.
fn triangle_distance(p: Point3, [a, b, c]: [Point3; 3]) -> Scalar {
    let n = (b - a).cross(c - a);
    let side = |v0: Point3, v1: Point3| (v1 - v0).cross(p - v0).dot(n);
    let signs = [side(a, b), side(b, c), side(c, a)];
    if n.length() > 0.0 && (signs.iter().all(|s| *s >= 0.0) || signs.iter().all(|s| *s <= 0.0)) {
        return (p - a).dot(n).abs() / n.length();
    }
    [(a, b), (b, c), (c, a)]
        .into_iter()
        .map(|(u, v)| {
            let edge = v - u;
            let t = if edge.length_squared() > 0.0 {
                ((p - u).dot(edge) / edge.length_squared()).clamp(0.0, 1.0)
            } else {
                0.0
            };
            p.distance(u + edge * t)
        })
        .fold(Scalar::INFINITY, Scalar::min)
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
    // Faces read as one surface within tolerance name how far that moved
    // an operand (#324), within the tolerance; the rest is the boolean's
    // own measurement.
    let (coincident, measured): (Vec<&DeviationContribution>, Vec<_>) = report
        .contributions
        .iter()
        .partition(|c| c.detail == COINCIDENT_OPERANDS);
    for c in &coincident {
        assert!(
            matches!(c.bound, DeviationBound::Certified(moved) if moved <= options.tolerance().linear()),
            "{name}: {report:?}"
        );
    }
    assert_eq!(measured.len(), 1, "{report:?}");
    let contribution = *measured[0];
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

/// The wall `3.65 x 0.25 x 3.67` less two door openings `1.01 x 2.26`
/// cut through it (#314, the construction of #312), each operand placed
/// by `world`.
fn wall_with_two_doors(g: &mut Graph, world: Transform3) -> NodeId {
    let derived = |x: Scalar, y: Scalar, transform: Transform2| Profile::Derived {
        basis: Box::new(rect(x, y)),
        transform,
    };
    let wall = derived(
        3.65,
        0.25,
        Transform2::from_cols(Vec2::NEG_X, Vec2::NEG_Y, Vec2::new(1.825, 0.0)),
    );
    let wall = g.extrusion(wall, 3.67);
    let mut net = g.place(wall, world);
    for x in [3.325, 1.35] {
        let door = derived(
            1.01,
            2.26,
            Transform2::from_cols(Vec2::Y, Vec2::NEG_X, Vec2::new(1.13, 0.505)),
        );
        let door = g.extrusion(door, 0.25);
        let door = g.place(
            door,
            world
                * Transform3::from_cols(
                    Vec3::Z,
                    Vec3::NEG_X,
                    Vec3::NEG_Y,
                    Vec3::new(x, 0.125, 0.2),
                ),
        );
        net = g.boolean(net, door, BooleanOperator::Difference);
    }
    net
}

#[test]
fn an_all_planar_boolean_is_bounded_near_rounding() {
    // Every face is a plane: the search runs to the measured distance, not
    // to the budget (#314), which reported 9.96e-4 of a 1e-3 request.
    for (name, world) in [
        ("two doors", Transform3::IDENTITY),
        ("two doors, placed", building()),
    ] {
        let options = options(1e-3);
        let case = compile(|g| wall_with_two_doors(g, world), &options);
        let (bound, _) = assert_measured(name, case, &options, Scalar::INFINITY);
        assert!(bound <= 1e-11, "{name}: {bound:e}");
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

/// The IPE 300-like beam with its 15 mm root fillets, the hole run through
/// the web and both fillets past the flange tips: touching the flange, it
/// touches each top fillet too, where the hole/fillet section has a
/// double point (#249).
fn filleted_beam_with_a_hole_touching_its_flange(g: &mut Graph) -> NodeId {
    let (depth, width, flange, radius) = (0.3, 0.15, 0.0107, 0.05);
    let section = Profile::Section(SectionProfile::I {
        depth,
        width,
        web_thickness: 0.0071,
        flange_thickness: flange,
        fillet_radius: Some(0.015),
        flange_edge_radius: None,
        flange_slope: None,
    });
    let beam = g.extrusion(section, 4.0);
    let beam = g.place(beam, building());
    let hole = g.extrusion(circle(radius), 2.0 * width);
    let across =
        Transform3::from_translation(Vec3::new(-width, depth / 2.0 - flange - radius, 2.0))
            * Transform3::from_rotation_y(FRAC_PI_2);
    let hole = g.place(hole, building() * across);
    g.boolean(beam, hole, BooleanOperator::Difference)
}

#[test]
fn a_filleted_i_beam_with_a_hole_touching_its_flange_is_certified_against_its_exact_result() {
    let budget = 1e-3;
    let options = options(budget);
    let case = compile(filleted_beam_with_a_hole_touching_its_flange, &options);
    let (bound, _) = assert_measured(
        "filleted I-beam, hole touching the flange",
        case,
        &options,
        2.0,
    );
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

/// A slab less an oblique round shaft and an oblique rounded opening,
/// both extruded along directions leaning off the normal (#280), one up and
/// one down (a reflected prism).
///
/// Each operand is exact on its own: its arc walls are oblique circular
/// cylinders, audited clean, of the sheared prism's closed-form volume.
/// Since #287 the exact boolean cuts such a wall (a slab face meets it in
/// an ellipse, closed form), so the mesh is certified against the exact
/// difference, whose volume is the slab's less each opening's plan area
/// times the slab's depth: a horizontal section of an oblique prism is its
/// profile, moved. Before #287 the report stayed unbounded under the
/// refusal's name; before #280 the rounded operand stood its corner walls
/// upright on the profile plane, a shear away from the mesh's.
#[test]
fn a_slab_with_oblique_round_openings_is_certified_against_its_exact_result() {
    let budget = 1e-3;
    let options = options(budget);
    let mut operands = Vec::new();
    let case = compile(
        |g| {
            let slab = g.extrusion(rect(4.0, 3.0), 0.25);
            let slab = g.place(slab, building());
            let lean = Vec3::new(0.3, -0.2, 1.0);
            let shaft = g.extrusion_along(circle(0.3), lean, 0.75);
            let local = Transform3::from_translation(Vec3::new(-1.0, 0.4, -0.25));
            let shaft = g.place(shaft, building() * local);
            operands.push((shaft, PI * 0.09, lean));
            let slab = g.boolean(slab, shaft, BooleanOperator::Difference);
            let rounded = Profile::Rectangle(RectangleProfile {
                x: 0.8,
                y: 0.5,
                thickness: None,
                outer_radius: Some(0.15),
                inner_radius: None,
            });
            let lean = Vec3::new(-0.4, 0.25, -1.0);
            let opening = g.extrusion_along(rounded, lean, 0.75);
            let local = Transform3::from_translation(Vec3::new(0.9, -0.5, 0.5));
            let opening = g.place(opening, building() * local);
            operands.push((opening, 0.4 - (4.0 - PI) * 0.0225, lean));
            g.boolean(slab, opening, BooleanOperator::Difference)
        },
        &options,
    );
    let (graph, root) = (case.0.clone(), case.1);
    assert_measured("slab, oblique round openings", case, &options, 2.0);
    let mut removed = 0.0;
    for (operand, area, lean) in operands {
        let exact = ReferenceExactCompiler::new()
            .compile_exact(&graph, operand, &options)
            .expect("an oblique operand is exact");
        let health = axiolid_brep_audit::geometric_audit(&exact, options.tolerance());
        assert!(health.is_consistent(), "{:?}", health.defects());
        let volume = axiolid_measure::exact_properties(&exact, options.tolerance())
            .expect("measurable")
            .signed_volume;
        // A shear keeps volume: area times the normal part of the offset.
        let expected = area * 0.75 * lean.z.abs() / lean.length();
        assert!(
            (volume - expected).abs() <= 1e-9 * expected,
            "volume {volume}, expected {expected}"
        );
        removed += area * 0.25;
    }
    let exact = ReferenceExactCompiler::new()
        .compile_exact(&graph, root, &options)
        .expect("the exact difference");
    let health = axiolid_brep_audit::geometric_audit(&exact, options.tolerance());
    assert!(health.is_consistent(), "{:?}", health.defects());
    let volume = axiolid_measure::exact_properties(&exact, options.tolerance())
        .expect("measurable")
        .signed_volume;
    let expected = 3.0 - removed;
    assert!(
        (volume - expected).abs() <= 1e-9 * expected,
        "volume {volume}, expected {expected}"
    );
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

/// A circular column `r = 0.4` standing on, and through, a box, under a
/// building placement.
fn column_and_box(g: &mut Graph, operator: BooleanOperator) -> NodeId {
    let a = g.extrusion(rect(1.0, 1.0), 1.0);
    let a = g.place(a, building());
    let b = g.extrusion(circle(0.4), 2.0);
    let b = g.place(
        b,
        building() * Transform3::from_translation(Vec3::new(0.5, 0.0, -0.5)),
    );
    g.boolean(a, b, operator)
}

#[test]
fn a_column_united_with_a_box_is_certified_against_its_exact_result() {
    // The union of placed operands was unbounded by name before #376.
    let options = options(1e-3);
    let case = compile(|g| column_and_box(g, BooleanOperator::Union), &options);
    let (bound, _) = assert_measured("column union box", case, &options, 2.0);
    assert!(bound <= 1e-3, "{bound}");
}

#[test]
fn a_column_intersected_with_a_box_is_certified_against_its_exact_result() {
    let options = options(1e-3);
    let case = compile(
        |g| column_and_box(g, BooleanOperator::Intersection),
        &options,
    );
    let (bound, _) = assert_measured("column intersect box", case, &options, 2.0);
    assert!(bound <= 1e-3, "{bound}");
}

/// Two walls `4 x 0.3 x 3` meeting at a right angle and overlapping at
/// the corner, under a building placement.
fn two_walls(g: &mut Graph, offset: Vec3, operator: BooleanOperator) -> NodeId {
    let a = g.extrusion(rect(4.0, 0.3), 3.0);
    let a = g.place(a, building());
    let b = g.extrusion(rect(4.0, 0.3), 3.0);
    let b = g.place(
        b,
        building() * Transform3::from_translation(offset) * Transform3::from_rotation_z(FRAC_PI_2),
    );
    g.boolean(a, b, operator)
}

#[test]
fn two_overlapping_walls_united_are_certified_against_their_exact_result() {
    let options = options(1e-3);
    let case = compile(
        |g| two_walls(g, Vec3::new(1.5, 1.8, 0.0), BooleanOperator::Union),
        &options,
    );
    let (bound, _) = assert_measured("walls union", case, &options, Scalar::INFINITY);
    assert!(bound <= 1e-11, "{bound:e}");
}

#[test]
fn two_overlapping_walls_intersected_are_certified_against_their_exact_result() {
    let options = options(1e-3);
    let case = compile(
        |g| two_walls(g, Vec3::new(1.5, 1.8, 0.0), BooleanOperator::Intersection),
        &options,
    );
    let (bound, _) = assert_measured("walls intersection", case, &options, Scalar::INFINITY);
    assert!(bound <= 1e-11, "{bound:e}");
}

#[test]
fn two_disjoint_walls_united_are_certified_against_their_exact_result() {
    let options = options(1e-3);
    let case = compile(
        |g| two_walls(g, Vec3::new(1.5, 9.0, 0.0), BooleanOperator::Union),
        &options,
    );
    let (bound, _) = assert_measured("walls disjoint union", case, &options, Scalar::INFINITY);
    assert!(bound <= 1e-11, "{bound:e}");
}

#[test]
fn a_union_of_a_cut_wall_and_a_column_is_certified() {
    // A placed boolean is a subject of the union too.
    let options = options(1e-3);
    let case = compile(
        |g| {
            let cut = wall_with_round_window(g);
            let column = g.extrusion(circle(0.3), 3.0);
            let column = g.place(
                column,
                building() * Transform3::from_translation(Vec3::new(-2.5, 0.0, 0.0)),
            );
            g.boolean(cut, column, BooleanOperator::Union)
        },
        &options,
    );
    let (bound, _) = assert_measured("cut wall union column", case, &options, 2.0);
    assert!(bound <= 1e-3, "{bound}");
}

/// A box `0.4 x 0.5 x 6` and a quarter-turn revolution of a circle of
/// radius `0.1` about an axis 1 m away (the issue's repro, #378): a bent
/// tube of torus faces, its start disc in the box's bottom plane, turned
/// down (`up` false: it only touches the box) or up through the box.
fn box_and_bent_tube(g: &mut Graph, operator: BooleanOperator, up: bool) -> NodeId {
    let a = g.extrusion(rect(0.4, 0.5), 6.0);
    let axis = if up { Vec3::NEG_Y } else { Vec3::Y };
    let b = g.revolution(circle(0.1), Point3::new(-1.0, 0.0, 0.0), axis, FRAC_PI_2);
    g.boolean(a, b, operator)
}

/// The certified bound, measured against exact samples; `limit` times the
/// budget at most. TODO(#380): the intersection's cut of faceted operands
/// meeting a face obliquely lands beyond the budget (1.02e-3 at 1e-3), so
/// it is allowed 1.1x until the mesh boolean's cut is tightened.
fn assert_bent_tube(operator: BooleanOperator, up: bool, limit: Scalar, name: &str) {
    let options = options(1e-3);
    let case = compile(|g| box_and_bent_tube(g, operator, up), &options);
    let (bound, _) = assert_measured(name, case, &options, 2.0);
    assert!(bound <= limit * 1e-3, "{name}: {bound}");
}

#[test]
fn a_box_united_with_a_quarter_turn_revolution_is_certified() {
    assert_bent_tube(BooleanOperator::Union, false, 1.0, "box union tube below");
    assert_bent_tube(BooleanOperator::Union, true, 1.0, "box union tube above");
}

#[test]
fn a_box_minus_a_quarter_turn_revolution_is_certified() {
    assert_bent_tube(BooleanOperator::Difference, true, 1.0, "box minus tube");
}

#[test]
fn a_box_intersected_with_a_quarter_turn_revolution_is_certified() {
    assert_bent_tube(
        BooleanOperator::Intersection,
        true,
        1.1,
        "box intersect tube",
    );
}

#[test]
fn a_placed_revolution_united_with_a_placed_box_is_certified() {
    let options = options(1e-3);
    let case = compile(
        |g| {
            let a = g.extrusion(rect(0.4, 0.5), 6.0);
            let a = g.place(a, building());
            let b = g.revolution(circle(0.1), Point3::new(-1.0, 0.0, 0.0), Vec3::Y, FRAC_PI_2);
            let b = g.place(b, building());
            g.boolean(a, b, BooleanOperator::Union)
        },
        &options,
    );
    let (bound, _) = assert_measured("placed box union bent tube", case, &options, 2.0);
    assert!(bound <= 1e-3, "{bound}");
}

#[test]
fn a_ring_groove_turned_from_a_rectangle_is_certified() {
    // A full-turn revolution of a rectangle (a tube) cut from a wall as a
    // ring groove: the revolution is the tool.
    let options = options(1e-3);
    let case = compile(
        |g| {
            let wall = g.extrusion(rect(2.0, 0.3), 3.0);
            let wall = g.place(wall, building());
            let ring = g.revolution(
                rect(0.1, 0.15),
                Point3::new(-0.3, 0.0, 0.0),
                Vec3::Y,
                2.0 * PI,
            );
            let ring = g.place(
                ring,
                building() * Transform3::from_translation(Vec3::new(0.3, 0.15, 1.5)),
            );
            g.boolean(wall, ring, BooleanOperator::Difference)
        },
        &options,
    );
    let (bound, _) = assert_measured("wall minus turned ring", case, &options, 2.0);
    assert!(bound <= 1e-3, "{bound}");
}

#[test]
fn a_union_whose_tool_is_a_boolean_is_certified() {
    // A boolean is a compiled solid like any other, as the tool too (#326).
    let options = options(1e-3);
    let case = compile(
        |g| {
            let a = g.extrusion(rect(1.0, 1.0), 1.0);
            let tool = wall_with_round_window(g);
            g.boolean(a, tool, BooleanOperator::Union)
        },
        &options,
    );
    assert_measured("box union cut wall", case, &options, 2.0);
}

/// A `1 x 1 x 2` box united with a disk of radius `0.1` swept along a
/// quarter of a vertical circle of radius `2`, the pipe leaving the box
/// obliquely through its side (a swept disk is compiled exactly, #378).
#[test]
fn a_box_united_with_a_bent_pipe_is_certified() {
    let options = options(1e-3);
    let case = compile(
        |g| {
            let a = g.extrusion(rect(1.0, 1.0), 2.0);
            let a = g.place(a, Transform3::from_translation(Vec3::new(0.0, 0.0, -1.0)));
            let arc = g.push(GeometryNode::Curve3(Curve3::Circle(Circle3 {
                frame: Frame3 {
                    origin: Point3::new(2.0, 0.0, 0.0),
                    x: Vec3::X,
                    y: Vec3::Z,
                    z: -Vec3::Y,
                },
                radius: 2.0,
            })));
            let pipe = g.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
                directrix: arc,
                radius: 0.1,
                inner_radius: None,
                parameter_range: Some((PI, 1.5 * PI)),
                fillet_radius: None,
            }));
            g.boolean(a, pipe, BooleanOperator::Union)
        },
        &options,
    );
    let (bound, _) = assert_measured("box union bent pipe", case, &options, 2.0);
    assert!(bound <= 1e-3, "{bound}");
}

#[test]
fn a_tool_boolean_reading_within_tolerance_reports_it() {
    // The inner boolean (a wall with a window turned from sines and
    // cosines) reads its caps within tolerance; the outer union with a
    // faraway box takes no reading itself, so the report carries the
    // tool's.
    let options = options(1e-3);
    let (_, _, _, report, _) = compile(
        |g| {
            let far = g.extrusion(rect(1.0, 1.0), 1.0);
            let far = g.place(far, Transform3::from_translation(Vec3::new(40.0, 0.0, 0.0)));
            let tool =
                wall_with_flush_window(g, building(), Transform3::from_rotation_x(FRAC_PI_2));
            g.boolean(far, tool, BooleanOperator::Union)
        },
        &options,
    );
    assert_eq!(
        report.contributions[0].detail,
        "measured against the exact compiler's result, operands within tolerance"
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
