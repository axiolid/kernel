//! Chains of half-space intersections and the #194 pinch check (#328).
//!
//! A wall `[0, L] x [-T/2, T/2] x [0, H]`, placed by a rigid `world`
//! transform, is intersected by unbounded half-spaces in turn, each plane
//! through an edge an earlier cut left. The intersection of convex solids is
//! convex and cannot touch itself, and the mesh compiler does not refuse it:
//! every step compiles, closed and with no welded edge of more than two
//! faces, to its analytic volume.
//!
//! What openbimrs/ifc met (#328) is a chain over a NON-convex subject: a
//! wall a sloped clip had already notched, cut by a horizontal plane through
//! the notch's corner. That intersection is two solids touching along the
//! corner line, a true contact the check refuses by name. The same set cut
//! with the planes in another order never passes through the contact and
//! compiles to its analytic volume.

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Plane3, Point3, Tolerance, Transform3, Vec3};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{Profile, RectangleProfile};

/// The wall of openbimrs/ifc's report: length, thickness and height (m).
const L: f64 = 10.125;
const T: f64 = 0.35;
const H: f64 = 16.165033064627303;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

/// A rotation by `angle` about z, then a move by `(x, y)`.
fn turned(angle: f64, x: f64, y: f64) -> Transform3 {
    Transform3::from_cols(
        Vec3::new(angle.cos(), angle.sin(), 0.0),
        Vec3::new(-angle.sin(), angle.cos(), 0.0),
        Vec3::Z,
        Vec3::new(x, y, 0.0),
    )
}

/// Placements the chains are compiled under: none, a quarter turn (its
/// cosine a rounding residue), and turns placed away from the origin.
fn worlds() -> [Transform3; 4] {
    [
        Transform3::IDENTITY,
        turned(core::f64::consts::FRAC_PI_2, -8.071, -0.722),
        turned(0.959, 512.3, 259.7),
        turned(2.0, 1000.0, -2000.0),
    ]
}

struct Graph {
    builder: GeometryGraphBuilder,
    world: Transform3,
}

impl Graph {
    fn new(world: Transform3) -> Self {
        Self {
            builder: GeometryGraphBuilder::new(),
            world,
        }
    }

    fn push(&mut self, node: GeometryNode) -> NodeId {
        self.builder.push(node).expect("a valid node")
    }

    fn wall(&mut self) -> NodeId {
        let rect = self.push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: L,
                y: T,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )));
        let wall = self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: rect,
            direction: Vec3::Z,
            depth: H,
        }));
        let transform = self.world * Transform3::from_translation(Vec3::new(L / 2.0, 0.0, 0.0));
        self.push(GeometryNode::Instance(Instance {
            source: wall,
            transform,
        }))
    }

    /// `left` with `operator` and the half-space on the normal side of the
    /// plane through `origin` with `normal`, both in the wall's frame.
    fn clip(
        &mut self,
        left: NodeId,
        (origin, normal): (Point3, Vec3),
        operator: BooleanOperator,
    ) -> NodeId {
        let right = self.push(GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 {
                origin: self.world.transform_point3(origin),
                normal: self.world.transform_vector3(normal),
            },
            agreement: true,
        }));
        self.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator,
        }))
    }

    fn boolean(&mut self, left: NodeId, right: NodeId, operator: BooleanOperator) -> NodeId {
        self.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator,
        }))
    }

    fn finish(self, roots: Vec<NodeId>) -> GeometryGraph {
        self.builder.finish(roots).expect("a valid graph")
    }
}

fn compile(graph: &GeometryGraph, root: NodeId) -> Result<TriMesh, GeomError> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(graph, root, &options())
}

fn volume(mesh: &TriMesh) -> f64 {
    volume_properties(mesh, Tolerance::new(1e-7, 1e-9).unwrap())
        .expect("a closed solid")
        .signed_volume
}

/// Edges shared by more than two triangles once coincident positions are
/// welded, as a consumer that validates by position sees them.
fn welded_non_manifold(mesh: &TriMesh) -> usize {
    use std::collections::HashMap;
    let key = |p: Point3| {
        (
            (p.x + 0.0).to_bits(),
            (p.y + 0.0).to_bits(),
            (p.z + 0.0).to_bits(),
        )
    };
    let mut id: HashMap<(u64, u64, u64), u32> = HashMap::new();
    let welded: Vec<u32> = mesh
        .positions
        .iter()
        .map(|p| {
            let n = id.len() as u32;
            *id.entry(key(*p)).or_insert(n)
        })
        .collect();
    let mut edges: HashMap<(u32, u32), usize> = HashMap::new();
    for t in mesh.indices.chunks_exact(3) {
        for k in 0..3 {
            let (a, b) = (welded[t[k] as usize], welded[t[(k + 1) % 3] as usize]);
            *edges.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    edges.values().filter(|&&c| c > 2).count()
}

/// A line in the wall's xz plane: a point on it and its normal toward the
/// side kept.
type Line = ((f64, f64), (f64, f64));

/// The side of the line through `a` and `b` (in the wall's xz plane) that
/// holds the wall's point `(5, 5)`, as a plane containing the y axis.
fn side(a: (f64, f64), b: (f64, f64)) -> Line {
    let normal = (a.1 - b.1, b.0 - a.0);
    let inward = (5.0 - a.0) * normal.0 + (5.0 - a.1) * normal.1;
    (
        a,
        if inward < 0.0 {
            (-normal.0, -normal.1)
        } else {
            normal
        },
    )
}

fn plane(((x, z), (nx, nz)): ((f64, f64), (f64, f64))) -> (Point3, Vec3) {
    (
        Point3::new(x, -T / 2.0, z),
        Vec3::new(nx, 0.0, nz).normalize(),
    )
}

/// `polygon` clipped to the side of `line`, exactly as the planes cut the
/// wall's cross-section.
fn clip_polygon(
    polygon: &[(f64, f64)],
    ((ox, oz), (nx, nz)): ((f64, f64), (f64, f64)),
) -> Vec<(f64, f64)> {
    let at = |p: (f64, f64)| (p.0 - ox) * nx + (p.1 - oz) * nz;
    let mut out = Vec::new();
    for (i, &a) in polygon.iter().enumerate() {
        let b = polygon[(i + 1) % polygon.len()];
        let (sa, sb) = (at(a), at(b));
        if sa >= 0.0 {
            out.push(a);
        }
        if (sa > 0.0 && sb < 0.0) || (sa < 0.0 && sb > 0.0) {
            let t = sa / (sa - sb);
            out.push((a.0 + t * (b.0 - a.0), a.1 + t * (b.1 - a.1)));
        }
    }
    out
}

fn area(polygon: &[(f64, f64)]) -> f64 {
    let twice: f64 = (0..polygon.len())
        .map(|i| {
            let (a, b) = (polygon[i], polygon[(i + 1) % polygon.len()]);
            a.0 * b.1 - a.1 * b.0
        })
        .sum();
    twice / 2.0
}

/// Corner `(X0, Z0)` of the strip clip in openbimrs/ifc's wall.
const X0: f64 = 9.88645064727072;
const Z0: f64 = 11.404747431413968;

#[test]
fn convex_chains_through_earlier_cut_edges_compile_to_their_analytic_volume() {
    let top = |x: f64| (x, H);
    let cases: [(&str, Vec<Line>); 3] = [
        // A gable whose two roofs meet on the top face, then eaves through
        // the corners the roofs left on the ends.
        (
            "gable and eaves",
            vec![
                side((L, 12.0), top(5.0625)),
                side(top(5.0625), (0.0, 12.0)),
                side((L, 12.0), (9.0, 0.0)),
                side((0.0, 12.0), (1.0, 0.0)),
            ],
        ),
        // Three planes through the edge one clip left on the top face, the
        // last the third again from another of its points: flush with the
        // face it cut.
        (
            "a fan through one edge",
            vec![
                side((L, 10.0), top(7.0)),
                side(top(7.0), (0.0, 16.0)),
                side(top(7.0), (L, 9.0)),
                side((L, 9.0), top(7.0)),
            ],
        ),
        // The convex form of #328's corner: a steep plane, then the
        // horizontal and the vertical plane through the corner it passes.
        (
            "planes meeting on a sloped edge",
            vec![
                side((X0, Z0), (X0 - 0.475_693_878_096_870_3, Z0 + 1.7272497)),
                side((X0, Z0), (0.0, Z0)),
                side((X0, Z0), (X0, 0.0)),
            ],
        ),
    ];
    for (name, lines) in cases {
        for world in worlds() {
            let what = format!("{name} under {world:?}");
            let mut g = Graph::new(world);
            let mut node = g.wall();
            let mut steps = Vec::new();
            let mut polygon = vec![(0.0, 0.0), (L, 0.0), (L, H), (0.0, H)];
            for &line in &lines {
                node = g.clip(node, plane(line), BooleanOperator::Intersection);
                polygon = clip_polygon(&polygon, line);
                steps.push((node, area(&polygon) * T));
            }
            let graph = g.finish(steps.iter().map(|s| s.0).collect());
            for (k, (node, want)) in steps.into_iter().enumerate() {
                let mesh =
                    compile(&graph, node).unwrap_or_else(|e| panic!("{what}, step {k}: {e}"));
                assert_eq!(welded_non_manifold(&mesh), 0, "{what}, step {k}");
                let got = volume(&mesh);
                assert!(
                    (got - want).abs() <= 1e-9 * want,
                    "{what}, step {k}: {got} vs {want}"
                );
            }
        }
    }
}

/// A small deterministic generator, so the sweep is the same every run.
struct Lcg(u64);

impl Lcg {
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Chains of three to seven half-spaces, each plane through an edge (both
/// its ends, to rounding) or a corner of the mesh the chain has reached,
/// half of them vertical as a polygonal prism's walls are, under the
/// placements above: none is refused, and each step is a closed solid.
#[test]
fn chains_through_the_edges_they_leave_are_never_refused() {
    let mut rng = Lcg(328);
    for case in 0..48 {
        let world = worlds()[case % 4];
        let mut planes: Vec<(Point3, Vec3)> = Vec::new();
        for step in 0..3 + case % 5 {
            let mut g = Graph::new(world);
            let mut node = g.wall();
            for &p in &planes {
                node = g.clip(node, p, BooleanOperator::Intersection);
            }
            let graph = g.finish(vec![node]);
            let what = format!("case {case}, step {step}, planes {planes:?}, under {world:?}");
            let mesh = compile(&graph, node).unwrap_or_else(|e| panic!("{what}: {e}"));
            assert_eq!(welded_non_manifold(&mesh), 0, "{what}");
            assert!(volume(&mesh) > 0.0, "{what}");
            // The next plane, through an edge or a corner of this mesh, in
            // the wall's frame.
            let back = world.inverse();
            let triangle = (rng.unit() * (mesh.indices.len() / 3) as f64) as usize;
            let k = (rng.unit() * 3.0) as usize;
            let corner = |k: usize| {
                back.transform_point3(mesh.positions[mesh.indices[3 * triangle + k % 3] as usize])
            };
            let (a, b) = (corner(k), corner(k + 1));
            let centre = mesh
                .positions
                .iter()
                .fold(Vec3::ZERO, |s, p| s + back.transform_point3(*p))
                / mesh.positions.len() as f64;
            let r = Vec3::new(rng.unit() - 0.5, rng.unit() - 0.5, rng.unit() - 0.5);
            let mut normal = if rng.unit() < 0.8 {
                (b - a).cross(r)
            } else {
                r
            };
            if rng.unit() < 0.5 {
                normal.z = 0.0;
            }
            let normal = normal.normalize();
            planes.push((
                a,
                if normal.dot(centre - a) < 0.0 {
                    -normal
                } else {
                    normal
                },
            ));
        }
    }
}

/// openbimrs/ifc's wall (#328): a steep clip removes the wall above the
/// plane through `(X0, Z0)` over the strip `[X0 - W4, X0]`, leaving a notch
/// whose corner is that point; the strip clip then removes `[X0, L]` above
/// `Z0`. Both are lowered as openbimrs/ifc#518 does, `S - (S ∩ H ∩ W...)`,
/// `H` the clip's plane and `W` its prism's walls.
struct NotchedWall {
    graph: GeometryGraph,
    /// The notched wall cut by the strip clip's plane alone.
    through_the_corner: NodeId,
    /// The strip removed with the prism's walls cut first.
    walls_first: NodeId,
    /// Its closed-form volume.
    want: f64,
}

fn notched_wall() -> NotchedWall {
    // Clip 4 of the export, its normal as written.
    let n4 = Vec3::new(-0.964_105_691_204_325_9, 0.0, -0.265_518_768_051_204_1);
    let w4 = 0.475_693_878_096_870_3;
    let mut g = Graph::new(Transform3::IDENTITY);
    let wall = g.wall();
    // `walls` are the prism's sides over `[x0, x1]`, the wall's whole
    // thickness, all through the clip plane's origin height.
    let prism = |x0: f64, x1: f64| {
        [
            (Point3::new(x0, -T / 2.0, Z0), Vec3::Y),
            (Point3::new(x1, -T / 2.0, Z0), -Vec3::X),
            (Point3::new(x1, T / 2.0, Z0), -Vec3::Y),
            (Point3::new(x0, T / 2.0, Z0), Vec3::X),
        ]
    };
    let corner = Point3::new(X0, -T / 2.0, Z0);
    let mut tool = g.clip(wall, (corner, -n4), BooleanOperator::Intersection);
    for p in prism(X0 - w4, X0) {
        tool = g.clip(tool, p, BooleanOperator::Intersection);
    }
    let notched = g.boolean(wall, tool, BooleanOperator::Difference);
    let through_the_corner = g.clip(notched, (corner, Vec3::Z), BooleanOperator::Intersection);
    let mut strip = notched;
    for p in prism(X0, L) {
        strip = g.clip(strip, p, BooleanOperator::Intersection);
    }
    let strip = g.clip(strip, (corner, Vec3::Z), BooleanOperator::Intersection);
    let walls_first = g.boolean(notched, strip, BooleanOperator::Difference);
    // Clip 4 removes, over `[X0 - W4, X0]`, the wall above the line rising
    // from `(X0, Z0)` by `n4.x / n4.z` per unit leftwards; the strip clip
    // removes `[X0, L] x [Z0, H]`.
    let rise = n4.x / n4.z;
    let notch = w4 * (H - Z0) - rise * w4 * w4 / 2.0;
    let want = (L * H - notch - (L - X0) * (H - Z0)) * T;
    NotchedWall {
        graph: g.finish(vec![through_the_corner, walls_first]),
        through_the_corner,
        walls_first,
        want,
    }
}

#[test]
fn a_notched_wall_cut_through_its_notch_corner_is_refused_naming_the_contact() {
    // Above `Z0` the notched wall is the strip's column and the wall left
    // of the notch, which meet only along the notch's corner line: two
    // solids touching along an edge, which the check refuses by name.
    let wall = notched_wall();
    match compile(&wall.graph, wall.through_the_corner) {
        Err(GeomError::Degenerate(message)) => {
            assert!(
                message.contains("touches itself along the edge"),
                "{message}"
            );
            let ends = [
                format!("{:?}", Point3::new(X0, -T / 2.0, Z0)),
                format!("{:?}", Point3::new(X0, T / 2.0, Z0)),
            ];
            assert!(ends.iter().all(|end| message.contains(end)), "{message}");
        }
        Ok(mesh) => panic!(
            "a result touching itself was returned, {} welded edges with more than two faces",
            welded_non_manifold(&mesh)
        ),
        Err(other) => panic!("{other:?}"),
    }
}

#[test]
fn the_notched_wall_cut_walls_first_compiles_to_its_analytic_volume() {
    // The same strip, its prism's walls cut before its plane: no step
    // passes through the notch's corner alone.
    let wall = notched_wall();
    let mesh = compile(&wall.graph, wall.walls_first).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(welded_non_manifold(&mesh), 0);
    let got = volume(&mesh);
    assert!(
        (got - wall.want).abs() <= 1e-9 * wall.want,
        "{got} vs {}",
        wall.want
    );
}
