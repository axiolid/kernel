//! Curved B-rep faces whose edges carry only 3D curves (#310).
//!
//! Exchange formats such as IFC give an edge its 3D curve and no p-curve.
//! The tessellator inverts the edge's shared 3D samples onto the face's
//! surface instead. Each case compiles such a face with
//! `compile_mesh_with_deviation`, samples the exact surface and the exact
//! edges densely (in closed form or through the reference evaluator, never
//! read back from the compiler), and asserts that no sample lies further
//! from the mesh than the reported bound. The refusals pin the singular
//! cases by name.

use std::collections::HashMap;
use std::f64::consts::TAU;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Frame3, Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::{BSplineCurve3, Circle3, Curve3, Ellipse3, KnotSpec};
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationReport, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshClosure;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId};
use axiolid_surface::{BSplineSurface, Cone, Cylinder, Plane, Sphere, Surface};
use axiolid_topology::{
    BRep, Edge, EdgeId, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex, VertexId,
};

fn options(budget: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(budget.min(1e-3), 1e-9).unwrap())
        .with_chord_error(budget)
        .unwrap()
}

fn frame_at(z: Scalar) -> Frame3 {
    Frame3 {
        origin: Point3::new(0.0, 0.0, z),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

/// A B-rep of curved faces, each one ring of edge uses, none with a
/// p-curve. `use_edge` adds to the ring under construction; `close_face`
/// makes it a face, and `compile` closes the last one.
struct OneFace {
    builder: GeometryGraphBuilder,
    brep: BRep<NodeId>,
    uses: Vec<EdgeUse<NodeId>>,
    faces: Vec<axiolid_topology::FaceId>,
}

impl OneFace {
    fn new() -> Self {
        Self {
            builder: GeometryGraphBuilder::new(),
            brep: BRep::default(),
            uses: Vec::new(),
            faces: Vec::new(),
        }
    }

    fn close_face(&mut self, surface: &Surface) {
        self.close_oriented_face(surface, Orientation::Forward);
    }

    /// Close the ring as a face whose normal runs with (`Forward`) or
    /// against (`Reversed`) its surface's.
    fn close_oriented_face(&mut self, surface: &Surface, orientation: Orientation) {
        let surface = self
            .builder
            .push(GeometryNode::Surface(surface.clone()))
            .unwrap();
        let loop_id = self.brep.add_loop(Loop {
            edges: std::mem::take(&mut self.uses),
        });
        self.faces.push(self.brep.add_face(Face {
            surface: Some(surface),
            bounds: vec![FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation,
        }));
    }

    fn vertex(&mut self, position: Point3) -> VertexId {
        self.brep.add_vertex(Vertex { position })
    }

    fn edge(&mut self, start: VertexId, end: VertexId, curve: Option<Curve3>) -> EdgeId {
        let curve = curve.map(|c| self.builder.push(GeometryNode::Curve3(c)).unwrap());
        self.brep.add_edge(Edge { start, end, curve })
    }

    fn use_edge(&mut self, edge: EdgeId, orientation: Orientation) {
        self.uses.push(EdgeUse {
            edge,
            orientation,
            pcurve: None,
        });
    }

    fn compile(
        mut self,
        surface: &Surface,
        options: &ExecutionOptions,
    ) -> axiolid_contracts::GeomResult<(TriMesh, DeviationReport)> {
        if !self.uses.is_empty() {
            self.close_face(surface);
        }
        self.compile_shell(false, options)
            .map(|(mesh, _, report)| (mesh, report))
    }

    /// Compile the closed faces as one solid's shell, `closed` as declared.
    fn compile_shell(
        mut self,
        closed: bool,
        options: &ExecutionOptions,
    ) -> axiolid_contracts::GeomResult<(TriMesh, MeshClosure, DeviationReport)> {
        let shell = self.brep.add_shell(Shell {
            faces: self
                .faces
                .iter()
                .map(|&face| (face, Orientation::Forward))
                .collect(),
            closed,
        });
        self.brep.add_solid(Solid {
            outer: shell,
            voids: Vec::new(),
        });
        let root = self.builder.push(GeometryNode::BRep(self.brep)).unwrap();
        let graph = self.builder.finish(vec![root]).unwrap();
        let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
            .compile_mesh_with_deviation(&graph, root, options)?;
        Ok((outcome.mesh, outcome.closure, report))
    }
}

/// Triangles binned on a uniform grid, each in every cell its box touches
/// once grown by `reach`.
struct TriangleGrid {
    triangles: Vec<[Point3; 3]>,
    cells: HashMap<[i64; 3], Vec<usize>>,
    cell: Scalar,
    reach: Scalar,
}

fn key(p: Point3, cell: Scalar) -> [i64; 3] {
    [p.x, p.y, p.z].map(|v| (v / cell).floor() as i64)
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

/// Assert no sample is further from the mesh than the reported bound and
/// that the bound is within `ceiling` times the measured maximum.
fn assert_bounded(
    name: &str,
    mesh: &TriMesh,
    samples: &[Point3],
    report: &DeviationReport,
    ceiling: Scalar,
) {
    let bound = report
        .bound
        .unwrap_or_else(|| panic!("{name}: expected a certified bound, got {report:?}"));
    let grid = TriangleGrid::new(mesh, 2.0 * bound);
    let worst = samples
        .iter()
        .map(|p| grid.distance(*p))
        .fold(0.0, Scalar::max);
    let ratio = bound / worst;
    eprintln!(
        "{name}: reported {bound:e}, measured {worst:e}, ratio {ratio:.2}, {} triangles, {} samples",
        mesh.indices.len() / 3,
        samples.len()
    );
    assert!(
        worst <= bound,
        "{name}: a surface point lies {worst:e} from the mesh, above the reported {bound:e}"
    );
    assert!(
        ratio <= ceiling,
        "{name}: bound {bound:e} is {ratio:.1}x the measured {worst:e}"
    );
}

fn span(a: Scalar, b: Scalar, n: usize) -> impl Iterator<Item = Scalar> + Clone {
    (0..=n).map(move |i| a + (b - a) * (i as Scalar) / n as Scalar)
}

fn surface_samples(
    surface: &Surface,
    u: (Scalar, Scalar),
    v: (Scalar, Scalar),
    n: usize,
) -> Vec<Point3> {
    let mut out = Vec::with_capacity((n + 1) * (n + 1));
    for uu in span(u.0, u.1, n) {
        for vv in span(v.0, v.1, n) {
            out.push(axiolid_reference::surface::evaluate(surface, uu, vv).unwrap());
        }
    }
    out
}

fn curve_samples(curve: &Curve3, n: usize) -> Vec<Point3> {
    let domain = axiolid_reference::curve::domain3(curve);
    span(domain.start, domain.end, n)
        .map(|t| axiolid_reference::curve::evaluate3(curve, t).unwrap())
        .collect()
}

fn circle(z: Scalar, radius: Scalar) -> Curve3 {
    Curve3::Circle(Circle3 {
        frame: frame_at(z),
        radius,
    })
}

/// A whole cylinder band the way IFC writes one: a bottom and a top circle
/// and one straight seam, used once each way, none with a p-curve. The seam
/// is inverted at `u = 0` both times and unwrapped to `2 pi` on its way up.
#[test]
fn a_cylinder_band_with_only_3d_edges_compiles_within_its_reported_bound() {
    let (radius, height) = (0.3, 0.5);
    let surface = Surface::Cylinder(Cylinder {
        frame: frame_at(0.0),
        radius,
    });
    for budget in [1e-3, 1e-4] {
        let mut face = OneFace::new();
        let bottom = face.vertex(Point3::new(radius, 0.0, 0.0));
        let top = face.vertex(Point3::new(radius, 0.0, height));
        let rim0 = face.edge(bottom, bottom, Some(circle(0.0, radius)));
        let rim1 = face.edge(top, top, Some(circle(height, radius)));
        let seam = face.edge(bottom, top, None);
        face.use_edge(rim0, Orientation::Forward);
        face.use_edge(seam, Orientation::Forward);
        face.use_edge(rim1, Orientation::Reversed);
        face.use_edge(seam, Orientation::Reversed);
        let (mesh, report) = face
            .compile(&surface, &options(budget))
            .expect("a cylinder face with only 3D edges compiles");

        // Watertight along the seam: each of its samples is one vertex,
        // shared by both of its uses.
        let on_seam: Vec<_> = mesh
            .positions
            .iter()
            .filter(|p| (p.x - radius).abs() < 1e-12 && p.y.abs() < 1e-12)
            .collect();
        let mut heights: Vec<Scalar> = on_seam.iter().map(|p| p.z).collect();
        heights.sort_by(Scalar::total_cmp);
        heights.dedup();
        assert_eq!(heights.len(), on_seam.len(), "a seam sample is duplicated");
        assert!(on_seam.len() >= 2);

        let mut samples = surface_samples(&surface, (0.0, TAU), (0.0, height), 400);
        samples.extend(curve_samples(&circle(0.0, radius), 20_000));
        samples.extend(curve_samples(&circle(height, radius), 20_000));
        assert_bounded(
            &format!("cylinder band at {budget:e}"),
            &mesh,
            &samples,
            &report,
            4.0,
        );
    }
}

/// Half a circle of `radius` at height `z`, from angle `start` (0 or pi),
/// as two rational quadratic quarter arcs.
fn half_circle(z: Scalar, radius: Scalar, start: Scalar) -> Curve3 {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let s = if start == 0.0 { 1.0 } else { -1.0 };
    let p = |x: Scalar, y: Scalar| Point3::new(s * x * radius, s * y * radius, z);
    Curve3::BSpline(BSplineCurve3 {
        degree: 2,
        control_points: vec![
            p(1.0, 0.0),
            p(1.0, 1.0),
            p(0.0, 1.0),
            p(-1.0, 1.0),
            p(-1.0, 0.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![3, 2, 3],
        weights: Some(vec![1.0, w, 1.0, w, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// A tube of two half-cylinder faces sharing their two straight rulings,
/// no edge with a p-curve: the second face inverts the rulings the first
/// sampled, so the tube is closed along them by index, open only at its
/// rims.
#[test]
fn two_faces_share_the_samples_of_an_edge_without_a_pcurve() {
    let (radius, height) = (0.3, 0.5);
    let surface = Surface::Cylinder(Cylinder {
        frame: frame_at(0.0),
        radius,
    });
    let mut faces = OneFace::new();
    let a = faces.vertex(Point3::new(radius, 0.0, 0.0));
    let b = faces.vertex(Point3::new(-radius, 0.0, 0.0));
    let c = faces.vertex(Point3::new(-radius, 0.0, height));
    let d = faces.vertex(Point3::new(radius, 0.0, height));
    let front_low = faces.edge(a, b, Some(half_circle(0.0, radius, 0.0)));
    let front_high = faces.edge(d, c, Some(half_circle(height, radius, 0.0)));
    let back_low = faces.edge(b, a, Some(half_circle(0.0, radius, TAU / 2.0)));
    let back_high = faces.edge(c, d, Some(half_circle(height, radius, TAU / 2.0)));
    let ruling0 = faces.edge(a, d, None);
    let ruling1 = faces.edge(b, c, None);
    faces.use_edge(front_low, Orientation::Forward);
    faces.use_edge(ruling1, Orientation::Forward);
    faces.use_edge(front_high, Orientation::Reversed);
    faces.use_edge(ruling0, Orientation::Reversed);
    faces.close_face(&surface);
    faces.use_edge(back_low, Orientation::Forward);
    faces.use_edge(ruling0, Orientation::Forward);
    faces.use_edge(back_high, Orientation::Reversed);
    faces.use_edge(ruling1, Orientation::Reversed);
    let (mesh, report) = faces
        .compile(&surface, &options(1e-4))
        .expect("two cylinder faces with only 3D edges compile");

    // Every mesh edge off the rims is used by exactly two triangles.
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for t in mesh.indices.chunks_exact(3) {
        for (x, y) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *uses.entry((x.min(y), x.max(y))).or_default() += 1;
        }
    }
    let on_rim = |i: u32| {
        let z = mesh.positions[i as usize].z;
        z.abs() < 1e-12 || (z - height).abs() < 1e-12
    };
    for (&(x, y), &count) in &uses {
        if on_rim(x) && on_rim(y) {
            assert!(count <= 2, "rim edge used {count} times");
        } else {
            assert_eq!(
                count, 2,
                "an edge off the rims is used {count} times: the tube is open"
            );
        }
    }

    let mut samples = surface_samples(&surface, (0.0, TAU), (0.0, height), 400);
    for curve in [
        half_circle(0.0, radius, 0.0),
        half_circle(height, radius, 0.0),
        half_circle(0.0, radius, TAU / 2.0),
        half_circle(height, radius, TAU / 2.0),
    ] {
        samples.extend(curve_samples(&curve, 10_000));
    }
    assert_bounded("half-cylinder pair", &mesh, &samples, &report, 4.0);
}

/// A quarter of a cylinder between two rational quarter arcs and two
/// straight rulings: a patch that does not wrap.
#[test]
fn a_cylinder_patch_bounded_by_rational_arcs_is_within_its_reported_bound() {
    let (radius, height) = (0.3, 0.5);
    let surface = Surface::Cylinder(Cylinder {
        frame: frame_at(0.0),
        radius,
    });
    // The exact quarter circle as a rational quadratic.
    let arc = |z: Scalar| {
        Curve3::BSpline(BSplineCurve3 {
            degree: 2,
            control_points: vec![
                Point3::new(radius, 0.0, z),
                Point3::new(radius, radius, z),
                Point3::new(0.0, radius, z),
            ],
            knots: vec![0.0, 1.0],
            multiplicities: vec![3, 3],
            weights: Some(vec![1.0, std::f64::consts::FRAC_1_SQRT_2, 1.0]),
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        })
    };
    let mut face = OneFace::new();
    let corners = [
        Point3::new(radius, 0.0, 0.0),
        Point3::new(0.0, radius, 0.0),
        Point3::new(0.0, radius, height),
        Point3::new(radius, 0.0, height),
    ];
    let v: Vec<_> = corners.iter().map(|&p| face.vertex(p)).collect();
    let low = face.edge(v[0], v[1], Some(arc(0.0)));
    let left = face.edge(v[1], v[2], None);
    let high = face.edge(v[3], v[2], Some(arc(height)));
    let right = face.edge(v[0], v[3], None);
    face.use_edge(low, Orientation::Forward);
    face.use_edge(left, Orientation::Forward);
    face.use_edge(high, Orientation::Reversed);
    face.use_edge(right, Orientation::Reversed);
    let (mesh, report) = face
        .compile(&surface, &options(1e-4))
        .expect("a cylinder patch with only 3D edges compiles");
    let mut samples = surface_samples(&surface, (0.0, TAU / 4.0), (0.0, height), 400);
    samples.extend(curve_samples(&arc(0.0), 20_000));
    samples.extend(curve_samples(&arc(height), 20_000));
    assert_bounded("cylinder patch", &mesh, &samples, &report, 4.0);
}

/// The bicubic bump of the deviation report, its four boundary edges given
/// as the 3D B-spline curves they are (the net's border rows), with no
/// p-curves: each sample is inverted by the certified global search.
#[test]
fn a_bspline_face_with_only_3d_edges_compiles_within_its_reported_bound() {
    let net: Vec<Vec<Point3>> = (0..5)
        .map(|i| {
            (0..5)
                .map(|j| {
                    let (x, y) = (i as Scalar * 0.25, j as Scalar * 0.25);
                    let bump = if (1..4).contains(&i) && (1..4).contains(&j) {
                        0.12
                    } else {
                        0.0
                    };
                    Point3::new(x, y, bump + 0.05 * x * y)
                })
                .collect()
        })
        .collect();
    let surface = Surface::BSpline(BSplineSurface {
        u_degree: 3,
        v_degree: 3,
        control_points: net.clone(),
        u_knots: vec![0.0, 0.5, 1.0],
        u_multiplicities: vec![4, 1, 4],
        v_knots: vec![0.0, 0.5, 1.0],
        v_multiplicities: vec![4, 1, 4],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    });
    let border = |points: Vec<Point3>| {
        Curve3::BSpline(BSplineCurve3 {
            degree: 3,
            control_points: points,
            knots: vec![0.0, 0.5, 1.0],
            multiplicities: vec![4, 1, 4],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        })
    };
    // Along u at v = 0 and v = 1, along v at u = 0 and u = 1.
    let v0 = border((0..5).map(|i| net[i][0]).collect());
    let v1 = border((0..5).map(|i| net[i][4]).collect());
    let u0 = border(net[0].clone());
    let u1 = border(net[4].clone());

    let mut face = OneFace::new();
    let c00 = face.vertex(net[0][0]);
    let c10 = face.vertex(net[4][0]);
    let c11 = face.vertex(net[4][4]);
    let c01 = face.vertex(net[0][4]);
    let e0 = face.edge(c00, c10, Some(v0.clone()));
    let e1 = face.edge(c10, c11, Some(u1.clone()));
    let e2 = face.edge(c01, c11, Some(v1.clone()));
    let e3 = face.edge(c00, c01, Some(u0.clone()));
    face.use_edge(e0, Orientation::Forward);
    face.use_edge(e1, Orientation::Forward);
    face.use_edge(e2, Orientation::Reversed);
    face.use_edge(e3, Orientation::Reversed);
    let (mesh, report) = face
        .compile(&surface, &options(1e-3))
        .expect("a B-spline face with only 3D edges compiles");
    let mut samples = surface_samples(&surface, (0.0, 1.0), (0.0, 1.0), 400);
    for curve in [&v0, &v1, &u0, &u1] {
        samples.extend(curve_samples(curve, 20_000));
    }
    assert_bounded("B-spline face", &mesh, &samples, &report, 5.0);
}

/// A whole cone up to its apex, bounded by its base circle and one ruling
/// used once each way. At the apex every `u` names the point, so the
/// ruling's parameters there are not determined: refused by name, never
/// guessed.
#[test]
fn an_edge_through_a_cone_apex_is_refused() {
    let (radius, semi_angle) = (0.4, -0.4_f64);
    let surface = Surface::Cone(Cone {
        frame: frame_at(0.0),
        radius,
        semi_angle,
    });
    let apex_height = -radius / semi_angle.tan();
    let mut face = OneFace::new();
    let base = face.vertex(Point3::new(radius, 0.0, 0.0));
    let apex = face.vertex(Point3::new(0.0, 0.0, apex_height));
    let rim = face.edge(base, base, Some(circle(0.0, radius)));
    let ruling = face.edge(base, apex, None);
    face.use_edge(rim, Orientation::Forward);
    face.use_edge(ruling, Orientation::Forward);
    face.use_edge(ruling, Orientation::Reversed);
    let error = face
        .compile(&surface, &options(1e-3))
        .expect_err("an edge into the apex has no parameters there");
    assert!(
        matches!(&error, GeomError::Degenerate(reason) if reason.contains("cone apex")),
        "got {error:?}"
    );
}

/// A meridian from the equator over the north pole: refused by name.
#[test]
fn an_edge_through_a_sphere_pole_is_refused() {
    let radius = 0.5;
    let surface = Surface::Sphere(Sphere {
        frame: frame_at(0.0),
        radius,
    });
    // Half a great circle in the x-z plane, from (r, 0, 0) over the pole
    // to (-r, 0, 0), as a rational quadratic pair of quarter arcs.
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let meridian = Curve3::BSpline(BSplineCurve3 {
        degree: 2,
        control_points: vec![
            Point3::new(radius, 0.0, 0.0),
            Point3::new(radius, 0.0, radius),
            Point3::new(0.0, 0.0, radius),
            Point3::new(-radius, 0.0, radius),
            Point3::new(-radius, 0.0, 0.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![3, 2, 3],
        weights: Some(vec![1.0, w, 1.0, w, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    let equator = Curve3::BSpline(BSplineCurve3 {
        degree: 2,
        control_points: vec![
            Point3::new(-radius, 0.0, 0.0),
            Point3::new(-radius, radius, 0.0),
            Point3::new(0.0, radius, 0.0),
            Point3::new(radius, radius, 0.0),
            Point3::new(radius, 0.0, 0.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![3, 2, 3],
        weights: Some(vec![1.0, w, 1.0, w, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    let mut face = OneFace::new();
    let east = face.vertex(Point3::new(radius, 0.0, 0.0));
    let west = face.vertex(Point3::new(-radius, 0.0, 0.0));
    let over = face.edge(east, west, Some(meridian));
    let back = face.edge(west, east, Some(equator));
    face.use_edge(over, Orientation::Forward);
    face.use_edge(back, Orientation::Forward);
    let error = face
        .compile(&surface, &options(1e-3))
        .expect_err("an edge over the pole has no parameters there");
    assert!(
        matches!(&error, GeomError::Degenerate(reason) if reason.contains("sphere pole")),
        "got {error:?}"
    );
}

/// An edge whose 3D curve leaves the face's surface is refused, not
/// projected onto it.
#[test]
fn an_edge_off_its_surface_is_refused() {
    let (radius, height) = (0.3, 0.5);
    let surface = Surface::Cylinder(Cylinder {
        frame: frame_at(0.0),
        radius,
    });
    let mut face = OneFace::new();
    let bottom = face.vertex(Point3::new(radius, 0.0, 0.0));
    let top = face.vertex(Point3::new(radius, 0.0, height));
    // The top rim is a little too wide for the cylinder.
    let rim0 = face.edge(bottom, bottom, Some(circle(0.0, radius)));
    let wide = face.vertex(Point3::new(radius * 1.01, 0.0, height));
    let rim1 = face.edge(wide, wide, Some(circle(height, radius * 1.01)));
    let seam = face.edge(bottom, top, None);
    let step = face.edge(top, wide, None);
    face.use_edge(rim0, Orientation::Forward);
    face.use_edge(seam, Orientation::Forward);
    face.use_edge(step, Orientation::Forward);
    face.use_edge(rim1, Orientation::Reversed);
    face.use_edge(step, Orientation::Reversed);
    face.use_edge(seam, Orientation::Reversed);
    let error = face
        .compile(&surface, &options(1e-3))
        .expect_err("an edge off the surface is refused");
    assert!(
        matches!(&error, GeomError::Degenerate(reason) if reason.contains("from the surface")),
        "got {error:?}"
    );
}

/// A closed curve is trimmed at its vertex (#319): a circle that starts a
/// quarter turn away from the edge's single vertex is the whole turn from
/// that vertex, not refused.
#[test]
fn a_full_circle_is_the_whole_turn_from_its_one_vertex() {
    let (radius, height) = (0.3, 0.5);
    let surface = Surface::Cylinder(Cylinder {
        frame: frame_at(0.0),
        radius,
    });
    let turned = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::Y,
            y: -Vec3::X,
            z: Vec3::Z,
        },
        radius,
    });
    let mut face = OneFace::new();
    let bottom = face.vertex(Point3::new(radius, 0.0, 0.0));
    let top = face.vertex(Point3::new(radius, 0.0, height));
    let rim0 = face.edge(bottom, bottom, Some(turned.clone()));
    let rim1 = face.edge(top, top, Some(circle(height, radius)));
    let seam = face.edge(bottom, top, None);
    face.use_edge(rim0, Orientation::Forward);
    face.use_edge(seam, Orientation::Forward);
    face.use_edge(rim1, Orientation::Reversed);
    face.use_edge(seam, Orientation::Reversed);
    let (mesh, report) = face
        .compile(&surface, &options(1e-4))
        .expect("a circle trimmed at its one vertex compiles");
    let mut samples = surface_samples(&surface, (0.0, TAU), (0.0, height), 400);
    samples.extend(curve_samples(&turned, 20_000));
    samples.extend(curve_samples(&circle(height, radius), 20_000));
    assert_bounded("circle from its vertex", &mesh, &samples, &report, 4.0);
}

/// The IFC cylinder of openbimrs/ifc's `advanced_brep_3d_edges` fixture
/// with its seam an `IfcLine` whose vector is a unit length (#319): the
/// line's parameter runs past `[0, 1]` to the top vertex at `t = 2`, and
/// the seam is the segment between its vertices' parameters.
#[test]
fn a_seam_line_with_a_unit_direction_is_trimmed_at_its_vertices() {
    let (radius, height) = (0.5, 2.0);
    let line = Curve3::Line(axiolid_curve::Line3 {
        origin: Point3::new(radius, 0.0, 0.0),
        direction: Vec3::Z,
    });
    for budget in [1e-3, 1e-4] {
        let (mesh, closure, report) = capped_cylinder(
            radius,
            height,
            frame_at(0.0),
            circle(0.0, radius),
            Some(line.clone()),
            budget,
            false,
        );
        let name = format!("unit-direction seam at {budget:e}");
        assert_eq!(closure, MeshClosure::Solid, "{name}");
        assert_watertight(&name, &mesh);
        let top = mesh.positions.iter().map(|p| p.z).fold(0.0, Scalar::max);
        assert!(top <= height + 1e-12, "{name}: the seam runs to {top}");
        let area = TAU * radius * height + TAU * radius * radius;
        assert_volume(
            &name,
            &mesh,
            &report,
            0.5 * TAU * radius * radius * height,
            area,
        );
        let wall = Surface::Cylinder(Cylinder {
            frame: frame_at(0.0),
            radius,
        });
        let mut samples = surface_samples(&wall, (0.0, TAU), (0.0, height), 400);
        samples.extend(span(0.0, height, 4_000).map(|z| Point3::new(radius, 0.0, z)));
        assert_bounded(&name, &mesh, &samples, &report, 6.0);
    }
}

/// A quarter of a cylinder whose arcs are whole circles trimmed at the
/// patch's corners, its rulings lines with unit directions, none with a
/// p-curve. `clockwise` gives the arcs circles running the other way about
/// the axis, so each arc's edge runs from the quarter's far corner back to
/// its near one and the patch uses it reversed.
fn quarter_patch(clockwise: bool) -> (TriMesh, DeviationReport, Vec<Point3>) {
    let (radius, height) = (0.3, 0.5);
    let surface = Surface::Cylinder(Cylinder {
        frame: frame_at(0.0),
        radius,
    });
    let arc = |z: Scalar| {
        let frame = if clockwise {
            Frame3 {
                origin: Point3::new(0.0, 0.0, z),
                x: Vec3::X,
                y: -Vec3::Y,
                z: -Vec3::Z,
            }
        } else {
            frame_at(z)
        };
        Curve3::Circle(Circle3 { frame, radius })
    };
    let ruling = |x: Scalar, y: Scalar| {
        Curve3::Line(axiolid_curve::Line3 {
            origin: Point3::new(x, y, 0.0),
            direction: Vec3::Z,
        })
    };
    let mut face = OneFace::new();
    let corners = [
        Point3::new(radius, 0.0, 0.0),
        Point3::new(0.0, radius, 0.0),
        Point3::new(0.0, radius, height),
        Point3::new(radius, 0.0, height),
    ];
    let v: Vec<_> = corners.iter().map(|&p| face.vertex(p)).collect();
    let (low, high, sense) = if clockwise {
        (
            face.edge(v[1], v[0], Some(arc(0.0))),
            face.edge(v[2], v[3], Some(arc(height))),
            [Orientation::Reversed, Orientation::Forward],
        )
    } else {
        (
            face.edge(v[0], v[1], Some(arc(0.0))),
            face.edge(v[3], v[2], Some(arc(height))),
            [Orientation::Forward, Orientation::Reversed],
        )
    };
    let left = face.edge(v[1], v[2], Some(ruling(0.0, radius)));
    let right = face.edge(v[0], v[3], Some(ruling(radius, 0.0)));
    face.use_edge(low, sense[0]);
    face.use_edge(left, Orientation::Forward);
    face.use_edge(high, sense[1]);
    face.use_edge(right, Orientation::Reversed);
    let (mesh, report) = face
        .compile(&surface, &options(1e-4))
        .expect("a quarter patch bounded by trimmed circles compiles");
    let mut samples = surface_samples(&surface, (0.0, TAU / 4.0), (0.0, height), 400);
    for z in [0.0, height] {
        samples.extend(
            span(0.0, TAU / 4.0, 20_000)
                .map(|t| Point3::new(radius * t.cos(), radius * t.sin(), z)),
        );
    }
    (mesh, report, samples)
}

/// Every mesh vertex of a quarter patch lies in its quarter: the arcs were
/// trimmed at the corners, not sampled whole.
fn assert_in_quarter(name: &str, mesh: &TriMesh) {
    for p in &mesh.positions {
        assert!(
            p.x >= -1e-12 && p.y >= -1e-12,
            "{name}: a mesh vertex {p:?} lies outside the quarter"
        );
    }
}

#[test]
fn a_quarter_circle_edge_is_trimmed_at_its_vertices() {
    let (mesh, report, samples) = quarter_patch(false);
    assert_in_quarter("quarter circle", &mesh);
    assert_bounded("quarter circle", &mesh, &samples, &report, 4.0);
}

/// The same quarter on circles running clockwise: each arc runs forward
/// along its own circle from its start vertex, and the patch uses it
/// reversed.
#[test]
fn a_reversed_sense_arc_is_trimmed_at_its_vertices() {
    let (mesh, report, samples) = quarter_patch(true);
    assert_in_quarter("reversed-sense arc", &mesh);
    assert_bounded("reversed-sense arc", &mesh, &samples, &report, 4.0);
}

/// The refusal of a cylinder face whose bottom rim is the edge `curve`
/// (the bottom circle by default) from `(r, 0, 0)` to a second vertex at
/// `end`, closed back to its start by a straight edge.
fn refused_rim(end: Point3, curve: Option<Curve3>) -> GeomError {
    let (radius, height) = (0.3, 0.5);
    let surface = Surface::Cylinder(Cylinder {
        frame: frame_at(0.0),
        radius,
    });
    let mut face = OneFace::new();
    let bottom = face.vertex(Point3::new(radius, 0.0, 0.0));
    let other = face.vertex(end);
    let top = face.vertex(Point3::new(radius, 0.0, height));
    let curve = curve.unwrap_or_else(|| circle(0.0, radius));
    let rim0 = face.edge(bottom, other, Some(curve));
    let back = face.edge(other, bottom, None);
    let rim1 = face.edge(top, top, Some(circle(height, radius)));
    let seam = face.edge(bottom, top, None);
    face.use_edge(rim0, Orientation::Forward);
    face.use_edge(back, Orientation::Forward);
    face.use_edge(seam, Orientation::Forward);
    face.use_edge(rim1, Orientation::Reversed);
    face.use_edge(seam, Orientation::Reversed);
    face.compile(&surface, &options(1e-3))
        .expect_err("the rim is refused")
}

/// A vertex a little off its edge's circle is not projected onto it.
#[test]
fn a_vertex_off_its_edge_curve_is_refused_by_name() {
    let error = refused_rim(Point3::new(0.0, 0.303, 0.0), None);
    assert!(
        matches!(&error, GeomError::InvalidInput(reason)
            if reason.contains("end vertex of an edge without a p-curve is not on its curve")),
        "got {error:?}"
    );
}

/// Two vertices at one point of a circle bound either no arc or the whole
/// turn: refused as ambiguous rather than read as either.
#[test]
fn two_vertices_at_one_point_of_a_circle_are_refused_as_ambiguous() {
    let error = refused_rim(Point3::new(0.3, 0.0, 0.0), None);
    assert!(
        matches!(&error, GeomError::InvalidInput(reason) if reason.contains("ambiguous")),
        "got {error:?}"
    );
}

/// A vertex off a B-spline edge curve is refused by name (#333): the
/// certified projection's distance lower bound exceeds the tolerance.
#[test]
fn a_vertex_off_its_bspline_edge_curve_is_refused_by_name() {
    let error = refused_rim(
        Point3::new(0.0, 0.303, 0.0),
        Some(half_circle(0.0, 0.3, 0.0)),
    );
    assert!(
        matches!(&error, GeomError::InvalidInput(reason)
            if reason.contains("end vertex of an edge without a p-curve is not on its curve")),
        "got {error:?}"
    );
}

/// A section curve's inversion takes the first of several turns naming a
/// point, so an edge on one that does not run its whole domain stays
/// refused by name (#333): here the bottom rim as the section of the
/// cylinder at height zero.
#[test]
fn a_section_edge_not_ending_at_its_vertices_is_refused_by_name() {
    let section = Curve3::RuledSection(axiolid_curve::RuledSection3 {
        carrier: axiolid_curve::RuledCarrier {
            frame: frame_at(0.0),
            x_radius: 0.3,
            y_radius: 0.3,
            slope: 0.0,
        },
        graph: axiolid_curve::QuadraticGraph2 {
            a: axiolid_curve::Trig2 {
                constant: 1.0,
                ..axiolid_curve::Trig2::default()
            },
            b: axiolid_curve::Trig2::default(),
            c: axiolid_curve::Trig2::default(),
            branch: axiolid_curve::Branch::Plus,
        },
    });
    let error = refused_rim(Point3::new(0.0, 0.3, 0.0), Some(section));
    assert!(
        matches!(&error, GeomError::InvalidInput(reason)
            if reason.contains("neither runs its whole domain")),
        "got {error:?}"
    );
}

/// A whole circle as one closed rational B-spline of four quarter arcs,
/// starting at angle zero.
fn closed_circle_spline(z: Scalar, radius: Scalar) -> Curve3 {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let p = |x: Scalar, y: Scalar| Point3::new(x * radius, y * radius, z);
    Curve3::BSpline(BSplineCurve3 {
        degree: 2,
        control_points: vec![
            p(1.0, 0.0),
            p(1.0, 1.0),
            p(0.0, 1.0),
            p(-1.0, 1.0),
            p(-1.0, 0.0),
            p(-1.0, -1.0),
            p(0.0, -1.0),
            p(1.0, -1.0),
            p(1.0, 0.0),
        ],
        knots: vec![0.0, 0.25, 0.5, 0.75, 1.0],
        multiplicities: vec![3, 2, 2, 2, 3],
        weights: Some(vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0]),
        closed: true,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// The patch of a cylinder between angles `from` and `to` (`from < to`),
/// its arcs `low` and `high` trimmed at the patch's corners, its rulings
/// straight, no p-curve anywhere; with dense samples of the exact patch
/// and arcs.
fn trimmed_patch(
    low: Curve3,
    high: Curve3,
    (from, to): (Scalar, Scalar),
) -> (TriMesh, DeviationReport, Vec<Point3>) {
    let (radius, height) = (0.3, 0.5);
    let surface = Surface::Cylinder(Cylinder {
        frame: frame_at(0.0),
        radius,
    });
    let at = |t: Scalar, z: Scalar| Point3::new(radius * t.cos(), radius * t.sin(), z);
    let mut face = OneFace::new();
    let corners = [at(from, 0.0), at(to, 0.0), at(to, height), at(from, height)];
    let v: Vec<_> = corners.iter().map(|&p| face.vertex(p)).collect();
    let bottom = face.edge(v[0], v[1], Some(low));
    let top = face.edge(v[3], v[2], Some(high));
    let left = face.edge(v[1], v[2], None);
    let right = face.edge(v[0], v[3], None);
    face.use_edge(bottom, Orientation::Forward);
    face.use_edge(left, Orientation::Forward);
    face.use_edge(top, Orientation::Reversed);
    face.use_edge(right, Orientation::Reversed);
    let (mesh, report) = face
        .compile(&surface, &options(1e-4))
        .expect("a patch bounded by trimmed B-spline arcs compiles");
    let mut samples = surface_samples(&surface, (from, to), (0.0, height), 400);
    for z in [0.0, height] {
        samples.extend(span(from, to, 20_000).map(|t| at(t, z)));
    }
    (mesh, report, samples)
}

/// Every mesh vertex lies between the angles `from` and `to`: the arcs
/// were trimmed at the corners, not sampled whole.
fn assert_within_angles(name: &str, mesh: &TriMesh, (from, to): (Scalar, Scalar)) {
    let middle = 0.5 * (from + to);
    for p in &mesh.positions {
        let off = (p.y.atan2(p.x) - middle + 0.5 * TAU).rem_euclid(TAU) - 0.5 * TAU;
        assert!(
            off.abs() <= 0.5 * (to - from) + 1e-9,
            "{name}: a mesh vertex {p:?} lies outside the patch"
        );
    }
}

/// An open rational B-spline half circle trimmed at two interior vertices
/// (#333), on either side of its double knot: the span is located by the
/// certified projection and its bounds are split at the knot.
#[test]
fn an_open_bspline_edge_is_trimmed_at_interior_vertices() {
    let angles = (TAU / 12.0, TAU / 3.0);
    let (mesh, report, samples) = trimmed_patch(
        half_circle(0.0, 0.3, 0.0),
        half_circle(0.5, 0.3, 0.0),
        angles,
    );
    assert_within_angles("trimmed open B-spline", &mesh, angles);
    assert_bounded("trimmed open B-spline", &mesh, &samples, &report, 4.0);
}

/// A closed B-spline circle trimmed at two vertices either side of its seam
/// (#333): the span runs forward from the start vertex across the seam.
#[test]
fn a_closed_bspline_edge_is_trimmed_across_its_seam() {
    let angles = (-TAU / 8.0, TAU / 8.0);
    let (mesh, report, samples) = trimmed_patch(
        closed_circle_spline(0.0, 0.3),
        closed_circle_spline(0.5, 0.3),
        angles,
    );
    assert_within_angles("closed B-spline across its seam", &mesh, angles);
    assert_bounded(
        "closed B-spline across its seam",
        &mesh,
        &samples,
        &report,
        4.0,
    );
}

/// A closed B-spline circle whose one vertex lies a quarter turn from its
/// seam is the whole curve from that vertex (#333).
#[test]
fn a_closed_bspline_edge_is_the_whole_curve_from_its_one_vertex() {
    let (radius, height) = (0.3, 0.5);
    let surface = Surface::Cylinder(Cylinder {
        frame: frame_at(0.0),
        radius,
    });
    let rim = closed_circle_spline(0.0, radius);
    let mut face = OneFace::new();
    let bottom = face.vertex(Point3::new(0.0, radius, 0.0));
    let top = face.vertex(Point3::new(0.0, radius, height));
    let turned = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::new(0.0, 0.0, height),
            x: Vec3::Y,
            y: -Vec3::X,
            z: Vec3::Z,
        },
        radius,
    });
    let rim0 = face.edge(bottom, bottom, Some(rim.clone()));
    let rim1 = face.edge(top, top, Some(turned.clone()));
    let seam = face.edge(bottom, top, None);
    face.use_edge(rim0, Orientation::Forward);
    face.use_edge(seam, Orientation::Forward);
    face.use_edge(rim1, Orientation::Reversed);
    face.use_edge(seam, Orientation::Reversed);
    let (mesh, report) = face
        .compile(&surface, &options(1e-4))
        .expect("a closed B-spline trimmed at its one vertex compiles");
    let mut samples = surface_samples(&surface, (0.0, TAU), (0.0, height), 400);
    samples.extend(curve_samples(&rim, 20_000));
    samples.extend(curve_samples(&turned, 20_000));
    assert_bounded(
        "closed B-spline from its vertex",
        &mesh,
        &samples,
        &report,
        4.0,
    );
}

/// The area of a triangle mesh.
fn mesh_area(mesh: &TriMesh) -> Scalar {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            0.5 * (b - a).cross(c - a).length()
        })
        .sum()
}

/// Points every `step` along the closed polygon `corners`.
fn polygon_samples(corners: &[Point3], per_side: usize) -> Vec<Point3> {
    let mut out = Vec::new();
    for (i, &a) in corners.iter().enumerate() {
        let b = corners[(i + 1) % corners.len()];
        out.extend(span(0.0, 1.0, per_side).map(|s| a + (b - a) * s));
    }
    out
}

/// A planar face whose top edge is a zigzag polyline trimmed at two
/// vertices inside its end segments (#333), closed by straight edges: the
/// vertices are located on their segments exactly, and the edge is the
/// zigzag between them, its corners included within its bound.
#[test]
fn an_open_polyline_edge_is_trimmed_at_interior_vertices() {
    let p = |x: Scalar, y: Scalar| Point3::new(x, y, 0.0);
    let zigzag = Curve3::Polyline(axiolid_curve::Polyline3 {
        points: vec![
            p(-0.2, 0.0),
            p(0.25, 0.3),
            p(0.5, 0.0),
            p(0.75, 0.3),
            p(1.2, 0.0),
        ],
        closed: false,
    });
    let (a, b) = (p(0.025, 0.15), p(0.975, 0.15));
    let (c, d) = (p(0.975, -0.5), p(0.025, -0.5));
    let mut face = OneFace::new();
    let [va, vb, vc, vd] = [a, b, c, d].map(|q| face.vertex(q));
    let top = face.edge(va, vb, Some(zigzag));
    let down = face.edge(va, vd, None);
    let bottom = face.edge(vd, vc, None);
    let up = face.edge(vc, vb, None);
    face.use_edge(down, Orientation::Forward);
    face.use_edge(bottom, Orientation::Forward);
    face.use_edge(up, Orientation::Forward);
    face.use_edge(top, Orientation::Reversed);
    let plane = Surface::Plane(Plane {
        frame: frame_at(0.0),
    });
    let (mesh, report) = face
        .compile(&plane, &options(1e-4))
        .expect("a face bounded by a trimmed polyline compiles");
    let exact = [a, d, c, b, p(0.75, 0.3), p(0.5, 0.0), p(0.25, 0.3)];
    let shoelace: Scalar = (0..exact.len())
        .map(|i| {
            let (u, v) = (exact[i], exact[(i + 1) % exact.len()]);
            0.5 * (u.x * v.y - v.x * u.y)
        })
        .sum();
    let bound = report.bound.expect("a certified bound");
    let perimeter: Scalar = (0..exact.len())
        .map(|i| (exact[(i + 1) % exact.len()] - exact[i]).length())
        .sum();
    let area = mesh_area(&mesh);
    eprintln!("trimmed polyline: area {area}, exact {shoelace}, bound {bound:e}");
    assert!(
        (area - shoelace).abs() <= bound * perimeter,
        "area {area} misses {shoelace} by more than the bound allows"
    );
    for q in &mesh.positions {
        assert!(
            q.x >= 0.025 - 1e-12 && q.x <= 0.975 + 1e-12,
            "{q:?} past a vertex"
        );
    }
    let samples = polygon_samples(&exact, 4_000);
    assert_bounded("trimmed polyline", &mesh, &samples, &report, 6.0);
}

/// A planar face bounded by one closed polyline edge whose one vertex lies
/// inside a side (#333): the edge is the whole square from it, its corners
/// between samples and within its bound.
#[test]
fn a_closed_polyline_edge_is_the_whole_loop_from_its_one_vertex() {
    let p = |x: Scalar, y: Scalar| Point3::new(x, y, 0.0);
    let corners = [p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)];
    let square = Curve3::Polyline(axiolid_curve::Polyline3 {
        points: corners.to_vec(),
        closed: true,
    });
    let mut face = OneFace::new();
    let start = face.vertex(p(0.3, 0.0));
    let edge = face.edge(start, start, Some(square));
    face.use_edge(edge, Orientation::Forward);
    let plane = Surface::Plane(Plane {
        frame: frame_at(0.0),
    });
    let (mesh, report) = face
        .compile(&plane, &options(1e-3))
        .expect("a face bounded by a closed polyline compiles");
    let area = mesh_area(&mesh);
    let bound = report.bound.expect("a certified bound");
    eprintln!("closed polyline: area {area}, bound {bound:e}");
    assert!((area - 1.0).abs() <= 4.0 * bound, "area {area}");
    assert_bounded(
        "closed polyline",
        &mesh,
        &polygon_samples(&corners, 4_000),
        &report,
        6.0,
    );
}

/// A flat B-spline patch trimmed by one closed circle given only as its 3D
/// curve: the surface does not bend, so the triangles' own bound is zero
/// and the rim's bulge past its chords is all the report has to carry.
#[test]
fn the_edge_terms_are_carried_by_the_report() {
    let surface = Surface::BSpline(BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        control_points: vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
            vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ],
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![2, 2],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![2, 2],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    });
    let rim = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::new(0.5, 0.5, 0.0),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius: 0.2,
    });
    let mut face = OneFace::new();
    let start = face.vertex(Point3::new(0.7, 0.5, 0.0));
    let edge = face.edge(start, start, Some(rim.clone()));
    face.use_edge(edge, Orientation::Forward);
    let (mesh, report) = face
        .compile(&surface, &options(1e-3))
        .expect("a flat patch trimmed by a 3D circle compiles");
    let samples = curve_samples(&rim, 200_000);
    assert_bounded("rim on a flat patch", &mesh, &samples, &report, 4.0);
}

/// A closed cylinder the way an IFC advanced B-rep writes it: a wall with
/// one straight seam and two planar caps bounded by the wall's rims, no
/// pcurve anywhere. The bottom rim is `bottom` on the plane `floor`, the
/// top a circle at `height`; both start on the seam at `u = 0`.
fn capped_cylinder(
    radius: Scalar,
    height: Scalar,
    floor: Frame3,
    bottom: Curve3,
    seam: Option<Curve3>,
    budget: Scalar,
    caps_first: bool,
) -> (TriMesh, MeshClosure, DeviationReport) {
    let wall = Surface::Cylinder(Cylinder {
        frame: frame_at(0.0),
        radius,
    });
    let start = axiolid_reference::curve::evaluate3(&bottom, 0.0).unwrap();
    let mut solid = OneFace::new();
    let low = solid.vertex(start);
    let high = solid.vertex(Point3::new(radius, 0.0, height));
    let rim0 = solid.edge(low, low, Some(bottom));
    let rim1 = solid.edge(high, high, Some(circle(height, radius)));
    let seam = solid.edge(low, high, seam);
    solid.use_edge(rim0, Orientation::Forward);
    solid.use_edge(seam, Orientation::Forward);
    solid.use_edge(rim1, Orientation::Reversed);
    solid.use_edge(seam, Orientation::Reversed);
    solid.close_face(&wall);
    // The floor's frame normal points into the solid: the face runs
    // against it, and its loop runs clockwise about it.
    solid.use_edge(rim0, Orientation::Reversed);
    solid.close_oriented_face(
        &Surface::Plane(Plane { frame: floor }),
        Orientation::Reversed,
    );
    solid.use_edge(rim1, Orientation::Forward);
    solid.close_face(&Surface::Plane(Plane {
        frame: frame_at(height),
    }));
    if caps_first {
        // The caps then sample the rims, and the wall inverts their samples.
        solid.faces.rotate_left(1);
    }
    solid
        .compile_shell(true, &options(budget))
        .expect("a capped cylinder with only 3D edges compiles")
}

/// Every mesh edge is used by exactly two triangles, once each way.
fn assert_watertight(name: &str, mesh: &TriMesh) {
    let mut directed: HashMap<(u32, u32), usize> = HashMap::new();
    for t in mesh.indices.chunks_exact(3) {
        for (x, y) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *directed.entry((x, y)).or_default() += 1;
        }
    }
    for (&(x, y), &count) in &directed {
        assert_eq!(
            count, 1,
            "{name}: edge {x}-{y} is used {count} times one way"
        );
        assert_eq!(
            directed.get(&(y, x)),
            Some(&1),
            "{name}: edge {x}-{y} has no opposite use: the mesh is open there"
        );
    }
}

/// The enclosed volume is within the bound's effect, the bound times the
/// surface's area, of the exact one.
fn assert_volume(
    name: &str,
    mesh: &TriMesh,
    report: &DeviationReport,
    exact: Scalar,
    area: Scalar,
) {
    let bound = report
        .bound
        .unwrap_or_else(|| panic!("{name}: expected a certified bound, got {report:?}"));
    let volume = axiolid_inspect::enclosed_volume(mesh).expect("a closed mesh");
    let slack = bound * area;
    eprintln!(
        "{name}: volume in [{:e}, {:e}], exact {exact:e}, bound {bound:e}, slack {slack:e}",
        volume.lower, volume.upper
    );
    assert!(
        volume.lower - slack <= exact && exact <= volume.upper + slack,
        "{name}: exact volume {exact} outside [{}, {}] grown by {slack}",
        volume.lower,
        volume.upper
    );
}

#[test]
fn a_closed_cylinder_with_circular_caps_compiles_as_a_watertight_solid() {
    let (radius, height) = (0.3, 0.5);
    for (budget, caps_first) in [(1e-3, false), (1e-4, false), (1e-3, true)] {
        let (mesh, closure, report) = capped_cylinder(
            radius,
            height,
            frame_at(0.0),
            circle(0.0, radius),
            None,
            budget,
            caps_first,
        );
        let name = format!("capped cylinder at {budget:e}, caps first {caps_first}");
        assert_eq!(closure, MeshClosure::Solid, "{name}");
        assert_watertight(&name, &mesh);
        let area = TAU * radius * height + TAU * radius * radius;
        assert_volume(
            &name,
            &mesh,
            &report,
            0.5 * TAU * radius * radius * height,
            area,
        );
        let wall = Surface::Cylinder(Cylinder {
            frame: frame_at(0.0),
            radius,
        });
        let mut samples = surface_samples(&wall, (0.0, TAU), (0.0, height), 300);
        for z in [0.0, height] {
            for rr in span(0.0, radius, 60) {
                for t in span(0.0, TAU, 240) {
                    samples.push(Point3::new(rr * t.cos(), rr * t.sin(), z));
                }
            }
        }
        assert_bounded(&name, &mesh, &samples, &report, 6.0);
    }
}

/// The floor cut obliquely: the plane `z = k x` meets the wall in an
/// ellipse, the floor cap's only edge. The volume is still `pi r^2 h`,
/// `k x` averaging to zero over the disk.
#[test]
fn a_cylinder_cut_by_an_oblique_plane_compiles_as_a_watertight_solid() {
    let (radius, height, k) = (0.3, 0.5, 0.5);
    let tilt = Vec3::new(1.0, 0.0, k).normalize();
    let floor = Frame3 {
        origin: Point3::ZERO,
        x: tilt,
        y: Vec3::Y,
        z: tilt.cross(Vec3::Y),
    };
    let ellipse = Curve3::Ellipse(Ellipse3 {
        frame: floor,
        semi_axis_x: radius * (1.0 + k * k).sqrt(),
        semi_axis_y: radius,
    });
    let (mesh, closure, report) =
        capped_cylinder(radius, height, floor, ellipse.clone(), None, 1e-4, false);
    assert_eq!(closure, MeshClosure::Solid);
    assert_watertight("oblique cylinder", &mesh);
    let floor_area = 0.5 * TAU * radius * radius * (1.0 + k * k).sqrt();
    let area = TAU * radius * height + 0.5 * TAU * radius * radius + floor_area;
    assert_volume(
        "oblique cylinder",
        &mesh,
        &report,
        0.5 * TAU * radius * radius * height,
        area,
    );
    let mut samples = curve_samples(&ellipse, 20_000);
    for rr in span(0.0, 1.0, 60) {
        for t in span(0.0, TAU, 240) {
            let (x, y) = (radius * rr * t.cos(), radius * rr * t.sin());
            samples.push(Point3::new(x, y, k * x));
        }
    }
    assert_bounded("oblique cylinder", &mesh, &samples, &report, 6.0);
}

/// A D-shaped half cylinder the way an IFC advanced B-rep writes it (#336):
/// a half-cylinder wall and a planar face across its diameter sharing two
/// straight rulings, and two half-disk caps, no p-curve anywhere. The
/// faces are listed rotated left by `rotate` and then, with `reverse`,
/// backwards, so every face reaches the shared edges first in some order.
fn d_shaped_half_cylinder(
    radius: Scalar,
    height: Scalar,
    budget: Scalar,
    rotate: usize,
    reverse: bool,
) -> (TriMesh, MeshClosure, DeviationReport) {
    let wall = Surface::Cylinder(Cylinder {
        frame: frame_at(0.0),
        radius,
    });
    let mut solid = OneFace::new();
    let a = solid.vertex(Point3::new(0.0, -radius, 0.0));
    let b = solid.vertex(Point3::new(0.0, radius, 0.0));
    let c = solid.vertex(Point3::new(0.0, radius, height));
    let d = solid.vertex(Point3::new(0.0, -radius, height));
    // The arcs run forward about the axis from the start vertex at -pi/2
    // to the end vertex at pi/2, through x = r (#319).
    let arc_low = solid.edge(a, b, Some(circle(0.0, radius)));
    let arc_high = solid.edge(d, c, Some(circle(height, radius)));
    let ruling_a = solid.edge(a, d, None);
    let ruling_b = solid.edge(b, c, None);
    let diameter_low = solid.edge(a, b, None);
    let diameter_high = solid.edge(d, c, None);
    // The wall, counter-clockwise in (u, v): outward.
    solid.use_edge(arc_low, Orientation::Forward);
    solid.use_edge(ruling_b, Orientation::Forward);
    solid.use_edge(arc_high, Orientation::Reversed);
    solid.use_edge(ruling_a, Orientation::Reversed);
    solid.close_face(&wall);
    // The flat side at x = 0, facing -x.
    solid.use_edge(ruling_a, Orientation::Forward);
    solid.use_edge(diameter_high, Orientation::Forward);
    solid.use_edge(ruling_b, Orientation::Reversed);
    solid.use_edge(diameter_low, Orientation::Reversed);
    solid.close_face(&Surface::Plane(Plane {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::Z,
            y: Vec3::Y,
            z: -Vec3::X,
        },
    }));
    // The floor faces -z, against its plane's normal.
    solid.use_edge(diameter_low, Orientation::Forward);
    solid.use_edge(arc_low, Orientation::Reversed);
    solid.close_oriented_face(
        &Surface::Plane(Plane {
            frame: frame_at(0.0),
        }),
        Orientation::Reversed,
    );
    solid.use_edge(arc_high, Orientation::Forward);
    solid.use_edge(diameter_high, Orientation::Reversed);
    solid.close_face(&Surface::Plane(Plane {
        frame: frame_at(height),
    }));
    solid.faces.rotate_left(rotate);
    if reverse {
        solid.faces.reverse();
    }
    solid
        .compile_shell(true, &options(budget))
        .expect("a D-shaped half cylinder with only 3D edges compiles")
}

/// The wall samples each ruling in pieces sized by its curvature; the flat
/// side takes the same pieces, whichever face comes first (#336).
#[test]
fn a_d_shaped_half_cylinder_compiles_as_a_watertight_solid_in_any_face_order() {
    let (radius, height) = (0.3, 0.5);
    let wall = Surface::Cylinder(Cylinder {
        frame: frame_at(0.0),
        radius,
    });
    let quarter = 0.25 * TAU;
    let mut samples = surface_samples(&wall, (-quarter, quarter), (0.0, height), 300);
    for z in [0.0, height] {
        for rr in span(0.0, radius, 60) {
            for t in span(-quarter, quarter, 120) {
                samples.push(Point3::new(rr * t.cos(), rr * t.sin(), z));
            }
        }
    }
    for y in span(-radius, radius, 60) {
        for z in span(0.0, height, 60) {
            samples.push(Point3::new(0.0, y, z));
        }
    }
    let area = 0.5 * TAU * radius * height + 2.0 * radius * height + 0.5 * TAU * radius * radius;
    let mut counts = Vec::new();
    for budget in [1e-3, 1e-4] {
        for rotate in 0..4 {
            for reverse in [false, true] {
                let (mesh, closure, report) =
                    d_shaped_half_cylinder(radius, height, budget, rotate, reverse);
                let name = format!("D-shape at {budget:e}, rotated {rotate}, reversed {reverse}");
                assert_eq!(closure, MeshClosure::Solid, "{name}");
                assert_watertight(&name, &mesh);
                assert_volume(
                    &name,
                    &mesh,
                    &report,
                    0.25 * TAU * radius * radius * height,
                    area,
                );
                assert_bounded(&name, &mesh, &samples, &report, 6.0);
                counts.push((budget, mesh.positions.len(), mesh.indices.len()));
            }
        }
    }
    // The order changes nothing the mesh is made of.
    for &(budget, positions, indices) in &counts {
        let first = counts.iter().find(|c| c.0 == budget).unwrap();
        assert_eq!((positions, indices), (first.1, first.2), "{budget:e}");
    }
}
