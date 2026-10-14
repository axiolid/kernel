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
use axiolid_curve::{BSplineCurve3, Circle3, Curve3, KnotSpec};
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationReport, ReferenceMeshCompiler};
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId};
use axiolid_surface::{BSplineSurface, Cone, Cylinder, Sphere, Surface};
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
            orientation: Orientation::Forward,
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
        let shell = self.brep.add_shell(Shell {
            faces: self
                .faces
                .iter()
                .map(|&face| (face, Orientation::Forward))
                .collect(),
            closed: false,
        });
        self.brep.add_solid(Solid {
            outer: shell,
            voids: Vec::new(),
        });
        let root = self.builder.push(GeometryNode::BRep(self.brep)).unwrap();
        let graph = self.builder.finish(vec![root]).unwrap();
        let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
            .compile_mesh_with_deviation(&graph, root, options)?;
        Ok((outcome.mesh, report))
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

/// Without a pcurve or an edge interval the curve's whole domain is the
/// edge: a circle that starts a quarter turn away from its vertex is
/// refused, not trimmed by a guess.
#[test]
fn an_edge_curve_that_does_not_start_at_its_vertex_is_refused() {
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
    let rim0 = face.edge(bottom, bottom, Some(turned));
    let rim1 = face.edge(top, top, Some(circle(height, radius)));
    let seam = face.edge(bottom, top, None);
    face.use_edge(rim0, Orientation::Forward);
    face.use_edge(seam, Orientation::Forward);
    face.use_edge(rim1, Orientation::Reversed);
    face.use_edge(seam, Orientation::Reversed);
    let error = face
        .compile(&surface, &options(1e-3))
        .expect_err("a curve that misses its vertex is refused");
    assert!(
        matches!(&error, GeomError::InvalidInput(reason) if reason.contains("whole domain")),
        "got {error:?}"
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
