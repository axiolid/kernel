//! B-rep faces on swept and ruled surface relations (#411, ADR 0103).
//!
//! A face over a `Revolution`, `LinearExtrusion` or `Ruled` relation is
//! meshed when every edge use carries a pcurve, with a bound certified from
//! the relation's closed-form derivative bounds. Closed solids are a wall
//! on the relation with two planar caps, each checked against the same
//! solid on the atomic surface it coincides with, for the fixture's sake.
//! Open sheets are checked against a closed-form point of the surface.
//! Every case reads the mesh back: closure, watertightness, volume and a
//! dense sample of the exact surface against the reported bound, which is
//! certified within the budget or the case fails.

use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve3, Circle2, Circle3, Curve2, Curve3, KnotSpec, Line3, Polyline2, Polyline3,
};
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, DeviationReport, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshClosure;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId, SurfaceRelation};
use axiolid_surface::{Cone, Cylinder, Plane, Sphere, Surface, Torus};
use axiolid_topology::{
    BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};

const BUDGET: Scalar = 1e-3;
const NO_PCURVE: &str = "an edge use without a pcurve on a face over a surface relation";

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(BUDGET, 1e-9).unwrap())
        .with_chord_error(BUDGET)
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

/// Pushes the support node of a face and returns it.
type Support<'a> = &'a dyn Fn(&mut GeometryGraphBuilder) -> NodeId;

fn push_curve3(builder: &mut GeometryGraphBuilder, curve: Curve3) -> NodeId {
    builder.push(GeometryNode::Curve3(curve)).unwrap()
}

fn line3(origin: Point3, direction: Vec3) -> Curve3 {
    Curve3::Line(Line3 { origin, direction })
}

fn polyline(points: &[(Scalar, Scalar)]) -> Curve2 {
    Curve2::Polyline(Polyline2 {
        points: points.iter().map(|&(u, v)| Point2::new(u, v)).collect(),
        closed: false,
    })
}

/// A circle about the vertical line through `(x, y)`, at height `z`.
fn horizontal_circle(x: Scalar, y: Scalar, z: Scalar, radius: Scalar) -> Curve3 {
    Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::new(x, y, z),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius,
    })
}

/// The circle of `radius` in the plane of `x` and `y` about `origin`.
fn circle_in(origin: Point3, x: Vec3, y: Vec3, radius: Scalar) -> Curve3 {
    Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin,
            x,
            y,
            z: x.cross(y),
        },
        radius,
    })
}

fn revolution(
    builder: &mut GeometryGraphBuilder,
    curve: Curve3,
    axis_origin: Point3,
    axis_direction: Vec3,
) -> NodeId {
    let swept_curve = push_curve3(builder, curve);
    builder
        .push(GeometryNode::SurfaceRelation(SurfaceRelation::Revolution {
            swept_curve,
            axis_origin,
            axis_direction,
        }))
        .unwrap()
}

fn atomic(builder: &mut GeometryGraphBuilder, surface: Surface) -> NodeId {
    builder.push(GeometryNode::Surface(surface)).unwrap()
}

/// A closed solid: a wall on the relation, closed by two planar caps.
struct Wall<'a> {
    support: Support<'a>,
    /// The wall's `u` at the start and the end of the turn.
    u: (Scalar, Scalar),
    /// The wall's `v` at the bottom and the top rim.
    v: (Scalar, Scalar),
    /// The rims' heights and radii.
    heights: (Scalar, Scalar),
    radii: (Scalar, Scalar),
    /// The seam's 3D curve, when it is not straight.
    seam: Option<Curve3>,
}

/// A cap's circle about its plane's origin, anticlockwise or clockwise.
fn cap_circle(radius: Scalar, anticlockwise: bool) -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: if anticlockwise { Vec2::Y } else { -Vec2::Y },
        },
        radius,
    })
}

/// What a compilation returns.
type Compiled = (TriMesh, MeshClosure, DeviationReport);

fn compile(mut builder: GeometryGraphBuilder, brep: BRep<NodeId>) -> GeomResult<Compiled> {
    let root = builder.push(GeometryNode::BRep(brep)).unwrap();
    let graph = builder.finish(vec![root]).unwrap();
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options())?;
    Ok((outcome.mesh, outcome.closure, report))
}

/// A closed solid, the way an exchange format with pcurves writes it: the
/// wall with one seam and two planar caps, every edge's every use with its
/// own pcurve.
fn closed_wall(wall: &Wall<'_>) -> GeomResult<Compiled> {
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let (u0, u1) = wall.u;
    let (vb, vt) = wall.v;
    let ((zb, zt), (rb, rt)) = (wall.heights, wall.radii);
    let pc = |builder: &mut GeometryGraphBuilder, curve| {
        builder.push(GeometryNode::Curve2(curve)).unwrap()
    };
    let bottom = pc(&mut builder, polyline(&[(u0, vb), (u1, vb)]));
    let top = pc(&mut builder, polyline(&[(u1, vt), (u0, vt)]));
    let seam_up = pc(&mut builder, polyline(&[(u1, vb), (u1, vt)]));
    let seam_down = pc(&mut builder, polyline(&[(u0, vt), (u0, vb)]));
    let floor_circle = pc(&mut builder, cap_circle(rb, false));
    let lid_circle = pc(&mut builder, cap_circle(rt, true));
    let wall_support = (wall.support)(&mut builder);
    let floor = atomic(
        &mut builder,
        Surface::Plane(Plane {
            frame: frame_at(zb),
        }),
    );
    let lid = atomic(
        &mut builder,
        Surface::Plane(Plane {
            frame: frame_at(zt),
        }),
    );
    let low = brep.add_vertex(Vertex {
        position: Point3::new(rb, 0.0, zb),
    });
    let high = brep.add_vertex(Vertex {
        position: Point3::new(rt, 0.0, zt),
    });
    let rims = [
        horizontal_circle(0.0, 0.0, zb, rb),
        horizontal_circle(0.0, 0.0, zt, rt),
    ];
    let mut edge = |start, end, curve: Option<Curve3>| {
        let curve = curve.map(|c| builder.push(GeometryNode::Curve3(c)).unwrap());
        brep.add_edge(Edge { start, end, curve })
    };
    let [bottom_curve, top_curve] = rims;
    let (rim0, rim1, seam) = (
        edge(low, low, Some(bottom_curve)),
        edge(high, high, Some(top_curve)),
        edge(low, high, wall.seam.clone()),
    );
    let used = |edge, orientation, pcurve| EdgeUse {
        edge,
        orientation,
        pcurve: Some(pcurve),
    };
    let wall_loop = brep.add_loop(Loop {
        edges: vec![
            used(rim0, Orientation::Forward, bottom),
            used(seam, Orientation::Forward, seam_up),
            used(rim1, Orientation::Reversed, top),
            used(seam, Orientation::Reversed, seam_down),
        ],
    });
    let floor_loop = brep.add_loop(Loop {
        edges: vec![used(rim0, Orientation::Reversed, floor_circle)],
    });
    let lid_loop = brep.add_loop(Loop {
        edges: vec![used(rim1, Orientation::Forward, lid_circle)],
    });
    let mut face = |surface, loop_id, orientation| {
        brep.add_face(Face {
            surface: Some(surface),
            bounds: vec![FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation,
        })
    };
    // The floor's frame normal points into the solid: the face runs against
    // it, and its loop runs clockwise about it.
    let faces = [
        face(wall_support, wall_loop, Orientation::Forward),
        face(floor, floor_loop, Orientation::Reversed),
        face(lid, lid_loop, Orientation::Forward),
    ];
    let shell = brep.add_shell(Shell {
        faces: faces.map(|face| (face, Orientation::Forward)).to_vec(),
        closed: true,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    compile(builder, brep)
}

/// A hole in a sheet: a closed edge whose pcurve is a circle in `(u, v)`.
struct Hole {
    centre: (Scalar, Scalar),
    radius: Scalar,
}

/// One face over the rectangle `u` by `v`, an open sheet. Its edges run
/// bottom (`v` low), right (`u` high), top and left, each from the lower
/// corner; the top and the left are used reversed.
struct Sheet<'a> {
    u: (Scalar, Scalar),
    v: (Scalar, Scalar),
    /// The exact surface point.
    point: &'a dyn Fn(Scalar, Scalar) -> Point3,
    /// The edges' 3D curves in the order above, `None` for a straight one.
    curves: [Option<Curve3>; 4],
    /// Whether the edge uses carry pcurves.
    pcurves: bool,
    hole: Option<Hole>,
}

fn open_sheet(support: Support<'_>, sheet: &Sheet<'_>) -> GeomResult<Compiled> {
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let ((u0, u1), (v0, v1)) = (sheet.u, sheet.v);
    let support = support(&mut builder);
    let corners = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)].map(|(u, v)| {
        brep.add_vertex(Vertex {
            position: (sheet.point)(u, v),
        })
    });
    // Each edge from the lower corner, with the pcurve in the use's direction.
    let runs = [
        (
            corners[0],
            corners[1],
            polyline(&[(u0, v0), (u1, v0)]),
            Orientation::Forward,
        ),
        (
            corners[1],
            corners[2],
            polyline(&[(u1, v0), (u1, v1)]),
            Orientation::Forward,
        ),
        (
            corners[3],
            corners[2],
            polyline(&[(u1, v1), (u0, v1)]),
            Orientation::Reversed,
        ),
        (
            corners[0],
            corners[3],
            polyline(&[(u0, v1), (u0, v0)]),
            Orientation::Reversed,
        ),
    ];
    let mut uses = Vec::new();
    for ((start, end, pcurve, orientation), curve) in runs.into_iter().zip(&sheet.curves) {
        let curve = curve
            .clone()
            .map(|c| builder.push(GeometryNode::Curve3(c)).unwrap());
        let edge = brep.add_edge(Edge { start, end, curve });
        let pcurve = sheet
            .pcurves
            .then(|| builder.push(GeometryNode::Curve2(pcurve)).unwrap());
        uses.push(EdgeUse {
            edge,
            orientation,
            pcurve,
        });
    }
    let outer = brep.add_loop(Loop { edges: uses });
    let mut bounds = vec![FaceBound {
        loop_id: outer,
        orientation: Orientation::Forward,
        outer: true,
    }];
    if let Some(Hole {
        centre: (cu, cv),
        radius,
    }) = &sheet.hole
    {
        let vertex = brep.add_vertex(Vertex {
            position: (sheet.point)(cu + radius, *cv),
        });
        let edge = brep.add_edge(Edge {
            start: vertex,
            end: vertex,
            curve: None,
        });
        let circle = builder
            .push(GeometryNode::Curve2(Curve2::Circle(Circle2 {
                frame: Frame2 {
                    origin: Point2::new(*cu, *cv),
                    x: Vec2::X,
                    y: -Vec2::Y,
                },
                radius: *radius,
            })))
            .unwrap();
        let inner = brep.add_loop(Loop {
            edges: vec![EdgeUse {
                edge,
                orientation: Orientation::Forward,
                pcurve: Some(circle),
            }],
        });
        bounds.push(FaceBound {
            loop_id: inner,
            orientation: Orientation::Forward,
            outer: false,
        });
    }
    let face = brep.add_face(Face {
        surface: Some(support),
        bounds,
        orientation: Orientation::Forward,
    });
    brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    compile(builder, brep)
}

/// Every mesh edge is used by exactly two triangles, once each way.
fn assert_watertight(name: &str, mesh: &TriMesh) {
    let mut directed: HashMap<(u32, u32), usize> = HashMap::new();
    for t in mesh.indices.as_chunks::<3>().0 {
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

/// The farthest any of `samples` lies from the mesh.
fn farthest(mesh: &TriMesh, samples: &[Point3]) -> Scalar {
    let triangles: Vec<[Point3; 3]> = mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| [0, 1, 2].map(|i| mesh.positions[t[i] as usize]))
        .collect();
    samples
        .iter()
        .map(|&p| {
            triangles
                .iter()
                .filter_map(|&t| closest_point_on_triangle(p, t).ok())
                .map(|q| p.distance(q))
                .fold(Scalar::INFINITY, Scalar::min)
        })
        .fold(0.0, Scalar::max)
}

fn signed_volume(mesh: &TriMesh) -> Scalar {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            a.dot(b.cross(c)) / 6.0
        })
        .sum()
}

/// The reported bound, certified and within the budget, with a certified
/// contribution of the wall's family.
fn certified(name: &str, report: &DeviationReport, family: &str) -> Scalar {
    let bound = certified_bound(name, report, family);
    assert!(
        report.meets_requested(),
        "{name}: bound {bound:e} above the budget"
    );
    bound
}

/// The reported bound, certified, whether or not it meets the budget.
fn certified_bound(name: &str, report: &DeviationReport, family: &str) -> Scalar {
    let bound = report
        .bound
        .unwrap_or_else(|| panic!("{name}: expected a certified bound, got {report:?}"));
    assert!(
        report
            .contributions
            .iter()
            .any(|c| c.path == DeviationPath::BRepFace
                && c.detail == family
                && matches!(c.bound, DeviationBound::Certified(_))),
        "{name}: no certified {family} contribution in {report:?}"
    );
    bound
}

/// What a closed wall must come to.
struct Solid3<'a> {
    /// The exact wall point at the turn fraction `t` and height fraction `s`.
    wall: &'a dyn Fn(Scalar, Scalar) -> Point3,
    volume: Scalar,
    area: Scalar,
    family: &'static str,
}

fn assert_closed_wall(name: &str, wall: &Wall<'_>, expect: &Solid3<'_>) {
    let (mesh, closure, report) = closed_wall(wall).unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_eq!(closure, MeshClosure::Solid, "{name}");
    assert_watertight(name, &mesh);
    let bound = certified(name, &report, expect.family);
    // A signed sum, not `enclosed_volume`: a rim whose boundary vertices are
    // collinear in the wall's parameters leaves ears of exactly no area.
    let volume = signed_volume(&mesh);
    let slack = bound * expect.area;
    assert!(
        (volume - expect.volume).abs() <= slack,
        "{name}: mesh volume {volume} against the exact {}, slack {slack}",
        expect.volume
    );
    let ((zb, zt), (rb, rt)) = (wall.heights, wall.radii);
    let mut samples = Vec::new();
    for k in 0..720 {
        let angle = TAU * k as Scalar / 720.0;
        let (c, s) = (angle.cos(), angle.sin());
        for (z, r) in [(zb, rb), (zt, rt)] {
            for f in [0.25, 0.5, 0.75, 1.0] {
                samples.push(Point3::new(f * r * c, f * r * s, z));
            }
        }
        for h in 0..=10 {
            samples.push((expect.wall)(angle / TAU, h as Scalar / 10.0));
        }
    }
    let worst = farthest(&mesh, &samples);
    assert!(
        worst <= bound,
        "{name}: a surface point lies {worst:e} from the mesh, above the reported {bound:e}"
    );
}

fn assert_sheet(name: &str, support: Support<'_>, sheet: &Sheet<'_>, family: &str) -> Compiled {
    check_sheet(name, support, sheet, family, true)
}

/// `assert_sheet`; `budgeted` false lets the certified bound exceed the budget.
fn check_sheet(
    name: &str,
    support: Support<'_>,
    sheet: &Sheet<'_>,
    family: &str,
    budgeted: bool,
) -> Compiled {
    let result = open_sheet(support, sheet).unwrap_or_else(|e| panic!("{name}: {e}"));
    let (mesh, closure, report) = &result;
    assert_eq!(*closure, MeshClosure::Surface, "{name}");
    let bound = if budgeted {
        certified(name, report, family)
    } else {
        certified_bound(name, report, family)
    };
    let ((u0, u1), (v0, v1)) = (sheet.u, sheet.v);
    let mut samples = Vec::new();
    let (nu, nv) = (61, 21);
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * (i as Scalar + 0.37) / (nu as Scalar + 0.74);
            let v = v0 + (v1 - v0) * (j as Scalar + 0.41) / (nv as Scalar + 0.82);
            let (u, v) = (u.clamp(u0, u1), v.clamp(v0, v1));
            let skip = sheet.hole.as_ref().is_some_and(|h| {
                ((u - h.centre.0).powi(2) + (v - h.centre.1).powi(2)).sqrt() < 1.2 * h.radius
            });
            if !skip {
                samples.push((sheet.point)(u, v));
            }
        }
    }
    for k in 0..=200 {
        let t = k as Scalar / 200.0;
        samples.push((sheet.point)(u0 + (u1 - u0) * t, v0));
        samples.push((sheet.point)(u0 + (u1 - u0) * t, v1));
        samples.push((sheet.point)(u0, v0 + (v1 - v0) * t));
        samples.push((sheet.point)(u1, v0 + (v1 - v0) * t));
    }
    let worst = farthest(mesh, &samples);
    assert!(
        worst <= bound * (1.0 + 1e-9),
        "{name}: a surface point lies {worst:e} from the mesh, above the reported {bound:e}"
    );
    result
}

// --- the fixture, on atomic surfaces ------------------------------------------

const RADIUS: Scalar = 0.3;
const HEIGHT: Scalar = 0.5;

fn cylinder_wall(point: Support<'_>) -> Wall<'_> {
    Wall {
        support: point,
        u: (0.0, TAU),
        v: (0.0, HEIGHT),
        heights: (0.0, HEIGHT),
        radii: (RADIUS, RADIUS),
        seam: None,
    }
}

fn cylinder_point(t: Scalar, s: Scalar) -> Point3 {
    let angle = TAU * t;
    Point3::new(RADIUS * angle.cos(), RADIUS * angle.sin(), HEIGHT * s)
}

fn cylinder_solid() -> Solid3<'static> {
    Solid3 {
        wall: &cylinder_point,
        volume: 0.1413716694115407,
        area: 1.5079644737231008,
        family: "surface of revolution",
    }
}

/// The frustum of radii 0.2 at the floor and 0.3 at the lid.
fn frustum_wall(support: Support<'_>) -> Wall<'_> {
    Wall {
        support,
        u: (0.0, TAU),
        v: (0.0, HEIGHT),
        heights: (0.0, HEIGHT),
        radii: (0.2, 0.3),
        seam: None,
    }
}

fn frustum_point(t: Scalar, s: Scalar) -> Point3 {
    let (angle, r) = (TAU * t, 0.2 + 0.1 * s);
    Point3::new(r * angle.cos(), r * angle.sin(), HEIGHT * s)
}

fn frustum_solid() -> Solid3<'static> {
    Solid3 {
        wall: &frustum_point,
        volume: 0.09948376736367677,
        area: 1.2093591571873779,
        family: "surface of revolution",
    }
}

/// The zone of the sphere of radius 0.3 between latitudes -45 and 45 degrees.
fn zone_circle() -> Curve3 {
    circle_in(Point3::ZERO, Vec3::X, Vec3::Z, 0.3)
}

fn zone_wall(support: Support<'_>) -> Wall<'_> {
    let h = 0.3 * FRAC_PI_4.sin();
    Wall {
        support,
        u: (0.0, TAU),
        v: (-FRAC_PI_4, FRAC_PI_4),
        heights: (-h, h),
        radii: (h, h),
        seam: Some(zone_circle()),
    }
}

fn zone_point(t: Scalar, s: Scalar) -> Point3 {
    let (angle, latitude) = (TAU * t, -FRAC_PI_4 + FRAC_PI_2 * s);
    Point3::new(
        0.3 * latitude.cos() * angle.cos(),
        0.3 * latitude.cos() * angle.sin(),
        0.3 * latitude.sin(),
    )
}

fn zone_solid() -> Solid3<'static> {
    Solid3 {
        wall: &zone_point,
        volume: 0.09996486610856323,
        area: 1.0824622676915872,
        family: "surface of revolution",
    }
}

#[test]
fn atomic_controls_compile_certified() {
    let cylinder = |b: &mut GeometryGraphBuilder| {
        atomic(
            b,
            Surface::Cylinder(Cylinder {
                frame: frame_at(0.0),
                radius: RADIUS,
            }),
        )
    };
    let cone = |b: &mut GeometryGraphBuilder| {
        atomic(
            b,
            Surface::Cone(Cone {
                frame: frame_at(0.0),
                radius: 0.2,
                semi_angle: 0.2_f64.atan(),
            }),
        )
    };
    let sphere = |b: &mut GeometryGraphBuilder| {
        atomic(
            b,
            Surface::Sphere(Sphere {
                frame: frame_at(0.0),
                radius: 0.3,
            }),
        )
    };
    let family = |name| Solid3 {
        family: name,
        ..cylinder_solid()
    };
    assert_closed_wall("cylinder", &cylinder_wall(&cylinder), &family("cylinder"));
    assert_closed_wall(
        "cone",
        &frustum_wall(&cone),
        &Solid3 {
            family: "cone",
            ..frustum_solid()
        },
    );
    assert_closed_wall(
        "sphere",
        &zone_wall(&sphere),
        &Solid3 {
            family: "sphere",
            ..zone_solid()
        },
    );
}

// --- revolution: closed solids --------------------------------------------------

#[test]
fn a_line_revolved_about_a_parallel_axis_meshes_as_its_cylinder() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(
            b,
            line3(Point3::new(0.3, 0.0, 0.0), Vec3::Z),
            Point3::ZERO,
            Vec3::Z,
        )
    };
    assert_closed_wall("revolved line", &cylinder_wall(&support), &cylinder_solid());
}

#[test]
fn a_slanted_line_revolved_meshes_as_its_frustum() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(
            b,
            line3(Point3::new(0.2, 0.0, 0.0), Vec3::new(0.2, 0.0, 1.0)),
            Point3::ZERO,
            Vec3::Z,
        )
    };
    assert_closed_wall(
        "revolved slanted line",
        &frustum_wall(&support),
        &frustum_solid(),
    );
}

#[test]
fn an_arc_revolved_meshes_as_a_spherical_zone() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(
            b,
            circle_in(Point3::ZERO, Vec3::X, Vec3::Z, 0.3),
            Point3::ZERO,
            Vec3::Z,
        )
    };
    assert_closed_wall("revolved arc", &zone_wall(&support), &zone_solid());
}

// --- revolution: open sheets ----------------------------------------------------

/// The patch of a torus a Revit surface of revolution would carry: an arc
/// profile about a vertical axis through `(1, 2)`.
fn torus_profile() -> Curve3 {
    circle_in(Point3::new(1.4, 2.0, 0.1), Vec3::X, Vec3::Z, 0.1)
}

fn torus_point(u: Scalar, v: Scalar) -> Point3 {
    let reach = 0.4 + 0.1 * v.cos();
    Point3::new(
        1.0 + reach * u.cos(),
        2.0 + reach * u.sin(),
        0.1 + 0.1 * v.sin(),
    )
}

fn torus_patch<'a>(point: &'a dyn Fn(Scalar, Scalar) -> Point3) -> Sheet<'a> {
    let (u0, u1) = (0.3, 1.5);
    let meridian = |u: Scalar| {
        circle_in(
            Point3::new(1.0 + 0.4 * u.cos(), 2.0 + 0.4 * u.sin(), 0.1),
            Vec3::new(u.cos(), u.sin(), 0.0),
            Vec3::Z,
            0.1,
        )
    };
    Sheet {
        u: (u0, u1),
        v: (-FRAC_PI_2, 0.0),
        point,
        curves: [
            Some(horizontal_circle(1.0, 2.0, 0.0, 0.4)),
            Some(meridian(u1)),
            Some(horizontal_circle(1.0, 2.0, 0.1, 0.5)),
            Some(meridian(u0)),
        ],
        pcurves: true,
        hole: None,
    }
}

#[test]
fn a_revolved_arc_patch_like_revit_s_meshes_within_its_bound() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(b, torus_profile(), Point3::new(1.0, 2.0, 0.0), Vec3::Z)
    };
    let sheet = torus_patch(&torus_point);
    assert_sheet(
        "revolved arc patch",
        &support,
        &sheet,
        "surface of revolution",
    );
    let control = |b: &mut GeometryGraphBuilder| {
        atomic(
            b,
            Surface::Torus(Torus {
                frame: Frame3 {
                    origin: Point3::new(1.0, 2.0, 0.1),
                    x: Vec3::X,
                    y: Vec3::Y,
                    z: Vec3::Z,
                },
                major_radius: 0.4,
                minor_radius: 0.1,
            }),
        )
    };
    assert_sheet("torus control", &control, &sheet, "torus");
}

fn bezier_points() -> [Point3; 4] {
    [
        Point3::new(0.3, 0.0, 0.0),
        Point3::new(0.35, 0.0, 0.1),
        Point3::new(0.25, 0.0, 0.2),
        Point3::new(0.3, 0.0, 0.3),
    ]
}

fn bezier(points: [Point3; 4]) -> Curve3 {
    Curve3::BSpline(BSplineCurve3 {
        degree: 3,
        control_points: points.to_vec(),
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// The profile at `v` by de Casteljau, turned by `u` about the vertical axis.
fn bezier_point(u: Scalar, v: Scalar) -> Point3 {
    let mut p = bezier_points().to_vec();
    while p.len() > 1 {
        p = p.windows(2).map(|w| w[0].lerp(w[1], v)).collect();
    }
    Point3::new(p[0].x * u.cos(), p[0].x * u.sin(), p[0].z)
}

#[test]
fn a_bezier_profile_revolved_like_a_revit_hermite_meshes_within_its_bound() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(b, bezier(bezier_points()), Point3::ZERO, Vec3::Z)
    };
    let turned = bezier_points().map(|p| Point3::new(-p.x, -p.y, p.z));
    let sheet = Sheet {
        u: (0.0, PI),
        v: (0.0, 1.0),
        point: &bezier_point,
        curves: [
            Some(horizontal_circle(0.0, 0.0, 0.0, 0.3)),
            Some(bezier(turned)),
            Some(horizontal_circle(0.0, 0.0, 0.3, 0.3)),
            Some(bezier(bezier_points())),
        ],
        pcurves: true,
        hole: None,
    };
    assert_sheet("revolved Bezier", &support, &sheet, "surface of revolution");
}

fn cylinder_sheet_point(u: Scalar, v: Scalar) -> Point3 {
    Point3::new(RADIUS * u.cos(), RADIUS * u.sin(), v)
}

fn cylinder_sheet<'a>(pcurves: bool, hole: Option<Hole>) -> Sheet<'a> {
    Sheet {
        u: (0.0, PI),
        v: (0.0, HEIGHT),
        point: &cylinder_sheet_point,
        curves: [
            Some(horizontal_circle(0.0, 0.0, 0.0, RADIUS)),
            None,
            Some(horizontal_circle(0.0, 0.0, HEIGHT, RADIUS)),
            None,
        ],
        pcurves,
        hole,
    }
}

fn cylinder_revolution(b: &mut GeometryGraphBuilder) -> NodeId {
    revolution(
        b,
        line3(Point3::new(0.3, 0.0, 0.0), Vec3::Z),
        Point3::ZERO,
        Vec3::Z,
    )
}

/// A hole on a curved face is bounded above the budget today, on an atomic
/// surface as much as on a relation: the relation's bound is the atomic
/// one's, which this checks instead of the budget.
#[test]
fn a_revolved_sheet_with_a_hole_meshes_round_it() {
    let hole = || Hole {
        centre: (FRAC_PI_2, 0.25),
        radius: 0.1,
    };
    let (mesh, _, report) = check_sheet(
        "revolved sheet with a hole",
        &cylinder_revolution,
        &cylinder_sheet(true, Some(hole())),
        "surface of revolution",
        false,
    );
    let control = |b: &mut GeometryGraphBuilder| {
        atomic(
            b,
            Surface::Cylinder(Cylinder {
                frame: frame_at(0.0),
                radius: RADIUS,
            }),
        )
    };
    let (_, _, atomic_report) = check_sheet(
        "cylinder sheet with a hole",
        &control,
        &cylinder_sheet(true, Some(hole())),
        "cylinder",
        false,
    );
    let (bound, atomic_bound) = (report.bound.unwrap(), atomic_report.bound.unwrap());
    assert!(
        (bound - atomic_bound).abs() <= 1e-9 * atomic_bound,
        "{bound:e} against the atomic cylinder's {atomic_bound:e}"
    );
    // The hole is 0.1 radians wide, which is 0.03 along the wall: inside it,
    // shrunk by the budget, no triangle has its centre.
    for t in mesh.indices.as_chunks::<3>().0 {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        let m = (a + b + c) / 3.0;
        let (u, z) = (m.y.atan2(m.x), m.z);
        let inside = ((u - FRAC_PI_2) * RADIUS / (0.1 * RADIUS - BUDGET)).powi(2)
            + ((z - 0.25) / (0.1 - BUDGET)).powi(2);
        assert!(inside > 1.0, "a triangle centred in the hole at ({u}, {z})");
    }
}

// --- refusals ---------------------------------------------------------------------

fn refused(result: GeomResult<Compiled>) -> GeomError {
    match result {
        Ok((_, _, report)) => panic!("expected a refusal, got {report:?}"),
        Err(error) => error,
    }
}

#[test]
fn a_revolved_face_whose_edges_have_no_pcurve_is_refused_by_name() {
    let error = refused(open_sheet(
        &cylinder_revolution,
        &cylinder_sheet(false, None),
    ));
    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                operation: Operation::Tessellation,
                input: NO_PCURVE,
                ..
            }
        ),
        "{error:?}"
    );
}

#[test]
fn a_revolution_over_a_polyline_is_refused_by_name() {
    let support = |b: &mut GeometryGraphBuilder| {
        let polyline = Curve3::Polyline(Polyline3 {
            points: vec![Point3::new(0.3, 0.0, 0.0), Point3::new(0.3, 0.0, 1.0)],
            closed: false,
        });
        revolution(b, polyline, Point3::ZERO, Vec3::Z)
    };
    let error = refused(open_sheet(&support, &cylinder_sheet(true, None)));
    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                operation: Operation::SurfaceEvaluation,
                input: "a swept curve other than a line, circle, ellipse or B-spline",
                ..
            }
        ),
        "{error:?}"
    );
}

#[test]
fn a_revolution_about_a_zero_axis_is_refused() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(
            b,
            line3(Point3::new(0.3, 0.0, 0.0), Vec3::Z),
            Point3::ZERO,
            Vec3::ZERO,
        )
    };
    let error = refused(open_sheet(&support, &cylinder_sheet(true, None)));
    assert!(matches!(error, GeomError::InvalidInput(_)), "{error:?}");
}
