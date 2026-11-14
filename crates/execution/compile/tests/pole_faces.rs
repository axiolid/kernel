//! B-rep faces that reach a pole of their surface (#463, ADR 0112).
//!
//! A pole is a vertex, with no edge along the pole line (ADR 0059): the
//! tessellator closes the face's parameter ring there with a row of samples
//! welded to the pole's vertex, takes an edge's `u` at the pole from its
//! neighbouring sample, and bounds a grid's pole row with one fan per cell.
//! Each case is built on atomic surfaces and on the relations Revit and IFC
//! write (a revolution reaching its axis, a ruled surface with an apex),
//! with and without pcurves, and read back: closure, watertightness, the
//! volume against its closed form and a dense sample of the exact surface
//! against the reported bound. Faces without pcurves may report above the
//! budget, as ADR 0110 records for any face without pcurves, but never
//! below a measured deviation.

mod common;

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame3, Interval, Point3, Scalar, Vec3};
use axiolid_curve::{Curve3, KnotSpec};
use axiolid_mesh::TriMesh;
use axiolid_mesh_compile::{DeviationBound, DeviationReport};
use axiolid_mesh_compile_contract::MeshClosure;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId, RuledSide, SurfaceRelation};
use axiolid_surface::{BSplineSurface, Cone, Plane, Sphere, Surface};
use axiolid_topology::{
    BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};

use common::*;

const R: Scalar = 0.5;
const H: Scalar = 0.8;
const REVOLUTION: &str = "surface of revolution";

/// The reported bound and the triangle count of a closed solid that must
/// come to `volume`, and every point of `samples` within the bound. The
/// bound is held within the budget when `budgeted`, else within twice it.
fn solid(
    name: &str,
    result: GeomResult<Compiled>,
    family: &str,
    exact: (Scalar, Scalar),
    samples: &[Point3],
    budgeted: bool,
) -> (Scalar, usize) {
    solid_within(name, result, family, exact, (samples, 0.0), budgeted)
}

/// `solid`, each sample allowed `rounding` for its own: a sample computed
/// in floating point lies off a plane its exact point is on by an ulp or
/// so, which an exact face's bound of rounding size cannot absorb.
fn solid_within(
    name: &str,
    result: GeomResult<Compiled>,
    family: &str,
    (volume, area): (Scalar, Scalar),
    (samples, rounding): (&[Point3], Scalar),
    budgeted: bool,
) -> (Scalar, usize) {
    let (mesh, closure, report) = result.unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_eq!(closure, MeshClosure::Solid, "{name}");
    assert_watertight(name, &mesh);
    let bound = held(name, &report, family, budgeted);
    let got = signed_volume(&mesh);
    assert!(
        (got - volume).abs() <= bound * area,
        "{name}: volume {got} against {volume}, slack {}",
        bound * area
    );
    let worst = farthest(&mesh, samples);
    assert!(
        worst <= bound * (1.0 + 1e-9) + rounding,
        "{name}: a surface point lies {worst:e} from the mesh, above the reported {bound:e}"
    );
    let triangles = mesh.indices.len() / 3;
    eprintln!(
        "{name}: {triangles} triangles, {:.2?}, certified {bound:e}, measured {worst:e}",
        LAST_COMPILE.get()
    );
    (bound, triangles)
}

fn held(name: &str, report: &DeviationReport, family: &str, budgeted: bool) -> Scalar {
    if budgeted {
        certified(name, report, family)
    } else {
        let bound = certified_bound(name, report, family);
        assert!(
            bound <= 2.0 * BUDGET,
            "{name}: bound {bound:e} past twice the budget"
        );
        bound
    }
}

/// An open sheet's bound, every point of `samples` within it.
fn sheet(
    name: &str,
    result: GeomResult<Compiled>,
    family: &str,
    samples: &[Point3],
    budgeted: bool,
) -> (Scalar, TriMesh) {
    let (mesh, closure, report) = result.unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_eq!(closure, MeshClosure::Surface, "{name}");
    assert_manifold(name, &mesh);
    let bound = held(name, &report, family, budgeted);
    let worst = farthest(&mesh, samples);
    assert!(
        worst <= bound * (1.0 + 1e-9),
        "{name}: a surface point lies {worst:e} from the mesh, above the reported {bound:e}"
    );
    eprintln!(
        "{name}: {} triangles, {:.2?}, certified {bound:e}, measured {worst:e}",
        mesh.indices.len() / 3,
        LAST_COMPILE.get()
    );
    (bound, mesh)
}

/// No mesh edge of a sheet is used twice the same way: an open sheet is a
/// two-manifold with boundary, consistently wound.
fn assert_manifold(name: &str, mesh: &TriMesh) {
    let mut directed = std::collections::HashSet::new();
    for t in mesh.indices.as_chunks::<3>().0 {
        for (x, y) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            assert!(
                directed.insert((x, y)),
                "{name}: edge {x}-{y} used twice one way"
            );
        }
    }
}

/// A dense sample of `point` over `[0, 1]^2`, its border included.
fn dense(point: &dyn Fn(Scalar, Scalar) -> Point3) -> Vec<Point3> {
    let mut samples = Vec::new();
    let (n, m) = (97, 33);
    for i in 0..=n {
        for j in 0..=m {
            samples.push(point(
                (i as Scalar + 0.37) / (n as Scalar + 0.74),
                (j as Scalar + 0.41) / (m as Scalar + 0.82),
            ));
        }
    }
    for k in 0..=200 {
        let t = k as Scalar / 200.0;
        for s in [0.0, 1.0] {
            samples.push(point(t, s));
            samples.push(point(s, t));
        }
    }
    samples
}

/// A dense sample of the disc of `radius` at height `z` about the axis.
fn disc(radius: Scalar, z: Scalar) -> Vec<Point3> {
    dense(&|t, s| {
        let (angle, r) = (TAU * t, radius * s);
        Point3::new(r * angle.cos(), r * angle.sin(), z)
    })
}

fn relation(builder: &mut GeometryGraphBuilder, node: SurfaceRelation) -> NodeId {
    builder.push(GeometryNode::SurfaceRelation(node)).unwrap()
}

fn sphere(builder: &mut GeometryGraphBuilder) -> NodeId {
    atomic(
        builder,
        Surface::Sphere(Sphere {
            frame: frame_at(0.0),
            radius: R,
        }),
    )
}

fn arc_revolution(builder: &mut GeometryGraphBuilder) -> NodeId {
    revolution(builder, meridian_arc(R), Point3::ZERO, Vec3::Z)
}

/// The quarter of the meridian from the equator to the north pole, one
/// rational quadratic over `[0, 1]`: no knot, so no tangent break of the
/// profile (#498) reaches the hemisphere's edges.
fn quarter_revolution(builder: &mut GeometryGraphBuilder) -> NodeId {
    let w = core::f64::consts::FRAC_1_SQRT_2;
    let quarter = Curve3::BSpline(axiolid_curve::BSplineCurve3 {
        degree: 2,
        control_points: vec![
            Point3::new(R, 0.0, 0.0),
            Point3::new(R, 0.0, R),
            Point3::new(0.0, 0.0, R),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![3, 3],
        weights: Some(vec![1.0, w, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    revolution(builder, quarter, Point3::ZERO, Vec3::Z)
}

/// The quarter meridian's point at its ends, where the vertices lie.
fn quarter_point(u: Scalar, v: Scalar) -> Point3 {
    sphere_point(R, u, v * FRAC_PI_2)
}

/// The latitude of the B-spline meridian's parameter at its knots, where
/// the fixtures place their vertices.
fn arc_point(u: Scalar, v: Scalar) -> Point3 {
    sphere_point(R, u, (v - 0.5) * PI)
}

fn latitude_point(u: Scalar, v: Scalar) -> Point3 {
    sphere_point(R, u, v)
}

/// A full turn from the equator to the north pole, its seam the meridian
/// at `u = 0` and its floor a disc.
fn hemisphere<'a>(
    v: (Scalar, Scalar),
    point: &'a dyn Fn(Scalar, Scalar) -> Point3,
    pcurves: bool,
) -> PoleSheet<'a> {
    PoleSheet {
        u: (0.0, TAU),
        v,
        point,
        poles: [false, true],
        curves: [
            Some(horizontal_circle(0.0, 0.0, 0.0, R)),
            Some(circle_in(Point3::ZERO, Vec3::X, Vec3::Z, R)),
            None,
            None,
        ],
        pcurves,
        seam: true,
    }
}

fn hemisphere_samples() -> Vec<Point3> {
    let mut samples = dense(&|t, s| sphere_point(R, TAU * t, FRAC_PI_2 * s));
    samples.extend(disc(R, 0.0));
    samples
}

const HEMISPHERE: (Scalar, Scalar) = (2.0 / 3.0 * PI * R * R * R, 3.0 * PI * R * R);

#[test]
fn a_hemisphere_closes_at_its_pole() {
    for pcurves in [true, false] {
        let name = format!("hemisphere by revolution, pcurves {pcurves}");
        let result = pole_solid(
            &quarter_revolution,
            &hemisphere((0.0, 1.0), &quarter_point, pcurves),
        );
        solid(
            &name,
            result,
            REVOLUTION,
            HEMISPHERE,
            &hemisphere_samples(),
            false,
        );
        let name = format!("hemisphere on a sphere, pcurves {pcurves}");
        let result = pole_solid(
            &sphere,
            &hemisphere((0.0, FRAC_PI_2), &latitude_point, pcurves),
        );
        solid(
            &name,
            result,
            "sphere",
            HEMISPHERE,
            &hemisphere_samples(),
            pcurves,
        );
    }
}

/// A whole sphere with one seam: two poles, the full-turn rule (I3).
#[test]
fn a_sphere_with_one_seam_closes_at_both_poles() {
    let samples = dense(&|t, s| sphere_point(R, TAU * t, PI * (s - 0.5)));
    let exact = (4.0 / 3.0 * PI * R * R * R, 4.0 * PI * R * R);
    for pcurves in [true, false] {
        for (name, support, v, point) in [
            (
                "sphere by revolution",
                &arc_revolution as Support<'_>,
                (0.0, 1.0),
                &arc_point as &dyn Fn(Scalar, Scalar) -> Point3,
            ),
            ("sphere", &sphere, (-FRAC_PI_2, FRAC_PI_2), &latitude_point),
        ] {
            let sheet = PoleSheet {
                u: (0.0, TAU),
                v,
                point,
                poles: [true, true],
                curves: [
                    None,
                    Some(circle_in(Point3::ZERO, Vec3::X, Vec3::Z, R)),
                    None,
                    None,
                ],
                pcurves,
                seam: true,
            };
            let family = if name == "sphere" {
                "sphere"
            } else {
                REVOLUTION
            };
            let name = format!("{name}, pcurves {pcurves}");
            let result = pole_solid(support, &sheet);
            if family == REVOLUTION && !pcurves {
                // The semicircle's knot at the equator is a tangent break of
                // the profile, which the seam crosses: no span across it is
                // bounded (#498), so the seam stays coarse and its spans
                // into the poles cannot be shown regular. A solid, unbounded
                // by name.
                let (mesh, closure, report) = result.unwrap_or_else(|e| panic!("{name}: {e}"));
                assert_eq!(closure, MeshClosure::Solid, "{name}");
                assert_watertight(&name, &mesh);
                assert!(
                    report.contributions.iter().any(|c| matches!(c.bound,
                        DeviationBound::Unbounded(reason) if reason.contains("into a pole")
                            || reason.contains("not certifiably invertible"))),
                    "{name}: {report:?}"
                );
                continue;
            }
            solid(&name, result, family, exact, &samples, pcurves);
        }
    }
}

/// Two halves of a sphere, each between the two meridians: with pcurves,
/// the second half's far meridian sits at `u = 2 pi`, which names its side;
/// without, the side is not determined and the face is refused (P2).
fn lunes(support: Support<'_>, v: (Scalar, Scalar), pcurves: bool) -> GeomResult<Compiled> {
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let support = support(&mut builder);
    let south = brep.add_vertex(Vertex {
        position: Point3::new(0.0, 0.0, -R),
    });
    let north = brep.add_vertex(Vertex {
        position: Point3::new(0.0, 0.0, R),
    });
    let mut meridian = |x: Vec3| {
        let curve = push_curve3(&mut builder, circle_in(Point3::ZERO, x, Vec3::Z, R));
        brep.add_edge(Edge {
            start: south,
            end: north,
            curve: Some(curve),
        })
    };
    let (m0, m1) = (meridian(Vec3::X), meridian(-Vec3::X));
    let (v0, v1) = v;
    let mut half = |(u0, u1): (Scalar, Scalar), right, left| {
        let mut pc = |points: &[(Scalar, Scalar)]| {
            pcurves.then(|| {
                builder
                    .push(GeometryNode::Curve2(polyline(points)))
                    .unwrap()
            })
        };
        let edges = vec![
            EdgeUse {
                edge: right,
                orientation: Orientation::Forward,
                pcurve: pc(&[(u1, v0), (u1, v1)]),
            },
            EdgeUse {
                edge: left,
                orientation: Orientation::Reversed,
                pcurve: pc(&[(u0, v1), (u0, v0)]),
            },
        ];
        let loop_id = brep.add_loop(Loop { edges });
        brep.add_face(Face {
            surface: Some(support),
            bounds: vec![FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation: Orientation::Forward,
        })
    };
    let faces = [half((0.0, PI), m1, m0), half((PI, TAU), m0, m1)];
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

#[test]
fn two_halves_of_a_sphere_take_their_side_from_their_pcurves() {
    let samples = dense(&|t, s| sphere_point(R, TAU * t, PI * (s - 0.5)));
    let exact = (4.0 / 3.0 * PI * R * R * R, 4.0 * PI * R * R);
    for (name, support, v) in [
        (
            "lunes by revolution",
            &arc_revolution as Support<'_>,
            (0.0, 1.0),
        ),
        ("lunes on a sphere", &sphere, (-FRAC_PI_2, FRAC_PI_2)),
    ] {
        let family = if name.ends_with("sphere") {
            "sphere"
        } else {
            REVOLUTION
        };
        solid(name, lunes(support, v, true), family, exact, &samples, true);
        match lunes(support, v, false) {
            Err(GeomError::UnsupportedInput { input, .. }) => assert!(
                input.contains("between two poles whose edge uses carry no pcurves"),
                "{name}: {input}"
            ),
            other => panic!("{name}: expected P2, got {:?}", other.map(|c| c.2)),
        }
    }
}

fn cone_point(u: Scalar, v: Scalar) -> Point3 {
    let r = R * (1.0 - v);
    Point3::new(r * u.cos(), r * u.sin(), H * v)
}

/// The atomic cone's `v` is its height.
fn atomic_cone_point(u: Scalar, v: Scalar) -> Point3 {
    cone_point(u, v / H)
}

fn cone(builder: &mut GeometryGraphBuilder) -> NodeId {
    atomic(
        builder,
        Surface::Cone(Cone {
            frame: frame_at(0.0),
            radius: R,
            semi_angle: -(R / H).atan(),
        }),
    )
}

fn line_revolution(builder: &mut GeometryGraphBuilder) -> NodeId {
    revolution(
        builder,
        line3(Point3::new(R, 0.0, 0.0), Vec3::new(-R, 0.0, H)),
        Point3::ZERO,
        Vec3::Z,
    )
}

#[test]
fn a_cone_closes_at_its_apex() {
    let mut samples = dense(&|t, s| cone_point(TAU * t, s));
    samples.extend(disc(R, 0.0));
    let slant = (R * R + H * H).sqrt();
    let exact = (PI * R * R * H / 3.0, PI * R * slant + PI * R * R);
    for pcurves in [true, false] {
        for (name, support, v, point) in [
            (
                "cone by revolution",
                &line_revolution as Support<'_>,
                1.0,
                &cone_point as &dyn Fn(Scalar, Scalar) -> Point3,
            ),
            ("cone", &cone, H, &atomic_cone_point),
        ] {
            let sheet = PoleSheet {
                u: (0.0, TAU),
                v: (0.0, v),
                point,
                poles: [false, true],
                curves: [Some(horizontal_circle(0.0, 0.0, 0.0, R)), None, None, None],
                pcurves,
                seam: true,
            };
            let family = if name == "cone" { "cone" } else { REVOLUTION };
            let name = format!("{name}, pcurves {pcurves}");
            solid(
                &name,
                pole_solid(support, &sheet),
                family,
                exact,
                &samples,
                pcurves,
            );
        }
    }
}

/// A cone's wall one cell high is one fan per rim span: no triangle more,
/// and its bound is the rim's sagitta bound, `R du^2 / 8` (I7).
#[test]
fn a_cone_wall_one_cell_high_is_one_fan_per_rim_span() {
    let sheet = PoleSheet {
        u: (0.0, TAU),
        v: (0.0, H),
        point: &atomic_cone_point,
        poles: [false, true],
        curves: [Some(horizontal_circle(0.0, 0.0, 0.0, R)), None, None, None],
        pcurves: true,
        seam: true,
    };
    let samples = dense(&|t, s| cone_point(TAU * t, s));
    let (bound, mesh) = sheet_case("cone wall", open_pole_sheet(&cone, &sheet), &samples);
    let spans = mesh.indices.len() / 3;
    let du = TAU / spans as Scalar;
    let sagitta = R * du * du / 8.0;
    assert!(
        bound <= sagitta * (1.0 + 1e-9) && bound >= 0.5 * sagitta,
        "bound {bound:e} against the rim's {sagitta:e} over {spans} spans"
    );
    // The rim's own samples, no more: one fan each, all through the apex.
    let apex = Point3::new(0.0, 0.0, H);
    assert!(mesh.positions.contains(&apex));
    assert_eq!(mesh.positions.len(), spans + 1, "{spans} fans");
}

fn sheet_case(name: &str, result: GeomResult<Compiled>, samples: &[Point3]) -> (Scalar, TriMesh) {
    sheet(name, result, "cone", samples, true)
}

fn ruled_apex(
    builder: &mut GeometryGraphBuilder,
    curve: Curve3,
    interval: (Scalar, Scalar),
) -> NodeId {
    let curve = push_curve3(builder, curve);
    relation(
        builder,
        SurfaceRelation::Ruled {
            first: RuledSide::Curve {
                curve,
                interval: Interval {
                    start: interval.0,
                    end: interval.1,
                },
            },
            second: RuledSide::Apex(Point3::new(0.0, 0.0, H)),
        },
    )
}

/// Three quarters of a cone, ruled from a circle's arc to an apex: an
/// open sheet whose top is the apex.
#[test]
fn a_ruled_cone_closes_at_its_apex() {
    let turn = 1.5 * PI;
    let support = |b: &mut GeometryGraphBuilder| {
        ruled_apex(b, horizontal_circle(0.0, 0.0, 0.0, R), (0.0, turn))
    };
    let point = |u: Scalar, v: Scalar| cone_point(turn * u, v);
    let samples = dense(&|t, s| cone_point(turn * t, s));
    for pcurves in [true, false] {
        let sheet_ = PoleSheet {
            u: (0.0, 1.0),
            v: (0.0, 1.0),
            point: &point,
            poles: [false, true],
            curves: [Some(horizontal_circle(0.0, 0.0, 0.0, R)), None, None, None],
            pcurves,
            seam: false,
        };
        let name = format!("ruled cone, pcurves {pcurves}");
        sheet(
            &name,
            open_pole_sheet(&support, &sheet_),
            RULED,
            &samples,
            pcurves,
        );
    }
}

/// A flat ruled face from a segment to an apex is its own triangle: one
/// fan, bound 0.
#[test]
fn a_flat_ruled_face_into_an_apex_is_one_triangle() {
    // In the plane `y = 0`, so every sample lies in its triangle's plane
    // exactly and is measured without rounding.
    let (a, b) = (Point3::new(-0.4, 0.0, 0.0), Point3::new(0.5, 0.0, 0.0));
    let support =
        |builder: &mut GeometryGraphBuilder| ruled_apex(builder, line3(a, b - a), (0.0, 1.0));
    let point = |u: Scalar, v: Scalar| Point3::new((a.x + (b.x - a.x) * u) * (1.0 - v), 0.0, H * v);
    for pcurves in [true, false] {
        let sheet_ = PoleSheet {
            u: (0.0, 1.0),
            v: (0.0, 1.0),
            point: &point,
            poles: [false, true],
            curves: [None, None, None, None],
            pcurves,
            seam: false,
        };
        let name = format!("flat ruled fan, pcurves {pcurves}");
        // Inside the triangle, off its edges, where the in-plane distance is
        // read exactly.
        let samples: Vec<Point3> = (1..60)
            .flat_map(|i| (1..60).map(move |j| (i as Scalar / 60.0, j as Scalar / 60.0)))
            .map(|(u, v)| point(u, v))
            .collect();
        let (bound, mesh) = sheet(
            &name,
            open_pole_sheet(&support, &sheet_),
            RULED,
            &samples,
            true,
        );
        assert_eq!(bound, 0.0, "{name}");
        assert_eq!(mesh.indices.len(), 3, "{name}: one triangle");
    }
}

/// A pyramid whose walls are ruled from each side of its base to its apex,
/// the base a star polygon (pleated, as a luminaire's shade), or a regular
/// one with each side an arc bulging outward. The flat walls are exact,
/// one triangle each.
fn pleated(corners: &[Point3], arcs: Option<Scalar>, pcurves: bool) -> GeomResult<Compiled> {
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let n = corners.len();
    let apex_point = Point3::new(0.0, 0.0, H);
    let apex = brep.add_vertex(Vertex {
        position: apex_point,
    });
    let base: Vec<_> = corners
        .iter()
        .map(|&position| brep.add_vertex(Vertex { position }))
        .collect();
    // The side from corner i to i + 1, as a curve over `[0, 1]` (a line) or
    // an arc of its circle.
    let side = |i: usize| -> (Curve3, (Scalar, Scalar)) {
        let (p, q) = (corners[i], corners[(i + 1) % n]);
        match arcs {
            None => (line3(p, q - p), (0.0, 1.0)),
            Some(sagitta) => {
                let (chord, middle) = ((q - p).length(), (p + q) * 0.5);
                let out = Vec3::new(middle.x, middle.y, 0.0).normalize();
                let radius = (chord * chord / 4.0 + sagitta * sagitta) / (2.0 * sagitta);
                let centre = middle + out * (sagitta - radius);
                let x = (p - centre).normalize();
                let y = Vec3::Z.cross(x);
                let sweep = 2.0 * (0.5 * chord / radius).asin();
                (circle_in(centre, x, y, radius), (0.0, sweep))
            }
        }
    };
    let sides: Vec<_> = (0..n)
        .map(|i| {
            let (curve, _) = side(i);
            let curve = push_curve3(&mut builder, curve);
            brep.add_edge(Edge {
                start: base[i],
                end: base[(i + 1) % n],
                curve: Some(curve),
            })
        })
        .collect();
    let rulings: Vec<_> = (0..n)
        .map(|i| {
            brep.add_edge(Edge {
                start: base[i],
                end: apex,
                curve: None,
            })
        })
        .collect();
    let mut faces = Vec::new();
    for i in 0..n {
        let (curve, interval) = side(i);
        let support = ruled_apex(&mut builder, curve, interval);
        let mut pc = |points: &[(Scalar, Scalar)]| {
            pcurves.then(|| {
                builder
                    .push(GeometryNode::Curve2(polyline(points)))
                    .unwrap()
            })
        };
        let edges = vec![
            EdgeUse {
                edge: sides[i],
                orientation: Orientation::Forward,
                pcurve: pc(&[(0.0, 0.0), (1.0, 0.0)]),
            },
            EdgeUse {
                edge: rulings[(i + 1) % n],
                orientation: Orientation::Forward,
                pcurve: pc(&[(1.0, 0.0), (1.0, 1.0)]),
            },
            EdgeUse {
                edge: rulings[i],
                orientation: Orientation::Reversed,
                pcurve: pc(&[(0.0, 1.0), (0.0, 0.0)]),
            },
        ];
        let loop_id = brep.add_loop(Loop { edges });
        faces.push(brep.add_face(Face {
            surface: Some(support),
            bounds: vec![FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation: Orientation::Forward,
        }));
    }
    // The base faces down: its loop runs anticlockwise about its plane's
    // normal, `-z`.
    let floor = atomic(
        &mut builder,
        Surface::Plane(Plane {
            frame: Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: -Vec3::Y,
                z: -Vec3::Z,
            },
        }),
    );
    let edges = (0..n)
        .rev()
        .map(|i| EdgeUse {
            edge: sides[i],
            orientation: Orientation::Reversed,
            pcurve: None,
        })
        .collect();
    let loop_id = brep.add_loop(Loop { edges });
    faces.push(brep.add_face(Face {
        surface: Some(floor),
        bounds: vec![FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    }));
    let shell = brep.add_shell(Shell {
        faces: faces
            .into_iter()
            .map(|f| (f, Orientation::Forward))
            .collect(),
        closed: true,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    compile(builder, brep)
}

fn polygon_area(corners: &[Point3]) -> Scalar {
    let n = corners.len();
    0.5 * (0..n)
        .map(|i| {
            let (p, q) = (corners[i], corners[(i + 1) % n]);
            p.x * q.y - q.x * p.y
        })
        .sum::<Scalar>()
}

#[test]
fn a_pleated_shade_ruled_to_its_apex_is_exact() {
    let star: Vec<Point3> = (0..10)
        .map(|k| {
            let (angle, r) = (
                TAU * k as Scalar / 10.0,
                if k % 2 == 0 { R } else { 0.45 * R },
            );
            Point3::new(r * angle.cos(), r * angle.sin(), 0.0)
        })
        .collect();
    let apex = Point3::new(0.0, 0.0, H);
    let mut samples = Vec::new();
    for i in 0..star.len() {
        let (p, q) = (star[i], star[(i + 1) % star.len()]);
        samples.extend(dense(&|t, s| (p + (q - p) * t) * (1.0 - s) + apex * s));
    }
    let area = polygon_area(&star);
    let wall: Scalar = (0..star.len())
        .map(|i| {
            let (p, q) = (star[i], star[(i + 1) % star.len()]);
            0.5 * (q - p).cross(apex - p).length()
        })
        .sum();
    for pcurves in [true, false] {
        let name = format!("pleated shade, pcurves {pcurves}");
        let rounding = 4.0 * Scalar::EPSILON;
        let (bound, triangles) = solid_within(
            &name,
            pleated(&star, None, pcurves),
            RULED,
            (area * H / 3.0, area + wall),
            (&samples, rounding),
            true,
        );
        // Exact up to the rounding of the ruled surface's own points.
        assert!(bound <= rounding, "{name}: {bound:e}");
        // One triangle per wall, and the star's ten corners cut into eight.
        assert_eq!(triangles, 10 + 8, "{name}");
    }
}

#[test]
fn a_shade_with_arc_pleats_closes_at_its_apex() {
    let n = 6;
    let hexagon: Vec<Point3> = (0..n)
        .map(|k| {
            let angle = TAU * k as Scalar / n as Scalar;
            Point3::new(R * angle.cos(), R * angle.sin(), 0.0)
        })
        .collect();
    let sagitta = 0.06;
    let chord = R;
    let radius = (chord * chord / 4.0 + sagitta * sagitta) / (2.0 * sagitta);
    let theta = 2.0 * (0.5 * chord / radius).asin();
    let segment = 0.5 * radius * radius * (theta - theta.sin());
    let base = polygon_area(&hexagon) + n as Scalar * segment;
    let apex = Point3::new(0.0, 0.0, H);
    let mut samples = Vec::new();
    for i in 0..n {
        let (p, q) = (hexagon[i], hexagon[(i + 1) % n]);
        let middle = (p + q) * 0.5;
        let out = Vec3::new(middle.x, middle.y, 0.0).normalize();
        let centre = middle + out * (sagitta - radius);
        let x = (p - centre).normalize();
        let y = Vec3::Z.cross(x);
        let arc =
            move |t: Scalar| centre + (x * (theta * t).cos() + y * (theta * t).sin()) * radius;
        samples.extend(dense(&|t, s| arc(t) * (1.0 - s) + apex * s));
    }
    // The base's area bounds the walls' from below; the slack only scales.
    let area = 3.0 * base;
    solid(
        "arc pleats",
        pleated(&hexagon, Some(sagitta), false),
        RULED,
        (base * H / 3.0, area),
        &samples,
        false,
    );
    // With pcurves on the walls only, the base takes its arcs' samples
    // through them and has no chord bound on them (ADR 0086): the walls
    // are certified, the base is unbounded by that name.
    let (mesh, closure, report) = pleated(&hexagon, Some(sagitta), true).unwrap();
    assert_eq!(closure, MeshClosure::Solid);
    assert_watertight("arc pleats with pcurves", &mesh);
    for contribution in &report.contributions {
        match contribution.detail {
            RULED => assert!(
                matches!(contribution.bound, DeviationBound::Certified(b) if b <= BUDGET),
                "{contribution:?}"
            ),
            _ => assert!(
                matches!(contribution.bound, DeviationBound::Unbounded(reason)
                    if reason.contains("sampled through a neighbour's pcurve")),
                "{contribution:?}"
            ),
        }
    }
}

/// A canopy profile, a cubic from the axis at the bottom to the axis at
/// the top, bulging out.
fn canopy_points() -> [Point3; 4] {
    [
        Point3::ZERO,
        Point3::new(0.6, 0.0, 0.1),
        Point3::new(0.5, 0.0, 0.9),
        Point3::new(0.0, 0.0, 1.0),
    ]
}

fn casteljau(points: &[Point3], t: Scalar) -> Point3 {
    let mut p = points.to_vec();
    while p.len() > 1 {
        p = p.windows(2).map(|w| w[0].lerp(w[1], t)).collect();
    }
    p[0]
}

fn turned(points: &[Point3], u: Scalar, v: Scalar) -> Point3 {
    let p = casteljau(points, v);
    Point3::new(p.x * u.cos(), p.x * u.sin(), p.z)
}

/// The volume and area of the profile turned once, by Simpson's rule.
fn turned_measure(points: &[Point3]) -> (Scalar, Scalar) {
    let n = 4000;
    let (mut volume, mut area) = (0.0, 0.0);
    for i in 0..=n {
        let t = i as Scalar / n as Scalar;
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let h = 1e-6;
        let d = (casteljau(points, (t + h).min(1.0)) - casteljau(points, (t - h).max(0.0)))
            / ((t + h).min(1.0) - (t - h).max(0.0));
        let p = casteljau(points, t);
        volume += w * PI * p.x * p.x * d.z;
        area += w * TAU * p.x * d.length();
    }
    (volume / (3.0 * n as Scalar), area / (3.0 * n as Scalar))
}

#[test]
fn a_canopy_from_the_axis_to_the_axis_closes_at_both_poles() {
    let points = canopy_points();
    let support =
        |b: &mut GeometryGraphBuilder| revolution(b, bezier(points), Point3::ZERO, Vec3::Z);
    let point = |u: Scalar, v: Scalar| turned(&points, u, v);
    let samples = dense(&|t, s| turned(&points, TAU * t, s));
    let exact = turned_measure(&points);
    for pcurves in [true, false] {
        let sheet_ = PoleSheet {
            u: (0.0, TAU),
            v: (0.0, 1.0),
            point: &point,
            poles: [true, true],
            curves: [None, Some(bezier(points)), None, None],
            pcurves,
            seam: true,
        };
        let name = format!("canopy, one seam, pcurves {pcurves}");
        solid(
            &name,
            pole_solid(&support, &sheet_),
            REVOLUTION,
            exact,
            &samples,
            false,
        );
    }
    // Two half turns between the profile and its mirror, with pcurves.
    let mirror = points.map(|p| Point3::new(-p.x, p.y, p.z));
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let support = revolution(&mut builder, bezier(points), Point3::ZERO, Vec3::Z);
    let [bottom, top] = [points[0], points[3]].map(|position| brep.add_vertex(Vertex { position }));
    let mut meridian = |curve| {
        let curve = push_curve3(&mut builder, curve);
        brep.add_edge(Edge {
            start: bottom,
            end: top,
            curve: Some(curve),
        })
    };
    let (m0, m1) = (meridian(bezier(points)), meridian(bezier(mirror)));
    let mut half = |(u0, u1): (Scalar, Scalar), right, left| {
        let mut pc = |points: &[(Scalar, Scalar)]| {
            builder
                .push(GeometryNode::Curve2(polyline(points)))
                .unwrap()
        };
        let edges = vec![
            EdgeUse {
                edge: right,
                orientation: Orientation::Forward,
                pcurve: Some(pc(&[(u1, 0.0), (u1, 1.0)])),
            },
            EdgeUse {
                edge: left,
                orientation: Orientation::Reversed,
                pcurve: Some(pc(&[(u0, 1.0), (u0, 0.0)])),
            },
        ];
        let loop_id = brep.add_loop(Loop { edges });
        brep.add_face(Face {
            surface: Some(support),
            bounds: vec![FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation: Orientation::Forward,
        })
    };
    let faces = [half((0.0, PI), m1, m0), half((PI, TAU), m0, m1)];
    let shell = brep.add_shell(Shell {
        faces: faces.map(|face| (face, Orientation::Forward)).to_vec(),
        closed: true,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    solid(
        "canopy, two half turns",
        compile(builder, brep),
        REVOLUTION,
        exact,
        &samples,
        true,
    );
}

#[test]
fn a_canopy_from_a_rim_to_the_axis_closes_with_its_floor() {
    let points = [
        Point3::new(0.5, 0.0, 0.0),
        Point3::new(0.7, 0.0, 0.4),
        Point3::new(0.4, 0.0, 0.9),
        Point3::new(0.0, 0.0, 1.0),
    ];
    let support =
        |b: &mut GeometryGraphBuilder| revolution(b, bezier(points), Point3::ZERO, Vec3::Z);
    let point = |u: Scalar, v: Scalar| turned(&points, u, v);
    let mut samples = dense(&|t, s| turned(&points, TAU * t, s));
    samples.extend(disc(0.5, 0.0));
    let (volume, area) = turned_measure(&points);
    for pcurves in [true, false] {
        let sheet_ = PoleSheet {
            u: (0.0, TAU),
            v: (0.0, 1.0),
            point: &point,
            poles: [false, true],
            curves: [
                Some(horizontal_circle(0.0, 0.0, 0.0, 0.5)),
                Some(bezier(points)),
                None,
                None,
            ],
            pcurves,
            seam: true,
        };
        let name = format!("canopy on a floor, pcurves {pcurves}");
        solid(
            &name,
            pole_solid(&support, &sheet_),
            REVOLUTION,
            (volume, area + PI * 0.25),
            &samples,
            false,
        );
    }
}

/// A spherical triangle with its tip at the north pole, its far side a
/// great circle's arc: not a grid, so its pole row is meshed by the
/// clipper and bounded by the nominal triangles (I8). Reported as it is,
/// never below the measured deviation.
#[test]
fn a_teardrop_with_its_tip_at_a_pole_is_bounded_honestly() {
    let latitude = PI / 6.0;
    let (a, b) = (
        sphere_point(R, 0.0, latitude),
        sphere_point(R, FRAC_PI_2, latitude),
    );
    let north = Point3::new(0.0, 0.0, R);
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let support = sphere(&mut builder);
    let [va, vb, vn] = [a, b, north].map(|position| brep.add_vertex(Vertex { position }));
    let mut edge = |start, end, curve| {
        let curve = push_curve3(&mut builder, curve);
        brep.add_edge(Edge {
            start,
            end,
            curve: Some(curve),
        })
    };
    let x = a.normalize();
    let y = (b - x * b.dot(x)).normalize();
    let across = edge(va, vb, circle_in(Point3::ZERO, x, y, R));
    let up = edge(vb, vn, circle_in(Point3::ZERO, Vec3::Y, Vec3::Z, R));
    let down = edge(va, vn, circle_in(Point3::ZERO, Vec3::X, Vec3::Z, R));
    let uses = [
        (across, Orientation::Forward),
        (up, Orientation::Forward),
        (down, Orientation::Reversed),
    ]
    .map(|(edge, orientation)| EdgeUse {
        edge,
        orientation,
        pcurve: None,
    });
    let loop_id = brep.add_loop(Loop {
        edges: uses.to_vec(),
    });
    let face = brep.add_face(Face {
        surface: Some(support),
        bounds: vec![FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    let samples = dense(&|t, s| {
        let w = [(1.0 - t) * (1.0 - s), t * (1.0 - s), s];
        (a * w[0] + b * w[1] + north * w[2]).normalize() * R
    });
    let (mesh, _, report) = compile(builder, brep).expect("a teardrop at a pole meshes");
    assert_manifold("teardrop", &mesh);
    let bound = certified_bound("teardrop", &report, "sphere");
    let worst = farthest(&mesh, &samples);
    eprintln!(
        "teardrop: {} triangles, {:.2?}, certified {bound:e}, measured {worst:e}",
        mesh.indices.len() / 3,
        LAST_COMPILE.get()
    );
    assert!(
        worst <= bound,
        "teardrop: measured {worst:e} above {bound:e}"
    );
}

/// A face between a meridian and a small circle through the north pole,
/// which reaches the pole obliquely: the small circle's last span into the
/// pole is bounded along the meridians (I6), however far its preimage
/// swings in `u` there, and no surface point of the face lies farther from
/// the mesh than the report says.
#[test]
fn an_oblique_small_circle_into_a_pole_is_bounded_along_the_meridians() {
    let north = Point3::new(0.0, 0.0, R);
    let a = sphere_point(R, 0.0, PI / 6.0);
    let b = sphere_point(R, 0.9, 50f64.to_radians());
    let normal = (a - north).cross(b - north).normalize();
    let centre = normal * normal.dot(north);
    let radius = (R * R - centre.length_squared()).sqrt();
    let x = (a - centre).normalize();
    let angle = |y: Vec3, p: Point3| {
        let q = p - centre;
        q.dot(y).atan2(q.dot(x)).rem_euclid(TAU)
    };
    let y = if angle(normal.cross(x), b) < angle(normal.cross(x), north) {
        normal.cross(x)
    } else {
        x.cross(normal)
    };
    let (end, middle) = (angle(y, north), angle(y, b));
    assert!(middle < end, "{middle} {end}");
    let circle = |t: Scalar| centre + (x * t.cos() + y * t.sin()) * radius;
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let support = sphere(&mut builder);
    let [va, vn] = [a, north].map(|position| brep.add_vertex(Vertex { position }));
    let small = push_curve3(&mut builder, circle_in(centre, x, y, radius));
    let meridian = push_curve3(&mut builder, circle_in(Point3::ZERO, Vec3::X, Vec3::Z, R));
    let small = brep.add_edge(Edge {
        start: va,
        end: vn,
        curve: Some(small),
    });
    let meridian = brep.add_edge(Edge {
        start: va,
        end: vn,
        curve: Some(meridian),
    });
    let loop_id = brep.add_loop(Loop {
        edges: vec![
            EdgeUse {
                edge: small,
                orientation: Orientation::Forward,
                pcurve: None,
            },
            EdgeUse {
                edge: meridian,
                orientation: Orientation::Reversed,
                pcurve: None,
            },
        ],
    });
    let face = brep.add_face(Face {
        surface: Some(support),
        bounds: vec![FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    // Each point of the small arc and the parallel from it back to the
    // meridian: the face.
    let samples = dense(&|t, s| {
        let p = circle(end * t);
        let latitude = (p.z / R).clamp(-1.0, 1.0).asin();
        let longitude = p.y.atan2(p.x);
        sphere_point(R, s * longitude, latitude)
    });
    sheet(
        "oblique small circle",
        compile(builder, brep),
        "sphere",
        &samples,
        false,
    );
}

/// One edge use of a face: the edge, its orientation and its pcurve's ends.
type Use = (usize, Orientation, Option<[(Scalar, Scalar); 2]>);

/// A closed shell on the sphere: each face a loop of `(edge, orientation,
/// pcurve)` over edges between `vertices`, each edge a curve.
fn sphere_shell(
    vertices: &[Point3],
    edges: &[(usize, usize, Curve3)],
    faces: &[&[Use]],
) -> GeomResult<Compiled> {
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let support = sphere(&mut builder);
    let vertices: Vec<_> = vertices
        .iter()
        .map(|&position| brep.add_vertex(Vertex { position }))
        .collect();
    let edges: Vec<_> = edges
        .iter()
        .map(|(start, end, curve)| {
            let curve = push_curve3(&mut builder, curve.clone());
            brep.add_edge(Edge {
                start: vertices[*start],
                end: vertices[*end],
                curve: Some(curve),
            })
        })
        .collect();
    let faces: Vec<_> = faces
        .iter()
        .map(|uses| {
            let edges = uses
                .iter()
                .map(|&(edge, orientation, pcurve)| EdgeUse {
                    edge: edges[edge],
                    orientation,
                    pcurve: pcurve
                        .map(|ends| builder.push(GeometryNode::Curve2(polyline(&ends))).unwrap()),
                })
                .collect();
            let loop_id = brep.add_loop(Loop { edges });
            brep.add_face(Face {
                surface: Some(support),
                bounds: vec![FaceBound {
                    loop_id,
                    orientation: Orientation::Forward,
                    outer: true,
                }],
                orientation: Orientation::Forward,
            })
        })
        .collect();
    let shell = brep.add_shell(Shell {
        faces: faces
            .into_iter()
            .map(|f| (f, Orientation::Forward))
            .collect(),
        closed: true,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    compile(builder, brep)
}

/// Two distinct edges between the poles are no seam, even where they leave
/// both poles along each other (here they coincide): without pcurves the
/// side of the meridians each face covers is not determined (P2). Read as
/// a full turn, each face would be the whole sphere plus a sliver.
#[test]
fn two_edges_between_the_poles_are_no_full_turn() {
    let poles = [Point3::new(0.0, 0.0, -R), Point3::new(0.0, 0.0, R)];
    let meridian = || circle_in(Point3::ZERO, Vec3::X, Vec3::Z, R);
    let (f, r) = (Orientation::Forward, Orientation::Reversed);
    let result = sphere_shell(
        &poles,
        &[(0, 1, meridian()), (0, 1, meridian())],
        &[&[(1, f, None), (0, r, None)], &[(0, f, None), (1, r, None)]],
    );
    match result {
        Err(GeomError::UnsupportedInput { input, .. }) => assert!(
            input.contains("between two poles whose edge uses carry no pcurves"),
            "{input}"
        ),
        other => panic!("expected P2, got {:?}", other.map(|c| c.2)),
    }
}

/// One seam used once each way is a full turn whether or not each use
/// carries its pcurve. The use without one is inverted at the samples the
/// other's pcurve set, so it has no chord bound on its curve and the face
/// is unbounded by that name (ADR 0086), but a solid of the sphere's
/// volume.
#[test]
fn a_full_turn_takes_one_seam_use_with_a_pcurve_and_one_without() {
    let poles = [Point3::new(0.0, 0.0, -R), Point3::new(0.0, 0.0, R)];
    let samples = dense(&|t, s| sphere_point(R, TAU * t, PI * (s - 0.5)));
    let exact = (4.0 / 3.0 * PI * R * R * R, 4.0 * PI * R * R);
    let (f, r) = (Orientation::Forward, Orientation::Reversed);
    let up = Some([(TAU, -FRAC_PI_2), (TAU, FRAC_PI_2)]);
    let down = Some([(0.0, FRAC_PI_2), (0.0, -FRAC_PI_2)]);
    for (name, uses) in [
        ("pcurve up only", [(0, f, up), (0, r, None)]),
        ("pcurve down only", [(0, f, None), (0, r, down)]),
    ] {
        let result = sphere_shell(
            &poles,
            &[(0, 1, circle_in(Point3::ZERO, Vec3::X, Vec3::Z, R))],
            &[&uses],
        );
        let (mesh, closure, report) = result.unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(closure, MeshClosure::Solid, "{name}");
        assert_watertight(name, &mesh);
        assert!(
            report.contributions.iter().any(|c| matches!(c.bound,
                DeviationBound::Unbounded(reason) if reason.contains("sampled through a neighbour's pcurve"))),
            "{name}: {report:?}"
        );
        let volume = signed_volume(&mesh);
        assert!(
            (volume - exact.0).abs() <= 2.0 * BUDGET * exact.1,
            "{name}: volume {volume}"
        );
        let worst = farthest(&mesh, &samples);
        assert!(worst <= 2.0 * BUDGET, "{name}: {worst:e}");
    }
}

/// A face whose true region is the complement of a geodesic triangle with
/// a vertex at the north pole, with no seam, is invalid input (ADR 0059: a
/// face round a pole needs a seam), and its ring is the triangle's: the
/// complement is meshed as the triangle, wound as the face says, so the
/// shell of the two is the triangle twice the same way round, which the
/// closure check reports open. It used to be refused at the pole; now its
/// side is not read (#500 owns reading it from the loop's winding).
#[test]
fn a_geodesic_triangle_s_complement_without_a_seam_meshes_as_the_triangle() {
    let latitude = PI / 6.0;
    let (a, b) = (
        sphere_point(R, 0.0, latitude),
        sphere_point(R, FRAC_PI_2, latitude),
    );
    let north = Point3::new(0.0, 0.0, R);
    let x = a.normalize();
    let y = (b - x * b.dot(x)).normalize();
    let (f, r) = (Orientation::Forward, Orientation::Reversed);
    let (mesh, closure, _) = sphere_shell(
        &[a, b, north],
        &[
            (0, 1, circle_in(Point3::ZERO, x, y, R)),
            (1, 2, circle_in(Point3::ZERO, Vec3::Y, Vec3::Z, R)),
            (0, 2, circle_in(Point3::ZERO, Vec3::X, Vec3::Z, R)),
        ],
        &[
            &[(0, f, None), (1, f, None), (2, r, None)],
            &[(2, f, None), (1, r, None), (0, r, None)],
        ],
    )
    .expect("the complement is not detected");
    assert_eq!(closure, MeshClosure::OpenSolid);
    // Both faces lie on the triangle, whose great-circle side bends toward
    // the pole: no vertex lies below its corners.
    let floor = mesh
        .positions
        .iter()
        .map(|p| p.z)
        .fold(Scalar::INFINITY, Scalar::min);
    assert!(floor >= R * latitude.sin() - 1e-12, "{floor}");
}

/// The trims the tessellator refuses at a pole, each by name.
fn refusal(result: GeomResult<Compiled>) -> String {
    refused(result).to_string()
}

/// A hole whose trim reaches a pole (P4): the outer ring is a patch of the
/// sphere, the hole three edges round the north pole.
#[test]
fn a_hole_reaching_a_pole_is_refused() {
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let support = sphere(&mut builder);
    let mut pc = |points: &[(Scalar, Scalar)]| {
        builder
            .push(GeometryNode::Curve2(polyline(points)))
            .unwrap()
    };
    // Each edge's pcurve; the edges run from one pcurve's start to the next.
    let mut ring = |pcurves: &[[(Scalar, Scalar); 2]]| {
        let vertices: Vec<_> = pcurves
            .iter()
            .map(|&[(u, v), _]| {
                brep.add_vertex(Vertex {
                    position: sphere_point(R, u, v),
                })
            })
            .collect();
        let n = pcurves.len();
        let edges = (0..n)
            .map(|i| EdgeUse {
                edge: brep.add_edge(Edge {
                    start: vertices[i],
                    end: vertices[(i + 1) % n],
                    curve: None,
                }),
                orientation: Orientation::Forward,
                pcurve: Some(pc(&pcurves[i])),
            })
            .collect();
        brep.add_loop(Loop { edges })
    };
    let outer = ring(&[
        [(0.0, 0.2), (3.0, 0.2)],
        [(3.0, 0.2), (3.0, 1.2)],
        [(3.0, 1.2), (0.0, 1.2)],
        [(0.0, 1.2), (0.0, 0.2)],
    ]);
    let hole = ring(&[
        [(1.0, 1.3), (1.0, FRAC_PI_2)],
        [(2.0, FRAC_PI_2), (2.0, 1.3)],
        [(2.0, 1.3), (1.0, 1.3)],
    ]);
    let face = brep.add_face(Face {
        surface: Some(support),
        bounds: vec![
            FaceBound {
                loop_id: outer,
                orientation: Orientation::Forward,
                outer: true,
            },
            FaceBound {
                loop_id: hole,
                orientation: Orientation::Forward,
                outer: false,
            },
        ],
        orientation: Orientation::Forward,
    });
    brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    let text = refusal(compile(builder, brep));
    assert!(
        text.contains("a hole of a curved B-rep face whose trim reaches a pole"),
        "{text}"
    );
}

/// A B-spline surface whose top row of control points is one point: its
/// collapsed edge is no pole the tessellator closes (P1).
#[test]
fn a_collapsed_edge_of_a_spline_surface_is_refused() {
    let apex = Point3::new(0.0, 0.0, H);
    let net: Vec<Vec<Point3>> = (0..3)
        .map(|i| {
            let x = i as Scalar * 0.25;
            vec![Point3::new(x, 0.0, 0.0), Point3::new(x, 0.2, 0.3), apex]
        })
        .collect();
    let surface = Surface::BSpline(BSplineSurface {
        u_degree: 2,
        v_degree: 2,
        control_points: net,
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![3, 3],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![3, 3],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    });
    let point = |u: Scalar, v: Scalar| {
        if v == 1.0 {
            apex
        } else {
            Point3::new(0.5 * u, 0.0, 0.0)
        }
    };
    let support = |b: &mut GeometryGraphBuilder| atomic(b, surface.clone());
    let sheet_ = PoleSheet {
        u: (0.0, 1.0),
        v: (0.0, 1.0),
        point: &point,
        poles: [false, true],
        curves: [None, None, None, None],
        pcurves: true,
        seam: false,
    };
    let text = refusal(open_pole_sheet(&support, &sheet_));
    assert!(
        text.contains("closes across a collapsed edge of a surface that is not a sphere"),
        "{text}"
    );
}

/// A spherical cap bounded by one circle and no seam: its trim circles the
/// pole without reaching it. Refused as before (P8).
#[test]
fn a_cap_circling_a_pole_without_a_seam_is_refused_as_before() {
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let support = sphere(&mut builder);
    let latitude = PI / 4.0;
    let start = brep.add_vertex(Vertex {
        position: sphere_point(R, 0.0, latitude),
    });
    let curve = push_curve3(
        &mut builder,
        horizontal_circle(0.0, 0.0, R * latitude.sin(), R * latitude.cos()),
    );
    let rim = brep.add_edge(Edge {
        start,
        end: start,
        curve: Some(curve),
    });
    let loop_id = brep.add_loop(Loop {
        edges: vec![EdgeUse {
            edge: rim,
            orientation: Orientation::Forward,
            pcurve: None,
        }],
    });
    let face = brep.add_face(Face {
        surface: Some(support),
        bounds: vec![FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    let text = refusal(compile(builder, brep));
    eprintln!("cap without a seam: {text}");
    assert!(text.contains("curved"), "{text}");
}

/// An edge whose interior passes over a pole is refused there (P5, or P6
/// where no sample lands within the tolerance of it): the producer splits
/// it at the pole's vertex.
#[test]
fn an_edge_over_a_pole_off_its_samples_is_refused() {
    // From 50 degrees at u = 0 over the pole to 70 degrees at u = pi: the
    // pole lies at two thirds of the edge, on no dyadic sample.
    let (a, b) = (
        sphere_point(R, 0.0, 50f64.to_radians()),
        sphere_point(R, PI, 70f64.to_radians()),
    );
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let support = sphere(&mut builder);
    let [va, vb] = [a, b].map(|position| brep.add_vertex(Vertex { position }));
    let over = push_curve3(&mut builder, circle_in(Point3::ZERO, Vec3::X, Vec3::Z, R));
    let e0 = brep.add_edge(Edge {
        start: va,
        end: vb,
        curve: Some(over),
    });
    let loop_id = brep.add_loop(Loop {
        edges: vec![
            EdgeUse {
                edge: e0,
                orientation: Orientation::Forward,
                pcurve: None,
            },
            EdgeUse {
                edge: e0,
                orientation: Orientation::Reversed,
                pcurve: None,
            },
        ],
    });
    let face = brep.add_face(Face {
        surface: Some(support),
        bounds: vec![FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    let text = refusal(compile(builder, brep));
    assert!(
        text.contains("sphere pole") || text.contains("turns a quarter period"),
        "{text}"
    );
}
