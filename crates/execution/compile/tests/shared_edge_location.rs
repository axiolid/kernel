//! Shared edge samples located on each face's pcurve (#413).
//!
//! The B-rep edge pass samples an edge once, through one pcurve use. A
//! neighbouring use runs its own pcurve over the same edge, and nothing
//! makes it walk the edge at the same pace: a cap whose circle starts
//! elsewhere than the wall's rim, a polyline pcurve with an uneven vertex
//! spacing. Each case compiles a closed cylinder whose every edge carries a
//! pcurve on both its faces, and reads the mesh back: watertight, the
//! enclosed volume, and a dense sample of the rims and the wall against the
//! reported bound, which is certified or the case fails.

use std::collections::HashMap;
use std::f64::consts::TAU;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{BSplineCurve, Circle2, Circle3, Curve2, Curve3, KnotSpec, Polyline2};
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationReport, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshClosure;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId};
use axiolid_surface::{Cylinder, Plane, Surface};
use axiolid_topology::{
    BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};

const RADIUS: Scalar = 0.3;
const HEIGHT: Scalar = 0.5;
const BUDGET: Scalar = 1e-3;

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

/// How the rim and the caps of the closed cylinder are authored.
#[derive(Clone, Copy)]
struct Rim {
    /// The angle the rim's one vertex, the seam and both cap circles start at.
    start: Scalar,
    /// Where a rim's polyline pcurve has its middle vertex, as a fraction of
    /// the turn: `0.5` walks it at the caps' pace, anything else does not.
    middle: Scalar,
    /// How far above the rim the middle vertex of the bottom polyline sits.
    bump: Scalar,
    /// How far the caps' circles are wider than the wall.
    cap_growth: Scalar,
    /// Whether the caps' circles are pushed to the graph before the wall's
    /// pcurves, which makes a cap the use the edge is sampled through.
    caps_first: bool,
    /// The bottom rim as one quadratic B-spline along the rim, its middle
    /// control point at `middle` of the turn: collinear, and at a pace of
    /// its own, where the polyline bends.
    spline: bool,
}

impl Rim {
    const EVEN: Self = Self {
        start: 0.0,
        middle: 0.5,
        bump: 0.0,
        cap_growth: 0.0,
        caps_first: false,
        spline: false,
    };
}

fn polyline(points: &[(Scalar, Scalar)]) -> Curve2 {
    Curve2::Polyline(Polyline2 {
        points: points.iter().map(|&(u, v)| Point2::new(u, v)).collect(),
        closed: false,
    })
}

/// A cap's circle, centred on its plane's origin, from angle `start`,
/// anticlockwise or clockwise about the plane's normal.
fn cap_circle(start: Scalar, radius: Scalar, anticlockwise: bool) -> Curve2 {
    let x = Vec2::new(start.cos(), start.sin());
    let along = Vec2::new(-start.sin(), start.cos());
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x,
            y: if anticlockwise { along } else { -along },
        },
        radius,
    })
}

/// A closed cylinder, the way an exchange format with pcurves writes it: a
/// wall with one seam and two planar caps, every edge's every use with its
/// own pcurve.
fn closed_cylinder(rim: Rim) -> GeomResult<(TriMesh, MeshClosure, DeviationReport)> {
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let Rim {
        start,
        middle,
        bump,
        cap_growth,
        caps_first,
        spline,
    } = rim;
    let (u0, u1) = (start, start + TAU);
    let at_height = |z: Scalar| Point3::new(RADIUS * start.cos(), RADIUS * start.sin(), z);
    let caps = |builder: &mut GeometryGraphBuilder| {
        [
            push_curve(builder, cap_circle(start, RADIUS + cap_growth, false)),
            push_curve(builder, cap_circle(start, RADIUS + cap_growth, true)),
        ]
    };
    let early = caps_first.then(|| caps(&mut builder));
    // The wall's pcurves run in the direction of their use.
    let bottom_points = [(u0, 0.0), (u0 + middle * TAU, bump), (u1, 0.0)];
    let bottom = if spline {
        push_curve(
            &mut builder,
            Curve2::BSpline(BSplineCurve {
                degree: 2,
                control_points: bottom_points.map(|(u, v)| Point2::new(u, v)).to_vec(),
                knots: vec![0.0, 1.0],
                multiplicities: vec![3, 3],
                weights: None,
                closed: false,
                self_intersect: None,
                knot_spec: KnotSpec::Unspecified,
            }),
        )
    } else {
        push_curve(&mut builder, polyline(&bottom_points))
    };
    let top = push_curve(
        &mut builder,
        polyline(&[(u1, HEIGHT), (u1 - middle * TAU, HEIGHT), (u0, HEIGHT)]),
    );
    let seam_up = push_curve(&mut builder, polyline(&[(u1, 0.0), (u1, HEIGHT)]));
    let seam_down = push_curve(&mut builder, polyline(&[(u0, HEIGHT), (u0, 0.0)]));
    let [floor_circle, lid_circle] = early.unwrap_or_else(|| caps(&mut builder));
    let wall = push_surface(
        &mut builder,
        Surface::Cylinder(Cylinder {
            frame: frame_at(0.0),
            radius: RADIUS,
        }),
    );
    let floor = push_surface(
        &mut builder,
        Surface::Plane(Plane {
            frame: frame_at(0.0),
        }),
    );
    let lid = push_surface(
        &mut builder,
        Surface::Plane(Plane {
            frame: frame_at(HEIGHT),
        }),
    );
    let low = brep.add_vertex(Vertex {
        position: at_height(0.0),
    });
    let high = brep.add_vertex(Vertex {
        position: at_height(HEIGHT),
    });
    let rim_circle = |z: Scalar| {
        Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: Point3::new(0.0, 0.0, z),
                x: Vec3::new(start.cos(), start.sin(), 0.0),
                y: Vec3::new(-start.sin(), start.cos(), 0.0),
                z: Vec3::Z,
            },
            radius: RADIUS,
        })
    };
    let rim_curves = [rim_circle(0.0), rim_circle(HEIGHT)];
    let mut edge = |start, end, curve: Option<Curve3>| {
        let curve = curve.map(|c| builder.push(GeometryNode::Curve3(c)).unwrap());
        brep.add_edge(Edge { start, end, curve })
    };
    let [bottom_curve, top_curve] = rim_curves;
    let (rim0, rim1, seam) = (
        edge(low, low, Some(bottom_curve)),
        edge(high, high, Some(top_curve)),
        edge(low, high, None),
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
        face(wall, wall_loop, Orientation::Forward),
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
    let root = builder.push(GeometryNode::BRep(brep)).unwrap();
    let graph = builder.finish(vec![root]).unwrap();
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options())?;
    Ok((outcome.mesh, outcome.closure, report))
}

fn push_curve(builder: &mut GeometryGraphBuilder, curve: Curve2) -> NodeId {
    builder.push(GeometryNode::Curve2(curve)).unwrap()
}

fn push_surface(builder: &mut GeometryGraphBuilder, surface: Surface) -> NodeId {
    builder.push(GeometryNode::Surface(surface)).unwrap()
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

fn assert_closed_cylinder(name: &str, rim: Rim) {
    let (mesh, closure, report) = closed_cylinder(rim).unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_eq!(closure, MeshClosure::Solid, "{name}");
    assert_watertight(name, &mesh);
    let bound = report
        .bound
        .unwrap_or_else(|| panic!("{name}: expected a certified bound, got {report:?}"));
    // A signed sum, not `enclosed_volume`: a rim whose boundary vertices are
    // collinear in the wall's parameters leaves ears of exactly no area, as
    // before, which that audit refuses at zero tolerance.
    let volume: Scalar = mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            a.dot(b.cross(c)) / 6.0
        })
        .sum();
    let exact = 0.5 * TAU * RADIUS * RADIUS * HEIGHT;
    let area = TAU * RADIUS * (HEIGHT + RADIUS);
    let slack = bound * area;
    assert!(
        (volume - exact).abs() <= slack,
        "{name}: mesh volume {volume} against the exact {exact}, slack {slack}"
    );
    let mut samples = Vec::new();
    for k in 0..720 {
        let angle = TAU * k as Scalar / 720.0;
        for z in [0.0, HEIGHT] {
            samples.push(Point3::new(RADIUS * angle.cos(), RADIUS * angle.sin(), z));
            for r in [0.25, 0.5, 0.75] {
                samples.push(Point3::new(
                    r * RADIUS * angle.cos(),
                    r * RADIUS * angle.sin(),
                    z,
                ));
            }
        }
        for z in [0.1, 0.25, 0.4] {
            samples.push(Point3::new(RADIUS * angle.cos(), RADIUS * angle.sin(), z));
        }
    }
    let worst = farthest(&mesh, &samples);
    assert!(
        worst <= bound,
        "{name}: a surface point lies {worst:e} from the mesh, above the reported {bound:e}"
    );
}

#[test]
fn caps_and_wall_walking_the_rim_at_one_pace_compile_certified() {
    for caps_first in [false, true] {
        let name = format!("even pace, caps first {caps_first}");
        assert_closed_cylinder(
            &name,
            Rim {
                caps_first,
                ..Rim::EVEN
            },
        );
    }
}

#[test]
fn a_cap_circle_starting_at_another_angle_than_the_origin_compiles_certified() {
    for (start, caps_first) in [(0.7, false), (2.5, true), (-1.3, false), (-1.3, true)] {
        let name = format!("circles starting at {start}, caps first {caps_first}");
        assert_closed_cylinder(
            &name,
            Rim {
                start,
                caps_first,
                ..Rim::EVEN
            },
        );
    }
}

#[test]
fn a_wall_rim_walked_at_another_pace_than_its_caps_compiles_certified() {
    for caps_first in [false, true] {
        for (start, middle) in [(0.0, 0.35), (0.7, 0.62), (2.5, 0.42)] {
            let name = format!(
                "rim polyline bent at {middle}, starting at {start}, caps first {caps_first}"
            );
            assert_closed_cylinder(
                &name,
                Rim {
                    start,
                    middle,
                    caps_first,
                    ..Rim::EVEN
                },
            );
        }
    }
}

#[test]
fn the_bound_over_a_rim_at_another_pace_reads_the_located_spans() {
    // The wall's rim is sampled at the caps' uniform angles, which the
    // polyline reaches at parameters of its own. Its chord bounds read over
    // the uniform parameters would put the corner in the wrong span and give
    // 6.1e-4 here; over the located ones the bound is 8.6e-4.
    let (_, _, report) = closed_cylinder(Rim {
        middle: 0.35,
        ..Rim::EVEN
    })
    .unwrap();
    let bound = report.bound.expect("a certified bound");
    assert!(bound >= 8e-4, "{bound:e}");
}

#[test]
fn a_rim_bent_at_the_caps_pace_with_a_bump_compiles_certified_whoever_samples_it() {
    // Both rims are corners at the middle of the turn, a sample of the caps
    // whichever use samples the edge: the located spans across a corner are
    // bounded piecewise, and the mesh stays closed.
    for caps_first in [false, true] {
        for bump in [1e-4, 3e-4, 6e-4, 9e-4] {
            let name = format!("bump {bump:e}, caps first {caps_first}");
            assert_closed_cylinder(
                &name,
                Rim {
                    bump,
                    caps_first,
                    ..Rim::EVEN
                },
            );
        }
    }
}

#[test]
fn pcurves_a_hundredth_apart_are_still_refused() {
    // The wall's bottom rim bulges by 1e-2 at its middle; both ends stay
    // on the cap's circle.
    let error = closed_cylinder(Rim {
        bump: 1e-2,
        ..Rim::EVEN
    })
    .map(|_| ())
    .expect_err("a rim a hundredth off its cap is refused");
    assert!(
        matches!(&error, GeomError::Degenerate(text)
            if text == "a shared edge sample is not located on a neighbouring pcurve"),
        "{error:?}"
    );
}

#[test]
fn a_cap_a_hundredth_wider_than_the_wall_is_refused() {
    let error = closed_cylinder(Rim {
        cap_growth: 1e-2,
        ..Rim::EVEN
    })
    .map(|_| ())
    .expect_err("a cap a hundredth off the wall is refused");
    assert!(matches!(&error, GeomError::Degenerate(_)), "{error:?}");
}

#[test]
fn a_wall_rim_on_an_uneven_spline_is_located_and_refined_when_the_caps_sample_it() {
    // The rim's spline walks the turn at nine times the pace at one end as
    // at the other. The caps' circles are the sampling use, and the located
    // spans of the spline outgrow the budget at the samples their circle
    // needs: the edge is sampled again at twice the count.
    for caps_first in [false, true] {
        let rim = Rim {
            middle: 0.05,
            spline: true,
            caps_first,
            ..Rim::EVEN
        };
        let (mesh, closure, report) = closed_cylinder(rim).unwrap();
        let name = format!("spline rim, caps first {caps_first}");
        assert_eq!(closure, MeshClosure::Solid, "{name}");
        assert_watertight(&name, &mesh);
        let bound = report.bound.expect("a certified bound");
        assert!(bound <= BUDGET, "{name}: bound {bound:e} misses the budget");
    }
}

#[test]
fn a_bump_within_the_budget_is_admitted_and_reported() {
    let at = |bump| {
        closed_cylinder(Rim { bump, ..Rim::EVEN })
            .unwrap()
            .2
            .bound
            .expect("a certified bound")
    };
    let (flat, bumped) = (at(0.0), at(5e-4));
    // The 5e-4 offset of the samples is in the triangles' bound once, and
    // not again in the trim terms.
    assert!(
        bumped >= flat + 4.5e-4,
        "a rim 5e-4 off its cap widens the bound from {flat:e} only to {bumped:e}"
    );
    assert!(
        bumped <= flat + 5.5e-4,
        "a rim 5e-4 off its cap widens the bound from {flat:e} to {bumped:e}"
    );
}
