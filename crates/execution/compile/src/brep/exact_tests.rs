//! `tessellate_exact` against the exact B-rep it reads (#391).
//!
//! Each case compares the mesh with its source: the volume from the mesh
//! and from the exact integral, the boundary Hausdorff distance against the
//! certified bound, a dense sample of the surface against the mesh, and the
//! closure.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_6, PI};

use axiolid_brep::ExactBRep;
use axiolid_brep::ExactBRepBuilder;
use axiolid_brep_boolean::ToleranceDecisionKind;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{
    BooleanOperator, Frame3, Interval, Point2, Point3, Tolerance, Transform3, Vec2, Vec3,
};
use axiolid_curve::{Circle3, Curve2, Curve3, Line2, Line3, Polyline2};
use axiolid_measure::{
    exact_properties, one_sided_boundary_hausdorff_with_budget, volume_properties,
};
use axiolid_mesh::TriMesh;
use axiolid_mesh_compile_contract::MeshClosure;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, Instance, SolidOperation};
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};
use axiolid_surface::{Cylinder, Surface};
use axiolid_topology::{Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Vertex};

use super::view::from_exact;
use super::{tessellate_exact, Deviation};
use crate::deviation::{DeviationBound, DeviationPath};
use crate::exact::polyhedral::from_triangles;
use crate::ReferenceExactCompiler;

const TOLERANCE: Tolerance = Tolerance::METRE;

/// The refinement budget of each one-sided Hausdorff measurement.
const SPLITS: usize = 4_000;

fn rectangle(x: f64, y: f64) -> Profile {
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

/// The mesh of `exact` at `chord`, which must be a closed solid.
fn mesh_of(exact: &ExactBRep, chord: f64) -> (TriMesh, Deviation) {
    let (mesh, closure, deviation) =
        tessellate_exact(exact, TOLERANCE, chord).expect("a tessellable exact B-rep");
    assert_eq!(closure, MeshClosure::Solid);
    (mesh, deviation)
}

fn mesh_volume(mesh: &TriMesh) -> f64 {
    volume_properties(mesh, TOLERANCE)
        .expect("a closed two-manifold")
        .signed_volume
}

fn exact_volume(exact: &ExactBRep) -> f64 {
    exact_properties(exact, TOLERANCE)
        .expect("measurable")
        .signed_volume
}

/// Both one-sided boundary distances' lower bounds, between the exact
/// B-rep and the polyhedron of its mesh, are within the certified bound.
fn assert_within_bound(exact: &ExactBRep, mesh: &TriMesh, deviation: &Deviation) {
    let bound = deviation.bound().expect("a certified bound");
    let polyhedron = from_triangles(mesh).expect("a closed mesh");
    // To the bound itself and within a fixed budget: a lower bound is a
    // witness however coarse the interval above it is.
    let measure = |from, to| {
        one_sided_boundary_hausdorff_with_budget(from, to, bound, TOLERANCE, SPLITS)
            .expect("measurable")
    };
    let sides = [
        ("forward", measure(exact, &polyhedron)),
        ("backward", measure(&polyhedron, exact)),
    ];
    for (name, side) in sides {
        assert!(
            side.lower <= bound,
            "{name}: lower bound {} above the certified {bound}",
            side.lower
        );
    }
}

#[test]
fn a_rectangle_prism_is_its_own_polyhedron() {
    let exact = extrude_profile_exact(&rectangle(2.0, 3.0), Vec3::Z, 4.0, TOLERANCE).unwrap();
    let (mesh, deviation) = mesh_of(&exact, 1e-2);
    assert_eq!(deviation.bound(), Some(0.0));
    assert!((mesh_volume(&mesh) - 24.0).abs() <= 1e-12);
    assert!((exact_volume(&exact) - 24.0).abs() <= 1e-12);
    assert_within_bound(&exact, &mesh, &deviation);
}

#[test]
fn a_cylinder_is_within_its_certified_bound() {
    let radius = 1.5;
    let exact = extrude_profile_exact(&circle(radius), Vec3::Z, 2.0, TOLERANCE).unwrap();
    let (mesh, deviation) = mesh_of(&exact, 1e-2);
    let bound = deviation.bound().expect("a certified bound");
    assert!(bound > 0.0 && bound <= 1e-2, "{bound}");
    // A chord-flattened wall removes at most bound x the wall's area, and
    // the caps are planar.
    let volume = exact_volume(&exact);
    assert!((volume - PI * radius * radius * 2.0).abs() <= 1e-9);
    let wall = 2.0 * PI * radius * 2.0;
    assert!((mesh_volume(&mesh) - volume).abs() <= bound * wall);
    assert_within_bound(&exact, &mesh, &deviation);
}

#[test]
fn a_quarter_torus_is_within_its_certified_bound() {
    let (major, minor) = (2.0, 0.25);
    let profile = circle(minor);
    // The profile circle is centred on the origin: move the axis to the
    // major radius.
    let exact = revolve_profile_exact(
        &profile,
        Point3::new(-major, 0.0, 0.0),
        Vec3::Y,
        FRAC_PI_2,
        TOLERANCE,
    )
    .unwrap();
    let (mesh, deviation) = mesh_of(&exact, 1e-2);
    let bound = deviation.bound().expect("a certified bound");
    assert!(bound <= 2e-2, "{bound}");
    let volume = exact_volume(&exact);
    let expected = 0.25 * 2.0 * PI * PI * major * minor * minor;
    assert!((volume - expected).abs() <= 1e-9 * expected, "{volume}");
    let area = exact_properties(&exact, TOLERANCE).unwrap().area;
    assert!((mesh_volume(&mesh) - volume).abs() <= bound * area);
    assert_within_bound(&exact, &mesh, &deviation);
}

/// A box and a dowel through its top face, in a graph compiled exactly.
fn box_and_dowel() -> ExactBRep {
    let mut builder = GeometryGraphBuilder::new();
    let mut push = |node| builder.push(node).expect("a valid node");
    let slab = push(GeometryNode::Profile(rectangle(4.0, 4.0)));
    let slab = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile: slab,
        direction: Vec3::Z,
        depth: 1.0,
    }));
    let pin = push(GeometryNode::Profile(circle(0.5)));
    let pin = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile: pin,
        direction: Vec3::Z,
        depth: 2.0,
    }));
    // Standing on the slab's top face, off-centre: one shell through the
    // contact.
    let pin = push(GeometryNode::Instance(Instance {
        source: pin,
        transform: Transform3::from_translation(Vec3::new(0.5, 0.25, 1.0)),
    }));
    let root = push(GeometryNode::SolidOperation(SolidOperation::Boolean {
        left: slab,
        right: pin,
        operator: BooleanOperator::Union,
    }));
    let graph = builder.finish(vec![root]).expect("a valid graph");
    let options = ExecutionOptions::new(TOLERANCE);
    let mut solids = ReferenceExactCompiler::new()
        .compile_exact_batch_with_reports(&graph, &[root], &options)
        .expect("exact union");
    let (union, report) = solids.remove(0);
    assert!(!report.contains(ToleranceDecisionKind::CoincidentSupports));
    union
}

#[test]
fn an_exact_union_of_a_box_and_a_dowel_is_one_solid_within_its_bound() {
    let exact = box_and_dowel();
    let (mesh, deviation) = mesh_of(&exact, 1e-2);
    let bound = deviation.bound().expect("a certified bound");
    let volume = exact_volume(&exact);
    let expected = 16.0 + PI * 0.25 * 2.0;
    assert!((volume - expected).abs() <= 1e-9, "{volume}");
    let area = exact_properties(&exact, TOLERANCE).unwrap().area;
    assert!((mesh_volume(&mesh) - volume).abs() <= bound * area);
    assert_within_bound(&exact, &mesh, &deviation);
}

/// The cylinder `x^2 + y^2 = 1` between `u = 0` and `u = pi / 2` and
/// `z = 0` to `1`, as one face of four edges walked anticlockwise in
/// `(u, v)`. The two edges the face walks against their own direction run
/// their pcurves over a reversed span. The bottom pcurve is a line, or a
/// polyline through its middle (a family with no chord bound) when
/// `kinked`.
fn cylinder_patch(kinked: bool) -> ExactBRep {
    let mut b = ExactBRepBuilder::default();
    let world = Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    };
    let at = |x, y, z| Vertex {
        position: Point3::new(x, y, z),
    };
    let topology = b.topology_mut();
    let v = [
        topology.add_vertex(at(1.0, 0.0, 0.0)),
        topology.add_vertex(at(0.0, 1.0, 0.0)),
        topology.add_vertex(at(0.0, 1.0, 1.0)),
        topology.add_vertex(at(1.0, 0.0, 1.0)),
    ];
    let arc = |height| {
        Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: Point3::new(0.0, 0.0, height),
                ..world
            },
            radius: 1.0,
        })
    };
    let rule = |x, y| {
        Curve3::Line(Line3 {
            origin: Point3::new(x, y, 0.0),
            direction: Vec3::Z,
        })
    };
    let curves = [arc(0.0), rule(0.0, 1.0), arc(1.0), rule(1.0, 0.0)];
    // Bottom, right, top (from `v3` to `v2`), left (from `v0` to `v3`).
    let ends = [(v[0], v[1]), (v[1], v[2]), (v[3], v[2]), (v[0], v[3])];
    let spans = [
        Interval::new(0.0, FRAC_PI_2),
        Interval::new(0.0, 1.0),
        Interval::new(0.0, FRAC_PI_2),
        Interval::new(0.0, 1.0),
    ];
    let mut edges = Vec::new();
    for ((curve, (start, end)), span) in curves.into_iter().zip(ends).zip(spans) {
        let curve = b.add_curve3(curve);
        let edge = b.topology_mut().add_edge(Edge {
            start,
            end,
            curve: Some(curve),
        });
        b.set_edge_interval(edge, span);
        edges.push(edge);
    }
    let line = |x, y, dx, dy| {
        Curve2::Line(Line2 {
            origin: Vec2::new(x, y),
            direction: Vec2::new(dx, dy),
        })
    };
    let (bottom, bottom_span) = if kinked {
        (
            Curve2::Polyline(Polyline2 {
                points: vec![
                    Point2::new(0.0, 0.0),
                    Point2::new(FRAC_PI_6, 0.0),
                    Point2::new(FRAC_PI_3, 0.0),
                    Point2::new(FRAC_PI_2, 0.0),
                ],
                closed: false,
            }),
            Interval::new(0.0, 3.0),
        )
    } else {
        (line(0.0, 0.0, 1.0, 0.0), Interval::new(0.0, FRAC_PI_2))
    };
    // (pcurve, span in loop order, orientation of the use)
    let uses = [
        (bottom, bottom_span, Orientation::Forward),
        (
            line(FRAC_PI_2, 0.0, 0.0, 1.0),
            Interval::new(0.0, 1.0),
            Orientation::Forward,
        ),
        (
            line(0.0, 1.0, 1.0, 0.0),
            Interval::new(FRAC_PI_2, 0.0),
            Orientation::Reversed,
        ),
        (
            line(0.0, 0.0, 0.0, 1.0),
            Interval::new(1.0, 0.0),
            Orientation::Reversed,
        ),
    ];
    let mut loop_uses = Vec::new();
    let mut intervals = Vec::new();
    for (position, ((pcurve, span, orientation), edge)) in uses.into_iter().zip(&edges).enumerate()
    {
        let pcurve = b.add_curve2(pcurve);
        loop_uses.push(EdgeUse {
            edge: *edge,
            orientation,
            pcurve: Some(pcurve),
        });
        intervals.push((position, span));
    }
    let wire = b.topology_mut().add_loop(Loop { edges: loop_uses });
    for (position, span) in intervals {
        b.set_pcurve_interval(wire, position, span);
    }
    let surface = b.add_surface(Surface::Cylinder(Cylinder {
        frame: world,
        radius: 1.0,
    }));
    b.topology_mut().add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id: wire,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    b.finish().expect("a valid patch")
}

#[test]
fn a_patch_whose_uses_run_against_their_edges_is_meshed_within_its_bound() {
    let patch = cylinder_patch(false);
    let (mesh, closure, deviation) = tessellate_exact(&patch, TOLERANCE, 1e-2).unwrap();
    assert_eq!(closure, MeshClosure::Surface);
    let bound = deviation.bound().expect("a certified bound");
    // The spans of the two walked-against uses run backward.
    let view = from_exact(&patch).unwrap();
    let spans: Vec<Interval> = view.loops()[0]
        .edges
        .iter()
        .map(|use_| use_.pcurve.unwrap().span)
        .collect();
    assert_eq!(spans[2], Interval::new(FRAC_PI_2, 0.0));
    assert_eq!(spans[3], Interval::new(1.0, 0.0));
    // Every sampled point of the exact surface is within the bound of the
    // mesh.
    let mut worst = 0.0_f64;
    for i in 0..=60 {
        for j in 0..=10 {
            let u = FRAC_PI_2 * f64::from(i) / 60.0;
            let point = Point3::new(u.cos(), u.sin(), f64::from(j) / 10.0);
            let probe = TriMesh::new(vec![point], vec![0, 0, 0]);
            let gap = axiolid_inspect::surface_gap(&probe, &mesh, 1.0).expect("near the mesh");
            worst = worst.max(gap);
        }
    }
    assert!(worst <= bound, "{worst} above the certified {bound}");
    // The mesh covers the patch: its area is the quarter wall's.
    let area: f64 = mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [t[0], t[1], t[2]].map(|i| mesh.positions[i as usize]);
            0.5 * (b - a).cross(c - a).length()
        })
        .sum();
    assert!((area - FRAC_PI_2).abs() <= bound * 4.0, "{area}");
}

#[test]
fn a_pcurve_family_with_no_chord_bound_is_unbounded_by_name() {
    let (_, _, deviation) = tessellate_exact(&cylinder_patch(true), TOLERANCE, 1e-2).unwrap();
    assert_eq!(deviation.bound(), None);
    let report = deviation.report(1e-2);
    let named: Vec<_> = report.unbounded().map(|c| (c.path, c.bound)).collect();
    assert_eq!(
        named,
        [(
            DeviationPath::BRepFace,
            DeviationBound::Unbounded("trim pcurve family with no chord bound")
        )]
    );
}
