//! Exact compilation under rigid placements, and exact swept disks (#223).
//!
//! A consumer builds a body's exact B-rep from the same graph its mesh
//! comes from: an extrusion, revolution or swept disk under one or more
//! `Instance` placements. Each case here compiles such a graph and checks
//! the certified distance to a slab whose top is the plane `z = 0` against
//! a closed form (see `axiolid-construct`'s `placement` tests for the
//! geometry), and the exact volume against its formula.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use axiolid_brep::ExactBRep;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Frame3, Point3, Tolerance, Transform3, Vec3};
use axiolid_curve::{Circle3, Curve3, Line3, Polyline3};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::{boundary_distance, exact_properties};
use axiolid_mesh_compile::ReferenceExactCompiler;
use axiolid_model::{
    CurveRelation, GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId,
    SolidOperation, TrimSelector, TrimmingPreference,
};
use axiolid_profile::{CircleProfile, EllipseProfile, Profile, RectangleProfile};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

struct Graph {
    builder: GeometryGraphBuilder,
}

impl Graph {
    fn new() -> Self {
        Self {
            builder: GeometryGraphBuilder::new(),
        }
    }

    fn push(&mut self, node: GeometryNode) -> NodeId {
        self.builder.push(node).expect("a valid node")
    }

    fn extrusion(&mut self, profile: Profile, depth: f64) -> NodeId {
        let profile = self.push(GeometryNode::Profile(profile));
        self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
    }

    fn instance(&mut self, source: NodeId, transform: Transform3) -> NodeId {
        self.push(GeometryNode::Instance(Instance { source, transform }))
    }

    fn swept_disk(
        &mut self,
        directrix: NodeId,
        radius: f64,
        inner_radius: Option<f64>,
        parameter_range: Option<(f64, f64)>,
    ) -> NodeId {
        self.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
            directrix,
            radius,
            inner_radius,
            parameter_range,
            fillet_radius: None,
        }))
    }

    /// The slab: a 40 x 40 x 1 block moved down so its top is `z = 0`.
    fn slab(&mut self) -> NodeId {
        let block = self.extrusion(
            Profile::Rectangle(RectangleProfile {
                x: 40.0,
                y: 40.0,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            }),
            1.0,
        );
        self.instance(block, Transform3::from_translation(-Vec3::Z))
    }

    fn finish(self, roots: Vec<NodeId>) -> GeometryGraph {
        self.builder.finish(roots).expect("a valid graph")
    }
}

fn compile(graph: &GeometryGraph, roots: &[NodeId]) -> Result<Vec<ExactBRep>, GeomError> {
    let mut out = Vec::new();
    ReferenceExactCompiler::new()
        .compile_exact_batch_into(graph, roots, &options(), &mut out)
        .map(|()| out)
}

/// Compile `body` and the slab, check the volume and the certified
/// distance between them.
fn check(graph: Graph, body: NodeId, slab: NodeId, volume: f64, distance: f64) {
    let graph = graph.finish(vec![body, slab]);
    let solids = compile(&graph, &[body, slab]).expect("exact");
    let measured = exact_properties(&solids[0], Tolerance::METRE)
        .expect("measurable")
        .signed_volume;
    assert!(
        (measured - volume).abs() <= 1e-9 * volume,
        "volume {measured}, expected {volume}"
    );
    let bounds =
        boundary_distance(&solids[0], &solids[1], 1e-9, Tolerance::METRE).expect("bounded");
    assert!(
        bounds.lower <= distance + 1e-12 && distance <= bounds.upper + 1e-12,
        "[{}, {}] must contain {distance}",
        bounds.lower,
        bounds.upper
    );
    assert!(bounds.upper - bounds.lower <= 1e-9);
}

/// Tilt by `tilt` about x, turn by `turn` about z, then shift to `at`.
fn placement(tilt: f64, turn: f64, at: Vec3) -> Transform3 {
    Transform3::from_translation(at)
        * Transform3::from_rotation_z(turn)
        * Transform3::from_rotation_x(tilt)
}

#[test]
fn an_instanced_circular_extrusion_compiles_exactly_under_a_tilt() {
    let mut graph = Graph::new();
    let (r, h, tilt) = (0.3, 3.0, 0.4);
    let column = graph.extrusion(
        Profile::Circle(CircleProfile {
            radius: r,
            thickness: None,
        }),
        h,
    );
    // Two nested placements compose: the inner tilt, then the outer turn
    // and shift.
    let tilted = graph.instance(column, Transform3::from_rotation_x(tilt));
    let placed = graph.instance(
        tilted,
        Transform3::from_translation(Vec3::new(2.0, -1.0, 1.5)) * Transform3::from_rotation_z(0.7),
    );
    let slab = graph.slab();
    check(graph, placed, slab, PI * r * r * h, 1.5 - r * tilt.sin());
}

#[test]
fn an_instanced_elliptical_extrusion_compiles_exactly_under_a_tilt() {
    let mut graph = Graph::new();
    let (a, b, h, tilt) = (0.5, 0.2, 2.0, 0.6);
    let column = graph.extrusion(
        Profile::Ellipse(EllipseProfile {
            semi_axis_x: a,
            semi_axis_y: b,
        }),
        h,
    );
    let placed = graph.instance(column, placement(tilt, -1.1, Vec3::new(-3.0, 0.5, 1.2)));
    let slab = graph.slab();
    check(graph, placed, slab, PI * a * b * h, 1.2 - b * tilt.sin());
}

#[test]
fn an_instanced_full_turn_revolution_compiles_exactly_under_a_tilt() {
    let mut graph = Graph::new();
    let (major, minor, tilt) = (1.0, 0.2, 0.5);
    let profile = graph.push(GeometryNode::Profile(Profile::Circle(CircleProfile {
        radius: minor,
        thickness: None,
    })));
    // Centred on the origin, axis along y.
    let torus = graph.push(GeometryNode::SolidOperation(SolidOperation::Revolution {
        profile,
        axis_origin: Point3::new(-major, 0.0, 0.0),
        axis_direction: Vec3::Y,
        angle: TAU,
    }));
    let centred = graph.instance(torus, Transform3::from_translation(Vec3::X * major));
    let placed = graph.instance(centred, placement(tilt, 0.9, Vec3::new(0.5, -0.5, 2.0)));
    let slab = graph.slab();
    check(
        graph,
        placed,
        slab,
        2.0 * PI * PI * major * minor * minor,
        2.0 - major * tilt.cos() - minor,
    );
}

#[test]
fn a_disk_swept_along_a_line_compiles_exactly() {
    let mut graph = Graph::new();
    let line = graph.push(GeometryNode::Curve3(Curve3::Line(Line3 {
        // Parameter 0.5 is (0, 0, 3), parameter 1.5 is (1, 2, 1.5).
        origin: Point3::new(-0.5, -1.0, 3.75),
        direction: Vec3::new(1.0, 2.0, -1.5),
    })));
    let r = 0.2;
    let pipe = graph.swept_disk(line, r, None, Some((0.5, 1.5)));
    let slab = graph.slab();
    let d = Vec3::new(1.0, 2.0, -1.5);
    let length = d.length();
    let dz = d.z / length;
    check(
        graph,
        pipe,
        slab,
        PI * r * r * length,
        1.5 - r * (1.0 - dz * dz).sqrt(),
    );
}

#[test]
fn a_disk_swept_along_a_two_point_polyline_under_a_placement_compiles_exactly() {
    let mut graph = Graph::new();
    let polyline = graph.push(GeometryNode::Curve3(Curve3::Polyline(Polyline3 {
        points: vec![Point3::ZERO, Point3::new(0.0, 0.0, 2.0)],
        closed: false,
    })));
    let (r, ri, tilt) = (0.3, 0.2, 0.4);
    let pipe = graph.swept_disk(polyline, r, Some(ri), None);
    let placed = graph.instance(pipe, placement(tilt, 0.3, Vec3::new(0.0, 0.0, 1.5)));
    let slab = graph.slab();
    check(
        graph,
        placed,
        slab,
        PI * (r * r - ri * ri) * 2.0,
        1.5 - r * tilt.sin(),
    );
}

fn vertical_arc(radius: f64) -> Circle3 {
    // In the plane y = 0.5, centre at height 3; t = 3 pi / 2 is lowest.
    Circle3 {
        frame: Frame3 {
            origin: Point3::new(0.2, 0.5, 3.0),
            x: Vec3::X,
            y: Vec3::Z,
            z: -Vec3::Y,
        },
        radius,
    }
}

#[test]
fn a_disk_swept_along_an_arc_compiles_exactly() {
    let mut graph = Graph::new();
    let circle = graph.push(GeometryNode::Curve3(Curve3::Circle(vertical_arc(2.0))));
    let r = 0.25;
    let (start, end) = (PI, 2.0 * PI - 0.2);
    let pipe = graph.swept_disk(circle, r, None, Some((start, end)));
    let slab = graph.slab();
    check(
        graph,
        pipe,
        slab,
        PI * r * r * 2.0 * (end - start),
        3.0 - 2.0 - r,
    );
}

#[test]
fn a_disk_swept_along_a_trimmed_arc_across_the_seam_compiles_exactly() {
    // A trim from 3 pi / 2 + 0.5 forward to 0.5 crosses the seam: the
    // swept angle is pi / 2, not its complement, and the arc's lowest
    // point is just before it, so the nearest point is the start cap's
    // rim.
    let mut graph = Graph::new();
    let basis = graph.push(GeometryNode::Curve3(Curve3::Circle(vertical_arc(2.0))));
    let (a, b) = (1.5 * PI + 0.5, 0.5);
    let trimmed = graph.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
        basis,
        start: vec![TrimSelector::Parameter(a)],
        end: vec![TrimSelector::Parameter(b)],
        sense_agreement: true,
        preference: TrimmingPreference::Parameter,
    }));
    let r = 0.25;
    let pipe = graph.swept_disk(trimmed, r, None, None);
    let slab = graph.slab();
    // The start cap is centred at angle a on the circle, at height
    // 3 + 2 sin(a); its disk is perpendicular to the tangent, whose
    // vertical part is cos(a), so it dips r sqrt(1 - cos^2 a) below that.
    let centre_z = 3.0 + 2.0 * a.sin();
    let tangent_z = a.cos();
    let dip = r * (1.0 - tangent_z * tangent_z).sqrt();
    check(
        graph,
        pipe,
        slab,
        PI * r * r * 2.0 * FRAC_PI_2,
        centre_z - dip,
    );
}

#[test]
fn a_scaled_instance_is_refused_by_name() {
    let mut graph = Graph::new();
    let column = graph.extrusion(
        Profile::Circle(CircleProfile {
            radius: 0.3,
            thickness: None,
        }),
        1.0,
    );
    let scaled = graph.instance(column, Transform3::from_scale(Vec3::new(1.0, 1.0, 2.0)));
    let graph = graph.finish(vec![scaled]);
    let error = compile(&graph, &[scaled]).unwrap_err();
    assert!(
        matches!(error, GeomError::UnsupportedInput { input, .. } if input.contains("scaled")),
        "{error:?}"
    );
}

#[test]
fn a_swept_disk_the_exact_path_cannot_follow_is_refused_by_name() {
    let mut graph = Graph::new();
    let bent = graph.push(GeometryNode::Curve3(Curve3::Polyline(Polyline3 {
        points: vec![Point3::ZERO, Point3::X, Point3::new(1.0, 1.0, 0.0)],
        closed: false,
    })));
    let cornered = graph.swept_disk(bent, 0.1, None, None);
    let line = graph.push(GeometryNode::Curve3(Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::X,
    })));
    let unbounded = graph.swept_disk(line, 0.1, None, None);
    let circle = graph.push(GeometryNode::Curve3(Curve3::Circle(vertical_arc(0.1))));
    let fat = graph.swept_disk(circle, 0.2, None, None);
    let graph = graph.finish(vec![cornered, unbounded, fat]);
    for root in [cornered, unbounded, fat] {
        let error = compile(&graph, &[root]).unwrap_err();
        assert!(
            matches!(error, GeomError::UnsupportedInput { .. }),
            "{error:?}"
        );
    }
}
