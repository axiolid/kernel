//! Curved profiles that change between stations (#361): profiles that pair
//! segment for segment are flattened together, so the run is certified and
//! a dense sampling of the *exact* interpolated surface (each section the
//! equal-fraction blend of the two profiles' curves, written out in closed
//! form below, independent of the compiler's rings) never strays from the
//! mesh by more than the bound.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Vec2};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, DeviationReport, ReferenceMeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, SectionAtStation, SolidOperation,
    Station, StationFrame, StationOffsets,
};
use axiolid_profile::{CircleProfile, EllipseProfile, Profile, RectangleProfile};
use axiolid_reference::station::{station_section2, SectionFrame};

/// The chord budget of every run here, in metres.
const BUDGET: Scalar = 0.01;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(BUDGET)
        .unwrap()
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn disk(b: &mut GeometryGraphBuilder, radius: Scalar) -> NodeId {
    b.push_value(Profile::Circle(CircleProfile {
        radius,
        thickness: None,
    }))
    .unwrap()
}

fn rounded(b: &mut GeometryGraphBuilder, radius: Scalar) -> NodeId {
    b.push_value(Profile::Rectangle(RectangleProfile {
        x: 2.0,
        y: 1.4,
        thickness: None,
        outer_radius: Some(radius),
        inner_radius: None,
    }))
    .unwrap()
}

fn section(profile: NodeId, distance: Scalar) -> SectionAtStation {
    SectionAtStation::new(profile, Station::new(distance, StationOffsets::ZERO))
}

fn spine(
    b: &mut GeometryGraphBuilder,
    directrix: NodeId,
    sections: Vec<SectionAtStation>,
) -> NodeId {
    b.push(GeometryNode::SolidOperation(
        SolidOperation::SectionsAtStations {
            directrix,
            sections,
            frame: StationFrame::Section,
        },
    ))
    .unwrap()
}

fn line() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::new(2.0, 1.0),
        direction: Vec2::new(0.6, 0.8),
    })
}

fn arc(radius: Scalar) -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius,
    })
}

/// Points of the exact surface: the exact section at `steps` + 1 fractions
/// of the run, placed on the frame law, `exact(u)` its boundary samples.
fn exact_surface(
    frame_at: &dyn Fn(Scalar) -> SectionFrame,
    (s0, s1): (Scalar, Scalar),
    exact: &dyn Fn(Scalar) -> Vec<Point2>,
    steps: usize,
) -> Vec<Point3> {
    let mut out = Vec::new();
    for a in 0..=steps {
        let u = a as Scalar / steps as Scalar;
        let frame = frame_at(s0 + (s1 - s0) * u);
        let origin = frame.place(0.0, 0.0, 0.0);
        out.extend(
            exact(u)
                .into_iter()
                .map(|p| origin + p.x * frame.lateral + p.y * frame.up),
        );
    }
    out
}

fn point_triangle(p: Point3, [a, b, c]: [Point3; 3]) -> Scalar {
    axiolid_measure::proximity::closest_point_on_triangle(p, [a, b, c])
        .map(|q| p.distance(q))
        .unwrap_or(Scalar::INFINITY)
}

/// The largest distance from `points` to the mesh's triangles.
fn measured(mesh: &TriMesh, points: &[Point3]) -> Scalar {
    let triangles: Vec<[Point3; 3]> = mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| [0, 1, 2].map(|k| mesh.positions[t[k] as usize]))
        .collect();
    points
        .iter()
        .map(|&p| {
            triangles
                .iter()
                .map(|t| point_triangle(p, *t))
                .fold(Scalar::INFINITY, Scalar::min)
        })
        .fold(0.0, Scalar::max)
}

fn certified(report: &DeviationReport) -> Scalar {
    let contribution = report
        .contributions
        .iter()
        .find(|c| c.path == DeviationPath::StationedSpine)
        .expect("the run reports its path");
    match contribution.bound {
        DeviationBound::Certified(value) => value,
        other => panic!("not certified: {other:?}"),
    }
}

/// Compile `root`, check it is certified within the budget, and that the
/// exact surface strays from the mesh by no more than the bound.
fn assert_holds(graph: &GeometryGraph, root: NodeId, points: &[Point3]) -> Scalar {
    let (outcome, report) = compiler()
        .compile_mesh_with_deviation(graph, root, &options())
        .unwrap();
    assert!(report.meets_requested(), "{report:?}");
    let bound = certified(&report);
    let far = measured(&outcome.mesh, points);
    assert!(
        far <= bound + 1e-9,
        "the exact surface strays {far} from the mesh, past the bound {bound}"
    );
    assert!(bound > 0.0, "a curved run has a bound");
    bound
}

/// A circle of `radius`, `count` points.
fn circle_points(radius: Scalar, count: usize) -> Vec<Point2> {
    (0..count)
        .map(|k| {
            let t = std::f64::consts::TAU * k as Scalar / count as Scalar;
            Point2::new(radius * t.cos(), radius * t.sin())
        })
        .collect()
}

fn disk_taper(basis: Curve2, run: (Scalar, Scalar)) {
    let mut b = GeometryGraphBuilder::new();
    let directrix = b.push_value(basis.clone()).unwrap();
    let (small, large) = (disk(&mut b, 0.5), disk(&mut b, 0.8));
    let root = spine(
        &mut b,
        directrix,
        vec![section(small, run.0), section(large, run.1)],
    );
    let graph = b.finish(vec![root]).unwrap();
    let frame_at = move |s: Scalar| station_section2(&basis, s).unwrap();
    // The equal-fraction blend of two circles is the circle of the blended
    // radius.
    let points = exact_surface(&frame_at, run, &|u| circle_points(0.5 + 0.3 * u, 360), 40);
    assert_holds(&graph, root, &points);
}

#[test]
fn a_disk_growing_along_a_line_is_certified() {
    disk_taper(line(), (0.0, 10.0));
}

#[test]
fn a_disk_growing_along_an_arc_is_certified() {
    disk_taper(arc(12.0), (0.0, 20.0));
}

#[test]
fn a_rounded_rectangle_with_growing_corners_is_certified() {
    let mut b = GeometryGraphBuilder::new();
    let basis = line();
    let directrix = b.push_value(basis.clone()).unwrap();
    let (tight, round) = (rounded(&mut b, 0.1), rounded(&mut b, 0.3));
    let root = spine(
        &mut b,
        directrix,
        vec![section(tight, 0.0), section(round, 10.0)],
    );
    let graph = b.finish(vec![root]).unwrap();
    let frame_at = move |s: Scalar| station_section2(&basis, s).unwrap();
    // Blending the corner arcs and the edges at equal fractions gives the
    // rounded rectangle of the blended radius.
    let exact = |u: Scalar| {
        let (hx, hy, r) = (1.0, 0.7, 0.1 + 0.2 * u);
        let mut points = Vec::new();
        for (cx, cy, from) in [
            (hx - r, hy - r, 0.0),
            (-(hx - r), hy - r, 0.5),
            (-(hx - r), -(hy - r), 1.0),
            (hx - r, -(hy - r), 1.5),
        ] {
            for k in 0..=40 {
                let t = (from + 0.5 * k as Scalar / 40.0) * std::f64::consts::PI;
                points.push(Point2::new(cx + r * t.cos(), cy + r * t.sin()));
            }
        }
        // The straight edges between the arcs.
        let corners = [
            (hx, hy - r, hx, -(hy - r)),
            (-(hx - r), hy, hx - r, hy),
            (-hx, hy - r, -hx, -(hy - r)),
            (-(hx - r), -hy, hx - r, -hy),
        ];
        for (x0, y0, x1, y1) in corners {
            for k in 0..=20 {
                let v = k as Scalar / 20.0;
                points.push(Point2::new(x0 + (x1 - x0) * v, y0 + (y1 - y0) * v));
            }
        }
        points
    };
    let points = exact_surface(&frame_at, (0.0, 10.0), &exact, 40);
    assert_holds(&graph, root, &points);
}

#[test]
fn an_ellipse_becoming_a_disk_is_certified() {
    let mut b = GeometryGraphBuilder::new();
    let basis = line();
    let directrix = b.push_value(basis.clone()).unwrap();
    let ellipse = b
        .push_value(Profile::Ellipse(EllipseProfile {
            semi_axis_x: 1.0,
            semi_axis_y: 0.5,
        }))
        .unwrap();
    let round = disk(&mut b, 0.8);
    let root = spine(
        &mut b,
        directrix,
        vec![section(ellipse, 0.0), section(round, 10.0)],
    );
    let graph = b.finish(vec![root]).unwrap();
    let frame_at = move |s: Scalar| station_section2(&basis, s).unwrap();
    // The blend at equal parameters is the ellipse of the blended axes.
    let exact = |u: Scalar| {
        let (a, c) = (1.0 + (0.8 - 1.0) * u, 0.5 + (0.8 - 0.5) * u);
        (0..360)
            .map(|k| {
                let t = std::f64::consts::TAU * k as Scalar / 360.0;
                Point2::new(a * t.cos(), c * t.sin())
            })
            .collect()
    };
    let points = exact_surface(&frame_at, (0.0, 10.0), &exact, 40);
    assert_holds(&graph, root, &points);
}

#[test]
fn a_mirrored_section_against_an_unmirrored_one_is_refused() {
    // The mirrored square's ring is reversed, so matching vertices by index
    // pairs opposite corners and pinches the loft; the run is refused, not
    // certified with a bound of zero.
    let mut b = GeometryGraphBuilder::new();
    let directrix = b.push_value(line()).unwrap();
    let square = |size: Scalar| {
        Profile::Rectangle(RectangleProfile {
            x: size,
            y: size,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        })
    };
    let plain = b.push_value(square(2.0)).unwrap();
    let mirrored = b
        .push_value(Profile::Derived {
            basis: Box::new(square(1.0)),
            transform: axiolid_core::Transform2::from_scale(Vec2::new(-1.0, 1.0)),
        })
        .unwrap();
    let root = spine(
        &mut b,
        directrix,
        vec![section(mirrored, 0.0), section(plain, 10.0)],
    );
    let graph = b.finish(vec![root]).unwrap();
    let err = compiler()
        .compile_mesh_with_deviation(&graph, root, &options())
        .expect_err("opposite handedness is refused");
    assert!(
        matches!(&err, axiolid_contracts::GeomError::InvalidInput(d) if d.contains("mirrored")),
        "{err:?}"
    );
}
