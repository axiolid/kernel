//! A bounded half-space's slab is sized by the subject it meets (#325).
//!
//! A bounded half-space is unbounded on its kept side. The mesh compiler
//! used to build it standalone, a slab twice its boundary's extent deep,
//! so a subject reaching further than that along the sweep was cut too
//! little by a difference, or kept too little by an intersection, with no
//! error. Its difference and intersection operands are now built from the
//! subject, through any instances above them, as an unbounded half-space's
//! are; with no subject (a root, a union, a collection member) it is
//! refused by name.
//!
//! The case: a 2 x 2 x 20 column whose top 19 lie above the plane `z = 1`,
//! cut by the prism of a 0.5 x 0.5 square, along either axis, both
//! compilers, both operators, plain and placed by an instance.

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Plane3, Point2, Point3, Scalar, Tolerance, Transform3, Vec3};
use axiolid_curve::{Curve2, Polyline2};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::exact_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{ReferenceExactCompiler, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{Profile, RectangleProfile};

/// Which variant states the clip: the prism along the plane's normal
/// (`BoundedHalfSpace`) or along the placement's `z`
/// (`BoundedHalfSpaceAlongZ`, #321).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PrismAxis {
    PlaneNormal,
    PlacementZ,
}

impl PrismAxis {
    fn operation(
        self,
        half_space: NodeId,
        boundary: NodeId,
        placement: Transform3,
    ) -> SolidOperation {
        match self {
            Self::PlaneNormal => SolidOperation::BoundedHalfSpace {
                half_space,
                boundary,
                placement,
            },
            Self::PlacementZ => SolidOperation::BoundedHalfSpaceAlongZ {
                half_space,
                boundary,
                placement,
            },
        }
    }
}

/// The column's height and the plane's.
const HEIGHT: Scalar = 20.0;
const PLANE: Scalar = 1.0;
/// The square's side: a quarter of the column's plan.
const SIDE: Scalar = 0.5;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn general() -> Transform3 {
    Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7)
}

struct Graph {
    builder: GeometryGraphBuilder,
}

impl Graph {
    fn push(&mut self, node: GeometryNode) -> NodeId {
        self.builder.push(node).expect("a valid node")
    }

    /// The column under `placement`.
    fn column(&mut self, placement: Transform3) -> NodeId {
        let profile = self.push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: 2.0,
                y: 2.0,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )));
        let column = self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: HEIGHT,
        }));
        self.push(GeometryNode::Instance(Instance {
            source: column,
            transform: placement,
        }))
    }

    /// The side above the plane within the square's prism along `axis`,
    /// under an instance of `placement`: stated with the normal up and its
    /// side kept, or with the normal down and the other side kept.
    fn square_clip(&mut self, axis: PrismAxis, placement: Transform3, up: bool) -> NodeId {
        let half_space = self.push(GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 {
                origin: Point3::new(0.0, 0.0, PLANE),
                normal: if up { Vec3::Z } else { -Vec3::Z },
            },
            agreement: up,
        }));
        let h = SIDE / 2.0;
        let boundary = self.push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![
                Point2::new(-h, -h),
                Point2::new(h, -h),
                Point2::new(h, h),
                Point2::new(-h, h),
            ],
            closed: true,
        })));
        let bounded = self.push(GeometryNode::SolidOperation(axis.operation(
            half_space,
            boundary,
            Transform3::IDENTITY,
        )));
        self.push(GeometryNode::Instance(Instance {
            source: bounded,
            transform: placement,
        }))
    }

    fn boolean(&mut self, left: NodeId, right: NodeId, operator: BooleanOperator) -> NodeId {
        self.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator,
        }))
    }
}

fn graph(
    operator: BooleanOperator,
    axis: PrismAxis,
    placement: Transform3,
    up: bool,
) -> (GeometryGraph, NodeId) {
    let mut g = Graph {
        builder: GeometryGraphBuilder::new(),
    };
    let column = g.column(placement);
    let clip = g.square_clip(axis, placement, up);
    let root = g.boolean(column, clip, operator);
    (g.builder.finish(vec![root]).expect("a valid graph"), root)
}

fn mesh_volume(mesh: &TriMesh) -> Scalar {
    let base = mesh.positions[0];
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize] - base);
            a.dot(b.cross(c))
        })
        .sum::<Scalar>()
        / 6.0
}

/// The square's prism above the plane through the whole column: 19 high,
/// where a slab sized by the square would have reached 0.5.
#[test]
fn a_small_boundary_cuts_a_tall_subject_to_its_full_height() {
    let cut = SIDE * SIDE * (HEIGHT - PLANE);
    for (placement, up) in [
        (Transform3::IDENTITY, true),
        (Transform3::IDENTITY, false),
        (general(), true),
    ] {
        for axis in [PrismAxis::PlaneNormal, PrismAxis::PlacementZ] {
            for (operator, expected) in [
                (BooleanOperator::Difference, 4.0 * HEIGHT - cut),
                (BooleanOperator::Intersection, cut),
            ] {
                let (graph, root) = graph(operator, axis, placement, up);
                let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
                    .compile_mesh(&graph, root, &options())
                    .unwrap_or_else(|error| panic!("{axis:?} {operator:?}: {error:?}"));
                let measured = mesh_volume(&mesh);
                assert!(
                    (measured - expected).abs() < 1e-9 * HEIGHT,
                    "{axis:?} {operator:?}: mesh {measured}, expected {expected}"
                );
                let exact = ReferenceExactCompiler::new()
                    .compile_exact(&graph, root, &options())
                    .unwrap_or_else(|error| panic!("{axis:?} {operator:?}: {error:?}"));
                let exact = exact_properties(&exact, Tolerance::METRE)
                    .expect("measurable")
                    .signed_volume;
                assert!(
                    (exact - expected).abs() < 1e-9 * HEIGHT,
                    "{axis:?} {operator:?}: exact {exact}, expected {expected}"
                );
            }
        }
    }
}

fn refused_without_a_subject(graph: &GeometryGraph, root: NodeId, case: &str) {
    let error = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(graph, root, &options())
        .expect_err("refused");
    assert!(
        matches!(&error, GeomError::UnsupportedInput { input, .. } if input.contains("without a subject")),
        "{case}: {error:?}"
    );
}

/// Where nothing sizes the slab, the mesh compiler refuses by name, as the
/// exact compiler does: a bounded half-space as a root, through an
/// instance, and in a union.
#[test]
fn a_bounded_half_space_without_a_subject_is_refused_by_name() {
    for axis in [PrismAxis::PlaneNormal, PrismAxis::PlacementZ] {
        let mut g = Graph {
            builder: GeometryGraphBuilder::new(),
        };
        let clip = g.square_clip(axis, general(), true);
        let column = g.column(Transform3::IDENTITY);
        let union = g.boolean(column, clip, BooleanOperator::Union);
        let graph = g.builder.finish(vec![clip, union]).expect("a valid graph");
        refused_without_a_subject(&graph, clip, "a placed root");
        refused_without_a_subject(&graph, union, "a union");
        let exact = ReferenceExactCompiler::new()
            .compile_exact(&graph, clip, &options())
            .expect_err("refused");
        assert!(
            matches!(&exact, GeomError::UnsupportedInput { input, .. } if input.contains("half-space")),
            "{exact:?}"
        );
    }
}

/// openbimrs/ifc's repro (`tall_wall_small_bounded_clip_ifc2x3.ifc`, its
/// numbers as written, in millimetres): a 16.165 m wall less what lies above
/// `z = 11404.75` within an aligned 238.5 x 350 boundary at its end
/// (`AgreementFlag .F.`: the side the normal points to is the half-space).
/// Along the plane normal, so through `BoundedHalfSpace`. Sized by its
/// boundary, the slab was 0.7 m deep and the wall stayed 4.06 m too tall
/// there: 57.226391331354 m³, a certified bound of 0.175 m that missed the
/// request. Sized by the wall, it is the exact and analytic
/// 56.887388852960 m³, certified within the request.
#[test]
fn the_tall_wall_under_a_small_aligned_clip_is_cut_to_its_full_height() {
    const EXPECTED: Scalar = 56.887_388_852_960;
    let mut g = Graph {
        builder: GeometryGraphBuilder::new(),
    };
    let profile = g.push(GeometryNode::Profile(Profile::Rectangle(
        RectangleProfile {
            x: 10125.0,
            y: 350.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        },
    )));
    let wall = g.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: 16165.033064627303,
    }));
    // `#32`: centred at (5062.5, 0), its x along -X.
    let wall = g.push(GeometryNode::Instance(Instance {
        source: wall,
        transform: Transform3::from_cols(-Vec3::X, -Vec3::Y, Vec3::Z, Vec3::new(5062.5, 0.0, 0.0)),
    }));
    let origin = Point3::new(9886.45064727072, -175.0, 11404.747431413967);
    let half_space = g.push(GeometryNode::HalfSpace(HalfSpace {
        boundary: Plane3 {
            origin,
            normal: Vec3::Z,
        },
        agreement: true,
    }));
    let boundary = g.push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(238.54935272928, 0.0),
            Point2::new(238.54935272928, 350.0),
            Point2::new(0.0, 350.0),
            Point2::new(0.0, 0.0),
        ],
        closed: true,
    })));
    let clip = g.push(GeometryNode::SolidOperation(
        SolidOperation::BoundedHalfSpace {
            half_space,
            boundary,
            placement: Transform3::from_translation(origin),
        },
    ));
    let root = g.boolean(wall, clip, BooleanOperator::Difference);
    let graph = g.builder.finish(vec![root]).expect("a valid graph");
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options)
        .expect("the clipped wall compiles");
    let mesh = mesh_volume(&outcome.mesh) * 1e-9;
    let exact = ReferenceExactCompiler::new()
        .compile_exact(&graph, root, &options)
        .expect("the exact result");
    let exact = exact_properties(&exact, Tolerance::MILLIMETRE)
        .expect("measurable")
        .signed_volume
        * 1e-9;
    eprintln!("tall wall: mesh {mesh:.12} m³, exact {exact:.12} m³, {report:?}");
    for (path, volume) in [("mesh", mesh), ("exact", exact)] {
        assert!(
            (volume - EXPECTED).abs() < 1e-9,
            "{path}: {volume} m³, expected {EXPECTED}"
        );
    }
    assert!(report.meets_requested(), "{report:?}");
}
