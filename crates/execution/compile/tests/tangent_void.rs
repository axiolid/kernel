//! A void tangent to a face of its host (#194).
//!
//! A beam [-2, 2] x [-0.1, 0.1] x [0, 0.3] with a cylindrical hole of
//! radius 0.05 across its width, axis along y at x = 0, z = 0.25: the
//! hole's top touches the beam's top face along a line.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Tolerance, Transform3, Vec3};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, Instance, SolidOperation};
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

const CHORD: f64 = 1e-4;

fn compile(hole_z: f64, radius: f64) -> Result<TriMesh, axiolid_contracts::GeomError> {
    compile_span(hole_z, radius, CHORD, 0.0, 0.4)
}

/// The beam with a hole of `radius` at height `hole_z`, chorded to `chord`,
/// its polygon turned by `turn` about the axis, `span` long across the
/// beam's 0.2 width (0.2: its caps flush with the beam's sides).
fn compile_span(
    hole_z: f64,
    radius: f64,
    chord: f64,
    turn: f64,
    span: f64,
) -> Result<TriMesh, axiolid_contracts::GeomError> {
    let mut b = GeometryGraphBuilder::new();
    let rect = b
        .push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: 4.0,
                y: 0.2,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )))
        .unwrap();
    let beam = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: rect,
            direction: Vec3::Z,
            depth: 0.3,
        }))
        .unwrap();
    let circle = b
        .push(GeometryNode::Profile(Profile::Circle(CircleProfile {
            radius,
            thickness: None,
        })))
        .unwrap();
    let rod = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: circle,
            direction: Vec3::Z,
            depth: span,
        }))
        .unwrap();
    // The rod's axis z -> y, its x -> x, its y -> -z: from y = -0.2 to 0.2.
    let hole = b
        .push(GeometryNode::Instance(Instance {
            source: rod,
            transform: Transform3::from_cols(
                Vec3::new(turn.cos(), 0.0, turn.sin()),
                Vec3::new(turn.sin(), 0.0, -turn.cos()),
                Vec3::Y,
                Vec3::new(0.0, -span / 2.0, hole_z),
            ),
        }))
        .unwrap();
    let cut = b
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            operator: BooleanOperator::Difference,
            left: beam,
            right: hole,
        }))
        .unwrap();
    let graph = b.finish(vec![cut]).unwrap();
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(chord)
        .expect("a chord budget");
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(&graph, cut, &options)
}

#[test]
fn a_hole_tangent_to_the_top_face_meshes_closed() {
    let r = 0.05;
    let mesh = compile(0.3 - r, r).expect("the beam meshes");
    let volume = volume_properties(&mesh, Tolerance::new(1e-7, 1e-9).unwrap())
        .expect("a closed solid")
        .signed_volume;
    let exact = 4.0 * 0.2 * 0.3 - core::f64::consts::PI * r * r * 0.2;
    // The hole's polygon lies inside its circle by at most the chord
    // budget all round.
    let bound = core::f64::consts::TAU * r * CHORD * 0.2 + 1e-12;
    assert!((volume - exact).abs() <= bound, "{volume} vs {exact}");
}

/// Over chord budgets, hole sizes, turned polygons and holes flush with
/// the beam's sides, the result is a closed solid with no pinch, or a
/// refusal naming the contact -- never a pinched mesh.
#[test]
fn tangent_holes_never_come_out_pinched() {
    for chord in [1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 1e-4, 1e-5] {
        for r in [0.05, 0.04, 0.1] {
            for (turn, span) in [
                (0.0, 0.4),
                (0.1, 0.4),
                (0.0, 0.2),
                (0.1, 0.2),
                (core::f64::consts::FRAC_PI_4, 0.2),
            ] {
                let case = format!("chord {chord} r {r} turn {turn} span {span}");
                match compile_span(0.3 - r, r, chord, turn, span) {
                    Ok(mesh) => {
                        volume_properties(&mesh, Tolerance::new(1e-7, 1e-9).unwrap())
                            .unwrap_or_else(|e| panic!("{case}: {e}"));
                        assert_eq!(welded_non_manifold(&mesh), 0, "{case}");
                    }
                    Err(axiolid_contracts::GeomError::Degenerate(_)) => {}
                    Err(e) => panic!("{case}: {e}"),
                }
            }
        }
    }
}

/// Edges shared by more than two triangles once coincident positions are
/// welded, as a consumer that validates by position sees them.
fn welded_non_manifold(mesh: &TriMesh) -> usize {
    use std::collections::HashMap;
    let key = |p: axiolid_core::Point3| (p.x.to_bits(), p.y.to_bits(), p.z.to_bits());
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

#[test]
fn the_tangent_hole_leaves_no_pinch() {
    // Whatever the chord budget or the hole's size, the chorded circle
    // keeps a sliver of material under the face it is tangent to.
    for chord in [1e-2, 1e-3, 1e-4] {
        for r in [0.04, 0.05, 0.1] {
            let mesh = compile_span(0.3 - r, r, chord, 0.0, 0.2).expect("the beam meshes");
            assert_eq!(welded_non_manifold(&mesh), 0, "chord {chord} r {r}");
        }
    }
}

/// A void whose corner touches the top face exactly: a diamond across the
/// beam, its top corner at z = 0.3. No chord slack can part them, so the
/// pinch is refused with the contact named.
#[test]
fn an_exact_pinch_is_refused_naming_the_contact() {
    use axiolid_core::{Interval, Point2, Vec2};
    use axiolid_curve::{Curve2, Line2};
    use axiolid_profile::{Contour, ContourProfile, ProfileSegment};
    let corners = [
        Point2::new(0.0, -0.05),
        Point2::new(0.05, 0.0),
        Point2::new(0.0, 0.05),
        Point2::new(-0.05, 0.0),
    ];
    let segments = (0..4)
        .map(|k| {
            let (a, b) = (corners[k], corners[(k + 1) % 4]);
            ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: a,
                    direction: Vec2::new(b.x - a.x, b.y - a.y),
                }),
                domain: Interval {
                    start: 0.0,
                    end: 1.0,
                },
                same_sense: true,
            }
        })
        .collect();
    let diamond = Profile::Contour(ContourProfile {
        outer: Contour { segments },
        holes: Vec::new(),
    });
    let mut b = GeometryGraphBuilder::new();
    let rect = b
        .push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: 4.0,
                y: 0.2,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )))
        .unwrap();
    let beam = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: rect,
            direction: Vec3::Z,
            depth: 0.3,
        }))
        .unwrap();
    let profile = b.push(GeometryNode::Profile(diamond)).unwrap();
    let prism = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 0.4,
        }))
        .unwrap();
    // Axis along y; the diamond's local -y corner up, at z = 0.25 + 0.05.
    let hole = b
        .push(GeometryNode::Instance(Instance {
            source: prism,
            transform: Transform3::from_cols(
                Vec3::X,
                Vec3::new(0.0, 0.0, -1.0),
                Vec3::Y,
                Vec3::new(0.0, -0.2, 0.25),
            ),
        }))
        .unwrap();
    let cut = b
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            operator: BooleanOperator::Difference,
            left: beam,
            right: hole,
        }))
        .unwrap();
    let graph = b.finish(vec![cut]).unwrap();
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(CHORD)
        .unwrap();
    match ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(&graph, cut, &options) {
        Err(axiolid_contracts::GeomError::Degenerate(message)) => {
            assert!(message.contains("touches itself"), "{message}");
        }
        Ok(mesh) => panic!(
            "a pinched result was returned, {} welded edges with more than two faces",
            welded_non_manifold(&mesh)
        ),
        Err(other) => panic!("{other:?}"),
    }
}

/// A void touching the face at one point: a pyramid whose apex meets the
/// top face. Refused naming the point, or returned with no pinch.
#[test]
fn a_void_touching_at_a_point_is_not_returned_pinched() {
    use axiolid_primitive::Primitive;
    let mut b = GeometryGraphBuilder::new();
    let rect = b
        .push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: 4.0,
                y: 0.2,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )))
        .unwrap();
    let beam = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: rect,
            direction: Vec3::Z,
            depth: 0.3,
        }))
        .unwrap();
    let pyramid = b
        .push(GeometryNode::Primitive(Primitive::Pyramid {
            x: 0.1,
            y: 0.1,
            height: 0.1,
        }))
        .unwrap();
    let hole = b
        .push(GeometryNode::Instance(Instance {
            source: pyramid,
            transform: Transform3::from_translation(Vec3::new(0.5, 0.0, 0.2)),
        }))
        .unwrap();
    let cut = b
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            operator: BooleanOperator::Difference,
            left: beam,
            right: hole,
        }))
        .unwrap();
    let graph = b.finish(vec![cut]).unwrap();
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(CHORD)
        .unwrap();
    match ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(&graph, cut, &options) {
        Err(axiolid_contracts::GeomError::Degenerate(message)) => {
            assert!(message.contains("touches itself"), "{message}");
            eprintln!("refused: {message}");
        }
        Ok(mesh) => {
            eprintln!("meshed");
            assert_eq!(welded_non_manifold(&mesh), 0);
        }
        Err(other) => panic!("{other:?}"),
    }
}
