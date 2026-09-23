//! Channels and normals survive `Instance` and `Collection` (#115).
//!
//! Before #115 both nodes rebuilt the mesh from positions and indices only,
//! so a textured item lost its `uv` channel the moment a product had a second
//! item or was instanced -- silently, with no report. These tests pin each
//! path through the public compiler, and pin the report that now says what
//! happened to every channel.
//!
//! The fixture is a unit tetrahedron whose `uv` channel is corner-indexed
//! with a DIFFERENT value on every corner (twelve values, four triangles), so
//! any shuffle of corners -- a mirror that swaps winding without swapping
//! the channel, an offset applied to the wrong member -- reads a wrong value
//! rather than a coincidentally right one.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Point3, Tolerance, Transform3, Vec3};
use axiolid_mesh::{AttributeChannel, AttributeFate, Blend, DropReason, NormalAttribute, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId};

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

/// Outward-wound unit tetrahedron, no channels.
fn tetrahedron() -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ],
        vec![0, 2, 1, 0, 1, 3, 1, 2, 3, 0, 3, 2],
    )
}

/// The value corner `c` carries: `(c, 100 + c)`, unique per corner.
fn corner_uv(c: usize) -> [f64; 2] {
    [c as f64, 100.0 + c as f64]
}

/// Tetrahedron with a corner-indexed `uv` that differs on every corner.
fn textured() -> TriMesh {
    let mut mesh = tetrahedron();
    let values = (0..12).flat_map(corner_uv).collect();
    mesh.attributes.push(AttributeChannel::corner_indexed(
        "uv",
        values,
        2,
        Blend::Linear,
        (0..12).collect(),
    ));
    mesh
}

/// The `uv` value at corner `c` of `mesh`.
fn uv_at(mesh: &TriMesh, c: usize) -> Option<[f64; 2]> {
    let channel = mesh.attributes.iter().find(|ch| ch.name == "uv")?;
    let v = channel.at_corner(&mesh.indices, c)?;
    Some([v[0], v[1]])
}

fn graph(
    nodes: Vec<GeometryNode>,
    root: impl Fn(&[NodeId]) -> GeometryNode,
) -> (GeometryGraph, NodeId) {
    let mut builder = GeometryGraphBuilder::new();
    let ids: Vec<NodeId> = nodes
        .into_iter()
        .map(|n| builder.push(n).expect("push"))
        .collect();
    let root = builder.push(root(&ids)).expect("push root");
    (builder.finish(vec![root]).expect("finish"), root)
}

// ---- Instance ------------------------------------------------------------

/// A plain translation keeps every value on its own corner.
#[test]
fn an_instance_keeps_the_channel_on_every_corner() {
    let (g, root) = graph(vec![GeometryNode::TriMesh(textured())], |ids| {
        GeometryNode::Instance(Instance {
            source: ids[0],
            transform: Transform3::from_translation(Vec3::new(5.0, 0.0, 0.0)),
        })
    });
    let out = compiler()
        .compile_mesh(&g, root, &options())
        .expect("compiles");
    out.validate_structure().expect("valid");
    for c in 0..12 {
        assert_eq!(uv_at(&out, c), Some(corner_uv(c)), "corner {c}");
    }
    assert_eq!(out.positions[0], Point3::new(5.0, 0.0, 0.0), "moved");
}

/// A mirror swaps corners 1 and 2 of every triangle to stay outward. The
/// channel must swap with them, or corner 1 reads corner 2's value.
#[test]
fn a_mirroring_instance_swaps_channel_corners_with_the_triangle() {
    let mirror = Transform3::from_scale(Vec3::new(-1.0, 1.0, 1.0));
    let (g, root) = graph(vec![GeometryNode::TriMesh(textured())], |ids| {
        GeometryNode::Instance(Instance {
            source: ids[0],
            transform: mirror,
        })
    });
    let out = compiler()
        .compile_mesh(&g, root, &options())
        .expect("compiles");
    out.validate_structure().expect("valid");
    let source = textured();
    for t in 0..4 {
        for (corner, from) in [(0, 0), (1, 2), (2, 1)] {
            let (c, f) = (3 * t + corner, 3 * t + from);
            assert_eq!(
                out.indices[c], source.indices[f],
                "triangle {t}: winding flipped"
            );
            assert_eq!(
                uv_at(&out, c),
                Some(corner_uv(f)),
                "triangle {t} corner {corner}"
            );
        }
    }
}

/// Normals were dropped on every instance before #115. They take the
/// inverse transpose, so under non-uniform scale they stay perpendicular to
/// the transformed surface instead of skewing with it.
#[test]
fn an_instance_transforms_normals_by_the_inverse_transpose() {
    let mut mesh = tetrahedron();
    // A normal along (1,1,0)/sqrt2, one per vertex.
    let n = Vec3::new(1.0, 1.0, 0.0).normalize();
    mesh.normals = Some(NormalAttribute {
        values: vec![n; 4],
        indices: None,
    });
    let stretch = Transform3::from_scale(Vec3::new(2.0, 1.0, 1.0));
    let (g, root) = graph(vec![GeometryNode::TriMesh(mesh)], |ids| {
        GeometryNode::Instance(Instance {
            source: ids[0],
            transform: stretch,
        })
    });
    let out = compiler()
        .compile_mesh(&g, root, &options())
        .expect("compiles");
    let got = out.normals.expect("normals kept").values[0];
    // Inverse transpose of diag(2,1,1) is diag(0.5,1,1).
    let want = Vec3::new(0.5, 1.0, 0.0).normalize();
    assert!((got - want).length() < 1e-12, "got {got:?}, want {want:?}");
}

// ---- Collection ----------------------------------------------------------

/// The reproduction from #115: a textured item plus an untextured one. The
/// channel stays; the untextured item's triangles are explicitly unmapped.
#[test]
fn a_collection_keeps_a_channel_only_some_members_have() {
    let (g, root) = graph(
        vec![
            GeometryNode::TriMesh(textured()),
            GeometryNode::TriMesh(tetrahedron()),
        ],
        |ids| GeometryNode::Collection(ids.to_vec()),
    );
    let out = compiler()
        .compile_mesh(&g, root, &options())
        .expect("compiles");
    out.validate_structure().expect("valid");
    assert_eq!(out.indices.len(), 24, "both members' triangles");
    for c in 0..12 {
        assert_eq!(uv_at(&out, c), Some(corner_uv(c)), "textured corner {c}");
    }
    for c in 12..24 {
        assert_eq!(uv_at(&out, c), None, "untextured corner {c} is unmapped");
    }
}

/// Both members textured: the second member's corners must address the
/// second member's values, i.e. its corner indices are offset into the
/// merged pool. Without the offset it would read the first member's.
#[test]
fn a_collection_offsets_the_second_members_corners_into_the_merged_pool() {
    let mut second = textured();
    let uv = second
        .attributes
        .iter_mut()
        .find(|c| c.name == "uv")
        .unwrap();
    for v in &mut uv.values {
        *v += 1000.0;
    }
    let (g, root) = graph(
        vec![
            GeometryNode::TriMesh(textured()),
            GeometryNode::TriMesh(second),
        ],
        |ids| GeometryNode::Collection(ids.to_vec()),
    );
    let out = compiler()
        .compile_mesh(&g, root, &options())
        .expect("compiles");
    out.validate_structure().expect("valid");
    for c in 0..12 {
        let [s, t] = corner_uv(c);
        assert_eq!(
            uv_at(&out, 12 + c),
            Some([s + 1000.0, t + 1000.0]),
            "corner {c}"
        );
    }
}

/// Per-vertex in every member stays per-vertex: no needless conversion.
#[test]
fn a_collection_of_per_vertex_channels_stays_per_vertex() {
    let per_vertex = |base: f64| {
        let mut mesh = tetrahedron();
        mesh.attributes.push(AttributeChannel::new(
            "temperature",
            vec![base, base + 1.0, base + 2.0, base + 3.0],
            1,
            Blend::Linear,
        ));
        mesh
    };
    let (g, root) = graph(
        vec![
            GeometryNode::TriMesh(per_vertex(0.0)),
            GeometryNode::TriMesh(per_vertex(10.0)),
        ],
        |ids| GeometryNode::Collection(ids.to_vec()),
    );
    let out = compiler()
        .compile_mesh(&g, root, &options())
        .expect("compiles");
    out.validate_structure().expect("valid");
    let channel = &out.attributes[0];
    assert!(!channel.is_corner_indexed());
    assert_eq!(
        channel.values,
        vec![0.0, 1.0, 2.0, 3.0, 10.0, 11.0, 12.0, 13.0]
    );
}

/// Two members defining `uv` with different widths cannot share one
/// channel. It is dropped -- and, the #84 contract, reported.
#[test]
fn incompatible_widths_are_dropped_and_reported() {
    let mut wide = tetrahedron();
    wide.attributes
        .push(AttributeChannel::new("uv", vec![0.0; 12], 3, Blend::Linear));
    let (g, root) = graph(
        vec![
            GeometryNode::TriMesh(textured()),
            GeometryNode::TriMesh(wide),
        ],
        |ids| GeometryNode::Collection(ids.to_vec()),
    );
    let out = compiler()
        .compile_mesh_reported(&g, root, &options())
        .expect("compiles");
    assert!(
        out.mesh.attributes.iter().all(|c| c.name != "uv"),
        "dropped"
    );
    assert_eq!(
        out.fate("uv"),
        Some(&AttributeFate::Dropped(DropReason::IncompatibleChannels))
    );
}

// ---- Nesting and the report ---------------------------------------------

/// The IFC shape: a mapped item (Instance) inside a multi-item product
/// (Collection), mirrored. Every value must land on its own corner.
#[test]
fn a_mirrored_instance_inside_a_collection_keeps_every_value() {
    let (g, root) = {
        let mut b = GeometryGraphBuilder::new();
        let src = b.push(GeometryNode::TriMesh(textured())).unwrap();
        let mirrored = b
            .push(GeometryNode::Instance(Instance {
                source: src,
                transform: Transform3::from_scale(Vec3::new(1.0, -1.0, 1.0)),
            }))
            .unwrap();
        let plain = b.push(GeometryNode::TriMesh(tetrahedron())).unwrap();
        let root = b
            .push(GeometryNode::Collection(vec![plain, mirrored]))
            .unwrap();
        (b.finish(vec![root]).unwrap(), root)
    };
    let out = compiler()
        .compile_mesh(&g, root, &options())
        .expect("compiles");
    out.validate_structure().expect("valid");
    for c in 0..12 {
        assert_eq!(uv_at(&out, c), None, "plain member unmapped");
    }
    for t in 0..4 {
        for (corner, from) in [(0, 0), (1, 2), (2, 1)] {
            let c = 12 + 3 * t + corner;
            assert_eq!(
                uv_at(&out, c),
                Some(corner_uv(3 * t + from)),
                "t{t} c{corner}"
            );
        }
    }
}

/// Every channel that reached the root unchanged is reported Preserved.
/// The report used not to exist; a caller could not tell dropped from
/// never-there.
#[test]
fn the_report_names_every_preserved_channel() {
    let (g, root) = graph(
        vec![
            GeometryNode::TriMesh(textured()),
            GeometryNode::TriMesh(tetrahedron()),
        ],
        |ids| GeometryNode::Collection(ids.to_vec()),
    );
    let out = compiler()
        .compile_mesh_reported(&g, root, &options())
        .expect("compiles");
    assert_eq!(out.fate("uv"), Some(&AttributeFate::Preserved));
    assert_eq!(out.attribute_fates.as_ref().map(Vec::len), Some(1));
}

/// A compiler that does not override `compile_mesh_reported` says it does
/// not know, rather than claiming nothing was lost.
#[test]
fn the_default_report_is_unknown_not_empty() {
    #[derive(Debug)]
    struct Plain(ReferenceMeshCompiler<BoolmeshBoolean>);
    impl axiolid_contracts::Backend for Plain {
        fn descriptor(&self) -> axiolid_contracts::BackendDescriptor {
            self.0.descriptor()
        }
    }
    impl MeshCompiler for Plain {
        fn compile_mesh(
            &self,
            graph: &GeometryGraph,
            root: NodeId,
            options: &ExecutionOptions,
        ) -> axiolid_contracts::GeomResult<TriMesh> {
            self.0.compile_mesh(graph, root, options)
        }
    }
    let (g, root) = graph(vec![GeometryNode::TriMesh(textured())], |ids| {
        GeometryNode::Collection(ids.to_vec())
    });
    let out = Plain(compiler())
        .compile_mesh_reported(&g, root, &options())
        .expect("compiles");
    assert_eq!(out.attribute_fates, None);
    assert!(out.mesh.attributes.iter().any(|c| c.name == "uv"));
}

// ---- Composition with the boolean ---------------------------------------

/// Axis-aligned box `[lo, hi]`, outward-wound.
fn cuboid(lo: [f64; 3], hi: [f64; 3]) -> TriMesh {
    let p = |x: usize, y: usize, z: usize| {
        Point3::new([lo[0], hi[0]][x], [lo[1], hi[1]][y], [lo[2], hi[2]][z])
    };
    let positions = vec![
        p(0, 0, 0),
        p(1, 0, 0),
        p(1, 1, 0),
        p(0, 1, 0),
        p(0, 0, 1),
        p(1, 0, 1),
        p(1, 1, 1),
        p(0, 1, 1),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, // z = lo
        4, 5, 6, 4, 6, 7, // z = hi
        0, 1, 5, 0, 5, 4, // y = lo
        2, 3, 7, 2, 7, 6, // y = hi
        1, 2, 6, 1, 6, 5, // x = hi
        0, 4, 7, 0, 7, 3, // x = lo
    ];
    TriMesh::new(positions, indices)
}

/// A boolean downstream of an Instance: the boolean's fate is composed onto
/// the Instance's, not replacing it and not replaced by it. The tag channel
/// forbids derivation (`Blend::None`), the cut creates vertices, so the
/// boolean drops it `NotBlendable` -- and that must be what the root
/// reports, even though the Instance preserved it. (A Collection is not a
/// legal boolean operand; the graph rejects it.)
#[test]
fn a_boolean_after_an_instance_reports_the_boolean_drop() {
    let mut tagged = cuboid([0.0; 3], [2.0; 3]);
    tagged
        .attributes
        .push(AttributeChannel::new("tag", vec![7.0; 8], 1, Blend::None));
    let (g, root) = {
        let mut b = GeometryGraphBuilder::new();
        let a = b.push(GeometryNode::TriMesh(tagged)).unwrap();
        let subject = b
            .push(GeometryNode::Instance(Instance {
                source: a,
                transform: Transform3::IDENTITY,
            }))
            .unwrap();
        let tool = b
            .push(GeometryNode::TriMesh(cuboid([1.0; 3], [3.0; 3])))
            .unwrap();
        let root = b
            .push(GeometryNode::SolidOperation(
                axiolid_model::SolidOperation::Boolean {
                    left: subject,
                    right: tool,
                    operator: axiolid_core::BooleanOperator::Difference,
                },
            ))
            .unwrap();
        (b.finish(vec![root]).unwrap(), root)
    };
    let out = compiler()
        .compile_mesh_reported(&g, root, &ExecutionOptions::new(Tolerance::METRE))
        .expect("compiles");
    assert_eq!(
        out.fate("tag"),
        Some(&AttributeFate::Dropped(DropReason::NotBlendable))
    );
}
