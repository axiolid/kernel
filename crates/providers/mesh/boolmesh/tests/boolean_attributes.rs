//! A boolean carries attribute channels, or says why it could not (#116).
//!
//! Every result triangle lies inside ONE operand triangle, and the CSG core
//! records which. So a channel's value at any result corner is fixed by the
//! source triangle: copied where the corner is a source corner, derived
//! under the channel's own `Blend` where it lies on a cut.
//!
//! The oracle here is an independent one: a LINEAR field over space. Its
//! exact value at any point is known in closed form, and barycentric
//! interpolation inside a planar triangle reproduces a linear field
//! exactly. So every result corner can be checked against the formula,
//! not against the implementation's own arithmetic.

mod support;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Tolerance};
use axiolid_mesh::{AttributeChannel, AttributeFate, Blend, DropReason, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
use support::boxx;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn subject() -> TriMesh {
    boxx(0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 0.0)
}

fn tool() -> TriMesh {
    boxx(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 0.0)
}

fn run(
    subject: &TriMesh,
    tool: &TriMesh,
    op: BooleanOperator,
) -> axiolid_mesh_boolean_contract::BooleanOutcome {
    BoolmeshBoolean::new()
        .boolean(subject, tool, op, &options())
        .expect("boolean")
}

/// A linear field over space: exact under barycentric interpolation.
fn field(p: Point3) -> [f64; 2] {
    [0.5 * p.x - 0.25 * p.y + 2.0, p.z + 0.125 * p.x]
}

/// Attach `field` as a per-vertex channel.
fn with_field(mut mesh: TriMesh, name: &str) -> TriMesh {
    let values = mesh.positions.iter().flat_map(|&p| field(p)).collect();
    mesh.attributes
        .push(AttributeChannel::new(name, values, 2, Blend::Linear));
    mesh
}

/// Every mapped corner of `channel` on `mesh` equals `field` at its position.
fn assert_matches_field(mesh: &TriMesh, channel: &AttributeChannel) {
    let mut mapped = 0;
    for (corner, &vertex) in mesh.indices.iter().enumerate() {
        let Some(value) = channel.at_corner(&mesh.indices, corner) else {
            continue;
        };
        mapped += 1;
        let want = field(mesh.positions[vertex as usize]);
        for d in 0..2 {
            assert!(
                (value[d] - want[d]).abs() < 1e-12,
                "corner {corner}: {value:?} vs field {want:?}"
            );
        }
    }
    assert!(mapped > 0, "no corner was mapped, so nothing was checked");
}

/// Both operands carry the field: every result corner, cut or not, reads it.
#[test]
fn a_linear_channel_is_interpolated_exactly_across_every_cut() {
    for op in [
        BooleanOperator::Difference,
        BooleanOperator::Union,
        BooleanOperator::Intersection,
    ] {
        let a = with_field(subject(), "uv");
        let b = with_field(tool(), "uv");
        let outcome = run(&a, &b, op);
        let mesh = &outcome.mesh;
        mesh.validate_structure().expect("channel is well formed");
        assert_eq!(
            outcome.evidence.attribute_fates,
            vec![("uv".to_owned(), AttributeFate::Interpolated)],
            "{op:?}"
        );
        let channel = &mesh.attributes[0];
        assert!(
            channel.is_corner_indexed(),
            "{op:?}: cut vertices are shared"
        );
        let unmapped = channel.corner_indices.as_ref().map_or(0, |c| {
            c.iter()
                .filter(|&&i| i == AttributeChannel::UNMAPPED)
                .count()
        });
        assert_eq!(unmapped, 0, "{op:?}: both operands carry it");
        assert_matches_field(mesh, channel);
    }
}

/// Faces from a tool WITHOUT the channel are unmapped, never zero-filled.
#[test]
fn tool_faces_without_the_channel_are_unmapped() {
    let a = with_field(subject(), "uv");
    let outcome = run(&a, &tool(), BooleanOperator::Difference);
    let mesh = &outcome.mesh;
    mesh.validate_structure()
        .expect("unmapped triangles are whole");
    let channel = &mesh.attributes[0];
    let corners = channel.corner_indices.as_ref().expect("corner-indexed");
    let unmapped = corners
        .iter()
        .filter(|&&i| i == AttributeChannel::UNMAPPED)
        .count();
    // The notch's walls come from the tool, so some faces MUST be unmapped.
    assert!(unmapped > 0, "tool faces must be unmapped");
    assert_eq!(unmapped % 3, 0, "whole triangles only");
    assert_matches_field(mesh, channel);
}

/// A UV seam: one position, a different value on each side. It must come
/// through exactly: each face keeps its own side's value, copied.
#[test]
fn a_corner_indexed_seam_survives_the_cut_exactly() {
    let mut a = subject();
    // Every corner its own tag, so every vertex is a seam: each face meeting
    // there carries a different value at it.
    let n = a.indices.len() as u32;
    a.attributes.push(AttributeChannel::corner_indexed(
        "tag",
        (0..n).map(f64::from).collect(),
        1,
        Blend::Nearest,
        (0..n).collect(),
    ));
    a.validate_structure().expect("valid operand");

    let outcome = run(&a, &tool(), BooleanOperator::Difference);
    let mesh = &outcome.mesh;
    mesh.validate_structure().expect("valid result");
    let channel = &mesh.attributes[0];

    // (-1,-1,0) is a subject corner the tool never reaches.
    let far = Point3::new(-1.0, -1.0, 0.0);
    let mut tags = Vec::new();
    for (corner, &vertex) in mesh.indices.iter().enumerate() {
        if mesh.positions[vertex as usize] != far {
            continue;
        }
        let tag = channel
            .at_corner(&mesh.indices, corner)
            .expect("subject face");
        let source_corner = tag[0] as usize;
        // Copied, not borrowed from a neighbour: the tag names a source
        // corner at exactly this position.
        assert_eq!(
            a.positions[a.indices[source_corner] as usize], far,
            "corner {corner} took tag {source_corner} from elsewhere"
        );
        tags.push(source_corner);
    }
    tags.sort_unstable();
    tags.dedup();
    // The seam itself: several faces meet here, each with its own value.
    // Flattening to one value per position would leave exactly one tag.
    assert!(tags.len() >= 2, "seam collapsed to {tags:?}");
}

/// `Blend::None` forbids deriving a value, and a cut needs one.
#[test]
fn a_non_blendable_channel_is_dropped_when_the_cut_needs_a_value() {
    let mut a = subject();
    let n = a.positions.len();
    a.attributes.push(AttributeChannel::new(
        "opaque",
        (0..n).map(|i| i as f64).collect(),
        1,
        Blend::None,
    ));
    let outcome = run(&a, &tool(), BooleanOperator::Difference);
    assert_eq!(
        outcome.evidence.attribute_fates,
        vec![(
            "opaque".to_owned(),
            AttributeFate::Dropped(DropReason::NotBlendable)
        )]
    );
    assert!(outcome.mesh.attributes.is_empty(), "dropped means absent");
}

/// A disjoint tool leaves the subject untouched: nothing derived, so the
/// channel is `Preserved` and every value is the input's own.
#[test]
fn an_untouched_subject_preserves_its_channel() {
    let a = with_field(subject(), "uv");
    let far = boxx(10.0, 10.0, 0.0, 1.0, 1.0, 1.0, 0.0);
    let outcome = run(&a, &far, BooleanOperator::Difference);
    assert_eq!(
        outcome.evidence.attribute_fates,
        vec![("uv".to_owned(), AttributeFate::Preserved)]
    );
    assert_matches_field(&outcome.mesh, &outcome.mesh.attributes[0]);
}

/// A tool channel with the same name but a different width is not the same
/// channel: tool faces stay unmapped rather than being read wrongly.
#[test]
fn a_mismatched_tool_channel_is_not_mixed_in() {
    let a = with_field(subject(), "uv");
    let mut b = tool();
    let n = b.positions.len();
    b.attributes
        .push(AttributeChannel::new("uv", vec![9.0; n], 1, Blend::Linear));
    let outcome = run(&a, &b, BooleanOperator::Difference);
    outcome.mesh.validate_structure().expect("valid");
    let channel = &outcome.mesh.attributes[0];
    assert_eq!(channel.width, 2);
    assert!(channel
        .corner_indices
        .as_ref()
        .expect("corner-indexed")
        .contains(&AttributeChannel::UNMAPPED));
    assert_matches_field(&outcome.mesh, channel);
}
