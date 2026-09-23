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

/// Tag every corner with its own triangle and barycentric identity:
/// `[triangle + base, e0, e1, e2]`. Blended at a result corner this reads
/// `[source, w0, w1, w2]` -- the weights the sampler actually used.
fn tagged(mut m: TriMesh, base: f64) -> TriMesh {
    let n = m.indices.len();
    let mut v = Vec::with_capacity(n * 4);
    for c in 0..n {
        let mut e = [0.0; 3];
        e[c % 3] = 1.0;
        v.extend_from_slice(&[base + (c / 3) as f64, e[0], e[1], e[2]]);
    }
    let ix = (0..n as u32).collect();
    m.attributes.push(AttributeChannel::corner_indexed(
        "tag",
        v,
        4,
        Blend::Linear,
        ix,
    ));
    m
}

/// Every mapped corner's weights are those of a point INSIDE its source.
///
/// The linear-field oracle cannot see a wrong source triangle: a linear
/// field extrapolates to the right value from any coplanar triangle. This
/// catches it: a corner sampled outside its triangle has a negative weight,
/// and piecewise data (texture atlases, per-face ids) then reads a value
/// belonging to a different region.
fn assert_inside(mesh: &TriMesh, what: &str) -> usize {
    let ch = mesh
        .attributes
        .iter()
        .find(|c| c.name == "tag")
        .expect("tag");
    let mut checked = 0;
    for c in 0..mesh.indices.len() {
        let Some(t) = ch.at_corner(&mesh.indices, c) else {
            continue;
        };
        let worst = t[1].min(t[2]).min(t[3]);
        assert!(
            worst >= -1e-9,
            "{what}: corner {c} sampled outside its source (weight {worst})"
        );
        checked += 1;
    }
    checked
}

/// A rotated tool: the cut crosses the subject's face diagonals, so a
/// simplified result face spans BOTH triangles of a face. The regression
/// #116's first cut missed (weights to -0.65 on this very layout).
#[test]
fn a_corner_is_sampled_inside_its_true_source_triangle() {
    let s = tagged(subject(), 0.0);
    let t = tagged(boxx(0.5, 0.3, 0.5, 1.5, 1.5, 3.0, 0.4), 1000.0);
    for op in [
        BooleanOperator::Difference,
        BooleanOperator::Union,
        BooleanOperator::Intersection,
    ] {
        let out = run(&s, &t, op);
        out.mesh.validate_structure().expect("valid result");
        let n = assert_inside(&out.mesh, &format!("{op:?}"));
        assert!(n > 0, "{op:?}: nothing was mapped, so nothing was checked");
    }
}

fn cutters() -> Vec<TriMesh> {
    [-0.6, 0.0, 0.6]
        .iter()
        .map(|&x| tagged(boxx(x, 0.0, 0.5, 0.3, 3.0, 1.0, 0.0), 1000.0))
        .collect()
}

/// `subtract_many` (grouped: disjoint cutters fused into one tool) carries
/// channels, samples inside, and matches the one-at-a-time result.
#[test]
fn subtract_many_carries_channels_like_the_sequential_path() {
    let s = tagged(subject(), 0.0);
    let tools = cutters();
    let batch = BoolmeshBoolean::new()
        .subtract_many(&s, &tools, &options())
        .expect("subtract_many");
    batch.mesh.validate_structure().expect("valid");
    assert!(assert_inside(&batch.mesh, "batch") > 0);
    assert_eq!(
        batch.evidence.attribute_fates,
        vec![("tag".to_owned(), AttributeFate::Interpolated)]
    );
    // Same data as subtracting one at a time. Compared on a CONTINUOUS
    // field: the two paths simplify differently, so their faces differ,
    // but every corner of both must read the field exactly.
    let s2 = with_field(subject(), "uv");
    let tools2: Vec<TriMesh> = tools.iter().map(|t| with_field(t.clone(), "uv")).collect();
    let batch2 = BoolmeshBoolean::new()
        .subtract_many(&s2, &tools2, &options())
        .expect("batch");
    let mut seq = s2.clone();
    for t in &tools2 {
        seq = run(&seq, t, BooleanOperator::Difference).mesh;
    }
    let uv = |m: &TriMesh| {
        m.attributes
            .iter()
            .find(|c| c.name == "uv")
            .cloned()
            .expect("uv")
    };
    assert_matches_field(&batch2.mesh, &uv(&batch2.mesh));
    assert_matches_field(&seq, &uv(&seq));
}

/// Order-independent fingerprint of which source each mapped corner took:
/// the multiset of `(source id, corner position)`, rounded. Two results
/// with the same fingerprint put the same source data at the same places.
fn value_sum(mesh: &TriMesh) -> Vec<(i64, [i64; 3])> {
    let ch = mesh
        .attributes
        .iter()
        .find(|c| c.name == "tag")
        .expect("tag");
    let r = |x: f64| (x * 1e6).round() as i64;
    let mut out: Vec<_> = (0..mesh.indices.len())
        .filter_map(|c| {
            let t = ch.at_corner(&mesh.indices, c)?;
            let p = mesh.positions[mesh.indices[c] as usize];
            Some((r(t[0]), [r(p.x), r(p.y), r(p.z)]))
        })
        .collect();
    out.sort_unstable();
    out
}

/// The analytic box path carries channels: plane lookup recovers sources.
#[test]
fn the_analytic_box_path_carries_channels() {
    let s = tagged(subject(), 0.0);
    let tools = cutters();
    let out = BoolmeshBoolean::new()
        .subtract_boxes_analytic(&s, &tools, &options(), 1_000_000)
        .expect("ok")
        .expect("boxes are recognised");
    assert!(out.evidence.analytic_path);
    out.mesh.validate_structure().expect("valid");
    let ch = out
        .mesh
        .attributes
        .iter()
        .find(|c| c.name == "tag")
        .expect("carried");
    let unmapped = ch.corner_indices.as_ref().map_or(0, |c| {
        c.iter()
            .filter(|&&i| i == AttributeChannel::UNMAPPED)
            .count()
    });
    assert_eq!(unmapped, 0, "every face lies in a host or cutter face");
    assert!(assert_inside(&out.mesh, "analytic") > 0);
    // The general path puts the same source data at the same places.
    let general = BoolmeshBoolean::new()
        .subtract_many(&s, &tools, &options())
        .expect("general");
    let mut a: Vec<_> = value_sum(&out.mesh).into_iter().map(|(s, _)| s).collect();
    let mut g: Vec<_> = value_sum(&general.mesh)
        .into_iter()
        .map(|(s, _)| s)
        .collect();
    a.dedup();
    g.dedup();
    assert_eq!(a, g, "the same source triangles contribute");
}

/// `union_many` (balanced tree) carries channels; a solid WITHOUT the
/// channel contributes unmapped faces, not zeros, and not a lost channel.
#[test]
fn union_many_carries_channels_and_unmaps_a_solid_without_them() {
    let mut row: Vec<TriMesh> = (0..5)
        .map(|i| {
            tagged(
                boxx(i as f64 * 0.7, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0),
                1000.0 * i as f64,
            )
        })
        .collect();
    row[3].attributes.clear();
    let out = BoolmeshBoolean::new()
        .union_many(&row, &options())
        .expect("union");
    out.mesh.validate_structure().expect("valid");
    assert_eq!(
        out.evidence.attribute_fates,
        vec![("tag".to_owned(), AttributeFate::Interpolated)]
    );
    assert!(assert_inside(&out.mesh, "union_many") > 0);
    let ch = &out.mesh.attributes[0];
    let unmapped = ch.corner_indices.as_ref().expect("corner").iter();
    let unmapped = unmapped
        .filter(|&&i| i == AttributeChannel::UNMAPPED)
        .count();
    assert!(
        unmapped > 0,
        "solid 3 had no channel: its faces are unmapped"
    );
    let sources: Vec<i64> = value_sum(&out.mesh)
        .into_iter()
        .map(|(s, _)| (s / 1_000_000) / 1000)
        .collect();
    assert!(!sources.contains(&3), "no value was invented for solid 3");
}

/// Symmetric difference composes union then difference: the fate reflects
/// both, not just the last step.
#[test]
fn symmetric_difference_reports_the_composed_fate() {
    let s = tagged(subject(), 0.0);
    let t = tagged(tool(), 1000.0);
    let out = run(&s, &t, BooleanOperator::SymmetricDifference);
    out.mesh.validate_structure().expect("valid");
    assert_eq!(out.evidence.sub_operations, 3);
    assert_eq!(
        out.evidence.attribute_fates,
        vec![("tag".to_owned(), AttributeFate::Interpolated)]
    );
    assert!(assert_inside(&out.mesh, "symdiff") > 0);
}

/// An empty result keeps the channel, empty, and says nothing was lost.
#[test]
fn an_empty_result_reports_the_channel_preserved() {
    let s = tagged(subject(), 0.0);
    let far = tagged(boxx(50.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0), 1000.0);
    let out = run(&s, &far, BooleanOperator::Intersection);
    assert!(out.mesh.indices.is_empty());
    out.mesh.validate_structure().expect("valid");
    assert_eq!(
        out.evidence.attribute_fates,
        vec![("tag".to_owned(), AttributeFate::Preserved)]
    );
    assert_eq!(out.mesh.attributes.len(), 1);
}

/// A batch whose LAST step touches nothing: its own fate is `Preserved`,
/// but the result carries values the first step derived. The batch must
/// report the composition, not the last step.
#[test]
fn a_batch_reports_every_step_not_just_the_last() {
    let s = tagged(subject(), 0.0);
    // Overlapping cutters so they land in separate groups (two steps); the
    // second is fully inside the first, so its step cuts nothing new.
    let tools = vec![
        tagged(boxx(0.0, 0.0, 0.5, 0.8, 3.0, 1.0, 0.0), 1000.0),
        tagged(boxx(0.0, 0.0, 0.7, 0.4, 2.8, 0.6, 0.0), 2000.0),
    ];
    let out = BoolmeshBoolean::new()
        .subtract_many(&s, &tools, &options())
        .expect("batch");
    assert_eq!(out.evidence.sub_operations, 2, "two groups, two steps");
    // Precondition the test relies on: the last step alone is a no-op.
    let first = run(&s, &tools[0], BooleanOperator::Difference);
    let last = run(&first.mesh, &tools[1], BooleanOperator::Difference);
    assert_eq!(
        last.evidence.attribute_fates,
        vec![("tag".to_owned(), AttributeFate::Preserved)]
    );
    assert_eq!(
        out.evidence.attribute_fates,
        vec![("tag".to_owned(), AttributeFate::Interpolated)]
    );
}

/// A channel-less solid that the tree places as a LEFT operand: pairing
/// reads fates from the left, so without conforming it the channel would
/// vanish from that subtree and the batch would call it dropped.
#[test]
fn union_many_keeps_a_channel_a_left_operand_lacks() {
    let mut row: Vec<TriMesh> = (0..4)
        .map(|i| {
            tagged(
                boxx(i as f64 * 0.7, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0),
                1000.0 * i as f64,
            )
        })
        .collect();
    // Level 0 pairs (0,1) and (2,3): solid 2 is a left operand.
    row[2].attributes.clear();
    let out = BoolmeshBoolean::new()
        .union_many(&row, &options())
        .expect("union");
    out.mesh.validate_structure().expect("valid");
    assert_eq!(
        out.evidence.attribute_fates,
        vec![("tag".to_owned(), AttributeFate::Interpolated)]
    );
    assert!(out.mesh.attributes.iter().any(|c| c.name == "tag"));
}

/// Seam side: every subject face is split by its diagonal into two
/// triangles tagged apart (a UV atlas seam along the diagonal). A cut
/// vertex landing ON that diagonal is shared by result faces on both
/// sides. Each face must read ITS side's triangle: located by a point just
/// inside the face, not by the corner (which is in both). Nearest-blend so
/// a wrong side reads the other triangle's id outright.
#[test]
fn a_corner_on_a_seam_reads_the_side_its_face_is_on() {
    let mut s = subject();
    let tris = s.triangle_count();
    let ids: Vec<f64> = (0..tris).flat_map(|t| [t as f64; 3]).collect();
    let ix = (0..s.indices.len() as u32).collect();
    s.attributes.push(AttributeChannel::corner_indexed(
        "face",
        ids,
        1,
        Blend::Nearest,
        ix,
    ));
    // Rotated so the cut crosses every face's diagonal.
    let t = boxx(0.5, 0.3, 0.5, 1.5, 1.5, 3.0, 0.4);
    let out = run(&s, &t, BooleanOperator::Difference);
    out.mesh.validate_structure().expect("valid");
    let ch = &out.mesh.attributes[0];
    let mut checked = 0;
    for f in 0..out.mesh.triangle_count() {
        let q = |k: usize| out.mesh.positions[out.mesh.indices[f * 3 + k] as usize];
        let c = (q(0) + q(1) + q(2)) / 3.0;
        // Where this face is at corner k: a point just inside it from there.
        // The id read must name a subject triangle containing that point --
        // the corner's own side of any seam through it.
        for k in 0..3 {
            let Some(v) = ch.at_corner(&out.mesh.indices, f * 3 + k) else {
                continue;
            };
            let near = q(k) + (c - q(k)) * 1e-4;
            let sides: Vec<usize> = (0..tris).filter(|&t| closed_in(&s, t, near)).collect();
            assert!(
                sides.contains(&(v[0] as usize)),
                "face {f} corner {k}: read {} not {sides:?}",
                v[0]
            );
            checked += 1;
        }
    }
    assert!(checked > 0, "no subject face was checked");
}

/// Whether `p` lies in subject triangle `t`, edges included.
fn closed_in(m: &TriMesh, t: usize, p: Point3) -> bool {
    let q = |k: usize| m.positions[m.indices[t * 3 + k] as usize];
    let n = (q(1) - q(0)).cross(q(2) - q(0));
    (p - q(0)).dot(n).abs() <= 1e-9 * n.length()
        && (0..3).all(|k| (q((k + 1) % 3) - q(k)).cross(p - q(k)).dot(n) >= -1e-12)
}
