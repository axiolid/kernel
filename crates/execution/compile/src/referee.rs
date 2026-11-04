//! The exact referee for a union the mesh cut cannot settle (#391).
//!
//! The mesh boolean has no tolerance (ADR 0047). Operands whose faces meet
//! in exact arithmetic but were placed by linear parts that differ in a bit
//! (a turn computed two ways) come out of it as two bodies, a pinch, or one
//! body with sub-rounding wedge cracks (ADR 0096, D6). The exact compiler
//! reads such faces as one surface within the tolerance and fuses them
//! (ADR 0080), so it decides: [`union`] compiles the union node exactly and
//! tessellates the exact result ([`crate::brep::tessellate_exact`]).
//!
//! The result is accepted only as a closed solid that meets the same guards
//! as a mesh union: no pinch (#194) and no two bodies within the tolerance
//! (ADR 0096, D2). Anything else is a refusal with its reason, and the
//! caller keeps what the mesh cut gave (`ReferenceMeshCompiler::build_boolean`):
//! the referee is never worse than the mesh path.
//!
//! The mesh is the tessellation of the exact union. Its deviation is the
//! tessellation's certified bound under [`REFEREED_UNION_EXACT`], or under
//! [`REFEREED_UNION`] when the exact compiler read a feature within the
//! tolerance, with the move named as [`COINCIDENT_OPERANDS`]. It carries no
//! channel: the operands' are `Dropped(ProviderLimitation)`, as no triangle
//! of the result descends from a triangle of an operand.
//!
//! Mesh operands are not refereed: the exact compiler refuses them, and
//! flattening them is reserved for the deviation report (ADR 0080, #266).

use axiolid_contracts::ExecutionOptions;
use axiolid_mesh::{AttributeFate, DropReason};
use axiolid_mesh_compile_contract::MeshClosure;
use axiolid_model::{GeometryGraph, NodeId};

use crate::channels::{self, Built};
use crate::deviation::{with_readings, REFEREED_UNION, REFEREED_UNION_EXACT};
use crate::ReferenceExactCompiler;

/// The union node `id` of `graph`, compiled exactly and tessellated, or
/// the reason it is refused.
pub(crate) fn union(
    graph: &GeometryGraph,
    id: NodeId,
    options: &ExecutionOptions,
    subject: &Built,
    tool: &Built,
) -> Result<Built, String> {
    let (exact, report) = ReferenceExactCompiler::new()
        .compile_exact_with_report(graph, id, options)
        .map_err(|error| error.to_string())?;
    let eps = options.tolerance().linear();
    let (mesh, closure, tessellation) = crate::brep::tessellate_exact(
        &exact,
        options.tolerance(),
        crate::compiler::chord_error(options),
    )
    .map_err(|error| error.to_string())?;
    if closure != MeshClosure::Solid {
        return Err("the exact union does not tessellate to a closed solid".to_owned());
    }
    if crate::pinch::find(&mesh).is_some() {
        return Err("the tessellated exact union touches itself".to_owned());
    }
    if crate::contact::unmerged(&mesh, eps).is_some() {
        return Err("the tessellated exact union holds two bodies within the tolerance".to_owned());
    }
    let dropped = |built: &Built| {
        built
            .fates
            .clone()
            .into_vec()
            .into_iter()
            .map(|(name, _)| (name, AttributeFate::Dropped(DropReason::ProviderLimitation)))
    };
    let provider = dropped(subject).chain(dropped(tool)).collect();
    let mut built = channels::after_boolean(mesh, &subject.fates, Some(&tool.fates), provider);
    let detail = if report.is_exact() {
        REFEREED_UNION_EXACT
    } else {
        REFEREED_UNION
    };
    let mut deviation = with_readings(detail, tessellation.worst_bound(), &report);
    deviation.carry_perturbations(&subject.deviation);
    deviation.carry_perturbations(&tool.deviation);
    built.deviation = deviation;
    Ok(built)
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::{BooleanOperator, Point3, Tolerance, Transform3, Vec3};
    use axiolid_mesh::{AttributeChannel, Blend, TriMesh};
    use axiolid_model::{GeometryGraphBuilder, GeometryNode, Instance, SolidOperation};
    use axiolid_profile::{Profile, RectangleProfile};

    /// A tetrahedron with a corner-indexed `uv` channel: the stand-in for
    /// what an operand's history carried.
    fn textured() -> TriMesh {
        let mut mesh = TriMesh::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 0.0, 1.0),
            ],
            vec![0, 2, 1, 0, 1, 3, 1, 2, 3, 0, 3, 2],
        );
        mesh.attributes.push(AttributeChannel::corner_indexed(
            "uv",
            (0..24).map(f64::from).collect(),
            2,
            Blend::Linear,
            (0..12).collect(),
        ));
        mesh
    }

    #[test]
    fn an_operands_channel_is_dropped_by_name_and_reason() {
        let mut builder = GeometryGraphBuilder::new();
        let mut push = |node| builder.push(node).expect("a valid node");
        let profile = push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: 1.0,
                y: 1.0,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )));
        let solid = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 1.0,
        }));
        let at = |x| Transform3::from_translation(Vec3::new(x, 0.0, 0.0));
        let first = push(GeometryNode::Instance(Instance {
            source: solid,
            transform: at(0.0),
        }));
        let second = push(GeometryNode::Instance(Instance {
            source: solid,
            transform: at(1.0),
        }));
        let root = push(GeometryNode::SolidOperation(
            axiolid_model::SolidOperation::Boolean {
                left: first,
                right: second,
                operator: BooleanOperator::Union,
            },
        ));
        let graph = builder.finish(vec![root]).expect("a valid graph");
        let options = ExecutionOptions::new(Tolerance::METRE);
        // The operands as the mesh compiler would have built them: the
        // subject carried a channel.
        let subject = Built::leaf(textured());
        let tool = Built::leaf(TriMesh::new(Vec::new(), Vec::new()));
        let built = union(&graph, root, &options, &subject, &tool).expect("fused");
        assert!(built.mesh.attributes.is_empty());
        assert_eq!(
            built.fates.into_vec(),
            vec![(
                "uv".to_owned(),
                AttributeFate::Dropped(DropReason::ProviderLimitation)
            )]
        );
        assert!(built.local.is_none());
    }
}
