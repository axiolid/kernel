//! The exact referee for a boolean the mesh cut cannot settle (#391, #423).
//!
//! The mesh boolean has no tolerance (ADR 0047). Operands whose faces meet
//! in exact arithmetic but were placed by linear parts that differ in a bit
//! (a turn computed two ways) come out of it as two bodies, a pinch, or one
//! body with sub-rounding wedge cracks (ADR 0096, D6). The exact compiler
//! reads such faces as one surface within the tolerance and fuses them
//! (ADR 0080), so it decides: [`boolean`] compiles the boolean node exactly
//! and tessellates the exact result ([`crate::brep::tessellate_exact`]).
//! A difference or an intersection is refereed like a union (#423, which
//! lifts the union-only scope of ADR 0096, D7): a flush cutter a residue off
//! its host's face leaves the void it cut out, skins a rounding thick, as a
//! shell of its own.
//!
//! The result is accepted only as a closed solid, or an empty one, that
//! meets the same guards as a mesh union: no pinch (#194) and no two bodies
//! within the tolerance (ADR 0096, D2). Anything else is a refusal with its
//! reason, and the
//! caller keeps what the mesh cut gave (`ReferenceMeshCompiler::build_boolean`):
//! the referee is never worse than the mesh path.
//!
//! The mesh is the tessellation of the exact union. Its deviation is the
//! tessellation's certified bound under [`REFEREED_UNION_EXACT`], or under
//! [`REFEREED_UNION`] (the difference's and the intersection's details name
//! their operator) when the exact compiler read a feature within the
//! tolerance, with the move named as [`COINCIDENT_OPERANDS`]. It carries no
//! channel: the operands' are `Dropped(ProviderLimitation)`, as no triangle
//! of the result descends from a triangle of an operand.
//!
//! Mesh operands are not refereed: the exact compiler refuses them, and
//! flattening them is reserved for the deviation report (ADR 0080, #266).

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Scalar};
use axiolid_mesh::{AttributeFate, DropReason, TriMesh};
use axiolid_mesh_compile_contract::MeshClosure;
use axiolid_model::{GeometryGraph, NodeId};

use crate::channels::{self, Built};
use crate::deviation::{
    with_readings, REFEREED_DIFFERENCE, REFEREED_DIFFERENCE_EXACT, REFEREED_INTERSECTION,
    REFEREED_INTERSECTION_EXACT, REFEREED_UNION, REFEREED_UNION_EXACT,
};
use crate::ReferenceExactCompiler;

/// What a difference or an intersection must show before the exact result
/// replaces its mesh (#423): the mesh holds a shell the exact result lacks,
/// and the exact compiler read the operands' faces no further than `reach`
/// (see [`within_reach`]).
pub(crate) struct Against<'a> {
    /// The mesh result of the cut.
    pub mesh: &'a TriMesh,
    /// How far a reading may move an operand's face and still be rounding.
    pub reach: Scalar,
}

/// Whether a reading of `linear` is within the rounding `reach`
/// ([`crate::snap::reach`]): the exact compiler then reads faces that agree
/// up to rounding, and not an authored thin feature (ADR 0096, D7).
pub(crate) fn within_reach(linear: Scalar, reach: Scalar) -> bool {
    linear <= reach
}

/// The boolean node `id` of `graph`, an `operator`, compiled exactly and
/// tessellated, or the reason it is refused. A difference or an
/// intersection is decided only [`Against`] its mesh.
pub(crate) fn boolean(
    graph: &GeometryGraph,
    id: NodeId,
    operator: BooleanOperator,
    options: &ExecutionOptions,
    (subject, tool): (&Built, Option<&Built>),
    against: Option<Against<'_>>,
    rounding: Option<Scalar>,
) -> Result<Built, String> {
    let (exact, report) = ReferenceExactCompiler::new()
        .compile_exact_with_report(graph, id, options)
        .map_err(|error| error.to_string())?;
    if let Some(against) = against {
        if !crate::deviation::extra_shells(&exact, against.mesh) {
            return Err("the mesh result has no more shells than the exact result".to_owned());
        }
        if !within_reach(report.linear(), against.reach) {
            return Err(
                "the exact compiler read the operands' faces further than their rounding: \
                 an authored thin feature stays"
                    .to_owned(),
            );
        }
    }
    let name = format!("{operator:?}").to_lowercase();
    // A skin the guard refused (#400, ADR 0102) is decided only where the
    // exact compiler read nothing above the rounding reach as one: a skin
    // above the reach but below the tolerance is authored (ADR 0096, D7).
    if let Some(reach) = rounding {
        if !within_reach(report.linear(), reach) {
            return Err(format!(
                "the exact {name} reads features {:e} apart as one, above the rounding \
                 reach {reach:e}: a skin thinner than the tolerance is authored (ADR 0096, D7)",
                report.linear()
            ));
        }
    }
    let eps = options.tolerance().linear();
    let (mesh, closure, tessellation) = crate::brep::tessellate_exact(
        &exact,
        options.tolerance(),
        crate::compiler::chord_error(options),
        // The result closes only within the reading of its booleans.
        report.linear().min(eps),
    )
    .map_err(|error| error.to_string())?;
    if closure != MeshClosure::Solid && !mesh.indices.is_empty() {
        return Err(format!(
            "the exact {name} does not tessellate to a closed solid"
        ));
    }
    if crate::pinch::find(&mesh).is_some() {
        return Err(format!("the tessellated exact {name} touches itself"));
    }
    if crate::contact::unmerged(&mesh, eps).is_some() {
        return Err(format!(
            "the tessellated exact {name} holds two bodies within the tolerance"
        ));
    }
    if let Some(reach) = rounding {
        if axiolid_inspect::thin_sheet(&mesh, reach).is_some() {
            return Err(format!(
                "the tessellated exact {name} holds a skin within rounding"
            ));
        }
    }
    let dropped = |built: &Built| {
        built
            .fates
            .clone()
            .into_vec()
            .into_iter()
            .map(|(name, _)| (name, AttributeFate::Dropped(DropReason::ProviderLimitation)))
    };
    let provider = dropped(subject)
        .chain(tool.into_iter().flat_map(dropped))
        .collect();
    let mut built =
        channels::after_boolean(mesh, &subject.fates, tool.map(|tool| &tool.fates), provider);
    let (within, exact_detail) = match operator {
        BooleanOperator::Union => (REFEREED_UNION, REFEREED_UNION_EXACT),
        BooleanOperator::Difference => (REFEREED_DIFFERENCE, REFEREED_DIFFERENCE_EXACT),
        BooleanOperator::Intersection => (REFEREED_INTERSECTION, REFEREED_INTERSECTION_EXACT),
        other => return Err(format!("no exact referee for the {other:?} operator")),
    };
    let detail = if report.is_exact() {
        exact_detail
    } else {
        within
    };
    let mut deviation = with_readings(detail, tessellation.worst_bound(), &report);
    deviation.carry_perturbations(&subject.deviation);
    if let Some(tool) = tool {
        deviation.carry_perturbations(&tool.deviation);
    }
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
        let built = boolean(
            &graph,
            root,
            BooleanOperator::Union,
            &options,
            (&subject, Some(&tool)),
            None,
            None,
        )
        .expect("fused");
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

    /// Two unit boxes a half apart and their intersection.
    fn overlapping_boxes() -> (GeometryGraph, NodeId) {
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
            transform: at(0.5),
        }));
        let root = push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: first,
            right: second,
            operator: BooleanOperator::Intersection,
        }));
        (builder.finish(vec![root]).expect("a valid graph"), root)
    }

    #[test]
    fn a_difference_or_intersection_is_decided_only_against_an_extra_shell() {
        let (graph, root) = overlapping_boxes();
        let options = ExecutionOptions::new(Tolerance::METRE);
        let operands = (
            &Built::leaf(TriMesh::new(Vec::new(), Vec::new())),
            Some(&Built::leaf(TriMesh::new(Vec::new(), Vec::new()))),
        );
        let tetra = textured();
        // One shell, as the exact intersection has: the mesh stands.
        let refused = boolean(
            &graph,
            root,
            BooleanOperator::Intersection,
            &options,
            operands,
            Some(Against {
                mesh: &tetra,
                reach: 1e-9,
            }),
            None,
        )
        .expect_err("no more shells than the exact result");
        assert!(refused.contains("no more shells"), "{refused}");
        // Two shells: the exact result decides, under the operator's detail.
        let mut two = tetra.clone();
        two.attributes.clear();
        let base = two.positions.len() as u32;
        two.positions.extend(
            tetra
                .positions
                .iter()
                .map(|p| *p + Vec3::new(5.0, 0.0, 0.0)),
        );
        two.indices.extend(tetra.indices.iter().map(|i| i + base));
        let built = boolean(
            &graph,
            root,
            BooleanOperator::Intersection,
            &options,
            operands,
            Some(Against {
                mesh: &two,
                reach: 1e-9,
            }),
            None,
        )
        .expect("decided");
        let report = format!("{:?}", built.deviation);
        assert!(
            report.contains(crate::deviation::REFEREED_INTERSECTION_EXACT),
            "{report}"
        );
        assert_eq!(axiolid_mesh::component_count(&built.mesh), 1);
    }
}
