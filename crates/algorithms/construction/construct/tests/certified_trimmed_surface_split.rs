use std::collections::BTreeSet;

use axiolid_construct::{
    split_surface_pair_certified, CertifiedSurfacePairSplit3, CertifiedSurfacePairSplitOptions,
    CertifiedTrimmedSurfacePair3, SurfacePairMember, SurfacePairSplitUnresolvedReason,
};
use axiolid_core::Point3;
use axiolid_curve::KnotSpec;
use axiolid_nurbs::{
    CertifiedSurfaceSurfaceIntersection3, CertifiedSurfaceSurfaceIntersectionOptions,
};
use axiolid_surface::BSplineSurface;
use axiolid_topology::{audit_brep, Orientation};

fn xy_plane(x_min: f64, x_max: f64) -> BSplineSurface {
    BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        control_points: vec![
            vec![Point3::new(x_min, -1.0, 0.0), Point3::new(x_min, 1.0, 0.0)],
            vec![Point3::new(x_max, -1.0, 0.0), Point3::new(x_max, 1.0, 0.0)],
        ],
        u_knots: vec![10.0, 12.0],
        u_multiplicities: vec![2, 2],
        v_knots: vec![-4.0, 0.0],
        v_multiplicities: vec![2, 2],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    }
}

fn xz_plane(x_min: f64, x_max: f64) -> BSplineSurface {
    BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        control_points: vec![
            vec![Point3::new(x_min, 0.0, -1.0), Point3::new(x_min, 0.0, 1.0)],
            vec![Point3::new(x_max, 0.0, -1.0), Point3::new(x_max, 0.0, 1.0)],
        ],
        u_knots: vec![20.0, 24.0],
        u_multiplicities: vec![2, 2],
        v_knots: vec![3.0, 7.0],
        v_multiplicities: vec![2, 2],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    }
}

fn options() -> CertifiedSurfacePairSplitOptions {
    CertifiedSurfacePairSplitOptions::new(
        CertifiedSurfaceSurfaceIntersectionOptions::new(1.0e-7, 100_000, 100_000, 64)
            .expect("valid intersection policy"),
        1.0e-6,
    )
    .expect("valid split policy")
}

fn assert_all_generated_entities_are_reachable_and_oriented(split: &CertifiedTrimmedSurfacePair3) {
    let topology = split.brep.topology();
    let expected_edges: BTreeSet<_> = (0..topology.edges().len()).collect();
    let expected_vertices: BTreeSet<_> = (0..topology.vertices().len()).collect();
    let expected_loops: BTreeSet<_> = (0..topology.loops().len()).collect();
    let expected_faces: BTreeSet<_> = (0..topology.faces().len()).collect();

    let mut used_edges = BTreeSet::new();
    for loop_ in topology.loops() {
        for (index, edge_use) in loop_.edges.iter().enumerate() {
            used_edges.insert(edge_use.edge.index());
            let edge = &topology.edges()[edge_use.edge.index()];
            let head = match edge_use.orientation {
                Orientation::Forward => edge.end,
                Orientation::Reversed => edge.start,
            };
            let next_use = &loop_.edges[(index + 1) % loop_.edges.len()];
            let next_edge = &topology.edges()[next_use.edge.index()];
            let next_tail = match next_use.orientation {
                Orientation::Forward => next_edge.start,
                Orientation::Reversed => next_edge.end,
            };
            assert_eq!(
                head, next_tail,
                "loop edge uses must close in traversal order"
            );
        }
    }
    assert_eq!(
        used_edges, expected_edges,
        "every generated edge must bound a loop"
    );

    let mut used_vertices = BTreeSet::new();
    let mut endpoint_pairs = BTreeSet::new();
    for edge in topology.edges() {
        used_vertices.insert(edge.start.index());
        used_vertices.insert(edge.end.index());
        let pair = if edge.start < edge.end {
            (edge.start.index(), edge.end.index())
        } else {
            (edge.end.index(), edge.start.index())
        };
        assert!(
            endpoint_pairs.insert(pair),
            "duplicate topological edge endpoints"
        );
    }
    assert_eq!(
        used_vertices, expected_vertices,
        "every generated vertex must bound an edge"
    );

    let mut used_loops = BTreeSet::new();
    for face in topology.faces() {
        for bound in &face.bounds {
            assert!(
                used_loops.insert(bound.loop_id.index()),
                "generated loops must have exactly one owning face"
            );
        }
    }
    assert_eq!(
        used_loops, expected_loops,
        "every generated loop must be owned"
    );

    let result_faces = BTreeSet::from([
        split.split_faces[0].index(),
        split.split_faces[1].index(),
        split.unsplit_face.index(),
    ]);
    assert_eq!(
        result_faces, expected_faces,
        "the result must own every generated face exactly once"
    );

    assert_all_generated_geometry_is_referenced(split);
}

fn assert_all_generated_geometry_is_referenced(split: &CertifiedTrimmedSurfacePair3) {
    let topology = split.brep.topology();
    let used_curve3: BTreeSet<_> = topology
        .edges()
        .iter()
        .map(|edge| edge.curve.expect("strict edge carrier").index())
        .collect();
    assert_eq!(used_curve3, (0..split.brep.curves3().len()).collect());

    let used_surfaces: BTreeSet<_> = topology
        .faces()
        .iter()
        .map(|face| face.surface.expect("strict face support").index())
        .collect();
    assert_eq!(used_surfaces, (0..split.brep.surfaces().len()).collect());

    let mut used_pcurves = BTreeSet::new();
    for loop_ in topology.loops() {
        for edge_use in &loop_.edges {
            assert!(
                used_pcurves.insert(edge_use.pcurve.expect("strict pcurve").index()),
                "each face-local edge use must own a distinct pcurve"
            );
        }
    }
    assert!(used_pcurves.insert(split.embedded_curve.pcurve.index()));
    assert_eq!(used_pcurves, (0..split.brep.curves2().len()).collect());
}

#[test]
fn splits_boundary_owned_face_and_embeds_same_edge_in_containing_face() {
    let first = xy_plane(-1.0, 1.0);
    let second = xz_plane(-0.5, 0.5);
    let result = split_surface_pair_certified(&first, &second, options())
        .expect("certified topology integration succeeds");

    let split = match result {
        CertifiedSurfacePairSplit3::Split(split) => split,
        other => panic!("expected a complete trimmed split, got {other:?}"),
    };
    assert_eq!(split.visited_patch_pairs, 1);
    assert_eq!(split.boundary_queries, 8);
    assert!(split.max_surface_residual_upper_bound <= 1.0e-6);
    assert_eq!(split.brep.surfaces().len(), 2);
    assert_eq!(split.brep.topology().faces().len(), 3);
    assert_eq!(split.brep.topology().loops().len(), 3);
    assert_eq!(split.split_faces.len(), 2);
    assert_eq!(
        split.split_surface,
        axiolid_construct::SurfacePairMember::Second
    );
    assert_eq!(split.embedded_curve.face, split.unsplit_face);
    assert_eq!(split.embedded_curve.edge, split.intersection_edge);
    assert_eq!(split.embedded_curve.interval, axiolid_core::Interval::UNIT);
    assert!(split
        .brep
        .curves2()
        .get(split.embedded_curve.pcurve.index())
        .is_some());

    let health = audit_brep(split.brep.topology());
    assert!(
        health.is_tessellable(),
        "invalid split topology: {health:?}"
    );
    assert_eq!(health.dangling_references, 0);
    assert_eq!(health.open_loops, 0);
    assert_all_generated_entities_are_reachable_and_oriented(&split);

    let mut shared_uses = 0;
    let mut forward_uses = 0;
    let mut reverse_uses = 0;
    for (loop_index, loop_) in split.brep.topology().loops().iter().enumerate() {
        let loop_id = split
            .brep
            .topology()
            .loop_id_at(loop_index)
            .expect("enumerated loop exists");
        for (use_index, edge_use) in loop_.edges.iter().enumerate() {
            assert!(edge_use.pcurve.is_some());
            assert!(split.brep.pcurve_interval(loop_id, use_index).is_some());
            if edge_use.edge == split.intersection_edge {
                shared_uses += 1;
                match edge_use.orientation {
                    Orientation::Forward => forward_uses += 1,
                    Orientation::Reversed => reverse_uses += 1,
                }
            }
        }
    }
    assert_eq!(shared_uses, 2);
    assert_eq!(forward_uses, 1);
    assert_eq!(reverse_uses, 1);
    assert_eq!(
        split.brep.edge_interval(split.intersection_edge),
        Some(axiolid_core::Interval::UNIT)
    );

    for endpoint in [&split.trace.start, &split.trace.end] {
        for interval in [
            endpoint.parameters.first_u,
            endpoint.parameters.first_v,
            endpoint.parameters.second_u,
            endpoint.parameters.second_v,
        ] {
            assert!(interval.start.is_finite());
            assert!(interval.end.is_finite());
            assert!(interval.start <= interval.end);
        }
        assert!(endpoint.parameters.first_u.start >= 10.0);
        assert!(endpoint.parameters.first_u.end <= 12.0);
        assert!(endpoint.parameters.first_v.start >= -4.0);
        assert!(endpoint.parameters.first_v.end <= 0.0);
        assert!(endpoint.parameters.second_u.start >= 20.0);
        assert!(endpoint.parameters.second_u.end <= 24.0);
        assert!(endpoint.parameters.second_v.start >= 3.0);
        assert!(endpoint.parameters.second_v.end <= 7.0);
    }
}

#[test]
fn disjoint_surfaces_produce_complete_empty_split_without_brep() {
    let first = xy_plane(-1.0, 1.0);
    let mut second = xy_plane(-1.0, 1.0);
    for row in &mut second.control_points {
        for point in row {
            point.z = 2.0;
        }
    }
    let result = split_surface_pair_certified(&first, &second, options())
        .expect("disjoint certified query succeeds");
    assert!(matches!(
        result,
        CertifiedSurfacePairSplit3::Empty {
            visited_patch_pairs: 1,
            boundary_queries: 0,
        }
    ));
}

#[test]
fn swapped_inputs_preserve_geometric_owner_and_reverse_member_label() {
    let containing = xy_plane(-1.0, 1.0);
    let owner = xz_plane(-0.5, 0.5);
    let result = split_surface_pair_certified(&owner, &containing, options())
        .expect("swapped certified query succeeds");
    let CertifiedSurfacePairSplit3::Split(split) = result else {
        panic!("expected swapped split result");
    };
    assert_eq!(split.split_surface, SurfacePairMember::First);
    assert_eq!(split.brep.topology().faces().len(), 3);
    assert_eq!(split.embedded_curve.face, split.unsplit_face);
}

#[test]
fn certified_trace_above_residual_policy_stays_unresolved() {
    let first = xy_plane(-1.0, 1.0);
    let second = xz_plane(-0.5, 0.5);
    let options = CertifiedSurfacePairSplitOptions::new(
        CertifiedSurfaceSurfaceIntersectionOptions::default(),
        f64::MIN_POSITIVE,
    )
    .expect("positive finite policy");
    let result =
        split_surface_pair_certified(&first, &second, options).expect("valid query terminates");
    assert!(matches!(
        result,
        CertifiedSurfacePairSplit3::Unresolved {
            reason: SurfacePairSplitUnresolvedReason::ResidualExceedsPolicy,
            ..
        }
    ));
}

/// A chord that partitions BOTH patches now splits both.
///
/// Two planes crossing corner-to-corner cut each other's rectangle in
/// half. Before certified boundary roots existed this returned
/// `IntersectionUnresolved`: every endpoint lies exactly ON a patch
/// domain edge, which ordinary Krawczyk containment can never prove.
#[test]
fn a_chord_partitioning_both_patches_splits_both() {
    let first = xy_plane(-1.0, 1.0);
    let second = xz_plane(-1.0, 1.0);
    let result = split_surface_pair_certified(&first, &second, options())
        .expect("valid dual-boundary query terminates");

    let CertifiedSurfacePairSplit3::DualSplit(split) = result else {
        panic!("expected a dual split, got {result:?}");
    };

    // Four closed faces: two per input surface.
    assert_eq!(split.brep.topology().faces().len(), 4);
    let faces: BTreeSet<_> = split
        .first_faces
        .iter()
        .chain(split.second_faces.iter())
        .collect();
    assert_eq!(faces.len(), 4, "the four faces must be distinct");

    // One shared edge is what makes this a single arrangement rather than
    // two that merely coincide in space.
    for bound in split.brep.topology().faces().iter().flat_map(|f| &f.bounds) {
        let uses = &split.brep.topology().loops()[bound.loop_id.index()].edges;
        assert!(
            uses.iter().any(|u| u.edge == split.intersection_edge),
            "every trimmed loop must use the shared intersection edge"
        );
    }

    let health = audit_brep(split.brep.topology());
    assert!(health.is_tessellable(), "invalid dual split: {health:?}");
    assert_eq!(health.dangling_references, 0);
    assert_eq!(health.open_loops, 0);
}
#[test]
fn mixed_endpoint_ownership_refuses_open_trim_chains() {
    // Each patch sees the trace run from an interior point to one boundary,
    // so each is slit rather than partitioned. That is provable, not a gap
    // in this crate: the refusal must say so.
    let first = xy_plane(-1.0, 0.5);
    let second = xz_plane(-0.5, 1.0);
    let result = split_surface_pair_certified(&first, &second, options())
        .expect("valid mixed-ownership query terminates");
    match result {
        CertifiedSurfacePairSplit3::Unresolved {
            reason: SurfacePairSplitUnresolvedReason::NoPartitionExists,
            intersection: CertifiedSurfaceSurfaceIntersection3::Complete { traces, .. },
        } => assert_eq!(traces.len(), 1),
        other => panic!("expected mixed-ownership refusal, got {other:?}"),
    }
}

#[test]
fn split_policy_rejects_non_finite_or_non_positive_residual_limit() {
    let intersection = CertifiedSurfaceSurfaceIntersectionOptions::default();
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(CertifiedSurfacePairSplitOptions::new(intersection, invalid).is_err());
    }
}

/// A proven no-partition refusal does not change with the residual policy.
///
/// This is what makes the reason worth its own variant: the caller can
/// stop. If a tighter or looser policy could flip the verdict, it would be
/// an unimplemented case wearing a proof's clothing.
#[test]
fn no_partition_is_terminal_across_residual_policies() {
    let first = xy_plane(-1.0, 0.5);
    let second = xz_plane(-0.5, 1.0);
    for residual in [1e-12, 1e-9, 1e-6, 1e-3, 1.0] {
        let policy = CertifiedSurfacePairSplitOptions::new(
            CertifiedSurfaceSurfaceIntersectionOptions::default(),
            residual,
        )
        .expect("a positive finite residual limit is valid");
        let result = split_surface_pair_certified(&first, &second, policy)
            .expect("the query terminates at every policy");
        assert!(
            matches!(
                result,
                CertifiedSurfacePairSplit3::Unresolved {
                    reason: SurfacePairSplitUnresolvedReason::NoPartitionExists,
                    ..
                }
            ),
            "residual {residual} must not change a proven refusal, got {result:?}"
        );
    }
}
