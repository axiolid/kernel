//! Outcome checks for `MeshPlaneSectionRegistry::section_request` (ADR 0091).
//!
//! A provider result is checked against the documented output invariants of
//! `axiolid_mesh_section_contract::SectionChain` and `InPlaneRegion`: the
//! frame, cardinalities, finiteness, adjacent duplicates, canonical ring
//! starts and orders, the budgets and the evidence counts. Ring orientation
//! is not checked: a floating-point area sign is unreliable far from the
//! origin (#274); the conformance suite checks it on fixtures.
//!
//! Cut regions (ADR 0092) are checked when, and only when, the request asked
//! for them: a provider that answers a request that did not ask violates the
//! contract. One that ignores a request (`cut_regions` is `None`) predates the
//! feature; the registry (`section.rs`) treats it as `Unsupported` and tries
//! the next provider, so this check never sees that case. Their rings follow the in-plane region checks,
//! and a cut request is refused for a sheet before any provider runs.
//!
//! The source topology depends on the request's role: a solid is closed and
//! oriented, a sheet may have boundary edges but not degenerate, non-manifold
//! or inconsistently wound ones. The evidence repeats the role and the
//! audited number of boundary edges, and a mesh without boundary edges (a
//! solid, or a closed sheet) has no chain end on the mesh boundary.

use std::cmp::Ordering;

use axiolid_contracts::{BackendId, ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::Point2;
use axiolid_mesh::{try_audit_mesh, TriMesh};
use axiolid_mesh_section_contract::{
    ChainEnd, CutRegion, InPlaneRegion, InPlaneSide, PlaneSection, PlaneSectionEvidence,
    SectionChain, SectionRequest, SectionRole,
};

use crate::section::validate_source_topology;

/// Check the source against the request's role; return its number of boundary
/// edges (zero for a solid, which has none).
///
/// A sheet follows `hidden_lines.rs::validate_topology`: boundary edges are
/// allowed, every other defect is refused.
pub(crate) fn validate_request_topology(
    mesh: &TriMesh,
    role: SectionRole,
    options: &ExecutionOptions,
) -> GeomResult<usize> {
    match role {
        SectionRole::Solid => {
            validate_source_topology(mesh, options)?;
            Ok(0)
        }
        SectionRole::Sheet => {
            let health = try_audit_mesh(mesh, options.tolerance())
                .map_err(|_| GeomError::BudgetExceeded { resource: "memory" })?;
            if let Some((triangle, index)) = health.first_invalid_index {
                return Err(GeomError::InvalidInput(format!(
                    "triangle {triangle} references position {index}, but the mesh has {}",
                    health.positions
                )));
            }
            if let Some(position) = health.first_non_finite_position {
                return Err(GeomError::InvalidInput(format!(
                    "position {position} is not finite"
                )));
            }
            if health.degenerate_triangles > 0 {
                return Err(GeomError::Degenerate(format!(
                    "{} degenerate triangle(s)",
                    health.degenerate_triangles
                )));
            }
            if health.non_manifold_edges > 0 || health.inconsistent_winding_edges > 0 {
                return Err(GeomError::NotManifold(format!(
                    "a sheet must be a consistently wound two-manifold; non_manifold={}, \
                     inconsistent_winding={}",
                    health.non_manifold_edges, health.inconsistent_winding_edges
                )));
            }
            Ok(health.boundary_edges)
        }
        _ => Err(GeomError::Unsupported {
            backend: BackendId::new("mesh-section-registry"),
            operation: Operation::MeshPlaneSection,
        }),
    }
}

/// Lexicographic order of plane-local points: `x`, then `y`, by `total_cmp`.
fn point_cmp(left: &Point2, right: &Point2) -> Ordering {
    left.x
        .total_cmp(&right.x)
        .then_with(|| left.y.total_cmp(&right.y))
}

/// Lexicographic order of point sequences; a prefix sorts first.
fn sequence_cmp(left: &[Point2], right: &[Point2]) -> Ordering {
    left.iter()
        .zip(right)
        .map(|(a, b)| point_cmp(a, b))
        .find(|order| order.is_ne())
        .unwrap_or_else(|| left.len().cmp(&right.len()))
}

/// Whether `points` has an adjacent duplicate; a ring also compares its last
/// point with its first.
fn has_adjacent_duplicate(points: &[Point2], ring: bool) -> bool {
    let pairs = if ring {
        points.len()
    } else {
        points.len().saturating_sub(1)
    };
    (0..pairs).any(|index| points[index] == points[(index + 1) % points.len()])
}

/// Whether a ring starts at its lexicographically smallest point.
fn starts_at_smallest(points: &[Point2]) -> bool {
    points
        .iter()
        .all(|point| point_cmp(&points[0], point).is_le())
}

/// Order of chains: point sequence, a closed chain before an open one.
fn chain_cmp(left: &SectionChain, right: &SectionChain) -> Ordering {
    sequence_cmp(&left.points, &right.points).then_with(|| right.is_closed().cmp(&left.is_closed()))
}

fn violation(backend: BackendId, detail: &str) -> GeomError {
    GeomError::BackendContractViolation {
        backend,
        detail: detail.into(),
    }
}

/// Refuse a cut-region request for a sheet: a sheet has no inside.
pub(crate) fn validate_cut_request(request: &SectionRequest) -> GeomResult<()> {
    if request.cut_regions.is_some() && request.role == SectionRole::Sheet {
        return Err(GeomError::InvalidInput(
            "cut regions need SectionRole::Solid: a sheet has no inside".into(),
        ));
    }
    Ok(())
}

/// Check one region ring and add its points to the running vertex total.
///
/// `label` names the kind of region in the messages and `resource` the
/// vertex budget that runs out.
fn check_ring(
    backend: BackendId,
    label: &str,
    resource: &'static str,
    max_vertices: usize,
    ring: &[Point2],
    vertices: &mut usize,
) -> GeomResult<()> {
    if ring.len() < 3 {
        return Err(violation(
            backend,
            &format!("{label} ring needs at least three points"),
        ));
    }
    if !ring.iter().all(|point| point.is_finite()) {
        return Err(violation(
            backend,
            &format!("{label} contains a non-finite point"),
        ));
    }
    if has_adjacent_duplicate(ring, true) {
        return Err(violation(
            backend,
            &format!("{label} ring contains adjacent duplicate points"),
        ));
    }
    if !starts_at_smallest(ring) {
        return Err(violation(
            backend,
            &format!("{label} ring does not start at its smallest point"),
        ));
    }
    *vertices = vertices
        .checked_add(ring.len())
        .ok_or(GeomError::BudgetExceeded { resource })?;
    if *vertices > max_vertices {
        return Err(GeomError::BudgetExceeded { resource });
    }
    Ok(())
}

/// Check the cut regions of a section against the request that asked for them
/// (or did not).
fn check_cut_regions(
    backend: BackendId,
    request: &SectionRequest,
    section: &PlaneSection,
) -> GeomResult<()> {
    match (request.cut_regions, &section.cut_regions) {
        (None, None) => Ok(()),
        // The registry falls back past a provider without cut regions.
        (Some(_), None) => Err(violation(
            backend,
            "provider did not trace the requested cut regions",
        )),
        (None, Some(_)) => Err(violation(
            backend,
            "provider returned cut regions that were not requested",
        )),
        (Some(limits), Some(cut_regions)) => {
            if section
                .regions
                .iter()
                .any(|region| region.side == InPlaneSide::Both)
            {
                return Err(violation(
                    backend,
                    "cut regions were returned beside in-plane faces of both sides",
                ));
            }
            if cut_regions.len() > limits.max_regions {
                return Err(GeomError::BudgetExceeded {
                    resource: "section cut regions",
                });
            }
            let mut vertices = 0usize;
            let mut previous: Option<&CutRegion> = None;
            for region in cut_regions {
                for ring in std::iter::once(&region.outer).chain(&region.holes) {
                    check_ring(
                        backend,
                        "cut region",
                        "section cut region vertices",
                        limits.max_region_vertices,
                        ring,
                        &mut vertices,
                    )?;
                }
                if region
                    .holes
                    .windows(2)
                    .any(|pair| sequence_cmp(&pair[0], &pair[1]).is_gt())
                {
                    return Err(violation(
                        backend,
                        "cut region holes are not in canonical order",
                    ));
                }
                if previous.is_some_and(|last| sequence_cmp(&last.outer, &region.outer).is_gt()) {
                    return Err(violation(backend, "cut regions are not in canonical order"));
                }
                previous = Some(region);
            }
            if section.evidence.cut_pairs > limits.max_pairs {
                return Err(GeomError::BudgetExceeded {
                    resource: "section cut pairs",
                });
            }
            Ok(())
        }
    }
}

pub(crate) fn validate_plane_section(
    backend: BackendId,
    mesh: &TriMesh,
    request: &SectionRequest,
    boundary_edges: usize,
    section: &PlaneSection,
) -> GeomResult<()> {
    let violation = |detail: &str| violation(backend, detail);
    let limits = request.limits;
    let region_limits = request.region_limits;
    if section.frame != request.frame {
        return Err(violation("provider changed the requested section frame"));
    }

    if section.chains.len() > limits.max_contours {
        return Err(GeomError::BudgetExceeded {
            resource: "section contours",
        });
    }
    let mut chain_vertices = 0usize;
    let mut previous: Option<&SectionChain> = None;
    for chain in &section.chains {
        let minimum = if chain.is_closed() { 3 } else { 2 };
        if chain.points.len() < minimum {
            return Err(violation(
                "section chain has too few points: three if closed, two if open",
            ));
        }
        if !chain.points.iter().all(|point| point.is_finite()) {
            return Err(violation("section chain contains a non-finite point"));
        }
        if has_adjacent_duplicate(&chain.points, chain.is_closed()) {
            return Err(violation(
                "section chain contains adjacent duplicate points",
            ));
        }
        let canonical = if chain.is_closed() {
            starts_at_smallest(&chain.points)
        } else {
            chain
                .points
                .first()
                .zip(chain.points.last())
                .is_some_and(|(first, last)| point_cmp(first, last).is_lt())
        };
        if boundary_edges == 0
            && chain
                .ends
                .is_some_and(|ends| ends.contains(&ChainEnd::MeshBoundary))
        {
            return Err(violation(
                "a section chain ends on the mesh boundary of a mesh without one",
            ));
        }
        if !canonical {
            return Err(violation(
                "section chain is not in canonical start or direction",
            ));
        }
        if previous.is_some_and(|last| chain_cmp(last, chain).is_gt()) {
            return Err(violation("section chains are not in canonical order"));
        }
        previous = Some(chain);
        chain_vertices =
            chain_vertices
                .checked_add(chain.points.len())
                .ok_or(GeomError::BudgetExceeded {
                    resource: "section output vertices",
                })?;
        if chain_vertices > limits.max_output_vertices {
            return Err(GeomError::BudgetExceeded {
                resource: "section output vertices",
            });
        }
    }

    if section.regions.len() > region_limits.max_regions {
        return Err(GeomError::BudgetExceeded {
            resource: "section in-plane regions",
        });
    }
    let mut region_vertices = 0usize;
    let mut previous: Option<&InPlaneRegion> = None;
    for region in &section.regions {
        for ring in std::iter::once(&region.outer).chain(&region.holes) {
            check_ring(
                backend,
                "in-plane region",
                "section in-plane region vertices",
                region_limits.max_region_vertices,
                ring,
                &mut region_vertices,
            )?;
        }
        if region
            .holes
            .windows(2)
            .any(|pair| sequence_cmp(&pair[0], &pair[1]).is_gt())
        {
            return Err(violation(
                "in-plane region holes are not in canonical order",
            ));
        }
        if previous.is_some_and(|last| {
            last.side
                .cmp(&region.side)
                .then_with(|| sequence_cmp(&last.outer, &region.outer))
                .is_gt()
        }) {
            return Err(violation("in-plane regions are not in canonical order"));
        }
        previous = Some(region);
    }

    check_cut_regions(backend, request, section)?;

    let evidence = section.evidence;
    if evidence.in_plane_pairs > region_limits.max_in_plane_pairs {
        return Err(GeomError::BudgetExceeded {
            resource: "section in-plane pairs",
        });
    }
    let derived = PlaneSectionEvidence::input_mesh(mesh.triangle_count())
        .with_source_counts(
            evidence.in_plane_triangles,
            evidence.touching_edges,
            evidence.in_plane_pairs,
        )
        .with_role(request.role, boundary_edges)
        .with_output(&section.chains, &section.regions)
        .with_cut_pairs(if request.cut_regions.is_some() {
            evidence.cut_pairs
        } else {
            0
        })
        .with_cut_output(section.cut_regions.as_deref().unwrap_or_default());
    if evidence != derived
        || !evidence.is_derived_from_input_mesh()
        || evidence.in_plane_triangles > mesh.triangle_count()
    {
        return Err(violation(
            "section evidence does not match input/output counts",
        ));
    }
    Ok(())
}
