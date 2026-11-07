//! Request checks shared by the registries of the orthographic view
//! operations (hidden lines and visible faces).
//!
//! Each takes the pieces of the request it needs, so both request types use
//! the same code and the same refusals. `label` names the operation in a
//! message ("hidden-line", "visible-face").

use std::collections::BTreeSet;

use axiolid_contracts::{ExecutionOptions, GeomError, ScratchRequirement};
use axiolid_core::Aabb2;
use axiolid_mesh::{audit_mesh_scratch_bytes, try_audit_mesh};
use axiolid_mesh_hidden_lines_contract::{
    FaceGroups, HiddenLineError, HiddenLineMesh, MeshRole, View, ViewProjection,
};

use crate::frame::validate_frame;

/// The limits both operations share, and the resource names their refusals
/// carry.
pub(crate) struct SharedLimits {
    pub label: &'static str,
    pub max_meshes: usize,
    pub max_source_vertices: usize,
    pub max_source_triangles: usize,
    pub meshes_resource: &'static str,
    pub vertices_resource: &'static str,
    pub triangles_resource: &'static str,
}

/// The fault of a `PerTriangle` group slice that does not have one id per
/// triangle.
pub(crate) fn group_length_fault(entry: &HiddenLineMesh<'_>, label: &str) -> Option<GeomError> {
    match entry.groups {
        FaceGroups::PerTriangle(groups) if groups.len() != entry.mesh.triangle_count() => {
            Some(GeomError::InvalidInput(format!(
                "{label} mesh has {} face group(s) for {} triangle(s)",
                groups.len(),
                entry.mesh.triangle_count()
            )))
        }
        _ => None,
    }
}

/// A `PerTriangle` group slice must have one id per triangle.
pub(crate) fn validate_group_lengths(
    meshes: &[HiddenLineMesh<'_>],
    label: &str,
) -> Result<(), HiddenLineError> {
    for entry in meshes {
        if let Some(error) = group_length_fault(entry, label) {
            return Err(HiddenLineError::mesh(entry.id, error));
        }
    }
    Ok(())
}

pub(crate) fn invalid(detail: impl Into<String>) -> HiddenLineError {
    GeomError::InvalidInput(detail.into()).into()
}

/// Total positions and triangles; `None` on overflow.
pub(crate) fn source_totals(meshes: &[HiddenLineMesh<'_>]) -> (Option<usize>, Option<usize>) {
    meshes
        .iter()
        .fold((Some(0usize), Some(0usize)), |(v, t), m| {
            (
                v.and_then(|v| v.checked_add(m.mesh.positions.len())),
                t.and_then(|t| t.checked_add(m.mesh.triangle_count())),
            )
        })
}

/// Mesh list, counts against the limits, unique ids and whole triangles.
pub(crate) fn validate_meshes(
    meshes: &[HiddenLineMesh<'_>],
    limits: &SharedLimits,
) -> Result<(), HiddenLineError> {
    let label = limits.label;
    if meshes.is_empty() {
        return Err(invalid(format!("{label} request has no meshes")));
    }
    if meshes.len() > limits.max_meshes {
        return Err(GeomError::BudgetExceeded {
            resource: limits.meshes_resource,
        }
        .into());
    }
    let (vertices, triangles) = source_totals(meshes);
    if vertices.is_none_or(|n| n > limits.max_source_vertices) {
        return Err(GeomError::BudgetExceeded {
            resource: limits.vertices_resource,
        }
        .into());
    }
    if triangles.is_none_or(|n| n > limits.max_source_triangles) {
        return Err(GeomError::BudgetExceeded {
            resource: limits.triangles_resource,
        }
        .into());
    }

    let mut ids = BTreeSet::new();
    for entry in meshes {
        if !ids.insert(entry.id) {
            return Err(HiddenLineError::mesh(
                entry.id,
                GeomError::InvalidInput("duplicate mesh id".into()),
            ));
        }
        if entry.mesh.indices.is_empty() || !entry.mesh.indices.len().is_multiple_of(3) {
            return Err(HiddenLineError::mesh(
                entry.id,
                GeomError::InvalidInput(format!(
                    "{label} mesh has no triangles or a partial triangle"
                )),
            ));
        }
    }
    Ok(())
}

/// Frame and depth range.
pub(crate) fn validate_view(
    view: &View,
    angular_tolerance: f64,
    label: &str,
) -> Result<(), HiddenLineError> {
    validate_frame(
        view.frame,
        angular_tolerance,
        &format!("{label} view frame"),
    )?;
    for (name, bound) in [("near", view.near), ("far", view.far)] {
        if bound.is_some_and(|value| !value.is_finite()) {
            return Err(invalid(format!("{label} {name} depth is not finite")));
        }
    }
    match view.projection {
        ViewProjection::Orthographic => {}
        ViewProjection::Perspective => {
            if !view.near.is_some_and(|near| near < 0.0) {
                return Err(invalid(format!(
                    "{label} perspective view needs a near depth below zero"
                )));
            }
        }
        _ => return Err(invalid(format!("{label} projection is not known"))),
    }
    if let (Some(near), Some(far)) = (view.near, view.far) {
        if far >= near {
            return Err(invalid(format!(
                "{label} far depth must be smaller than the near depth"
            )));
        }
    }
    Ok(())
}

/// The optional image crop.
pub(crate) fn validate_region(region: Option<Aabb2>, label: &str) -> Result<(), HiddenLineError> {
    if let Some(region) = region {
        let finite = region.min.is_finite() && region.max.is_finite();
        if !finite || region.min.x > region.max.x || region.min.y > region.max.y {
            return Err(invalid(format!(
                "{label} region must be finite with min not above max"
            )));
        }
    }
    Ok(())
}

/// Audit every mesh: indices, finiteness, degeneracy, manifoldness and the
/// closedness of solids. Allocating; call it after the memory budget admitted
/// a provider.
pub(crate) fn validate_topology(
    meshes: &[HiddenLineMesh<'_>],
    groups: Option<&str>,
    options: &ExecutionOptions,
) -> Result<(), HiddenLineError> {
    for entry in meshes {
        options.check_cancelled()?;
        if let Some(error) = mesh_verdict(entry, groups, options)? {
            return Err(HiddenLineError::mesh(entry.id, error));
        }
    }
    Ok(())
}

/// Why one mesh is refused, or `None` when it passes the audit. `groups`
/// names the operation when the mesh's face groups are read, so a short
/// `PerTriangle` slice is a fault of the mesh. Shared by the
/// registries' `validate_topology` and the hidden-line registry's `audit`, so
/// both give the same verdict per mesh.
pub(crate) fn mesh_verdict(
    entry: &HiddenLineMesh<'_>,
    groups: Option<&str>,
    options: &ExecutionOptions,
) -> Result<Option<GeomError>, HiddenLineError> {
    if let Some(fault) = groups.and_then(|label| group_length_fault(entry, label)) {
        return Ok(Some(fault));
    }
    let health = try_audit_mesh(entry.mesh, options.tolerance())
        .map_err(|_| GeomError::BudgetExceeded { resource: "memory" })?;
    if let Some((triangle, index)) = health.first_invalid_index {
        return Ok(Some(GeomError::InvalidInput(format!(
            "triangle {triangle} references position {index}, but the mesh has {}",
            health.positions
        ))));
    }
    if let Some(position) = health.first_non_finite_position {
        return Ok(Some(GeomError::InvalidInput(format!(
            "position {position} is not finite"
        ))));
    }
    if health.degenerate_triangles > 0 {
        return Ok(Some(GeomError::Degenerate(format!(
            "{} degenerate triangle(s)",
            health.degenerate_triangles
        ))));
    }
    if health.non_manifold_edges > 0 || health.inconsistent_winding_edges > 0 {
        return Ok(Some(GeomError::NotManifold(format!(
            "non_manifold={}, inconsistent_winding={}",
            health.non_manifold_edges, health.inconsistent_winding_edges
        ))));
    }
    if entry.role == MeshRole::Solid && health.boundary_edges > 0 {
        return Ok(Some(GeomError::NotManifold(format!(
            "a solid must be closed; boundary={}",
            health.boundary_edges
        ))));
    }
    Ok(None)
}

/// Whether a provider's scratch fits the memory budget next to the audit.
pub(crate) fn scratch_fits(
    provider: ScratchRequirement,
    meshes: &[HiddenLineMesh<'_>],
    options: &ExecutionOptions,
) -> bool {
    let Some(budget) = options.memory_budget_bytes() else {
        return true;
    };
    let (vertices, triangles) = source_totals(meshes);
    let (Some(vertices), Some(_)) = (vertices, triangles) else {
        return false;
    };
    // Meshes are audited one after another, so the audit peak is the largest.
    let mut audit_bytes = 0usize;
    for entry in meshes {
        let Some(bytes) = audit_mesh_scratch_bytes(entry.mesh.triangle_count()) else {
            return false;
        };
        audit_bytes = audit_bytes.max(bytes);
    }
    let elements = vertices.max(triangles.unwrap_or(0));
    let Some(provider_bytes) = provider.upper_bound_bytes(elements) else {
        return false;
    };
    // Audit and provider phases are sequential, so their peak is the maximum,
    // not their sum. Output storage is outside the scratch budget by contract.
    audit_bytes.max(provider_bytes) <= budget
}
