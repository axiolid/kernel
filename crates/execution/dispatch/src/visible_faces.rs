//! Provider registration, validation, fallback, and budget policy for
//! visible face regions of meshes in an orthographic view.
//!
//! The shape, topology and memory checks are those of the hidden-line
//! registry (`view_checks`). The outcome check has no geometric overlap test:
//! ring points are rounded, so only the structure of the answer is checked.

use core::cmp::Ordering;
use std::sync::Arc;

use axiolid_contracts::{BackendId, ExecutionOptions, GeomError, Operation};
use axiolid_core::Point2;
use axiolid_mesh_hidden_lines_contract::{
    conformance::faces::{FaceConformanceReport, FaceConformanceSuite},
    FaceGroups, FaceSource, HiddenLineError, HiddenLineMesh, MeshRole, MeshVisibleFaces,
    VisibleFaceOutcome, VisibleFaceRequest,
};

use crate::device::matches_device;
use crate::view_checks::{
    scratch_fits, validate_meshes, validate_region, validate_topology, validate_view, SharedLimits,
};

const LABEL: &str = "visible-face";

#[derive(Debug, Clone)]
struct RegisteredVisibleFaces {
    priority: i32,
    provider: Arc<dyn MeshVisibleFaces>,
}

/// Ordered executable visible-face providers.
///
/// Only unsupported/unavailable providers permit fallback. Geometry,
/// degeneracy, budget, cancellation, and contract errors fail immediately.
/// A refusal caused by one input mesh names it
/// ([`HiddenLineError::mesh`]); a refusal of the whole request does not.
#[derive(Debug, Clone, Default)]
pub struct MeshVisibleFacesRegistry {
    providers: Vec<RegisteredVisibleFaces>,
}

impl MeshVisibleFacesRegistry {
    /// Empty registry.
    pub const fn new() -> Self {
        Self {
            providers: Vec::new(),
        }
    }

    /// Register an implementation. Higher priorities run first.
    pub fn register<P>(&mut self, priority: i32, provider: P)
    where
        P: MeshVisibleFaces + 'static,
    {
        self.register_arc(priority, Arc::new(provider));
    }

    /// Validate a provider with the shared conformance suite before registration.
    pub fn register_conformant<P>(
        &mut self,
        priority: i32,
        provider: P,
    ) -> Result<(), Box<FaceConformanceReport>>
    where
        P: MeshVisibleFaces + 'static,
    {
        let report = FaceConformanceSuite::run(&provider);
        if !report.is_success() {
            return Err(Box::new(report));
        }
        self.register(priority, provider);
        Ok(())
    }

    /// Register a shared implementation. Higher priorities run first.
    pub fn register_arc(&mut self, priority: i32, provider: Arc<dyn MeshVisibleFaces>) {
        self.providers
            .push(RegisteredVisibleFaces { priority, provider });
        self.providers
            .sort_by_key(|entry| std::cmp::Reverse(entry.priority));
    }

    /// Registered providers in dispatch order.
    pub fn providers(&self) -> impl Iterator<Item = &dyn MeshVisibleFaces> {
        self.providers.iter().map(|entry| entry.provider.as_ref())
    }

    /// Validate and execute one visible-face request.
    pub fn visible_faces(
        &self,
        request: &VisibleFaceRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<VisibleFaceOutcome, HiddenLineError> {
        options.check_cancelled()?;
        validate_request_shape(request, options)?;

        let mut has_matching_provider = false;
        let mut has_budgeted_provider = false;
        for entry in &self.providers {
            let descriptor = entry.provider.descriptor();
            if matches_device(options.device(), descriptor.id, descriptor.target) {
                has_matching_provider = true;
                has_budgeted_provider |= scratch_fits(
                    entry.provider.scratch_requirement(),
                    request.meshes,
                    options,
                );
            }
        }
        if !has_matching_provider {
            return Err(unsupported());
        }
        if !has_budgeted_provider {
            return Err(GeomError::BudgetExceeded { resource: "memory" }.into());
        }

        // The allocating topology audit happens only after at least one provider
        // and the audit itself have been admitted by the memory budget.
        validate_topology(request.meshes, options)?;

        let mut last_retryable = None;
        for entry in &self.providers {
            let descriptor = entry.provider.descriptor();
            if !matches_device(options.device(), descriptor.id, descriptor.target)
                || !scratch_fits(
                    entry.provider.scratch_requirement(),
                    request.meshes,
                    options,
                )
            {
                continue;
            }
            let result = entry.provider.visible_faces(request, options);
            if let Err(HiddenLineError { mesh: Some(id), .. }) = &result {
                if !request.meshes.iter().any(|entry| entry.id == *id) {
                    return Err(contract_violation(
                        descriptor.id,
                        "refusal names a mesh that is not in the request",
                    ));
                }
            }
            match result {
                Ok(outcome) => {
                    validate_outcome(descriptor.id, request, &outcome)?;
                    return Ok(outcome);
                }
                Err(
                    error @ HiddenLineError {
                        error: GeomError::Unsupported { .. } | GeomError::Unavailable { .. },
                        ..
                    },
                ) => {
                    last_retryable = Some(error);
                }
                Err(error) => return Err(error),
            }
        }

        Err(last_retryable.unwrap_or_else(unsupported))
    }
}

fn unsupported() -> HiddenLineError {
    GeomError::Unsupported {
        backend: BackendId::new("mesh-visible-faces-registry"),
        operation: Operation::MeshVisibleFaces,
    }
    .into()
}

fn contract_violation(backend: BackendId, detail: &str) -> HiddenLineError {
    GeomError::BackendContractViolation {
        backend,
        detail: detail.into(),
    }
    .into()
}

fn validate_request_shape(
    request: &VisibleFaceRequest<'_>,
    options: &ExecutionOptions,
) -> Result<(), HiddenLineError> {
    let limits = request.limits;
    validate_meshes(
        request.meshes,
        &SharedLimits {
            label: LABEL,
            max_meshes: limits.max_meshes,
            max_source_vertices: limits.max_source_vertices,
            max_source_triangles: limits.max_source_triangles,
            meshes_resource: "visible-face meshes",
            vertices_resource: "visible-face source vertices",
            triangles_resource: "visible-face source triangles",
        },
    )?;
    validate_view(&request.view, options.tolerance().angular(), LABEL)?;
    validate_region(request.region, LABEL)?;
    for entry in request.meshes {
        if let FaceGroups::PerTriangle(groups) = entry.groups {
            if groups.len() != entry.mesh.triangle_count() {
                return Err(HiddenLineError::mesh(
                    entry.id,
                    GeomError::InvalidInput(format!(
                        "{LABEL} mesh has {} face group(s) for {} triangle(s)",
                        groups.len(),
                        entry.mesh.triangle_count()
                    )),
                ));
            }
        }
    }
    Ok(())
}

fn point_order(a: &Point2, b: &Point2) -> Ordering {
    a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y))
}

fn sequence_order(a: &[Point2], b: &[Point2]) -> Ordering {
    a.iter()
        .zip(b)
        .map(|(p, q)| point_order(p, q))
        .find(|o| o.is_ne())
        .unwrap_or_else(|| a.len().cmp(&b.len()))
}

fn validate_outcome(
    backend: BackendId,
    request: &VisibleFaceRequest<'_>,
    outcome: &VisibleFaceOutcome,
) -> Result<(), HiddenLineError> {
    let violation = |detail: &str| contract_violation(backend, detail);
    let limits = request.limits;
    if outcome.view != request.view {
        return Err(violation("provider changed the requested view"));
    }
    if outcome.regions.len() > limits.max_regions {
        return Err(GeomError::BudgetExceeded {
            resource: "visible-face regions",
        }
        .into());
    }
    let mut vertices = 0usize;
    let mut caps = 0usize;
    let mut previous: Option<(usize, FaceSource, &[Point2])> = None;
    for region in &outcome.regions {
        let Some(position) = request
            .meshes
            .iter()
            .position(|entry| entry.id == region.mesh)
        else {
            return Err(violation("region names a mesh that is not in the request"));
        };
        let entry: &HiddenLineMesh<'_> = &request.meshes[position];
        for ring in std::iter::once(&region.outer).chain(&region.holes) {
            if ring.len() < 3 {
                return Err(violation("region ring has fewer than three points"));
            }
            if ring.iter().any(|p| !p.is_finite()) {
                return Err(violation("region ring contains a non-finite point"));
            }
            if ring
                .iter()
                .zip(ring.iter().cycle().skip(1))
                .any(|(a, b)| a == b)
            {
                return Err(violation("region ring repeats a point in a row"));
            }
            vertices = vertices.saturating_add(ring.len());
        }
        let count = entry.mesh.triangle_count();
        if region.triangles.windows(2).any(|w| w[0] >= w[1]) {
            return Err(violation("region triangles are not strictly ascending"));
        }
        if region.triangles.iter().any(|&t| t as usize >= count) {
            return Err(violation("region triangle is not in its mesh"));
        }
        match region.source {
            FaceSource::Group(group) => {
                if region
                    .triangles
                    .iter()
                    .any(|&t| entry.groups.of(t as usize) != group)
                {
                    return Err(violation("region triangle is not in the region's group"));
                }
            }
            FaceSource::Cap => {
                if entry.role != MeshRole::Solid || request.view.near.is_none() {
                    return Err(violation(
                        "cap region of a mesh that is not a solid cut by a near plane",
                    ));
                }
                if !region.triangles.is_empty() {
                    return Err(violation("cap region lists triangles"));
                }
                caps += 1;
            }
            _ => return Err(violation("region has an unknown face source")),
        }
        if let Some((last_position, last_source, last_outer)) = previous {
            let order = last_position
                .cmp(&position)
                .then(last_source.cmp(&region.source))
                .then_with(|| sequence_order(last_outer, &region.outer));
            if order.is_ge() {
                return Err(violation(
                    "regions are not ordered by mesh, source and outer ring",
                ));
            }
        }
        previous = Some((position, region.source, &region.outer));
    }
    if vertices > limits.max_region_vertices {
        return Err(GeomError::BudgetExceeded {
            resource: "visible-face region vertices",
        }
        .into());
    }
    let evidence = &outcome.evidence;
    if evidence.candidate_pairs > limits.max_candidate_pairs {
        return Err(GeomError::BudgetExceeded {
            resource: "visible-face candidate pairs",
        }
        .into());
    }
    let triangles = request
        .meshes
        .iter()
        .try_fold(0usize, |t, m| t.checked_add(m.mesh.triangle_count()));
    if !evidence.is_derived_from_input_mesh()
        || evidence.source_meshes != request.meshes.len()
        || Some(evidence.source_triangles) != triangles
        || evidence.regions != outcome.regions.len()
        || evidence.cap_regions != caps
    {
        return Err(violation(
            "visible-face evidence does not match input/output counts",
        ));
    }
    Ok(())
}
