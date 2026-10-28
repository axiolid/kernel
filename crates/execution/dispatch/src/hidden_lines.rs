//! Provider registration, validation, fallback, and budget policy for
//! orthographic hidden-line drawings.

use std::collections::BTreeSet;
use std::sync::Arc;

use axiolid_contracts::{BackendId, ExecutionOptions, GeomError, Operation, ScratchRequirement};
use axiolid_mesh::{audit_mesh_scratch_bytes, try_audit_mesh};
use axiolid_mesh_hidden_lines_contract::{
    conformance::{ConformanceReport, ConformanceSuite},
    EdgeKind, HiddenLineError, HiddenLineMesh, HiddenLineOutcome, HiddenLineRequest,
    MeshHiddenLines, MeshRole, Visibility,
};

use crate::device::matches_device;
use crate::frame::validate_frame;

#[derive(Debug, Clone)]
struct RegisteredHiddenLines {
    priority: i32,
    provider: Arc<dyn MeshHiddenLines>,
}

/// Ordered executable hidden-line providers.
///
/// Only unsupported/unavailable providers permit fallback. Geometry,
/// degeneracy, budget, cancellation, and contract errors fail immediately.
/// A refusal caused by one input mesh names it
/// ([`HiddenLineError::mesh`]); a refusal of the whole request does not.
#[derive(Debug, Clone, Default)]
pub struct MeshHiddenLinesRegistry {
    providers: Vec<RegisteredHiddenLines>,
}

impl MeshHiddenLinesRegistry {
    /// Empty registry.
    pub const fn new() -> Self {
        Self {
            providers: Vec::new(),
        }
    }

    /// Register an implementation. Higher priorities run first.
    pub fn register<P>(&mut self, priority: i32, provider: P)
    where
        P: MeshHiddenLines + 'static,
    {
        self.register_arc(priority, Arc::new(provider));
    }

    /// Validate a provider with the shared conformance suite before registration.
    pub fn register_conformant<P>(
        &mut self,
        priority: i32,
        provider: P,
    ) -> Result<(), Box<ConformanceReport>>
    where
        P: MeshHiddenLines + 'static,
    {
        let report = ConformanceSuite::run(&provider);
        if !report.is_success() {
            return Err(Box::new(report));
        }
        self.register(priority, provider);
        Ok(())
    }

    /// Register a shared implementation. Higher priorities run first.
    pub fn register_arc(&mut self, priority: i32, provider: Arc<dyn MeshHiddenLines>) {
        self.providers
            .push(RegisteredHiddenLines { priority, provider });
        self.providers
            .sort_by_key(|entry| std::cmp::Reverse(entry.priority));
    }

    /// Registered providers in dispatch order.
    pub fn providers(&self) -> impl Iterator<Item = &dyn MeshHiddenLines> {
        self.providers.iter().map(|entry| entry.provider.as_ref())
    }

    /// Validate and execute one hidden-line request.
    pub fn hidden_lines(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<HiddenLineOutcome, HiddenLineError> {
        options.check_cancelled()?;
        validate_request_shape(request, options)?;

        let mut has_matching_provider = false;
        let mut has_budgeted_provider = false;
        for entry in &self.providers {
            let descriptor = entry.provider.descriptor();
            if matches_device(options.device(), descriptor.id, descriptor.target) {
                has_matching_provider = true;
                has_budgeted_provider |=
                    scratch_fits(entry.provider.scratch_requirement(), request, options);
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
        validate_topology(request, options)?;

        let mut last_retryable = None;
        for entry in &self.providers {
            let descriptor = entry.provider.descriptor();
            if !matches_device(options.device(), descriptor.id, descriptor.target)
                || !scratch_fits(entry.provider.scratch_requirement(), request, options)
            {
                continue;
            }
            let result = entry.provider.hidden_lines(request, options);
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
        backend: BackendId::new("mesh-hidden-lines-registry"),
        operation: Operation::MeshHiddenLines,
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

fn invalid(detail: impl Into<String>) -> HiddenLineError {
    GeomError::InvalidInput(detail.into()).into()
}

fn validate_request_shape(
    request: &HiddenLineRequest<'_>,
    options: &ExecutionOptions,
) -> Result<(), HiddenLineError> {
    let limits = request.limits;
    if request.meshes.is_empty() {
        return Err(invalid("hidden-line request has no meshes"));
    }
    if request.meshes.len() > limits.max_meshes {
        return Err(GeomError::BudgetExceeded {
            resource: "hidden-line meshes",
        }
        .into());
    }
    let (vertices, triangles) = source_totals(request.meshes);
    if vertices.is_none_or(|n| n > limits.max_source_vertices) {
        return Err(GeomError::BudgetExceeded {
            resource: "hidden-line source vertices",
        }
        .into());
    }
    if triangles.is_none_or(|n| n > limits.max_source_triangles) {
        return Err(GeomError::BudgetExceeded {
            resource: "hidden-line source triangles",
        }
        .into());
    }

    let mut ids = BTreeSet::new();
    for entry in request.meshes {
        if !ids.insert(entry.id) {
            return Err(HiddenLineError::mesh(
                entry.id,
                GeomError::InvalidInput("duplicate mesh id".into()),
            ));
        }
        if entry.mesh.indices.is_empty() || !entry.mesh.indices.len().is_multiple_of(3) {
            return Err(HiddenLineError::mesh(
                entry.id,
                GeomError::InvalidInput(
                    "hidden-line mesh has no triangles or a partial triangle".into(),
                ),
            ));
        }
    }

    let view = &request.view;
    validate_frame(
        view.frame,
        options.tolerance().angular(),
        "hidden-line view frame",
    )?;
    for (name, bound) in [("near", view.near), ("far", view.far)] {
        if bound.is_some_and(|value| !value.is_finite()) {
            return Err(invalid(format!("hidden-line {name} depth is not finite")));
        }
    }
    if let (Some(near), Some(far)) = (view.near, view.far) {
        if far >= near {
            return Err(invalid(
                "hidden-line far depth must be smaller than the near depth",
            ));
        }
    }
    if request
        .edges
        .crease_cosine
        .is_some_and(|cosine| !(-1.0..=1.0).contains(&cosine))
    {
        return Err(invalid("hidden-line crease cosine must lie in [-1, 1]"));
    }
    if let Some(region) = request.region {
        let finite = region.min.is_finite() && region.max.is_finite();
        if !finite || region.min.x > region.max.x || region.min.y > region.max.y {
            return Err(invalid(
                "hidden-line region must be finite with min not above max",
            ));
        }
    }
    Ok(())
}

/// Total positions and triangles; `None` on overflow.
fn source_totals(meshes: &[HiddenLineMesh<'_>]) -> (Option<usize>, Option<usize>) {
    meshes
        .iter()
        .fold((Some(0usize), Some(0usize)), |(v, t), m| {
            (
                v.and_then(|v| v.checked_add(m.mesh.positions.len())),
                t.and_then(|t| t.checked_add(m.mesh.triangle_count())),
            )
        })
}

fn validate_topology(
    request: &HiddenLineRequest<'_>,
    options: &ExecutionOptions,
) -> Result<(), HiddenLineError> {
    for entry in request.meshes {
        options.check_cancelled()?;
        let refuse = |error: GeomError| HiddenLineError::mesh(entry.id, error);
        let health = try_audit_mesh(entry.mesh, options.tolerance())
            .map_err(|_| GeomError::BudgetExceeded { resource: "memory" })?;
        if let Some((triangle, index)) = health.first_invalid_index {
            return Err(refuse(GeomError::InvalidInput(format!(
                "triangle {triangle} references position {index}, but the mesh has {}",
                health.positions
            ))));
        }
        if let Some(position) = health.first_non_finite_position {
            return Err(refuse(GeomError::InvalidInput(format!(
                "position {position} is not finite"
            ))));
        }
        if health.degenerate_triangles > 0 {
            return Err(refuse(GeomError::Degenerate(format!(
                "{} degenerate triangle(s)",
                health.degenerate_triangles
            ))));
        }
        if health.non_manifold_edges > 0 || health.inconsistent_winding_edges > 0 {
            return Err(refuse(GeomError::NotManifold(format!(
                "non_manifold={}, inconsistent_winding={}",
                health.non_manifold_edges, health.inconsistent_winding_edges
            ))));
        }
        if entry.role == MeshRole::Solid && health.boundary_edges > 0 {
            return Err(refuse(GeomError::NotManifold(format!(
                "a solid must be closed; boundary={}",
                health.boundary_edges
            ))));
        }
    }
    Ok(())
}

fn scratch_fits(
    provider: ScratchRequirement,
    request: &HiddenLineRequest<'_>,
    options: &ExecutionOptions,
) -> bool {
    let Some(budget) = options.memory_budget_bytes() else {
        return true;
    };
    let (vertices, triangles) = source_totals(request.meshes);
    let (Some(vertices), Some(_)) = (vertices, triangles) else {
        return false;
    };
    // Meshes are audited one after another, so the audit peak is the largest.
    let mut audit_bytes = 0usize;
    for entry in request.meshes {
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

fn validate_outcome(
    backend: BackendId,
    request: &HiddenLineRequest<'_>,
    outcome: &HiddenLineOutcome,
) -> Result<(), HiddenLineError> {
    let violation = |detail: &str| contract_violation(backend, detail);
    if outcome.view != request.view {
        return Err(violation("provider changed the requested view"));
    }
    if outcome.segments.len() > request.limits.max_output_segments {
        return Err(GeomError::BudgetExceeded {
            resource: "hidden-line output segments",
        }
        .into());
    }
    let selected = |kind: EdgeKind| match kind {
        EdgeKind::Boundary => request.edges.boundary,
        EdgeKind::Silhouette => request.edges.silhouette,
        EdgeKind::Crease => request.edges.crease_cosine.is_some(),
        _ => false,
    };
    let mut previous: Option<(usize, [u32; 2])> = None;
    let mut previous_kind = None;
    let mut visible = 0usize;
    let mut hidden = 0usize;
    for segment in &outcome.segments {
        let Some(position) = request
            .meshes
            .iter()
            .position(|entry| entry.id == segment.source.mesh)
        else {
            return Err(violation("segment names a mesh that is not in the request"));
        };
        let [a, b] = segment.source.vertices;
        let positions = request.meshes[position].mesh.positions.len();
        if a >= b || b as usize >= positions {
            return Err(violation(
                "segment source edge must be ascending and within the mesh",
            ));
        }
        if !segment.start.is_finite() || !segment.end.is_finite() {
            return Err(violation("segment contains a non-finite point"));
        }
        if segment.start == segment.end {
            return Err(violation("segment has zero length"));
        }
        if !selected(segment.kind) {
            return Err(violation("segment has an edge kind that was not selected"));
        }
        match segment.visibility {
            Visibility::Visible => visible += 1,
            Visibility::Hidden if request.edges.hidden => hidden += 1,
            _ => return Err(violation("hidden segment returned although not requested")),
        }
        let key = (position, segment.source.vertices);
        if previous.is_some_and(|last| last > key) {
            return Err(violation(
                "segments are not ordered by mesh and source edge",
            ));
        }
        if previous == Some(key) && previous_kind != Some(segment.kind) {
            return Err(violation("pieces of one source edge differ in edge kind"));
        }
        previous = Some(key);
        previous_kind = Some(segment.kind);
    }
    let triangles = source_totals(request.meshes).1;
    let evidence = &outcome.evidence;
    if !evidence.is_derived_from_input_mesh()
        || evidence.source_meshes != request.meshes.len()
        || Some(evidence.source_triangles) != triangles
        || evidence.visible_segments != visible
        || evidence.hidden_segments != hidden
    {
        return Err(violation(
            "hidden-line evidence does not match input/output counts",
        ));
    }
    Ok(())
}
