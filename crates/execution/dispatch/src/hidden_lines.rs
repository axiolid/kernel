//! Provider registration, validation, fallback, and budget policy for
//! orthographic hidden-line drawings.

use std::sync::Arc;

use axiolid_contracts::{BackendId, ExecutionOptions, GeomError, Operation};
use axiolid_mesh::audit_mesh_scratch_bytes;
use axiolid_mesh_hidden_lines_contract::{
    conformance::{
        perspective::{
            PerspectiveConformanceFailure, PerspectiveConformanceReport,
            PerspectiveConformanceSuite,
        },
        ConformanceReport, ConformanceSuite,
    },
    EdgeKind, HiddenLineError, HiddenLineOutcome, HiddenLineRequest, MeshHiddenLines, RefusedMesh,
    Visibility,
};

use crate::device::matches_device;
use crate::view_checks::{
    invalid, mesh_verdict, scratch_fits, source_totals, validate_meshes, validate_region,
    validate_topology, validate_view, SharedLimits,
};

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
/// [`hidden_lines`](Self::hidden_lines) names the first refused mesh only; a
/// caller that drops refused meshes and draws the rest learns all of them in
/// one call from [`audit`](Self::audit).
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

    /// Validate a provider with the shared suite and the perspective suite
    /// before registration (ADR 0107). The claim is the provider's
    /// `supports_projection(Projection::Perspective)`: a provider that does
    /// not claim it must refuse perspective views by name, one that does must
    /// pass the perspective fixtures. Failures of either suite are in the
    /// report.
    pub fn register_perspective_conformant<P>(
        &mut self,
        priority: i32,
        provider: P,
    ) -> Result<(), Box<PerspectiveConformanceReport>>
    where
        P: MeshHiddenLines + 'static,
    {
        let mut report = PerspectiveConformanceSuite::run(&provider);
        report.failures.extend(
            ConformanceSuite::run(&provider)
                .failures
                .into_iter()
                .map(PerspectiveConformanceFailure::HiddenLinesBase),
        );
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

    /// Every input mesh that [`hidden_lines`](Self::hidden_lines) would
    /// refuse, in request order, without searching.
    ///
    /// Runs only the per-mesh audit of shape and topology, through the code
    /// `hidden_lines` uses, so a mesh is listed exactly when `hidden_lines`
    /// would refuse a request holding it. An empty list means every mesh
    /// passes the audit (a provider can still refuse). A caller that leaves
    /// the listed meshes out and calls `hidden_lines` again draws the rest;
    /// no substitute drawing is made. A refusal of the whole request (limits,
    /// frame, memory, cancellation, an empty request, a duplicate id or a mesh
    /// without triangles) is returned as `Err`, as it is by `hidden_lines`.
    pub fn audit(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<Vec<RefusedMesh>, HiddenLineError> {
        options.check_cancelled()?;
        validate_request_shape(request, options)?;
        if let Some(budget) = options.memory_budget_bytes() {
            let fits = request.meshes.iter().all(|entry| {
                audit_mesh_scratch_bytes(entry.mesh.triangle_count())
                    .is_some_and(|bytes| bytes <= budget)
            });
            if !fits {
                return Err(GeomError::BudgetExceeded { resource: "memory" }.into());
            }
        }
        let groups = request.edges.group_boundary.then_some(LABEL);
        let mut refused = Vec::new();
        for entry in request.meshes {
            options.check_cancelled()?;
            if let Some(error) = mesh_verdict(entry, groups, options)? {
                refused.push(RefusedMesh::new(entry.id, error));
            }
        }
        Ok(refused)
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
            if matches_device(options.device(), descriptor.id, descriptor.target)
                && entry.provider.supports_projection(request.view.projection)
            {
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
        validate_topology(
            request.meshes,
            request.edges.group_boundary.then_some(LABEL),
            options,
        )?;

        let mut last_retryable = None;
        for entry in &self.providers {
            let descriptor = entry.provider.descriptor();
            if !matches_device(options.device(), descriptor.id, descriptor.target)
                || !entry.provider.supports_projection(request.view.projection)
                || !scratch_fits(
                    entry.provider.scratch_requirement(),
                    request.meshes,
                    options,
                )
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

const LABEL: &str = "hidden-line";

fn validate_request_shape(
    request: &HiddenLineRequest<'_>,
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
            meshes_resource: "hidden-line meshes",
            vertices_resource: "hidden-line source vertices",
            triangles_resource: "hidden-line source triangles",
        },
    )?;
    validate_view(&request.view, options.tolerance().angular(), LABEL)?;
    if request
        .edges
        .crease_cosine
        .is_some_and(|cosine| !(-1.0..=1.0).contains(&cosine))
    {
        return Err(invalid("hidden-line crease cosine must lie in [-1, 1]"));
    }
    validate_region(request.region, LABEL)
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
        EdgeKind::GroupBoundary => request.edges.group_boundary,
        _ => false,
    };
    let mut previous: Option<(usize, [u32; 2])> = None;
    let mut previous_kind = None;
    let mut previous_visibility = None;
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
        // With hidden pieces omitted, two visible pieces may be neighbours
        // along the edge with a hidden one between; otherwise a piece is
        // maximal, so neighbours differ.
        if request.edges.hidden
            && previous == Some(key)
            && previous_visibility == Some(segment.visibility)
        {
            return Err(violation(
                "adjacent pieces of one source edge have the same visibility",
            ));
        }
        previous = Some(key);
        previous_kind = Some(segment.kind);
        previous_visibility = Some(segment.visibility);
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
