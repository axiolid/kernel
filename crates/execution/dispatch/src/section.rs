//! Provider registration, validation, fallback, and budget policy for mesh sections.

use std::sync::Arc;

use axiolid_contracts::{
    BackendId, ExecutionOptions, GeomError, GeomResult, Operation, ScratchRequirement,
};
use axiolid_core::{Frame3, SpaceFrame, Tolerance};
use axiolid_mesh::{audit_mesh_scratch_bytes, try_audit_mesh, TriMesh};
use axiolid_mesh_contracts::SolidRequirements;
use axiolid_mesh_section_contract::{
    conformance::{ConformanceReport, ConformanceSuite},
    MeshPlaneSection, SectionLimits, SectionOutcome,
};

use crate::device::matches_device;

#[derive(Debug, Clone)]
struct RegisteredSection {
    priority: i32,
    provider: Arc<dyn MeshPlaneSection>,
}

/// Ordered executable mesh plane-section providers.
///
/// Only unsupported/unavailable providers permit fallback. Geometry,
/// degeneracy, budget, cancellation, and contract errors fail immediately.
#[derive(Debug, Clone, Default)]
pub struct MeshPlaneSectionRegistry {
    providers: Vec<RegisteredSection>,
}

impl MeshPlaneSectionRegistry {
    /// Empty registry.
    pub const fn new() -> Self {
        Self {
            providers: Vec::new(),
        }
    }

    /// Register an implementation. Higher priorities run first.
    pub fn register<P>(&mut self, priority: i32, provider: P)
    where
        P: MeshPlaneSection + 'static,
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
        P: MeshPlaneSection + 'static,
    {
        let report = ConformanceSuite::run(&provider);
        if !report.is_success() {
            return Err(Box::new(report));
        }
        self.register(priority, provider);
        Ok(())
    }

    /// Register a shared implementation. Higher priorities run first.
    pub fn register_arc(&mut self, priority: i32, provider: Arc<dyn MeshPlaneSection>) {
        self.providers
            .push(RegisteredSection { priority, provider });
        self.providers
            .sort_by_key(|entry| std::cmp::Reverse(entry.priority));
    }

    /// Registered providers in dispatch order.
    pub fn providers(&self) -> impl Iterator<Item = &dyn MeshPlaneSection> {
        self.providers.iter().map(|entry| entry.provider.as_ref())
    }

    /// Validate and execute one plane-section request.
    pub fn section(
        &self,
        mesh: &TriMesh,
        frame: Frame3,
        limits: SectionLimits,
        options: &ExecutionOptions,
    ) -> GeomResult<SectionOutcome> {
        options.check_cancelled()?;
        validate_source_shape(mesh, frame, limits, options)?;

        let mut has_matching_provider = false;
        let mut has_budgeted_provider = false;
        for entry in &self.providers {
            let descriptor = entry.provider.descriptor();
            if matches_device(options.device(), descriptor.id, descriptor.target) {
                has_matching_provider = true;
                has_budgeted_provider |=
                    section_scratch_fits(entry.provider.scratch_requirement(), mesh, options);
            }
        }
        if !has_matching_provider {
            return Err(GeomError::Unsupported {
                backend: BackendId::new("mesh-section-registry"),
                operation: Operation::MeshPlaneSection,
            });
        }
        if !has_budgeted_provider {
            return Err(GeomError::BudgetExceeded { resource: "memory" });
        }

        // The allocating topology audit happens only after at least one provider
        // and the audit itself have been admitted by the memory budget.
        validate_source_topology(mesh, options)?;

        let mut last_retryable = None;
        for entry in &self.providers {
            let descriptor = entry.provider.descriptor();
            if !matches_device(options.device(), descriptor.id, descriptor.target)
                || !section_scratch_fits(entry.provider.scratch_requirement(), mesh, options)
            {
                continue;
            }
            match entry.provider.section(mesh, frame, limits, options) {
                Ok(outcome) => {
                    validate_outcome(descriptor.id, mesh, frame, limits, &outcome)?;
                    return Ok(outcome);
                }
                Err(error @ (GeomError::Unsupported { .. } | GeomError::Unavailable { .. })) => {
                    last_retryable = Some(error);
                }
                Err(error) => return Err(error),
            }
        }

        Err(last_retryable.unwrap_or(GeomError::Unsupported {
            backend: BackendId::new("mesh-section-registry"),
            operation: Operation::MeshPlaneSection,
        }))
    }
}

fn validate_source_shape(
    mesh: &TriMesh,
    frame: Frame3,
    limits: SectionLimits,
    options: &ExecutionOptions,
) -> GeomResult<()> {
    if mesh.positions.len() > limits.max_source_vertices {
        return Err(GeomError::BudgetExceeded {
            resource: "section source vertices",
        });
    }
    if mesh.triangle_count() > limits.max_source_triangles {
        return Err(GeomError::BudgetExceeded {
            resource: "section source triangles",
        });
    }
    if mesh.indices.is_empty() {
        return Err(GeomError::InvalidInput(
            "section source mesh has no triangles".into(),
        ));
    }
    validate_frame(frame, options.tolerance().angular())?;

    Ok(())
}

fn validate_source_topology(mesh: &TriMesh, options: &ExecutionOptions) -> GeomResult<()> {
    let health = try_audit_mesh(mesh, options.tolerance())
        .map_err(|_| GeomError::BudgetExceeded { resource: "memory" })?;
    if !health.is_closed_two_manifold() {
        return Err(GeomError::NotManifold(format!(
            "plane section requires a closed consistently wound two-manifold mesh; \
             boundary={}, non_manifold={}, inconsistent_winding={}, degenerate={}",
            health.boundary_edges,
            health.non_manifold_edges,
            health.inconsistent_winding_edges,
            health.degenerate_triangles
        )));
    }
    SolidRequirements::Oriented.validate(mesh, "section source")?;
    Ok(())
}

fn section_scratch_fits(
    provider: ScratchRequirement,
    mesh: &TriMesh,
    options: &ExecutionOptions,
) -> bool {
    let Some(budget) = options.memory_budget_bytes() else {
        return true;
    };
    let Some(audit_bytes) = audit_mesh_scratch_bytes(mesh.triangle_count()) else {
        return false;
    };
    let elements = mesh.positions.len().max(mesh.triangle_count());
    let Some(provider_bytes) = provider.upper_bound_bytes(elements) else {
        return false;
    };
    // Audit and provider phases are sequential, so their peak is the maximum,
    // not their sum. Output storage is outside the scratch budget by contract.
    audit_bytes.max(provider_bytes) <= budget
}

fn validate_frame(frame: Frame3, angular_tolerance: f64) -> GeomResult<()> {
    // Structural frame validity cannot be relaxed into zero/scaled axes by a
    // coarse geometric tolerance, so the caller's angular tolerance is capped
    // here before it reaches the shared validator. Tighter caller tolerances
    // still apply; looser ones do not.
    let limit = angular_tolerance.clamp(128.0 * f64::EPSILON, 1.0e-6);
    let tolerance = Tolerance::new(limit, limit)
        .map_err(|_| GeomError::InvalidInput("section frame tolerance is invalid".into()))?;
    SpaceFrame::new(frame.origin, frame.x, frame.y, frame.z, tolerance)
        .map(|_| ())
        .map_err(|error| GeomError::InvalidInput(format!("section frame is invalid: {error}")))
}

fn validate_outcome(
    backend: BackendId,
    mesh: &TriMesh,
    frame: Frame3,
    limits: SectionLimits,
    outcome: &SectionOutcome,
) -> GeomResult<()> {
    let violation = |detail: &str| GeomError::BackendContractViolation {
        backend,
        detail: detail.into(),
    };
    if outcome.frame != frame {
        return Err(violation("provider changed the requested section frame"));
    }
    if outcome.contours.len() > limits.max_contours {
        return Err(GeomError::BudgetExceeded {
            resource: "section contours",
        });
    }
    let mut vertices = 0usize;
    for contour in &outcome.contours {
        if contour.points.len() < 3 {
            return Err(violation(
                "section contour must be closed with at least three points",
            ));
        }
        if !contour.points.iter().all(|point| point.is_finite()) {
            return Err(violation("section contour contains a non-finite point"));
        }
        for index in 0..contour.points.len() {
            if contour.points[index] == contour.points[(index + 1) % contour.points.len()] {
                return Err(violation(
                    "section contour contains adjacent duplicate points",
                ));
            }
        }
        vertices = vertices
            .checked_add(contour.points.len())
            .ok_or(GeomError::BudgetExceeded {
                resource: "section output vertices",
            })?;
        if vertices > limits.max_output_vertices {
            return Err(GeomError::BudgetExceeded {
                resource: "section output vertices",
            });
        }
    }
    if outcome.evidence.source_triangles != mesh.triangle_count()
        || outcome.evidence.output_vertices != vertices
        || outcome.evidence.output_contours != outcome.contours.len()
        || !outcome.evidence.is_derived_from_input_mesh()
    {
        return Err(violation(
            "section evidence does not match input/output counts",
        ));
    }
    Ok(())
}
