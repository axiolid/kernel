//! Request validation and scene preparation, shared by the providers of the
//! orthographic view operations (hidden lines now, visible faces next).
//!
//! `prepare` checks the shape of a request (counts against the limits, ids,
//! frame, depth range, crop), audits every mesh, then builds what both
//! operations start from: the exact view-space vertices, the triangles with
//! their exact image orientation, and the cross-sections of the solids a near
//! plane cuts. Edge collection is left to `edges`, which `draw` runs on the
//! result. The request types differ per operation, so the pieces they share
//! come in as plain values: the meshes, the view, the crop, the three common
//! limits and the resource names a budget refusal carries.

use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::ops::Range;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Aabb2, Frame3, SpaceFrame, Tolerance};
use axiolid_mesh::try_audit_mesh;
use axiolid_mesh_hidden_lines_contract::{HiddenLineError, HiddenLineMesh, MeshRole, View};

use crate::caps::Cap;
use crate::decide::decide;
use crate::edges::{push_triangles, Tri};
use crate::recipe::{Constraint, Kind, ValueSign};
use crate::view::Vert;

/// The limits both operations share.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SharedLimits {
    pub max_meshes: usize,
    pub max_source_vertices: usize,
    pub max_source_triangles: usize,
}

/// How an operation names itself in a refusal.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ResourceNames {
    /// Prefix of a message: "hidden-line".
    pub label: &'static str,
    /// Budget resource of the mesh count.
    pub meshes: &'static str,
    /// Budget resource of the vertex total.
    pub vertices: &'static str,
    /// Budget resource of the triangle total.
    pub triangles: &'static str,
}

impl ResourceNames {
    pub(crate) const HIDDEN_LINES: Self = Self {
        label: "hidden-line",
        meshes: "hidden-line meshes",
        vertices: "hidden-line source vertices",
        triangles: "hidden-line source triangles",
    };

    pub(crate) const VISIBLE_FACES: Self = Self {
        label: "visible-face",
        meshes: "visible-face meshes",
        vertices: "visible-face source vertices",
        triangles: "visible-face source triangles",
    };
}

/// Where one mesh sits in the prepared lists.
#[derive(Debug, Clone)]
pub(crate) struct MeshRange {
    /// Scene vertex id of the mesh's first position.
    pub base: u32,
    /// The mesh's triangles in [`Prepared::tris`].
    pub triangles: Range<usize>,
}

/// The common start of both operations.
pub(crate) struct Prepared {
    /// View-space vertices of every mesh, one after another.
    pub verts: Vec<Vert>,
    /// Triangles of every mesh, one after another.
    pub tris: Vec<Tri>,
    /// The solids the near plane cuts.
    pub caps: Vec<Cap>,
    /// One range per request mesh.
    pub ranges: Vec<MeshRange>,
    /// Exact-tier decisions taken so far.
    pub escalations: usize,
    pub total_triangles: usize,
}

/// Validate the request pieces, audit the meshes and build the scene.
pub(crate) fn prepare(
    meshes: &[HiddenLineMesh<'_>],
    view: &View,
    region: Option<Aabb2>,
    limits: SharedLimits,
    names: &ResourceNames,
    options: &ExecutionOptions,
) -> Result<Prepared, HiddenLineError> {
    options.check_cancelled()?;
    let (total_vertices, total_triangles) =
        validate_shape(meshes, view, region, limits, names, options)?;
    for entry in meshes {
        options.check_cancelled()?;
        audit(entry, options)?;
    }
    let memory = || HiddenLineError::from(GeomError::BudgetExceeded { resource: "memory" });
    let mut verts: Vec<Vert> = Vec::new();
    verts
        .try_reserve_exact(total_vertices)
        .map_err(|_| memory())?;
    let mut tris: Vec<Tri> = Vec::new();
    let mut caps: Vec<Cap> = Vec::new();
    let mut ranges: Vec<MeshRange> = Vec::new();
    ranges
        .try_reserve_exact(meshes.len())
        .map_err(|_| memory())?;
    let mut escalations = 0;
    for entry in meshes {
        let base = u32::try_from(verts.len()).map_err(|_| {
            HiddenLineError::from(GeomError::BudgetExceeded {
                resource: names.vertices,
            })
        })?;
        verts.extend(
            entry
                .mesh
                .positions
                .iter()
                .map(|p| Vert::new(*p, &view.frame)),
        );
        let cut = match (view.near, entry.role) {
            (Some(near), MeshRole::Solid) => {
                cut_by_near(&verts[base as usize..], near, &mut escalations)
                    .map_err(|error| HiddenLineError::mesh(entry.id, error))?
            }
            _ => false,
        };
        let first_triangle = tris.len();
        push_triangles(
            entry.id,
            entry.mesh,
            base,
            &verts,
            options,
            &mut escalations,
            &mut tris,
        )?;
        let triangles = first_triangle..tris.len();
        if cut {
            caps.push(Cap {
                tris: triangles.clone(),
            });
        }
        ranges.push(MeshRange { base, triangles });
    }
    debug_assert_eq!(tris.len(), total_triangles);
    Ok(Prepared {
        verts,
        tris,
        caps,
        ranges,
        escalations,
        total_triangles,
    })
}

/// Whether the plane `W = near` cuts the vertices: one has `W > near` and one
/// has `W <= near`, each decided exactly.
fn cut_by_near(verts: &[Vert], near: f64, escalations: &mut usize) -> GeomResult<bool> {
    let (mut beyond, mut kept) = (false, false);
    for at in 0..verts.len() as u32 {
        let constraint = Constraint::new(Kind::BelowNear(near), false);
        let side = decide(
            &ValueSign {
                verts,
                constraint: &constraint,
                at,
            },
            escalations,
        )?;
        // `near - W`: negative means beyond the plane.
        if side == Ordering::Less {
            beyond = true;
        } else {
            kept = true;
        }
        if beyond && kept {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Request-level checks. Returns the total positions and triangles.
fn validate_shape(
    meshes: &[HiddenLineMesh<'_>],
    view: &View,
    region: Option<Aabb2>,
    limits: SharedLimits,
    names: &ResourceNames,
    options: &ExecutionOptions,
) -> Result<(usize, usize), HiddenLineError> {
    let label = names.label;
    let invalid = |detail: String| HiddenLineError::from(GeomError::InvalidInput(detail));
    let budget = |resource| HiddenLineError::from(GeomError::BudgetExceeded { resource });
    if meshes.is_empty() {
        return Err(invalid(format!("{label} request has no meshes")));
    }
    if meshes.len() > limits.max_meshes {
        return Err(budget(names.meshes));
    }
    let mut vertices = 0usize;
    let mut triangles = 0usize;
    for entry in meshes {
        vertices = vertices
            .checked_add(entry.mesh.positions.len())
            .ok_or_else(|| budget(names.vertices))?;
        triangles = triangles
            .checked_add(entry.mesh.triangle_count())
            .ok_or_else(|| budget(names.triangles))?;
    }
    if vertices > limits.max_source_vertices {
        return Err(budget(names.vertices));
    }
    if triangles > limits.max_source_triangles {
        return Err(budget(names.triangles));
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

    validate_frame(&view.frame, options.tolerance(), label)?;
    for (name, bound) in [("near", view.near), ("far", view.far)] {
        if bound.is_some_and(|value| !value.is_finite()) {
            return Err(invalid(format!("{label} {name} depth is not finite")));
        }
    }
    if let (Some(near), Some(far)) = (view.near, view.far) {
        if far >= near {
            return Err(invalid(format!(
                "{label} far depth must be smaller than the near depth"
            )));
        }
    }
    if let Some(region) = region {
        let finite = region.min.is_finite() && region.max.is_finite();
        if !finite || region.min.x > region.max.x || region.min.y > region.max.y {
            return Err(invalid(format!(
                "{label} region must be finite with min not above max"
            )));
        }
    }
    Ok((vertices, triangles))
}

/// A finite right-handed orthonormal frame, as the registry requires.
fn validate_frame(
    frame: &Frame3,
    tolerance: Tolerance,
    label: &str,
) -> Result<(), HiddenLineError> {
    // A coarse caller tolerance must not admit zero or scaled axes.
    let limit = tolerance.angular().clamp(128.0 * f64::EPSILON, 1.0e-6);
    let tolerance = Tolerance::new(limit, limit).map_err(|_| {
        HiddenLineError::from(GeomError::InvalidInput(format!(
            "{label} view frame tolerance is invalid"
        )))
    })?;
    SpaceFrame::new(frame.origin, frame.x, frame.y, frame.z, tolerance)
        .map(|_| ())
        .map_err(|error| {
            HiddenLineError::from(GeomError::InvalidInput(format!(
                "{label} view frame is invalid: {error}"
            )))
        })
}

/// Per-mesh audit: indices, positions, degenerate triangles, topology.
fn audit(entry: &HiddenLineMesh<'_>, options: &ExecutionOptions) -> Result<(), HiddenLineError> {
    let refuse = |error: GeomError| HiddenLineError::mesh(entry.id, error);
    let health = try_audit_mesh(entry.mesh, options.tolerance())
        .map_err(|_| HiddenLineError::from(GeomError::BudgetExceeded { resource: "memory" }))?;
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
    Ok(())
}
