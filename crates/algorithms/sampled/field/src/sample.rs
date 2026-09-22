//! Deterministic scalar CPU triangle coverage.
//!
//! Coverage emits [`SurfaceHit`]s only. A triangle has no thickness, so it can
//! never produce an occupancy span here; occupancy is a separate, explicit
//! construction over a closed shell (see [`LayeredField::derive_occupancy`]).

use axiolid_core::{Point3, Scalar, Vec3};

use crate::{
    FieldConfig, FieldEvidence, LayeredCell, LayeredField, LayeredFieldError, SurfaceFacing,
    SurfaceHit,
};

// Re-exported rather than redefined: this crate used to carry its own
// `Triangle3`, which made the most reusable 3D primitive in the kernel
// reachable only by depending on field sampling. The type now lives in
// `axiolid_core` beside the other primitives, and this alias keeps the
// existing `axiolid_field_ops::Triangle3` path working for consumers.
pub use axiolid_core::Triangle3;

/// The scalar reference provider for triangle coverage.
///
/// Traversal order is row-major over `(x, y)` and input order within a cell, so
/// repeated runs on identical input produce byte-identical fields.
#[derive(Debug, Clone, Copy, Default)]
pub struct CpuCoverageProvider;

impl CpuCoverageProvider {
    /// Construct the provider. It holds no state, no pool, and no global config.
    pub const fn new() -> Self {
        Self
    }

    /// Sample every cell centre against every triangle.
    pub fn sample(
        &self,
        config: &FieldConfig,
        triangles: &[Triangle3],
    ) -> Result<LayeredField, LayeredFieldError> {
        sample_triangles_cpu(config, triangles)
    }
}

/// Deterministic scalar triangle coverage over a validated configuration.
pub fn sample_triangles_cpu(
    config: &FieldConfig,
    triangles: &[Triangle3],
) -> Result<LayeredField, LayeredFieldError> {
    if triangles
        .iter()
        .any(|t| !t.a.is_finite() || !t.b.is_finite() || !t.c.is_finite())
    {
        return Err(LayeredFieldError::NonFiniteGeometry);
    }

    let tolerance = config.tolerance();
    let linear = tolerance.linear();
    let direction = config.frame().z;
    let span = config.bounds().normal_span();
    let (w_low, w_high) = (span.start, span.end);
    let budget = config.budget();

    // Precompute per-triangle plane data once; the inner loop is per cell.
    let mut prepared = Vec::with_capacity(triangles.len());
    let mut evidence = FieldEvidence::default();
    for triangle in triangles {
        let normal = triangle.normal();
        // Area scales as |normal| / 2; reject slivers using the linear tolerance.
        if !normal.is_finite() || normal.length() <= linear * linear {
            evidence.degenerate_triangles += 1;
            continue;
        }
        let denominator = normal.dot(direction);
        if denominator.abs() <= linear.max(Scalar::EPSILON) {
            evidence.parallel_triangles_skipped += 1;
            continue;
        }
        prepared.push(Prepared {
            triangle: *triangle,
            normal,
            denominator,
        });
    }

    let empty = LayeredField::with_config(config)?;
    let (width, height) = config.dimensions();
    let mut cells = empty.cells().to_vec();
    let mut stored = 0usize;

    for y in 0..height {
        for x in 0..width {
            let origin = config.cell_center(x, y);
            let mut hits: Vec<SurfaceHit> = Vec::new();
            for item in &prepared {
                let offset = item.triangle.a - origin;
                let w = item.normal.dot(offset) / item.denominator;
                if !w.is_finite() {
                    continue;
                }
                if w < w_low - linear || w > w_high + linear {
                    evidence.out_of_bounds_hits += 1;
                    continue;
                }
                let point = origin + direction * w;
                match classify(&item.triangle, item.normal, point, linear) {
                    Containment::Outside => continue,
                    Containment::Boundary => evidence.boundary_contacts += 1,
                    Containment::Interior => {}
                }
                hits.push(SurfaceHit::new(w, facing_of(item.denominator)));
            }
            evidence.cells_sampled += 1;
            // Two facets sharing an edge both report a crossing when the
            // sampling line passes through that edge. Collapse same-facing
            // crossings that agree within the linear tolerance: they describe
            // one surface. Distinct facings are never merged, because an
            // enter/exit pair at the same coordinate is a real thin feature.
            hits.sort_by(|left, right| {
                left.w()
                    .total_cmp(&right.w())
                    .then_with(|| left.facing().cmp(&right.facing()))
            });
            let before = hits.len();
            hits.dedup_by(|right, left| {
                left.facing() == right.facing() && (right.w() - left.w()).abs() <= linear
            });
            evidence.coincident_hits_merged += before - hits.len();
            evidence.surface_hits += hits.len();
            if hits.is_empty() {
                evidence.empty_cells += 1;
            } else if hits.len() > 1 {
                evidence.multi_layer_cells += 1;
            }
            stored = stored
                .checked_add(hits.len())
                .ok_or(LayeredFieldError::SampleBudgetExceeded)?;
            if stored > budget.max_intervals {
                return Err(LayeredFieldError::SampleBudgetExceeded);
            }
            let index = empty
                .linear_index(x, y)
                .ok_or(LayeredFieldError::NodeOutsideField)?;
            cells[index] = LayeredCell::with_layers(hits, Vec::new())?;
        }
    }

    LayeredField::from_cells(width, height, cells, evidence)
}

struct Prepared {
    triangle: Triangle3,
    normal: Vec3,
    denominator: Scalar,
}

enum Containment {
    Interior,
    Boundary,
    Outside,
}

/// A crossing whose triangle normal opposes the sampling direction enters the
/// solid; the opposite sign exits it. This is winding, not semantics.
fn facing_of(denominator: Scalar) -> SurfaceFacing {
    if denominator < 0.0 {
        SurfaceFacing::AgainstNormal
    } else {
        SurfaceFacing::WithNormal
    }
}

/// Edge-function containment scaled by the triangle's own size, so the
/// tolerance means the same thing for large and small triangles.
fn classify(triangle: &Triangle3, normal: Vec3, point: Point3, linear: Scalar) -> Containment {
    let length = normal.length();
    let scale = if length > 0.0 { length } else { 1.0 };
    let edges = [
        (triangle.b - triangle.a)
            .cross(point - triangle.a)
            .dot(normal)
            / scale,
        (triangle.c - triangle.b)
            .cross(point - triangle.b)
            .dot(normal)
            / scale,
        (triangle.a - triangle.c)
            .cross(point - triangle.c)
            .dot(normal)
            / scale,
    ];
    let limit = linear.max(Scalar::EPSILON) * scale.max(1.0);
    if edges.iter().any(|value| *value < -limit) {
        return Containment::Outside;
    }
    if edges.iter().any(|value| value.abs() <= limit) {
        return Containment::Boundary;
    }
    Containment::Interior
}
