//! Topology-aware integration of certified affine surface intersections.
//!
//! The supported arrangement splits one rectangular patch whose certified trace
//! is a boundary-to-boundary chord. The same edge is retained as an explicit
//! embedded pcurve on the other patch when both endpoints are strictly interior
//! there. This avoids pretending that an interior segment partitions a disk.

use axiolid_contracts::GeomResult;
use axiolid_nurbs::{intersect_surface_surface_certified, CertifiedSurfaceSurfaceIntersection3};
use axiolid_surface::BSplineSurface;

pub use crate::trimmed_intersection_types::{
    CertifiedDualTrimmedSurfacePair3, CertifiedSurfacePairSplit3, CertifiedSurfacePairSplitOptions,
    CertifiedTrimmedSurfacePair3, EmbeddedFaceCurve, SurfacePairMember,
    SurfacePairSplitUnresolvedReason,
};

use crate::trimmed_intersection_assembly::{assemble, assemble_dual};
use crate::trimmed_intersection_classify::{classify, Classification};

/// Intersect two supported affine patches and construct a topology-aware trimmed arrangement.
///
/// `Split` is returned only for a single certified trace whose endpoints form a
/// boundary-to-boundary chord on exactly one patch and are both interior to the
/// other patch. The returned analytic B-rep has two closed trimmed faces for the
/// split patch and one closed rectangular face for the other patch. The wrapper's
/// embedded-curve relation attaches the same intersection edge to the unsplit
/// face without inventing a dangling trim loop.
///
/// Empty intersections remain `Empty`. Incomplete tracing, dual ownership,
/// corners, mixed ownership, coincidences, and residual-policy misses remain
/// `Unresolved` with the original intersection evidence.
pub fn split_surface_pair_certified(
    first: &BSplineSurface,
    second: &BSplineSurface,
    options: CertifiedSurfacePairSplitOptions,
) -> GeomResult<CertifiedSurfacePairSplit3> {
    let intersection =
        intersect_surface_surface_certified(first, second, options.intersection_options())?;
    let (mut traces, visited_patch_pairs, boundary_queries) = match intersection {
        CertifiedSurfaceSurfaceIntersection3::Unresolved { .. } => {
            return Ok(CertifiedSurfacePairSplit3::Unresolved {
                intersection,
                reason: SurfacePairSplitUnresolvedReason::IntersectionUnresolved,
            });
        }
        CertifiedSurfaceSurfaceIntersection3::Complete {
            traces,
            visited_patch_pairs,
            boundary_queries,
        } => (traces, visited_patch_pairs, boundary_queries),
        other => {
            return Ok(CertifiedSurfacePairSplit3::Unresolved {
                intersection: other,
                reason: SurfacePairSplitUnresolvedReason::IntersectionUnresolved,
            });
        }
    };
    if traces.is_empty() {
        return Ok(CertifiedSurfacePairSplit3::Empty {
            visited_patch_pairs,
            boundary_queries: u32::from(boundary_queries),
        });
    }
    if traces.len() != 1 {
        return Ok(unresolved_complete(
            traces,
            visited_patch_pairs,
            boundary_queries,
            SurfacePairSplitUnresolvedReason::UnsupportedTraceCount,
        ));
    }
    let Some(trace_ref) = traces.first() else {
        return Ok(CertifiedSurfacePairSplit3::Empty {
            visited_patch_pairs,
            boundary_queries: u32::from(boundary_queries),
        });
    };
    let residual_upper_bound = trace_ref
        .start
        .residual_upper_bound
        .max(trace_ref.end.residual_upper_bound);
    if !residual_upper_bound.is_finite() || residual_upper_bound > options.max_surface_residual() {
        return Ok(unresolved_complete(
            traces,
            visited_patch_pairs,
            boundary_queries,
            SurfacePairSplitUnresolvedReason::ResidualExceedsPolicy,
        ));
    }
    // The classifier distinguishes a proven terminal refusal from an
    // unimplemented combination; pass its verdict through unchanged rather
    // than collapsing both to one reason.
    let classification = match classify(first, second, trace_ref) {
        Ok(classification) => classification,
        Err(reason) => {
            return Ok(unresolved_complete(
                traces,
                visited_patch_pairs,
                boundary_queries,
                reason,
            ));
        }
    };
    let Some(trace) = traces.pop() else {
        return Ok(CertifiedSurfacePairSplit3::Empty {
            visited_patch_pairs,
            boundary_queries: u32::from(boundary_queries),
        });
    };
    // Both shapes assemble from the same certified trace; they differ only
    // in how many faces the chord bounds.
    match classification {
        Classification::Single(single) => {
            let split = assemble(
                first,
                second,
                trace,
                single,
                residual_upper_bound,
                visited_patch_pairs,
                boundary_queries,
            )?;
            Ok(CertifiedSurfacePairSplit3::Split(split))
        }
        Classification::Dual(dual) => {
            let split = assemble_dual(
                first,
                second,
                trace,
                dual,
                residual_upper_bound,
                visited_patch_pairs,
                boundary_queries,
            )?;
            Ok(CertifiedSurfacePairSplit3::DualSplit(split))
        }
    }
}

fn unresolved_complete(
    traces: Vec<axiolid_nurbs::TransverseSurfaceSurfaceTrace3>,
    visited_patch_pairs: u32,
    boundary_queries: u8,
    reason: SurfacePairSplitUnresolvedReason,
) -> CertifiedSurfacePairSplit3 {
    CertifiedSurfacePairSplit3::Unresolved {
        intersection: CertifiedSurfaceSurfaceIntersection3::Complete {
            traces,
            visited_patch_pairs,
            boundary_queries,
        },
        reason,
    }
}
