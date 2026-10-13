//! Memory-budget checks for one boolean and for the batch paths (#226, #259).
//!
//! [`MeshBoolean::scratch_requirement`] bounds one boolean in the input
//! triangles of its two operands and the workers it runs on. A batch is a
//! sequence of such booleans, and between them it holds state the
//! declaration does not see:
//!
//! - Intermediates. `union_many` feeds each level's results into the next,
//!   and `subtract_many` cuts a subject that grows with every group. A
//!   union of overlapping solids has more triangles than its operands, so no
//!   bound in the batch's *input* triangles holds for the later steps.
//! - Concurrency. With `parallel-batch` every pair of a `union_many` level
//!   can be in flight at once. With `parallel` as well, a worker waiting in
//!   one boolean's `rayon::join` steals another pair's boolean and runs it on
//!   top, so more booleans are live than there are workers: up to every
//!   pair of the level, not one per worker.
//!
//! Neither is known before the batch starts, so the batch paths charge each
//! step when its operands exist and refuse with
//! [`GeomError::BudgetExceeded`] before the step whose scratch would not fit.
//! A step's charge is the bytes the batch already holds (its own copies and
//! intermediates, by [`heap_bytes`]) plus the declaration applied to every
//! boolean that can be live during it ([`concurrent_bound`]). A finished
//! boolean's result stays live as an intermediate; it is covered by that
//! boolean's term, because the measured peak of one boolean includes its
//! result. `tests/scratch_bound.rs` measures both batch paths under forced
//! worst-case schedules and checks the peak against the smallest budget
//! each admits.
//!
//! [`MeshBoolean::scratch_requirement`]: axiolid_mesh_boolean_contract::MeshBoolean::scratch_requirement

use std::mem::size_of;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_mesh::{AttributeChannel, AttributeFate, TriMesh};

use crate::provider::{SCRATCH_BASE_BYTES, SCRATCH_BYTES_PER_TRIANGLE, SCRATCH_BYTES_PER_WORKER};

/// Worker threads this call can run on: the rayon pool it is in, or the
/// calling thread alone in a build without rayon.
pub(crate) fn running_workers() -> usize {
    #[cfg(any(feature = "parallel", feature = "parallel-batch"))]
    {
        rayon::current_num_threads()
    }
    #[cfg(not(any(feature = "parallel", feature = "parallel-batch")))]
    {
        1
    }
}

/// The declared scratch of `booleans` booleans live at once over `elements`
/// input triangles in total, on `workers` threads, or `None` on overflow.
///
/// Each boolean pays the fixed term and its triangles; the per-worker term
/// is paid once, because a worker's bookkeeping exists once however many
/// booleans are suspended on its stack, and the work stolen onto it is the
/// one it is running. No boolean, no scratch.
pub(crate) fn concurrent_bound(booleans: usize, elements: usize, workers: usize) -> Option<usize> {
    if booleans == 0 {
        return Some(0);
    }
    let base = SCRATCH_BASE_BYTES.checked_mul(booleans)?;
    let per_elements = SCRATCH_BYTES_PER_TRIANGLE.checked_mul(elements)?;
    let per_workers = SCRATCH_BYTES_PER_WORKER.checked_mul(workers)?;
    base.checked_add(per_elements)?.checked_add(per_workers)
}

/// Refuse unless `held` bytes of batch state plus `booleans` concurrent
/// booleans over `elements` triangles fit `options`' memory budget.
///
/// Without a budget nothing is computed, so an unbudgeted batch does not
/// start rayon's global pool early.
pub(crate) fn charge_step(
    options: &ExecutionOptions,
    held: usize,
    booleans: usize,
    elements: usize,
) -> GeomResult<()> {
    let Some(budget) = options.memory_budget_bytes() else {
        return Ok(());
    };
    match concurrent_bound(booleans, elements, running_workers())
        .and_then(|bound| bound.checked_add(held))
    {
        Some(needed) if needed <= budget => Ok(()),
        _ => Err(GeomError::BudgetExceeded { resource: "memory" }),
    }
}

/// Heap bytes `mesh` owns, by capacity.
pub(crate) fn heap_bytes(mesh: &TriMesh) -> usize {
    let normals = mesh.normals.as_ref().map_or(0, |normals| {
        normals.values.capacity() * size_of::<Vec3>()
            + normals.indices.as_ref().map_or(0, Vec::capacity) * size_of::<u32>()
    });
    let channels: usize = mesh
        .attributes
        .iter()
        .map(|channel| {
            channel.name.capacity()
                + channel.values.capacity() * size_of::<Scalar>()
                + channel.corner_indices.as_ref().map_or(0, Vec::capacity) * size_of::<u32>()
        })
        .sum();
    mesh.positions.capacity() * size_of::<Point3>()
        + mesh.indices.capacity() * size_of::<u32>()
        + normals
        + mesh.attributes.capacity() * size_of::<AttributeChannel>()
        + channels
}

/// Heap bytes of one node's attribute history.
pub(crate) fn fates_bytes(fates: &[(String, AttributeFate)]) -> usize {
    fates
        .iter()
        .map(|(name, _)| size_of::<(String, AttributeFate)>() + name.capacity())
        .sum()
}

/// Upper bound on the heap bytes of
/// [`fuse_with_channels`](crate::grouping::fuse_with_channels)`(members,
/// template)` plus the normals of a single member, before building it.
///
/// Positions and indices are sized exactly; the channel vectors grow by
/// pushing, so each is charged twice its length plus the smallest
/// allocation, which is what doubling growth can reach.
pub(crate) fn fused_bound(members: &[&TriMesh], template: &[AttributeChannel]) -> usize {
    let positions: usize = members.iter().map(|m| m.positions.len()).sum();
    let corners: usize = members.iter().map(|m| m.indices.len()).sum();
    let normals: usize = match members {
        [only] => only.normals.as_ref().map_or(0, |normals| {
            normals.values.len() * size_of::<Vec3>()
                + normals.indices.as_ref().map_or(0, Vec::len) * size_of::<u32>()
        }),
        _ => 0,
    };
    let channels: usize = template
        .iter()
        .map(|want| {
            want.name.len()
                + grown(corners * want.width) * size_of::<Scalar>()
                + corners * size_of::<u32>()
        })
        .sum();
    positions * size_of::<Point3>()
        + corners * size_of::<u32>()
        + normals
        + grown(template.len()) * size_of::<AttributeChannel>()
        + channels
}

/// Capacity a vector grown one push at a time can reach holding `len`.
fn grown(len: usize) -> usize {
    if len == 0 {
        0
    } else {
        2 * len + 8
    }
}
