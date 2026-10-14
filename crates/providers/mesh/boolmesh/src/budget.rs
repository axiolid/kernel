//! Memory-budget checks for one boolean and for the batch paths (#226, #259, #308).
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
//! boolean that can be live during it ([`concurrent_bound`]). Each boolean
//! then runs with the budget the batch does not hold ([`remaining`]), or an
//! even share of it among the booleans of a `parallel-batch` level
//! (`concurrent_share`), and charges its own result against that (#308,
//! below). A finished boolean's result stays live as an intermediate; it is
//! covered by that boolean's share, because the measured peak of one
//! boolean includes its result. `tests/scratch_bound.rs` measures both
//! batch paths under forced worst-case schedules and checks the peak
//! against the smallest budget each admits.
//!
//! # The result is charged inside the solve (#308)
//!
//! The declaration counts input triangles, so it bounds a boolean whose
//! result is about the size of its operands. Operands that cross many times
//! have a far larger result: two perpendicular combs of 32 teeth, 1,048
//! triangles, cross 17,716 times, and their union peaked at 8.8x the
//! declaration (31x for the symmetric difference). No bound in the input
//! alone holds for that, and one in broad-phase candidate pairs is
//! quadratic, so the declaration stays the admission bound and the result
//! is charged once its size is known.
//!
//! [`ResultBudget`] carries what the budget leaves after the declaration.
//! The intersection pass counts crossings, and `boolean45` then computes
//! the result's half-edge count exactly, before allocating the result. The
//! declaration covers [`RESULT_HALF_EDGES_PER_TRIANGLE`] result half-edges
//! per input triangle; each one beyond that costs
//! [`SCRATCH_BYTES_PER_EXCESS_HALF_EDGE`]. A result that does not fit is
//! refused with [`GeomError::BudgetExceeded`] there. Every crossing puts at
//! least three half-edges into the result, so the intersection pass stops
//! as soon as its crossings alone could not fit, before its own tables grow
//! further.
//!
//! [`MeshBoolean::scratch_requirement`]: axiolid_mesh_boolean_contract::MeshBoolean::scratch_requirement

use std::mem::size_of;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_mesh::{AttributeChannel, AttributeFate, TriMesh};

use crate::provider::{SCRATCH_BASE_BYTES, SCRATCH_BYTES_PER_TRIANGLE, SCRATCH_BYTES_PER_WORKER};

/// Result half-edges per input triangle that the declared per-triangle term
/// covers (#308).
///
/// Measured with a counting allocator (`src/bin/scratch_probe/`, one
/// worker): over combs, crossing slab soups and pierced slabs, a boolean's
/// peak is within the declaration plus 245 bytes per result half-edge beyond
/// three per input triangle. Ordinary booleans stay below three: two
/// overlapping boxes subdivided three times reach 2.7 in their union, and
/// disjoint operands, whose faces pass through whole, exactly three. Those
/// are charged nothing beyond the declaration.
pub(crate) const RESULT_HALF_EDGES_PER_TRIANGLE: usize = 3;

/// Scratch per result half-edge beyond [`RESULT_HALF_EDGES_PER_TRIANGLE`]
/// per input triangle (#308): 1.8x the 245 bytes measured.
///
/// The half-edges are those of the result's polygonal faces before
/// triangulation, which `boolean45` counts exactly. The slope covers what
/// grows with them: the intersection tables (at most a third as many
/// crossings), the result's vertices and half-edges, its triangulation and
/// simplification, the returned mesh, and the channels carried onto it.
pub(crate) const SCRATCH_BYTES_PER_EXCESS_HALF_EDGE: usize = 448;

/// What one boolean may allocate for its result, once its operands are
/// admitted on the declaration (#308).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ResultBudget {
    /// Bytes the budget leaves after the declaration on the operands and
    /// the caller's held bytes.
    headroom: usize,
    /// Result half-edges the declaration already covers.
    covered: usize,
}

impl ResultBudget {
    /// Admit one boolean over `elements` input triangles on the pool it runs
    /// on, with `budget` bytes for it and nothing to spare for anything
    /// else. `None` without a budget.
    ///
    /// Refused with `BudgetExceeded` when the declaration does not fit: the
    /// check before dispatch charges the per-worker term for the width the
    /// options or the dispatcher name, and rayon runs the call on whatever
    /// pool it is in, which can be wider (#226).
    pub(crate) fn admit(budget: Option<usize>, elements: usize) -> GeomResult<Option<Self>> {
        let Some(budget) = budget else {
            return Ok(None);
        };
        match concurrent_bound(1, elements, running_workers()) {
            Some(declared) if declared <= budget => Ok(Some(Self {
                headroom: budget - declared,
                covered: elements.saturating_mul(RESULT_HALF_EDGES_PER_TRIANGLE),
            })),
            _ => Err(GeomError::BudgetExceeded { resource: "memory" }),
        }
    }

    /// Whether a result of `half_edges` polygon half-edges fits.
    pub(crate) fn fits(&self, half_edges: usize) -> bool {
        half_edges
            .saturating_sub(self.covered)
            .checked_mul(SCRATCH_BYTES_PER_EXCESS_HALF_EDGE)
            .is_some_and(|bytes| bytes <= self.headroom)
    }

    /// The most crossings whose result can still fit: each crossing adds at
    /// least three half-edges to the result (one to each face it cuts and
    /// one to the face it pierces), so past this [`Self::fits`] fails
    /// whatever else the result holds.
    pub(crate) fn crossing_limit(&self) -> usize {
        let half_edges = self
            .covered
            .saturating_add(self.headroom / SCRATCH_BYTES_PER_EXCESS_HALF_EDGE);
        half_edges / 3
    }
}

/// The budget left for a step after `held` bytes, or `None` without one.
pub(crate) fn remaining(options: &ExecutionOptions, held: usize) -> Option<usize> {
    options
        .memory_budget_bytes()
        .map(|budget| budget.saturating_sub(held))
}

/// The budget of one of `booleans` booleans that can run at once over
/// `elements` input triangles in all, `own` of them its operands', while the
/// batch holds `held` bytes (#308). `None` without a budget.
///
/// Each gets its own declaration and an equal share of what the step's
/// charge ([`charge_step`]) leaves, and charges its result against that.
/// Every share includes the per-worker term, which the step charges once:
/// the workers' bookkeeping exists once however many booleans run on them,
/// so what the booleans can allocate together still fits the budget.
#[cfg(feature = "parallel-batch")]
pub(crate) fn concurrent_share(
    options: &ExecutionOptions,
    held: usize,
    booleans: usize,
    elements: usize,
    own: usize,
) -> Option<usize> {
    let budget = options.memory_budget_bytes()?;
    let workers = running_workers();
    let step = concurrent_bound(booleans, elements, workers).and_then(|b| b.checked_add(held));
    let surplus = step.map_or(0, |step| budget.saturating_sub(step));
    // A step that does not fit was refused by `charge_step` already; a
    // share of nothing then admits nothing.
    let own = concurrent_bound(1, own, workers).unwrap_or(usize::MAX);
    Some(own.saturating_add(surplus / booleans.max(1)))
}

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

#[cfg(test)]
mod tests {
    use super::*;

    fn declared(elements: usize) -> usize {
        concurrent_bound(1, elements, running_workers()).expect("finite")
    }

    fn admitted(budget: usize, elements: usize) -> ResultBudget {
        ResultBudget::admit(Some(budget), elements)
            .expect("the declaration fits")
            .expect("a budget")
    }

    /// Admitted on the declaration; the result is free up to the allowance
    /// and costs its slope beyond it, to the byte (#308).
    #[test]
    fn a_result_is_charged_beyond_what_the_declaration_covers() {
        assert_eq!(ResultBudget::admit(None, 100), Ok(None));
        assert_eq!(
            ResultBudget::admit(Some(declared(100) - 1), 100),
            Err(GeomError::BudgetExceeded { resource: "memory" })
        );
        let exact = admitted(declared(100), 100);
        assert!(exact.fits(300), "three half-edges per triangle are covered");
        assert!(!exact.fits(301), "the next one is not");
        let room = admitted(declared(100) + 10 * SCRATCH_BYTES_PER_EXCESS_HALF_EDGE, 100);
        assert!(room.fits(310));
        assert!(!room.fits(311));
        // Three half-edges per crossing: 310 / 3.
        assert_eq!(room.crossing_limit(), 103);
        assert!(!room.fits(3 * (room.crossing_limit() + 1)));
    }

    /// Booleans running at once get shares that, less the per-worker term
    /// each repeats, add up to what the step leaves (#308).
    #[cfg(feature = "parallel-batch")]
    #[test]
    fn concurrent_shares_add_up_to_the_budget() {
        use axiolid_core::Tolerance;
        let workers = running_workers();
        let held = 1000;
        let step = concurrent_bound(3, 300, workers).expect("finite") + held;
        let options = ExecutionOptions::new(Tolerance::METRE).with_memory_budget(step + 3 * 4096);
        let shares: Vec<usize> = [50, 100, 150]
            .iter()
            .map(|&own| concurrent_share(&options, held, 3, 300, own).expect("a budget"))
            .collect();
        for (share, own) in shares.iter().zip([50, 100, 150]) {
            assert_eq!(*share, declared(own) + 4096, "its declaration and a third");
        }
        let repeated = 2 * SCRATCH_BYTES_PER_WORKER * workers;
        assert_eq!(
            shares.iter().sum::<usize>() - repeated + held,
            step + 3 * 4096,
            "the shares spend the budget exactly"
        );
    }
}
