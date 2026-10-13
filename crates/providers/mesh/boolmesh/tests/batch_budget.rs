//! What the batch paths charge a memory budget, step by step (#259).
//!
//! `tests/scratch_bound.rs` measures that a batch stays within the budget it
//! admits. These pin what it admits: each step is charged the declared bound
//! of every boolean that can be live during it, plus the copies and
//! intermediates the batch holds, and is refused with `BudgetExceeded`
//! before it starts when that does not fit. A charge that dropped a term
//! would admit a budget below the one computed here.
//!
//! Disjoint boxes keep the arithmetic exact: their unions concatenate, so
//! the triangle counts of every level are known in advance.

mod support;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, ScratchRequirement};
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::{BooleanOutcome, MeshBoolean};
use support::boxx;

/// Workers the batch runs on: two with rayon, fewer than the pairs of a
/// level, so a charge per worker instead of per pair would show.
const WORKERS: usize = if cfg!(any(feature = "parallel", feature = "parallel-batch")) {
    2
} else {
    1
};

/// Run `call` on [`WORKERS`] threads.
fn on_workers<R: Send>(call: impl FnOnce() -> R + Send) -> R {
    #[cfg(any(feature = "parallel", feature = "parallel-batch"))]
    {
        rayon::ThreadPoolBuilder::new()
            .num_threads(WORKERS)
            .build()
            .expect("pool")
            .install(call)
    }
    #[cfg(not(any(feature = "parallel", feature = "parallel-batch")))]
    {
        call()
    }
}

/// The declared terms: fixed, per triangle, per worker.
fn terms() -> (usize, usize, usize) {
    match BoolmeshBoolean::new().scratch_requirement() {
        ScratchRequirement::Affine {
            base_bytes,
            bytes_per_element,
            bytes_per_worker,
        } => (base_bytes, bytes_per_element, bytes_per_worker),
        other => panic!("boolmesh declares an affine bound, not {other:?}"),
    }
}

/// The declared bound of one boolean over `triangles` on [`WORKERS`].
fn one_boolean(triangles: usize) -> usize {
    let (base, per_triangle, per_worker) = terms();
    base + per_triangle * triangles + per_worker * WORKERS
}

/// `n` disjoint unit boxes of 12 triangles each.
fn disjoint_boxes(n: usize) -> Vec<TriMesh> {
    (0..n)
        .map(|i| boxx(2.0 * i as f64, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0))
        .collect()
}

fn budget(bytes: usize) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE).with_memory_budget(bytes)
}

fn refused(outcome: GeomResult<BooleanOutcome>, what: &str) {
    match outcome {
        Err(GeomError::BudgetExceeded { resource: "memory" }) => {}
        other => panic!(
            "{what} must be refused for memory, got {:?}",
            other.map(drop)
        ),
    }
}

fn union_many(solids: &[TriMesh], bytes: usize) -> GeomResult<BooleanOutcome> {
    on_workers(|| BoolmeshBoolean::new().union_many(solids, &budget(bytes)))
}

/// Bookkeeping a batch holds beyond its meshes: node and index lists, a few
/// hundred bytes per solid, generously.
const BOOKKEEPING: usize = 16 * 1024;

/// With `parallel-batch` every pair of a level is charged at once: the
/// eight pairs of the first level run concurrently, and stolen work can put
/// all eight in flight on two workers. That, and the batch's bookkeeping,
/// is the smallest budget it admits; every later level has fewer pairs.
#[cfg(feature = "parallel-batch")]
#[test]
fn union_many_charges_every_pair_of_a_level_at_once() {
    let solids = disjoint_boxes(16);
    let (base, _, _) = terms();
    let first_level = one_boolean(16 * 12) + 7 * base;
    refused(
        union_many(&solids, first_level - 1),
        "a budget one byte short of eight concurrent booleans",
    );
    union_many(&solids, first_level + BOOKKEEPING).expect("eight concurrent booleans fit");
}

/// Without `parallel-batch` the pairs run one at a time, but the level and
/// the results so far stay live: the last union, of two halves of 96
/// triangles each, is charged their copies on top of its own bound, so the
/// bound alone is refused. The copies are a few KB, not another boolean.
#[cfg(not(feature = "parallel-batch"))]
#[test]
fn union_many_charges_the_intermediates_it_holds() {
    let solids = disjoint_boxes(16);
    let last_union = one_boolean(16 * 12);
    refused(
        union_many(&solids, last_union),
        "the last union's own bound, without the halves it holds",
    );
    union_many(&solids, last_union + 2 * BOOKKEEPING).expect("the halves are a few KB");
}

/// One solid is no boolean: only its bookkeeping is charged, not a
/// boolean's terms (with rayon, not a worker's either).
#[test]
fn union_many_of_one_solid_runs_no_boolean() {
    union_many(&disjoint_boxes(1), 1024).expect("no boolean runs");
}

/// `subtract_many` fuses disjoint tools into one, and that copy is held
/// through the boolean: eight tools of eight vertices and 36 corners make
/// 8 * (8 * 24 + 36 * 4) = 2,688 bytes beyond the boolean's own bound and
/// the batch's bookkeeping.
#[test]
fn subtract_many_charges_the_fused_tool() {
    let subject = boxx(8.0, 0.5, 0.0, 20.0, 3.0, 1.0, 0.0);
    let tools: Vec<TriMesh> = (0..8)
        .map(|i| boxx(2.0 * i as f64 + 1.0, 0.5, -0.5, 1.0, 1.0, 2.0, 0.0))
        .collect();
    let subtract = |bytes| {
        on_workers(|| BoolmeshBoolean::new().subtract_many(&subject, &tools, &budget(bytes)))
    };
    let step = one_boolean(12 + 8 * 12);
    refused(subtract(step + 2_687), "the boolean without its fused tool");
    let outcome = subtract(step + 2_688 + BOOKKEEPING).expect("the boolean and its fused tool");
    assert_eq!(outcome.evidence.sub_operations, 1, "one fused group");
}

/// Split every triangle into four, `levels` times.
fn subdivide(mesh: &TriMesh, levels: usize) -> TriMesh {
    let mut mesh = mesh.clone();
    for _ in 0..levels {
        let mut positions = mesh.positions.clone();
        let mut indices = Vec::with_capacity(mesh.indices.len() * 4);
        for corner in mesh.indices.chunks_exact(3) {
            let [a, b, c] = [corner[0], corner[1], corner[2]];
            let mid =
                |i: u32, j: u32| (mesh.positions[i as usize] + mesh.positions[j as usize]) * 0.5;
            let next = positions.len() as u32;
            positions.extend([mid(a, b), mid(b, c), mid(c, a)]);
            let (ab, bc, ca) = (next, next + 1, next + 2);
            indices.extend([a, ab, ca, ab, b, bc, ca, bc, c, ab, bc, ca]);
        }
        mesh = TriMesh::new(positions, indices);
    }
    mesh
}

/// After the first group the subject is an intermediate the batch holds,
/// larger than the input it came from (#259). Two groups: four notches,
/// then four tools overlapping them. The second step is charged the cut
/// subject on top of its own bound and the fused tools, so that bound, the
/// fused tools and the bookkeeping alone are refused. The subject is
/// subdivided so the cut copy is clearly larger than the bookkeeping.
#[test]
fn subtract_many_charges_the_subject_it_has_cut() {
    let subject = subdivide(&boxx(8.0, 0.5, 0.0, 20.0, 3.0, 1.0, 0.0), 3);
    let notch = |x: f64| boxx(x, 0.5, -0.5, 1.0, 1.0, 2.0, 0.0);
    let first: Vec<TriMesh> = (0..4).map(|i| notch(4.0 * i as f64 + 1.0)).collect();
    let second: Vec<TriMesh> = (0..4).map(|i| notch(4.0 * i as f64 + 1.5)).collect();
    // Interleaved, so first-fit grouping puts each set in a group of its own.
    let tools: Vec<TriMesh> = first
        .iter()
        .zip(&second)
        .flat_map(|(a, b)| [a.clone(), b.clone()])
        .collect();
    let subtract = |tools: &[TriMesh], bytes| {
        on_workers(|| BoolmeshBoolean::new().subtract_many(&subject, tools, &budget(bytes)))
    };
    let cut = subtract(&first, usize::MAX).expect("first group").mesh;
    assert!(
        cut.triangle_count() > subject.triangle_count(),
        "the cut subject grew"
    );
    // Four tools of eight vertices and 36 corners, fused.
    let fused = 4 * (8 * 24 + 36 * 4);
    let second_step = one_boolean(cut.triangle_count() + 4 * 12) + fused;
    refused(
        subtract(&tools, second_step + BOOKKEEPING),
        "the second step without the cut subject it holds",
    );
    // The batch holds exactly the mesh the first group returns.
    let held = cut.positions.capacity() * std::mem::size_of_val(&cut.positions[0])
        + cut.indices.capacity() * std::mem::size_of::<u32>();
    assert!(
        held > BOOKKEEPING,
        "the cut subject outweighs the bookkeeping"
    );
    let outcome =
        subtract(&tools, second_step + held + BOOKKEEPING).expect("the step and the cut subject");
    assert_eq!(outcome.evidence.sub_operations, 2, "two groups");
}

/// A budget given to the provider directly, without dispatch, is enforced in
/// every build, not only with `parallel` (#259).
#[test]
fn a_direct_boolean_honours_its_budget() {
    let subject = boxx(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0);
    let tool = boxx(0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 0.0);
    let run = |bytes| {
        on_workers(|| {
            BoolmeshBoolean::new().boolean(&subject, &tool, BooleanOperator::Union, &budget(bytes))
        })
    };
    let bound = one_boolean(24);
    refused(run(bound - 1), "one byte under the declared bound");
    run(bound).expect("the declared bound");
}
