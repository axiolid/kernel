//! The declared scratch bound holds for every measured boolean (#110, #226),
//! and a batch never exceeds a budget it admits (#259).
//!
//! `ScratchRequirement::fits_budget` admits work on the provider's word, so a
//! declared bound below a real peak makes a memory budget look enforced when
//! it is not. This measures the same workloads as the `scratch_probe` binary
//! and fails if any peak exceeds the declaration for the worker count it ran
//! on.
//!
//! With rayon the measurement forces the worst-case schedule instead of
//! waiting for machine load to produce it (#226): pools of 1 to 64 workers,
//! started inside the measured window, join halves pushed onto other workers.
//! An idle machine and a loaded one therefore measure the same costs, which
//! is what makes this deterministic.
//!
//! The batch paths have no bound in their input triangles, so they charge
//! each step and refuse the first that does not fit (`src/budget.rs`). For
//! each batch workload this finds the smallest budget it admits and checks
//! that the run at that budget, and the refused run a byte below it, stay
//! within it, under the same forced schedules, and that a budget below its
//! copies and bookkeeping is refused before they are made. On 16 and 64
//! workers every
//! pair of every `union_many` level can run at once, which is the most
//! booleans that can be live: stolen work can put more of them in flight
//! than there are workers, never more than there are pairs.
//!
//! No test harness (`harness = false`): the counting allocator is global, and
//! the harness's own threads would allocate while a boolean is measured.

#[path = "../src/bin/scratch_probe/batch.rs"]
mod batch;
#[path = "../src/bin/scratch_probe/scratch.rs"]
mod scratch;

use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;

fn main() {
    let requirement = BoolmeshBoolean::new().scratch_requirement();
    let mut measured = 0;
    let mut smallest_peak = Vec::new();
    for schedule in scratch::Schedule::all() {
        let workers = schedule.workers();
        let samples = scratch::measure(schedule);
        assert_eq!(samples.len(), 16, "four sizes, four operations");
        for sample in &samples {
            let declared = requirement
                .upper_bound_bytes_on(sample.elements, workers)
                .expect("boolmesh declares a finite scratch bound");
            assert!(
                sample.peak <= declared,
                "{:?} on {} triangles and {} workers peaked at {} bytes, above the declared {} \
                 ({:?})",
                sample.operation,
                sample.elements,
                workers,
                sample.peak,
                declared,
                requirement
            );
        }
        measured += samples.len();
        let smallest = samples[..4].iter().map(|s| s.peak).max().expect("samples");
        smallest_peak.push((workers, smallest));
    }

    // The forcing is itself checked, or a harness that stopped charging
    // worker start-up would pass any declaration that ignores workers. Each
    // rayon worker allocates several KB of its own; at least 2 KiB per
    // extra worker must show up in the smallest boolean's peak.
    if let (Some(&(1, one)), Some(&(widest, wide))) = (smallest_peak.first(), smallest_peak.last())
    {
        assert!(
            wide >= one + 2048 * (widest - 1),
            "24 triangles peaked at {wide} bytes on {widest} workers and {one} on one: \
             worker start-up is not inside the measured window"
        );
    }

    let mut batches = 0;
    for schedule in scratch::Schedule::all() {
        let workers = schedule.workers();
        for sample in batch::measure(schedule) {
            assert!(
                sample.peak <= sample.admitted,
                "{} ({} triangles) on {workers} workers peaked at {} bytes, above the {} it \
                 admitted",
                sample.name,
                sample.elements,
                sample.peak,
                sample.admitted
            );
            assert!(
                sample.refused_peak < sample.admitted,
                "{} on {workers} workers was refused {} bytes only after allocating {}",
                sample.name,
                sample.admitted - 1,
                sample.refused_peak
            );
            assert!(
                sample.small_peak <= batch::SMALL_BUDGET,
                "{} on {workers} workers was refused {} bytes only after allocating {}",
                sample.name,
                batch::SMALL_BUDGET,
                sample.small_peak
            );
            batches += 1;
        }
    }
    println!(
        "scratch_bound: {measured} booleans within {requirement:?}; {batches} batches within \
         the budgets they admit"
    );
}
