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
//! A single boolean is admitted on its declaration and charges its result
//! inside the solve (#308). Booleans whose result outgrows their operands
//! (crossing combs of 2 to 32 teeth, crossing slab soups, pierced slabs)
//! must stay within the smallest budget each admits, be refused a byte
//! below it before exceeding it, and either fit their operands' declared
//! bound or be refused within it. The combs with more than two teeth must
//! be refused there: that is the case the declaration alone missed.
//!
//! No test harness (`harness = false`): the counting allocator is global, and
//! the harness's own threads would allocate while a boolean is measured.

#[path = "../src/bin/scratch_probe/batch.rs"]
mod batch;
// The table with no budget is the probe's alone.
#[allow(dead_code)]
#[path = "../src/bin/scratch_probe/crossing.rs"]
mod crossing;
#[path = "../src/bin/scratch_probe/scratch.rs"]
mod scratch;
#[path = "../src/bin/scratch_probe/shapes.rs"]
mod shapes;

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
    let mut crossings = 0;
    // Every workload on the calling thread (a build without rayon). The
    // forced schedule stalls every allocation, which makes the largest
    // results take minutes per pool, so on rayon pools this runs the combs
    // up to 8 teeth and a few soups and slabs; `scratch_probe crossing`
    // measures every workload on every pool (`src/budget.rs`).
    let all = crossing::workloads();
    let forced: Vec<crossing::Crossing> = crossing::workloads()
        .into_iter()
        .filter(|w| {
            [
                "combs of 2 ",
                "combs of 4 ",
                "combs of 8 ",
                "soup 0,",
                "soup 1,",
                "slab pierced by 4x4",
            ]
            .iter()
            .any(|prefix| w.name.starts_with(prefix))
        })
        .collect();
    for schedule in scratch::Schedule::all() {
        let workers = schedule.workers();
        let workloads = match schedule {
            scratch::Schedule::Inline => &all,
            scratch::Schedule::ColdPool { .. } => &forced,
        };
        for sample in crossing::budgeted(workloads, schedule) {
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
                sample.declared_peak <= sample.declared,
                "{} on {workers} workers {} its declared {} bytes after allocating {}",
                sample.name,
                if sample.declared_refused {
                    "was refused"
                } else {
                    "ran within"
                },
                sample.declared,
                sample.declared_peak
            );
            if sample.name.starts_with("combs of") && !sample.name.starts_with("combs of 2 ") {
                assert!(
                    sample.declared_refused,
                    "{} on {workers} workers fit its operands' declared bound",
                    sample.name
                );
            }
            crossings += 1;
        }
    }
    println!(
        "scratch_bound: {measured} booleans within {requirement:?}; {batches} batches within \
         the budgets they admit; {crossings} crossing booleans within the budgets they admit"
    );
}
