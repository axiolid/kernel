//! Print `boolmesh`'s measured peak scratch against its declared bound.
//!
//! ```text
//! cargo run --release -p axiolid-mesh-boolean-boolmesh --bin scratch_probe --all-features
//! ```
//!
//! Each row of the first table is the worst of the four operations at one
//! size and worker count, measured under the forced worst-case schedule of
//! `scratch.rs`. The second table measures the batch paths against the
//! smallest budget each admits (`batch.rs`, #259). The third measures
//! booleans whose result outgrows their operands (`crossing.rs`, #308):
//! with no budget, peak against the operands' declared bound; with one, the
//! smallest budget each admits and what a budget of the declared bound does.
//!
//! `scratch_probe combs` prints only the crossing workloads with no budget,
//! the table of #308; `scratch_probe crossing` only the third table.

mod batch;
mod crossing;
mod scratch;
mod shapes;

use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;

fn main() {
    match std::env::args().nth(1).as_deref() {
        Some("combs") => return crossing_unbudgeted(),
        Some("crossing") => return crossing_budgeted(),
        _ => {}
    }
    let requirement = BoolmeshBoolean::new().scratch_requirement();
    println!("declared: {requirement:?}");
    println!(
        "{:>8}  {:>10}  {:>12}  {:>12}  {:>6}  worst operation",
        "workers", "triangles", "peak bytes", "declared", "used"
    );
    let mut worst_used = 0.0f64;
    for schedule in scratch::Schedule::all() {
        let workers = schedule.workers();
        for size in scratch::measure(schedule).chunks(4) {
            let worst = size.iter().max_by_key(|s| s.peak).expect("four samples");
            let declared = requirement
                .upper_bound_bytes_on(worst.elements, workers)
                .expect("boolmesh declares a finite bound");
            let used = worst.peak as f64 / declared as f64;
            worst_used = worst_used.max(used);
            println!(
                "{:>8}  {:>10}  {:>12}  {:>12}  {:>5.0}%  {:?}",
                workers,
                worst.elements,
                worst.peak,
                declared,
                used * 100.0,
                worst.operation
            );
        }
    }
    println!();
    println!("worst peak / declared bound: {:.0}%", worst_used * 100.0);
    println!("Above 100% the declaration is not an upper bound: raise it, never the other way.");

    println!();
    println!(
        "{:>8}  {:>10}  {:>12}  {:>12}  {:>6}  {:>12}  {:>6}  batch",
        "workers", "triangles", "peak bytes", "admitted", "used", "refused peak", "small"
    );
    let mut worst_batch = 0.0f64;
    for schedule in scratch::Schedule::all() {
        for sample in batch::measure(schedule) {
            let used = sample.peak as f64 / sample.admitted as f64;
            worst_batch = worst_batch.max(used);
            println!(
                "{:>8}  {:>10}  {:>12}  {:>12}  {:>5.0}%  {:>12}  {:>6}  {}",
                schedule.workers(),
                sample.elements,
                sample.peak,
                sample.admitted,
                used * 100.0,
                sample.refused_peak,
                sample.small_peak,
                sample.name
            );
        }
    }
    println!();
    println!(
        "worst batch peak / smallest admitted budget: {:.0}%",
        worst_batch * 100.0
    );
    println!();
    crossing_budgeted();
}

/// Booleans whose result outgrows their operands, against the smallest
/// budget each admits and against their operands' declared bound (#308).
fn crossing_budgeted() {
    println!(
        "{:>8}  {:>10}  {:>12}  {:>12}  {:>6}  {:>12}  {:>12}  {:>12}  workload",
        "workers",
        "triangles",
        "peak bytes",
        "admitted",
        "used",
        "refused peak",
        "declared",
        "at declared"
    );
    let workloads = crossing::workloads();
    let mut worst = 0.0f64;
    for schedule in scratch::Schedule::all() {
        for sample in crossing::budgeted(&workloads, schedule) {
            let used = sample.peak as f64 / sample.admitted as f64;
            worst = worst.max(used);
            let at_declared = if sample.declared_refused {
                format!("refused {}", sample.declared_peak)
            } else {
                format!("{}", sample.declared_peak)
            };
            println!(
                "{:>8}  {:>10}  {:>12}  {:>12}  {:>5.0}%  {:>12}  {:>12}  {:>12}  {}",
                schedule.workers(),
                sample.elements,
                sample.peak,
                sample.admitted,
                used * 100.0,
                sample.refused_peak,
                sample.declared,
                at_declared,
                sample.name
            );
        }
    }
    println!();
    println!(
        "worst crossing peak / smallest admitted budget: {:.0}%",
        worst * 100.0
    );
}

/// The table of #308: peak with no budget against the operands' declared
/// bound.
fn crossing_unbudgeted() {
    let requirement = BoolmeshBoolean::new().scratch_requirement();
    println!(
        "{:>8}  {:>10}  {:>10}  {:>12}  {:>12}  {:>7}  workload",
        "workers", "triangles", "output", "peak bytes", "declared", "ratio"
    );
    let workloads = crossing::workloads();
    for schedule in scratch::Schedule::all() {
        for sample in crossing::unbudgeted(&workloads, schedule) {
            let declared = requirement
                .upper_bound_bytes_on(sample.elements, schedule.workers())
                .expect("boolmesh declares a finite bound");
            println!(
                "{:>8}  {:>10}  {:>10}  {:>12}  {:>12}  {:>6.2}x  {}",
                schedule.workers(),
                sample.elements,
                sample.output,
                sample.peak,
                declared,
                sample.peak as f64 / declared as f64,
                sample.name
            );
        }
    }
}
