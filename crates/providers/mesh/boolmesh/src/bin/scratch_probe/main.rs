//! Print `boolmesh`'s measured peak scratch against its declared bound.
//!
//! ```text
//! cargo run --release -p axiolid-mesh-boolean-boolmesh --bin scratch_probe --all-features
//! ```
//!
//! Each row of the first table is the worst of the four operations at one
//! size and worker count, measured under the forced worst-case schedule of
//! `scratch.rs`. The second table measures the batch paths against the
//! smallest budget each admits (`batch.rs`, #259).

mod batch;
mod scratch;

use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;

fn main() {
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
}
