//! Print `boolmesh`'s measured peak scratch against its declared bound.
//!
//! ```text
//! cargo run --release -p axiolid-mesh-boolean-boolmesh --bin scratch_probe --all-features
//! ```
//!
//! Each row is the worst of the four operations at one size and worker
//! count, measured under the forced worst-case schedule of `scratch.rs`.

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
}
