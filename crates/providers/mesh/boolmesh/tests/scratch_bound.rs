//! The declared scratch bound holds for every measured boolean (#110).
//!
//! `ScratchRequirement::fits_budget` admits work on the provider's word, so a
//! declared `PerElement` bound below a real peak makes a memory budget look
//! enforced when it is not. This measures the same workloads as the
//! `scratch_probe` binary and fails if any peak exceeds the declaration.
//!
//! No test harness (`harness = false`): the counting allocator is global, and
//! the harness's own threads would allocate while a boolean is measured.

#[path = "../src/bin/scratch_probe/scratch.rs"]
mod scratch;

use axiolid_contracts::ScratchRequirement;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;

fn main() {
    let ScratchRequirement::PerElement { bytes_per_element } =
        BoolmeshBoolean::new().scratch_requirement()
    else {
        panic!("boolmesh declares a per-triangle scratch bound");
    };
    let samples = scratch::measure();
    assert_eq!(samples.len(), 16, "four sizes, four operations");
    for sample in &samples {
        let declared = bytes_per_element * sample.elements;
        assert!(
            sample.peak <= declared,
            "{:?} on {} triangles peaked at {} bytes, above the declared {} ({} per triangle)",
            sample.operation,
            sample.elements,
            sample.peak,
            declared,
            bytes_per_element
        );
    }
    println!("scratch_bound: 16 booleans within {bytes_per_element} bytes per triangle");
}
