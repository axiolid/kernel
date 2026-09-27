//! Print `boolmesh`'s measured peak scratch per input triangle.
//!
//! ```text
//! cargo run --release -p axiolid-mesh-boolean-boolmesh --bin scratch_probe --all-features
//! ```

mod scratch;

fn main() {
    println!(
        "{:>10}  {:>14}  {:>18}",
        "triangles", "peak bytes", "bytes/triangle"
    );
    let mut worst_per_triangle = 0usize;
    for sample in scratch::measure() {
        let per_triangle = sample.peak / sample.elements.max(1);
        worst_per_triangle = worst_per_triangle.max(per_triangle);
        println!(
            "{:>10}  {:>14}  {:>18}  {:?}",
            sample.elements, sample.peak, per_triangle, sample.operation
        );
    }
    println!();
    println!("worst observed bytes/triangle: {worst_per_triangle}");
    println!("Declare PerElement with headroom above this, never below.");
}
