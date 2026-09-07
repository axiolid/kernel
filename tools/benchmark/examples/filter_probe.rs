//! Which orient3d inputs actually escalate past the floating-point filter.
//!
//! Written because a benchmark row labelled "near-degenerate" measured
//! identically to the separated case: the input was not escalating at all, so
//! the row reported the filter's cost under a name that implied exact
//! arithmetic. Guessing at an offset twice produced two wrong answers, so this
//! asks the filter directly.

use axiolid_core::Point3;
use axiolid_guarantees::Certified;
use axiolid_predicates::orient3d_filter;

fn main() {
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1.0, 0.0, 0.0);
    let c = Point3::new(0.0, 1.0, 0.0);

    println!("{:>12}  {:>10}  outcome", "z-offset", "certain?");
    for exponent in [0i32, -1, -8, -14, -15, -16, -17, -18, -20, -30, -300] {
        let z = if exponent == 0 {
            0.0
        } else {
            10f64.powi(exponent)
        };
        let d = Point3::new(0.5, 0.5, z);
        let outcome = orient3d_filter(a, b, c, d);
        let certain = matches!(outcome, Certified::Certain { .. });
        println!("{z:>12e}  {certain:>10}  {outcome:?}");
    }
}
