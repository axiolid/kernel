//! The narrowest supported consumer: core values only.
//!
//! This profile exists to prove `axiolid-core` is independently usable.
//! Points, vectors, frames, tolerances, and intervals must not require a
//! representation, an algorithm, or a provider to be constructed.

use axiolid_core::{Interval, Point3, Tolerance, Vec3};

fn main() {
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(3.0, 4.0, 0.0);
    assert_eq!((b - a).length(), 5.0);

    let axis = Vec3::Z;
    assert_eq!(axis.dot(Vec3::X), 0.0);

    let span = Interval::new(0.0, 1.0);
    assert!(span.end > span.start);

    // A tolerance is a core value: comparing at a stated scale needs no
    // geometry package.
    let tol = Tolerance::MILLIMETRE;
    assert!(tol.linear() > 0.0);

    println!("core-only ok");
}
