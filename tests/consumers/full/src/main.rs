//! The maximal supported closure: every facade feature enabled.
//!
//! This profile is the UPPER bound. The narrow profiles assert what a
//! focused application avoids compiling; this one records what the
//! project costs when nothing is avoided, so growth at the top is
//! visible rather than silent.

fn main() {
    // Touch a core value so the binary is not optimised into nothing.
    let p = axiolid::core::Point3::new(1.0, 2.0, 2.0);
    assert_eq!((p - axiolid::core::Point3::new(0.0, 0.0, 0.0)).length(), 3.0);
    println!("full ok");
}
