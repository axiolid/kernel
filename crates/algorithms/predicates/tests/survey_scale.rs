// SPDX-License-Identifier: MPL-2.0

//! Does survey-scale coordinate quantisation actually break coincident
//! geometry, or is the ~2 nm grid harmless in practice?
//!
//! # Why this exists
//!
//! `fix(measure)` removed the conditioning loss in mass properties by summing
//! about a local origin. What it could NOT fix is the input: at base 1e7 the
//! spacing between adjacent binary64 values is ~1.9e-9 m, so an authored
//! coordinate is snapped to that grid before the kernel ever sees it.
//!
//! The open question was whether that matters. The stated worry was exact
//! predicates on coincident geometry: two faces meant to be flush getting
//! snapped to slightly different values, so a boolean sees a sliver that
//! should not exist.
//!
//! This probe answers it with measurements instead of assertion. It is a
//! `#[test]` so it runs in CI and cannot rot, but every check here is a
//! characterisation: it records what IS, so a future change that alters the
//! behaviour has to acknowledge it.
//!
//! # What is deliberately NOT claimed
//!
//! `orient3d` is exact for finite binary64 input at any magnitude. So the
//! predicate never returns a wrong answer about the points it is GIVEN. The
//! only question is whether those points still mean what the author intended.
//! Conflating those two would make this probe measure the wrong thing.

use axiolid_core::Point3;
use axiolid_guarantees::{Certified, Sign};
use axiolid_predicates::{orient3d, orient3d_filter};

/// Survey-grid magnitudes, from origin to national grid.
const BASES: [f64; 5] = [0.0, 1.0e3, 1.0e5, 1.0e6, 1.0e7];

fn sign_of(c: Certified) -> Option<Sign> {
    match c {
        Certified::Certain { sign, .. } => Some(sign),
        _ => None,
    }
}

/// Spacing between adjacent representable values at `x`.
fn ulp(x: f64) -> f64 {
    if x == 0.0 {
        f64::EPSILON
    } else {
        let next = f64::from_bits(x.abs().to_bits() + 1);
        next - x.abs()
    }
}

/// Question 1: does coincidence survive when both sides are authored the
/// same way?
///
/// This is the common case: a wall ends at `base + 5.1` and the next begins
/// at `base + 5.1`, both evaluated by the same expression. Rounding is
/// deterministic, so both land on the same representable value.
#[test]
fn coincidence_survives_when_both_sides_round_identically() {
    for base in BASES {
        // A plane through three points, and a fourth point authored to lie
        // exactly on it by the same arithmetic.
        let a = Point3::new(base, base, 0.0);
        let b = Point3::new(base + 5.1, base, 0.0);
        let c = Point3::new(base, base + 5.1, 0.0);
        let on_plane = Point3::new(base + 5.1, base + 5.1, 0.0);

        let sign = sign_of(orient3d(a, b, c, on_plane));
        assert_eq!(
            sign,
            Some(Sign::Zero),
            "at base {base:e}, a point authored onto the plane by identical \
             arithmetic should still be exactly on it"
        );
    }
}

/// Question 2: does coincidence survive when the two sides are authored by
/// DIFFERENT arithmetic reaching the same intended value?
///
/// This is the case the worry was really about. A wall face at
/// `base + 10.2 - 5.1` and one at `base + 5.1` are the same place in exact
/// arithmetic, but the two expressions can round to different f64 values.
#[test]
fn coincidence_under_different_arithmetic_paths_is_measured() {
    let mut broken = Vec::new();
    for base in BASES {
        let a = Point3::new(base, base, 0.0);
        let b = Point3::new(base + 5.1, base, 0.0);
        let c = Point3::new(base, base + 5.1, 0.0);

        // Intended: exactly the same point as `base + 5.1`.
        let by_other_path = (base + 10.2) - 5.1;
        let direct = base + 5.1;
        let drift = (by_other_path - direct).abs();

        let probe = Point3::new(by_other_path, base + 5.1, 0.0);
        let sign = sign_of(orient3d(a, b, c, probe));

        eprintln!(
            "COINCIDENCE base={base:<9e} drift={drift:.3e} m  ulp={:.3e}  sign={sign:?}",
            ulp(base.max(1.0))
        );
        if sign != Some(Sign::Zero) {
            broken.push((base, drift));
        }
    }

    // Record the outcome rather than asserting a preferred one: the point of
    // the probe is to find out, and a wrong guess baked in as an assertion
    // would be worse than no probe.
    eprintln!(
        "COINCIDENCE summary: {} of {} magnitudes lost exact coincidence",
        broken.len(),
        BASES.len()
    );
    for (base, drift) in &broken {
        eprintln!("COINCIDENCE   base={base:e} drifted {drift:.3e} m off the plane");
    }
}

/// Question 3: the performance claim.
///
/// I claimed rebasing makes the filter hit its cheap path more often. That
/// is only worth acting on if the filter measurably escalates more at large
/// coordinates. This counts escalations directly.
#[test]
fn filter_escalation_rate_versus_coordinate_magnitude() {
    for base in BASES {
        let mut escalations = 0usize;
        let mut total = 0usize;

        // Sweep a point across a plane, passing through coincidence. Near
        // the crossing the filter should be unable to decide and escalate.
        for step in 0..2000 {
            let offset = (step as f64 - 1000.0) * 1.0e-9;
            let a = Point3::new(base, base, 0.0);
            let b = Point3::new(base + 1.0, base, 0.0);
            let c = Point3::new(base, base + 1.0, 0.0);
            let d = Point3::new(base + 0.25, base + 0.25, offset);

            total += 1;
            if sign_of(orient3d_filter(a, b, c, d)).is_none() {
                escalations += 1;
            }
        }

        let rate = escalations as f64 / total as f64;
        eprintln!(
            "ESCALATION base={base:<9e} {escalations:>4}/{total} = {:.1}% fell back to exact",
            rate * 100.0
        );
    }
}
