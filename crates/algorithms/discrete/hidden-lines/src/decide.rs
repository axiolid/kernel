//! Sign decisions: the interval filter first, the exact tier when it cannot
//! decide. Nothing is ever read as zero because it was undecided.

use std::cmp::Ordering;

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_exact::{Dyadic, Interval, SignExpr};
use axiolid_guarantees::Sign;

/// The sign of `expr` as an [`Ordering`] against zero.
///
/// `escalations` counts the decisions the interval filter could not settle.
/// A sign the exact tier cannot produce is a refusal.
pub(crate) fn decide<E: SignExpr>(expr: &E, escalations: &mut usize) -> GeomResult<Ordering> {
    if let Some(sign) = expr.sign_in::<Interval>() {
        return ordering(sign);
    }
    *escalations += 1;
    match expr.sign_in::<Dyadic>() {
        Some(sign) => ordering(sign),
        None => Err(GeomError::Degenerate(
            "a visibility sign could not be decided exactly".into(),
        )),
    }
}

fn ordering(sign: Sign) -> GeomResult<Ordering> {
    match sign {
        Sign::Negative => Ok(Ordering::Less),
        Sign::Zero => Ok(Ordering::Equal),
        Sign::Positive => Ok(Ordering::Greater),
        _ => Err(GeomError::Degenerate(
            "an unknown sign cannot decide visibility".into(),
        )),
    }
}
