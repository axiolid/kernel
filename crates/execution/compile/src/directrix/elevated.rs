//! Elevated and banked directrices (#252).
//!
//! A [`Curve3::Elevated`] (what a gradient curve lowers to) and a
//! [`Curve3::Banked`] are read by PLAN distance, their own parameter
//! (ADR 0082), over the plan's length: a sweep range on one is a span of
//! plan distance. A line plan has no end, so a sweep along it needs a
//! range. The curve is flattened by `flatten3` against its certified 3D
//! chord bound, the plan's and the profile's composed in quadrature
//! (`axiolid_reference::elevated`).
//!
//! It is ONE smooth curve -- end tangents reported, refined until a sweep's
//! walls fit -- unless the grade jumps at a seam of its profile (or of a
//! banked curve's pivot) inside the swept span by more than
//! [`CORNER_ANGLE`]. A chain plan is tangent-continuous at its joins by
//! construction. A curve with such a corner is swept as sampled and its
//! deviation is unbounded by name. A banked curve's point path is its
//! rotation point's: a disk swept along it is the same tube whatever the
//! roll, and the rolled section frames are what station-placed sections
//! stand in (`crate::station`).
//!
//! A straight plan under a profile of degree at most one is a straight 3D
//! segment, so it resolves as one for the exact swept disk and the
//! deviation report, like a line.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Scalar, Vec3};
use axiolid_curve::{Curve2, Curve3, ElevationLaw};

use super::ExactDirectrix;

/// Largest grade-angle jump (radians) at a seam that still counts as
/// tangent-continuous. Data authored tangent-continuous carries rounding
/// in its stored start grades; a real break is far larger. It only decides
/// whether the sweep refines and reports end tangents: the tube's
/// deviation is certified against the exact curve either way.
pub(crate) const CORNER_ANGLE: Scalar = 1e-9;

/// Whether `curve` is read here.
pub(crate) fn is_elevated(curve: &Curve3) -> bool {
    matches!(curve, Curve3::Elevated(_) | Curve3::Banked(_))
}

/// Whether the grade breaks strictly inside `[lo, hi]`. A seam whose laws
/// cannot be read there is a break: nothing is assumed smooth.
pub(crate) fn has_corner(curve: &Curve3, lo: Scalar, hi: Scalar) -> bool {
    let (lo, hi) = (lo.min(hi), lo.max(hi));
    axiolid_reference::elevated::grade_corners3(curve, CORNER_ANGLE)
        .map_or(true, |corners| corners.iter().any(|&c| c > lo && c < hi))
}

/// Unit tangents at both ends of `[start, end]`, in increasing plan
/// distance, when the span has no corner; `None` otherwise.
pub(crate) fn smooth_ends(curve: &Curve3, start: Scalar, end: Scalar) -> Option<[Vec3; 2]> {
    if has_corner(curve, start, end) {
        return None;
    }
    let unit = |d| {
        axiolid_reference::curve::derivative3(curve, d)
            .ok()
            .map(Vec3::normalize_or_zero)
            .filter(|v| *v != Vec3::ZERO && v.is_finite())
    };
    Some([unit(start)?, unit(end)?])
}

/// Refuse sampling an elevated curve with no end.
pub(crate) fn check_bounded(span: axiolid_core::Interval) -> GeomResult<()> {
    if span.start.is_finite() && span.end.is_finite() {
        Ok(())
    } else {
        Err(GeomError::InvalidInput(
            "an elevated directrix over a line plan has no end: give the sweep a \
             parameter range in plan distance"
                .into(),
        ))
    }
}

/// A straight 3D segment: a line plan under a profile of degree at most
/// one (a constant grade or a level), for the exact swept disk; `None`
/// for any other elevated curve.
pub(crate) fn exact_segment(
    curve: &Curve3,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
    unsupported: fn(&'static str) -> GeomError,
) -> Option<GeomResult<ExactDirectrix>> {
    let Curve3::Elevated(elevated) = curve else {
        return None;
    };
    let Curve2::Line(_) = elevated.plan.as_ref() else {
        return None;
    };
    let ElevationLaw::Polynomial { coefficients } = &elevated.elevation else {
        return None;
    };
    if coefficients.iter().skip(2).any(|c| *c != 0.0) {
        return None;
    }
    Some((|| {
        let range = range.ok_or_else(|| {
            unsupported("exact swept disk along an elevated curve over an unbounded line plan")
        })?;
        let (start, end) = super::finite_range(range)?;
        if start.min(end) < -options.tolerance().linear() {
            return Err(GeomError::InvalidInput(
                "sweep parameter range falls outside curve domain".into(),
            ));
        }
        let (start, end) = (start.max(0.0), end.max(0.0));
        Ok(ExactDirectrix::Segment(
            axiolid_reference::curve::evaluate3(curve, start)?,
            axiolid_reference::curve::evaluate3(curve, end)?,
        ))
    })())
}
