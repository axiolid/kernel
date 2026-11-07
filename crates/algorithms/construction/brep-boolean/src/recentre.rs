//! Exact recentring of operands far from the origin (#418, #419).
//!
//! A coordinate a million metres out has an ulp of a nanometre, and every
//! certified decision of the pipeline (a section's roots, a point's
//! distance to a face boundary, the sewing slacks) works at a scale of the
//! operands' own size, not of their distance from the origin. Placed far
//! out, the same solids lose most of their digits and the pipeline refuses
//! cases it decides at the origin.
//!
//! So operands far from the origin, relative to their own size, are moved
//! by `-c` first, cut there, and the result moved back by `+c`:
//!
//! - `c` is the centre of the operands' common vertex box, rounded on every
//!   axis to a multiple of the ulp of the largest coordinate, so that
//!   `x - c` is a multiple of `x`'s own ulp and, for a coordinate within
//!   the box, fits in 53 bits. An axis whose centre lies within four
//!   diagonals of zero keeps `c = 0` there, so a box far out in `x` and `y`
//!   alone is moved in those, and a solid standing on a plane near zero
//!   keeps its heights;
//! - the move is accepted only if it is *proved* exact: every vertex
//!   position, frame origin, curve origin and control point is `x - c` with
//!   no rounding ([`ExactBRep::translated_exactly`]). The recentred
//!   operands are then the same solids, shifted, and a vertex of an operand
//!   that survives in the result returns to its own coordinates (`(x - c) +
//!   c == x`). An operand whose move is not exact is not recentred, nor is
//!   one with a support a translation cannot carry;
//! - if the cut of the recentred operands is refused, they are cut again as
//!   given, in a session as it was opened: planar faces charted in world
//!   coordinates (#227, a contour with absolute coordinates under an
//!   identity placement) keep parameters as far out as before the move, and
//!   the move can be what the classifier cannot decide. The recentring never
//!   turns a boolean that works into a refusal;
//! - new vertices of the result, computed near the origin, are rounded once
//!   when moved back, to the ulp of their far coordinates: the least any
//!   result at that distance can carry.
//!
//! Operands within four diagonals of the origin on every axis keep their
//! own path unchanged, bit for bit.

use axiolid_brep::ExactBRep;
use axiolid_core::{Scalar, Transform3, Vec3};

use crate::BooleanError;

/// The centre is recentred on only beyond this many diagonals of the
/// operands' common box from the origin.
const DIAGONALS: Scalar = 4.0;

/// The operands moved by `-c`, and `c`, when they lie far from the origin
/// and move exactly.
pub(crate) fn operands(a: &ExactBRep, b: &ExactBRep) -> Option<(ExactBRep, ExactBRep, Vec3)> {
    let centre = centre(&[a, b])?;
    let near_a = a.translated_exactly(-centre).ok()??;
    let near_b = b.translated_exactly(-centre).ok()??;
    Some((near_a, near_b, centre))
}

/// A result moved back by `centre`.
///
/// # Errors
///
/// [`BooleanError::Evaluation`] where the result cannot be moved: it holds a
/// curve or surface family `ExactBRep::transformed` does not carry.
pub(crate) fn restore(result: &ExactBRep, centre: Vec3) -> Result<ExactBRep, BooleanError> {
    result
        .transformed(&Transform3::from_translation(centre))
        .map_err(|_| BooleanError::Evaluation)
}

/// The recentring offset of the operands, if they lie far enough from the
/// origin for it to matter: the centre of their vertex box, rounded per
/// axis to a multiple of the ulp of the largest coordinate, and zero on an
/// axis where it lies within [`DIAGONALS`] diagonals of zero.
pub(crate) fn centre(operands: &[&ExactBRep]) -> Option<Vec3> {
    let mut lo = Vec3::splat(Scalar::INFINITY);
    let mut hi = Vec3::splat(Scalar::NEG_INFINITY);
    for operand in operands {
        for vertex in operand.topology().vertices() {
            lo = lo.min(vertex.position);
            hi = hi.max(vertex.position);
        }
    }
    let largest = lo.abs().max(hi.abs()).max_element();
    if !(lo.is_finite() && hi.is_finite()) || !(largest > 0.0 && largest.is_finite()) {
        return None;
    }
    let quantum = ulp(largest)?;
    let mid = (lo + hi) * 0.5;
    // Per axis: an axis whose centre lies within a few diagonals of zero is
    // left alone, so a box far out in `x` and `y` only is moved in those
    // (its `z` near zero may hold coordinates `x - c` would round).
    let near = DIAGONALS * (hi - lo).length();
    let round = |mid: Scalar| {
        if mid.abs() <= near {
            0.0
        } else {
            (mid / quantum).round() * quantum
        }
    };
    let rounded = Vec3::new(round(mid.x), round(mid.y), round(mid.z));
    if !rounded.is_finite() || rounded == Vec3::ZERO {
        return None;
    }
    Some(rounded)
}

/// The ulp of `x`, a power of two.
fn ulp(x: Scalar) -> Option<Scalar> {
    let exponent = (x.to_bits() >> 52) & 0x7ff;
    // Normal numbers whose ulp is itself normal.
    (exponent > 52 && exponent < 0x7ff).then(|| Scalar::from_bits((exponent - 52) << 52))
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_construct::extrude::extrude_profile_exact;
    use axiolid_core::Tolerance;
    use axiolid_profile::{Profile, RectangleProfile};

    fn block(site: Vec3) -> ExactBRep {
        let profile = Profile::Rectangle(RectangleProfile {
            x: 0.4,
            y: 0.5,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        });
        extrude_profile_exact(&profile, Vec3::Z, 6.0, Tolerance::METRE)
            .expect("a block")
            .transformed(&Transform3::from_translation(site))
            .expect("placed")
    }

    #[test]
    fn the_ulp_is_the_spacing_above_a_number() {
        for x in [1.0_f64, 1.5, 6e5, 5.6e6, 123_456.789] {
            let spacing = Scalar::from_bits(x.to_bits() + 1) - x;
            assert_eq!(ulp(x), Some(spacing), "{x}");
        }
        assert_eq!(ulp(0.0), None);
        assert_eq!(ulp(Scalar::NAN), None);
    }

    #[test]
    fn operands_near_the_origin_are_not_recentred() {
        for site in [
            Vec3::ZERO,
            Vec3::new(1.0, 2.0, 0.0),
            Vec3::new(0.0, 0.0, 10.0),
        ] {
            let a = block(site);
            assert_eq!(centre(&[&a]), None, "{site:?}");
            assert!(operands(&a, &a).is_none(), "{site:?}");
        }
    }

    #[test]
    fn a_far_centre_is_a_multiple_of_the_ulp_of_the_largest_coordinate() {
        let site = Vec3::new(6e5 + 0.1, 5.6e6 + 0.3, 200.7);
        let a = block(site);
        let c = centre(&[&a]).expect("far from the origin");
        let q = ulp(5.6e6 + 3.0).expect("a normal number");
        for axis in [c.x, c.y, c.z] {
            assert_eq!(axis % q, 0.0, "{axis}");
        }
        assert!((c - site - Vec3::new(0.0, 0.0, 3.0)).length() <= 1e-6);
    }

    #[test]
    fn recentred_operands_come_back_to_their_own_coordinates() {
        let a = block(Vec3::new(6e5, 5.6e6, 200.0));
        let b = block(Vec3::new(6e5, 5.6e6, 212.0));
        let (near_a, near_b, c) = operands(&a, &b).expect("far and exact");
        // The centre of both, not of the first alone (z = 203).
        assert_eq!(c.z, 209.0);
        assert!(near_a
            .topology()
            .vertices()
            .iter()
            .all(|v| v.position.length() < 10.0));
        assert!(near_b
            .topology()
            .vertices()
            .iter()
            .all(|v| v.position.length() < 10.0));
        assert_eq!(restore(&near_a, c).expect("moved back"), a);
        assert_eq!(restore(&near_b, c).expect("moved back"), b);
    }

    #[test]
    fn an_axis_near_zero_keeps_its_own_coordinates() {
        // x and y far out, z in [0.3, 6.3]: moved in x and y alone.
        let a = block(Vec3::new(6e5, 5.6e6, 0.3));
        let c = centre(&[&a]).expect("far out in x and y");
        assert_eq!((c.x, c.y, c.z), (6e5, 5.6e6, 0.0));
        let (near, _, c) = operands(&a, &a).expect("exact");
        assert_eq!(restore(&near, c).expect("moved back"), a);
    }
}
