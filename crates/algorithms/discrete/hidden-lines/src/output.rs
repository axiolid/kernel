//! Output endpoints: the exact rational point, rounded once.
//!
//! A piece's endpoint is the parameter `t` of the edge `A -> B`; each image
//! coordinate `U_A + t·(U_B - U_A)` is evaluated exactly and rounded once to
//! the nearest binary64 (ties to even). For a root `t = a / (a - b)` this is
//! `(a·U_B - b·U_A) / (a - b)`. The result depends on the exact point alone,
//! so equal points reached through different constraints are bit-identical.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Point2;
use axiolid_exact::{Arith, Dyadic};
use axiolid_exact_arrangement::round_ratio;

use crate::param::Param;
use crate::view::Vert;

/// The image point of the edge `start -> end` at `param`.
pub(crate) fn endpoint(verts: &[Vert], start: u32, end: u32, param: Param) -> GeomResult<Point2> {
    // Homogeneous: the point of vertex `v` is `exact[axis] / exact[3]`, and the
    // point of the edge where `L` vanishes is `(a q - b p) / (a h_q - b h_p)`.
    let coord = |vertex: u32, axis: usize| &verts[vertex as usize].exact[axis];
    let weight = |vertex: u32| &verts[vertex as usize].exact[3];
    let mut out = [0.0; 2];
    for (axis, slot) in out.iter_mut().enumerate() {
        let rounded = match param {
            Param::Zero => round_ratio(coord(start, axis), weight(start)),
            Param::One => round_ratio(coord(end, axis), weight(end)),
            Param::Root(root) => {
                let a = root.constraint.value::<Dyadic>(verts, start);
                let b = root.constraint.value::<Dyadic>(verts, end);
                let numerator = a.mul(coord(end, axis)).sub(&b.mul(coord(start, axis)));
                let denominator = a.mul(weight(end)).sub(&b.mul(weight(start)));
                round_ratio(&numerator, &denominator)
            }
        };
        *slot = rounded.ok_or_else(|| {
            GeomError::Degenerate("an image coordinate is outside the binary64 range".into())
        })?;
    }
    Ok(Point2::new(out[0], out[1]))
}
