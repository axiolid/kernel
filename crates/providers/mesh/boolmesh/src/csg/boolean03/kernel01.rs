//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

use crate::csg::{Half, Real, Vec2, Vec3, Vec4};
// These two functions (Interpolate and Intersect) are the only places where
// floating-point operations take place in the whole Boolean function. These
// are carefully designed to minimize rounding error and to remove it at edge
// cases to ensure consistency.

pub fn interpolate(pl: Vec3, pr: Vec3, x: Real) -> Vec2 {
    let dx_l = x - pl.x;
    let dx_r = x - pr.x;
    let diff = pr - pl;
    let use_l = dx_l.abs() < dx_r.abs();
    let lambda = if use_l { dx_l } else { dx_r } / diff.x;

    if lambda.is_infinite()
        || lambda.is_nan()
        || diff.y.is_infinite()
        || diff.y.is_nan()
        || diff.z.is_infinite()
        || diff.z.is_nan()
    {
        return Vec2::new(pl.y, pl.z);
    }

    Vec2::new(
        lambda * diff.y + if use_l { pl.y } else { pr.y },
        lambda * diff.z + if use_l { pl.z } else { pr.z },
    )
}

pub fn intersect(pl: Vec3, pr: Vec3, ql: Vec3, qr: Vec3) -> Vec4 {
    let dy_l = ql.y - pl.y;
    let dy_r = qr.y - pr.y;
    assert!(dy_l * dy_r <= 0., "Boolean manifold error: no intersection");
    let use_l = dy_l.abs() < dy_r.abs();
    let dx = pr.x - pl.x;
    let mut lambda = if use_l { dy_l } else { dy_r } / (dy_l - dy_r);
    if lambda.is_infinite() || lambda.is_nan() {
        lambda = 0.;
    }
    let p_dy = pr.y - pl.y;
    let q_dy = qr.y - ql.y;
    let use_p = p_dy.abs() < q_dy.abs();
    // A 2x2 choice on (use_l, use_p): pick the endpoint, then the operand.
    // Nested if/else obscured that it is one selection, not two decisions.
    let y_base = match (use_l, use_p) {
        (true, true) => pl.y,
        (true, false) => ql.y,
        (false, true) => pr.y,
        (false, false) => qr.y,
    };
    // Built in one expression rather than `default()` + field writes, so the
    // value is never observable in a partly-initialised state.
    Vec4::new(
        lambda * dx + if use_l { pl.x } else { pr.x },
        lambda * if use_p { p_dy } else { q_dy } + y_base,
        lambda * (pr.z - pl.z) + if use_l { pl.z } else { pr.z },
        lambda * (qr.z - ql.z) + if use_l { ql.z } else { qr.z },
    )
}

pub fn shadows(p: Real, q: Real, dir: Real) -> bool {
    if p == q {
        dir < 0.
    } else {
        p < q
    }
}

// This is equivalent to Kernel01 or X01 in the thesis.
// Expand represents the sign of the normal.
/// The two operands as `shadows01` sees them, in caller-chosen order.
///
/// `shadows01` is called with P and Q swapped depending on which side is
/// casting, so this is deliberately not "P and Q" but "caster and
/// receiver": kernel11 calls it once each way within a single operation.
/// Grouping the four slices keeps that symmetry visible at the call site.
pub struct ShadowOperands<'a> {
    /// Positions of the operand supplying the point.
    pub ps_p: &'a [Vec3],
    /// Positions of the operand supplying the edge.
    pub ps_q: &'a [Vec3],
    /// Halfedges of the operand supplying the edge.
    pub hs_q: &'a [Half],
    /// Vertex normals used to expand the shadow test.
    pub ns: &'a [Vec3],
}

pub fn shadows01(
    p0: usize,
    q1: usize,
    ops: &ShadowOperands,
    expand: Real,
    reverse: bool,
) -> Option<(i32, Vec2)> {
    let q1s = ops.hs_q[q1].tail;
    let q1e = ops.hs_q[q1].head;
    let p0x = ops.ps_p[p0].x;
    let q1sx = ops.ps_q[q1s].x;
    let q1ex = ops.ps_q[q1e].x;

    // check weather the vert is in between the half from the x-axis point of view
    let mut s01 = if reverse {
        let a = if shadows(q1sx, p0x, expand * ops.ns[q1s].x) {
            1
        } else {
            0
        };
        let b = if shadows(q1ex, p0x, expand * ops.ns[q1e].x) {
            1
        } else {
            0
        };
        a - b
    } else {
        let a = if shadows(p0x, q1ex, expand * ops.ns[p0].x) {
            1
        } else {
            0
        };
        let b = if shadows(p0x, q1sx, expand * ops.ns[p0].x) {
            1
        } else {
            0
        };
        a - b
    };

    // if in between...
    if s01 != 0 {
        let yz01 = interpolate(ops.ps_q[q1s], ops.ps_q[q1e], ops.ps_p[p0].x);
        if reverse {
            let d1 = ops.ps_q[q1s] - ops.ps_p[p0];
            let d2 = ops.ps_q[q1e] - ops.ps_p[p0];
            let sta2 = d1.length_squared();
            let end2 = d2.length_squared();
            let dir = if sta2 < end2 {
                ops.ns[q1s].y
            } else {
                ops.ns[q1e].y
            };
            if !shadows(yz01[0], ops.ps_p[p0].y, expand * dir) {
                s01 = 0;
            }
        } else {
            // return sign as 0 if vert from mfd_p is above
            if !shadows(ops.ps_p[p0].y, yz01[0], expand * ops.ns[p0].y) {
                s01 = 0;
            }
        }
        return Some((s01, yz01));
    }
    None
}
