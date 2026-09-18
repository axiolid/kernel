//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

use super::{
    collapse_triangle, form_loops, head_of, hids_of, next_of, pair_of, remove_if_folded,
    update_vid_around_star,
};
use crate::csg::{is_ccw_3d, Half, Real, Tref, Vec3};

// Check around a halfedges from the same tail vertex.
// If they consist of only two tris, then their edge is collapsable.
fn record_if_collinear(hs: &[Half], ns: &[Vec3], rs: &[Tref], hid: usize, nv: usize) -> bool {
    let h = &hs[hid];
    if h.pair().is_none() || (h.tail < nv) {
        return false;
    }

    let cw_next = |i: usize| next_of(hs[i].pair);

    let bgn = hid;
    let mut cur = cw_next(bgn);
    let t0 = bgn / 3;
    let mut t1 = cur / 3;
    let mut same = is_coplanar_at(ns, t0, t1, rs);
    while cur != bgn {
        cur = cw_next(cur);
        let t2 = cur / 3;
        if !is_coplanar_at(ns, t2, t0, rs) && !is_coplanar_at(ns, t2, t1, rs) {
            if same {
                t1 = t2;
                same = false;
            } else {
                return false;
            }
        }
    }
    true
}

fn record_if_short(hs: &[Half], ps: &[Vec3], hid: usize, nv: usize, ep: Real) -> bool {
    let h = &hs[hid];
    if h.pair().is_none() || (h.tail < nv && h.head < nv) {
        return false;
    }
    (ps[hs[hid].head] - ps[hs[hid].tail]).length_squared() < ep.powi(2)
}

pub fn collapse_edge(
    hs: &mut [Half],
    ps: &mut Vec<Vec3>,
    ns: &mut [Vec3],
    rs: &mut [Tref],
    hid: usize,
    eps: Real,
    store: &mut Vec<usize>, // storing the halfedge data for form_loops
) -> bool {
    let to_rmv = &hs[hid];
    if to_rmv.pair().is_none() {
        return false;
    }

    let vid_keep = to_rmv.head;
    let vid_delt = to_rmv.tail;
    let pos_keep = ps[vid_keep];
    let pos_delt = ps[vid_delt];

    let t0 = hids_of(hid);
    let t1 = hids_of(to_rmv.pair);

    let mut bgn = pair_of(hs, t1.1); // the bgn half heading delt vert
    let end = t0.2; // the end half heading delt vert

    // check validity by orbiting start vert ccw order
    if (pos_keep - pos_delt).length_squared() >= eps.powi(2) {
        let mut cur = bgn;
        let mut tri0 = to_rmv.pair / 3;
        let mut p_prev = ps[head_of(hs, t1.1)];
        while cur != to_rmv.pair {
            cur = next_of(cur); // incoming half around delt vert
            let p_next = ps[head_of(hs, cur)];
            let tri_curr = cur / 3;
            let n_curr = &ns[cur / 3];
            let n_pair = &ns[to_rmv.pair / 3];
            let ccw = |p0, p1, p2| is_ccw_3d(p0, p1, p2, n_curr, eps);
            if !is_coplanar_at(ns, tri_curr, tri0, rs) {
                let tri2 = tri0;
                tri0 = hid / 3;
                if !is_coplanar_at(ns, tri_curr, tri0, rs) {
                    return false;
                }
                if rs[tri0].mid != rs[tri2].mid || n_pair.dot(*n_curr) < -0.5 {
                    // Restrict collapse to co-linear edges when the edge separates faces or the edge is sharp.
                    // This ensures large shifts are not introduced parallel to the tangent plane.
                    if ccw(&p_prev, &pos_delt, &pos_keep) != 0 {
                        return false;
                    }
                }
            }

            // Don't collapse edge if it would cause a triangle to invert
            if ccw(&p_next, &p_prev, &pos_keep) < 0 {
                return false;
            }

            p_prev = p_next;
            cur = pair_of(hs, cur); // outgoing half around delt vert
        }
    }

    // find a candidate by orbiting end verts ccw order
    let mut cur = pair_of(hs, t0.1);
    while cur != t1.2 {
        cur = next_of(cur);
        store.push(cur); // storing outgoing half here
        cur = pair_of(hs, cur);
    }

    ps[to_rmv.tail] = Vec3::NAN;
    collapse_triangle(hs, &t1);

    let mut cur = bgn;
    while cur != end {
        cur = next_of(cur);
        let pair = pair_of(hs, cur);
        let head = head_of(hs, cur);
        if let Some((i, &v)) = store
            .iter()
            .enumerate()
            .find(|&(_, &s)| head_of(hs, s) == head)
        {
            form_loops(hs, ps, v, cur);
            bgn = pair;
            store.truncate(i);
        }
        cur = pair;
    }

    // do collapse
    update_vid_around_star(hs, bgn, end, vid_keep);
    collapse_triangle(hs, &t0);
    remove_if_folded(hs, ps, bgn);
    true
}

pub fn collapse_collinear_edges(
    hs: &mut [Half],
    ps: &mut Vec<Vec3>,
    ns: &mut [Vec3],
    rs: &mut [Tref],
    nv: usize,
    ep: Real,
) {
    let mut _flag = 0;
    let rec = (0..hs.len())
        .filter(|&hid| record_if_collinear(hs, ns, rs, hid, nv))
        .collect::<Vec<_>>();
    for hid in rec {
        if collapse_edge(hs, ps, ns, rs, hid, ep, &mut vec![]) {
            _flag += 1;
        }
    }
}

pub fn collapse_short_edges(
    hs: &mut [Half],
    ps: &mut Vec<Vec3>,
    ns: &mut [Vec3],
    rs: &mut [Tref],
    nv: usize,
    ep: Real,
) {
    loop {
        let mut flag = 0;
        let rec = (0..hs.len())
            .filter(|&hid| record_if_short(hs, ps, hid, nv, ep))
            .collect::<Vec<_>>();
        for hid in rec {
            if collapse_edge(hs, ps, ns, rs, hid, ep, &mut vec![]) {
                flag += 1;
            }
        }
        if flag == 0 {
            break;
        }
    }
}

/// Do two triangles lie on the same plane?
///
/// Provenance (`mid`/`pid`) is the fast path: faces from one operand's own
/// coplanar-index entry are known coplanar without arithmetic.
///
/// But provenance is only a PROXY for the geometric question. `pid` indexes
/// the coplanar set WITHIN one operand, and `mid` says which operand a face
/// came from, so two faces that are genuinely on one plane compare unequal
/// whenever they arrive from different operands -- or from the same operand
/// under a later call, where the index was recomputed. A reconstructed seam
/// is exactly that case: the halves of a split face carry different ids, so
/// the seam never collapses and the result stays two shells (kernel#100).
///
/// Falling back to the normals answers the real question. Being a fallback,
/// it cannot make previously-merged faces stop merging: it only admits pairs
/// provenance rejected.
#[inline]
fn is_coplanar_at(ns: &[Vec3], t0: usize, t1: usize, rs: &[Tref]) -> bool {
    if rs[t0].mid == rs[t1].mid && rs[t0].pid == rs[t1].pid {
        return true;
    }
    // Same plane means parallel normals. Antiparallel counts too: a seam
    // between two solids has opposed windings by construction.
    let (n0, n1) = (ns[t0], ns[t1]);
    n0.cross(n1).length_squared() <= COPLANAR_EPS && n0.dot(n1).abs() > 0.0
}

/// Squared sine of the angle two normals may differ by and still count as
/// one plane. Normals here are unit length, so this is an angle bound, not
/// a scale-dependent distance.
const COPLANAR_EPS: Real = 1e-20;
