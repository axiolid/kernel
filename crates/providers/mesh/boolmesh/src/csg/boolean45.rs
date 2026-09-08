//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

use crate::csg::boolean03::Boolean03;
use crate::csg::bounds::BBox;
use crate::csg::OpType;
use crate::csg::{face_of, Half, Manifold, Real, Tref, Vec3};
use std::collections::HashMap;
use std::mem;

fn duplicate_verts(inc: &[i32], vt_r: &[i32], ps_p: &[Vec3], ps_r: &mut [Vec3], vid: usize) {
    let n = inc[vid].unsigned_abs() as usize;
    for i in 0..n {
        ps_r[vt_r[vid] as usize + i] = ps_p[vid];
    }
}

fn inclusive_scan(input: &[i32], output: &mut [i32], offset: i32) {
    if input.is_empty() || output.is_empty() {
        return;
    }
    let mut sum = offset;
    for (i, &v) in input.iter().enumerate() {
        sum += v;
        if i < output.len() {
            output[i] = sum;
        }
    }
}

fn exclusive_scan(input: &[i32], output: &mut [i32], offset: i32) {
    if input.is_empty() || output.is_empty() {
        return;
    }
    let mut sum = offset;
    output[0] = sum;
    for i in 1..input.len() {
        sum += input[i - 1];
        if i < output.len() {
            output[i] = sum;
        }
    }
}

/// The result mesh under construction.
///
/// `append_partial_edges`, `append_new_edges`, and `append_whole_edges`
/// each fill the same three parallel arrays, and every call site was
/// threading them individually. Bundling them names the thing being
/// built and drops three arguments from each signature.
struct ResultEdges<'a> {
    /// Halfedge data of the result, filled in as edges are appended.
    hs: &'a mut [Half],
    /// Maps a result halfedge back to the triangle it came from.
    rs: &'a mut [Tref],
    /// Write cursor per result face.
    face_ptr: &'a mut [i32],
}

/// One operand's contribution to the result, with its index maps.
///
/// `append_partial_edges` and `append_whole_edges` are each called twice,
/// once per operand, and every call threads the same halfedges, winding
/// numbers, and the two P-to-R maps. Grouping them names "the side being
/// appended" and makes the two calls visibly symmetric.
struct SourceSide<'a> {
    /// Winding contribution per vertex of this operand.
    i03: &'a [i32],
    /// Halfedges of this operand.
    hs: &'a [Half],
    /// Maps a vertex of this operand to its result vertex.
    vid2r: &'a [i32],
    /// Maps a face of this operand to its result face.
    fid2r: &'a [i32],
    /// Whether this operand is the forward (P) side.
    fwd: bool,
}

/// Winding contributions, adjusted for the operation being performed.
///
/// All four are derived together from the `Boolean03` result by the same
/// three coefficients (which encode union vs difference vs intersection),
/// and every consumer needs the set rather than any one of them. Deriving
/// them in one place keeps the operation's sign convention in one place
/// too.
struct Windings {
    /// Per-vertex winding for the P operand.
    i03: Vec<i32>,
    /// Per-vertex winding for the Q operand.
    i30: Vec<i32>,
    /// Per-intersection winding, P edge against Q face.
    i12: Vec<i32>,
    /// Per-intersection winding, Q edge against P face.
    i21: Vec<i32>,
}

impl Windings {
    /// Apply the operation's coefficients to the raw `boolean03` result.
    fn new(b03: &Boolean03, c1: i32, c2: i32, c3: i32) -> Self {
        Self {
            i12: b03.x12.iter().map(|v| c3 * v).collect(),
            i21: b03.x21.iter().map(|v| c3 * v).collect(),
            i03: b03.w03.iter().map(|v| c1 + c3 * v).collect(),
            i30: b03.w30.iter().map(|v| c2 + c3 * v).collect(),
        }
    }
}

fn size_output(
    mp: &Manifold,
    mq: &Manifold,
    w: &Windings,
    b03: &Boolean03,
    fns: &mut Vec<Vec3>,
    inv: bool, // whether to invert mesh of q
) -> (Vec<i32>, Vec<i32>) {
    let mut side_p = vec![0; mp.nf];
    let mut side_q = vec![0; mq.nf];

    // equivalent to CountVerts
    for (i, h) in mp.hs.iter().enumerate() {
        side_p[face_of(i)] += w.i03[h.tail].abs();
    }
    for (i, h) in mq.hs.iter().enumerate() {
        side_q[face_of(i)] += w.i30[h.tail].abs();
    }

    // equivalent to CountNewVerts
    for i in 0..w.i12.len() {
        let hid0 = b03.p1q2[i][0];
        let hid1 = mp.hs[hid0].pair;
        let inc = w.i12[i].abs();
        side_p[face_of(hid0)] += inc;
        side_p[face_of(hid1)] += inc;
        side_q[b03.p1q2[i][1]] += inc;
    }

    for i in 0..w.i21.len() {
        let hid0 = b03.p2q1[i][1];
        let hid1 = mq.hs[hid0].pair;
        let inc = w.i21[i].abs();
        side_q[face_of(hid0)] += inc;
        side_q[face_of(hid1)] += inc;
        side_p[b03.p2q1[i][0]] += inc;
    }

    // a map from face_p and face_q to face_r
    let mut face_pq2r = vec![0; mp.nf + mq.nf + 1];
    let side_pq = [&side_p[..], &side_q[..]].concat();
    let keep_fs = side_pq
        .iter()
        .map(|&x| if x > 0 { 1 } else { 0 })
        .collect::<Vec<_>>();

    inclusive_scan(&keep_fs, &mut face_pq2r[1..], 0);
    let nf_r = *face_pq2r.last().unwrap() as usize;
    face_pq2r.truncate(mp.nf + mq.nf);
    fns.resize(nf_r, Vec3::ZERO);

    let mut fid_r = 0;
    for (i, n) in mp.face_normals.iter().enumerate() {
        if side_p[i] > 0 {
            fns[fid_r] = *n;
            fid_r += 1;
        }
    }
    for (i, n) in mq.face_normals.iter().enumerate() {
        if side_q[i] > 0 {
            fns[fid_r] = *n * if inv { -1. } else { 1. };
            fid_r += 1;
        }
    }

    let truncated = side_pq
        .iter()
        .filter(|s| **s > 0)
        .copied()
        .collect::<Vec<_>>();
    let mut ih_per_f = vec![0; truncated.len()];

    inclusive_scan(&truncated, &mut ih_per_f, 0);
    ih_per_f.insert(0, 0);

    (ih_per_f, face_pq2r)
}

// Sort of intermediate data store for halfedge creation
#[derive(Clone, Debug)]
struct EdgePt {
    val: Real,     // dot value of edge
    vid: usize,    // vertex id
    cid: usize,    // collision id
    is_tail: bool, //
}

/// Where new intersection vertices accumulate before triangulation.
///
/// `pt_old` keys by the halfedge a vertex splits; `pt_new` keys by the face
/// pair that created it. They are filled together on every call, so passing
/// them as one value keeps the two maps' relationship explicit.
struct EdgePoints<'a> {
    /// New vertices lying on an existing halfedge.
    on_edge: &'a mut HashMap<usize, Vec<EdgePt>>,
    /// New vertices interior to a face pair.
    on_face: &'a mut HashMap<(usize, usize), Vec<EdgePt>>,
}

fn add_new_edge_verts(
    p1q2: &[[usize; 2]],
    i12: &[i32],
    v12_r: &[i32],
    hs_p: &[Half],
    pts: &mut EdgePoints,
    fwd: bool,
    oft: usize,
) {
    for i in 0..p1q2.len() {
        let hid_p = p1q2[i][if fwd { 0 } else { 1 }];
        let fid_q = p1q2[i][if fwd { 1 } else { 0 }];
        let vid_r = v12_r[i] as usize;
        let inc = i12[i];
        let hid0 = hid_p;
        let hid1 = hs_p[hid_p].pair;
        let key_l = if fwd {
            (face_of(hid0), fid_q)
        } else {
            (fid_q, face_of(hid0))
        };
        let key_r = if fwd {
            (face_of(hid1), fid_q)
        } else {
            (fid_q, face_of(hid1))
        };
        let dir = inc < 0;
        pts.on_edge.entry(hid_p).or_default();
        pts.on_face.entry(key_l).or_default();
        pts.on_face.entry(key_r).or_default();
        let dir0 = dir ^ !fwd;
        let dir1 = dir ^ fwd;
        let inc_ = inc.unsigned_abs() as usize;
        for j in 0..inc_ {
            pts.on_edge.get_mut(&hid_p).unwrap().push(EdgePt {
                val: 0.,
                vid: vid_r + j,
                cid: i + oft,
                is_tail: dir,
            });
        }
        for j in 0..inc_ {
            pts.on_face.get_mut(&key_r).unwrap().push(EdgePt {
                val: 0.,
                vid: vid_r + j,
                cid: i + oft,
                is_tail: dir0,
            });
        }
        for j in 0..inc_ {
            pts.on_face.get_mut(&key_l).unwrap().push(EdgePt {
                val: 0.,
                vid: vid_r + j,
                cid: i + oft,
                is_tail: dir1,
            });
        }
    }
}

// Creating a partial halfedges from a list of positions.
// It's very confusing, but it's not aiming to pair twins (pair is -1).
// It's more likely to say pairing sta-end vertex and make a halfedge
fn pair_up(pts: &mut [EdgePt]) -> Vec<Half> {
    assert_eq!(pts.len() % 2, 0);
    let nh = pts.len() / 2;
    let mid_idx = {
        let mut sta_idx = 0;
        let mut end_idx = pts.len();

        while sta_idx < end_idx {
            if pts[sta_idx].is_tail {
                sta_idx += 1;
            } else {
                end_idx -= 1;
                pts.swap(sta_idx, end_idx);
            }
        }
        sta_idx
    };

    let cmp = |a: &EdgePt, b: &EdgePt| {
        a.val
            .partial_cmp(&b.val)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.cid.cmp(&b.cid))
    };
    // The comparator breaks ties on `cid`, so the order is total and a
    // stable sort has nothing left to preserve.
    pts[..mid_idx].sort_unstable_by(cmp);
    pts[mid_idx..].sort_unstable_by(cmp);

    let mut edges = Vec::with_capacity(nh);
    for i in 0..nh {
        edges.push(Half::new_without_pair(pts[i].vid, pts[i + nh].vid));
    }
    edges
}

fn append_partial_edges(
    side: &SourceSide,
    ps_p: &[Vec3],
    ps_r: &[Vec3], // the vert pos of mfd_r, already fulfilled so far
    out: &mut ResultEdges,
    pt_p: &mut HashMap<usize, Vec<EdgePt>>, //
    whole_flag: &mut [bool], // a flag to find out a halfedge from mfd_p is entirely usable in mfd_r
) {
    for (hid_p, pt) in pt_p {
        let hpos_p = pt;
        let h = &side.hs[*hid_p];
        whole_flag[*hid_p] = false;
        whole_flag[h.pair] = false;

        // assigning 0-1 value to hpos_p
        let dif = ps_p[h.head] - ps_p[h.tail];
        for p in hpos_p.iter_mut() {
            p.val = dif.dot(ps_r[p.vid]);
        }

        let i_tail = side.i03[h.tail]; // mostly 0 or 1
        let i_head = side.i03[h.head]; // mostly 0 or 1
        let p_tail = ps_r[side.vid2r[h.tail] as usize];
        let p_head = ps_r[side.vid2r[h.head] as usize];

        for i in 0..i_tail.unsigned_abs() as usize {
            hpos_p.push(EdgePt {
                val: p_tail.dot(dif),
                vid: side.vid2r[h.tail] as usize + i,
                cid: usize::MAX,
                is_tail: i_tail > 0,
            });
        }
        for i in 0..i_head.unsigned_abs() as usize {
            hpos_p.push(EdgePt {
                val: p_head.dot(dif),
                vid: side.vid2r[h.head] as usize + i,
                cid: usize::MAX,
                is_tail: i_head < 0,
            });
        }

        let mut half_seq = pair_up(hpos_p);
        let fp_l = face_of(*hid_p);
        let fp_r = face_of(h.pair);
        let fid_l = side.fid2r[fp_l] as usize;
        let fid_r = side.fid2r[fp_r] as usize;

        // Negative inclusion means the halfedges are reversed, which means our
        // reference is now to the head instead of the tail, which is one
        // position advanced CCW. This is only valid if this is a retained vert;
        // it will be ignored later if the vert is new.

        let fw_tri = Tref {
            mid: if side.fwd { 0 } else { 1 },
            fid: fp_l,
            ..Default::default()
        };
        let bk_tri = Tref {
            mid: if side.fwd { 0 } else { 1 },
            fid: fp_r,
            ..Default::default()
        };

        for h in half_seq.iter_mut() {
            let fw_edge = out.face_ptr[fid_l] as usize;
            let bk_edge = out.face_ptr[fid_r] as usize;
            out.face_ptr[fid_l] += 1;
            out.face_ptr[fid_r] += 1;
            out.hs[fw_edge] = Half::new(h.tail, h.head, bk_edge);
            out.hs[bk_edge] = Half::new(h.head, h.tail, fw_edge);
            out.rs[fw_edge] = fw_tri;
            out.rs[bk_edge] = bk_tri;
        }
    }
}

fn append_new_edges(
    ps_r: &[Vec3],    // the vert pos of mfd_r, already fulfilled so far
    fid_pq2r: &[i32], //
    nf_p: usize,      // num of faces in mfd_p
    pt_new: &mut HashMap<(usize, usize), Vec<EdgePt>>,
    out: &mut ResultEdges,
) {
    for ((fid_p, fid_q), pt_init) in pt_new.iter_mut() {
        let pt = pt_init;
        let mut bb = BBox::default();
        for p in pt.iter() {
            bb.union(&ps_r[p.vid]);
        }

        let d = bb.longest_dim();
        for p in pt.iter_mut() {
            p.val = ps_r[p.vid][d];
        }

        let mut half_seq = pair_up(pt);
        let fid_l = fid_pq2r[*fid_p] as usize;
        let fid_r = fid_pq2r[*fid_q + nf_p] as usize;
        let fw_ref = Tref {
            mid: 0,
            fid: *fid_p,
            ..Default::default()
        };
        let bk_ref = Tref {
            mid: 1,
            fid: *fid_q,
            ..Default::default()
        };

        for h in half_seq.iter_mut() {
            let fw_edge = out.face_ptr[fid_l] as usize;
            let bk_edge = out.face_ptr[fid_r] as usize;
            out.face_ptr[fid_l] += 1;
            out.face_ptr[fid_r] += 1;
            out.hs[fw_edge] = Half::new(h.tail, h.head, bk_edge);
            out.hs[bk_edge] = Half::new(h.head, h.tail, fw_edge);
            out.rs[fw_edge] = fw_ref;
            out.rs[bk_edge] = bk_ref;
        }
    }
}

fn append_whole_edges(side: &SourceSide, whole_flag: &[bool], out: &mut ResultEdges) {
    for (i, hp) in side.hs.iter().enumerate() {
        if !whole_flag[i] || !hp.is_forward() {
            continue;
        }

        let mut h = hp.clone();
        let inc = side.i03[h.tail];
        if inc == 0 {
            continue;
        }
        if inc < 0 {
            mem::swap(&mut h.tail, &mut h.head);
        }

        h.tail = side.vid2r[h.tail] as usize;
        h.head = side.vid2r[h.head] as usize;

        let fp_l = face_of(i);
        let fp_r = face_of(hp.pair);
        let fid_l = side.fid2r[fp_l] as usize;
        let fid_r = side.fid2r[fp_r] as usize;
        let fw_ref = Tref {
            mid: if side.fwd { 0 } else { 1 },
            fid: fp_l,
            ..Default::default()
        };
        let bk_ref = Tref {
            mid: if side.fwd { 0 } else { 1 },
            fid: fp_r,
            ..Default::default()
        };

        for _ in 0..inc.unsigned_abs() as usize {
            let fw_edge = out.face_ptr[fid_l] as usize;
            let bk_edge = out.face_ptr[fid_r] as usize;
            out.face_ptr[fid_l] += 1;
            out.face_ptr[fid_r] += 1;
            out.hs[fw_edge] = Half::new(h.tail, h.head, bk_edge);
            out.hs[bk_edge] = Half::new(h.head, h.tail, fw_edge);
            out.rs[fw_edge] = fw_ref;
            out.rs[bk_edge] = bk_ref;
            h.tail += 1;
            h.head += 1;
        }
    }
}

pub struct Boolean45 {
    pub ps: Vec<Vec3>,
    pub ns: Vec<Vec3>,
    pub hs: Vec<Half>,
    pub rs: Vec<Tref>,
    pub hid_per_f: Vec<i32>,
    pub nv_from_p: usize,
    pub nv_from_q: usize,
}

pub fn boolean45(mp: &Manifold, mq: &Manifold, b03: &Boolean03, op: &OpType) -> Boolean45 {
    let c1 = if op == &OpType::Intersect { 0 } else { 1 };
    let c2 = if op == &OpType::Add { 1 } else { 0 };
    let c3 = if op == &OpType::Intersect { 1 } else { -1 };
    let w = Windings::new(b03, c1, c2, c3);
    let mut nv = 0;
    let mut vid_p2r = vec![0; mp.nv];
    let mut vid_q2r = vec![0; mq.nv];
    let mut vid_12r = vec![0; b03.v12.len()];
    let mut vid_21r = vec![0; b03.v21.len()];

    exclusive_scan(
        &w.i03.iter().map(|i| i.abs()).collect::<Vec<_>>(),
        &mut vid_p2r,
        nv,
    );
    nv = (*vid_p2r.last().unwrap()).abs() + w.i03.last().unwrap().abs();
    let nv_rp = nv;

    exclusive_scan(
        &w.i30.iter().map(|i| i.abs()).collect::<Vec<_>>(),
        &mut vid_q2r,
        nv,
    );
    nv = (*vid_q2r.last().unwrap()).abs() + w.i30.last().unwrap().abs();
    let nv_rq = nv - nv_rp;

    if !b03.v12.is_empty() {
        exclusive_scan(
            &w.i12.iter().map(|i| i.abs()).collect::<Vec<_>>(),
            &mut vid_12r,
            nv,
        );
        nv = (*vid_12r.last().unwrap()).abs() + w.i12.last().unwrap().abs();
    }
    let nv_12 = nv - nv_rp - nv_rq;

    if !b03.v21.is_empty() {
        exclusive_scan(
            &w.i21.iter().map(|i| i.abs()).collect::<Vec<_>>(),
            &mut vid_21r,
            nv,
        );
        nv = (*vid_21r.last().unwrap()).abs() + w.i21.last().unwrap().abs();
    }
    let nv_21 = nv - nv_rp - nv_rq - nv_12;

    let mut ps_r = vec![Vec3::ZERO; nv as usize];

    for i in 0..mp.nv {
        duplicate_verts(&w.i03, &vid_p2r, &mp.ps, &mut ps_r, i);
    }
    for i in 0..mq.nv {
        duplicate_verts(&w.i30, &vid_q2r, &mq.ps, &mut ps_r, i);
    }
    for i in 0..nv_12 {
        duplicate_verts(&w.i12, &vid_12r, &b03.v12, &mut ps_r, i as usize);
    }
    for i in 0..nv_21 {
        duplicate_verts(&w.i21, &vid_21r, &b03.v21, &mut ps_r, i as usize);
    }

    let mut pt_p = HashMap::new();
    let mut pt_q = HashMap::new();
    let mut pt_new = HashMap::new();
    add_new_edge_verts(
        &b03.p1q2,
        &w.i12,
        &vid_12r,
        &mp.hs,
        &mut EdgePoints {
            on_edge: &mut pt_p,
            on_face: &mut pt_new,
        },
        true,
        0,
    );
    add_new_edge_verts(
        &b03.p2q1,
        &w.i21,
        &vid_21r,
        &mq.hs,
        &mut EdgePoints {
            on_edge: &mut pt_q,
            on_face: &mut pt_new,
        },
        false,
        b03.p1q2.len(),
    );

    let mut ns_r = vec![];
    let inv = op == &OpType::Subtract;
    let (hid_per_f, fid_pq2r) = size_output(mp, mq, &w, b03, &mut ns_r, inv);

    let nh = *hid_per_f.last().unwrap() as usize;
    let mut face_ptr_r = hid_per_f.clone();
    let mut whole_flag_p = vec![true; mp.nh];
    let mut whole_flag_q = vec![true; mq.nh];
    let mut rs_r = vec![Tref::default(); nh];
    let mut hs_r = vec![Half::default(); nh];
    let fid_p2r = &fid_pq2r[0..mp.nf];
    let fid_q2r = &fid_pq2r[mp.nf..];

    let mut out = ResultEdges {
        hs: &mut hs_r,
        rs: &mut rs_r,
        face_ptr: &mut face_ptr_r,
    };

    let side_p = SourceSide {
        i03: &w.i03,
        hs: &mp.hs,
        vid2r: &vid_p2r,
        fid2r: fid_p2r,
        fwd: true,
    };
    let side_q = SourceSide {
        i03: &w.i30,
        hs: &mq.hs,
        vid2r: &vid_q2r,
        fid2r: fid_q2r,
        fwd: false,
    };

    append_partial_edges(
        &side_p,
        &mp.ps,
        &ps_r,
        &mut out,
        &mut pt_p,
        &mut whole_flag_p,
    );
    append_partial_edges(
        &side_q,
        &mq.ps,
        &ps_r,
        &mut out,
        &mut pt_q,
        &mut whole_flag_q,
    );

    append_new_edges(&ps_r, &fid_pq2r, mp.nf, &mut pt_new, &mut out);

    append_whole_edges(&side_p, &whole_flag_p, &mut out);
    append_whole_edges(&side_q, &whole_flag_q, &mut out);

    Boolean45 {
        ps: ps_r,
        ns: ns_r,
        hs: hs_r,
        rs: rs_r,
        hid_per_f,
        nv_from_p: nv_rp as usize,
        nv_from_q: nv_rq as usize,
    }
}
