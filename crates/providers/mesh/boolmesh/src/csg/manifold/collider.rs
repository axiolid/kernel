//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

use crate::csg::bounds::{union_bbs, Aabb, BBox, BPos, QueryShape};
use crate::csg::{Real, Vec3};

pub const K_NO_CODE: u32 = 0xFFFFFFFF;
const K_INITIAL_LENGTH: i32 = 128;
const K_LENGTH_MULTIPLE: i32 = 4;
const K_ROOT: i32 = 1;

fn spread_bits_3(v: u32) -> u32 {
    assert!(v <= 1023);
    let mut v = v;
    v = 0xFF0000FFu32 & v.wrapping_mul(0x00010001u32);
    v = 0x0F00F00Fu32 & v.wrapping_mul(0x00000101u32);
    v = 0xC30C30C3u32 & v.wrapping_mul(0x00000011u32);
    v = 0x49249249u32 & v.wrapping_mul(0x00000005u32);
    v
}

pub fn morton_code(p: &Vec3, bb: &BBox) -> u32 {
    if p.x.is_nan() {
        return K_NO_CODE;
    }
    let mut xyz = (p - bb.min) / (bb.max - bb.min);
    xyz = (1024. * xyz)
        .max(Vec3::ZERO)
        .min(Vec3::new(1023., 1023., 1023.));
    let x = spread_bits_3(xyz.x as u32);
    let y = spread_bits_3(xyz.y as u32);
    let z = spread_bits_3(xyz.z as u32);
    x * 4 + y * 2 + z
}

fn node2intl(node: i32) -> Option<i32> {
    if node % 2 == 1 {
        Some((node - 1) / 2)
    } else {
        None
    }
}
fn node2leaf(node: i32) -> Option<i32> {
    if node % 2 == 0 {
        Some(node / 2)
    } else {
        None
    }
}
fn intl2node(intl: i32) -> i32 {
    intl * 2 + 1
}
fn leaf2node(leaf: i32) -> i32 {
    leaf * 2
}
fn prefix_length(a: u32, b: u32) -> u32 {
    (a ^ b).leading_zeros()
}

struct RadixTree<'a> {
    parent: &'a mut [i32],
    children: &'a mut [(i32, i32)],
    leaf_morton: &'a [u32],
}

impl<'a> RadixTree<'a> {
    fn prefix_length(&self, i: i32, j: i32) -> i32 {
        if j < 0 || j >= self.leaf_morton.len() as i32 {
            return -1;
        }
        let lmi = self.leaf_morton[i as usize];
        let lmj = self.leaf_morton[j as usize];
        if lmi == lmj {
            return 32 + prefix_length(i as u32, j as u32) as i32;
        }
        prefix_length(lmi, lmj) as i32
    }

    fn range_end(&self, i: i32) -> i32 {
        let mut dir = self.prefix_length(i, i + 1) - self.prefix_length(i, i - 1);
        dir = if dir > 0 {
            1
        } else if dir < 0 {
            -1
        } else {
            0
        };

        let common = self.prefix_length(i, i - dir);
        let mut max = K_INITIAL_LENGTH;
        while self.prefix_length(i, i + dir * max) > common {
            max *= K_LENGTH_MULTIPLE;
        }

        // compute precise range length with binary search
        let mut len = 0;
        let mut stp = max / 2;
        while stp > 0 {
            if self.prefix_length(i, i + dir * (len + stp)) > common {
                len += stp;
            }
            stp /= 2;
        }
        i + dir * len
    }

    fn find_split(&self, bgn: i32, end: i32) -> i32 {
        let common = self.prefix_length(bgn, end);
        // Find the furthest object that shares more than common_prefix bits
        // with the first one, using binary search.
        let mut split = bgn;
        let mut step = end - bgn;

        loop {
            step = (step + 1) >> 1; // divide by 2, rounding up
            let new_split = split + step;
            if new_split < end && self.prefix_length(bgn, new_split) > common {
                split = new_split;
            }
            if step <= 1 {
                break;
            }
        }

        split
    }

    fn op(&mut self, intl: i32) {
        let mut bgn = intl;
        let mut end = self.range_end(bgn);
        if bgn > end {
            std::mem::swap(&mut bgn, &mut end);
        }

        let mut s = self.find_split(bgn, end);
        let child1 = if s == bgn { leaf2node(s) } else { intl2node(s) };
        s += 1;
        let child2 = if s == end { leaf2node(s) } else { intl2node(s) };

        self.children[intl as usize] = (child1, child2);
        self.parent[child1 as usize] = intl2node(intl);
        self.parent[child2 as usize] = intl2node(intl);
    }
}

fn build_internal_boxes(
    node_bb: &mut [Aabb],
    counter: &mut [i32],
    node_parent: &[i32],
    intl_children: &[(i32, i32)],
    leaf: i32,
) {
    let mut node = leaf2node(leaf);
    let mut flag = false;
    loop {
        if flag && node == K_ROOT {
            return;
        }
        node = node_parent[node as usize];
        let intl_idx = node2intl(node).unwrap() as usize;
        let c = counter[intl_idx];
        counter[intl_idx] += 1;
        if c == 0 {
            return;
        }
        node_bb[node as usize] = union_bbs(
            &node_bb[intl_children[intl_idx].0 as usize],
            &node_bb[intl_children[intl_idx].1 as usize],
        );
        flag = true;
    }
}

/// A uniform grid over x/y, for queries that only test x and y.
///
/// `winding03`'s point-in-polygon test reads only x/y (see
/// `BPos::overlaps_node`); the earlier attempt at this problem built a
/// second BVH sorted the same way and found that its O(n log n)
/// construction cost as much as the traversal it saved, because the tree
/// is built once and used for exactly one pass of queries -- there is no
/// second call to amortize it against.
///
/// A grid is O(n) to build: bucket each face's xy footprint into the
/// cells it overlaps, no sort and no tree. It trades the BVH's adaptivity
/// for a resolution tuned to the expected density -- appropriate here
/// because a subdivided icosphere's faces are close to uniform in size,
/// which is exactly the case a uniform grid is suited to and a tree
/// derives no extra benefit from adapting to.
#[derive(Clone, Debug)]
pub struct PlanarGrid {
    min: Vec3,
    cell: Real,
    dim: u32,
    /// CSR row offsets: cell `(cx, cy)` at index `cy * dim + cx` owns
    /// `bucket[start[i]..start[i + 1]]`.
    start: Vec<u32>,
    bucket: Vec<u32>,
    /// The exact box the query must still confirm against: cell
    /// membership is a broad-phase over-approximation for faces whose
    /// bbox straddles a cell boundary, not the collision test itself.
    face_bb: Vec<BBox>,
}

impl PlanarGrid {
    pub fn new(face_bb: &[BBox], bb: &BBox) -> Self {
        let n = face_bb.len();
        // Target ~2 faces per cell: dense enough that a query's candidate
        // list stays short, sparse enough that few faces straddle more
        // than one cell.
        let dim = ((n as f64 / 2.0).sqrt().ceil() as u32).max(1);
        let sx = (bb.max.x - bb.min.x).max(Real::EPSILON);
        let sy = (bb.max.y - bb.min.y).max(Real::EPSILON);
        let cell = (sx.max(sy)) / dim as Real;
        let min = bb.min;

        let cell_of = |x: Real, y: Real| -> (u32, u32) {
            let cx = (((x - min.x) / cell) as i64).clamp(0, dim as i64 - 1) as u32;
            let cy = (((y - min.y) / cell) as i64).clamp(0, dim as i64 - 1) as u32;
            (cx, cy)
        };

        let n_cells = (dim * dim) as usize;
        let mut count = vec![0u32; n_cells + 1];
        // Pass 1: how many (face, cell) pairs exist, per cell -- a face
        // whose bbox spans several cells is counted once per cell it
        // touches, matching how it will be looked up from any of them.
        let mut ranges = Vec::with_capacity(n);
        for fb in face_bb {
            let (cx0, cy0) = cell_of(fb.min.x, fb.min.y);
            let (cx1, cy1) = cell_of(fb.max.x, fb.max.y);
            for cy in cy0..=cy1 {
                for cx in cx0..=cx1 {
                    count[(cy * dim + cx) as usize + 1] += 1;
                }
            }
            ranges.push((cx0, cy0, cx1, cy1));
        }
        for i in 0..n_cells {
            count[i + 1] += count[i];
        }
        let start = count.clone();

        let mut bucket = vec![0u32; *count.last().unwrap() as usize];
        let mut cursor = start.clone();
        for (i, &(cx0, cy0, cx1, cy1)) in ranges.iter().enumerate() {
            for cy in cy0..=cy1 {
                for cx in cx0..=cx1 {
                    let c = (cy * dim + cx) as usize;
                    bucket[cursor[c] as usize] = i as u32;
                    cursor[c] += 1;
                }
            }
        }

        PlanarGrid {
            min,
            cell,
            dim,
            start,
            bucket,
            face_bb: face_bb.to_vec(),
        }
    }

    /// Test every point query against the grid.
    ///
    /// Specific to `BPos` rather than generic over `QueryShape`: this is a
    /// spatial-hash lookup, not a tree descent, and the only caller
    /// (`winding03`) only ever has points. Forcing it through the shared
    /// trait would buy genericity nothing has asked for.
    pub fn collision<F>(&self, queries: &[BPos], record: &mut F)
    where
        F: FnMut(usize, usize),
    {
        let hi_x = self.min.x + self.cell * self.dim as Real;
        let hi_y = self.min.y + self.cell * self.dim as Real;
        for q in queries {
            let Some(qid) = q.id else { continue };
            let px = q.pos.x as Real;
            let py = q.pos.y as Real;
            if px < self.min.x || py < self.min.y || px > hi_x || py > hi_y {
                continue;
            }
            let cx = (((px - self.min.x) / self.cell) as i64).clamp(0, self.dim as i64 - 1) as u32;
            let cy = (((py - self.min.y) / self.cell) as i64).clamp(0, self.dim as i64 - 1) as u32;
            let c = (cy * self.dim + cx) as usize;
            for &f in &self.bucket[self.start[c] as usize..self.start[c + 1] as usize] {
                let fb = &self.face_bb[f as usize];
                if q.overlaps_node(&fb.into()) {
                    record(qid, f as usize);
                }
            }
        }
    }
}

#[derive(Clone, Debug)]
pub struct MortonCollider {
    pub node_bb: Vec<Aabb>,
    pub node_parent: Vec<i32>,
    pub intl_children: Vec<(i32, i32)>,
}

impl MortonCollider {
    fn num_intl(&self) -> usize {
        self.intl_children.len()
    }
    fn num_leaf(&self) -> usize {
        if self.intl_children.is_empty() {
            0
        } else {
            self.num_intl() + 1
        }
    }

    fn update_boxes(&mut self, leaf_bb: &[BBox]) {
        for (i, box_val) in leaf_bb.iter().enumerate() {
            self.node_bb[i * 2] = box_val.into();
        }
        let mut counter: Vec<i32> = vec![0; self.num_intl()];
        for i in 0..self.num_leaf() {
            build_internal_boxes(
                &mut self.node_bb,
                &mut counter,
                &self.node_parent,
                &self.intl_children,
                i as i32,
            );
        }
    }

    pub fn new(leaf_bb: &[BBox], leaf_morton: &[u32]) -> Self {
        let n_intl = leaf_bb.len() - 1;
        let n_node = 2 * leaf_bb.len() - 1;
        let mut node_parent = vec![-1; n_node];
        let mut intl_children = vec![(0, 0); n_intl];
        let mut tree = RadixTree {
            parent: &mut node_parent,
            children: &mut intl_children,
            leaf_morton,
        };

        for i in 0..n_intl {
            tree.op(i as i32);
        }

        let mut res = MortonCollider {
            node_bb: vec![Aabb::empty(); n_node],
            node_parent,
            intl_children,
        };

        res.update_boxes(leaf_bb);
        res
    }

    /// Test every query against the tree, recording (query id, leaf) hits.
    ///
    /// Generic over the query shape rather than taking a `Query` enum: each
    /// caller passes a homogeneous slice, so this monomorphises into one
    /// traversal with the overlap test inlined and no per-node branch on
    /// the variant.
    pub fn collision<Q, F>(&self, queries: &[Q], record: &mut F)
    where
        Q: QueryShape,
        F: FnMut(usize, usize),
    {
        for (i, query) in queries.iter().enumerate() {
            find_collisions(query, i, &self.node_bb, &self.intl_children, record, false)
        }
    }
}

fn find_collisions<Q, F>(
    query: &Q,
    query_idx: usize,
    node_bb: &[Aabb],
    children: &[(i32, i32)],
    record: &mut F,
    self_collision: bool,
) where
    Q: QueryShape,
    F: FnMut(usize, usize),
{
    // The query's id is loop-invariant: upstream re-matched the enum and
    // re-read `id` at every leaf hit. Nothing to record without one, so
    // an id-less query can skip the descent entirely.
    let Some(query_id) = query.id() else {
        return;
    };

    // depth-first search
    let mut stack = [0i32; 64];
    let mut top = -1i32;
    let mut node = K_ROOT;

    loop {
        let intl = node2intl(node).unwrap();
        let (c1, c2) = children[intl as usize];

        let mut visit = |n: i32| -> bool {
            if !query.overlaps_node(&node_bb[n as usize]) {
                return false;
            }
            match node2leaf(n) {
                Some(leaf) => {
                    if !self_collision || leaf != query_idx as i32 {
                        record(query_id, leaf as usize);
                    }
                    // A leaf is never descended into.
                    false
                }
                // Overlapping and internal: traverse.
                None => true,
            }
        };

        let traverse1 = visit(c1);
        let traverse2 = visit(c2);

        if !traverse1 && !traverse2 {
            if top < 0 {
                break;
            }
            node = stack[top as usize];
            top -= 1;
        } else {
            node = if traverse1 { c1 } else { c2 };
            if traverse1 && traverse2 {
                top += 1;
                stack[top as usize] = c2;
            }
        }
    }
}
