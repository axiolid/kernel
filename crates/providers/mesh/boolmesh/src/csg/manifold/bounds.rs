//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

use crate::csg::{Real, Vec2, Vec3};

#[derive(Clone, Copy, Debug)]
pub struct BBox {
    pub id: Option<usize>,
    pub min: Vec3,
    pub max: Vec3,
}

#[derive(Clone, Debug)]
pub struct BPos {
    pub id: Option<usize>,
    pub pos: Vec2,
}

impl BBox {
    pub fn default() -> Self {
        BBox {
            id: None,
            min: Vec3::MAX,
            max: Vec3::MIN,
        }
    }

    pub fn new(id: Option<usize>, pts: &[Vec3]) -> Self {
        let mut b = BBox {
            id,
            min: Vec3::MAX,
            max: Vec3::MIN,
        };
        for pt in pts {
            b.union(pt);
        }
        b
    }

    pub fn size(&self) -> Vec3 {
        self.max - self.min
    }

    pub fn scale(&self) -> Real {
        let s = self.size();
        s.x.abs().max(s.y.abs()).max(s.z.abs())
    }

    pub fn union(&mut self, p: &Vec3) {
        if p.x.is_nan() {
            return;
        }
        self.min = self.min.min(*p);
        self.max = self.max.max(*p);
    }

    pub fn longest_dim(&self) -> usize {
        let s = self.size();
        if s.x > s.y && s.x > s.z {
            0
        } else if s.y > s.z {
            1
        } else {
            2
        }
    }
}

pub fn union_bbs(b0: &Aabb, b1: &Aabb) -> Aabb {
    Aabb {
        min: b0.min.min(b1.min),
        max: b0.max.max(b1.max),
    }
}

/// An axis-aligned box with no identity attached.
///
/// The BVH stores one of these per node, and `node_bb` is the array the
/// traversal walks. `BBox` carries an `Option<usize>` id so a QUERY can
/// name itself to the recorder; a tree node never uses it, but paid 16
/// bytes of every cache line for it -- a third of the struct, always
/// `None`. At 81920 triangles per operand that is 163839 nodes, so the
/// node array shrinks from 10.0 MiB to 7.5 MiB. The traversal is
/// memory-bound, so those bytes are the cost.
#[derive(Clone, Copy, Debug)]
pub struct Aabb {
    pub min: Vec3,
    pub max: Vec3,
}

impl Aabb {
    /// An empty box that absorbs any point on the first `union`.
    pub fn empty() -> Self {
        Aabb {
            min: Vec3::MAX,
            max: Vec3::MIN,
        }
    }
}

impl From<&BBox> for Aabb {
    fn from(b: &BBox) -> Self {
        Aabb {
            min: b.min,
            max: b.max,
        }
    }
}

/// A query shape, resolved once per traversal instead of once per node.
///
/// `Query` is an enum, so `BBox::overlaps` had to branch on the variant
/// at every node of every descent -- for a value that cannot change
/// during a traversal, since each `collision` call passes a homogeneous
/// slice. Making the shape a type parameter lets the compiler emit one
/// specialised traversal per kind: the branch disappears, the overlap
/// test inlines, and the query's id is read once rather than per hit.
pub trait QueryShape {
    /// Whether this query overlaps a node's bounding box.
    fn overlaps_node(&self, bb: &Aabb) -> bool;
    /// The caller's index for this query, if it carries one.
    fn id(&self) -> Option<usize>;
}

impl QueryShape for BBox {
    #[inline]
    fn overlaps_node(&self, bb: &Aabb) -> bool {
        bb.min.cmple(self.max).all() && bb.max.cmpge(self.min).all()
    }

    #[inline]
    fn id(&self) -> Option<usize> {
        self.id
    }
}

impl QueryShape for BPos {
    #[inline]
    fn overlaps_node(&self, bb: &Aabb) -> bool {
        // Only the xy axes are evaluated, matching the enum arm this
        // replaces: a point query is a 2-D test against a 3-D node box.
        bb.min.x <= self.pos.x
            && bb.min.y <= self.pos.y
            && bb.max.x >= self.pos.x
            && bb.max.y >= self.pos.y
    }

    #[inline]
    fn id(&self) -> Option<usize> {
        self.id
    }
}
