//! Welding a shell through the corners its faces insert into shared edges
//! (#265).
//!
//! The certified clipper triangulates a planar face under
//! [`PinchPolicy::Accept`](axiolid_construct::profile::PinchPolicy): a
//! corner of the face lying inside one of the face's own ring edges (a
//! pocket rim touching the face's outer edge, two holes touching) is
//! inserted into that edge, so the face's triangles run along the edge in
//! two pieces. In a closed shell the edge is shared with a neighbouring
//! face, which knows nothing of the corner and keeps the edge whole: the
//! mesh has a T-junction, two boundary edges on one side and one on the
//! other, and is no longer closed.
//!
//! So the corner is inserted into every face's copy of the edge, before
//! anything is triangulated:
//!
//! 1. Each planar face is projected exactly as it will be triangulated,
//!    and the corners its clipper would insert are found by the clipper's
//!    own exact tests ([`crate::planar::edge_touches`]). Each is recorded
//!    against its edge, keyed by the edge's two shared corners regardless
//!    of direction.
//! 2. Every face then has each recorded corner inserted into each ring
//!    edge it uses, in order along the edge, and is triangulated with it.
//!    The touching face now sees the corner twice at one point, a pinch it
//!    accepts and triangulates the same way; the neighbour sees one more
//!    corner on its straight edge, an ordinary ring corner.
//!
//! Both faces then run along the edge in the same pieces, each piece used
//! once from each side, as a closed two-manifold needs. Corners are keyed
//! by their shared identity (a B-rep's topological vertex, an authored
//! polygon's position index), so the weld introduces no new position and
//! moves none: the corner was already a corner of the shell. The volume
//! changes by nothing more than rounding, since the neighbour gains a
//! corner on its own edge.
//!
//! A face that cannot take the corner keeps the edge whole: a curved B-rep
//! face, which samples its edges itself. The mesh is then still open
//! there, and the B-rep path's closure check reports it
//! ([`axiolid_mesh_compile_contract::MeshClosure::OpenSolid`]) rather than
//! claiming a solid.

use std::collections::HashMap;
use std::hash::Hash;

use axiolid_core::{Scalar, Vec3};

use crate::planar::EdgeTouch;

/// The corners recorded inside each shared edge, keyed by the edge's
/// corners in ascending order.
#[derive(Debug, Clone)]
pub(crate) struct EdgeSplits<K> {
    inside: HashMap<(K, K), Vec<(K, Vec3)>>,
}

impl<K: Copy + Eq + Hash + Ord> EdgeSplits<K> {
    pub(crate) fn new() -> Self {
        Self {
            inside: HashMap::new(),
        }
    }

    /// Whether no corner lies inside any edge.
    pub(crate) fn is_empty(&self) -> bool {
        self.inside.is_empty()
    }

    /// Record the touches of one face whose corners, in the order its
    /// projection lists them, are `corners`.
    pub(crate) fn record(&mut self, corners: &[(K, Vec3)], touches: &[EdgeTouch]) {
        for touch in touches {
            let (a, b) = (corners[touch.from].0, corners[touch.to].0);
            let corner = corners[touch.corner];
            if a == b || corner.0 == a || corner.0 == b {
                continue;
            }
            let list = self.inside.entry((a.min(b), a.max(b))).or_default();
            if !list.iter().any(|&(k, _)| k == corner.0) {
                list.push(corner);
            }
        }
    }

    /// Insert every recorded corner into each ring edge that has one, in
    /// order from the edge's start to its end.
    pub(crate) fn split<T: Copy>(
        &self,
        rings: &mut [Vec<T>],
        corner: impl Fn(T) -> (K, Vec3),
        make: impl Fn(K, Vec3) -> T,
    ) {
        if self.is_empty() {
            return;
        }
        for ring in rings.iter_mut() {
            let len = ring.len();
            let mut out: Vec<T> = Vec::with_capacity(len);
            for i in 0..len {
                let (from, to) = (corner(ring[i]), corner(ring[(i + 1) % len]));
                out.push(ring[i]);
                let Some(list) = self.inside.get(&(from.0.min(to.0), from.0.max(to.0))) else {
                    continue;
                };
                let direction = to.1 - from.1;
                let mut along: Vec<(Scalar, K, Vec3)> = list
                    .iter()
                    .filter(|&&(k, _)| k != from.0 && k != to.0)
                    .map(|&(k, p)| ((p - from.1).dot(direction), k, p))
                    .collect();
                along.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
                out.extend(along.into_iter().map(|(_, k, p)| make(k, p)));
            }
            *ring = out;
        }
    }
}
