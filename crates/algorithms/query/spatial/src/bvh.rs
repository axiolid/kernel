//! Deterministic median-split bounding-volume hierarchy.
//!
//! The tree stores only caller-owned keys and axis-aligned bounds. It is a
//! broad-phase structure: overlap and ray results are candidates, never an
//! assertion about exact geometry. Build and queries are serial. The
//! immutable representation is deliberately
//! provider-neutral; a parallel or GPU builder can implement the same
//! [`crate::SpatialIndex`] contract later without exposing hardware concepts.

use core::cmp::Ordering;
use core::ops::ControlFlow;
use std::collections::BinaryHeap;

use axiolid_core::{Aabb, Ray3, Scalar};

use crate::{RayHit, SpatialIndex, SpatialItem};

const LEAF_SIZE: usize = 8;

/// Observable cost counters for an allocating candidate-pair query.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SpatialQueryStats {
    /// Tree nodes whose envelope was examined.
    pub visited_nodes: usize,
    /// Leaf-item bounds tested after envelope pruning.
    pub tested_items: usize,
}

/// One conservative pair emitted by a broad-phase query.
#[derive(Debug, Clone, PartialEq)]
pub struct CandidatePair<K> {
    /// First key, ordered by accepted input position.
    pub a: K,
    /// Second key, ordered by accepted input position.
    pub b: K,
    /// Nonnegative AABB lower bound; zero means touching or overlap.
    pub lower_bound: Scalar,
}

/// Deterministically ordered conservative pair candidates plus execution cost.
#[derive(Debug, Clone, PartialEq)]
pub struct PairCandidates<K> {
    /// Pairs in accepted input order `(i, j)` where `i < j`.
    pub pairs: Vec<CandidatePair<K>>,
    /// Broad-phase work performed to obtain `pairs`.
    pub stats: SpatialQueryStats,
}

/// Nearest accepted key according to AABB lower-bound distance.
#[derive(Debug, Clone, PartialEq)]
pub struct NearestCandidate<K> {
    /// Accepted caller key.
    pub key: K,
    /// Nonnegative AABB lower bound.
    pub lower_bound: Scalar,
    /// Broad-phase work performed to find this key.
    pub stats: SpatialQueryStats,
}

#[derive(Debug)]
enum NodeKind {
    Leaf(Vec<usize>),
    Branch { left: usize, right: usize },
}

#[derive(Debug)]
struct Node {
    bounds: Aabb,
    kind: NodeKind,
}

/// Immutable median-split AABB hierarchy over opaque caller keys.
///
/// Invalid input bounds are not silently indexed: empty and non-finite boxes
/// are rejected during construction and counted by [`Self::rejected_items`].
/// Accepted keys keep their input position, which makes pair output stable even
/// though node layout is optimized for pruning.
#[derive(Debug)]
pub struct Bvh<K> {
    items: Vec<SpatialItem<K>>,
    nodes: Vec<Node>,
    root: Option<usize>,
    rejected_items: usize,
}

impl<K> Bvh<K> {
    /// Build a deterministic median-split hierarchy.
    pub fn build(items: impl IntoIterator<Item = SpatialItem<K>>) -> Self {
        let mut rejected_items = 0;
        let items = items
            .into_iter()
            .filter(|item| {
                let accepted = item.bounds.is_finite() && !item.bounds.is_empty();
                rejected_items += usize::from(!accepted);
                accepted
            })
            .collect();
        let mut tree = Self {
            items,
            nodes: Vec::new(),
            root: None,
            rejected_items,
        };
        if !tree.items.is_empty() {
            let indices = (0..tree.items.len()).collect();
            tree.root = Some(tree.build_node(indices));
        }
        tree
    }

    /// Number of accepted items.
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// Whether the hierarchy has no accepted items.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// Number of source items rejected for empty or non-finite bounds.
    pub fn rejected_items(&self) -> usize {
        self.rejected_items
    }

    /// Access an accepted item by its stable input position.
    pub fn item(&self, index: usize) -> Option<&SpatialItem<K>> {
        self.items.get(index)
    }

    fn build_node(&mut self, mut indices: Vec<usize>) -> usize {
        let bounds = union_bounds(indices.iter().map(|&index| self.items[index].bounds));
        let node_index = self.nodes.len();
        self.nodes.push(Node {
            bounds,
            kind: NodeKind::Leaf(Vec::new()),
        });
        if indices.len() <= LEAF_SIZE {
            self.nodes[node_index].kind = NodeKind::Leaf(indices);
            return node_index;
        }

        let extent = bounds.diagonal();
        let axis = if extent.x >= extent.y && extent.x >= extent.z {
            0
        } else if extent.y >= extent.z {
            1
        } else {
            2
        };
        indices.sort_unstable_by(|&left, &right| {
            component(self.items[left].bounds.center(), axis)
                .total_cmp(&component(self.items[right].bounds.center(), axis))
                .then_with(|| left.cmp(&right))
        });
        let right_indices = indices.split_off(indices.len() / 2);
        let left = self.build_node(indices);
        let right = self.build_node(right_indices);
        self.nodes[node_index].kind = NodeKind::Branch { left, right };
        node_index
    }
}

impl<K: Clone> Bvh<K> {
    /// Return pairs whose AABBs overlap with at least `min_penetration` on all
    /// axes. A zero threshold includes touching boxes.
    ///
    /// # Panics
    /// Panics when `min_penetration` is negative or non-finite; an invalid query
    /// must not masquerade as an evaluated empty result.
    /// Keys whose bounds overlap `probe`.
    ///
    /// The existing pair queries are self-joins within one tree. A clash
    /// test joins TWO meshes, so it probes one tree with the other's boxes.
    /// Without this a caller has no option but a quadratic scan.
    pub fn query_aabb(&self, probe: &Aabb, out: &mut Vec<usize>) {
        out.clear();
        let Some(root) = self.root else {
            return;
        };
        let mut stack = vec![root];
        while let Some(index) = stack.pop() {
            let node = &self.nodes[index];
            if !node.bounds.intersects(probe) {
                continue;
            }
            match &node.kind {
                NodeKind::Leaf(items) => out.extend_from_slice(items),
                NodeKind::Branch { left, right } => {
                    stack.push(*left);
                    stack.push(*right);
                }
            }
        }
    }

    pub fn overlap_pairs(&self, min_penetration: Scalar) -> PairCandidates<K> {
        assert!(
            min_penetration.is_finite() && min_penetration >= 0.0,
            "minimum penetration must be finite and non-negative"
        );
        self.collect_pairs(|left, right| {
            penetrates(left, right, min_penetration).then_some(left.gap(right))
        })
    }

    /// Return pairs whose AABB lower-bound distance is at most `max_distance`.
    ///
    /// # Panics
    /// Panics when `max_distance` is negative or non-finite; an invalid query
    /// must not masquerade as an evaluated empty result.
    pub fn pairs_within_distance(&self, max_distance: Scalar) -> PairCandidates<K> {
        assert!(
            max_distance.is_finite() && max_distance >= 0.0,
            "maximum distance must be finite and non-negative"
        );
        self.collect_pairs(|left, right| {
            let gap = left.gap(right);
            (gap <= max_distance).then_some(gap)
        })
    }

    /// Find the accepted key with the smallest AABB lower-bound distance.
    ///
    /// Equal distances resolve to the earliest accepted input item. Invalid query
    /// bounds panic rather than being misreported as an evaluated empty result.
    pub fn nearest_to(
        &self,
        query: &Aabb,
        accept: impl Fn(&K) -> bool,
    ) -> Option<NearestCandidate<K>> {
        assert!(
            query.is_finite() && !query.is_empty(),
            "nearest-neighbour query bounds must be finite and non-empty"
        );
        let root = self.root?;
        let mut stats = SpatialQueryStats::default();
        let mut pending = BinaryHeap::new();
        pending.push(NearestQueueEntry::new(
            query.gap(&self.nodes[root].bounds),
            root,
        ));
        let mut best = None;

        while let Some(entry) = pending.pop() {
            stats.visited_nodes += 1;
            if best.is_some_and(|(distance, _)| entry.distance > distance) {
                break;
            }
            match &self.nodes[entry.node].kind {
                NodeKind::Leaf(indices) => {
                    for &index in indices {
                        if !accept(&self.items[index].key) {
                            continue;
                        }
                        stats.tested_items += 1;
                        let distance = query.gap(&self.items[index].bounds);
                        if best.is_none_or(|(current, current_index)| {
                            distance < current || (distance == current && index < current_index)
                        }) {
                            best = Some((distance, index));
                        }
                    }
                }
                NodeKind::Branch { left, right } => {
                    for child in [*left, *right] {
                        let distance = query.gap(&self.nodes[child].bounds);
                        if best.is_none_or(|(current, _)| distance <= current) {
                            pending.push(NearestQueueEntry::new(distance, child));
                        }
                    }
                }
            }
        }

        best.map(|(lower_bound, index)| NearestCandidate {
            key: self.items[index].key.clone(),
            lower_bound,
            stats,
        })
    }

    fn collect_pairs(&self, matches: impl Fn(&Aabb, &Aabb) -> Option<Scalar>) -> PairCandidates<K> {
        let mut pairs = Vec::new();
        let mut stats = SpatialQueryStats::default();
        let Some(root) = self.root else {
            return PairCandidates { pairs, stats };
        };

        for index in 0..self.items.len() {
            let bounds = &self.items[index].bounds;
            let mut stack = vec![root];
            let mut matches_for_item = Vec::new();
            while let Some(node_index) = stack.pop() {
                stats.visited_nodes += 1;
                let node = &self.nodes[node_index];
                if matches(bounds, &node.bounds).is_none() {
                    continue;
                }
                match &node.kind {
                    NodeKind::Leaf(indices) => {
                        for &other_index in indices {
                            if other_index <= index {
                                continue;
                            }
                            stats.tested_items += 1;
                            if let Some(lower_bound) =
                                matches(bounds, &self.items[other_index].bounds)
                            {
                                matches_for_item.push((other_index, lower_bound));
                            }
                        }
                    }
                    NodeKind::Branch { left, right } => {
                        stack.push(*right);
                        stack.push(*left);
                    }
                }
            }
            matches_for_item.sort_unstable_by_key(|(other_index, _)| *other_index);
            pairs.extend(
                matches_for_item
                    .into_iter()
                    .map(|(other_index, lower_bound)| CandidatePair {
                        a: self.items[index].key.clone(),
                        b: self.items[other_index].key.clone(),
                        lower_bound,
                    }),
            );
        }
        PairCandidates { pairs, stats }
    }
}

impl<K> SpatialIndex<K> for Bvh<K>
where
    K: core::fmt::Debug + Send + Sync,
{
    fn visit_aabb(&self, query: &Aabb, visitor: &mut dyn FnMut(&K) -> ControlFlow<()>) {
        if query.is_empty() || !query.is_finite() {
            return;
        }
        let Some(root) = self.root else {
            return;
        };
        let mut stack = vec![root];
        while let Some(node_index) = stack.pop() {
            let node = &self.nodes[node_index];
            if !query.intersects(&node.bounds) {
                continue;
            }
            match &node.kind {
                NodeKind::Leaf(indices) => {
                    for &item_index in indices {
                        let item = &self.items[item_index];
                        if query.intersects(&item.bounds) && visitor(&item.key).is_break() {
                            return;
                        }
                    }
                }
                NodeKind::Branch { left, right } => {
                    stack.push(*right);
                    stack.push(*left);
                }
            }
        }
    }

    fn visit_ray(&self, ray: &Ray3, visitor: &mut dyn FnMut(RayHit<&K>) -> ControlFlow<()>) {
        let Some(root) = self.root else {
            return;
        };
        let Some(root_distance) = ray_aabb_entry(ray, &self.nodes[root].bounds) else {
            return;
        };
        let mut pending = BinaryHeap::new();
        pending.push(RayQueueEntry::node(root_distance, root));
        while let Some(entry) = pending.pop() {
            match entry.kind {
                RayQueueKind::Node(node_index) => match &self.nodes[node_index].kind {
                    NodeKind::Leaf(indices) => {
                        for &item_index in indices {
                            if let Some(distance) =
                                ray_aabb_entry(ray, &self.items[item_index].bounds)
                            {
                                pending.push(RayQueueEntry::item(distance, item_index));
                            }
                        }
                    }
                    NodeKind::Branch { left, right } => {
                        for child in [*left, *right] {
                            if let Some(distance) = ray_aabb_entry(ray, &self.nodes[child].bounds) {
                                pending.push(RayQueueEntry::node(distance, child));
                            }
                        }
                    }
                },
                RayQueueKind::Item(item_index) => {
                    if visitor(RayHit {
                        key: &self.items[item_index].key,
                        distance: entry.distance,
                    })
                    .is_break()
                    {
                        return;
                    }
                }
            }
        }
    }

    fn len(&self) -> usize {
        self.len()
    }
}

#[derive(Debug, Clone, Copy)]
struct NearestQueueEntry {
    distance: Scalar,
    node: usize,
}

impl NearestQueueEntry {
    const fn new(distance: Scalar, node: usize) -> Self {
        Self { distance, node }
    }
}

impl PartialEq for NearestQueueEntry {
    fn eq(&self, other: &Self) -> bool {
        self.distance == other.distance && self.node == other.node
    }
}
impl Eq for NearestQueueEntry {}
impl PartialOrd for NearestQueueEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for NearestQueueEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .distance
            .total_cmp(&self.distance)
            .then_with(|| other.node.cmp(&self.node))
    }
}

#[derive(Debug, Clone, Copy)]
enum RayQueueKind {
    Node(usize),
    Item(usize),
}

#[derive(Debug, Clone, Copy)]
struct RayQueueEntry {
    distance: Scalar,
    kind: RayQueueKind,
}

impl RayQueueEntry {
    fn node(distance: Scalar, index: usize) -> Self {
        Self {
            distance,
            kind: RayQueueKind::Node(index),
        }
    }

    fn item(distance: Scalar, index: usize) -> Self {
        Self {
            distance,
            kind: RayQueueKind::Item(index),
        }
    }

    fn order_key(self) -> (u8, usize) {
        match self.kind {
            // Nodes at a shared distance expand before items, so all same-distance
            // candidates enter the heap before stable item-index ordering applies.
            RayQueueKind::Node(index) => (0, index),
            RayQueueKind::Item(index) => (1, index),
        }
    }
}

impl PartialEq for RayQueueEntry {
    fn eq(&self, other: &Self) -> bool {
        self.distance == other.distance && self.order_key() == other.order_key()
    }
}
impl Eq for RayQueueEntry {}
impl PartialOrd for RayQueueEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for RayQueueEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .distance
            .total_cmp(&self.distance)
            .then_with(|| other.order_key().cmp(&self.order_key()))
    }
}

fn union_bounds(bounds: impl IntoIterator<Item = Aabb>) -> Aabb {
    let mut union = Aabb::empty();
    for bounds in bounds {
        union.union(&bounds);
    }
    union
}

fn component(point: axiolid_core::Point3, axis: usize) -> Scalar {
    match axis {
        0 => point.x,
        1 => point.y,
        _ => point.z,
    }
}

fn penetrates(left: &Aabb, right: &Aabb, minimum: Scalar) -> bool {
    let overlap = left.max.min(right.max) - left.min.max(right.min);
    overlap.x >= minimum && overlap.y >= minimum && overlap.z >= minimum
}

fn ray_aabb_entry(ray: &Ray3, bounds: &Aabb) -> Option<Scalar> {
    if bounds.is_empty()
        || !bounds.is_finite()
        || !ray.origin.is_finite()
        || !ray.direction.is_finite()
    {
        return None;
    }
    let mut entry = Scalar::NEG_INFINITY;
    let mut exit = Scalar::INFINITY;
    for (origin, direction, min, max) in [
        (ray.origin.x, ray.direction.x, bounds.min.x, bounds.max.x),
        (ray.origin.y, ray.direction.y, bounds.min.y, bounds.max.y),
        (ray.origin.z, ray.direction.z, bounds.min.z, bounds.max.z),
    ] {
        if direction == 0.0 {
            if origin < min || origin > max {
                return None;
            }
            continue;
        }
        let first = (min - origin) / direction;
        let second = (max - origin) / direction;
        entry = entry.max(first.min(second));
        exit = exit.min(first.max(second));
        if exit < entry {
            return None;
        }
    }
    (exit >= 0.0).then_some(entry.max(0.0))
}
