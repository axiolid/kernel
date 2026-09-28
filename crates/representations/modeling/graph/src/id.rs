//! Typed graph identity.

use core::fmt;
use core::sync::atomic::{AtomicU64, Ordering};

static NEXT_GRAPH_ID: AtomicU64 = AtomicU64::new(1);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct GraphId(u64);

impl GraphId {
    pub(crate) fn fresh() -> Self {
        let value = NEXT_GRAPH_ID
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |current| {
                current.checked_add(1)
            })
            .expect("geometry graph identity space exhausted");
        Self(value)
    }
}

/// Stable index owned by one immutable [`crate::GeometryGraph`].
///
/// A handle carries its graph's brand. Another builder refuses it
/// ([`crate::GraphError::ForeignReference`]) and another graph's
/// [`crate::GeometryGraph::get`] returns `None`, so a graph that references
/// a node it does not own cannot be built at all rather than being detected
/// later.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeId {
    graph: GraphId,
    index: u32,
}

impl NodeId {
    pub(crate) fn from_index(graph: GraphId, index: usize) -> Self {
        Self {
            graph,
            index: u32::try_from(index).expect("geometry graph exceeds u32 capacity"),
        }
    }

    pub(crate) fn belongs_to(self, graph: GraphId) -> bool {
        self.graph == graph
    }

    /// Zero-based graph index.
    pub const fn index(self) -> usize {
        self.index as usize
    }
}

impl fmt::Display for NodeId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "geometry#{}", self.index)
    }
}
