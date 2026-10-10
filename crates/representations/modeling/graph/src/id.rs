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

#[cfg(feature = "serde")]
std::thread_local! {
    static WIRE_BRAND: core::cell::Cell<Option<GraphId>> = const { core::cell::Cell::new(None) };
}

/// The brand node references take while a wire payload is read (ADR 0085).
///
/// A [`NodeId`] on the wire is only an index; its brand is the graph being
/// rebuilt from the payload, set for the duration of the read by
/// [`WireBrand::enter`]. Outside that scope there is no graph to brand a
/// reference with, so deserialising a `NodeId` refuses rather than inventing
/// one.
#[cfg(feature = "serde")]
pub(crate) struct WireBrand {
    previous: Option<GraphId>,
}

#[cfg(feature = "serde")]
impl WireBrand {
    /// Brand every `NodeId` read on this thread with `graph` until the
    /// returned guard drops.
    pub(crate) fn enter(graph: GraphId) -> Self {
        Self {
            previous: WIRE_BRAND.with(|brand| brand.replace(Some(graph))),
        }
    }
}

#[cfg(feature = "serde")]
impl Drop for WireBrand {
    fn drop(&mut self) {
        WIRE_BRAND.with(|brand| brand.set(self.previous));
    }
}

/// A node reference is written as its zero-based index (ADR 0085).
#[cfg(feature = "serde")]
impl serde::Serialize for NodeId {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_u32(self.index)
    }
}

/// A node reference reads as an index into the graph payload being read.
///
/// It is refused outside a graph payload: a handle carries the brand of the
/// graph that owns it, and only the graph reader knows that graph.
#[cfg(feature = "serde")]
impl<'de> serde::Deserialize<'de> for NodeId {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let index = <u32 as serde::Deserialize>::deserialize(deserializer)?;
        match WIRE_BRAND.with(core::cell::Cell::get) {
            Some(graph) => Ok(Self { graph, index }),
            None => Err(serde::de::Error::custom(
                "a node reference is only meaningful inside a geometry graph payload",
            )),
        }
    }
}
