//! Typed arena handles prevent mixing vertex, edge, face, and shell indices.

use core::fmt;

macro_rules! topology_id {
    ($name:ident, $label:literal) => {
        #[doc = concat!("Stable handle into a B-rep ", $label, " arena.")]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        #[cfg_attr(
            feature = "serde",
            derive(serde::Serialize, serde::Deserialize),
            serde(transparent)
        )]
        pub struct $name(u32);

        impl $name {
            pub(crate) fn from_index(index: usize) -> Self {
                Self(u32::try_from(index).expect("topology arena exceeds u32 capacity"))
            }

            /// Zero-based arena index.
            pub const fn index(self) -> usize {
                self.0 as usize
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, concat!($label, "#{}"), self.0)
            }
        }
    };
}

topology_id!(VertexId, "vertex");
topology_id!(EdgeId, "edge");
topology_id!(LoopId, "loop");
topology_id!(FaceId, "face");
topology_id!(ShellId, "shell");
topology_id!(SolidId, "solid");
