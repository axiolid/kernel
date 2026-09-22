// SPDX-License-Identifier: MPL-2.0

//! Typed arena handles for the planar arrangement.
//!
//! Mirrors the convention in `axiolid-topology`: a bare `usize` cannot be
//! passed where a vertex, half-edge, or face is expected, because the three
//! are different types. The arrangement's handles are deliberately NOT the
//! same types as the B-rep's -- a planar face and a B-rep face are different
//! concepts, and letting one stand in for the other would be a bug the
//! compiler should catch.

use core::fmt;

macro_rules! arrangement_id {
    ($name:ident, $label:literal) => {
        #[doc = concat!("Stable handle into the arrangement's ", $label, " arena.")]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(u32);

        impl $name {
            pub(crate) fn from_index(index: usize) -> Self {
                Self(u32::try_from(index).expect("arrangement arena exceeds u32 capacity"))
            }

            /// Zero-based arena index.
            #[must_use]
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

arrangement_id!(VertexId, "vertex");
arrangement_id!(HalfEdgeId, "halfedge");
arrangement_id!(FaceId, "face");

impl FaceId {
    /// The unbounded face is always arena slot zero, created with the
    /// arrangement itself. Fixing it by construction means `outer_face()`
    /// needs no bookkeeping and cannot drift.
    pub(crate) const fn outer() -> Self {
        Self(0)
    }
}
