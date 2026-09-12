//! Persistent names for exact B-rep faces and edges.
//!
//! [`crate::SurfaceId`] and [`axiolid_topology::FaceId`] are arena positions.
//! They answer "which slot" and are only meaningful inside one assembled
//! value: rebuild the catalog and slot 7 is a different face. That is fine
//! for validation, which never outlives assembly, and useless for anything
//! that must refer to a face *across* an operation -- selecting an edge to
//! fillet, carrying a material through a boolean, or re-applying a feature
//! after an upstream edit.
//!
//! A [`FaceName`] instead records WHERE THE FACE CAME FROM. It is derived
//! from the generating inputs, so the same construction names the same face
//! no matter how the arenas are packed, and an operation that splits a face
//! can say which face each fragment came from.
//!
//! # What a name is not
//!
//! A name is not a guarantee that the face still exists, and not a claim
//! that two equal names are geometrically identical. It is a statement about
//! provenance: "this face was produced by that part of that input". Callers
//! that need existence must look the name up and handle absence.

use core::fmt;

/// Which side of a swept profile a face came from.
///
/// The cap variants carry no index because a sweep has exactly one of each,
/// while side walls are per profile edge and need to say which.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SweptFace {
    /// The face closing the sweep at its start.
    StartCap,
    /// The face closing the sweep at its end.
    EndCap,
    /// The wall swept from one profile edge, indexed in profile order.
    Side(u32),
}

impl fmt::Display for SweptFace {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::StartCap => f.write_str("start-cap"),
            Self::EndCap => f.write_str("end-cap"),
            Self::Side(index) => write!(f, "side[{index}]"),
        }
    }
}

/// Which operand of a binary operation a fragment came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Operand {
    /// The left-hand operand.
    Subject,
    /// The right-hand operand.
    Tool,
}

impl fmt::Display for Operand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Subject => f.write_str("subject"),
            Self::Tool => f.write_str("tool"),
        }
    }
}

/// A persistent, structural name for one exact B-rep face.
///
/// Names compose: a boolean fragment of a filleted extrusion wall carries the
/// whole chain, so a caller can ask both "which operand" and "which original
/// profile edge" without consulting the operation that produced it.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FaceName {
    /// A face of a swept solid, named by its role in the sweep.
    Swept(SweptFace),
    /// A blend face introduced by filleting or chamfering a named edge.
    ///
    /// The name is the *edge that was blended*, not a fresh anonymous id, so
    /// re-running the feature after an edit still names the same blend.
    Blend(Box<EdgeName>),
    /// A fragment of a face of one operand of a boolean.
    ///
    /// Boolean output faces are always pieces of input faces -- the operation
    /// never invents a surface -- so every fragment can say which input face
    /// it is part of.
    Fragment {
        /// Which operand contributed the face this fragment came from.
        operand: Operand,
        /// The name that face carried in its own solid.
        source: Box<FaceName>,
    },
    /// A face whose provenance was not tracked by the producing operation.
    ///
    /// Deliberately explicit: an operation that cannot name its output says
    /// so, instead of fabricating a name that would collide with a real one.
    Anonymous(u32),
}

impl FaceName {
    /// Name a swept cap or wall.
    pub const fn swept(part: SweptFace) -> Self {
        Self::Swept(part)
    }

    /// Name the blend that replaced `edge`.
    pub fn blend(edge: EdgeName) -> Self {
        Self::Blend(Box::new(edge))
    }

    /// Name a fragment of `self` cut out by a boolean.
    pub fn fragment(self, operand: Operand) -> Self {
        Self::Fragment {
            operand,
            source: Box::new(self),
        }
    }

    /// The original face name with all boolean fragment layers stripped.
    ///
    /// Answers "what was this before any boolean touched it", which is what a
    /// material or finish assignment needs.
    pub fn origin(&self) -> &Self {
        match self {
            Self::Fragment { source, .. } => source.origin(),
            other => other,
        }
    }

    /// Whether this name, or anything it came from, is anonymous.
    ///
    /// A caller that requires full provenance checks this rather than pattern
    /// matching the outermost layer, which a fragment would hide.
    pub fn is_anonymous(&self) -> bool {
        matches!(self.origin(), Self::Anonymous(_))
    }
}

impl fmt::Display for FaceName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Swept(part) => write!(f, "{part}"),
            Self::Blend(edge) => write!(f, "blend({edge})"),
            Self::Fragment { operand, source } => write!(f, "{operand}/{source}"),
            Self::Anonymous(index) => write!(f, "anon#{index}"),
        }
    }
}

/// A persistent, structural name for one exact B-rep edge.
///
/// An edge is named by the faces that meet there rather than by its own
/// arena slot, because that is what survives a rebuild: the intersection of
/// two named faces is the same edge however the arenas are ordered.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EdgeName {
    /// The edge where two named faces meet.
    ///
    /// The pair is stored in a canonical order so that naming the edge from
    /// either side produces the same name.
    Between(Box<FaceName>, Box<FaceName>),
    /// An edge whose provenance was not tracked by the producing operation.
    Anonymous(u32),
}

impl EdgeName {
    /// Name the edge shared by two faces, order-independently.
    ///
    /// `between(a, b)` and `between(b, a)` are equal, because "the edge where
    /// these two faces meet" does not depend on which face is mentioned
    /// first. Without this an edge would have two names and a fillet applied
    /// from the other side would miss.
    pub fn between(first: FaceName, second: FaceName) -> Self {
        let (low, high) = if first <= second {
            (first, second)
        } else {
            (second, first)
        };
        Self::Between(Box::new(low), Box::new(high))
    }

    /// Whether this name, or either face it names, is anonymous.
    pub fn is_anonymous(&self) -> bool {
        match self {
            Self::Between(first, second) => first.is_anonymous() || second.is_anonymous(),
            Self::Anonymous(_) => true,
        }
    }
}

impl fmt::Display for EdgeName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Between(first, second) => write!(f, "{first}|{second}"),
            Self::Anonymous(index) => write!(f, "anon#{index}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A fragment remembers the face it was cut from, however deep.
    ///
    /// This is what carries a material or finish through a boolean: the
    /// caller asks the fragment what it originally was.
    #[test]
    fn a_fragment_reports_the_face_it_came_from() {
        let wall = FaceName::swept(SweptFace::Side(2));
        let once = wall.clone().fragment(Operand::Subject);
        let twice = once.clone().fragment(Operand::Tool);

        assert_eq!(once.origin(), &wall);
        assert_eq!(
            twice.origin(),
            &wall,
            "a fragment of a fragment still came from the original wall"
        );
    }

    /// Anonymity is inherited, so a wrapper cannot launder it.
    #[test]
    fn anonymity_survives_being_wrapped() {
        let unnamed = FaceName::Anonymous(3);
        assert!(unnamed.is_anonymous());
        assert!(
            unnamed.fragment(Operand::Subject).is_anonymous(),
            "wrapping an unnamed face must not make it look named"
        );
    }

    /// An edge naming an anonymous face is itself not fully named.
    #[test]
    fn an_edge_is_anonymous_when_either_side_is() {
        let known = FaceName::swept(SweptFace::StartCap);
        let unknown = FaceName::Anonymous(0);
        assert!(EdgeName::between(known.clone(), unknown).is_anonymous());
        assert!(!EdgeName::between(known.clone(), known).is_anonymous());
    }

    /// Distinct provenance must not collide into one name.
    #[test]
    fn different_origins_are_different_names() {
        assert_ne!(
            FaceName::swept(SweptFace::Side(1)),
            FaceName::swept(SweptFace::Side(2))
        );
        assert_ne!(
            FaceName::swept(SweptFace::StartCap),
            FaceName::swept(SweptFace::EndCap)
        );
        let wall = FaceName::swept(SweptFace::Side(1));
        assert_ne!(
            wall.clone().fragment(Operand::Subject),
            wall.fragment(Operand::Tool),
            "the same wall cut from each operand yields distinct fragments"
        );
    }
}
