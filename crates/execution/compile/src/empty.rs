//! What a boolean with an empty operand is (#326).
//!
//! A half-space clip can leave nothing of its subject, and a boolean of such
//! an empty intermediate with another operand is decided by an identity, not
//! by a cut: `S - 0 = S`, `0 - T = 0`, `S & 0 = 0`, `S | 0 = S` and
//! `S ^ 0 = S` (and `0 | T = T`, `0 ^ T = T`). Both compilers read this one
//! table, so the empty case cannot drift between them: the mesh compiler
//! holds an empty mesh, the exact compiler an empty body (`ExactBRep` itself
//! cannot be empty, so only an empty root is refused there).
//!
//! The table covers the four operators the model has. An unknown operator
//! has no identity, and its caller refuses it as it does when neither
//! operand is empty.

use axiolid_core::BooleanOperator;

/// What a boolean with an empty operand is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Identity {
    /// The subject, unchanged.
    Subject,
    /// The tool, unchanged.
    Tool,
    /// Nothing.
    Empty,
}

/// The identity of `operator` when an operand is empty, or `None` when
/// neither is, or the operator is not one of the four.
pub(crate) fn identity(
    operator: BooleanOperator,
    subject_empty: bool,
    tool_empty: bool,
) -> Option<Identity> {
    if !subject_empty && !tool_empty {
        return None;
    }
    match operator {
        BooleanOperator::Difference => Some(if subject_empty {
            Identity::Empty
        } else {
            Identity::Subject
        }),
        BooleanOperator::Intersection => Some(Identity::Empty),
        BooleanOperator::Union | BooleanOperator::SymmetricDifference => {
            Some(if subject_empty && tool_empty {
                Identity::Empty
            } else if subject_empty {
                Identity::Tool
            } else {
                Identity::Subject
            })
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::{identity, Identity};
    use axiolid_core::BooleanOperator::{Difference, Intersection, SymmetricDifference, Union};

    #[test]
    fn the_table_is_every_operator_against_every_empty_operand() {
        use Identity::{Empty, Subject, Tool};
        // (operator, subject empty, tool empty) -> identity.
        let table = [
            (Difference, true, false, Empty),
            (Difference, false, true, Subject),
            (Difference, true, true, Empty),
            (Intersection, true, false, Empty),
            (Intersection, false, true, Empty),
            (Intersection, true, true, Empty),
            (Union, true, false, Tool),
            (Union, false, true, Subject),
            (Union, true, true, Empty),
            (SymmetricDifference, true, false, Tool),
            (SymmetricDifference, false, true, Subject),
            (SymmetricDifference, true, true, Empty),
        ];
        for (operator, subject, tool, expected) in table {
            assert_eq!(
                identity(operator, subject, tool),
                Some(expected),
                "{operator:?} with subject empty {subject}, tool empty {tool}"
            );
        }
    }

    #[test]
    fn two_solid_operands_have_no_identity() {
        for operator in [Difference, Intersection, Union, SymmetricDifference] {
            assert_eq!(identity(operator, false, false), None, "{operator:?}");
        }
    }
}
