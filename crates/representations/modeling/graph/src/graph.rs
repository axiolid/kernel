//! Immutable geometry DAG and append-only builder.

use core::fmt;

use crate::{id::GraphId, BuiltInNode, GeometryNode, NodeId};

/// Invalid graph construction.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GraphError {
    /// A node handle belongs to a different graph builder.
    ForeignReference { reference: NodeId },
    /// A node referenced itself or a later/not-yet-inserted node.
    NonPriorReference { node: NodeId, reference: NodeId },
    /// A reference resolves locally but points to the wrong node family.
    InvalidReferenceType {
        /// Existing node whose family is invalid for this edge.
        reference: NodeId,
        /// Human-readable family accepted by the edge.
        expected: &'static str,
        /// Human-readable family of the referenced node.
        actual: &'static str,
    },
    /// A requested root does not exist.
    UnknownRoot { root: NodeId, node_count: usize },
    /// A master representation names a side the relation does not have.
    ContradictoryMaster {
        /// What the contradiction is, in the caller's terms.
        detail: &'static str,
    },
    /// A station or a list of stations is malformed (#241): a non-finite
    /// or negative distance, a non-finite offset, too few stations,
    /// distances that do not increase, or inconsistent tags.
    InvalidStation {
        /// What is wrong, in the caller's terms.
        detail: &'static str,
    },
}

impl fmt::Display for GraphError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ForeignReference { reference } => {
                write!(f, "{reference} belongs to another geometry graph")
            }
            Self::NonPriorReference { node, reference } => {
                write!(f, "{node} references non-prior {reference}")
            }
            Self::InvalidReferenceType {
                reference,
                expected,
                actual,
            } => write!(f, "{reference} has node type {actual}; expected {expected}"),
            Self::ContradictoryMaster { detail } => {
                write!(f, "contradictory surface-curve master: {detail}")
            }
            Self::InvalidStation { detail } => write!(f, "invalid station: {detail}"),
            Self::UnknownRoot { root, node_count } => {
                write!(f, "root {root} exceeds graph size {node_count}")
            }
        }
    }
}

impl std::error::Error for GraphError {}

/// Immutable acyclic geometry graph with one or more roots.
#[derive(Debug, Clone, PartialEq)]
pub struct GeometryGraph {
    owner: GraphId,
    nodes: Vec<GeometryNode>,
    roots: Vec<NodeId>,
}

impl Default for GeometryGraph {
    fn default() -> Self {
        Self {
            owner: GraphId::fresh(),
            nodes: Vec::new(),
            roots: Vec::new(),
        }
    }
}

impl GeometryGraph {
    /// Number of nodes.
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Whether the graph has no nodes.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// Read a node by typed handle. A handle owned by another graph returns `None`.
    pub fn get(&self, id: NodeId) -> Option<&GeometryNode> {
        if !id.belongs_to(self.owner) {
            return None;
        }
        self.nodes.get(id.index())
    }

    /// Output roots in caller-specified order.
    pub fn roots(&self) -> &[NodeId] {
        &self.roots
    }

    /// All nodes in stable topological insertion order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (NodeId, &GeometryNode)> {
        self.nodes
            .iter()
            .enumerate()
            .map(|(index, node)| (NodeId::from_index(self.owner, index), node))
    }
}

/// Append-only builder that makes cycles and dangling references unrepresentable.
#[derive(Debug)]
pub struct GeometryGraphBuilder {
    owner: Option<GraphId>,
    nodes: Vec<GeometryNode>,
}

impl Default for GeometryGraphBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl GeometryGraphBuilder {
    /// Create an empty builder.
    pub const fn new() -> Self {
        Self {
            owner: None,
            nodes: Vec::new(),
        }
    }

    /// Insert a node. Every reference must be to an earlier node.
    pub fn push(&mut self, node: GeometryNode) -> Result<NodeId, GraphError> {
        let owner = *self.owner.get_or_insert_with(GraphId::fresh);
        let id = NodeId::from_index(owner, self.nodes.len());
        let references = node.references();
        if let Some(&reference) = references
            .iter()
            .find(|reference| !reference.belongs_to(owner))
        {
            return Err(GraphError::ForeignReference { reference });
        }
        if let Some(&reference) = references
            .iter()
            .find(|reference| reference.index() >= id.index())
        {
            return Err(GraphError::NonPriorReference {
                node: id,
                reference,
            });
        }
        crate::validation::validate_reference_types(&node, &self.nodes)?;
        self.nodes.push(node);
        Ok(id)
    }

    /// Insert one canonical representation without spelling its enum variant.
    ///
    /// The accepted set is deliberately sealed; adapters must translate custom
    /// values into a built-in neutral representation before graph construction.
    pub fn push_value<T>(&mut self, value: T) -> Result<NodeId, GraphError>
    where
        T: BuiltInNode,
    {
        self.push(value.into())
    }

    /// Freeze the graph after validating roots.
    pub fn finish(self, roots: Vec<NodeId>) -> Result<GeometryGraph, GraphError> {
        let owner = self.owner.unwrap_or_else(GraphId::fresh);
        if let Some(&reference) = roots.iter().find(|root| !root.belongs_to(owner)) {
            return Err(GraphError::ForeignReference { reference });
        }
        if let Some(&root) = roots.iter().find(|root| root.index() >= self.nodes.len()) {
            return Err(GraphError::UnknownRoot {
                root,
                node_count: self.nodes.len(),
            });
        }
        Ok(GeometryGraph {
            owner,
            nodes: self.nodes,
            roots,
        })
    }
}

#[cfg(test)]
mod tests {
    use axiolid_core::Vec3;

    use super::*;
    use crate::Instance;

    const EMPTY_BUILDER: GeometryGraphBuilder = GeometryGraphBuilder::new();

    #[test]
    fn const_constructor_remains_source_compatible() {
        assert!(EMPTY_BUILDER.finish(Vec::new()).unwrap().is_empty());
    }

    #[test]
    fn insertion_order_is_topological_order() {
        let mut builder = GeometryGraphBuilder::new();
        let source = builder.push(GeometryNode::Point3(Vec3::ZERO)).unwrap();
        let instance = builder
            .push(GeometryNode::Instance(Instance {
                source,
                transform: axiolid_core::Transform3::IDENTITY,
            }))
            .unwrap();
        let graph = builder.finish(vec![instance]).unwrap();
        assert_eq!(graph.len(), 2);
        assert_eq!(graph.roots(), &[instance]);
    }

    #[test]
    fn sealed_built_in_values_have_an_ergonomic_builder_path() {
        let mut builder = GeometryGraphBuilder::new();
        let sphere = builder
            .push_value(axiolid_primitive::Primitive::Sphere { radius: 1.0 })
            .unwrap();
        let graph = builder.finish(vec![sphere]).unwrap();
        assert!(matches!(
            graph.get(sphere),
            Some(GeometryNode::Primitive(
                axiolid_primitive::Primitive::Sphere { radius: 1.0 }
            ))
        ));
    }

    #[test]
    fn handles_from_another_builder_cannot_alias_local_nodes() {
        let mut foreign_builder = GeometryGraphBuilder::new();
        let foreign = foreign_builder
            .push(GeometryNode::Point3(Vec3::ZERO))
            .unwrap();

        let mut builder = GeometryGraphBuilder::new();
        let local = builder.push(GeometryNode::Point3(Vec3::ZERO)).unwrap();
        let error = builder
            .push(GeometryNode::Instance(Instance {
                source: foreign,
                transform: axiolid_core::Transform3::IDENTITY,
            }))
            .unwrap_err();
        assert!(matches!(
            error,
            GraphError::ForeignReference { reference } if reference == foreign
        ));
        let error = builder.finish(vec![foreign]).unwrap_err();
        assert!(matches!(
            error,
            GraphError::ForeignReference { reference } if reference == foreign
        ));

        let graph = foreign_builder.finish(vec![foreign]).unwrap();
        assert!(graph.get(local).is_none());
    }

    #[test]
    fn semantic_reference_types_are_validated_before_insertion() {
        let mut builder = GeometryGraphBuilder::new();
        let point = builder.push(GeometryNode::Point3(Vec3::ZERO)).unwrap();
        let error = builder
            .push(GeometryNode::SolidOperation(
                crate::SolidOperation::Extrusion {
                    profile: point,
                    direction: Vec3::Z,
                    depth: 1.0,
                },
            ))
            .unwrap_err();
        assert!(matches!(
            error,
            GraphError::InvalidReferenceType {
                reference,
                expected: "profile",
                actual: "point3",
            } if reference == point
        ));
    }

    #[test]
    fn instance_nodes_preserve_their_source_reference_family() {
        let mut builder = GeometryGraphBuilder::new();
        let point = builder.push(GeometryNode::Point3(Vec3::ZERO)).unwrap();
        let instance = builder
            .push(GeometryNode::Instance(Instance {
                source: point,
                transform: axiolid_core::Transform3::IDENTITY,
            }))
            .unwrap();

        let error = builder
            .push(GeometryNode::SolidOperation(
                crate::SolidOperation::Boolean {
                    left: instance,
                    right: instance,
                    operator: axiolid_core::BooleanOperator::Union,
                },
            ))
            .unwrap_err();
        let GraphError::InvalidReferenceType {
            reference,
            expected,
            actual,
        } = error
        else {
            panic!("unexpected graph error: {error:?}");
        };
        assert_eq!(reference, instance);
        assert_eq!(expected, "solid");
        assert_eq!(actual, "instance");

        let solid = builder
            .push(GeometryNode::Primitive(
                axiolid_primitive::Primitive::Sphere { radius: 1.0 },
            ))
            .unwrap();
        let solid_instance = builder
            .push(GeometryNode::Instance(Instance {
                source: solid,
                transform: axiolid_core::Transform3::IDENTITY,
            }))
            .unwrap();
        assert!(builder
            .push(GeometryNode::SolidOperation(
                crate::SolidOperation::Boolean {
                    left: solid_instance,
                    right: solid_instance,
                    operator: axiolid_core::BooleanOperator::Union,
                },
            ))
            .is_ok());
    }

    #[test]
    fn surface_curve_requires_a_three_dimensional_basis() {
        let mut builder = GeometryGraphBuilder::new();
        let curve_2d = builder
            .push(GeometryNode::Curve2(axiolid_curve::Curve2::Line(
                axiolid_curve::Line2 {
                    origin: axiolid_core::Vec2::ZERO,
                    direction: axiolid_core::Vec2::X,
                },
            )))
            .unwrap();
        let plane = builder
            .push(GeometryNode::Surface(axiolid_surface::Surface::Plane(
                axiolid_surface::Plane {
                    frame: axiolid_core::Frame3 {
                        origin: Vec3::ZERO,
                        x: Vec3::X,
                        y: Vec3::Y,
                        z: Vec3::Z,
                    },
                },
            )))
            .unwrap();
        let error = builder
            .push(GeometryNode::CurveRelation(
                crate::CurveRelation::SurfaceCurve {
                    curve_3d: curve_2d,
                    sides: crate::SurfaceSides::one(plane, curve_2d),
                    master: crate::MasterRepresentation::Curve3d,
                },
            ))
            .unwrap_err();
        assert!(matches!(
            error,
            GraphError::InvalidReferenceType {
                reference,
                expected: "curve3",
                actual: "curve2",
            } if reference == curve_2d
        ));

        let curve_3d = builder
            .push(GeometryNode::Curve3(axiolid_curve::Curve3::Line(
                axiolid_curve::Line3 {
                    origin: Vec3::ZERO,
                    direction: Vec3::X,
                },
            )))
            .unwrap();
        let trimmed_3d = builder
            .push(GeometryNode::CurveRelation(crate::CurveRelation::Trimmed {
                basis: curve_3d,
                start: Vec::new(),
                end: Vec::new(),
                sense_agreement: true,
                preference: crate::TrimmingPreference::Unspecified,
            }))
            .unwrap();
        assert!(builder
            .push(GeometryNode::CurveRelation(
                crate::CurveRelation::SurfaceCurve {
                    curve_3d: trimmed_3d,
                    sides: crate::SurfaceSides::one(plane, curve_2d),
                    master: crate::MasterRepresentation::Curve3d,
                },
            ))
            .is_ok());
    }

    #[test]
    fn parameter_curve_requires_a_two_dimensional_reference() {
        let mut builder = GeometryGraphBuilder::new();
        let surface = builder
            .push(GeometryNode::Surface(axiolid_surface::Surface::Plane(
                axiolid_surface::Plane {
                    frame: axiolid_core::Frame3 {
                        origin: Vec3::ZERO,
                        x: Vec3::X,
                        y: Vec3::Y,
                        z: Vec3::Z,
                    },
                },
            )))
            .unwrap();
        let curve_3d = builder
            .push(GeometryNode::Curve3(axiolid_curve::Curve3::Line(
                axiolid_curve::Line3 {
                    origin: Vec3::ZERO,
                    direction: Vec3::X,
                },
            )))
            .unwrap();
        let error = builder
            .push(GeometryNode::CurveRelation(
                crate::CurveRelation::ParameterCurve {
                    basis_surface: surface,
                    reference_curve: curve_3d,
                },
            ))
            .unwrap_err();
        assert!(matches!(
            error,
            GraphError::InvalidReferenceType {
                reference,
                expected: "curve2",
                actual: "curve3",
            } if reference == curve_3d
        ));
    }
}
