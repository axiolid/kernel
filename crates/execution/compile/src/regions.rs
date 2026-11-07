//! Which nodes are regions, decided on demand and remembered (#326).
//!
//! A region is a half-space, an instance of a region, or an intersection of
//! regions (the tool of a clip, see `exact::clip`). Both compilers ask the
//! question of the operands of a boolean, and an intersection of n
//! half-spaces is n nodes deep, so it must be neither walked again at each
//! link nor decided for the whole graph up front: a batch of independent
//! roots would pay the graph once per root, and one exact compilation per
//! boolean node the same.
//!
//! [`RegionMemo::has`] walks only what the question reaches, following
//! instance sources and the operands of an intersection, with an explicit
//! stack (a nest is as deep as the file makes it), and remembers every node
//! it decides. The memo lives as long as the compilation or the batch that
//! owns it, so over a chain of n clips the total work is O(n). An
//! intersection with one operand that is no region is decided at the first,
//! without reading the other.

use std::cell::RefCell;
use std::collections::HashMap;

use axiolid_core::BooleanOperator;
use axiolid_model::{GeometryGraph, GeometryNode, NodeId, SolidOperation};

/// What is known of the nodes of one graph.
#[derive(Debug, Default)]
pub(crate) struct RegionMemo {
    known: RefCell<HashMap<NodeId, bool>>,
}

#[cfg(test)]
thread_local! {
    /// Nodes decided by [`RegionMemo::has`] on this thread, for the tests
    /// that the work of a batch is linear in its size.
    pub(crate) static DECIDED: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

impl RegionMemo {
    /// Whether `id` is a region, `leaf` saying whether a node that is
    /// neither an instance nor an intersection is a half-space.
    pub(crate) fn has(
        &self,
        graph: &GeometryGraph,
        id: NodeId,
        leaf: impl Fn(NodeId) -> bool,
    ) -> bool {
        let mut known = self.known.borrow_mut();
        let mut stack = vec![id];
        while let Some(&top) = stack.last() {
            if known.contains_key(&top) {
                stack.pop();
                continue;
            }
            // `Some` once the node is decided, else the operand to decide
            // first.
            let decided: Result<bool, NodeId> = match graph.get(top) {
                Some(GeometryNode::Instance(instance)) => match known.get(&instance.source) {
                    Some(&is) => Ok(is),
                    None => Err(instance.source),
                },
                Some(GeometryNode::SolidOperation(SolidOperation::Boolean {
                    left,
                    right,
                    operator: BooleanOperator::Intersection,
                })) => match (known.get(left), known.get(right)) {
                    (Some(false), _) | (_, Some(false)) => Ok(false),
                    (Some(true), Some(true)) => Ok(true),
                    (None, _) => Err(*left),
                    (_, None) => Err(*right),
                },
                Some(_) => Ok(leaf(top)),
                None => Ok(false),
            };
            match decided {
                Ok(is) => {
                    #[cfg(test)]
                    DECIDED.with(|n| n.set(n.get() + 1));
                    known.insert(top, is);
                    stack.pop();
                }
                Err(operand) => stack.push(operand),
            }
        }
        known[&id]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_contracts::ExecutionOptions;
    use axiolid_core::{Plane3, Point3, Tolerance, Transform3, Vec3};
    use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
    use axiolid_mesh_compile_contract::MeshCompiler;
    use axiolid_model::{GeometryGraphBuilder, Instance};
    use axiolid_primitive::HalfSpace;
    use axiolid_profile::{Profile, RectangleProfile};

    fn half_space(b: &mut GeometryGraphBuilder, z: f64) -> NodeId {
        b.push(GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 {
                origin: Point3::new(0.0, 0.0, z),
                normal: Vec3::Z,
            },
            agreement: true,
        }))
        .unwrap()
    }

    fn decided() -> usize {
        DECIDED.with(std::cell::Cell::get)
    }

    /// A batch of N independent clips decides O(N) nodes, not N times the
    /// graph: the memo is the batch's, and a root asks only about its own
    /// operands.
    #[test]
    fn a_batch_of_roots_decides_a_number_of_nodes_linear_in_its_size() {
        let n = 200;
        let mut b = GeometryGraphBuilder::new();
        let profile = b
            .push(GeometryNode::Profile(Profile::Rectangle(
                RectangleProfile {
                    x: 1.0,
                    y: 1.0,
                    thickness: None,
                    outer_radius: None,
                    inner_radius: None,
                },
            )))
            .unwrap();
        let block = b
            .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
                profile,
                direction: Vec3::Z,
                depth: 1.0,
            }))
            .unwrap();
        let mut roots = Vec::new();
        for k in 0..n {
            let placed = b
                .push(GeometryNode::Instance(Instance {
                    source: block,
                    transform: Transform3::from_translation(Vec3::new(3.0 * k as f64, 0.0, 0.0)),
                }))
                .unwrap();
            let side = half_space(&mut b, 0.5);
            roots.push(
                b.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                    left: placed,
                    right: side,
                    operator: BooleanOperator::Difference,
                }))
                .unwrap(),
            );
        }
        let graph = b.finish(roots.clone()).unwrap();
        let compiler = crate::ReferenceMeshCompiler::new(BoolmeshBoolean::new());
        DECIDED.with(|d| d.set(0));
        let mut meshes = Vec::new();
        compiler
            .compile_mesh_batch_into(
                &graph,
                &roots,
                &ExecutionOptions::new(Tolerance::METRE),
                &mut meshes,
            )
            .unwrap();
        assert_eq!(meshes.len(), n);
        // Per root: its half-space, its placed block and its source, the
        // last shared. A walk of the graph per root would be 3 N^2.
        assert!(
            decided() <= 4 * n,
            "{} nodes decided for {n} roots",
            decided()
        );
    }

    /// Asking about each link of a nest of intersections, as a chain of
    /// clips does, decides each node once; and the walk is iterative.
    #[test]
    fn a_nest_is_decided_once_over_all_its_links() {
        let n = 5_000;
        let mut b = GeometryGraphBuilder::new();
        let mut links = Vec::new();
        let mut so_far = half_space(&mut b, 0.0);
        links.push(so_far);
        for k in 1..n {
            let side = half_space(&mut b, -(k as f64));
            so_far = b
                .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                    left: so_far,
                    right: side,
                    operator: BooleanOperator::Intersection,
                }))
                .unwrap();
            links.push(so_far);
        }
        let graph = b.finish(vec![so_far]).unwrap();
        let memo = RegionMemo::default();
        DECIDED.with(|d| d.set(0));
        for &link in &links {
            assert!(memo.has(&graph, link, |_| true));
        }
        assert!(decided() <= graph.len(), "{} decided", decided());
        // The deepest link first, on a fresh memo.
        assert!(RegionMemo::default().has(&graph, so_far, |_| true));
    }
}
