//! The B-rep tessellator's resolved support view.
//!
//! The tessellator reads a B-rep whose edges, uses and faces carry geometry
//! handles. Resolving them once, here, gives it a [`ResolvedBRep`] that
//! needs no graph lookups: an edge's 3D curve, a use's pcurve with the
//! parameter span the use runs over, and a face's surface are borrowed
//! values, so one tessellator serves any source of supports (#412, #391).
//!
//! [`from_graph`] resolves handles into a [`GeometryGraph`]. The copy adds
//! every entity in index order, so each id of the source B-rep names the
//! same entity in the view. A handle that does not resolve to a node of its
//! kind is refused by name here, up front, whether or not a meshed face
//! reaches it.
//!
//! [`from_exact`] resolves the typed handles of an [`ExactBRep`] into its own
//! catalogs, and a use's span is the interval the B-rep records for it
//! (`pcurve_interval`), which may run against the pcurve's own direction. A
//! pcurve use with no recorded interval, or a handle outside its catalog, is
//! refused (#391).
//!
//! A graph use's span is the parameter domain of its pcurve; sampling always runs
//! from `span.start` to `span.end`. `key` ranks a pcurve among the distinct
//! pcurve handles of the B-rep, so an edge shared by several uses keeps the
//! pcurve of the smallest handle ([`super::edges`]).

use super::support::Support;
use axiolid_brep::{Curve2Id, ExactBRep};
use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::Interval;
use axiolid_curve::{Curve2, Curve3};
use axiolid_model::{GeometryGraph, GeometryNode, NodeId};
use axiolid_topology::{BRep, Edge, EdgeUse, Face, Loop};
use std::collections::{BTreeSet, HashMap};

/// A B-rep whose supports are borrowed values.
pub(crate) type ResolvedBRep<'a> = BRep<&'a Curve3, Pcurve<'a>, Support<'a>>;

/// A pcurve with the parameter span one edge use runs over.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Pcurve<'a> {
    pub curve: &'a Curve2,
    /// The parameters sampled, from the use's start to its end.
    pub span: Interval,
    /// Rank among the B-rep's distinct pcurve handles, in handle order.
    pub key: u32,
}

/// Resolve `brep`'s handles in `graph`.
pub(crate) fn from_graph<'a>(
    brep: &'a BRep<NodeId>,
    graph: &'a GeometryGraph,
) -> GeomResult<ResolvedBRep<'a>> {
    let handles: BTreeSet<NodeId> = brep
        .loops()
        .iter()
        .flat_map(|wire| wire.edges.iter())
        .filter_map(|use_| use_.pcurve)
        .collect();
    let ranks: HashMap<NodeId, u32> = handles
        .into_iter()
        .enumerate()
        .map(|(rank, id)| (id, rank as u32))
        .collect();
    let mut view = ResolvedBRep::default();
    for vertex in brep.vertices() {
        view.add_vertex(*vertex);
    }
    for edge in brep.edges() {
        let curve = edge.curve.map(|id| edge_curve(graph, id)).transpose()?;
        view.add_edge(Edge {
            start: edge.start,
            end: edge.end,
            curve,
        });
    }
    for wire in brep.loops() {
        let edges = wire
            .edges
            .iter()
            .map(|use_| {
                let pcurve = use_
                    .pcurve
                    .map(|id| {
                        let curve = use_pcurve(graph, id)?;
                        Ok(Pcurve {
                            curve,
                            span: axiolid_reference::curve::domain2(curve),
                            key: ranks[&id],
                        })
                    })
                    .transpose()?;
                Ok(EdgeUse {
                    edge: use_.edge,
                    orientation: use_.orientation,
                    pcurve,
                })
            })
            .collect::<GeomResult<Vec<_>>>()?;
        view.add_loop(Loop { edges });
    }
    for face in brep.faces() {
        view.add_face(Face {
            surface: face_surface(graph, face)?,
            bounds: face.bounds.clone(),
            orientation: face.orientation,
        });
    }
    for shell in brep.shells() {
        view.add_shell(shell.clone());
    }
    for solid in brep.solids() {
        view.add_solid(solid.clone());
    }
    Ok(view)
}

/// Resolve `exact`'s handles in its own catalogs.
pub(crate) fn from_exact(exact: &ExactBRep) -> GeomResult<ResolvedBRep<'_>> {
    let brep = exact.topology();
    let handles: BTreeSet<Curve2Id> = brep
        .loops()
        .iter()
        .flat_map(|wire| wire.edges.iter())
        .filter_map(|use_| use_.pcurve)
        .collect();
    let ranks: HashMap<Curve2Id, u32> = handles
        .into_iter()
        .enumerate()
        .map(|(rank, id)| (id, rank as u32))
        .collect();
    let missing =
        |what: &str| GeomError::InvalidInput(format!("exact B-rep {what} is not in its catalog"));
    let mut view = ResolvedBRep::default();
    for vertex in brep.vertices() {
        view.add_vertex(*vertex);
    }
    for edge in brep.edges() {
        let curve = edge
            .curve
            .map(|id| {
                exact
                    .curves3()
                    .get(id.index())
                    .ok_or_else(|| missing("edge curve"))
            })
            .transpose()?;
        view.add_edge(Edge {
            start: edge.start,
            end: edge.end,
            curve,
        });
    }
    for (index, wire) in brep.loops().iter().enumerate() {
        let loop_id = brep.loop_id_at(index).ok_or_else(|| missing("loop"))?;
        let edges = wire
            .edges
            .iter()
            .enumerate()
            .map(|(position, use_)| {
                let pcurve = use_
                    .pcurve
                    .map(|id| {
                        let curve = exact
                            .curves2()
                            .get(id.index())
                            .ok_or_else(|| missing("pcurve"))?;
                        let span = exact
                            .pcurve_interval(loop_id, position)
                            .ok_or_else(|| missing("pcurve interval"))?;
                        Ok(Pcurve {
                            curve,
                            span,
                            key: ranks[&id],
                        })
                    })
                    .transpose()?;
                Ok(EdgeUse {
                    edge: use_.edge,
                    orientation: use_.orientation,
                    pcurve,
                })
            })
            .collect::<GeomResult<Vec<_>>>()?;
        view.add_loop(Loop { edges });
    }
    for face in brep.faces() {
        let surface = face
            .surface
            .map(|id| {
                exact
                    .surfaces()
                    .get(id.index())
                    .map(Support::Surface)
                    .ok_or_else(|| missing("face surface"))
            })
            .transpose()?;
        view.add_face(Face {
            surface,
            bounds: face.bounds.clone(),
            orientation: face.orientation,
        });
    }
    for shell in brep.shells() {
        view.add_shell(shell.clone());
    }
    for solid in brep.solids() {
        view.add_solid(solid.clone());
    }
    Ok(view)
}

fn edge_curve(graph: &GeometryGraph, id: NodeId) -> GeomResult<&Curve3> {
    match graph.get(id) {
        Some(GeometryNode::Curve3(curve)) => Ok(curve),
        _ => Err(GeomError::InvalidInput(
            "edge curve must reference a Curve3 node".to_owned(),
        )),
    }
}

fn use_pcurve(graph: &GeometryGraph, id: NodeId) -> GeomResult<&Curve2> {
    match graph.get(id) {
        Some(GeometryNode::Curve2(curve)) => Ok(curve),
        _ => Err(GeomError::InvalidInput(
            "edge pcurve must reference a Curve2 node".to_owned(),
        )),
    }
}

/// Resolve a face's support surface from the graph, if it declares one.
///
/// A face may legitimately omit its surface: a planar polygon's boundary
/// already determines its plane. A declared handle that does not resolve to
/// a Surface node is a graph error, not an absent surface; one that resolves
/// to a surface relation is refused by name.
fn face_surface<'g>(
    graph: &'g GeometryGraph,
    face: &Face<NodeId>,
) -> GeomResult<Option<Support<'g>>> {
    let Some(id) = face.surface else {
        return Ok(None);
    };
    match graph.get(id) {
        Some(GeometryNode::Surface(s)) => Ok(Some(Support::Surface(s))),
        // A relation (a ruled or swept surface) is a model node with no
        // `Surface` value to tessellate a face over; a refusal by name, not
        // a malformed graph (#402).
        Some(GeometryNode::SurfaceRelation(_)) => Err(GeomError::UnsupportedInput {
            backend: crate::BACKEND_ID,
            operation: Operation::Tessellation,
            input: "a B-rep face on a surface relation",
        }),
        Some(_) => Err(GeomError::InvalidInput(format!(
            "face support {id:?} is not a Surface node"
        ))),
        None => Err(GeomError::InvalidInput(format!(
            "face support {id:?} does not belong to this graph"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::{Frame3, Point2, Point3, Vec3};
    use axiolid_curve::{Polyline2, Polyline3};
    use axiolid_model::GeometryGraphBuilder;
    use axiolid_surface::{Plane, Surface};
    use axiolid_topology::{FaceBound, Orientation, Shell, Vertex};

    fn pcurve(length: f64) -> Polyline2 {
        Polyline2 {
            points: vec![Point2::new(0.0, 0.0), Point2::new(length, 0.0)],
            closed: false,
        }
    }

    fn line() -> Polyline3 {
        Polyline3 {
            points: vec![Point3::ZERO, Point3::new(1.0, 0.0, 0.0)],
            closed: false,
        }
    }

    fn plane() -> Plane {
        Plane {
            frame: Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
        }
    }

    /// A shelled face whose outer loop uses the later-pushed pcurve and
    /// whose inner loop the earlier one, and a loose face on no shell whose
    /// loop uses the handle `stray` picks.
    struct Fixture {
        brep: BRep<NodeId>,
        graph: GeometryGraph,
        pcurves: [NodeId; 2],
        curve: NodeId,
        surface: NodeId,
    }

    fn fixture(stray: impl FnOnce(&[NodeId; 2], NodeId) -> NodeId) -> Fixture {
        let mut builder = GeometryGraphBuilder::default();
        let first = builder.push_value(Curve2::Polyline(pcurve(1.0))).unwrap();
        let second = builder.push_value(Curve2::Polyline(pcurve(2.5))).unwrap();
        let curve = builder.push_value(Curve3::Polyline(line())).unwrap();
        let surface = builder.push_value(Surface::Plane(plane())).unwrap();
        let graph = builder.finish(vec![surface]).unwrap();
        let pcurves = [first, second];
        let mut brep = BRep::<NodeId>::default();
        let at = |x| Vertex {
            position: Point3::new(x, 0.0, 0.0),
        };
        let (a, b) = (brep.add_vertex(at(0.0)), brep.add_vertex(at(1.0)));
        let edge = brep.add_edge(Edge {
            start: a,
            end: b,
            curve: Some(curve),
        });
        let bare = brep.add_edge(Edge {
            start: b,
            end: a,
            curve: None,
        });
        let used = |pcurve| EdgeUse {
            edge,
            orientation: Orientation::Reversed,
            pcurve,
        };
        let outer = brep.add_loop(Loop {
            edges: vec![
                used(Some(second)),
                EdgeUse {
                    edge: bare,
                    orientation: Orientation::Forward,
                    pcurve: None,
                },
            ],
        });
        let inner = brep.add_loop(Loop {
            edges: vec![used(Some(first))],
        });
        let bound = |loop_id, outer| FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer,
        };
        let face = brep.add_face(Face {
            surface: Some(surface),
            bounds: vec![bound(outer, true), bound(inner, false)],
            orientation: Orientation::Reversed,
        });
        brep.add_shell(Shell {
            faces: vec![(face, Orientation::Reversed)],
            closed: false,
        });
        let lone = brep.add_loop(Loop {
            edges: vec![used(Some(stray(&pcurves, curve)))],
        });
        brep.add_face(Face {
            surface: None,
            bounds: vec![bound(lone, true)],
            orientation: Orientation::Forward,
        });
        Fixture {
            brep,
            graph,
            pcurves,
            curve,
            surface,
        }
    }

    #[test]
    fn the_view_keeps_every_id_and_order_of_the_source() {
        let f = fixture(|pcurves, _| pcurves[0]);
        let view = from_graph(&f.brep, &f.graph).unwrap();
        assert_eq!(view.vertices(), f.brep.vertices());
        assert_eq!(view.shells(), f.brep.shells());
        assert_eq!(view.solids(), f.brep.solids());
        assert_eq!(view.edges().len(), f.brep.edges().len());
        for (seen, source) in view.edges().iter().zip(f.brep.edges()) {
            assert_eq!((seen.start, seen.end), (source.start, source.end));
            assert_eq!(seen.curve.is_some(), source.curve.is_some());
        }
        assert_eq!(view.loops().len(), f.brep.loops().len());
        for (seen, source) in view.loops().iter().zip(f.brep.loops()) {
            assert_eq!(seen.edges.len(), source.edges.len());
            for (a, b) in seen.edges.iter().zip(&source.edges) {
                assert_eq!(
                    (a.edge, a.orientation, a.pcurve.is_some()),
                    (b.edge, b.orientation, b.pcurve.is_some())
                );
            }
        }
        assert_eq!(view.faces().len(), f.brep.faces().len());
        for (seen, source) in view.faces().iter().zip(f.brep.faces()) {
            assert_eq!(seen.bounds, source.bounds);
            assert_eq!(seen.orientation, source.orientation);
            assert_eq!(seen.surface.is_some(), source.surface.is_some());
        }
    }

    #[test]
    fn a_pcurve_carries_its_domain_as_the_span() {
        let f = fixture(|pcurves, _| pcurves[0]);
        let view = from_graph(&f.brep, &f.graph).unwrap();
        let second = view.loops()[0].edges[0].pcurve.unwrap();
        let Some(GeometryNode::Curve2(curve)) = f.graph.get(f.pcurves[1]) else {
            panic!("pcurve node");
        };
        assert_eq!(second.span, axiolid_reference::curve::domain2(curve));
        assert!(std::ptr::eq(second.curve, curve));
    }

    #[test]
    fn ranks_follow_the_handle_order_not_the_use_order() {
        let f = fixture(|pcurves, _| pcurves[0]);
        let view = from_graph(&f.brep, &f.graph).unwrap();
        // The outer loop uses the later handle, the inner loop and the
        // loose face the earlier one: one handle, one rank.
        let key = |index: usize| view.loops()[index].edges[0].pcurve.unwrap().key;
        assert_eq!((key(0), key(1), key(2)), (1, 0, 0));
    }

    #[test]
    fn a_mistyped_pcurve_handle_no_meshed_face_reaches_is_refused_by_name() {
        // The loose face is on no shell, so no meshed face reads its loop.
        let f = fixture(|_, curve| curve);
        let error = from_graph(&f.brep, &f.graph).expect_err("refused");
        assert!(
            matches!(&error, GeomError::InvalidInput(text)
                if text == "edge pcurve must reference a Curve2 node"),
            "{error:?}"
        );
    }

    #[test]
    fn a_mistyped_edge_curve_or_face_support_is_refused_by_name() {
        let f = fixture(|pcurves, _| pcurves[0]);
        let mut wrong_curve = f.brep.clone();
        let (start, end) = (wrong_curve.edges()[0].start, wrong_curve.edges()[0].end);
        wrong_curve.add_edge(Edge {
            start,
            end,
            curve: Some(f.surface),
        });
        let error = from_graph(&wrong_curve, &f.graph).expect_err("refused");
        assert!(
            matches!(&error, GeomError::InvalidInput(text)
                if text == "edge curve must reference a Curve3 node"),
            "{error:?}"
        );
        let mut wrong_face = f.brep.clone();
        wrong_face.add_face(Face {
            surface: Some(f.curve),
            bounds: Vec::new(),
            orientation: Orientation::Forward,
        });
        let error = from_graph(&wrong_face, &f.graph).expect_err("refused");
        assert!(
            matches!(&error, GeomError::InvalidInput(text) if text.ends_with("is not a Surface node")),
            "{error:?}"
        );
    }
}
