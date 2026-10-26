//! A closed triangle mesh as an exact polyhedral B-rep (#266).
//!
//! A boolean whose operand is an authored mesh, or a B-rep whose faces are
//! all planes or warped faces without a surface, is measured against the
//! exact boolean of that operand's triangle mesh: the operand flattened
//! within its slab width (`crate::planar`), the perturbed operand of
//! ADR 0080. This module turns that mesh into the exact B-rep the general
//! boolean takes. It is only reached from the measuring path of
//! `crate::deviation::of_boolean`; the public exact compiler keeps refusing
//! meshes by name.
//!
//! Every triangle is one planar face: its plane through its first corner,
//! `x` along its first edge, `z` along its normal, so every corner lies in
//! it up to the rounding of that frame, and an axis-aligned triangle's
//! frame is exact. Every edge is a line between two mesh vertices, shared
//! by the two triangles that use it, with a line pcurve in each face's
//! parameters. Nothing is merged or moved: the B-rep is the mesh, so the
//! exact boolean of it is the exact boolean of the operand as meshed.
//!
//! Refused by name: a mesh that is not a closed, consistently wound
//! two-manifold by its corner indices, a triangle without area, and a
//! shell wound inside out (non-positive enclosed volume).

use std::collections::HashMap;

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame3, Interval, Point2, Point3, Scalar};
use axiolid_curve::{Curve2, Curve3, Line2, Line3};
use axiolid_mesh::TriMesh;
use axiolid_surface::{Plane, Surface};
use axiolid_topology::{
    audit_brep, Edge, EdgeId, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};

use super::unsupported;

/// The refusal of a mesh operand that is not a closed, consistently wound
/// two-manifold by its corner indices.
const NOT_CLOSED: &str =
    "exact boolean over a mesh operand that is not a closed, consistently wound two-manifold";

/// The exact polyhedral B-rep of a closed triangle mesh, its triangles its
/// faces.
pub(super) fn from_triangles(mesh: &TriMesh) -> GeomResult<ExactBRep> {
    // A triangle repeating a corner index covers no edge pair; the mesh
    // compiler's closure check skips it too.
    let triangles: Vec<[u32; 3]> = mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| [t[0], t[1], t[2]])
        .filter(|[a, b, c]| a != b && b != c && c != a)
        .collect();
    if triangles.len() < 4 {
        return Err(unsupported(NOT_CLOSED));
    }
    let position = |index: u32| -> GeomResult<Point3> {
        mesh.positions
            .get(index as usize)
            .copied()
            .filter(|p| p.is_finite())
            .ok_or_else(|| unsupported(NOT_CLOSED))
    };

    // Each undirected edge once each way: a closed, consistently wound
    // two-manifold, the mesh compiler's own closure test.
    let mut directed: HashMap<(u32, u32), usize> = HashMap::new();
    for t in &triangles {
        for k in 0..3 {
            *directed.entry((t[k], t[(k + 1) % 3])).or_default() += 1;
        }
    }
    let manifold = directed
        .iter()
        .all(|(&(a, b), &count)| count == 1 && directed.get(&(b, a)) == Some(&1));
    if !manifold {
        return Err(unsupported(NOT_CLOSED));
    }

    // Outward: six times the enclosed volume is positive.
    let mut volume: Scalar = 0.0;
    for t in &triangles {
        let [a, b, c] = [position(t[0])?, position(t[1])?, position(t[2])?];
        volume += a.dot(b.cross(c));
    }
    if volume.is_nan() || volume <= 0.0 {
        return Err(unsupported(
            "exact boolean over a mesh operand wound inside out",
        ));
    }

    let mut builder = ExactBRepBuilder::default();
    let mut vertices = HashMap::new();
    let mut edges: HashMap<(u32, u32), EdgeId> = HashMap::new();
    let mut faces = Vec::with_capacity(triangles.len());
    for t in &triangles {
        let corners = [position(t[0])?, position(t[1])?, position(t[2])?];
        let frame = triangle_frame(corners)
            .ok_or_else(|| unsupported("exact boolean over a mesh triangle without area"))?;
        let surface = builder.add_surface(Surface::Plane(Plane { frame }));
        let uv = |p: Point3| {
            let offset = p - frame.origin;
            Point2::new(offset.dot(frame.x), offset.dot(frame.y))
        };
        let mut uses = Vec::with_capacity(3);
        let mut intervals = Vec::with_capacity(3);
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            // The edge runs from its lower corner index to its higher; the
            // triangle walking it the other way uses it reversed.
            let (start, end) = (a.min(b), a.max(b));
            let edge = match edges.get(&(start, end)) {
                Some(&edge) => edge,
                None => {
                    let mut vertex = |index: u32| -> GeomResult<_> {
                        let at = position(index)?;
                        Ok(*vertices.entry(index).or_insert_with(|| {
                            builder.topology_mut().add_vertex(Vertex { position: at })
                        }))
                    };
                    let (from, to) = (vertex(start)?, vertex(end)?);
                    let (p, q) = (position(start)?, position(end)?);
                    let curve = builder.add_curve3(Curve3::Line(Line3 {
                        origin: p,
                        direction: q - p,
                    }));
                    let edge = builder.topology_mut().add_edge(Edge {
                        start: from,
                        end: to,
                        curve: Some(curve),
                    });
                    builder.set_edge_interval(edge, Interval::UNIT);
                    edges.insert((start, end), edge);
                    edge
                }
            };
            let (p, q) = (uv(position(start)?), uv(position(end)?));
            let pcurve = builder.add_curve2(Curve2::Line(Line2 {
                origin: p,
                direction: q - p,
            }));
            let forward = a == start;
            uses.push(EdgeUse {
                edge,
                orientation: if forward {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                },
                pcurve: Some(pcurve),
            });
            intervals.push(if forward {
                Interval::UNIT
            } else {
                Interval::new(1.0, 0.0)
            });
        }
        let loop_id = builder.topology_mut().add_loop(Loop { edges: uses });
        for (index, interval) in intervals.into_iter().enumerate() {
            builder.set_pcurve_interval(loop_id, index, interval);
        }
        faces.push(builder.topology_mut().add_face(Face {
            surface: Some(surface),
            bounds: vec![FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation: Orientation::Forward,
        }));
    }
    let shell = builder.topology_mut().add_shell(Shell {
        faces: faces
            .into_iter()
            .map(|face| (face, Orientation::Forward))
            .collect(),
        closed: true,
    });
    builder.topology_mut().add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    let exact = builder.finish().map_err(|error| {
        GeomError::Degenerate(format!(
            "a mesh operand did not assemble as an exact polyhedron: {error}"
        ))
    })?;
    if !audit_brep(exact.topology()).is_closed_manifold() {
        return Err(unsupported(NOT_CLOSED));
    }
    Ok(exact)
}

/// The plane of a triangle: origin at its first corner, `x` along its
/// first edge, `z` its normal, or `None` when it has no area.
fn triangle_frame([a, b, c]: [Point3; 3]) -> Option<Frame3> {
    let x = (b - a).normalize_or_zero();
    let z = (b - a).cross(c - a).normalize_or_zero();
    let y = z.cross(x);
    (x != Point3::ZERO && z != Point3::ZERO && y.is_finite()).then_some(Frame3 {
        origin: a,
        x,
        y,
        z,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The unit cube as twelve outward triangles.
    fn cube() -> TriMesh {
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 1.0),
        ];
        let quads = [
            [0, 3, 2, 1],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
        ];
        let indices = quads
            .iter()
            .flat_map(|q| [q[0], q[1], q[2], q[0], q[2], q[3]])
            .collect();
        TriMesh::new(positions, indices)
    }

    #[test]
    fn a_closed_cube_becomes_twelve_planar_faces() {
        let exact = from_triangles(&cube()).unwrap();
        assert_eq!(exact.topology().faces().len(), 12);
        assert_eq!(exact.topology().edges().len(), 18);
        assert_eq!(exact.topology().vertices().len(), 8);
        assert!(audit_brep(exact.topology()).is_closed_manifold());
    }

    #[test]
    fn an_open_or_inverted_mesh_is_refused_by_name() {
        let mut open = cube();
        open.indices.truncate(33);
        let named = |error: GeomError, needle: &str| matches!(error, GeomError::UnsupportedInput { input, .. } if input.contains(needle));
        assert!(named(from_triangles(&open).unwrap_err(), "closed"));
        let mut inverted = cube();
        for t in inverted.indices.as_chunks_mut::<3>().0 {
            t.swap(1, 2);
        }
        assert!(named(from_triangles(&inverted).unwrap_err(), "inside out"));
    }
}
