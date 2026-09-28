//! Pinches in a boolean's result (#194).
//!
//! Where two parts of a solid's boundary touch -- a void tangent to its
//! host's face, two bodies meeting along an edge -- the solid itself is not
//! a two-manifold there. A mesh boolean represents that by keeping two
//! copies of the vertices along the contact: closed and manifold by index,
//! but two surfaces at the same positions. A consumer that validates by
//! position sees an edge with four faces, or a vertex with two separate
//! fans, and rejects the body. Such a result is refused with the contact
//! named rather than returned.

use std::collections::HashMap;

use axiolid_core::Point3;
use axiolid_mesh::TriMesh;

/// A contact where the result touches itself, welded by exact position.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Pinch {
    /// Along a line: an edge shared by more than two faces.
    Edge { from: Point3, to: Point3 },
    /// At a point: a vertex whose faces form separate fans.
    Vertex { at: Point3 },
}

/// The first pinch in `mesh`, welding coincident positions; `None` when
/// every welded edge has at most two faces and every welded vertex one fan.
pub(crate) fn find(mesh: &TriMesh) -> Option<Pinch> {
    let key = |p: Point3| (p.x.to_bits(), p.y.to_bits(), p.z.to_bits());
    let mut ids: HashMap<(u64, u64, u64), u32> = HashMap::new();
    let mut at: Vec<Point3> = Vec::new();
    let welded: Vec<u32> = mesh
        .positions
        .iter()
        .map(|p| {
            *ids.entry(key(*p)).or_insert_with(|| {
                at.push(*p);
                (at.len() - 1) as u32
            })
        })
        .collect();
    let triangles: Vec<[u32; 3]> = mesh
        .indices
        .chunks_exact(3)
        .map(|t| [0, 1, 2].map(|k| welded[t[k] as usize]))
        .filter(|t| t[0] != t[1] && t[1] != t[2] && t[2] != t[0])
        .collect();
    let mut edges: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    for (index, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            edges.entry((a.min(b), a.max(b))).or_default().push(index);
        }
    }
    let mut sorted: Vec<(&(u32, u32), &Vec<usize>)> = edges.iter().collect();
    sorted.sort_unstable_by_key(|(edge, _)| **edge);
    if let Some((&(a, b), _)) = sorted.iter().find(|(_, faces)| faces.len() > 2) {
        return Some(Pinch::Edge {
            from: at[a as usize],
            to: at[b as usize],
        });
    }
    // Each vertex's faces, joined across the edges through the vertex.
    let mut around: Vec<Vec<usize>> = vec![Vec::new(); at.len()];
    for (index, t) in triangles.iter().enumerate() {
        for &v in t {
            around[v as usize].push(index);
        }
    }
    for (v, faces) in around.iter().enumerate() {
        if faces.len() < 2 {
            continue;
        }
        let v = v as u32;
        let mut parent: Vec<usize> = (0..faces.len()).collect();
        fn root(parent: &mut [usize], mut i: usize) -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        }
        let local: HashMap<usize, usize> = faces.iter().enumerate().map(|(i, f)| (*f, i)).collect();
        for &f in faces {
            for &w in &triangles[f] {
                if w == v {
                    continue;
                }
                for &g in &edges[&(v.min(w), v.max(w))] {
                    let (i, j) = (root(&mut parent, local[&f]), root(&mut parent, local[&g]));
                    parent[i] = j;
                }
            }
        }
        let first = root(&mut parent, 0);
        if (1..faces.len()).any(|i| root(&mut parent, i) != first) {
            return Some(Pinch::Vertex { at: at[v as usize] });
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tetrahedron's outward triangles over four corners.
    fn tetra(corners: [Point3; 4]) -> (Vec<Point3>, Vec<u32>) {
        (corners.to_vec(), vec![0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3])
    }

    /// Two meshes as one, the second's positions duplicated as a boolean
    /// keeps them: separate indices, same points.
    fn join(a: (Vec<Point3>, Vec<u32>), b: (Vec<Point3>, Vec<u32>)) -> TriMesh {
        let offset = a.0.len() as u32;
        let mut positions = a.0;
        positions.extend(b.0);
        let mut indices = a.1;
        indices.extend(b.1.into_iter().map(|i| i + offset));
        TriMesh::new(positions, indices)
    }

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn a_single_solid_has_no_pinch() {
        let (positions, indices) = tetra([
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, 0.0, 1.0),
        ]);
        assert_eq!(find(&TriMesh::new(positions, indices)), None);
    }

    #[test]
    fn two_solids_meeting_at_a_corner_pinch_there() {
        let a = tetra([
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, 0.0, 1.0),
        ]);
        let b = tetra([
            p(0.0, 0.0, 0.0),
            p(-1.0, 0.0, 0.0),
            p(0.0, -1.0, 0.0),
            p(0.0, 0.0, -1.0),
        ]);
        assert_eq!(
            find(&join(a, b)),
            Some(Pinch::Vertex {
                at: p(0.0, 0.0, 0.0)
            })
        );
    }

    #[test]
    fn two_solids_meeting_along_an_edge_pinch_there() {
        let a = tetra([
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, 0.0, 1.0),
        ]);
        let b = tetra([
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, -1.0, 0.0),
            p(0.0, 0.0, -1.0),
        ]);
        match find(&join(a, b)) {
            Some(Pinch::Edge { from, to }) => {
                let mut ends = [from, to];
                ends.sort_by(|u, v| u.x.total_cmp(&v.x));
                assert_eq!(ends, [p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0)]);
            }
            other => panic!("{other:?}"),
        }
    }
}
