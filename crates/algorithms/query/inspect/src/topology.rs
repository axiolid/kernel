//! Topology of a triangle mesh, per connected component (#144).
//!
//! [`genus`](fn@crate::genus) answers one number for one closed orientable
//! surface and refuses everything else. This reports every component of a
//! two-manifold mesh, with or without boundary and whether or not it is
//! orientable: its counts, Euler characteristic, boundary loops,
//! orientability and the surface it is (genus `g` orientable, or `k`
//! crosscaps). For a closed orientable component it also gives a basis of
//! its first homology: `2g` closed edge loops, by the tree-cotree
//! construction (Eppstein, "Dynamic generators of topologically embedded
//! graphs", 2003).
//!
//! Everything is combinatorial: positions are never read, so the answers
//! are exact. A mesh that is not a two-manifold -- an edge on three or more
//! triangles, or a vertex whose triangles form more than one fan -- is
//! refused, since the classification of surfaces does not describe it.
//!
//! Not provided: homology generators of surfaces with boundary or of
//! non-orientable surfaces, and homotopy questions (whether a given closed
//! path is contractible, or two paths homotopic).

use std::collections::BTreeMap;

use axiolid_mesh::TriMesh;
use thiserror::Error;

/// Why a mesh's topology could not be classified.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum TopologyError {
    /// A triangle names a vertex the mesh does not have.
    #[error("triangle {triangle} names vertex {vertex}, but the mesh has {vertices}")]
    IndexOutOfRange {
        /// The offending triangle.
        triangle: usize,
        /// The index it names.
        vertex: u32,
        /// How many positions the mesh has.
        vertices: usize,
    },
    /// A triangle repeats a vertex, so it has no edges of its own.
    #[error("triangle {triangle} repeats a vertex")]
    DegenerateTriangle {
        /// The offending triangle.
        triangle: usize,
    },
    /// The mesh is not a two-manifold.
    #[error(
        "the mesh is not a two-manifold: {edges} edges on three or more \
         triangles, {vertices} vertices whose triangles form several fans"
    )]
    NonManifold {
        /// Edges used by three or more triangles.
        edges: usize,
        /// Vertices whose incident triangles form more than one fan.
        vertices: usize,
    },
}

/// Which closed surface a component is, with its boundary loops filled in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SurfaceKind {
    /// An orientable surface of genus `genus`: a sphere is 0, a torus 1.
    Orientable {
        /// Number of handles.
        genus: u32,
    },
    /// A non-orientable surface with `crosscaps` crosscaps: a projective
    /// plane or Möbius strip is 1, a Klein bottle 2.
    NonOrientable {
        /// Number of crosscaps.
        crosscaps: u32,
    },
}

/// A closed loop of mesh edges, as its vertices in order: consecutive
/// vertices, and the last and the first, are joined by an edge.
pub type EdgeLoop = Vec<u32>;

/// The topology of one connected component.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct ComponentTopology {
    /// The component's triangles, as indices into the mesh's triangles.
    pub triangles: Vec<usize>,
    /// Vertices the component's triangles use.
    pub vertices: usize,
    /// Distinct edges of the component's triangles.
    pub edges: usize,
    /// `vertices - edges + triangles`.
    pub euler_characteristic: i64,
    /// Closed loops of boundary edges (edges on one triangle).
    pub boundary_loops: usize,
    /// Whether the triangles can be wound consistently.
    pub orientable: bool,
    /// Whether they are: every interior edge is used once in each
    /// direction. Only an orientable component can be.
    pub consistently_oriented: bool,
    /// The surface, from `euler_characteristic = 2 - 2g - b` (orientable)
    /// or `2 - k - b` (non-orientable).
    pub surface: SurfaceKind,
    /// For a closed orientable component, `2g` closed edge loops forming a
    /// basis of its first homology; cutting along all of them leaves a
    /// disc. `None` for a component with boundary or a non-orientable one.
    pub homology_basis: Option<Vec<EdgeLoop>>,
}

/// The topology of a two-manifold triangle mesh.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct MeshTopology {
    /// Connected components (triangles joined through shared edges), in
    /// order of their first triangle.
    pub components: Vec<ComponentTopology>,
}

/// One use of an edge by a triangle.
#[derive(Debug, Clone, Copy)]
struct Use {
    triangle: usize,
    /// Whether the triangle runs the edge from its smaller vertex to its
    /// larger one.
    forward: bool,
}

/// Classify every connected component of `mesh`.
///
/// Positions not used by any triangle are ignored.
///
/// # Errors
///
/// A triangle naming a missing vertex or repeating one, and a mesh that is
/// not a two-manifold ([`TopologyError::NonManifold`]).
pub fn topology(mesh: &TriMesh) -> Result<MeshTopology, TopologyError> {
    let triangles: Vec<[u32; 3]> = mesh
        .indices
        .chunks_exact(3)
        .map(|t| [t[0], t[1], t[2]])
        .collect();
    let vertex_count = mesh.positions.len();
    for (triangle, t) in triangles.iter().enumerate() {
        if let Some(&vertex) = t.iter().find(|&&v| v as usize >= vertex_count) {
            return Err(TopologyError::IndexOutOfRange {
                triangle,
                vertex,
                vertices: vertex_count,
            });
        }
        if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] {
            return Err(TopologyError::DegenerateTriangle { triangle });
        }
    }

    let mut edges: BTreeMap<(u32, u32), Vec<Use>> = BTreeMap::new();
    for (triangle, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            edges.entry((a.min(b), a.max(b))).or_default().push(Use {
                triangle,
                forward: a < b,
            });
        }
    }
    let non_manifold_edges = edges.values().filter(|uses| uses.len() > 2).count();
    let non_manifold_vertices = count_split_vertices(&triangles, &edges, vertex_count);
    if non_manifold_edges > 0 || non_manifold_vertices > 0 {
        return Err(TopologyError::NonManifold {
            edges: non_manifold_edges,
            vertices: non_manifold_vertices,
        });
    }

    // Components through shared edges, and a relative orientation of each
    // triangle found on the way: `flip[t]` says whether `t` must be reversed
    // to agree with the first triangle of its component.
    let mut neighbours: Vec<Vec<(usize, bool)>> = vec![Vec::new(); triangles.len()];
    for uses in edges.values() {
        if let [a, b] = uses[..] {
            // Consistent neighbours run their shared edge in opposite
            // directions; the same direction means one must flip.
            let flips = a.forward == b.forward;
            neighbours[a.triangle].push((b.triangle, flips));
            neighbours[b.triangle].push((a.triangle, flips));
        }
    }
    let mut component_of = vec![usize::MAX; triangles.len()];
    let mut flip = vec![false; triangles.len()];
    let mut components = Vec::new();
    for seed in 0..triangles.len() {
        if component_of[seed] != usize::MAX {
            continue;
        }
        let id = components.len();
        component_of[seed] = id;
        let mut members = vec![seed];
        let mut orientable = true;
        let mut consistent = true;
        let mut next = 0;
        while next < members.len() {
            let t = members[next];
            next += 1;
            for &(n, flips) in &neighbours[t] {
                consistent &= !flips;
                let wanted = flip[t] ^ flips;
                if component_of[n] == usize::MAX {
                    component_of[n] = id;
                    flip[n] = wanted;
                    members.push(n);
                } else if flip[n] != wanted {
                    orientable = false;
                }
            }
        }
        members.sort_unstable();
        components.push((members, orientable, consistent));
    }

    let mut out = Vec::with_capacity(components.len());
    for (members, orientable, consistent) in components {
        let component = component_of[members[0]];
        let mut used: Vec<u32> = members.iter().flat_map(|&t| triangles[t]).collect();
        used.sort_unstable();
        used.dedup();
        let own_edges: Vec<(&(u32, u32), &Vec<Use>)> = edges
            .iter()
            .filter(|(_, uses)| component_of[uses[0].triangle] == component)
            .collect();
        let boundary: Vec<(u32, u32)> = own_edges
            .iter()
            .filter(|(_, uses)| uses.len() == 1)
            .map(|(&edge, _)| edge)
            .collect();
        let boundary_loops = count_loops(&boundary);
        let characteristic = used.len() as i64 - own_edges.len() as i64 + members.len() as i64;
        // chi = 2 - 2g - b or 2 - k - b; both deficits are non-negative on
        // a connected two-manifold.
        let deficit = 2 - characteristic - boundary_loops as i64;
        let surface = if orientable {
            SurfaceKind::Orientable {
                genus: u32::try_from(deficit / 2).unwrap_or(0),
            }
        } else {
            SurfaceKind::NonOrientable {
                crosscaps: u32::try_from(deficit).unwrap_or(0),
            }
        };
        let homology_basis =
            (orientable && boundary.is_empty()).then(|| tree_cotree(&members, &own_edges));
        out.push(ComponentTopology {
            triangles: members,
            vertices: used.len(),
            edges: own_edges.len(),
            euler_characteristic: characteristic,
            boundary_loops,
            orientable,
            consistently_oriented: consistent,
            surface,
            homology_basis,
        });
    }
    Ok(MeshTopology { components: out })
}

/// Vertices whose incident triangles form more than one fan (joined
/// through edges at the vertex), as at the tip where two cones touch.
fn count_split_vertices(
    triangles: &[[u32; 3]],
    edges: &BTreeMap<(u32, u32), Vec<Use>>,
    vertex_count: usize,
) -> usize {
    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); vertex_count];
    for (t, tri) in triangles.iter().enumerate() {
        for &v in tri {
            incident[v as usize].push(t);
        }
    }
    let mut split = 0;
    for (v, around) in incident.iter().enumerate() {
        if around.len() < 2 {
            continue;
        }
        let v = v as u32;
        let slot = |t: usize| around.iter().position(|&u| u == t);
        let mut parent: Vec<usize> = (0..around.len()).collect();
        for (i, &t) in around.iter().enumerate() {
            for &w in &triangles[t] {
                if w == v {
                    continue;
                }
                for other in &edges[&(v.min(w), v.max(w))] {
                    if let Some(j) = slot(other.triangle) {
                        union(&mut parent, i, j);
                    }
                }
            }
        }
        let fans = (0..around.len())
            .filter(|&i| find(&mut parent, i) == i)
            .count();
        if fans > 1 {
            split += 1;
        }
    }
    split
}

fn find(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

fn union(parent: &mut [usize], a: usize, b: usize) {
    let (a, b) = (find(parent, a), find(parent, b));
    if a != b {
        parent[a.max(b)] = a.min(b);
    }
}

/// Boundary loops: on a two-manifold every boundary vertex has exactly two
/// boundary edges, so the loops are the components of the boundary edges.
fn count_loops(boundary: &[(u32, u32)]) -> usize {
    let mut index: BTreeMap<u32, usize> = BTreeMap::new();
    for &(a, b) in boundary {
        let n = index.len();
        index.entry(a).or_insert(n);
        let n = index.len();
        index.entry(b).or_insert(n);
    }
    let mut parent: Vec<usize> = (0..index.len()).collect();
    for &(a, b) in boundary {
        union(&mut parent, index[&a], index[&b]);
    }
    (0..parent.len())
        .filter(|&i| find(&mut parent, i) == i)
        .count()
}

/// A homology basis of a closed orientable component by tree-cotree: a
/// spanning tree `T` of its vertices, a spanning tree of its triangles
/// across edges not in `T`, and one loop for each edge in neither -- the
/// edge closed through `T`. There are `E - (V - 1) - (F - 1) = 2g` of them.
fn tree_cotree(members: &[usize], own_edges: &[(&(u32, u32), &Vec<Use>)]) -> Vec<EdgeLoop> {
    // Primal spanning tree, breadth first from the smallest vertex, so the
    // loops come out short and the result is deterministic.
    let mut adjacent: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for (&(a, b), _) in own_edges {
        adjacent.entry(a).or_default().push(b);
        adjacent.entry(b).or_default().push(a);
    }
    let root = *adjacent.keys().next().expect("a component has vertices");
    let mut parent: BTreeMap<u32, u32> = BTreeMap::new();
    let mut depth: BTreeMap<u32, usize> = BTreeMap::new();
    parent.insert(root, root);
    depth.insert(root, 0);
    let mut queue = std::collections::VecDeque::from([root]);
    while let Some(v) = queue.pop_front() {
        for &w in &adjacent[&v] {
            if let std::collections::btree_map::Entry::Vacant(slot) = parent.entry(w) {
                slot.insert(v);
                depth.insert(w, depth[&v] + 1);
                queue.push_back(w);
            }
        }
    }
    // The root is its own parent, and no edge joins a vertex to itself.
    let in_tree = |a: u32, b: u32| parent[&a] == b || parent[&b] == a;

    // Dual spanning tree over the edges left: union-find on triangles.
    let slot: BTreeMap<usize, usize> = members.iter().enumerate().map(|(i, &t)| (t, i)).collect();
    let mut dual: Vec<usize> = (0..members.len()).collect();
    let mut leftover = Vec::new();
    for (&(a, b), uses) in own_edges {
        if in_tree(a, b) {
            continue;
        }
        let (s, t) = (slot[&uses[0].triangle], slot[&uses[1].triangle]);
        if find(&mut dual, s) == find(&mut dual, t) {
            leftover.push((a, b));
        } else {
            union(&mut dual, s, t);
        }
    }
    leftover
        .into_iter()
        .map(|(a, b)| {
            // Walk both ends up to their lowest common ancestor.
            let (mut up_a, mut up_b) = (vec![a], vec![b]);
            let (mut x, mut y) = (a, b);
            while depth[&x] > depth[&y] {
                x = parent[&x];
                up_a.push(x);
            }
            while depth[&y] > depth[&x] {
                y = parent[&y];
                up_b.push(y);
            }
            while x != y {
                x = parent[&x];
                y = parent[&y];
                up_a.push(x);
                up_b.push(y);
            }
            // a .. lca, then back down to b; the loop closes over (b, a).
            up_b.pop();
            up_a.extend(up_b.into_iter().rev());
            up_a
        })
        .collect()
}
