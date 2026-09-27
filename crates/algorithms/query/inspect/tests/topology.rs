//! Mesh topology per component, and homology bases (#144).
//!
//! Oracles are the classification of surfaces: every fixture is built so
//! its genus, crosscaps and boundary loops are known by construction.

use std::collections::{BTreeMap, BTreeSet};

use axiolid_core::Point3;
use axiolid_inspect::{topology, ComponentTopology, SurfaceKind, TopologyError};
use axiolid_mesh::TriMesh;

/// Welds integer points into indexed positions.
#[derive(Default)]
struct Builder {
    index: BTreeMap<[i64; 3], u32>,
    positions: Vec<Point3>,
    indices: Vec<u32>,
}

impl Builder {
    fn vertex(&mut self, p: [i64; 3]) -> u32 {
        let next = self.positions.len() as u32;
        *self.index.entry(p).or_insert_with(|| {
            self.positions
                .push(Point3::new(p[0] as f64, p[1] as f64, p[2] as f64));
            next
        })
    }

    /// A quad as two triangles, wound `a b c d`.
    fn quad(&mut self, corners: [[i64; 3]; 4]) {
        let [a, b, c, d] = corners.map(|p| self.vertex(p));
        self.indices.extend_from_slice(&[a, b, c, a, c, d]);
    }

    fn finish(self) -> TriMesh {
        TriMesh::new(self.positions, self.indices)
    }
}

/// The outward boundary of a set of unit voxels.
fn voxels(cells: &[[i64; 3]]) -> TriMesh {
    let set: BTreeSet<[i64; 3]> = cells.iter().copied().collect();
    let mut b = Builder::default();
    for &cell in cells {
        for axis in 0..3 {
            let (u, w) = ((axis + 1) % 3, (axis + 2) % 3);
            for high in [false, true] {
                let mut neighbour = cell;
                neighbour[axis] += if high { 1 } else { -1 };
                if set.contains(&neighbour) {
                    continue;
                }
                let mut base = cell;
                if high {
                    base[axis] += 1;
                }
                let at = |du: i64, dw: i64| {
                    let mut p = base;
                    p[u] += du;
                    p[w] += dw;
                    p
                };
                // e_u x e_w = e_axis, so this order faces +axis.
                let mut corners = [at(0, 0), at(1, 0), at(1, 1), at(0, 1)];
                if !high {
                    corners.reverse();
                }
                b.quad(corners);
            }
        }
    }
    b.finish()
}

/// A flat mesh over unit squares of the plane z = 0.
fn pixels(cells: &[[i64; 2]]) -> TriMesh {
    let mut b = Builder::default();
    for &[x, y] in cells {
        b.quad([[x, y, 0], [x + 1, y, 0], [x + 1, y + 1, 0], [x, y + 1, 0]]);
    }
    b.finish()
}

/// A grid of `n x m` quads with the sides glued: `twist_u` glues the right
/// edge to the left one reversed (a Möbius band), `close_v` glues the top
/// to the bottom. Positions are placeholders; only the gluing matters.
fn glued_grid(n: u32, m: u32, twist_u: bool, close_v: bool) -> TriMesh {
    let rows = if close_v { m } else { m + 1 };
    let vertex = |i: u32, j: u32| -> u32 {
        let (mut i, mut j) = (i, j);
        if close_v {
            j %= m;
        }
        if i == n {
            i = 0;
            if twist_u {
                j = if close_v { (m - j) % m } else { m - j };
            }
        }
        j * n + i
    };
    let positions = (0..n * rows)
        .map(|k| Point3::new(f64::from(k % n), f64::from(k / n), 0.0))
        .collect();
    let mut indices = Vec::new();
    for i in 0..n {
        for j in 0..m {
            let (a, b, c, d) = (
                vertex(i, j),
                vertex(i + 1, j),
                vertex(i + 1, j + 1),
                vertex(i, j + 1),
            );
            indices.extend_from_slice(&[a, b, c, a, c, d]);
        }
    }
    TriMesh::new(positions, indices)
}

fn only(mesh: &TriMesh) -> ComponentTopology {
    let report = topology(mesh).expect("a two-manifold");
    assert_eq!(report.components.len(), 1, "one component");
    report.components.into_iter().next().unwrap()
}

/// A ring of voxels around `holes` missing cells of a `w x 3` slab.
fn slab_with_holes(w: i64, holes: &[i64]) -> TriMesh {
    let mut cells = Vec::new();
    for x in 0..w {
        for y in 0..3 {
            if !(y == 1 && holes.contains(&x)) {
                cells.push([x, y, 0]);
            }
        }
    }
    voxels(&cells)
}

/// Checks a claimed homology basis of a closed orientable component:
/// `2g` simple closed edge loops, independent as edge chains mod 2, whose
/// union does not separate the surface. A sum of loops bounding a set of
/// triangles `S` would cut `S` off from the rest, so together these make
/// the loops independent in homology, and `2g` of them a basis.
fn assert_basis(mesh: &TriMesh, component: &ComponentTopology, genus: u32) {
    let loops = component.homology_basis.as_ref().expect("a basis");
    assert_eq!(loops.len(), 2 * genus as usize);
    let mut edges: BTreeMap<(u32, u32), Vec<usize>> = BTreeMap::new();
    for (t, tri) in mesh.indices.chunks_exact(3).enumerate() {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            edges.entry((a.min(b), a.max(b))).or_default().push(t);
        }
    }
    let ids: BTreeMap<(u32, u32), usize> = edges.keys().enumerate().map(|(i, &e)| (e, i)).collect();
    let mut rows: Vec<Vec<bool>> = Vec::new();
    let mut cut = BTreeSet::new();
    for l in loops {
        assert!(l.len() >= 3, "a loop has at least three edges: {l:?}");
        let distinct: BTreeSet<u32> = l.iter().copied().collect();
        assert_eq!(distinct.len(), l.len(), "a simple loop: {l:?}");
        let mut row = vec![false; ids.len()];
        for k in 0..l.len() {
            let (a, b) = (l[k], l[(k + 1) % l.len()]);
            let edge = ids
                .get(&(a.min(b), a.max(b)))
                .unwrap_or_else(|| panic!("{a}-{b} is not a mesh edge"));
            row[*edge] ^= true;
            cut.insert((a.min(b), a.max(b)));
        }
        rows.push(row);
    }
    // Independent over GF(2).
    let mut rank = 0;
    let columns = ids.len();
    for col in 0..columns {
        if let Some(pivot) = (rank..rows.len()).find(|&r| rows[r][col]) {
            rows.swap(rank, pivot);
            for r in 0..rows.len() {
                if r != rank && rows[r][col] {
                    let pivot_row = rows[rank].clone();
                    for (x, y) in rows[r].iter_mut().zip(pivot_row) {
                        *x ^= y;
                    }
                }
            }
            rank += 1;
        }
    }
    assert_eq!(rank, loops.len(), "loops independent as chains");
    // Non-separating: triangles stay connected across uncut edges.
    let triangles = mesh.indices.len() / 3;
    let mut seen = vec![false; triangles];
    let mut stack = vec![0];
    seen[0] = true;
    while let Some(t) = stack.pop() {
        for (edge, uses) in &edges {
            if uses.contains(&t) && !cut.contains(edge) {
                for &n in uses {
                    if !seen[n] {
                        seen[n] = true;
                        stack.push(n);
                    }
                }
            }
        }
    }
    assert!(seen.iter().all(|&s| s), "the loops separate the surface");
}

#[test]
fn a_cube_is_a_sphere_with_an_empty_basis() {
    let mesh = voxels(&[[0, 0, 0]]);
    let c = only(&mesh);
    assert_eq!((c.vertices, c.edges, c.triangles.len()), (8, 18, 12));
    assert_eq!(c.euler_characteristic, 2);
    assert_eq!(c.boundary_loops, 0);
    assert!(c.orientable && c.consistently_oriented);
    assert_eq!(c.surface, SurfaceKind::Orientable { genus: 0 });
    assert_basis(&mesh, &c, 0);
}

#[test]
fn a_ring_of_voxels_is_a_torus_with_two_generators() {
    let mesh = slab_with_holes(3, &[1]);
    let c = only(&mesh);
    assert_eq!(c.euler_characteristic, 0);
    assert_eq!(c.surface, SurfaceKind::Orientable { genus: 1 });
    assert_basis(&mesh, &c, 1);
}

#[test]
fn a_slab_with_two_holes_is_a_double_torus_with_four_generators() {
    let mesh = slab_with_holes(5, &[1, 3]);
    let c = only(&mesh);
    assert_eq!(c.euler_characteristic, -2);
    assert_eq!(c.surface, SurfaceKind::Orientable { genus: 2 });
    assert_basis(&mesh, &c, 2);
    let mesh = slab_with_holes(7, &[1, 3, 5]);
    let c = only(&mesh);
    assert_eq!(c.surface, SurfaceKind::Orientable { genus: 3 });
    assert_basis(&mesh, &c, 3);
}

#[test]
fn a_triangle_is_a_disc() {
    let mesh = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        vec![0, 1, 2],
    );
    let c = only(&mesh);
    assert_eq!((c.euler_characteristic, c.boundary_loops), (1, 1));
    assert_eq!(c.surface, SurfaceKind::Orientable { genus: 0 });
    assert_eq!(c.homology_basis, None, "no basis with boundary");
}

#[test]
fn a_square_ring_is_an_annulus() {
    let cells: Vec<[i64; 2]> = (0..3)
        .flat_map(|x| (0..3).map(move |y| [x, y]))
        .filter(|&c| c != [1, 1])
        .collect();
    let c = only(&pixels(&cells));
    assert_eq!((c.euler_characteristic, c.boundary_loops), (0, 2));
    assert_eq!(c.surface, SurfaceKind::Orientable { genus: 0 });
}

#[test]
fn a_twisted_strip_is_a_mobius_band() {
    let c = only(&glued_grid(6, 1, true, false));
    assert_eq!((c.euler_characteristic, c.boundary_loops), (0, 1));
    assert!(!c.orientable && !c.consistently_oriented);
    assert_eq!(c.surface, SurfaceKind::NonOrientable { crosscaps: 1 });
    assert_eq!(c.homology_basis, None);
    // The same strip untwisted is an annulus.
    let c = only(&glued_grid(6, 1, false, false));
    assert_eq!(c.boundary_loops, 2);
    assert_eq!(c.surface, SurfaceKind::Orientable { genus: 0 });
}

#[test]
fn a_twisted_torus_is_a_klein_bottle() {
    let c = only(&glued_grid(6, 6, true, true));
    assert_eq!((c.euler_characteristic, c.boundary_loops), (0, 0));
    assert_eq!(c.surface, SurfaceKind::NonOrientable { crosscaps: 2 });
    let mesh = glued_grid(6, 6, false, true);
    let c = only(&mesh);
    assert_eq!(c.surface, SurfaceKind::Orientable { genus: 1 });
    assert_basis(&mesh, &c, 1);
}

#[test]
fn an_inconsistently_wound_cube_is_still_orientable() {
    let mut mesh = voxels(&[[0, 0, 0]]);
    mesh.indices.swap(0, 1);
    let c = only(&mesh);
    assert!(c.orientable);
    assert!(!c.consistently_oriented);
    assert_eq!(c.surface, SurfaceKind::Orientable { genus: 0 });
}

#[test]
fn components_are_reported_separately() {
    // A cube, a torus and a triangle, with an unused position between.
    let mut mesh = voxels(&[[0, 0, 0]]);
    let torus = slab_with_holes(3, &[1]);
    mesh.positions.push(Point3::new(9.0, 9.0, 9.0));
    let base = mesh.positions.len() as u32;
    mesh.positions.extend(torus.positions.iter().copied());
    mesh.indices.extend(torus.indices.iter().map(|i| i + base));
    let base = mesh.positions.len() as u32;
    mesh.positions.extend([
        Point3::new(20.0, 0.0, 0.0),
        Point3::new(21.0, 0.0, 0.0),
        Point3::new(20.0, 1.0, 0.0),
    ]);
    mesh.indices.extend([base, base + 1, base + 2]);
    let report = topology(&mesh).unwrap();
    let kinds: Vec<(SurfaceKind, usize)> = report
        .components
        .iter()
        .map(|c| (c.surface, c.boundary_loops))
        .collect();
    assert_eq!(
        kinds,
        [
            (SurfaceKind::Orientable { genus: 0 }, 0),
            (SurfaceKind::Orientable { genus: 1 }, 0),
            (SurfaceKind::Orientable { genus: 0 }, 1),
        ]
    );
    let total: usize = report.components.iter().map(|c| c.triangles.len()).sum();
    assert_eq!(total, mesh.indices.len() / 3);
}

#[test]
fn non_manifold_meshes_are_refused() {
    // Two voxels sharing only an edge: that edge is on four triangles.
    assert!(matches!(
        topology(&voxels(&[[0, 0, 0], [1, 1, 0]])),
        Err(TopologyError::NonManifold { edges: 1, .. })
    ));
    // Two voxels sharing only a corner: the corner has two fans.
    assert_eq!(
        topology(&voxels(&[[0, 0, 0], [1, 1, 1]])),
        Err(TopologyError::NonManifold {
            edges: 0,
            vertices: 1
        })
    );
    // Three triangles on one edge.
    let fin = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ],
        vec![0, 1, 2, 1, 0, 3, 0, 1, 4],
    );
    assert!(matches!(
        topology(&fin),
        Err(TopologyError::NonManifold { edges: 1, .. })
    ));
}

#[test]
fn malformed_triangles_are_refused() {
    let p = vec![Point3::new(0.0, 0.0, 0.0); 3];
    assert_eq!(
        topology(&TriMesh::new(p.clone(), vec![0, 1, 1])),
        Err(TopologyError::DegenerateTriangle { triangle: 0 })
    );
    assert!(matches!(
        topology(&TriMesh::new(p, vec![0, 1, 3])),
        Err(TopologyError::IndexOutOfRange { vertex: 3, .. })
    ));
}
