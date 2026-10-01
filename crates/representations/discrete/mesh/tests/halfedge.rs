//! Halfedge mesh (#140, ledger row D1).
//!
//! The suite pins the structure's contract for the algorithms that will
//! edit it in place: exact round-trips through `TriMesh`, O(1) navigation
//! that agrees with the input, counter-clockwise circulation, boundary
//! loops, refusal of non-manifold input by name, and local edits that keep
//! every invariant and the Euler characteristic. The random-edit tests
//! replay fixed-seed sequences and, after every step, check the invariants
//! twice: with `validate`, and independently by converting to a `TriMesh`
//! and rebuilding, which refuses anything non-manifold.

use axiolid_core::Point3;
use axiolid_mesh::{
    EdgeId, FaceId, HalfedgeBuildError, HalfedgeEditError, HalfedgeId, HalfedgeMesh, TriMesh,
    VertexId,
};

// ---- fixtures ---------------------------------------------------------------

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// Outward unit tetrahedron.
fn tetrahedron() -> TriMesh {
    TriMesh::new(
        vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, 0.0, 1.0),
        ],
        vec![0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3],
    )
}

/// Tetrahedron without its base: three triangles around an apex.
fn open_tetrahedron() -> TriMesh {
    TriMesh::new(
        vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, 0.0, 1.0),
        ],
        vec![0, 1, 3, 0, 3, 2, 1, 2, 3],
    )
}

fn lone_triangle() -> TriMesh {
    TriMesh::new(
        vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0)],
        vec![0, 1, 2],
    )
}

/// Octahedron, outward wound.
fn octahedron() -> TriMesh {
    let positions = vec![
        p(1.0, 0.0, 0.0),
        p(0.0, 1.0, 0.0),
        p(-1.0, 0.0, 0.0),
        p(0.0, -1.0, 0.0),
        p(0.0, 0.0, 1.0),
        p(0.0, 0.0, -1.0),
    ];
    let mut indices = Vec::new();
    for i in 0..4u32 {
        let j = (i + 1) % 4;
        indices.extend([i, j, 4, j, i, 5]);
    }
    TriMesh::new(positions, indices)
}

/// `nx` by `ny` cells in the xy plane, two counter-clockwise triangles each.
fn grid(nx: u32, ny: u32) -> TriMesh {
    let mut positions = Vec::new();
    for y in 0..=ny {
        for x in 0..=nx {
            positions.push(p(f64::from(x), f64::from(y), 0.0));
        }
    }
    let at = |x: u32, y: u32| y * (nx + 1) + x;
    let mut indices = Vec::new();
    for y in 0..ny {
        for x in 0..nx {
            let (a, b, c, d) = (at(x, y), at(x + 1, y), at(x + 1, y + 1), at(x, y + 1));
            indices.extend([a, b, c, a, c, d]);
        }
    }
    TriMesh::new(positions, indices)
}

/// Closed torus of `n` by `m` quads, two triangles each: genus one.
fn torus(n: u32, m: u32) -> TriMesh {
    let mut positions = Vec::new();
    for i in 0..n {
        for j in 0..m {
            let u = f64::from(i) / f64::from(n) * core::f64::consts::TAU;
            let v = f64::from(j) / f64::from(m) * core::f64::consts::TAU;
            let r = 2.0 + v.cos();
            positions.push(p(r * u.cos(), r * u.sin(), v.sin()));
        }
    }
    let at = |i: u32, j: u32| (i % n) * m + (j % m);
    let mut indices = Vec::new();
    for i in 0..n {
        for j in 0..m {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            indices.extend([a, b, c, a, c, d]);
        }
    }
    TriMesh::new(positions, indices)
}

/// Two triangles on the same three vertices, wound oppositely: a closed
/// sphere with two faces.
fn pillow() -> TriMesh {
    TriMesh::new(
        vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0)],
        vec![0, 1, 2, 0, 2, 1],
    )
}

fn build(mesh: &TriMesh) -> HalfedgeMesh {
    HalfedgeMesh::from_tri_mesh(mesh).expect("fixture is manifold")
}

/// Deterministic SplitMix64, so every run replays the same edit sequence.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// Every navigation identity the consumers rely on, checked from the
/// public API rather than from `validate`.
fn assert_navigation(mesh: &HalfedgeMesh) {
    let mut degree_sum = 0;
    for v in mesh.vertices() {
        for h in mesh.outgoing_halfedges(v) {
            assert_eq!(mesh.source(h), v);
            degree_sum += 1;
        }
        if mesh.is_boundary_vertex(v) {
            let first = mesh.outgoing_halfedges(v).next().unwrap();
            assert!(
                mesh.is_boundary_halfedge(first),
                "boundary vertex starts at its hole"
            );
            assert_eq!(
                mesh.outgoing_halfedges(v)
                    .filter(|&h| mesh.is_boundary_halfedge(h))
                    .count(),
                1
            );
        }
    }
    assert_eq!(degree_sum, mesh.halfedge_count());
    for h in mesh.halfedges() {
        assert_eq!(mesh.prev(mesh.next(h)), h);
        assert_eq!(mesh.next(mesh.prev(h)), h);
        assert_eq!(mesh.opposite(mesh.opposite(h)), h);
        assert_ne!(mesh.opposite(h), h);
        assert_eq!(mesh.source(mesh.next(h)), mesh.target(h));
        assert_eq!(mesh.target(mesh.opposite(h)), mesh.source(h));
        assert_eq!(mesh.face(mesh.next(h)), mesh.face(h));
        assert_eq!(mesh.edge(mesh.opposite(h)), mesh.edge(h));
        assert!(!(mesh.is_boundary_halfedge(h) && mesh.is_boundary_halfedge(mesh.opposite(h))));
    }
    let mut face_sides = 0;
    for f in mesh.faces() {
        assert_eq!(mesh.face(mesh.face_halfedge(f)), Some(f));
        face_sides += mesh.face_degree(f);
    }
    let boundary = mesh.boundary_halfedges().count();
    assert_eq!(face_sides + boundary, mesh.halfedge_count());
    let loop_sides: usize = mesh
        .boundary_loops()
        .into_iter()
        .map(|h| mesh.loop_halfedges(h).count())
        .sum();
    assert_eq!(loop_sides, boundary);
}

/// `validate`, the navigation identities, and an independent rebuild.
fn assert_sound(mesh: &HalfedgeMesh) {
    mesh.validate().expect("invariants hold");
    assert_navigation(mesh);
    if mesh.faces().all(|f| mesh.face_degree(f) == 3) {
        let tri = mesh.to_tri_mesh().expect("all triangles");
        let rebuilt = HalfedgeMesh::from_tri_mesh(&tri).expect("still a manifold");
        assert_eq!(rebuilt.euler_characteristic(), mesh.euler_characteristic());
        assert_eq!(rebuilt.boundary_loops().len(), mesh.boundary_loops().len());
    }
}

fn triangle_set(mesh: &TriMesh) -> Vec<[u32; 3]> {
    // Rotate each triangle to start at its smallest corner, keeping winding.
    let mut set: Vec<[u32; 3]> = mesh
        .triangles()
        .map(|t| {
            let k = (0..3).min_by_key(|&k| t[k]).unwrap();
            [t[k], t[(k + 1) % 3], t[(k + 2) % 3]]
        })
        .collect();
    set.sort_unstable();
    set
}

// ---- round trips --------------------------------------------------------

#[test]
fn every_fixture_round_trips_exactly() {
    for mesh in [
        tetrahedron(),
        open_tetrahedron(),
        lone_triangle(),
        octahedron(),
        grid(4, 3),
        torus(6, 5),
        pillow(),
    ] {
        let halfedge = build(&mesh);
        assert_sound(&halfedge);
        assert_eq!(halfedge.to_tri_mesh().unwrap(), mesh);
    }
}

#[test]
fn ids_follow_the_input_numbering() {
    let mesh = tetrahedron();
    let halfedge = build(&mesh);
    for (j, triangle) in mesh.triangles().enumerate() {
        let corners: Vec<u32> = halfedge
            .face_vertices(FaceId::new(j as u32))
            .map(VertexId::get)
            .collect();
        assert_eq!(corners, triangle);
    }
    for v in halfedge.vertices() {
        assert_eq!(halfedge.position(v), mesh.positions[v.index()]);
    }
}

#[test]
fn an_isolated_vertex_survives_the_round_trip() {
    let mut mesh = grid(1, 1);
    mesh.positions.push(p(9.0, 9.0, 9.0));
    let halfedge = build(&mesh);
    let lonely = VertexId::new(4);
    assert!(halfedge.is_isolated(lonely));
    assert!(!halfedge.is_boundary_vertex(lonely));
    assert_eq!(halfedge.outgoing_halfedges(lonely).count(), 0);
    assert_eq!(halfedge.vertex_count(), 5);
    // Five vertices, five edges, two faces.
    assert_eq!(halfedge.euler_characteristic(), 2);
    assert_eq!(halfedge.to_tri_mesh().unwrap(), mesh);
}

#[test]
fn polygon_faces_round_trip_after_triangulation_only() {
    let positions = vec![
        p(0.0, 0.0, 0.0),
        p(1.0, 0.0, 0.0),
        p(1.0, 1.0, 0.0),
        p(0.0, 1.0, 0.0),
    ];
    let mut quad = HalfedgeMesh::from_faces(positions, &[[0u32, 1, 2, 3]]).unwrap();
    assert_sound(&quad);
    assert_eq!(
        quad.to_tri_mesh(),
        Err(HalfedgeEditError::NotATriangle {
            face: FaceId::new(0),
            degree: 4,
        })
    );
    let h = quad.face_halfedge(FaceId::new(0)); // 0 -> 1
    let across = quad.next(quad.next(h)); // 2 -> 3
    let diagonal = quad.split_face_diagonal(h, across).unwrap();
    assert_eq!(quad.source(diagonal), VertexId::new(1));
    assert_eq!(quad.target(diagonal), VertexId::new(3));
    assert_sound(&quad);
    assert_eq!(quad.face_count(), 2);
    assert_eq!(
        triangle_set(&quad.to_tri_mesh().unwrap()),
        vec![[0, 1, 3], [1, 2, 3]]
    );
}

// ---- navigation and circulation ---------------------------------------------

#[test]
fn closed_meshes_have_no_boundary_and_their_euler_characteristic() {
    let tetra = build(&tetrahedron());
    assert_eq!(tetra.euler_characteristic(), 2);
    assert!(tetra.boundary_loops().is_empty());
    for v in tetra.vertices() {
        assert_eq!(tetra.degree(v), 3);
        assert_eq!(tetra.vertex_faces(v).count(), 3);
    }
    for h in tetra.halfedges() {
        assert_eq!(tetra.next(tetra.next(tetra.next(h))), h);
    }
    let torus = build(&torus(6, 5));
    assert_eq!(torus.euler_characteristic(), 0);
    assert!(torus.vertices().all(|v| torus.degree(v) == 6));
    assert_eq!(build(&octahedron()).euler_characteristic(), 2);
}

#[test]
fn outgoing_halfedges_turn_counter_clockwise() {
    let mesh = build(&grid(2, 2));
    let centre = VertexId::new(4);
    let out: Vec<HalfedgeId> = mesh.outgoing_halfedges(centre).collect();
    assert_eq!(out.len(), 6);
    for i in 0..out.len() {
        let a = mesh.position(mesh.target(out[i])) - mesh.position(centre);
        let b = mesh.position(mesh.target(out[(i + 1) % out.len()])) - mesh.position(centre);
        assert!(a.cross(b).z > 0.0, "step {i} turns clockwise");
    }
    // Every face around the vertex is reached once, in the same order.
    let faces: Vec<FaceId> = mesh.vertex_faces(centre).collect();
    assert_eq!(faces.len(), 6);
    let mut sorted = faces.clone();
    sorted.sort_unstable();
    sorted.dedup();
    assert_eq!(sorted.len(), 6);
}

#[test]
fn a_boundary_vertex_circulates_from_its_hole() {
    let mesh = build(&grid(2, 2));
    let corner = VertexId::new(0);
    let out: Vec<HalfedgeId> = mesh.outgoing_halfedges(corner).collect();
    assert!(mesh.is_boundary_halfedge(out[0]));
    assert!(out[1..].iter().all(|&h| !mesh.is_boundary_halfedge(h)));
    // Up the left side first (the boundary halfedge), then on round
    // counter-clockwise: across the hole to the bottom side, then the
    // diagonal between the two faces.
    let neighbours: Vec<u32> = mesh.vertex_vertices(corner).map(VertexId::get).collect();
    assert_eq!(neighbours, vec![3, 1, 4]);
    assert_eq!(mesh.vertex_faces(corner).count(), 2);
}

#[test]
fn boundary_loops_run_against_the_faces() {
    let triangle = build(&lone_triangle());
    let loops = triangle.boundary_loops();
    assert_eq!(loops.len(), 1);
    let hole: Vec<HalfedgeId> = triangle.loop_halfedges(loops[0]).collect();
    assert_eq!(hole.len(), 3);
    for &h in &hole {
        assert!(triangle.is_boundary_halfedge(h));
        assert_eq!(triangle.face(triangle.opposite(h)), Some(FaceId::new(0)));
        assert!(triangle.is_boundary_edge(triangle.edge(h)));
    }
    let face_faces: Vec<Option<FaceId>> = triangle.face_faces(FaceId::new(0)).collect();
    assert_eq!(face_faces, vec![None, None, None]);

    let mesh = build(&grid(3, 2));
    let loops = mesh.boundary_loops();
    assert_eq!(loops.len(), 1);
    assert_eq!(mesh.loop_halfedges(loops[0]).count(), 2 * (3 + 2));
    let interior = VertexId::new(5);
    assert!(!mesh.is_boundary_vertex(interior));
    assert!(mesh.is_boundary_vertex(VertexId::new(1)));
}

#[test]
fn find_halfedge_and_edge_vertices_agree() {
    let mesh = build(&octahedron());
    for e in mesh.edges() {
        let [a, b] = mesh.edge_vertices(e);
        let h = mesh.find_halfedge(a, b).unwrap();
        assert_eq!(mesh.edge(h), e);
        assert_eq!(mesh.find_halfedge(b, a), Some(mesh.opposite(h)));
    }
    assert_eq!(mesh.find_halfedge(VertexId::new(4), VertexId::new(5)), None);
}

// ---- refusals -----------------------------------------------------------------

fn refusal(positions: usize, indices: Vec<u32>) -> HalfedgeBuildError {
    let positions = (0..positions).map(|i| p(i as f64, 0.0, 0.0)).collect();
    HalfedgeMesh::from_tri_mesh(&TriMesh::new(positions, indices)).unwrap_err()
}

#[test]
fn an_edge_shared_by_three_faces_is_refused_by_name() {
    assert_eq!(
        refusal(5, vec![0, 1, 2, 1, 0, 3, 0, 1, 4]),
        HalfedgeBuildError::NonManifoldEdge {
            a: 0,
            b: 1,
            faces: 3
        }
    );
}

#[test]
fn faces_disagreeing_on_winding_are_refused_by_name() {
    assert_eq!(
        refusal(4, vec![0, 1, 2, 0, 1, 3]),
        HalfedgeBuildError::InconsistentOrientation { from: 0, to: 1 }
    );
    // The same face twice is the same disagreement.
    assert_eq!(
        refusal(3, vec![0, 1, 2, 0, 1, 2]),
        HalfedgeBuildError::InconsistentOrientation { from: 0, to: 1 }
    );
}

#[test]
fn a_bowtie_vertex_is_refused_by_name() {
    assert_eq!(
        refusal(5, vec![0, 1, 2, 0, 3, 4]),
        HalfedgeBuildError::NonManifoldVertex { vertex: 0 }
    );
}

#[test]
fn two_holes_meeting_at_a_vertex_are_refused_by_name() {
    // A 3x3 grid missing cells (0, 0) and (1, 1): the holes touch only at
    // grid vertex (1, 1).
    let mut mesh = grid(3, 3);
    let cell = |x: usize, y: usize| (y * 3 + x) * 6;
    for start in [cell(1, 1), cell(0, 0)] {
        mesh.indices.drain(start..start + 6);
    }
    assert_eq!(
        HalfedgeMesh::from_tri_mesh(&mesh).unwrap_err(),
        HalfedgeBuildError::NonManifoldVertex { vertex: 5 }
    );
    // Three fans at one vertex.
    assert_eq!(
        refusal(7, vec![0, 1, 2, 0, 3, 4, 0, 5, 6]),
        HalfedgeBuildError::NonManifoldVertex { vertex: 0 }
    );
}

#[test]
fn two_closed_fans_at_one_vertex_are_refused_by_name() {
    // Two tetrahedra touching at vertex 0: every edge is fine, the vertex
    // is not.
    let mut indices = vec![0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
    indices.extend([0, 5, 4, 0, 4, 6, 0, 6, 5, 4, 5, 6]);
    assert_eq!(
        refusal(7, indices),
        HalfedgeBuildError::NonManifoldVertex { vertex: 0 }
    );
}

#[test]
fn malformed_faces_are_refused_by_name() {
    assert_eq!(
        refusal(3, vec![0, 1, 2, 0]),
        HalfedgeBuildError::IncompleteTriangle { index_count: 4 }
    );
    assert_eq!(
        refusal(3, vec![0, 1, 7]),
        HalfedgeBuildError::IndexOutOfRange {
            face: 0,
            index: 7,
            position_count: 3,
        }
    );
    assert_eq!(
        refusal(3, vec![0, 1, 2, 2, 1, 1]),
        HalfedgeBuildError::DegenerateFace { face: 1, vertex: 1 }
    );
    let positions = vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0)];
    assert_eq!(
        HalfedgeMesh::from_faces(positions, &[vec![0u32, 1]]).unwrap_err(),
        HalfedgeBuildError::FaceTooSmall {
            face: 0,
            corners: 2
        }
    );
}

#[test]
fn the_empty_mesh_is_a_valid_mesh() {
    let empty = build(&TriMesh::default());
    assert_sound(&empty);
    assert_eq!(empty.euler_characteristic(), 0);
    assert_eq!(empty.to_tri_mesh().unwrap(), TriMesh::default());
}

// ---- flips --------------------------------------------------------------------

fn interior_edge(mesh: &HalfedgeMesh) -> EdgeId {
    mesh.edges().find(|&e| !mesh.is_boundary_edge(e)).unwrap()
}

#[test]
fn flipping_twice_restores_the_triangles() {
    let mut mesh = build(&grid(1, 1));
    let e = interior_edge(&mesh);
    let before = triangle_set(&mesh.to_tri_mesh().unwrap());
    assert_eq!(before, vec![[0, 1, 3], [0, 3, 2]]);
    mesh.flip_edge(e).unwrap();
    assert_sound(&mesh);
    let flipped = triangle_set(&mesh.to_tri_mesh().unwrap());
    assert_eq!(flipped, vec![[0, 1, 2], [1, 3, 2]]);
    let mut endpoints = mesh.edge_vertices(e).map(VertexId::get);
    endpoints.sort_unstable();
    assert_eq!(endpoints, [1, 2]);
    mesh.flip_edge(e).unwrap();
    assert_sound(&mesh);
    assert_eq!(triangle_set(&mesh.to_tri_mesh().unwrap()), before);
}

#[test]
fn flipping_every_interior_edge_keeps_the_structure() {
    let mut mesh = build(&grid(3, 3));
    let interior: Vec<EdgeId> = mesh
        .edges()
        .filter(|&e| !mesh.is_boundary_edge(e))
        .collect();
    for e in interior {
        if mesh.flip_edge(e).is_ok() {
            mesh.validate().unwrap();
        }
    }
    assert_sound(&mesh);
}

#[test]
fn flips_that_would_break_the_surface_are_refused() {
    let mut mesh = build(&grid(1, 1));
    let boundary = mesh.edges().find(|&e| mesh.is_boundary_edge(e)).unwrap();
    assert_eq!(
        mesh.flip_edge(boundary),
        Err(HalfedgeEditError::BoundaryEdge { edge: boundary })
    );

    let mut tetra = build(&tetrahedron());
    let snapshot = tetra.clone();
    let e = EdgeId::new(0);
    assert!(matches!(
        tetra.flip_edge(e),
        Err(HalfedgeEditError::EdgeExists { .. })
    ));
    assert_eq!(tetra, snapshot, "a refused flip changes nothing");

    let mut pillow = build(&pillow());
    assert_eq!(
        pillow.flip_edge(e),
        Err(HalfedgeEditError::WouldDegenerate { edge: e })
    );
}

// ---- splits -------------------------------------------------------------------

#[test]
fn splitting_an_interior_edge_keeps_triangles_and_euler() {
    let mut mesh = build(&grid(2, 2));
    let chi = mesh.euler_characteristic();
    let (v, e, f) = (mesh.vertex_count(), mesh.edge_count(), mesh.face_count());
    let edge = interior_edge(&mesh);
    let [a, b] = mesh.edge_vertices(edge);
    let mid = (mesh.position(a) + mesh.position(b)) * 0.5;
    let m = mesh.split_edge(edge, mid).unwrap();
    assert_sound(&mesh);
    assert_eq!(mesh.position(m), mid);
    assert_eq!(
        (mesh.vertex_count(), mesh.edge_count(), mesh.face_count()),
        (v + 1, e + 3, f + 2)
    );
    assert_eq!(mesh.euler_characteristic(), chi);
    assert_eq!(mesh.degree(m), 4);
    assert_eq!(mesh.edge_vertices(edge), [a, m]);
}

#[test]
fn splitting_a_boundary_edge_puts_the_new_vertex_on_the_boundary() {
    let mut mesh = build(&lone_triangle());
    let m = mesh.split_edge(EdgeId::new(0), p(0.5, 0.0, 0.0)).unwrap();
    assert_sound(&mesh);
    assert!(mesh.is_boundary_vertex(m));
    assert_eq!(mesh.degree(m), 3);
    assert_eq!(
        (mesh.vertex_count(), mesh.edge_count(), mesh.face_count()),
        (4, 5, 2)
    );
    assert_eq!(mesh.euler_characteristic(), 1);
    assert_eq!(mesh.loop_halfedges(mesh.boundary_loops()[0]).count(), 4);
}

#[test]
fn splitting_a_boundary_edge_works_from_either_side() {
    // A build always puts halfedge 2e on a face. Collapsing 1 into 4 moves
    // halfedge 4 -> 0 (the 2e of edge 0-4) onto the hole, so that edge has
    // its boundary on the other side.
    let mut mesh = build(&grid(2, 2));
    let h = mesh
        .find_halfedge(VertexId::new(1), VertexId::new(4))
        .unwrap();
    mesh.collapse_edge(h).unwrap();
    let flipped_side = mesh.edge(
        mesh.find_halfedge(VertexId::new(4), VertexId::new(0))
            .unwrap(),
    );
    assert!(mesh.is_boundary_halfedge(mesh.edge_halfedge(flipped_side, 0)));
    let usual_side = mesh.edge(
        mesh.find_halfedge(VertexId::new(2), VertexId::new(5))
            .unwrap(),
    );
    assert!(mesh.is_boundary_halfedge(mesh.edge_halfedge(usual_side, 1)));
    for e in [flipped_side, usual_side] {
        let m = mesh.split_edge(e, p(0.0, 0.0, 1.0)).unwrap();
        assert_sound(&mesh);
        assert!(mesh.is_boundary_vertex(m));
        assert_eq!(mesh.degree(m), 3);
    }
}

#[test]
fn splitting_a_polygon_edge_adds_a_side_instead_of_a_triangle() {
    let positions = vec![
        p(0.0, 0.0, 0.0),
        p(1.0, 0.0, 0.0),
        p(1.0, 1.0, 0.0),
        p(0.0, 1.0, 0.0),
    ];
    let mut quad = HalfedgeMesh::from_faces(positions, &[[0u32, 1, 2, 3]]).unwrap();
    quad.split_edge(EdgeId::new(0), p(0.5, 0.0, 0.0)).unwrap();
    assert_sound(&quad);
    assert_eq!(quad.face_count(), 1);
    assert_eq!(quad.face_degree(FaceId::new(0)), 5);
}

#[test]
fn splitting_a_face_fans_it_around_a_new_vertex() {
    let mut mesh = build(&tetrahedron());
    let c = mesh.split_face(FaceId::new(0), p(0.3, 0.3, 0.0)).unwrap();
    assert_sound(&mesh);
    assert_eq!(
        (mesh.vertex_count(), mesh.edge_count(), mesh.face_count()),
        (5, 9, 6)
    );
    assert_eq!(mesh.euler_characteristic(), 2);
    assert_eq!(mesh.degree(c), 3);

    let positions = vec![
        p(0.0, 0.0, 0.0),
        p(1.0, 0.0, 0.0),
        p(1.0, 1.0, 0.0),
        p(0.0, 1.0, 0.0),
    ];
    let mut quad = HalfedgeMesh::from_faces(positions, &[[0u32, 1, 2, 3]]).unwrap();
    let c = quad.split_face(FaceId::new(0), p(0.5, 0.5, 0.0)).unwrap();
    assert_sound(&quad);
    assert_eq!(quad.face_count(), 4);
    assert_eq!(quad.euler_characteristic(), 1);
    assert_eq!(
        triangle_set(&quad.to_tri_mesh().unwrap()),
        vec![[0, 1, 4], [0, 4, 3], [1, 2, 4], [2, 3, 4]]
    );
    assert_eq!(c, VertexId::new(4));
}

#[test]
fn face_diagonals_that_would_duplicate_a_side_are_refused() {
    let positions = (0..5).map(|i| p(f64::from(i), 0.0, 0.0)).collect();
    let mut pentagon = HalfedgeMesh::from_faces(positions, &[[0u32, 1, 2, 3, 4]]).unwrap();
    let h = pentagon.face_halfedge(FaceId::new(0));
    assert_eq!(
        pentagon.split_face_diagonal(h, h),
        Err(HalfedgeEditError::AdjacentCorners { a: h, b: h })
    );
    let n = pentagon.next(h);
    assert_eq!(
        pentagon.split_face_diagonal(h, n),
        Err(HalfedgeEditError::AdjacentCorners { a: h, b: n })
    );
    assert_eq!(
        pentagon.split_face_diagonal(n, h),
        Err(HalfedgeEditError::AdjacentCorners { a: n, b: h })
    );
    let hole = pentagon.opposite(h);
    assert_eq!(
        pentagon.split_face_diagonal(h, hole),
        Err(HalfedgeEditError::NotSameFace { a: h, b: hole })
    );
    // 1 -> 3 splits off triangle (1, 2, 3); 3 -> 1 again would duplicate it.
    let d = pentagon
        .split_face_diagonal(h, pentagon.next(pentagon.next(h)))
        .unwrap();
    assert_sound(&pentagon);
    let rest = pentagon.face(pentagon.opposite(d)).unwrap();
    assert_eq!(pentagon.face_degree(rest), 3);
    let quad = pentagon.face(d).unwrap();
    assert_eq!(pentagon.face_degree(quad), 4);
}

#[test]
fn a_diagonal_through_another_face_is_refused() {
    // Quads (0, 1, 2, 3) and (3, 2, 4, 0) share sides 2-3 and 0-3. Joining
    // 0 to 2 across the first leaves the second unable to take the same
    // diagonal.
    let positions = (0..5)
        .map(|i| p(f64::from(i), f64::from(i * i), 0.0))
        .collect();
    let faces: Vec<Vec<u32>> = vec![vec![0, 1, 2, 3], vec![3, 2, 4, 0]];
    let mut mesh = HalfedgeMesh::from_faces(positions, &faces).unwrap();
    assert_sound(&mesh);
    let f0 = mesh.face_halfedge(FaceId::new(0)); // 0 -> 1
    let into_2 = mesh.next(f0); // 1 -> 2
    let into_0 = mesh.prev(f0); // 3 -> 0
    mesh.split_face_diagonal(into_0, into_2).unwrap();
    assert_sound(&mesh);
    let to_2 = mesh.face_halfedge(FaceId::new(1)); // 3 -> 2
    let to_0 = mesh.next(mesh.next(to_2)); // 4 -> 0
    assert_eq!(
        mesh.split_face_diagonal(to_2, to_0),
        Err(HalfedgeEditError::EdgeExists {
            a: VertexId::new(2),
            b: VertexId::new(0),
        })
    );
}

// ---- collapses ----------------------------------------------------------------

#[test]
fn collapsing_an_octahedron_edge_keeps_a_sphere() {
    let mut mesh = build(&octahedron());
    let h = mesh
        .find_halfedge(VertexId::new(4), VertexId::new(0))
        .unwrap();
    let kept = mesh.collapse_edge(h).unwrap();
    assert_eq!(kept, VertexId::new(0));
    assert!(!mesh.contains_vertex(VertexId::new(4)));
    assert_sound(&mesh);
    assert_eq!(
        (mesh.vertex_count(), mesh.edge_count(), mesh.face_count()),
        (5, 9, 6)
    );
    assert_eq!(mesh.euler_characteristic(), 2);
    assert_eq!(
        mesh.collapse_edge(h),
        Err(HalfedgeEditError::RemovedElement {
            element: "edge",
            index: h.get() / 2,
        })
    );
}

#[test]
fn collapsing_a_boundary_edge_shortens_the_hole() {
    let mut mesh = build(&grid(2, 2));
    let h = mesh
        .find_halfedge(VertexId::new(1), VertexId::new(0))
        .unwrap();
    let before = mesh.loop_halfedges(mesh.boundary_loops()[0]).count();
    mesh.collapse_edge(h).unwrap();
    assert_sound(&mesh);
    assert_eq!(mesh.euler_characteristic(), 1);
    assert_eq!(
        mesh.loop_halfedges(mesh.boundary_loops()[0]).count(),
        before - 1
    );
    assert!(mesh.is_boundary_vertex(VertexId::new(0)));
}

#[test]
fn collapsing_toward_the_boundary_moves_the_hole_onto_the_kept_vertex() {
    // Interior vertex 4 of a 2x2 grid collapsed onto boundary vertex 1 and
    // the reverse: both keep a valid boundary.
    for (from, to) in [(4, 1), (1, 4)] {
        let mut mesh = build(&grid(2, 2));
        let h = mesh
            .find_halfedge(VertexId::new(from), VertexId::new(to))
            .unwrap();
        let kept = mesh.collapse_edge(h).unwrap();
        assert_sound(&mesh);
        assert!(mesh.is_boundary_vertex(kept));
        assert_eq!(mesh.euler_characteristic(), 1);
    }
}

#[test]
fn a_collapse_violating_the_link_condition_is_refused() {
    // Open tetrahedron: boundary edge 0-1 has apex 3 opposite, but 0 and 1
    // also share boundary neighbour 2.
    let mut open = build(&open_tetrahedron());
    let snapshot = open.clone();
    let h = open
        .find_halfedge(VertexId::new(0), VertexId::new(1))
        .unwrap();
    let edge = open.edge(h);
    assert_eq!(
        open.collapse_edge(h),
        Err(HalfedgeEditError::LinkCondition { edge })
    );
    assert_eq!(open, snapshot);

    // A strip one cell wide: the diagonal joins two boundary vertices but
    // is interior, so collapsing it would pinch the strip.
    let mut strip = build(&grid(2, 1));
    let h = strip
        .find_halfedge(VertexId::new(0), VertexId::new(4))
        .unwrap();
    assert!(!strip.is_boundary_edge(strip.edge(h)));
    assert_eq!(
        strip.collapse_edge(h),
        Err(HalfedgeEditError::LinkCondition {
            edge: strip.edge(h)
        })
    );
}

#[test]
fn collapses_that_would_degenerate_are_refused() {
    for mesh in [tetrahedron(), lone_triangle(), pillow()] {
        let mut halfedge = build(&mesh);
        let snapshot = halfedge.clone();
        for h in snapshot.halfedges() {
            assert_eq!(
                halfedge.collapse_edge(h),
                Err(HalfedgeEditError::WouldDegenerate {
                    edge: snapshot.edge(h)
                })
            );
        }
        assert_eq!(halfedge, snapshot);
    }
}

#[test]
fn collapsing_next_to_a_polygon_is_refused_by_name() {
    let positions = vec![
        p(0.0, 0.0, 0.0),
        p(1.0, 0.0, 0.0),
        p(1.0, 1.0, 0.0),
        p(0.0, 1.0, 0.0),
        p(0.5, -1.0, 0.0),
    ];
    let faces: Vec<Vec<u32>> = vec![vec![0, 1, 2, 3], vec![0, 4, 1]];
    let mut mesh = HalfedgeMesh::from_faces(positions, &faces).unwrap();
    let h = mesh
        .find_halfedge(VertexId::new(0), VertexId::new(4))
        .unwrap();
    assert_eq!(
        mesh.collapse_edge(h),
        Err(HalfedgeEditError::NotATriangle {
            face: FaceId::new(0),
            degree: 4,
        })
    );
}

// ---- hole filling and compaction ----------------------------------------------

#[test]
fn filling_a_hole_closes_the_surface() {
    let mut mesh = build(&grid(2, 2));
    let hole = mesh.boundary_loops()[0];
    assert_eq!(
        mesh.fill_hole(mesh.opposite(interior_halfedge(&mesh))),
        Err(HalfedgeEditError::NotBoundary {
            halfedge: mesh.opposite(interior_halfedge(&mesh)),
        })
    );
    let face = mesh.fill_hole(hole).unwrap();
    assert_sound(&mesh);
    assert!(mesh.boundary_loops().is_empty());
    assert_eq!(mesh.euler_characteristic(), 2);
    assert_eq!(mesh.face_degree(face), 8);
    mesh.split_face(face, p(1.0, 1.0, -1.0)).unwrap();
    assert_sound(&mesh);
    let closed = mesh.to_tri_mesh().unwrap();
    assert_eq!(build(&closed).euler_characteristic(), 2);
}

fn interior_halfedge(mesh: &HalfedgeMesh) -> HalfedgeId {
    mesh.edge_halfedge(interior_edge(mesh), 0)
}

#[test]
fn compaction_renumbers_densely_and_keeps_the_mesh() {
    let mut mesh = build(&octahedron());
    let h = mesh
        .find_halfedge(VertexId::new(4), VertexId::new(0))
        .unwrap();
    mesh.collapse_edge(h).unwrap();
    let before = mesh.to_tri_mesh().unwrap();
    let remap = mesh.compact();
    assert_sound(&mesh);
    assert_eq!(remap.vertices.len(), 6);
    assert_eq!(remap.vertices[4], None);
    assert_eq!(remap.vertices[5], Some(VertexId::new(4)));
    assert_eq!(remap.edges.iter().filter(|e| e.is_none()).count(), 3);
    assert_eq!(remap.faces.iter().filter(|f| f.is_none()).count(), 2);
    assert_eq!(mesh.positions().len(), 5);
    assert_eq!(mesh.to_tri_mesh().unwrap(), before);
    assert!(mesh.vertices().map(VertexId::index).eq(0..5));
}

// ---- random edit sequences ------------------------------------------------------

#[derive(Default, Debug)]
struct Tally {
    accepted: [usize; 4],
    refused: Vec<HalfedgeEditError>,
}

/// Replay `steps` random edits; after every one, the mesh is sound, a
/// refused edit changed nothing, and the Euler characteristic and the
/// number of holes are those of the start.
fn replay(mut mesh: HalfedgeMesh, seed: u64, steps: usize) -> (HalfedgeMesh, Tally) {
    let mut rng = Rng(seed);
    let chi = mesh.euler_characteristic();
    let holes = mesh.boundary_loops().len();
    let mut tally = Tally::default();
    for _ in 0..steps {
        let before = mesh.clone();
        let edges: Vec<EdgeId> = mesh.edges().collect();
        let e = edges[rng.below(edges.len())];
        let [a, b] = mesh.edge_vertices(e);
        let mid = (mesh.position(a) + mesh.position(b)) * 0.5;
        // Collapses dominate so the mesh neither explodes nor vanishes.
        let roll = rng.below(10);
        let op = match roll {
            0..=2 => 0,
            3..=4 => 1,
            5 => 2,
            _ => 3,
        };
        let result = match op {
            0 => mesh.flip_edge(e),
            1 => mesh.split_edge(e, mid).map(|_| ()),
            2 => {
                let faces: Vec<FaceId> = mesh.faces().collect();
                let f = faces[rng.below(faces.len())];
                let centre = mesh
                    .face_vertices(f)
                    .map(|v| mesh.position(v))
                    .sum::<Point3>()
                    / 3.0;
                mesh.split_face(f, centre).map(|_| ())
            }
            _ => {
                let h = mesh.edge_halfedge(e, rng.below(2) as u32);
                mesh.collapse_edge(h).map(|_| ())
            }
        };
        match result {
            Ok(()) => tally.accepted[op] += 1,
            Err(error) => {
                assert_eq!(mesh, before, "refused {error:?} changed the mesh");
                tally.refused.push(error);
            }
        }
        assert_sound(&mesh);
        assert_eq!(mesh.euler_characteristic(), chi);
        assert_eq!(mesh.boundary_loops().len(), holes);
    }
    (mesh, tally)
}

#[test]
fn random_edits_on_a_torus_keep_every_invariant() {
    let (mut mesh, tally) = replay(build(&torus(6, 5)), 0x5EED_0140, 400);
    assert!(tally.accepted.iter().all(|&n| n > 0), "{tally:?}");
    assert!(
        tally
            .refused
            .iter()
            .any(|e| matches!(e, HalfedgeEditError::LinkCondition { .. })),
        "the link condition was exercised: {tally:?}"
    );
    mesh.compact();
    assert_sound(&mesh);
    assert_eq!(mesh.euler_characteristic(), 0);
}

#[test]
fn random_edits_on_an_open_grid_keep_every_invariant() {
    let (mesh, tally) = replay(build(&grid(5, 4)), 0xB0_0D_A7_E5, 400);
    assert!(tally.accepted.iter().all(|&n| n > 0), "{tally:?}");
    assert!(
        tally
            .refused
            .iter()
            .any(|e| matches!(e, HalfedgeEditError::BoundaryEdge { .. })),
        "{tally:?}"
    );
    assert_eq!(mesh.euler_characteristic(), 1);
}

#[test]
fn random_edits_on_a_sphere_keep_every_invariant() {
    let (mesh, tally) = replay(build(&octahedron()), 7, 300);
    assert!(tally.accepted.iter().all(|&n| n > 0), "{tally:?}");
    assert_eq!(mesh.euler_characteristic(), 2);
}

#[test]
fn random_edits_on_a_disc_with_holes_keep_every_invariant() {
    // A grid with two cells removed: three boundary loops.
    let mut mesh = grid(6, 5);
    let cell = |x: usize, y: usize| (y * 6 + x) * 6;
    for start in [cell(4, 3), cell(1, 1)] {
        mesh.indices.drain(start..start + 6);
    }
    let halfedge = build(&mesh);
    assert_eq!(halfedge.boundary_loops().len(), 3);
    let (mesh, tally) = replay(halfedge, 0xD15C, 400);
    assert!(tally.accepted.iter().all(|&n| n > 0), "{tally:?}");
    assert_eq!(mesh.euler_characteristic(), -1);
}
