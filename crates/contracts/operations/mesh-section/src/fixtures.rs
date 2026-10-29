//! Fixtures of the region conformance suite (ADR 0091).
//!
//! A solid case is a closed, consistently wound, outward oriented solid; a
//! sheet case is an open, consistently wound, non-degenerate surface whose
//! boundary edges were counted by hand (the tests below audit both). Every
//! expected value is derived by hand in the comment of its case. Frame `h(z)` has its origin at `(0, 0, z)` and
//! the model axes as `x`, `y`, `z`, so plane-local `(x, y)` equals model
//! `(x, y)`. A ring starts at its lexicographically smallest point; an outer
//! ring runs counter-clockwise and a hole clockwise.

use axiolid_core::{Frame3, Point2, Point3, Vec3};
use axiolid_mesh::TriMesh;

use crate::cut::CutRegion;
use crate::request::{
    ChainEnd, InPlaneRegion, InPlaneSide, SectionChain, SectionRole,
    REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES, REFUSED_OVERLAPPING_IN_PLANE_FACES,
    REFUSED_SADDLE_VERTEX,
};

/// What a fixture must produce.
pub(crate) enum Expected {
    /// This output, with these source counts.
    Output(Output),
    /// Either the refusal or this output: a provider may merge overlapping
    /// in-plane faces instead of refusing them.
    OutputOrRefused(&'static str, Output),
    /// Only the refusal.
    Refused(&'static str),
}

/// Expected chains, regions and source counts.
pub(crate) struct Output {
    pub chains: Vec<SectionChain>,
    pub regions: Vec<InPlaneRegion>,
    pub in_plane_triangles: usize,
    pub touching_edges: usize,
}

/// What a cut-region request must produce.
pub(crate) enum Cut {
    /// These cut regions; the output otherwise equals the plain request's,
    /// which may be a refusal the cut request must repeat.
    Regions(Vec<CutRegion>),
    /// Only this refusal.
    Refused(&'static str),
    /// `GeomError::InvalidInput`: a sheet has no cut.
    Invalid,
}

/// One named fixture.
pub(crate) struct Case {
    pub name: &'static str,
    pub mesh: TriMesh,
    pub frame: Frame3,
    pub role: SectionRole,
    /// Mesh edges used by one triangle only, counted by hand; zero for a solid.
    pub boundary_edges: usize,
    pub expected: Expected,
    /// What the same request must return with cut regions.
    pub cut: Cut,
}

fn h(z: f64) -> Frame3 {
    Frame3 {
        origin: Point3::new(0.0, 0.0, z),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

fn ring(points: &[(f64, f64)]) -> Vec<Point2> {
    points.iter().map(|&(x, y)| Point2::new(x, y)).collect()
}

fn output(
    chains: Vec<SectionChain>,
    regions: Vec<InPlaneRegion>,
    in_plane_triangles: usize,
    touching_edges: usize,
) -> Output {
    Output {
        chains,
        regions,
        in_plane_triangles,
        touching_edges,
    }
}

/// Mesh under construction; a vertex is shared by exact position.
#[derive(Default)]
struct Builder {
    positions: Vec<Point3>,
    indices: Vec<u32>,
}

impl Builder {
    fn vertex(&mut self, p: [f64; 3]) -> u32 {
        let point = Point3::new(p[0], p[1], p[2]);
        let index = self
            .positions
            .iter()
            .position(|&q| q == point)
            .unwrap_or_else(|| {
                self.positions.push(point);
                self.positions.len() - 1
            });
        u32::try_from(index).expect("fixture index fits u32")
    }

    fn tri(&mut self, p: [[f64; 3]; 3]) {
        for corner in p {
            let index = self.vertex(corner);
            self.indices.push(index);
        }
    }

    /// Two triangles `(a, b, c)` and `(a, c, d)`.
    fn quad(&mut self, p: [[f64; 3]; 4]) {
        self.tri([p[0], p[1], p[2]]);
        self.tri([p[0], p[2], p[3]]);
    }

    fn finish(self) -> TriMesh {
        TriMesh::new(self.positions, self.indices)
    }
}

type Span = (f64, f64);

// Rectangular faces with the outward normal named by the function. The point
// order is counter-clockwise seen from outside.
fn y_minus(b: &mut Builder, y: f64, (x0, x1): Span, (z0, z1): Span) {
    b.quad([[x0, y, z0], [x1, y, z0], [x1, y, z1], [x0, y, z1]]);
}
fn y_plus(b: &mut Builder, y: f64, (x0, x1): Span, (z0, z1): Span) {
    b.quad([[x0, y, z0], [x0, y, z1], [x1, y, z1], [x1, y, z0]]);
}
fn z_minus(b: &mut Builder, z: f64, (x0, x1): Span, (y0, y1): Span) {
    b.quad([[x0, y0, z], [x0, y1, z], [x1, y1, z], [x1, y0, z]]);
}
fn z_plus(b: &mut Builder, z: f64, (x0, x1): Span, (y0, y1): Span) {
    b.quad([[x0, y0, z], [x1, y0, z], [x1, y1, z], [x0, y1, z]]);
}
fn x_minus(b: &mut Builder, x: f64, (y0, y1): Span, (z0, z1): Span) {
    b.quad([[x, y0, z0], [x, y0, z1], [x, y1, z1], [x, y1, z0]]);
}
fn x_plus(b: &mut Builder, x: f64, (y0, y1): Span, (z0, z1): Span) {
    b.quad([[x, y0, z0], [x, y1, z0], [x, y1, z1], [x, y0, z1]]);
}

/// Six faces, twelve triangles.
fn add_box(b: &mut Builder, x: Span, y: Span, z: Span) {
    y_minus(b, y.0, x, z);
    y_plus(b, y.1, x, z);
    z_minus(b, z.0, x, y);
    z_plus(b, z.1, x, y);
    x_minus(b, x.0, y, z);
    x_plus(b, x.1, y, z);
}

fn boxes(parts: &[(Span, Span, Span)]) -> TriMesh {
    let mut b = Builder::default();
    for &(x, y, z) in parts {
        add_box(&mut b, x, y, z);
    }
    b.finish()
}

/// Square annulus prism: outer `[o0,o1]^2`, inner `[i0,i1]^2`, `z` in `[z0,z1]`.
///
/// 16 vertices and 32 triangles (8 per annulus face, 8 per wall ring): with
/// 48 edges `V - E + F = 0`, a torus.
fn add_annulus_prism(b: &mut Builder, (o0, o1): Span, (i0, i1): Span, (z0, z1): Span) {
    let outer = [(o0, o0), (o1, o0), (o1, o1), (o0, o1)];
    let inner = [(i0, i0), (i1, i0), (i1, i1), (i0, i1)];
    let at = |(x, y): (f64, f64), z: f64| [x, y, z];
    for k in 0..4 {
        let n = (k + 1) % 4;
        // Top annulus, normal +z: one trapezoid per side, split on o_k - i_n.
        b.tri([at(outer[k], z1), at(outer[n], z1), at(inner[n], z1)]);
        b.tri([at(outer[k], z1), at(inner[n], z1), at(inner[k], z1)]);
        // Bottom annulus, normal -z: the same triangles reversed.
        b.tri([at(outer[k], z0), at(inner[n], z0), at(outer[n], z0)]);
        b.tri([at(outer[k], z0), at(inner[k], z0), at(inner[n], z0)]);
        // Outer wall faces away from the axis.
        b.tri([at(outer[k], z0), at(outer[n], z0), at(outer[n], z1)]);
        b.tri([at(outer[k], z0), at(outer[n], z1), at(outer[k], z1)]);
        // Inner wall faces into the hole.
        b.tri([at(inner[k], z0), at(inner[n], z1), at(inner[n], z0)]);
        b.tri([at(inner[k], z0), at(inner[k], z1), at(inner[n], z1)]);
    }
}

/// Square annulus prism: outer `[0,3]^2`, inner `[1,2]^2`, `z` in `[0,1]`.
fn frame_prism() -> TriMesh {
    let mut b = Builder::default();
    add_annulus_prism(&mut b, (0.0, 3.0), (1.0, 2.0), (0.0, 1.0));
    b.finish()
}

/// Two nested tubes: `[0,7]^2` minus `[1,6]^2` and `[2,5]^2` minus `[3,4]^2`,
/// `z` in `[0,1]`. The second tube stands in the hole of the first without
/// touching it.
fn nested_tubes() -> TriMesh {
    let mut b = Builder::default();
    add_annulus_prism(&mut b, (0.0, 7.0), (1.0, 6.0), (0.0, 1.0));
    add_annulus_prism(&mut b, (2.0, 5.0), (3.0, 4.0), (0.0, 1.0));
    b.finish()
}

/// A prism over a profile in `(x, z)`, extruded over `y`. `profile` is
/// counter-clockwise in `(x, z)` and `triangles` triangulates it. The front
/// (`y` smallest) copy keeps the profile's winding, which faces `-y`; the
/// back copy is reversed; each profile edge `p -> q` becomes the outward side
/// quad `p@y0, p@y1, q@y1, q@y0`. With `n` profile vertices and `t`
/// triangles: `2n` vertices, `2t + 2n` triangles, `3t + 3n` edges.
fn extruded_profile(profile: &[(f64, f64)], triangles: &[[usize; 3]], (y0, y1): Span) -> TriMesh {
    let at = |index: usize, y: f64| [profile[index].0, y, profile[index].1];
    let mut b = Builder::default();
    for &[p, q, r] in triangles {
        b.tri([at(p, y0), at(q, y0), at(r, y0)]);
        b.tri([at(p, y1), at(r, y1), at(q, y1)]);
    }
    for from in 0..profile.len() {
        let to = (from + 1) % profile.len();
        b.quad([at(from, y0), at(from, y1), at(to, y1), at(to, y0)]);
    }
    b.finish()
}

/// A block `[0,3]^2 x [0,2]` with a pocket `[1,2]^2` from `z = 2` down to
/// `z = 1`: 16 vertices (8 block corners, 4 pocket rim, 4 pocket floor) and
/// 28 triangles (8 block walls, 2 bottom, 8 top annulus, 8 pocket walls, 2
/// floor); `V - E + F = 16 - 42 + 28 = 2`.
fn pocket_block() -> TriMesh {
    let (all, deep) = ((0.0, 3.0), (0.0, 2.0));
    let (pocket, wall) = ((1.0, 2.0), (1.0, 2.0));
    let mut b = Builder::default();
    y_minus(&mut b, 0.0, all, deep);
    y_plus(&mut b, 3.0, all, deep);
    x_minus(&mut b, 0.0, all, deep);
    x_plus(&mut b, 3.0, all, deep);
    z_minus(&mut b, 0.0, all, all);
    // Top annulus at z = 2, normal +z, around the pocket opening.
    let outer = [(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0)];
    let inner = [(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)];
    let at = |(x, y): (f64, f64)| [x, y, 2.0];
    for k in 0..4 {
        let n = (k + 1) % 4;
        b.tri([at(outer[k]), at(outer[n]), at(inner[n])]);
        b.tri([at(outer[k]), at(inner[n]), at(inner[k])]);
    }
    // Pocket walls face into the pocket; the floor faces up.
    x_plus(&mut b, 1.0, pocket, wall);
    x_minus(&mut b, 2.0, pocket, wall);
    y_plus(&mut b, 1.0, pocket, wall);
    y_minus(&mut b, 2.0, pocket, wall);
    z_plus(&mut b, 1.0, pocket, pocket);
    b.finish()
}

/// The two chains around the piers of the wall fixtures.
fn pier_chains() -> Vec<SectionChain> {
    let open = |points: &[(f64, f64)]| {
        SectionChain::open(
            ring(points),
            ChainEnd::InPlaneRegion,
            ChainEnd::InPlaneRegion,
        )
    };
    vec![
        open(&[(1.0, 0.0), (0.0, 0.0), (0.0, 0.25), (1.0, 0.25)]),
        open(&[(2.0, 0.0), (3.0, 0.0), (3.0, 0.25), (2.0, 0.25)]),
    ]
}

/// The two piers of the wall fixtures as cut regions.
fn pier_cuts() -> Vec<CutRegion> {
    let pier = |x0: f64, x1: f64| ring(&[(x0, 0.0), (x1, 0.0), (x1, 0.25), (x0, 0.25)]);
    vec![
        CutRegion::new(pier(0.0, 1.0), vec![]),
        CutRegion::new(pier(2.0, 3.0), vec![]),
    ]
}

/// A parapet `[0,3] x [0,0.1] x [0,1]` with a post `[1,1.1] x [0,0.1]`
/// rising to `z = 1.2`.
///
/// Front and back faces are split into the rectangles `[0,1]`, `[1,1.1]`,
/// `[1.1,3]` times `[0,1]` and `[1,1.1]` times `[1,1.2]`, so no T-junction
/// occurs. Vertices: `x` in {0, 1, 1.1, 3} and `z` in {0, 1} on both `y`
/// planes (16), plus the post top (1, 1.2) and (1.1, 1.2) on both (4): 20.
/// Triangles: front 4 rectangles and back 4 (16), bottom 3 rectangles (6),
/// top strips at `z = 1` left and right of the post (4), post top (2), post
/// sides (4), end caps (4): 36. `V - E + F = 20 - 54 + 36 = 2`.
fn parapet_with_post() -> TriMesh {
    let y = (0.0, 0.1);
    let mut b = Builder::default();
    for x in [(0.0, 1.0), (1.0, 1.1), (1.1, 3.0)] {
        y_minus(&mut b, y.0, x, (0.0, 1.0));
        y_plus(&mut b, y.1, x, (0.0, 1.0));
        z_minus(&mut b, 0.0, x, y);
    }
    y_minus(&mut b, y.0, (1.0, 1.1), (1.0, 1.2));
    y_plus(&mut b, y.1, (1.0, 1.1), (1.0, 1.2));
    z_plus(&mut b, 1.0, (0.0, 1.0), y);
    z_plus(&mut b, 1.0, (1.1, 3.0), y);
    z_plus(&mut b, 1.2, (1.0, 1.1), y);
    x_minus(&mut b, 1.0, y, (1.0, 1.2));
    x_plus(&mut b, 1.1, y, (1.0, 1.2));
    x_minus(&mut b, 0.0, y, (0.0, 1.0));
    x_plus(&mut b, 3.0, y, (0.0, 1.0));
    b.finish()
}

/// Slab `[0,3]^2 x [0,1]` with a triangular post over A = (0,0), F = (2,1),
/// E = (1,2) for `z` in `[1,2]`; A is also a corner of the slab.
///
/// The slab top is the square minus the post footprint, in five triangles
/// counter-clockwise from above: (A,B,F), (B,C,F), (C,E,F), (C,D,E) and
/// (D,A,E) with B = (3,0), C = (3,3), D = (0,3). Their areas are 1.5 each,
/// 7.5 in all, which is 9 minus the footprint's 1.5. The post is one closed
/// solid with the slab: its walls start on the footprint edges, which the
/// slab top borders in the opposite direction.
fn slab_with_corner_post() -> TriMesh {
    let (a, b, c, d) = ([0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0]);
    let (e, f) = ([1.0, 2.0], [2.0, 1.0]);
    let at = |p: [f64; 2], z: f64| [p[0], p[1], z];
    let mut m = Builder::default();
    y_minus(&mut m, 0.0, (0.0, 3.0), (0.0, 1.0));
    y_plus(&mut m, 3.0, (0.0, 3.0), (0.0, 1.0));
    x_minus(&mut m, 0.0, (0.0, 3.0), (0.0, 1.0));
    x_plus(&mut m, 3.0, (0.0, 3.0), (0.0, 1.0));
    z_minus(&mut m, 0.0, (0.0, 3.0), (0.0, 3.0));
    for [p, q, r] in [[a, b, f], [b, c, f], [c, e, f], [c, d, e], [d, a, e]] {
        m.tri([at(p, 1.0), at(q, 1.0), at(r, 1.0)]);
    }
    // The footprint A, F, E is counter-clockwise; each wall is outward.
    for (p, q) in [(a, f), (f, e), (e, a)] {
        m.quad([at(p, 1.0), at(q, 1.0), at(q, 2.0), at(p, 2.0)]);
    }
    m.tri([at(a, 2.0), at(f, 2.0), at(e, 2.0)]);
    m.finish()
}

/// Triangular prism whose lower edge is `v0 - v3` and whose top is the
/// quadrilateral `v1 v2 v5 v4` at `z = 1`.
fn prism() -> TriMesh {
    let p = |x, y, z| Point3::new(x, y, z);
    TriMesh::new(
        vec![
            p(0.0, 0.0, 0.0),
            p(0.0, 1.0, 1.0),
            p(0.0, -1.0, 1.0),
            p(2.0, 0.0, 0.0),
            p(2.0, 1.0, 1.0),
            p(2.0, -1.0, 1.0),
        ],
        vec![
            0, 2, 1, 3, 4, 5, 0, 1, 4, 0, 4, 3, 0, 3, 5, 0, 5, 2, 1, 2, 5, 1, 5, 4,
        ],
    )
}

/// Octahedron with vertices `(+-1,0,0)`, `(0,+-1,0)`, `(0,0,+-1)`.
fn octahedron() -> TriMesh {
    let mut b = Builder::default();
    let ring = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    ];
    for k in 0..4 {
        let n = (k + 1) % 4;
        b.tri([ring[k], ring[n], [0.0, 0.0, 1.0]]);
        b.tri([ring[n], ring[k], [0.0, 0.0, -1.0]]);
    }
    b.finish()
}

/// Frame on the plane `x = x0` with local `(x, y)` equal to model `(y, z)`;
/// its normal `Y x Z` is model `+x`.
fn plane_x(x0: f64) -> Frame3 {
    Frame3 {
        origin: Point3::new(x0, 0.0, 0.0),
        x: Vec3::Y,
        y: Vec3::Z,
        z: Vec3::X,
    }
}

/// Terrain over a 3 x 3 grid: vertex `v(i, j) = i + 3 j` lies at
/// `(i, j, h(i, j))` with the heights `h` below, and cell `(i, j)` is the two
/// triangles `(v(i,j), v(i+1,j), v(i+1,j+1))` and `(v(i,j), v(i+1,j+1),
/// v(i,j+1))`, both facing `+z`. Eight triangles; the boundary is the
/// perimeter of the grid: 4 sides of 2 edges, 8 boundary edges.
fn terrain() -> TriMesh {
    // Rows j = 0, 1, 2; columns i = 0, 1, 2.
    let heights = [[0.0, 2.0, 1.0], [2.0, 0.0, 1.0], [0.0, 3.0, 1.0]];
    let mut positions = Vec::new();
    for (j, row) in heights.iter().enumerate() {
        for (i, &z) in row.iter().enumerate() {
            positions.push(Point3::new(i as f64, j as f64, z));
        }
    }
    let v = |i: u32, j: u32| i + 3 * j;
    let mut indices = Vec::new();
    for j in 0..2 {
        for i in 0..2 {
            indices.extend([v(i, j), v(i + 1, j), v(i + 1, j + 1)]);
            indices.extend([v(i, j), v(i + 1, j + 1), v(i, j + 1)]);
        }
    }
    TriMesh::new(positions, indices)
}

/// The square `(0,0) (2,0) (2,2) (0,2)` at `z = 0` as two triangles; the
/// flipped one winds clockwise seen from `+z`. Four boundary edges.
fn square_sheet(flipped: bool) -> TriMesh {
    let p = |x, y| Point3::new(x, y, 0.0);
    let positions = vec![p(0.0, 0.0), p(2.0, 0.0), p(2.0, 2.0), p(0.0, 2.0)];
    let indices = if flipped {
        vec![0, 2, 1, 0, 3, 2]
    } else {
        vec![0, 1, 2, 0, 2, 3]
    };
    TriMesh::new(positions, indices)
}

/// Square annulus sheet at `z = 0`: outer `v0..v3 = (0,0) (3,0) (3,3) (0,3)`,
/// inner `v4..v7 = (1,1) (2,1) (2,2) (1,2)`, eight triangles facing `+z`
/// (each has a positive cross product). Eight boundary edges: four outer and
/// four inner.
fn annulus_sheet() -> TriMesh {
    let p = |x, y| Point3::new(x, y, 0.0);
    TriMesh::new(
        vec![
            p(0.0, 0.0),
            p(3.0, 0.0),
            p(3.0, 3.0),
            p(0.0, 3.0),
            p(1.0, 1.0),
            p(2.0, 1.0),
            p(2.0, 2.0),
            p(1.0, 2.0),
        ],
        vec![
            0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0, 4, 3, 4, 7,
        ],
    )
}

/// A four-triangle fan around `(0,0,0)` over the ring `(1,0,1) (0,1,-1)
/// (-1,0,1) (0,-1,-1)`: a saddle. Four boundary edges, the ring.
fn saddle_sheet() -> TriMesh {
    let p = |x, y, z| Point3::new(x, y, z);
    TriMesh::new(
        vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 1.0),
            p(0.0, 1.0, -1.0),
            p(-1.0, 0.0, 1.0),
            p(0.0, -1.0, -1.0),
        ],
        vec![0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 1],
    )
}

/// The unit square ring starting at the origin.
fn unit_square() -> Vec<Point2> {
    ring(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
}

/// All fixtures; the first decides whether a provider claims the method.
pub(crate) fn cases() -> Vec<Case> {
    use InPlaneSide::{Above, Below, Both};
    let region = InPlaneRegion::new;
    let square =
        |x0: f64, y0: f64, x1: f64, y1: f64| ring(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]);
    vec![
        // Cube top at z = 1: its two triangles face +z, so the material is
        // below. The diagonal is interior; the four corners remain. The
        // side triangles reach the plane in an edge or a vertex, neither is a
        // crossing, and each plane edge borders the top face: no touching.
        Case {
            name: "cube through its top face",
            role: SectionRole::Solid,
            boundary_edges: 0,
            cut: Cut::Regions(vec![]),
            mesh: crate::conformance::unit_cube(),
            frame: h(1.0),
            expected: Expected::Output(output(
                vec![],
                vec![region(Below, unit_square(), vec![])],
                2,
                0,
            )),
        },
        // Cube bottom at z = 0: the faces face -z, the material is above.
        Case {
            name: "cube through its bottom face",
            role: SectionRole::Solid,
            boundary_edges: 0,
            cut: Cut::Regions(vec![]),
            mesh: crate::conformance::unit_cube(),
            frame: h(0.0),
            expected: Expected::Output(output(
                vec![],
                vec![region(Above, unit_square(), vec![])],
                2,
                0,
            )),
        },
        // Annulus prism top at z = 1: 8 triangles facing +z form one
        // component with an outer boundary [0,3]^2 and an inner [1,2]^2.
        // The hole runs clockwise from its smallest point (1,1).
        Case {
            name: "frame top with a hole",
            role: SectionRole::Solid,
            boundary_edges: 0,
            cut: Cut::Regions(vec![]),
            mesh: frame_prism(),
            frame: h(1.0),
            expected: Expected::Output(output(
                vec![],
                vec![region(
                    Below,
                    square(0.0, 0.0, 3.0, 3.0),
                    vec![ring(&[(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)])],
                )],
                8,
                0,
            )),
        },
        // Parapet top at z = 1: the strips [0,1] and [1.1,3] are in-plane,
        // facing +z (two triangles each: 4). The post is a separate
        // crossing: its front edge (1,0,1)-(1.1,0,1) lies in the plane
        // between a triangle whose third vertex is at z = 0 and one whose
        // third vertex is at z = 1.2, so it is a cut segment; likewise at
        // y = 0.1. Its endpoints (1,0), (1.1,0) are vertices of the strips,
        // so each chain is open and ends on a region. The post walls at
        // x = 1 and 1.1 border the strips: region edges, not cuts.
        Case {
            name: "parapet with a post",
            role: SectionRole::Solid,
            boundary_edges: 0,
            cut: Cut::Regions(vec![CutRegion::new(square(1.0, 0.0, 1.1, 0.1), vec![])]),
            mesh: parapet_with_post(),
            frame: h(1.0),
            expected: Expected::Output(output(
                vec![
                    SectionChain::open(
                        ring(&[(1.0, 0.0), (1.1, 0.0)]),
                        ChainEnd::InPlaneRegion,
                        ChainEnd::InPlaneRegion,
                    ),
                    SectionChain::open(
                        ring(&[(1.0, 0.1), (1.1, 0.1)]),
                        ChainEnd::InPlaneRegion,
                        ChainEnd::InPlaneRegion,
                    ),
                ],
                vec![
                    region(Below, square(0.0, 0.0, 1.0, 0.1), vec![]),
                    region(Below, square(1.1, 0.0, 3.0, 0.1), vec![]),
                ],
                4,
                0,
            )),
        },
        // The same parapet seen with the frame x = -X, y = Y, z = -Z (a half
        // turn about Y) at origin (0,0,1): local (x, y) is model (-x, y). The
        // strips now face against the frame's +z: Above. Their rings stay
        // counter-clockwise from the smallest point: [-1,0] x [0,0.1] and
        // [-3,-1.1] x [0,0.1], the latter sorting first. The post chains run
        // from x = -1 to x = -1.1 in mesh order and are reversed so the
        // smaller end comes first; each end is a vertex of a strip.
        Case {
            name: "parapet with a post seen from behind",
            role: SectionRole::Solid,
            boundary_edges: 0,
            cut: Cut::Regions(vec![CutRegion::new(square(-1.1, 0.0, -1.0, 0.1), vec![])]),
            mesh: parapet_with_post(),
            frame: Frame3 {
                origin: Point3::new(0.0, 0.0, 1.0),
                x: -Vec3::X,
                y: Vec3::Y,
                z: -Vec3::Z,
            },
            expected: Expected::Output(output(
                vec![
                    SectionChain::open(
                        ring(&[(-1.1, 0.0), (-1.0, 0.0)]),
                        ChainEnd::InPlaneRegion,
                        ChainEnd::InPlaneRegion,
                    ),
                    SectionChain::open(
                        ring(&[(-1.1, 0.1), (-1.0, 0.1)]),
                        ChainEnd::InPlaneRegion,
                        ChainEnd::InPlaneRegion,
                    ),
                ],
                vec![
                    region(Above, square(-3.0, 0.0, -1.1, 0.1), vec![]),
                    region(Above, square(-1.0, 0.0, 0.0, 0.1), vec![]),
                ],
                4,
                0,
            )),
        },
        // The slab top at z = 1 is the five triangles above, all facing +z:
        // one component. Its boundary is the square A B C D and the footprint
        // A E F, which meet at A, so the boundary is one weakly simple ring
        // visiting A twice: A B C D A E F. A is the smallest point in every
        // ascending order, so the outer/hole decision falls to a descending
        // order. The two rotations starting at A are [A B C D A E F] and
        // [A E F A B C D]; (1,2) sorts before (3,0), so the second is
        // canonical. No chains: the post walls border the top.
        Case {
            name: "slab with a post touching its corner",
            role: SectionRole::Solid,
            boundary_edges: 0,
            cut: Cut::Regions(vec![CutRegion::new(
                ring(&[(0.0, 0.0), (2.0, 1.0), (1.0, 2.0)]),
                vec![],
            )]),
            mesh: slab_with_corner_post(),
            frame: h(1.0),
            expected: Expected::Output(output(
                vec![],
                vec![region(
                    Below,
                    ring(&[
                        (0.0, 0.0),
                        (1.0, 2.0),
                        (2.0, 1.0),
                        (0.0, 0.0),
                        (3.0, 0.0),
                        (3.0, 3.0),
                        (0.0, 3.0),
                    ]),
                    vec![],
                )],
                5,
                0,
            )),
        },
        // The plane z = 0 meets the prism only along its lower edge v0-v3.
        // Both triangles at that edge rise to z = 1 on the same side, so it
        // touches: nothing is emitted, one touching edge is counted.
        Case {
            name: "prism touching along an edge",
            role: SectionRole::Solid,
            boundary_edges: 0,
            cut: Cut::Regions(vec![]),
            mesh: prism(),
            frame: h(0.0),
            expected: Expected::Output(output(vec![], vec![], 0, 1)),
        },
        // The prism top at z = 1 faces +z: the quadrilateral v2 v5 v4 v1,
        // (0,-1) (2,-1) (2,1) (0,1), counter-clockwise from (0,-1).
        Case {
            name: "prism top",
            role: SectionRole::Solid,
            boundary_edges: 0,
            cut: Cut::Regions(vec![]),
            mesh: prism(),
            frame: h(1.0),
            expected: Expected::Output(output(
                vec![],
                vec![region(Below, square(0.0, -1.0, 2.0, 1.0), vec![])],
                2,
                0,
            )),
        },
        // The octahedron's apex (0,0,1) is its only point at z = 1: a point
        // touch, neither a segment nor an edge.
        Case {
            name: "point touch",
            role: SectionRole::Solid,
            boundary_edges: 0,
            cut: Cut::Regions(vec![]),
            mesh: octahedron(),
            frame: h(1.0),
            expected: Expected::Output(output(vec![], vec![], 0, 0)),
        },
        // Slab [0,4]^2 x [0,1] with a box [1,2]^2 x [1,2] standing on it,
        // two shells in one mesh. At z = 1 the slab top (+z, Below) covers
        // the box bottom (-z, Above): faces of opposite sides overlap. A
        // provider may refuse, or report the slab top minus the footprint as
        // a Below region with a hole and the footprint as Both.
        Case {
            name: "box on a slab",
            role: SectionRole::Solid,
            boundary_edges: 0,
            cut: Cut::Refused(REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES),
            mesh: boxes(&[
                ((0.0, 4.0), (0.0, 4.0), (0.0, 1.0)),
                ((1.0, 2.0), (1.0, 2.0), (1.0, 2.0)),
            ]),
            frame: h(1.0),
            expected: Expected::OutputOrRefused(
                REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES,
                output(
                    vec![],
                    vec![
                        region(
                            Below,
                            square(0.0, 0.0, 4.0, 4.0),
                            vec![ring(&[(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)])],
                        ),
                        region(Both, square(1.0, 1.0, 2.0, 2.0), vec![]),
                    ],
                    4,
                    0,
                ),
            ),
        },
        // Boxes [0,2] x [0,1] x [0,1] and [1,3] x [0,1] x [0,1] overlap in
        // [1,2] x [0,1]. At z = 1 both tops face +z: overlapping faces of one
        // side. A provider may refuse, or merge them into [0,3] x [0,1].
        Case {
            name: "overlapping boxes",
            role: SectionRole::Solid,
            boundary_edges: 0,
            cut: Cut::Regions(vec![]),
            mesh: boxes(&[
                ((0.0, 2.0), (0.0, 1.0), (0.0, 1.0)),
                ((1.0, 3.0), (0.0, 1.0), (0.0, 1.0)),
            ]),
            frame: h(1.0),
            expected: Expected::OutputOrRefused(
                REFUSED_OVERLAPPING_IN_PLANE_FACES,
                output(
                    vec![],
                    vec![region(Below, square(0.0, 0.0, 3.0, 1.0), vec![])],
                    4,
                    0,
                ),
            ),
        },
        // Cut regions (ADR 0092). The L-step is the profile (0,0) (2,0) (2,1)
        // (1,1) (1,2) (0,2) in (x, z), extruded over y in [0,1]: a column
        // x in [0,1] up to z = 2 beside a ledge x in [1,2] topped at z = 1.
        // At z = 1 the ledge top faces +z: a Below region [1,2] x [0,1], two
        // triangles. The column is cut: its front and back faces cross the
        // plane along y = 0 and y = 1 (the mid-face crossing on the diagonal
        // is collinear and goes), and its wall x = 0 crosses at x = 0. One
        // open chain (1,0) (0,0) (0,1) (1,1) runs from the region's corner
        // (1,0) to its corner (1,1). The ledge's jamb x = 1 rises from the
        // region edge (1,0)-(1,1), away from the material below, so the cut
        // lies across that edge: the cut ring is the chain closed along it,
        // [0,1] x [0,1], counter-clockwise from (0,0).
        Case {
            name: "L-step through a ledge top",
            role: SectionRole::Solid,
            boundary_edges: 0,
            mesh: extruded_profile(
                &[
                    (0.0, 0.0),
                    (2.0, 0.0),
                    (2.0, 1.0),
                    (1.0, 1.0),
                    (1.0, 2.0),
                    (0.0, 2.0),
                ],
                &[[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5]],
                (0.0, 1.0),
            ),
            frame: h(1.0),
            expected: Expected::Output(output(
                vec![SectionChain::open(
                    ring(&[(1.0, 0.0), (0.0, 0.0), (0.0, 1.0), (1.0, 1.0)]),
                    ChainEnd::InPlaneRegion,
                    ChainEnd::InPlaneRegion,
                )],
                vec![region(Below, square(1.0, 0.0, 2.0, 1.0), vec![])],
                2,
                0,
            )),
            cut: Cut::Regions(vec![CutRegion::new(unit_square(), vec![])]),
        },
        // A wall of thickness 0.25 (y in [0,0.25]) with a sill: the U profile
        // (0,0) (3,0) (3,2) (2,2) (2,1) (1,1) (1,2) (0,2) in (x, z). At z = 1
        // the sill top between the piers faces +z: a Below region
        // [1,2] x [0,0.25]. Each pier is cut across: the chains (1,0) (0,0)
        // (0,0.25) (1,0.25) and (2,0) (3,0) (3,0.25) (2,0.25). The jambs x = 1
        // and x = 2 rise from the region edges, so the cuts are the two piers
        // [0,1] x [0,0.25] and [2,3] x [0,0.25]; the region lies between them.
        Case {
            name: "wall at a sill",
            role: SectionRole::Solid,
            boundary_edges: 0,
            mesh: extruded_profile(
                &[
                    (0.0, 0.0),
                    (3.0, 0.0),
                    (3.0, 2.0),
                    (2.0, 2.0),
                    (2.0, 1.0),
                    (1.0, 1.0),
                    (1.0, 2.0),
                    (0.0, 2.0),
                ],
                &[
                    [0, 1, 4],
                    [1, 2, 3],
                    [1, 3, 4],
                    [0, 4, 5],
                    [0, 5, 6],
                    [0, 6, 7],
                ],
                (0.0, 0.25),
            ),
            frame: h(1.0),
            expected: Expected::Output(output(
                pier_chains(),
                vec![region(Below, square(1.0, 0.0, 2.0, 0.25), vec![])],
                2,
                0,
            )),
            cut: Cut::Regions(pier_cuts()),
        },
        // The same wall with a door: the profile (0,0) (1,0) (1,1) (2,1)
        // (2,0) (3,0) (3,2) (0,2). At z = 1 the lintel underside between the
        // jambs faces -z: an Above region [1,2] x [0,0.25]. The jambs fall
        // from its edges, so the chains and the cut piers are those of the
        // sill.
        Case {
            name: "wall at a head",
            role: SectionRole::Solid,
            boundary_edges: 0,
            mesh: extruded_profile(
                &[
                    (0.0, 0.0),
                    (1.0, 0.0),
                    (1.0, 1.0),
                    (2.0, 1.0),
                    (2.0, 0.0),
                    (3.0, 0.0),
                    (3.0, 2.0),
                    (0.0, 2.0),
                ],
                &[
                    [0, 1, 2],
                    [0, 2, 7],
                    [2, 3, 7],
                    [3, 6, 7],
                    [3, 4, 5],
                    [3, 5, 6],
                ],
                (0.0, 0.25),
            ),
            frame: h(1.0),
            expected: Expected::Output(output(
                pier_chains(),
                vec![region(Above, square(1.0, 0.0, 2.0, 0.25), vec![])],
                2,
                0,
            )),
            cut: Cut::Regions(pier_cuts()),
        },
        // A stepped slab: the profile (0,0) (3,0) (3,0.5) (5,0.5) (5,1) (2,1)
        // (2,0.5) (0,0.5), y in [0,1], cut at z = 0.5. The lower part's top
        // [0,2] x [0,1] faces +z (Below, two triangles) and the upper part's
        // underside [3,5] x [0,1] faces -z (Above, two). Between them, x in
        // [2,3], the slab is solid across the plane: the cut. Its front and
        // back faces hold the mesh edge (2,0.5)-(3,0.5) in the plane between
        // a triangle below and one above: the chains (2,0) (3,0) and
        // (2,1) (3,1), each from one region to the other. The riser x = 2
        // rises from the Below edge and the riser x = 3 falls from the Above
        // edge, so both bound the cut [2,3] x [0,1].
        Case {
            name: "stepped slab",
            role: SectionRole::Solid,
            boundary_edges: 0,
            mesh: extruded_profile(
                &[
                    (0.0, 0.0),
                    (3.0, 0.0),
                    (3.0, 0.5),
                    (5.0, 0.5),
                    (5.0, 1.0),
                    (2.0, 1.0),
                    (2.0, 0.5),
                    (0.0, 0.5),
                ],
                &[
                    [0, 1, 2],
                    [0, 2, 6],
                    [0, 6, 7],
                    [2, 3, 4],
                    [2, 4, 5],
                    [2, 5, 6],
                ],
                (0.0, 1.0),
            ),
            frame: h(0.5),
            expected: Expected::Output(output(
                vec![
                    SectionChain::open(
                        ring(&[(2.0, 0.0), (3.0, 0.0)]),
                        ChainEnd::InPlaneRegion,
                        ChainEnd::InPlaneRegion,
                    ),
                    SectionChain::open(
                        ring(&[(2.0, 1.0), (3.0, 1.0)]),
                        ChainEnd::InPlaneRegion,
                        ChainEnd::InPlaneRegion,
                    ),
                ],
                vec![
                    region(Below, square(0.0, 0.0, 2.0, 1.0), vec![]),
                    region(Above, square(3.0, 0.0, 5.0, 1.0), vec![]),
                ],
                4,
                0,
            )),
            cut: Cut::Regions(vec![CutRegion::new(square(2.0, 0.0, 3.0, 1.0), vec![])]),
        },
        // A block [0,3]^2 x [0,2] with a pocket [1,2]^2 from z = 2 down to
        // z = 1 (16 vertices, 28 triangles). At z = 1 the pocket floor faces
        // +z: a Below region [1,2]^2, two triangles. The block's walls cross
        // the plane along the closed chain [0,3]^2. The pocket walls rise
        // from the floor edges, so the cut is the square with the floor as a
        // hole, clockwise from (1,1).
        Case {
            name: "pocket floor",
            role: SectionRole::Solid,
            boundary_edges: 0,
            mesh: pocket_block(),
            frame: h(1.0),
            expected: Expected::Output(output(
                vec![SectionChain::closed(square(0.0, 0.0, 3.0, 3.0))],
                vec![region(Below, square(1.0, 1.0, 2.0, 2.0), vec![])],
                2,
                0,
            )),
            cut: Cut::Regions(vec![CutRegion::new(
                square(0.0, 0.0, 3.0, 3.0),
                vec![ring(&[(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)])],
            )]),
        },
        // Two nested tubes, annulus prisms [0,7]^2 minus [1,6]^2 and
        // [2,5]^2 minus [3,4]^2, z in [0,1], cut at z = 0.5. Four closed
        // chains: the outer walls run counter-clockwise and the inner walls
        // clockwise, keeping the material on the left. The cut is two
        // regions, each a ring with a hole; the inner tube lies in the hole
        // of the outer one.
        Case {
            name: "nested tubes",
            role: SectionRole::Solid,
            boundary_edges: 0,
            mesh: nested_tubes(),
            frame: h(0.5),
            expected: Expected::Output(output(
                vec![
                    SectionChain::closed(square(0.0, 0.0, 7.0, 7.0)),
                    SectionChain::closed(ring(&[(1.0, 1.0), (1.0, 6.0), (6.0, 6.0), (6.0, 1.0)])),
                    SectionChain::closed(square(2.0, 2.0, 5.0, 5.0)),
                    SectionChain::closed(ring(&[(3.0, 3.0), (3.0, 4.0), (4.0, 4.0), (4.0, 3.0)])),
                ],
                vec![],
                0,
                0,
            )),
            cut: Cut::Regions(vec![
                CutRegion::new(
                    square(0.0, 0.0, 7.0, 7.0),
                    vec![ring(&[(1.0, 1.0), (1.0, 6.0), (6.0, 6.0), (6.0, 1.0)])],
                ),
                CutRegion::new(
                    square(2.0, 2.0, 5.0, 5.0),
                    vec![ring(&[(3.0, 3.0), (3.0, 4.0), (4.0, 4.0), (4.0, 3.0)])],
                ),
            ]),
        },
        // Sheets (open meshes). Terrain on the plane x = 0.5 (local (x, y) =
        // model (y, z)): only cell column i = 0 is cut. Its edges crossing
        // x = 0.5 are v00-v10 at (y, z) = (0, 1), the diagonal v00-v11 at
        // (0.5, 0), v01-v11 at (1, 1), v01-v12 at (1.5, 2.5) and v12-v02 at
        // (2, 1.5). The four triangles give the segments (0,1)-(0.5,0),
        // (0.5,0)-(1,1), (1,1)-(1.5,2.5), (1.5,2.5)-(2,1.5): one open chain,
        // no node collinear. Its first node lies on the boundary edge
        // v00-v10 (y = 0) and its last on v12-v02 (y = 2): MeshBoundary at
        // both ends. 8 boundary edges.
        Case {
            name: "terrain cut by a vertical plane",
            mesh: terrain(),
            frame: plane_x(0.5),
            role: SectionRole::Sheet,
            cut: Cut::Invalid,
            boundary_edges: 8,
            expected: Expected::Output(output(
                vec![SectionChain::open(
                    ring(&[(0.0, 1.0), (0.5, 0.0), (1.0, 1.0), (1.5, 2.5), (2.0, 1.5)]),
                    ChainEnd::MeshBoundary,
                    ChainEnd::MeshBoundary,
                )],
                vec![],
                0,
                0,
            )),
        },
        // The plane x = 1 holds v10 = (y, z) (0, 2), v11 (1, 0), v12 (2, 3).
        // Edge v10-v11 lies between (v00,v10,v11), whose third vertex is at
        // x = 0, and (v10,v21,v11), whose third vertex is at x = 2: opposite
        // sides, a cut segment. Likewise v11-v12 between (v01,v11,v12) and
        // (v11,v22,v12). The triangles with only one vertex on the plane have
        // the others on one side and add nothing. One open chain (0,2) (1,0)
        // (2,3), not collinear; v10 and v12 lie on the boundary y = 0 and
        // y = 2 and touch no in-plane triangle: MeshBoundary.
        Case {
            name: "terrain cut through vertices",
            mesh: terrain(),
            frame: plane_x(1.0),
            role: SectionRole::Sheet,
            cut: Cut::Invalid,
            boundary_edges: 8,
            expected: Expected::Output(output(
                vec![SectionChain::open(
                    ring(&[(0.0, 2.0), (1.0, 0.0), (2.0, 3.0)]),
                    ChainEnd::MeshBoundary,
                    ChainEnd::MeshBoundary,
                )],
                vec![],
                0,
                0,
            )),
        },
        // Two triangles (0,1,2), (0,2,3) face +z: one Below region, the
        // square counter-clockwise from (0,0); the diagonal is interior. Each
        // square edge has one in-plane triangle: a region edge, no touching.
        Case {
            name: "sheet lying in the plane",
            mesh: square_sheet(false),
            frame: h(0.0),
            role: SectionRole::Sheet,
            cut: Cut::Invalid,
            boundary_edges: 4,
            expected: Expected::Output(output(
                vec![],
                vec![region(
                    Below,
                    ring(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]),
                    vec![],
                )],
                2,
                0,
            )),
        },
        // The same square wound the other way faces -z: the region is Above,
        // its ring still counter-clockwise in the frame.
        Case {
            name: "sheet lying in the plane, flipped",
            mesh: square_sheet(true),
            frame: h(0.0),
            role: SectionRole::Sheet,
            cut: Cut::Invalid,
            boundary_edges: 4,
            expected: Expected::Output(output(
                vec![],
                vec![region(
                    Above,
                    ring(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]),
                    vec![],
                )],
                2,
                0,
            )),
        },
        // Eight triangles facing +z form one component. Its boundary edges
        // are the outer square and the inner square [1,2]^2, whose directed
        // edges run 4 -> 7 -> 6 -> 5 -> 4 = (1,1) (1,2) (2,2) (2,1):
        // clockwise, a hole. No chains; 8 boundary edges.
        Case {
            name: "annulus sheet",
            mesh: annulus_sheet(),
            frame: h(0.0),
            role: SectionRole::Sheet,
            cut: Cut::Invalid,
            boundary_edges: 8,
            expected: Expected::Output(output(
                vec![],
                vec![region(
                    Below,
                    square(0.0, 0.0, 3.0, 3.0),
                    vec![ring(&[(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)])],
                )],
                8,
                0,
            )),
        },
        // The centre lies on z = 0 and the ring alternates +1, -1, +1, -1:
        // every one of the four triangles cuts from the centre to the middle
        // of its ring edge, so four segments meet at the centre.
        Case {
            name: "saddle vertex",
            mesh: saddle_sheet(),
            frame: h(0.0),
            role: SectionRole::Sheet,
            cut: Cut::Invalid,
            boundary_edges: 4,
            expected: Expected::Refused(REFUSED_SADDLE_VERTEX),
        },
    ]
}

#[cfg(test)]
mod tests {
    use axiolid_core::Tolerance;
    use axiolid_mesh::try_audit_mesh;
    use axiolid_mesh_contracts::SolidRequirements;

    use super::*;

    #[test]
    fn every_solid_is_closed_and_oriented_and_every_sheet_a_consistent_surface() {
        let tolerance = Tolerance::new(1e-9, 1e-9).expect("valid tolerance");
        for case in cases() {
            let health = try_audit_mesh(&case.mesh, tolerance).expect("audit");
            assert_eq!(health.boundary_edges, case.boundary_edges, "{}", case.name);
            match case.role {
                SectionRole::Solid => {
                    assert!(health.is_closed_two_manifold(), "{}: {health:?}", case.name);
                    SolidRequirements::Oriented
                        .validate(&case.mesh, case.name)
                        .unwrap_or_else(|error| panic!("{}: {error}", case.name));
                }
                SectionRole::Sheet => {
                    assert!(case.boundary_edges > 0, "{}", case.name);
                    assert_eq!(health.degenerate_triangles, 0, "{}", case.name);
                    assert_eq!(health.non_manifold_edges, 0, "{}", case.name);
                    assert_eq!(health.inconsistent_winding_edges, 0, "{}", case.name);
                }
            }
        }
    }

    #[test]
    fn parapet_with_post_has_the_derived_counts() {
        let mesh = parapet_with_post();
        assert_eq!(mesh.positions.len(), 20);
        assert_eq!(mesh.triangle_count(), 36);
    }

    #[test]
    fn frame_prism_is_a_torus() {
        let mesh = frame_prism();
        assert_eq!(mesh.positions.len(), 16);
        assert_eq!(mesh.triangle_count(), 32);
    }
}
