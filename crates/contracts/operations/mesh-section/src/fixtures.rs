//! Fixtures of the region conformance suite (ADR 0091).
//!
//! Every mesh is a closed, consistently wound, outward oriented solid (the
//! tests below audit that), and every expected value is derived by hand in
//! the comment of its case. Frame `h(z)` has its origin at `(0, 0, z)` and
//! the model axes as `x`, `y`, `z`, so plane-local `(x, y)` equals model
//! `(x, y)`. A ring starts at its lexicographically smallest point; an outer
//! ring runs counter-clockwise and a hole clockwise.

use axiolid_core::{Frame3, Point2, Point3, Vec3};
use axiolid_mesh::TriMesh;

use crate::request::{
    ChainEnd, InPlaneRegion, InPlaneSide, SectionChain, REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES,
    REFUSED_OVERLAPPING_IN_PLANE_FACES,
};

/// What a fixture must produce.
pub(crate) enum Expected {
    /// This output, with these source counts.
    Output(Output),
    /// Either the refusal or this output: a provider may merge overlapping
    /// in-plane faces instead of refusing them.
    OutputOrRefused(&'static str, Output),
}

/// Expected chains, regions and source counts.
pub(crate) struct Output {
    pub chains: Vec<SectionChain>,
    pub regions: Vec<InPlaneRegion>,
    pub in_plane_triangles: usize,
    pub touching_edges: usize,
}

/// One named fixture.
pub(crate) struct Case {
    pub name: &'static str,
    pub mesh: TriMesh,
    pub frame: Frame3,
    pub expected: Expected,
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

/// Square annulus prism: outer `[0,3]^2`, inner `[1,2]^2`, `z` in `[0,1]`.
///
/// 16 vertices and 32 triangles (8 per annulus face, 8 per wall ring): with
/// 48 edges `V - E + F = 0`, a torus.
fn frame_prism() -> TriMesh {
    let outer = [(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0)];
    let inner = [(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)];
    let at = |(x, y): (f64, f64), z: f64| [x, y, z];
    let mut b = Builder::default();
    for k in 0..4 {
        let n = (k + 1) % 4;
        // Top annulus, normal +z: one trapezoid per side, split on o_k - i_n.
        b.tri([at(outer[k], 1.0), at(outer[n], 1.0), at(inner[n], 1.0)]);
        b.tri([at(outer[k], 1.0), at(inner[n], 1.0), at(inner[k], 1.0)]);
        // Bottom annulus, normal -z: the same triangles reversed.
        b.tri([at(outer[k], 0.0), at(inner[n], 0.0), at(outer[n], 0.0)]);
        b.tri([at(outer[k], 0.0), at(inner[k], 0.0), at(inner[n], 0.0)]);
        // Outer wall faces away from the axis.
        b.tri([at(outer[k], 0.0), at(outer[n], 0.0), at(outer[n], 1.0)]);
        b.tri([at(outer[k], 0.0), at(outer[n], 1.0), at(outer[k], 1.0)]);
        // Inner wall faces into the hole.
        b.tri([at(inner[k], 0.0), at(inner[n], 1.0), at(inner[n], 0.0)]);
        b.tri([at(inner[k], 0.0), at(inner[k], 1.0), at(inner[n], 1.0)]);
    }
    b.finish()
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
            mesh: prism(),
            frame: h(0.0),
            expected: Expected::Output(output(vec![], vec![], 0, 1)),
        },
        // The prism top at z = 1 faces +z: the quadrilateral v2 v5 v4 v1,
        // (0,-1) (2,-1) (2,1) (0,1), counter-clockwise from (0,-1).
        Case {
            name: "prism top",
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
    ]
}

#[cfg(test)]
mod tests {
    use axiolid_core::Tolerance;
    use axiolid_mesh::try_audit_mesh;
    use axiolid_mesh_contracts::SolidRequirements;

    use super::*;

    #[test]
    fn every_fixture_is_a_closed_oriented_two_manifold() {
        let tolerance = Tolerance::new(1e-9, 1e-9).expect("valid tolerance");
        for case in cases() {
            let health = try_audit_mesh(&case.mesh, tolerance).expect("audit");
            assert!(health.is_closed_two_manifold(), "{}: {health:?}", case.name);
            SolidRequirements::Oriented
                .validate(&case.mesh, case.name)
                .unwrap_or_else(|error| panic!("{}: {error}", case.name));
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
