//! Self-contained scenes and expected values for the conformance suite.
//!
//! Every number here is derived by hand in the comment beside it and checked
//! again by the unit tests at the bottom of this file, so a wrong fixture
//! fails here before it can fail a provider.

use axiolid_contracts::GeomError;
use axiolid_core::{Aabb2, Frame3, Point2, Point3, Vec3};
use axiolid_mesh::TriMesh;

use crate::{
    EdgeKind, EdgeSelection, HiddenLineLimits, MeshId, MeshRole, OrthographicView, SourceEdge,
    Visibility,
};

/// One mesh of a scene.
pub(crate) struct Fixture {
    pub id: MeshId,
    pub role: MeshRole,
    pub mesh: TriMesh,
}

/// Expected count and total length of the segments of one kind and visibility.
/// A kind and visibility without a row must have no segments.
pub(crate) struct Row {
    pub kind: EdgeKind,
    pub visibility: Visibility,
    pub count: usize,
    pub length: f64,
}

/// One expected piece of a source edge.
pub(crate) struct Piece {
    pub visibility: Visibility,
    pub start: [f64; 2],
    pub end: [f64; 2],
}

/// The expected pieces of one source edge, in order.
/// All pieces of one edge share one kind.
pub(crate) struct EdgePieces {
    pub edge: SourceEdge,
    pub kind: EdgeKind,
    pub pieces: Vec<Piece>,
}

#[derive(Default)]
pub(crate) struct Expect {
    /// Per kind and visibility; `None` checks only the listed edges.
    pub rows: Option<Vec<Row>>,
    pub edges: Vec<EdgePieces>,
    pub selected_edges: Option<usize>,
    pub point_edges: Option<usize>,
    pub clipped_edges: Option<usize>,
}

/// A scene with a drawing the provider must reproduce.
pub(crate) struct Case {
    pub name: &'static str,
    pub meshes: Vec<Fixture>,
    pub view: OrthographicView,
    pub edges: EdgeSelection,
    pub region: Option<Aabb2>,
    pub limits: HiddenLineLimits,
    pub expect: Expect,
    /// Also run with every triangle reversed and require identical output.
    pub orientation_independent: bool,
}

/// A scene the provider must refuse.
pub(crate) struct Refusal {
    pub name: &'static str,
    pub meshes: Vec<Fixture>,
    pub view: OrthographicView,
    pub edges: EdgeSelection,
    pub limits: HiddenLineLimits,
    /// The refusing variant.
    pub accepts: fn(&GeomError) -> bool,
    /// The mesh the refusal must name; `None` for the whole request.
    pub mesh: Option<MeshId>,
}

pub(crate) const LENGTH_TOLERANCE: f64 = 1e-12;

pub(crate) fn limits() -> HiddenLineLimits {
    HiddenLineLimits::new(8, 100_000, 100_000, 10_000_000, 100_000)
}

/// Axis-aligned box numbered and wound like the section cube: vertices 0..4
/// are the bottom face counter-clockwise from `min`, 4..8 the same at the top.
pub(crate) fn box_mesh(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    TriMesh::new(
        vec![
            Point3::new(x0, y0, z0),
            Point3::new(x1, y0, z0),
            Point3::new(x1, y1, z0),
            Point3::new(x0, y1, z0),
            Point3::new(x0, y0, z1),
            Point3::new(x1, y0, z1),
            Point3::new(x1, y1, z1),
            Point3::new(x0, y1, z1),
        ],
        vec![
            0, 2, 1, 0, 3, 2, // bottom
            4, 5, 6, 4, 6, 7, // top
            0, 1, 5, 0, 5, 4, // front
            1, 2, 6, 1, 6, 5, // right
            2, 3, 7, 2, 7, 6, // back
            3, 0, 4, 3, 4, 7, // left
        ],
    )
}

/// Horizontal quad at height `z`, counter-clockwise seen from `+z`:
/// vertices `(x0,y0) (x1,y0) (x1,y1) (x0,y1)`.
pub(crate) fn sheet_quad(x0: f64, y0: f64, x1: f64, y1: f64, z: f64) -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(x0, y0, z),
            Point3::new(x1, y0, z),
            Point3::new(x1, y1, z),
            Point3::new(x0, y1, z),
        ],
        vec![0, 1, 2, 0, 2, 3],
    )
}

/// Quad in the plane x=1 over y in [0,1], z in [0,2], seen edge-on in plan:
/// vertices `(1,0,0) (1,1,0) (1,1,2) (1,0,2)`.
pub(crate) fn vertical_sheet() -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 2.0),
            Point3::new(1.0, 0.0, 2.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    )
}

/// Tetrahedron on four points, wound outward whatever their order.
pub(crate) fn tetrahedron(points: [Point3; 4]) -> TriMesh {
    let mut indices = vec![0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
    let mesh = TriMesh::new(points.to_vec(), indices.clone());
    if signed_volume(&mesh) < 0.0 {
        for triangle in indices.as_chunks_mut::<3>().0.iter_mut() {
            triangle.swap(1, 2);
        }
    }
    TriMesh::new(points.to_vec(), indices)
}

/// The same mesh with every triangle reversed.
pub(crate) fn reversed(mesh: &TriMesh) -> TriMesh {
    let mut indices = mesh.indices.clone();
    for triangle in indices.as_chunks_mut::<3>().0.iter_mut() {
        triangle.swap(1, 2);
    }
    TriMesh::new(mesh.positions.clone(), indices)
}

/// Signed enclosed volume of a closed mesh; positive when wound outward.
pub(crate) fn signed_volume(mesh: &TriMesh) -> f64 {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [t[0], t[1], t[2]].map(|i| mesh.positions[i as usize]);
            a.dot(b.cross(c)) / 6.0
        })
        .sum()
}

/// Plan view: the world axes, looking down `-z`.
pub(crate) fn plan_frame() -> Frame3 {
    Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

/// Isometric view of the first octant: the viewer sits on `(1,1,1)`.
pub(crate) fn isometric_frame(origin: Point3) -> Frame3 {
    let s2 = core::f64::consts::FRAC_1_SQRT_2;
    let s6 = 1.0 / 6.0_f64.sqrt();
    let s3 = 1.0 / 3.0_f64.sqrt();
    Frame3 {
        origin,
        x: Vec3::new(s2, -s2, 0.0),
        y: Vec3::new(s6, s6, -2.0 * s6),
        z: Vec3::new(s3, s3, s3),
    }
}

/// Image length of a unit cube edge in the isometric view: `sqrt(1 - 1/3)`.
pub(crate) fn iso_edge() -> f64 {
    (2.0_f64 / 3.0).sqrt()
}

fn solid(id: u64, mesh: TriMesh) -> Fixture {
    Fixture {
        id: MeshId(id),
        role: MeshRole::Solid,
        mesh,
    }
}

fn sheet(id: u64, mesh: TriMesh) -> Fixture {
    Fixture {
        id: MeshId(id),
        role: MeshRole::Sheet,
        mesh,
    }
}

fn edge(mesh: u64, a: u32, b: u32) -> SourceEdge {
    SourceEdge {
        mesh: MeshId(mesh),
        vertices: [a, b],
    }
}

fn piece(visibility: Visibility, start: [f64; 2], end: [f64; 2]) -> Piece {
    Piece {
        visibility,
        start,
        end,
    }
}

fn row(kind: EdgeKind, visibility: Visibility, count: usize, length: f64) -> Row {
    Row {
        kind,
        visibility,
        count,
        length,
    }
}

fn unit_cube() -> TriMesh {
    box_mesh([0.0; 3], [1.0; 3])
}

fn plan(view: Frame3) -> OrthographicView {
    OrthographicView::new(view)
}

fn crease_half() -> EdgeSelection {
    EdgeSelection::new(Some(0.5))
}

fn case(
    name: &'static str,
    meshes: Vec<Fixture>,
    view: OrthographicView,
    edges: EdgeSelection,
    expect: Expect,
) -> Case {
    Case {
        name,
        meshes,
        view,
        edges,
        region: None,
        limits: limits(),
        expect,
        orientation_independent: false,
    }
}

/// The isometric cube: the near corner is vertex 6, the far corner vertex 0.
/// The six edges of the hexagonal outline are silhouettes; the three edges
/// at each corner are creases, visible at the near corner and hidden at the
/// far one. Every cube edge has image length `sqrt(2/3)`.
fn isometric_cube(
    name: &'static str,
    origin: Point3,
    edges: EdgeSelection,
    rows: Vec<Row>,
    selected: usize,
) -> Case {
    let expect = Expect {
        rows: Some(rows),
        selected_edges: Some(selected),
        point_edges: Some(0),
        clipped_edges: Some(0),
        ..Expect::default()
    };
    let mut case = case(
        name,
        vec![solid(1, unit_cube())],
        plan(isometric_frame(origin)),
        edges,
        expect,
    );
    case.orientation_independent = true;
    case
}

pub(crate) fn cases() -> Vec<Case> {
    use EdgeKind::{Boundary, Crease, Silhouette};
    use Visibility::{Hidden, Visible};
    let l = iso_edge();
    let mut cases = Vec::new();

    let full = || {
        vec![
            row(Silhouette, Visible, 6, 6.0 * l),
            row(Crease, Visible, 3, 3.0 * l),
            row(Crease, Hidden, 3, 3.0 * l),
        ]
    };
    cases.push(isometric_cube(
        "isometric cube",
        Point3::ZERO,
        crease_half(),
        full(),
        12,
    ));
    cases.push(isometric_cube(
        "isometric cube, offset origin",
        Point3::new(0.1, 0.2, 0.3),
        crease_half(),
        full(),
        12,
    ));
    cases.push(isometric_cube(
        "isometric cube, no creases",
        Point3::ZERO,
        EdgeSelection::new(None),
        vec![row(Silhouette, Visible, 6, 6.0 * l)],
        6,
    ));
    cases.push(isometric_cube(
        "isometric cube, visible only",
        Point3::ZERO,
        crease_half().with_hidden(false),
        vec![
            row(Silhouette, Visible, 6, 6.0 * l),
            row(Crease, Visible, 3, 3.0 * l),
        ],
        12,
    ));

    // Box behind box in plan. F = [0,2]^2 x [2,3] hides part of
    // K = [1,3]^2 x [0,1]. Silhouettes are the top and bottom ring edges
    // (a side face is edge-on, its own orientation); the 4 vertical edges of
    // each box are points in plan: 8 point edges, 24 selected edges.
    //   visible: F top ring 4 x 2 = 8; K top ring outside F: the edges y=1
    //   and x=1 keep [2,3] (1 each), the edges y=3 and x=3 keep 2 each: 6.
    //     total 14 in 8 pieces.
    //   hidden: F bottom ring 8, K top ring inside F 2, K bottom ring 8
    //     (behind K's own top): 18 in 10 pieces.
    cases.push(case(
        "box behind box",
        vec![
            solid(1, box_mesh([0.0, 0.0, 2.0], [2.0, 2.0, 3.0])),
            solid(2, box_mesh([1.0, 1.0, 0.0], [3.0, 3.0, 1.0])),
        ],
        plan(plan_frame()),
        crease_half(),
        Expect {
            rows: Some(vec![
                row(Silhouette, Visible, 8, 14.0),
                row(Silhouette, Hidden, 10, 18.0),
            ]),
            edges: vec![
                EdgePieces {
                    edge: edge(2, 4, 5),
                    kind: Silhouette,
                    pieces: vec![
                        piece(Hidden, [1.0, 1.0], [2.0, 1.0]),
                        piece(Visible, [2.0, 1.0], [3.0, 1.0]),
                    ],
                },
                EdgePieces {
                    edge: edge(2, 4, 7),
                    kind: Silhouette,
                    pieces: vec![
                        piece(Hidden, [1.0, 1.0], [1.0, 2.0]),
                        piece(Visible, [1.0, 2.0], [1.0, 3.0]),
                    ],
                },
            ],
            selected_edges: Some(24),
            point_edges: Some(8),
            clipped_edges: Some(0),
        },
    ));

    // A tetrahedron ridge R0 -> R1 pierces a plate [0,2]^2 x [1,2]. On the
    // ridge z = 4t and x = -1 + 4t: it enters the plate's top image at x=0
    // (z=1, under the top) and rises through the top plane z=2 at x=1.
    cases.push(case(
        "interpenetration",
        vec![
            solid(1, box_mesh([0.0, 0.0, 1.0], [2.0, 2.0, 2.0])),
            solid(2, tetrahedron(ridge_points())),
        ],
        plan(plan_frame()),
        crease_half(),
        Expect {
            edges: vec![EdgePieces {
                edge: edge(2, 0, 1),
                kind: Crease,
                pieces: vec![
                    piece(Visible, [-1.0, 1.0], [0.0, 1.0]),
                    piece(Hidden, [0.0, 1.0], [1.0, 1.0]),
                    piece(Visible, [1.0, 1.0], [3.0, 1.0]),
                ],
            }],
            ..Expect::default()
        },
    ));

    // An open sheet far above a cube: its 4 boundary edges are visible
    // (12 in all), the cube's 8 horizontal edges are hidden under it.
    let mut sheet_over_cube = case(
        "sheet over solid",
        vec![
            sheet(1, sheet_quad(-1.0, -1.0, 2.0, 2.0, 5.0)),
            solid(2, unit_cube()),
        ],
        plan(plan_frame()),
        crease_half(),
        Expect {
            rows: Some(vec![
                row(Boundary, Visible, 4, 12.0),
                row(Silhouette, Hidden, 8, 8.0),
            ]),
            selected_edges: Some(16),
            point_edges: Some(4),
            clipped_edges: Some(0),
            ..Expect::default()
        },
    );
    sheet_over_cube.orientation_independent = true;
    cases.push(sheet_over_cube);

    // The same scene with one edge kind deselected: the sheet's 4 boundary
    // edges alone (12 visible) or the cube's 8 horizontal silhouettes alone
    // (8 hidden). The cube's 4 vertical creases are points in plan, so they
    // are selected and counted as point edges whenever creases are on:
    // 4 + 8 = 12 without the boundary; 4 + 4 = 8 without silhouettes.
    for (name, edges, rows, selected) in [
        (
            "sheet over solid, no boundary",
            crease_half().with_boundary(false),
            vec![row(Silhouette, Hidden, 8, 8.0)],
            12,
        ),
        (
            "sheet over solid, no silhouettes",
            crease_half().with_silhouette(false),
            vec![row(Boundary, Visible, 4, 12.0)],
            8,
        ),
    ] {
        cases.push(case(
            name,
            vec![
                sheet(1, sheet_quad(-1.0, -1.0, 2.0, 2.0, 5.0)),
                solid(2, unit_cube()),
            ],
            plan(plan_frame()),
            edges,
            Expect {
                rows: Some(rows),
                selected_edges: Some(selected),
                point_edges: Some(4),
                clipped_edges: Some(0),
                ..Expect::default()
            },
        ));
    }

    // Touching: equal depth never occludes. A's edge (5,6) and B's edge
    // (4,7) are the same segment x=1, y in [0,1], at z=1, under a sheet at
    // z=1; all of it stays visible.
    let sheet_edge = |a: u32, b: u32, start: [f64; 2], end: [f64; 2]| EdgePieces {
        edge: edge(3, a, b),
        kind: Boundary,
        pieces: vec![piece(Visible, start, end)],
    };
    cases.push(case(
        "touching at equal depth",
        vec![
            solid(1, unit_cube()),
            solid(2, box_mesh([1.0, 0.0, 0.0], [2.0, 1.0, 1.0])),
            sheet(3, sheet_quad(0.5, 0.25, 1.5, 0.75, 1.0)),
        ],
        plan(plan_frame()),
        crease_half(),
        Expect {
            edges: vec![
                EdgePieces {
                    edge: edge(1, 5, 6),
                    kind: Silhouette,
                    pieces: vec![piece(Visible, [1.0, 0.0], [1.0, 1.0])],
                },
                EdgePieces {
                    edge: edge(2, 4, 7),
                    kind: Silhouette,
                    pieces: vec![piece(Visible, [1.0, 0.0], [1.0, 1.0])],
                },
                sheet_edge(0, 1, [0.5, 0.25], [1.5, 0.25]),
                sheet_edge(1, 2, [1.5, 0.25], [1.5, 0.75]),
                sheet_edge(2, 3, [1.5, 0.75], [0.5, 0.75]),
                sheet_edge(0, 3, [0.5, 0.25], [0.5, 0.75]),
            ],
            ..Expect::default()
        },
    ));

    // A tall box's wall x=1 is edge-on in plan and spans z in [0,2]; it
    // hides the edge of the short neighbour at z=1 through its depth range.
    cases.push(case(
        "edge-on wall hides through its depth",
        vec![
            solid(1, box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 2.0])),
            solid(2, box_mesh([1.0, 0.0, 0.0], [2.0, 1.0, 1.0])),
        ],
        plan(plan_frame()),
        crease_half(),
        Expect {
            edges: vec![EdgePieces {
                edge: edge(2, 4, 7),
                kind: Silhouette,
                pieces: vec![piece(Hidden, [1.0, 0.0], [1.0, 1.0])],
            }],
            ..Expect::default()
        },
    ));

    // A vertical sheet in the plane x=1 (y in [0,1], z in [0,2]) is edge-on
    // in plan; it hides B's edge (4,7) at z=1 only through its fibre, since
    // B is a solid whose own faces would not (the edge lies on B's
    // boundary at equal depth). The fibre reaches z=2 > 1: hidden over the
    // whole edge. Its own edge (0,1) at z=0 is hidden by the same fibre; its
    // top edge (2,3) at z=2 is not strictly under anything: visible. Its
    // other two boundary edges and B's 4 verticals are points: 6 point
    // edges; 4 + 12 = 16 selected edges.
    cases.push(case(
        "edge-on sheet hides through its fibre",
        vec![
            sheet(1, vertical_sheet()),
            solid(2, box_mesh([1.0, 0.0, 0.0], [2.0, 1.0, 1.0])),
        ],
        plan(plan_frame()),
        crease_half(),
        Expect {
            edges: vec![
                EdgePieces {
                    edge: edge(2, 4, 7),
                    kind: Silhouette,
                    pieces: vec![piece(Hidden, [1.0, 0.0], [1.0, 1.0])],
                },
                EdgePieces {
                    edge: edge(1, 0, 1),
                    kind: Boundary,
                    pieces: vec![piece(Hidden, [1.0, 0.0], [1.0, 1.0])],
                },
                EdgePieces {
                    edge: edge(1, 2, 3),
                    kind: Boundary,
                    pieces: vec![piece(Visible, [1.0, 1.0], [1.0, 0.0])],
                },
            ],
            selected_edges: Some(16),
            point_edges: Some(6),
            clipped_edges: Some(0),
            ..Expect::default()
        },
    ));

    // Far clipping at 0.5 removes the 4 bottom-ring edges (z=0) entirely;
    // the 4 verticals are points; the top ring (z=1) stays visible.
    cases.push(case(
        "far clipping",
        vec![solid(1, unit_cube())],
        OrthographicView::new(plan_frame()).with_far(0.5),
        crease_half(),
        Expect {
            rows: Some(vec![row(Silhouette, Visible, 4, 4.0)]),
            selected_edges: Some(12),
            point_edges: Some(4),
            clipped_edges: Some(4),
            ..Expect::default()
        },
    ));

    // Near clipping at 1 removes the sheet at z=3 (4 boundary edges, so it
    // neither draws nor occludes); the cube, declared a sheet, lies at or
    // behind the near plane, which is a closed bound.
    cases.push(case(
        "near clipping",
        vec![
            sheet(1, sheet_quad(-1.0, -1.0, 2.0, 2.0, 3.0)),
            sheet(2, unit_cube()),
        ],
        OrthographicView::new(plan_frame()).with_near(1.0),
        crease_half(),
        Expect {
            rows: Some(vec![
                row(Silhouette, Visible, 4, 4.0),
                row(Silhouette, Hidden, 4, 4.0),
            ]),
            selected_edges: Some(16),
            point_edges: Some(4),
            clipped_edges: Some(4),
            ..Expect::default()
        },
    ));

    // A crop to x in [0.5, 2]: of the top ring the edges y=0 and y=1 keep
    // 0.5 each and x=1 keeps 1 (2 in all); x=0 is outside, on both rings
    // (2 clipped edges); the bottom ring is hidden the same way.
    let mut cropped = case(
        "region crop",
        vec![solid(1, unit_cube())],
        plan(plan_frame()),
        crease_half(),
        Expect {
            rows: Some(vec![
                row(Silhouette, Visible, 3, 2.0),
                row(Silhouette, Hidden, 3, 2.0),
            ]),
            selected_edges: Some(12),
            point_edges: Some(4),
            clipped_edges: Some(2),
            ..Expect::default()
        },
    );
    cropped.region = Some(Aabb2 {
        min: Point2::new(0.5, -1.0),
        max: Point2::new(2.0, 2.0),
    });
    cases.push(cropped);

    cases
}

/// The ridge of the interpenetration scene: R0, R1 and the two low points.
pub(crate) fn ridge_points() -> [Point3; 4] {
    [
        Point3::new(-1.0, 1.0, 0.0),
        Point3::new(3.0, 1.0, 4.0),
        Point3::new(1.0, -1.0, -2.0),
        Point3::new(1.0, 3.0, -2.0),
    ]
}

fn soup(positions: Vec<Point3>, indices: Vec<u32>) -> TriMesh {
    TriMesh::new(positions, indices)
}

fn triangle_soup_points() -> Vec<Point3> {
    vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, -1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
    ]
}

fn refusal(name: &'static str, offender: Fixture, accepts: fn(&GeomError) -> bool) -> Refusal {
    Refusal {
        name,
        mesh: Some(offender.id),
        meshes: vec![solid(1, unit_cube()), offender],
        view: plan(plan_frame()),
        edges: crease_half(),
        limits: limits(),
        accepts,
    }
}

pub(crate) fn refusals() -> Vec<Refusal> {
    let not_manifold: fn(&GeomError) -> bool = |e| matches!(e, GeomError::NotManifold(_));
    let mut refusals = vec![
        // Edge (0,1) carries three triangles.
        refusal(
            "three triangles on one edge",
            sheet(
                2,
                soup(triangle_soup_points(), vec![0, 1, 2, 1, 0, 3, 0, 1, 4]),
            ),
            not_manifold,
        ),
        // Both triangles traverse edge (0,1) from 0 to 1.
        refusal(
            "inconsistent winding",
            sheet(2, soup(triangle_soup_points(), vec![0, 1, 2, 0, 1, 3])),
            not_manifold,
        ),
        // An open quad is no solid.
        refusal(
            "sheet declared solid",
            solid(2, sheet_quad(0.0, 0.0, 1.0, 1.0, 0.0)),
            not_manifold,
        ),
    ];
    let mut output = Refusal {
        name: "output limit",
        meshes: vec![solid(1, unit_cube())],
        view: plan(isometric_frame(Point3::ZERO)),
        edges: crease_half(),
        limits: limits(),
        accepts: |e| {
            matches!(
                e,
                GeomError::BudgetExceeded {
                    resource: "hidden-line output segments"
                }
            )
        },
        mesh: None,
    };
    output.limits.max_output_segments = 1;
    refusals.push(output);
    let mut pairs = Refusal {
        name: "candidate limit",
        meshes: vec![solid(1, unit_cube())],
        view: plan(isometric_frame(Point3::ZERO)),
        edges: crease_half(),
        limits: limits(),
        accepts: |e| {
            matches!(
                e,
                GeomError::BudgetExceeded {
                    resource: "hidden-line candidate pairs"
                }
            )
        },
        mesh: None,
    };
    pairs.limits.max_candidate_pairs = 0;
    refusals.push(pairs);
    refusals
}

#[cfg(test)]
mod tests {
    use axiolid_core::Tolerance;
    use axiolid_mesh::audit_mesh;

    use super::*;

    const EPS: f64 = 1e-12;

    fn assert_closed_manifold(mesh: &TriMesh) {
        let health = audit_mesh(mesh, Tolerance::ZERO);
        assert!(health.is_closed_two_manifold(), "{health:?}");
    }

    fn assert_orthonormal_right_handed(frame: &Frame3) {
        for axis in [frame.x, frame.y, frame.z] {
            assert!((axis.length() - 1.0).abs() < EPS);
        }
        assert!(frame.x.dot(frame.y).abs() < EPS);
        assert!(frame.y.dot(frame.z).abs() < EPS);
        assert!(frame.z.dot(frame.x).abs() < EPS);
        assert!((frame.x.cross(frame.y) - frame.z).length() < EPS);
    }

    #[test]
    fn boxes_are_closed_outward_and_have_the_box_volume() {
        for (min, max) in [
            ([0.0; 3], [1.0; 3]),
            ([0.0, 0.0, 2.0], [2.0, 2.0, 3.0]),
            ([1.0, 1.0, 0.0], [3.0, 3.0, 1.0]),
            ([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]),
        ] {
            let mesh = box_mesh(min, max);
            assert_closed_manifold(&mesh);
            let volume: f64 = (0..3).map(|a| max[a] - min[a]).product();
            assert!((signed_volume(&mesh) - volume).abs() < EPS);
            assert!((signed_volume(&reversed(&mesh)) + volume).abs() < EPS);
        }
    }

    #[test]
    fn tetrahedron_is_outward_in_either_point_order() {
        let points = ridge_points();
        let mesh = tetrahedron(points);
        assert_closed_manifold(&mesh);
        assert!(signed_volume(&mesh) > 0.0);
        // The same points in the opposite parity must also come out outward.
        let swapped = tetrahedron([points[1], points[0], points[2], points[3]]);
        assert!(signed_volume(&swapped) > 0.0);
        // Edge-vector triple product: |det| / 6.
        let [r0, r1, q0, q1] = points;
        let det = (r1 - r0).dot((q0 - r0).cross(q1 - r0));
        assert!((signed_volume(&mesh) - det.abs() / 6.0).abs() < EPS);
    }

    #[test]
    fn sheet_quad_faces_up_and_is_open() {
        let mesh = sheet_quad(0.0, 0.0, 2.0, 3.0, 1.0);
        let health = audit_mesh(&mesh, Tolerance::ZERO);
        assert_eq!(health.boundary_edges, 4);
        assert_eq!(health.non_manifold_edges, 0);
        assert_eq!(health.inconsistent_winding_edges, 0);
        let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[i]);
        assert!((b - a).cross(c - a).z > 0.0);
    }

    #[test]
    fn frames_are_orthonormal_and_right_handed() {
        assert_orthonormal_right_handed(&plan_frame());
        assert_orthonormal_right_handed(&isometric_frame(Point3::ZERO));
        assert_orthonormal_right_handed(&isometric_frame(Point3::new(0.1, 0.2, 0.3)));
    }

    #[test]
    fn isometric_cube_edges_project_to_sqrt_two_thirds() {
        let f = isometric_frame(Point3::ZERO);
        for axis in [Vec3::X, Vec3::Y, Vec3::Z] {
            let image = (axis.dot(f.x).powi(2) + axis.dot(f.y).powi(2)).sqrt();
            assert!((image - iso_edge()).abs() < EPS);
            // The depth gained along any cube axis is 1/sqrt(3).
            assert!((axis.dot(f.z) - 1.0 / 3.0_f64.sqrt()).abs() < EPS);
        }
        // Corner 6 is nearest and corner 0 farthest; both sit at the image
        // centre, so the hexagon is the outline of the other six corners.
        let mesh = unit_cube();
        let depth = |i: usize| mesh.positions[i].dot(f.z);
        for i in 0..8 {
            assert!(depth(i) <= depth(6) + EPS && depth(i) >= depth(0) - EPS);
        }
        let image = |i: usize| (mesh.positions[i].dot(f.x), mesh.positions[i].dot(f.y));
        assert!((image(0).0 - image(6).0).abs() < EPS && (image(0).1 - image(6).1).abs() < EPS);
    }

    /// Length of `[lo, hi]` outside the closed interval `[a, b]`.
    fn outside(lo: f64, hi: f64, a: f64, b: f64) -> f64 {
        (hi - lo) - (hi.min(b) - lo.max(a)).max(0.0)
    }

    #[test]
    fn box_behind_box_totals_are_the_analytic_lengths() {
        // F = [0,2]^2 x [2,3], K = [1,3]^2 x [0,1]; F covers [0,2]^2.
        let f_top = 4.0 * 2.0;
        // K's top ring: x in [1,3] at y=1 and y=3, y in [1,3] at x=1 and x=3.
        // Along y=1 and x=1 the covered part is [1,2]; along y=3 and x=3 the
        // line lies outside F's square (3 > 2).
        let k_visible = outside(1.0, 3.0, 0.0, 2.0) * 2.0 + 2.0 * 2.0;
        assert!((f_top + k_visible - 14.0).abs() < EPS);
        let k_hidden_top = (3.0 - 1.0) - outside(1.0, 3.0, 0.0, 2.0);
        let hidden = f_top + 2.0 * k_hidden_top + 4.0 * 2.0;
        assert!((hidden - 18.0).abs() < EPS);
        // Visible pieces: F's 4 top edges, one on each of K's 4 top edges.
        // Hidden pieces: F's 4 bottom edges, K's 2 covered top edges, K's
        // 4 bottom edges.
        let visible_pieces = 4 + 4;
        let hidden_pieces = 4 + 2 + 4;
        assert_eq!((visible_pieces, hidden_pieces), (8, 10));
    }

    #[test]
    fn ridge_crosses_the_plate_where_the_fixture_says() {
        let [r0, r1, q0, q1] = ridge_points();
        let at = |t: f64| r0 + (r1 - r0) * t;
        // z = 4t, x = -1 + 4t, y = 1.
        assert!((at(0.25).x - 0.0).abs() < EPS && (at(0.25).z - 1.0).abs() < EPS);
        assert!((at(0.5).x - 1.0).abs() < EPS && (at(0.5).z - 2.0).abs() < EPS);
        assert!((at(0.5).y - 1.0).abs() < EPS);
        // Both faces at the ridge face the viewer and meet at a crease.
        let mesh = tetrahedron([r0, r1, q0, q1]);
        let normal = |t: usize| {
            let [a, b, c] = [0, 1, 2].map(|k| mesh.positions[mesh.indices[3 * t + k] as usize]);
            (b - a).cross(c - a).normalize()
        };
        let ridge_faces: Vec<usize> = (0..4)
            .filter(|&t| {
                let tri = &mesh.indices[3 * t..3 * t + 3];
                tri.contains(&0) && tri.contains(&1)
            })
            .collect();
        assert_eq!(ridge_faces.len(), 2);
        let (n1, n2) = (normal(ridge_faces[0]), normal(ridge_faces[1]));
        assert!(n1.z > 0.0 && n2.z > 0.0);
        assert!((n1.dot(n2) + 1.0 / 3.0).abs() < 1e-12);
        assert!(n1.dot(n2) < 0.5);
        // The plate [0,2]^2 x [1,2] contains the ridge's image y=1.
        let plate = box_mesh([0.0, 0.0, 1.0], [2.0, 2.0, 2.0]);
        assert!(plate.positions[0].y < r0.y && r0.y < plate.positions[2].y);
    }

    #[test]
    fn touching_fixture_shares_the_edge_and_the_depth() {
        let a = box_mesh([0.0; 3], [1.0; 3]);
        let b = box_mesh([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        let s = sheet_quad(0.5, 0.25, 1.5, 0.75, 1.0);
        assert_eq!(a.positions[5], b.positions[4]);
        assert_eq!(a.positions[6], b.positions[7]);
        assert_eq!(a.positions[5].z, s.positions[0].z);
        assert_eq!(b.positions[4].z, s.positions[0].z);
        // The vertical sheet lies in the plane x=1 of B edge (4,7), over
        // its image, and reaches above it.
        let wall = vertical_sheet();
        assert!(wall.positions.iter().all(|p| p.x == b.positions[4].x));
        assert!(wall.positions.iter().any(|p| p.z > b.positions[4].z));
        assert!(wall.positions.iter().any(|p| p.z < b.positions[4].z));
        assert_eq!(
            (wall.positions[0].y, wall.positions[1].y),
            (b.positions[4].y, b.positions[7].y)
        );
        let tall = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 2.0]);
        assert_eq!(tall.positions[5].x, b.positions[4].x);
        assert!(tall.positions[5].z > b.positions[4].z);
    }

    #[test]
    fn clipping_fixtures_match_their_arithmetic() {
        // Far 0.5: bottom ring at z=0 is out, top ring at z=1 is in.
        let cube = unit_cube();
        assert!(cube.positions[..4].iter().all(|p| p.z < 0.5));
        assert!(cube.positions[4..].iter().all(|p| p.z >= 0.5));
        // Near 1: the cube's top z=1 is kept (closed), the sheet at z=3 not.
        assert!(cube.positions.iter().all(|p| p.z <= 1.0));
        assert!(sheet_quad(-1.0, -1.0, 2.0, 2.0, 3.0).positions[0].z > 1.0);
        // Region x in [0.5, 2]: kept lengths of the cube's top ring.
        let kept = |a: f64, b: f64| (b.min(2.0) - a.max(0.5)).max(0.0);
        assert!((kept(0.0, 1.0) * 2.0 + 1.0 - 2.0).abs() < EPS);
        assert_eq!(kept(0.0, 0.0), 0.0);
    }

    #[test]
    fn every_case_is_well_formed() {
        let mut names = std::collections::BTreeSet::new();
        for case in cases() {
            assert!(names.insert(case.name), "duplicate case {}", case.name);
            let ids: std::collections::BTreeSet<_> = case.meshes.iter().map(|m| m.id).collect();
            assert_eq!(ids.len(), case.meshes.len(), "{}", case.name);
            for fixture in &case.meshes {
                let health = audit_mesh(&fixture.mesh, Tolerance::ZERO);
                assert_eq!(health.degenerate_triangles, 0, "{}", case.name);
                assert_eq!(health.non_manifold_edges, 0, "{}", case.name);
                assert_eq!(health.inconsistent_winding_edges, 0, "{}", case.name);
                if fixture.role == MeshRole::Solid {
                    assert!(health.is_closed_two_manifold(), "{}", case.name);
                }
            }
            if let Some(rows) = &case.expect.rows {
                let mut seen = std::collections::BTreeSet::new();
                for r in rows {
                    let key = (r.kind as u8, r.visibility as u8);
                    assert!(seen.insert(key), "{}: duplicate row", case.name);
                    assert!(r.count > 0 && r.length > 0.0);
                }
            }
            for pieces in &case.expect.edges {
                assert!(pieces.edge.vertices[0] < pieces.edge.vertices[1]);
                assert!(case.meshes.iter().any(|m| m.id == pieces.edge.mesh));
                for p in &pieces.pieces {
                    assert_ne!(p.start, p.end);
                }
            }
        }
        assert_eq!(names.len(), 15);
    }

    #[test]
    fn every_refusal_names_a_mesh_in_its_request() {
        for refusal in refusals() {
            if let Some(id) = refusal.mesh {
                assert!(
                    refusal.meshes.iter().any(|m| m.id == id),
                    "{}",
                    refusal.name
                );
                // The offender is the second mesh, so a provider that blames
                // the first or the request fails the check.
                assert_eq!(refusal.meshes[1].id, id);
            }
        }
    }
}
