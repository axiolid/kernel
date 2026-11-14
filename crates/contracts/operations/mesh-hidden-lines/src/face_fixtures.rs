//! Self-contained scenes and expected values for the visible-face
//! conformance suite.
//!
//! Every ring here is derived by hand in the comment beside it. Image
//! coordinates in the plan view are the world `x` and `y`, and the depth is
//! `z`, so most scenes read off a drawing in the plane. The unit tests at the
//! bottom check each expected ring (orientation, canonical start, area) and
//! every scene's mesh, so a wrong fixture fails here before it can fail a
//! provider.

use axiolid_contracts::GeomError;
use axiolid_core::{Aabb2, Frame3, Point2, Point3, Vec3};
use axiolid_mesh::TriMesh;

mod snap;

pub(crate) use snap::snap_cases;

use crate::fixtures::{box_mesh, isometric_frame, plan_frame, refusals, sheet_quad};
use crate::MeshRole::{Sheet, Solid};
use crate::{
    BehindSurface, FaceGroup, FaceSource, MeshId, MeshRole, Occlusion, View, VisibleFaceLimits,
    VisibleRegion,
};

/// One mesh of a scene. An empty `groups` means every triangle is in
/// `FaceGroup(0)`. A `planar` fixture is requested with
/// `FaceGroups::Planar`; its `groups` are then the groups that request must
/// resolve to, used by the checks and not sent. An empty `occlusion` sends no
/// flags; otherwise it is sent `PerTriangle`.
pub(crate) struct FaceFixture {
    pub id: MeshId,
    pub role: MeshRole,
    pub mesh: TriMesh,
    pub groups: Vec<FaceGroup>,
    pub planar: bool,
    pub occlusion: Vec<Occlusion>,
}

impl FaceFixture {
    fn new(id: u64, role: MeshRole, mesh: TriMesh) -> Self {
        Self {
            id: MeshId(id),
            role,
            mesh,
            groups: Vec::new(),
            planar: false,
            occlusion: Vec::new(),
        }
    }

    /// Every triangle see-through.
    fn see_through(mut self) -> Self {
        self.occlusion = vec![Occlusion::SeeThrough; self.mesh.triangle_count()];
        self
    }

    /// The same mesh with every triangle reversed; triangle order, groups and
    /// flags are kept.
    pub(crate) fn flipped(&self) -> Self {
        Self {
            id: self.id,
            role: self.role,
            mesh: crate::fixtures::reversed(&self.mesh),
            groups: self.groups.clone(),
            planar: self.planar,
            occlusion: self.occlusion.clone(),
        }
    }

    /// Request planar groups, which resolve to `resolved`.
    fn planar(mut self, resolved: &[u32]) -> Self {
        self.groups = resolved.iter().map(|&g| FaceGroup(g)).collect();
        self.planar = true;
        self
    }

    fn uniform(mut self, group: u32) -> Self {
        self.groups = vec![FaceGroup(group); self.mesh.triangle_count()];
        self
    }

    /// The group of triangle `triangle`.
    pub fn group(&self, triangle: usize) -> FaceGroup {
        self.groups.get(triangle).copied().unwrap_or(FaceGroup(0))
    }
}

/// What the evidence of a case must show beyond the generic checks.
#[derive(Default)]
pub(crate) struct EvidenceExpect {
    /// `intersection_segments` is positive.
    pub intersections: bool,
    /// `tied_sides` is positive.
    pub ties: bool,
    /// `edge_on_triangles` is exactly this.
    pub edge_on: Option<usize>,
}

/// A scene with its exact expected regions.
pub(crate) struct FaceCase {
    pub name: &'static str,
    pub meshes: Vec<FaceFixture>,
    pub view: View,
    pub region: Option<Aabb2>,
    pub limits: VisibleFaceLimits,
    pub regions: Vec<VisibleRegion>,
    pub evidence: EvidenceExpect,
}

/// A scene with the expected total area per source, for coordinates that are
/// not binary64.
pub(crate) struct AreaCase {
    pub name: &'static str,
    pub meshes: Vec<FaceFixture>,
    pub view: View,
    pub limits: VisibleFaceLimits,
    /// `(mesh, source, area)`, in output order. Every other source must have
    /// no region.
    pub areas: Vec<(MeshId, FaceSource, f64)>,
}

/// A request a provider must refuse.
pub(crate) struct FaceRefusal {
    pub name: &'static str,
    pub meshes: Vec<FaceFixture>,
    pub view: View,
    pub limits: VisibleFaceLimits,
    /// Whether the refusal's error is the expected kind.
    pub accepts: fn(&GeomError) -> bool,
    /// The mesh the refusal must name, or `None` for a request-level one.
    pub mesh: Option<MeshId>,
}

/// Relative tolerance of the area cases.
pub(crate) const AREA_TOLERANCE: f64 = 1e-12;

pub(crate) fn limits() -> VisibleFaceLimits {
    VisibleFaceLimits::new(8, 100_000, 100_000, 10_000_000, 100_000, 1_000_000)
}

/// Elevation looking along `+y` from the `-y` side: image `(x, z)`, depth
/// `-y`. It is the plan view rotated about the `x` axis: the plan point
/// `(x, y, z)` appears here as the world point `(x, -z, y)`.
pub(crate) fn elevation_frame() -> Frame3 {
    Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Z,
        z: Vec3::NEG_Y,
    }
}

fn plan() -> View {
    View::new(plan_frame())
}

fn solid(id: u64, mesh: TriMesh) -> FaceFixture {
    FaceFixture::new(id, Solid, mesh)
}

fn sheet(id: u64, mesh: TriMesh) -> FaceFixture {
    FaceFixture::new(id, Sheet, mesh)
}

fn pts(points: &[[f64; 2]]) -> Vec<Point2> {
    points.iter().map(|&[x, y]| Point2::new(x, y)).collect()
}

fn region(
    mesh: u64,
    source: FaceSource,
    outer: &[[f64; 2]],
    holes: &[&[[f64; 2]]],
    triangles: &[u32],
) -> VisibleRegion {
    VisibleRegion::new(
        MeshId(mesh),
        source,
        pts(outer),
        holes.iter().map(|h| pts(h)).collect(),
        triangles.to_vec(),
    )
}

fn group(mesh: u64, group: u32, outer: &[[f64; 2]], triangles: &[u32]) -> VisibleRegion {
    region(
        mesh,
        FaceSource::Group(FaceGroup(group)),
        outer,
        &[],
        triangles,
    )
}

fn case(
    name: &'static str,
    meshes: Vec<FaceFixture>,
    view: View,
    regions: Vec<VisibleRegion>,
) -> FaceCase {
    FaceCase {
        name,
        meshes,
        view,
        region: None,
        limits: limits(),
        regions,
        evidence: EvidenceExpect::default(),
    }
}

/// Front box F over the part of a back box B. In plan, F is `[0,2]x[0,2]`
/// at `z` in `[2,3]` and B is `[1,4]x[1,3]` at `z` in `[0,1]`; the top
/// faces are triangles 2 and 3 of a box. F is whole. B is its top minus the
/// overlap `[1,2]x[1,2]`; the vertices `(2,1)`, `(2,2)` and `(1,2)` are
/// shared bit for bit.
fn two_boxes_regions(triangles: &[u32]) -> Vec<VisibleRegion> {
    vec![
        group(
            1,
            0,
            &[
                [0.0, 0.0],
                [2.0, 0.0],
                [2.0, 1.0],
                [2.0, 2.0],
                [1.0, 2.0],
                [0.0, 2.0],
            ],
            triangles,
        ),
        group(
            2,
            0,
            &[
                [1.0, 2.0],
                [2.0, 2.0],
                [2.0, 1.0],
                [4.0, 1.0],
                [4.0, 3.0],
                [1.0, 3.0],
            ],
            triangles,
        ),
    ]
}

fn two_boxes_overlapping() -> FaceCase {
    case(
        "two_boxes_overlapping",
        vec![
            solid(1, box_mesh([0.0, 0.0, 2.0], [2.0, 2.0, 3.0])),
            solid(2, box_mesh([1.0, 1.0, 0.0], [4.0, 3.0, 1.0])),
        ],
        plan(),
        two_boxes_regions(&[2, 3]),
    )
}

/// A large low box B with a small high box F over its middle: B's top is
/// the whole `[0,4]^2` with a clockwise hole where F stands.
fn box_over_slab_regions() -> Vec<VisibleRegion> {
    vec![
        group(
            1,
            0,
            &[[1.0, 1.0], [2.0, 1.0], [2.0, 2.0], [1.0, 2.0]],
            &[2, 3],
        ),
        region(
            2,
            FaceSource::Group(FaceGroup(0)),
            &[[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]],
            &[&[[1.0, 1.0], [1.0, 2.0], [2.0, 2.0], [2.0, 1.0]]],
            &[2, 3],
        ),
    ]
}

fn box_over_slab(name: &'static str) -> FaceCase {
    case(
        name,
        vec![
            solid(1, box_mesh([1.0, 1.0, 2.0], [2.0, 2.0, 3.0])),
            solid(2, box_mesh([0.0, 0.0, 0.0], [4.0, 4.0, 1.0])),
        ],
        plan(),
        box_over_slab_regions(),
    )
}

/// The first case seen from the front: the plan point `(x, y, z)` is the
/// world point `(x, -z, y)` (a rotation, so winding is kept), which gives the
/// same image and the same depth, and so the same rings. The nearest face of
/// each box is now the one at its smallest `y`, the front face with
/// triangles 4 and 5.
fn two_boxes_elevation() -> FaceCase {
    case(
        "two_boxes_elevation",
        vec![
            solid(1, box_mesh([0.0, -3.0, 0.0], [2.0, -2.0, 2.0])),
            solid(2, box_mesh([1.0, -1.0, 1.0], [4.0, 0.0, 3.0])),
        ],
        View::new(elevation_frame()),
        two_boxes_regions(&[4, 5]),
    )
}

/// A column C `[0,2]^2 x [0,3]` (all group 7) on a slab S `[-1,3]^2 x
/// [-1,0]` (all group 9), with the near plane at 1.5. The column is cut: its
/// cap covers `[0,2]^2` and hides the column's own bottom face and the slab
/// under it. The slab's top is the square with a hole where the cap is.
fn cut_column_on_slab() -> FaceCase {
    case(
        "cut_column_on_slab",
        vec![
            solid(1, box_mesh([0.0, 0.0, 0.0], [2.0, 2.0, 3.0])).uniform(7),
            solid(2, box_mesh([-1.0, -1.0, -1.0], [3.0, 3.0, 0.0])).uniform(9),
        ],
        plan().with_near(1.5),
        vec![
            region(
                1,
                FaceSource::Cap,
                &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]],
                &[],
                &[],
            ),
            region(
                2,
                FaceSource::Group(FaceGroup(9)),
                &[[-1.0, -1.0], [3.0, -1.0], [3.0, 3.0], [-1.0, 3.0]],
                &[&[[0.0, 0.0], [0.0, 2.0], [2.0, 2.0], [2.0, 0.0]]],
                &[2, 3],
            ),
        ],
    )
}

/// The same column as a sheet is never capped: only its bottom face
/// (triangles 0 and 1, seen from inside) lies at depth below the near plane;
/// the top is beyond it and the walls are edge-on.
fn cut_sheet_has_no_cap() -> FaceCase {
    case(
        "cut_sheet_has_no_cap",
        vec![sheet(1, box_mesh([0.0, 0.0, 0.0], [2.0, 2.0, 3.0]))],
        plan().with_near(1.5),
        vec![group(
            1,
            0,
            &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]],
            &[0, 1],
        )],
    )
}

/// A horizontal sheet S at depth 1 over `[0,4]^2` and a ramp R over
/// `[1,3]x[0,4]` whose depth is `y / 2`. They meet where `y / 2 = 1`, that
/// is on `y = 2`: R is in front for `y > 2`. S is the square with a notch.
fn interpenetrating_sheets() -> FaceCase {
    let ramp = TriMesh::new(
        vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(3.0, 4.0, 2.0),
            Point3::new(1.0, 4.0, 2.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );
    let mut c = case(
        "interpenetrating_sheets",
        vec![
            sheet(1, sheet_quad(0.0, 0.0, 4.0, 4.0, 1.0)),
            sheet(2, ramp),
        ],
        plan(),
        vec![
            group(
                1,
                0,
                &[
                    [0.0, 0.0],
                    [4.0, 0.0],
                    [4.0, 4.0],
                    [3.0, 4.0],
                    [3.0, 2.0],
                    [1.0, 2.0],
                    [1.0, 4.0],
                    [0.0, 4.0],
                ],
                &[0, 1],
            ),
            group(
                2,
                0,
                &[[1.0, 2.0], [3.0, 2.0], [3.0, 4.0], [1.0, 4.0]],
                &[0, 1],
            ),
        ],
    );
    c.evidence.intersections = true;
    c
}

/// Two coplanar sheets overlap over `[1,2]^2`. Equal depth never occludes, so
/// both own the overlap and the regions overlap exactly there.
fn coplanar_tie() -> FaceCase {
    let mut c = case(
        "coplanar_tie",
        vec![
            sheet(1, sheet_quad(0.0, 0.0, 2.0, 2.0, 0.0)),
            sheet(2, sheet_quad(1.0, 1.0, 3.0, 3.0, 0.0)),
        ],
        plan(),
        vec![
            group(
                1,
                0,
                &[
                    [0.0, 0.0],
                    [2.0, 0.0],
                    [2.0, 1.0],
                    [2.0, 2.0],
                    [1.0, 2.0],
                    [0.0, 2.0],
                ],
                &[0, 1],
            ),
            group(
                2,
                0,
                &[
                    [1.0, 1.0],
                    [2.0, 1.0],
                    [3.0, 1.0],
                    [3.0, 3.0],
                    [1.0, 3.0],
                    [1.0, 2.0],
                ],
                &[0, 1],
            ),
        ],
    );
    c.evidence.ties = true;
    c
}

/// Two floors side by side and a wall standing on their common edge. In plan
/// the wall is edge-on: it makes no region, never ties and never splits a
/// floor.
fn edge_on_wall_between_floors() -> FaceCase {
    let wall = TriMesh::new(
        vec![
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 4.0, 0.0),
            Point3::new(2.0, 4.0, 3.0),
            Point3::new(2.0, 0.0, 3.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );
    let mut c = case(
        "edge_on_wall_between_floors",
        vec![
            sheet(1, sheet_quad(0.0, 0.0, 2.0, 4.0, 0.0)),
            sheet(2, sheet_quad(2.0, 0.0, 4.0, 4.0, 0.0)),
            sheet(3, wall),
        ],
        plan(),
        vec![
            group(
                1,
                0,
                &[[0.0, 0.0], [2.0, 0.0], [2.0, 4.0], [0.0, 4.0]],
                &[0, 1],
            ),
            group(
                2,
                0,
                &[[2.0, 0.0], [4.0, 0.0], [4.0, 4.0], [2.0, 4.0]],
                &[0, 1],
            ),
        ],
    );
    c.evidence.edge_on = Some(2);
    c
}

/// A box `[0,2]x[0,2]x[0,1]` whose top triangles 2 and 3 are in different
/// groups. Both are visible and meet along the diagonal.
fn two_groups_on_one_face() -> FaceCase {
    let mesh = box_mesh([0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
    let mut fixture = solid(1, mesh);
    fixture.groups = vec![FaceGroup(0); 12];
    fixture.groups[2] = FaceGroup(1);
    fixture.groups[3] = FaceGroup(2);
    case(
        "two_groups_on_one_face",
        vec![fixture],
        plan(),
        vec![
            group(1, 1, &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0]], &[2]),
            group(1, 2, &[[0.0, 0.0], [2.0, 2.0], [0.0, 2.0]], &[3]),
        ],
    )
}

/// Two unit squares side by side, as a sheet of four triangles in one plane,
/// requested with planar groups: one group, one region without the interior
/// vertices `(1,0)` and `(1,1)`.
fn flat_sheet_is_one_planar_group() -> FaceCase {
    let mesh = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
        ],
        vec![0, 1, 4, 0, 4, 3, 1, 2, 5, 1, 5, 4],
    );
    case(
        "flat_sheet_is_one_planar_group",
        vec![sheet(1, mesh).planar(&[0; 4])],
        plan(),
        vec![group(
            1,
            0,
            &[[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]],
            &[0, 1, 2, 3],
        )],
    )
}

/// A fold: triangle 0 faces `+z` and triangle 1 faces `-z` over the same
/// edge, so their normals are parallel but not same-facing and planar groups
/// keep them apart. Seen from above they tie in depth over triangle 1, so
/// each keeps its area; the vertex `(0.5, 0.5)` of triangle 1 stays on the
/// hypotenuse of triangle 0.
fn fold_is_two_planar_groups() -> FaceCase {
    let mesh = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.5, 0.5, 0.0),
        ],
        vec![0, 1, 2, 1, 0, 3],
    );
    case(
        "fold_is_two_planar_groups",
        vec![sheet(1, mesh).planar(&[0, 1])],
        plan(),
        vec![
            group(
                1,
                0,
                &[[0.0, 0.0], [1.0, 0.0], [0.5, 0.5], [0.0, 1.0]],
                &[0],
            ),
            group(1, 1, &[[0.0, 0.0], [1.0, 0.0], [0.5, 0.5]], &[1]),
        ],
    )
}

/// The box-over-slab scene cropped to `[1.5,3]^2`. F is cut to its corner
/// square `[1.5,2]^2` and B's top to the L around it.
fn crop() -> FaceCase {
    let mut c = box_over_slab("crop");
    c.region = Some(Aabb2 {
        min: Point2::new(1.5, 1.5),
        max: Point2::new(3.0, 3.0),
    });
    c.regions = vec![
        group(
            1,
            0,
            &[[1.5, 1.5], [2.0, 1.5], [2.0, 2.0], [1.5, 2.0]],
            &[2, 3],
        ),
        group(
            2,
            0,
            &[
                [1.5, 2.0],
                [2.0, 2.0],
                [2.0, 1.5],
                [3.0, 1.5],
                [3.0, 3.0],
                [1.5, 3.0],
            ],
            &[2, 3],
        ),
    ];
    c
}

/// The ramp of the interpenetration case (depth `y / 2` over `[1,3]x[0,4]`)
/// with `far = 1`: only `y >= 2` is kept, the half that reaches depth 1.
fn far_clip() -> FaceCase {
    let ramp = TriMesh::new(
        vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(3.0, 4.0, 2.0),
            Point3::new(1.0, 4.0, 2.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );
    case(
        "far_clip",
        vec![sheet(1, ramp)],
        plan().with_far(1.0),
        vec![group(
            1,
            0,
            &[[1.0, 2.0], [3.0, 2.0], [3.0, 4.0], [1.0, 4.0]],
            &[0, 1],
        )],
    )
}

/// One sheet in two components: the annulus `[0,4]^2 \ [1,3]^2`
/// (triangles 0..8, vertices 0..8) and the square `[1,3]^2` with its own
/// vertices (triangles 8 and 9). They share a boundary in space but no
/// vertex. The visible region is one square without a hole.
fn floating_cluster() -> FaceCase {
    let positions = [
        [0.0, 0.0],
        [4.0, 0.0],
        [4.0, 4.0],
        [0.0, 4.0],
        [1.0, 1.0],
        [3.0, 1.0],
        [3.0, 3.0],
        [1.0, 3.0],
        [1.0, 1.0],
        [3.0, 1.0],
        [3.0, 3.0],
        [1.0, 3.0],
    ]
    .map(|[x, y]| Point3::new(x, y, 0.0))
    .to_vec();
    let indices = vec![
        0, 1, 5, 0, 5, 4, // bottom strip
        1, 2, 6, 1, 6, 5, // right strip
        2, 3, 7, 2, 7, 6, // top strip
        3, 0, 4, 3, 4, 7, // left strip
        8, 9, 10, 8, 10, 11, // inner square
    ];
    case(
        "floating_cluster",
        vec![sheet(1, TriMesh::new(positions, indices))],
        plan(),
        vec![group(
            1,
            0,
            &[[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]],
            &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        )],
    )
}

pub(crate) fn cases() -> Vec<FaceCase> {
    vec![
        two_boxes_overlapping(),
        box_over_slab("box_over_slab_hole"),
        two_boxes_elevation(),
        cut_column_on_slab(),
        cut_sheet_has_no_cap(),
        interpenetrating_sheets(),
        coplanar_tie(),
        edge_on_wall_between_floors(),
        two_groups_on_one_face(),
        flat_sheet_is_one_planar_group(),
        fold_is_two_planar_groups(),
        crop(),
        far_clip(),
        floating_cluster(),
    ]
}

/// The isometric unit cube shows three faces, each a rhombus of area
/// `1 / sqrt(3)` (a unit face seen at `cos = 1/sqrt(3)`). The six faces of
/// the cube are six groups: bottom 0, top 1, front 2, right 3, back 4 and
/// left 5; the top, the right and the back face the viewer.
pub(crate) fn area_cases() -> Vec<AreaCase> {
    let mut cube = solid(1, box_mesh([0.0; 3], [1.0; 3]));
    cube.groups = (0..6u32)
        .flat_map(|face| [FaceGroup(face), FaceGroup(face)])
        .collect();
    // The same cube with planar groups: each face is named by its smaller
    // triangle, so face k is group 2k.
    let planar =
        solid(1, box_mesh([0.0; 3], [1.0; 3])).planar(&[0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10]);
    let face = 1.0 / 3.0_f64.sqrt();
    vec![
        AreaCase {
            name: "isometric_cube",
            meshes: vec![cube],
            view: View::new(isometric_frame(Point3::ZERO)),
            limits: limits(),
            areas: vec![
                (MeshId(1), FaceSource::Group(FaceGroup(1)), face),
                (MeshId(1), FaceSource::Group(FaceGroup(3)), face),
                (MeshId(1), FaceSource::Group(FaceGroup(4)), face),
            ],
        },
        AreaCase {
            name: "isometric_cube_planar",
            meshes: vec![planar],
            view: View::new(isometric_frame(Point3::ZERO)),
            limits: limits(),
            areas: vec![
                (MeshId(1), FaceSource::Group(FaceGroup(2)), face),
                (MeshId(1), FaceSource::Group(FaceGroup(6)), face),
                (MeshId(1), FaceSource::Group(FaceGroup(8)), face),
            ],
        },
    ]
}

fn refusal_of(
    name: &'static str,
    meshes: Vec<FaceFixture>,
    limits: VisibleFaceLimits,
    accepts: fn(&GeomError) -> bool,
    mesh: Option<MeshId>,
) -> FaceRefusal {
    FaceRefusal {
        name,
        meshes,
        view: plan(),
        limits,
        accepts,
        mesh,
    }
}

fn clone_fixture(f: &crate::fixtures::Fixture) -> FaceFixture {
    let mut clone = FaceFixture::new(f.id.0, f.role, f.mesh.clone());
    clone.groups = f.groups.clone();
    clone.planar = f.planar;
    clone.occlusion = f.occlusion.clone();
    clone
}

pub(crate) fn face_refusals() -> Vec<FaceRefusal> {
    // The input refusals of hidden lines, which name their mesh.
    let mut all: Vec<FaceRefusal> = refusals()
        .into_iter()
        .filter(|r| r.mesh.is_some())
        .map(|r| FaceRefusal {
            name: r.name,
            meshes: r.meshes.iter().map(clone_fixture).collect(),
            view: r.view,
            limits: limits(),
            accepts: r.accepts,
            mesh: r.mesh,
        })
        .collect();

    // Four group ids for a box of twelve triangles.
    let mut wrong = solid(2, box_mesh([0.0; 3], [1.0; 3]));
    wrong.groups = vec![FaceGroup(0); 4];
    all.push(refusal_of(
        "group count mismatch",
        vec![solid(1, box_mesh([5.0; 3], [6.0; 3])), wrong],
        limits(),
        |e| matches!(e, GeomError::InvalidInput(_)),
        Some(MeshId(2)),
    ));

    let two_boxes = || two_boxes_overlapping().meshes;
    let mut regions = limits();
    regions.max_regions = 1;
    all.push(refusal_of(
        "region limit",
        two_boxes(),
        regions,
        |e| {
            matches!(
                e,
                GeomError::BudgetExceeded {
                    resource: "visible-face regions"
                }
            )
        },
        None,
    ));
    let mut vertices = limits();
    vertices.max_region_vertices = 5;
    all.push(refusal_of(
        "region vertex limit",
        two_boxes(),
        vertices,
        |e| {
            matches!(
                e,
                GeomError::BudgetExceeded {
                    resource: "visible-face region vertices"
                }
            )
        },
        None,
    ));
    let mut pairs = limits();
    pairs.max_candidate_pairs = 1;
    all.push(refusal_of(
        "candidate limit",
        two_boxes(),
        pairs,
        |e| {
            matches!(
                e,
                GeomError::BudgetExceeded {
                    resource: "visible-face candidate pairs"
                }
            )
        },
        None,
    ));
    all
}

// --- see-through triangles (ADR 0111) ----------------------------------------

const G0: FaceSource = FaceSource::Group(FaceGroup(0));

/// `region` with `behind` entries `(mesh, source, layer)`.
fn seen(region: VisibleRegion, behind: &[(u64, FaceSource, u32)]) -> VisibleRegion {
    region.with_behind(
        behind
            .iter()
            .map(|&(mesh, source, layer)| BehindSurface::new(MeshId(mesh), source, layer))
            .collect(),
    )
}

fn see_through_case(
    name: &'static str,
    meshes: Vec<FaceFixture>,
    view: View,
    regions: Vec<VisibleRegion>,
) -> FaceCase {
    let mut c = case(name, meshes, view, regions);
    c.limits = limits().with_max_behind_surfaces(16);
    c
}

/// The square `[0,4]^2` at `z = 0`.
fn floor() -> FaceFixture {
    sheet(1, sheet_quad(0.0, 0.0, 4.0, 4.0, 0.0))
}

/// The square `[0,4]^2` with the hole `[1,3]^2`.
fn square_with_hole(mesh: u64, triangles: &[u32]) -> VisibleRegion {
    region(
        mesh,
        G0,
        &[[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]],
        &[&[[1.0, 1.0], [1.0, 3.0], [3.0, 3.0], [3.0, 1.0]]],
        triangles,
    )
}

fn middle(mesh: u64, source: FaceSource, triangles: &[u32]) -> VisibleRegion {
    region(
        mesh,
        source,
        &[[1.0, 1.0], [3.0, 1.0], [3.0, 3.0], [1.0, 3.0]],
        &[],
        triangles,
    )
}

/// A pane over a box. The box's top `[0,4]^2` at `z = 1` (triangles 2 and
/// 3) is frontmost outside the see-through pane `[1,3] x [1,5]` at `z = 2`:
/// a square with the notch `[1,3] x [1,4]`, area 10. The pane is frontmost
/// over all of its image: over the box it lists the box's top behind it,
/// beyond the box (`y` in `[4,5]`) nothing. The pane's diagonal from `(1,1)`
/// to `(3,5)` separates no class, and both of its triangles reach both
/// parts.
fn pane_over_box() -> FaceCase {
    see_through_case(
        "pane over a box",
        vec![
            solid(1, box_mesh([0.0, 0.0, 0.0], [4.0, 4.0, 1.0])),
            sheet(2, sheet_quad(1.0, 1.0, 3.0, 5.0, 2.0)).see_through(),
        ],
        plan(),
        vec![
            group(
                1,
                0,
                &[
                    [0.0, 0.0],
                    [4.0, 0.0],
                    [4.0, 4.0],
                    [3.0, 4.0],
                    [3.0, 1.0],
                    [1.0, 1.0],
                    [1.0, 4.0],
                    [0.0, 4.0],
                ],
                &[2, 3],
            ),
            seen(
                group(
                    2,
                    0,
                    &[[1.0, 1.0], [3.0, 1.0], [3.0, 4.0], [1.0, 4.0]],
                    &[0, 1],
                ),
                &[(1, G0, 1)],
            ),
            group(
                2,
                0,
                &[[1.0, 4.0], [3.0, 4.0], [3.0, 5.0], [1.0, 5.0]],
                &[0, 1],
            ),
        ],
    )
}

/// Two panes over the box of `pane_over_box`: A = `[1,3]^2` at `z = 2` and
/// B = `[2,5]^2` at `z = 3`, both see-through. The box's top is the square
/// minus `A ∪ B`. A is frontmost where B is not, the L `A \ B`, over the
/// box. B is frontmost over all of its image: over `[2,3]^2` it sees A and
/// then the box, over the rest of the box only the box, beyond the box
/// nothing. Every corner where three class boundaries meet is a ring vertex.
fn panes_over_a_box() -> FaceCase {
    see_through_case(
        "panes over a box",
        vec![
            solid(1, box_mesh([0.0, 0.0, 0.0], [4.0, 4.0, 1.0])),
            sheet(2, sheet_quad(1.0, 1.0, 3.0, 3.0, 2.0)).see_through(),
            sheet(3, sheet_quad(2.0, 2.0, 5.0, 5.0, 3.0)).see_through(),
        ],
        plan(),
        vec![
            group(
                1,
                0,
                &[
                    [0.0, 0.0],
                    [4.0, 0.0],
                    [4.0, 2.0],
                    [3.0, 2.0],
                    [3.0, 1.0],
                    [1.0, 1.0],
                    [1.0, 3.0],
                    [2.0, 3.0],
                    [2.0, 4.0],
                    [0.0, 4.0],
                ],
                &[2, 3],
            ),
            seen(
                group(
                    2,
                    0,
                    &[
                        [1.0, 1.0],
                        [3.0, 1.0],
                        [3.0, 2.0],
                        [2.0, 2.0],
                        [2.0, 3.0],
                        [1.0, 3.0],
                    ],
                    &[0, 1],
                ),
                &[(1, G0, 1)],
            ),
            seen(
                group(
                    3,
                    0,
                    &[[2.0, 2.0], [3.0, 2.0], [3.0, 3.0], [2.0, 3.0]],
                    &[0, 1],
                ),
                &[(2, G0, 1), (1, G0, 2)],
            ),
            seen(
                group(
                    3,
                    0,
                    &[
                        [2.0, 3.0],
                        [3.0, 3.0],
                        [3.0, 2.0],
                        [4.0, 2.0],
                        [4.0, 4.0],
                        [2.0, 4.0],
                    ],
                    &[0, 1],
                ),
                &[(1, G0, 1)],
            ),
            group(
                3,
                0,
                &[
                    [2.0, 4.0],
                    [4.0, 4.0],
                    [4.0, 2.0],
                    [5.0, 2.0],
                    [5.0, 5.0],
                    [2.0, 5.0],
                ],
                &[0, 1],
            ),
        ],
    )
}

/// Three see-through sheets in a depth cycle, each in front of the next
/// where two overlap and no three overlapping: A, the strip `[0,4] x [0,1]`
/// with `z = x`; B, the strip `[3,4] x [0,4]` with `z = 2y`; C, the band
/// between `y = x - 1` and `y = x` for `y` in `[0,3]` at `z = 3`. A is in
/// front of B over `[3,4] x [0,1]` (`x >= 3 > 2 >= 2y`), B of C over the
/// triangle `(3,2) (4,3) (3,3)` (`2y >= 4 > 3`) and C of A over the
/// parallelogram `(0,0) (1,0) (2,1) (1,1)` (`3 > 2 >= x`). A's diagonal
/// `y = x / 4` keeps its triangle 0 out of the corner left of C.
fn see_through_cycle() -> FaceCase {
    let ramp = |id: u64, corners: [[f64; 3]; 4]| {
        let positions = corners.map(|[x, y, z]| Point3::new(x, y, z)).to_vec();
        sheet(id, TriMesh::new(positions, vec![0, 1, 2, 0, 2, 3])).see_through()
    };
    see_through_case(
        "see-through cycle",
        vec![
            ramp(
                1,
                [
                    [0.0, 0.0, 0.0],
                    [4.0, 0.0, 4.0],
                    [4.0, 1.0, 4.0],
                    [0.0, 1.0, 0.0],
                ],
            ),
            ramp(
                2,
                [
                    [3.0, 0.0, 0.0],
                    [4.0, 0.0, 0.0],
                    [4.0, 4.0, 8.0],
                    [3.0, 4.0, 8.0],
                ],
            ),
            ramp(
                3,
                [
                    [0.0, 0.0, 3.0],
                    [1.0, 0.0, 3.0],
                    [4.0, 3.0, 3.0],
                    [3.0, 3.0, 3.0],
                ],
            ),
        ],
        plan(),
        vec![
            group(1, 0, &[[0.0, 0.0], [1.0, 1.0], [0.0, 1.0]], &[1]),
            group(
                1,
                0,
                &[[1.0, 0.0], [3.0, 0.0], [3.0, 1.0], [2.0, 1.0]],
                &[0, 1],
            ),
            seen(
                group(
                    1,
                    0,
                    &[[3.0, 0.0], [4.0, 0.0], [4.0, 1.0], [3.0, 1.0]],
                    &[0, 1],
                ),
                &[(2, G0, 1)],
            ),
            group(
                2,
                0,
                &[[3.0, 1.0], [4.0, 1.0], [4.0, 3.0], [3.0, 2.0]],
                &[0, 1],
            ),
            seen(
                group(2, 0, &[[3.0, 2.0], [4.0, 3.0], [3.0, 3.0]], &[0, 1]),
                &[(3, G0, 1)],
            ),
            group(
                2,
                0,
                &[[3.0, 3.0], [4.0, 3.0], [4.0, 4.0], [3.0, 4.0]],
                &[0, 1],
            ),
            seen(
                group(
                    3,
                    0,
                    &[[0.0, 0.0], [1.0, 0.0], [2.0, 1.0], [1.0, 1.0]],
                    &[0, 1],
                ),
                &[(1, G0, 1)],
            ),
            group(
                3,
                0,
                &[[1.0, 1.0], [2.0, 1.0], [3.0, 2.0], [3.0, 3.0]],
                &[0, 1],
            ),
        ],
    )
}

/// A see-through box `[1,3]^2 x [1,2]` over the floor. Over its image the
/// top (triangles 2 and 3) is frontmost; behind it lie its own bottom, the
/// same source at another depth, and then the floor, which ends the stack.
/// Its sides are seen edge-on.
fn see_through_box() -> FaceCase {
    let mut c = see_through_case(
        "see-through box over a floor",
        vec![
            floor(),
            solid(2, box_mesh([1.0, 1.0, 1.0], [3.0, 3.0, 2.0])).see_through(),
        ],
        plan(),
        vec![
            square_with_hole(1, &[0, 1]),
            seen(middle(2, G0, &[2, 3]), &[(2, G0, 1), (1, G0, 2)]),
        ],
    );
    c.evidence.edge_on = Some(8);
    c
}

/// A see-through pane in the floor's plane. The tie gives both the area;
/// the floor occludes, so it lists nothing, and the pane lists the floor in
/// its own layer 0.
fn see_through_tie() -> FaceCase {
    let mut c = see_through_case(
        "see-through pane tied with a floor",
        vec![
            floor(),
            sheet(2, sheet_quad(1.0, 1.0, 3.0, 3.0, 0.0)).see_through(),
        ],
        plan(),
        vec![
            group(
                1,
                0,
                &[[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]],
                &[0, 1],
            ),
            seen(middle(2, G0, &[0, 1]), &[(1, G0, 0)]),
        ],
    );
    c.evidence.ties = true;
    c
}

/// A see-through sheet over the whole floor with an opaque sheet
/// `[1,3]^2` between them: the stack stops at the first opaque surface, so
/// the middle lists the opaque sheet only and the rest the floor only. The
/// floor and the opaque sheet are frontmost nowhere.
fn stack_stops_at_an_occluder() -> FaceCase {
    see_through_case(
        "a stack stops at an opaque sheet",
        vec![
            floor(),
            sheet(2, sheet_quad(1.0, 1.0, 3.0, 3.0, 1.0)),
            sheet(3, sheet_quad(0.0, 0.0, 4.0, 4.0, 2.0)).see_through(),
        ],
        plan(),
        vec![
            seen(square_with_hole(3, &[0, 1]), &[(1, G0, 1)]),
            seen(middle(3, G0, &[0, 1]), &[(2, G0, 1)]),
        ],
    )
}

/// A see-through box `[1,3]^2 x [1,4]` cut by the near plane `z = 2`: its
/// cap is decided from all of its triangles and occludes, so the cap region
/// lists nothing behind it.
fn cap_of_a_see_through_solid() -> FaceCase {
    see_through_case(
        "cap of a see-through solid",
        vec![
            floor(),
            solid(2, box_mesh([1.0, 1.0, 1.0], [3.0, 3.0, 4.0])).see_through(),
        ],
        plan().with_near(2.0),
        vec![
            square_with_hole(1, &[0, 1]),
            middle(2, FaceSource::Cap, &[]),
        ],
    )
}

/// One sheet at `z = 1` over the floor `[0,4] x [0,2]`: the square
/// `[0,2]^2` (triangles 0 and 1) opaque, the square `[2,4] x [0,2]`
/// (triangles 2 and 3) see-through. Planar groups do not join triangles of
/// different occlusion, so the coplanar sheet is two groups, named 0 and 2.
fn planar_groups_split_by_occlusion() -> FaceFixture {
    let positions = [
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 2.0],
        [0.0, 2.0],
        [4.0, 0.0],
        [4.0, 2.0],
    ]
    .map(|[x, y]| Point3::new(x, y, 1.0))
    .to_vec();
    let mut glazed = sheet(
        2,
        TriMesh::new(positions, vec![0, 1, 2, 0, 2, 3, 1, 4, 5, 1, 5, 2]),
    );
    glazed.occlusion = vec![
        Occlusion::Opaque,
        Occlusion::Opaque,
        Occlusion::SeeThrough,
        Occlusion::SeeThrough,
    ];
    glazed
}

fn planar_glazing() -> FaceCase {
    see_through_case(
        "planar groups split by occlusion",
        vec![
            sheet(1, sheet_quad(0.0, 0.0, 4.0, 2.0, 0.0)),
            planar_groups_split_by_occlusion().planar(&[0, 0, 2, 2]),
        ],
        plan(),
        vec![
            group(
                2,
                0,
                &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]],
                &[0, 1],
            ),
            seen(
                group(
                    2,
                    2,
                    &[[2.0, 0.0], [4.0, 0.0], [4.0, 2.0], [2.0, 2.0]],
                    &[2, 3],
                ),
                &[(1, G0, 1)],
            ),
        ],
    )
}

/// The regions a provider that claims `Occlusion::SeeThrough` must
/// reproduce, in plan view; each is also run with every triangle reversed.
pub(crate) fn see_through_face_cases() -> Vec<FaceCase> {
    vec![
        pane_over_box(),
        panes_over_a_box(),
        see_through_cycle(),
        see_through_box(),
        see_through_tie(),
        stack_stops_at_an_occluder(),
        cap_of_a_see_through_solid(),
        planar_glazing(),
    ]
}

/// The base cases a see-through claim must answer unchanged when every
/// triangle is flagged opaque.
pub(crate) const NEUTRAL_FACE_CASES: [&str; 3] = [
    "two_boxes_overlapping",
    "cut_column_on_slab",
    "coplanar_tie",
];

/// The see-through requests a provider that claims them must refuse.
pub(crate) fn see_through_face_refusals() -> Vec<FaceRefusal> {
    let invalid = |e: &GeomError| matches!(e, GeomError::InvalidInput(_));
    let mut short = pane_over_box().meshes;
    short[1].occlusion.truncate(1);
    vec![
        refusal_of(
            "a face group mixing see-through and opaque triangles",
            vec![
                sheet(1, sheet_quad(0.0, 0.0, 4.0, 2.0, 0.0)),
                planar_groups_split_by_occlusion(),
            ],
            limits().with_max_behind_surfaces(16),
            invalid,
            Some(MeshId(2)),
        ),
        refusal_of(
            "occlusion flags shorter than the triangles",
            short,
            limits().with_max_behind_surfaces(16),
            invalid,
            Some(MeshId(2)),
        ),
        refusal_of(
            "behind surfaces over the limit",
            pane_over_box().meshes,
            limits(),
            |e| {
                matches!(
                    e,
                    GeomError::BudgetExceeded {
                        resource: "visible-face behind surfaces"
                    }
                )
            },
            None,
        ),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn area(ring: &[Point2]) -> f64 {
        let n = ring.len();
        (0..n)
            .map(|i| {
                let (a, b) = (ring[i], ring[(i + 1) % n]);
                a.x * b.y - b.x * a.y
            })
            .sum::<f64>()
            / 2.0
    }

    fn smallest_first(ring: &[Point2]) -> bool {
        let first = (ring[0].x, ring[0].y);
        ring.iter().all(|p| (p.x, p.y) >= first)
    }

    #[test]
    fn expected_rings_have_the_documented_orientation_and_start() {
        for case in cases().into_iter().chain(see_through_face_cases()) {
            for r in &case.regions {
                assert!(r.outer.len() >= 3, "{}", case.name);
                assert!(area(&r.outer) > 0.0, "{}: outer is clockwise", case.name);
                assert!(smallest_first(&r.outer), "{}: outer start", case.name);
                for hole in &r.holes {
                    assert!(hole.len() >= 3, "{}", case.name);
                    assert!(area(hole) < 0.0, "{}: hole is counter-clockwise", case.name);
                    assert!(smallest_first(hole), "{}: hole start", case.name);
                }
                assert!(r.triangles.windows(2).all(|w| w[0] < w[1]), "{}", case.name);
            }
        }
    }

    #[test]
    fn expected_regions_are_sorted_by_mesh_position_then_source() {
        for case in cases().into_iter().chain(see_through_face_cases()) {
            let position = |id: MeshId| case.meshes.iter().position(|m| m.id == id).unwrap();
            let keys: Vec<_> = case
                .regions
                .iter()
                .map(|r| (position(r.mesh), r.source))
                .collect();
            assert!(keys.windows(2).all(|w| w[0] <= w[1]), "{}", case.name);
        }
    }

    #[test]
    fn expected_areas_of_the_exact_cases() {
        let total = |name: &str, mesh: u64| -> f64 {
            let case = cases().into_iter().find(|c| c.name == name).unwrap();
            case.regions
                .iter()
                .filter(|r| r.mesh == MeshId(mesh))
                .map(|r| area(&r.outer) + r.holes.iter().map(|h| area(h)).sum::<f64>())
                .sum()
        };
        // F is the square of side 2; B's top is 3x2 minus the 1x1 overlap.
        assert_eq!(total("two_boxes_overlapping", 1), 4.0);
        assert_eq!(total("two_boxes_overlapping", 2), 5.0);
        assert_eq!(total("two_boxes_elevation", 1), 4.0);
        assert_eq!(total("two_boxes_elevation", 2), 5.0);
        // B's top is 4x4 minus the 1x1 hole.
        assert_eq!(total("box_over_slab_hole", 2), 15.0);
        assert_eq!(total("crop", 1), 0.25);
        assert_eq!(total("crop", 2), 2.25 - 0.25);
        // The slab is 4x4 minus the 2x2 hole; the cap is the hole.
        assert_eq!(total("cut_column_on_slab", 2), 12.0);
        assert_eq!(total("cut_column_on_slab", 1), 4.0);
        // The square with a notch of 2x2, and the half of the ramp.
        assert_eq!(total("interpenetrating_sheets", 1), 16.0 - 4.0);
        assert_eq!(total("interpenetrating_sheets", 2), 4.0);
        assert_eq!(total("far_clip", 1), 4.0);
        assert_eq!(total("floating_cluster", 1), 16.0);
        // The overlap [1,2]^2 is counted twice.
        assert_eq!(total("coplanar_tie", 1), 4.0);
        assert_eq!(total("coplanar_tie", 2), 4.0);
    }

    #[test]
    fn every_scene_mesh_has_whole_triangles_and_matching_groups() {
        let check = |name: &str, meshes: &[FaceFixture]| {
            for m in meshes {
                assert_eq!(m.mesh.indices.len() % 3, 0, "{name}");
                let n = m.mesh.positions.len() as u32;
                assert!(m.mesh.indices.iter().all(|&i| i < n), "{name}");
                assert!(
                    m.groups.is_empty() || m.groups.len() == m.mesh.triangle_count(),
                    "{name}: group count"
                );
            }
        };
        for case in cases().into_iter().chain(see_through_face_cases()) {
            check(case.name, &case.meshes);
        }
        for case in area_cases() {
            check(case.name, &case.meshes);
        }
    }

    #[test]
    fn expected_triangles_belong_to_the_region_group_and_exist() {
        for case in cases().into_iter().chain(see_through_face_cases()) {
            for r in &case.regions {
                let fixture = case.meshes.iter().find(|m| m.id == r.mesh).unwrap();
                for &t in &r.triangles {
                    assert!(
                        (t as usize) < fixture.mesh.triangle_count(),
                        "{}",
                        case.name
                    );
                    assert_eq!(
                        FaceSource::Group(fixture.group(t as usize)),
                        r.source,
                        "{}",
                        case.name
                    );
                }
                assert!(
                    r.source != FaceSource::Cap || r.triangles.is_empty(),
                    "{}: a cap has no triangles",
                    case.name
                );
            }
        }
    }

    #[test]
    fn the_elevation_scene_is_the_plan_scene_rotated() {
        // The plan point (x, y, z) is the world point (x, -z, y): image
        // (x, y) = (x, world z) and depth z = -world y. So each elevation
        // box is the plan box with its y range and z range exchanged.
        let bounds = |m: &FaceFixture| {
            let mut min = Point3::splat(f64::INFINITY);
            let mut max = Point3::splat(f64::NEG_INFINITY);
            for p in &m.mesh.positions {
                min = min.min(*p);
                max = max.max(*p);
            }
            (min, max)
        };
        let plan = two_boxes_overlapping();
        let elevation = two_boxes_elevation();
        for (p, e) in plan.meshes.iter().zip(&elevation.meshes) {
            let ((pmin, pmax), (emin, emax)) = (bounds(p), bounds(e));
            assert_eq!(emin, Point3::new(pmin.x, -pmax.z, pmin.y));
            assert_eq!(emax, Point3::new(pmax.x, -pmin.z, pmax.y));
        }
        for (p, e) in plan.regions.iter().zip(&elevation.regions) {
            assert_eq!((&p.outer, &p.holes), (&e.outer, &e.holes));
        }
    }

    #[test]
    fn the_isometric_cube_names_six_faces_and_three_visible() {
        let case = &area_cases()[0];
        assert_eq!(case.meshes[0].groups.len(), 12);
        assert_eq!(case.areas.len(), 3);
        // Top (z = 1) is the box's triangles 2 and 3.
        assert_eq!(case.meshes[0].group(2), FaceGroup(1));
        assert_eq!(case.meshes[0].group(3), FaceGroup(1));
        // Right (x = 1) is 6 and 7, back (y = 1) is 8 and 9.
        assert_eq!(case.meshes[0].group(6), FaceGroup(3));
        assert_eq!(case.meshes[0].group(8), FaceGroup(4));
    }

    #[test]
    fn the_planar_cube_relabels_each_face_by_its_smaller_triangle() {
        let planar = &area_cases()[1].meshes[0];
        let named = &area_cases()[0].meshes[0];
        assert!(planar.planar && !named.planar);
        for t in 0..12 {
            assert_eq!(planar.group(t).0, 2 * named.group(t).0, "triangle {t}");
        }
    }

    #[test]
    fn refusals_name_their_mesh_where_one_is_at_fault() {
        let refusals = face_refusals();
        assert_eq!(refusals.len(), 8);
        for r in &refusals {
            if let Some(id) = r.mesh {
                assert!(r.meshes.iter().any(|m| m.id == id), "{}", r.name);
            }
        }
        let mismatch = refusals
            .iter()
            .find(|r| r.name == "group count mismatch")
            .unwrap();
        let offender = mismatch.meshes.iter().find(|m| m.id == MeshId(2)).unwrap();
        assert_ne!(offender.groups.len(), offender.mesh.triangle_count());
    }

    #[test]
    fn see_through_cases_cover_their_frontmost_areas_and_stack_by_layer() {
        let total = |name: &str, mesh: u64| -> f64 {
            let case = see_through_face_cases()
                .into_iter()
                .find(|c| c.name == name)
                .unwrap();
            case.regions
                .iter()
                .filter(|r| r.mesh == MeshId(mesh))
                .map(|r| area(&r.outer) + r.holes.iter().map(|h| area(h)).sum::<f64>())
                .sum()
        };
        // The box's top 16 minus the pane over it (2 x 3); the pane 2 x 4.
        assert_eq!(total("pane over a box", 1), 10.0);
        assert_eq!(total("pane over a box", 2), 8.0);
        // The box's top minus A ∪ B (4 + 4 - 1); A minus B; all of B.
        assert_eq!(total("panes over a box", 1), 9.0);
        assert_eq!(total("panes over a box", 2), 3.0);
        assert_eq!(total("panes over a box", 3), 9.0);
        // A minus the parallelogram of C (area 1); B minus [3,4]x[0,1]; C
        // (area 3) minus B's triangle (area 1/2).
        assert_eq!(total("see-through cycle", 1), 3.0);
        assert_eq!(total("see-through cycle", 2), 3.0);
        assert_eq!(total("see-through cycle", 3), 2.5);
        assert_eq!(total("see-through pane tied with a floor", 1), 16.0);
        assert_eq!(total("a stack stops at an opaque sheet", 3), 16.0);
        assert_eq!(total("planar groups split by occlusion", 2), 8.0);

        for case in see_through_face_cases() {
            assert!(
                case.meshes
                    .iter()
                    .any(|m| m.occlusion.contains(&Occlusion::SeeThrough)),
                "{}",
                case.name
            );
            assert!(case.limits.max_behind_surfaces > 0, "{}", case.name);
            for m in &case.meshes {
                assert!(
                    m.occlusion.is_empty() || m.occlusion.len() == m.mesh.triangle_count(),
                    "{}",
                    case.name
                );
            }
            for r in &case.regions {
                let position = |id: MeshId| case.meshes.iter().position(|m| m.id == id).unwrap();
                let keys: Vec<_> = r
                    .behind
                    .iter()
                    .map(|b| (b.layer, position(b.mesh), b.source))
                    .collect();
                assert!(keys.windows(2).all(|w| w[0] < w[1]), "{}", case.name);
                assert!(
                    keys.first().is_none_or(|k| k.0 <= 1)
                        && keys.windows(2).all(|w| w[1].0 - w[0].0 <= 1),
                    "{}",
                    case.name
                );
                let own = case.meshes.iter().find(|m| m.id == r.mesh).unwrap();
                let see_through = r.triangles.first().is_some_and(|&t| {
                    own.occlusion.get(t as usize) == Some(&Occlusion::SeeThrough)
                });
                assert!(r.behind.is_empty() || see_through, "{}", case.name);
            }
        }
        for refusal in see_through_face_refusals() {
            if let Some(id) = refusal.mesh {
                assert_eq!(refusal.meshes[1].id, id, "{}", refusal.name);
            }
        }
        for name in NEUTRAL_FACE_CASES {
            assert!(cases().iter().any(|c| c.name == name), "{name}");
        }
    }
}
