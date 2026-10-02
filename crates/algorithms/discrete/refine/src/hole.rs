//! Hole filling: triangulate a boundary loop, refine the patch to the
//! density of its surroundings, and fair it.
//!
//! This is Liepa's pipeline (*Filling Holes in Meshes*, SGP 2003), the one
//! behind CGAL's `triangulate_refine_and_fair_hole`:
//!
//! 1. **Triangulate.** Among all triangulations of the boundary polygon that
//!    use only its own vertices, dynamic programming picks one whose largest
//!    dihedral angle (between neighbouring patch triangles and between the
//!    patch and the faces around the hole) is the minimum over all of them,
//!    and among those, by a greedy tie-break, a small total area. Liepa's
//!    and CGAL's `O(n^3)` recurrence is only a heuristic for that minimum;
//!    here each sub-polygon's state also names the triangle above it, which
//!    makes the angle exact at `O(n^4)` time and `O(n^3)` memory for `n`
//!    boundary vertices. A triangle with (numerically) collinear corners,
//!    and a diagonal joining two boundary vertices that the mesh already
//!    joins, are never used: either would break the manifold.
//! 2. **Refine.** Each boundary vertex gets a scale `sigma`, the mean
//!    length of its edges. A patch triangle is split at its centroid while
//!    `density * |centroid - corner|` exceeds both the centroid's
//!    interpolated scale and the corner's, for every corner; after each
//!    round of splits, interior patch edges are flipped wherever the two
//!    opposite angles sum to more than pi (the Delaunay criterion).
//! 3. **Fair.** The new interior vertices are placed by [`crate::fair()`]
//!    with the border fixed, so the patch continues the surrounding surface
//!    instead of spanning the hole flat.
//!
//! The patch is wound like the surrounding faces, since every patch face
//! comes from splitting the one face that closes the hole, and it shares
//! the hole's boundary edges, so a mesh with only this hole becomes closed.
//!
//! # Refusals
//!
//! A boundary loop whose projection onto its own mean plane (the Newell
//! normal's) is not a simple polygon is refused: such a loop either
//! crosses itself or folds over, and the minimum-weight triangulation of
//! it would fold too. This is conservative. A strongly curved hole, for
//! example a slot wrapping more than halfway round a cylinder, projects
//! onto itself and is refused although it could be filled. Intersections
//! between the finished patch and the rest of the mesh are not checked;
//! run `axiolid-heal`'s `self_intersections` on the result when that
//! matters.
//!
//! Every refusal leaves the mesh unchanged: the work happens on a copy that
//! replaces the mesh only on success.

use core::f64::consts::PI;

use axiolid_contracts::Sign;
use axiolid_core::{Point2, Point3, Scalar, Vec3};
use axiolid_mesh::halfedge::{
    EdgeId, FaceId, HalfedgeEditError, HalfedgeId, HalfedgeMesh, VertexId,
};
use axiolid_predicates::orient2d;

use crate::fair::{fair, FairError, FairOptions, FairReport};

/// Options of [`fill_hole`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HoleFillOptions {
    /// Refine the triangulation to the density of the surrounding mesh.
    /// Without refinement the patch uses only the boundary vertices and
    /// there is nothing to fair.
    pub refine: bool,
    /// Liepa's density factor `alpha`; larger gives smaller patch
    /// triangles. Must be positive and finite. Default `sqrt 2`.
    pub density: Scalar,
    /// Fair the new interior vertices; `None` leaves them where refinement
    /// put them (on the flat initial triangulation).
    pub fairing: Option<FairOptions>,
    /// Largest hole accepted, in boundary vertices. The triangulation costs
    /// `O(n^4)` time and `O(n^3)` memory: about 30 MB at the default 200.
    pub max_boundary_vertices: usize,
    /// Largest number of new vertices refinement may add.
    pub max_new_vertices: usize,
}

impl Default for HoleFillOptions {
    fn default() -> Self {
        Self {
            refine: true,
            density: core::f64::consts::SQRT_2,
            fairing: Some(FairOptions::default()),
            max_boundary_vertices: 200,
            max_new_vertices: 1_000_000,
        }
    }
}

/// What [`fill_hole`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct HoleFillReport {
    /// Vertices on the hole's boundary loop.
    pub boundary_vertices: usize,
    /// The patch faces, in id order.
    pub faces: Vec<FaceId>,
    /// The new interior vertices, in id order.
    pub vertices: Vec<VertexId>,
    /// Largest dihedral angle, in radians, of the initial triangulation:
    /// between neighbouring patch triangles and across the hole's boundary.
    pub max_dihedral: Scalar,
    /// The fairing solve, when one ran.
    pub fairing: Option<FairReport>,
}

/// Why [`fill_hole`] refused. On every refusal the mesh is unchanged.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum HoleFillError {
    /// The halfedge id names no live halfedge.
    #[error("halfedge {halfedge} is not a live halfedge")]
    RemovedHalfedge {
        /// The halfedge.
        halfedge: HalfedgeId,
    },
    /// The halfedge has a face, so it bounds no hole.
    #[error("halfedge {halfedge} has a face; it does not bound a hole")]
    NotBoundary {
        /// The halfedge.
        halfedge: HalfedgeId,
    },
    /// The hole has more boundary vertices than allowed.
    #[error("hole has {boundary_vertices} boundary vertices, over the limit of {limit}")]
    HoleTooLarge {
        /// Boundary vertices of the hole.
        boundary_vertices: usize,
        /// The limit in force.
        limit: usize,
    },
    /// The boundary loop spans no area (its Newell normal is zero or not
    /// finite), so it has no mean plane.
    #[error("hole boundary spans no area")]
    DegenerateBoundary,
    /// The boundary loop, projected onto its mean plane, crosses or touches
    /// itself. Positions are indices into the loop starting at the given
    /// halfedge.
    #[error("hole boundary edges {first} and {second} cross in projection")]
    SelfIntersectingBoundary {
        /// Position of one boundary edge in the loop.
        first: usize,
        /// Position of the other.
        second: usize,
    },
    /// Every triangulation of the loop needs a degenerate triangle or an
    /// edge the mesh already has.
    #[error("no triangulation of the hole avoids degenerate triangles and existing edges")]
    NoValidTriangulation,
    /// An option is outside its valid range.
    #[error("option {name} {reason}")]
    InvalidOption {
        /// The option.
        name: &'static str,
        /// What it must satisfy.
        reason: &'static str,
    },
    /// Refinement would add more vertices than allowed.
    #[error("refinement would add more than {limit} vertices")]
    BudgetExceeded {
        /// The limit in force.
        limit: usize,
    },
    /// Fairing refused.
    #[error("fairing refused: {0}")]
    Fair(#[from] FairError),
    /// A connectivity edit refused. Not expected on a valid mesh; reported
    /// rather than unwrapped.
    #[error("mesh edit refused: {0}")]
    Edit(#[from] HalfedgeEditError),
}

/// Fill the hole bounded by `boundary`, a boundary halfedge of `mesh`.
///
/// See the [module documentation](self) for the method. The result is
/// deterministic: the same mesh and halfedge give the same patch, ids
/// included.
///
/// # Errors
///
/// [`HoleFillError`]; the mesh is unchanged on every error.
pub fn fill_hole(
    mesh: &mut HalfedgeMesh,
    boundary: HalfedgeId,
    options: &HoleFillOptions,
) -> Result<HoleFillReport, HoleFillError> {
    if !(options.density.is_finite() && options.density > 0.0) {
        return Err(HoleFillError::InvalidOption {
            name: "density",
            reason: "must be positive and finite",
        });
    }
    if !mesh.contains_halfedge(boundary) {
        return Err(HoleFillError::RemovedHalfedge { halfedge: boundary });
    }
    if !mesh.is_boundary_halfedge(boundary) {
        return Err(HoleFillError::NotBoundary { halfedge: boundary });
    }
    let loop_halfedges: Vec<HalfedgeId> = mesh.loop_halfedges(boundary).collect();
    let n = loop_halfedges.len();
    if n > options.max_boundary_vertices {
        return Err(HoleFillError::HoleTooLarge {
            boundary_vertices: n,
            limit: options.max_boundary_vertices,
        });
    }
    let corners: Vec<VertexId> = loop_halfedges.iter().map(|&h| mesh.source(h)).collect();
    let points: Vec<Point3> = corners.iter().map(|&v| mesh.position(v)).collect();
    require_simple_projection(&points)?;

    let plan = triangulate(mesh, &loop_halfedges, &corners, &points)?;
    let max_dihedral = plan.max_dihedral();
    // Scale of each boundary vertex, read before the patch adds edges.
    let boundary_sigma: Vec<Scalar> = corners.iter().map(|&v| mean_edge_length(mesh, v)).collect();

    let mut work = mesh.clone();
    let first_face = work.fill_hole(boundary)?;
    let mut patch = PatchFaces::new(first_face);
    insert_diagonals(&mut work, &mut patch, &corners, &plan)?;

    let mut new_vertices = Vec::new();
    if options.refine {
        let mut sigma = vec![0.0; work.positions().len()];
        for (&v, &s) in corners.iter().zip(&boundary_sigma) {
            sigma[v.index()] = s;
        }
        refine_patch(
            &mut work,
            &mut patch,
            &mut sigma,
            &mut new_vertices,
            options,
        )?;
    }
    let fairing = match options.fairing {
        Some(fair_options) if !new_vertices.is_empty() => {
            Some(fair(&mut work, &new_vertices, &fair_options)?)
        }
        _ => None,
    };
    *mesh = work;
    new_vertices.sort_unstable();
    Ok(HoleFillReport {
        boundary_vertices: n,
        faces: patch.sorted(),
        vertices: new_vertices,
        max_dihedral,
        fairing,
    })
}

/// The patch's faces as a membership table that grows with the mesh.
struct PatchFaces {
    member: Vec<bool>,
}

impl PatchFaces {
    fn new(first: FaceId) -> Self {
        let mut out = Self { member: Vec::new() };
        out.insert(first);
        out
    }

    fn insert(&mut self, f: FaceId) {
        if self.member.len() <= f.index() {
            self.member.resize(f.index() + 1, false);
        }
        self.member[f.index()] = true;
    }

    fn contains(&self, f: Option<FaceId>) -> bool {
        f.is_some_and(|f| self.member.get(f.index()).copied().unwrap_or(false))
    }

    fn sorted(&self) -> Vec<FaceId> {
        self.member
            .iter()
            .enumerate()
            .filter(|(_, &m)| m)
            .map(|(i, _)| FaceId::new(i as u32))
            .collect()
    }
}

// ---- step 0: the loop must project to a simple polygon --------------------

fn require_simple_projection(points: &[Point3]) -> Result<(), HoleFillError> {
    let n = points.len();
    let centroid = points.iter().fold(Vec3::ZERO, |s, &p| s + p) / n as Scalar;
    let mut normal = Vec3::ZERO;
    for i in 0..n {
        let a = points[i] - centroid;
        let b = points[(i + 1) % n] - centroid;
        normal += a.cross(b);
    }
    let length = normal.length();
    if !(length.is_finite() && length > 0.0) {
        return Err(HoleFillError::DegenerateBoundary);
    }
    let normal = normal / length;
    let u = normal.any_orthonormal_vector();
    let v = normal.cross(u);
    let flat: Vec<Point2> = points
        .iter()
        .map(|&p| {
            let d = p - centroid;
            Point2::new(d.dot(u), d.dot(v))
        })
        .collect();
    for i in 0..n {
        let (a, b) = (flat[i], flat[(i + 1) % n]);
        for j in (i + 1)..n {
            let (c, d) = (flat[j], flat[(j + 1) % n]);
            let adjacent = j == i + 1 || (i == 0 && j == n - 1);
            let crossing = if adjacent {
                // Neighbouring edges share a corner; they intersect
                // elsewhere only by folding back along one line.
                let (shared, p, q) = if j == i + 1 { (b, a, d) } else { (a, b, c) };
                n > 3 && side(p, shared, q) == 0 && (p - shared).dot(q - shared) > 0.0
            } else {
                segments_meet(a, b, c, d)
            };
            if crossing {
                return Err(HoleFillError::SelfIntersectingBoundary {
                    first: i,
                    second: j,
                });
            }
        }
    }
    Ok(())
}

/// Exact orientation sign of `c` against `a -> b`: 1, -1 or 0.
fn side(a: Point2, b: Point2, c: Point2) -> i32 {
    match orient2d(a, b, c).sign() {
        Some(Sign::Positive) => 1,
        Some(Sign::Negative) => -1,
        _ => 0,
    }
}

/// Whether closed segments `ab` and `cd` share a point.
fn segments_meet(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    let (o1, o2, o3, o4) = (side(a, b, c), side(a, b, d), side(c, d, a), side(c, d, b));
    if o1 * o2 < 0 && o3 * o4 < 0 {
        return true;
    }
    let within = |p: Point2, q: Point2, r: Point2| {
        r.x >= p.x.min(q.x) && r.x <= p.x.max(q.x) && r.y >= p.y.min(q.y) && r.y <= p.y.max(q.y)
    };
    (o1 == 0 && within(a, b, c))
        || (o2 == 0 && within(a, b, d))
        || (o3 == 0 && within(c, d, a))
        || (o4 == 0 && within(c, d, b))
}

// ---- step 1: minimum-weight triangulation ---------------------------------
//
// Liepa's and CGAL's recurrence weighs a sub-polygon `i..=k` by its own
// best triangulation and only then compares its top triangle with the
// triangle on the other side of the diagonal `(i, k)`. That is not optimal:
// the best triangulation below a diagonal depends on what lies above it. So
// the state here also names the apex `j` of the triangle above the diagonal
// (`F(i, k | j)`), which makes the largest dihedral angle an exact minimum
// at `O(n^4)` time and `O(n^3)` memory.

/// Weight of a (partial) triangulation: the cosine of its largest dihedral
/// angle (larger is better) and its area (smaller is better).
#[derive(Debug, Clone, Copy)]
struct Weight {
    cos: Scalar,
    area: Scalar,
}

/// Cosines closer than this are treated as equal, so area decides between
/// triangulations whose dihedral angles differ only by rounding (a planar
/// hole).
const COS_TIE: Scalar = 1e-12;

impl Weight {
    fn better_than(self, other: Self) -> bool {
        if (self.cos - other.cos).abs() <= COS_TIE {
            self.area < other.area
        } else {
            self.cos > other.cos
        }
    }
}

/// Smallest `|sin|` of a triangle's angle for it to count as non-degenerate.
const MIN_TRIANGLE_SINE: Scalar = 1e-10;

/// One state `F(i, k | j)`: the best weight and the apex `m` of the
/// triangle on diagonal `(i, k)`; `m == u32::MAX` when infeasible.
#[derive(Debug, Clone, Copy)]
struct Cell {
    weight: Weight,
    apex: u32,
}

const INFEASIBLE: Cell = Cell {
    weight: Weight {
        cos: Scalar::NEG_INFINITY,
        area: Scalar::INFINITY,
    },
    apex: u32::MAX,
};

/// The table of states, `F(i, k | j)` for every diagonal `(i, k)` with
/// `k - i >= 2` and every apex `j` outside `i..=k`; the closing edge
/// `(0, n - 1)` has one state, for the face outside the hole.
struct Plan {
    n: usize,
    /// First cell of diagonal `(i, k)`, at `i * n + k`.
    offset: Vec<usize>,
    cells: Vec<Cell>,
}

impl Plan {
    fn slot(&self, i: usize, k: usize, j: usize) -> usize {
        let base = self.offset[i * self.n + k];
        if i == 0 && k == self.n - 1 {
            base
        } else if j < i {
            base + j
        } else {
            base + j - (k - i + 1)
        }
    }

    fn cell(&self, i: usize, k: usize, j: usize) -> Cell {
        self.cells[self.slot(i, k, j)]
    }

    fn max_dihedral(&self) -> Scalar {
        self.cell(0, self.n - 1, 0)
            .weight
            .cos
            .clamp(-1.0, 1.0)
            .acos()
    }
}

fn unit_normal(a: Point3, b: Point3, c: Point3) -> Option<Vec3> {
    let e1 = b - a;
    let e2 = c - a;
    let e3 = c - b;
    let cross = e1.cross(e2);
    let length = cross.length();
    // `|sin|` of every corner angle must clear the threshold.
    if !(length > MIN_TRIANGLE_SINE * e1.length() * e2.length()
        && length > MIN_TRIANGLE_SINE * e1.length() * e3.length()
        && length > MIN_TRIANGLE_SINE * e2.length() * e3.length())
    {
        return None;
    }
    Some(cross / length)
}

fn face_normal(mesh: &HalfedgeMesh, f: FaceId) -> Vec3 {
    // Newell's normal: well defined for any planar or near-planar polygon.
    let corners: Vec<Point3> = mesh.face_vertices(f).map(|v| mesh.position(v)).collect();
    let mut normal = Vec3::ZERO;
    for i in 0..corners.len() {
        normal += corners[i].cross(corners[(i + 1) % corners.len()]);
    }
    normal.normalize_or_zero()
}

fn triangulate(
    mesh: &HalfedgeMesh,
    loop_halfedges: &[HalfedgeId],
    corners: &[VertexId],
    points: &[Point3],
) -> Result<Plan, HoleFillError> {
    let n = corners.len();
    // Normal of the face across each boundary edge `i -> i + 1`.
    let outside: Vec<Vec3> = loop_halfedges
        .iter()
        .map(|&h| {
            let f = mesh
                .face(mesh.opposite(h))
                .expect("an edge never has a boundary on both sides");
            face_normal(mesh, f)
        })
        .collect();
    let mut offset = vec![0usize; n * n];
    let mut total = 0;
    for span in 2..n {
        for i in 0..n - span {
            let k = i + span;
            offset[i * n + k] = total;
            total += if span == n - 1 { 1 } else { n - span - 1 };
        }
    }
    let mut plan = Plan {
        n,
        offset,
        cells: vec![INFEASIBLE; total],
    };
    // Normals of the triangles above diagonal `(i, k)`, one per apex `j`
    // outside `i..=k` (`None` when degenerate), or the outside face's for
    // the closing edge.
    let mut above: Vec<(usize, Option<Vec3>)> = Vec::with_capacity(n);
    for span in 2..n {
        for i in 0..n - span {
            let k = i + span;
            let closing = i == 0 && k == n - 1;
            if !closing && mesh.find_halfedge(corners[i], corners[k]).is_some() {
                continue;
            }
            above.clear();
            if closing {
                above.push((0, Some(outside[n - 1])));
            } else {
                for j in (0..i).chain(k + 1..n) {
                    let normal = if j < i {
                        unit_normal(points[j], points[i], points[k])
                    } else {
                        unit_normal(points[i], points[k], points[j])
                    };
                    above.push((j, normal));
                }
            }
            for m in (i + 1)..k {
                let Some(tn) = unit_normal(points[i], points[m], points[k]) else {
                    continue;
                };
                // The best of each side given that triangle (i, m, k) lies
                // above it.
                let left = if m == i + 1 {
                    Weight {
                        cos: tn.dot(outside[i]),
                        area: 0.0,
                    }
                } else {
                    plan.cell(i, m, k).weight
                };
                let right = if k == m + 1 {
                    Weight {
                        cos: tn.dot(outside[m]),
                        area: 0.0,
                    }
                } else {
                    plan.cell(m, k, i).weight
                };
                if left.cos == Scalar::NEG_INFINITY || right.cos == Scalar::NEG_INFINITY {
                    continue;
                }
                let below = Weight {
                    cos: left.cos.min(right.cos),
                    area: left.area
                        + right.area
                        + 0.5
                            * (points[m] - points[i])
                                .cross(points[k] - points[i])
                                .length(),
                };
                for &(j, normal) in &above {
                    let Some(normal) = normal else { continue };
                    let candidate = Weight {
                        cos: below.cos.min(tn.dot(normal)),
                        area: below.area,
                    };
                    let slot = plan.slot(i, k, j);
                    let cell = &mut plan.cells[slot];
                    if cell.apex == u32::MAX || candidate.better_than(cell.weight) {
                        *cell = Cell {
                            weight: candidate,
                            apex: m as u32,
                        };
                    }
                }
            }
        }
    }
    if plan.cell(0, n - 1, 0).apex == u32::MAX {
        return Err(HoleFillError::NoValidTriangulation);
    }
    Ok(plan)
}

/// Turn the face that closes the hole into the planned triangles.
fn insert_diagonals(
    mesh: &mut HalfedgeMesh,
    patch: &mut PatchFaces,
    corners: &[VertexId],
    plan: &Plan,
) -> Result<(), HoleFillError> {
    let n = corners.len();
    // Diagonal and the apex above it.
    let mut stack = vec![(0usize, n - 1, 0usize)];
    while let Some((i, k, j)) = stack.pop() {
        let m = plan.cell(i, k, j).apex as usize;
        for (a, b, above) in [(i, m, k), (m, k, i)] {
            if b - a >= 2 {
                connect(mesh, patch, corners[a], corners[b])?;
                stack.push((a, b, above));
            }
        }
    }
    Ok(())
}

/// Join two corners of one patch face by a diagonal.
fn connect(
    mesh: &mut HalfedgeMesh,
    patch: &mut PatchFaces,
    a: VertexId,
    b: VertexId,
) -> Result<(), HoleFillError> {
    let mut found = None;
    for h in mesh.incoming_halfedges(a) {
        if !patch.contains(mesh.face(h)) {
            continue;
        }
        if let Some(g) = mesh.loop_halfedges(h).find(|&g| mesh.target(g) == b) {
            found = Some((h, g));
            break;
        }
    }
    let (h, g) = found.expect("the planned diagonal's corners share a patch face");
    let diagonal = mesh.split_face_diagonal(h, g)?;
    for side in [diagonal, mesh.opposite(diagonal)] {
        if let Some(f) = mesh.face(side) {
            patch.insert(f);
        }
    }
    Ok(())
}

// ---- step 2: refinement ---------------------------------------------------

fn mean_edge_length(mesh: &HalfedgeMesh, v: VertexId) -> Scalar {
    let p = mesh.position(v);
    let (sum, count) = mesh.vertex_vertices(v).fold((0.0, 0usize), |(s, c), u| {
        (s + (mesh.position(u) - p).length(), c + 1)
    });
    sum / count as Scalar
}

/// Rounds of edge flips after each round of splits. Each flip strictly
/// raises the smallest angle of its quadrilateral in the plane; the cap only
/// guards against cycling on strongly curved patches.
const MAX_RELAX_SWEEPS: usize = 100;

fn refine_patch(
    mesh: &mut HalfedgeMesh,
    patch: &mut PatchFaces,
    sigma: &mut Vec<Scalar>,
    new_vertices: &mut Vec<VertexId>,
    options: &HoleFillOptions,
) -> Result<(), HoleFillError> {
    let alpha = options.density;
    loop {
        let mut split_any = false;
        for f in patch.sorted() {
            let corners: Vec<VertexId> = mesh.face_vertices(f).collect();
            let positions: Vec<Point3> = corners.iter().map(|&v| mesh.position(v)).collect();
            let centroid = (positions[0] + positions[1] + positions[2]) / 3.0;
            let sigma_c = corners.iter().map(|&v| sigma[v.index()]).sum::<Scalar>() / 3.0;
            let split = corners.iter().zip(&positions).all(|(&v, &p)| {
                let d = alpha * (centroid - p).length();
                d > sigma_c && d > sigma[v.index()]
            });
            if !split {
                continue;
            }
            if new_vertices.len() >= options.max_new_vertices {
                return Err(HoleFillError::BudgetExceeded {
                    limit: options.max_new_vertices,
                });
            }
            let v = mesh.split_face(f, centroid)?;
            if sigma.len() <= v.index() {
                sigma.resize(v.index() + 1, 0.0);
            }
            sigma[v.index()] = sigma_c;
            for g in mesh.vertex_faces(v).collect::<Vec<_>>() {
                patch.insert(g);
            }
            new_vertices.push(v);
            split_any = true;
        }
        if !split_any {
            return Ok(());
        }
        relax(mesh, patch);
    }
}

/// Flip interior patch edges that violate the Delaunay angle criterion.
fn relax(mesh: &mut HalfedgeMesh, patch: &PatchFaces) {
    for _ in 0..MAX_RELAX_SWEEPS {
        let mut flipped = false;
        let edges: Vec<EdgeId> = mesh
            .edges()
            .filter(|&e| {
                let h = mesh.edge_halfedge(e, 0);
                patch.contains(mesh.face(h)) && patch.contains(mesh.face(mesh.opposite(h)))
            })
            .collect();
        for e in edges {
            if should_flip(mesh, e) && mesh.flip_edge(e).is_ok() {
                flipped = true;
            }
        }
        if !flipped {
            return;
        }
    }
}

/// Margin over pi before an edge is flipped, so co-circular quadrilaterals
/// do not flip back and forth.
const FLIP_MARGIN: Scalar = 1e-10;

fn should_flip(mesh: &HalfedgeMesh, e: EdgeId) -> bool {
    let h = mesh.edge_halfedge(e, 0);
    let g = mesh.opposite(h);
    let a = mesh.position(mesh.source(h));
    let b = mesh.position(mesh.target(h));
    let c = mesh.position(mesh.target(mesh.next(h)));
    let d = mesh.position(mesh.target(mesh.next(g)));
    let angle = |apex: Point3| (a - apex).angle_between(b - apex);
    if angle(c) + angle(d) <= PI + FLIP_MARGIN {
        return false;
    }
    // The flipped triangles must not be degenerate.
    unit_normal(c, d, b).is_some() && unit_normal(d, c, a).is_some()
}
