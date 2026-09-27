//! Certified line of sight from an eye point to a mesh past blockers
//! (#185).
//!
//! # Visible: one ray, checked exactly
//!
//! A ray from the eye through a point `p` is a witness when it crosses the
//! interior of a target triangle strictly before it meets any blocker
//! triangle -- touching a blocker's edge counts as meeting it. Where a ray
//! meets a plane is a ratio of two exact `orient3d` values, so which comes
//! first is an exact comparison. Candidate rays aim at points spread over
//! each target triangle; the first that checks out is returned.
//!
//! # Hidden: every ray, by cones
//!
//! Every target triangle is cut into sub-triangles whose corners are exact
//! dyadic points on it (midpoints of dyadic points are dyadic, so they stay
//! on the triangle exactly). A sub-triangle is hidden when one blocker
//! piece covers it: its three corners' rays pass strictly inside the
//! piece's cone from the eye -- so every ray through the sub-triangle does,
//! the cone being convex -- and a plane of the piece has the eye strictly
//! on one side and the sub-triangle strictly on the other, so every such
//! ray meets the piece before the target. A piece is a blocker triangle,
//! two coplanar triangles forming a convex quadrilateral (a wall), or a
//! whole closed convex blocker (a column). When every sub-triangle of
//! every target triangle is covered, the target is hidden, and the
//! blockers used are named.
//!
//! # Undecided
//!
//! Neither argument may be available: a target just grazed, or covered only
//! by several pieces together, so that some sub-cone straddles a seam
//! between them at every depth. Then the answer is undecided, never a
//! guess.

use axiolid_core::{Point3, Tolerance};
use axiolid_exact::{certify, Arith, Dyadic, SignExpr};
use axiolid_guarantees::Sign;
use axiolid_mesh::{audit_mesh, TriangleMeshView};

/// Sub-triangles the hidden argument examines at most.
pub const MAX_SIGHT_CELLS: usize = 50_000;

/// What can be proven about the view from an eye to a target.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Sight {
    /// Some part of the target is in view.
    Visible {
        /// The ray from the eye through this point is the witness.
        through: Point3,
        /// The target triangle it crosses, strictly inside, before any
        /// blocker.
        triangle: usize,
    },
    /// No part of the target is in view.
    Hidden {
        /// Indices of the blockers the argument used, ascending.
        occluders: Vec<usize>,
    },
    /// Neither could be proven within the budget.
    Undecided,
}

/// Why the question was not asked.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SightError {
    /// The eye or a vertex is not finite.
    NonFinite,
    /// The target has no triangles with area.
    EmptyTarget,
}

/// Whether any part of `target` is in view from `eye`, past `blockers`.
///
/// # Errors
///
/// [`SightError`] for non-finite input or an empty target.
pub fn line_of_sight<T, B>(eye: Point3, target: &T, blockers: &[&B]) -> Result<Sight, SightError>
where
    T: TriangleMeshView + ?Sized,
    B: TriangleMeshView + ?Sized,
{
    line_of_sight_within(eye, target, blockers, MAX_SIGHT_CELLS)
}

/// [`line_of_sight`] with a caller-chosen budget of sub-triangles.
///
/// # Errors
///
/// As [`line_of_sight`].
pub fn line_of_sight_within<T, B>(
    eye: Point3,
    target: &T,
    blockers: &[&B],
    budget: usize,
) -> Result<Sight, SightError>
where
    T: TriangleMeshView + ?Sized,
    B: TriangleMeshView + ?Sized,
{
    if !eye.is_finite() {
        return Err(SightError::NonFinite);
    }
    let targets = triangles(target)?;
    if targets.is_empty() {
        return Err(SightError::EmptyTarget);
    }
    let mut walls: Vec<[Point3; 3]> = Vec::new();
    let mut pieces: Vec<Piece> = Vec::new();
    for (index, mesh) in blockers.iter().enumerate() {
        let own = triangles(*mesh)?;
        pieces.extend(pieces_of(eye, index, &own, *mesh));
        walls.extend(own);
    }
    if let Some(sight) = witness(eye, &targets, &walls) {
        return Ok(sight);
    }
    Ok(hidden(eye, &targets, &pieces, budget)
        .map_or(Sight::Undecided, |occluders| Sight::Hidden { occluders }))
}

fn triangles<M: TriangleMeshView + ?Sized>(mesh: &M) -> Result<Vec<[Point3; 3]>, SightError> {
    let mut out = Vec::with_capacity(mesh.triangle_count());
    for t in 0..mesh.triangle_count() {
        let corners = mesh.triangle(t).map(|i| mesh.position(i as usize));
        if !corners.iter().all(|p| p.is_finite()) {
            return Err(SightError::NonFinite);
        }
        let [a, b, c] = corners;
        if !collinear(a, b, c) {
            out.push(corners);
        }
    }
    Ok(out)
}

/// Whether three points are collinear, exactly.
fn collinear(a: Point3, b: Point3, c: Point3) -> bool {
    let d = |p: Point3| [exact(p.x), exact(p.y), exact(p.z)];
    let (a, b, c) = (d(a), d(b), d(c));
    let u = sub(&b, &a);
    let v = sub(&c, &a);
    let n = cross(&u, &v);
    n.iter().all(|x| x.sign() == Some(Sign::Zero))
}

// ---------------------------------------------------------------------
// Visible
// ---------------------------------------------------------------------

/// A ray through a point spread over some target triangle that meets it
/// before every blocker.
fn witness(eye: Point3, targets: &[[Point3; 3]], walls: &[[Point3; 3]]) -> Option<Sight> {
    // Barycentric weights: the centroid first, then a grid inside.
    let mut weights: Vec<[f64; 3]> = vec![[1.0 / 3.0; 3]];
    let n = 6;
    for i in 1..n {
        for j in 1..n - i {
            let (u, v) = (i as f64 / n as f64, j as f64 / n as f64);
            weights.push([u, v, 1.0 - u - v]);
        }
    }
    for (index, t) in targets.iter().enumerate() {
        for w in &weights {
            let through = t[0] * w[0] + t[1] * w[1] + t[2] * w[2];
            if through == eye {
                continue;
            }
            if let Some(near) = strict_hit(eye, through, t) {
                if walls.iter().all(|wall| !blocks(eye, through, wall, &near)) {
                    return Some(Sight::Visible {
                        through,
                        triangle: index,
                    });
                }
            }
        }
    }
    None
}

/// Where the ray from `eye` through `p` meets a triangle's plane: the
/// parameter `t = O(eye) / (O(eye) - O(p))` along `eye + t (p - eye)`,
/// kept as numerator and denominator, exact.
#[derive(Debug, Clone)]
struct Param {
    num: Dyadic,
    den: Dyadic,
}

impl Param {
    fn of(eye: Point3, p: Point3, t: &[Point3; 3]) -> Option<Self> {
        let oe = orient(&t.map(dpoint), &dpoint(eye));
        let op = orient(&t.map(dpoint), &dpoint(p));
        let den = oe.sub(&op);
        if den.sign() == Some(Sign::Zero) {
            return None;
        }
        Some(Self { num: oe, den })
    }

    fn sign(&self) -> Sign {
        times(
            self.num.sign().unwrap_or(Sign::Zero),
            self.den.sign().unwrap_or(Sign::Zero),
        )
    }

    /// The sign of `self - other`.
    fn compare(&self, other: &Self) -> Sign {
        let gap = self.num.mul(&other.den).sub(&other.num.mul(&self.den));
        let d = times(
            self.den.sign().unwrap_or(Sign::Zero),
            other.den.sign().unwrap_or(Sign::Zero),
        );
        times(gap.sign().unwrap_or(Sign::Zero), d)
    }
}

/// The ray's parameter at a triangle it crosses strictly inside, ahead of
/// the eye.
fn strict_hit(eye: Point3, p: Point3, t: &[Point3; 3]) -> Option<Param> {
    let s = inside_signs(eye, p, t);
    let all_same = s[0] != Sign::Zero && s[0] == s[1] && s[1] == s[2];
    if !all_same {
        return None;
    }
    let param = Param::of(eye, p, t)?;
    (param.sign() == Sign::Positive).then_some(param)
}

/// Whether a blocker triangle meets the ray (edges included) no farther
/// than `near`.
fn blocks(eye: Point3, p: Point3, wall: &[Point3; 3], near: &Param) -> bool {
    let s = inside_signs(eye, p, wall);
    let has_pos = s.contains(&Sign::Positive);
    let has_neg = s.contains(&Sign::Negative);
    if has_pos && has_neg {
        return false;
    }
    let Some(param) = Param::of(eye, p, wall) else {
        // Parallel to the wall's plane: off it, the ray never meets the
        // wall; in it, count the wall as blocking, conservatively.
        let eye_side = orient(&wall.map(dpoint), &dpoint(eye)).sign();
        return eye_side == Some(Sign::Zero);
    };
    param.sign() != Sign::Negative && param.compare(near) != Sign::Positive
}

/// `orient3d(eye, p, t_i, t_{i+1})` for each edge: all one strict sign
/// when the line through the eye and `p` passes strictly inside.
fn inside_signs(eye: Point3, p: Point3, t: &[Point3; 3]) -> [Sign; 3] {
    [0, 1, 2].map(|i| {
        sign(&Orient3 {
            p: [dpoint(eye), dpoint(p), dpoint(t[i]), dpoint(t[(i + 1) % 3])],
        })
    })
}

// ---------------------------------------------------------------------
// Hidden
// ---------------------------------------------------------------------

type DPoint = [Dyadic; 3];

/// A blocker piece: a convex polygon or a convex solid, with the planes
/// bounding its cone from the eye and the planes that may separate it
/// from the target.
struct Piece {
    blocker: usize,
    /// Planes through the eye, as point triples, whose positive side is
    /// the cone's inside.
    cone: Vec<[DPoint; 3]>,
    /// Candidate separating planes, as point triples whose positive side
    /// holds the eye.
    faces: Vec<[DPoint; 3]>,
}

fn pieces_of<M: TriangleMeshView + ?Sized>(
    eye: Point3,
    blocker: usize,
    own: &[[Point3; 3]],
    mesh: &M,
) -> Vec<Piece> {
    let e = dpoint(eye);
    let mut out = Vec::new();
    let mut polygon = |ring: Vec<Point3>| {
        let ring: Vec<DPoint> = ring.into_iter().map(dpoint).collect();
        // The eye must be off the polygon's plane.
        let side = orient(&[ring[0].clone(), ring[1].clone(), ring[2].clone()], &e);
        let Some(s) = side.sign().filter(|s| *s != Sign::Zero) else {
            return;
        };
        let face = if s == Sign::Positive {
            [ring[0].clone(), ring[1].clone(), ring[2].clone()]
        } else {
            [ring[0].clone(), ring[2].clone(), ring[1].clone()]
        };
        // Seen from the eye the ring turns one way; orient the cone planes
        // so the inside is positive.
        let n = ring.len();
        let planes: Vec<[DPoint; 3]> = (0..n)
            .map(|i| [e.clone(), ring[i].clone(), ring[(i + 1) % n].clone()])
            .collect();
        let inward = orient(&planes[0], &ring[2 % n]).sign() == Some(Sign::Positive);
        let cone = planes
            .into_iter()
            .map(|[a, b, c]| if inward { [a, b, c] } else { [a, c, b] })
            .collect();
        out.push(Piece {
            blocker,
            cone,
            faces: vec![face],
        });
    };
    for t in own {
        polygon(t.to_vec());
    }
    // Coplanar edge neighbours forming a convex quadrilateral: a wall.
    for (i, a) in own.iter().enumerate() {
        for b in &own[i + 1..] {
            if let Some(quad) = convex_quad(a, b) {
                polygon(quad);
            }
        }
    }
    if let Some(solid) = convex_solid(eye, blocker, own, mesh) {
        out.push(solid);
    }
    out
}

/// Two triangles sharing an edge, coplanar, whose union is a convex
/// quadrilateral: its ring.
fn convex_quad(a: &[Point3; 3], b: &[Point3; 3]) -> Option<Vec<Point3>> {
    let shared: Vec<usize> = (0..3).filter(|&i| b.contains(&a[i])).collect();
    if shared.len() != 2 {
        return None;
    }
    let apex_a = (0..3).find(|i| !shared.contains(i))?;
    let apex_b = *b.iter().find(|p| !a.contains(p))?;
    // Ring: a's apex, then round a to the shared edge, b's apex between.
    let (u, v) = (a[(apex_a + 1) % 3], a[(apex_a + 2) % 3]);
    let ring = vec![a[apex_a], u, apex_b, v];
    let d: Vec<DPoint> = ring.iter().copied().map(dpoint).collect();
    if orient(&[d[0].clone(), d[1].clone(), d[2].clone()], &d[3]).sign() != Some(Sign::Zero) {
        return None;
    }
    // Convex: each corner turns the same way within the plane, judged by
    // the cross products against the plane's normal.
    let normal = cross(&sub(&d[1], &d[0]), &sub(&d[3], &d[0]));
    let turns: Vec<Option<Sign>> = (0..4)
        .map(|i| {
            let (p, q, r) = (&d[i], &d[(i + 1) % 4], &d[(i + 2) % 4]);
            dot(&cross(&sub(q, p), &sub(r, q)), &normal).sign()
        })
        .collect();
    let first = turns[0]?;
    (first != Sign::Zero && turns.iter().all(|t| *t == Some(first))).then_some(ring)
}

/// A closed convex blocker the eye is strictly outside of, as one piece:
/// its cone is bounded by the planes through the eye and two of its
/// vertices with every vertex on one side.
fn convex_solid<M: TriangleMeshView + ?Sized>(
    eye: Point3,
    blocker: usize,
    own: &[[Point3; 3]],
    mesh: &M,
) -> Option<Piece> {
    if !audit_mesh(mesh, Tolerance::ZERO).is_closed_two_manifold() {
        return None;
    }
    let mut vertices: Vec<Point3> = own.iter().flatten().copied().collect();
    vertices.sort_by(|a, b| {
        a.x.total_cmp(&b.x)
            .then(a.y.total_cmp(&b.y))
            .then(a.z.total_cmp(&b.z))
    });
    vertices.dedup();
    let dv: Vec<DPoint> = vertices.iter().copied().map(dpoint).collect();
    let e = dpoint(eye);
    // Convex, with each face's inner side: every vertex on one side of it.
    let mut faces = Vec::new();
    let mut outside = false;
    for t in own {
        let plane = t.map(dpoint);
        let mut side = Sign::Zero;
        for v in &dv {
            match orient(&plane, v).sign()? {
                Sign::Zero => {}
                s if side == Sign::Zero => side = s,
                s if s != side => return None,
                _ => {}
            }
        }
        let eye_side = orient(&plane, &e).sign()?;
        if eye_side != side {
            outside = true;
        }
        // Oriented so the eye's side is positive, for separation.
        if eye_side != Sign::Zero {
            faces.push(if eye_side == Sign::Positive {
                plane
            } else {
                [plane[0].clone(), plane[2].clone(), plane[1].clone()]
            });
        }
    }
    if !outside {
        return None;
    }
    let mut cone = Vec::new();
    for i in 0..dv.len() {
        for j in i + 1..dv.len() {
            let plane = [e.clone(), dv[i].clone(), dv[j].clone()];
            let mut side = Sign::Zero;
            let mut ok = true;
            for v in &dv {
                match orient(&plane, v).sign() {
                    Some(Sign::Zero) => {}
                    Some(s) if side == Sign::Zero => side = s,
                    Some(s) if s != side => {
                        ok = false;
                        break;
                    }
                    Some(_) => {}
                    None => return None,
                }
            }
            if ok && side != Sign::Zero {
                cone.push(if side == Sign::Positive {
                    plane
                } else {
                    [plane[0].clone(), plane[2].clone(), plane[1].clone()]
                });
            }
        }
    }
    Some(Piece {
        blocker,
        cone,
        faces,
    })
}

/// Cover every target triangle by pieces, subdividing; the blockers used,
/// or `None` when some part could not be covered within the budget.
fn hidden(
    eye: Point3,
    targets: &[[Point3; 3]],
    pieces: &[Piece],
    budget: usize,
) -> Option<Vec<usize>> {
    let _ = eye;
    let mut used: Vec<usize> = Vec::new();
    let mut stack: Vec<[DPoint; 3]> = targets.iter().map(|t| t.map(dpoint)).collect();
    let mut cells = 0usize;
    while let Some(cell) = stack.pop() {
        cells += 1;
        if cells > budget {
            return None;
        }
        if let Some(piece) = pieces.iter().find(|p| covers(p, &cell)) {
            if !used.contains(&piece.blocker) {
                used.push(piece.blocker);
            }
            continue;
        }
        let half = exact(0.5);
        let mid =
            |a: &DPoint, b: &DPoint| -> DPoint { [0, 1, 2].map(|k| a[k].add(&b[k]).mul(&half)) };
        let [a, b, c] = cell;
        let (ab, bc, ca) = (mid(&a, &b), mid(&b, &c), mid(&c, &a));
        stack.push([a, ab.clone(), ca.clone()]);
        stack.push([ab.clone(), b, bc.clone()]);
        stack.push([ca.clone(), bc.clone(), c]);
        stack.push([ab, bc, ca]);
    }
    used.sort_unstable();
    Some(used)
}

/// Whether the piece hides every point of the cell from the eye.
fn covers(piece: &Piece, cell: &[DPoint; 3]) -> bool {
    let in_cone = piece.cone.iter().all(|plane| {
        cell.iter().all(|g| {
            sign(&Orient3 {
                p: [
                    plane[0].clone(),
                    plane[1].clone(),
                    plane[2].clone(),
                    g.clone(),
                ],
            }) == Sign::Positive
        })
    });
    if !in_cone {
        return false;
    }
    piece.faces.iter().any(|face| {
        cell.iter().all(|g| {
            sign(&Orient3 {
                p: [face[0].clone(), face[1].clone(), face[2].clone(), g.clone()],
            }) == Sign::Negative
        })
    })
}

// ---------------------------------------------------------------------
// Exact helpers
// ---------------------------------------------------------------------

fn exact(x: f64) -> Dyadic {
    Dyadic::from_f64(x)
}

fn dpoint(p: Point3) -> DPoint {
    [exact(p.x), exact(p.y), exact(p.z)]
}

fn sub<T: Arith>(a: &[T; 3], b: &[T; 3]) -> [T; 3] {
    [a[0].sub(&b[0]), a[1].sub(&b[1]), a[2].sub(&b[2])]
}

fn cross<T: Arith>(u: &[T; 3], v: &[T; 3]) -> [T; 3] {
    [
        u[1].mul(&v[2]).sub(&u[2].mul(&v[1])),
        u[2].mul(&v[0]).sub(&u[0].mul(&v[2])),
        u[0].mul(&v[1]).sub(&u[1].mul(&v[0])),
    ]
}

fn dot<T: Arith>(u: &[T; 3], v: &[T; 3]) -> T {
    u[0].mul(&v[0]).add(&u[1].mul(&v[1])).add(&u[2].mul(&v[2]))
}

/// `orient3d(a, b, c, d)` exactly: positive when `d` lies on the side of
/// the plane through `a, b, c` that `(b - a) x (c - a)` points to.
fn orient(plane: &[DPoint; 3], d: &DPoint) -> Dyadic {
    let n = cross(&sub(&plane[1], &plane[0]), &sub(&plane[2], &plane[0]));
    dot(&n, &sub(d, &plane[0]))
}

/// [`orient`] as a filtered predicate: intervals first, dyadics if
/// undecided.
struct Orient3 {
    p: [DPoint; 4],
}

impl SignExpr for Orient3 {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let q = |i: usize| self.p[i].clone().map(|x| T::from_dyadic(&x));
        let (a, b, c, d) = (q(0), q(1), q(2), q(3));
        let n = cross(&sub(&b, &a), &sub(&c, &a));
        dot(&n, &sub(&d, &a)).sign()
    }
}

/// The inputs are finite, so the exact tier always decides.
fn sign<E: SignExpr>(e: &E) -> Sign {
    certify(e).unwrap_or(Sign::Zero)
}

fn times(a: Sign, b: Sign) -> Sign {
    match (a, b) {
        (Sign::Zero, _) | (_, Sign::Zero) => Sign::Zero,
        (x, y) if x == y => Sign::Positive,
        _ => Sign::Negative,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_mesh::TriMesh;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    /// An extrusion of a counter-clockwise outline in z, capped by a fan
    /// from the first vertex (which must see every other).
    fn extrusion(outline: &[(f64, f64)], z: (f64, f64)) -> TriMesh {
        let n = outline.len() as u32;
        let mut positions: Vec<Point3> = outline.iter().map(|&(x, y)| p(x, y, z.0)).collect();
        positions.extend(outline.iter().map(|&(x, y)| p(x, y, z.1)));
        let mut indices = Vec::new();
        for i in 1..n - 1 {
            indices.extend([0, i + 1, i]);
            indices.extend([n, n + i, n + i + 1]);
        }
        for i in 0..n {
            let j = (i + 1) % n;
            indices.extend([i, j, j + n, i, j + n, i + n]);
        }
        TriMesh::new(positions, indices)
    }

    #[test]
    fn only_convex_closed_blockers_are_solids() {
        let eye = p(-10.0, 0.5, 0.0);
        let cube = extrusion(
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
            (-1.0, 1.0),
        );
        let own = triangles(&cube).unwrap();
        assert!(convex_solid(eye, 0, &own, &cube).is_some());
        // An L, fanned from its reflex corner.
        let l = extrusion(
            &[(1.0, 1.0), (0.0, 2.0), (0.0, 0.0), (2.0, 0.0), (2.0, 1.0)],
            (-1.0, 1.0),
        );
        let own = triangles(&l).unwrap();
        assert!(audit_mesh(&l, Tolerance::ZERO).is_closed_two_manifold());
        assert!(convex_solid(eye, 0, &own, &l).is_none());
        // The eye inside a convex solid: not a solid piece.
        let own = triangles(&cube).unwrap();
        assert!(convex_solid(p(0.5, 0.5, 0.0), 0, &own, &cube).is_none());
    }
}
