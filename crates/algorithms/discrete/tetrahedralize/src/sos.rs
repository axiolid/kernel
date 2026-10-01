//! The signs the triangulation is built from, with symbolic perturbation.
//!
//! Every decision is an exact sign from `axiolid-predicates`. Where the exact
//! answer is zero -- five cospherical points, a point in the plane of a hull
//! face and on that face's circumcircle -- the tie is broken by *simulation
//! of simplicity* on the lifted points, following Devillers and Teillaud,
//! "Perturbations for Delaunay and weighted Delaunay 3D triangulations"
//! (CGAL's `Delaunay_triangulation_3::side_of_oriented_sphere`):
//!
//! Each point's lifted coordinate `|p|^2` is raised by `eps^(rank)`, where
//! the lexicographically largest point (`x`, then `y`, then `z`) gets the
//! largest perturbation. The perturbed in-sphere determinant is a polynomial
//! in `eps` whose coefficients are the orientations obtained by replacing one
//! tetrahedron vertex with the query point, so the first non-zero one, taken
//! in descending lexicographic order of the five points, decides. The order
//! depends only on coordinates, never on insertion order, so the perturbed
//! Delaunay triangulation is unique: the same point set gives the same
//! tetrahedra whatever order it is inserted in.

use std::cmp::Ordering;

use axiolid_core::Point3;
use axiolid_guarantees::{Certified, Sign};
use axiolid_predicates::{in_diametral_sphere, insphere, orient2d, orient3d};

use axiolid_core::Point2;

/// The proven sign of a certified predicate.
///
/// Input is range-checked on entry (see `crate::delaunay::check_point`), and
/// within that range every predicate used here certifies. `Zero` is the
/// conservative reading of anything else: it routes to the perturbation,
/// which always yields a definite answer.
fn decided(certified: Certified) -> Sign {
    match certified {
        Certified::Certain { sign, .. } => sign,
        _ => Sign::Zero,
    }
}

/// Sign of the signed volume `(b - a) . ((c - a) x (d - a))`.
///
/// Positive when `d` lies on the side of the plane `abc` that the normal
/// `(b - a) x (c - a)` points to: the right-handed convention, under which
/// the unit tetrahedron `0, e1, e2, e3` is positive. `orient3d` uses the
/// opposite sign, so its last two arguments are swapped here.
pub(crate) fn orientation(a: Point3, b: Point3, c: Point3, d: Point3) -> Sign {
    decided(orient3d(a, b, d, c))
}

/// Whether `a`, `b`, `c` are exactly collinear (or not all distinct).
pub(crate) fn collinear(a: Point3, b: Point3, c: Point3) -> bool {
    [project_xy, project_yz, project_zx]
        .iter()
        .all(|project| decided(orient2d(project(a), project(b), project(c))) == Sign::Zero)
}

fn project_xy(p: Point3) -> Point2 {
    Point2::new(p.x, p.y)
}

fn project_yz(p: Point3) -> Point2 {
    Point2::new(p.y, p.z)
}

fn project_zx(p: Point3) -> Point2 {
    Point2::new(p.z, p.x)
}

/// Lexicographic order on coordinates: the perturbation's ranking.
pub(crate) fn lexicographic(p: Point3, q: Point3) -> Ordering {
    p.x.total_cmp(&q.x)
        .then(p.y.total_cmp(&q.y))
        .then(p.z.total_cmp(&q.z))
}

/// Is point `e` inside the perturbed circumsphere of the positively oriented
/// tetrahedron `t`?
///
/// Never undecided: on an exact tie the perturbation answers. `e` must not
/// coincide with a vertex of `t`.
pub(crate) fn in_sphere(points: &[Point3], t: [u32; 4], e: u32) -> bool {
    let [a, b, c, d] = t.map(|v| points[v as usize]);
    let query = points[e as usize];
    // `insphere` is positive inside for `orient3d`-positive input, which is
    // our orientation with the last two vertices swapped.
    match decided(insphere(a, b, d, c, query)) {
        Sign::Positive => return true,
        Sign::Negative => return false,
        _ => {}
    }
    let mut order = [t[0], t[1], t[2], t[3], e];
    order.sort_unstable_by(|&p, &q| lexicographic(points[q as usize], points[p as usize]));
    for q in order {
        if q == e {
            // Raising the query point's own lift moves it outside.
            return false;
        }
        // Raising vertex q's lift pulls the sphere towards points on q's
        // side of the opposite face: the coefficient is the orientation of
        // `t` with q replaced by e.
        let replaced = t.map(|v| points[if v == q { e } else { v } as usize]);
        match orientation(replaced[0], replaced[1], replaced[2], replaced[3]) {
            Sign::Positive => return true,
            Sign::Negative => return false,
            _ => {}
        }
    }
    // Unreachable for distinct points with `t` non-degenerate: two vanishing
    // coefficients would put three collinear points on one sphere.
    false
}

/// Is point `e` in conflict with the hull face `face` -- strictly beyond it,
/// or in its plane and inside its perturbed circumcircle?
///
/// `face` is ordered so that `orientation(face, x) > 0` for `x` outside the
/// hull. This is the in-sphere test of the infinite cell over that face: the
/// limit of the spheres through the face and a point running off to
/// infinity beyond it.
pub(crate) fn beyond_face(points: &[Point3], face: [u32; 3], e: u32) -> bool {
    let [a, b, c] = face.map(|v| points[v as usize]);
    match orientation(a, b, c, points[e as usize]) {
        Sign::Positive => true,
        Sign::Negative => false,
        _ => in_coplanar_circle(points, face, e),
    }
}

/// Is `e`, coplanar with `face`, inside the face's perturbed circumcircle?
///
/// The perturbation agrees with [`in_sphere`] on every finite cell through
/// the face, which is what keeps the conflict region a topological ball.
fn in_coplanar_circle(points: &[Point3], face: [u32; 3], e: u32) -> bool {
    let [a, b, c] = face.map(|v| points[v as usize]);
    let query = points[e as usize];
    match decided(in_diametral_sphere(a, b, c, query)) {
        Sign::Positive => return true,
        Sign::Negative => return false,
        _ => {}
    }
    // Orientation inside the plane, read off a coordinate projection that
    // the plane maps onto bijectively; one projection is used for every
    // comparison, so their product is meaningful.
    let project = [project_xy, project_yz, project_zx]
        .into_iter()
        .find(|project| decided(orient2d(project(a), project(b), project(c))) != Sign::Zero)
        .unwrap_or(project_xy);
    let side = |r: u32, s: u32, x: u32| {
        decided(orient2d(
            project(points[r as usize]),
            project(points[s as usize]),
            project(points[x as usize]),
        ))
    };
    let mut order = [face[0], face[1], face[2], e];
    order.sort_unstable_by(|&p, &q| lexicographic(points[q as usize], points[p as usize]));
    for q in order {
        if q == e {
            return false;
        }
        let [r, s] = match face.iter().position(|&v| v == q) {
            Some(0) => [face[1], face[2]],
            Some(1) => [face[0], face[2]],
            _ => [face[0], face[1]],
        };
        let query_side = side(r, s, e);
        if query_side != Sign::Zero {
            return query_side == side(r, s, q);
        }
    }
    false
}
