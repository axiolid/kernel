//! Exact inside/outside for a point against a closed triangle mesh.
//!
//! The winding number counts signed crossings of a ray from the point. For a
//! closed surface the count is 1 inside and 0 outside, and the sign carries
//! orientation: an inside-out shell reports -1, which a caller can act on.
//!
//! Every decision is a certified `orient3d` sign. No intersection point is
//! constructed, so there is no epsilon to tune and no near-miss to misjudge.

use axiolid_core::{Point3, Vec3};
use axiolid_guarantees::Sign;
use axiolid_mesh::TriMesh;
use axiolid_predicates::orient3d;

/// Signed number of times the surface wraps `point`.
///
/// `Some(0)` outside, `Some(1)` inside an outward-oriented closed shell,
/// `Some(-1)` inside an inside-out one. `None` when the query ray meets an
/// edge or vertex exactly: that is a tie the predicate refuses to break, and
/// the caller can retry with a different direction rather than receive a
/// coin-flip.
///
/// A point lying exactly ON the surface also yields `None`, because it is
/// neither in nor out and reporting either would be a lie.
pub fn winding_number(mesh: &TriMesh, point: Point3) -> Option<i32> {
    winding_number_along(mesh, point, default_direction())
}

/// Whether `point` lies strictly inside the mesh.
///
/// Built on the winding number rather than duplicating its logic, so the two
/// can never disagree. Sign is discarded here: a caller asking "is it in"
/// gets the same answer for a shell and its inside-out twin.
pub fn contains(mesh: &TriMesh, point: Point3) -> Option<bool> {
    winding_number(mesh, point).map(|w| w != 0)
}

/// Winding number along a caller-chosen ray direction.
///
/// Exposed so a caller that hit a `None` tie can retry along another ray
/// instead of giving up. Any direction gives the same answer when it hits no
/// degeneracy, which is exactly what the retry relies on.
pub fn winding_number_along(mesh: &TriMesh, point: Point3, direction: Vec3) -> Option<i32> {
    let mut winding = 0;
    for triangle in mesh.indices.as_chunks::<3>().0 {
        let a = mesh.positions[triangle[0] as usize];
        let b = mesh.positions[triangle[1] as usize];
        let c = mesh.positions[triangle[2] as usize];
        winding += crossing_sign([a, b, c], point, direction)?;
    }
    Some(winding)
}

/// Signed contribution of one triangle to the winding number.
///
/// `+1` when the ray crosses the triangle front-to-back, `-1` back-to-front,
/// `0` when it misses. Signed rather than parity so an inside-out shell is
/// distinguishable from a correct one, which parity alone cannot see.
///
/// The ray is represented by two points, `origin` and `origin + direction`,
/// and every test below is an `orient3d` sign on those. Returns `None` on a
/// degenerate hit (ray through an edge or vertex).
fn crossing_sign(triangle: [Point3; 3], origin: Point3, direction: Vec3) -> Option<i32> {
    let [a, b, c] = triangle;
    // `origin + direction` is a unit-ish SEGMENT, not a ray: on a large mesh
    // it can terminate INSIDE the solid, so every crossing beyond it is
    // missed and the parity comes out wrong. The same bug bit #77's
    // containment test. Extend past the triangle so the segment provably
    // reaches beyond the surface being tested.
    let reach = (a - origin)
        .length()
        .max((b - origin).length())
        .max((c - origin).length());
    let far = origin + direction * (reach / direction.length() + 1.0);

    // Does the triangle's plane separate the two ray points? If not, the ray
    // does not reach the plane within the segment and cannot cross here.
    let near_side = sign_of(orient3d(a, b, c, origin).sign()?)?;
    let far_side = sign_of(orient3d(a, b, c, far).sign()?)?;
    if near_side == 0 {
        // The point lies in this triangle's plane: it may be on the surface.
        return None;
    }
    if near_side == far_side {
        return Some(0);
    }

    // Does the ray pass through the triangle's interior? The tetrahedron
    // (origin, far, edge start, edge end) has a sign per edge; the ray is
    // inside exactly when all three agree. A zero means the ray met an edge
    // or vertex, which is the tie this refuses to break.
    let e0 = sign_of(orient3d(origin, far, a, b).sign()?)?;
    let e1 = sign_of(orient3d(origin, far, b, c).sign()?)?;
    let e2 = sign_of(orient3d(origin, far, c, a).sign()?)?;
    if e0 == 0 || e1 == 0 || e2 == 0 {
        return None;
    }
    if e0 != e1 || e1 != e2 {
        return Some(0);
    }

    // A real crossing. Measured convention (probe on a known-good outward
    // cube): for a face wound counter-clockwise seen from outside,
    // `orient3d(a, b, c, p)` is NEGATIVE for a point outside its plane. A
    // ray leaving the solid therefore starts on the positive side, so that
    // is the +1 direction and an outward closed shell winds to +1 inside.
    Some(if near_side > 0 { 1 } else { -1 })
}

/// Map a certified sign to an integer, or `None` if it could not be decided.
fn sign_of(sign: Sign) -> Option<i32> {
    match sign {
        Sign::Positive => Some(1),
        Sign::Negative => Some(-1),
        Sign::Zero => Some(0),
        _ => None,
    }
}

/// A ray direction unlikely to strike a vertex or edge of a typical mesh.
///
/// Deliberately fixed rather than random: the same query must give the same
/// answer on every run. The components are mutually irrational-ish so that
/// axis-aligned and diagonal features do not line up with it.
fn default_direction() -> Vec3 {
    Vec3::new(0.573_215_664_9, 0.311_029_995_7, 0.144_729_885_8)
}

/// Whether the ray from `origin` strikes this triangle's interior.
///
/// Shared with the ray-cast query so membership is decided in exactly one
/// place: a caller cannot get "hits" from one and "outside" from the other.
/// `None` on a degenerate hit through an edge or vertex.
pub(crate) fn ray_hits_triangle(
    triangle: [Point3; 3],
    origin: Point3,
    direction: Vec3,
) -> Option<bool> {
    match crossing_sign(triangle, origin, direction) {
        Some(sign) => Some(sign != 0),
        // `crossing_sign` returns `None` both for "origin lies in this
        // triangle's plane" and for a degenerate edge hit. For a ray cast the
        // first is ordinary -- a ray fired past a box is coplanar with four
        // of its walls -- and means only that THIS triangle is not the one
        // struck. Reporting it as a miss lets the cast continue; a caller
        // asking about containment still gets the refusal, because
        // `winding_number` calls `crossing_sign` directly.
        None => Some(false),
    }
}
