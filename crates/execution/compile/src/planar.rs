//! Planar polygon triangulation shared by every compiler path that turns a
//! flat, possibly concave or holed polygon into triangles.
//!
//! B-rep faces and authored polygon faces (#160) go through the same three
//! steps, so the same polygon triangulates the same way whichever node it
//! arrived in:
//!
//! 1. a plane from the Newell normal of the outer ring, which is stable for
//!    concave rings and for rings whose first corners are collinear;
//! 2. projection onto orthonormal in-plane axes `(u, v)` with `u x v = n`,
//!    so a ring that winds counter-clockwise about `n` stays
//!    counter-clockwise in 2D;
//! 3. ear clipping (`earcut`) of the projected outer ring with its holes.

use axiolid_contracts::GeomError;
use axiolid_core::{Scalar, Vec3};

/// Newell normal of a closed ring: correct for concave rings and rings that
/// start with collinear corners. Its length is twice the ring's area.
///
/// A cross product of the first two edges fails when they are collinear,
/// which is common at the start of an exported ring.
pub(crate) fn newell_normal(ring: impl IntoIterator<Item = Vec3>) -> Vec3 {
    let mut normal = Vec3::ZERO;
    let mut iter = ring.into_iter();
    let Some(first) = iter.next() else {
        return normal;
    };
    let mut current = first;
    for next in iter.chain(std::iter::once(first)) {
        normal.x += (current.y - next.y) * (current.z + next.z);
        normal.y += (current.z - next.z) * (current.x + next.x);
        normal.z += (current.x - next.x) * (current.y + next.y);
        current = next;
    }
    normal
}

/// Orthonormal in-plane axes `(u, v)` with `u x v = n / |n|`, or `None` when
/// the normal is zero or not finite.
pub(crate) fn plane_axes(normal: Vec3) -> Option<(Vec3, Vec3)> {
    let length = normal.length();
    if !length.is_finite() || length <= f64::EPSILON {
        return None;
    }
    let n = normal / length;
    // Pick the axis least aligned with n so the cross product stays stable.
    let helper = if n.x.abs() <= n.y.abs() && n.x.abs() <= n.z.abs() {
        Vec3::X
    } else if n.y.abs() <= n.z.abs() {
        Vec3::Y
    } else {
        Vec3::Z
    };
    let u = n.cross(helper).normalize();
    let v = n.cross(u);
    Some((u, v))
}

/// Ear-clip a projected polygon: the outer ring first, each hole starting at
/// the matching entry of `hole_starts`. Returns three indices per triangle
/// into `flat`, or an empty list when earcut found no triangle.
pub(crate) fn earcut_projected(flat: &[[Scalar; 2]], hole_starts: &[usize]) -> Vec<usize> {
    let mut earcutter = earcut::Earcut::new();
    let mut indices: Vec<usize> = Vec::new();
    earcutter.earcut(flat.iter().copied(), hole_starts, &mut indices);
    indices
}

/// Twice the signed area of a projected triangle.
fn doubled_area(a: [Scalar; 2], b: [Scalar; 2], c: [Scalar; 2]) -> Scalar {
    (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
}

/// Twice the signed area of a projected ring (shoelace).
fn doubled_ring_area(ring: &[[Scalar; 2]]) -> Scalar {
    let mut sum = 0.0;
    for index in 0..ring.len() {
        let p = ring[index];
        let q = ring[(index + 1) % ring.len()];
        sum += p[0] * q[1] - q[0] * p[1];
    }
    sum
}

/// Why an authored polygon face could not be triangulated.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum PolygonRefusal {
    /// The outer ring has zero or non-finite area, so it has no plane.
    NoPlane,
    /// A corner lies further from the face plane than the tolerance; the
    /// value is that distance. A non-planar polygon has no unique
    /// triangulation, so choosing one would invent geometry.
    NotPlanar(Scalar),
    /// The triangles do not cover the polygon's area: the rings cross
    /// themselves or each other, or a hole is not inside the outer ring.
    AreaMismatch {
        /// Area of the outer ring minus its holes.
        polygon: Scalar,
        /// Area covered by the triangles.
        triangles: Scalar,
    },
}

impl core::fmt::Display for PolygonRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoPlane => f.write_str("its outer ring encloses no area"),
            Self::NotPlanar(distance) => write!(
                f,
                "it is not planar: a corner lies {distance:e} from the face plane, \
                 beyond the linear tolerance"
            ),
            Self::AreaMismatch { polygon, triangles } => write!(
                f,
                "its rings cross or a hole lies outside the outer ring \
                 (polygon area {polygon:e}, triangulated area {triangles:e})"
            ),
        }
    }
}

/// Triangulate one planar polygon given as an outer ring plus holes.
///
/// Returns three indices per triangle into the rings concatenated in order
/// (outer, then each hole). Every triangle winds counter-clockwise about the
/// outer ring's Newell normal, i.e. in the outer ring's authored sense.
///
/// # Refusals
///
/// Returns a [`PolygonRefusal`] rather than a triangulation that could be
/// wrong: no plane, a corner off the plane by more than `linear`, or
/// triangles whose area differs from the polygon's (earcut returns a partial
/// result on crossing rings instead of failing, so the area is checked).
pub(crate) fn triangulate_polygon(
    rings: &[&[Vec3]],
    linear: Scalar,
) -> Result<Vec<usize>, PolygonRefusal> {
    let Some(outer) = rings.first() else {
        return Err(PolygonRefusal::NoPlane);
    };
    let normal = newell_normal(outer.iter().copied());
    let (u, v) = plane_axes(normal).ok_or(PolygonRefusal::NoPlane)?;
    let n = normal.normalize();

    // The plane through the outer ring's centroid, with the Newell normal:
    // Newell's is the least-squares plane direction for a near-planar ring.
    let centroid = outer.iter().copied().fold(Vec3::ZERO, |sum, p| sum + p) / outer.len() as Scalar;
    let mut worst: Scalar = 0.0;
    for ring in rings {
        for &p in *ring {
            let distance = (p - centroid).dot(n).abs();
            if !distance.is_finite() {
                return Err(PolygonRefusal::NoPlane);
            }
            worst = worst.max(distance);
        }
    }
    if worst > linear {
        return Err(PolygonRefusal::NotPlanar(worst));
    }

    let mut flat: Vec<[Scalar; 2]> = Vec::new();
    let mut hole_starts: Vec<usize> = Vec::with_capacity(rings.len().saturating_sub(1));
    let mut expected = 0.0;
    for (index, ring) in rings.iter().enumerate() {
        if index > 0 {
            hole_starts.push(flat.len());
        }
        let start = flat.len();
        flat.extend(ring.iter().map(|&p| {
            let d = p - centroid;
            [d.dot(u), d.dot(v)]
        }));
        let area = doubled_ring_area(&flat[start..]).abs();
        expected += if index == 0 { area } else { -area };
    }

    let mut indices = earcut_projected(&flat, &hole_starts);
    if indices.len() % 3 != 0 {
        indices.clear();
    }
    let mut covered = 0.0;
    for triangle in indices.chunks_exact_mut(3) {
        let area = doubled_area(flat[triangle[0]], flat[triangle[1]], flat[triangle[2]]);
        // The outer ring is counter-clockwise in (u, v) by construction of
        // the axes, so every triangle must be too. earcut 0.4 already emits
        // them that way, so this never fires today; it keeps the authored
        // winding independent of earcut's output convention across
        // upgrades (an equivalent mutant in the probe, deliberately kept).
        if area < 0.0 {
            triangle.swap(1, 2);
        }
        covered += area.abs();
    }
    // Earcut only adds diagonals between existing corners, so a correct
    // result covers the polygon exactly up to rounding in the projection.
    let scale = flat
        .iter()
        .fold(0.0_f64, |m, p| m.max(p[0].abs()).max(p[1].abs()));
    let slack = 1e-9 * expected.abs().max(scale * scale) + f64::EPSILON;
    if expected <= slack || (covered - expected).abs() > slack {
        return Err(PolygonRefusal::AreaMismatch {
            polygon: expected / 2.0,
            triangles: covered / 2.0,
        });
    }
    Ok(indices)
}

/// A typed error naming the authored face that could not be triangulated.
pub(crate) fn face_error(face: usize, refusal: PolygonRefusal) -> GeomError {
    match refusal {
        PolygonRefusal::NoPlane => GeomError::Degenerate(format!(
            "authored polygon face {face} cannot be triangulated: {refusal}"
        )),
        _ => GeomError::InvalidInput(format!(
            "authored polygon face {face} cannot be triangulated: {refusal}"
        )),
    }
}
