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
///
/// earcut may leave out a corner that lies on a straight run of its ring.
/// Planar faces that share such a corner with a neighbour must use
/// [`earcut_planar_face`] instead, which puts those corners back.
pub(crate) fn earcut_projected(flat: &[[Scalar; 2]], hole_starts: &[usize]) -> Vec<usize> {
    let mut earcutter = earcut::Earcut::new();
    let mut indices: Vec<usize> = Vec::new();
    earcutter.earcut(flat.iter().copied(), hole_starts, &mut indices);
    indices
}

/// [`earcut_projected`] for a planar face whose corners are shared with
/// neighbouring faces: no triangle edge passes over a corner of the face.
///
/// earcut returns a correct triangulation of ONE polygon, but two of its
/// habits break a mesh built from several faces:
///
/// - it drops a corner that lies on a straight run of its ring (its
///   `filter_points` removes any corner with zero turn), and
/// - it may bridge collinear corners of different rings with one long
///   edge, e.g. a window head flush with a door head in the same wall face.
///
/// Either way a triangle edge of this face passes over a corner that the
/// neighbouring face still splits its edge at, so the two faces stop
/// sharing an edge and the mesh gets a T-junction crack. A closed authored
/// shell then fails its closure check and a boolean sees a non-manifold
/// operand. Real exports carry these corners wherever a face meets a finer
/// neighbour (a wall face against a split floor slab, a reveal against a
/// window lining, openings whose heads line up).
///
/// Curved faces keep [`earcut_projected`]: their parameter-space boundary
/// is refined afterwards against the surface, which this split would
/// disturb.
pub(crate) fn earcut_planar_face(
    flat: &[[Scalar; 2]],
    hole_starts: &[usize],
    linear: Scalar,
) -> Vec<usize> {
    let mut indices = earcut_projected(flat, hole_starts);
    if indices.len() % 3 == 0 && !indices.is_empty() {
        split_invented_edges(flat, hole_starts, collinear_band(linear), &mut indices);
    }
    indices
}

/// How far off a segment a corner may lie and still count as on it.
///
/// A thousandth of the model's linear tolerance: 1 um at
/// `Tolerance::MILLIMETRE`. Measured on 743,346 ring corners of two real
/// exports (EFH, ArchiCAD), float noise on shared corners stays below
/// 1e-7 m and genuine turns start at 1e-5 m, with only 28 corners in
/// between; 1e-6 m sits in that valley. It is tied to the declared
/// tolerance, not to the face's own size, so a long straight stair flight
/// whose corners wander 5e-8 m off their line is judged the same on every
/// face that shares them, and it stays three orders below anything the
/// model treats as a feature.
pub(crate) fn collinear_band(linear: Scalar) -> Scalar {
    1.0e-3 * linear
}

/// Split every edge earcut invented where it passes over a corner of the
/// face.
///
/// An authored ring edge is shared with the neighbouring face exactly as
/// written, so it is never split. Every other triangle edge is earcut's:
/// a diagonal, a hole bridge, or a shortcut over a corner it dropped from
/// a straight run. When such an edge passes over a corner `c` of the face,
/// the neighbour that meets the face at `c` has an edge ending there, so
/// the mesh cracks unless the edge is split at `c`.
///
/// Two steps:
///
/// 1. Slivers go. Along a straight run earcut can emit a zero-area
///    triangle `(a, m, b)` with `m` on segment `a-b`. Its three edges cancel
///    on that line once every edge there is split at every corner, so it
///    adds nothing but duplicate edges; it is dropped.
/// 2. The corners strictly inside each remaining invented edge are found
///    once per undirected edge, so the two triangles sharing it insert the
///    same corners and still pair up. A triangle whose sides gained corners
///    is re-triangulated by [`split_triangle`] into proper triangles only.
///
/// The result covers the same area with the same winding, and every
/// segment between consecutive corners on a line is used once from each
/// side, which is what a closed mesh needs.
///
/// Work is `O(triangles * corners)`; beyond `MAX_SPLIT_WORK` the face keeps
/// earcut's output unchanged, which is no worse than before this pass and
/// is reported downstream as an open mesh, not hidden.
fn split_invented_edges(
    flat: &[[Scalar; 2]],
    hole_starts: &[usize],
    noise: Scalar,
    indices: &mut Vec<usize>,
) {
    const MAX_SPLIT_WORK: usize = 1 << 26;
    let n = flat.len();
    let triangles = indices.len() / 3;
    if triangles == 0 || triangles.saturating_mul(3).saturating_mul(n) > MAX_SPLIT_WORK {
        return;
    }
    let mut authored: std::collections::HashSet<(usize, usize)> =
        std::collections::HashSet::with_capacity(n);
    let mut bounds: Vec<usize> = Vec::with_capacity(hole_starts.len() + 2);
    bounds.push(0);
    bounds.extend(hole_starts.iter().copied().filter(|&s| s > 0 && s < n));
    bounds.push(n);
    for ring in bounds.windows(2) {
        let (start, end) = (ring[0], ring[1]);
        for k in start..end {
            let next = if k + 1 == end { start } else { k + 1 };
            authored.insert((k.min(next), k.max(next)));
        }
    }

    // Step 1: drop slivers, triangles with one corner on the segment
    // between the other two, where that segment is an edge earcut
    // invented. Its two short sides then lie on the same straight run, and
    // splitting the invented long side at the corner cancels the triangle
    // against itself, so it carries no area and no edge the mesh needs.
    //
    // When the long side is an authored ring edge the triangle is real
    // geometry, however thin: a faceted B-rep face may itself be a 1 um
    // wide triangle (a Revit fastener has 1,096 of them), and its long
    // edge is shared with a neighbour exactly as written. Dropping it
    // would delete the face.
    let is_sliver = |t: &[usize]| {
        (0..3).any(|i| {
            let (a, b) = (t[(i + 1) % 3], t[(i + 2) % 3]);
            !authored.contains(&(a.min(b), a.max(b)))
                && strictly_between(flat[a], flat[t[i]], flat[b], noise)
        })
    };
    let kept: Vec<usize> = indices
        .chunks_exact(3)
        .filter(|t| !is_sliver(t))
        .flatten()
        .copied()
        .collect();

    // Step 2: corners strictly inside each invented edge, ordered from the
    // lower index to the higher, computed on the canonical (low, high) key
    // so both sides of an edge see the identical list.
    let mut on_edge: std::collections::HashMap<(usize, usize), Vec<usize>> =
        std::collections::HashMap::new();
    for t in kept.chunks_exact(3) {
        for i in 0..3 {
            let (p, q) = (t[i], t[(i + 1) % 3]);
            let key = (p.min(q), p.max(q));
            if authored.contains(&key) || on_edge.contains_key(&key) {
                continue;
            }
            let (a, b) = (flat[key.0], flat[key.1]);
            let mut found: Vec<(Scalar, usize)> = flat
                .iter()
                .enumerate()
                .filter(|&(c, _)| c != key.0 && c != key.1)
                .filter_map(|(c, &point)| along_if_between(a, point, b, noise).map(|s| (s, c)))
                .collect();
            found.sort_by(|x, y| x.0.total_cmp(&y.0));
            on_edge.insert(key, found.into_iter().map(|(_, c)| c).collect());
        }
    }
    if kept.len() == indices.len() && on_edge.values().all(Vec::is_empty) {
        return;
    }

    let mut out: Vec<usize> = Vec::with_capacity(kept.len());
    for t in kept.chunks_exact(3) {
        // Corners on side i run from t[i] to t[i+1].
        let side = |i: usize| -> Vec<usize> {
            let (p, q) = (t[i], t[(i + 1) % 3]);
            let key = (p.min(q), p.max(q));
            let mut run = on_edge.get(&key).cloned().unwrap_or_default();
            if p > q {
                run.reverse();
            }
            run
        };
        let sides = [side(0), side(1), side(2)];
        if sides.iter().all(Vec::is_empty) {
            out.extend_from_slice(t);
        } else {
            split_triangle([t[0], t[1], t[2]], &sides, &mut out);
        }
    }
    *indices = out;
}

/// Triangulate a proper triangle `t` whose side `i` (from `t[i]` to
/// `t[i+1]`) carries the corners `sides[i]` in order, into proper
/// triangles with the same winding.
///
/// The triangle plus its side corners is a convex polygon whose points all
/// lie on three lines. Collinearity is decided from which sides a point is
/// on, not from coordinates: `t[i]` is on sides `i` and `i-1`, a side
/// corner only on its own side, and three points are collinear exactly when
/// they share a side. Ears are clipped while one exists whose three points
/// are not collinear and whose removal leaves a polygon that is not
/// collinear either; clipping such an ear from a convex polygon leaves a
/// convex polygon, so the loop always finds the next one and ends at a
/// proper triangle. `O(k^2)` for `k` points, and `k` is small.
fn split_triangle(t: [usize; 3], sides: &[Vec<usize>; 3], out: &mut Vec<usize>) {
    // (index, side-membership bitmask)
    let mut ring: Vec<(usize, u8)> =
        Vec::with_capacity(3 + sides.iter().map(Vec::len).sum::<usize>());
    for i in 0..3 {
        ring.push((t[i], (1 << i) | (1 << ((i + 2) % 3))));
        ring.extend(sides[i].iter().map(|&c| (c, 1u8 << i)));
    }
    let collinear =
        |masks: &mut dyn Iterator<Item = u8>| masks.fold(0b111u8, |acc, m| acc & m) != 0;
    while ring.len() > 3 {
        let n = ring.len();
        // Clip the first corner whose removal leaves a ring that is not
        // all on one line. That alone keeps every emitted ear proper: the
        // ring is a triangle with corners on its sides, so the first such
        // corner always joins two different lines (checked exhaustively
        // for up to four corners per side; see
        // `scripts/probe_closed_mesh_mutants.py`).
        let ear = (0..n).find(|&k| {
            !collinear(
                &mut ring
                    .iter()
                    .enumerate()
                    .filter(|&(j, _)| j != k)
                    .map(|(_, p)| p.1),
            )
        });
        let Some(k) = ear else {
            // Only reachable if the input was not a proper triangle; keep
            // the triangle unsplit rather than emit a degenerate piece.
            out.extend_from_slice(&t);
            return;
        };
        let (prev, next) = ((k + n - 1) % n, (k + 1) % n);
        out.extend_from_slice(&[ring[prev].0, ring[k].0, ring[next].0]);
        ring.remove(k);
    }
    out.extend_from_slice(&[ring[0].0, ring[1].0, ring[2].0]);
}

/// Where `c` lies along `a-b` (as `along / length`, in (0, 1)) if it is on
/// the open segment, else `None`. See [`strictly_between`].
fn along_if_between(
    a: [Scalar; 2],
    c: [Scalar; 2],
    b: [Scalar; 2],
    noise: Scalar,
) -> Option<Scalar> {
    if !strictly_between(a, c, b, noise) {
        return None;
    }
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ac = [c[0] - a[0], c[1] - a[1]];
    Some((ab[0] * ac[0] + ab[1] * ac[1]) / (ab[0] * ab[0] + ab[1] * ab[1]))
}

/// Is `c` on the open segment `a-b`? Within `noise` (a distance) of the
/// line, or within 1e-9 of the edge length, and at a parameter strictly
/// inside (0, 1).
fn strictly_between(a: [Scalar; 2], c: [Scalar; 2], b: [Scalar; 2], noise: Scalar) -> bool {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ac = [c[0] - a[0], c[1] - a[1]];
    let length = ab[0] * ab[0] + ab[1] * ab[1];
    let along = ab[0] * ac[0] + ab[1] * ac[1];
    if !(along > 0.0 && along < length) {
        return false;
    }
    // |turn| / |ab| is the distance from c to the line through a and b.
    let turn = ab[0] * ac[1] - ab[1] * ac[0];
    turn.abs() <= (1e-9 * length).max(noise * length.sqrt())
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
    // After the area check, so it judges earcut's own cover; each split
    // keeps winding and moves area only by float noise.
    split_invented_edges(&flat, &hole_starts, collinear_band(linear), &mut indices);
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Undirected edge use counts of a triangle list.
    fn edge_uses(indices: &[usize]) -> std::collections::HashMap<(usize, usize), usize> {
        let mut uses = std::collections::HashMap::new();
        for t in indices.chunks_exact(3) {
            for i in 0..3 {
                let (p, q) = (t[i], t[(i + 1) % 3]);
                *uses.entry((p.min(q), p.max(q))).or_insert(0) += 1;
            }
        }
        uses
    }

    /// Two triangles sharing the invented diagonal 0-2 of a square, with
    /// corner 4 on that diagonal: both halves must split at 4 in the same
    /// order, so every interior edge is still used twice.
    #[test]
    fn a_shared_invented_edge_splits_identically_on_both_sides() {
        let flat = [[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0], [1.0, 1.0]];
        // Ring 0-1-2-3; 4 is a hole corner in real use, here any corner.
        let mut indices = vec![0, 1, 2, 0, 2, 3];
        split_invented_edges(&flat, &[4], collinear_band(1.0e-3), &mut indices);
        let uses = edge_uses(&indices);
        assert_eq!(uses.get(&(0, 4)), Some(&2), "{indices:?}");
        assert_eq!(uses.get(&(2, 4)), Some(&2), "{indices:?}");
        assert_eq!(uses.get(&(0, 2)), None, "the diagonal is gone: {indices:?}");
        let area: Scalar = indices
            .chunks_exact(3)
            .map(|t| doubled_area(flat[t[0]], flat[t[1]], flat[t[2]]))
            .sum();
        assert!((area - 8.0).abs() < 1e-12, "area {area}");
    }

    /// A zero-area triangle earcut left along a straight run: its apex
    /// lies on its opposite side. It must not survive with an unpaired
    /// long edge.
    #[test]
    fn a_zero_area_sliver_leaves_no_unpaired_edge() {
        // Straight run 0-1-2 along y = 0, closing through 3 at the top.
        let flat = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [1.0, 1.0]];
        // earcut-style output with the sliver (0, 1, 2) spanning 0-2.
        let mut indices = vec![0, 1, 2, 0, 2, 3];
        split_invented_edges(&flat, &[], collinear_band(1.0e-3), &mut indices);
        let uses = edge_uses(&indices);
        assert_eq!(uses.get(&(0, 2)), None, "{indices:?}");
        for t in indices.chunks_exact(3) {
            assert!(t[0] != t[1] && t[1] != t[2] && t[2] != t[0], "{indices:?}");
        }
    }

    /// The census slab (`arc_boe` #61210): a ring corner 5.9 nm off an
    /// invented 1 m edge on a face a few metres across. That is exporter
    /// noise, so the zero-width sliver through it must go and the edge must
    /// split there, or the edge is used four times and the closed slab
    /// reads as non-manifold.
    #[test]
    fn a_corner_nanometres_off_an_invented_edge_is_on_it() {
        // Ring 0-1-2-3-4; 0-2 is invented, 1 lies 5.9e-9 off it.
        let flat = [
            [0.0, 0.0],
            [0.24, 5.9e-9],
            [1.0, 0.0],
            [0.5, -1.0],
            [-2.0, 2.0],
        ];
        let mut indices = vec![0, 1, 2, 0, 2, 3];
        split_invented_edges(&flat, &[], collinear_band(1.0e-3), &mut indices);
        let uses = edge_uses(&indices);
        assert_eq!(uses.get(&(0, 2)), None, "{indices:?}");
        assert_eq!(uses.get(&(0, 1)), Some(&1), "{indices:?}");
        assert_eq!(uses.get(&(1, 2)), Some(&1), "{indices:?}");
        assert!(uses.values().all(|&n| n <= 2), "{indices:?}");
    }

    /// A real stair stringer (`arc_boe` #21009): a 19-corner sawtooth
    /// whose inner corners 7, 9, 11, 13, 15, 17 lie on one line up to
    /// 1e-8..5e-8 m of export noise. earcut bridges them with long
    /// diagonals; with a band scaled to the face's size those diagonals
    /// were split unevenly (one side saw a corner, the other did not) and
    /// three edges were left unpaired. With the tolerance-tied band every
    /// edge is used at most twice and the authored ring edges exactly once.
    #[test]
    fn a_stair_stringer_with_noisy_collinear_steps_stays_paired() {
        // (x, z) of the stringer's outer ring, projected onto its plane
        // along the stair direction (x, y) = (-0.4659, 0.8848) of the export.
        let raw: [(f64, f64, f64); 19] = [
            (19.055730261854272, -2.9397826126246684, -2.8),
            (18.838372983807588, -2.526678826402099, -2.8),
            (18.231554326511727, -1.3733745196303921, -1.906627452209237),
            (18.231554326511727, -1.3733745196303921, -1.754117659085965),
            (18.161708907836992, -1.2406280702631145, -1.754117659085965),
            (18.161708907836992, -1.2406280702631145, -1.654117588496447),
            (18.329337942590353, -1.5592196056363878, -1.654117588496447),
            (18.329337942590353, -1.5592196056363878, -1.832352844473842),
            (18.450403334222592, -1.789313449805038, -1.832352844473842),
            (18.450403334222592, -1.789313449805038, -2.010588133233825),
            (18.57146869532532, -2.019407235950066, -2.010588133233825),
            (18.57146869532532, -2.019407235950066, -2.18882338921122),
            (18.69253411748707, -2.2495011381423398, -2.18882338921122),
            (18.69253411748707, -2.2495011381423398, -2.36705871075379),
            (18.81359949385455, -2.4795949532991797, -2.36705871075379),
            (18.81359949385455, -2.4795949532991797, -2.545294065078947),
            (18.934664854957276, -2.7096887394442053, -2.545294065078947),
            (18.934664854957276, -2.7096887394442053, -2.723529321056342),
            (19.055730261854272, -2.9397826126246684, -2.723529321056342),
        ];
        let o = raw[0];
        let dir = {
            let (dx, dy) = (raw[1].0 - o.0, raw[1].1 - o.1);
            let l = (dx * dx + dy * dy).sqrt();
            (dx / l, dy / l)
        };
        let flat: Vec<[Scalar; 2]> = raw
            .iter()
            .map(|p| [(p.0 - o.0) * dir.0 + (p.1 - o.1) * dir.1, p.2 - o.2])
            .collect();
        let mut indices = earcut_projected(&flat, &[]);
        split_invented_edges(&flat, &[], collinear_band(1.0e-3), &mut indices);
        let uses = edge_uses(&indices);
        assert!(uses.values().all(|&n| n <= 2), "over-used edge: {uses:?}");
        for k in 0..flat.len() {
            let key = (k.min((k + 1) % flat.len()), k.max((k + 1) % flat.len()));
            assert_eq!(uses.get(&key), Some(&1), "ring edge {key:?}: {indices:?}");
        }
        for (key, &n) in &uses {
            let ring = (key.1 - key.0 == 1) || (key.0 == 0 && key.1 == flat.len() - 1);
            if !ring {
                assert_eq!(n, 2, "invented edge {key:?} unpaired: {indices:?}");
            }
        }
    }

    /// A face that is itself one thin triangle keeps it: every side is an
    /// authored ring edge. The numbers are a real Revit fastener face
    /// (`twp_sfw` #16626), 0.79 um from apex to base, inside the 1 mm
    /// collinearity band.
    #[test]
    fn a_thin_authored_triangle_face_is_kept() {
        let flat = [
            [1.198193270495232, -0.262123154757972],
            [1.2022602528833906, -0.267596768443262],
            [1.2022292140899142, -0.2675563179000615],
        ];
        let mut indices = earcut_projected(&flat, &[]);
        assert_eq!(indices.len(), 3, "earcut keeps the one triangle");
        split_invented_edges(&flat, &[], collinear_band(1.0e-3), &mut indices);
        assert_eq!(indices.len(), 3, "the face is not a sliver: {indices:?}");
    }

    /// An authored ring edge is shared with the neighbour exactly as
    /// written, so a corner that merely lies near its line never splits it.
    #[test]
    fn an_authored_ring_edge_is_never_split() {
        // Ring 0-1-2-3 with corner 4 on edge 0-1's line: 4 belongs to a
        // hole and 0-1 is authored, so 0-1 must stay whole.
        let flat = [[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0], [2.0, 0.0]];
        let mut indices = vec![0, 1, 2, 0, 2, 3];
        split_invented_edges(&flat, &[4], collinear_band(1.0e-3), &mut indices);
        assert_eq!(edge_uses(&indices).get(&(0, 1)), Some(&1), "{indices:?}");
    }
}
