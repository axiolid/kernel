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
//! 3. the certified ear clipper of `axiolid-construct` (ADR 0083, #260)
//!    on the projected outer ring with its holes, under
//!    [`PinchPolicy::Accept`]: a face is a surface patch, so rings touching
//!    at a single point (a hole whose corner sits on the outer ring or on
//!    another hole, an outer ring pinched at a vertex) bound a valid face,
//!    and every edge of the shell is still shared by two faces once a
//!    corner the clipper inserts into a shared edge is inserted into the
//!    neighbour's copy too ([`edge_touches`], [`crate::weld`], #265). Its
//!    certificate proves the triangles tile the projected face exactly
//!    once, every ring edge a triangle edge, so no triangle edge runs past
//!    a corner and the face's own edges are the edges its neighbours
//!    share; rings that do not bound a region (crossing, overlapping) are
//!    refused by name.
//!
//! # Warped authored faces (#254)
//!
//! An authored polygon face whose corners lie off its plane by more than
//! the linear tolerance has no unique surface: a warped quad has two
//! triangulations and neither is the face. [`triangulate_polygon`] refuses
//! it ([`PolygonRefusal::NotPlanar`]); the mesh compiler's authored polygon
//! path uses [`triangulate_authored_polygon`] instead, which triangulates
//! it through the same three steps and reports how far it is warped. The
//! exact compiler refuses polygon meshes altogether.
//!
//! The steps are unchanged, so a face within the tolerance of its plane
//! triangulates exactly as before. A warped face is projected onto its fit
//! plane `P` (the outer ring's centroid with its Newell normal `n`), ear
//! clipped there with the same area check, and the triangles are lifted
//! back to the authored corners, which are never moved. Rings that cross
//! in that projection, or enclose no area in it, are still refused.
//!
//! # The reported bound: the slab width (#261)
//!
//! Let `d_i` be the signed distance of corner `i` (outer ring and holes)
//! from `P` along `n`, and `S` the slab of points whose distance lies in
//! `[min d_i, max d_i]`. The reported bound is its width
//! `W = max d_i - min d_i`, plus the rounding of computing it.
//!
//! A warped face has no single true surface. Either diagonal of a quad,
//! the face flattened onto `P` (each corner replaced by its foot) and a
//! bilinear patch through the corners are all reasonable readings, and
//! they differ from each other by up to the whole spread of the corners
//! along `n`: with one corner of a square lifted by `h`, the two diagonal
//! triangulations are `h / 2` apart at the centre while every corner is
//! only `h / 4` from `P`. `W` bounds that difference for every reading:
//!
//! - Every reading lies in `S`. A triangulation through the corners or a
//!   bilinear patch lies in the corners' convex hull, and `S` is convex
//!   and holds every corner. `P` itself lies in `S`, so the flattened
//!   face does too: `P` passes through the outer ring's centroid, so the
//!   outer corners' `d_i` sum to zero and `min d_i <= 0 <= max d_i`.
//! - The mesh lies in `S` too, and it is a graph over the projected face:
//!   its triangles are those of a triangulation of the face's projection
//!   `R` onto `P` (outer ring minus holes, certified by the clipper), each lifted to
//!   the authored corners. Over every point `x'` of `R` there is exactly
//!   one mesh point `x = x' + h(x') n`, with `h(x')` in `[min d_i, max d_i]`.
//! - So for every point `y` of any reading in `S` whose projection `y'`
//!   onto `P` lies in `R`, the mesh point over `y'` differs from `y` only
//!   along `n`, by the difference of two distances in `[min d_i, max d_i]`:
//!   `y` is within `W` of the mesh.
//!
//! Readings that cover exactly `R` in projection get the bound both ways:
//! every triangulation through the corners that triangulates `R` in
//! projection, the flattened face (which is `R` lifted onto `P`), and a
//! bilinear patch over a quad whose projection is convex (projection is
//! linear, so the patch projects to the bilinear patch over the projected
//! corners, which covers exactly that convex quad). Over each point of `R`
//! such a reading and the mesh are then within `W` of each other along
//! `n`, so each is within `W` of the other. What `W` does not bound is a
//! reading that spills outside `R` in projection, over a concave notch or
//! a hole the face does not cover; that is a different face, not a reading
//! of this one. Every authored corner is a mesh vertex.
//!
//! `W` lies between the largest corner distance from `P` (the bound before
//! #261, which `min d_i <= 0 <= max d_i` keeps below it) and twice that.
//! A saddle whose corners alternate `+h` and `-h` reports `2h`. A square
//! with one corner lifted by `h` is `z = h x y`: its fit plane takes the
//! linear part and leaves a saddle of `+-h/4` about it, so it reports
//! `h/2` (along the tilted normal, a fraction of a percent less), exactly
//! the gap between its two diagonal triangulations at the centre; the
//! plane of the three unlifted corners would give a slab `h` wide, which
//! also holds every reading but is twice as loose. A face counts as warped,
//! as before, when a corner lies further than the linear tolerance from
//! `P`.

use axiolid_construct::profile::{ring_touches, triangulate_with, PinchPolicy, Rings};
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Scalar, Vec3};

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

/// Triangulate a projected polygon with the certified ear clipper (#260):
/// the outer ring first, each hole starting at the matching entry of
/// `hole_starts`. Returns three indices per triangle into `flat`, every
/// triangle counter-clockwise, every corner used.
///
/// A corner repeating the one before it exactly (an exporter's closing
/// point, a doubled corner) is dropped first and appears in no triangle;
/// earcut dropped it too. Rings that touch at single points are accepted
/// ([`PinchPolicy::Accept`], see the module notes).
///
/// # Errors
///
/// The clipper's refusal, by name: rings with fewer than three distinct
/// corners, non-finite, crossing, overlapping or folding back, a hole
/// outside the outer ring or inside another, or a triangulation that fails
/// its certificate.
pub(crate) fn clip_projected(
    flat: &[[Scalar; 2]],
    hole_starts: &[usize],
) -> GeomResult<Vec<usize>> {
    let (rings, kept) = projected_rings(flat, hole_starts);
    let (_, triangles) = triangulate_with(&rings, PinchPolicy::Accept)?;
    Ok(triangles
        .into_iter()
        .flatten()
        .map(|corner| kept[corner as usize])
        .collect())
}

/// The rings [`clip_projected`] hands the clipper, and for each of their
/// points (in `outer ++ holes` order) its index in `flat`: a corner
/// repeating the one before it exactly is dropped, as is a ring's closing
/// repeat of its first corner.
fn projected_rings(flat: &[[Scalar; 2]], hole_starts: &[usize]) -> (Rings, Vec<usize>) {
    let mut bounds = Vec::with_capacity(hole_starts.len() + 2);
    bounds.push(0);
    bounds.extend(hole_starts.iter().map(|&start| start.min(flat.len())));
    bounds.push(flat.len());
    let mut kept: Vec<usize> = Vec::with_capacity(flat.len());
    let mut rings: Vec<Vec<Point2>> = Vec::with_capacity(bounds.len() - 1);
    for span in bounds.windows(2) {
        let mut ring: Vec<Point2> = Vec::with_capacity(span[1].saturating_sub(span[0]));
        for (index, corner) in flat.iter().enumerate().take(span[1]).skip(span[0]) {
            let point = Point2::new(corner[0], corner[1]);
            if ring.last() != Some(&point) {
                ring.push(point);
                kept.push(index);
            }
        }
        while ring.len() > 1 && ring.first() == ring.last() {
            ring.pop();
            kept.pop();
        }
        rings.push(ring);
    }
    let outer = rings.remove(0);
    (
        Rings {
            outer,
            holes: rings,
        },
        kept,
    )
}

/// A corner of a face lying inside one of the face's ring edges (#265):
/// the edge runs from corner `from` to corner `to`, indices into `flat`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct EdgeTouch {
    pub from: usize,
    pub to: usize,
    pub corner: usize,
}

/// Every corner that [`clip_projected`] inserts into a ring edge of the
/// same face, found by the clipper's own exact tests
/// ([`axiolid_construct::profile::ring_touches`]) on exactly the rings it
/// triangulates. Empty when the rings do not bound a region: the
/// triangulation then refuses them by name.
pub(crate) fn edge_touches(flat: &[[Scalar; 2]], hole_starts: &[usize]) -> Vec<EdgeTouch> {
    let (rings, kept) = projected_rings(flat, hole_starts);
    let Ok(touches) = ring_touches(&rings) else {
        return Vec::new();
    };
    let mut starts = Vec::with_capacity(rings.holes.len() + 1);
    let mut start = 0;
    for ring in std::iter::once(&rings.outer).chain(&rings.holes) {
        starts.push((start, ring.len()));
        start += ring.len();
    }
    touches
        .into_iter()
        .map(|touch| {
            let (start, len) = starts[touch.ring];
            EdgeTouch {
                from: kept[start + touch.edge],
                to: kept[start + (touch.edge + 1) % len],
                corner: kept[touch.vertex],
            }
        })
        .collect()
}

/// [`clip_projected`] for a planar face whose corners are shared with
/// neighbouring faces, then the noise-band split of
/// [`split_invented_edges`].
///
/// The clipper never runs a triangle edge exactly past a corner, and every
/// ring edge is a triangle edge, so the face meets its neighbours along its
/// own edges. Export noise still puts corners a few nanometres off the
/// straight run they lie on, and the clipper may then cut a sliver across
/// them whose long diagonal another face, sharing those corners, cuts too:
/// that edge is used four times and the closed shell reads as
/// non-manifold. The split removes such slivers as it did for earcut.
pub(crate) fn clip_planar_face(
    flat: &[[Scalar; 2]],
    hole_starts: &[usize],
    linear: Scalar,
) -> GeomResult<Vec<usize>> {
    let mut indices = clip_projected(flat, hole_starts)?;
    split_invented_edges(flat, hole_starts, collinear_band(linear), &mut indices);
    Ok(indices)
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

/// Split every edge the triangulation invented where it passes over a
/// corner of the face, within the collinearity band.
///
/// An authored ring edge is shared with the neighbouring face exactly as
/// written, so it is never split. Every other triangle edge is invented:
/// a diagonal or a hole bridge (and, under earcut before #260, a shortcut
/// over a corner it dropped from a straight run; the certified clipper
/// passes exactly over no corner, so only corners within the band but off
/// the edge remain). When such an edge passes over a corner `c` of the face,
/// the neighbour that meets the face at `c` has an edge ending there, so
/// the mesh cracks unless the edge is split at `c`.
///
/// Two steps:
///
/// 1. Slivers go. Along a straight run the triangulation can emit a
///    triangle `(a, m, b)` with `m` on segment `a-b` within the band. Its three edges cancel
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
/// the clipper's output unchanged, which is certified on its own and at
/// worst reported downstream as an open mesh, not hidden.
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
    // between the other two, where that segment is an invented edge. Its two short sides then lie on the same straight run, and
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

/// Why an authored polygon face could not be triangulated.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum PolygonRefusal {
    /// The outer ring has zero or non-finite area, so it has no plane.
    NoPlane,
    /// A corner lies further from the face plane than the tolerance; the
    /// value is that distance. Only [`triangulate_polygon`] refuses it,
    /// whose callers mean the polygon to be planar (a curve-bounded plane);
    /// an authored face reports its warp instead (#254).
    NotPlanar(Scalar),
    /// The rings do not bound a region in the face's plane: the clipper's
    /// refusal (#260), naming the ring -- crossing or overlapping rings, a
    /// hole outside the outer ring or inside another, too few distinct
    /// corners -- or its certificate's.
    Rings(String),
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
            Self::Rings(reason) => write!(f, "its rings do not bound a region: {reason}"),
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
/// wrong: no plane, a corner off the plane by more than `linear`, or rings
/// the certified clipper refuses (crossing, overlapping, a hole outside).
pub(crate) fn triangulate_polygon(
    rings: &[&[Vec3]],
    linear: Scalar,
) -> Result<Vec<usize>, PolygonRefusal> {
    triangulate_in_fit_plane(rings, linear, true).map(|(indices, _)| indices)
}

/// An authored polygon face's triangles and how far it is warped (#254).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AuthoredTriangulation {
    /// Three indices per triangle, as [`triangulate_polygon`] returns them.
    pub indices: Vec<usize>,
    /// `None` when every corner lies within the linear tolerance of the fit
    /// plane: the face is planar and its triangles are exactly what
    /// [`triangulate_polygon`] returns. Otherwise `Some(w)`, the bound of
    /// the module notes, above the tolerance.
    pub warp: Option<Scalar>,
}

/// [`triangulate_polygon`] for an authored polygon face: a face warped
/// beyond the linear tolerance is triangulated too and its warp reported,
/// so [`PolygonRefusal::NotPlanar`] never comes back. The module notes say
/// what the warp bounds.
///
/// # Refusals
///
/// No plane, or rings that cross (or a hole outside the outer ring) in the
/// projection onto the fit plane, as [`triangulate_polygon`].
pub(crate) fn triangulate_authored_polygon(
    rings: &[&[Vec3]],
    linear: Scalar,
) -> Result<AuthoredTriangulation, PolygonRefusal> {
    triangulate_in_fit_plane(rings, linear, false)
        .map(|(indices, warp)| AuthoredTriangulation { indices, warp })
}

/// The shared triangulation. With `refuse_warp`, a corner off the fit
/// plane by more than `linear` is refused before anything is triangulated,
/// as it always was; without, the face is triangulated and its warp
/// returned.
fn triangulate_in_fit_plane(
    rings: &[&[Vec3]],
    linear: Scalar,
    refuse_warp: bool,
) -> Result<(Vec<usize>, Option<Scalar>), PolygonRefusal> {
    let Some(outer) = rings.first() else {
        return Err(PolygonRefusal::NoPlane);
    };
    let normal = newell_normal(outer.iter().copied());
    let (u, v) = plane_axes(normal).ok_or(PolygonRefusal::NoPlane)?;
    let n = normal.normalize();

    // The plane through the outer ring's centroid, with the Newell normal:
    // Newell's is the least-squares plane direction for a near-planar ring.
    let centroid = outer.iter().copied().fold(Vec3::ZERO, |sum, p| sum + p) / outer.len() as Scalar;
    let spread = CornerSpread::of(rings, centroid, n)?;
    let warp = spread.beyond(linear);
    if refuse_warp && warp.is_some() {
        return Err(PolygonRefusal::NotPlanar(spread.largest));
    }

    let (flat, hole_starts) = project_about(rings, centroid, u, v);

    // The outer ring is counter-clockwise in (u, v) by construction of the
    // axes, and the clipper's triangles always are, so they keep the
    // authored winding. Its certificate replaces the area cross-check that
    // earcut, which returned partial covers on crossing rings, needed.
    let indices = clip_planar_face(&flat, &hole_starts, linear).map_err(|error| {
        PolygonRefusal::Rings(match error {
            GeomError::InvalidInput(reason) | GeomError::Degenerate(reason) => reason,
            other => other.to_string(),
        })
    })?;
    Ok((indices, warp))
}

/// `rings` projected onto the axes `(u, v)` about `origin`, concatenated,
/// with the index at which each hole starts.
fn project_about(
    rings: &[&[Vec3]],
    origin: Vec3,
    u: Vec3,
    v: Vec3,
) -> (Vec<[Scalar; 2]>, Vec<usize>) {
    let mut flat: Vec<[Scalar; 2]> = Vec::new();
    let mut hole_starts: Vec<usize> = Vec::with_capacity(rings.len().saturating_sub(1));
    for (index, ring) in rings.iter().enumerate() {
        if index > 0 {
            hole_starts.push(flat.len());
        }
        flat.extend(ring.iter().map(|&p| {
            let d = p - origin;
            [d.dot(u), d.dot(v)]
        }));
    }
    (flat, hole_starts)
}

/// [`edge_touches`] of an authored polygon face, projected exactly as
/// [`triangulate_authored_polygon`] projects it: the corners its
/// triangulation inserts into one of its own ring edges, by index into the
/// rings concatenated (#265). Empty when the face has no plane.
pub(crate) fn authored_edge_touches(rings: &[&[Vec3]]) -> Vec<EdgeTouch> {
    let Some(outer) = rings.first().filter(|outer| !outer.is_empty()) else {
        return Vec::new();
    };
    let Some((u, v)) = plane_axes(newell_normal(outer.iter().copied())) else {
        return Vec::new();
    };
    let centroid = outer.iter().copied().fold(Vec3::ZERO, |sum, p| sum + p) / outer.len() as Scalar;
    let (flat, hole_starts) = project_about(rings, centroid, u, v);
    edge_touches(&flat, &hole_starts)
}

/// How the corners of a polygon spread about its fit plane (#254, #261):
/// the plane through `centroid` with unit normal `n`, as stored.
#[derive(Debug, Clone, Copy, PartialEq)]
struct CornerSpread {
    /// The largest corner distance from the plane, `max |d_i|`: what the
    /// planarity check compares with the linear tolerance.
    largest: Scalar,
    /// The largest signed distance `max d_i`.
    above: Scalar,
    /// The smallest signed distance `min d_i`.
    below: Scalar,
    /// The largest corner distance from `centroid`, which scales the
    /// rounding of every distance.
    reach: Scalar,
}

impl CornerSpread {
    /// Measure every corner of every ring, holes included. A non-finite
    /// distance has no plane.
    fn of(rings: &[&[Vec3]], centroid: Vec3, n: Vec3) -> Result<Self, PolygonRefusal> {
        let mut spread = Self {
            largest: 0.0,
            above: 0.0,
            below: 0.0,
            reach: 0.0,
        };
        for ring in rings {
            for &p in *ring {
                let offset = p - centroid;
                let distance = offset.dot(n);
                if !distance.is_finite() {
                    return Err(PolygonRefusal::NoPlane);
                }
                spread.largest = spread.largest.max(distance.abs());
                spread.above = spread.above.max(distance);
                spread.below = spread.below.min(distance);
                spread.reach = spread.reach.max(offset.length());
            }
        }
        Ok(spread)
    }

    /// The slab width `max d_i - min d_i` plus its rounding (#261): the
    /// bound the module notes derive.
    ///
    /// Starting `above` and `below` at zero changes nothing: the plane
    /// passes through the outer ring's centroid, so the outer corners'
    /// signed distances sum to zero and already straddle it.
    ///
    /// Each computed distance differs from the true distance to the plane
    /// through `centroid` along `n` (both as stored) by the rounding of a
    /// subtraction, a three-term dot product and the normalisation of `n`,
    /// within 16 ulps of `reach`; the difference of two such distances is
    /// within 32, and the subtraction and the final sum each round by at
    /// most another ulp of `2 reach`. The projection axes `(u, v)` are
    /// orthogonal to `n` up to a few ulps, which tilts the direction the
    /// mesh is a graph along by as much. 64 ulps of `reach` cover all of it.
    fn slab_width(self) -> Scalar {
        (self.above - self.below) + 64.0 * f64::EPSILON * self.reach
    }

    /// The warp a caller reports: `None` when every corner lies within
    /// `linear` of the plane (the face counts as planar, as it always
    /// did), else the [`Self::slab_width`].
    fn beyond(self, linear: Scalar) -> Option<Scalar> {
        (self.largest > linear).then(|| self.slab_width())
    }
}

/// How far a polygon given as an outer ring plus holes is warped, as
/// [`triangulate_authored_polygon`] reports it: `None` when every corner
/// lies within `linear` of the fit plane, else `Some` slab width of the
/// module notes (#261). For a caller that triangulates the polygon itself
/// projected along the same Newell normal, such as a B-rep face that
/// declares no surface (#257).
///
/// # Refusals
///
/// [`PolygonRefusal::NoPlane`] when the outer ring encloses no area or a
/// corner is not finite.
pub(crate) fn polygon_warp(
    rings: &[&[Vec3]],
    linear: Scalar,
) -> Result<Option<Scalar>, PolygonRefusal> {
    let Some(outer) = rings.first().filter(|outer| !outer.is_empty()) else {
        return Err(PolygonRefusal::NoPlane);
    };
    let normal = newell_normal(outer.iter().copied());
    plane_axes(normal).ok_or(PolygonRefusal::NoPlane)?;
    let n = normal.normalize();
    let centroid = outer.iter().copied().fold(Vec3::ZERO, |sum, p| sum + p) / outer.len() as Scalar;
    let spread = CornerSpread::of(rings, centroid, n)?;
    Ok(spread.beyond(linear))
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

    /// Twice the signed area of a projected triangle.
    fn doubled_area(a: [Scalar; 2], b: [Scalar; 2], c: [Scalar; 2]) -> Scalar {
        (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
    }

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
    /// edge is used at most twice and the authored ring edges exactly once,
    /// on the certified clipper's triangles as on earcut's (#260).
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
        let indices = clip_planar_face(&flat, &[], 1.0e-3).expect("the stringer triangulates");
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
        let indices = clip_planar_face(&flat, &[], 1.0e-3).expect("one triangle");
        assert_eq!(indices.len(), 3, "the face is not a sliver: {indices:?}");
    }

    /// A planar face (within the tolerance) triangulates bit for bit as
    /// [`triangulate_polygon`] does and reports no warp (#254); a warped one
    /// is still refused by [`triangulate_polygon`], whose callers mean a
    /// plane, and triangulated with its warp by the authored path.
    #[test]
    fn the_authored_path_is_the_planar_path_on_planar_faces() {
        let holed: [&[Vec3]; 2] = [
            &[
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(4.0, 0.0, 1.0),
                Vec3::new(4.0, 3.0, 1.0),
                Vec3::new(2.0, 1.5, 0.5),
                Vec3::new(0.0, 3.0, 0.0),
            ],
            &[
                Vec3::new(0.5, 0.5, 0.125),
                Vec3::new(0.5, 1.0, 0.125),
                Vec3::new(1.0, 1.0, 0.25),
                Vec3::new(1.0, 0.5, 0.25),
            ],
        ];
        let noisy: [&[Vec3]; 1] = [&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 9e-4),
            Vec3::new(0.0, 1.0, 0.0),
        ]];
        for rings in [&holed[..], &noisy[..]] {
            let planar = triangulate_polygon(rings, 1e-3).unwrap();
            let authored = triangulate_authored_polygon(rings, 1e-3).unwrap();
            assert_eq!(authored.indices, planar);
            assert_eq!(authored.warp, None);
        }
        let warped: [&[Vec3]; 1] = [&[
            Vec3::new(0.0, 0.0, 0.05),
            Vec3::new(1.0, 0.0, -0.05),
            Vec3::new(1.0, 1.0, 0.05),
            Vec3::new(0.0, 1.0, -0.05),
        ]];
        assert!(matches!(
            triangulate_polygon(&warped, 1e-3),
            Err(PolygonRefusal::NotPlanar(d)) if (d - 0.05).abs() < 1e-12
        ));
        let authored = triangulate_authored_polygon(&warped, 1e-3).unwrap();
        assert_eq!(authored.indices.len(), 6);
        // The +-5 cm saddle's slab is 10 cm wide (#261).
        let warp = authored.warp.unwrap();
        assert!((0.1..0.1 + 1e-12).contains(&warp), "{warp}");
    }

    /// The slab width (#261): a saddle reports its full spread, a lifted
    /// corner the spread about its fit plane, and holes count.
    #[test]
    fn the_warp_is_the_slab_width_of_every_corner() {
        let h = 0.05;
        let lifted: [&[Vec3]; 1] = [&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, h),
            Vec3::new(0.0, 1.0, 0.0),
        ]];
        let saddle: [&[Vec3]; 1] = [&[
            Vec3::new(0.0, 0.0, h),
            Vec3::new(1.0, 0.0, -h),
            Vec3::new(1.0, 1.0, h),
            Vec3::new(0.0, 1.0, -h),
        ]];
        // A flat outer ring whose hole has one corner below it.
        let holed: [&[Vec3]; 2] = [
            &[
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(4.0, 0.0, 0.0),
                Vec3::new(4.0, 4.0, 0.0),
                Vec3::new(0.0, 4.0, 0.0),
            ],
            &[
                Vec3::new(1.0, 1.0, 0.0),
                Vec3::new(1.0, 2.0, 0.0),
                Vec3::new(2.0, 2.0, -h),
                Vec3::new(2.0, 1.0, 0.0),
            ],
        ];
        // The lifted quad is `z = h x y`: its fit plane takes the linear
        // part, leaving a saddle of +-h/4 about it, so the slab is h/2
        // wide along that plane's normal (Newell normal (-h, -h, 2)), the
        // gap between the two diagonal triangulations at the centre.
        let tilt = 2.0 / (4.0 + 2.0 * h * h).sqrt();
        for (rings, want) in [
            (&lifted[..], 0.5 * h * tilt),
            (&saddle[..], 2.0 * h),
            (&holed[..], h),
        ] {
            let warp = triangulate_authored_polygon(rings, 1e-3)
                .unwrap()
                .warp
                .unwrap();
            assert!((want..want + 1e-12).contains(&warp), "{warp}, want {want}");
            // The B-rep path measures the same way (#257).
            assert_eq!(polygon_warp(rings, 1e-3), Ok(Some(warp)));
        }
        let mut planar = lifted[0].to_vec();
        planar[2].z = 0.0;
        assert_eq!(polygon_warp(&[&planar], 1e-3), Ok(None));
        assert_eq!(polygon_warp(&[], 1e-3), Err(PolygonRefusal::NoPlane));
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
