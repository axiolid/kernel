//! Simple-polygon triangulation by ear clipping, with hole support.
//!
//! Orientation is decided by the certified `orient2d` predicate rather than a
//! raw sign test, so a near-degenerate vertex cannot silently flip a triangle.

use axiolid_contracts::{GeomError, GeomResult, Sign};
use axiolid_core::{Point2, Scalar};

use crate::orient2d;

/// Twice the signed area of a closed ring. Positive means counter-clockwise.
pub fn signed_area2(ring: &[Point2]) -> Scalar {
    let mut acc = 0.0;
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        acc += a.x * b.y - b.x * a.y;
    }
    acc
}

/// Orientation of a ring, or `None` when it encloses no certifiable area.
pub fn ring_orientation(ring: &[Point2]) -> Option<Sign> {
    if ring.len() < 3 {
        return None;
    }
    match signed_area2(ring) {
        a if a > 0.0 => Some(Sign::Positive),
        a if a < 0.0 => Some(Sign::Negative),
        _ => None,
    }
}

/// Certified orientation sign of a finite triple.
///
/// `orient2d` is exact for every finite input and only a non-finite
/// coordinate leaves it uncertain, so [`triangulate_simple`] refuses such a
/// ring before any triple reaches this helper.
fn sign_of(a: Point2, b: Point2, c: Point2) -> Sign {
    orient2d(a, b, c)
        .sign()
        .expect("orient2d is certain for finite coordinates, which the caller checked")
}

/// Whether `p` lies in or on the boundary of counter-clockwise triangle
/// `a,b,c`, using certified signs.
///
/// A non-corner vertex on an ear boundary blocks that ear. Otherwise the clip
/// can leave a negatively wound remainder whose signed area merely cancels the
/// outside area it emitted.
fn blocks_ear(a: Point2, b: Point2, c: Point2, p: Point2) -> bool {
    let s1 = sign_of(a, b, p);
    let s2 = sign_of(b, c, p);
    let s3 = sign_of(c, a, p);
    s1 != Sign::Negative && s2 != Sign::Negative && s3 != Sign::Negative
}

/// Whether vertex `i` of `ring` is a clippable ear of a CCW polygon.
fn is_ear(ring: &[Point2], indices: &[usize], at: usize) -> bool {
    let n = indices.len();
    let (ia, ib, ic) = (
        indices[(at + n - 1) % n],
        indices[at],
        indices[(at + 1) % n],
    );
    let (a, b, c) = (ring[ia], ring[ib], ring[ic]);

    // Reflex or collinear vertices are never ears. Collinear is excluded
    // because a zero-area ear adds a degenerate triangle to the output.
    if sign_of(a, b, c) != Sign::Positive {
        return false;
    }

    // No other vertex may fall inside the candidate ear.
    //
    // Bridged rings contain DUPLICATE vertices by construction: a bridge visits
    // the same two points twice. Comparing by index alone would treat a
    // duplicate of an ear corner as an intruder and reject every ear, so
    // coincident positions are skipped as well.
    !indices
        .iter()
        .filter(|&&idx| idx != ia && idx != ib && idx != ic)
        .map(|&idx| ring[idx])
        .filter(|q| *q != a && *q != b && *q != c)
        .any(|q| blocks_ear(a, b, c, q))
}

/// Triangulate a simple CCW polygon by ear clipping.
///
/// Returns index triples into `ring`. Fails rather than emitting a partial fan:
/// a caller that receives triangles must be able to trust they cover the input.
///
/// # Errors
///
/// [`GeomError::InvalidInput`] for a NaN or infinite coordinate, which no
/// orientation predicate can decide; [`GeomError::Degenerate`] for a ring that
/// is too short, self-intersecting or without area.
pub fn triangulate_simple(ring: &[Point2]) -> GeomResult<Vec<[u32; 3]>> {
    if ring.iter().any(|p| !(p.x.is_finite() && p.y.is_finite())) {
        return Err(GeomError::InvalidInput(
            "ring has a non-finite coordinate".to_owned(),
        ));
    }
    if ring.len() < 3 {
        return Err(GeomError::Degenerate(format!(
            "ring has {} vertices, need at least 3",
            ring.len()
        )));
    }
    let mut indices: Vec<usize> = (0..ring.len()).collect();
    let mut out = Vec::with_capacity(ring.len().saturating_sub(2));

    // Each successful clip removes one vertex, so the loop is bounded. The
    // `guard` counts consecutive failures: a full pass with no ear found means
    // the polygon is not simple, which is a data fault, not a retry case.
    let mut at = 0usize;
    let mut guard = 0usize;
    while indices.len() > 3 {
        if guard > indices.len() {
            return Err(GeomError::Degenerate(
                "no ear found; polygon is self-intersecting or degenerate".to_owned(),
            ));
        }
        if is_ear(ring, &indices, at % indices.len()) {
            let n = indices.len();
            let cur = at % n;
            out.push([
                indices[(cur + n - 1) % n] as u32,
                indices[cur] as u32,
                indices[(cur + 1) % n] as u32,
            ]);
            indices.remove(cur);
            at = cur;
            guard = 0;
        } else {
            at += 1;
            guard += 1;
        }
    }
    // The final triangle must still enclose area; a collinear remainder means
    // the whole ring was degenerate and no valid fan exists.
    let (a, b, c) = (ring[indices[0]], ring[indices[1]], ring[indices[2]]);
    if sign_of(a, b, c) != Sign::Positive {
        return Err(GeomError::Degenerate(
            "final triangle is not positively wound; polygon is self-intersecting or degenerate"
                .to_owned(),
        ));
    }
    out.push([indices[0] as u32, indices[1] as u32, indices[2] as u32]);
    Ok(out)
}
