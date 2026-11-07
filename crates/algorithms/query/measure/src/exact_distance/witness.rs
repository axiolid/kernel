//! Certified witnesses beyond patch centres and edge midpoints, which close
//! a gap the lower bound already has right (#318).
//!
//! A lower bound that does not depend on the size of a patch -- along a
//! planar patch's own normal, against a parallel planar patch or a line edge
//! -- is exact on a whole region, so refinement never raises it. What never
//! closes is the upper bound: patch-centre witnesses line up by chance. A
//! boundary point nearest the other element's witness does close it:
//!
//! - the foot, on a planar face, of the other element's witness;
//! - the closest points of two line spans;
//! - for a point query, the foot of the point on a planar face, or its
//!   closest point on a line span.
//!
//! Every point is evaluated on its own boundary as any witness is: on the
//! edge's curve, or on the face's surface at parameters that lie in a patch
//! already certified inside the face or that the face's domain certifies. A
//! foot beyond the face is no witness.

use axiolid_core::{Frame3, Point2, Point3, Scalar};
use axiolid_curve::Curve3;
use axiolid_surface::Surface;

use super::touch::{line_span, on_curve, segment_parameters};
use super::{point_on, surface_of, Element, Shape, Side};
use crate::exact::ExactMeasureError;

/// A boundary point of each element, and the shapes holding them.
type Witnesses = (Point3, Point3, super::Pair);

/// The unclamped `(u, v)` of the foot of `p` on the plane of `frame`, or
/// `None` for a degenerate frame.
pub(super) fn plane_parameters(frame: &Frame3, p: Point3) -> Option<Point2> {
    let (gxx, gxy, gyy) = (
        frame.x.dot(frame.x),
        frame.x.dot(frame.y),
        frame.y.dot(frame.y),
    );
    let gram = gxx * gyy - gxy * gxy;
    if gram.is_nan() || gram <= 0.0 {
        return None;
    }
    let r = p - frame.origin;
    let (rx, ry) = (frame.x.dot(r), frame.y.dot(r));
    Some(Point2::new(
        (gyy * rx - gxy * ry) / gram,
        (gxx * ry - gxy * rx) / gram,
    ))
}

/// The best certified pair of boundary points of two elements, and the
/// shapes holding them: the closest points of two line spans, the foot of
/// `b`'s witness on `a`'s planar face, the foot of `a`'s on `b`'s -- the
/// nearest of them, the first on a tie.
pub(super) fn pair(
    side_a: &Side<'_>,
    a: &Element,
    side_b: &Side<'_>,
    b: &Element,
) -> Result<Option<Witnesses>, ExactMeasureError> {
    let mut best: Option<(Scalar, Point3, Point3, (Shape, Shape))> = None;
    let mut offer = |p: Point3, q: Point3, on: (Shape, Shape)| {
        let d = (p - q).length();
        if d.is_finite() && best.is_none_or(|(current, ..)| d < current) {
            best = Some((d, p, q, on));
        }
    };
    if let Some((p, q, on)) = spans(side_a, a, side_b, b)? {
        offer(p, q, on);
    }
    if let Some(w) = b.witness {
        if let Some((p, shape)) = toward_face(w, side_a, a)? {
            offer(p, w, (shape, b.shape));
        }
    }
    if let Some(w) = a.witness {
        if let Some((q, shape)) = toward_face(w, side_b, b)? {
            offer(w, q, (a.shape, shape));
        }
    }
    Ok(best.map(|(_, p, q, on)| (p, q, on)))
}

/// A certified boundary point of `element` nearest `p` by its family: the
/// foot of `p` on a planar face, its closest point on a line span.
pub(crate) fn toward(
    p: Point3,
    side: &Side<'_>,
    element: &Element,
) -> Result<Option<(Point3, Shape)>, ExactMeasureError> {
    if let Some(found) = toward_face(p, side, element)? {
        return Ok(Some(found));
    }
    let Some((curve, t0, t1)) = line_span(side, element)? else {
        return Ok(None);
    };
    let (Curve3::Line(line), Shape::Edge { edge, .. }) = (curve, element.shape) else {
        return Ok(None);
    };
    let along = line.direction.dot(line.direction);
    let t = ((p - line.origin).dot(line.direction) / along).clamp(t0, t1);
    if !(along > 0.0 && t.is_finite()) {
        return Ok(None);
    }
    Ok(Some((
        on_curve(curve, t)?,
        Shape::Edge { edge, t0: t, t1: t },
    )))
}

/// The closest points of two line spans, each evaluated on its own curve.
fn spans(
    side_a: &Side<'_>,
    a: &Element,
    side_b: &Side<'_>,
    b: &Element,
) -> Result<Option<Witnesses>, ExactMeasureError> {
    let (Some((ca, a0, a1)), Some((cb, b0, b1))) = (line_span(side_a, a)?, line_span(side_b, b)?)
    else {
        return Ok(None);
    };
    let (Shape::Edge { edge: ea, .. }, Shape::Edge { edge: eb, .. }) = (a.shape, b.shape) else {
        return Ok(None);
    };
    let (p0, p1) = (on_curve(ca, a0)?, on_curve(ca, a1)?);
    let (q0, q1) = (on_curve(cb, b0)?, on_curve(cb, b1)?);
    let (s, t) = segment_parameters(p0, p1, q0, q1);
    if !(s.is_finite() && t.is_finite()) {
        return Ok(None);
    }
    let (ta, tb) = (a0 + (a1 - a0) * s, b0 + (b1 - b0) * t);
    let p = on_curve(ca, ta)?;
    let q = on_curve(cb, tb)?;
    Ok(Some((
        p,
        q,
        (
            Shape::Edge {
                edge: ea,
                t0: ta,
                t1: ta,
            },
            Shape::Edge {
                edge: eb,
                t0: tb,
                t1: tb,
            },
        ),
    )))
}

/// The foot of `p` on a planar face patch, when the patch lies in the face
/// or the face's domain certifies the foot.
fn toward_face(
    p: Point3,
    side: &Side<'_>,
    element: &Element,
) -> Result<Option<(Point3, Shape)>, ExactMeasureError> {
    let Shape::Face {
        face,
        lo,
        hi,
        inside,
    } = element.shape
    else {
        return Ok(None);
    };
    let surface = surface_of(side.brep, side.brep.topology().faces()[face].surface)?;
    let Surface::Plane(plane) = surface else {
        return Ok(None);
    };
    let Some(at) = plane_parameters(&plane.frame, p) else {
        return Ok(None);
    };
    if !(at.x.is_finite() && at.y.is_finite()) {
        return Ok(None);
    }
    let in_patch = at.x >= lo.x.min(hi.x)
        && at.x <= lo.x.max(hi.x)
        && at.y >= lo.y.min(hi.y)
        && at.y <= lo.y.max(hi.y);
    let certified = (inside && in_patch)
        || match &side.domains[face] {
            Some(domain) => domain.contains(at)? == Some(true),
            None => false,
        };
    if !certified {
        return Ok(None);
    }
    Ok(Some((
        point_on(surface, at)?,
        Shape::Face {
            face,
            lo: at,
            hi: at,
            inside: true,
        },
    )))
}
