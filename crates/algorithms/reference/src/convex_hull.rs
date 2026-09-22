//! Deterministic 2D convex hulls and oriented bounding rectangles.
use crate::orient2d;
use axiolid_contracts::{GeomError, GeomResult, Sign};
use axiolid_core::{Point2, Rectangle2, Scalar};

/// The oriented rectangle this module used to define itself.
///
/// Kept as an alias rather than a distinct type: it stored four corners plus a
/// cached area and side pair, which is the same rectangle `Rectangle2` models
/// as an origin and two edge vectors, and having two spellings of one concept
/// forced every caller to convert. The edge-vector form is also the one that
/// cannot drift -- four independently stored corners can be edited into a
/// non-parallelogram, and the cached area can disagree with them.
pub type OrientedRectangle2 = Rectangle2;

/// Side lengths of an oriented rectangle, shortest edge first.
///
/// A free function rather than an inherent method because `Rectangle2` lives
/// in `axiolid-core`, which holds data and no algorithms. Ordering is fixed so
/// the result does not depend on which hull edge the caliper happened to stop
/// on.
#[must_use]
pub fn side_lengths(rectangle: &Rectangle2) -> [Scalar; 2] {
    let mut sides = [rectangle.x.length(), rectangle.y.length()];
    sides.sort_by(|left, right| left.total_cmp(right));
    sides
}

pub fn strict_convex_hull(points: &[Point2]) -> GeomResult<Vec<usize>> {
    for (index, point) in points.iter().enumerate() {
        if !point.is_finite() {
            return Err(GeomError::InvalidInput(format!(
                "point {index} is not finite"
            )));
        }
    }
    let mut ordered: Vec<usize> = (0..points.len()).collect();
    ordered.sort_by(|&a, &b| {
        points[a]
            .x
            .total_cmp(&points[b].x)
            .then_with(|| points[a].y.total_cmp(&points[b].y))
            .then(a.cmp(&b))
    });
    ordered.dedup_by(|a, b| points[*a] == points[*b]);
    if ordered.len() < 3 {
        return Err(GeomError::Degenerate("need three distinct points".into()));
    }
    let mut lower = Vec::new();
    for &index in &ordered {
        push_strict(&mut lower, index, points);
    }
    let mut upper = Vec::new();
    for &index in ordered.iter().rev() {
        push_strict(&mut upper, index, points);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    if lower.len() < 3 {
        return Err(GeomError::Degenerate("points are collinear".into()));
    }
    Ok(lower)
}

pub fn minimum_area_rectangle(points: &[Point2]) -> GeomResult<OrientedRectangle2> {
    let hull = strict_convex_hull(points)?;
    let mut best: Option<Rectangle2> = None;
    // Tracked beside the rectangle because `Rectangle2` deliberately stores no
    // cached area: recomputing it per candidate is a multiply, and a cached
    // field is one more thing that can disagree with the geometry.
    let mut best_area: Option<Scalar> = None;
    for i in 0..hull.len() {
        let a = points[hull[i]];
        let b = points[hull[(i + 1) % hull.len()]];
        let edge = b - a;
        let width = edge.length();
        let u = edge / width;
        let v = Point2::new(-u.y, u.x);
        let (mut ulo, mut uhi, mut vlo, mut vhi) = (
            Scalar::INFINITY,
            Scalar::NEG_INFINITY,
            Scalar::INFINITY,
            Scalar::NEG_INFINITY,
        );
        for &index in &hull {
            let point = points[index];
            let pu = point.dot(u);
            let pv = point.dot(v);
            ulo = ulo.min(pu);
            uhi = uhi.max(pu);
            vlo = vlo.min(pv);
            vhi = vhi.max(pv);
        }
        let sides = [uhi - ulo, vhi - vlo];
        let area = sides[0] * sides[1];
        // Origin plus two edge vectors, rather than four corners: the
        // parallelogram property then holds by construction instead of being
        // an invariant four separately stored points could violate.
        let rectangle = Rectangle2::new(u * ulo + v * vlo, u * sides[0], v * sides[1]);
        if best_area.is_none_or(|current| area < current) {
            best_area = Some(area);
            best = Some(rectangle);
        }
    }
    Ok(best.expect("non-empty strict hull has an edge"))
}

fn push_strict(hull: &mut Vec<usize>, index: usize, points: &[Point2]) {
    while hull.len() >= 2 {
        let n = hull.len();
        if sign(orient2d(
            points[hull[n - 2]],
            points[hull[n - 1]],
            points[index],
        )) == Sign::Positive
        {
            break;
        }
        hull.pop();
    }
    hull.push(index);
}
fn sign(value: axiolid_contracts::Certified) -> Sign {
    value.sign().expect("orient2d is total")
}
