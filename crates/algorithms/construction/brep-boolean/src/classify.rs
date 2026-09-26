//! Classification: is a split region inside the other operand?
//! (ADR 0075, step 4.)
//!
//! A region lies wholly inside or wholly outside the other solid -- its
//! boundary is where the two meet -- so one point strictly inside it
//! decides. That point is classified by ray parity: a ray from it is
//! intersected exactly with every face's surface
//! (`exact_curve_surface_intersection`), and a hit counts when the face's
//! certified domain contains it. A ray that grazes a surface (a touching
//! hit), passes too close to a face boundary to decide, or starts on the
//! other solid's boundary is not trusted: the next of a fixed set of
//! directions is tried, and if all fail the region is refused.

use axiolid_brep::ExactBRep;
use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{Curve3, Line3};
use axiolid_evaluate::curve::{derivative2, evaluate2};
use axiolid_evaluate::surface::{evaluate, invert};
use axiolid_measure::FaceDomain;
use axiolid_nurbs::{exact_curve_surface_intersection, ExactCurveIntersection};
use axiolid_surface::Surface;

use crate::split::{Piece, Region};
use crate::BooleanError;

/// Directions tried in turn, none axis-aligned.
const DIRECTIONS: [[Scalar; 3]; 6] = [
    [0.573, 0.341, 0.745],
    [-0.482, 0.613, 0.625],
    [0.297, -0.861, 0.413],
    [-0.664, -0.229, -0.712],
    [0.812, 0.463, -0.355],
    [0.123, 0.951, -0.284],
];

/// A solid, prepared for point classification.
pub(crate) struct Solid<'a> {
    brep: &'a ExactBRep,
    domains: Vec<FaceDomain<'a>>,
}

impl<'a> Solid<'a> {
    pub(crate) fn new(brep: &'a ExactBRep, tolerance: Tolerance) -> Result<Self, BooleanError> {
        let topology = brep.topology();
        let mut domains = Vec::with_capacity(topology.faces().len());
        for index in 0..topology.faces().len() {
            let face = topology
                .face_id_at(index)
                .ok_or(BooleanError::DanglingReference)?;
            domains.push(
                FaceDomain::new(brep, face, tolerance)
                    .map_err(BooleanError::Measure)?
                    .ok_or(BooleanError::UnsupportedTrim)?,
            );
        }
        Ok(Self { brep, domains })
    }

    fn surface(&self, face: usize) -> Result<&'a Surface, BooleanError> {
        self.brep.topology().faces()[face]
            .surface
            .and_then(|id| self.brep.surfaces().get(id.index()))
            .ok_or(BooleanError::DanglingReference)
    }

    /// Whether `point` lies inside the solid.
    pub(crate) fn contains(
        &self,
        point: Point3,
        tolerance: Tolerance,
    ) -> Result<bool, BooleanError> {
        'direction: for direction in DIRECTIONS {
            let direction = Vec3::from_array(direction).normalize();
            let ray = Curve3::Line(Line3 {
                origin: point,
                direction,
            });
            let mut crossings = 0usize;
            for face in 0..self.domains.len() {
                let surface = self.surface(face)?;
                let hits = match exact_curve_surface_intersection(&ray, surface) {
                    Ok(ExactCurveIntersection::Points(hits)) => hits,
                    // The ray lies in the surface: try another.
                    Ok(_) => continue 'direction,
                    Err(_) => return Err(BooleanError::UnsupportedTrim),
                };
                for hit in hits {
                    let t = hit.parameter.approx();
                    let at_start = t.abs() <= tolerance.linear().max(1e-9);
                    if t < 0.0 && !at_start {
                        continue;
                    }
                    let (u, v) = invert(surface, hit.point, tolerance)
                        .map_err(|_| BooleanError::Evaluation)?;
                    let on_face = self.domains[face]
                        .contains(Point2::new(u, v))
                        .map_err(BooleanError::Measure)?;
                    if at_start {
                        // The point sits on this face's support. Outside the
                        // face that is no crossing at all; on the face it is
                        // on the boundary, which no direction can settle.
                        match on_face {
                            Some(false) => continue,
                            _ => return Err(BooleanError::Undecided),
                        }
                    }
                    if hit.multiplicity != 1 {
                        continue 'direction;
                    }
                    match on_face {
                        Some(true) => crossings += 1,
                        Some(false) => {}
                        None => continue 'direction,
                    }
                }
            }
            return Ok(crossings % 2 == 1);
        }
        Err(BooleanError::Undecided)
    }
}

/// A point strictly inside a region, on its face's surface.
///
/// Stepped inwards from the middle of the longest boundary piece, a little
/// way to the loop's left (regions run anticlockwise), then halved until it
/// falls inside the region's sampled outline and outside its holes.
pub(crate) fn interior_point(region: &Region, surface: &Surface) -> Result<Point3, BooleanError> {
    let outline = |pieces: &[Piece]| -> Result<Vec<Point2>, BooleanError> {
        let mut out = Vec::new();
        for piece in pieces {
            for i in 0..64 {
                let p =
                    piece.pspan.start + (piece.pspan.end - piece.pspan.start) * i as Scalar / 64.0;
                out.push(evaluate2(&piece.pcurve, p).map_err(|_| BooleanError::Evaluation)?);
            }
        }
        Ok(out)
    };
    let outer = outline(&region.outer)?;
    let holes = region
        .holes
        .iter()
        .map(|h| outline(h))
        .collect::<Result<Vec<_>, _>>()?;
    let (mut lo, mut hi) = (outer[0], outer[0]);
    for p in &outer {
        lo = lo.min(*p);
        hi = hi.max(*p);
    }
    let size = (hi - lo).length();
    let inside = |p: Point2| {
        crate::split::inside_polygon(&outer, p)
            && holes.iter().all(|h| !crate::split::inside_polygon(h, p))
    };
    for piece in region.outer.iter().chain(region.holes.iter().flatten()) {
        let mid = 0.5 * (piece.pspan.start + piece.pspan.end);
        let at = evaluate2(&piece.pcurve, mid).map_err(|_| BooleanError::Evaluation)?;
        let mut tangent = derivative2(&piece.pcurve, mid).map_err(|_| BooleanError::Evaluation)?;
        if piece.pspan.end < piece.pspan.start {
            tangent = -tangent;
        }
        let length = tangent.length();
        if length == 0.0 {
            continue;
        }
        let left = Vec2::new(-tangent.y, tangent.x) / length;
        let mut step = 0.25 * size;
        for _ in 0..40 {
            let probe = at + left * step;
            if inside(probe) {
                return evaluate(surface, probe.x, probe.y).map_err(|_| BooleanError::Evaluation);
            }
            step *= 0.5;
        }
    }
    Err(BooleanError::Undecided)
}
