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
//!
//! Where a ray starts is read in two ways ([`Start`]). First within the
//! linear tolerance: a point that close to any face's support is not
//! decided, since a reading within tolerance may have moved that face onto
//! the point's own. Where no point of a region decides so, the point is
//! read exactly (#332, #327): it lies on a support only within the
//! rounding floor, but for the faces read as its own support, which keep
//! the tolerance. A point micrometres off a surface its own support is not
//! read as one with (a box's face tangent to a curved wall's inner
//! cylinder, away from the contact ruling; a side face of one opening just
//! inside its neighbour, their planes turned by more than the angular
//! tolerance) is then off it, and the ray's crossing of it is an exact
//! simple root that counts. Such a region lies within the tolerance of
//! that surface wherever its points were taken, so either side is within
//! the perturbation the readings already stand for.

use axiolid_brep::ExactBRep;
use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{Curve3, Line3};
use axiolid_evaluate::curve::{derivative2, evaluate2};
use axiolid_evaluate::surface::{evaluate, locate, normal};
use axiolid_measure::FaceDomain;
use axiolid_nurbs::{exact_curve_surface_intersection, ExactCurveIntersection};
use axiolid_surface::Surface;
use axiolid_topology::Orientation;

use crate::bounds::{face_box, Aabb};
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
    /// Per face, whether its loops wind clockwise in its parameters, so it
    /// faces opposite to its orientation flag (`split::face_against`).
    against: Vec<bool>,
    /// A sound box around each face, enlarged by the tolerance: a ray that
    /// misses it cannot cross the face (#228).
    boxes: Vec<Option<Aabb>>,
    /// The box around all of them, when every face has one.
    bounds: Option<Aabb>,
}

impl<'a> Solid<'a> {
    pub(crate) fn new(brep: &'a ExactBRep, tolerance: Tolerance) -> Result<Self, BooleanError> {
        let topology = brep.topology();
        let mut domains = Vec::with_capacity(topology.faces().len());
        let mut against = Vec::with_capacity(topology.faces().len());
        for index in 0..topology.faces().len() {
            let face = topology
                .face_id_at(index)
                .ok_or(BooleanError::DanglingReference)?;
            against.push(crate::split::face_against(brep, face)?);
            domains.push(
                FaceDomain::new(brep, face, tolerance)
                    .map_err(BooleanError::Measure)?
                    .ok_or(BooleanError::UnsupportedTrim)?,
            );
        }
        let mut boxes = Vec::with_capacity(domains.len());
        for (index, domain) in domains.iter().enumerate() {
            let surface = topology.faces()[index]
                .surface
                .and_then(|id| brep.surfaces().get(id.index()))
                .ok_or(BooleanError::DanglingReference)?;
            let (lo, hi) = domain.bounds();
            boxes.push(face_box(surface, lo, hi, tolerance));
        }
        let bounds = boxes
            .iter()
            .try_fold(None::<Aabb>, |all, b| {
                let b = (*b)?;
                Some(Some(all.map_or(b, |a| a.union(&b))))
            })
            .flatten();
        Ok(Self {
            brep,
            domains,
            against,
            boxes,
            bounds,
        })
    }

    fn surface(&self, face: usize) -> Result<&'a Surface, BooleanError> {
        self.brep.topology().faces()[face]
            .surface
            .and_then(|id| self.brep.surfaces().get(id.index()))
            .ok_or(BooleanError::DanglingReference)
    }

    /// The first of `candidates` (faces on the point's own support) whose
    /// face contains `point`, with that face's outward normal there.
    ///
    /// The outward normal is the surface's turned by the face's flag and
    /// by its winding: a face whose loops wind clockwise in its parameters
    /// faces against its flag, as a split region of it does (#324: a
    /// niche's inner face, an extrusion's face on its profile's inner arc,
    /// read with its normal turned round where it lies on a curved wall's
    /// inner face, so the difference kept the wall's patch on it).
    pub(crate) fn on_face(
        &self,
        point: Point3,
        candidates: &[usize],
        tolerance: Tolerance,
    ) -> Result<Option<Vec3>, BooleanError> {
        for &face in candidates {
            if self.boxes[face]
                .as_ref()
                .is_some_and(|b| !b.contains(point))
            {
                continue;
            }
            let surface = self.surface(face)?;
            let (u, v) = locate(surface, point, crate::report::floored(tolerance))
                .map_err(|_| BooleanError::Evaluation)?;
            match self.domains[face]
                .contains(Point2::new(u, v))
                .map_err(BooleanError::Measure)?
            {
                Some(true) => {
                    let n = normal(surface, u, v).map_err(|_| BooleanError::Evaluation)?;
                    let sign = match (
                        self.brep.topology().faces()[face].orientation,
                        self.against[face],
                    ) {
                        (Orientation::Forward, false) | (Orientation::Reversed, true) => 1.0,
                        _ => -1.0,
                    };
                    return Ok(Some(n * sign));
                }
                Some(false) => {}
                None => return Err(BooleanError::Undecided),
            }
        }
        Ok(None)
    }

    /// Whether `point` lies inside the solid, its ray's start read as
    /// `start` says.
    pub(crate) fn contains(
        &self,
        point: Point3,
        start: Start<'_>,
        tolerance: Tolerance,
    ) -> Result<bool, BooleanError> {
        // Outside the box around every face, the point is outside the solid.
        if self.bounds.as_ref().is_some_and(|b| !b.contains(point)) {
            return Ok(false);
        }
        let start_slack = tolerance.linear().max(crate::report::rounding());
        'direction: for direction in DIRECTIONS {
            let direction = Vec3::from_array(direction).normalize();
            let ray = Curve3::Line(Line3 {
                origin: point,
                direction,
            });
            let mut crossings = 0usize;
            for face in 0..self.domains.len() {
                if self.boxes[face]
                    .as_ref()
                    .is_some_and(|b| !b.met_by_ray(point, direction, start_slack))
                {
                    continue;
                }
                let surface = self.surface(face)?;
                let hits = match exact_curve_surface_intersection(&ray, surface) {
                    Ok(ExactCurveIntersection::Points(hits)) => hits,
                    // The ray lies in the surface: try another.
                    Ok(_) => continue 'direction,
                    Err(_) => return Err(BooleanError::UnsupportedTrim),
                };
                for hit in hits {
                    let t = hit.parameter.approx();
                    let near = t.abs() <= start_slack;
                    if t < 0.0 && !near {
                        continue;
                    }
                    let (u, v) = locate(surface, hit.point, crate::report::floored(tolerance))
                        .map_err(|_| BooleanError::Evaluation)?;
                    // A hit at the ray's start: the point lies on this
                    // face's support. Not a reading: it only skips a hit
                    // outside the face or refuses the point.
                    let at_start = near
                        && match start {
                            Start::Tolerance => true,
                            Start::Exact(own) if own.contains(&face) => true,
                            // Off the support by more than the rounding of
                            // one exact point, measured along the normal,
                            // not along a ray that may graze it.
                            Start::Exact(_) => {
                                let n =
                                    normal(surface, u, v).map_err(|_| BooleanError::Evaluation)?;
                                (t * direction.dot(n.normalize())).abs()
                                    <= crate::report::rounding()
                            }
                        };
                    if t < 0.0 && !at_start {
                        continue;
                    }
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

/// How a ray's start is read in [`Solid::contains`].
#[derive(Debug, Clone, Copy)]
pub(crate) enum Start<'c> {
    /// A hit within the linear tolerance of the start puts the point on
    /// that face's support, which no direction settles.
    Tolerance,
    /// Only a hit within the rounding floor does (#332), but on the faces
    /// listed (read as the point's own support within tolerance), which
    /// keep the linear tolerance: the decision is exact for the operands
    /// as given.
    Exact(&'c [usize]),
}

/// Points strictly inside a region, on its face's surface, best first.
///
/// Each is stepped inwards from a point on a boundary piece (at a few
/// fractions along it, off-centre so symmetric inputs do not put every
/// candidate on a line of contact), a little way to the loop's left
/// (regions run anticlockwise), then halved until it falls inside the
/// region's sampled outline and outside its holes. A caller whose
/// classification of one candidate is undecided tries the next.
pub(crate) fn interior_points(
    region: &Region,
    surface: &Surface,
) -> Result<Vec<Point3>, BooleanError> {
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
    let mut out = Vec::new();
    for fraction in [0.382, 0.618, 0.5, 0.25, 0.75] {
        for piece in region.outer.iter().chain(region.holes.iter().flatten()) {
            let at_t = piece.pspan.start + fraction * (piece.pspan.end - piece.pspan.start);
            let at = evaluate2(&piece.pcurve, at_t).map_err(|_| BooleanError::Evaluation)?;
            let mut tangent =
                derivative2(&piece.pcurve, at_t).map_err(|_| BooleanError::Evaluation)?;
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
                    out.push(
                        evaluate(surface, probe.x, probe.y)
                            .map_err(|_| BooleanError::Evaluation)?,
                    );
                    break;
                }
                step *= 0.5;
            }
            if out.len() >= 8 {
                return Ok(out);
            }
        }
    }
    if out.is_empty() {
        return Err(BooleanError::Undecided);
    }
    Ok(out)
}
