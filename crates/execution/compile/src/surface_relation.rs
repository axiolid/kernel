//! Swept and ruled surface relations as a B-rep face's support (#411).
//!
//! A [`SurfaceRelation`](axiolid_model::SurfaceRelation) is a graph node, not
//! an `axiolid_surface::Surface` value (ADR 0009), so the B-rep tessellator
//! reads it through a crate-private [`RelationSurface`]: a point at
//! `(u, v)`, the parameter periods, a family name, and a
//! [`RelationOracle`] returning the same
//! [`SurfaceBounds`](axiolid_reference::bound::SurfaceBounds) an atomic
//! surface does. Nothing here is public; ADR 0097 kept the evaluators of
//! its ruled surface crate-private for the same reason.
//!
//! # Parameterisation
//!
//! - **Revolution:** `S(u, v) = O + R_a(u) (c(v) - O)`, `c` the swept curve,
//!   `O` the axis origin, `a` the axis direction normalised (its length is
//!   not read). `R_a(u) w = (w.a) a + cos u (w - (w.a) a) + sin u (a x w)`
//!   turns right-handed about `a`; `u` is in radians, `u = 0` is the swept
//!   curve itself, and `v` is the curve's own parameter.
//!
//! - **LinearExtrusion:** `S(u, v) = c(u) + v d`, `d` the direction. `u` is
//!   the curve's own parameter and `v` is in units of the direction's
//!   length: it is not normalised.
//!
//! - **Ruled:** ADR 0097, `S(u, v) = (1 - v) A(u) + v B(u)` over `u, v` in
//!   `[0, 1]`, each side `A(u) = c(t0 + u (t1 - t0))` or an apex. A read
//!   outside the unit square beyond the slack is refused.
//!
//! The swept curve is a line, circle, ellipse or B-spline (the family
//! [`crate::ruled`] reads too); any other curve is refused by name. A
//! B-spline is read only inside its domain, within a slack of
//! `1e-12 max(1, span)` that is clamped; a conic is periodic and is read at
//! any parameter.
//!
//! # The bounds
//!
//! Rotations preserve length, so over a box `v0..v1` of the swept curve's
//! parameter, with `d` its certified derivative bounds
//! ([`axiolid_reference::bound::curve_derivative_bounds3`]): `|S_u| = |S_uu|
//! = dist(c(v), axis)`, `|S_v| = |c'|`, `|S_uv| = |c'_perp|` (the part of
//! `c'` normal to the axis) and `|S_vv| = |c''|`. The distance from the axis
//! is read at the box's midpoint and widened by the half width times the
//! perpendicular speed, which is `|dir - (dir.a) a|` for a line and
//! `d.first` otherwise. The bounds are smooth unless a tangent break of the
//! swept curve lies strictly inside the box. An extrusion over `u0..u1` has
//! `|S_u| = |c'|`, `|S_v| = |d|`, `|S_uu| = |c''|` and no mixed or `v`
//! second derivative.
//!
//! A ruled surface bounds each side over the box, `a1 = k sup|c'|` and `a2 =
//! k^2 sup|c''|` with `k = |t1 - t0|` (zero for an apex). With `um` the
//! box's middle in `u` and `h` its width, `|S_uv| <= |B'(um) - A'(um)| + h/2
//! (a2_A + a2_B)` and `|S_v| = |B - A| <= |B(um) - A(um)| + h/2 sup|S_uv|`;
//! `|S_u|` and `|S_uu|` are largest at an end of the `v` range, being a
//! blend of the sides' bounds, and `S_vv = 0`. The bounds are smooth unless
//! a side's tangent break, mapped to `u`, lies strictly inside.

use crate::ruled::{resolve, Side};
use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Interval, Point3, Scalar, Vec3};
use axiolid_curve::Curve3;
use axiolid_model::{GeometryGraph, GeometryNode, NodeId, RuledSide};
use axiolid_reference::bound::{
    continuity_breaks3, curve_derivative_bounds3, SurfaceBounds, ROUNDING,
};
use std::f64::consts::TAU;

/// A swept curve that is not an atomic 3D curve node.
pub(crate) const SWEPT_CURVE_NOT_ATOMIC: &str = "a swept curve that is not an atomic 3D curve";

/// A swept curve of a family that is not read.
pub(crate) const SWEPT_CURVE_FAMILY: &str =
    "a swept curve other than a line, circle, ellipse or B-spline";

/// A ruled surface read outside its unit parameter square.
pub(crate) const OUTSIDE_RULED_SQUARE: &str =
    "a ruled surface read outside its unit parameter square";

/// A swept B-spline read outside its domain.
pub(crate) const OUTSIDE_SWEPT_DOMAIN: &str =
    "a surface relation read outside its swept curve's domain";

/// How far past its domain a B-spline is read and clamped, relative to
/// `max(1, span)`.
const DOMAIN_SLACK: Scalar = 1e-12;

/// A swept curve and, for a B-spline, the domain it is read in.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SweptCurve<'g> {
    curve: &'g Curve3,
    domain: Option<Interval>,
}

impl<'g> SweptCurve<'g> {
    fn resolve(graph: &'g GeometryGraph, id: NodeId) -> GeomResult<Self> {
        let node = graph.get(id).ok_or_else(|| {
            GeomError::InvalidInput(format!(
                "surface relation swept curve {id:?} does not belong to this graph"
            ))
        })?;
        let GeometryNode::Curve3(curve) = node else {
            return Err(refused(SWEPT_CURVE_NOT_ATOMIC));
        };
        let domain = match curve {
            Curve3::Line(_) | Curve3::Circle(_) | Curve3::Ellipse(_) => None,
            Curve3::BSpline(_) => Some(axiolid_reference::curve::domain3(curve)),
            _ => return Err(refused(SWEPT_CURVE_FAMILY)),
        };
        Ok(Self { curve, domain })
    }

    /// The parameter to read the curve at: itself, or clamped into a
    /// B-spline's domain, refused beyond the slack.
    fn parameter(&self, t: Scalar) -> GeomResult<Scalar> {
        let Some(domain) = self.domain else {
            return Ok(t);
        };
        let slack = DOMAIN_SLACK * (domain.end - domain.start).abs().max(1.0);
        if !(t >= domain.start - slack && t <= domain.end + slack) {
            return Err(GeomError::InvalidInput(OUTSIDE_SWEPT_DOMAIN.into()));
        }
        Ok(t.clamp(domain.start, domain.end))
    }

    fn point(&self, t: Scalar) -> GeomResult<Point3> {
        axiolid_reference::curve::evaluate3(self.curve, self.parameter(t)?)
    }

    #[cfg(test)]
    fn derivatives(&self, t: Scalar) -> GeomResult<(Vec3, Vec3)> {
        let t = self.parameter(t)?;
        Ok((
            axiolid_reference::curve::derivative3(self.curve, t)?,
            axiolid_reference::curve::second_derivative3(self.curve, t)?,
        ))
    }

    /// The box clamped into the domain, or `None` when nothing is left.
    fn clamp(&self, (lo, hi): (Scalar, Scalar)) -> Option<(Scalar, Scalar)> {
        let (lo, hi) = (lo.min(hi), lo.max(hi));
        let (lo, hi) = match self.domain {
            Some(domain) => {
                let slack = DOMAIN_SLACK * (domain.end - domain.start).abs().max(1.0);
                if lo > domain.end + slack || hi < domain.start - slack {
                    return None;
                }
                (
                    lo.clamp(domain.start, domain.end),
                    hi.clamp(domain.start, domain.end),
                )
            }
            None => (lo, hi),
        };
        (lo.is_finite() && hi.is_finite() && lo <= hi).then_some((lo, hi))
    }

    fn is_conic(&self) -> bool {
        matches!(self.curve, Curve3::Circle(_) | Curve3::Ellipse(_))
    }

    fn breaks(&self) -> Vec<Scalar> {
        continuity_breaks3(self.curve, 1)
    }
}

/// A surface relation a B-rep face lies on.
#[derive(Debug, Clone, Copy)]
pub(crate) enum RelationSurface<'g> {
    /// `S(u, v) = O + R_a(u) (c(v) - O)`, `a` a unit axis.
    Revolution {
        curve: SweptCurve<'g>,
        origin: Point3,
        axis: Vec3,
    },
    /// `S(u, v) = c(u) + v d`.
    Extrusion {
        curve: SweptCurve<'g>,
        direction: Vec3,
    },
    /// `S(u, v) = (1 - v) A(u) + v B(u)` over the unit square.
    Ruled { sides: [Side<'g>; 2] },
}

/// A ruled-square parameter, clamped within the slack of `[0, 1]`.
fn unit(t: Scalar) -> GeomResult<Scalar> {
    if !(-DOMAIN_SLACK..=1.0 + DOMAIN_SLACK).contains(&t) {
        return Err(GeomError::InvalidInput(OUTSIDE_RULED_SQUARE.into()));
    }
    Ok(t.clamp(0.0, 1.0))
}

/// A parameter box clamped to `[0, 1]`, or `None` when nothing is left.
fn unit_box((lo, hi): (Scalar, Scalar)) -> Option<(Scalar, Scalar)> {
    let (lo, hi) = (lo.min(hi), lo.max(hi));
    if lo > 1.0 + DOMAIN_SLACK || hi < -DOMAIN_SLACK {
        return None;
    }
    let (lo, hi) = (lo.clamp(0.0, 1.0), hi.clamp(0.0, 1.0));
    (lo.is_finite() && hi.is_finite() && lo <= hi).then_some((lo, hi))
}

impl<'g> RelationSurface<'g> {
    /// A surface of revolution of `swept_curve` about the axis.
    pub(crate) fn revolution(
        graph: &'g GeometryGraph,
        swept_curve: NodeId,
        axis_origin: Point3,
        axis_direction: Vec3,
    ) -> GeomResult<Self> {
        if !axis_origin.is_finite() {
            return Err(GeomError::InvalidInput(
                "a revolution's axis origin is not finite".into(),
            ));
        }
        let axis = axis_direction.normalize_or_zero();
        if !axis_direction.is_finite() || axis == Vec3::ZERO {
            return Err(GeomError::InvalidInput(
                "a revolution's axis direction is zero or not finite".into(),
            ));
        }
        Ok(Self::Revolution {
            curve: SweptCurve::resolve(graph, swept_curve)?,
            origin: axis_origin,
            axis,
        })
    }

    /// A surface of linear extrusion of `swept_curve` along `direction`.
    pub(crate) fn extrusion(
        graph: &'g GeometryGraph,
        swept_curve: NodeId,
        direction: Vec3,
    ) -> GeomResult<Self> {
        if !direction.is_finite() || direction == Vec3::ZERO {
            return Err(GeomError::InvalidInput(
                "a linear extrusion's direction is zero or not finite".into(),
            ));
        }
        Ok(Self::Extrusion {
            curve: SweptCurve::resolve(graph, swept_curve)?,
            direction,
        })
    }

    /// The ruled surface between two sides.
    pub(crate) fn ruled(
        graph: &'g GeometryGraph,
        first: &RuledSide,
        second: &RuledSide,
    ) -> GeomResult<Self> {
        Ok(Self::Ruled {
            sides: [resolve(graph, first)?, resolve(graph, second)?],
        })
    }

    /// The point at `(u, v)`.
    pub(crate) fn evaluate(&self, u: Scalar, v: Scalar) -> GeomResult<Point3> {
        match self {
            Self::Revolution {
                curve,
                origin,
                axis,
            } => Ok(*origin + rotate(*axis, u, curve.point(v)? - *origin)),
            Self::Extrusion { curve, direction } => Ok(curve.point(u)? + *direction * v),
            Self::Ruled { sides } => {
                let (u, v) = (unit(u)?, unit(v)?);
                Ok(sides[0].at(u)? * (1.0 - v) + sides[1].at(u)? * v)
            }
        }
    }

    /// `S_u, S_v, S_uu, S_uv, S_vv`.
    #[cfg(test)]
    pub(crate) fn partials(&self, u: Scalar, v: Scalar) -> GeomResult<[Vec3; 5]> {
        match self {
            Self::Revolution {
                curve,
                origin,
                axis,
            } => {
                let a = *axis;
                let w = curve.point(v)? - *origin;
                let (d1, d2) = curve.derivatives(v)?;
                let spin = |x: Vec3| rotate(a, u, x);
                let (rw, r1) = (spin(w), spin(d1));
                Ok([a.cross(rw), r1, a.cross(a.cross(rw)), a.cross(r1), spin(d2)])
            }
            Self::Extrusion { curve, direction } => {
                let (d1, d2) = curve.derivatives(u)?;
                Ok([d1, *direction, d2, Vec3::ZERO, Vec3::ZERO])
            }
            Self::Ruled { sides } => {
                let (u, v) = (unit(u)?, unit(v)?);
                let [a, b] = sides;
                let (a1, a2) = (a.derivative(u)?, a.second_derivative(u)?);
                let (b1, b2) = (b.derivative(u)?, b.second_derivative(u)?);
                Ok([
                    a1 * (1.0 - v) + b1 * v,
                    b.at(u)? - a.at(u)?,
                    a2 * (1.0 - v) + b2 * v,
                    b1 - a1,
                    Vec3::ZERO,
                ])
            }
        }
    }

    /// The parameter periods.
    pub(crate) fn periods(&self) -> (Option<Scalar>, Option<Scalar>) {
        match self {
            Self::Revolution { curve, .. } => (Some(TAU), curve.is_conic().then_some(TAU)),
            Self::Extrusion { curve, .. } => (curve.is_conic().then_some(TAU), None),
            Self::Ruled { .. } => (None, None),
        }
    }

    /// The family's name, for a report.
    pub(crate) fn family(&self) -> &'static str {
        match self {
            Self::Revolution { .. } => "surface of revolution",
            Self::Extrusion { .. } => "surface of linear extrusion",
            Self::Ruled { .. } => "ruled surface",
        }
    }
}

/// `w` turned by `u` right-handed about the unit axis `a`.
fn rotate(a: Vec3, u: Scalar, w: Vec3) -> Vec3 {
    let along = a * w.dot(a);
    along + (w - along) * u.cos() + a.cross(w) * u.sin()
}

fn refused(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: crate::BACKEND_ID,
        operation: Operation::SurfaceEvaluation,
        input,
    }
}

/// The derivative bounds of a [`RelationSurface`] over parameter boxes.
#[derive(Debug, Clone)]
pub(crate) struct RelationOracle<'g> {
    surface: RelationSurface<'g>,
    /// The tangent breaks in `u` and in `v`.
    breaks: [Vec<Scalar>; 2],
}

impl<'g> RelationOracle<'g> {
    pub(crate) fn new(surface: RelationSurface<'g>) -> Self {
        let breaks = match &surface {
            RelationSurface::Revolution { curve, .. } => [Vec::new(), curve.breaks()],
            RelationSurface::Extrusion { curve, .. } => [curve.breaks(), Vec::new()],
            RelationSurface::Ruled { sides } => {
                let mut breaks: Vec<Scalar> = sides.iter().flat_map(Side::breaks_in_u).collect();
                breaks.sort_by(Scalar::total_cmp);
                [breaks, Vec::new()]
            }
        };
        Self { surface, breaks }
    }

    /// The bounds over the box, clamped to where the surface is read, or
    /// `None` when none can be had.
    pub(crate) fn bounds(&self, u: (Scalar, Scalar), v: (Scalar, Scalar)) -> Option<SurfaceBounds> {
        let smooth = |axis: usize, (lo, hi): (Scalar, Scalar)| {
            !self.breaks[axis].iter().any(|&b| b > lo && b < hi)
        };
        let bounds = match &self.surface {
            RelationSurface::Revolution {
                curve,
                origin,
                axis,
            } => {
                let (v0, v1) = curve.clamp(v)?;
                let d = curve_derivative_bounds3(curve.curve, v0, v1)?;
                let perp = match curve.curve {
                    Curve3::Line(line) => {
                        (line.direction - *axis * line.direction.dot(*axis)).length()
                    }
                    _ => d.first,
                };
                let w = curve.point(0.5 * (v0 + v1)).ok()? - *origin;
                let reach = (w - *axis * w.dot(*axis)).length()
                    + 0.5 * (v1 - v0) * perp
                    + ROUNDING * (w.length() + origin.length());
                SurfaceBounds {
                    du: reach,
                    dv: d.first,
                    duu: reach,
                    duv: perp,
                    dvv: d.second,
                    smooth: smooth(1, (v0, v1)),
                    interpolation: [reach, perp, d.second],
                }
            }
            RelationSurface::Ruled { sides } => {
                let (u0, u1) = unit_box(u)?;
                let (v0, v1) = unit_box(v)?;
                let (h, um) = (u1 - u0, 0.5 * (u0 + u1));
                let [(a1, a2), (b1, b2)] = [sides[0].speeds(u0, u1)?, sides[1].speeds(u0, u1)?];
                let (a, b) = (sides[0].at(um).ok()?, sides[1].at(um).ok()?);
                let (da, db) = (sides[0].derivative(um).ok()?, sides[1].derivative(um).ok()?);
                // Across a tangent break `A'` or `B'` jumps, so the Lipschitz
                // form no longer holds; `|B' - A'| <= a1 + b1` does.
                let smooth_box = smooth(0, (u0, u1));
                let lipschitz = (db - da).length()
                    + 0.5 * h * (a2 + b2)
                    + ROUNDING * (da.length() + db.length());
                let duv = if smooth_box {
                    lipschitz.min(a1 + b1)
                } else {
                    a1 + b1
                };
                let dv = (b - a).length() + 0.5 * h * duv + ROUNDING * (a.length() + b.length());
                let blend = |a: Scalar, b: Scalar| {
                    [v0, v1]
                        .into_iter()
                        .map(|v| (1.0 - v) * a + v * b)
                        .fold(0.0, Scalar::max)
                };
                let (du, duu) = (blend(a1, b1), blend(a2, b2));
                SurfaceBounds {
                    du,
                    dv,
                    duu,
                    duv,
                    dvv: 0.0,
                    smooth: smooth_box,
                    interpolation: [duu, duv, 0.0],
                }
            }
            RelationSurface::Extrusion { curve, direction } => {
                let (u0, u1) = curve.clamp(u)?;
                let d = curve_derivative_bounds3(curve.curve, u0, u1)?;
                SurfaceBounds {
                    du: d.first,
                    dv: direction.length(),
                    duu: d.second,
                    duv: 0.0,
                    dvv: 0.0,
                    smooth: smooth(0, (u0, u1)),
                    interpolation: [d.second, 0.0, 0.0],
                }
            }
        };
        finite(inflate(bounds))
    }
}

fn inflate(bounds: SurfaceBounds) -> SurfaceBounds {
    let up = |x: Scalar| x * (1.0 + ROUNDING);
    SurfaceBounds {
        du: up(bounds.du),
        dv: up(bounds.dv),
        duu: up(bounds.duu),
        duv: up(bounds.duv),
        dvv: up(bounds.dvv),
        smooth: bounds.smooth,
        interpolation: bounds.interpolation.map(up),
    }
}

fn finite(bounds: SurfaceBounds) -> Option<SurfaceBounds> {
    [bounds.du, bounds.dv, bounds.duu, bounds.duv, bounds.dvv]
        .into_iter()
        .chain(bounds.interpolation)
        .all(Scalar::is_finite)
        .then_some(bounds)
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::{Frame3, Point2, Vec2};
    use axiolid_curve::{
        BSplineCurve, Circle3, Curve2, Ellipse3, KnotSpec, Line2, Line3, Polyline3,
    };
    use axiolid_model::GeometryGraphBuilder;
    use axiolid_surface::{Cone, Cylinder, Sphere, Surface, Torus};
    use std::f64::consts::PI;

    fn frame_at(z: Scalar) -> Frame3 {
        Frame3 {
            origin: Point3::new(0.0, 0.0, z),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        }
    }

    fn graph_of(node: GeometryNode) -> (GeometryGraph, NodeId) {
        let mut builder = GeometryGraphBuilder::new();
        let id = builder.push(node).unwrap();
        (builder.finish(vec![id]).unwrap(), id)
    }

    fn curve_graph(curve: Curve3) -> (GeometryGraph, NodeId) {
        graph_of(GeometryNode::Curve3(curve))
    }

    fn line(origin: Point3, direction: Vec3) -> Curve3 {
        Curve3::Line(Line3 { origin, direction })
    }

    fn spline(degree: u16, points: &[Point3], weights: Option<Vec<Scalar>>) -> Curve3 {
        Curve3::BSpline(BSplineCurve {
            degree,
            control_points: points.to_vec(),
            knots: vec![0.0, 1.0],
            multiplicities: vec![u32::from(degree) + 1; 2],
            weights,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        })
    }

    fn about_z(graph: &GeometryGraph, id: NodeId) -> RelationSurface<'_> {
        RelationSurface::revolution(graph, id, Point3::ZERO, Vec3::Z).unwrap()
    }

    /// Deterministic values in `[0, 1)`.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> Scalar {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as Scalar / (1u64 << 53) as Scalar
        }
    }

    fn grid() -> impl Iterator<Item = (Scalar, Scalar)> {
        (0..7).flat_map(|i| (0..5).map(move |j| (-1.3 + 0.9 * i as Scalar, 0.1 * j as Scalar)))
    }

    fn assert_same(relation: &RelationSurface<'_>, surface: &Surface) {
        for (u, v) in grid() {
            let exact = axiolid_reference::surface::evaluate(surface, u, v).unwrap();
            let got = relation.evaluate(u, v).unwrap();
            assert!(
                (got - exact).length() <= 1e-12,
                "({u}, {v}): {got:?} vs {exact:?}"
            );
        }
    }

    #[test]
    fn a_line_revolved_about_a_parallel_axis_evaluates_to_its_cylinder() {
        let (graph, id) = curve_graph(line(Point3::new(0.3, 0.0, 0.0), Vec3::Z));
        let cylinder = Surface::Cylinder(Cylinder {
            frame: frame_at(0.0),
            radius: 0.3,
        });
        assert_same(&about_z(&graph, id), &cylinder);
    }

    #[test]
    fn a_slanted_line_revolved_evaluates_to_its_cone() {
        let (graph, id) = curve_graph(line(Point3::new(0.2, 0.0, 0.0), Vec3::new(0.2, 0.0, 1.0)));
        let cone = Surface::Cone(Cone {
            frame: frame_at(0.0),
            radius: 0.2,
            semi_angle: 0.2_f64.atan(),
        });
        assert_same(&about_z(&graph, id), &cone);
    }

    #[test]
    fn an_arc_revolved_about_its_diameter_evaluates_to_its_sphere() {
        let arc = Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Z,
                z: -Vec3::Y,
            },
            radius: 0.3,
        });
        let (graph, id) = curve_graph(arc);
        let sphere = Surface::Sphere(Sphere {
            frame: frame_at(0.0),
            radius: 0.3,
        });
        assert_same(&about_z(&graph, id), &sphere);
    }

    #[test]
    fn a_circle_revolved_off_its_axis_evaluates_to_its_torus() {
        let tube = Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: Point3::new(0.4, 0.0, 0.0),
                x: Vec3::X,
                y: Vec3::Z,
                z: -Vec3::Y,
            },
            radius: 0.1,
        });
        let (graph, id) = curve_graph(tube);
        let torus = Surface::Torus(Torus {
            frame: frame_at(0.0),
            major_radius: 0.4,
            minor_radius: 0.1,
        });
        assert_same(&about_z(&graph, id), &torus);
    }

    #[test]
    fn the_rotation_is_right_handed_about_the_axis() {
        let (graph, id) = curve_graph(line(Point3::new(0.3, 0.0, 0.0), Vec3::Z));
        let up = about_z(&graph, id);
        let quarter = core::f64::consts::FRAC_PI_2;
        assert!((up.evaluate(quarter, 0.0).unwrap() - Point3::new(0.0, 0.3, 0.0)).length() < 1e-12);
        let down = RelationSurface::revolution(&graph, id, Point3::ZERO, -Vec3::Z).unwrap();
        assert!(
            (down.evaluate(quarter, 0.0).unwrap() - Point3::new(0.0, -0.3, 0.0)).length() < 1e-12
        );
    }

    #[test]
    fn a_non_unit_axis_is_normalised() {
        let (graph, id) = curve_graph(line(Point3::new(0.3, 0.1, 0.0), Vec3::new(0.1, 0.0, 1.0)));
        let unit = about_z(&graph, id);
        let long = RelationSurface::revolution(&graph, id, Point3::ZERO, Vec3::new(0.0, 0.0, 2.0))
            .unwrap();
        for (u, v) in grid() {
            let (a, b) = (unit.evaluate(u, v).unwrap(), long.evaluate(u, v).unwrap());
            assert!((a - b).length() <= 1e-15);
        }
    }

    fn bounds_of(
        surface: RelationSurface<'_>,
        u: (Scalar, Scalar),
        v: (Scalar, Scalar),
    ) -> SurfaceBounds {
        RelationOracle::new(surface).bounds(u, v).unwrap()
    }

    #[test]
    fn the_bounds_of_a_cylinder_revolution_are_its_own() {
        let (graph, id) = curve_graph(line(Point3::new(0.3, 0.0, 0.0), Vec3::Z));
        let b = bounds_of(about_z(&graph, id), (0.0, 1.0), (0.0, 0.5));
        assert!((b.du - 0.3).abs() <= 1e-9 && (b.duu - 0.3).abs() <= 1e-9);
        assert!((b.dv - 1.0).abs() <= 1e-9);
        assert_eq!(b.duv, 0.0);
        assert_eq!(b.dvv, 0.0);
        assert!(b.smooth);
    }

    #[test]
    fn the_bounds_of_a_cone_revolution() {
        let (graph, id) = curve_graph(line(Point3::new(0.2, 0.0, 0.0), Vec3::new(0.2, 0.0, 1.0)));
        let b = bounds_of(about_z(&graph, id), (0.0, 1.0), (0.0, 0.5));
        // The radius is largest at the top of the box, 0.2 + 0.2 * 0.5.
        assert!((b.du - 0.3).abs() <= 1e-9, "{}", b.du);
        assert!((b.duv - 0.2).abs() <= 1e-9, "{}", b.duv);
        assert!((b.dv - 1.019803902718557).abs() <= 1e-9, "{}", b.dv);
        assert_eq!(b.dvv, 0.0);
    }

    /// Profiles of an ellipse, a rational quadratic B-spline and a skew line.
    fn profiles() -> Vec<(Curve3, (Scalar, Scalar))> {
        let ellipse = Curve3::Ellipse(Ellipse3 {
            frame: Frame3 {
                origin: Point3::new(0.5, 0.1, 0.2),
                x: Vec3::X,
                y: Vec3::new(0.0, 0.6, 0.8),
                z: Vec3::new(0.0, -0.8, 0.6),
            },
            semi_axis_x: 0.3,
            semi_axis_y: 0.15,
        });
        let rational = spline(
            2,
            &[
                Point3::new(0.3, 0.0, 0.0),
                Point3::new(0.6, 0.1, 0.4),
                Point3::new(0.2, -0.1, 0.9),
            ],
            Some(vec![1.0, 0.6, 1.0]),
        );
        let skew = line(Point3::new(0.3, 0.1, 0.2), Vec3::new(0.5, -0.2, 1.0));
        vec![
            (ellipse, (0.0, TAU)),
            (rational, (0.0, 1.0)),
            (skew, (-1.0, 1.0)),
        ]
    }

    fn skewed_axis() -> (Point3, Vec3) {
        (Point3::new(0.1, -0.2, 0.3), Vec3::new(0.2, 0.1, 1.0))
    }

    #[test]
    fn bounds_dominate_sampled_derivatives() {
        let mut random = Lcg(7);
        for (curve, (lo, hi)) in profiles() {
            let (graph, id) = curve_graph(curve);
            let (origin, direction) = skewed_axis();
            let surface = RelationSurface::revolution(&graph, id, origin, direction).unwrap();
            let oracle = RelationOracle::new(surface);
            for _ in 0..50 {
                let u0 = -3.0 + 6.0 * random.next();
                let u1 = u0 + 3.0 * random.next();
                let v0 = lo + (hi - lo) * random.next();
                let v1 = v0 + (hi - v0) * random.next();
                let b = oracle.bounds((u0, u1), (v0, v1)).unwrap();
                for _ in 0..30 {
                    let u = u0 + (u1 - u0) * random.next();
                    let v = v0 + (v1 - v0) * random.next();
                    let [su, sv, suu, suv, svv] =
                        surface.partials(u, v).unwrap().map(|x| x.length());
                    let limit = |x: Scalar, bound: Scalar| {
                        assert!(x <= bound * (1.0 + 1e-9), "{x} > {bound} at ({u}, {v})");
                    };
                    limit(su, b.du);
                    limit(sv, b.dv);
                    limit(suu, b.duu);
                    limit(suv, b.duv);
                    limit(svv, b.dvv);
                    limit(suu, b.interpolation[0]);
                    limit(suv, b.interpolation[1]);
                    limit(svv, b.interpolation[2]);
                }
            }
        }
    }

    #[test]
    fn the_partials_match_central_differences() {
        let h = 1e-5;
        for (curve, (lo, hi)) in profiles() {
            let (graph, id) = curve_graph(curve);
            let (origin, direction) = skewed_axis();
            let surface = RelationSurface::revolution(&graph, id, origin, direction).unwrap();
            let (u, v) = (0.7, lo + 0.4 * (hi - lo));
            let at = |u: Scalar, v: Scalar| surface.evaluate(u, v).unwrap();
            let [su, sv, suu, suv, svv] = surface.partials(u, v).unwrap();
            let first = |a: Vec3, b: Vec3| (a - b) / (2.0 * h);
            let second = |a: Vec3, b: Vec3, c: Vec3| (a - 2.0 * b + c) / (h * h);
            let close = |got: Vec3, want: Vec3, scale: Scalar| {
                assert!(
                    (got - want).length() <= scale * (1.0 + want.length()),
                    "{got:?} vs {want:?}"
                );
            };
            close(first(at(u + h, v), at(u - h, v)), su, 1e-6);
            close(first(at(u, v + h), at(u, v - h)), sv, 1e-6);
            close(second(at(u + h, v), at(u, v), at(u - h, v)), suu, 1e-3);
            close(second(at(u, v + h), at(u, v), at(u, v - h)), svv, 1e-3);
            let cross = (at(u + h, v + h) - at(u + h, v - h) - at(u - h, v + h) + at(u - h, v - h))
                / (4.0 * h * h);
            close(cross, suv, 1e-3);
        }
    }

    #[test]
    fn periods_are_a_turn_in_u_and_the_profile_s_own_in_v() {
        let (graph, id) = curve_graph(line(Point3::new(0.3, 0.0, 0.0), Vec3::Z));
        assert_eq!(about_z(&graph, id).periods(), (Some(TAU), None));
        let (graph, id) = curve_graph(profiles().remove(0).0);
        assert_eq!(about_z(&graph, id).periods(), (Some(TAU), Some(TAU)));
        let (graph, id) = curve_graph(profiles().remove(1).0);
        assert_eq!(about_z(&graph, id).periods(), (Some(TAU), None));
    }

    fn kinked() -> Curve3 {
        Curve3::BSpline(BSplineCurve {
            degree: 1,
            control_points: vec![
                Point3::new(0.3, 0.0, 0.0),
                Point3::new(0.3, 0.0, 0.5),
                Point3::new(0.4, 0.0, 1.0),
            ],
            knots: vec![0.0, 0.5, 1.0],
            multiplicities: vec![2, 1, 2],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        })
    }

    #[test]
    fn a_kink_inside_the_box_is_not_smooth() {
        let (graph, id) = curve_graph(kinked());
        let surface = about_z(&graph, id);
        assert!(!bounds_of(surface, (0.0, 1.0), (0.4, 0.6)).smooth);
        assert!(bounds_of(surface, (0.0, 1.0), (0.0, 0.4)).smooth);
    }

    #[test]
    fn a_spline_read_past_its_domain_is_refused() {
        let (graph, id) = curve_graph(kinked());
        let surface = about_z(&graph, id);
        let error = surface.evaluate(0.0, 1.0 + 1e-6).unwrap_err();
        assert!(matches!(&error, GeomError::InvalidInput(text) if text == OUTSIDE_SWEPT_DOMAIN));
        assert!(surface.evaluate(0.0, 1.0 + 1e-14).is_ok());
    }

    #[test]
    fn a_box_reaching_past_a_spline_domain_is_bounded_over_its_part_inside() {
        let (graph, id) = curve_graph(kinked());
        let oracle = RelationOracle::new(about_z(&graph, id));
        assert!(oracle.bounds((0.0, 1.0), (0.5, 2.0)).is_some());
        assert!(oracle.bounds((0.0, 1.0), (1.5, 2.0)).is_none());
    }

    fn unsupported(error: GeomError) -> &'static str {
        match error {
            GeomError::UnsupportedInput { input, .. } => input,
            other => panic!("not a named refusal: {other:?}"),
        }
    }

    #[test]
    fn swept_curves_that_are_not_read_are_refused_by_name() {
        let polyline = Curve3::Polyline(Polyline3 {
            points: vec![Point3::ZERO, Point3::X],
            closed: false,
        });
        let (graph, id) = curve_graph(polyline);
        let error = RelationSurface::revolution(&graph, id, Point3::ZERO, Vec3::Z).unwrap_err();
        assert_eq!(unsupported(error), SWEPT_CURVE_FAMILY);

        let planar = Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        });
        let (graph, id) = graph_of(GeometryNode::Curve2(planar));
        let error = RelationSurface::revolution(&graph, id, Point3::ZERO, Vec3::Z).unwrap_err();
        assert_eq!(unsupported(error), SWEPT_CURVE_NOT_ATOMIC);
    }

    #[test]
    fn a_bad_axis_is_invalid_input() {
        let (graph, id) = curve_graph(line(Point3::new(0.3, 0.0, 0.0), Vec3::Z));
        for (origin, direction) in [
            (Point3::ZERO, Vec3::ZERO),
            (Point3::ZERO, Vec3::new(0.0, 0.0, Scalar::NAN)),
            (Point3::new(Scalar::INFINITY, 0.0, 0.0), Vec3::Z),
        ] {
            let error = RelationSurface::revolution(&graph, id, origin, direction).unwrap_err();
            assert!(matches!(error, GeomError::InvalidInput(_)), "{error:?}");
        }
    }

    fn circle_profile() -> Curve3 {
        Curve3::Circle(Circle3 {
            frame: frame_at(0.0),
            radius: 0.3,
        })
    }

    #[test]
    fn a_circle_extruded_along_its_normal_evaluates_to_its_cylinder() {
        let (graph, id) = curve_graph(circle_profile());
        let surface = RelationSurface::extrusion(&graph, id, Vec3::Z).unwrap();
        let cylinder = Surface::Cylinder(Cylinder {
            frame: frame_at(0.0),
            radius: 0.3,
        });
        assert_same(&surface, &cylinder);
    }

    #[test]
    fn an_extrusion_reads_v_in_units_of_its_direction() {
        let (graph, id) = curve_graph(circle_profile());
        let surface = RelationSurface::extrusion(&graph, id, Vec3::new(0.0, 0.0, 2.0)).unwrap();
        let moved = surface.evaluate(0.0, 1.0).unwrap() - surface.evaluate(0.0, 0.0).unwrap();
        assert_eq!(moved, Vec3::new(0.0, 0.0, 2.0));
    }

    #[test]
    fn the_bounds_of_a_cylinder_extrusion_are_its_own() {
        let (graph, id) = curve_graph(circle_profile());
        let surface = RelationSurface::extrusion(&graph, id, Vec3::Z).unwrap();
        let b = bounds_of(surface, (0.0, 1.0), (0.0, 0.5));
        assert!((b.du - 0.3).abs() <= 1e-9 && (b.duu - 0.3).abs() <= 1e-9);
        assert!((b.dv - 1.0).abs() <= 1e-9);
        assert_eq!((b.duv, b.dvv), (0.0, 0.0));
        assert!(b.smooth);
    }

    fn oblique() -> Vec3 {
        Vec3::new(0.2, -0.1, 1.0)
    }

    #[test]
    fn extrusion_bounds_dominate_sampled_derivatives() {
        let mut random = Lcg(11);
        for (curve, (lo, hi)) in profiles() {
            let (graph, id) = curve_graph(curve);
            let surface = RelationSurface::extrusion(&graph, id, oblique()).unwrap();
            let oracle = RelationOracle::new(surface);
            for _ in 0..50 {
                let u0 = lo + (hi - lo) * random.next();
                let u1 = u0 + (hi - u0) * random.next();
                let v0 = -2.0 + 4.0 * random.next();
                let v1 = v0 + 2.0 * random.next();
                let b = oracle.bounds((u0, u1), (v0, v1)).unwrap();
                for _ in 0..30 {
                    let u = u0 + (u1 - u0) * random.next();
                    let v = v0 + (v1 - v0) * random.next();
                    let [su, sv, suu, suv, svv] =
                        surface.partials(u, v).unwrap().map(|x| x.length());
                    let limit = |x: Scalar, bound: Scalar| {
                        assert!(x <= bound * (1.0 + 1e-9), "{x} > {bound} at ({u}, {v})");
                    };
                    limit(su, b.du);
                    limit(sv, b.dv);
                    limit(suu, b.duu);
                    limit(suv, b.duv);
                    limit(svv, b.dvv);
                    limit(suu, b.interpolation[0]);
                }
            }
        }
    }

    #[test]
    fn the_extrusion_partials_match_central_differences() {
        let h = 1e-5;
        for (curve, (lo, hi)) in profiles() {
            let (graph, id) = curve_graph(curve);
            let surface = RelationSurface::extrusion(&graph, id, oblique()).unwrap();
            let (u, v) = (lo + 0.4 * (hi - lo), 0.7);
            let at = |u: Scalar, v: Scalar| surface.evaluate(u, v).unwrap();
            let [su, sv, suu, _, _] = surface.partials(u, v).unwrap();
            let close = |got: Vec3, want: Vec3, scale: Scalar| {
                assert!(
                    (got - want).length() <= scale * (1.0 + want.length()),
                    "{got:?} vs {want:?}"
                );
            };
            close((at(u + h, v) - at(u - h, v)) / (2.0 * h), su, 1e-6);
            close((at(u, v + h) - at(u, v - h)) / (2.0 * h), sv, 1e-6);
            close(
                (at(u + h, v) - 2.0 * at(u, v) + at(u - h, v)) / (h * h),
                suu,
                1e-3,
            );
        }
    }

    #[test]
    fn extrusion_periods_are_a_turn_in_u_for_a_conic() {
        let (graph, id) = curve_graph(circle_profile());
        let surface = RelationSurface::extrusion(&graph, id, Vec3::Z).unwrap();
        assert_eq!(surface.periods(), (Some(TAU), None));
        let (graph, id) = curve_graph(line(Point3::ZERO, Vec3::X));
        let surface = RelationSurface::extrusion(&graph, id, Vec3::Z).unwrap();
        assert_eq!(surface.periods(), (None, None));
        let (graph, id) = curve_graph(profiles().remove(1).0);
        let surface = RelationSurface::extrusion(&graph, id, Vec3::Z).unwrap();
        assert_eq!(surface.periods(), (None, None));
    }

    #[test]
    fn an_extrusion_over_a_spline_is_read_inside_its_domain_only() {
        let (graph, id) = curve_graph(kinked());
        let surface = RelationSurface::extrusion(&graph, id, Vec3::Z).unwrap();
        let error = surface.evaluate(1.0 + 1e-6, 0.0).unwrap_err();
        assert!(matches!(&error, GeomError::InvalidInput(text) if text == OUTSIDE_SWEPT_DOMAIN));
        assert!(surface.evaluate(1.0 + 1e-14, 0.0).is_ok());
        assert!(!bounds_of(surface, (0.4, 0.6), (0.0, 1.0)).smooth);
        assert!(bounds_of(surface, (0.0, 0.4), (0.0, 1.0)).smooth);
        let oracle = RelationOracle::new(surface);
        assert!(oracle.bounds((0.5, 2.0), (0.0, 1.0)).is_some());
        assert!(oracle.bounds((1.5, 2.0), (0.0, 1.0)).is_none());
    }

    #[test]
    fn a_bad_extrusion_is_refused() {
        let (graph, id) = curve_graph(circle_profile());
        for direction in [Vec3::ZERO, Vec3::new(0.0, Scalar::NAN, 1.0)] {
            let error = RelationSurface::extrusion(&graph, id, direction).unwrap_err();
            assert!(matches!(error, GeomError::InvalidInput(_)), "{error:?}");
        }
        let polyline = Curve3::Polyline(Polyline3 {
            points: vec![Point3::ZERO, Point3::X],
            closed: false,
        });
        let (graph, id) = curve_graph(polyline);
        let error = RelationSurface::extrusion(&graph, id, Vec3::Z).unwrap_err();
        assert_eq!(unsupported(error), SWEPT_CURVE_FAMILY);
    }

    enum Spec {
        Curve(Box<Curve3>, Scalar, Scalar),
        Apex(Point3),
    }

    impl Spec {
        fn curve(curve: Curve3, start: Scalar, end: Scalar) -> Self {
            Self::Curve(Box::new(curve), start, end)
        }
    }

    /// A graph holding the curves of two sides, and the sides.
    fn ruled_graph(specs: [Spec; 2]) -> (GeometryGraph, [RuledSide; 2]) {
        let mut builder = GeometryGraphBuilder::new();
        let mut first = None;
        let sides = specs.map(|spec| match spec {
            Spec::Curve(curve, start, end) => {
                let id = builder.push(GeometryNode::Curve3(*curve)).unwrap();
                first.get_or_insert(id);
                RuledSide::Curve {
                    curve: id,
                    interval: Interval { start, end },
                }
            }
            Spec::Apex(point) => RuledSide::Apex(point),
        });
        (builder.finish(vec![first.unwrap()]).unwrap(), sides)
    }

    fn unit_circle(z: Scalar, radius: Scalar) -> Curve3 {
        Curve3::Circle(Circle3 {
            frame: frame_at(z),
            radius,
        })
    }

    fn arc_and_line(reversed: bool) -> [Spec; 2] {
        let (a, b) = if reversed {
            ((core::f64::consts::FRAC_PI_2, 0.0), (1.0, 0.0))
        } else {
            ((0.0, core::f64::consts::FRAC_PI_2), (0.0, 1.0))
        };
        [
            Spec::curve(unit_circle(0.0, 1.0), a.0, a.1),
            Spec::curve(
                line(Point3::new(1.0, 0.0, 1.0), Vec3::new(-1.0, 1.0, 0.0)),
                b.0,
                b.1,
            ),
        ]
    }

    /// A degree-1 side with a sharp kink at its middle knot.
    fn sharp_kink() -> Curve3 {
        Curve3::BSpline(BSplineCurve {
            degree: 1,
            control_points: vec![
                Point3::ZERO,
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 30.0, 0.0),
            ],
            knots: vec![0.0, 0.5, 1.0],
            multiplicities: vec![2, 1, 2],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        })
    }

    fn ruled_cases() -> Vec<[Spec; 2]> {
        let rational = profiles().remove(1).0;
        vec![
            [
                Spec::curve(sharp_kink(), 0.0, 1.0),
                Spec::curve(line(Point3::new(0.0, 0.0, 1.0), Vec3::X), 0.0, 1.0),
            ],
            [
                Spec::curve(sharp_kink(), 0.9, 0.3),
                Spec::curve(line(Point3::new(0.0, 0.0, 1.0), Vec3::X), 1.0, 0.0),
            ],
            arc_and_line(false),
            arc_and_line(true),
            [
                Spec::curve(unit_circle(0.0, 1.0), 0.0, PI),
                Spec::Apex(Point3::new(0.0, 0.0, 1.0)),
            ],
            [
                Spec::Apex(Point3::new(0.0, 0.0, 1.0)),
                Spec::curve(unit_circle(0.0, 1.0), 0.0, PI),
            ],
            [
                Spec::curve(unit_circle(0.0, 1.0), 0.0, core::f64::consts::FRAC_PI_2),
                Spec::curve(unit_circle(1.0, 0.5), 0.0, core::f64::consts::FRAC_PI_3),
            ],
            [
                Spec::curve(rational, 0.9, 0.1),
                Spec::curve(profiles().remove(0).0, 0.2, 4.0),
            ],
        ]
    }

    #[test]
    fn ruled_evaluates_as_the_sheet_does() {
        let (graph, [first, second]) = ruled_graph(arc_and_line(false));
        let surface = RelationSurface::ruled(&graph, &first, &second).unwrap();
        for (u, v) in grid().map(|(u, v)| ((u + 1.3) / 5.4, v * 2.0)) {
            let angle = core::f64::consts::FRAC_PI_2 * u;
            let a = Point3::new(angle.cos(), angle.sin(), 0.0);
            let b = Point3::new(1.0 - u, u, 1.0);
            let want = a * (1.0 - v) + b * v;
            assert!((surface.evaluate(u, v).unwrap() - want).length() <= 1e-12);
        }
    }

    #[test]
    fn the_bounds_of_two_parallel_circles_have_no_twist() {
        let (graph, [first, second]) = ruled_graph([
            Spec::curve(unit_circle(0.0, 0.3), 0.0, TAU),
            Spec::curve(unit_circle(0.5, 0.3), 0.0, TAU),
        ]);
        let surface = RelationSurface::ruled(&graph, &first, &second).unwrap();
        let h = 0.1;
        let b = bounds_of(surface, (0.2, 0.2 + h), (0.0, 1.0));
        assert!(
            b.duv <= 0.5 * h * 2.0 * TAU * TAU * 0.3 * (1.0 + 1e-9),
            "{}",
            b.duv
        );
        assert_eq!(b.dvv, 0.0);
        assert!(b.smooth);
    }

    #[test]
    fn ruled_bounds_dominate_sampled_derivatives() {
        let mut random = Lcg(23);
        for specs in ruled_cases() {
            let (graph, [first, second]) = ruled_graph(specs);
            let surface = RelationSurface::ruled(&graph, &first, &second).unwrap();
            let oracle = RelationOracle::new(surface);
            for _ in 0..50 {
                let u0 = random.next();
                let u1 = u0 + (1.0 - u0) * random.next();
                let v0 = random.next();
                let v1 = v0 + (1.0 - v0) * random.next();
                let b = oracle.bounds((u0, u1), (v0, v1)).unwrap();
                for _ in 0..30 {
                    let u = u0 + (u1 - u0) * random.next();
                    let v = v0 + (v1 - v0) * random.next();
                    let [su, sv, suu, suv, svv] =
                        surface.partials(u, v).unwrap().map(|x| x.length());
                    let limit = |x: Scalar, bound: Scalar| {
                        assert!(x <= bound * (1.0 + 1e-9), "{x} > {bound} at ({u}, {v})");
                    };
                    limit(su, b.du);
                    limit(sv, b.dv);
                    limit(suu, b.duu);
                    limit(suv, b.duv);
                    limit(svv, b.dvv);
                    limit(suu, b.interpolation[0]);
                    limit(suv, b.interpolation[1]);
                }
            }
        }
    }

    #[test]
    fn the_ruled_partials_match_central_differences() {
        let h = 1e-5;
        for specs in ruled_cases() {
            let (graph, [first, second]) = ruled_graph(specs);
            let surface = RelationSurface::ruled(&graph, &first, &second).unwrap();
            let (u, v) = (0.4, 0.7);
            let at = |u: Scalar, v: Scalar| surface.evaluate(u, v).unwrap();
            let [su, sv, suu, suv, _] = surface.partials(u, v).unwrap();
            let close = |got: Vec3, want: Vec3, scale: Scalar| {
                assert!(
                    (got - want).length() <= scale * (1.0 + want.length()),
                    "{got:?} vs {want:?}"
                );
            };
            close((at(u + h, v) - at(u - h, v)) / (2.0 * h), su, 1e-6);
            close((at(u, v + h) - at(u, v - h)) / (2.0 * h), sv, 1e-6);
            close(
                (at(u + h, v) - 2.0 * at(u, v) + at(u - h, v)) / (h * h),
                suu,
                1e-3,
            );
            let cross = (at(u + h, v + h) - at(u + h, v - h) - at(u - h, v + h) + at(u - h, v - h))
                / (4.0 * h * h);
            close(cross, suv, 1e-3);
        }
    }

    #[test]
    fn ruled_breaks_make_the_box_not_smooth() {
        let (graph, [first, second]) = ruled_graph([
            Spec::curve(kinked(), 0.0, 1.0),
            Spec::curve(line(Point3::new(0.0, 1.0, 0.0), Vec3::Z), 0.0, 1.0),
        ]);
        let surface = RelationSurface::ruled(&graph, &first, &second).unwrap();
        assert!(!bounds_of(surface, (0.4, 0.6), (0.0, 1.0)).smooth);
        assert!(bounds_of(surface, (0.0, 0.4), (0.0, 1.0)).smooth);
    }

    #[test]
    fn a_kink_across_a_box_bounds_the_jump_of_the_sides_speed() {
        let (graph, [first, second]) = ruled_graph([
            Spec::curve(sharp_kink(), 0.0, 1.0),
            Spec::curve(line(Point3::new(0.0, 0.0, 1.0), Vec3::X), 0.0, 1.0),
        ]);
        let surface = RelationSurface::ruled(&graph, &first, &second).unwrap();
        let b = bounds_of(surface, (0.2, 0.6), (0.0, 1.0));
        assert!(!b.smooth);
        let [_, sv, _, suv, _] = surface.partials(0.6, 0.5).unwrap().map(|x| x.length());
        assert!(
            sv <= b.dv && suv <= b.duv,
            "{sv} {} / {suv} {}",
            b.dv,
            b.duv
        );
    }

    #[test]
    fn a_break_is_mapped_through_a_shifted_reversed_interval() {
        // The kink at t = 0.5 is at u = (0.5 - 0.9) / (0.3 - 0.9) = 2/3.
        let (graph, [first, second]) = ruled_graph([
            Spec::curve(sharp_kink(), 0.9, 0.3),
            Spec::curve(line(Point3::new(0.0, 0.0, 1.0), Vec3::X), 0.0, 1.0),
        ]);
        let surface = RelationSurface::ruled(&graph, &first, &second).unwrap();
        assert!(!bounds_of(surface, (0.6, 0.7), (0.0, 1.0)).smooth);
        assert!(bounds_of(surface, (0.4, 0.6), (0.0, 1.0)).smooth);
    }

    #[test]
    fn a_box_wholly_within_the_slack_past_the_end_is_clamped_to_it() {
        let (graph, id) = curve_graph(kinked());
        let past = 1.0 + 2.0 * Scalar::EPSILON;
        let oracle = RelationOracle::new(about_z(&graph, id));
        assert!(oracle.bounds((0.0, 1.0), (past, past)).is_some());
        let (graph, [first, second]) = ruled_graph(arc_and_line(false));
        let oracle = RelationOracle::new(RelationSurface::ruled(&graph, &first, &second).unwrap());
        assert!(oracle.bounds((past, past), (0.0, 1.0)).is_some());
    }

    #[test]
    fn reading_outside_the_unit_square_is_refused() {
        let (graph, [first, second]) = ruled_graph(arc_and_line(false));
        let surface = RelationSurface::ruled(&graph, &first, &second).unwrap();
        for (u, v) in [
            (1.0 + 1e-6, 0.5),
            (0.5, 1.0 + 1e-6),
            (-1e-6, 0.5),
            (0.5, -1e-6),
        ] {
            let error = surface.evaluate(u, v).unwrap_err();
            assert!(
                matches!(&error, GeomError::InvalidInput(text) if text == OUTSIDE_RULED_SQUARE),
                "{error:?}"
            );
        }
        assert!(surface.evaluate(1.0 + 1e-14, 0.5).is_ok());
        let oracle = RelationOracle::new(surface);
        assert!(oracle.bounds((0.5, 2.0), (0.0, 1.0)).is_some());
        assert!(oracle.bounds((1.5, 2.0), (0.0, 1.0)).is_none());
    }

    #[test]
    fn a_ruled_side_that_is_not_read_is_refused() {
        let polyline = Curve3::Polyline(Polyline3 {
            points: vec![Point3::ZERO, Point3::X],
            closed: false,
        });
        let (graph, [first, second]) =
            ruled_graph([Spec::curve(polyline, 0.0, 1.0), Spec::Apex(Point3::Z)]);
        assert!(RelationSurface::ruled(&graph, &first, &second).is_err());
    }
}
