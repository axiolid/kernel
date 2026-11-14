//! The surface a B-rep face lies on, as the tessellator reads it.
//!
//! The tessellator never matches a surface itself. It asks a [`Support`]
//! for a point, its parameter periods, its family name and whether two
//! faces share it, and a [`SupportOracle`] for the derivative bounds over
//! a parameter box. A support is an atomic [`axiolid_surface::Surface`] or
//! a surface relation read through a [`RelationSurface`] (#411, ADR 0103);
//! the view keeps the tessellator closed over one `Copy` type, so a further
//! support lands as one more arm here.

use crate::ruled::Side;
use crate::surface_relation::{RelationOracle, RelationSurface};
use axiolid_contracts::GeomResult;
use axiolid_core::{Point2, Point3, Scalar, Vec3};
use axiolid_model::SurfaceRelation;
use axiolid_reference::bound::{SurfaceBoundOracle, SurfaceBounds};
use axiolid_surface::Surface;

/// How a support collapses a line of its parameters to one point, a pole
/// (#463, ADR 0112): along `v = v_pole` every `u` names the same point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PoleKind {
    /// The support turns about its pole: a sphere's pole, a cone's apex, a
    /// revolution's profile on its axis. Its meridians (`u` fixed) are
    /// curves, everywhere orthogonal to its parallels, and `|S_v|` does not
    /// depend on `u`.
    Rotational,
    /// A ruled surface's side that is an apex: its rulings (`u` fixed) are
    /// straight.
    Apex,
}

/// What a face lies on.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Support<'a> {
    /// An atomic surface value.
    Surface(&'a Surface),
    /// A surface relation, with the node it was read from.
    Relation {
        node: &'a SurfaceRelation,
        surface: RelationSurface<'a>,
    },
}

impl<'a> Support<'a> {
    /// The atomic surface, when the support is one.
    pub(crate) fn atomic(self) -> Option<&'a Surface> {
        match self {
            Self::Surface(surface) => Some(surface),
            Self::Relation { .. } => None,
        }
    }

    /// Whether the support has poles of a kind the tessellator closes with
    /// a row of samples (#463): a sphere, a cone and a surface of revolution
    /// turn about theirs, a ruled surface with an apex side has straight
    /// rulings into it. `None` for every other support.
    pub(crate) fn pole_kind(self) -> Option<PoleKind> {
        match self {
            Self::Surface(Surface::Sphere(_) | Surface::Cone(_)) => Some(PoleKind::Rotational),
            Self::Surface(_) => None,
            Self::Relation { surface, .. } => match surface {
                RelationSurface::Revolution { .. } => Some(PoleKind::Rotational),
                RelationSurface::Ruled { sides } => sides
                    .iter()
                    .any(|side| matches!(side, Side::Apex(_)))
                    .then_some(PoleKind::Apex),
                RelationSurface::Extrusion { .. } => None,
            },
        }
    }

    /// Whether `p` is a pole of the support within `tau`, as the edge
    /// inverter reads one: a rotational support's `|S_u(p)|`, its distance
    /// from the axis per radian, at most `tau`; a ruled surface's point
    /// within `tau` of its apex.
    pub(crate) fn is_pole(self, p: Point2, tau: Scalar) -> bool {
        match self.pole_kind() {
            None => false,
            Some(PoleKind::Rotational) => self
                .partials(p.x, p.y)
                .is_ok_and(|(su, _)| su.length() <= tau),
            Some(PoleKind::Apex) => {
                let Self::Relation {
                    surface: RelationSurface::Ruled { sides },
                    ..
                } = self
                else {
                    return false;
                };
                let Ok(point) = self.evaluate(p.x, p.y) else {
                    return false;
                };
                sides.iter().any(|side| match side {
                    Side::Apex(apex) => (point - *apex).length() <= tau,
                    Side::Curve { .. } => false,
                })
            }
        }
    }

    /// The support's point at `(u, v)`.
    pub(crate) fn evaluate(self, u: Scalar, v: Scalar) -> GeomResult<Point3> {
        match self {
            Self::Surface(surface) => axiolid_reference::surface::evaluate(surface, u, v),
            Self::Relation { surface, .. } => surface.evaluate(u, v),
        }
    }

    /// The support's tangents `S_u` and `S_v` at `(u, v)`.
    pub(crate) fn partials(self, u: Scalar, v: Scalar) -> GeomResult<(Vec3, Vec3)> {
        match self {
            Self::Surface(surface) => axiolid_reference::surface::partials(surface, u, v),
            Self::Relation { surface, .. } => surface.tangents(u, v),
        }
    }

    /// The parameter range of each open axis, which a preimage never
    /// leaves; unbounded elsewhere.
    pub(crate) fn open_domain(self) -> ((Scalar, Scalar), (Scalar, Scalar)) {
        match self {
            Self::Surface(surface) => super::inversion::parameter_domain(surface),
            Self::Relation { surface, .. } => surface.open_domain(),
        }
    }

    /// How far `point` lies from the support's poles: from the axis a
    /// rotational support turns about, or from a ruled surface's apex.
    /// `None` without a [`PoleKind`]. It changes no faster than the point
    /// moves, so an edge of length `L` between points at `d_a` and `d_b`
    /// keeps off the poles when `d_a + d_b > L`.
    pub(crate) fn pole_distance(self, point: Point3) -> Option<Scalar> {
        let from_axis = |origin: Point3, axis: Vec3| {
            let (offset, axis) = (point - origin, axis.normalize_or_zero());
            (offset - axis * offset.dot(axis)).length()
        };
        match self {
            Self::Surface(Surface::Sphere(sphere)) => {
                Some(from_axis(sphere.frame.origin, sphere.frame.z))
            }
            Self::Surface(Surface::Cone(cone)) => Some(from_axis(cone.frame.origin, cone.frame.z)),
            Self::Surface(_) => None,
            Self::Relation { surface, .. } => match surface {
                RelationSurface::Revolution { origin, axis, .. } => Some(from_axis(origin, axis)),
                RelationSurface::Ruled { sides } => sides
                    .iter()
                    .filter_map(|side| match side {
                        Side::Apex(apex) => Some((point - *apex).length()),
                        Side::Curve { .. } => None,
                    })
                    .reduce(Scalar::min),
                RelationSurface::Extrusion { .. } => None,
            },
        }
    }

    /// The range of `v` the support's points are inverted into, which an
    /// edge's preimage never leaves: a sphere's latitudes, a cone's side of
    /// its apex, an open axis's domain; unbounded elsewhere (#463).
    pub(crate) fn v_range(self) -> (Scalar, Scalar) {
        let unbounded = (Scalar::NEG_INFINITY, Scalar::INFINITY);
        match self {
            Self::Surface(Surface::Sphere(_)) => {
                let half = core::f64::consts::FRAC_PI_2;
                (-half, half)
            }
            Self::Surface(Surface::Cone(cone)) => {
                let slope = cone.semi_angle.tan();
                if !(slope.is_finite() && slope != 0.0) {
                    return unbounded;
                }
                let apex = -cone.radius / slope;
                if slope > 0.0 {
                    (apex, Scalar::INFINITY)
                } else {
                    (Scalar::NEG_INFINITY, apex)
                }
            }
            _ => self.open_domain().1,
        }
    }

    /// The support's parameter periods.
    pub(crate) fn periods(self) -> (Option<Scalar>, Option<Scalar>) {
        match self {
            Self::Surface(surface) => super::surface_periods(surface),
            Self::Relation { surface, .. } => surface.periods(),
        }
    }

    /// The support's family name, for a report.
    pub(crate) fn family(self) -> &'static str {
        match self {
            Self::Surface(surface) => super::surface_family(surface),
            Self::Relation { surface, .. } => surface.family(),
        }
    }

    /// Whether both are the same node, by identity.
    pub(crate) fn same(self, other: Self) -> bool {
        match (self, other) {
            (Self::Surface(a), Self::Surface(b)) => std::ptr::eq(a, b),
            (Self::Relation { node: a, .. }, Self::Relation { node: b, .. }) => std::ptr::eq(a, b),
            _ => false,
        }
    }
}

/// The derivative bounds of a [`Support`] over parameter boxes.
#[derive(Debug, Clone)]
pub(crate) enum SupportOracle<'a> {
    /// An atomic surface's oracle.
    Surface(SurfaceBoundOracle<'a>),
    /// A surface relation's oracle.
    Relation(RelationOracle<'a>),
}

impl<'a> SupportOracle<'a> {
    /// The oracle of `support`, or `None` when it has no derivative bound.
    pub(crate) fn new(support: Support<'a>) -> Option<Self> {
        match support {
            Support::Surface(surface) => SurfaceBoundOracle::new(surface).map(Self::Surface),
            Support::Relation { surface, .. } => Some(Self::Relation(RelationOracle::new(surface))),
        }
    }

    /// The bounds over the parameter box, or `None` when it has none.
    pub(crate) fn bounds(&self, u: (Scalar, Scalar), v: (Scalar, Scalar)) -> Option<SurfaceBounds> {
        match self {
            Self::Surface(oracle) => oracle.bounds(u, v),
            Self::Relation(oracle) => oracle.bounds(u, v),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::{Frame3, Interval};
    use axiolid_curve::{Circle3, Curve3, Line3};
    use axiolid_model::{GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, RuledSide};
    use axiolid_surface::{Cone, Cylinder, Sphere, Torus};
    use std::f64::consts::{FRAC_PI_2, PI};

    const TAU: Scalar = 1e-3;

    fn frame() -> Frame3 {
        Frame3 {
            origin: Point3::new(0.2, -0.1, 0.3),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        }
    }

    fn graph_of(curves: Vec<Curve3>) -> (GeometryGraph, Vec<NodeId>) {
        let mut builder = GeometryGraphBuilder::new();
        let ids: Vec<NodeId> = curves
            .into_iter()
            .map(|curve| builder.push(GeometryNode::Curve3(curve)).unwrap())
            .collect();
        (builder.finish(vec![ids[0]]).unwrap(), ids)
    }

    fn line(origin: Point3, direction: Vec3) -> Curve3 {
        Curve3::Line(Line3 { origin, direction })
    }

    fn circle() -> Curve3 {
        Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius: 1.0,
        })
    }

    /// Any node: a support compares nodes by identity only.
    fn node(id: NodeId) -> SurfaceRelation {
        SurfaceRelation::LinearExtrusion {
            swept_curve: id,
            direction: Vec3::Z,
        }
    }

    fn relation<'a>(node: &'a SurfaceRelation, surface: RelationSurface<'a>) -> Support<'a> {
        Support::Relation { node, surface }
    }

    #[test]
    fn a_sphere_s_poles_and_a_cone_s_apex_are_poles() {
        let radius = 0.5;
        let sphere = Surface::Sphere(Sphere {
            frame: frame(),
            radius,
        });
        let support = Support::Surface(&sphere);
        assert_eq!(support.pole_kind(), Some(PoleKind::Rotational));
        for v in [FRAC_PI_2, -FRAC_PI_2] {
            assert!(support.is_pole(Point2::new(0.7, v), TAU), "{v}");
        }
        assert!(!support.is_pole(Point2::new(0.7, 0.0), TAU));
        // Half the tolerance from the axis is a pole, one and a half is not.
        assert!(support.is_pole(Point2::new(0.7, (0.5 * TAU / radius).acos()), TAU));
        assert!(!support.is_pole(Point2::new(0.7, (1.5 * TAU / radius).acos()), TAU));

        let (radius, semi_angle) = (0.4, -0.4_f64);
        let cone = Surface::Cone(Cone {
            frame: frame(),
            radius,
            semi_angle,
        });
        let support = Support::Surface(&cone);
        assert_eq!(support.pole_kind(), Some(PoleKind::Rotational));
        let apex = -radius / semi_angle.tan();
        assert!(support.is_pole(Point2::new(2.0, apex), TAU));
        assert!(!support.is_pole(Point2::new(2.0, 0.0), TAU));
    }

    #[test]
    fn a_revolved_line_meeting_its_axis_has_a_pole_there() {
        let (graph, ids) = graph_of(vec![line(Point3::ZERO, Vec3::new(0.2, 0.0, 1.0))]);
        let node = node(ids[0]);
        let surface = RelationSurface::revolution(&graph, ids[0], Point3::ZERO, Vec3::Z).unwrap();
        let support = relation(&node, surface);
        assert_eq!(support.pole_kind(), Some(PoleKind::Rotational));
        assert!(support.is_pole(Point2::new(1.0, 0.0), TAU));
        assert!(!support.is_pole(Point2::new(1.0, 0.5), TAU));
    }

    #[test]
    fn a_ruled_surface_has_a_pole_at_an_apex_on_either_side() {
        let (graph, ids) = graph_of(vec![circle()]);
        let apex = Point3::new(0.0, 0.0, 1.0);
        let curve = RuledSide::Curve {
            curve: ids[0],
            interval: Interval {
                start: 0.0,
                end: PI,
            },
        };
        let node = node(ids[0]);
        for (apex_v, first, second) in [
            (0.0, RuledSide::Apex(apex), curve),
            (1.0, curve, RuledSide::Apex(apex)),
        ] {
            let surface = RelationSurface::ruled(&graph, &first, &second).unwrap();
            let support = relation(&node, surface);
            assert_eq!(support.pole_kind(), Some(PoleKind::Apex), "{apex_v}");
            assert!(support.is_pole(Point2::new(0.3, apex_v), TAU), "{apex_v}");
            assert!(!support.is_pole(Point2::new(0.3, 0.5), TAU), "{apex_v}");
            assert!(
                !support.is_pole(Point2::new(0.3, 1.0 - apex_v), TAU),
                "{apex_v}"
            );
        }
    }

    #[test]
    fn surfaces_without_a_pole_have_none() {
        let cylinder = Surface::Cylinder(Cylinder {
            frame: frame(),
            radius: 0.3,
        });
        let torus = Surface::Torus(Torus {
            frame: frame(),
            major_radius: 1.0,
            minor_radius: 0.2,
        });
        for surface in [&cylinder, &torus] {
            let support = Support::Surface(surface);
            assert_eq!(support.pole_kind(), None);
            assert!(!support.is_pole(Point2::new(0.1, 0.1), TAU));
        }
        let (graph, ids) = graph_of(vec![circle(), line(Point3::new(0.0, 0.0, 1.0), Vec3::X)]);
        let node = node(ids[0]);
        let extrusion = RelationSurface::extrusion(&graph, ids[0], Vec3::Z).unwrap();
        assert_eq!(relation(&node, extrusion).pole_kind(), None);
        let side = |i: usize, end: Scalar| RuledSide::Curve {
            curve: ids[i],
            interval: Interval { start: 0.0, end },
        };
        let ruled = RelationSurface::ruled(&graph, &side(0, PI), &side(1, 1.0)).unwrap();
        let support = relation(&node, ruled);
        assert_eq!(support.pole_kind(), None);
        assert!(!support.is_pole(Point2::new(0.5, 0.0), TAU));
    }
}
