//! Half-space clipping of placed solids (#234).
//!
//! A building model clips a wall by a roof plane: the wall minus the side of
//! the plane above the roof, and a gable wall is clipped by two. The clip is
//! a boolean whose tool is a half-space (`GeometryNode::HalfSpace`), or a
//! half-space bounded by a polygon or a profile in its plane
//! (`SolidOperation::BoundedHalfSpace`), either read through rigid
//! placements. Both keep the mesh compiler's semantics: `agreement`
//! selects the normal side of the plane as the half-space; a bounded one is
//! that side intersected with the infinite prism of its boundary, swept
//! along the plane normal, the boundary framed by the operation's placement
//! projected into the plane, its scale, shear and handedness kept (#322), or
//! along the placement's `z` (#321, see below). A frame that scales or
//! mirrors a profile boundary is applied as a derived profile, the part of
//! the map that keeps handedness, and the prism reflected where it does
//! not, so a similarity's arcs stay exact cylinder walls; any other map
//! makes them ellipses, refused by name.
//!
//! # Prisms along the placement's `z` (#321)
//!
//! `SolidOperation::BoundedHalfSpaceAlongZ` sweeps the boundary, as written in the
//! placement's `xy` plane, along the placement's `z`, at any angle to the
//! clip plane. Its tool is built in two steps
//! (`ExactCompilation::placement_z_tool`): the prism in the boundary's own
//! plane, the boundary as authored there (an arc a circle; a scaled or
//! sheared frame applied as a derived profile, which keeps a similarity's
//! arcs arcs and refuses the ellipses of any other map by name), extruded
//! along the axis from twice the margin below the envelope to twice the
//! margin above it; then its intersection with the unbounded half-space's
//! prism over the envelope, through the general boolean, whose report is
//! merged into the clip's. A vertical prism under a sloped plane, the IFC
//! gable, keeps right cylinder walls, each cut by the plane in an ellipse.
//! A `z` leaning off the boundary's plane sweeps elliptical cylinders
//! (#287), which the general boolean may refuse to sew, by name.
//!
//! # Boundaries with arcs (#277)
//!
//! A boundary of line and circular-arc segments is a profile node, and the
//! prism is the exact extrusion of that profile
//! ([`extrude_profile_exact`]): a line becomes a plane wall, an arc a right
//! circular cylinder wall, exactly, never chords. The sweep is along the
//! plane normal, perpendicular to the profile, so no arc wall is oblique.
//! The profile is checked first, as the mesh compiler checks it
//! (`crate::half_space_boundary`): open, self-crossing or zero-radius
//! contours are refused by name. Like the polygon, the profile is mirrored
//! in its `x` axis when the kept side is opposite the normal, so that the
//! prism's frame stays a rotation: its prism, extruded as authored, is
//! reflected in the profile's `x` axis by the same placement
//! (`ExactBRep::transformed`), and the general boolean cuts by the
//! reflected cylinder wall as by one built directly (#288).
//!
//! # A finite tool through the general boolean
//!
//! The half-space is replaced by a finite prism, and the clip is the general
//! exact boolean (ADR 0075, ADR 0080) of the subject and that prism: the
//! subject's faces against one plane face, under the same tolerance contract
//! as a placed opening. A dedicated plane-clip of an exact B-rep was not
//! written: it would be a second implementation of section, split and
//! classification to keep in agreement with the general one, and a clip
//! must compose with the openings the general boolean already cuts (a wall
//! minus its windows, clipped by its roof, in either order).
//!
//! The prism stands on the clip plane and covers the subject's envelope (a
//! sound box around it, see below) plus a margin `m`. For an unbounded
//! half-space its footprint is the envelope's projection into the plane,
//! widened by `m`; for a bounded one it is the boundary. It reaches
//! `m` past the envelope's farthest point on the kept side of the plane.
//!
//! **Why the result does not depend on the margin.** Write the half-space as
//! `H` and the prism as `P = H ∩ C`, where `C` is the footprint's prism along
//! the normal cut off `m` beyond the envelope `E`. The subject `S` lies in
//! `E`, and `E ∩ H` lies in `C` minus its far cap and (unbounded case) minus
//! its sides, so `S ∩ P = S ∩ H` and `S - P = S - H` as sets, for every
//! `m > 0`. The faces `m` controls (the far cap, and the sides of the
//! unbounded prism) stay at least `m` from `S`. With `m` above the linear
//! tolerance `eps` no within-tolerance decision of the general boolean
//! involves them (its bounding boxes, enlarged by `eps`, already skip
//! them), so every decision it takes is about the plane face (and, bounded,
//! the boundary's side faces) against the subject's faces: the same for
//! every margin. `m` is a quarter of the envelope's diagonal plus `4 eps`.
//!
//! **The envelope.** The subject is a placed extrusion, or a difference or
//! clip of one; neither removes anything outside its first operand, so the
//! envelope of a difference or clip is its subject's, of any other boolean
//! of placed operands both operands', and of an instance its source's box
//! mapped by the placement. An extrusion's faces are its caps and walls
//! ruled along the extrusion, so it lies in the box of its edges: lines by
//! their vertices, circles and ellipses in closed form, polylines and
//! B-splines by their points and control points.
//!
//! When the envelope lies wholly off the half-space the clip is decided
//! there, without a boolean: a difference keeps the subject unchanged and
//! an intersection is empty. An unbounded half-space containing the whole
//! envelope likewise empties a difference and keeps the subject of an
//! intersection.
//!
//! # Refused, by name
//!
//! A union with a half-space (unbounded), a half-space as the subject, a
//! subject not built from placed extrusions, a bounded half-space whose
//! boundary curve is not a polyline, whose profile does not bound a region
//! or cannot be extruded exactly, or whose plane is placed by an instance,
//! one along its placement's `z` whose axis lies within the angular
//! tolerance of the plane or whose placement is degenerate, a scaled or
//! sheared placement above the half-space, an envelope with an edge family
//! it cannot bound, and every refusal of the general boolean. An emptied
//! subject is `GeomError::Degenerate`. The model's half-space carries a
//! plane, so a curved base surface cannot reach this path.

use axiolid_brep::{ExactBRep, TransformError};
use axiolid_brep_boolean::{boolean_with_report, BooleanError, BooleanReport};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::half_space::BoundaryPrism;
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{
    BooleanOperator, Interval, Mat3, Point2, Point3, Scalar, Transform2, Transform3, Vec2, Vec3,
};
use axiolid_curve::{Curve2, Curve3, Line2};
use axiolid_model::{GeometryNode, NodeId, SolidOperation};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile};

use super::boolean::remap_boolean_error;
use super::{remap_construction_error, unsupported, ExactCompilation};
use crate::half_space_boundary::PrismAxis;

/// A boolean tool that is a half-space, read through its placements.
#[derive(Debug, Clone)]
pub(super) struct ClipTool {
    clip: Clip,
    /// The placements above the half-space, composed.
    placement: Transform3,
}

#[derive(Debug, Clone)]
enum Clip {
    /// One side of a plane.
    Unbounded(HalfSpace),
    /// One side of a plane within the prism of a polyline in the plane.
    Bounded {
        half_space: HalfSpace,
        boundary: NodeId,
        /// The boundary's own frame.
        frame: Transform3,
        /// Which way its prism runs (#321).
        axis: PrismAxis,
    },
}

/// A bounded half-space's boundary, in its own frame.
#[derive(Debug, Clone)]
enum Footprint {
    /// A closed polyline's distinct points, as wound.
    Polygon(Vec<Point2>),
    /// A profile: lines and circular arcs, exactly (#277).
    Profile(Profile),
}

/// What the envelope alone decides.
enum Decided {
    /// The subject, unchanged.
    Subject,
    /// Nothing is left.
    Empty,
}

/// A sound axis-aligned box.
#[derive(Debug, Clone, Copy)]
struct Envelope {
    lo: Point3,
    hi: Point3,
}

impl Envelope {
    const EMPTY: Self = Self {
        lo: Vec3::splat(Scalar::INFINITY),
        hi: Vec3::splat(Scalar::NEG_INFINITY),
    };

    fn include(&mut self, point: Point3) {
        self.lo = self.lo.min(point);
        self.hi = self.hi.max(point);
    }

    fn union(self, other: Self) -> Self {
        Self {
            lo: self.lo.min(other.lo),
            hi: self.hi.max(other.hi),
        }
    }

    fn corners(&self) -> [Point3; 8] {
        let (a, b) = (self.lo, self.hi);
        [
            Vec3::new(a.x, a.y, a.z),
            Vec3::new(b.x, a.y, a.z),
            Vec3::new(a.x, b.y, a.z),
            Vec3::new(b.x, b.y, a.z),
            Vec3::new(a.x, a.y, b.z),
            Vec3::new(b.x, a.y, b.z),
            Vec3::new(a.x, b.y, b.z),
            Vec3::new(b.x, b.y, b.z),
        ]
    }

    /// The box around this one's corners under `transform`.
    fn mapped(&self, transform: &Transform3) -> Self {
        let mut out = Self::EMPTY;
        for corner in self.corners() {
            out.include(transform.transform_point3(corner));
        }
        out
    }
}

impl ExactCompilation<'_> {
    /// The half-space `id` is, through its placements, or `None` when it is
    /// not one.
    pub(super) fn clip_tool(&self, mut id: NodeId) -> GeomResult<Option<ClipTool>> {
        let mut placement = Transform3::IDENTITY;
        for _ in 0..=self.graph.len() {
            let node = self.graph.get(id).ok_or_else(|| {
                GeomError::InvalidInput(format!("operand {id:?} does not belong to this graph"))
            })?;
            let clip = match node {
                GeometryNode::Instance(instance) => {
                    placement *= instance.transform;
                    id = instance.source;
                    continue;
                }
                GeometryNode::HalfSpace(half_space) => Clip::Unbounded(*half_space),
                GeometryNode::SolidOperation(operation) => {
                    let Some((half_space, boundary, frame, axis)) =
                        crate::half_space_boundary::bounded(operation)
                    else {
                        return Ok(None);
                    };
                    let Some(GeometryNode::HalfSpace(half_space)) = self.graph.get(half_space)
                    else {
                        return Err(unsupported(
                            "exact bounded half-space whose plane is placed by an instance",
                        ));
                    };
                    Clip::Bounded {
                        half_space: *half_space,
                        boundary,
                        frame,
                        axis,
                    }
                }
                _ => return Ok(None),
            };
            return Ok(Some(ClipTool { clip, placement }));
        }
        Err(GeomError::InvalidInput(
            "exact boolean operand placements form a cycle".to_owned(),
        ))
    }

    /// `left` clipped by the half-space `tool`.
    pub(super) fn compile_clip(
        &mut self,
        left: NodeId,
        tool: ClipTool,
        operator: BooleanOperator,
    ) -> GeomResult<ExactBRep> {
        if operator == BooleanOperator::Union {
            return Err(unsupported(
                "exact union with a half-space, which is unbounded",
            ));
        }
        if self.clip_tool(left)?.is_some() {
            return Err(unsupported("exact boolean whose subject is a half-space"));
        }
        let subject = self.compile(left)?;
        let envelope = self.envelope(left)?;
        let tolerance = self.options.tolerance();
        let (half_space, footprint) = match tool.clip {
            Clip::Unbounded(half_space) => (half_space, None),
            Clip::Bounded {
                half_space,
                boundary,
                frame,
                axis,
            } => (half_space, Some((self.boundary(boundary)?, frame, axis))),
        };
        let (solid, tool_report) = match self.clip_solid(
            &half_space,
            footprint.as_ref(),
            &tool.placement,
            envelope,
            operator,
        )? {
            Ok(solid) => solid,
            Err(Decided::Subject) => return Ok(subject),
            Err(Decided::Empty) => return Err(emptied()),
        };
        let (body, report) = boolean_with_report(&subject, &solid, operator, tolerance).map_err(
            |error| match error {
                BooleanError::EmptyResult => emptied(),
                error => remap_boolean_error(error),
            },
        )?;
        // What it read within tolerance, for the compiler's report (#236),
        // building the tool included.
        let report = match tool_report {
            Some(tool) => tool.merged(&report),
            None => report,
        };
        self.boolean_report = Some(report);
        Ok(body)
    }

    /// The finite prism standing in for the half-space over `envelope`, or
    /// what the envelope decides alone; with what building it read within
    /// tolerance, when that took a boolean.
    fn clip_solid(
        &self,
        half_space: &HalfSpace,
        footprint: Option<&(Footprint, Transform3, PrismAxis)>,
        placement: &Transform3,
        envelope: Envelope,
        operator: BooleanOperator,
    ) -> GeomResult<Result<(ExactBRep, Option<BooleanReport>), Decided>> {
        let tolerance = self.options.tolerance();
        let origin = half_space.boundary.origin;
        let normal = half_space.boundary.normal.normalize_or_zero();
        if !origin.is_finite() || !normal.is_finite() || normal == Vec3::ZERO {
            return Err(GeomError::InvalidInput(
                "half-space boundary plane needs a finite point and a non-zero normal".to_owned(),
            ));
        }
        // Everything is built in the half-space's own frame; the subject's
        // envelope is brought into it.
        let inverse = placement.inverse();
        if !inverse.matrix3.is_finite() || !inverse.translation.is_finite() {
            return Err(unsupported(
                "exact half-space under a scaled or sheared transform",
            ));
        }
        let corners = envelope.mapped(&inverse).corners();
        let side = if half_space.agreement {
            normal
        } else {
            -normal
        };
        let (near, far) = corners.iter().fold(
            (Scalar::INFINITY, Scalar::NEG_INFINITY),
            |(near, far), &corner| {
                let distance = (corner - origin).dot(side);
                (near.min(distance), far.max(distance))
            },
        );
        let keeps = operator == BooleanOperator::Intersection;
        if far < 0.0 {
            // The subject lies wholly off the half-space.
            return Ok(Err(if keeps {
                Decided::Empty
            } else {
                Decided::Subject
            }));
        }
        if footprint.is_none() && near > 0.0 {
            // The subject lies wholly in the half-space.
            return Ok(Err(if keeps {
                Decided::Subject
            } else {
                Decided::Empty
            }));
        }
        let margin = 0.25 * (envelope.hi - envelope.lo).length() + 4.0 * tolerance.linear();
        let footprint = match footprint {
            Some((footprint, frame, PrismAxis::PlacementZ)) => {
                return self.placement_z_tool(
                    half_space,
                    (footprint, *frame),
                    placement,
                    envelope.mapped(&inverse),
                    operator,
                );
            }
            Some((footprint, frame, PrismAxis::PlaneNormal)) => Some((footprint, frame)),
            None => None,
        };

        // The in-plane frame: `x` and `y = normal x x` as the mesh compiler
        // frames the boundary, then `y` turned with the kept side so that
        // `(x, y_side, side)` is a rotation and the prism runs along `+z`.
        // A bounded half-space's placement is read whole, as the mesh
        // compiler reads it (#322): `shape` maps the boundary into
        // `(x, normal x x)`, the placement's scale, shear and handedness kept.
        let (anchor, x, shape) = match footprint {
            None => {
                let reference = if normal.x.abs() < 0.9 {
                    Vec3::X
                } else {
                    Vec3::Y
                };
                (
                    origin,
                    reference - normal * normal.dot(reference),
                    Transform2::IDENTITY,
                )
            }
            Some((_, &frame)) => {
                let prism =
                    BoundaryPrism::along_plane_normal(frame, half_space.boundary, tolerance)?;
                let x = prism.x.normalize();
                let shape = Transform2::from_cols(
                    Vec2::new(prism.x.dot(x), 0.0),
                    Vec2::new(prism.y.dot(x), prism.y.dot(normal.cross(x))),
                    Vec2::ZERO,
                );
                let shape = if is_identity(&shape) {
                    Transform2::IDENTITY
                } else {
                    shape
                };
                (prism.origin, x, shape)
            }
        };
        let x = x.normalize_or_zero();
        if x == Vec3::ZERO {
            return Err(GeomError::InvalidInput(
                "authored boundary frame x axis is parallel to the clip plane normal, \
                 so it fixes no in-plane direction"
                    .to_owned(),
            ));
        }
        let flip = if half_space.agreement { 1.0 } else { -1.0 };
        let y = normal.cross(x) * flip;
        // The boundary into `(x, y)`: its shape, then `y` turned with the
        // side. A profile is extruded under the part of that map that keeps
        // its handedness, and its prism reflected in the profile's `x` axis
        // where the map reflects (#288), so its arcs stay arcs.
        let into = Transform2::from_cols(Vec2::X, Vec2::new(0.0, flip), Vec2::ZERO) * shape;
        let reflects = into.matrix2.determinant() < 0.0;
        let kept = if reflects {
            Transform2::from_cols(Vec2::X, -Vec2::Y, Vec2::ZERO) * into
        } else {
            into
        };
        let (profile, centre) = match footprint {
            None => {
                let (mut lo, mut hi) = (
                    Vec2::splat(Scalar::INFINITY),
                    Vec2::splat(-Scalar::INFINITY),
                );
                for corner in corners {
                    let at = Vec2::new((corner - anchor).dot(x), (corner - anchor).dot(y));
                    lo = lo.min(at);
                    hi = hi.max(at);
                }
                let size = hi - lo + Vec2::splat(2.0 * margin);
                let profile = Profile::Rectangle(RectangleProfile {
                    x: size.x,
                    y: size.y,
                    thickness: None,
                    outer_radius: None,
                    inner_radius: None,
                });
                (profile, (lo + hi) * 0.5)
            }
            Some((Footprint::Polygon(points), _)) => {
                let points: Vec<Point2> = points
                    .iter()
                    .map(|&p| {
                        let p = shape.transform_point2(p);
                        Point2::new(p.x, p.y * flip)
                    })
                    .collect();
                (polygon(&points), Vec2::ZERO)
            }
            // A profile (#277) is extruded as authored, under its kept map.
            // Every arc of it becomes a right circular cylinder wall: the
            // sweep is along the plane normal, never oblique. A similarity
            // keeps it an arc; any other map makes it an ellipse, which the
            // exact extruder refuses by name.
            Some((Footprint::Profile(profile), _)) if is_identity(&kept) => {
                (profile.clone(), Vec2::ZERO)
            }
            Some((Footprint::Profile(profile), _)) => (
                Profile::Derived {
                    basis: Box::new(profile.clone()),
                    transform: kept,
                },
                Vec2::ZERO,
            ),
        };
        let prism = extrude_profile_exact(&profile, Vec3::Z, far + margin, tolerance)
            .map_err(remap_construction_error)?;
        let base = anchor + x * centre.x + y * centre.y;
        let frame = Transform3::from_cols(x, y, side, base);
        // A profile turns over with `y` as the polygon does, by reflecting
        // its prism in the profile's `x` axis (#288).
        let turned = match footprint {
            Some((Footprint::Profile(_), _)) if reflects => {
                frame * Transform3::from_mat3(Mat3::from_diagonal(Vec3::new(1.0, -1.0, 1.0)))
            }
            _ => frame,
        };
        prism
            .transformed(&(*placement * turned))
            .map(|solid| Ok((solid, None)))
            .map_err(transform_error)
    }

    /// The tool of a half-space bounded along its placement's `z` (#321),
    /// built in the half-space's own frame over `local`, the envelope
    /// brought into it, then placed.
    ///
    /// The prism need not stand perpendicular to the plane, so it is not
    /// stood on the plane: its cap there would be the boundary carried
    /// obliquely onto the plane, each arc an ellipse no exact extrusion
    /// sweeps along a leaning axis. It is built instead in the boundary's
    /// own plane, the boundary as authored (an arc a circle), extruded
    /// along the axis -- a right prism, or an oblique one under a sheared
    /// placement, arcs sweeping its elliptical cylinders (#287) -- from a
    /// margin below the envelope to a margin above it, measured across that
    /// plane. Its intersection with the unbounded half-space's tool over
    /// the same envelope is then the tool: on the envelope, the side of the
    /// plane within the prism. The margin argument of the module docs holds
    /// for every face either adds. Where the envelope lies wholly in the
    /// half-space, the prism alone is the tool; where the prism misses the
    /// side within the margin, the envelope decides.
    fn placement_z_tool(
        &self,
        half_space: &HalfSpace,
        (footprint, frame): (&Footprint, Transform3),
        placement: &Transform3,
        local: Envelope,
        operator: BooleanOperator,
    ) -> GeomResult<Result<(ExactBRep, Option<BooleanReport>), Decided>> {
        let tolerance = self.options.tolerance();
        let prism = BoundaryPrism::along_placement_z(frame, tolerance)?;
        prism.crossing(half_space.boundary, tolerance)?;
        // A rigid frame on the boundary's plane, `z` across it on the side
        // `x` and `y` turn to, so a mirrored frame keeps its boundary
        // unmirrored in the plane and the frame a rotation.
        let x = prism.x.normalize();
        let z = prism.x.cross(prism.y).normalize();
        let y = z.cross(x);
        // The boundary in that frame: scaled, sheared, never reflected.
        let shape = Transform2::from_cols(
            Vec2::new(prism.x.dot(x), 0.0),
            Vec2::new(prism.y.dot(x), prism.y.dot(y)),
            Vec2::ZERO,
        );
        let profile = match footprint {
            Footprint::Polygon(points) => polygon(
                &points
                    .iter()
                    .map(|&p| shape.transform_point2(p))
                    .collect::<Vec<_>>(),
            ),
            // An arc stays an arc under a similarity; under any other map
            // it is an ellipse, which the exact extruder refuses by name.
            Footprint::Profile(profile) if is_identity(&shape) => profile.clone(),
            Footprint::Profile(profile) => Profile::Derived {
                basis: Box::new(profile.clone()),
                transform: shape,
            },
        };
        // The axis in that frame, turned across the plane to `+z`; a lean
        // below the rounding of the frame's own axes is none.
        let axis = prism.axis.normalize();
        let mut along = Vec3::new(axis.dot(x), axis.dot(y), axis.dot(z));
        if along.z < 0.0 {
            along = -along;
        }
        if along.truncate().length() <= 8.0 * Scalar::EPSILON {
            along = Vec3::Z;
        }
        // Twice the side's margin, so the prism's caps pass the side's far
        // cap and sides rather than meet them in one plane, a coincidence
        // the boolean would have to read for nothing.
        let margin = 0.5 * (local.hi - local.lo).length() + 8.0 * tolerance.linear();
        let (lo, hi) = local.corners().iter().fold(
            (Scalar::INFINITY, Scalar::NEG_INFINITY),
            |(lo, hi), &corner| {
                let height = (corner - prism.origin).dot(z);
                (lo.min(height), hi.max(height))
            },
        );
        let (lo, hi) = (lo - margin, hi + margin);
        let long = extrude_profile_exact(&profile, along, (hi - lo) / along.z, tolerance)
            .map_err(remap_construction_error)?;
        let start = prism.origin + (x * along.x + y * along.y + z * along.z) * (lo / along.z);
        let long = long
            .transformed(&Transform3::from_cols(x, y, z, start))
            .map_err(transform_error)?;
        // The side of the plane over the same envelope, unless it holds the
        // envelope whole.
        let side =
            match self.clip_solid(half_space, None, &Transform3::IDENTITY, local, operator)? {
                Ok((side, _)) => side,
                Err(_) => {
                    return long
                        .transformed(placement)
                        .map(|solid| Ok((solid, None)))
                        .map_err(transform_error)
                }
            };
        let (tool, report) =
            match boolean_with_report(&long, &side, BooleanOperator::Intersection, tolerance) {
                Ok(tool) => tool,
                // The prism misses the plane's side over the envelope.
                Err(BooleanError::EmptyResult) => {
                    return Ok(Err(if operator == BooleanOperator::Intersection {
                        Decided::Empty
                    } else {
                        Decided::Subject
                    }))
                }
                Err(error) => return Err(remap_boolean_error(error)),
            };
        tool.transformed(placement)
            .map(|solid| Ok((solid, Some(report))))
            .map_err(transform_error)
    }

    /// A bounded half-space's boundary, read as the mesh compiler reads it.
    fn boundary(&self, id: NodeId) -> GeomResult<Footprint> {
        let curve = match self.graph.get(id) {
            Some(GeometryNode::Curve2(curve)) => curve,
            // Lines and circular arcs, exactly (#277); whatever the exact
            // extruder refuses of it is refused by its name.
            Some(GeometryNode::Profile(profile)) => {
                crate::half_space_boundary::check(profile, self.options.tolerance())?;
                return Ok(Footprint::Profile(profile.clone()));
            }
            _ => {
                return Err(GeomError::InvalidInput(format!(
                    "half-space boundary {id:?} is neither a Curve2 nor a Profile node"
                )))
            }
        };
        let Curve2::Polyline(polyline) = curve else {
            return Err(unsupported(
                "exact bounded half-space whose boundary is not a polyline",
            ));
        };
        let mut points = polyline.points.clone();
        // A closed polyline may or may not repeat its first point.
        if points.len() >= 2 && points[0] == points[points.len() - 1] {
            points.pop();
        }
        if points.len() < 3 {
            return Err(GeomError::InvalidInput(
                "half-space boundary needs at least 3 distinct points".to_owned(),
            ));
        }
        Ok(Footprint::Polygon(points))
    }

    /// A sound box around what `id` compiles to.
    fn envelope(&mut self, id: NodeId) -> GeomResult<Envelope> {
        let node = self.graph.get(id).ok_or_else(|| {
            GeomError::InvalidInput(format!("operand {id:?} does not belong to this graph"))
        })?;
        match node {
            GeometryNode::Instance(instance) => {
                let instance = *instance;
                Ok(self.envelope(instance.source)?.mapped(&instance.transform))
            }
            GeometryNode::SolidOperation(SolidOperation::Extrusion { .. }) => {
                let solid = self.compile(id)?;
                edge_box(&solid)
            }
            GeometryNode::SolidOperation(SolidOperation::Boolean {
                left,
                right,
                operator,
            }) => {
                let (left, right, operator) = (*left, *right, *operator);
                // A difference or a clip removes nothing outside its subject.
                if operator == BooleanOperator::Difference || self.clip_tool(right)?.is_some() {
                    self.envelope(left)
                } else {
                    Ok(self.envelope(left)?.union(self.envelope(right)?))
                }
            }
            // A mesh flattened for a boolean's report (#266): every edge of
            // its polyhedron is a segment, so its edges' box holds it.
            GeometryNode::PolygonMesh(_) | GeometryNode::TriMesh(_) | GeometryNode::BRep(_)
                if self.flatten =>
            {
                let solid = self.compile(id)?;
                edge_box(&solid)
            }
            _ => Err(unsupported(
                "exact half-space clip of a subject not built from placed extrusions",
            )),
        }
    }
}

/// The box of a solid whose faces are bounded by their edges.
fn edge_box(solid: &ExactBRep) -> GeomResult<Envelope> {
    let mut envelope = Envelope::EMPTY;
    for vertex in solid.topology().vertices() {
        envelope.include(vertex.position);
    }
    // `|a cos t + b sin t| <= sqrt(a^2 + b^2)` per axis.
    let reach = |a: Vec3, b: Vec3| (a * a + b * b).powf(0.5);
    for curve in solid.curves3() {
        match curve {
            // Bounded by its edges' vertices.
            Curve3::Line(_) => {}
            Curve3::Circle(circle) => {
                let r = reach(circle.frame.x, circle.frame.y) * circle.radius;
                envelope.include(circle.frame.origin - r);
                envelope.include(circle.frame.origin + r);
            }
            Curve3::Ellipse(ellipse) => {
                let r = reach(
                    ellipse.frame.x * ellipse.semi_axis_x,
                    ellipse.frame.y * ellipse.semi_axis_y,
                );
                envelope.include(ellipse.frame.origin - r);
                envelope.include(ellipse.frame.origin + r);
            }
            Curve3::Polyline(polyline) => {
                for &point in &polyline.points {
                    envelope.include(point);
                }
            }
            // In the hull of its control points, its weights positive.
            Curve3::BSpline(spline)
                if spline
                    .weights
                    .as_ref()
                    .is_none_or(|weights| weights.iter().all(|&w| w > 0.0)) =>
            {
                for &point in &spline.control_points {
                    envelope.include(point);
                }
            }
            _ => {
                return Err(unsupported(
                    "exact half-space clip of a subject with an edge it cannot bound",
                ))
            }
        }
    }
    if !(envelope.lo.is_finite() && envelope.hi.is_finite()) {
        return Err(GeomError::InvalidInput(
            "exact half-space clip of a subject with no finite extent".to_owned(),
        ));
    }
    Ok(envelope)
}

/// A closed polygon as a contour of line segments.
fn polygon(points: &[Point2]) -> Profile {
    let segments = (0..points.len())
        .map(|i| {
            let (a, b) = (points[i], points[(i + 1) % points.len()]);
            ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: a,
                    direction: b - a,
                }),
                domain: Interval::new(0.0, 1.0),
                same_sense: true,
            }
        })
        .collect();
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

/// A placement's refusal, named.
fn transform_error(error: TransformError) -> GeomError {
    match error {
        TransformError::NotRigid => {
            unsupported("exact half-space under a scaled or sheared transform")
        }
        TransformError::Unsupported(what) => unsupported(what),
        other => GeomError::InvalidInput(format!("half-space transform: {other}")),
    }
}

/// Whether a boundary map is the identity to the rounding of its entries.
fn is_identity(map: &Transform2) -> bool {
    let floor = 8.0 * Scalar::EPSILON;
    (map.matrix2.x_axis - Vec2::X).abs().max_element() <= floor
        && (map.matrix2.y_axis - Vec2::Y).abs().max_element() <= floor
}

fn emptied() -> GeomError {
    GeomError::Degenerate("exact half-space clip leaves nothing of the subject".to_owned())
}
