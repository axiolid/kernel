//! Extrude an arc-capable cross-section into an exact B-rep (ADR 0050).
//!
//! # Why this is a second extruder
//!
//! `extrude_polygon_rings` builds a planar wall per edge. An arc edge
//! sweeps a CYLINDRICAL wall instead, so the surface, the 3D curve and the
//! parameter-space curves all differ. Rather than thread a branch through
//! the polygon extruder, the arc case is built here and the polygon path
//! is left untouched.
//!
//! # Exactness
//!
//! An arc edge becomes a `Circle3` bound to the swept `Cylinder` surface:
//! no sampling, no tessellation. The bulge is converted to centre, radius
//! and sweep by closed form, verified against the defining circle before
//! any B-rep is built.

use axiolid_brep::{ExactBRep, ExactBRepBuilder, FaceName, SweptFace};
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Frame3, Interval, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    Circle2, Circle3, Curve2, Curve3, Ellipse2, Ellipse3, Line2, Line3, Sinusoid2,
};
use axiolid_nurbs::exact_surface_intersection;
use axiolid_overlay::ArcRing;
use axiolid_surface::{Cylinder, Plane, Surface};
use axiolid_topology::{
    Edge, EdgeId, EdgeUse, Face, FaceBound, FaceId, LoopId, Orientation, Vertex, VertexId,
};

use crate::extrude_exact::{
    add_line_edge, add_loop, add_single_bound_face, finish_closed, identity_frame3, reserve,
};

/// The circle an arc edge lies on.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ArcGeometry {
    /// Circle centre in the cross-section plane.
    pub(crate) centre: Point2,
    /// Circle radius.
    pub(crate) radius: Scalar,
    /// Signed sweep, positive counter-clockwise.
    pub(crate) sweep: Scalar,
}

/// Recover the circle of an arc edge from its bulge.
///
/// The bulge is `tan(theta / 4)` for included angle `theta`, so
/// `theta = 4 atan(bulge)` carries both magnitude and direction. The
/// centre sits on the chord's left normal at `R cos(theta / 2)`, which
/// goes negative for a major arc and moves the centre to the other
/// side -- exactly the wanted behaviour, so no case split is needed.
pub(crate) fn arc_geometry(from: Point2, to: Point2, bulge: Scalar) -> GeomResult<ArcGeometry> {
    let chord = to - from;
    let length = chord.length();
    if length == 0.0 || !length.is_finite() || bulge == 0.0 {
        return Err(GeomError::Degenerate(
            "arc edge has no chord or no bulge".to_owned(),
        ));
    }
    let sweep = 4.0 * bulge.atan();
    let half = 0.5 * sweep;
    let sine = half.abs().sin();
    if sine == 0.0 {
        return Err(GeomError::Degenerate(
            "arc edge sweeps a full turn".to_owned(),
        ));
    }
    let radius = length / (2.0 * sine);
    let unit = chord / length;
    let normal = Vec2::new(-unit.y, unit.x);
    // The apothem must carry the sweep's SIGN. `cos` is even, so using it
    // unsigned puts +bulge and -bulge on the same centre: every arc would
    // bulge the same way regardless of its stated direction. A positive
    // sweep turns left, so its centre sits left of the chord.
    let centre = from + chord * 0.5 + normal * (radius * half.cos()) * sweep.signum();
    if !centre.is_finite() || !radius.is_finite() {
        return Err(GeomError::Degenerate(
            "arc edge produced a non-finite circle".to_owned(),
        ));
    }
    Ok(ArcGeometry {
        centre,
        radius,
        sweep,
    })
}

/// Height as an affine function of plan position:
/// `z = height + gradient . (x, y)`.
///
/// A flat level (zero gradient) is an ordinary cap plane. A sloped level is
/// a plane cut across the extrusion: over an arc edge it meets the swept
/// cylinder in an ellipse, whose pcurve on that cylinder is a
/// [`Sinusoid2`] (ADR 0071).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Level {
    /// Height at the plan origin.
    pub(crate) height: Scalar,
    /// Rise per unit `x` and per unit `y`.
    pub(crate) gradient: Vec2,
}

impl Level {
    /// A horizontal level at `height`.
    pub(crate) fn flat(height: Scalar) -> Self {
        Self {
            height,
            gradient: Vec2::ZERO,
        }
    }

    /// Height of the level above plan point `p`.
    pub(crate) fn at(&self, p: Point2) -> Scalar {
        self.height + self.gradient.dot(p)
    }

    pub(crate) fn is_flat(&self) -> bool {
        self.gradient == Vec2::ZERO
    }

    /// Orthonormal frame of a sloped level, origin above the plan origin.
    ///
    /// `x` climbs along plan `x`, so a cap pcurve's first coordinate still
    /// grows with plan `x`; `z` is the upward normal.
    pub(crate) fn sloped_frame(&self) -> GeomResult<Frame3> {
        let normal = Vec3::new(-self.gradient.x, -self.gradient.y, 1.0).normalize();
        let x = Vec3::new(1.0, 0.0, self.gradient.x).normalize();
        let y = normal.cross(x);
        let frame = Frame3 {
            origin: Vec3::new(0.0, 0.0, self.height),
            x,
            y,
            z: normal,
        };
        if !(frame.origin.is_finite() && x.is_finite() && y.is_finite() && normal.is_finite()) {
            return Err(GeomError::Degenerate(
                "sloped cap produced a non-finite plane".to_owned(),
            ));
        }
        Ok(frame)
    }
}

/// The two caps of an arc extrusion.
///
/// `shear` moves the top horizontally (an oblique extrusion) and is only
/// meaningful with two flat levels; the sloped path requires it to be
/// zero. `named` says which caps are the solid's own start and end caps:
/// a cap cut by another operand is not, and is left unnamed rather than
/// mislabelled.
#[derive(Debug, Clone, Copy)]
struct Span {
    bottom: Level,
    top: Level,
    shear: Vec2,
    named: (bool, bool),
}

impl Span {
    fn level(&self, is_top: bool) -> Level {
        if is_top {
            self.top
        } else {
            self.bottom
        }
    }
}

/// Per-edge topology of an extruded arc ring.
struct RingTopology {
    bottom_edges: Vec<EdgeId>,
    top_edges: Vec<EdgeId>,
    vertical_edges: Vec<EdgeId>,
    bottom_points: Vec<Point3>,
    top_points: Vec<Point3>,
    /// The 3D support and native span of each bottom edge, kept so a sloped
    /// cap can derive its pcurves from the edges it actually uses.
    bottom_curves: Vec<(Curve3, Interval)>,
    top_curves: Vec<(Curve3, Interval)>,
}

/// [`extrude_arc_rings`] from the plane `z = base` instead of `z = 0`.
///
/// Every absolute position (vertices, circle centres, cap planes, wall
/// frames) derives from the bottom points and the cap levels, so lifting
/// those is the whole change; surface parameters are relative and stay as
/// they are.
pub(crate) fn extrude_arc_rings_from(
    rings: &[ArcRing],
    base: Scalar,
    offset: Vec3,
) -> GeomResult<ExactBRep> {
    build_arc_rings(rings, flat_span(base, offset))
}

/// Extrude an arc-capable section with holes along `offset`.
///
/// Ring 0 is the outer boundary; any further rings are through-holes. Each
/// ring contributes its own walls and its own cap LOOP, and the two cap faces
/// carry one bound per ring -- that is what makes a hole a hole rather than a
/// second disconnected outline.
///
/// Hole rings must be wound CLOCKWISE, matching the polygon path: a ring's
/// wall normals follow its winding, so a counter-clockwise hole would face
/// its walls outward and produce a solid that is inside-out along the
/// passage.
pub(crate) fn extrude_arc_rings(rings: &[ArcRing], offset: Vec3) -> GeomResult<ExactBRep> {
    build_arc_rings(rings, flat_span(0.0, offset))
}

/// Extrude an arc-capable section straight up between two levels, either of
/// which may be sloped.
///
/// The caller guarantees `top.at(p) > bottom.at(p)` over the whole section;
/// a level that crosses the other inside the section is not a prism with
/// two caps and must be refused before this is called. `named` marks which
/// caps are the solid's own start and end caps.
pub(crate) fn extrude_arc_rings_between(
    rings: &[ArcRing],
    bottom: Level,
    top: Level,
    named: (bool, bool),
) -> GeomResult<ExactBRep> {
    build_arc_rings(
        rings,
        Span {
            bottom,
            top,
            shear: Vec2::ZERO,
            named,
        },
    )
}

fn flat_span(base: Scalar, offset: Vec3) -> Span {
    Span {
        bottom: Level::flat(base),
        top: Level::flat(base + offset.z),
        shear: Vec2::new(offset.x, offset.y),
        named: (true, true),
    }
}

fn build_arc_rings(rings: &[ArcRing], span: Span) -> GeomResult<ExactBRep> {
    if rings.is_empty() {
        return Err(GeomError::Degenerate(
            "arc extrusion needs at least one ring".to_owned(),
        ));
    }
    for ring in rings {
        if ring.vertices.len() < 2 {
            return Err(GeomError::Degenerate(
                "arc ring needs at least two vertices".to_owned(),
            ));
        }
    }
    if span.shear != Vec2::ZERO && !(span.bottom.is_flat() && span.top.is_flat()) {
        return Err(GeomError::Degenerate(
            "an oblique arc extrusion cannot have a sloped cap".to_owned(),
        ));
    }
    let total: usize = rings.iter().map(|ring| ring.vertices.len()).sum();
    let ring_count = rings.len();

    let mut builder = ExactBRepBuilder::default();
    reserve(
        &mut builder,
        total * 2,
        total * 3,
        2 * ring_count + total,
        2 + total,
        total * 3,
        total * 8,
        2 + total,
    )?;

    let topologies = rings
        .iter()
        .map(|ring| add_arc_ring(&mut builder, ring, &span))
        .collect::<GeomResult<Vec<_>>>()?;

    let bottom_frame = cap_frame(span.bottom, Vec2::ZERO)?;
    let top_frame = cap_frame(span.top, span.shear)?;
    let bottom_surface = builder.add_surface(Surface::Plane(Plane {
        frame: bottom_frame,
    }));
    let top_surface = builder.add_surface(Surface::Plane(Plane { frame: top_frame }));

    let mut bottom_bounds = Vec::with_capacity(ring_count);
    let mut top_bounds = Vec::with_capacity(ring_count);
    for (index, (ring, topology)) in rings.iter().zip(&topologies).enumerate() {
        let bottom_loop = if span.bottom.is_flat() {
            arc_cap_loop(&mut builder, ring, topology, false)
        } else {
            sloped_cap_loop(&mut builder, topology, false, &bottom_frame)?
        };
        let top_loop = if span.top.is_flat() {
            arc_cap_loop(&mut builder, ring, topology, true)
        } else {
            sloped_cap_loop(&mut builder, topology, true, &top_frame)?
        };
        bottom_bounds.push(FaceBound {
            loop_id: bottom_loop,
            orientation: Orientation::Forward,
            outer: index == 0,
        });
        top_bounds.push(FaceBound {
            loop_id: top_loop,
            orientation: Orientation::Forward,
            outer: index == 0,
        });
    }

    let bottom_face = builder.topology_mut().add_face(Face {
        surface: Some(bottom_surface),
        bounds: bottom_bounds,
        orientation: Orientation::Reversed,
    });
    let top_face = builder.topology_mut().add_face(Face {
        surface: Some(top_surface),
        bounds: top_bounds,
        orientation: Orientation::Forward,
    });
    if span.named.0 {
        builder.set_face_name(bottom_face, FaceName::swept(SweptFace::StartCap));
    }
    if span.named.1 {
        builder.set_face_name(top_face, FaceName::swept(SweptFace::EndCap));
    }

    let mut faces = vec![bottom_face, top_face];
    // Wall ordinals run across all rings so each side face keeps a distinct
    // name; restarting per ring would collide the outer and hole walls.
    let mut ordinal: u32 = 0;
    for (ring, topology) in rings.iter().zip(&topologies) {
        for index in 0..ring.vertices.len() {
            let bulge = ring.vertices[index].bulge;
            let face = if bulge == 0.0 {
                add_planar_wall(&mut builder, ring, topology, index, &span)?
            } else {
                add_cylindrical_wall(&mut builder, ring, topology, index, &span)?
            };
            builder.set_face_name(face, FaceName::swept(SweptFace::Side(ordinal)));
            ordinal = ordinal.checked_add(1).ok_or_else(|| {
                GeomError::Degenerate("profile edge count exceeds u32 capacity".to_owned())
            })?;
            faces.push(face);
        }
    }
    finish_closed(builder, faces)
}

/// Frame of a cap plane. A flat cap keeps the identity axes (so its pcurves
/// are plain plan coordinates); a sloped cap uses [`Level::sloped_frame`].
fn cap_frame(level: Level, shear: Vec2) -> GeomResult<Frame3> {
    if level.is_flat() {
        Ok(identity_frame3(Vec3::new(shear.x, shear.y, level.height)))
    } else {
        level.sloped_frame()
    }
}

/// Build the vertices and the three edge families of an arc ring.
fn add_arc_ring(
    builder: &mut ExactBRepBuilder,
    ring: &ArcRing,
    span: &Span,
) -> GeomResult<RingTopology> {
    let count = ring.vertices.len();
    let bottom_points: Vec<Point3> = ring
        .vertices
        .iter()
        .map(|vertex| Point3::new(vertex.point.x, vertex.point.y, span.bottom.at(vertex.point)))
        .collect();
    let top_points: Vec<Point3> = ring
        .vertices
        .iter()
        .map(|vertex| {
            Point3::new(
                vertex.point.x + span.shear.x,
                vertex.point.y + span.shear.y,
                span.top.at(vertex.point),
            )
        })
        .collect();
    let bottom_vertices: Vec<VertexId> = bottom_points
        .iter()
        .map(|position| {
            builder.topology_mut().add_vertex(Vertex {
                position: *position,
            })
        })
        .collect();
    let top_vertices: Vec<VertexId> = top_points
        .iter()
        .map(|position| {
            builder.topology_mut().add_vertex(Vertex {
                position: *position,
            })
        })
        .collect();

    let mut bottom_edges = Vec::with_capacity(count);
    let mut top_edges = Vec::with_capacity(count);
    let mut bottom_curves = Vec::with_capacity(count);
    let mut top_curves = Vec::with_capacity(count);
    for index in 0..count {
        let next = (index + 1) % count;
        let from = ring.vertices[index].point;
        let to = ring.vertices[next].point;
        let bulge = ring.vertices[index].bulge;
        // Iterate by flag, not by height: a zero-height offset would make
        // the two levels equal and put every edge in the bottom family. The
        // caller refuses zero height, but a silent mis-binding here would be
        // far harder to see than a refusal.
        for (is_top, vertices) in [(false, &bottom_vertices), (true, &top_vertices)] {
            let level = span.level(is_top);
            let (curve, interval) = if bulge == 0.0 {
                // A straight edge at either level is the chord between its
                // two vertices at that level, flat or sloped. The shear is
                // deliberately not applied to the edge support: the flat
                // path has always built edges over the plan section.
                let start = Point3::new(from.x, from.y, level.at(from));
                let end = Point3::new(to.x, to.y, level.at(to));
                (
                    Curve3::Line(Line3 {
                        origin: start,
                        direction: end - start,
                    }),
                    Interval::UNIT,
                )
            } else {
                let arc = arc_geometry(from, to, bulge)?;
                if level.is_flat() {
                    let start = Point3::new(from.x, from.y, level.height);
                    (
                        Curve3::Circle(circle_of(&arc, level.height, start)?),
                        Interval::new(0.0, arc.sweep),
                    )
                } else {
                    sloped_arc_edge(&arc, from, level, span.bottom.height)?
                }
            };
            let curve_id = builder.add_curve3(curve.clone());
            let edge = builder.topology_mut().add_edge(Edge {
                start: vertices[index],
                end: vertices[next],
                curve: Some(curve_id),
            });
            builder.set_edge_interval(edge, interval);
            if is_top {
                top_edges.push(edge);
                top_curves.push((curve, interval));
            } else {
                bottom_edges.push(edge);
                bottom_curves.push((curve, interval));
            }
        }
    }

    let mut vertical_edges = Vec::with_capacity(count);
    for index in 0..count {
        vertical_edges.push(add_line_edge(
            builder,
            bottom_vertices[index],
            top_vertices[index],
            bottom_points[index],
            top_points[index] - bottom_points[index],
        ));
    }

    Ok(RingTopology {
        bottom_edges,
        top_edges,
        vertical_edges,
        bottom_points,
        top_points,
        bottom_curves,
        top_curves,
    })
}

/// The edge where a sloped level cuts the cylinder swept by an arc edge.
///
/// The ellipse comes from the exact cylinder/plane intersection (#119), not
/// from a formula repeated here. Its native parameter is linear in the
/// cylinder angle -- the minor semi-axis equals the radius, and the major
/// axis's plan projection has the radius as its length -- so the edge spans
/// exactly `sweep` of it, starting at the parameter of the arc's start
/// point and running the way the arc does.
pub(crate) fn sloped_arc_edge(
    arc: &ArcGeometry,
    from: Point2,
    level: Level,
    anchor: Scalar,
) -> GeomResult<(Curve3, Interval)> {
    let circle = circle_of(arc, anchor, Point3::new(from.x, from.y, anchor))?;
    let cylinder = Surface::Cylinder(Cylinder {
        frame: circle.frame,
        radius: arc.radius,
    });
    let plane = Surface::Plane(Plane {
        frame: level.sloped_frame()?,
    });
    let cut = exact_surface_intersection(&cylinder, &plane).map_err(|refusal| {
        GeomError::Degenerate(format!("sloped cap cannot cut an arc wall: {refusal:?}"))
    })?;
    let ellipse = match cut.branches.as_slice() {
        [Curve3::Ellipse(ellipse)] => *ellipse,
        [Curve3::Circle(circle)] => Ellipse3 {
            frame: circle.frame,
            semi_axis_x: circle.radius,
            semi_axis_y: circle.radius,
        },
        _ => {
            return Err(GeomError::Degenerate(
                "a sloped cap cut an arc wall in something other than one ellipse".to_owned(),
            ))
        }
    };
    let start = Point3::new(from.x, from.y, level.at(from)) - ellipse.frame.origin;
    let phase = (start.dot(ellipse.frame.y) / ellipse.semi_axis_y)
        .atan2(start.dot(ellipse.frame.x) / ellipse.semi_axis_x);
    // Which way the ellipse parameter runs relative to the arc: compare its
    // tangent at the start with the arc's own counter-clockwise direction.
    // Today both frames point up, so this is always +1 and the sweep's sign
    // alone orients the edge. The check stays as a guard: if the
    // intersection ever returns a downward frame, edges still run from
    // their start vertex to their end vertex (pinned by `clip_arc_prism`).
    let (sin, cos) = phase.sin_cos();
    let tangent = ellipse.frame.x * (-ellipse.semi_axis_x * sin)
        + ellipse.frame.y * (ellipse.semi_axis_y * cos);
    let turn = if tangent.dot(circle.frame.y) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    if !phase.is_finite() {
        return Err(GeomError::Degenerate(
            "sloped cap edge has no finite start parameter".to_owned(),
        ));
    }
    Ok((
        Curve3::Ellipse(ellipse),
        Interval::new(phase, phase + turn * arc.sweep),
    ))
}

/// The `Circle3` for an arc at height `level`, oriented so its frame
/// x-axis points at the arc start.
///
/// Anchoring x at the start vertex makes the surface parameter u the angle
/// measured from that vertex, so an edge interval of `0..sweep` is
/// literally the arc and needs no offset term.
pub(crate) fn circle_of(arc: &ArcGeometry, level: Scalar, start: Point3) -> GeomResult<Circle3> {
    let centre = Point3::new(arc.centre.x, arc.centre.y, level);
    let radial = start - centre;
    let length = radial.length();
    if length == 0.0 || !length.is_finite() {
        return Err(GeomError::Degenerate(
            "arc edge start lies on its own centre".to_owned(),
        ));
    }
    let x = radial / length;
    let z = Vec3::Z;
    let y = z.cross(x);
    Ok(Circle3 {
        frame: Frame3 {
            origin: centre,
            x,
            y,
            z,
        },
        radius: arc.radius,
    })
}

/// A cylindrical wall swept by an arc edge.
///
/// The surface is parameterised by angle `u` from the arc start and height
/// `v` above the bottom level's anchor height. A flat rim is the straight
/// pcurve `v = const`; a sloped rim is the [`Sinusoid2`] the plane traces
/// across the cylinder, so the trim stays exact rather than approximated.
fn add_cylindrical_wall(
    builder: &mut ExactBRepBuilder,
    ring: &ArcRing,
    topology: &RingTopology,
    index: usize,
    span: &Span,
) -> GeomResult<FaceId> {
    let count = ring.vertices.len();
    let next = (index + 1) % count;
    let from = ring.vertices[index].point;
    let to = ring.vertices[next].point;
    let arc = arc_geometry(from, to, ring.vertices[index].bulge)?;
    let anchor = span.bottom.height;
    let circle = circle_of(&arc, anchor, Point3::new(from.x, from.y, anchor))?;
    let surface = builder.add_surface(Surface::Cylinder(Cylinder {
        frame: circle.frame,
        radius: arc.radius,
    }));

    // A rim pcurve at `level`, with its native span.
    let rim = |level: Level| -> (Curve2, Interval) {
        if level.is_flat() {
            (
                Curve2::Line(Line2 {
                    origin: Vec2::new(0.0, level.height - anchor),
                    direction: Vec2::new(arc.sweep, 0.0),
                }),
                Interval::UNIT,
            )
        } else {
            // Above the point at angle u the level sits at
            //   h + g.c + r (g.x_h cos u + g.y_h sin u),
            // with c the centre and x_h, y_h the frame's plan axes.
            let x = Vec2::new(circle.frame.x.x, circle.frame.x.y);
            let y = Vec2::new(circle.frame.y.x, circle.frame.y.y);
            (
                Curve2::Sinusoid(Sinusoid2 {
                    mean: level.at(arc.centre) - anchor,
                    cosine: arc.radius * level.gradient.dot(x),
                    sine: arc.radius * level.gradient.dot(y),
                }),
                Interval::new(0.0, arc.sweep),
            )
        }
    };
    // A vertical pcurve at angle `u` between the two levels above `p`.
    let rise = |u: Scalar, p: Point2| -> (Curve2, Interval) {
        let low = span.bottom.at(p) - anchor;
        let high = span.top.at(p) - anchor;
        (
            Curve2::Line(Line2 {
                origin: Vec2::new(u, low),
                direction: Vec2::new(0.0, high - low),
            }),
            Interval::UNIT,
        )
    };
    let uses = [
        (
            topology.bottom_edges[index],
            Orientation::Forward,
            rim(span.bottom),
        ),
        (
            topology.vertical_edges[next],
            Orientation::Forward,
            rise(arc.sweep, to),
        ),
        (
            topology.top_edges[index],
            Orientation::Reversed,
            rim(span.top),
        ),
        (
            topology.vertical_edges[index],
            Orientation::Reversed,
            rise(0.0, from),
        ),
    ];
    Ok(wall_face(builder, surface, uses))
}

/// A planar wall swept by a straight edge.
///
/// Every pcurve is the straight line between its two corners, projected
/// into the wall plane. That covers flat and sloped caps with one rule: a
/// sloped cap only moves the corners along the wall, never off it.
fn add_planar_wall(
    builder: &mut ExactBRepBuilder,
    ring: &ArcRing,
    topology: &RingTopology,
    index: usize,
    span: &Span,
) -> GeomResult<FaceId> {
    let count = ring.vertices.len();
    let next = (index + 1) % count;
    let start = topology.bottom_points[index];
    let from = ring.vertices[index].point;
    let to = ring.vertices[next].point;
    let plan = Vec3::new(to.x - from.x, to.y - from.y, 0.0);
    let length = plan.length();
    if length == 0.0 || !length.is_finite() {
        return Err(GeomError::Degenerate(
            "straight edge has zero length".to_owned(),
        ));
    }
    let x = plan / length;
    // The wall contains the sweep direction. With a shear that is the flat
    // offset itself; a sloped span has no shear and sweeps straight up.
    let sweep = if span.shear == Vec2::ZERO {
        Vec3::Z
    } else {
        Vec3::new(
            span.shear.x,
            span.shear.y,
            span.top.height - span.bottom.height,
        )
    };
    let z = x.cross(sweep).normalize();
    if !z.is_finite() {
        return Err(GeomError::Degenerate(
            "straight wall produced a degenerate plane".to_owned(),
        ));
    }
    let y = z.cross(x);
    let surface = builder.add_surface(Surface::Plane(Plane {
        frame: Frame3 {
            origin: start,
            x,
            y,
            z,
        },
    }));
    let uv = |p: Point3| -> Vec2 {
        let d = p - start;
        Vec2::new(d.dot(x), d.dot(y))
    };
    let line = |a: Point3, b: Point3| -> (Curve2, Interval) {
        let (a, b) = (uv(a), uv(b));
        (
            Curve2::Line(Line2 {
                origin: a,
                direction: b - a,
            }),
            Interval::UNIT,
        )
    };
    let (b0, b1) = (topology.bottom_points[index], topology.bottom_points[next]);
    let (t0, t1) = (topology.top_points[index], topology.top_points[next]);
    let uses = [
        (
            topology.bottom_edges[index],
            Orientation::Forward,
            line(b0, b1),
        ),
        (
            topology.vertical_edges[next],
            Orientation::Forward,
            line(b1, t1),
        ),
        (
            topology.top_edges[index],
            Orientation::Reversed,
            line(t0, t1),
        ),
        (
            topology.vertical_edges[index],
            Orientation::Reversed,
            line(b0, t0),
        ),
    ];
    Ok(wall_face(builder, surface, uses))
}

/// Assemble a four-sided wall face from its edge uses.
///
/// Each pcurve is given in its EDGE's direction with its native span; a
/// reversed use walks that span backwards, which is how the loop's
/// traversal order is stated without building a second, mirrored curve.
fn wall_face(
    builder: &mut ExactBRepBuilder,
    surface: axiolid_brep::SurfaceId,
    edge_uses: [(EdgeId, Orientation, (Curve2, Interval)); 4],
) -> FaceId {
    let mut uses = Vec::with_capacity(4);
    let mut intervals = Vec::with_capacity(4);
    for (edge, orientation, (pcurve, native)) in edge_uses {
        let pcurve = builder.add_curve2(pcurve);
        uses.push(EdgeUse {
            edge,
            orientation,
            pcurve: Some(pcurve),
        });
        intervals.push(match orientation {
            Orientation::Forward => native,
            Orientation::Reversed => Interval::new(native.end, native.start),
        });
    }
    let loop_id = add_loop(builder, uses, intervals);
    add_single_bound_face(builder, surface, loop_id, Orientation::Forward)
}

/// The cap loop of a ring on a SLOPED cap.
///
/// Each pcurve is the cap edge's own 3D support expressed in the cap
/// plane's frame: a line stays a line, and an ellipse lying in the plane
/// stays the same ellipse with its frame projected. The edge's native span
/// carries over unchanged, because projecting into the plane it lies in
/// does not reparameterise it.
fn sloped_cap_loop(
    builder: &mut ExactBRepBuilder,
    topology: &RingTopology,
    top: bool,
    frame: &Frame3,
) -> GeomResult<LoopId> {
    let (edges, curves) = if top {
        (&topology.top_edges, &topology.top_curves)
    } else {
        (&topology.bottom_edges, &topology.bottom_curves)
    };
    let point = |p: Point3| {
        let d = p - frame.origin;
        Vec2::new(d.dot(frame.x), d.dot(frame.y))
    };
    let vector = |v: Vec3| Vec2::new(v.dot(frame.x), v.dot(frame.y));
    let mut uses = Vec::with_capacity(edges.len());
    let mut intervals = Vec::with_capacity(edges.len());
    for (edge, (curve, interval)) in edges.iter().zip(curves) {
        let pcurve = match curve {
            Curve3::Line(line) => Curve2::Line(Line2 {
                origin: point(line.origin),
                direction: vector(line.direction),
            }),
            Curve3::Ellipse(ellipse) => Curve2::Ellipse(Ellipse2 {
                frame: Frame2 {
                    origin: point(ellipse.frame.origin),
                    x: vector(ellipse.frame.x),
                    y: vector(ellipse.frame.y),
                },
                semi_axis_x: ellipse.semi_axis_x,
                semi_axis_y: ellipse.semi_axis_y,
            }),
            _ => {
                return Err(GeomError::Degenerate(
                    "a sloped cap edge is neither a line nor an ellipse".to_owned(),
                ))
            }
        };
        let pcurve = builder.add_curve2(pcurve);
        uses.push(EdgeUse {
            edge: *edge,
            orientation: Orientation::Forward,
            pcurve: Some(pcurve),
        });
        intervals.push(*interval);
    }
    Ok(add_loop(builder, uses, intervals))
}

/// The cap loop of an arc ring, in the cross-section plane.
///
/// Each pcurve is the cross-section curve itself: a `Line2` for a straight
/// edge, a circle for an arc. The cap surface is the z-plane, so the
/// parameter space IS the cross-section plane and no mapping is needed.
fn arc_cap_loop(
    builder: &mut ExactBRepBuilder,
    ring: &ArcRing,
    topology: &RingTopology,
    top: bool,
) -> LoopId {
    let count = ring.vertices.len();
    let edges = if top {
        &topology.top_edges
    } else {
        &topology.bottom_edges
    };
    let mut uses = Vec::with_capacity(count);
    let mut intervals = Vec::with_capacity(count);
    for (index, edge) in edges.iter().enumerate().take(count) {
        let next = (index + 1) % count;
        let from = ring.vertices[index].point;
        let to = ring.vertices[next].point;
        let bulge = ring.vertices[index].bulge;
        // The cap pcurve must be the cross-section curve itself. Using a
        // chord for an arc edge would make the cap boundary disagree with
        // the wall boundary along the same edge -- a silent inconsistency
        // that no area check would reveal.
        let (curve, interval) = if bulge == 0.0 {
            (
                Curve2::Line(Line2 {
                    origin: Vec2::new(from.x, from.y),
                    direction: Vec2::new(to.x - from.x, to.y - from.y),
                }),
                Interval::UNIT,
            )
        } else {
            match arc_geometry(from, to, bulge) {
                Ok(arc) => (
                    Curve2::Circle(Circle2 {
                        frame: circle2_frame(&arc, from),
                        radius: arc.radius,
                    }),
                    Interval::new(0.0, arc.sweep),
                ),
                // A malformed arc was already refused when its 3D edge was
                // built, so this branch is unreachable in practice; falling
                // back to the chord keeps the loop total rather than
                // panicking on an impossible state.
                Err(_) => (
                    Curve2::Line(Line2 {
                        origin: Vec2::new(from.x, from.y),
                        direction: Vec2::new(to.x - from.x, to.y - from.y),
                    }),
                    Interval::UNIT,
                ),
            }
        };
        let pcurve = builder.add_curve2(curve);
        uses.push(EdgeUse {
            edge: *edge,
            orientation: Orientation::Forward,
            pcurve: Some(pcurve),
        });
        intervals.push(interval);
    }
    add_loop(builder, uses, intervals)
}

/// A 2D frame whose x-axis points from the arc centre at its start.
pub(crate) fn circle2_frame(arc: &ArcGeometry, start: Point2) -> Frame2 {
    let radial = start - arc.centre;
    let length = radial.length();
    let x = if length == 0.0 {
        Vec2::X
    } else {
        radial / length
    };
    Frame2 {
        origin: Vec2::new(arc.centre.x, arc.centre.y),
        x,
        y: Vec2::new(-x.y, x.x),
    }
}
