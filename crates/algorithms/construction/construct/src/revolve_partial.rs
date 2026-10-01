//! Exact partial-turn revolution (#172).
//!
//! A full turn closes on itself, so its walls need seams and its caps are
//! annuli between whole circles (`revolve_contour`). A partial turn through
//! `0 < |angle| < 2 pi` is a different topology: every profile station sweeps
//! a circular ARC, every profile segment appears twice -- once in the start
//! wall and once in the end wall -- and the two planar end caps are the
//! profile itself, at the start and end angles. Nothing closes round the
//! axis, so no seam edge is needed.
//!
//! # Topology
//!
//! For each station `i` (radius `r`, height `h` in the profile plane) there
//! is a start vertex `S_i` and an end vertex `E_i`, joined by the arc `A_i`.
//! Each segment `i -> i+1` yields a start-cap edge `P_i` (`S_i -> S_{i+1}`),
//! an end-cap edge `Q_i` (`E_i -> E_{i+1}`) and one wall face bounded by
//! `A_i`, `Q_i`, `A_{i+1}` reversed and `P_i` reversed: in the wall's
//! `(angle, v)` parameters the trim is the rectangle `[u0, u1] x [v_i,
//! v_{i+1}]`, the same rectangle the full-turn seam loop spans, only
//! narrower.
//!
//! # Touching the axis
//!
//! A station ON the axis sweeps no arc: `S_i = E_i` is one vertex and the
//! adjacent walls lose that edge (a cone wall closes at its apex, a planar
//! wall becomes a sector). A straight segment lying ON the axis sweeps no
//! wall at all; its start and end copies are the same edge, shared by the
//! two caps. A section CROSSING the axis would sweep through itself and is
//! refused, as is an arc whose circle reaches the axis (a spindle or horn
//! torus), exactly as the full turn refuses it.
//!
//! # Orientation
//!
//! Every ring is walked so the material is on its left in `(r, h)`: the
//! outer ring anticlockwise, each hole clockwise. Walls then follow the
//! full-turn convention (face forward, loop anticlockwise in `(u, v)` for a
//! rising segment). The start cap's plane normal `radial(u0) x axis` points
//! against the sweep, out of the material, so it is used forward; the end
//! cap's points into the material, so it is used reversed.

use axiolid_brep::{ExactBRep, ExactBRepBuilder, SurfaceId};
use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Frame3, Interval, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Circle3, Curve2, Curve3, Line2, Line3};
use axiolid_overlay::{ArcRing, ArcVertex};
use axiolid_surface::{Cone, Cylinder, Plane, Surface, Torus};
use axiolid_topology::{
    Edge, EdgeId, EdgeUse, Face, FaceBound, FaceId, Loop, Orientation, Shell, Solid, Vertex,
    VertexId,
};

use crate::contour_lower::arc_ring_signed_area;
use crate::revolve_contour::frame_at;
use crate::BACKEND_ID;

const TAU: Scalar = core::f64::consts::TAU;
const PI: Scalar = core::f64::consts::PI;

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BACKEND_ID,
        operation: Operation::Sweep,
        input,
    }
}

/// One section vertex in the axis frame: distance from the axis, height
/// along it, and the bulge of the segment leaving it.
#[derive(Debug, Clone, Copy)]
struct Station {
    radius: Scalar,
    height: Scalar,
    bulge: Scalar,
}

/// The topology one revolved ring contributes, kept for its two cap loops.
struct RingEdges {
    /// Start-cap edges `P_i` with their cap pcurves.
    start: Vec<(EdgeId, Curve2, Interval)>,
    /// End-cap edges `Q_i` with their cap pcurves.
    end: Vec<(EdgeId, Curve2, Interval)>,
}

/// Revolve one section -- an outer ring and its holes, in the profile plane
/// `z = 0` -- about the axis through `axis_origin` along world `+y`, through
/// `angle` radians by the right-hand rule, into an exact closed B-rep.
///
/// `0 < |angle| < 2 pi`. A negative angle sweeps the other way. The section
/// may lie on either side of the axis and may touch it along straight
/// segments or at vertices; it must not cross it. Each hole becomes a
/// tunnel running from the start cap to the end cap, not a void: a partial
/// turn does not close round it.
///
/// # Errors
///
/// [`GeomError::InvalidInput`] for a zero, non-finite or full-turn angle,
/// [`GeomError::Degenerate`] for a section without area or an angle whose
/// sweep is below the tolerance, and [`GeomError::UnsupportedInput`] for a
/// section crossing the axis or an arc whose circle reaches it.
pub fn revolve_section_partial(
    outer: &ArcRing,
    holes: &[ArcRing],
    axis_origin: Point3,
    angle: Scalar,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    if !angle.is_finite() || angle == 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "partial revolution angle must be finite and non-zero, got {angle}"
        )));
    }
    let epsilon = tolerance.linear();
    if angle.abs() >= TAU - epsilon {
        return Err(GeomError::InvalidInput(format!(
            "a partial revolution must turn less than a full turn, got {angle}"
        )));
    }
    if !axis_origin.is_finite() {
        return Err(GeomError::InvalidInput(
            "revolution axis origin must be finite".to_owned(),
        ));
    }

    // Which side of the axis the section lies on. A full turn makes the two
    // sides the same solid; a partial turn does not, so the side becomes the
    // base angle the sweep starts from rather than a mirror.
    let offsets = || {
        core::iter::once(outer)
            .chain(holes)
            .flat_map(|ring| ring.vertices.iter())
            .map(|vertex| vertex.point.x - axis_origin.x)
    };
    let min = offsets().fold(Scalar::INFINITY, Scalar::min);
    let max = offsets().fold(Scalar::NEG_INFINITY, Scalar::max);
    if !min.is_finite() || !max.is_finite() {
        return Err(GeomError::InvalidInput(
            "revolved section must have finite vertices".to_owned(),
        ));
    }
    if min < -epsilon && max > epsilon {
        return Err(unsupported(
            "exact revolution of a section crossing the axis",
        ));
    }
    let side = if max > epsilon { 1.0 } else { -1.0 };
    if max.abs().max(min.abs()) <= epsilon {
        return Err(GeomError::Degenerate(
            "revolved section lies on the axis".to_owned(),
        ));
    }
    let base = if side > 0.0 { 0.0 } else { PI };
    let (u0, u1) = (base + angle.min(0.0), base + angle.max(0.0));

    let outer_stations = stations(outer, axis_origin.x, side, true, epsilon)?;
    let hole_stations = holes
        .iter()
        .map(|hole| stations(hole, axis_origin.x, side, false, epsilon))
        .collect::<GeomResult<Vec<_>>>()?;

    let reach = outer_stations
        .iter()
        .map(|station| station.radius)
        .fold(0.0, Scalar::max);
    if reach * angle.abs() <= epsilon {
        return Err(GeomError::Degenerate(format!(
            "revolution through {angle} rad sweeps less than the tolerance"
        )));
    }

    let sweep = Sweep {
        axis_x: axis_origin.x,
        u0,
        u1,
        epsilon,
    };
    let mut builder = ExactBRepBuilder::default();
    let mut faces = Vec::new();
    let mut rings = Vec::with_capacity(1 + holes.len());
    for ring in core::iter::once(&outer_stations).chain(&hole_stations) {
        rings.push(sweep.ring(&mut builder, ring, &mut faces)?);
    }

    // The two caps: the profile itself, at the start and end angles.
    let start_plane = builder.add_surface(Surface::Plane(Plane {
        frame: sweep.cap_frame(u0),
    }));
    let end_plane = builder.add_surface(Surface::Plane(Plane {
        frame: sweep.cap_frame(u1),
    }));
    faces.push(cap(
        &mut builder,
        start_plane,
        rings.iter().map(|ring| ring.start.as_slice()),
        Orientation::Forward,
    ));
    faces.push(cap(
        &mut builder,
        end_plane,
        rings.iter().map(|ring| ring.end.as_slice()),
        Orientation::Reversed,
    ));

    let shell = builder.topology_mut().add_shell(Shell {
        faces: faces
            .into_iter()
            .map(|face| (face, Orientation::Forward))
            .collect(),
        closed: true,
    });
    builder.topology_mut().add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    let exact = builder
        .finish()
        .map_err(|error| GeomError::BackendContractViolation {
            backend: BACKEND_ID,
            detail: format!("partial revolution assembly failed: {error}"),
        })?;
    let health = axiolid_topology::audit_brep(exact.topology());
    if !health.is_closed_manifold() {
        return Err(GeomError::BackendContractViolation {
            backend: BACKEND_ID,
            detail: format!("partial revolution is not a closed manifold: {health:?}"),
        });
    }
    Ok(exact)
}

/// A ring's stations in the axis frame, walked with the material on the
/// left: anticlockwise in `(r, h)` for the outer ring, clockwise for a hole.
fn stations(
    ring: &ArcRing,
    axis_x: Scalar,
    side: Scalar,
    outer: bool,
    epsilon: Scalar,
) -> GeomResult<Vec<Station>> {
    // Measured from the axis on the section's own side; mirroring the far
    // side flips every arc's turning direction with it.
    let mut mapped: Vec<ArcVertex> = Vec::with_capacity(ring.vertices.len());
    for vertex in &ring.vertices {
        let radius = side * (vertex.point.x - axis_x);
        // A vertex within tolerance of the axis is ON it: it sweeps no arc.
        let radius = if radius <= epsilon { 0.0 } else { radius };
        let point = Point2::new(radius, vertex.point.y);
        let bulge = side * vertex.bulge;
        // A straight run of zero length states nothing.
        if let Some(last) = mapped.last() {
            if last.bulge == 0.0 && (last.point - point).length() <= epsilon {
                if let Some(last) = mapped.last_mut() {
                    last.bulge = bulge;
                }
                continue;
            }
        }
        mapped.push(ArcVertex { point, bulge });
    }
    while mapped.len() > 1 {
        let (first, last) = (mapped[0], mapped[mapped.len() - 1]);
        if last.bulge == 0.0 && (last.point - first.point).length() <= epsilon {
            mapped.pop();
        } else {
            break;
        }
    }
    let count = mapped.len();
    if count < 3 && mapped.iter().all(|vertex| vertex.bulge == 0.0) {
        return Err(GeomError::Degenerate(format!(
            "a revolved section ring needs at least three vertices, got {count}"
        )));
    }
    if count < 2 {
        return Err(GeomError::Degenerate(
            "a revolved section ring needs at least two vertices".to_owned(),
        ));
    }
    let area = arc_ring_signed_area(&ArcRing {
        vertices: mapped.clone(),
    });
    if area == 0.0 || !area.is_finite() {
        return Err(GeomError::Degenerate(
            "revolved section ring encloses no area".to_owned(),
        ));
    }
    let keep = (area > 0.0) == outer;
    Ok((0..count)
        .map(|index| {
            let source = if keep { index } else { count - 1 - index };
            let bulge = if keep {
                mapped[source].bulge
            } else {
                // Reversing the walk moves each bulge onto the segment that
                // now leaves this station, and flips its side.
                -mapped[(source + count - 1) % count].bulge
            };
            Station {
                radius: mapped[source].point.x,
                height: mapped[source].point.y,
                bulge,
            }
        })
        .collect())
}

/// The swept angle range and axis, shared by every face built.
struct Sweep {
    axis_x: Scalar,
    u0: Scalar,
    u1: Scalar,
    epsilon: Scalar,
}

impl Sweep {
    /// Unit direction from the axis at angle `u`, matching `frame_at`'s
    /// `(x, y) = (X, -Z)`: `u` increases by the right-hand rule about `+y`.
    fn radial(u: Scalar) -> Vec3 {
        let (sin, cos) = u.sin_cos();
        Vec3::new(cos, 0.0, -sin)
    }

    /// A profile point `(r, h)` placed at angle `u`.
    fn place(&self, u: Scalar, radius: Scalar, height: Scalar) -> Point3 {
        Point3::new(self.axis_x, height, 0.0) + Self::radial(u) * radius
    }

    /// The plane holding the profile at angle `u`: `x` radial, `y` along the
    /// axis, so a cap pcurve is just the profile's own `(r, h)`.
    fn cap_frame(&self, u: Scalar) -> Frame3 {
        let x = Self::radial(u);
        Frame3 {
            origin: Point3::new(self.axis_x, 0.0, 0.0),
            x,
            y: Vec3::Y,
            z: x.cross(Vec3::Y),
        }
    }

    /// Build one ring's vertices, arcs, cap edges and walls.
    fn ring(
        &self,
        builder: &mut ExactBRepBuilder,
        stations: &[Station],
        faces: &mut Vec<FaceId>,
    ) -> GeomResult<RingEdges> {
        let count = stations.len();
        let (u0, u1) = (self.u0, self.u1);

        // Vertices: one per end, or one shared on the axis.
        let mut starts: Vec<VertexId> = Vec::with_capacity(count);
        let mut ends: Vec<VertexId> = Vec::with_capacity(count);
        let mut arcs: Vec<Option<EdgeId>> = Vec::with_capacity(count);
        for station in stations {
            let start = builder.topology_mut().add_vertex(Vertex {
                position: self.place(u0, station.radius, station.height),
            });
            if station.radius == 0.0 {
                starts.push(start);
                ends.push(start);
                arcs.push(None);
                continue;
            }
            let end = builder.topology_mut().add_vertex(Vertex {
                position: self.place(u1, station.radius, station.height),
            });
            let curve = builder.add_curve3(Curve3::Circle(Circle3 {
                frame: frame_at(self.axis_x, station.height),
                radius: station.radius,
            }));
            let arc = builder.topology_mut().add_edge(Edge {
                start,
                end,
                curve: Some(curve),
            });
            builder.set_edge_interval(arc, Interval::new(u0, u1));
            starts.push(start);
            ends.push(end);
            arcs.push(Some(arc));
        }

        let mut edges = RingEdges {
            start: Vec::with_capacity(count),
            end: Vec::with_capacity(count),
        };
        for index in 0..count {
            let next = (index + 1) % count;
            let (here, there) = (stations[index], stations[next]);
            let segment = Segment::of(here, there)?;

            let start = self.profile_edge(builder, u0, &segment, starts[index], starts[next]);
            edges.start.push(start.clone());
            if here.bulge == 0.0 && here.radius == 0.0 && there.radius == 0.0 {
                // On the axis: no wall, and the end copy is the same edge.
                edges.end.push(start);
                continue;
            }
            let end = self.profile_edge(builder, u1, &segment, ends[index], ends[next]);
            edges.end.push(end.clone());

            let wall = self.wall(builder, here, there, &segment)?;
            let face = self.wall_face(
                builder,
                &wall,
                (arcs[index], arcs[next]),
                (start.0, end.0),
                (here, there),
            );
            faces.push(face);
        }
        Ok(edges)
    }

    /// A profile segment placed at angle `u`, with its pcurve on that cap.
    fn profile_edge(
        &self,
        builder: &mut ExactBRepBuilder,
        u: Scalar,
        segment: &Segment,
        from: VertexId,
        to: VertexId,
    ) -> (EdgeId, Curve2, Interval) {
        let (curve3, curve2, interval) = match *segment {
            Segment::Line { from: a, to: b } => {
                let origin = self.place(u, a.x, a.y);
                let curve3 = Curve3::Line(Line3 {
                    origin,
                    direction: self.place(u, b.x, b.y) - origin,
                });
                let curve2 = Curve2::Line(Line2 {
                    origin: a,
                    direction: b - a,
                });
                (curve3, curve2, Interval::UNIT)
            }
            Segment::Arc {
                centre,
                radius,
                start,
                end,
            } => {
                let frame = self.cap_frame(u);
                let curve3 = Curve3::Circle(Circle3 {
                    frame: Frame3 {
                        origin: self.place(u, centre.x, centre.y),
                        ..frame
                    },
                    radius,
                });
                let curve2 = Curve2::Circle(Circle2 {
                    frame: Frame2 {
                        origin: centre,
                        x: Vec2::X,
                        y: Vec2::Y,
                    },
                    radius,
                });
                (curve3, curve2, Interval::new(start, end))
            }
        };
        let curve = builder.add_curve3(curve3);
        let edge = builder.topology_mut().add_edge(Edge {
            start: from,
            end: to,
            curve: Some(curve),
        });
        builder.set_edge_interval(edge, interval);
        (edge, curve2, interval)
    }

    /// The surface one segment sweeps, and where its two ends sit in `v`.
    fn wall(
        &self,
        builder: &mut ExactBRepBuilder,
        here: Station,
        there: Station,
        segment: &Segment,
    ) -> GeomResult<Wall> {
        let epsilon = self.epsilon;
        let wall = match *segment {
            Segment::Arc {
                centre,
                radius,
                start,
                end,
            } => {
                // As in the full turn: a tube reaching the axis is a spindle
                // or horn torus, which self-intersects there.
                if radius >= centre.x.abs() - Scalar::EPSILON {
                    return Err(unsupported(
                        "revolved arc whose tube reaches the axis (spindle torus)",
                    ));
                }
                let surface = builder.add_surface(Surface::Torus(Torus {
                    frame: frame_at(self.axis_x, centre.y),
                    major_radius: centre.x,
                    minor_radius: radius,
                }));
                Wall::Ruled {
                    surface,
                    from: start,
                    to: end,
                }
            }
            Segment::Line { .. } if (here.radius - there.radius).abs() <= epsilon => {
                let surface = builder.add_surface(Surface::Cylinder(Cylinder {
                    frame: frame_at(self.axis_x, here.height),
                    radius: here.radius,
                }));
                Wall::Ruled {
                    surface,
                    from: 0.0,
                    to: there.height - here.height,
                }
            }
            Segment::Line { .. } if (here.height - there.height).abs() <= epsilon => {
                let surface = builder.add_surface(Surface::Plane(Plane {
                    frame: frame_at(self.axis_x, here.height),
                }));
                Wall::Planar { surface }
            }
            Segment::Line { .. } => {
                // A cone. Its frame sits at the apex when one end is on the
                // axis, so the apex is `v = 0` exactly rather than a radius
                // rounded a hair below zero, which the surface refuses.
                let basis = if there.radius == 0.0 { there } else { here };
                let surface = builder.add_surface(Surface::Cone(Cone {
                    frame: frame_at(self.axis_x, basis.height),
                    radius: basis.radius,
                    semi_angle: ((there.radius - here.radius) / (there.height - here.height))
                        .atan(),
                }));
                Wall::Ruled {
                    surface,
                    from: here.height - basis.height,
                    to: there.height - basis.height,
                }
            }
        };
        Ok(wall)
    }

    /// The wall face: `A_i`, `Q_i`, `A_{i+1}` reversed, `P_i` reversed.
    fn wall_face(
        &self,
        builder: &mut ExactBRepBuilder,
        wall: &Wall,
        (lower, upper): (Option<EdgeId>, Option<EdgeId>),
        (start_edge, end_edge): (EdgeId, EdgeId),
        (here, there): (Station, Station),
    ) -> FaceId {
        let (u0, u1) = (self.u0, self.u1);
        // Each use: edge, orientation, pcurve, pcurve interval.
        let mut uses: Vec<(EdgeId, Orientation, Curve2, Interval)> = Vec::with_capacity(4);
        let surface = match *wall {
            Wall::Ruled { surface, from, to } => {
                let along = |v: Scalar| {
                    Curve2::Line(Line2 {
                        origin: Vec2::new(0.0, v),
                        direction: Vec2::X,
                    })
                };
                let across = |u: Scalar| {
                    Curve2::Line(Line2 {
                        origin: Vec2::new(u, from),
                        direction: Vec2::new(0.0, to - from),
                    })
                };
                if let Some(edge) = lower {
                    uses.push((
                        edge,
                        Orientation::Forward,
                        along(from),
                        Interval::new(u0, u1),
                    ));
                }
                uses.push((end_edge, Orientation::Forward, across(u1), Interval::UNIT));
                if let Some(edge) = upper {
                    uses.push((
                        edge,
                        Orientation::Reversed,
                        along(to),
                        Interval::new(u1, u0),
                    ));
                }
                uses.push((
                    start_edge,
                    Orientation::Reversed,
                    across(u0),
                    Interval::new(1.0, 0.0),
                ));
                surface
            }
            Wall::Planar { surface } => {
                let circle = |radius: Scalar| {
                    Curve2::Circle(Circle2 {
                        frame: Frame2 {
                            origin: Vec2::ZERO,
                            x: Vec2::X,
                            y: Vec2::Y,
                        },
                        radius,
                    })
                };
                // `frame_at` puts plane x on world X and plane y on -Z, so
                // the radial at angle `u` is `(cos u, sin u)` here.
                let radial = |u: Scalar| {
                    let (sin, cos) = u.sin_cos();
                    let direction = Vec2::new(cos, sin);
                    Curve2::Line(Line2 {
                        origin: direction * here.radius,
                        direction: direction * (there.radius - here.radius),
                    })
                };
                if let Some(edge) = lower {
                    uses.push((
                        edge,
                        Orientation::Forward,
                        circle(here.radius),
                        Interval::new(u0, u1),
                    ));
                }
                uses.push((end_edge, Orientation::Forward, radial(u1), Interval::UNIT));
                if let Some(edge) = upper {
                    uses.push((
                        edge,
                        Orientation::Reversed,
                        circle(there.radius),
                        Interval::new(u1, u0),
                    ));
                }
                uses.push((
                    start_edge,
                    Orientation::Reversed,
                    radial(u0),
                    Interval::new(1.0, 0.0),
                ));
                surface
            }
        };
        let loop_id = add_loop(builder, uses);
        builder.topology_mut().add_face(Face {
            surface: Some(surface),
            bounds: vec![FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation: Orientation::Forward,
        })
    }
}

/// A profile segment in `(r, h)`.
enum Segment {
    Line {
        from: Point2,
        to: Point2,
    },
    /// Circle centre and radius, and the angles of its two ends, `end`
    /// reached from `start` in the arc's own direction.
    Arc {
        centre: Point2,
        radius: Scalar,
        start: Scalar,
        end: Scalar,
    },
}

impl Segment {
    fn of(here: Station, there: Station) -> GeomResult<Self> {
        let from = Point2::new(here.radius, here.height);
        let to = Point2::new(there.radius, there.height);
        if here.bulge == 0.0 {
            return Ok(Self::Line { from, to });
        }
        let arc = crate::extrude_arc::arc_geometry(from, to, here.bulge)?;
        if arc.centre.x <= 0.0 {
            return Err(unsupported(
                "revolved arc centred on the far side of the axis",
            ));
        }
        let start = (from.y - arc.centre.y).atan2(from.x - arc.centre.x);
        Ok(Self::Arc {
            centre: arc.centre,
            radius: arc.radius,
            start,
            end: start + arc.sweep,
        })
    }
}

/// The surface a segment sweeps.
enum Wall {
    /// Cylinder, cone or torus: `(angle, v)` with the segment running from
    /// `v = from` to `v = to`.
    Ruled {
        surface: SurfaceId,
        from: Scalar,
        to: Scalar,
    },
    /// A planar sector at the segment's height, in Cartesian coordinates.
    Planar { surface: SurfaceId },
}

fn add_loop(
    builder: &mut ExactBRepBuilder,
    uses: Vec<(EdgeId, Orientation, Curve2, Interval)>,
) -> axiolid_topology::LoopId {
    let mut edges = Vec::with_capacity(uses.len());
    let mut intervals = Vec::with_capacity(uses.len());
    for (edge, orientation, pcurve, interval) in uses {
        let pcurve = builder.add_curve2(pcurve);
        edges.push(EdgeUse {
            edge,
            orientation,
            pcurve: Some(pcurve),
        });
        intervals.push(interval);
    }
    let loop_id = builder.topology_mut().add_loop(Loop { edges });
    for (index, interval) in intervals.into_iter().enumerate() {
        builder.set_pcurve_interval(loop_id, index, interval);
    }
    loop_id
}

/// A planar end cap: one loop per ring, the first outer.
fn cap<'a>(
    builder: &mut ExactBRepBuilder,
    surface: SurfaceId,
    rings: impl Iterator<Item = &'a [(EdgeId, Curve2, Interval)]>,
    orientation: Orientation,
) -> FaceId {
    let mut bounds = Vec::new();
    for (index, ring) in rings.enumerate() {
        let uses = ring
            .iter()
            .map(|(edge, pcurve, interval)| {
                (*edge, Orientation::Forward, pcurve.clone(), *interval)
            })
            .collect();
        let loop_id = add_loop(builder, uses);
        bounds.push(FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer: index == 0,
        });
    }
    builder.topology_mut().add_face(Face {
        surface: Some(surface),
        bounds,
        orientation,
    })
}
