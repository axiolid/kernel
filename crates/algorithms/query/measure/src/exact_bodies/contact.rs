//! The boundary of a union of items that share patches of face (#229).
//!
//! Two items in exact face contact -- planar faces lying in one plane on
//! opposite sides of it -- share the patch where those faces overlap. That
//! patch is interior to the union and not on its boundary. The union's
//! boundary is the closure of what is left of every face once the faces of
//! other items lying against it are cut away (see [`super`] for the
//! argument), so it is itself a set of faces: the faces without contact as
//! they are, and each face in contact replaced by its free region.
//!
//! # Exact
//!
//! Contact is cut only where it is exact in the B-rep's own numbers: both
//! faces lie on one plane normal to a coordinate axis, with every vertex,
//! line, circle and the plane itself at the very same coordinate along that
//! axis (which a placement turning about that axis keeps). Dropping that
//! coordinate is then an exact map onto the plane, the faces' boundaries
//! (straight edges between their vertices, circular arcs) keep their
//! coordinates, and the plane is cut by all of them at once by
//! [`ArcArrangement`], whose decisions are exact (ADR 0070): which pieces
//! cross, coincide, and which face contains what. Two items abutting on
//! their common vertices leave no sliver between them. The free regions
//! become planar faces on that plane, charted by the two kept coordinates.
//!
//! A face bounded by anything but lines and circles is not cut: the pair
//! is refused by the caller.

use std::collections::HashMap;

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_core::{Frame2, Frame3, Interval, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Circle3, Curve2, Curve3, Line2, Line3};
use axiolid_overlay::{ArcArrangement, ArcRing, ArcVertex};
use axiolid_surface::{Plane, Surface};
use axiolid_topology::{Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Vertex};

/// A plane normal to coordinate axis `axis`, at `level` along it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct AxisPlane {
    pub(super) axis: usize,
    pub(super) level: Scalar,
}

impl AxisPlane {
    /// The two kept coordinates, ordered so that a counter-clockwise turn
    /// in them is a positive turn about the axis.
    fn project(self, p: Point3) -> Point2 {
        match self.axis {
            0 => Point2::new(p.y, p.z),
            1 => Point2::new(p.z, p.x),
            _ => Point2::new(p.x, p.y),
        }
    }

    fn project_vec(self, v: Vec3) -> Vec2 {
        match self.axis {
            0 => Vec2::new(v.y, v.z),
            1 => Vec2::new(v.z, v.x),
            _ => Vec2::new(v.x, v.y),
        }
    }

    /// The point of the plane over `q`; exact.
    fn lift(self, q: Point2) -> Point3 {
        match self.axis {
            0 => Point3::new(self.level, q.x, q.y),
            1 => Point3::new(q.y, self.level, q.x),
            _ => Point3::new(q.x, q.y, self.level),
        }
    }

    fn lift_vec(self, v: Vec2) -> Vec3 {
        match self.axis {
            0 => Vec3::new(0.0, v.x, v.y),
            1 => Vec3::new(v.y, 0.0, v.x),
            _ => Vec3::new(v.x, v.y, 0.0),
        }
    }

    /// The plane charted by its kept coordinates: `S(u, v)` is the point
    /// over `(u, v)`, exactly.
    fn frame(self) -> Frame3 {
        let unit = |k: usize| {
            let mut v = Vec3::ZERO;
            v[k] = 1.0;
            v
        };
        let (x, y) = match self.axis {
            0 => (unit(1), unit(2)),
            1 => (unit(2), unit(0)),
            _ => (unit(0), unit(1)),
        };
        Frame3 {
            origin: self.lift(Point2::ZERO),
            x,
            y,
            z: unit(self.axis),
        }
    }
}

/// Why a face in contact cannot be cut.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Uncut {
    /// The face's plane is not normal to a coordinate axis.
    NotAxisNormal,
    /// The face's plane is normal to an axis, but its boundary is not at
    /// the plane's coordinate along it in its own numbers.
    OffLevel,
    /// An edge of the face is neither a line nor a circle.
    Edge,
}

/// The axis plane a planar face lies on in its own numbers: its plane and
/// every vertex and edge of its boundary at the same coordinate along one
/// axis, each edge a line or a circle.
pub(super) fn axis_plane(brep: &ExactBRep, face: usize) -> Result<AxisPlane, Uncut> {
    let normal = brep
        .topology()
        .faces()
        .get(face)
        .and_then(|f| f.surface)
        .and_then(|id| brep.surfaces().get(id.index()))
        .and_then(|surface| match surface {
            Surface::Plane(plane) => {
                let frame = plane.frame;
                (0..3).find(|&k| frame.x[k] == 0.0 && frame.y[k] == 0.0)
            }
            _ => None,
        });
    if normal.is_none() {
        return Err(Uncut::NotAxisNormal);
    }
    let found = on_axis_plane(brep, face).ok_or(Uncut::OffLevel)?;
    let topology = brep.topology();
    for bound in &topology.faces()[face].bounds {
        for use_ in &topology.loops()[bound.loop_id.index()].edges {
            let curve = topology.edges()[use_.edge.index()]
                .curve
                .and_then(|id| brep.curves3().get(id.index()));
            if !matches!(curve, Some(Curve3::Line(_) | Curve3::Circle(_))) {
                return Err(Uncut::Edge);
            }
        }
    }
    Ok(found)
}

fn on_axis_plane(brep: &ExactBRep, face: usize) -> Option<AxisPlane> {
    let topology = brep.topology();
    let face = topology.faces().get(face)?;
    let Surface::Plane(plane) = brep.surfaces().get(face.surface?.index())? else {
        return None;
    };
    let frame = plane.frame;
    let axis = (0..3).find(|&k| frame.x[k] == 0.0 && frame.y[k] == 0.0)?;
    let level = frame.origin[axis];
    let at = |p: Point3| p[axis] == level;
    let flat = |v: Vec3| v[axis] == 0.0;
    for bound in &face.bounds {
        let wire = topology.loops().get(bound.loop_id.index())?;
        for use_ in &wire.edges {
            let edge = topology.edges().get(use_.edge.index())?;
            for vertex in [edge.start, edge.end] {
                if !at(topology.vertices().get(vertex.index())?.position) {
                    return None;
                }
            }
            match brep.curves3().get(edge.curve?.index())? {
                Curve3::Line(line) if at(line.origin) && flat(line.direction) => {}
                Curve3::Circle(circle)
                    if at(circle.frame.origin) && flat(circle.frame.x) && flat(circle.frame.y) => {}
                Curve3::Line(_) | Curve3::Circle(_) => return None,
                // Judged by `axis_plane`: not cut, whatever its plane.
                _ => {}
            }
        }
    }
    Some(AxisPlane { axis, level })
}

/// A face's boundary loops as rings on the plane, each with whether it is
/// the outer one. `None` for an edge the arrangement cannot take.
fn rings(brep: &ExactBRep, face: usize, plane: AxisPlane) -> Option<Vec<(ArcRing, bool)>> {
    let topology = brep.topology();
    let face = topology.faces().get(face)?;
    let mut out = Vec::with_capacity(face.bounds.len());
    for bound in &face.bounds {
        let wire = topology.loops().get(bound.loop_id.index())?;
        let mut vertices = Vec::with_capacity(wire.edges.len());
        for use_ in &wire.edges {
            let edge = topology.edges().get(use_.edge.index())?;
            let forward = use_.orientation == Orientation::Forward;
            let start = if forward { edge.start } else { edge.end };
            let point = plane.project(topology.vertices().get(start.index())?.position);
            let bulge = match brep.curves3().get(edge.curve?.index())? {
                Curve3::Line(_) => 0.0,
                Curve3::Circle(circle) => {
                    let span = brep.edge_interval(use_.edge)?;
                    let sweep = if forward {
                        span.end - span.start
                    } else {
                        span.start - span.end
                    };
                    let (x, y) = (
                        plane.project_vec(circle.frame.x),
                        plane.project_vec(circle.frame.y),
                    );
                    let turn = (x.x * y.y - x.y * y.x).signum() * sweep.signum();
                    if sweep.abs() >= core::f64::consts::TAU - 1e-12 {
                        // A whole circle: as two half turns, exactly.
                        if wire.edges.len() != 1 {
                            return None;
                        }
                        let centre = plane.project(circle.frame.origin);
                        // A frame turned in the plane keeps its axes of
                        // unit length to rounding.
                        let unit = |v: Vec3| (v.length() - 1.0).abs() <= 1e-12;
                        if !unit(circle.frame.x) || !unit(circle.frame.y) {
                            return None;
                        }
                        vertices = ArcRing::circle(centre, circle.radius).vertices;
                        break;
                    }
                    turn * (sweep.abs() / 4.0).tan()
                }
                _ => return None,
            };
            vertices.push(ArcVertex::bulged(point, bulge));
        }
        out.push((ArcRing::new(vertices), bound.outer));
    }
    Some(out)
}

/// What is left of face `face` of `brep` once the faces `partners` (of
/// other items, on the same plane) are cut away: `Ok(None)` when they do
/// not overlap it, else its free regions as rings on the plane, each an
/// outer ring and its holes. `Err` when a boundary cannot be cut.
pub(super) fn free_regions(
    brep: &ExactBRep,
    face: usize,
    partners: &[(&ExactBRep, usize)],
    plane: AxisPlane,
) -> Result<Option<Vec<Vec<ArcRing>>>, ()> {
    let own = rings(brep, face, plane).ok_or(())?;
    let mut all: Vec<ArcRing> = own.iter().map(|(ring, _)| ring.clone()).collect();
    // Ring indices of the face, and of each partner: (outer, holes).
    let split = |rings: &[(ArcRing, bool)], first: usize| -> (Vec<usize>, Vec<usize>) {
        let mut outer = Vec::new();
        let mut holes = Vec::new();
        for (k, (_, is_outer)) in rings.iter().enumerate() {
            if *is_outer {
                outer.push(first + k);
            } else {
                holes.push(first + k);
            }
        }
        (outer, holes)
    };
    let mine = split(&own, 0);
    let mut theirs = Vec::with_capacity(partners.len());
    for &(other, other_face) in partners {
        let found = rings(other, other_face, plane).ok_or(())?;
        theirs.push(split(&found, all.len()));
        all.extend(found.into_iter().map(|(ring, _)| ring));
    }
    let arrangement = ArcArrangement::new(&all, Tolerance::ZERO).map_err(|_| ())?;
    let inside = |flags: &[bool], (outer, holes): &(Vec<usize>, Vec<usize>)| {
        outer.iter().any(|&k| flags[k]) && !holes.iter().any(|&k| flags[k])
    };
    let covered = |flags: &[bool]| theirs.iter().any(|part| inside(flags, part));
    let mut overlap = false;
    for edge in arrangement.edges() {
        let sides = (0..2).map(|side| {
            (0..all.len())
                .map(|ring| {
                    if side == 0 {
                        edge.inside_left(ring)
                    } else {
                        edge.inside_right(ring)
                    }
                })
                .collect::<Vec<bool>>()
        });
        for flags in sides {
            if inside(&flags, &mine) && covered(&flags) {
                overlap = true;
            }
        }
    }
    if !overlap {
        return Ok(None);
    }
    let regions = arrangement
        .regions(|flags| inside(flags, &mine) && !covered(flags))
        .map_err(|_| ())?;
    let mut out: Vec<Vec<ArcRing>> = regions
        .iter()
        .map(|region| {
            let mut rings = vec![canonical(arrangement.ring(&region.outer))];
            let mut holes: Vec<ArcRing> = region
                .holes
                .iter()
                .map(|hole| canonical(arrangement.ring(hole)))
                .collect();
            holes.sort_by(|p, q| first_key(p).total_cmp(&first_key(q)));
            rings.extend(holes);
            rings
        })
        .collect();
    // Regions in a fixed order, so a translated copy lists them alike.
    out.sort_by(|p, q| first_key(&p[0]).total_cmp(&first_key(&q[0])));
    Ok(Some(out))
}

/// A ring started at its lowest vertex (by `x`, then `y`), so a translated
/// copy is traversed alike.
fn canonical(ring: ArcRing) -> ArcRing {
    let n = ring.vertices.len();
    let start = (0..n)
        .min_by(|&i, &j| {
            let (p, q) = (ring.vertices[i].point, ring.vertices[j].point);
            p.x.total_cmp(&q.x).then(p.y.total_cmp(&q.y))
        })
        .unwrap_or(0);
    ArcRing::new((0..n).map(|k| ring.vertices[(start + k) % n]).collect())
}

fn first_key(ring: &ArcRing) -> Scalar {
    ring.vertices
        .first()
        .map_or(0.0, |v| v.point.x + 1e-3 * v.point.y)
}

/// Copies faces of the items into one B-rep, with their edges, and adds
/// free regions as new planar faces.
#[derive(Default)]
pub(super) struct Assembler {
    builder: ExactBRepBuilder,
    surfaces: HashMap<usize, axiolid_brep::SurfaceId>,
    curves3: HashMap<usize, axiolid_brep::Curve3Id>,
    curves2: HashMap<usize, axiolid_brep::Curve2Id>,
    vertices: HashMap<usize, axiolid_topology::VertexId>,
    edges: HashMap<usize, axiolid_topology::EdgeId>,
}

impl Assembler {
    /// Start a new item: its catalogs are its own.
    pub(super) fn next_item(&mut self) {
        self.surfaces.clear();
        self.curves3.clear();
        self.curves2.clear();
        self.vertices.clear();
        self.edges.clear();
    }

    pub(super) fn counts(&mut self) -> (usize, usize) {
        let topology = self.builder.topology_mut();
        (topology.faces().len(), topology.edges().len())
    }

    /// Copy face `face` of `brep` as it is.
    pub(super) fn copy_face(&mut self, brep: &ExactBRep, face: usize) -> Option<()> {
        let topology = brep.topology();
        let source = topology.faces().get(face)?;
        let surface_index = source.surface?.index();
        let surface = match self.surfaces.get(&surface_index) {
            Some(id) => *id,
            None => {
                let id = self
                    .builder
                    .add_surface(brep.surfaces().get(surface_index)?.clone());
                self.surfaces.insert(surface_index, id);
                id
            }
        };
        let mut bounds = Vec::with_capacity(source.bounds.len());
        for bound in &source.bounds {
            let wire = topology.loops().get(bound.loop_id.index())?;
            let mut uses = Vec::with_capacity(wire.edges.len());
            for use_ in &wire.edges {
                let edge = self.copy_edge(brep, use_.edge)?;
                let pcurve = match use_.pcurve {
                    Some(id) => Some(match self.curves2.get(&id.index()) {
                        Some(copied) => *copied,
                        None => {
                            let copied = self
                                .builder
                                .add_curve2(brep.curves2().get(id.index())?.clone());
                            self.curves2.insert(id.index(), copied);
                            copied
                        }
                    }),
                    None => None,
                };
                uses.push(EdgeUse {
                    edge,
                    orientation: use_.orientation,
                    pcurve,
                });
            }
            let loop_id = self.builder.topology_mut().add_loop(Loop { edges: uses });
            for index in 0..wire.edges.len() {
                if let Some(interval) = brep.pcurve_interval(bound.loop_id, index) {
                    self.builder.set_pcurve_interval(loop_id, index, interval);
                }
            }
            bounds.push(FaceBound { loop_id, ..*bound });
        }
        self.builder.topology_mut().add_face(Face {
            surface: Some(surface),
            bounds,
            orientation: source.orientation,
        });
        Some(())
    }

    fn copy_edge(
        &mut self,
        brep: &ExactBRep,
        edge: axiolid_topology::EdgeId,
    ) -> Option<axiolid_topology::EdgeId> {
        if let Some(copied) = self.edges.get(&edge.index()) {
            return Some(*copied);
        }
        let topology = brep.topology();
        let source = topology.edges().get(edge.index())?;
        let vertex = |id: axiolid_topology::VertexId, this: &mut Self| -> Option<_> {
            if let Some(copied) = this.vertices.get(&id.index()) {
                return Some(*copied);
            }
            let position = topology.vertices().get(id.index())?.position;
            let copied = this.builder.topology_mut().add_vertex(Vertex { position });
            this.vertices.insert(id.index(), copied);
            Some(copied)
        };
        let start = vertex(source.start, self)?;
        let end = vertex(source.end, self)?;
        let curve = match source.curve {
            Some(id) => Some(match self.curves3.get(&id.index()) {
                Some(copied) => *copied,
                None => {
                    let copied = self
                        .builder
                        .add_curve3(brep.curves3().get(id.index())?.clone());
                    self.curves3.insert(id.index(), copied);
                    copied
                }
            }),
            None => None,
        };
        let copied = self
            .builder
            .topology_mut()
            .add_edge(Edge { start, end, curve });
        if let Some(interval) = brep.edge_interval(edge) {
            self.builder.set_edge_interval(copied, interval);
        }
        self.edges.insert(edge.index(), copied);
        Some(copied)
    }

    /// Add a free region -- an outer ring and its holes on `plane` -- as a
    /// planar face charted by the plane's kept coordinates.
    pub(super) fn add_region(&mut self, rings: &[ArcRing], plane: AxisPlane) {
        let surface = self.builder.add_surface(Surface::Plane(Plane {
            frame: plane.frame(),
        }));
        let mut bounds = Vec::with_capacity(rings.len());
        for (k, ring) in rings.iter().enumerate() {
            let n = ring.vertices.len();
            let vertices: Vec<_> = ring
                .vertices
                .iter()
                .map(|v| {
                    self.builder.topology_mut().add_vertex(Vertex {
                        position: plane.lift(v.point),
                    })
                })
                .collect();
            let mut uses = Vec::with_capacity(n);
            let mut intervals = Vec::with_capacity(n);
            for i in 0..n {
                let (from, to) = (ring.vertices[i], ring.vertices[(i + 1) % n]);
                let (curve3, curve2, interval) = piece(from, to, plane);
                let curve3 = self.builder.add_curve3(curve3);
                let curve2 = self.builder.add_curve2(curve2);
                let edge = self.builder.topology_mut().add_edge(Edge {
                    start: vertices[i],
                    end: vertices[(i + 1) % n],
                    curve: Some(curve3),
                });
                self.builder.set_edge_interval(edge, interval);
                uses.push(EdgeUse {
                    edge,
                    orientation: Orientation::Forward,
                    pcurve: Some(curve2),
                });
                intervals.push(interval);
            }
            let loop_id = self.builder.topology_mut().add_loop(Loop { edges: uses });
            for (index, interval) in intervals.into_iter().enumerate() {
                self.builder.set_pcurve_interval(loop_id, index, interval);
            }
            bounds.push(FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: k == 0,
            });
        }
        self.builder.topology_mut().add_face(Face {
            surface: Some(surface),
            bounds,
            orientation: Orientation::Forward,
        });
    }

    pub(super) fn finish(self) -> Option<ExactBRep> {
        self.builder.finish().ok()
    }
}

/// One edge of a ring, from `from` to `to`, as a curve in space, its
/// pcurve on the plane's chart, and the parameter span both share.
fn piece(from: ArcVertex, to: ArcVertex, plane: AxisPlane) -> (Curve3, Curve2, Interval) {
    let (p, q) = (from.point, to.point);
    let chord = Vec2::new(q.x - p.x, q.y - p.y);
    let bulge = from.bulge;
    if bulge == 0.0 {
        return (
            Curve3::Line(Line3 {
                origin: plane.lift(p),
                direction: plane.lift_vec(chord),
            }),
            Curve2::Line(Line2 {
                origin: p,
                direction: chord,
            }),
            Interval::UNIT,
        );
    }
    // The arc's centre, from the chord and the signed included angle
    // `4 atan(bulge)`, positive counter-clockwise.
    let sweep = 4.0 * bulge.atan();
    let length = chord.length();
    let radius = length * (1.0 + bulge * bulge) / (4.0 * bulge.abs());
    let left = Vec2::new(-chord.y, chord.x) / length;
    let sagitta = 0.5 * bulge.abs() * length;
    let mid = Point2::new(0.5 * (p.x + q.x), 0.5 * (p.y + q.y));
    let offset = bulge.signum() * (radius - sagitta);
    let centre = Point2::new(mid.x + left.x * offset, mid.y + left.y * offset);
    let start = (p.y - centre.y).atan2(p.x - centre.x);
    let interval = Interval {
        start,
        end: start + sweep,
    };
    let frame = plane.frame();
    (
        Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: plane.lift(centre),
                ..frame
            },
            radius,
        }),
        Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: centre,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius,
        }),
        interval,
    )
}
