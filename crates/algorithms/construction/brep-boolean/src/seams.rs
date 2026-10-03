//! Seam edges for faces that wind round their surface without one.
//!
//! A face on a cylinder, cone, sphere or torus may be bounded by loops that
//! go all the way round the surface: a dome's single rim circle, a band
//! between two circles. In the face's parameters such a loop does not close
//! -- it ends a whole turn from where it starts -- so splitting the face in
//! its parameters (`split`) cannot trace regions from it. Before the
//! boolean runs, each such face gets a seam edge along the iso-curve where
//! the loops' parameters wrap, used once each way, and its winding loops
//! become one closed loop in parameters: the lower loop, up the seam, the
//! upper loop (or across a pole or apex), down the seam. Nothing changes
//! for faces that already close.

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_core::{Frame3, Interval, Point2, Point3, Scalar, Tolerance, Vec2};
use axiolid_curve::{Circle3, Curve2, Curve3, Line2, Line3};
use axiolid_evaluate::curve::evaluate2;
use axiolid_evaluate::surface::{evaluate, partials};
use axiolid_surface::Surface;
use axiolid_topology::{Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Vertex, VertexId};
use core::f64::consts::{FRAC_PI_2, TAU};

use crate::support::periods;
use crate::BooleanError;

/// One use of a loop, in traversal order, with its pcurve moved by whole
/// periods where the loop needs it.
#[derive(Debug, Clone)]
struct Use {
    edge: axiolid_topology::EdgeId,
    orientation: Orientation,
    pcurve: Curve2,
    interval: Interval,
    shift: Vec2,
}

impl Use {
    fn at(&self, start: bool) -> Result<Point2, BooleanError> {
        let t = if start {
            self.interval.start
        } else {
            self.interval.end
        };
        Ok(evaluate2(&self.pcurve, t).map_err(|_| BooleanError::Evaluation)? + self.shift)
    }
}

/// A pcurve moved by `d` in parameters: the curve itself for the families
/// that are plain geometry, or its interval for the graphs whose parameter
/// is their first coordinate.
fn moved(curve: &Curve2, interval: Interval, d: Vec2) -> Option<(Curve2, Interval)> {
    if d == Vec2::ZERO {
        return Some((curve.clone(), interval));
    }
    Some(match curve {
        Curve2::Line(l) => (
            Curve2::Line(Line2 {
                origin: l.origin + d,
                ..*l
            }),
            interval,
        ),
        Curve2::Circle(c) => {
            let mut c = *c;
            c.frame.origin += d;
            (Curve2::Circle(c), interval)
        }
        Curve2::Ellipse(e) => {
            let mut e = *e;
            e.frame.origin += d;
            (Curve2::Ellipse(e), interval)
        }
        Curve2::Polyline(p) => {
            let mut p = p.clone();
            for q in &mut p.points {
                *q += d;
            }
            (Curve2::Polyline(p), interval)
        }
        Curve2::BSpline(b) => {
            let mut b = b.clone();
            for q in &mut b.control_points {
                *q += d;
            }
            (Curve2::BSpline(b), interval)
        }
        Curve2::Implicit(c) => (Curve2::Implicit(c.shifted(d.x, d.y)), interval),
        Curve2::Lifted(l) => {
            let mut l = l.clone();
            for g in &mut l.guide {
                *g += d;
            }
            (Curve2::Lifted(l), interval)
        }
        // Graphs over their first coordinate, periodic in it: the same
        // curve, read a whole turn along.
        Curve2::Sinusoid(_) | Curve2::QuadraticGraph(_) if d.y == 0.0 => (
            curve.clone(),
            Interval::new(interval.start + d.x, interval.end + d.x),
        ),
        _ => return None,
    })
}

/// The iso-curve of `surface` at the first parameter `s`, parameterised by
/// the second, and at the second parameter for a torus (`along_v`).
fn iso_curve(surface: &Surface, s: Scalar, along_v: bool) -> Option<Curve3> {
    let (cs, ss) = (s.cos(), s.sin());
    Some(match surface {
        Surface::Cylinder(_) | Surface::Cone(_) | Surface::EllipticalCylinder(_) if !along_v => {
            // Linear in v along a ruling.
            let origin = evaluate(surface, s, 0.0).ok()?;
            let (_, dv) = partials(surface, s, 0.0).ok()?;
            Curve3::Line(Line3 {
                origin,
                direction: dv,
            })
        }
        Surface::Sphere(sp) if !along_v => {
            let f = &sp.frame;
            let x = f.x * cs + f.y * ss;
            Curve3::Circle(Circle3 {
                frame: Frame3 {
                    origin: f.origin,
                    x,
                    y: f.z,
                    z: x.cross(f.z),
                },
                radius: sp.radius,
            })
        }
        Surface::Torus(t) if !along_v => {
            let f = &t.frame;
            let x = f.x * cs + f.y * ss;
            Curve3::Circle(Circle3 {
                frame: Frame3 {
                    origin: f.origin + x * t.major_radius,
                    x,
                    y: f.z,
                    z: x.cross(f.z),
                },
                radius: t.minor_radius,
            })
        }
        Surface::Torus(t) => {
            let f = &t.frame;
            Curve3::Circle(Circle3 {
                frame: Frame3 {
                    origin: f.origin + f.z * (t.minor_radius * ss),
                    ..*f
                },
                radius: t.major_radius + t.minor_radius * cs,
            })
        }
        _ => return None,
    })
}

/// The pole a face's domain reaches past its last winding loop, in the
/// wound direction: a sphere's `v = +-pi/2`, a cone's apex.
fn pole(surface: &Surface, upward: bool) -> Option<Scalar> {
    match surface {
        Surface::Sphere(_) => Some(if upward { FRAC_PI_2 } else { -FRAC_PI_2 }),
        Surface::Cone(c) => {
            let slope = c.semi_angle.tan();
            let apex = -c.radius / slope;
            // The apex is on the side the radius shrinks towards.
            ((slope > 0.0) != upward).then_some(apex)
        }
        _ => None,
    }
}

/// `brep` with a seam edge in every face whose loops wind round its surface
/// without one; `None` when no face needs one.
pub(crate) fn with_seams(
    brep: &ExactBRep,
    tolerance: Tolerance,
) -> Result<Option<ExactBRep>, BooleanError> {
    let topology = brep.topology();
    let eps = tolerance.linear();
    // Per face: its loops as uses, unwrapped, and each loop's winding.
    let mut plans: Vec<Option<Plan>> = Vec::with_capacity(topology.faces().len());
    for face in topology.faces() {
        let surface = face
            .surface
            .and_then(|id| brep.surfaces().get(id.index()))
            .ok_or(BooleanError::DanglingReference)?;
        let (pu, pv) = periods(surface);
        if !pu && !pv {
            plans.push(None);
            continue;
        }
        let mut loops = Vec::new();
        let mut any = false;
        for bound in &face.bounds {
            let wire = topology
                .loops()
                .get(bound.loop_id.index())
                .ok_or(BooleanError::DanglingReference)?;
            let mut uses = Vec::with_capacity(wire.edges.len());
            for (k, u) in wire.edges.iter().enumerate() {
                uses.push(Use {
                    edge: u.edge,
                    orientation: u.orientation,
                    pcurve: u
                        .pcurve
                        .and_then(|id| brep.curves2().get(id.index()))
                        .ok_or(BooleanError::DanglingReference)?
                        .clone(),
                    interval: brep
                        .pcurve_interval(bound.loop_id, k)
                        .ok_or(BooleanError::DanglingReference)?,
                    shift: Vec2::ZERO,
                });
            }
            if bound.orientation == Orientation::Reversed {
                uses.reverse();
                for u in &mut uses {
                    u.orientation = match u.orientation {
                        Orientation::Forward => Orientation::Reversed,
                        Orientation::Reversed => Orientation::Forward,
                    };
                    u.interval = Interval::new(u.interval.end, u.interval.start);
                }
            }
            let winding = unwrap(&mut uses, pu, pv)?;
            any |= winding != Vec2::ZERO;
            loops.push((uses, winding, bound.outer));
        }
        plans.push(any.then_some(Plan { loops }));
    }
    if plans.iter().all(Option::is_none) {
        return Ok(None);
    }
    rebuild(brep, &plans, eps).map(Some)
}

/// A winding face's loops.
struct Plan {
    loops: Vec<(Vec<Use>, Vec2, bool)>,
}

/// Move each use by whole periods so the loop runs on continuously in
/// parameters (stepping over poles), and return the whole periods it gains
/// round the loop: zero for a loop that closes.
fn unwrap(uses: &mut [Use], pu: bool, pv: bool) -> Result<Vec2, BooleanError> {
    let snap = |x: Scalar, periodic: bool| {
        if periodic {
            (x / TAU).round() * TAU
        } else {
            0.0
        }
    };
    for k in 1..uses.len() {
        let end = uses[k - 1].at(false)?;
        let start = uses[k].at(true)?;
        let d = end - start;
        uses[k].shift += Vec2::new(snap(d.x, pu), snap(d.y, pv));
    }
    let (first, last) = (uses[0].at(true)?, uses[uses.len() - 1].at(false)?);
    let d = last - first;
    Ok(Vec2::new(snap(d.x, pu), snap(d.y, pv)))
}

/// A rebuilt copy of `brep` with seams where `plans` asks.
fn rebuild(
    brep: &ExactBRep,
    plans: &[Option<Plan>],
    eps: Scalar,
) -> Result<ExactBRep, BooleanError> {
    let topology = brep.topology();
    let mut b = ExactBRepBuilder::default();
    for c in brep.curves3() {
        b.add_curve3(c.clone());
    }
    for c in brep.curves2() {
        b.add_curve2(c.clone());
    }
    for s in brep.surfaces() {
        b.add_surface(s.clone());
    }
    let mut vertices: Vec<VertexId> = Vec::new();
    for v in topology.vertices() {
        vertices.push(b.topology_mut().add_vertex(*v));
    }
    for (index, e) in topology.edges().iter().enumerate() {
        let id = b.topology_mut().add_edge(e.clone());
        let old = topology
            .edge_id_at(index)
            .ok_or(BooleanError::DanglingReference)?;
        b.set_edge_interval(
            id,
            brep.edge_interval(old)
                .ok_or(BooleanError::DanglingReference)?,
        );
    }
    let vertex_at =
        |b: &mut ExactBRepBuilder, vertices: &mut Vec<VertexId>, p: Point3| -> VertexId {
            let found = b
                .topology_mut()
                .vertices()
                .iter()
                .position(|v| (v.position - p).length() <= eps);
            match found {
                Some(i) => vertices[i],
                None => {
                    let id = b.topology_mut().add_vertex(Vertex { position: p });
                    vertices.push(id);
                    id
                }
            }
        };
    // Loops, face by face; faces keep their order so shells stay valid.
    let mut faces = Vec::with_capacity(topology.faces().len());
    for (fi, face) in topology.faces().iter().enumerate() {
        let surface = face
            .surface
            .and_then(|id| brep.surfaces().get(id.index()))
            .ok_or(BooleanError::DanglingReference)?;
        let add_loop = |b: &mut ExactBRepBuilder,
                        uses: &[Use]|
         -> Result<axiolid_topology::LoopId, BooleanError> {
            let mut edges = Vec::with_capacity(uses.len());
            let mut spans = Vec::with_capacity(uses.len());
            for u in uses {
                let (pcurve, interval) =
                    moved(&u.pcurve, u.interval, u.shift).ok_or(BooleanError::UnclosedSplit)?;
                let pc = b.add_curve2(pcurve);
                edges.push(EdgeUse {
                    edge: u.edge,
                    orientation: u.orientation,
                    pcurve: Some(pc),
                });
                spans.push(interval);
            }
            let id = b.topology_mut().add_loop(Loop { edges });
            for (k, s) in spans.into_iter().enumerate() {
                b.set_pcurve_interval(id, k, s);
            }
            Ok(id)
        };
        let Some(plan) = &plans[fi] else {
            // Copy the face's loops as they are.
            let mut bounds = Vec::new();
            for bound in &face.bounds {
                let wire = topology
                    .loops()
                    .get(bound.loop_id.index())
                    .ok_or(BooleanError::DanglingReference)?;
                let id = b.topology_mut().add_loop(wire.clone());
                for k in 0..wire.edges.len() {
                    b.set_pcurve_interval(
                        id,
                        k,
                        brep.pcurve_interval(bound.loop_id, k)
                            .ok_or(BooleanError::DanglingReference)?,
                    );
                }
                bounds.push(FaceBound {
                    loop_id: id,
                    ..*bound
                });
            }
            faces.push(Face {
                surface: face.surface,
                bounds,
                orientation: face.orientation,
            });
            continue;
        };
        // Winding loops: which way round, where they start.
        let (pu, _) = periods(surface);
        let winding: Vec<usize> = plan
            .loops
            .iter()
            .enumerate()
            .filter(|(_, (_, w, _))| *w != Vec2::ZERO)
            .map(|(i, _)| i)
            .collect();
        let along_v = plan.loops[winding[0]].1.x == 0.0;
        if along_v && pu && plan.loops.iter().any(|(_, w, _)| w.x != 0.0) {
            return Err(BooleanError::UnclosedSplit);
        }
        // Coordinates: `a` the wound parameter, `c` the other.
        let wound = |p: Point2| if along_v { p.y } else { p.x };
        let across = |p: Point2| if along_v { p.x } else { p.y };
        let make = |a: Scalar, c: Scalar| {
            if along_v {
                Point2::new(c, a)
            } else {
                Point2::new(a, c)
            }
        };
        // Every winding loop must wrap at the same seam angle: rotate each to
        // start at a vertex there.
        let mut paths: Vec<(Vec<Use>, Scalar)> = Vec::new();
        let seam = wound(plan.loops[winding[0]].0[0].at(true)?);
        for &i in &winding {
            let (uses, w, _) = &plan.loops[i];
            let turns = |x: Scalar| (x - seam) / TAU;
            let start = (0..uses.len())
                .find(|&k| {
                    uses[k]
                        .at(true)
                        .is_ok_and(|p| (turns(wound(p)) - turns(wound(p)).round()).abs() <= 1e-9)
                })
                .ok_or(BooleanError::UnclosedSplit)?;
            let mut rotated: Vec<Use> = uses[start..].to_vec();
            rotated.extend(uses[..start].iter().cloned());
            // Re-unwrap from the new start so the path starts at the seam.
            let (pu2, pv2) = periods(surface);
            unwrap(&mut rotated, pu2, pv2)?;
            let first = rotated[0].at(true)?;
            let to_seam = seam - wound(first);
            let base = make((to_seam / TAU).round() * TAU, 0.0);
            for u in &mut rotated {
                u.shift += base;
            }
            paths.push((rotated, if along_v { w.y } else { w.x }));
        }
        // Lower loops run the wound way (+), upper ones back (-).
        let lower: Vec<&(Vec<Use>, Scalar)> = paths.iter().filter(|(_, w)| *w > 0.0).collect();
        let upper: Vec<&(Vec<Use>, Scalar)> = paths.iter().filter(|(_, w)| *w < 0.0).collect();
        let at_start = |p: &(Vec<Use>, Scalar)| -> Result<Scalar, BooleanError> {
            Ok(across(p.0[0].at(true)?))
        };
        let iso = |c0: Scalar| -> Result<Point3, BooleanError> {
            let p = make(seam, c0);
            evaluate(surface, p.x, p.y).map_err(|_| BooleanError::Evaluation)
        };
        let curve = iso_curve(surface, seam, along_v).ok_or(BooleanError::UnclosedSplit)?;
        let mut outer: Vec<Use> = Vec::new();
        let (from, to, top_path): (Scalar, Scalar, Option<&(Vec<Use>, Scalar)>) =
            match (lower.as_slice(), upper.as_slice()) {
                ([low], [high]) => {
                    let (c0, mut c1) = (at_start(low)?, at_start(high)?);
                    if c1 < c0 && matches!(surface, Surface::Torus(_)) {
                        c1 += TAU;
                    }
                    if c1 <= c0 {
                        return Err(BooleanError::UnclosedSplit);
                    }
                    (c0, c1, Some(*high))
                }
                ([low], []) => (
                    at_start(low)?,
                    pole(surface, true).ok_or(BooleanError::UnclosedSplit)?,
                    None,
                ),
                ([], [high]) => (
                    pole(surface, false).ok_or(BooleanError::UnclosedSplit)?,
                    at_start(high)?,
                    None,
                ),
                _ => return Err(BooleanError::UnclosedSplit),
            };
        // The seam edge, laid from `from` to `to` along the iso-curve.
        let (va, vb) = (
            vertex_at(&mut b, &mut vertices, iso(from)?),
            vertex_at(&mut b, &mut vertices, iso(to)?),
        );
        let c3 = b.add_curve3(curve);
        let seam_edge = b.topology_mut().add_edge(Edge {
            start: va,
            end: vb,
            curve: Some(c3),
        });
        b.set_edge_interval(seam_edge, Interval::new(from, to));
        let up_line = |at: Scalar| {
            let o = make(at, 0.0);
            let d = make(0.0, 1.0);
            Curve2::Line(Line2 {
                origin: o,
                direction: d,
            })
        };
        let seam_use = |at: Scalar, up: bool| Use {
            edge: seam_edge,
            orientation: if up {
                Orientation::Forward
            } else {
                Orientation::Reversed
            },
            pcurve: up_line(at),
            interval: if up {
                Interval::new(from, to)
            } else {
                Interval::new(to, from)
            },
            shift: Vec2::ZERO,
        };
        match (lower.first(), top_path) {
            (Some(low), Some(high)) => {
                // Lower loop (seam -> seam + 2 pi), up the seam, upper loop
                // (seam + 2 pi -> seam), down the seam.
                outer.extend(low.0.iter().cloned());
                outer.push(seam_use(seam + TAU, true));
                let lift = make(TAU, 0.0)
                    + make(
                        0.0,
                        if to - at_start(high)? > 1e-9 {
                            TAU
                        } else {
                            0.0
                        },
                    );
                outer.extend(high.0.iter().cloned().map(|mut u| {
                    u.shift += lift;
                    u
                }));
                outer.push(seam_use(seam, false));
            }
            (Some(low), None) => {
                // Lower loop, up the seam to the pole, back down.
                outer.extend(low.0.iter().cloned());
                outer.push(seam_use(seam + TAU, true));
                outer.push(seam_use(seam, false));
            }
            (None, _) => {
                // From the pole below: up the seam at seam + 2 pi, along the
                // upper loop back to the seam (lifted a turn), down the seam.
                outer.push(seam_use(seam + TAU, true));
                outer.extend(upper[0].0.iter().cloned().map(|mut u| {
                    u.shift += make(TAU, 0.0);
                    u
                }));
                outer.push(seam_use(seam, false));
            }
        }
        let mut bounds = vec![FaceBound {
            loop_id: add_loop(&mut b, &outer)?,
            orientation: Orientation::Forward,
            outer: true,
        }];
        for (uses, w, _) in &plan.loops {
            if *w == Vec2::ZERO {
                bounds.push(FaceBound {
                    loop_id: add_loop(&mut b, uses)?,
                    orientation: Orientation::Forward,
                    outer: false,
                });
            }
        }
        faces.push(Face {
            surface: face.surface,
            bounds,
            orientation: face.orientation,
        });
    }
    for f in faces {
        b.topology_mut().add_face(f);
    }
    for s in topology.shells() {
        b.topology_mut().add_shell(s.clone());
    }
    for s in topology.solids() {
        b.topology_mut().add_solid(s.clone());
    }
    b.finish().map_err(|_| BooleanError::Assembly)
}
