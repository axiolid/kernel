//! Rigid placement of an exact B-rep (#223).
//!
//! # What a rigid transform preserves
//!
//! A rotation, a translation and a reflection map every support family the
//! catalog carries onto the same family: a plane stays a plane, a cylinder
//! a cylinder of the same radius, a circle a circle, a B-spline a B-spline
//! over the same knots and weights with its control points moved. So
//! [`ExactBRep::transformed`] is a new exact B-rep of the same topology and
//! the same families, never a tessellation or a refit. Any other affine map
//! -- a scale, a shear -- turns circles into ellipses or worse and is
//! refused, not approximated.
//!
//! Coordinates are mapped in `f64`: each frame axis, origin, control point
//! and vertex is the transform applied once, rounded once. The result is the
//! exact B-rep of those numbers, as every constructed B-rep is the exact
//! B-rep of the numbers it was built from; the rounding is a few units in
//! the last place of the coordinates, not a change of representation.
//!
//! # Proper rotations
//!
//! A rotation moves every frame rigidly. Every point keeps its parameters,
//! so edge intervals, pcurves and pcurve intervals are unchanged, and the
//! curves that carry their own surface (sections, traced and lifted
//! curves) move with their carriers.
//!
//! # Reflections
//!
//! A reflection `M` (determinant -1) keeps a frame orthonormal but makes it
//! left-handed, which the evaluators refuse. So each frame is re-chosen
//! right-handed, with its parameterisation adjusted so every point keeps
//! its place:
//!
//! - a plane, circle or ellipse takes `(M x, M y, -M z)`: points keep their
//!   parameters;
//! - a cylinder, elliptical cylinder, cone, sphere and torus take
//!   `(M x, -M y, M z)`, which reads the angle about the axis backwards, so
//!   a point at `u` sits at `2 pi - u` and each pcurve on the face is
//!   reflected in `u` to match;
//! - a B-spline keeps its parameters; its control points are reflected.
//!
//! A reflection reverses the winding of every loop about the mirrored
//! outward normal, so every face is flipped: the solid stays outward
//! oriented and measures a positive volume. A pcurve family that cannot be
//! reflected in `u` in closed form here (graphs of quadric and torus
//! sections, traced and lifted curves, intrinsic curves) on a curved face,
//! and a curve that carries such a surface, are refused by name.
//!
//! # Refused for every transform
//!
//! An elevated alignment curve is a plan curve with a height law: only a
//! motion that keeps the vertical axis keeps that form, and re-expressing
//! the law is not implemented. An unknown support family is refused too.

use std::f64::consts::TAU;
use std::fmt;

use axiolid_core::{Frame2, Frame3, Interval, Point2, Scalar, Transform3, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, BSplineSurface, Carrier, Circle2, Circle3, Curve2, Curve3, Ellipse2, Ellipse3,
    ImplicitSection3, Intrinsic3, LiftedCurve2, Line2, Line3, PairNode, PairSection3, Polyline,
    RuledCarrier, RuledSection3, Sinusoid2, TorusCarrier, TorusSection3,
};
use axiolid_surface::{Cone, Cylinder, EllipticalCylinder, Plane, Sphere, Surface, Torus};
use axiolid_topology::{Edge, Face, Loop, Orientation, Vertex};

use crate::{ExactBRep, ExactTopology};

/// How far the columns of a transform's linear part may stray from an
/// orthonormal basis, in each dot product, and still be read as rigid.
///
/// A rotation assembled from unit axes or from sines and cosines is
/// orthonormal to a few units in the last place; a matrix off by more than
/// this carries a scale or shear that no rigid copy of the solid represents.
pub const RIGID_TOLERANCE: Scalar = 1e-12;

/// Why an exact B-rep could not be placed.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TransformError {
    /// The transform has a non-finite entry.
    NonFinite,
    /// The linear part is not orthonormal within [`RIGID_TOLERANCE`]: a
    /// scale or shear, which no rigid copy of the solid represents.
    NotRigid,
    /// A support the transform cannot carry exactly, named.
    Unsupported(&'static str),
}

impl fmt::Display for TransformError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite => f.write_str("the transform has a non-finite entry"),
            Self::NotRigid => f.write_str(
                "the transform is not rigid: its linear part is not orthonormal \
                 (a scale or shear)",
            ),
            Self::Unsupported(what) => write!(f, "{what} cannot be placed exactly"),
        }
    }
}

impl std::error::Error for TransformError {}

/// A validated rigid motion.
#[derive(Debug, Clone, Copy)]
struct Rigid {
    transform: Transform3,
    mirror: bool,
}

impl Rigid {
    fn new(transform: &Transform3) -> Result<Self, TransformError> {
        let m = transform.matrix3;
        let columns = [m.x_axis, m.y_axis, m.z_axis];
        if !columns.iter().all(|c| c.is_finite()) || !transform.translation.is_finite() {
            return Err(TransformError::NonFinite);
        }
        for i in 0..3 {
            for j in i..3 {
                let want = if i == j { 1.0 } else { 0.0 };
                if (columns[i].dot(columns[j]) - want).abs() > RIGID_TOLERANCE {
                    return Err(TransformError::NotRigid);
                }
            }
        }
        Ok(Self {
            transform: *transform,
            mirror: m.determinant() < 0.0,
        })
    }

    fn point(&self, p: Vec3) -> Vec3 {
        self.transform.transform_point3(p)
    }

    fn vector(&self, v: Vec3) -> Vec3 {
        self.transform.transform_vector3(v)
    }

    /// A frame whose points keep their parameters in the plane of `x, y`:
    /// `(M x, M y, +-M z)`, right-handed either way.
    fn planar_frame(&self, frame: &Frame3) -> Frame3 {
        let z = self.vector(frame.z);
        Frame3 {
            origin: self.point(frame.origin),
            x: self.vector(frame.x),
            y: self.vector(frame.y),
            z: if self.mirror { -z } else { z },
        }
    }

    /// A frame about an axis: `(M x, -+M y, M z)`. Under a reflection the
    /// angle about `z` reads backwards, `u -> 2 pi - u`.
    fn axial_frame(&self, frame: &Frame3) -> Frame3 {
        let y = self.vector(frame.y);
        Frame3 {
            origin: self.point(frame.origin),
            x: self.vector(frame.x),
            y: if self.mirror { -y } else { y },
            z: self.vector(frame.z),
        }
    }

    fn unsupported_mirrored(&self, what: &'static str) -> Result<(), TransformError> {
        if self.mirror {
            Err(TransformError::Unsupported(what))
        } else {
            Ok(())
        }
    }
}

/// How a transformed curve's parameter relates to the original's.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Reparam {
    /// The same point at the same parameter.
    Same,
    /// The point at `t` is now at `-t`.
    Negated,
    /// The point at `t` is now at `2 pi - t`: a graph over the angle `u`,
    /// whose parameter is `u` itself, reflected in `u`.
    Reflected,
}

impl Reparam {
    fn interval(self, interval: Interval) -> Interval {
        match self {
            Self::Same => interval,
            Self::Negated => Interval::new(-interval.start, -interval.end),
            Self::Reflected => Interval::new(TAU - interval.start, TAU - interval.end),
        }
    }
}

impl ExactBRep {
    /// This B-rep moved by the rigid motion `transform` (#223).
    ///
    /// Every vertex, curve and surface is mapped onto the same family; the
    /// topology, structural names and parameter intervals are kept, apart
    /// from the reparameterisation a reflection needs (see the module
    /// docs). A reflection flips every face, so a solid stays outward
    /// oriented. Coordinates are mapped once in `f64`; nothing is
    /// tessellated or refitted.
    ///
    /// # Errors
    ///
    /// [`TransformError::NonFinite`] for a non-finite transform,
    /// [`TransformError::NotRigid`] for a linear part that is not
    /// orthonormal within [`RIGID_TOLERANCE`], and
    /// [`TransformError::Unsupported`] naming a support the transform
    /// cannot carry exactly.
    pub fn transformed(&self, transform: &Transform3) -> Result<ExactBRep, TransformError> {
        let rigid = Rigid::new(transform)?;
        let source = &self.topology;

        let mut curves3 = Vec::with_capacity(self.curves3.len());
        let mut curve_reparam = Vec::with_capacity(self.curves3.len());
        for curve in &self.curves3 {
            let (curve, reparam) = curve3(curve, &rigid)?;
            curves3.push(curve);
            curve_reparam.push(reparam);
        }
        let mut surfaces = Vec::with_capacity(self.surfaces.len());
        let mut axial = Vec::with_capacity(self.surfaces.len());
        for surface in &self.surfaces {
            let (surface, reflected) = surface_of(surface, &rigid)?;
            surfaces.push(surface);
            axial.push(reflected);
        }

        // A pcurve is reflected in `u` when its face's surface reads its
        // angle backwards; a pcurve shared by faces that disagree cannot
        // be both, and is refused.
        let mut pcurve_reflect: Vec<Option<bool>> = vec![None; self.curves2.len()];
        for face in source.faces() {
            let Some(surface) = face.surface else {
                continue;
            };
            let reflect = axial.get(surface.index()).copied().unwrap_or(false);
            for bound in &face.bounds {
                let Some(wire) = source.loops().get(bound.loop_id.index()) else {
                    continue;
                };
                for use_ in &wire.edges {
                    let Some(pcurve) = use_.pcurve else {
                        continue;
                    };
                    let Some(slot) = pcurve_reflect.get_mut(pcurve.index()) else {
                        continue;
                    };
                    match slot {
                        Some(existing) if *existing != reflect => {
                            return Err(TransformError::Unsupported(
                                "a pcurve shared by faces of different families under a reflection",
                            ));
                        }
                        _ => *slot = Some(reflect),
                    }
                }
            }
        }
        let mut curves2 = Vec::with_capacity(self.curves2.len());
        let mut pcurve_reparam = Vec::with_capacity(self.curves2.len());
        for (index, curve) in self.curves2.iter().enumerate() {
            let (curve, reparam) = curve2(curve, &rigid, pcurve_reflect[index].unwrap_or(false))?;
            curves2.push(curve);
            pcurve_reparam.push(reparam);
        }

        let mut topology = ExactTopology::default();
        for vertex in source.vertices() {
            topology.add_vertex(Vertex {
                position: rigid.point(vertex.position),
            });
        }
        let mut edge_intervals = self.edge_intervals.clone();
        for (index, edge) in source.edges().iter().enumerate() {
            let id = topology.add_edge(Edge {
                start: edge.start,
                end: edge.end,
                curve: edge.curve,
            });
            let reparam = edge
                .curve
                .and_then(|curve| curve_reparam.get(curve.index()).copied())
                .unwrap_or(Reparam::Same);
            debug_assert_eq!(source.edge_id_at(index), Some(id));
            if let Some(interval) = edge_intervals.get_mut(&id) {
                *interval = reparam.interval(*interval);
            }
        }
        let mut pcurve_intervals = self.pcurve_intervals.clone();
        for (index, wire) in source.loops().iter().enumerate() {
            let id = topology.add_loop(Loop {
                edges: wire.edges.clone(),
            });
            debug_assert_eq!(source.loop_id_at(index), Some(id));
            for (use_index, use_) in wire.edges.iter().enumerate() {
                let reparam = use_
                    .pcurve
                    .and_then(|curve| pcurve_reparam.get(curve.index()).copied())
                    .unwrap_or(Reparam::Same);
                if let Some(interval) = pcurve_intervals.get_mut(&(id, use_index)) {
                    *interval = reparam.interval(*interval);
                }
            }
        }
        for face in source.faces() {
            let orientation = match (face.orientation, rigid.mirror) {
                (orientation, false) => orientation,
                (Orientation::Forward, true) => Orientation::Reversed,
                (Orientation::Reversed, true) => Orientation::Forward,
            };
            topology.add_face(Face {
                surface: face.surface,
                bounds: face.bounds.clone(),
                orientation,
            });
        }
        for shell in source.shells() {
            topology.add_shell(shell.clone());
        }
        for solid in source.solids() {
            topology.add_solid(solid.clone());
        }

        Ok(ExactBRep {
            topology,
            curves3,
            curves2,
            surfaces,
            edge_intervals,
            pcurve_intervals,
            face_names: self.face_names.clone(),
            edge_names: self.edge_names.clone(),
        })
    }
}

/// The surface moved, and whether its angle now reads backwards.
fn surface_of(surface: &Surface, rigid: &Rigid) -> Result<(Surface, bool), TransformError> {
    let axial = rigid.mirror;
    Ok(match surface {
        Surface::Plane(p) => (
            Surface::Plane(Plane {
                frame: rigid.planar_frame(&p.frame),
            }),
            false,
        ),
        Surface::Cylinder(c) => (
            Surface::Cylinder(Cylinder {
                frame: rigid.axial_frame(&c.frame),
                ..*c
            }),
            axial,
        ),
        Surface::EllipticalCylinder(c) => (
            Surface::EllipticalCylinder(EllipticalCylinder {
                frame: rigid.axial_frame(&c.frame),
                ..*c
            }),
            axial,
        ),
        Surface::Cone(c) => (
            Surface::Cone(Cone {
                frame: rigid.axial_frame(&c.frame),
                ..*c
            }),
            axial,
        ),
        Surface::Sphere(s) => (
            Surface::Sphere(Sphere {
                frame: rigid.axial_frame(&s.frame),
                ..*s
            }),
            axial,
        ),
        Surface::Torus(t) => (
            Surface::Torus(Torus {
                frame: rigid.axial_frame(&t.frame),
                ..*t
            }),
            axial,
        ),
        Surface::BSpline(b) => (Surface::BSpline(spline_surface(b, rigid)), false),
        _ => return Err(TransformError::Unsupported("an unknown surface family")),
    })
}

fn spline_surface(surface: &BSplineSurface, rigid: &Rigid) -> BSplineSurface {
    BSplineSurface {
        control_points: surface
            .control_points
            .iter()
            .map(|row| row.iter().map(|p| rigid.point(*p)).collect())
            .collect(),
        ..surface.clone()
    }
}

fn carrier(value: &Carrier, rigid: &Rigid) -> Result<Carrier, TransformError> {
    Ok(match value {
        Carrier::Plane(frame) => Carrier::Plane(rigid.planar_frame(frame)),
        Carrier::Ruled(ruled) => {
            rigid.unsupported_mirrored("a curve on a ruled carrier under a reflection")?;
            Carrier::Ruled(RuledCarrier {
                frame: rigid.axial_frame(&ruled.frame),
                ..*ruled
            })
        }
        Carrier::Sphere { frame, radius } => {
            rigid.unsupported_mirrored("a curve on a spherical carrier under a reflection")?;
            Carrier::Sphere {
                frame: rigid.axial_frame(frame),
                radius: *radius,
            }
        }
        Carrier::Torus(torus) => {
            rigid.unsupported_mirrored("a curve on a toroidal carrier under a reflection")?;
            Carrier::Torus(TorusCarrier {
                frame: rigid.axial_frame(&torus.frame),
                ..*torus
            })
        }
        Carrier::Spline(surface) => Carrier::Spline(Box::new(spline_surface(surface, rigid))),
    })
}

fn curve3(curve: &Curve3, rigid: &Rigid) -> Result<(Curve3, Reparam), TransformError> {
    Ok(match curve {
        Curve3::Line(line) => (
            Curve3::Line(Line3 {
                origin: rigid.point(line.origin),
                direction: rigid.vector(line.direction),
            }),
            Reparam::Same,
        ),
        Curve3::Circle(circle) => (
            Curve3::Circle(Circle3 {
                frame: rigid.planar_frame(&circle.frame),
                ..*circle
            }),
            Reparam::Same,
        ),
        Curve3::Ellipse(ellipse) => (
            Curve3::Ellipse(Ellipse3 {
                frame: rigid.planar_frame(&ellipse.frame),
                ..*ellipse
            }),
            Reparam::Same,
        ),
        Curve3::Polyline(polyline) => (
            Curve3::Polyline(Polyline {
                points: polyline.points.iter().map(|p| rigid.point(*p)).collect(),
                closed: polyline.closed,
            }),
            Reparam::Same,
        ),
        Curve3::BSpline(spline) => (
            Curve3::BSpline(BSplineCurve {
                control_points: spline
                    .control_points
                    .iter()
                    .map(|p| rigid.point(*p))
                    .collect(),
                ..spline.clone()
            }),
            Reparam::Same,
        ),
        Curve3::Intrinsic(intrinsic) => {
            // A reflection reverses the sense of the screw: the torsion law
            // would have to change sign.
            rigid.unsupported_mirrored("an intrinsic space curve under a reflection")?;
            (
                Curve3::Intrinsic(Intrinsic3 {
                    start: rigid.planar_frame(&intrinsic.start),
                    ..intrinsic.clone()
                }),
                Reparam::Same,
            )
        }
        Curve3::RuledSection(section) => {
            rigid.unsupported_mirrored("a ruled section under a reflection")?;
            (
                Curve3::RuledSection(RuledSection3 {
                    carrier: RuledCarrier {
                        frame: rigid.axial_frame(&section.carrier.frame),
                        ..section.carrier
                    },
                    graph: section.graph,
                }),
                Reparam::Same,
            )
        }
        Curve3::TorusSection(section) => {
            rigid.unsupported_mirrored("a torus section under a reflection")?;
            (
                Curve3::TorusSection(TorusSection3 {
                    torus: TorusCarrier {
                        frame: rigid.axial_frame(&section.torus.frame),
                        ..section.torus
                    },
                    graph: section.graph,
                }),
                Reparam::Same,
            )
        }
        Curve3::ImplicitSection(section) => (
            Curve3::ImplicitSection(ImplicitSection3 {
                carrier: carrier(&section.carrier, rigid)?,
                curve: section.curve.clone(),
            }),
            Reparam::Same,
        ),
        Curve3::PairSection(section) => (
            Curve3::PairSection(PairSection3 {
                first: carrier(&section.first, rigid)?,
                second: carrier(&section.second, rigid)?,
                nodes: section
                    .nodes
                    .iter()
                    .map(|node| PairNode {
                        point: rigid.point(node.point),
                        ..*node
                    })
                    .collect(),
            }),
            Reparam::Same,
        ),
        Curve3::Elevated(_) => {
            return Err(TransformError::Unsupported(
                "an elevated alignment curve (plan and height law)",
            ))
        }
        _ => return Err(TransformError::Unsupported("an unknown curve family")),
    })
}

/// A parameter-plane point reflected in `u`: `(2 pi - u, v)`.
fn reflect_point(p: Point2) -> Point2 {
    Point2::new(TAU - p.x, p.y)
}

/// A parameter-plane direction reflected in `u`.
fn reflect_vector(v: Vec2) -> Vec2 {
    Vec2::new(-v.x, v.y)
}

/// A conic frame reflected in `u`, kept right-handed by reversing `y`; the
/// point at `t` is then at `-t`.
fn reflect_conic_frame(frame: &Frame2) -> Frame2 {
    Frame2 {
        origin: reflect_point(frame.origin),
        x: reflect_vector(frame.x),
        y: -reflect_vector(frame.y),
    }
}

/// The pcurve under the motion. On a face whose angle now reads backwards
/// (`reflect`), it is reflected in `u`; otherwise its points keep their
/// parameters, and only a pcurve holding geometry in space moves.
fn curve2(
    curve: &Curve2,
    rigid: &Rigid,
    reflect: bool,
) -> Result<(Curve2, Reparam), TransformError> {
    if let Curve2::Lifted(lifted) = curve {
        // The pcurve shares its space curve's parameter: it moves with it,
        // and only where nothing is reparameterised.
        rigid.unsupported_mirrored("a lifted pcurve under a reflection")?;
        let (space, _) = curve3(&lifted.curve, rigid)?;
        return Ok((
            Curve2::Lifted(LiftedCurve2 {
                curve: Box::new(space),
                carrier: carrier(&lifted.carrier, rigid)?,
                ..lifted.clone()
            }),
            Reparam::Same,
        ));
    }
    if !reflect {
        return Ok((curve.clone(), Reparam::Same));
    }
    Ok(match curve {
        Curve2::Line(line) => (
            Curve2::Line(Line2 {
                origin: reflect_point(line.origin),
                direction: reflect_vector(line.direction),
            }),
            Reparam::Same,
        ),
        Curve2::Circle(circle) => (
            Curve2::Circle(Circle2 {
                frame: reflect_conic_frame(&circle.frame),
                ..*circle
            }),
            Reparam::Negated,
        ),
        Curve2::Ellipse(ellipse) => (
            Curve2::Ellipse(Ellipse2 {
                frame: reflect_conic_frame(&ellipse.frame),
                ..*ellipse
            }),
            Reparam::Negated,
        ),
        Curve2::Polyline(polyline) => (
            Curve2::Polyline(Polyline {
                points: polyline.points.iter().map(|p| reflect_point(*p)).collect(),
                closed: polyline.closed,
            }),
            Reparam::Same,
        ),
        Curve2::BSpline(spline) => (
            Curve2::BSpline(BSplineCurve {
                control_points: spline
                    .control_points
                    .iter()
                    .map(|p| reflect_point(*p))
                    .collect(),
                ..spline.clone()
            }),
            Reparam::Same,
        ),
        // The point at `t` is `(t, v(t))`; reflected it is `(s, v)` with
        // `s = 2 pi - t`, where `cos` keeps its sign and `sin` changes it.
        Curve2::Sinusoid(wave) => (
            Curve2::Sinusoid(Sinusoid2 {
                sine: -wave.sine,
                ..*wave
            }),
            Reparam::Reflected,
        ),
        _ => {
            return Err(TransformError::Unsupported(
                "a pcurve family on a curved face under a reflection",
            ))
        }
    })
}
