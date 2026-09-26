//! Fixtures the constructors do not build yet: a whole sphere, and any
//! solid moved by a rigid motion.

#![allow(dead_code)]

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_core::{Frame3, Interval, Point3, Vec2, Vec3};
use axiolid_curve::{Circle3, Curve2, Curve3, Line2};
use axiolid_surface::{Sphere, Surface};
use axiolid_topology::{Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

/// A whole sphere: one face, its seam a meridian from the south pole to
/// the north pole used once each way, the poles its only vertices.
pub fn sphere(centre: Point3, radius: f64) -> ExactBRep {
    let frame = Frame3 {
        origin: centre,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    };
    let mut b = ExactBRepBuilder::default();
    let south = b.topology_mut().add_vertex(Vertex {
        position: centre - Vec3::Z * radius,
    });
    let north = b.topology_mut().add_vertex(Vertex {
        position: centre + Vec3::Z * radius,
    });
    let meridian = b.add_curve3(Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: centre,
            x: Vec3::X,
            y: Vec3::Z,
            z: -Vec3::Y,
        },
        radius,
    }));
    let seam = b.topology_mut().add_edge(Edge {
        start: south,
        end: north,
        curve: Some(meridian),
    });
    b.set_edge_interval(seam, Interval::new(-FRAC_PI_2, FRAC_PI_2));
    let right = b.add_curve2(Curve2::Line(Line2 {
        origin: Vec2::new(TAU, -FRAC_PI_2),
        direction: Vec2::Y,
    }));
    let left = b.add_curve2(Curve2::Line(Line2 {
        origin: Vec2::new(0.0, -FRAC_PI_2),
        direction: Vec2::Y,
    }));
    let ring = b.topology_mut().add_loop(Loop {
        edges: vec![
            EdgeUse {
                edge: seam,
                orientation: Orientation::Forward,
                pcurve: Some(right),
            },
            EdgeUse {
                edge: seam,
                orientation: Orientation::Reversed,
                pcurve: Some(left),
            },
        ],
    });
    b.set_pcurve_interval(ring, 0, Interval::new(0.0, PI));
    b.set_pcurve_interval(ring, 1, Interval::new(PI, 0.0));
    let surface = b.add_surface(Surface::Sphere(Sphere { frame, radius }));
    let face = b.topology_mut().add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id: ring,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let shell = b.topology_mut().add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: true,
    });
    b.topology_mut().add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    b.finish().expect("a sphere")
}

/// A rotation about a unit axis by an angle, then a translation.
#[derive(Clone, Copy)]
pub struct Motion {
    pub axis: Vec3,
    pub angle: f64,
    pub shift: Vec3,
}

impl Motion {
    /// A quarter turn about y taking z onto x, exactly: `(x, y, z)` to
    /// `(z, y, -x)`.
    pub const QUARTER_Y: f64 = f64::NAN;

    fn turn(&self, v: Vec3) -> Vec3 {
        if self.angle.is_nan() {
            return Vec3::new(v.z, v.y, -v.x);
        }
        let k = self.axis.normalize();
        let (s, c) = self.angle.sin_cos();
        v * c + k.cross(v) * s + k * (k.dot(v) * (1.0 - c))
    }

    fn point(&self, p: Point3) -> Point3 {
        self.turn(p) + self.shift
    }

    fn frame(&self, f: Frame3) -> Frame3 {
        Frame3 {
            origin: self.point(f.origin),
            x: self.turn(f.x),
            y: self.turn(f.y),
            z: self.turn(f.z),
        }
    }
}

/// `brep` moved rigidly; parameters, pcurves and topology are unchanged.
pub fn moved(brep: &ExactBRep, m: Motion) -> ExactBRep {
    let mut b = ExactBRepBuilder::default();
    for curve in brep.curves3() {
        b.add_curve3(match curve {
            Curve3::Line(l) => Curve3::Line(axiolid_curve::Line3 {
                origin: m.point(l.origin),
                direction: m.turn(l.direction),
            }),
            Curve3::Circle(c) => Curve3::Circle(Circle3 {
                frame: m.frame(c.frame),
                ..*c
            }),
            Curve3::Ellipse(e) => Curve3::Ellipse(axiolid_curve::Ellipse3 {
                frame: m.frame(e.frame),
                ..*e
            }),
            other => panic!("no motion for {other:?}"),
        });
    }
    for curve in brep.curves2() {
        b.add_curve2(curve.clone());
    }
    for surface in brep.surfaces() {
        b.add_surface(match surface {
            Surface::Plane(p) => Surface::Plane(axiolid_surface::Plane {
                frame: m.frame(p.frame),
            }),
            Surface::Cylinder(c) => Surface::Cylinder(axiolid_surface::Cylinder {
                frame: m.frame(c.frame),
                ..*c
            }),
            Surface::Cone(c) => Surface::Cone(axiolid_surface::Cone {
                frame: m.frame(c.frame),
                ..*c
            }),
            Surface::Sphere(s) => Surface::Sphere(Sphere {
                frame: m.frame(s.frame),
                ..*s
            }),
            Surface::Torus(t) => Surface::Torus(axiolid_surface::Torus {
                frame: m.frame(t.frame),
                ..*t
            }),
            other => panic!("no motion for {other:?}"),
        });
    }
    let t = brep.topology();
    for v in t.vertices() {
        b.topology_mut().add_vertex(Vertex {
            position: m.point(v.position),
        });
    }
    for (index, e) in t.edges().iter().enumerate() {
        let id = b.topology_mut().add_edge(e.clone());
        let old = t.edge_id_at(index).expect("edge");
        b.set_edge_interval(id, brep.edge_interval(old).expect("interval"));
    }
    for (index, l) in t.loops().iter().enumerate() {
        let id = b.topology_mut().add_loop(l.clone());
        let old = t.loop_id_at(index).expect("loop");
        for k in 0..l.edges.len() {
            b.set_pcurve_interval(
                id,
                k,
                brep.pcurve_interval(old, k).expect("pcurve interval"),
            );
        }
    }
    for f in t.faces() {
        b.topology_mut().add_face(f.clone());
    }
    for s in t.shells() {
        b.topology_mut().add_shell(s.clone());
    }
    for s in t.solids() {
        b.topology_mut().add_solid(s.clone());
    }
    b.finish().expect("a moved solid")
}
