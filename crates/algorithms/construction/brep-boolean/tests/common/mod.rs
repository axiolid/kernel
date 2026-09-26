//! Fixtures the constructors do not build yet: a whole sphere, and any
//! solid moved by a rigid motion.

#![allow(dead_code, clippy::needless_range_loop)]

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

/// A box over `[0, 2] x [0, 2]` from `z = 0` up to a bi-quadratic Bezier
/// roof with control heights `h[i][j]` over `x = i`, `y = j`: the roof is
/// `z(x, y) = sum h[i][j] B_i(x / 2) B_j(y / 2)`, its boundary four quadratic
/// B-spline edges shared with the vertical walls.
pub fn spline_box(h: [[f64; 3]; 3]) -> ExactBRep {
    use axiolid_curve::{BSplineCurve2, BSplineCurve3, BSplineSurface, KnotSpec, Line3};
    use axiolid_surface::Plane;
    let roof_point = |i: usize, j: usize| Point3::new(i as f64, j as f64, h[i][j]);
    let roof = BSplineSurface {
        u_degree: 2,
        v_degree: 2,
        control_points: (0..3)
            .map(|i| (0..3).map(|j| roof_point(i, j)).collect())
            .collect(),
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![3, 3],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![3, 3],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::PiecewiseBezier,
        self_intersect: None,
    };
    let curve3 = |pts: Vec<Point3>| {
        Curve3::BSpline(BSplineCurve3 {
            degree: 2,
            control_points: pts,
            knots: vec![0.0, 1.0],
            multiplicities: vec![3, 3],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::PiecewiseBezier,
        })
    };
    let curve2 = |pts: Vec<Vec2>| {
        Curve2::BSpline(BSplineCurve2 {
            degree: 2,
            control_points: pts,
            knots: vec![0.0, 1.0],
            multiplicities: vec![3, 3],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::PiecewiseBezier,
        })
    };
    let mut b = ExactBRepBuilder::default();
    // Vertices: bottom corners 0..4, top corners 4..8, anticlockwise from
    // (0, 0).
    let corners = [(0usize, 0usize), (2, 0), (2, 2), (0, 2)];
    let mut v = Vec::new();
    for (i, j) in corners {
        v.push(b.topology_mut().add_vertex(Vertex {
            position: Point3::new(i as f64, j as f64, 0.0),
        }));
    }
    for (i, j) in corners {
        v.push(b.topology_mut().add_vertex(Vertex {
            position: roof_point(i, j),
        }));
    }
    let line = |b: &mut ExactBRepBuilder, from: Point3, to: Point3| {
        b.add_curve3(Curve3::Line(Line3 {
            origin: from,
            direction: to - from,
        }))
    };
    let pos = |k: usize| {
        let (i, j) = corners[k % 4];
        if k < 4 {
            Point3::new(i as f64, j as f64, 0.0)
        } else {
            roof_point(i, j)
        }
    };
    let edge = |b: &mut ExactBRepBuilder, s: usize, e: usize, curve| {
        let id = b.topology_mut().add_edge(Edge {
            start: v[s],
            end: v[e],
            curve: Some(curve),
        });
        b.set_edge_interval(id, Interval::new(0.0, 1.0));
        id
    };
    // Bottom edges 0->1->2->3->0, verticals k -> k+4, roof edges along the
    // patch boundary: v = 0 (y = 0), u = 1 (x = 2), v = 1 (y = 2), u = 0.
    let mut bottom = Vec::new();
    for k in 0..4 {
        let c = line(&mut b, pos(k), pos((k + 1) % 4));
        bottom.push(edge(&mut b, k, (k + 1) % 4, c));
    }
    let mut vertical = Vec::new();
    for k in 0..4 {
        let c = line(&mut b, pos(k), pos(k + 4));
        vertical.push(edge(&mut b, k, k + 4, c));
    }
    let rows = [
        (0..3).map(|i| roof_point(i, 0)).collect::<Vec<_>>(),
        (0..3).map(|j| roof_point(2, j)).collect(),
        (0..3).map(|i| roof_point(2 - i, 2)).collect(),
        (0..3).map(|j| roof_point(0, 2 - j)).collect(),
    ];
    let mut top = Vec::new();
    for k in 0..4 {
        let c = b.add_curve3(curve3(rows[k].clone()));
        top.push(edge(&mut b, 4 + k, 4 + (k + 1) % 4, c));
    }
    let add_loop =
        |b: &mut ExactBRepBuilder,
         uses: Vec<(axiolid_topology::EdgeId, Orientation, Curve2, Interval)>| {
            let mut edges = Vec::new();
            let mut spans = Vec::new();
            for (edge, orientation, pcurve, span) in uses {
                let pc = b.add_curve2(pcurve);
                edges.push(EdgeUse {
                    edge,
                    orientation,
                    pcurve: Some(pc),
                });
                spans.push(span);
            }
            let id = b.topology_mut().add_loop(Loop { edges });
            for (k, span) in spans.into_iter().enumerate() {
                b.set_pcurve_interval(id, k, span);
            }
            id
        };
    let fwd = Interval::new(0.0, 1.0);
    let rev = Interval::new(1.0, 0.0);
    let seg = |a: Vec2, bb: Vec2| {
        Curve2::Line(Line2 {
            origin: a,
            direction: bb - a,
        })
    };
    let mut faces = Vec::new();
    let mut face = |b: &mut ExactBRepBuilder, surface: Surface, ring, orientation| {
        let s = b.add_surface(surface);
        let f = b.topology_mut().add_face(Face {
            surface: Some(s),
            bounds: vec![FaceBound {
                loop_id: ring,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation,
        });
        faces.push((f, Orientation::Forward));
    };
    // Bottom: plane z = 0 in world (x, y); loop anticlockwise, facing down.
    let world = Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    };
    let c2 = |k: usize| {
        let (i, j) = corners[k];
        Vec2::new(i as f64, j as f64)
    };
    let ring = add_loop(
        &mut b,
        (0..4)
            .map(|k| {
                (
                    bottom[k],
                    Orientation::Forward,
                    seg(c2(k), c2((k + 1) % 4)),
                    fwd,
                )
            })
            .collect(),
    );
    face(
        &mut b,
        Surface::Plane(Plane { frame: world }),
        ring,
        Orientation::Reversed,
    );
    // Roof: the patch, (u, v) in [0, 1]^2, anticlockwise, facing up.
    let uv = [
        Vec2::new(0.0, 0.0),
        Vec2::new(1.0, 0.0),
        Vec2::new(1.0, 1.0),
        Vec2::new(0.0, 1.0),
    ];
    let ring = add_loop(
        &mut b,
        (0..4)
            .map(|k| {
                (
                    top[k],
                    Orientation::Forward,
                    seg(uv[k], uv[(k + 1) % 4]),
                    fwd,
                )
            })
            .collect(),
    );
    face(&mut b, Surface::BSpline(roof), ring, Orientation::Forward);
    // Walls: wall k runs from corner k to corner k + 1, outward. Its plane
    // frame: x along the bottom edge, y up, z = x cross y outward.
    for k in 0..4 {
        let (a, bb) = (pos(k), pos((k + 1) % 4));
        let x = (bb - a).normalize();
        let frame = Frame3 {
            origin: a,
            x,
            y: Vec3::Z,
            z: x.cross(Vec3::Z),
        };
        let local = |p: Point3| Vec2::new((p - a).dot(x), p.z);
        let length = (bb - a).length();
        let top_pts: Vec<Vec2> = rows[k].iter().map(|p| local(*p)).collect();
        let ring = add_loop(
            &mut b,
            vec![
                (
                    bottom[k],
                    Orientation::Forward,
                    seg(Vec2::ZERO, Vec2::new(length, 0.0)),
                    fwd,
                ),
                (
                    vertical[(k + 1) % 4],
                    Orientation::Forward,
                    seg(Vec2::new(length, 0.0), local(pos(4 + (k + 1) % 4))),
                    fwd,
                ),
                (top[k], Orientation::Reversed, curve2(top_pts), rev),
                (
                    vertical[k],
                    Orientation::Reversed,
                    seg(Vec2::ZERO, local(pos(4 + k))),
                    rev,
                ),
            ],
        );
        face(
            &mut b,
            Surface::Plane(Plane { frame }),
            ring,
            Orientation::Forward,
        );
    }
    let shell = b.topology_mut().add_shell(Shell {
        faces: faces.clone(),
        closed: true,
    });
    b.topology_mut().add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    b.finish().expect("a spline-roofed box")
}

/// A solid bounded by `wall` (a surface of revolution about world z whose
/// face winds round it with no seam) and flat discs at the given heights:
/// the wall's loops are circles with no seam edge, as a dome or a can read
/// from a file might be. `rims` are `(height, radius, wall v)` from bottom
/// to top; a single rim closes the wall at a pole above it.
pub fn seamless(wall: Surface, rims: &[(f64, f64, f64)]) -> ExactBRep {
    use axiolid_core::Frame2;
    use axiolid_curve::Circle2;
    use axiolid_surface::Plane;
    let mut b = ExactBRepBuilder::default();
    let mut edges = Vec::new();
    for &(z, r, _) in rims {
        let v = b.topology_mut().add_vertex(Vertex {
            position: Point3::new(r, 0.0, z),
        });
        let c = b.add_curve3(Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: Point3::new(0.0, 0.0, z),
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius: r,
        }));
        let e = b.topology_mut().add_edge(Edge {
            start: v,
            end: v,
            curve: Some(c),
        });
        b.set_edge_interval(e, Interval::new(0.0, TAU));
        edges.push(e);
    }
    let wall_id = b.add_surface(wall);
    // Wall loops: the lowest rim runs +u (the wall above it), the highest -u.
    let mut bounds = Vec::new();
    let n = rims.len();
    for (k, &(_, _, v)) in rims.iter().enumerate() {
        let up = k == 0;
        let pc = b.add_curve2(Curve2::Line(Line2 {
            origin: Vec2::new(0.0, v),
            direction: Vec2::X,
        }));
        let ring = b.topology_mut().add_loop(Loop {
            edges: vec![EdgeUse {
                edge: edges[k],
                orientation: if up {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                },
                pcurve: Some(pc),
            }],
        });
        b.set_pcurve_interval(
            ring,
            0,
            if up {
                Interval::new(0.0, TAU)
            } else {
                Interval::new(TAU, 0.0)
            },
        );
        bounds.push(FaceBound {
            loop_id: ring,
            orientation: Orientation::Forward,
            outer: k == 0,
        });
    }
    let mut faces = vec![(
        b.topology_mut().add_face(Face {
            surface: Some(wall_id),
            bounds,
            orientation: Orientation::Forward,
        }),
        Orientation::Forward,
    )];
    // Discs: the bottom one faces down, the top one (when there are two
    // rims) faces up.
    for (k, &(z, r, _)) in rims.iter().enumerate() {
        let bottom = k == 0;
        if !bottom && k != n - 1 {
            continue;
        }
        let plane = b.add_surface(Surface::Plane(Plane {
            frame: Frame3 {
                origin: Point3::new(0.0, 0.0, z),
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
        }));
        let pc = b.add_curve2(Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Vec2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: r,
        }));
        let ring = b.topology_mut().add_loop(Loop {
            edges: vec![EdgeUse {
                edge: edges[k],
                orientation: if bottom {
                    Orientation::Reversed
                } else {
                    Orientation::Forward
                },
                pcurve: Some(pc),
            }],
        });
        b.set_pcurve_interval(
            ring,
            0,
            if bottom {
                Interval::new(TAU, 0.0)
            } else {
                Interval::new(0.0, TAU)
            },
        );
        faces.push((
            b.topology_mut().add_face(Face {
                surface: Some(plane),
                bounds: vec![FaceBound {
                    loop_id: ring,
                    orientation: Orientation::Forward,
                    outer: true,
                }],
                orientation: Orientation::Forward,
            }),
            Orientation::Forward,
        ));
    }
    let shell = b.topology_mut().add_shell(Shell {
        faces,
        closed: true,
    });
    b.topology_mut().add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    b.finish().expect("a seamless solid")
}
