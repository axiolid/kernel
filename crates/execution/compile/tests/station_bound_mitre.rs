//! Cells beside a mitre (#335, ADR 0082): the certified bound of a run of
//! sections across a seam whose tangents differ, against a dense sampling
//! of the exact mitred solid.
//!
//! The exact solid at a mitre is each side extruded past its end along its
//! end tangent and cut by the mitre plane through the seam point moved by
//! the longitudinal offset along the plane's normal (each side keeps its
//! own half space), plus the planar faces of that plane inside one side's
//! cut but not the other's. The oracle samples its walls, its extrusions
//! and its faces; every test asserts the measured deviation of that solid
//! from the mesh stays within the reported bound.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::{Curve3, Polyline3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, ReferenceMeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, SeamSide, SectionAtStation,
    SolidOperation, Station, StationFrame, StationOffsets, StationOrientation,
};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_reference::station::{station_section3_on, Mitre, SectionFrame};

/// A run compiled, with what its report says of the stationed spine.
struct Compiled {
    mesh: TriMesh,
    bound: DeviationBound,
    meets: bool,
}

impl Compiled {
    fn certified(&self) -> Scalar {
        match self.bound {
            DeviationBound::Certified(value) => value,
            other => panic!("not certified: {other:?}"),
        }
    }
}

fn compile(graph: &GeometryGraph, root: NodeId, chord: Scalar) -> Compiled {
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(chord)
        .unwrap();
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(graph, root, &options)
        .unwrap();
    let bound = report
        .contributions
        .iter()
        .find(|c| c.path == DeviationPath::StationedSpine)
        .expect("the run reports its path")
        .bound;
    Compiled {
        mesh: outcome.mesh,
        bound,
        meets: report.meets_requested(),
    }
}

/// A rectangle `x` wide and one high.
fn rectangle(b: &mut GeometryGraphBuilder, x: Scalar) -> NodeId {
    b.push_value(Profile::Rectangle(RectangleProfile {
        x,
        y: 1.0,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    }))
    .unwrap()
}

/// A spine of rectangles `widths` wide at the two stations `at`, with the
/// longitudinal offset `long` and, when given, an orientation, along
/// `directrix`.
fn spine(
    b: &mut GeometryGraphBuilder,
    directrix: NodeId,
    at: [Scalar; 2],
    widths: [Scalar; 2],
    long: Scalar,
    turn: Option<StationOrientation>,
) -> NodeId {
    let sections = at
        .iter()
        .zip(widths)
        .map(|(&distance, width)| {
            let profile = rectangle(b, width);
            let section = SectionAtStation::new(
                profile,
                Station::new(distance, StationOffsets::new(0.0, 0.0, long)),
            );
            match turn {
                Some(orientation) => section.with_orientation(orientation),
                None => section,
            }
        })
        .collect();
    b.push(GeometryNode::SolidOperation(
        SolidOperation::SectionsAtStations {
            directrix,
            sections,
            frame: StationFrame::Section,
        },
    ))
    .unwrap()
}

/// The exact mitred solid across one corner, sampled against a mesh.
struct Mitred<'a> {
    /// The section frames of the piece ending at the knot and of the one
    /// starting there.
    incoming: &'a dyn Fn(Scalar) -> SectionFrame,
    outgoing: &'a dyn Fn(Scalar) -> SectionFrame,
    knot: Scalar,
    at: [Scalar; 2],
    widths: [Scalar; 2],
    long: Scalar,
    turn: Option<StationOrientation>,
}

/// How far the exact solid strays from a mesh: its walls (with their
/// extrusions to the plane) and its faces in the plane.
#[derive(Debug, Clone, Copy)]
struct Strays {
    walls: Scalar,
    faces: Scalar,
}

impl Strays {
    fn worst(&self) -> Scalar {
        self.walls.max(self.faces)
    }
}

/// The wall edges of a rectangle section, as corner pairs of `(x, y)`.
const CORNERS: [(Scalar, Scalar); 4] = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];

fn lerp(a: Scalar, b: Scalar, u: Scalar) -> Scalar {
    a + (b - a) * u
}

fn nearest(triangles: &[[Point3; 3]], p: Point3) -> Scalar {
    triangles
        .iter()
        .map(|t| {
            axiolid_measure::proximity::closest_point_on_triangle(p, *t)
                .map_or(Scalar::INFINITY, |q| p.distance(q))
        })
        .fold(Scalar::INFINITY, Scalar::min)
}

/// Whether `q` lies inside the polygon `ring` (crossing number).
fn inside(ring: &[(Scalar, Scalar)], q: (Scalar, Scalar)) -> bool {
    let mut in_ring = false;
    for (i, a) in ring.iter().enumerate() {
        let b = ring[(i + 1) % ring.len()];
        if (a.1 > q.1) != (b.1 > q.1) && q.0 < (b.0 - a.0) * (q.1 - a.1) / (b.1 - a.1) + a.0 {
            in_ring = !in_ring;
        }
    }
    in_ring
}

impl Mitred<'_> {
    /// The profile corner `(x, y)` of the section at `s` placed in the
    /// frame `f`: half the width and half the height out, turned.
    fn corner(&self, s: Scalar, (x, y): (Scalar, Scalar), f: &SectionFrame) -> Point3 {
        let u = (s - self.at[0]) / (self.at[1] - self.at[0]);
        let half = 0.5 * lerp(self.widths[0], self.widths[1], u);
        let turned = match self.turn {
            Some(o) => {
                let (axis, reference) = o.unit_axes().unwrap();
                f.oriented(Some(axis), Some(reference)).unwrap()
            }
            None => *f,
        };
        f.point + f.tangent * self.long + turned.lateral * (x * half) + turned.up * (y * 0.5)
    }

    fn measure(&self, mesh: &TriMesh) -> Strays {
        let triangles: Vec<[Point3; 3]> = mesh
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .map(|t| [0, 1, 2].map(|k| mesh.positions[t[k] as usize]))
            .collect();
        let mitre = Mitre::between(&(self.incoming)(self.knot), &(self.outgoing)(self.knot))
            .unwrap()
            .unwrap();
        let n = mitre.normal;
        // The plane the solid is cut at: the seam point moved along n.
        let q0 = mitre.incoming.point + n * self.long;
        let height = |p: Point3| (p - q0).dot(n);
        let (t_in, t_out) = (mitre.incoming.tangent, mitre.outgoing.tangent);
        let (steps, per_edge) = (600, 24);
        let mut walls: Scalar = 0.0;
        for e in 0..4 {
            let (a, c) = (CORNERS[e], CORNERS[(e + 1) % 4]);
            for iv in 0..=per_edge {
                let v = iv as Scalar / per_edge as Scalar;
                let xy = (lerp(a.0, c.0, v), lerp(a.1, c.1, v));
                for j in 0..=steps {
                    let u = j as Scalar / steps as Scalar;
                    // Each wall on its own side of the plane.
                    let s = lerp(self.at[0], self.knot, u);
                    let p = self.corner(s, xy, &(self.incoming)(s));
                    if height(p) <= 0.0 {
                        walls = walls.max(nearest(&triangles, p));
                    }
                    let s = lerp(self.knot, self.at[1], u);
                    let p = self.corner(s, xy, &(self.outgoing)(s));
                    if height(p) >= 0.0 {
                        walls = walls.max(nearest(&triangles, p));
                    }
                }
                // Each side's end section extruded to the plane.
                let c_in = self.corner(self.knot, xy, &(self.incoming)(self.knot));
                if height(c_in) < 0.0 {
                    let reach = -height(c_in) / t_in.dot(n);
                    for j in 0..=40 {
                        let at = c_in + t_in * (reach * j as Scalar / 40.0);
                        walls = walls.max(nearest(&triangles, at));
                    }
                }
                let c_out = self.corner(self.knot, xy, &(self.outgoing)(self.knot));
                if height(c_out) > 0.0 {
                    let reach = height(c_out) / t_out.dot(n);
                    for j in 0..=40 {
                        let at = c_out - t_out * (reach * j as Scalar / 40.0);
                        walls = walls.max(nearest(&triangles, at));
                    }
                }
            }
        }
        // The faces: coordinates in the plane through q0.
        let helper = if n.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
        let e1 = n.cross(helper).normalize();
        let e2 = n.cross(e1);
        let to_plane = |p: Point3| ((p - q0).dot(e1), (p - q0).dot(e2));
        // Where a side's section, cut by the plane, crosses it.
        let cut = |incoming: bool, xy: (Scalar, Scalar)| -> Point3 {
            let (frame, sign, tangent, from, to): (&dyn Fn(Scalar) -> SectionFrame, _, _, _, _) =
                if incoming {
                    (self.incoming, 1.0, t_in, self.at[0], self.knot)
                } else {
                    (self.outgoing, -1.0, t_out, self.knot, self.at[1])
                };
            let at_knot = self.corner(self.knot, xy, &frame(self.knot));
            let h = height(at_knot);
            if h * sign <= 0.0 {
                return at_knot - tangent * (h / tangent.dot(n));
            }
            // The wall itself crosses the plane before the knot.
            let (mut near, mut far) = if incoming { (from, to) } else { (to, from) };
            for _ in 0..100 {
                let mid = 0.5 * (near + far);
                if height(self.corner(mid, xy, &frame(mid))) * sign > 0.0 {
                    far = mid;
                } else {
                    near = mid;
                }
            }
            let s = 0.5 * (near + far);
            self.corner(s, xy, &frame(s))
        };
        let ring = |incoming: bool| -> Vec<(Scalar, Scalar)> {
            let mut out = Vec::new();
            for e in 0..4 {
                let (a, c) = (CORNERS[e], CORNERS[(e + 1) % 4]);
                for i in 0..64 {
                    let v = i as Scalar / 64.0;
                    out.push(to_plane(cut(
                        incoming,
                        (lerp(a.0, c.0, v), lerp(a.1, c.1, v)),
                    )));
                }
            }
            out
        };
        let (ring_in, ring_out) = (ring(true), ring(false));
        let (mut x0, mut x1, mut y0, mut y1) = (Scalar::MAX, Scalar::MIN, Scalar::MAX, Scalar::MIN);
        for p in ring_in.iter().chain(&ring_out) {
            x0 = x0.min(p.0);
            x1 = x1.max(p.0);
            y0 = y0.min(p.1);
            y1 = y1.max(p.1);
        }
        let mut faces: Scalar = 0.0;
        for i in 0..=400 {
            for j in 0..=40 {
                let q = (
                    lerp(x0, x1, i as Scalar / 400.0),
                    lerp(y0, y1, j as Scalar / 40.0),
                );
                if inside(&ring_in, q) != inside(&ring_out, q) {
                    // On the shifted plane, where the faces stand.
                    let p = q0 + e1 * q.0 + e2 * q.1;
                    faces = faces.max(nearest(&triangles, p));
                }
            }
        }
        Strays { walls, faces }
    }
}

/// A polyline of two legs and its frames on each side of the corner.
fn legs(
    points: [Point3; 3],
) -> (
    Curve3,
    impl Fn(Scalar) -> SectionFrame,
    impl Fn(Scalar) -> SectionFrame,
) {
    let curve = Curve3::Polyline(Polyline3 {
        points: points.to_vec(),
        closed: false,
    });
    let (a, b) = (curve.clone(), curve.clone());
    (
        curve,
        move |s| station_section3_on(&a, s, SeamSide::Incoming).unwrap(),
        move |s| station_section3_on(&b, s, SeamSide::Outgoing).unwrap(),
    )
}

/// An L of two legs of ten metres, over a rise of `rise` on the second,
/// with `long` as the longitudinal offset, compiled at a 0.01 chord with
/// stations at 2 and 18; the sampled strays beside it.
fn ell(rise: Scalar, long: Scalar) -> (Compiled, Strays) {
    let points = [
        Point3::ZERO,
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(10.0, 10.0, rise),
    ];
    let (curve, incoming, outgoing) = legs(points);
    let mut b = GeometryGraphBuilder::new();
    let directrix = b.push_value(curve).unwrap();
    let at = [2.0, 18.0];
    let root = spine(&mut b, directrix, at, [2.0, 2.0], long, None);
    let graph = b.finish(vec![root]).unwrap();
    let compiled = compile(&graph, root, 0.01);
    let strays = Mitred {
        incoming: &incoming,
        outgoing: &outgoing,
        knot: 10.0,
        at,
        widths: [2.0, 2.0],
        long,
        turn: None,
    }
    .measure(&compiled.mesh);
    (compiled, strays)
}

fn assert_holds(compiled: &Compiled, strays: &Strays) {
    let bound = compiled.certified();
    assert!(
        strays.worst() <= bound,
        "the exact solid strays {strays:?} from the mesh, past the bound {bound}"
    );
}

#[test]
fn a_turn_on_a_grade_reports_the_gap_between_the_planes_cuts() {
    // The mesh's mitre row is the midpoint of the two sides' projections;
    // on a grade the two sides' cuts differ, and each lies 0.12 m from it.
    let (compiled, strays) = ell(2.0, 0.0);
    assert_holds(&compiled, &strays);
    let bound = compiled.certified();
    assert!((bound - 0.22356947).abs() <= 1e-8, "{bound}");
    assert_eq!(compiled.mesh.positions.len(), 12);
    assert!(!compiled.meets);
    assert!(strays.walls >= 0.1006, "{strays:?}");
}

#[test]
fn a_longitudinal_offset_at_a_mitre_reports_its_shift() {
    // The mesh's mitre row is moved 0.5 along the plane's normal; the
    // exact solid is cut at the moved plane, and the walls end on it.
    let (compiled, strays) = ell(0.0, 0.5);
    assert_holds(&compiled, &strays);
    let bound = compiled.certified();
    assert!((bound - 0.5).abs() <= 1e-8, "{bound}");
    assert_eq!(compiled.mesh.positions.len(), 12);
    assert!((strays.walls - 0.5).abs() <= 1e-3, "{strays:?}");
}

#[test]
fn a_flat_turn_without_an_offset_is_exact() {
    let (compiled, strays) = ell(0.0, 0.0);
    assert_holds(&compiled, &strays);
    assert!(compiled.certified() <= 1e-9, "{}", compiled.certified());
    assert_eq!(compiled.mesh.positions.len(), 12);
    assert!(compiled.meets);
}
