//! Cells beside a mitre (#335, ADR 0082): the certified bound of a run of
//! sections across a seam whose tangents differ, against a dense sampling
//! of the exact mitred solid.
//!
//! The exact solid at a mitre is each side extruded past its end along its
//! end tangent and cut by the mitre plane through the seam point moved by
//! the longitudinal offset along the plane's normal (each side keeps its
//! own half space), plus the planar faces of that plane inside one side's
//! cut but not the other's. The oracle samples its walls, its extrusions
//! and its faces, piece by piece between the knots; every test asserts the
//! measured deviation of that solid from the mesh stays within the
//! reported bound, and the value the bound takes.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{
    Circle2, Curve2, Curve3, Elevated3, ElevationLaw, Line2, Polyline2, Polyline3,
};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, ReferenceMeshCompiler};
use axiolid_model::{
    CurveRelation, CurveSegment, GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId,
    OpenProfile, SeamSide, SectionAtStation, SolidOperation, Station, StationFrame, StationOffsets,
    StationOrientation, SurfaceRelation, Transition, TrimSelector, TrimmingPreference,
};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_reference::station::{
    station_section2, station_section2_on, station_section3_on, Mitre, SectionFrame,
};

const PI: Scalar = core::f64::consts::PI;

/// A run compiled, with what its report says of its path.
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

    fn vertices(&self) -> usize {
        self.mesh.positions.len()
    }
}

fn compile_on(graph: &GeometryGraph, root: NodeId, chord: Scalar, path: DeviationPath) -> Compiled {
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(chord)
        .unwrap();
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(graph, root, &options)
        .unwrap();
    let bound = report
        .contributions
        .iter()
        .find(|c| c.path == path)
        .expect("the run reports its path")
        .bound;
    Compiled {
        mesh: outcome.mesh,
        bound,
        meets: report.meets_requested(),
    }
}

fn compile(graph: &GeometryGraph, root: NodeId, chord: Scalar) -> Compiled {
    compile_on(graph, root, chord, DeviationPath::StationedSpine)
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

/// A section frame as a function of the distance along a piece.
type Frame = Box<dyn Fn(Scalar) -> SectionFrame>;

/// One piece of the basis between two knots: its frames over `[from, to]`.
struct Piece {
    frame: Frame,
    from: Scalar,
    to: Scalar,
}

/// A 2D curve's frames, reading `side` of a seam.
fn on2(curve: Curve2, side: SeamSide) -> Frame {
    Box::new(move |s| station_section2_on(&curve, s, side).unwrap())
}

/// A 3D curve's frames, reading `side` of a seam.
fn on3(curve: Curve3, side: SeamSide) -> Frame {
    Box::new(move |s| station_section3_on(&curve, s, side).unwrap())
}

/// A curve's frames, `shift` added to the distance asked for.
fn shifted(curve: Curve2, shift: Scalar) -> Frame {
    Box::new(move |s| station_section2(&curve, s + shift).unwrap())
}

/// The profile of a section for the oracle: `locals` at the two stations
/// (linear between), joined in a ring or as an open polyline.
struct Section {
    locals: [Vec<(Scalar, Scalar)>; 2],
    closed: bool,
}

impl Section {
    /// A rectangle `widths` wide and one high.
    fn rectangle(widths: [Scalar; 2]) -> Self {
        let ring = |w: Scalar| {
            vec![
                (-0.5 * w, -0.5),
                (0.5 * w, -0.5),
                (0.5 * w, 0.5),
                (-0.5 * w, 0.5),
            ]
        };
        Self {
            locals: [ring(widths[0]), ring(widths[1])],
            closed: true,
        }
    }

    fn edges(&self) -> Vec<(usize, usize)> {
        let n = self.locals[0].len();
        let count = if self.closed { n } else { n - 1 };
        (0..count).map(|k| (k, (k + 1) % n)).collect()
    }
}

/// The exact mitred solid along pieces joined at mitres, sampled against
/// a mesh.
struct Mitred {
    pieces: Vec<Piece>,
    at: [Scalar; 2],
    section: Section,
    long: Scalar,
    turn: Option<StationOrientation>,
}

/// How far the exact solid strays from a mesh: its walls (with their
/// extrusions to the planes) and its faces in the planes.
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

/// A mitre plane of the exact solid.
struct Cut {
    /// The seam point moved by the longitudinal offset along the normal.
    origin: Point3,
    normal: Vec3,
    /// The tangents of the piece ending there and of the one starting.
    incoming: Vec3,
    outgoing: Vec3,
}

impl Cut {
    fn height(&self, p: Point3) -> Scalar {
        (p - self.origin).dot(self.normal)
    }
}

impl Mitred {
    /// The profile vertex `k` at `s`.
    fn local(&self, s: Scalar, k: usize) -> (Scalar, Scalar) {
        let u = (s - self.at[0]) / (self.at[1] - self.at[0]);
        let (a, b) = (self.section.locals[0][k], self.section.locals[1][k]);
        (lerp(a.0, b.0, u), lerp(a.1, b.1, u))
    }

    /// The point `v` of the way along the profile edge `(k, l)` at `s`,
    /// placed in the frame `f` (turned, with the longitudinal offset).
    fn point(&self, s: Scalar, (k, l): (usize, usize), v: Scalar, f: &SectionFrame) -> Point3 {
        let (a, b) = (self.local(s, k), self.local(s, l));
        let (x, y) = (lerp(a.0, b.0, v), lerp(a.1, b.1, v));
        let turned = match self.turn {
            Some(o) => {
                let (axis, reference) = o.unit_axes().unwrap();
                f.oriented(Some(axis), Some(reference)).unwrap()
            }
            None => *f,
        };
        f.point + f.tangent * self.long + turned.lateral * x + turned.up * y
    }

    fn cuts(&self) -> Vec<Cut> {
        self.pieces
            .windows(2)
            .map(|pair| {
                let (a, b) = (&pair[0], &pair[1]);
                let mitre = Mitre::between(&(a.frame)(a.to), &(b.frame)(b.from))
                    .unwrap()
                    .unwrap();
                Cut {
                    origin: mitre.incoming.point + mitre.normal * self.long,
                    normal: mitre.normal,
                    incoming: mitre.incoming.tangent,
                    outgoing: mitre.outgoing.tangent,
                }
            })
            .collect()
    }

    /// Where piece `i`'s profile point standing at its end (`ending`) or
    /// its start meets the plane `cut`: the end section carried along its
    /// tangent, or, where the wall crosses the plane first, the wall's
    /// crossing point.
    fn cut_point(
        &self,
        i: usize,
        cut: &Cut,
        ending: bool,
        edge: (usize, usize),
        v: Scalar,
    ) -> Point3 {
        let piece = &self.pieces[i];
        let (knot, tangent, sign) = if ending {
            (piece.to, cut.incoming, 1.0)
        } else {
            (piece.from, cut.outgoing, -1.0)
        };
        let at_knot = self.point(knot, edge, v, &(piece.frame)(knot));
        let h = cut.height(at_knot);
        if h * sign <= 0.0 {
            return at_knot - tangent * (h / tangent.dot(cut.normal));
        }
        // The wall itself crosses the plane before the knot.
        let (mut near, mut far) = if ending {
            (piece.from, piece.to)
        } else {
            (piece.to, piece.from)
        };
        for _ in 0..100 {
            let mid = 0.5 * (near + far);
            if cut.height(self.point(mid, edge, v, &(piece.frame)(mid))) * sign > 0.0 {
                far = mid;
            } else {
                near = mid;
            }
        }
        let s = 0.5 * (near + far);
        self.point(s, edge, v, &(piece.frame)(s))
    }

    /// The largest distance of the exact solid from `mesh`, with `steps`
    /// distances along each piece and `per_edge` points along each edge.
    fn measure(&self, mesh: &TriMesh, steps: usize, per_edge: usize) -> Strays {
        let triangles: Vec<[Point3; 3]> = mesh
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .map(|t| [0, 1, 2].map(|k| mesh.positions[t[k] as usize]))
            .collect();
        let cuts = self.cuts();
        let last = self.pieces.len() - 1;
        let edges = self.section.edges();
        let mut walls: Scalar = 0.0;
        for (i, piece) in self.pieces.iter().enumerate() {
            // The plane ahead (the piece ends on it) and behind (it starts
            // on it).
            let ahead = (i < last).then(|| &cuts[i]);
            let behind = (i > 0).then(|| &cuts[i - 1]);
            for &edge in &edges {
                for iv in 0..=per_edge {
                    let v = iv as Scalar / per_edge as Scalar;
                    for j in 0..=steps {
                        let s = lerp(piece.from, piece.to, j as Scalar / steps as Scalar);
                        let p = self.point(s, edge, v, &(piece.frame)(s));
                        if ahead.is_none_or(|c| c.height(p) <= 0.0)
                            && behind.is_none_or(|c| c.height(p) >= 0.0)
                        {
                            walls = walls.max(nearest(&triangles, p));
                        }
                    }
                    // The end sections extruded to the planes.
                    if let Some(c) = ahead {
                        let end = self.point(piece.to, edge, v, &(piece.frame)(piece.to));
                        if c.height(end) < 0.0 {
                            let reach = -c.height(end) / c.incoming.dot(c.normal);
                            for j in 0..=40 {
                                let p = end + c.incoming * (reach * j as Scalar / 40.0);
                                walls = walls.max(nearest(&triangles, p));
                            }
                        }
                    }
                    if let Some(c) = behind {
                        let start = self.point(piece.from, edge, v, &(piece.frame)(piece.from));
                        if c.height(start) > 0.0 {
                            let reach = c.height(start) / c.outgoing.dot(c.normal);
                            for j in 0..=40 {
                                let p = start - c.outgoing * (reach * j as Scalar / 40.0);
                                walls = walls.max(nearest(&triangles, p));
                            }
                        }
                    }
                }
            }
        }
        let mut faces: Scalar = 0.0;
        if self.section.closed {
            for (m, cut) in cuts.iter().enumerate() {
                faces = faces.max(self.faces(&triangles, m, cut, &edges));
            }
        }
        Strays { walls, faces }
    }

    /// How far the planar faces of the plane `cut` between pieces `m` and
    /// `m + 1` stray from the mesh: the points inside one side's cut ring
    /// but not the other's, on the shifted plane.
    fn faces(
        &self,
        triangles: &[[Point3; 3]],
        m: usize,
        cut: &Cut,
        edges: &[(usize, usize)],
    ) -> Scalar {
        let n = cut.normal;
        let helper = if n.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
        let e1 = n.cross(helper).normalize();
        let e2 = n.cross(e1);
        let to_plane = |p: Point3| ((p - cut.origin).dot(e1), (p - cut.origin).dot(e2));
        let ring = |ending: bool| -> Vec<(Scalar, Scalar)> {
            let i = if ending { m } else { m + 1 };
            let mut out = Vec::new();
            for &edge in edges {
                for q in 0..64 {
                    let v = q as Scalar / 64.0;
                    out.push(to_plane(self.cut_point(i, cut, ending, edge, v)));
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
                    let p = cut.origin + e1 * q.0 + e2 * q.1;
                    faces = faces.max(nearest(triangles, p));
                }
            }
        }
        faces
    }
}

fn assert_holds(compiled: &Compiled, strays: &Strays) {
    let bound = compiled.certified();
    assert!(
        strays.worst() <= bound,
        "the exact solid strays {strays:?} from the mesh, past the bound {bound}"
    );
}

fn assert_value(compiled: &Compiled, expected: Scalar, tolerance: Scalar) {
    let bound = compiled.certified();
    assert!(
        (bound - expected).abs() <= tolerance,
        "{bound} != {expected} (+-{tolerance})"
    );
}

/// An L of two legs of ten metres, over a rise of `rise` on the second,
/// with `long` as the longitudinal offset, compiled at a 0.01 chord with
/// stations at 2 and 18; the sampled strays beside it.
fn ell(rise: Scalar, long: Scalar) -> (Compiled, Strays) {
    let curve = Curve3::Polyline(Polyline3 {
        points: vec![
            Point3::ZERO,
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, rise),
        ],
        closed: false,
    });
    let mut b = GeometryGraphBuilder::new();
    let directrix = b.push_value(curve.clone()).unwrap();
    let at = [2.0, 18.0];
    let root = spine(&mut b, directrix, at, [2.0, 2.0], long, None);
    let graph = b.finish(vec![root]).unwrap();
    let compiled = compile(&graph, root, 0.01);
    let strays = Mitred {
        pieces: vec![
            Piece {
                frame: on3(curve.clone(), SeamSide::Incoming),
                from: at[0],
                to: 10.0,
            },
            Piece {
                frame: on3(curve, SeamSide::Outgoing),
                from: 10.0,
                to: at[1],
            },
        ],
        at,
        section: Section::rectangle([2.0, 2.0]),
        long,
        turn: None,
    }
    .measure(&compiled.mesh, 600, 24);
    (compiled, strays)
}

#[test]
fn a_turn_on_a_grade_reports_the_gap_between_the_planes_cuts() {
    // The mesh's mitre row is the midpoint of the two sides' projections;
    // on a grade the two sides' cuts differ, and each lies 0.12 m from it.
    let (compiled, strays) = ell(2.0, 0.0);
    assert_holds(&compiled, &strays);
    assert_value(&compiled, 0.22356947, 1e-8);
    assert_eq!(compiled.vertices(), 12);
    assert!(!compiled.meets);
    assert!(strays.walls >= 0.1006, "{strays:?}");
}

#[test]
fn a_longitudinal_offset_at_a_mitre_reports_its_shift() {
    // The mesh's mitre row is moved 0.5 along the plane's normal; the
    // exact solid is cut at the moved plane, and the walls end on it.
    let (compiled, strays) = ell(0.0, 0.5);
    assert_holds(&compiled, &strays);
    assert_value(&compiled, 0.5, 1e-8);
    assert_eq!(compiled.vertices(), 12);
    assert!((strays.walls - 0.5).abs() <= 1e-3, "{strays:?}");
}

#[test]
fn a_flat_turn_without_an_offset_is_exact() {
    let (compiled, strays) = ell(0.0, 0.0);
    assert_holds(&compiled, &strays);
    assert!(compiled.certified() <= 1e-9, "{}", compiled.certified());
    assert_eq!(compiled.vertices(), 12);
    assert!(compiled.meets);
}

// --- curved pieces, changing sections ---------------------------------------

fn circle(cx: Scalar, cy: Scalar, radius: Scalar) -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(cx, cy),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius,
    })
}

fn line_along_x() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    })
}

/// A composite directrix of pieces, each trimmed by parameter.
fn composite(b: &mut GeometryGraphBuilder, pieces: &[(Curve2, Scalar, Scalar)]) -> NodeId {
    let segments = pieces
        .iter()
        .map(|(curve, start, end)| {
            let basis = b.push_value(curve.clone()).unwrap();
            let curve = b
                .push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
                    basis,
                    start: vec![TrimSelector::Parameter(*start)],
                    end: vec![TrimSelector::Parameter(*end)],
                    sense_agreement: true,
                    preference: TrimmingPreference::Parameter,
                }))
                .unwrap();
            CurveSegment {
                curve,
                same_sense: true,
                transition: Transition::Continuous,
            }
        })
        .collect();
    b.push(GeometryNode::CurveRelation(CurveRelation::Composite {
        segments,
    }))
    .unwrap()
}

/// Two arcs of one radius joined at a corner.
struct TwoArcs {
    directrix: NodeId,
    first: Curve2,
    second: Curve2,
    /// The radius, the second circle's start angle and the first arc's
    /// length.
    radius: Scalar,
    second_start: Scalar,
    length: Scalar,
}

/// An arc of `sweep` radians turning left from the origin along +X, a
/// corner of `corner` degrees turning left after it, then a second arc of
/// `sweep2`.
fn two_arcs(
    b: &mut GeometryGraphBuilder,
    r: Scalar,
    sweep: Scalar,
    corner: Scalar,
    sweep2: Scalar,
) -> TwoArcs {
    let end = (sweep.sin() * r, r - sweep.cos() * r);
    let corner = corner.to_radians();
    let tangent = (sweep.cos(), sweep.sin());
    let (tx, ty) = (
        tangent.0 * corner.cos() - tangent.1 * corner.sin(),
        tangent.0 * corner.sin() + tangent.1 * corner.cos(),
    );
    let centre = (end.0 - ty * r, end.1 + tx * r);
    let mut second_start = (end.1 - centre.1).atan2(end.0 - centre.0);
    if second_start < 0.0 {
        second_start += 2.0 * PI;
    }
    let (first, second) = (circle(0.0, r, r), circle(centre.0, centre.1, r));
    let directrix = composite(
        b,
        &[
            (first.clone(), 1.5 * PI, 1.5 * PI + sweep),
            (second.clone(), second_start, second_start + sweep2),
        ],
    );
    TwoArcs {
        directrix,
        first,
        second,
        radius: r,
        second_start,
        length: r * sweep,
    }
}

/// The oracle's pieces of the two arcs, over the stations `at`.
fn arc_pieces(arcs: &TwoArcs, at: [Scalar; 2]) -> Vec<Piece> {
    vec![
        Piece {
            frame: shifted(arcs.first.clone(), arcs.radius * 1.5 * PI),
            from: at[0],
            to: arcs.length,
        },
        Piece {
            frame: shifted(
                arcs.second.clone(),
                arcs.radius * arcs.second_start - arcs.length,
            ),
            from: arcs.length,
            to: at[1],
        },
    ]
}

/// Two radius-20 arcs of 0.3 rad with a 10 degree corner, stations 1 and
/// 11, a 2 by 1 section, compiled at a 0.01 chord.
fn two_arcs_run(long: Scalar, turn: Option<StationOrientation>) -> (Compiled, Strays) {
    let mut b = GeometryGraphBuilder::new();
    let arcs = two_arcs(&mut b, 20.0, 0.3, 10.0, 0.3);
    let at = [1.0, 11.0];
    let root = spine(&mut b, arcs.directrix, at, [2.0, 2.0], long, turn);
    let graph = b.finish(vec![root]).unwrap();
    let compiled = compile(&graph, root, 0.01);
    let strays = Mitred {
        pieces: arc_pieces(&arcs, at),
        at,
        section: Section::rectangle([2.0, 2.0]),
        long,
        turn,
    }
    .measure(&compiled.mesh, 400, 12);
    (compiled, strays)
}

#[test]
fn two_arcs_beside_a_small_corner_are_certified_within_the_request() {
    let (compiled, strays) = two_arcs_run(0.0, None);
    assert_holds(&compiled, &strays);
    assert_value(&compiled, 2.5635e-3, 1e-7);
    assert!(compiled.meets);
}

#[test]
fn a_turned_section_across_a_mitre_on_arcs_is_certified() {
    let turn = StationOrientation::new(
        Some(Vec3::new(0.0, 0.15, 1.0)),
        Some(Vec3::new(1.0, 0.0, 0.1)),
    );
    let (compiled, strays) = two_arcs_run(0.0, Some(turn));
    assert_holds(&compiled, &strays);
    assert_value(&compiled, 2.7691e-3, 1e-7);
    assert!(compiled.meets);
}

#[test]
fn a_longitudinal_offset_across_a_mitre_on_arcs_reports_the_gap() {
    let (compiled, strays) = two_arcs_run(0.2, None);
    assert_holds(&compiled, &strays);
    assert_value(&compiled, 0.025014, 1e-6);
    assert!(!compiled.meets);
}

#[test]
fn a_line_then_a_tight_arc_across_a_right_angle_is_certified() {
    let mut b = GeometryGraphBuilder::new();
    let arc = circle(8.0, 0.0, 5.0);
    let directrix = composite(
        &mut b,
        &[(line_along_x(), 0.0, 3.0), (arc.clone(), PI, 1.5 * PI)],
    );
    let at = [1.0, 6.0];
    let root = spine(&mut b, directrix, at, [2.0, 2.0], 0.0, None);
    let graph = b.finish(vec![root]).unwrap();
    let compiled = compile(&graph, root, 0.01);
    let strays = Mitred {
        pieces: vec![
            Piece {
                frame: on2(line_along_x(), SeamSide::Incoming),
                from: at[0],
                to: 3.0,
            },
            Piece {
                frame: shifted(arc, 5.0 * PI - 3.0),
                from: 3.0,
                to: at[1],
            },
        ],
        at,
        section: Section::rectangle([2.0, 2.0]),
        long: 0.0,
        turn: None,
    }
    .measure(&compiled.mesh, 400, 12);
    assert_holds(&compiled, &strays);
    assert_value(&compiled, 0.3486, 5e-5);
    assert_eq!(compiled.vertices(), 28);
    assert!(!compiled.meets);
}

#[test]
fn a_widening_section_across_a_corner_is_certified() {
    let curve = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::ZERO,
            Point2::new(10.0, 0.0),
            Point2::new(10.0, 10.0),
        ],
        closed: false,
    });
    let mut b = GeometryGraphBuilder::new();
    let directrix = b.push_value(curve.clone()).unwrap();
    let (at, widths) = ([2.0, 18.0], [2.0, 3.0]);
    let root = spine(&mut b, directrix, at, widths, 0.0, None);
    let graph = b.finish(vec![root]).unwrap();
    let compiled = compile(&graph, root, 0.01);
    let strays = Mitred {
        pieces: vec![
            Piece {
                frame: on2(curve.clone(), SeamSide::Incoming),
                from: at[0],
                to: 10.0,
            },
            Piece {
                frame: on2(curve, SeamSide::Outgoing),
                from: 10.0,
                to: at[1],
            },
        ],
        at,
        section: Section::rectangle(widths),
        long: 0.0,
        turn: None,
    }
    .measure(&compiled.mesh, 400, 12);
    assert_holds(&compiled, &strays);
    assert_value(&compiled, 0.0570247, 1e-6);
    assert!(!compiled.meets);
}

#[test]
fn a_grade_break_over_a_circular_plan_is_certified() {
    let crest = Curve3::Elevated(Elevated3::new(
        circle(0.0, 0.0, 80.0),
        ElevationLaw::Piecewise {
            breaks: vec![20.0],
            laws: vec![
                ElevationLaw::constant_grade(0.0, 0.0),
                ElevationLaw::constant_grade(0.0, 0.08),
            ],
        },
    ));
    let mut b = GeometryGraphBuilder::new();
    let directrix = b.push_value(crest.clone()).unwrap();
    let at = [0.0, 40.0];
    let root = spine(&mut b, directrix, at, [2.0, 2.0], 0.0, None);
    let graph = b.finish(vec![root]).unwrap();
    let compiled = compile(&graph, root, 0.01);
    let strays = Mitred {
        pieces: vec![
            Piece {
                frame: on3(crest.clone(), SeamSide::Incoming),
                from: at[0],
                to: 20.0,
            },
            Piece {
                frame: on3(crest, SeamSide::Outgoing),
                from: 20.0,
                to: at[1],
            },
        ],
        at,
        section: Section::rectangle([2.0, 2.0]),
        long: 0.0,
        turn: None,
    }
    .measure(&compiled.mesh, 120, 6);
    assert_holds(&compiled, &strays);
    assert_value(&compiled, 4.87e-3, 1e-5);
    assert!(compiled.meets);
    assert_eq!(compiled.vertices(), 140);
}

/// A line of five metres, an arc of radius 500 and `length` metres, a line
/// of four, each corner 10 degrees to the left, stations at 1 and the end:
/// the arc's cells lie between two mitres.
fn line_arc_line(length: Scalar) -> (Compiled, Strays) {
    let r = 500.0;
    let corner = 10.0_f64.to_radians();
    let centre = (5.0 - corner.sin() * r, corner.cos() * r);
    let mut start = (0.0 - centre.1).atan2(5.0 - centre.0);
    if start < 0.0 {
        start += 2.0 * PI;
    }
    let sweep = length / r;
    let end_angle = start + sweep;
    let end = (
        centre.0 + r * end_angle.cos(),
        centre.1 + r * end_angle.sin(),
    );
    let heading = corner + sweep + corner;
    let second_line = Curve2::Line(Line2 {
        origin: Point2::new(end.0, end.1),
        direction: Vec2::new(heading.cos(), heading.sin()),
    });
    let arc = circle(centre.0, centre.1, r);
    let mut b = GeometryGraphBuilder::new();
    let directrix = composite(
        &mut b,
        &[
            (line_along_x(), 0.0, 5.0),
            (arc.clone(), start, end_angle),
            (second_line.clone(), 0.0, 4.0),
        ],
    );
    let knots = [5.0, 5.0 + length];
    let at = [1.0, knots[1] + 4.0];
    let root = spine(&mut b, directrix, at, [2.0, 2.0], 0.0, None);
    let graph = b.finish(vec![root]).unwrap();
    let compiled = compile(&graph, root, 0.01);
    let strays = Mitred {
        pieces: vec![
            Piece {
                frame: on2(line_along_x(), SeamSide::Incoming),
                from: at[0],
                to: knots[0],
            },
            Piece {
                frame: shifted(arc, r * start - knots[0]),
                from: knots[0],
                to: knots[1],
            },
            Piece {
                frame: shifted(second_line, -knots[1]),
                from: knots[1],
                to: at[1],
            },
        ],
        at,
        section: Section::rectangle([2.0, 2.0]),
        long: 0.0,
        turn: None,
    }
    .measure(&compiled.mesh, 200, 8);
    (compiled, strays)
}

#[test]
fn a_short_arc_between_two_mitres_is_refined_until_certified() {
    // The arc's cells lie between two mitres, and are split.
    let (compiled, strays) = line_arc_line(4.0);
    assert_holds(&compiled, &strays);
    assert_value(&compiled, 2.184e-3, 1e-6);
    assert!(compiled.meets);
    assert_eq!(compiled.vertices(), 20);
}

#[test]
fn a_cell_between_two_close_mitres_is_split_though_it_could_be_bounded() {
    // The 0.4 m arc between the corners is one cell of the sampling: named
    // and split, whatever a single mitre's bound of it would say.
    let (compiled, strays) = line_arc_line(0.4);
    assert_holds(&compiled, &strays);
    assert_value(&compiled, 3.7663e-5, 1e-9);
    assert!(compiled.meets);
    assert_eq!(compiled.vertices(), 20);
}

#[test]
fn a_seam_inside_a_cell_beside_a_mitre_is_bounded() {
    // A line, then a gentle tangent-continuous arc of radius 20 km, then
    // a 10 degree corner onto a line: the cell next to the mitre holds the
    // joint of the line and the arc.
    let r = 20_000.0;
    let corner = 10.0_f64.to_radians();
    let sweep = 10.0 / r;
    let centre = (20.0, r);
    let end = (
        centre.0 + r * (1.5 * PI + sweep).cos(),
        centre.1 + r * (1.5 * PI + sweep).sin(),
    );
    let heading = sweep + corner;
    let second_line = Curve2::Line(Line2 {
        origin: Point2::new(end.0, end.1),
        direction: Vec2::new(heading.cos(), heading.sin()),
    });
    let arc = circle(centre.0, centre.1, r);
    let mut b = GeometryGraphBuilder::new();
    let directrix = composite(
        &mut b,
        &[
            (line_along_x(), 0.0, 20.0),
            (arc.clone(), 1.5 * PI, 1.5 * PI + sweep),
            (second_line.clone(), 0.0, 10.0),
        ],
    );
    let at = [2.0, 38.0];
    let root = spine(&mut b, directrix, at, [2.0, 2.0], 0.0, None);
    let graph = b.finish(vec![root]).unwrap();
    let compiled = compile(&graph, root, 0.01);
    let line = line_along_x();
    let before: Frame = Box::new(move |s| {
        if s <= 20.0 {
            station_section2_on(&line, s, SeamSide::Incoming).unwrap()
        } else {
            station_section2(&arc, r * 1.5 * PI + (s - 20.0)).unwrap()
        }
    });
    let strays = Mitred {
        pieces: vec![
            Piece {
                frame: before,
                from: at[0],
                to: 30.0,
            },
            Piece {
                frame: shifted(second_line, -30.0),
                from: 30.0,
                to: at[1],
            },
        ],
        at,
        section: Section::rectangle([2.0, 2.0]),
        long: 0.0,
        turn: None,
    }
    .measure(&compiled.mesh, 400, 12);
    assert_holds(&compiled, &strays);
    assert_value(&compiled, 4.5736e-3, 1e-7);
    assert!(compiled.meets);
    assert_eq!(compiled.vertices(), 12);
}

#[test]
fn a_seam_inside_a_cell_after_a_mitre_is_bounded() {
    // A line, a 10 degree corner onto a line, and from its end a gentle
    // tangent-continuous arc of radius 20 km: the cell next to the mitre
    // holds the joint of the line and the arc.
    let r = 20_000.0;
    let heading = 10.0_f64.to_radians();
    let corner_end = (10.0 + 10.0 * heading.cos(), 10.0 * heading.sin());
    let second_line = Curve2::Line(Line2 {
        origin: Point2::new(10.0, 0.0),
        direction: Vec2::new(heading.cos(), heading.sin()),
    });
    let centre = (
        corner_end.0 - r * heading.sin(),
        corner_end.1 + r * heading.cos(),
    );
    let start = heading - 0.5 * PI + 2.0 * PI;
    let sweep = 18.0 / r;
    let arc = circle(centre.0, centre.1, r);
    let mut b = GeometryGraphBuilder::new();
    let directrix = composite(
        &mut b,
        &[
            (line_along_x(), 0.0, 10.0),
            (second_line.clone(), 0.0, 10.0),
            (arc.clone(), start, start + sweep),
        ],
    );
    let at = [2.0, 38.0];
    let root = spine(&mut b, directrix, at, [2.0, 2.0], 0.0, None);
    let graph = b.finish(vec![root]).unwrap();
    let compiled = compile(&graph, root, 0.01);
    let after: Frame = Box::new(move |s| {
        if s <= 20.0 {
            station_section2(&second_line, s - 10.0).unwrap()
        } else {
            station_section2(&arc, r * start + (s - 20.0)).unwrap()
        }
    });
    let strays = Mitred {
        pieces: vec![
            Piece {
                frame: on2(line_along_x(), SeamSide::Incoming),
                from: at[0],
                to: 10.0,
            },
            Piece {
                frame: after,
                from: 10.0,
                to: at[1],
            },
        ],
        at,
        section: Section::rectangle([2.0, 2.0]),
        long: 0.0,
        turn: None,
    }
    .measure(&compiled.mesh, 400, 12);
    assert_holds(&compiled, &strays);
    assert_value(&compiled, 1.2253e-3, 1e-7);
    assert!(compiled.meets);
    assert_eq!(compiled.vertices(), 16);
}

#[test]
fn a_corner_too_sharp_for_its_cells_is_unbounded_by_name() {
    // Two radius-5 arcs, 120 degrees apart, a 1 by 1 section.
    let mut b = GeometryGraphBuilder::new();
    let arcs = two_arcs(&mut b, 5.0, 0.6, 120.0, 0.9);
    let at = [0.2, arcs.length + 0.85 * 5.0 * 0.9];
    let root = spine(&mut b, arcs.directrix, at, [1.0, 1.0], 0.0, None);
    let graph = b.finish(vec![root]).unwrap();
    let compiled = compile(&graph, root, 0.01);
    match compiled.bound {
        DeviationBound::Unbounded(reason) => {
            assert!(
                reason.contains("does not advance across the mitre plane"),
                "{reason}"
            );
        }
        other => panic!("not unbounded by name: {other:?}"),
    }
}

#[test]
fn a_sheet_across_a_mitre_on_arcs_is_certified() {
    let mut b = GeometryGraphBuilder::new();
    let arcs = two_arcs(&mut b, 20.0, 0.3, 10.0, 0.3);
    let points = vec![
        Point2::new(-1.0, -0.04),
        Point2::ZERO,
        Point2::new(1.0, -0.04),
    ];
    let path = b
        .push_value(Curve2::Polyline(Polyline2 {
            points: points.clone(),
            closed: false,
        }))
        .unwrap();
    let open = b
        .push(GeometryNode::OpenProfile(OpenProfile { path }))
        .unwrap();
    let at = [1.0, 11.0];
    let sections = at
        .iter()
        .map(|&d| SectionAtStation::new(open, Station::new(d, StationOffsets::ZERO)))
        .collect();
    let root = b
        .push(GeometryNode::SurfaceRelation(
            SurfaceRelation::OpenSectionsAtStations {
                directrix: arcs.directrix,
                sections,
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![root]).unwrap();
    let compiled = compile_on(&graph, root, 0.01, DeviationPath::SectionedSurface);
    let local: Vec<(Scalar, Scalar)> = points.iter().map(|p| (p.x, p.y)).collect();
    let strays = Mitred {
        pieces: arc_pieces(&arcs, at),
        at,
        section: Section {
            locals: [local.clone(), local],
            closed: false,
        },
        long: 0.0,
        turn: None,
    }
    .measure(&compiled.mesh, 400, 12);
    assert_holds(&compiled, &strays);
    assert!(compiled.meets, "{:?}", compiled.bound);
}
