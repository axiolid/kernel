//! Face splitting: a face cut along its section edges (ADR 0075, step 3).
//!
//! Everything happens in the face's own parameters, where its boundary uses
//! already have pcurves:
//!
//! 1. Each section edge on the face gets an exact pcurve, derived from the
//!    two surfaces: on a plane a line, circle or ellipse maps to its own
//!    family in the plane's coordinates; on a cylinder a ruling is a vertical
//!    line, a circle about the axis a horizontal one, and a plane's oblique
//!    cut a `Sinusoid2` (ADR 0071). Every other section on every analytic
//!    face -- a sphere, cone or torus, or a curved section on a plane or
//!    cylinder -- gets the implicit pcurve of ADR 0077: the other surface's
//!    equation read in the face's parameters, traced once per surface over
//!    the face's parameter box, and the stretch between the section's ends
//!    cut out of it.
//! 2. Boundary uses are split wherever a section edge ends on them.
//! 3. The pieces form a graph in `(u, v)`: boundary pieces are walked the
//!    way their loop runs, section pieces both ways. Starting from each
//!    unused half-edge, the walk turns at every vertex to the first
//!    outgoing half-edge clockwise from the one it arrived on, which keeps
//!    one region on its left.
//! 4. Loops running anticlockwise bound regions; clockwise ones are holes,
//!    each given to the smallest region around it.
//!
//! A seam is two vertices in `(u, v)`, so a periodic face splits like any
//! other. Two outgoing pieces leaving a vertex in the same direction
//! (tangent sections) are refused, not ordered by guesswork.

use axiolid_brep::ExactBRep;
use axiolid_core::{Frame2, Interval, Point2, Point3, Scalar, Tolerance, Vec2};
use axiolid_curve::{Curve2, Curve3, Ellipse2, Line2, Sinusoid2};
use axiolid_evaluate::curve::{derivative2, evaluate2, locate2, locate3, second_derivative2};
use axiolid_evaluate::evaluate3;
use axiolid_evaluate::surface::locate;
use axiolid_measure::FaceDomain;
use axiolid_surface::Surface;
use axiolid_topology::{EdgeId, FaceId, Orientation};
use core::f64::consts::{PI, TAU};

use crate::section::SectionEdge;
use crate::support::{periods, window};
use crate::BooleanError;

/// Where a piece of a split face's boundary came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PieceSource {
    /// Part of one of the face's own edges.
    Boundary(EdgeId),
    /// Part of the section edge with this index in the list given.
    Section(usize),
    /// A stretch of parameters that is one point in space -- a sphere's
    /// pole, a cone's apex -- closing a loop in the face's parameters. It
    /// bounds regions while the face is split but is no edge of the result.
    Collapsed,
}

/// One stretch of a split face's boundary, in the direction the region's
/// loop runs.
#[derive(Debug, Clone, PartialEq)]
pub struct Piece {
    /// The 3D curve it lies on.
    pub curve: Curve3,
    /// Its span on `curve`, from where the loop enters it to where it
    /// leaves.
    pub span: Interval,
    /// Its pcurve on the face.
    pub pcurve: Curve2,
    /// Its span on `pcurve`, in the same direction.
    pub pspan: Interval,
    /// Where it came from.
    pub source: PieceSource,
}

/// One region of a split face: an outer loop and its holes.
#[derive(Debug, Clone, PartialEq)]
pub struct Region {
    /// Outer loop, anticlockwise in the face's parameters.
    pub outer: Vec<Piece>,
    /// Holes, clockwise.
    pub holes: Vec<Vec<Piece>>,
    /// Whether the face's own loops wind clockwise in its parameters, so
    /// the region (always anticlockwise) faces opposite to the face's
    /// orientation flag.
    pub against: bool,
}

/// Split `face` of `brep` along the section edges lying on it.
///
/// `sections` are the section edges on this face. `first` says
/// whether `brep` is the first operand of the sections, which picks the
/// side whose `along_*` flag applies: a section running along one of the
/// face's own edges splits nothing but the edge, at its ends. `cuts` are
/// every point where the operand's edges are cut, by sections on any face:
/// an edge is cut the same way in both faces that share it.
///
/// # Errors
///
/// A face or section this stage cannot give exact pcurves to, tangent
/// pieces at a vertex, or a boundary that does not close in parameters.
pub fn split_face(
    brep: &ExactBRep,
    face: FaceId,
    sections: &[SectionEdge],
    first: bool,
    cuts: &[Point3],
    tolerance: Tolerance,
) -> Result<Vec<Region>, BooleanError> {
    let topology = brep.topology();
    let record = topology
        .faces()
        .get(face.index())
        .ok_or(BooleanError::DanglingReference)?;
    let surface = record
        .surface
        .and_then(|id| brep.surfaces().get(id.index()))
        .ok_or(BooleanError::DanglingReference)?;
    let domain = FaceDomain::new(brep, face, tolerance)
        .map_err(BooleanError::Measure)?
        .ok_or(BooleanError::UnsupportedTrim)?;
    let (lo, hi) = domain.bounds();

    // Section pieces with exact pcurves, their starts placed in the face's
    // own parameter range.
    let mut pieces: Vec<(Piece, bool)> = Vec::new();
    let mut ends: Vec<Point3> = cuts.to_vec();
    // One stretch of curve can reach a face from several face pairs (a
    // shared patch's edge is also where the neighbouring faces meet it).
    let mut seen: Vec<&SectionEdge> = Vec::new();
    let mut traces = Traces::default();
    for (index, section) in sections.iter().enumerate() {
        ends.push(section.start);
        ends.push(section.end);
        let (along, other) = if first {
            (section.along_a, &section.other_a)
        } else {
            (section.along_b, &section.other_b)
        };
        if along
            || seen
                .iter()
                .any(|other| same_stretch(other, section, tolerance))
        {
            continue;
        }
        seen.push(section);
        let piece = section_piece(
            surface,
            other,
            section,
            index,
            lo,
            hi,
            &mut traces,
            tolerance,
        )?;
        pieces.push((piece, true));
    }

    // Boundary uses, split where a section ends on them.
    for bound in &record.bounds {
        let wire = topology
            .loops()
            .get(bound.loop_id.index())
            .ok_or(BooleanError::DanglingReference)?;
        let mut uses: Vec<Piece> = Vec::with_capacity(wire.edges.len());
        for (index, use_) in wire.edges.iter().enumerate() {
            let edge = &topology.edges()[use_.edge.index()];
            let curve = edge
                .curve
                .and_then(|id| brep.curves3().get(id.index()))
                .ok_or(BooleanError::DanglingReference)?
                .clone();
            let span = brep
                .edge_interval(use_.edge)
                .ok_or(BooleanError::DanglingReference)?;
            let span = match use_.orientation {
                Orientation::Forward => span,
                Orientation::Reversed => Interval::new(span.end, span.start),
            };
            let pcurve = use_
                .pcurve
                .and_then(|id| brep.curves2().get(id.index()))
                .ok_or(BooleanError::DanglingReference)?
                .clone();
            let pspan = brep
                .pcurve_interval(bound.loop_id, index)
                .ok_or(BooleanError::DanglingReference)?;
            uses.push(Piece {
                curve,
                span,
                pcurve,
                pspan,
                source: PieceSource::Boundary(use_.edge),
            });
        }
        if bound.orientation == Orientation::Reversed {
            uses.reverse();
            for piece in &mut uses {
                piece.span = Interval::new(piece.span.end, piece.span.start);
                piece.pspan = Interval::new(piece.pspan.end, piece.pspan.start);
            }
        }
        let closed = close_poles(surface, uses, tolerance)?;
        for piece in closed {
            for part in split_use(surface, piece, &ends, tolerance)? {
                pieces.push((part, false));
            }
        }
    }

    // Sections ending at a pole meet the collapsed piece there: split it at
    // those points in parameters, so the graph has a vertex for them.
    let mut pole_ends: Vec<Point2> = Vec::new();
    for (piece, section) in &pieces {
        if *section {
            for t in [piece.pspan.start, piece.pspan.end] {
                pole_ends.push(evaluate2(&piece.pcurve, t).map_err(|_| BooleanError::Evaluation)?);
            }
        }
    }
    let mut split_pieces = Vec::with_capacity(pieces.len());
    for (piece, section) in pieces {
        if piece.source != PieceSource::Collapsed {
            split_pieces.push((piece, section));
            continue;
        }
        let (a, b) = (
            evaluate2(&piece.pcurve, piece.pspan.start).map_err(|_| BooleanError::Evaluation)?,
            evaluate2(&piece.pcurve, piece.pspan.end).map_err(|_| BooleanError::Evaluation)?,
        );
        let slack = 1e-9 * (1.0 + a.x.abs().max(b.x.abs()));
        let mut cuts: Vec<Scalar> = pole_ends
            .iter()
            .filter(|p| (p.y - a.y).abs() <= slack)
            .filter(|p| p.x > a.x.min(b.x) + slack && p.x < a.x.max(b.x) - slack)
            .map(|p| (p.x - a.x) / (b.x - a.x))
            .collect();
        cuts.sort_by(Scalar::total_cmp);
        cuts.dedup_by(|x, y| (*x - *y).abs() <= 1e-12);
        let mut from = piece.pspan.start;
        for c in cuts.into_iter().chain(std::iter::once(piece.pspan.end)) {
            split_pieces.push((
                Piece {
                    pspan: Interval::new(from, c),
                    ..piece.clone()
                },
                false,
            ));
            from = c;
        }
    }
    let mut pieces = split_pieces;

    // A face whose loops wind clockwise (a region an earlier boolean turned
    // over) is traced with its boundary reversed, so regions come out
    // anticlockwise; they then face against the face's flag.
    let mut swept = 0.0;
    for (piece, section) in &pieces {
        if !section {
            swept += sweep(piece)?;
        }
    }
    let against = swept < 0.0;
    if against {
        for (piece, section) in &mut pieces {
            if !*section {
                *piece = directed(piece, true);
            }
        }
    }
    let mut regions = trace(&pieces)?;
    for region in &mut regions {
        region.against = against;
    }
    Ok(regions)
}

/// `1/2 int (u dv - v du)` along a piece: summed over a face's loops, its
/// sign is the loops' winding, even where a seam leaves them open.
fn sweep(piece: &Piece) -> Result<Scalar, BooleanError> {
    let n = 64;
    let mut total = 0.0;
    let mut previous =
        evaluate2(&piece.pcurve, piece.pspan.start).map_err(|_| BooleanError::Evaluation)?;
    for i in 1..=n {
        let p =
            piece.pspan.start + (piece.pspan.end - piece.pspan.start) * i as Scalar / n as Scalar;
        let q = evaluate2(&piece.pcurve, p).map_err(|_| BooleanError::Evaluation)?;
        total += 0.5 * (previous.x * q.y - q.x * previous.y);
        previous = q;
    }
    Ok(total)
}

/// A loop's uses with a [`PieceSource::Collapsed`] piece wherever one use
/// ends and the next starts at the same point in space but not in
/// parameters, at a pole: a sphere's seam reaches its north pole at
/// `(2 pi, pi/2)` and leaves it at `(0, pi/2)`.
fn close_poles(
    surface: &Surface,
    uses: Vec<Piece>,
    tolerance: Tolerance,
) -> Result<Vec<Piece>, BooleanError> {
    let n = uses.len();
    let mut out = Vec::with_capacity(n + 2);
    for i in 0..n {
        let next = &uses[(i + 1) % n];
        let a =
            evaluate2(&uses[i].pcurve, uses[i].pspan.end).map_err(|_| BooleanError::Evaluation)?;
        let b = evaluate2(&next.pcurve, next.pspan.start).map_err(|_| BooleanError::Evaluation)?;
        out.push(uses[i].clone());
        let slack = 1e-7 * (1.0 + a.x.abs().max(a.y.abs()));
        if (a - b).length() <= slack {
            continue;
        }
        let pa = axiolid_evaluate::surface::evaluate(surface, a.x, a.y)
            .map_err(|_| BooleanError::Evaluation)?;
        let pb = axiolid_evaluate::surface::evaluate(surface, b.x, b.y)
            .map_err(|_| BooleanError::Evaluation)?;
        let (su, sv) = axiolid_evaluate::surface::partials(surface, a.x, a.y)
            .map_err(|_| BooleanError::Evaluation)?;
        let scale = 1.0 + sv.length();
        let pole =
            (pa - pb).length() <= tolerance.linear().max(1e-9) && su.length() <= 1e-9 * scale;
        if !pole {
            // A gap that is not a pole: the loop winds round a seam with no
            // seam edge, which this split does not close.
            return Err(BooleanError::UnclosedSplit);
        }
        out.push(Piece {
            curve: Curve3::Line(axiolid_curve::Line3 {
                origin: pa,
                direction: axiolid_core::Vec3::ZERO,
            }),
            span: Interval::new(0.0, 1.0),
            pcurve: Curve2::Line(Line2 {
                origin: a,
                direction: b - a,
            }),
            pspan: Interval::new(0.0, 1.0),
            source: PieceSource::Collapsed,
        });
    }
    Ok(out)
}

/// Whether two section edges are the same stretch of curve, either way.
fn same_stretch(a: &SectionEdge, b: &SectionEdge, tolerance: Tolerance) -> bool {
    let eps = tolerance.linear().max(1e-9);
    let near = |p: Point3, q: Point3| (p - q).length() <= eps;
    let mid = |e: &SectionEdge| evaluate3(&e.curve, 0.5 * (e.span.start + e.span.end));
    let ends = (near(a.start, b.start) && near(a.end, b.end))
        || (near(a.start, b.end) && near(a.end, b.start));
    ends && matches!((mid(a), mid(b)), (Ok(p), Ok(q)) if near(p, q))
}

/// Traced sections of the face's surface, one set per other surface, or
/// the fact that the trace could not be done (so it is not tried again).
#[derive(Default)]
struct Traces {
    done: Vec<(Surface, Option<Vec<axiolid_curve::ImplicitCurve2>>)>,
}

impl Traces {
    fn of(
        &mut self,
        surface: &Surface,
        other: &Surface,
        lo: Point2,
        hi: Point2,
    ) -> Result<&[axiolid_curve::ImplicitCurve2], BooleanError> {
        let index = match self.done.iter().position(|(s, _)| s == other) {
            Some(index) => index,
            None => {
                let curves =
                    axiolid_nurbs::trace_section_pcurves(surface, other, window(surface, lo, hi))
                        .ok();
                self.done.push((other.clone(), curves));
                self.done.len() - 1
            }
        };
        self.done[index]
            .1
            .as_deref()
            .ok_or(BooleanError::UnsupportedSplit)
    }
}

/// A section edge with its exact pcurve on `surface`.
///
/// On a seam the start's angle is ambiguous (`0` and `2 pi` name one
/// point), so each placement in the face's range is tried and the one whose
/// piece runs inside the range is kept.
#[allow(clippy::too_many_arguments)]
fn section_piece(
    surface: &Surface,
    other: &Surface,
    section: &SectionEdge,
    index: usize,
    lo: Point2,
    hi: Point2,
    traces: &mut Traces,
    tolerance: Tolerance,
) -> Result<Piece, BooleanError> {
    let closed_form = iso_curve(surface, &section.curve, tolerance)
        || matches!(
            (surface, &section.curve),
            (
                Surface::Plane(_),
                Curve3::Line(_) | Curve3::Circle(_) | Curve3::Ellipse(_)
            ) | (
                Surface::Cylinder(_),
                Curve3::Line(_) | Curve3::Circle(_) | Curve3::Ellipse(_)
            )
        );
    if !closed_form {
        return implicit_piece(surface, other, section, index, lo, hi, traces, tolerance);
    }
    // A start at a pole has no angle: read the piece a little way in.
    let first = match place(surface, section.start, lo, hi, tolerance) {
        Ok(p) => p,
        Err(_) => {
            let t = section.span.start + 0.01 * (section.span.end - section.span.start);
            let p = evaluate3(&section.curve, t).map_err(|_| BooleanError::Evaluation)?;
            place(surface, p, lo, hi, tolerance)?
        }
    };
    let slack = 1e-9 * (1.0 + lo.x.abs().max(hi.x.abs()));
    let mut candidates = vec![first];
    if periods(surface).0 {
        for shift in [TAU, -TAU] {
            let other_turn = Point2::new(first.x + shift, first.y);
            if other_turn.x >= lo.x - slack && other_turn.x <= hi.x + slack {
                candidates.push(other_turn);
            }
        }
    }
    let mut fallback = None;
    for start_uv in candidates {
        let piece = section_piece_from(surface, section, index, start_uv, tolerance)?;
        let mid = evaluate2(&piece.pcurve, 0.5 * (piece.pspan.start + piece.pspan.end))
            .map_err(|_| BooleanError::Evaluation)?;
        if mid.x >= lo.x - slack && mid.x <= hi.x + slack {
            return Ok(piece);
        }
        fallback.get_or_insert(piece);
    }
    fallback.ok_or(BooleanError::Evaluation)
}

/// A section edge's pcurve on an analytic face as its space curve read in
/// the face's parameters (`Curve2::Lifted`), sharing the edge's parameter;
/// the guide unwraps its angles into the face's range.
fn lifted_piece(
    surface: &Surface,
    section: &SectionEdge,
    index: usize,
    lo: Point2,
    hi: Point2,
    tolerance: Tolerance,
) -> Result<Piece, BooleanError> {
    let carrier = match surface {
        Surface::Plane(p) => axiolid_curve::Carrier::Plane(p.frame),
        Surface::Cylinder(c) => axiolid_curve::Carrier::Ruled(axiolid_curve::RuledCarrier {
            frame: c.frame,
            x_radius: c.radius,
            y_radius: c.radius,
            slope: 0.0,
        }),
        Surface::EllipticalCylinder(c) => {
            axiolid_curve::Carrier::Ruled(axiolid_curve::RuledCarrier {
                frame: c.frame,
                x_radius: c.semi_axis_x,
                y_radius: c.semi_axis_y,
                slope: 0.0,
            })
        }
        Surface::Cone(c) => axiolid_curve::Carrier::Ruled(axiolid_curve::RuledCarrier {
            frame: c.frame,
            x_radius: c.radius,
            y_radius: c.radius,
            slope: c.semi_angle.tan(),
        }),
        Surface::Sphere(s) => axiolid_curve::Carrier::Sphere {
            frame: s.frame,
            radius: s.radius,
        },
        Surface::Torus(t) => axiolid_curve::Carrier::Torus(axiolid_curve::TorusCarrier {
            frame: t.frame,
            major_radius: t.major_radius,
            minor_radius: t.minor_radius,
        }),
        // A B-spline carries a section of two B-splines in its own
        // parameters (ADR 0077).
        Surface::BSpline(b) => axiolid_curve::Carrier::Spline(Box::new(b.clone())),
        _ => return Err(BooleanError::UnsupportedSplit),
    };
    if let Curve3::PairSection(pair) = &section.curve {
        let first = pair.side(&carrier).ok_or(BooleanError::UnsupportedSplit)?;
        let (t0, t1) = (section.span.start, section.span.end);
        let n = 4 * pair.nodes.len();
        let mut guide = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let (a, b, _) = pair
                .solve(t0 + (t1 - t0) * i as Scalar / n as Scalar)
                .ok_or(BooleanError::Evaluation)?;
            guide.push(if first { a } else { b });
        }
        return Ok(Piece {
            curve: section.curve.clone(),
            span: section.span,
            pcurve: Curve2::Lifted(axiolid_curve::LiftedCurve2 {
                curve: Box::new(section.curve.clone()),
                carrier,
                start: t0,
                end: t1,
                guide,
            }),
            pspan: section.span,
            source: PieceSource::Section(index),
        });
    }
    let (pu, pv) = periods(surface);
    let n = 128;
    let (t0, t1) = (section.span.start, section.span.end);
    let mut guide: Vec<Point2> = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let t = t0 + (t1 - t0) * i as Scalar / n as Scalar;
        let p = evaluate3(&section.curve, t).map_err(|_| BooleanError::Evaluation)?;
        let (mut u, mut v) = locate(surface, p, tolerance).map_err(|_| BooleanError::Evaluation)?;
        if let Some(last) = guide.last() {
            if pu {
                u += ((last.x - u) / TAU).round() * TAU;
            }
            if pv {
                v += ((last.y - v) / TAU).round() * TAU;
            }
        }
        guide.push(Point2::new(u, v));
    }
    // Whole turns into the face's range, judged at the middle.
    let middle = guide[n / 2];
    let into = |x: Scalar, a: Scalar, b: Scalar, periodic: bool| {
        if !periodic || (x >= a - 1e-9 && x <= b + 1e-9) {
            0.0
        } else {
            ((0.5 * (a + b) - x) / TAU).round() * TAU
        }
    };
    let shift = Vec2::new(
        into(middle.x, lo.x, hi.x, pu),
        into(middle.y, lo.y, hi.y, pv),
    );
    for p in &mut guide {
        *p += shift;
    }
    Ok(Piece {
        curve: section.curve.clone(),
        span: section.span,
        pcurve: Curve2::Lifted(axiolid_curve::LiftedCurve2 {
            curve: Box::new(section.curve.clone()),
            carrier,
            start: t0,
            end: t1,
            guide,
        }),
        pspan: section.span,
        source: PieceSource::Section(index),
    })
}

/// A section edge's implicit pcurve (ADR 0077): the stretch of the other
/// surface's traced equation, in this face's parameters, from the edge's
/// start through its middle to its end.
#[allow(clippy::too_many_arguments)]
fn implicit_piece(
    surface: &Surface,
    other: &Surface,
    section: &SectionEdge,
    index: usize,
    lo: Point2,
    hi: Point2,
    traces: &mut Traces,
    tolerance: Tolerance,
) -> Result<Piece, BooleanError> {
    let uv = |p: Point3| -> Result<Point2, BooleanError> {
        let (u, v) = locate(surface, p, tolerance).map_err(|_| BooleanError::Evaluation)?;
        Ok(Point2::new(u, v))
    };
    let at = |f: Scalar| {
        evaluate3(
            &section.curve,
            section.span.start + f * (section.span.end - section.span.start),
        )
        .map_err(|_| BooleanError::Evaluation)
    };
    let closed = (section.start - section.end).length() <= tolerance.linear().max(1e-9);
    let curves = match traces.of(surface, other, lo, hi) {
        Ok(curves) => curves,
        // The other surface has no equation to read here (a B-spline), or
        // its trace here cannot be decided (branches nearly touching, read
        // far from the parameters' origin): the section's own space curve,
        // read on this face through its closed-form inverse.
        Err(BooleanError::UnsupportedSplit) => {
            return lifted_piece(surface, section, index, lo, hi, tolerance);
        }
        Err(e) => return Err(e),
    };
    let stretch = axiolid_nurbs::extract_stretch(
        curves,
        periods(surface),
        uv(section.start)?,
        [uv(at(1.0 / 3.0)?)?, uv(at(2.0 / 3.0)?)?],
        uv(section.end)?,
        closed,
    )
    .ok_or(BooleanError::UnsupportedSplit)?;
    // Whole turns into the face's own range, judged at the middle.
    let (pu, pv) = periods(surface);
    let middle = stretch
        .point(0.5 * stretch.end())
        .ok_or(BooleanError::Evaluation)?;
    let into = |x: Scalar, a: Scalar, b: Scalar, periodic: bool| {
        if !periodic || (x >= a - 1e-9 && x <= b + 1e-9) {
            return 0.0;
        }
        let k = ((0.5 * (a + b) - x) / TAU).round();
        k * TAU
    };
    let stretch = stretch.shifted(
        into(middle.x, lo.x, hi.x, pu),
        into(middle.y, lo.y, hi.y, pv),
    );
    let end = stretch.end();
    Ok(Piece {
        curve: section.curve.clone(),
        span: section.span,
        pcurve: Curve2::Implicit(stretch),
        pspan: Interval::new(0.0, end),
        source: PieceSource::Section(index),
    })
}

/// Whether a curve is an iso-curve of the surface -- a sphere's meridian or
/// latitude, a cone's ruling or circle about its axis, a torus's tube or
/// ring circle -- whose pcurve is a straight line, affine in the curve's own
/// parameter.
fn iso_curve(surface: &Surface, curve: &Curve3, tolerance: Tolerance) -> bool {
    let eps = tolerance.linear().max(1e-9);
    let parallel = |a: axiolid_core::Vec3, b: axiolid_core::Vec3| {
        a.normalize().cross(b.normalize()).length() <= 1e-12
    };
    let on_axis =
        |p: Point3, o: Point3, z: axiolid_core::Vec3| (p - o).cross(z.normalize()).length() <= eps;
    match (surface, curve) {
        (Surface::Sphere(sp), Curve3::Circle(c)) => {
            let n = c.frame.x.cross(c.frame.y);
            let meridian = (c.frame.origin - sp.frame.origin).length() <= eps
                && (c.radius - sp.radius).abs() <= eps
                && n.normalize().dot(sp.frame.z.normalize()).abs() <= 1e-12;
            let latitude =
                parallel(n, sp.frame.z) && on_axis(c.frame.origin, sp.frame.origin, sp.frame.z);
            meridian || latitude
        }
        (Surface::Cone(k), Curve3::Line(l)) => {
            let slope = k.semi_angle.tan();
            let apex = k.frame.origin - k.frame.z.normalize() * (k.radius / slope);
            let d = l.direction.normalize();
            let through = (apex - l.origin).cross(d).length() <= eps;
            let axis = k.frame.z.normalize();
            through && (d.dot(axis).abs() - k.semi_angle.cos().abs()).abs() <= 1e-12
        }
        (Surface::Cone(k), Curve3::Circle(c)) => {
            let n = c.frame.x.cross(c.frame.y);
            parallel(n, k.frame.z) && on_axis(c.frame.origin, k.frame.origin, k.frame.z)
        }
        (Surface::Torus(t), Curve3::Circle(c)) => {
            let n = c.frame.x.cross(c.frame.y);
            let z = t.frame.z.normalize();
            let ring = parallel(n, z) && on_axis(c.frame.origin, t.frame.origin, z);
            let d = c.frame.origin - t.frame.origin;
            let tube = n.normalize().dot(z).abs() <= 1e-12
                && d.dot(z).abs() <= eps
                && (d.length() - t.major_radius).abs() <= eps
                && (c.radius - t.minor_radius).abs() <= eps;
            ring || tube
        }
        _ => false,
    }
}

/// The straight pcurve of an iso-curve piece, affine in the curve's
/// parameter: read at two interior points (clear of a pole at an end,
/// where the angle has no value) and extended linearly.
fn affine_pcurve(
    surface: &Surface,
    section: &SectionEdge,
    start_uv: Point2,
    tolerance: Tolerance,
) -> Result<(Curve2, Interval), BooleanError> {
    let (t0, t1) = (section.span.start, section.span.end);
    let (ta, tb) = (t0 + 0.25 * (t1 - t0), t0 + 0.75 * (t1 - t0));
    let (pu, pv) = periods(surface);
    let uv = |t: Scalar, near: Point2| -> Result<Point2, BooleanError> {
        let p = evaluate3(&section.curve, t).map_err(|_| BooleanError::Evaluation)?;
        let (mut u, mut v) = locate(surface, p, tolerance).map_err(|_| BooleanError::Evaluation)?;
        if pu {
            u += ((near.x - u) / TAU).round() * TAU;
        }
        if pv {
            v += ((near.y - v) / TAU).round() * TAU;
        }
        Ok(Point2::new(u, v))
    };
    let a = uv(ta, start_uv)?;
    let b = uv(tb, a)?;
    let direction = (b - a) / (tb - ta);
    Ok((
        Curve2::Line(Line2 {
            origin: a - direction * ta,
            direction,
        }),
        section.span,
    ))
}

/// [`section_piece`] with the start's parameters given.
fn section_piece_from(
    surface: &Surface,
    section: &SectionEdge,
    index: usize,
    start_uv: Point2,
    tolerance: Tolerance,
) -> Result<Piece, BooleanError> {
    let (pcurve, pspan) = match (surface, &section.curve) {
        _ if iso_curve(surface, &section.curve, tolerance) => {
            affine_pcurve(surface, section, start_uv, tolerance)?
        }
        (Surface::Plane(p), curve) => {
            let f = p.frame;
            let local = |q: Point3| Point2::new((q - f.origin).dot(f.x), (q - f.origin).dot(f.y));
            let dir = |d: axiolid_core::Vec3| Vec2::new(d.dot(f.x), d.dot(f.y));
            let pcurve = match curve {
                Curve3::Line(l) => Curve2::Line(Line2 {
                    origin: local(l.origin),
                    direction: dir(l.direction),
                }),
                Curve3::Circle(c) => Curve2::Ellipse(Ellipse2 {
                    frame: Frame2 {
                        origin: local(c.frame.origin),
                        x: dir(c.frame.x),
                        y: dir(c.frame.y),
                    },
                    semi_axis_x: c.radius,
                    semi_axis_y: c.radius,
                }),
                Curve3::Ellipse(e) => Curve2::Ellipse(Ellipse2 {
                    frame: Frame2 {
                        origin: local(e.frame.origin),
                        x: dir(e.frame.x),
                        y: dir(e.frame.y),
                    },
                    semi_axis_x: e.semi_axis_x,
                    semi_axis_y: e.semi_axis_y,
                }),
                _ => return Err(BooleanError::UnsupportedSplit),
            };
            // On a plane the pcurve shares the curve's parameter.
            (pcurve, section.span)
        }
        (Surface::Cylinder(c), Curve3::Line(l)) => {
            // A ruling: constant angle, height linear in the parameter.
            let pcurve = Curve2::Line(Line2 {
                origin: start_uv - Vec2::new(0.0, l.direction.dot(c.frame.z)) * section.span.start,
                direction: Vec2::new(0.0, l.direction.dot(c.frame.z)),
            });
            (pcurve, section.span)
        }
        (Surface::Cylinder(c), Curve3::Circle(circle)) => {
            // About the axis: constant height, angle turning with the
            // parameter one way or the other.
            let turn = circle.frame.x.cross(circle.frame.y).dot(c.frame.z).signum();
            let pcurve = Curve2::Line(Line2 {
                origin: start_uv - Vec2::new(turn, 0.0) * section.span.start,
                direction: Vec2::new(turn, 0.0),
            });
            (pcurve, section.span)
        }
        (Surface::Cylinder(c), Curve3::Ellipse(ellipse)) => {
            // A plane's oblique cut: v = mean + a cos u + b sin u, with the
            // angle itself as parameter. The plane is the ellipse's own.
            let n = ellipse.frame.x.cross(ellipse.frame.y);
            let nz = n.dot(c.frame.z);
            if nz == 0.0 {
                return Err(BooleanError::UnsupportedSplit);
            }
            let wave = Sinusoid2 {
                mean: n.dot(ellipse.frame.origin - c.frame.origin) / nz,
                cosine: -c.radius * n.dot(c.frame.x) / nz,
                sine: -c.radius * n.dot(c.frame.y) / nz,
            };
            let pcurve = Curve2::Sinusoid(wave);
            let span = angle_span(surface, section, start_uv.x, tolerance)?;
            (pcurve, span)
        }
        _ => return Err(BooleanError::UnsupportedSplit),
    };
    Ok(Piece {
        curve: section.curve.clone(),
        span: section.span,
        pcurve,
        pspan,
        source: PieceSource::Section(index),
    })
}

/// The angle span a section covers on a cylinder, unwrapped through its
/// midpoint from `start`.
fn angle_span(
    surface: &Surface,
    section: &SectionEdge,
    start: Scalar,
    tolerance: Tolerance,
) -> Result<Interval, BooleanError> {
    let at = |t: Scalar| -> Result<Scalar, BooleanError> {
        let p = evaluate3(&section.curve, t).map_err(|_| BooleanError::Evaluation)?;
        Ok(locate(surface, p, tolerance)
            .map_err(|_| BooleanError::Evaluation)?
            .0)
    };
    let near = |raw: Scalar, reference: Scalar| raw + TAU * ((reference - raw) / TAU).round();
    let (t0, t1) = (section.span.start, section.span.end);
    let full = (t1 - t0).abs() >= TAU - 1e-12;
    let q1 = near(at(t0 + 0.25 * (t1 - t0))?, start);
    let q2 = near(at(t0 + 0.5 * (t1 - t0))?, q1);
    let q3 = near(at(t0 + 0.75 * (t1 - t0))?, q2);
    let end = if full {
        start + TAU * (q1 - start).signum()
    } else {
        near(at(t1)?, q3)
    };
    Ok(Interval::new(start, end))
}

/// A point's parameters on the surface, the angle shifted by whole turns
/// into the face's own range.
fn place(
    surface: &Surface,
    point: Point3,
    lo: Point2,
    hi: Point2,
    tolerance: Tolerance,
) -> Result<Point2, BooleanError> {
    let (mut u, mut v) = locate(surface, point, tolerance).map_err(|_| BooleanError::Evaluation)?;
    let (pu, pv) = periods(surface);
    let slack = 1e-9;
    let wrap = |x: &mut Scalar, a: Scalar, b: Scalar| {
        while *x < a - slack {
            *x += TAU;
        }
        while *x > b + slack {
            *x -= TAU;
        }
    };
    if pu {
        wrap(&mut u, lo.x, hi.x);
    }
    if pv {
        wrap(&mut v, lo.y, hi.y);
    }
    Ok(Point2::new(u, v))
}

/// Split one boundary use wherever a section ends strictly inside it.
fn split_use(
    surface: &Surface,
    piece: Piece,
    ends: &[Point3],
    tolerance: Tolerance,
) -> Result<Vec<Piece>, BooleanError> {
    let mut cuts: Vec<(Scalar, Scalar)> = Vec::new();
    let (lo, hi) = (
        piece.span.start.min(piece.span.end),
        piece.span.start.max(piece.span.end),
    );
    let slack = 1e-9 * (1.0 + lo.abs().max(hi.abs()));
    for &point in ends {
        let Ok(t) = locate3(&piece.curve, point, tolerance) else {
            continue;
        };
        let t = if matches!(piece.curve, Curve3::Circle(_) | Curve3::Ellipse(_)) {
            [t - TAU, t, t + TAU, t + 2.0 * TAU]
                .into_iter()
                .find(|c| *c > lo + slack && *c < hi - slack)
        } else {
            (t > lo + slack && t < hi - slack).then_some(t)
        };
        let Some(t) = t else { continue };
        let on = evaluate3(&piece.curve, t).map_err(|_| BooleanError::Evaluation)?;
        if (on - point).length() > tolerance.linear().max(1e-9) {
            continue;
        }
        // The same point on the pcurve, through the surface's parameters.
        let (u, v) = locate(surface, point, tolerance).map_err(|_| BooleanError::Evaluation)?;
        let (pu, pv) = periods(surface);
        let turns: &[Scalar] = &[0.0, TAU, -TAU, 2.0 * TAU, -2.0 * TAU];
        let mut shifts = Vec::new();
        for &a in if pu { turns } else { &turns[..1] } {
            for &b in if pv { turns } else { &turns[..1] } {
                shifts.push((a, b));
            }
        }
        let mut found = None;
        for (du, dv) in shifts {
            if let Ok(p) = locate2(&piece.pcurve, Point2::new(u + du, v + dv), tolerance) {
                let (plo, phi) = (
                    piece.pspan.start.min(piece.pspan.end),
                    piece.pspan.start.max(piece.pspan.end),
                );
                let candidates = [p - TAU, p, p + TAU];
                if let Some(p) = candidates
                    .into_iter()
                    .find(|c| *c >= plo - slack && *c <= phi + slack)
                {
                    found = Some(p);
                    break;
                }
            }
        }
        let p = found.ok_or(BooleanError::Evaluation)?;
        if !cuts.iter().any(|(c, _)| (c - t).abs() <= slack) {
            cuts.push((t, p));
        }
    }
    if cuts.is_empty() {
        return Ok(vec![piece]);
    }
    let forward = piece.span.end >= piece.span.start;
    cuts.sort_by(|a, b| {
        let order = a.0.total_cmp(&b.0);
        if forward {
            order
        } else {
            order.reverse()
        }
    });
    let mut out = Vec::with_capacity(cuts.len() + 1);
    let (mut t0, mut p0) = (piece.span.start, piece.pspan.start);
    for (t, p) in cuts {
        out.push(Piece {
            span: Interval::new(t0, t),
            pspan: Interval::new(p0, p),
            ..piece.clone()
        });
        (t0, p0) = (t, p);
    }
    out.push(Piece {
        span: Interval::new(t0, piece.span.end),
        pspan: Interval::new(p0, piece.pspan.end),
        ..piece
    });
    Ok(out)
}

/// A directed use of a piece in the graph.
#[derive(Debug, Clone, Copy)]
struct Half {
    piece: usize,
    reversed: bool,
    from: usize,
    to: usize,
    /// Direction leaving `from`, and direction arriving at `to`.
    leave: Scalar,
    arrive: Scalar,
    /// Signed curvature in the running direction at `from` and at `to`,
    /// which orders pieces leaving a vertex in the same direction.
    bend_leave: Scalar,
    bend_arrive: Scalar,
}

fn directed(piece: &Piece, reversed: bool) -> Piece {
    if reversed {
        Piece {
            span: Interval::new(piece.span.end, piece.span.start),
            pspan: Interval::new(piece.pspan.end, piece.pspan.start),
            ..piece.clone()
        }
    } else {
        piece.clone()
    }
}

fn angle_of(v: Vec2) -> Scalar {
    v.y.atan2(v.x)
}

/// Tangent direction of a piece in parameters at one of its ends, in the
/// piece's own running direction.
fn tangent(piece: &Piece, at_start: bool) -> Result<Vec2, BooleanError> {
    let p = if at_start {
        piece.pspan.start
    } else {
        piece.pspan.end
    };
    let d = derivative2(&piece.pcurve, p).map_err(|_| BooleanError::Evaluation)?;
    let sign = if piece.pspan.end >= piece.pspan.start {
        1.0
    } else {
        -1.0
    };
    Ok(d * sign)
}

/// Signed curvature of a piece in parameters at one of its ends, in the
/// piece's running direction.
fn bend(piece: &Piece, at_start: bool) -> Result<Scalar, BooleanError> {
    let p = if at_start {
        piece.pspan.start
    } else {
        piece.pspan.end
    };
    let d = derivative2(&piece.pcurve, p).map_err(|_| BooleanError::Evaluation)?;
    let dd = second_derivative2(&piece.pcurve, p).map_err(|_| BooleanError::Evaluation)?;
    let speed = d.length();
    if speed == 0.0 {
        return Err(BooleanError::Evaluation);
    }
    // Reversing the parameter keeps d x dd's magnitude but flips its sign.
    let sign = if piece.pspan.end >= piece.pspan.start {
        1.0
    } else {
        -1.0
    };
    Ok(sign * (d.x * dd.y - d.y * dd.x) / (speed * speed * speed))
}

fn trace(pieces: &[(Piece, bool)]) -> Result<Vec<Region>, BooleanError> {
    // Weld piece ends in parameters.
    let mut vertices: Vec<Point2> = Vec::new();
    let mut vertex = |p: Point2| -> usize {
        let slack = 1e-7 * (1.0 + p.x.abs().max(p.y.abs()));
        if let Some(i) = vertices.iter().position(|q| (*q - p).length() <= slack) {
            return i;
        }
        vertices.push(p);
        vertices.len() - 1
    };
    let mut halves = Vec::new();
    for (index, (piece, both_ways)) in pieces.iter().enumerate() {
        let a =
            evaluate2(&piece.pcurve, piece.pspan.start).map_err(|_| BooleanError::Evaluation)?;
        let b = evaluate2(&piece.pcurve, piece.pspan.end).map_err(|_| BooleanError::Evaluation)?;
        let (va, vb) = (vertex(a), vertex(b));
        let leave = angle_of(tangent(piece, true)?);
        let arrive = angle_of(tangent(piece, false)?);
        let (bend_leave, bend_arrive) = (bend(piece, true)?, bend(piece, false)?);
        halves.push(Half {
            piece: index,
            reversed: false,
            from: va,
            to: vb,
            leave,
            arrive,
            bend_leave,
            bend_arrive,
        });
        if *both_ways {
            halves.push(Half {
                piece: index,
                reversed: true,
                from: vb,
                to: va,
                leave: arrive + PI,
                arrive: leave + PI,
                bend_leave: -bend_arrive,
                bend_arrive: -bend_leave,
            });
        }
    }

    // Walk every loop.
    let mut used = vec![false; halves.len()];
    let mut loops: Vec<Vec<usize>> = Vec::new();
    for first in 0..halves.len() {
        if used[first] {
            continue;
        }
        let mut walk = Vec::new();
        let mut current = first;
        loop {
            if used[current] {
                if current == first {
                    break;
                }
                return Err(BooleanError::UnclosedSplit);
            }
            used[current] = true;
            walk.push(current);
            let here = halves[current];
            // Back along the arriving piece, then the first outgoing half
            // clockwise from there. Pieces leaving in one direction are
            // ordered a little way out, by how they bend: a piece bending
            // left of another lies clockwise-first from `back`.
            let back = here.arrive + PI;
            let back_bend = -here.bend_arrive;
            let mut keyed: Vec<((Scalar, Scalar), usize)> = Vec::new();
            for (index, candidate) in halves.iter().enumerate() {
                if candidate.from != here.to {
                    continue;
                }
                let is_twin = candidate.piece == here.piece && candidate.reversed != here.reversed;
                let turn = (back - candidate.leave).rem_euclid(TAU);
                // Second order: the turn a small step out is about
                // `turn + (back_bend - bend) * step`.
                let delta = back_bend - candidate.bend_leave;
                let key = if is_twin {
                    (TAU, Scalar::INFINITY)
                } else if !(1e-9..=TAU - 1e-9).contains(&turn) {
                    if delta.abs() <= 1e-9 * (1.0 + back_bend.abs()) {
                        // Leaving back along the arriving piece itself.
                        return Err(BooleanError::TangentSplit);
                    }
                    if delta > 0.0 {
                        (0.0, delta)
                    } else {
                        (TAU, delta)
                    }
                } else {
                    (turn, delta)
                };
                keyed.push((key, index));
            }
            let same_turn = |x: Scalar, y: Scalar| (x - y).abs() < 1e-9;
            let order = |x: &(Scalar, Scalar), y: &(Scalar, Scalar)| {
                if same_turn(x.0, y.0) {
                    x.1.total_cmp(&y.1)
                } else {
                    x.0.total_cmp(&y.0)
                }
            };
            keyed.sort_by(|x, y| order(&x.0, &y.0));
            if let [(first_key, _), (second_key, _), ..] = keyed.as_slice() {
                let same_bend = (first_key.1 - second_key.1).abs()
                    <= 1e-9 * (1.0 + first_key.1.abs().min(1e12));
                if same_turn(first_key.0, second_key.0) && same_bend && first_key.0 < TAU {
                    return Err(BooleanError::TangentSplit);
                }
            }
            let best = keyed.first().copied();
            current = best.ok_or(BooleanError::UnclosedSplit)?.1;
        }
        loops.push(walk);
    }

    // Loops as pieces, with signed areas.
    let as_pieces = |walk: &[usize]| -> Vec<Piece> {
        walk.iter()
            .map(|&h| directed(&pieces[halves[h].piece].0, halves[h].reversed))
            .collect()
    };
    let mut outers: Vec<(Vec<Piece>, Scalar, Vec<Point2>)> = Vec::new();
    let mut holes: Vec<(Vec<Piece>, Vec<Point2>)> = Vec::new();
    for walk in &loops {
        let loop_pieces = as_pieces(walk);
        let polygon = sample(&loop_pieces)?;
        let area = shoelace(&polygon);
        if area > 0.0 {
            outers.push((loop_pieces, area, polygon));
        } else {
            holes.push((loop_pieces, polygon));
        }
    }
    let mut regions: Vec<Region> = outers
        .iter()
        .map(|(pieces, _, _)| Region {
            outer: pieces.clone(),
            holes: Vec::new(),
            against: false,
        })
        .collect();
    for (hole, polygon) in holes {
        // A point just left of the hole's boundary, on the owner's side (a
        // hole runs clockwise): a point on the boundary itself would also
        // read as inside the region the same curve bounds from within.
        let piece = &hole[0];
        let t = 0.5 * (piece.pspan.start + piece.pspan.end);
        let at = evaluate2(&piece.pcurve, t).map_err(|_| BooleanError::Evaluation)?;
        let mut d = derivative2(&piece.pcurve, t).map_err(|_| BooleanError::Evaluation)?;
        if piece.pspan.end < piece.pspan.start {
            d = -d;
        }
        let size = {
            let (mut lo, mut hi) = (polygon[0], polygon[0]);
            for p in &polygon {
                lo = lo.min(*p);
                hi = hi.max(*p);
            }
            (hi - lo).length()
        };
        let left = Vec2::new(-d.y, d.x).normalize_or_zero();
        let probe = at + left * (1e-6 * size.max(1e-9));
        let owner = outers
            .iter()
            .enumerate()
            .filter(|(_, (_, _, outer))| inside_polygon(outer, probe))
            .min_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
            .map(|(index, _)| index)
            .ok_or(BooleanError::UnclosedSplit)?;
        regions[owner].holes.push(hole);
    }
    Ok(regions)
}

/// Dense points around a loop, for its area sign and containment only.
fn sample(pieces: &[Piece]) -> Result<Vec<Point2>, BooleanError> {
    let mut out = Vec::new();
    for piece in pieces {
        for i in 0..32 {
            let p = piece.pspan.start + (piece.pspan.end - piece.pspan.start) * i as Scalar / 32.0;
            out.push(evaluate2(&piece.pcurve, p).map_err(|_| BooleanError::Evaluation)?);
        }
    }
    Ok(out)
}

fn shoelace(polygon: &[Point2]) -> Scalar {
    let mut area = 0.0;
    for i in 0..polygon.len() {
        let (p, q) = (polygon[i], polygon[(i + 1) % polygon.len()]);
        area += p.x * q.y - q.x * p.y;
    }
    0.5 * area
}

pub(crate) fn inside_polygon(polygon: &[Point2], p: Point2) -> bool {
    let mut inside = false;
    for i in 0..polygon.len() {
        let (a, b) = (polygon[i], polygon[(i + 1) % polygon.len()]);
        if (a.y > p.y) != (b.y > p.y) {
            let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
            if x > p.x {
                inside = !inside;
            }
        }
    }
    inside
}
