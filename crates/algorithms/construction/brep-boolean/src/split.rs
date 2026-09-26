//! Face splitting: a face cut along its section edges (ADR 0075, step 3).
//!
//! Everything happens in the face's own parameters, where its boundary uses
//! already have pcurves:
//!
//! 1. Each section edge on the face gets an exact pcurve, derived from the
//!    two surfaces: on a plane a line, circle or ellipse maps to its own
//!    family in the plane's coordinates; on a cylinder a ruling is a vertical
//!    line, a circle about the axis a horizontal one, and a plane's oblique
//!    cut a `Sinusoid2` (ADR 0071). Other faces are refused in this stage.
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
use axiolid_evaluate::curve::{derivative2, evaluate2, invert2, invert3};
use axiolid_evaluate::evaluate3;
use axiolid_evaluate::surface::invert;
use axiolid_measure::FaceDomain;
use axiolid_surface::Surface;
use axiolid_topology::{EdgeId, FaceId, Orientation};
use core::f64::consts::{PI, TAU};

use crate::section::SectionEdge;
use crate::BooleanError;

/// Where a piece of a split face's boundary came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PieceSource {
    /// Part of one of the face's own edges.
    Boundary(EdgeId),
    /// Part of the section edge with this index in the list given.
    Section(usize),
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
}

/// Split `face` of `brep` along the section edges lying on it.
///
/// `sections` are the section edges on this face; `others` gives, for each,
/// the support surface of the face on the other operand (needed for exact
/// pcurves on a cylinder).
///
/// # Errors
///
/// A face or section this stage cannot give exact pcurves to, tangent
/// pieces at a vertex, or a boundary that does not close in parameters.
pub fn split_face(
    brep: &ExactBRep,
    face: FaceId,
    sections: &[SectionEdge],
    others: &[Surface],
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
    if !matches!(surface, Surface::Plane(_) | Surface::Cylinder(_)) {
        return Err(BooleanError::UnsupportedSplit);
    }
    let domain = FaceDomain::new(brep, face, tolerance)
        .map_err(BooleanError::Measure)?
        .ok_or(BooleanError::UnsupportedTrim)?;
    let (lo, hi) = domain.bounds();

    // Section pieces with exact pcurves, their starts placed in the face's
    // own parameter range.
    let mut pieces: Vec<(Piece, bool)> = Vec::new();
    let mut ends: Vec<Point3> = Vec::new();
    for (index, (section, other)) in sections.iter().zip(others).enumerate() {
        let piece = section_piece(surface, other, section, index, lo, hi, tolerance)?;
        ends.push(section.start);
        ends.push(section.end);
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
        for piece in uses {
            for part in split_use(surface, piece, &ends, tolerance)? {
                pieces.push((part, false));
            }
        }
    }

    trace(&pieces)
}

/// A section edge with its exact pcurve on `surface`.
///
/// On a seam the start's angle is ambiguous (`0` and `2 pi` name one
/// point), so each placement in the face's range is tried and the one whose
/// piece runs inside the range is kept.
fn section_piece(
    surface: &Surface,
    other: &Surface,
    section: &SectionEdge,
    index: usize,
    lo: Point2,
    hi: Point2,
    tolerance: Tolerance,
) -> Result<Piece, BooleanError> {
    let first = place(surface, section.start, lo, hi, tolerance)?;
    let slack = 1e-9 * (1.0 + lo.x.abs().max(hi.x.abs()));
    let mut candidates = vec![first];
    if matches!(surface, Surface::Cylinder(_)) {
        for shift in [TAU, -TAU] {
            let other_turn = Point2::new(first.x + shift, first.y);
            if other_turn.x >= lo.x - slack && other_turn.x <= hi.x + slack {
                candidates.push(other_turn);
            }
        }
    }
    let mut fallback = None;
    for start_uv in candidates {
        let piece = section_piece_from(surface, other, section, index, start_uv, tolerance)?;
        let mid = evaluate2(&piece.pcurve, 0.5 * (piece.pspan.start + piece.pspan.end))
            .map_err(|_| BooleanError::Evaluation)?;
        if mid.x >= lo.x - slack && mid.x <= hi.x + slack {
            return Ok(piece);
        }
        fallback.get_or_insert(piece);
    }
    fallback.ok_or(BooleanError::Evaluation)
}

/// [`section_piece`] with the start's parameters given.
fn section_piece_from(
    surface: &Surface,
    other: &Surface,
    section: &SectionEdge,
    index: usize,
    start_uv: Point2,
    tolerance: Tolerance,
) -> Result<Piece, BooleanError> {
    let (pcurve, pspan) = match (surface, &section.curve) {
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
        (Surface::Cylinder(c), Curve3::Ellipse(_)) => {
            // A plane's oblique cut: v = mean + a cos u + b sin u, with the
            // angle itself as parameter.
            let Surface::Plane(plane) = other else {
                return Err(BooleanError::UnsupportedSplit);
            };
            let n = plane.frame.z;
            let nz = n.dot(c.frame.z);
            if nz == 0.0 {
                return Err(BooleanError::UnsupportedSplit);
            }
            let wave = Sinusoid2 {
                mean: n.dot(plane.frame.origin - c.frame.origin) / nz,
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
        Ok(invert(surface, p, tolerance)
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
    let (mut u, v) = invert(surface, point, tolerance).map_err(|_| BooleanError::Evaluation)?;
    if matches!(surface, Surface::Cylinder(_)) {
        let slack = 1e-9;
        while u < lo.x - slack {
            u += TAU;
        }
        while u > hi.x + slack {
            u -= TAU;
        }
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
        let Ok(t) = invert3(&piece.curve, point, tolerance) else {
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
        let (u, v) = invert(surface, point, tolerance).map_err(|_| BooleanError::Evaluation)?;
        let mut found = None;
        for shift in [0.0, TAU, -TAU, 2.0 * TAU, -2.0 * TAU] {
            if let Ok(p) = invert2(&piece.pcurve, Point2::new(u + shift, v), tolerance) {
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
        halves.push(Half {
            piece: index,
            reversed: false,
            from: va,
            to: vb,
            leave,
            arrive,
        });
        if *both_ways {
            halves.push(Half {
                piece: index,
                reversed: true,
                from: vb,
                to: va,
                leave: arrive + PI,
                arrive: leave + PI,
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
            // clockwise from there.
            let back = here.arrive + PI;
            let mut best: Option<(Scalar, usize)> = None;
            for (index, candidate) in halves.iter().enumerate() {
                if candidate.from != here.to {
                    continue;
                }
                let is_twin = candidate.piece == here.piece && candidate.reversed != here.reversed;
                let mut turn = (back - candidate.leave).rem_euclid(TAU);
                if is_twin || turn < 1e-12 {
                    turn = TAU;
                }
                if let Some((best_turn, _)) = best {
                    if (turn - best_turn).abs() < 1e-9 && turn < TAU {
                        return Err(BooleanError::TangentSplit);
                    }
                }
                if best.is_none_or(|(best_turn, _)| turn < best_turn) {
                    best = Some((turn, index));
                }
            }
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
        })
        .collect();
    for (hole, polygon) in holes {
        let probe = polygon[0];
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
