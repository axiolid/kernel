//! Selection and sewing: the kept regions become the result's faces
//! (ADR 0075, steps 5 and 6).
//!
//! Every kept region becomes a face on its original surface, with its
//! original orientation. A region of the second operand kept by a
//! difference bounds the result from the other side: its loops are
//! reversed, which turns the face without touching its surface.
//!
//! Vertices are welded by position within the linear tolerance; a piece
//! that bounds two kept faces (a section edge, or part of an original edge)
//! becomes one edge used by both, in opposite directions. Faces joined by
//! edges form shells; a shell enclosing positive volume is a solid, one
//! enclosing negative volume a cavity of the solid around it.

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_core::{Interval, Point3, Tolerance};
use axiolid_evaluate::{derivative3, evaluate3};
use axiolid_surface::Surface;
use axiolid_topology::{
    Edge, EdgeId, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex, VertexId,
};

use crate::split::{Piece, Region};
use crate::BooleanError;

/// A kept region with what it needs to become a face.
pub(crate) struct Kept {
    pub(crate) surface: Surface,
    pub(crate) orientation: Orientation,
    pub(crate) region: Region,
    /// Whether the region bounds the result from its other side.
    pub(crate) flip: bool,
}

fn reversed(pieces: &[Piece]) -> Vec<Piece> {
    pieces
        .iter()
        .rev()
        .map(|piece| Piece {
            span: Interval::new(piece.span.end, piece.span.start),
            pspan: Interval::new(piece.pspan.end, piece.pspan.start),
            ..piece.clone()
        })
        .collect()
}

/// The loops of a kept region, outer first, turned if the region flips.
fn loops_of(kept: &Kept) -> Vec<Vec<Piece>> {
    let mut out = vec![kept.region.outer.clone()];
    out.extend(kept.region.holes.iter().cloned());
    if kept.flip {
        out = out.iter().map(|l| reversed(l)).collect();
    }
    out
}

/// Positions welded into vertices, and pieces matched into edges.
struct Welder {
    points: Vec<Point3>,
    /// Per edge: its two vertices, its midpoint, and its curve tangent at
    /// the midpoint in the edge's own direction.
    edges: Vec<(usize, usize, Point3)>,
    eps: f64,
}

impl Welder {
    fn vertex(&mut self, p: Point3) -> usize {
        if let Some(i) = self
            .points
            .iter()
            .position(|q| (*q - p).length() <= self.eps)
        {
            return i;
        }
        self.points.push(p);
        self.points.len() - 1
    }

    /// The edge a piece lies on, and whether the piece runs with it.
    fn edge(&mut self, piece: &Piece) -> Result<(usize, bool), BooleanError> {
        let eval = |t: f64| evaluate3(&piece.curve, t).map_err(|_| BooleanError::Evaluation);
        let (a, b) = (eval(piece.span.start)?, eval(piece.span.end)?);
        let mid = eval(0.5 * (piece.span.start + piece.span.end))?;
        let (va, vb) = (self.vertex(a), self.vertex(b));
        for (index, (s, e, m)) in self.edges.iter().enumerate() {
            if (*m - mid).length() > self.eps {
                continue;
            }
            if (*s, *e) == (va, vb) && va != vb {
                return Ok((index, true));
            }
            if (*s, *e) == (vb, va) && va != vb {
                return Ok((index, false));
            }
            if va == vb && *s == va && *e == vb {
                // A closed piece: compare directions at the midpoint.
                return Ok((index, true));
            }
        }
        self.edges.push((va, vb, mid));
        Ok((self.edges.len() - 1, true))
    }
}

/// Sew kept regions into an exact B-rep.
pub(crate) fn assemble(kept: &[Kept], tolerance: Tolerance) -> Result<ExactBRep, BooleanError> {
    if kept.is_empty() {
        return Err(BooleanError::EmptyResult);
    }
    let eps = tolerance.linear().max(1e-9);

    // Pass 1: which faces share edges, so shells can be told apart.
    let mut welder = Welder {
        points: Vec::new(),
        edges: Vec::new(),
        eps,
    };
    let mut face_edges: Vec<Vec<usize>> = Vec::with_capacity(kept.len());
    for k in kept {
        let mut edges = Vec::new();
        for l in loops_of(k) {
            for piece in &l {
                edges.push(welder.edge(piece)?.0);
            }
        }
        face_edges.push(edges);
    }
    let mut parent: Vec<usize> = (0..kept.len()).collect();
    fn find(parent: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while parent[r] != r {
            r = parent[r];
        }
        let mut y = x;
        while parent[y] != r {
            let next = parent[y];
            parent[y] = r;
            y = next;
        }
        r
    }
    let mut owner: Vec<Option<usize>> = vec![None; welder.edges.len()];
    for (face, edges) in face_edges.iter().enumerate() {
        for &e in edges {
            match owner[e] {
                None => owner[e] = Some(face),
                Some(other) => {
                    let (x, y) = (find(&mut parent, face), find(&mut parent, other));
                    parent[x] = y;
                }
            }
        }
    }
    let mut components: Vec<Vec<usize>> = Vec::new();
    for face in 0..kept.len() {
        let root = find(&mut parent, face);
        if let Some(index) = components
            .iter()
            .position(|c| find(&mut parent, c[0]) == root)
        {
            components[index].push(face);
        } else {
            components.push(vec![face]);
        }
    }

    // Pass 2: each shell on its own, measured to tell solids from cavities.
    let mut shells = Vec::with_capacity(components.len());
    for faces in &components {
        let brep = build(kept, faces, tolerance)?;
        let volume = axiolid_measure::exact_properties(&brep, tolerance)
            .map_err(BooleanError::Measure)?
            .signed_volume;
        shells.push((brep, volume));
    }
    let outers: Vec<usize> = (0..shells.len()).filter(|&i| shells[i].1 > 0.0).collect();
    let voids: Vec<usize> = (0..shells.len()).filter(|&i| shells[i].1 <= 0.0).collect();
    if !voids.is_empty() && outers.len() != 1 {
        return Err(BooleanError::AmbiguousCavity);
    }
    let mut builder = ExactBRepBuilder::default();
    for &outer in &outers {
        let shell = builder.append(&shells[outer].0, false);
        let mut cavities = Vec::new();
        if outers.len() == 1 {
            for &void in &voids {
                cavities.extend(builder.append(&shells[void].0, false));
            }
        }
        builder.topology_mut().add_solid(Solid {
            outer: shell[0],
            voids: cavities,
        });
    }
    builder.finish().map_err(|_| BooleanError::Assembly)
}

/// One shell's faces as an exact B-rep with a single solid.
fn build(kept: &[Kept], faces: &[usize], tolerance: Tolerance) -> Result<ExactBRep, BooleanError> {
    let eps = tolerance.linear().max(1e-9);
    let mut builder = ExactBRepBuilder::default();
    let mut welder = Welder {
        points: Vec::new(),
        edges: Vec::new(),
        eps,
    };
    let mut vertex_ids: Vec<VertexId> = Vec::new();
    // Per edge: its id, curve and increasing span.
    let mut edge_ids: Vec<(EdgeId, axiolid_curve::Curve3)> = Vec::new();
    let mut face_ids = Vec::with_capacity(faces.len());
    for &index in faces {
        let k = &kept[index];
        let surface = builder.add_surface(k.surface.clone());
        let mut bounds = Vec::new();
        for (loop_index, pieces) in loops_of(k).iter().enumerate() {
            let mut uses = Vec::with_capacity(pieces.len());
            let mut intervals = Vec::with_capacity(pieces.len());
            for piece in pieces {
                let before = welder.edges.len();
                let (edge, _) = welder.edge(piece)?;
                while vertex_ids.len() < welder.points.len() {
                    let position = welder.points[vertex_ids.len()];
                    vertex_ids.push(builder.topology_mut().add_vertex(Vertex { position }));
                }
                let ascending = piece.span.end > piece.span.start;
                let along = if edge == before {
                    // A new edge, laid along its increasing span.
                    let (s, e, _) = welder.edges[edge];
                    let (start, end, span) = if ascending {
                        (s, e, piece.span)
                    } else {
                        (e, s, Interval::new(piece.span.end, piece.span.start))
                    };
                    let curve = builder.add_curve3(piece.curve.clone());
                    let id = builder.topology_mut().add_edge(Edge {
                        start: vertex_ids[start],
                        end: vertex_ids[end],
                        curve: Some(curve),
                    });
                    builder.set_edge_interval(id, span);
                    edge_ids.push((id, piece.curve.clone()));
                    ascending
                } else {
                    // An edge already laid: compare directions where the
                    // piece is halfway along.
                    along_edge(&edge_ids[edge].1, piece, tolerance)?
                };
                let pcurve = builder.add_curve2(piece.pcurve.clone());
                uses.push(EdgeUse {
                    edge: edge_ids[edge].0,
                    orientation: if along {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    },
                    pcurve: Some(pcurve),
                });
                intervals.push(piece.pspan);
            }
            let loop_id = builder.topology_mut().add_loop(Loop { edges: uses });
            for (use_index, interval) in intervals.into_iter().enumerate() {
                builder.set_pcurve_interval(loop_id, use_index, interval);
            }
            bounds.push(FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: loop_index == 0,
            });
        }
        face_ids.push((
            builder.topology_mut().add_face(Face {
                surface: Some(surface),
                bounds,
                orientation: k.orientation,
            }),
            Orientation::Forward,
        ));
    }
    let shell = builder.topology_mut().add_shell(Shell {
        faces: face_ids,
        closed: true,
    });
    builder.topology_mut().add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    builder.finish().map_err(|_| BooleanError::Assembly)
}

/// Whether a piece runs the way an edge's curve parameter increases, by
/// the tangents at the piece's midpoint.
fn along_edge(
    edge_curve: &axiolid_curve::Curve3,
    piece: &Piece,
    tolerance: Tolerance,
) -> Result<bool, BooleanError> {
    let t = 0.5 * (piece.span.start + piece.span.end);
    let mid = evaluate3(&piece.curve, t).map_err(|_| BooleanError::Evaluation)?;
    let mut d = derivative3(&piece.curve, t).map_err(|_| BooleanError::Evaluation)?;
    if piece.span.end < piece.span.start {
        d = -d;
    }
    let s = axiolid_evaluate::curve::invert3(edge_curve, mid, tolerance)
        .map_err(|_| BooleanError::Evaluation)?;
    let e = derivative3(edge_curve, s).map_err(|_| BooleanError::Evaluation)?;
    Ok(d.dot(e) > 0.0)
}
