//! Selection and sewing: the kept regions become the result's faces
//! (ADR 0075, steps 5 and 6).
//!
//! Every kept region becomes a face on its original surface, its loops
//! anticlockwise in the surface's parameters. A region facing the other way
//! -- one of the second operand kept by a difference, or one whose face
//! wound clockwise -- gets the opposite orientation flag, which turns the
//! face (and the direction its edges count in) without touching its
//! surface or its loops.
//!
//! Vertices are welded by position within the linear tolerance; a piece
//! that bounds two kept faces (a section edge, or part of an original edge)
//! becomes one edge used by both, in opposite directions. Where solids of
//! the result touch along an edge, four or more faces meet there: they are
//! sorted by angle about the edge and paired across the wedges of material,
//! so each solid keeps its own copy of the edge and the result stays
//! manifold. Faces joined by edges form shells; a shell enclosing positive
//! volume is a solid, one enclosing negative volume a cavity of the
//! smallest solid around it.

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_core::{Interval, Point3, Tolerance};
use axiolid_evaluate::surface::locate;
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

/// The loops of a kept region, outer first, without the collapsed pieces
/// that closed them at poles in parameters (no edge lies there).
fn loops_of(kept: &Kept) -> Vec<Vec<Piece>> {
    let real = |pieces: &Vec<Piece>| -> Vec<Piece> {
        pieces
            .iter()
            .filter(|p| p.source != crate::split::PieceSource::Collapsed)
            .cloned()
            .collect()
    };
    let mut out = vec![real(&kept.region.outer)];
    out.extend(kept.region.holes.iter().map(real));
    out
}

/// The orientation flag a kept region's face gets.
fn orientation_of(kept: &Kept) -> Orientation {
    if kept.flip ^ kept.region.against {
        match kept.orientation {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
        }
    } else {
        kept.orientation
    }
}

/// Where one piece of one kept face's loops sits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Use {
    face: usize,
    ring: usize,
    index: usize,
}

/// Every piece of every kept face matched to the edge it becomes.
struct Sewing {
    /// Welded vertex positions.
    points: Vec<Point3>,
    /// Per edge: its first use, which lays it.
    edges: Vec<Use>,
    /// Per face, per loop, per piece: the edge it lies on.
    edge_of: Vec<Vec<Vec<usize>>>,
    /// Per face, per loop, per piece: its welded start and end vertices.
    ends_of: Vec<Vec<Vec<(usize, usize)>>>,
}

fn sew(kept: &[Kept], tolerance: Tolerance) -> Result<Sewing, BooleanError> {
    let eps = tolerance.linear();
    let mut points: Vec<Point3> = Vec::new();
    let mut vertex = |p: Point3| -> usize {
        if let Some(i) = points.iter().position(|q| (*q - p).length() <= eps) {
            return i;
        }
        points.push(p);
        points.len() - 1
    };
    // Group uses by the stretch of curve they cover: both end vertices and
    // the midpoint.
    let mut groups: Vec<((usize, usize, Point3), Vec<Use>)> = Vec::new();
    let mut ends_of = Vec::with_capacity(kept.len());
    for (face, k) in kept.iter().enumerate() {
        let mut rings = Vec::new();
        for (ring, pieces) in loops_of(k).iter().enumerate() {
            let mut ends = Vec::with_capacity(pieces.len());
            for (index, piece) in pieces.iter().enumerate() {
                let eval =
                    |t: f64| evaluate3(&piece.curve, t).map_err(|_| BooleanError::Evaluation);
                let (a, b) = (eval(piece.span.start)?, eval(piece.span.end)?);
                let mid = eval(0.5 * (piece.span.start + piece.span.end))?;
                let (va, vb) = (vertex(a), vertex(b));
                ends.push((va, vb));
                let key = (va.min(vb), va.max(vb));
                let at = Use { face, ring, index };
                match groups
                    .iter_mut()
                    .find(|((x, y, m), _)| (*x, *y) == key && (*m - mid).length() <= eps)
                {
                    Some((_, uses)) => uses.push(at),
                    None => groups.push(((key.0, key.1, mid), vec![at])),
                }
            }
            rings.push(ends);
        }
        ends_of.push(rings);
    }

    let mut edge_of: Vec<Vec<Vec<usize>>> = ends_of
        .iter()
        .map(|rings| rings.iter().map(|r| vec![usize::MAX; r.len()]).collect())
        .collect();
    let mut edges = Vec::new();
    for (_, uses) in &groups {
        let pairs = match uses.len() {
            2 => vec![(uses[0], uses[1])],
            n if n > 2 && n % 2 == 0 => radial_pairs(kept, uses, tolerance)?,
            _ => return Err(BooleanError::Assembly),
        };
        for (x, y) in pairs {
            let id = edges.len();
            edges.push(x);
            edge_of[x.face][x.ring][x.index] = id;
            edge_of[y.face][y.ring][y.index] = id;
        }
    }
    Ok(Sewing {
        points,
        edges,
        edge_of,
        ends_of,
    })
}

/// Pair the faces meeting at one edge across the wedges of material
/// between them, by their angle about the edge.
fn radial_pairs(
    kept: &[Kept],
    uses: &[Use],
    tolerance: Tolerance,
) -> Result<Vec<(Use, Use)>, BooleanError> {
    let piece_of = |u: &Use| loops_of(&kept[u.face])[u.ring][u.index].clone();
    let first = piece_of(&uses[0]);
    let t_mid = 0.5 * (first.span.start + first.span.end);
    let m = evaluate3(&first.curve, t_mid).map_err(|_| BooleanError::Evaluation)?;
    // The edge's axis, from the first use's running direction.
    let axis = {
        let mut d = derivative3(&first.curve, t_mid).map_err(|_| BooleanError::Evaluation)?;
        if first.span.end < first.span.start {
            d = -d;
        }
        d.normalize()
    };
    // Per use: the direction into its face, and whether the solid lies on
    // the side of increasing angle.
    let mut spokes = Vec::with_capacity(uses.len());
    for u in uses {
        let piece = piece_of(u);
        let k = &kept[u.face];
        let t = axiolid_evaluate::curve::locate3(&piece.curve, m, tolerance)
            .map_err(|_| BooleanError::Evaluation)?;
        let mut along = derivative3(&piece.curve, t).map_err(|_| BooleanError::Evaluation)?;
        if piece.span.end < piece.span.start {
            along = -along;
        }
        let (su, sv) = locate(&k.surface, m, tolerance).map_err(|_| BooleanError::Evaluation)?;
        let n_uv = axiolid_evaluate::surface::normal(&k.surface, su, sv)
            .map_err(|_| BooleanError::Evaluation)?;
        // Loops run anticlockwise in parameters: the face lies to the left
        // of each piece, seen from the parameter normal.
        let into = n_uv.cross(along);
        let into = (into - axis * into.dot(axis)).normalize();
        let outward = match orientation_of(k) {
            Orientation::Forward => n_uv,
            Orientation::Reversed => -n_uv,
        };
        let increasing = axis.cross(into);
        spokes.push((into, (-outward).dot(increasing) > 0.0));
    }
    let reference = spokes[0].0;
    let across = axis.cross(reference);
    let mut order: Vec<(f64, usize)> = spokes
        .iter()
        .enumerate()
        .map(|(i, (d, _))| (d.dot(across).atan2(d.dot(reference)), i))
        .collect();
    order.sort_by(|a, b| a.0.total_cmp(&b.0));
    let n = order.len();
    let start = (0..n)
        .find(|&i| spokes[order[i].1].1)
        .ok_or(BooleanError::Assembly)?;
    let mut pairs = Vec::with_capacity(n / 2);
    for step in 0..n / 2 {
        let (i, j) = (
            order[(start + 2 * step) % n].1,
            order[(start + 2 * step + 1) % n].1,
        );
        // The next face bounds the same wedge from the other side.
        if !spokes[i].1 || spokes[j].1 {
            return Err(BooleanError::Assembly);
        }
        pairs.push((uses[i], uses[j]));
    }
    Ok(pairs)
}

/// Sew kept regions into an exact B-rep.
pub(crate) fn assemble(kept: &[Kept], tolerance: Tolerance) -> Result<ExactBRep, BooleanError> {
    if kept.is_empty() {
        return Err(BooleanError::EmptyResult);
    }
    let sewing = sew(kept, tolerance)?;

    // Faces joined by edges form shells.
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
    let mut owner: Vec<Option<usize>> = vec![None; sewing.edges.len()];
    for (face, rings) in sewing.edge_of.iter().enumerate() {
        for &e in rings.iter().flatten() {
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

    // Each shell on its own, measured to tell solids from cavities.
    let mut shells = Vec::with_capacity(components.len());
    for faces in &components {
        let brep = build(kept, &sewing, faces, tolerance)?;
        let volume = axiolid_measure::exact_properties(&brep, tolerance)
            .map_err(BooleanError::Measure)?
            .signed_volume;
        shells.push((brep, volume));
    }
    let outers: Vec<usize> = (0..shells.len()).filter(|&i| shells[i].1 > 0.0).collect();
    let voids: Vec<usize> = (0..shells.len()).filter(|&i| shells[i].1 <= 0.0).collect();
    // Each cavity belongs to the smallest solid around it. A void shell
    // never touches an outer one, so any point on it is strictly inside or
    // strictly outside each solid.
    let mut owners: Vec<Vec<usize>> = vec![Vec::new(); shells.len()];
    if !voids.is_empty() {
        let classifiers = outers
            .iter()
            .map(|&outer| crate::classify::Solid::new(&shells[outer].0, tolerance))
            .collect::<Result<Vec<_>, _>>()?;
        for &void in &voids {
            let face = components[void][0];
            let point =
                crate::classify::interior_points(&kept[face].region, &kept[face].surface)?[0];
            let mut owner: Option<usize> = None;
            for (slot, classifier) in classifiers.iter().enumerate() {
                if classifier.contains(point, tolerance)?
                    && owner.is_none_or(|best| shells[outers[slot]].1 < shells[outers[best]].1)
                {
                    owner = Some(slot);
                }
            }
            let slot = owner.ok_or(BooleanError::AmbiguousCavity)?;
            owners[outers[slot]].push(void);
        }
    }
    let mut builder = ExactBRepBuilder::default();
    for &outer in &outers {
        let shell = builder.append(&shells[outer].0, false);
        let mut cavities = Vec::new();
        for &void in &owners[outer] {
            cavities.extend(builder.append(&shells[void].0, false));
        }
        builder.topology_mut().add_solid(Solid {
            outer: shell[0],
            voids: cavities,
        });
    }
    builder.finish().map_err(|_| BooleanError::Assembly)
}

/// One shell's faces as an exact B-rep with a single solid.
fn build(
    kept: &[Kept],
    sewing: &Sewing,
    faces: &[usize],
    tolerance: Tolerance,
) -> Result<ExactBRep, BooleanError> {
    let mut builder = ExactBRepBuilder::default();
    let mut vertex_ids: Vec<Option<VertexId>> = vec![None; sewing.points.len()];
    let mut edge_ids: Vec<Option<(EdgeId, axiolid_curve::Curve3)>> = vec![None; sewing.edges.len()];
    let mut face_ids = Vec::with_capacity(faces.len());
    for &index in faces {
        let k = &kept[index];
        let surface = builder.add_surface(k.surface.clone());
        let mut bounds = Vec::new();
        for (ring, pieces) in loops_of(k).iter().enumerate() {
            let mut uses = Vec::with_capacity(pieces.len());
            let mut intervals = Vec::with_capacity(pieces.len());
            for (at, piece) in pieces.iter().enumerate() {
                let edge = sewing.edge_of[index][ring][at];
                let (va, vb) = sewing.ends_of[index][ring][at];
                for v in [va, vb] {
                    if vertex_ids[v].is_none() {
                        vertex_ids[v] = Some(builder.topology_mut().add_vertex(Vertex {
                            position: sewing.points[v],
                        }));
                    }
                }
                let along = if let Some((_, curve)) = &edge_ids[edge] {
                    // An edge already laid: compare directions where the
                    // piece is halfway along.
                    along_edge(curve, piece, tolerance)?
                } else {
                    // A new edge, laid along its increasing span.
                    let ascending = piece.span.end > piece.span.start;
                    let (start, end, span) = if ascending {
                        (va, vb, piece.span)
                    } else {
                        (vb, va, Interval::new(piece.span.end, piece.span.start))
                    };
                    let curve = builder.add_curve3(piece.curve.clone());
                    let id = builder.topology_mut().add_edge(Edge {
                        start: vertex_ids[start].ok_or(BooleanError::Assembly)?,
                        end: vertex_ids[end].ok_or(BooleanError::Assembly)?,
                        curve: Some(curve),
                    });
                    builder.set_edge_interval(id, span);
                    edge_ids[edge] = Some((id, piece.curve.clone()));
                    ascending
                };
                let pcurve = builder.add_curve2(piece.pcurve.clone());
                uses.push(EdgeUse {
                    edge: edge_ids[edge].as_ref().ok_or(BooleanError::Assembly)?.0,
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
                outer: ring == 0,
            });
        }
        face_ids.push((
            builder.topology_mut().add_face(Face {
                surface: Some(surface),
                bounds,
                orientation: orientation_of(k),
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
    let s = axiolid_evaluate::curve::locate3(edge_curve, mid, tolerance)
        .map_err(|_| BooleanError::Evaluation)?;
    let e = derivative3(edge_curve, s).map_err(|_| BooleanError::Evaluation)?;
    Ok(d.dot(e) > 0.0)
}
