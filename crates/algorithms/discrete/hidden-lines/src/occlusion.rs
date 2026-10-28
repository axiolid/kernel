//! Visibility of one edge.
//!
//! The edge `A -> B` is cut by every kept triangle that can hide part of it.
//! A point `P` of the edge is hidden exactly when a kept triangle point lies
//! strictly in front of it at the same image point. Triangles are closed, so
//! touching an image boundary covers; equal depth never occludes; a triangle
//! seen edge-on occludes through the depth range it spans above its image.
//!
//! Each condition is an affine constraint in the edge parameter (`recipe`), so
//! the hidden part of the edge for one triangle is an interval or a short
//! union of intervals whose ends are roots, ordered exactly (`param`).

use std::cmp::Ordering;

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Aabb2;

use crate::broad::{Broad, Query};
use crate::caps::{hidden_by_caps, Cap};
use crate::decide::decide;
use crate::edges::{DrawEdge, Tri};
use crate::param::{Param, Solver, Span};
use crate::recipe::{Constraint, Delta, Kind, Lift, Orient2};
use crate::view::Vert;

/// The view-space scene one edge is resolved against.
pub(crate) struct Scene<'a> {
    pub verts: &'a [Vert],
    pub tris: &'a [Tri],
    pub near: Option<f64>,
    pub far: Option<f64>,
    pub region: Option<Aabb2>,
    /// The solids the near plane cuts.
    pub caps: &'a [Cap],
    /// The broad phase over `tris`.
    pub broad: &'a Broad,
}

/// A piece of an edge with one visibility.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Piece {
    pub lo: Param,
    pub hi: Param,
    pub hidden: bool,
}

pub(crate) enum Resolved {
    /// The edge's image is a point.
    Point,
    /// Depth clipping or the region removed the whole edge.
    Clipped,
    /// The visible and hidden pieces in parameter order.
    Pieces(Vec<Piece>),
}

/// Work counters shared by all edges.
#[derive(Debug, Default)]
pub(crate) struct Work {
    pub candidate_pairs: usize,
    pub escalations: usize,
}

/// Resolve `edge` against `scene`.
pub(crate) fn resolve(
    scene: &Scene<'_>,
    edge: &DrawEdge,
    max_pairs: usize,
    work: &mut Work,
) -> GeomResult<Resolved> {
    let [start, end] = edge.ends;
    let mut solver = Solver::new(scene.verts, start, end);
    let result = resolve_with(
        &mut solver,
        scene,
        edge,
        max_pairs,
        &mut work.candidate_pairs,
    );
    work.escalations += solver.escalations;
    result
}

fn resolve_with(
    solver: &mut Solver<'_>,
    scene: &Scene<'_>,
    edge: &DrawEdge,
    max_pairs: usize,
    pairs: &mut usize,
) -> GeomResult<Resolved> {
    let [start, end] = edge.ends;
    let verts = scene.verts;
    for axis in 0..2 {
        let delta = Delta {
            verts,
            from: start,
            to: end,
            axis,
        };
        if decide(&delta, &mut solver.escalations)? != Ordering::Equal {
            break;
        }
        if axis == 1 {
            return Ok(Resolved::Point);
        }
    }

    let mut kept = Some(Span::FULL);
    for constraint in kept_constraints(scene) {
        kept = match kept {
            Some(span) => solver.clip(span, &constraint)?,
            None => break,
        };
    }
    let Some(kept) = kept else {
        return Ok(Resolved::Clipped);
    };

    let mut hidden: Vec<Span> = Vec::new();
    let query = Query {
        image: edge.bounds,
        w_lo: edge.bounds.lo[2],
        w_hi: scene.near.unwrap_or(f64::INFINITY),
    };
    for t in scene.broad.candidates(scene.tris, &query) {
        let tri = &scene.tris[t];
        if *pairs >= max_pairs {
            return Err(GeomError::BudgetExceeded {
                resource: "hidden-line candidate pairs",
            });
        }
        *pairs += 1;
        if tri.sigma == Ordering::Equal {
            for span in edge_on(solver, scene, edge, tri, kept)? {
                solver.insert(&mut hidden, span)?;
            }
        } else if let Some(span) = facing(solver, scene, tri, kept)? {
            solver.insert(&mut hidden, span)?;
        }
    }

    if let Some(near) = scene.near {
        for span in hidden_by_caps(solver, scene, edge, near, kept, pairs, max_pairs)? {
            solver.insert(&mut hidden, span)?;
        }
    }

    let mut pieces = Vec::with_capacity(2 * hidden.len() + 1);
    let mut cursor = kept.lo;
    for span in hidden {
        if solver.cmp(cursor, span.lo)? == Ordering::Less {
            pieces.push(Piece {
                lo: cursor,
                hi: span.lo,
                hidden: false,
            });
        }
        pieces.push(Piece {
            lo: span.lo,
            hi: span.hi,
            hidden: true,
        });
        cursor = span.hi;
    }
    if solver.cmp(cursor, kept.hi)? == Ordering::Less {
        pieces.push(Piece {
            lo: cursor,
            hi: kept.hi,
            hidden: false,
        });
    }
    Ok(Resolved::Pieces(pieces))
}

/// The constraints that keep a point: depth range and region.
fn kept_constraints(scene: &Scene<'_>) -> Vec<Constraint> {
    let mut out = Vec::new();
    if let Some(near) = scene.near {
        out.push(Constraint::new(Kind::BelowNear(near), false));
    }
    if let Some(far) = scene.far {
        out.push(Constraint::new(Kind::AboveFar(far), false));
    }
    if let Some(region) = scene.region {
        for (axis, lower, upper) in [
            (0, region.min.x, region.max.x),
            (1, region.min.y, region.max.y),
        ] {
            for (upper_side, bound) in [(false, lower), (true, upper)] {
                out.push(Constraint::new(
                    Kind::Region {
                        axis,
                        upper: upper_side,
                        bound,
                    },
                    false,
                ));
            }
        }
    }
    out
}

/// The part of `kept` hidden by a triangle with a nonzero image area: where
/// the point is inside the triangle's image and the triangle's plane is
/// strictly in front of it, at a kept depth.
fn facing(
    solver: &mut Solver<'_>,
    scene: &Scene<'_>,
    tri: &Tri,
    kept: Span,
) -> GeomResult<Option<Span>> {
    let [v0, v1, v2] = tri.v;
    let clockwise = tri.sigma == Ordering::Less;
    let mut span = kept;
    for (p, q) in [(v0, v1), (v1, v2), (v2, v0)] {
        let side = Constraint::new(
            Kind::ImageSide {
                p,
                q,
                neg: clockwise,
            },
            false,
        );
        match solver.clip(span, &side)? {
            Some(clipped) => span = clipped,
            None => return Ok(None),
        }
    }
    let in_front = Constraint::new(
        Kind::PlaneAbove {
            t: tri.v,
            neg: !clockwise,
            lift: Lift::Own,
        },
        true,
    );
    let Some(mut span) = solver.clip(span, &in_front)? else {
        return Ok(None);
    };
    if let Some(near) = scene.near {
        // The triangle's point must not be removed by the near plane.
        let not_cut = Constraint::new(
            Kind::PlaneAbove {
                t: tri.v,
                neg: clockwise,
                lift: Lift::Near(near),
            },
            false,
        );
        match solver.clip(span, &not_cut)? {
            Some(clipped) => span = clipped,
            None => return Ok(None),
        }
    }
    Ok(Some(span))
}

/// The parts of `kept` hidden by a triangle seen edge-on.
///
/// Its image is a segment; only an edge lying along that segment can be
/// covered (a crossing is a single point). Over each image point the
/// triangle's fibre is the depth range between its extreme edges, so it hides
/// the point when some edge of the triangle over that image point is strictly
/// in front, and the fibre has a point at or beyond the near plane's kept side.
fn edge_on(
    solver: &mut Solver<'_>,
    scene: &Scene<'_>,
    edge: &DrawEdge,
    tri: &Tri,
    kept: Span,
) -> GeomResult<Vec<Span>> {
    let [start, end] = edge.ends;
    let verts = scene.verts;
    for corner in tri.v {
        let side = Orient2 {
            verts,
            a: start,
            b: end,
            c: corner,
        };
        if decide(&side, &mut solver.escalations)? != Ordering::Equal {
            return Ok(Vec::new());
        }
    }

    let [v0, v1, v2] = tri.v;
    let mut in_front: Vec<Span> = Vec::new();
    let mut reaches_kept: Vec<Span> = Vec::new();
    for (from, to) in [(v0, v1), (v1, v2), (v2, v0)] {
        let mut point_image = true;
        for axis in 0..2 {
            let delta = Delta {
                verts,
                from,
                to,
                axis,
            };
            if decide(&delta, &mut solver.escalations)? != Ordering::Equal {
                point_image = false;
                break;
            }
        }
        if point_image {
            continue;
        }
        let mut within = Some(kept);
        for kind in [Kind::ExtentFrom { from, to }, Kind::ExtentTo { from, to }] {
            within = match within {
                Some(span) => solver.clip(span, &Constraint::new(kind, false))?,
                None => break,
            };
        }
        let Some(within) = within else { continue };
        let above = Constraint::new(Kind::SegmentAbove { from, to }, true);
        if let Some(span) = solver.clip(within, &above)? {
            solver.insert(&mut in_front, span)?;
        }
        if let Some(near) = scene.near {
            let at_or_below = Constraint::new(Kind::SegmentAtOrBelow { from, to, near }, false);
            if let Some(span) = solver.clip(within, &at_or_below)? {
                solver.insert(&mut reaches_kept, span)?;
            }
        }
    }
    let Some(near) = scene.near else {
        return Ok(in_front);
    };
    let both = solver.intersect(&in_front, &reaches_kept)?;
    solver.clip_all(both, &Constraint::new(Kind::BelowNear(near), true))
}

#[cfg(test)]
mod tests {
    use axiolid_core::{Frame3, Point3, Vec3};
    use axiolid_mesh::TriMesh;

    use super::*;
    use crate::view::{Bounds, Vert};

    fn frame() -> Frame3 {
        Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        }
    }

    /// A triangle and an edge below it; returns the pieces.
    fn pieces(triangle: [Point3; 3], a: Point3, b: Point3) -> Vec<(bool, bool)> {
        let mesh = TriMesh::new(
            vec![triangle[0], triangle[1], triangle[2], a, b],
            vec![0, 1, 2],
        );
        let verts: Vec<Vert> = mesh
            .positions
            .iter()
            .map(|p| Vert::new(*p, &frame()))
            .collect();
        let v = [0u32, 1, 2];
        let sigma = decide(
            &Orient2 {
                verts: &verts,
                a: 0,
                b: 1,
                c: 2,
            },
            &mut 0,
        )
        .expect("decided");
        let tris = [Tri {
            v,
            sigma,
            bounds: Bounds::of(&verts, &v),
        }];
        let ends = [3, 4];
        let edge = DrawEdge {
            mesh: 0,
            local: [3, 4],
            ends,
            kind: axiolid_mesh_hidden_lines_contract::EdgeKind::Crease,
            bounds: Bounds::of(&verts, &ends),
        };
        let broad = Broad::build(&tris);
        let scene = Scene {
            verts: &verts,
            tris: &tris,
            near: None,
            far: None,
            region: None,
            caps: &[],
            broad: &broad,
        };
        match resolve(&scene, &edge, usize::MAX, &mut Work::default()).expect("resolved") {
            Resolved::Pieces(pieces) => pieces
                .iter()
                .map(|p| {
                    (
                        p.hidden,
                        matches!(p.lo, Param::Zero) && matches!(p.hi, Param::One),
                    )
                })
                .collect(),
            _ => panic!("expected pieces"),
        }
    }

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn a_triangle_strictly_in_front_hides_the_whole_edge_it_covers() {
        let tri = [p(0.0, 0.0, 1.0), p(4.0, 0.0, 1.0), p(0.0, 4.0, 1.0)];
        assert_eq!(
            pieces(tri, p(1.0, 1.0, 0.0), p(2.0, 1.0, 0.0)),
            [(true, true)]
        );
    }

    #[test]
    fn equal_depth_and_a_triangle_behind_do_not_hide() {
        let tri = [p(0.0, 0.0, 0.0), p(4.0, 0.0, 0.0), p(0.0, 4.0, 0.0)];
        assert_eq!(
            pieces(tri, p(1.0, 1.0, 0.0), p(2.0, 1.0, 0.0)),
            [(false, true)]
        );
        let behind = [p(0.0, 0.0, -1.0), p(4.0, 0.0, -1.0), p(0.0, 4.0, -1.0)];
        assert_eq!(
            pieces(behind, p(1.0, 1.0, 0.0), p(2.0, 1.0, 0.0)),
            [(false, true)]
        );
    }

    #[test]
    fn a_closed_image_boundary_covers_an_edge_running_along_it() {
        // The edge lies on the triangle's image boundary y = 0, below it.
        let tri = [p(0.0, 0.0, 1.0), p(4.0, 0.0, 1.0), p(0.0, 4.0, 1.0)];
        assert_eq!(
            pieces(tri, p(1.0, 0.0, 0.0), p(2.0, 0.0, 0.0)),
            [(true, true)]
        );
    }

    #[test]
    fn an_edge_leaving_the_image_is_split_at_the_boundary() {
        let tri = [p(0.0, 0.0, 1.0), p(4.0, 0.0, 1.0), p(0.0, 4.0, 1.0)];
        let split = pieces(tri, p(1.0, 1.0, 0.0), p(5.0, 1.0, 0.0));
        assert_eq!(split.iter().map(|s| s.0).collect::<Vec<_>>(), [true, false]);
    }
}
