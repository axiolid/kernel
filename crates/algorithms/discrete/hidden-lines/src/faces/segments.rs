//! Segments of the arrangement that come from edges, clipping and the crop.

use std::cmp::Ordering;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::Point2;
use axiolid_exact::Dyadic;
use axiolid_exact_arrangement::{PointId, Segment};
use axiolid_mesh_hidden_lines_contract::EdgeKind;

use super::points::Points;
use super::{Error, Faces, Work};
use crate::decide::decide;
use crate::edges::DrawEdge;
use crate::occlusion::{resolve, Resolved, Work as Resolver};
use crate::recipe::{Constraint, Kind, ValueSign};
use crate::view::Bounds;

/// The resource a pair-count refusal names.
pub(crate) const PAIRS: &str = "visible-face candidate pairs";

/// The visible pieces of the edges of every triangle that covers area.
///
/// With `prune` off every piece becomes a segment, hidden or not.
pub(crate) fn edge_segments(
    faces: &Faces<'_>,
    points: &mut Points,
    prune: bool,
    work: &mut Work,
    options: &ExecutionOptions,
    out: &mut Vec<Segment>,
) -> Result<(), Error> {
    let mut edges: Vec<(u32, u32, usize)> = Vec::new();
    for (t, tri) in faces.tris.iter().enumerate() {
        if !faces.covers(t) {
            continue;
        }
        for (i, j) in [(0, 1), (1, 2), (2, 0)] {
            let (a, b) = (tri.v[i], tri.v[j]);
            edges.push((a.min(b), a.max(b), t));
        }
    }
    edges.sort_unstable_by_key(|&(a, b, _)| (a, b));
    edges.dedup_by_key(|&mut (a, b, _)| (a, b));
    let scene = faces.scene();
    let mut resolver = Resolver {
        candidate_pairs: work.budget.used,
        escalations: 0,
    };
    for (a, b, t) in edges {
        options.check_cancelled()?;
        let mesh = faces.mesh_of[t];
        let ends = [a, b];
        let edge = DrawEdge {
            mesh,
            local: ends.map(|v| v - faces.bases[mesh]),
            ends,
            kind: EdgeKind::Boundary,
            bounds: Bounds::of(faces.verts, &ends),
        };
        let resolved = resolve(&scene, &edge, work.budget.max, &mut resolver)
            .map_err(|error| faces.blame(t, rename(error)))?;
        let Resolved::Pieces(pieces) = resolved else {
            continue;
        };
        for piece in pieces {
            if piece.hidden && prune {
                continue;
            }
            let [from, to] =
                [piece.lo, piece.hi].map(|param| points.at_param(faces.verts, a, b, param));
            out.push(Segment {
                ends: [
                    from.map_err(|error| faces.blame(t, error))?,
                    to.map_err(|error| faces.blame(t, error))?,
                ],
            });
        }
    }
    work.budget.used = resolver.candidate_pairs;
    work.escalations += resolver.escalations;
    Ok(())
}

/// The refusal of the resolver names its own resource; faces name theirs.
fn rename(error: GeomError) -> GeomError {
    match error {
        GeomError::BudgetExceeded {
            resource: "hidden-line candidate pairs",
        } => GeomError::BudgetExceeded { resource: PAIRS },
        other => other,
    }
}

/// The plane points of a triangle against a plane, as `(point, is_vertex)`:
/// the vertices on the plane and the crossings of edges whose ends lie
/// strictly on opposite sides. At most two for a triangle that is not
/// contained in the plane.
pub(crate) fn plane_points(
    points: &mut Points,
    verts: &[crate::view::Vert],
    corners: [u32; 3],
    signs: [Ordering; 3],
    value: impl Fn(u32) -> Dyadic,
) -> GeomResult<Vec<(PointId, bool)>> {
    let mut found = Vec::with_capacity(2);
    for k in 0..3 {
        if signs[k] == Ordering::Equal {
            found.push((points.vertex(verts, corners[k])?, true));
        }
    }
    for (i, j) in [(0, 1), (1, 2), (2, 0)] {
        if signs[i] != Ordering::Equal && signs[j] != Ordering::Equal && signs[i] != signs[j] {
            let ends = [corners[i], corners[j]];
            let values = ends.map(&value);
            found.push((
                points.crossing(verts, ends, [&values[0], &values[1]])?,
                false,
            ));
        }
    }
    Ok(found)
}

/// The segments where a triangle meets the near or far plane.
///
/// A triangle needs one with a vertex strictly on the removed side and two
/// plane points. Triangles that cover area, and every triangle of a solid the
/// near plane cuts (seen edge-on or not, since a wall seen edge-on carries
/// the outline of the cap), contribute to the near plane.
pub(crate) fn clip_segments(
    faces: &Faces<'_>,
    points: &mut Points,
    work: &mut Work,
    options: &ExecutionOptions,
    out: &mut Vec<Segment>,
) -> Result<(), Error> {
    let mut in_cap = vec![false; faces.tris.len()];
    for cap in faces.caps {
        in_cap[cap.tris.clone()].fill(true);
    }
    for (t, tri) in faces.tris.iter().enumerate() {
        if t % 1024 == 0 {
            options.check_cancelled()?;
        }
        if !(faces.covers(t) || in_cap[t]) {
            continue;
        }
        // (kind, whether the removed side is a strictly negative value).
        let mut kinds = Vec::with_capacity(2);
        if let Some(near) = faces.near {
            kinds.push(Kind::BelowNear(near));
        }
        if let (Some(far), true) = (faces.far, faces.covers(t)) {
            kinds.push(Kind::AboveFar(far));
        }
        for kind in kinds {
            let constraint = Constraint::new(kind, false);
            let found = (|| -> GeomResult<Option<[PointId; 2]>> {
                let mut signs = [Ordering::Equal; 3];
                for (k, sign) in signs.iter_mut().enumerate() {
                    let value = ValueSign {
                        verts: faces.verts,
                        constraint: &constraint,
                        at: tri.v[k],
                    };
                    *sign = decide(&value, &mut work.escalations)?;
                }
                // The removed side has a strictly negative value.
                if !signs.contains(&Ordering::Less) {
                    return Ok(None);
                }
                let value = |v: u32| constraint.value::<Dyadic>(faces.verts, v);
                let found = plane_points(points, faces.verts, tri.v, signs, value)?;
                Ok(match found[..] {
                    [(a, _), (b, _)] => Some([a, b]),
                    _ => None,
                })
            })()
            .map_err(|error| faces.blame(t, error))?;
            if let Some(ends) = found {
                out.push(Segment { ends });
            }
        }
    }
    Ok(())
}

/// The four sides of the crop box.
pub(crate) fn crop_segments(
    faces: &Faces<'_>,
    points: &mut Points,
    out: &mut Vec<Segment>,
) -> Result<(), Error> {
    let Some(region) = faces.region else {
        return Ok(());
    };
    let (min, max) = (region.min, region.max);
    let corners = [
        Point2::new(min.x, min.y),
        Point2::new(max.x, min.y),
        Point2::new(max.x, max.y),
        Point2::new(min.x, max.y),
    ];
    let mut ids = Vec::with_capacity(4);
    for corner in corners {
        ids.push(points.store.push_flat(corner)?);
    }
    for k in 0..4 {
        let next = (k + 1) % 4;
        if corners[k] != corners[next] {
            out.push(Segment {
                ends: [ids[k], ids[next]],
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use axiolid_core::{Aabb2, Point2};
    use axiolid_mesh::TriMesh;
    use axiolid_mesh_hidden_lines_contract::MeshRole::{Sheet, Solid};

    use super::*;
    use crate::faces::testing::{flat, options, p, plan, sheet, solid, Setup};

    fn edges(setup: &Setup, prune: bool) -> (Vec<Segment>, Points) {
        let mut work = setup.work(usize::MAX);
        let faces = setup.faces(&mut work);
        let mut points = setup.points();
        let mut out = Vec::new();
        edge_segments(&faces, &mut points, prune, &mut work, &options(), &mut out).expect("edges");
        (out, points)
    }

    #[test]
    fn a_fully_hidden_edge_yields_no_segment() {
        let under = sheet(0.5, 0.5, 1.5, 1.5, -1.0);
        let alone = Setup::new(&[(Solid, solid([0.0; 3], [2.0, 2.0, 1.0]))], plan(), None);
        let both = Setup::new(
            &[(Solid, solid([0.0; 3], [2.0, 2.0, 1.0])), (Sheet, under)],
            plan(),
            None,
        );
        assert_eq!(edges(&alone, true).0.len(), edges(&both, true).0.len());
        // Unpruned, the sheet's five edges stay in.
        assert_eq!(
            edges(&both, false).0.len(),
            edges(&alone, false).0.len() + 5
        );
    }

    #[test]
    fn an_edge_on_triangle_contributes_no_edge_of_its_own() {
        let wall = TriMesh::new(
            vec![
                p(0.0, 0.0, 0.0),
                p(2.0, 0.0, 0.0),
                p(2.0, 0.0, 3.0),
                p(0.0, 0.0, 3.0),
            ],
            vec![0, 1, 2, 0, 2, 3],
        );
        let setup = Setup::new(&[(Sheet, wall)], plan(), None);
        assert!(edges(&setup, false).0.is_empty());
    }

    #[test]
    fn an_edge_on_an_edge_on_walls_line_is_kept() {
        // The wall stands on the floor's bottom edge and is edge-on in plan.
        let wall = TriMesh::new(
            vec![
                p(0.0, 0.0, 0.0),
                p(2.0, 0.0, 0.0),
                p(2.0, 0.0, 3.0),
                p(0.0, 0.0, 3.0),
            ],
            vec![0, 1, 2, 0, 2, 3],
        );
        let setup = Setup::new(
            &[(Sheet, sheet(0.0, 0.0, 2.0, 2.0, 0.0)), (Sheet, wall)],
            plan(),
            None,
        );
        let (segments, points) = edges(&setup, true);
        let found = segments
            .iter()
            .any(|s| flat(&points, s) == [(0.0, 0.0), (2.0, 0.0)]);
        assert!(found, "the floor edge under the wall is a segment");
    }

    #[test]
    fn the_near_clip_of_a_cut_boxes_edge_on_face_is_present() {
        let setup = Setup::new(
            &[(Solid, solid([0.0; 3], [2.0, 2.0, 2.0]))],
            plan().with_near(1.0),
            None,
        );
        let mut work = setup.work(usize::MAX);
        let faces = setup.faces(&mut work);
        assert_eq!(faces.caps.len(), 1);
        let mut points = setup.points();
        let mut out = Vec::new();
        clip_segments(&faces, &mut points, &mut work, &options(), &mut out).expect("clip");
        // The wall at y = 0 is edge-on; its two triangles cut at z = 1 join
        // (0,0), (1,0) and (2,0).
        let on_wall: Vec<[(f64, f64); 2]> = out
            .iter()
            .map(|s| flat(&points, s))
            .filter(|ends| ends.iter().all(|e| e.1 == 0.0))
            .collect();
        let touches = |x: f64| on_wall.iter().any(|ends| ends.contains(&(x, 0.0)));
        assert!(touches(0.0) && touches(1.0) && touches(2.0), "{on_wall:?}");
    }

    #[test]
    fn the_crop_gives_the_four_sides_and_skips_degenerate_ones() {
        let region = Aabb2 {
            min: Point2::new(0.0, 0.0),
            max: Point2::new(2.0, 1.0),
        };
        let setup = Setup::new(
            &[(Sheet, sheet(0.0, 0.0, 4.0, 4.0, 0.0))],
            plan(),
            Some(region),
        );
        let mut work = setup.work(usize::MAX);
        let faces = setup.faces(&mut work);
        let mut points = setup.points();
        let mut out = Vec::new();
        crop_segments(&faces, &mut points, &mut out).expect("crop");
        assert_eq!(out.len(), 4);
        let flat_region = Aabb2 {
            min: Point2::new(0.0, 0.0),
            max: Point2::new(2.0, 0.0),
        };
        let setup = Setup::new(
            &[(Sheet, sheet(0.0, 0.0, 4.0, 4.0, 0.0))],
            plan(),
            Some(flat_region),
        );
        let mut work = setup.work(usize::MAX);
        let faces = setup.faces(&mut work);
        let mut points = setup.points();
        let mut out = Vec::new();
        crop_segments(&faces, &mut points, &mut out).expect("crop");
        assert_eq!(out.len(), 2, "the two vertical sides have no length");
    }

    #[test]
    fn the_pair_budget_is_enforced_for_edges() {
        let setup = Setup::new(&[(Solid, solid([0.0; 3], [2.0, 2.0, 1.0]))], plan(), None);
        let mut work = setup.work(0);
        let faces = setup.faces(&mut work);
        let mut points = setup.points();
        let error = edge_segments(
            &faces,
            &mut points,
            true,
            &mut work,
            &options(),
            &mut Vec::new(),
        )
        .expect_err("over budget");
        assert_eq!(error.error, GeomError::BudgetExceeded { resource: PAIRS });
    }
}
