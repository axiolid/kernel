//! Cross-section segments of interpenetrating triangles.
//!
//! For two triangles `a < b` that cover area, let `D_b` be the signed volume
//! against the plane of `b` (`det[b1-b0, b2-b0, x-b0]`). The chord of `a` is
//! where its edges meet that plane, the chord of `b` where its edges meet the
//! plane of `a`; both lie on the line where the planes meet, and the
//! triangles meet where the chords overlap. A chord that is an edge of its
//! triangle is already an edge segment, and coplanar or edge-sharing pairs
//! have no cross-section.

use std::cmp::Ordering;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_exact::Dyadic;
use axiolid_exact_arrangement::{PointId, Segment, CROSSING_SCRATCH_BYTES};

use super::points::Points;
use super::segments::plane_points;
use super::{Error, Faces, Work};
use crate::broad::Query;
use crate::decide::decide;
use crate::recipe::{Constraint, Kind, Lift, ValueSign};

/// A chord of one triangle: its two image points.
/// Bytes charged for one intersection segment: its two exact end points.
const SEGMENT_SCRATCH_BYTES: usize = 2 * CROSSING_SCRATCH_BYTES;

struct Chord {
    ends: [PointId; 2],
    is_edge: bool,
}

/// The chord of the triangle `own` against the plane of `other`, when it is
/// two points that differ in the image.
fn chord(
    faces: &Faces<'_>,
    points: &mut Points,
    own: usize,
    other: usize,
    work: &mut Work,
) -> GeomResult<Option<Chord>> {
    let corners = faces.tris[own].v;
    let constraint = Constraint::new(
        Kind::PlaneAbove {
            t: faces.tris[other].v,
            neg: false,
            lift: Lift::Own,
        },
        false,
    );
    let mut signs = [Ordering::Equal; 3];
    for (k, sign) in signs.iter_mut().enumerate() {
        let value = ValueSign {
            verts: faces.verts,
            constraint: &constraint,
            at: corners[k],
        };
        *sign = decide(&value, &mut work.escalations)?;
    }
    let value = |v: u32| constraint.value::<Dyadic>(faces.verts, v);
    let found = plane_points(points, faces.verts, corners, signs, value)?;
    let [(a, va), (b, vb)] = found[..] else {
        return Ok(None);
    };
    if points.store.same(a, b)? {
        return Ok(None);
    }
    Ok(Some(Chord {
        ends: [a, b],
        is_edge: va && vb,
    }))
}

/// The cross-section segments of every pair of triangles that cover area.
pub(crate) fn intersection_segments(
    faces: &Faces<'_>,
    points: &mut Points,
    work: &mut Work,
    declared: usize,
    options: &ExecutionOptions,
    out: &mut Vec<Segment>,
) -> Result<(), Error> {
    let first = out.len();
    for a in 0..faces.tris.len() {
        if !faces.covers(a) {
            continue;
        }
        options.check_cancelled()?;
        let bounds = faces.tris[a].bounds;
        let query = Query {
            image: bounds,
            w_lo: bounds.lo[2],
            w_hi: bounds.hi[2],
        };
        for b in faces.broad.candidates(faces.tris, &query) {
            if b <= a || !faces.covers(b) {
                continue;
            }
            let shared = faces.tris[a]
                .v
                .iter()
                .filter(|v| faces.tris[b].v.contains(v))
                .count();
            if shared >= 2 {
                continue;
            }
            work.budget.count().map_err(|e| faces.blame(a, e))?;
            if let Some(segment) =
                pair(faces, points, a, b, work).map_err(|error| faces.blame(a, error))?
            {
                out.push(segment);
                // Each segment keeps two exact points and may be split by
                // every other segment it meets; charge it like a crossing.
                let charged = SEGMENT_SCRATCH_BYTES.saturating_mul(out.len() - first);
                if let Some(budget) = options.memory_budget_bytes() {
                    if declared.saturating_add(charged) > budget {
                        return Err(Error::from(GeomError::BudgetExceeded {
                            resource: "memory",
                        }));
                    }
                }
            }
        }
    }
    Ok(())
}

fn pair(
    faces: &Faces<'_>,
    points: &mut Points,
    a: usize,
    b: usize,
    work: &mut Work,
) -> GeomResult<Option<Segment>> {
    let Some(first) = chord(faces, points, a, b, work)? else {
        return Ok(None);
    };
    let Some(second) = chord(faces, points, b, a, work)? else {
        return Ok(None);
    };
    if first.is_edge || second.is_edge {
        return Ok(None);
    }
    let order = |chord: &Chord, points: &Points| -> GeomResult<[PointId; 2]> {
        let [p, q] = chord.ends;
        Ok(if points.store.lexicographic(p, q)? == Ordering::Greater {
            [q, p]
        } else {
            [p, q]
        })
    };
    let (first, second) = (order(&first, points)?, order(&second, points)?);
    let start = if points.store.lexicographic(first[0], second[0])? == Ordering::Less {
        second[0]
    } else {
        first[0]
    };
    let end = if points.store.lexicographic(first[1], second[1])? == Ordering::Less {
        first[1]
    } else {
        second[1]
    };
    Ok((points.store.lexicographic(start, end)? == Ordering::Less)
        .then_some(Segment { ends: [start, end] }))
}

#[cfg(test)]
mod tests {
    use axiolid_contracts::GeomError;
    use axiolid_mesh::TriMesh;
    use axiolid_mesh_hidden_lines_contract::MeshRole::Sheet;

    use super::*;
    use crate::faces::segments::PAIRS;
    use crate::faces::testing::{flat, options, p, plan, sheet, Setup};

    /// A ramp over `x in [1,3]`, `y in [0,4]` rising from `z = -1` to `1`: it
    /// crosses the plane `z = 0` along `y = 2`.
    fn ramp() -> TriMesh {
        TriMesh::new(
            vec![
                p(1.0, 0.0, -1.0),
                p(3.0, 0.0, -1.0),
                p(3.0, 4.0, 1.0),
                p(1.0, 4.0, 1.0),
            ],
            vec![0, 1, 2, 0, 2, 3],
        )
    }

    fn found(parts: &[TriMesh], max_pairs: usize) -> Result<Vec<[(f64, f64); 2]>, Error> {
        let parts: Vec<_> = parts.iter().map(|m| (Sheet, m.clone())).collect();
        let setup = Setup::new(&parts, plan(), None);
        let mut work = setup.work(max_pairs);
        let faces = setup.faces(&mut work);
        let mut points = setup.points();
        let mut out = Vec::new();
        intersection_segments(&faces, &mut points, &mut work, 0, &options(), &mut out)?;
        Ok(out.iter().map(|s| flat(&points, s)).collect())
    }

    #[test]
    fn a_ramp_through_a_sheet_cuts_it_along_one_line() {
        let segments = found(&[sheet(0.0, 0.0, 4.0, 4.0, 0.0), ramp()], usize::MAX).expect("found");
        // Split at the sheets' diagonals, the two ends join at (2, 2).
        let mut merged: Vec<(f64, f64)> = segments.concat();
        merged.sort_by(|a, b| a.partial_cmp(b).expect("ordered"));
        merged.dedup();
        assert_eq!(merged, vec![(1.0, 2.0), (2.0, 2.0), (3.0, 2.0)]);
        assert!(segments.iter().all(|s| s[0].1 == 2.0 && s[1].1 == 2.0));
    }

    #[test]
    fn coplanar_and_edge_sharing_pairs_yield_none() {
        // One quad: its two triangles share an edge.
        assert!(found(&[sheet(0.0, 0.0, 2.0, 2.0, 0.0)], usize::MAX)
            .expect("found")
            .is_empty());
        // Two overlapping coplanar sheets.
        let both = [
            sheet(0.0, 0.0, 2.0, 2.0, 0.0),
            sheet(1.0, 1.0, 3.0, 3.0, 0.0),
        ];
        assert!(found(&both, usize::MAX).expect("found").is_empty());
    }

    #[test]
    fn planes_that_cross_outside_a_triangle_yield_none() {
        // The ramp meets `z = 0` along `y = 2`, outside this small sheet.
        let small = sheet(0.0, 0.0, 4.0, 1.0, 0.0);
        assert!(found(&[small, ramp()], usize::MAX)
            .expect("found")
            .is_empty());
    }

    #[test]
    fn a_segment_is_charged_against_the_memory_budget_when_it_is_found() {
        let parts: Vec<_> = [sheet(0.0, 0.0, 4.0, 4.0, 0.0), ramp()]
            .into_iter()
            .map(|m| (Sheet, m))
            .collect();
        let setup = Setup::new(&parts, plan(), None);
        let run = |bytes: usize| {
            let mut work = setup.work(usize::MAX);
            let faces = setup.faces(&mut work);
            let mut points = setup.points();
            let mut out = Vec::new();
            let options = options().with_memory_budget(bytes);
            intersection_segments(&faces, &mut points, &mut work, 0, &options, &mut out)
                .map(|()| out.len())
        };
        let count = run(usize::MAX).expect("found");
        assert!(count > 0);
        run(count * SEGMENT_SCRATCH_BYTES).expect("exactly enough");
        let error = run(count * SEGMENT_SCRATCH_BYTES - 1).expect_err("one byte short");
        assert_eq!(
            error.error,
            GeomError::BudgetExceeded { resource: "memory" }
        );
    }

    #[test]
    fn the_pair_budget_is_enforced() {
        let parts = [sheet(0.0, 0.0, 4.0, 4.0, 0.0), ramp()];
        let error = found(&parts, 0).expect_err("over budget");
        assert_eq!(error.error, GeomError::BudgetExceeded { resource: PAIRS });
    }
}
