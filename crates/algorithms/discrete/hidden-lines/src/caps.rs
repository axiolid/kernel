//! Near-plane caps of closed solids.
//!
//! A near plane that cuts a [`MeshRole::Solid`](axiolid_mesh_hidden_lines_contract::MeshRole)
//! removes what lies beyond it (`W > near`, nearer the viewer) and leaves
//! the solid's cross-section there: a cap that occludes. A point `P` of an
//! edge, with `W_P < near` (farther from the viewer than the plane), is hidden
//! by the cap of `M` exactly when the point of the plane above it,
//! `(U_P, V_P, near)`, lies in `M`. That is decided by ray parity: a ray from
//! it towards the viewer crosses the triangles of `M` that lie beyond the
//! plane an odd number of times if and only if it starts inside `M`.
//!
//! The ray runs through the image point perturbed by `δ = (ε, ε²)`, so it
//! never grazes an image edge: a point on the diagonal shared by two
//! triangles lies in exactly one of them. Triangles seen edge-on are skipped
//! for the same reason (they cover no area).
//!
//! A triangle `T` of `M` whose plane is exactly at `near` over the edge's
//! image is not a crossing (the strict clip is empty), and needs no
//! refusal. Over the image of `T` every point `P` with `W_P < near` is
//! already hidden by `T` itself: its point over `P'` has depth `near > W_P`
//! and is not removed. Outside that image `T` contributes no toggles in
//! either reading. So deciding the tie as "not crossed" cannot change a
//! drawing.

use std::ops::Range;

use axiolid_contracts::{GeomError, GeomResult};

use crate::edges::DrawEdge;
use crate::occlusion::Scene;
use crate::param::{Param, Solver, Span};
use crate::recipe::{Constraint, Kind, Lift};
use std::cmp::Ordering;

/// A solid the near plane cuts.
#[derive(Debug, Clone)]
pub(crate) struct Cap {
    /// Its triangles in the scene's triangle list.
    pub tris: Range<usize>,
}

/// The parts of `kept` hidden by the caps of `scene`, sorted and disjoint.
pub(crate) fn hidden_by_caps(
    solver: &mut Solver<'_>,
    scene: &Scene<'_>,
    edge: &DrawEdge,
    near: f64,
    kept: Span,
    pairs: &mut usize,
    max_pairs: usize,
) -> GeomResult<Vec<Span>> {
    let mut hidden: Vec<Span> = Vec::new();
    for cap in scene.caps {
        let mut ends: Vec<Param> = Vec::new();
        for t in cap.tris.clone() {
            let tri = &scene.tris[t];
            if tri.sigma == Ordering::Equal
                || tri.bounds.hi[2] < near
                || !tri.bounds.image_overlaps(&edge.bounds)
            {
                continue;
            }
            if *pairs >= max_pairs {
                return Err(GeomError::BudgetExceeded {
                    resource: "hidden-line candidate pairs",
                });
            }
            *pairs += 1;
            let clockwise = tri.sigma == Ordering::Less;
            let [v0, v1, v2] = tri.v;
            let mut span = kept;
            let mut covered = true;
            for (p, q) in [(v0, v1), (v1, v2), (v2, v0)] {
                let side = Constraint::new(
                    Kind::ImageSide {
                        p,
                        q,
                        neg: clockwise,
                    },
                    false,
                );
                match solver.clip_perturbed(span, &side)? {
                    Some(clipped) => span = clipped,
                    None => {
                        covered = false;
                        break;
                    }
                }
            }
            if !covered {
                continue;
            }
            let above = Constraint::new(
                Kind::PlaneAbove {
                    t: tri.v,
                    neg: !clockwise,
                    lift: Lift::Near(near),
                },
                true,
            );
            if let Some(crossing) = solver.clip(span, &above)? {
                sorted_insert(solver, &mut ends, crossing.lo)?;
                sorted_insert(solver, &mut ends, crossing.hi)?;
            }
        }
        // Odd parity between consecutive ends is inside the solid.
        let mut inside = false;
        let mut odd: Vec<Span> = Vec::new();
        for (i, &end) in ends.iter().enumerate() {
            inside = !inside;
            if let (true, Some(&next)) = (inside, ends.get(i + 1)) {
                if solver.cmp(end, next)? == Ordering::Less {
                    odd.push(Span { lo: end, hi: next });
                }
            }
        }
        // Equal depth does not occlude: not on the plane itself.
        let below = Constraint::new(Kind::BelowNear(near), true);
        for span in solver.clip_all(odd, &below)? {
            solver.insert(&mut hidden, span)?;
        }
    }
    Ok(hidden)
}

/// Insert `param` into the sorted `list`, after its equals.
fn sorted_insert(solver: &mut Solver<'_>, list: &mut Vec<Param>, param: Param) -> GeomResult<()> {
    let (mut low, mut high) = (0, list.len());
    while low < high {
        let mid = low + (high - low) / 2;
        if solver.cmp(list[mid], param)? == Ordering::Greater {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    list.insert(low, param);
    Ok(())
}
