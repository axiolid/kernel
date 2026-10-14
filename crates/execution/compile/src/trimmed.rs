//! A surface trimmed by a parameter rectangle, meshed as an open sheet
//! (#311).
//!
//! An atomic `Surface` node is unbounded (a plane, a cylinder along its
//! whole axis) and is still refused; a
//! [`SurfaceRelation::RectangularTrimmed`](axiolid_model::SurfaceRelation::RectangularTrimmed)
//! bounds one to `[u1, u2] x [v1, v2]`, so it has an area to mesh. The
//! mesh is a surface ([`axiolid_mesh_compile_contract::MeshClosure::Surface`]):
//! nothing downstream reads a volume from it, and a boolean refuses it.
//!
//! # The parameter rectangle
//!
//! Each sense says whether the trimmed surface's parameter runs with the
//! basis's (ISO 10303-42 `rectangular_trimmed_surface`): with it from `u1`
//! up to `u2`, against it from `u1` down to `u2`. On a periodic parameter
//! (round a cylinder, a cone, a sphere or a torus) a trim that runs past
//! the seam wraps by one period; on any other it is refused by name, as is
//! an empty range and one longer than a period. The mesh covers the
//! rectangle the trim spans, and faces the trimmed surface's normal: the
//! basis's, turned round when exactly one sense is against it.
//!
//! # The bound
//!
//! The rectangle is cut into a uniform grid of cells, each split into two
//! triangles whose corners are exact surface points. Over a parameter
//! triangle with corners `y_i` and a point `x = sum_i l_i y_i` in it, every
//! surface point `S(x)` is within `1/2 sum_i l_i q(y_i - x)` of the
//! triangle, `q(d) = A d_u^2 + 2 B |d_u d_v| + C d_v^2 <= (A + B) d_u^2 +
//! (C + B) d_v^2` with `A, B, C` certified bounds of the surface over the
//! cell ([`axiolid_reference::bound::SurfaceBounds::interpolation`]: the
//! second partials, or a rational B-spline's homogeneous form). `sum_i l_i
//! (y_i - x)_u^2` is the variance of the corners' `u` under the weights
//! `l_i`, and the corners' `u` lie in an interval of length `h_u`, so it is
//! at most `h_u^2 / 4` (Popoviciu); likewise in `v`. Each cell is
//! therefore within `1/8 ((A + B) h_u^2 + (C + B) h_v^2)` of its
//! triangles, the exact sagitta `r h_u^2 / 8` on a cylinder. The cell
//! counts are chosen from the bounds over the whole rectangle so that each
//! half fits half the chord budget, and every cell is then certified with
//! its own bounds; the worst is the reported [`DeviationBound::Certified`].
//! A cell a B-spline knot line of full multiplicity crosses has no
//! second-derivative bound across it and is reported
//! [`DeviationBound::Unbounded`] by name.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::Scalar;
use axiolid_mesh::TriMesh;
use axiolid_model::{GeometryGraph, GeometryNode, NodeId};
use axiolid_reference::bound::SurfaceBoundOracle;

use crate::deviation::DeviationBound;

/// Cells along one parameter.
const MAX_CELLS_PER_DIRECTION: usize = 1 << 12;

/// Cells of one patch.
const MAX_CELLS: usize = 1 << 20;

/// A trim's parameter rectangle and which way it faces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Trim {
    pub(crate) u: (Scalar, Scalar),
    pub(crate) v: (Scalar, Scalar),
    pub(crate) u_sense: bool,
    pub(crate) v_sense: bool,
}

/// Mesh `basis` over the trim within the options' chord budget: the mesh,
/// the basis's surface family, and the certified bound.
pub(crate) fn rectangular(
    graph: &GeometryGraph,
    basis: NodeId,
    trim: Trim,
    options: &ExecutionOptions,
) -> GeomResult<(TriMesh, &'static str, DeviationBound)> {
    let surface = match graph.get(basis) {
        Some(GeometryNode::Surface(surface)) => surface,
        Some(GeometryNode::SurfaceRelation(_)) => {
            return Err(refused(
                "a rectangular trim of a surface relation: only an atomic surface is trimmed",
            ))
        }
        Some(_) => {
            return Err(GeomError::InvalidInput(format!(
                "rectangular trim basis {basis:?} is not a surface"
            )))
        }
        None => {
            return Err(GeomError::InvalidInput(format!(
                "rectangular trim basis {basis:?} does not belong to this graph"
            )))
        }
    };
    let (u_period, v_period) = crate::brep::surface_periods(surface);
    let (u0, u1) = span(trim.u, trim.u_sense, u_period, "u")?;
    let (v0, v1) = span(trim.v, trim.v_sense, v_period, "v")?;
    let oracle = SurfaceBoundOracle::new(surface).ok_or_else(|| {
        refused("a rectangular trim of a surface with no certified derivative bounds")
    })?;
    let whole = oracle.bounds((u0, u1), (v0, v1)).ok_or_else(|| {
        GeomError::InvalidInput(
            "the trimmed surface is not defined over its whole rectangle (a cone past its apex)"
                .into(),
        )
    })?;
    let chord = crate::compiler::chord_error(options);
    let [a, b, c] = whole.interpolation;
    let nu = cells(u1 - u0, a + b, chord)?;
    let nv = cells(v1 - v0, c + b, chord)?;
    if nu.saturating_mul(nv) > MAX_CELLS {
        return Err(GeomError::BudgetExceeded {
            resource: "trimmed surface cells",
        });
    }
    let us = grid(u0, u1, nu);
    let vs = grid(v0, v1, nv);

    let mut positions = Vec::with_capacity((nu + 1) * (nv + 1));
    for &u in &us {
        for &v in &vs {
            positions.push(axiolid_reference::surface::evaluate(surface, u, v)?);
        }
    }
    let flip = trim.u_sense != trim.v_sense;
    let at = |i: usize, j: usize| u32::try_from(i * (nv + 1) + j).expect("cells are capped");
    let mut indices = Vec::with_capacity(nu * nv * 6);
    let mut worst: Scalar = 0.0;
    let mut unbounded = None;
    for i in 0..nu {
        for j in 0..nv {
            match oracle.bounds((us[i], us[i + 1]), (vs[j], vs[j + 1])) {
                Some(cell) if cell.smooth => {
                    let [a, b, c] = cell.interpolation;
                    let (hu, hv) = (us[i + 1] - us[i], vs[j + 1] - vs[j]);
                    let bound = 0.125 * ((a + b) * hu * hu + (c + b) * hv * hv);
                    worst = worst.max(bound * (1.0 + 1e-12));
                }
                Some(_) => {
                    unbounded = Some("a B-spline knot line of full multiplicity crosses a cell")
                }
                None => unbounded = Some("the surface has no certified bound over a cell"),
            }
            let corners = [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)];
            for [p, q, r] in [
                [corners[0], corners[1], corners[2]],
                [corners[0], corners[2], corners[3]],
            ] {
                let [pp, qp, rp] = [p, q, r].map(|k| positions[k as usize]);
                // A cell row on a pole (a sphere's) collapses a triangle to
                // the diagonal its neighbour already has.
                if pp == qp || qp == rp || rp == pp {
                    continue;
                }
                if flip {
                    indices.extend([p, r, q]);
                } else {
                    indices.extend([p, q, r]);
                }
            }
        }
    }
    let bound = match unbounded {
        Some(reason) => DeviationBound::Unbounded(reason),
        None => DeviationBound::Certified(worst),
    };
    Ok((
        TriMesh::new(positions, indices),
        crate::brep::surface_family(surface),
        bound,
    ))
}

fn refused(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: crate::BACKEND_ID,
        operation: Operation::SurfaceEvaluation,
        input,
    }
}

/// The parameter interval a trim covers, low end first.
fn span(
    (first, second): (Scalar, Scalar),
    sense: bool,
    period: Option<Scalar>,
    label: &str,
) -> GeomResult<(Scalar, Scalar)> {
    if !(first.is_finite() && second.is_finite()) {
        return Err(GeomError::InvalidInput(format!(
            "rectangular trim {label} range is not finite"
        )));
    }
    if first == second {
        return Err(GeomError::InvalidInput(format!(
            "rectangular trim {label} range is empty"
        )));
    }
    let (lo, hi) = match (sense, period) {
        (true, _) if second > first => (first, second),
        (false, _) if second < first => (second, first),
        // Past the seam: one more period in the trim's direction.
        (true, Some(period)) => (first, second + period),
        (false, Some(period)) => (second - period, first),
        (_, None) => {
            return Err(GeomError::InvalidInput(format!(
                "rectangular trim {label} range runs against its sense on a parameter that is \
                 not periodic"
            )))
        }
    };
    if period.is_some_and(|period| hi - lo > period) {
        return Err(GeomError::InvalidInput(format!(
            "rectangular trim {label} range is longer than one period"
        )));
    }
    Ok((lo, hi))
}

/// Cells along a parameter of length `length` so that `1/8 k h^2` is at
/// most half the chord budget: `h <= 2 sqrt(chord / k)`.
fn cells(length: Scalar, k: Scalar, chord: Scalar) -> GeomResult<usize> {
    if k <= 0.0 {
        return Ok(1);
    }
    let count = (0.5 * length * (k / chord).sqrt()).ceil();
    if !(count.is_finite() && count <= MAX_CELLS_PER_DIRECTION as Scalar) {
        return Err(GeomError::BudgetExceeded {
            resource: "trimmed surface cells",
        });
    }
    Ok((count as usize).max(1))
}

/// `n + 1` parameters from `lo` to `hi`, both ends exact.
fn grid(lo: Scalar, hi: Scalar, n: usize) -> Vec<Scalar> {
    (0..=n)
        .map(|i| {
            if i == n {
                hi
            } else {
                lo + (hi - lo) * (i as Scalar / n as Scalar)
            }
        })
        .collect()
}
