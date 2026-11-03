//! A surface ruled between two sides, meshed as an open sheet (#402).
//!
//! [`SurfaceRelation::Ruled`](axiolid_model::SurfaceRelation::Ruled) is
//! `S(u, v) = (1 - v) A(u) + v B(u)` over `u, v` in `[0, 1]`, where a side
//! is an atomic curve over an oriented interval, `A(u) = c(t0 + u (t1 -
//! t0))`, or an apex point. The mesh is a surface
//! ([`axiolid_mesh_compile_contract::MeshClosure::Surface`]): nothing
//! downstream reads a volume from it, and a boolean refuses it. It faces
//! `S_u x S_v`. A side is a line, circle, ellipse or B-spline (rational
//! too); a conic runs at most one turn and a B-spline stays in its domain,
//! so no side is read where its curve is not defined. The kernel does not
//! convert a conic to a NURBS (its angle parameter is not rational), and a
//! ruled surface is a graph relation, not an `axiolid_surface::Surface`
//! value (ADR 0009).
//!
//! # The bound
//!
//! The mesh is a grid of columns `u_0 < ... < u_n` and `m` rows `v_j = j /
//! m`. Every vertex is an exact surface point, `P_ij = (1 - v_j) a_i + v_j
//! b_i` with `a_i = A(u_i)` and `b_i = B(u_i)`; an apex row is one shared
//! vertex. A cell `[u_i, u_{i+1}] x [v_j, v_{j+1}]` is two triangles split
//! along the `(i, j)`-`(i+1, j+1)` diagonal. Its distance to the surface
//! has two terms.
//!
//! **Curve term.** Let `l_A` and `l_B` be the linear interpolants of the
//! sides between the cell's columns and `L = (1 - v) l_A + v l_B`. Then
//! `S - L = (1 - v)(A - l_A) + v (B - l_B)`, a convex combination, and a
//! function with `|A''| <= K` is within `K h^2 / 8` of its chord's linear
//! interpolant (the Green's-function bound of the interpolation error,
//! tight for a parabola). With `h = h_u` the cell's width and `A''(u) =
//! (t1 - t0)^2 c''`, `|S - L| <= h_u^2 / 8 max_sides (t1 - t0)^2 sup|c''|`,
//! `sup|c''|` certified over the cell's parameters by
//! [`axiolid_reference::bound::curve_derivative_bounds3`]. This is the
//! Taylor form: a circle's sagitta bounds the distance to the chord's line,
//! not the gap between equal-parameter points, which is what `L` needs, so
//! `chord_bound3` is not used. An apex contributes 0.
//!
//! **Twist term.** `L` is bilinear over a sub-cell with corners `p00`,
//! `p10`, `p11`, `p01`: `L = p00 + s (p10 - p00) + w (p01 - p00) + s w T`
//! with the twist `T = p11 - p10 - p01 + p00`. Both triangles carry the
//! corners exactly. In the one under the diagonal (`w <= s`) the mesh is
//! `p00 + s (p10 - p00) + w (p01 - p00) + w T`, so it differs from `L` by
//! `w (1 - s) T`; in the other by `s (1 - w) T`. Each factor is at most
//! `1/4`. Here `T = h_v T_i` with `T_i = (b_{i+1} - b_i) - (a_{i+1} -
//! a_i)`, independent of the row, so a cell is within `h_v |T_i| / 4` of
//! `L`. With an apex `L` is the planar triangle `a_i, a_{i+1}, apex` (one
//! row), which the mesh covers: the term is 0.
//!
//! A cell is within `h_u^2 K_i / 8 + h_v |T_i| / 4` of the surface, and the
//! worst over the cells, inflated by `1 + 1e-12`, is the reported
//! [`DeviationBound::Certified`]. Each term is budgeted half the chord.
//!
//! # The grid
//!
//! Columns are chosen in four steps, all deterministic. A curve side's
//! tangent breaks ([`axiolid_reference::bound::continuity_breaks3`]) mapped
//! into `u` are forced columns, kept with the side's exact parameter, so no
//! cell crosses a kink. Each piece between forced columns is cut into `n =
//! ceil(len_u sqrt(K / chord) / 2)` cells from the piece's `K`, so the curve
//! term fits its half. A cell whose two end rulings both vanish (the sides
//! meet at both ends) would be a segment, hiding the lens between the
//! sides: it is halved until a midpoint ruling does not vanish. If the
//! midpoint vanishes too and the certified gap `h^2 / 8 (K_A + K_B)` is at
//! most half the chord (or both sides are straight), the sides are the same
//! curve there and the sheet is refused as `Degenerate`. The corners are
//! then read and every cell split into `ceil(sqrt(tau_0 / (2 chord)))`:
//! `T_i` is first order in the cell's width (about `h_u |B' - A'|`), so
//! `m` parts cut it `m`-fold, and `m` about `sqrt(tau_0 / (2 chord))` makes
//! the columns and the rows about equally many, `tau_0` the largest
//! `|T_i|`. The rows are `ceil(tau_1 / (2 chord))` from the largest `|T_i|`
//! of the final columns, so the twist term fits its half.
//!
//! # Collapse
//!
//! A column whose ruling has no length (`a_i == b_i`) is one vertex, as an
//! apex is, so the triangle that would be a sliver repeats an index and is
//! skipped while the other still covers the strip. A triangle of three
//! distinct vertices and exactly zero area (rulings crossing at a row) is a
//! collapse; no mesh is returned with one, the sheet is `Degenerate`.
//! Refusal of coincident sides is guaranteed where the rulings vanish
//! bitwise or a triangle collapses; a near-coincident sheet otherwise
//! meshes, with a bound that is sound.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_curve::Curve3;
use axiolid_mesh::TriMesh;
use axiolid_model::{GeometryGraph, GeometryNode, RuledSide};
use std::f64::consts::TAU;

use crate::deviation::DeviationBound;

/// Cells along one parameter.
const MAX_CELLS_PER_DIRECTION: usize = 1 << 12;

/// Cells of one sheet.
const MAX_CELLS: usize = 1 << 20;

const BUDGET: &str = "ruled surface cells";

/// A resolved side.
enum Side<'g> {
    Curve {
        curve: &'g Curve3,
        start: Scalar,
        end: Scalar,
    },
    Apex(Point3),
}

/// A column: `u`, and the exact curve parameter of each side where a break
/// forced it.
#[derive(Clone, Copy)]
struct Column {
    u: Scalar,
    exact: [Option<Scalar>; 2],
}

impl Side<'_> {
    /// The curve parameter at `u`, exact at both ends and at a break.
    fn parameter(&self, u: Scalar, exact: Option<Scalar>) -> Scalar {
        match *self {
            Self::Curve { start, end, .. } => match exact {
                Some(t) => t,
                None if u <= 0.0 => start,
                None if u >= 1.0 => end,
                None => start + (end - start) * u,
            },
            Self::Apex(_) => 0.0,
        }
    }

    fn point(&self, u: Scalar, exact: Option<Scalar>) -> GeomResult<Point3> {
        match self {
            Self::Curve { curve, .. } => {
                axiolid_reference::curve::evaluate3(curve, self.parameter(u, exact))
            }
            Self::Apex(point) => Ok(*point),
        }
    }

    /// `(t1 - t0)^2 sup|c''|` over the curve parameters `t_a` to `t_b`.
    fn second(&self, t_a: Scalar, t_b: Scalar) -> GeomResult<Scalar> {
        match self {
            Self::Curve {
                curve, start, end, ..
            } => {
                let bounds = axiolid_reference::bound::curve_derivative_bounds3(curve, t_a, t_b)
                    .ok_or_else(|| {
                        refused("a ruled surface side with no certified derivative bounds")
                    })?;
                Ok((end - start) * (end - start) * bounds.second)
            }
            Self::Apex(_) => Ok(0.0),
        }
    }
}

/// Mesh the ruled surface of `first` and `second` within the options' chord
/// budget: the mesh and the certified bound.
pub(crate) fn ruled(
    graph: &GeometryGraph,
    first: &RuledSide,
    second: &RuledSide,
    options: &ExecutionOptions,
) -> GeomResult<(TriMesh, DeviationBound)> {
    let sides = [resolve(graph, first)?, resolve(graph, second)?];
    let chord = crate::compiler::chord_error(options);
    let apex = sides.iter().any(|side| matches!(side, Side::Apex(_)));

    let mut columns = settle(&sides, coarse_columns(&sides, chord)?, chord)?;
    let mut rims = rim_points(&sides, &columns)?;
    if !apex {
        let split = ratio(twist(&rims), 2.0 * chord, true)?;
        if split > 1 {
            columns = settle(&sides, split_columns(&columns, split)?, chord)?;
            rims = rim_points(&sides, &columns)?;
        }
    }
    let rows = if apex {
        1
    } else {
        ratio(twist(&rims), 2.0 * chord, false)?
    };
    let cells = columns.len() - 1;
    if cells.saturating_mul(rows) > MAX_CELLS {
        return Err(GeomError::BudgetExceeded { resource: BUDGET });
    }

    let mut worst: Scalar = 0.0;
    for i in 0..cells {
        let (lo, hi) = (columns[i], columns[i + 1]);
        let width = hi.u - lo.u;
        let mut curve_term: Scalar = 0.0;
        for (s, side) in sides.iter().enumerate() {
            let second = side.second(
                side.parameter(lo.u, lo.exact[s]),
                side.parameter(hi.u, hi.exact[s]),
            )?;
            curve_term = curve_term.max(0.125 * width * width * second);
        }
        let twist = if apex {
            0.0
        } else {
            0.25 * twist_at(&rims, i) / rows as Scalar
        };
        worst = worst.max((curve_term + twist) * (1.0 + 1e-12));
    }
    if !worst.is_finite() {
        return Err(GeomError::InvalidInput(
            "the ruled surface has no finite deviation bound".into(),
        ));
    }

    Ok((
        build(&sides, &rims, rows)?,
        DeviationBound::Certified(worst),
    ))
}

fn refused(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: crate::BACKEND_ID,
        operation: Operation::SurfaceEvaluation,
        input,
    }
}

/// A side's curve and interval, refused by name where they cannot be read.
fn resolve<'g>(graph: &'g GeometryGraph, side: &RuledSide) -> GeomResult<Side<'g>> {
    match side {
        RuledSide::Apex(point) => {
            if !point.is_finite() {
                return Err(GeomError::InvalidInput(
                    "a ruled surface apex is not finite".into(),
                ));
            }
            Ok(Side::Apex(*point))
        }
        RuledSide::Curve { curve, interval } => {
            let node = graph.get(*curve).ok_or_else(|| {
                GeomError::InvalidInput(format!(
                    "ruled surface side {curve:?} does not belong to this graph"
                ))
            })?;
            let GeometryNode::Curve3(value) = node else {
                return Err(GeomError::InvalidInput(format!(
                    "ruled surface side {curve:?} is not a 3D curve"
                )));
            };
            let (start, end) = (interval.start, interval.end);
            if !(start.is_finite() && end.is_finite()) || start == end {
                return Err(GeomError::InvalidInput(
                    "ruled surface side interval is empty or not finite".into(),
                ));
            }
            match value {
                Curve3::Line(_) => {}
                Curve3::Circle(_) | Curve3::Ellipse(_) => {
                    if (end - start).abs() > TAU {
                        return Err(GeomError::InvalidInput(
                            "ruled surface side runs over more than one turn of its conic".into(),
                        ));
                    }
                }
                Curve3::BSpline(_) => {
                    let domain = axiolid_reference::curve::domain3(value);
                    let inside = |t: Scalar| t >= domain.start && t <= domain.end;
                    if !(inside(start) && inside(end)) {
                        return Err(GeomError::InvalidInput(
                            "ruled surface side interval leaves its curve's domain".into(),
                        ));
                    }
                }
                _ => {
                    return Err(refused(
                        "a ruled surface side on a curve other than a line, circle, ellipse or \
                         B-spline",
                    ))
                }
            }
            Ok(Side::Curve {
                curve: value,
                start,
                end,
            })
        }
        _ => Err(refused("a ruled surface side of an unknown kind")),
    }
}

/// The columns that fit the curve term: forced lines at the sides' tangent
/// breaks, each piece between them cut by its own bound.
fn coarse_columns(sides: &[Side; 2], chord: Scalar) -> GeomResult<Vec<Column>> {
    let mut forced: Vec<Column> = Vec::new();
    for (s, side) in sides.iter().enumerate() {
        let Side::Curve { curve, start, end } = side else {
            continue;
        };
        for t in axiolid_reference::bound::continuity_breaks3(curve, 1) {
            let u = (t - start) / (end - start);
            if !(u > 0.0 && u < 1.0) {
                continue;
            }
            match forced.iter_mut().find(|column| column.u == u) {
                Some(column) => column.exact[s] = Some(t),
                None => {
                    let mut exact = [None; 2];
                    exact[s] = Some(t);
                    forced.push(Column { u, exact });
                }
            }
        }
    }
    forced.sort_by(|a, b| a.u.total_cmp(&b.u));
    let mut edges = vec![Column {
        u: 0.0,
        exact: [None; 2],
    }];
    edges.extend(forced);
    edges.push(Column {
        u: 1.0,
        exact: [None; 2],
    });

    let mut columns = vec![edges[0]];
    for pair in edges.windows(2) {
        let (lo, hi) = (pair[0], pair[1]);
        let mut k: Scalar = 0.0;
        for (s, side) in sides.iter().enumerate() {
            k = k.max(side.second(
                side.parameter(lo.u, lo.exact[s]),
                side.parameter(hi.u, hi.exact[s]),
            )?);
        }
        let n = if k > 0.0 {
            ratio((hi.u - lo.u) * (k / chord).sqrt(), 2.0, false)?
        } else {
            1
        };
        if columns.len() + n > MAX_CELLS_PER_DIRECTION + 1 {
            return Err(GeomError::BudgetExceeded { resource: BUDGET });
        }
        for step in 1..n {
            columns.push(Column {
                u: lo.u + (hi.u - lo.u) * (step as Scalar / n as Scalar),
                exact: [None; 2],
            });
        }
        columns.push(hi);
    }
    Ok(columns)
}

/// Halve every cell whose end rulings both vanish (`a_i == b_i` at both
/// columns) until a midpoint ruling does not. A cell is equal, and the sheet
/// refused as `Degenerate`, when the sides are straight (the gap `|A - B|`
/// is linear and vanishes at both ends), or the midpoint ruling vanishes
/// too and the certified gap `h^2 / 8 (K_A + K_B)` is at most half the
/// chord: the sides are the same curve to the budget. Otherwise a vanishing
/// cell is halved, and its gap falls fourfold.
fn settle(sides: &[Side; 2], mut columns: Vec<Column>, chord: Scalar) -> GeomResult<Vec<Column>> {
    let vanishes = |column: &Column| -> GeomResult<bool> {
        Ok(sides[0].point(column.u, column.exact[0])?
            == sides[1].point(column.u, column.exact[1])?)
    };
    let mut vanish = columns
        .iter()
        .map(vanishes)
        .collect::<GeomResult<Vec<bool>>>()?;
    let mut equal = 0;
    let mut i = 0;
    while i + 1 < columns.len() {
        if !(vanish[i] && vanish[i + 1]) {
            i += 1;
            continue;
        }
        let (lo, hi) = (columns[i], columns[i + 1]);
        let mut k: Scalar = 0.0;
        for (s, side) in sides.iter().enumerate() {
            k += side.second(
                side.parameter(lo.u, lo.exact[s]),
                side.parameter(hi.u, hi.exact[s]),
            )?;
        }
        let width = hi.u - lo.u;
        let mid = Column {
            u: lo.u + 0.5 * (hi.u - lo.u),
            exact: [None; 2],
        };
        let mid_vanishes = vanishes(&mid)?;
        if k == 0.0 || (mid_vanishes && 0.125 * width * width * k <= 0.5 * chord) {
            equal += 1;
            i += 1;
            continue;
        }
        if columns.len() > MAX_CELLS_PER_DIRECTION {
            return Err(GeomError::BudgetExceeded { resource: BUDGET });
        }
        vanish.insert(i + 1, mid_vanishes);
        columns.insert(i + 1, mid);
    }
    if equal == 0 {
        Ok(columns)
    } else if equal == columns.len() - 1 {
        Err(GeomError::Degenerate(
            "the ruled surface has no area: its sides coincide".into(),
        ))
    } else {
        Err(GeomError::Degenerate(
            "the ruled surface's sides coincide over part of their span".into(),
        ))
    }
}

/// Every cell cut into `split` equal ones in `u`.
fn split_columns(columns: &[Column], split: usize) -> GeomResult<Vec<Column>> {
    if (columns.len() - 1).saturating_mul(split) > MAX_CELLS_PER_DIRECTION {
        return Err(GeomError::BudgetExceeded { resource: BUDGET });
    }
    let mut out = vec![columns[0]];
    for pair in columns.windows(2) {
        for step in 1..split {
            out.push(Column {
                u: pair[0].u + (pair[1].u - pair[0].u) * (step as Scalar / split as Scalar),
                exact: [None; 2],
            });
        }
        out.push(pair[1]);
    }
    Ok(out)
}

/// `ceil(value / divisor)` at least 1 (or, `root`, `ceil(sqrt(value /
/// divisor))`), refused past a direction's cap.
fn ratio(value: Scalar, divisor: Scalar, root: bool) -> GeomResult<usize> {
    let mut quotient = value / divisor;
    if root {
        quotient = quotient.sqrt();
    }
    let count = quotient.ceil();
    if !(count.is_finite() && count <= MAX_CELLS_PER_DIRECTION as Scalar) {
        return Err(GeomError::BudgetExceeded { resource: BUDGET });
    }
    Ok((count as usize).max(1))
}

/// The two sides' points at every column.
fn rim_points(sides: &[Side; 2], columns: &[Column]) -> GeomResult<Vec<[Point3; 2]>> {
    columns
        .iter()
        .map(|column| {
            Ok([
                sides[0].point(column.u, column.exact[0])?,
                sides[1].point(column.u, column.exact[1])?,
            ])
        })
        .collect()
}

/// `|T_i|` of cell `i`: how far the rulings at its ends are from parallel
/// and equal in length.
fn twist_at(rim: &[[Point3; 2]], i: usize) -> Scalar {
    let [a0, b0] = rim[i];
    let [a1, b1] = rim[i + 1];
    ((b1 - b0) - (a1 - a0)).length()
}

/// The largest `|T_i|`.
fn twist(rim: &[[Point3; 2]]) -> Scalar {
    (0..rim.len() - 1)
        .map(|i| twist_at(rim, i))
        .fold(0.0, Scalar::max)
}

/// The sheet: `rows` rows of vertices over the columns, an apex row shared,
/// and one vertex for a column whose ruling vanishes. A triangle of three
/// distinct vertices and no area is a collapse the sheet is refused for.
fn build(sides: &[Side; 2], rim: &[[Point3; 2]], rows: usize) -> GeomResult<TriMesh> {
    let mut positions: Vec<Point3> = Vec::new();
    let mut index = vec![0u32; rim.len() * (rows + 1)];
    let mut shared: [Option<u32>; 2] = [None; 2];
    for (i, [a, b]) in rim.iter().enumerate() {
        if a == b {
            // The ruling has no length: one vertex, shared with the apex
            // when a side is one.
            let apex = sides.iter().position(|side| matches!(side, Side::Apex(_)));
            let vertex = match apex {
                Some(s) => *shared[s].get_or_insert_with(|| {
                    positions.push(*a);
                    u32::try_from(positions.len() - 1).expect("cells are capped")
                }),
                None => {
                    positions.push(*a);
                    u32::try_from(positions.len() - 1).expect("cells are capped")
                }
            };
            index[i * (rows + 1)..(i + 1) * (rows + 1)].fill(vertex);
            continue;
        }
        for j in 0..=rows {
            let end = if j == 0 {
                Some(0)
            } else if j == rows {
                Some(1)
            } else {
                None
            };
            let vertex = match end {
                Some(s) if matches!(sides[s], Side::Apex(_)) => {
                    *shared[s].get_or_insert_with(|| {
                        positions.push(if s == 0 { *a } else { *b });
                        u32::try_from(positions.len() - 1).expect("cells are capped")
                    })
                }
                _ => {
                    let v = j as Scalar / rows as Scalar;
                    positions.push(*a * (1.0 - v) + *b * v);
                    u32::try_from(positions.len() - 1).expect("cells are capped")
                }
            };
            index[i * (rows + 1) + j] = vertex;
        }
    }
    let at = |i: usize, j: usize| index[i * (rows + 1) + j];
    let mut indices = Vec::with_capacity((rim.len() - 1) * rows * 6);
    for i in 0..rim.len() - 1 {
        for j in 0..rows {
            let corners = [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)];
            for [p, q, r] in [
                [corners[0], corners[1], corners[2]],
                [corners[0], corners[2], corners[3]],
            ] {
                if p == q || q == r || r == p {
                    continue;
                }
                let [pp, qp, rp] = [p, q, r].map(|k| positions[k as usize]);
                if (qp - pp).cross(rp - pp) == Vec3::ZERO {
                    return Err(GeomError::Degenerate(format!(
                        "the ruled surface collapses to zero area over its cell {i}"
                    )));
                }
                indices.extend([p, q, r]);
            }
        }
    }
    Ok(TriMesh::new(positions, indices))
}
