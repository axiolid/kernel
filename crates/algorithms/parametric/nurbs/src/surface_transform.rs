//! Exact tensor-product NURBS surface transformations.

use axiolid_contracts::GeomResult;
use axiolid_core::{Point3, Scalar};
use axiolid_curve::BSplineCurve;
use axiolid_evaluate::surface::bspline_jet;
use axiolid_surface::BSplineSurface;

use crate::{axis::reverse_axis, transform::insert_knot3};

/// Insert one U-axis knot into every homogeneous control-net column.
///
/// The represented tensor-product surface is unchanged.
pub fn insert_surface_knot_u(
    surface: &BSplineSurface,
    parameter: Scalar,
) -> GeomResult<BSplineSurface> {
    bspline_jet(surface, parameter, 0.0)?;
    let columns = (0..surface.control_points[0].len())
        .map(|v| insert_knot3(&u_curve(surface, v), parameter))
        .collect::<GeomResult<Vec<_>>>()?;
    let u_count = columns[0].control_points.len();
    let control_points = (0..u_count)
        .map(|u| {
            columns
                .iter()
                .map(|column| column.control_points[u])
                .collect()
        })
        .collect();
    let weights = surface.weights.as_ref().map(|_| {
        (0..u_count)
            .map(|u| {
                columns
                    .iter()
                    .map(|column| column.weights.as_ref().expect("source was rational")[u])
                    .collect()
            })
            .collect()
    });
    Ok(BSplineSurface {
        u_degree: surface.u_degree,
        v_degree: surface.v_degree,
        control_points,
        u_knots: columns[0].knots.clone(),
        u_multiplicities: columns[0].multiplicities.clone(),
        v_knots: surface.v_knots.clone(),
        v_multiplicities: surface.v_multiplicities.clone(),
        weights,
        knot_spec: surface.knot_spec,
        u_closed: surface.u_closed,
        v_closed: surface.v_closed,
        self_intersect: surface.self_intersect,
    })
}

/// Insert one V-axis knot into every homogeneous control-net row.
///
/// The represented tensor-product surface is unchanged.
pub fn insert_surface_knot_v(
    surface: &BSplineSurface,
    parameter: Scalar,
) -> GeomResult<BSplineSurface> {
    bspline_jet(surface, 0.0, parameter)?;
    let rows = (0..surface.control_points.len())
        .map(|u| insert_knot3(&v_curve(surface, u), parameter))
        .collect::<GeomResult<Vec<_>>>()?;
    let control_points = rows.iter().map(|row| row.control_points.clone()).collect();
    let weights = surface.weights.as_ref().map(|_| {
        rows.iter()
            .map(|row| row.weights.as_ref().expect("source was rational").clone())
            .collect()
    });
    Ok(BSplineSurface {
        u_degree: surface.u_degree,
        v_degree: surface.v_degree,
        control_points,
        u_knots: surface.u_knots.clone(),
        u_multiplicities: surface.u_multiplicities.clone(),
        v_knots: rows[0].knots.clone(),
        v_multiplicities: rows[0].multiplicities.clone(),
        weights,
        knot_spec: surface.knot_spec,
        u_closed: surface.u_closed,
        v_closed: surface.v_closed,
        self_intersect: surface.self_intersect,
    })
}

/// Reverse the surface's U parameter without changing its image.
pub fn reverse_surface_u(surface: &BSplineSurface) -> GeomResult<BSplineSurface> {
    bspline_jet(surface, 0.0, 0.0)?;
    let mut result = surface.clone();
    result.control_points.reverse();
    if let Some(weights) = &mut result.weights {
        weights.reverse();
    }
    (result.u_knots, result.u_multiplicities) =
        reverse_axis(&surface.u_knots, &surface.u_multiplicities)?;
    Ok(result)
}

/// Reverse the surface's V parameter without changing its image.
pub fn reverse_surface_v(surface: &BSplineSurface) -> GeomResult<BSplineSurface> {
    bspline_jet(surface, 0.0, 0.0)?;
    let mut result = surface.clone();
    for row in &mut result.control_points {
        row.reverse();
    }
    if let Some(weights) = &mut result.weights {
        for row in weights {
            row.reverse();
        }
    }
    (result.v_knots, result.v_multiplicities) =
        reverse_axis(&surface.v_knots, &surface.v_multiplicities)?;
    Ok(result)
}

pub(crate) fn u_curve(surface: &BSplineSurface, v: usize) -> BSplineCurve<Point3> {
    BSplineCurve {
        degree: surface.u_degree,
        control_points: surface.control_points.iter().map(|row| row[v]).collect(),
        knots: surface.u_knots.clone(),
        multiplicities: surface.u_multiplicities.clone(),
        weights: surface
            .weights
            .as_ref()
            .map(|rows| rows.iter().map(|row| row[v]).collect()),
        knot_spec: surface.knot_spec,
        closed: surface.u_closed,
        self_intersect: surface.self_intersect,
    }
}

fn v_curve(surface: &BSplineSurface, u: usize) -> BSplineCurve<Point3> {
    BSplineCurve {
        degree: surface.v_degree,
        control_points: surface.control_points[u].clone(),
        knots: surface.v_knots.clone(),
        multiplicities: surface.v_multiplicities.clone(),
        weights: surface.weights.as_ref().map(|rows| rows[u].clone()),
        knot_spec: surface.knot_spec,
        closed: surface.v_closed,
        self_intersect: surface.self_intersect,
    }
}
