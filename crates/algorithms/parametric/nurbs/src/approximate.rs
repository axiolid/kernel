//! Least-squares curve and surface approximation (#150, ledger B11).
//!
//! # Approximation, not interpolation
//!
//! `fit.rs` interpolates: the curve passes through every input point.
//! This module approximates: the curve or surface minimises the sum of
//! squared deviations from scattered points, using fewer control points
//! than data points. The two are different operations with different
//! uses -- a caller fairing noisy survey points wants the first kind of
//! fit, not the second.
//!
//! # Method
//!
//! Parameters are assigned per point (uniform, chord-length or
//! centripetal, Piegl and Tiller section 9.2), the knot vector is placed
//! by the averaging formula for approximation (Piegl and Tiller eq.
//! 9.68, distinct from the interpolation averaging in `fit.rs` because
//! here the control count is smaller than the data count), and the
//! resulting Vandermonde-like system is solved by `axiolid-numeric`'s
//! column-pivoted QR rather than by normal equations: normal equations
//! square the condition number, which is already poor for high degree
//! or clustered parameters.
//!
//! Endpoint interpolation, when requested, is an equality constraint
//! (`P_0 = Q_0`, `P_n = Q_m`) solved by `constrained_least_squares`,
//! not a free consequence of the clamped basis: a plain least-squares
//! solve does not generally reproduce the first and last points exactly.
//!
//! Every reported deviation is the actual Euclidean distance between an
//! input point and the fitted curve or surface evaluated at that
//! point's own parameter -- computed from the solved control points,
//! never estimated from the residual norm the solver reports (which is
//! in the `f64` sense, not the geometric one, once a smoothing term is
//! mixed into the same system).

use axiolid_contracts::GeomError;
use axiolid_core::{Point3, Scalar};
use axiolid_curve::{BSplineCurve3, KnotSpec};
use axiolid_numeric::{constrained_least_squares, least_squares, Matrix, NumericError};
use axiolid_surface::BSplineSurface;

use crate::fit::{basis_at, chord_parameters, collapse, span_of};

/// Result of a geometry solve: `GeomError` on refusal.
pub type FitResult<T> = Result<T, GeomError>;

/// How a point sequence is mapped onto its curve parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Parameterisation {
    /// Equal parameter spacing, independent of point spacing.
    Uniform,
    /// Spacing proportional to the distance between consecutive points
    /// (Piegl and Tiller eq. 9.4).
    ChordLength,
    /// Spacing proportional to the square root of the distance (Piegl
    /// and Tiller eq. 9.6): tends to track curvature changes better than
    /// chord length when the input has sharp corners.
    Centripetal,
}

/// Options for `fit_curve3`.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveFitOptions {
    /// Degree of the fitted curve. Must be at least 1.
    pub degree: u16,
    /// Number of control points to solve for. Must be strictly between
    /// `degree` and the number of input points (inclusive of equality,
    /// which degenerates to interpolation).
    pub control_point_count: usize,
    /// Parameter assignment.
    pub parameterisation: Parameterisation,
    /// Constrain the first and last control points to equal the first
    /// and last input points exactly, by equality-constrained least
    /// squares rather than as an incidental consequence of the basis.
    pub interpolate_endpoints: bool,
    /// Weight of an optional second-difference smoothing term on the
    /// control polygon (a discrete fairness energy). Zero disables it.
    /// Must be finite and non-negative.
    pub smoothing: Scalar,
    /// Newton parameter-correction passes (Hoschek and Lasser, section
    /// 9.3) run after the initial solve: each pass projects every point
    /// onto the current curve and refits with the corrected parameters.
    /// Zero skips correction and uses the assigned parameters as-is.
    pub parameter_correction_iterations: usize,
}

impl Default for CurveFitOptions {
    fn default() -> Self {
        Self {
            degree: 3,
            control_point_count: 0,
            parameterisation: Parameterisation::ChordLength,
            interpolate_endpoints: false,
            smoothing: 0.0,
            parameter_correction_iterations: 0,
        }
    }
}

/// A fitted curve together with its certified deviation from the input.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveFit3 {
    /// The fitted curve.
    pub curve: BSplineCurve3,
    /// Largest Euclidean distance from an input point to the curve at
    /// that point's own parameter.
    pub max_deviation: Scalar,
    /// Root-mean-square of those distances.
    pub rms_deviation: Scalar,
}

/// Options for `fit_curve3_to_tolerance`.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveFitTolerance {
    /// Degree of the fitted curve.
    pub degree: u16,
    /// Parameter assignment.
    pub parameterisation: Parameterisation,
    /// Constrain the endpoints, as in `CurveFitOptions`.
    pub interpolate_endpoints: bool,
    /// Maximum deviation (see `CurveFit3::max_deviation`) the caller
    /// will accept.
    pub max_deviation: Scalar,
    /// Largest control-point count the search may try before refusing.
    pub max_control_points: usize,
}

fn validate_finite(points: &[Point3]) -> FitResult<()> {
    if !points.iter().all(|p| p.is_finite()) {
        return Err(GeomError::InvalidInput(
            "fitting points must be finite".to_owned(),
        ));
    }
    Ok(())
}

fn validate_degree(degree: u16) -> FitResult<usize> {
    if degree == 0 {
        return Err(GeomError::InvalidInput(
            "fitting degree must be at least 1".to_owned(),
        ));
    }
    Ok(usize::from(degree))
}

/// Parameters in `[0, 1]` for `points` under `mode`.
///
/// Uniform spacing never inspects point positions, so it cannot detect
/// coincident points; both chord-length and centripetal spacing refuse
/// coincident consecutive points, since the step they divide by is zero.
pub(crate) fn parameters_for(points: &[Point3], mode: Parameterisation) -> FitResult<Vec<Scalar>> {
    match mode {
        Parameterisation::Uniform => {
            let last = points.len() - 1;
            Ok((0..points.len())
                .map(|k| k as Scalar / last as Scalar)
                .collect())
        }
        Parameterisation::ChordLength => chord_parameters(points),
        Parameterisation::Centripetal => centripetal_parameters(points),
    }
}

fn centripetal_parameters(points: &[Point3]) -> FitResult<Vec<Scalar>> {
    let mut distances = Vec::with_capacity(points.len());
    distances.push(0.0);
    let mut total = 0.0;
    for pair in points.windows(2) {
        let step = (pair[1] - pair[0]).length();
        if step <= 0.0 {
            return Err(GeomError::Degenerate(
                "consecutive fitting points coincide, so centripetal spacing is undefined"
                    .to_owned(),
            ));
        }
        total += step.sqrt();
        distances.push(total);
    }
    Ok(distances.into_iter().map(|d| d / total).collect())
}

/// The clamped, expanded approximation knot vector (Piegl and Tiller eq.
/// 9.68): `control_count - degree - 1` interior knots placed by
/// averaging `degree` consecutive parameters, scaled by the ratio of
/// data points to control points.
pub(crate) fn approximation_knots(
    parameters: &[Scalar],
    degree: usize,
    control_count: usize,
) -> Vec<Scalar> {
    let m = parameters.len() - 1;
    let n = control_count - 1;
    let d = (m + 1) as Scalar / (n - degree + 1) as Scalar;
    let mut knots = vec![0.0; degree + 1];
    for j in 1..=(n - degree) {
        let temp = j as Scalar * d;
        let i = temp.floor() as usize;
        let alpha = temp - i as Scalar;
        let value = if i == 0 {
            parameters[0]
        } else {
            (1.0 - alpha) * parameters[i - 1] + alpha * parameters[i]
        };
        knots.push(value);
    }
    knots.extend(core::iter::repeat_n(1.0, degree + 1));
    knots
}

/// The design matrix: row `k` holds the `degree + 1` non-zero basis
/// values at `parameters[k]`, placed at their control-point columns.
fn design_matrix(
    parameters: &[Scalar],
    expanded_knots: &[Scalar],
    degree: usize,
    control_count: usize,
) -> Matrix {
    let n = control_count - 1;
    let mut dense = vec![0.0; parameters.len() * control_count];
    for (row, &t) in parameters.iter().enumerate() {
        let span = span_of(expanded_knots, n, degree, t);
        let basis = basis_at(span, t, degree, expanded_knots);
        for (offset, value) in basis.iter().enumerate() {
            dense[row * control_count + span - degree + offset] = *value;
        }
    }
    Matrix::from_fn(parameters.len(), control_count, |i, j| {
        dense[i * control_count + j]
    })
}

/// Evaluate a clamped B-spline curve described by its expanded knot
/// vector and control points at `t`.
fn evaluate(
    control_points: &[Point3],
    expanded_knots: &[Scalar],
    degree: usize,
    t: Scalar,
) -> Point3 {
    let n = control_points.len() - 1;
    let span = span_of(expanded_knots, n, degree, t);
    let basis = basis_at(span, t, degree, expanded_knots);
    let mut acc = Point3::new(0.0, 0.0, 0.0);
    for (offset, value) in basis.iter().enumerate() {
        let cp = control_points[span - degree + offset];
        acc = Point3::new(
            acc.x + value * cp.x,
            acc.y + value * cp.y,
            acc.z + value * cp.z,
        );
    }
    acc
}

/// First derivative of a clamped B-spline curve at `t`, by finite
/// difference of neighbouring control points (Piegl and Tiller eq. 3.3),
/// used only to take a Newton step in parameter correction.
fn evaluate_derivative(
    control_points: &[Point3],
    expanded_knots: &[Scalar],
    degree: usize,
    t: Scalar,
) -> axiolid_core::Vec3 {
    if degree == 0 {
        return axiolid_core::Vec3::new(0.0, 0.0, 0.0);
    }
    let n = control_points.len() - 1;
    let span = span_of(expanded_knots, n, degree, t);
    let basis = basis_at(span, t, degree - 1, expanded_knots);
    let mut acc = axiolid_core::Vec3::new(0.0, 0.0, 0.0);
    let p = degree as Scalar;
    for (offset, value) in basis.iter().enumerate() {
        let i = span - degree + offset;
        let denom = expanded_knots[i + degree + 1] - expanded_knots[i + 1];
        if denom == 0.0 {
            continue;
        }
        let coeff = p / denom * value;
        acc += (control_points[i + 1] - control_points[i]) * coeff;
    }
    acc
}

fn from_numeric(err: NumericError) -> GeomError {
    match err {
        NumericError::RankDeficient {
            name,
            rank,
            required,
        } => GeomError::Degenerate(format!(
            "{name} has numerical rank {rank}, needs {required} for a unique least-squares fit"
        )),
        other => GeomError::Degenerate(other.to_string()),
    }
}

/// Solve the least-squares (optionally endpoint-constrained, optionally
/// smoothed) system for one coordinate.
#[allow(clippy::too_many_arguments)]
fn solve_coordinate(
    design: &Matrix,
    smoothing_rows: Option<&Matrix>,
    values: &[Scalar],
    control_count: usize,
    interpolate_endpoints: bool,
    first: Scalar,
    last: Scalar,
) -> FitResult<Vec<Scalar>> {
    let (a, b): (Matrix, Vec<Scalar>) = if let Some(extra) = smoothing_rows {
        let rows = design.rows() + extra.rows();
        let a = Matrix::from_fn(rows, control_count, |i, j| {
            if i < design.rows() {
                design[(i, j)]
            } else {
                extra[(i - design.rows(), j)]
            }
        });
        let mut b = values.to_vec();
        b.resize(rows, 0.0);
        (a, b)
    } else {
        (design.clone(), values.to_vec())
    };

    if interpolate_endpoints {
        let mut c = Matrix::zeros(2, control_count);
        c[(0, 0)] = 1.0;
        c[(1, control_count - 1)] = 1.0;
        let d = [first, last];
        let solution = constrained_least_squares(&a, &b, &c, &d).map_err(from_numeric)?;
        Ok(solution.x)
    } else {
        let solution = least_squares(&a, &b).map_err(from_numeric)?;
        Ok(solution.x)
    }
}

/// Second-difference smoothing rows, scaled by `sqrt(weight)`: row `k`
/// is `weight^(1/2) * (P_{k-1} - 2 P_k + P_{k+1})`, driven to zero.
fn smoothing_rows(control_count: usize, weight: Scalar) -> Option<Matrix> {
    if weight <= 0.0 || control_count < 3 {
        return None;
    }
    let scale = weight.sqrt();
    let rows = control_count - 2;
    Some(Matrix::from_fn(rows, control_count, |i, j| {
        if j == i {
            scale
        } else if j == i + 1 {
            -2.0 * scale
        } else if j == i + 2 {
            scale
        } else {
            0.0
        }
    }))
}

fn assemble_curve(
    control_points: Vec<Point3>,
    expanded_knots: &[Scalar],
    degree: usize,
) -> FitResult<BSplineCurve3> {
    let (knots, multiplicities) = collapse(expanded_knots);
    Ok(BSplineCurve3 {
        degree: u16::try_from(degree)
            .map_err(|_| GeomError::InvalidInput("degree overflows".to_owned()))?,
        control_points,
        knots,
        multiplicities,
        weights: None,
        knot_spec: KnotSpec::Unspecified,
        closed: false,
        self_intersect: None,
    })
}

fn deviation(
    points: &[Point3],
    curve: &[Point3],
    knots: &[Scalar],
    degree: usize,
    t: &[Scalar],
) -> (Scalar, Scalar) {
    let mut max = 0.0;
    let mut sum_sq = 0.0;
    for (point, &tk) in points.iter().zip(t) {
        let fitted = evaluate(curve, knots, degree, tk);
        let d = (fitted - *point).length();
        max = Scalar::max(max, d);
        sum_sq += d * d;
    }
    let rms = (sum_sq / points.len() as Scalar).sqrt();
    (max, rms)
}

/// Fit a B-spline curve of `options.degree` through `points` by least
/// squares, with `options.control_point_count` unknowns.
///
/// # Errors
///
/// Refuses non-finite points, fewer than two points, an invalid degree,
/// a control-point count outside `degree + 1 ..= points.len()`,
/// coincident consecutive points under chord-length or centripetal
/// parameterisation, a non-finite or negative smoothing weight, and a
/// rank-deficient least-squares system (reported with the numeric rank
/// and the rank required).
pub fn fit_curve3(points: &[Point3], options: &CurveFitOptions) -> FitResult<CurveFit3> {
    if points.len() < 2 {
        return Err(GeomError::InvalidInput(
            "fitting needs at least two points".to_owned(),
        ));
    }
    validate_finite(points)?;
    let degree = validate_degree(options.degree)?;
    let control_count = options.control_point_count;
    if control_count < degree + 1 || control_count > points.len() {
        return Err(GeomError::InvalidInput(format!(
            "control point count {control_count} must be between {} and {} for degree {degree} \
             and {} points",
            degree + 1,
            points.len(),
            points.len()
        )));
    }
    if !options.smoothing.is_finite() || options.smoothing < 0.0 {
        return Err(GeomError::InvalidInput(
            "smoothing weight must be finite and non-negative".to_owned(),
        ));
    }

    let mut parameters = parameters_for(points, options.parameterisation)?;
    let mut expanded = approximation_knots(&parameters, degree, control_count);
    let mut design = design_matrix(&parameters, &expanded, degree, control_count);
    let extra = smoothing_rows(control_count, options.smoothing);

    let xs: Vec<Scalar> = points.iter().map(|p| p.x).collect();
    let ys: Vec<Scalar> = points.iter().map(|p| p.y).collect();
    let zs: Vec<Scalar> = points.iter().map(|p| p.z).collect();

    let mut control_points = solve_xyz(
        &design,
        extra.as_ref(),
        &xs,
        &ys,
        &zs,
        control_count,
        options.interpolate_endpoints,
        points[0],
        points[points.len() - 1],
    )?;

    for _ in 0..options.parameter_correction_iterations {
        parameters = correct_parameters(points, &control_points, &expanded, degree, &parameters);
        expanded = approximation_knots(&parameters, degree, control_count);
        design = design_matrix(&parameters, &expanded, degree, control_count);
        let extra = smoothing_rows(control_count, options.smoothing);
        control_points = solve_xyz(
            &design,
            extra.as_ref(),
            &xs,
            &ys,
            &zs,
            control_count,
            options.interpolate_endpoints,
            points[0],
            points[points.len() - 1],
        )?;
    }

    let (max_deviation, rms_deviation) =
        deviation(points, &control_points, &expanded, degree, &parameters);
    let curve = assemble_curve(control_points, &expanded, degree)?;
    Ok(CurveFit3 {
        curve,
        max_deviation,
        rms_deviation,
    })
}

#[allow(clippy::too_many_arguments)]
fn solve_xyz(
    design: &Matrix,
    extra: Option<&Matrix>,
    xs: &[Scalar],
    ys: &[Scalar],
    zs: &[Scalar],
    control_count: usize,
    interpolate_endpoints: bool,
    first: Point3,
    last: Point3,
) -> FitResult<Vec<Point3>> {
    let x = solve_coordinate(
        design,
        extra,
        xs,
        control_count,
        interpolate_endpoints,
        first.x,
        last.x,
    )?;
    let y = solve_coordinate(
        design,
        extra,
        ys,
        control_count,
        interpolate_endpoints,
        first.y,
        last.y,
    )?;
    let z = solve_coordinate(
        design,
        extra,
        zs,
        control_count,
        interpolate_endpoints,
        first.z,
        last.z,
    )?;
    Ok((0..control_count)
        .map(|i| Point3::new(x[i], y[i], z[i]))
        .collect())
}

/// One Newton step of point inversion per point (Hoschek and Lasser),
/// clamped to stay inside `[0, 1]`: the new parameter for point `q` is
/// `t - (C(t) - q) . C'(t) / |C'(t)|^2`.
fn correct_parameters(
    points: &[Point3],
    control_points: &[Point3],
    expanded_knots: &[Scalar],
    degree: usize,
    parameters: &[Scalar],
) -> Vec<Scalar> {
    let last = parameters.len() - 1;
    let mut out = Vec::with_capacity(points.len());
    for (k, (point, &t)) in points.iter().zip(parameters).enumerate() {
        if k == 0 || k == last {
            out.push(t);
            continue;
        }
        let c = evaluate(control_points, expanded_knots, degree, t);
        let d = evaluate_derivative(control_points, expanded_knots, degree, t);
        let speed_sq = d.x * d.x + d.y * d.y + d.z * d.z;
        let corrected = if speed_sq > 0.0 {
            let diff = c - *point;
            let num = diff.x * d.x + diff.y * d.y + diff.z * d.z;
            (t - num / speed_sq).clamp(0.0, 1.0)
        } else {
            t
        };
        out.push(corrected);
    }
    // Keep correction monotone: a parameter overtaking its neighbour
    // would break span search. Falling back to the previous value is
    // conservative rather than silently reordering the sequence.
    for i in 1..out.len() {
        if out[i] <= out[i - 1] {
            out[i] = parameters[i];
        }
    }
    out
}

/// Fit a curve whose maximum deviation does not exceed `tolerance.max_deviation`,
/// adding control points (starting from `degree + 1`) until it does.
///
/// # Errors
///
/// As `fit_curve3`, plus `GeomError::BudgetExceeded` when
/// `tolerance.max_control_points` is reached without meeting the
/// tolerance: the caller must widen the budget rather than receive a fit
/// that silently misses the requested accuracy.
pub fn fit_curve3_to_tolerance(
    points: &[Point3],
    tolerance: &CurveFitTolerance,
) -> FitResult<CurveFit3> {
    if !tolerance.max_deviation.is_finite() || tolerance.max_deviation < 0.0 {
        return Err(GeomError::InvalidInput(
            "tolerance must be finite and non-negative".to_owned(),
        ));
    }
    let degree = validate_degree(tolerance.degree)?;
    let mut control_count = degree + 1;
    loop {
        if control_count > tolerance.max_control_points || control_count > points.len() {
            return Err(GeomError::BudgetExceeded {
                resource: "control points",
            });
        }
        let options = CurveFitOptions {
            degree: tolerance.degree,
            control_point_count: control_count,
            parameterisation: tolerance.parameterisation,
            interpolate_endpoints: tolerance.interpolate_endpoints,
            smoothing: 0.0,
            parameter_correction_iterations: 0,
        };
        let fit = fit_curve3(points, &options)?;
        if fit.max_deviation <= tolerance.max_deviation {
            return Ok(fit);
        }
        // The next iteration's budget check (above) catches running out of
        // either the caller's ceiling or the data itself, so there is
        // nothing further to check here before trying one more point.
        control_count += 1;
    }
}

/// Options for `fit_surface_grid`.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceFitOptions {
    /// Degree along the first (row, `u`) axis.
    pub u_degree: u16,
    /// Degree along the second (column, `v`) axis.
    pub v_degree: u16,
    /// Control points along `u`.
    pub u_control_count: usize,
    /// Control points along `v`.
    pub v_control_count: usize,
    /// Parameter assignment, applied independently along each axis.
    pub parameterisation: Parameterisation,
}

/// A fitted surface together with its certified deviation from the grid.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceFit {
    /// The fitted surface.
    pub surface: BSplineSurface,
    /// Largest Euclidean distance from a grid point to the surface at
    /// that point's own `(u, v)` parameter.
    pub max_deviation: Scalar,
    /// Root-mean-square of those distances.
    pub rms_deviation: Scalar,
}

/// Fit a tensor-product B-spline surface to a rectangular grid of
/// points, `points[i][j]` at row `i` (the `u` direction) and column `j`
/// (the `v` direction), every row the same length.
///
/// The fit is separable (Piegl and Tiller section 9.4): each row is
/// least-squares fit along `u` sharing one design matrix and one global
/// `u` knot vector (parameters averaged across rows), producing an
/// intermediate `u_control_count x cols` net; each column of that net is
/// then least-squares fit along `v` the same way. The result is not the
/// global two-dimensional optimum a single joint solve would give, but
/// it is deterministic, reuses the already-certified curve fit, and
/// matches how `axiolid-nurbs` already builds surfaces from curve data
/// (`loft_surface`).
///
/// # Errors
///
/// Refuses fewer than two rows or columns, ragged rows, non-finite
/// points, an invalid degree, a control count outside
/// `degree + 1 ..= point count` on either axis, coincident consecutive
/// rows or columns under chord-length or centripetal parameterisation,
/// and a rank-deficient system on either axis.
pub fn fit_surface_grid(
    points: &[Vec<Point3>],
    options: &SurfaceFitOptions,
) -> FitResult<SurfaceFit> {
    let rows = points.len();
    if rows < 2 {
        return Err(GeomError::InvalidInput(
            "surface fitting needs at least two rows".to_owned(),
        ));
    }
    let cols = points[0].len();
    if cols < 2 {
        return Err(GeomError::InvalidInput(
            "surface fitting needs at least two columns".to_owned(),
        ));
    }
    for (i, row) in points.iter().enumerate() {
        if row.len() != cols {
            return Err(GeomError::InvalidInput(format!(
                "row {i} has {} points but row 0 has {cols}; the grid must be rectangular",
                row.len()
            )));
        }
    }
    let u_degree = validate_degree(options.u_degree)?;
    let v_degree = validate_degree(options.v_degree)?;
    if options.u_control_count < u_degree + 1 || options.u_control_count > rows {
        return Err(GeomError::InvalidInput(format!(
            "u control point count {} must be between {} and {rows}",
            options.u_control_count,
            u_degree + 1
        )));
    }
    if options.v_control_count < v_degree + 1 || options.v_control_count > cols {
        return Err(GeomError::InvalidInput(format!(
            "v control point count {} must be between {} and {cols}",
            options.v_control_count,
            v_degree + 1
        )));
    }

    // Average u-parameters across rows (one chord-length/centripetal/
    // uniform assignment per row, then averaged column-wise), as Piegl
    // and Tiller eq. 9.71-9.72 do for interpolation; the same averaging
    // keeps every row sharing one design matrix for approximation.
    let mut u_params = vec![0.0; rows];
    for row in points {
        let per_row = parameters_for(row, options.parameterisation)?;
        for (acc, value) in u_params.iter_mut().zip(&per_row) {
            *acc += value / cols as Scalar;
        }
    }
    u_params[0] = 0.0;
    *u_params.last_mut().expect("rows checked non-empty") = 1.0;

    let mut v_params = vec![0.0; cols];
    for j in 0..cols {
        let column: Vec<Point3> = points.iter().map(|row| row[j]).collect();
        let per_col = parameters_for(&column, options.parameterisation)?;
        for (acc, value) in v_params.iter_mut().zip(&per_col) {
            *acc += value / rows as Scalar;
        }
    }
    v_params[0] = 0.0;
    *v_params.last_mut().expect("columns checked non-empty") = 1.0;

    let u_knots = approximation_knots(&u_params, u_degree, options.u_control_count);
    let v_knots = approximation_knots(&v_params, v_degree, options.v_control_count);

    // Fit each row along u: one design matrix shared by every row and
    // every coordinate.
    let u_design = design_matrix(&u_params, &u_knots, u_degree, options.u_control_count);
    let mut intermediate: Vec<Vec<Point3>> =
        vec![Vec::with_capacity(cols); options.u_control_count];
    for j in 0..cols {
        let column: Vec<Point3> = points.iter().map(|row| row[j]).collect();
        let xs: Vec<Scalar> = column.iter().map(|p| p.x).collect();
        let ys: Vec<Scalar> = column.iter().map(|p| p.y).collect();
        let zs: Vec<Scalar> = column.iter().map(|p| p.z).collect();
        let solved = solve_xyz(
            &u_design,
            None,
            &xs,
            &ys,
            &zs,
            options.u_control_count,
            false,
            column[0],
            column[rows - 1],
        )?;
        for (row, point) in intermediate.iter_mut().zip(solved) {
            row.push(point);
        }
    }

    // Fit each row of the intermediate net along v.
    let v_design = design_matrix(&v_params, &v_knots, v_degree, options.v_control_count);
    let mut net: Vec<Vec<Point3>> = Vec::with_capacity(options.u_control_count);
    for row in &intermediate {
        let xs: Vec<Scalar> = row.iter().map(|p| p.x).collect();
        let ys: Vec<Scalar> = row.iter().map(|p| p.y).collect();
        let zs: Vec<Scalar> = row.iter().map(|p| p.z).collect();
        let solved = solve_xyz(
            &v_design,
            None,
            &xs,
            &ys,
            &zs,
            options.v_control_count,
            false,
            row[0],
            row[cols - 1],
        )?;
        net.push(solved);
    }

    // Deviation: evaluate the tensor-product surface at each grid
    // point's own (u, v) parameter.
    let mut max_deviation: Scalar = 0.0;
    let mut sum_sq = 0.0;
    for (i, row) in points.iter().enumerate() {
        for (j, point) in row.iter().enumerate() {
            let fitted = evaluate_surface(
                &net,
                &u_knots,
                u_degree,
                &v_knots,
                v_degree,
                u_params[i],
                v_params[j],
            );
            let d = (fitted - *point).length();
            max_deviation = Scalar::max(max_deviation, d);
            sum_sq += d * d;
        }
    }
    let rms_deviation = (sum_sq / (rows * cols) as Scalar).sqrt();

    let (u_knots_distinct, u_multiplicities) = collapse(&u_knots);
    let (v_knots_distinct, v_multiplicities) = collapse(&v_knots);
    let surface = BSplineSurface {
        u_degree: options.u_degree,
        v_degree: options.v_degree,
        control_points: net,
        u_knots: u_knots_distinct,
        u_multiplicities,
        v_knots: v_knots_distinct,
        v_multiplicities,
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    };
    Ok(SurfaceFit {
        surface,
        max_deviation,
        rms_deviation,
    })
}

/// Evaluate a tensor-product B-spline surface at `(u, v)` from its
/// expanded knot vectors and control net.
fn evaluate_surface(
    net: &[Vec<Point3>],
    u_knots: &[Scalar],
    u_degree: usize,
    v_knots: &[Scalar],
    v_degree: usize,
    u: Scalar,
    v: Scalar,
) -> Point3 {
    let rows = net.len();
    let cols = net[0].len();
    let u_span = span_of(u_knots, rows - 1, u_degree, u);
    let u_basis = basis_at(u_span, u, u_degree, u_knots);
    let v_span = span_of(v_knots, cols - 1, v_degree, v);
    let v_basis = basis_at(v_span, v, v_degree, v_knots);
    let mut acc = Point3::new(0.0, 0.0, 0.0);
    for (ui, &uw) in u_basis.iter().enumerate() {
        let row = u_span - u_degree + ui;
        for (vi, &vw) in v_basis.iter().enumerate() {
            let col = v_span - v_degree + vi;
            let weight = uw * vw;
            let cp = net[row][col];
            acc = Point3::new(
                acc.x + weight * cp.x,
                acc.y + weight * cp.y,
                acc.z + weight * cp.z,
            );
        }
    }
    acc
}
