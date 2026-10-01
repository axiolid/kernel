//! Dense linear algebra for small and medium systems.
//!
//! - [`Lu`]: Gaussian elimination with partial pivoting, for square
//!   systems (OCCT `math_Gauss`).
//! - [`Cholesky`]: `L L^T` for symmetric positive definite systems, such as
//!   normal equations and fairing energies (OCCT `math_Crout` covers the
//!   symmetric case there).
//! - [`Qr`]: Householder QR with column pivoting, for linear least squares
//!   and rank detection (OCCT `math_Householder`, LAPACK `xGEQP3`).
//! - [`constrained_least_squares`]: least squares subject to linear
//!   equality constraints, by the null-space method (Golub and Van Loan,
//!   *Matrix Computations*, 6.2).
//!
//! Every solve reports a condition estimate (Hager's 1-norm estimator as
//! refined by Higham, LAPACK `xLACON`) and a relative forward error estimate
//! derived from it. Singular and rank-deficient systems are refused by name
//! rather than solved into infinities.

// Triangular solves read clearest with the textbook indices.
#![allow(clippy::needless_range_loop)]

use std::ops::{Index, IndexMut};

use crate::error::{all_finite, NumericError, NumericResult};

const U: f64 = f64::EPSILON / 2.0;

/// A dense row-major matrix of `f64`.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl Matrix {
    /// A `rows x cols` matrix of zeros.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    /// The `n x n` identity.
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m[(i, i)] = 1.0;
        }
        m
    }

    /// A matrix whose entry `(i, j)` is `f(i, j)`.
    pub fn from_fn(rows: usize, cols: usize, mut f: impl FnMut(usize, usize) -> f64) -> Self {
        let mut m = Self::zeros(rows, cols);
        for i in 0..rows {
            for j in 0..cols {
                m[(i, j)] = f(i, j);
            }
        }
        m
    }

    /// A matrix from rows of equal length.
    ///
    /// # Errors
    ///
    /// Refuses rows of different lengths.
    pub fn from_rows(rows: &[&[f64]]) -> NumericResult<Self> {
        let cols = rows.first().map_or(0, |r| r.len());
        let mut data = Vec::with_capacity(rows.len() * cols);
        for r in rows {
            if r.len() != cols {
                return Err(NumericError::DimensionMismatch {
                    name: "matrix row",
                    expected: cols,
                    found: r.len(),
                });
            }
            data.extend_from_slice(r);
        }
        Ok(Self {
            rows: rows.len(),
            cols,
            data,
        })
    }

    /// Number of rows.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Row `i` as a slice.
    pub fn row(&self, i: usize) -> &[f64] {
        &self.data[i * self.cols..(i + 1) * self.cols]
    }

    /// The transpose.
    pub fn transpose(&self) -> Self {
        Self::from_fn(self.cols, self.rows, |i, j| self[(j, i)])
    }

    /// `A x`.
    ///
    /// # Errors
    ///
    /// Refuses `x` whose length is not the column count.
    pub fn mul_vec(&self, x: &[f64]) -> NumericResult<Vec<f64>> {
        if x.len() != self.cols {
            return Err(NumericError::DimensionMismatch {
                name: "vector",
                expected: self.cols,
                found: x.len(),
            });
        }
        Ok((0..self.rows)
            .map(|i| self.row(i).iter().zip(x).map(|(a, b)| a * b).sum())
            .collect())
    }

    /// The 1-norm: the largest column sum of absolute values.
    pub fn norm1(&self) -> f64 {
        (0..self.cols)
            .map(|j| (0..self.rows).map(|i| self[(i, j)].abs()).sum::<f64>())
            .fold(0.0, f64::max)
    }

    fn require_finite(&self, name: &'static str) -> NumericResult<()> {
        all_finite(&self.data, name)
    }

    fn require_square(&self, name: &'static str) -> NumericResult<()> {
        if self.rows == self.cols {
            Ok(())
        } else {
            Err(NumericError::DimensionMismatch {
                name,
                expected: self.rows,
                found: self.cols,
            })
        }
    }
}

impl Index<(usize, usize)> for Matrix {
    type Output = f64;
    fn index(&self, (i, j): (usize, usize)) -> &f64 {
        assert!(i < self.rows && j < self.cols, "matrix index out of range");
        &self.data[i * self.cols + j]
    }
}

impl IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut f64 {
        assert!(i < self.rows && j < self.cols, "matrix index out of range");
        &mut self.data[i * self.cols + j]
    }
}

/// The solution of a square linear system.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct LinearSolution {
    /// The solution vector.
    pub x: Vec<f64>,
    /// `||b - A x||_1`, computed in `f64`.
    pub residual_norm: f64,
    /// Estimate of `||x - x_exact||_1 / ||x_exact||_1`: twice the condition
    /// estimate times the normwise backward error, the latter enlarged by
    /// a bound on the rounding of the computed residual itself.
    pub relative_error_estimate: f64,
}

/// The solution of a linear least-squares problem.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct LeastSquaresSolution {
    /// The minimiser of `||A x - b||_2`.
    pub x: Vec<f64>,
    /// `||A x - b||_2` at the solution.
    pub residual_norm: f64,
    /// Estimate of `||x - x_exact|| / ||x_exact||` from the standard
    /// first-order perturbation bound for Householder least squares,
    /// `m n u (2 kappa + kappa^2 ||r|| / (||A|| ||x||))`.
    pub relative_error_estimate: f64,
}

/// Estimate `||M^-1||_1` given solves with `M` and `M^T` (Hager 1984, with
/// Higham's 1988 alternative vector).
fn inverse_norm1_estimate(
    n: usize,
    mut solve: impl FnMut(&mut [f64]),
    mut solve_transpose: impl FnMut(&mut [f64]),
) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let mut x = vec![1.0 / n as f64; n];
    let mut estimate = 0.0f64;
    let mut last_j = usize::MAX;
    for _ in 0..5 {
        let mut y = x.clone();
        solve(&mut y);
        estimate = estimate.max(y.iter().map(|v| v.abs()).sum());
        let mut z: Vec<f64> = y
            .iter()
            .map(|v| if *v >= 0.0 { 1.0 } else { -1.0 })
            .collect();
        solve_transpose(&mut z);
        let (j, zmax) = z
            .iter()
            .enumerate()
            .map(|(i, v)| (i, v.abs()))
            .fold((0, -1.0), |best, c| if c.1 > best.1 { c } else { best });
        let ztx: f64 = z.iter().zip(&x).map(|(a, b)| a * b).sum();
        if zmax <= ztx || j == last_j {
            break;
        }
        last_j = j;
        x = vec![0.0; n];
        x[j] = 1.0;
    }
    // Higham's alternating vector catches cases the power-like loop misses.
    let mut alt: Vec<f64> = (0..n)
        .map(|i| {
            let s = if i % 2 == 0 { 1.0 } else { -1.0 };
            s * (1.0 + i as f64 / (n.max(2) - 1) as f64)
        })
        .collect();
    solve(&mut alt);
    let alt_estimate = 2.0 * alt.iter().map(|v| v.abs()).sum::<f64>() / (3.0 * n as f64);
    estimate.max(alt_estimate)
}

fn forward_error_estimate(a: &Matrix, x: &[f64], b: &[f64], condition: f64) -> (f64, f64) {
    let n = a.cols;
    let mut r_norm = 0.0;
    let mut scale = 0.0;
    for (i, bi) in b.iter().enumerate() {
        let row = a.row(i);
        let ax: f64 = row.iter().zip(x).map(|(p, q)| p * q).sum();
        let mag: f64 = row.iter().zip(x).map(|(p, q)| (p * q).abs()).sum::<f64>() + bi.abs();
        r_norm += (bi - ax).abs();
        scale += mag;
    }
    let x_norm: f64 = x.iter().map(|v| v.abs()).sum();
    let b_norm: f64 = b.iter().map(|v| v.abs()).sum();
    let denominator = a.norm1() * x_norm + b_norm;
    if denominator == 0.0 {
        return (r_norm, 0.0);
    }
    // The computed residual carries its own rounding, up to
    // gamma(n + 1) * sum(|A||x| + |b|).
    let k = (n + 1) as f64;
    let residual_rounding = k * U / (1.0 - k * U) * scale;
    let eta = (r_norm + residual_rounding) / denominator;
    let kappa_eta = condition * eta;
    let estimate = if kappa_eta < 0.5 {
        2.0 * kappa_eta / (1.0 - kappa_eta)
    } else {
        f64::INFINITY
    };
    (r_norm, estimate)
}

/// `P A = L U` with partial pivoting.
#[derive(Debug, Clone, PartialEq)]
pub struct Lu {
    a: Matrix,
    lu: Matrix,
    perm: Vec<usize>,
    condition: f64,
}

impl Lu {
    /// Factor a square matrix.
    ///
    /// # Errors
    ///
    /// Refuses a non-square or non-finite matrix, and a matrix singular to
    /// working precision (a zero pivot, or a condition estimate above
    /// `1 / eps`).
    pub fn new(a: &Matrix) -> NumericResult<Self> {
        a.require_square("system matrix")?;
        a.require_finite("system matrix")?;
        let n = a.rows;
        let mut lu = a.clone();
        let mut perm: Vec<usize> = (0..n).collect();
        for k in 0..n {
            let p = (k..n)
                .max_by(|&i, &j| lu[(i, k)].abs().total_cmp(&lu[(j, k)].abs()))
                .unwrap_or(k);
            if lu[(p, k)] == 0.0 {
                return Err(NumericError::Singular {
                    name: "system matrix",
                });
            }
            if p != k {
                for j in 0..n {
                    lu.data.swap(k * n + j, p * n + j);
                }
                perm.swap(k, p);
            }
            let pivot = lu[(k, k)];
            for i in (k + 1)..n {
                let factor = lu[(i, k)] / pivot;
                lu[(i, k)] = factor;
                if factor != 0.0 {
                    for j in (k + 1)..n {
                        let v = lu[(k, j)];
                        lu[(i, j)] -= factor * v;
                    }
                }
            }
        }
        let mut out = Self {
            a: a.clone(),
            lu,
            perm,
            condition: 0.0,
        };
        let inverse_norm = inverse_norm1_estimate(
            n,
            |v| out.solve_in_place(v),
            |v| out.solve_transpose_in_place(v),
        );
        out.condition = a.norm1() * inverse_norm;
        if !out.condition.is_finite() || out.condition * f64::EPSILON >= 1.0 {
            return Err(NumericError::Singular {
                name: "system matrix",
            });
        }
        Ok(out)
    }

    /// Estimate of the 1-norm condition number `||A||_1 ||A^-1||_1`. It is
    /// a lower bound, in practice within a factor of 3.
    pub fn condition_estimate(&self) -> f64 {
        self.condition
    }

    /// Solve `A x = b`.
    ///
    /// # Errors
    ///
    /// Refuses `b` of the wrong length or with non-finite entries.
    pub fn solve(&self, b: &[f64]) -> NumericResult<LinearSolution> {
        check_rhs(b, self.a.rows)?;
        let mut x = b.to_vec();
        self.solve_in_place(&mut x);
        let (residual_norm, relative_error_estimate) =
            forward_error_estimate(&self.a, &x, b, self.condition);
        Ok(LinearSolution {
            x,
            residual_norm,
            relative_error_estimate,
        })
    }

    fn solve_in_place(&self, b: &mut [f64]) {
        let n = self.a.rows;
        let mut y: Vec<f64> = self.perm.iter().map(|&p| b[p]).collect();
        for i in 0..n {
            let mut s = y[i];
            for j in 0..i {
                s -= self.lu[(i, j)] * y[j];
            }
            y[i] = s;
        }
        for i in (0..n).rev() {
            let mut s = y[i];
            for j in (i + 1)..n {
                s -= self.lu[(i, j)] * y[j];
            }
            y[i] = s / self.lu[(i, i)];
        }
        b.copy_from_slice(&y);
    }

    fn solve_transpose_in_place(&self, c: &mut [f64]) {
        let n = self.a.rows;
        // A^T = U^T L^T P: solve U^T w = c, L^T v = w, then z[perm[i]] = v[i].
        let mut w = c.to_vec();
        for i in 0..n {
            let mut s = w[i];
            for j in 0..i {
                s -= self.lu[(j, i)] * w[j];
            }
            w[i] = s / self.lu[(i, i)];
        }
        for i in (0..n).rev() {
            let mut s = w[i];
            for j in (i + 1)..n {
                s -= self.lu[(j, i)] * w[j];
            }
            w[i] = s;
        }
        for (i, &p) in self.perm.iter().enumerate() {
            c[p] = w[i];
        }
    }
}

fn check_rhs(b: &[f64], expected: usize) -> NumericResult<()> {
    if b.len() != expected {
        return Err(NumericError::DimensionMismatch {
            name: "right-hand side",
            expected,
            found: b.len(),
        });
    }
    all_finite(b, "right-hand side")
}

/// `A = L L^T` for a symmetric positive definite matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct Cholesky {
    a: Matrix,
    l: Matrix,
    condition: f64,
}

impl Cholesky {
    /// Factor a symmetric positive definite matrix.
    ///
    /// Symmetry is checked to `n eps max|a_ij|`; only the lower triangle is
    /// used.
    ///
    /// # Errors
    ///
    /// Refuses a non-square or non-finite matrix, an asymmetric one
    /// ([`NumericError::NotSymmetric`]), one with a non-positive pivot
    /// ([`NumericError::NotPositiveDefinite`]), and one singular to working
    /// precision.
    pub fn new(a: &Matrix) -> NumericResult<Self> {
        a.require_square("symmetric matrix")?;
        a.require_finite("symmetric matrix")?;
        let n = a.rows;
        let max_abs = a.data.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let symmetry_tolerance = n as f64 * f64::EPSILON * max_abs;
        for i in 0..n {
            for j in 0..i {
                if (a[(i, j)] - a[(j, i)]).abs() > symmetry_tolerance {
                    return Err(NumericError::NotSymmetric { row: i, column: j });
                }
            }
        }
        let mut l = Matrix::zeros(n, n);
        for j in 0..n {
            let mut d = a[(j, j)];
            for k in 0..j {
                d -= l[(j, k)] * l[(j, k)];
            }
            if d <= 0.0 || !d.is_finite() {
                return Err(NumericError::NotPositiveDefinite { pivot: j });
            }
            let ljj = d.sqrt();
            l[(j, j)] = ljj;
            for i in (j + 1)..n {
                let mut s = a[(i, j)];
                for k in 0..j {
                    s -= l[(i, k)] * l[(j, k)];
                }
                l[(i, j)] = s / ljj;
            }
        }
        let mut out = Self {
            a: a.clone(),
            l,
            condition: 0.0,
        };
        let inverse_norm =
            inverse_norm1_estimate(n, |v| out.solve_in_place(v), |v| out.solve_in_place(v));
        out.condition = a.norm1() * inverse_norm;
        if !out.condition.is_finite() || out.condition * f64::EPSILON >= 1.0 {
            return Err(NumericError::Singular {
                name: "symmetric matrix",
            });
        }
        Ok(out)
    }

    /// Estimate of the 1-norm condition number.
    pub fn condition_estimate(&self) -> f64 {
        self.condition
    }

    /// The factor `L` (lower triangular, positive diagonal).
    pub fn factor(&self) -> &Matrix {
        &self.l
    }

    /// Solve `A x = b`.
    ///
    /// # Errors
    ///
    /// Refuses `b` of the wrong length or with non-finite entries.
    pub fn solve(&self, b: &[f64]) -> NumericResult<LinearSolution> {
        check_rhs(b, self.a.rows)?;
        let mut x = b.to_vec();
        self.solve_in_place(&mut x);
        let (residual_norm, relative_error_estimate) =
            forward_error_estimate(&self.a, &x, b, self.condition);
        Ok(LinearSolution {
            x,
            residual_norm,
            relative_error_estimate,
        })
    }

    fn solve_in_place(&self, b: &mut [f64]) {
        let n = self.a.rows;
        for i in 0..n {
            let mut s = b[i];
            for k in 0..i {
                s -= self.l[(i, k)] * b[k];
            }
            b[i] = s / self.l[(i, i)];
        }
        for i in (0..n).rev() {
            let mut s = b[i];
            for k in (i + 1)..n {
                s -= self.l[(k, i)] * b[k];
            }
            b[i] = s / self.l[(i, i)];
        }
    }
}

/// Householder QR with column pivoting: `A P = Q R`.
#[derive(Debug, Clone, PartialEq)]
pub struct Qr {
    a: Matrix,
    /// Upper triangle holds `R`; the rest is scratch.
    r: Matrix,
    /// Householder vectors `v_k` (length `m - k`) and their `beta_k`.
    reflectors: Vec<(Vec<f64>, f64)>,
    perm: Vec<usize>,
    rank: usize,
    condition: f64,
}

impl Qr {
    /// Factor an `m x n` matrix and determine its numerical rank: a pivot
    /// counts when it exceeds `max(m, n) eps |r_00|`.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite matrix.
    pub fn new(a: &Matrix) -> NumericResult<Self> {
        a.require_finite("least-squares matrix")?;
        let (m, n) = (a.rows, a.cols);
        let mut r = a.clone();
        let mut perm: Vec<usize> = (0..n).collect();
        let mut reflectors = Vec::new();
        let steps = m.min(n);
        let mut rank = steps;
        let mut first_pivot = 0.0;
        for k in 0..steps {
            // Pivot: the remaining column of largest norm.
            let norms: Vec<f64> = (k..n)
                .map(|j| (k..m).map(|i| r[(i, j)] * r[(i, j)]).sum::<f64>().sqrt())
                .collect();
            let (offset, &best) = norms
                .iter()
                .enumerate()
                .max_by(|x, y| x.1.total_cmp(y.1))
                .expect("k < n");
            if k == 0 {
                first_pivot = best;
            }
            if best <= (m.max(n) as f64) * f64::EPSILON * first_pivot || best == 0.0 {
                rank = k;
                break;
            }
            let p = k + offset;
            if p != k {
                for i in 0..m {
                    r.data.swap(i * n + k, i * n + p);
                }
                perm.swap(k, p);
            }
            let x0 = r[(k, k)];
            let alpha = if x0 >= 0.0 { -best } else { best };
            let mut v: Vec<f64> = (k..m).map(|i| r[(i, k)]).collect();
            v[0] -= alpha;
            let vtv: f64 = v.iter().map(|t| t * t).sum();
            let beta = if vtv == 0.0 { 0.0 } else { 2.0 / vtv };
            for j in k..n {
                let dot: f64 = (k..m).map(|i| v[i - k] * r[(i, j)]).sum();
                let s = beta * dot;
                for i in k..m {
                    r[(i, j)] -= s * v[i - k];
                }
            }
            r[(k, k)] = alpha;
            for i in (k + 1)..m {
                r[(i, k)] = 0.0;
            }
            reflectors.push((v, beta));
        }
        let mut out = Self {
            a: a.clone(),
            r,
            reflectors,
            perm,
            rank,
            condition: 0.0,
        };
        let k = out.rank;
        let inverse_norm =
            inverse_norm1_estimate(k, |v| out.solve_r(v), |v| out.solve_r_transpose(v));
        let r_norm = (0..k)
            .map(|j| (0..=j).map(|i| out.r[(i, j)].abs()).sum::<f64>())
            .fold(0.0, f64::max);
        out.condition = r_norm * inverse_norm;
        Ok(out)
    }

    /// Numerical rank.
    pub fn rank(&self) -> usize {
        self.rank
    }

    /// Estimate of the 1-norm condition number of the leading `rank x rank`
    /// block of `R`, which equals `A`'s 2-norm condition number up to a
    /// factor of `n`.
    pub fn condition_estimate(&self) -> f64 {
        self.condition
    }

    /// Minimise `||A x - b||_2`.
    ///
    /// # Errors
    ///
    /// Refuses `b` of the wrong length or with non-finite entries, and a
    /// matrix whose rank is below its column count
    /// ([`NumericError::RankDeficient`]): the minimiser is then not unique.
    pub fn solve_least_squares(&self, b: &[f64]) -> NumericResult<LeastSquaresSolution> {
        let (m, n) = (self.a.rows, self.a.cols);
        check_rhs(b, m)?;
        if self.rank < n {
            return Err(NumericError::RankDeficient {
                name: "least-squares matrix",
                rank: self.rank,
                required: n,
            });
        }
        let mut c = b.to_vec();
        self.apply_qt(&mut c);
        let mut z = c[..n].to_vec();
        self.solve_r(&mut z);
        let mut x = vec![0.0; n];
        for (k, &p) in self.perm.iter().enumerate() {
            x[p] = z[k];
        }
        let ax = self.a.mul_vec(&x)?;
        let residual_norm = ax
            .iter()
            .zip(b)
            .map(|(p, q)| (p - q) * (p - q))
            .sum::<f64>()
            .sqrt();
        let a_norm = self.a.data.iter().map(|v| v * v).sum::<f64>().sqrt();
        let x_norm = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        let kappa = self.condition;
        let relative_error_estimate = if x_norm == 0.0 {
            0.0
        } else {
            (m * n) as f64 * U * (2.0 * kappa + kappa * kappa * residual_norm / (a_norm * x_norm))
        };
        Ok(LeastSquaresSolution {
            x,
            residual_norm,
            relative_error_estimate,
        })
    }

    /// `c <- Q^T c` (length `m`).
    fn apply_qt(&self, c: &mut [f64]) {
        for (k, (v, beta)) in self.reflectors.iter().enumerate() {
            reflect(&mut c[k..], v, *beta);
        }
    }

    /// `c <- Q c` (length `m`).
    fn apply_q(&self, c: &mut [f64]) {
        for (k, (v, beta)) in self.reflectors.iter().enumerate().rev() {
            reflect(&mut c[k..], v, *beta);
        }
    }

    /// Solve `R_11 z = c` for the leading `rank` entries.
    fn solve_r(&self, c: &mut [f64]) {
        let k = self.rank;
        for i in (0..k).rev() {
            let mut s = c[i];
            for j in (i + 1)..k {
                s -= self.r[(i, j)] * c[j];
            }
            c[i] = s / self.r[(i, i)];
        }
    }

    /// Solve `R_11^T z = c`.
    fn solve_r_transpose(&self, c: &mut [f64]) {
        let k = self.rank;
        for i in 0..k {
            let mut s = c[i];
            for j in 0..i {
                s -= self.r[(j, i)] * c[j];
            }
            c[i] = s / self.r[(i, i)];
        }
    }
}

fn reflect(c: &mut [f64], v: &[f64], beta: f64) {
    let dot: f64 = v.iter().zip(c.iter()).map(|(a, b)| a * b).sum();
    let s = beta * dot;
    for (ci, vi) in c.iter_mut().zip(v) {
        *ci -= s * vi;
    }
}

/// Minimise `||A x - b||_2`: a convenience for [`Qr::new`] followed by
/// [`Qr::solve_least_squares`]. Factor once with [`Qr`] when several
/// right-hand sides share `A` (one per coordinate when fitting points).
///
/// # Errors
///
/// As [`Qr::new`] and [`Qr::solve_least_squares`].
pub fn least_squares(a: &Matrix, b: &[f64]) -> NumericResult<LeastSquaresSolution> {
    Qr::new(a)?.solve_least_squares(b)
}

/// The solution of an equality-constrained least-squares problem.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct ConstrainedSolution {
    /// The minimiser of `||A x - b||_2` subject to `C x = d`.
    pub x: Vec<f64>,
    /// `||A x - b||_2` at the solution.
    pub residual_norm: f64,
    /// `||C x - d||_2` at the solution: zero up to rounding.
    pub constraint_residual_norm: f64,
    /// Sum of the relative error estimates of the two stages (constraint
    /// solve and reduced least squares).
    pub relative_error_estimate: f64,
}

/// Minimise `||A x - b||_2` subject to `C x = d` by the null-space method.
///
/// `C` (`p x n`) must have full row rank `p <= n`, and `A` must have full
/// column rank on the null space of `C`, so the minimiser is unique. Use
/// this to fit with interpolated points or prescribed end tangents.
///
/// # Errors
///
/// Refuses mismatched sizes, non-finite input, a rank-deficient constraint
/// matrix (redundant or conflicting constraints) and an objective that does
/// not determine the solution on the constraints' null space.
pub fn constrained_least_squares(
    a: &Matrix,
    b: &[f64],
    c: &Matrix,
    d: &[f64],
) -> NumericResult<ConstrainedSolution> {
    let n = a.cols;
    if c.cols != n {
        return Err(NumericError::DimensionMismatch {
            name: "constraint matrix columns",
            expected: n,
            found: c.cols,
        });
    }
    check_rhs(b, a.rows)?;
    if d.len() != c.rows {
        return Err(NumericError::DimensionMismatch {
            name: "constraint right-hand side",
            expected: c.rows,
            found: d.len(),
        });
    }
    all_finite(d, "constraint right-hand side")?;
    a.require_finite("least-squares matrix")?;
    c.require_finite("constraint matrix")?;
    let p = c.rows;

    // C^T Pi = Q [R; 0].
    let ct = Qr::new(&c.transpose())?;
    if ct.rank < p {
        return Err(NumericError::RankDeficient {
            name: "constraint matrix",
            rank: ct.rank,
            required: p,
        });
    }
    // R^T y1 = Pi^T d.
    let mut y1: Vec<f64> = ct.perm.iter().map(|&i| d[i]).collect();
    ct.solve_r_transpose(&mut y1);

    // A Q, row by row: (A Q)_i = (Q^T a_i^T)^T.
    let mut aq = Matrix::zeros(a.rows, n);
    for i in 0..a.rows {
        let mut row = a.row(i).to_vec();
        ct.apply_qt(&mut row);
        aq.data[i * n..(i + 1) * n].copy_from_slice(&row);
    }
    let mut rhs = b.to_vec();
    for (i, value) in rhs.iter_mut().enumerate() {
        let row = aq.row(i);
        *value -= row[..p].iter().zip(&y1).map(|(s, t)| s * t).sum::<f64>();
    }
    let mut y = y1.clone();
    let mut reduced_estimate = 0.0;
    if p < n {
        let a2 = Matrix::from_fn(a.rows, n - p, |i, j| aq[(i, p + j)]);
        let qr = Qr::new(&a2)?;
        if qr.rank < n - p {
            return Err(NumericError::RankDeficient {
                name: "least-squares matrix on the constraints' null space",
                rank: qr.rank,
                required: n - p,
            });
        }
        let reduced = qr.solve_least_squares(&rhs)?;
        reduced_estimate = reduced.relative_error_estimate;
        y.extend_from_slice(&reduced.x);
    }
    let mut x = y;
    ct.apply_q(&mut x);

    let residual = |m: &Matrix, rhs: &[f64]| -> NumericResult<f64> {
        Ok(m.mul_vec(&x)?
            .iter()
            .zip(rhs)
            .map(|(s, t)| (s - t) * (s - t))
            .sum::<f64>()
            .sqrt())
    };
    let residual_norm = residual(a, b)?;
    let constraint_residual_norm = residual(c, d)?;
    let constraint_estimate = (n * p) as f64 * U * 2.0 * ct.condition;
    Ok(ConstrainedSolution {
        x,
        residual_norm,
        constraint_residual_norm,
        relative_error_estimate: constraint_estimate + reduced_estimate,
    })
}
