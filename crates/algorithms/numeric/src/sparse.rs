//! Sparse matrices and the preconditioned conjugate-gradient solve.
//!
//! Mesh operators (graph and cotangent Laplacians, their squares, mass
//! matrices) have a handful of non-zeros per row, so a dense [`Matrix`]
//! of a few thousand mesh vertices already costs more memory than the mesh
//! itself and an `O(n^3)` factorisation more time than the rest of the
//! algorithm. [`SparseMatrix`] stores only the non-zeros in compressed rows,
//! and [`conjugate_gradient`] solves a symmetric positive definite system
//! with it in `O(iterations * non-zeros)`.
//!
//! The solve follows Hestenes and Stiefel (1952) with a Jacobi (diagonal)
//! preconditioner, the textbook form in Saad, *Iterative Methods for Sparse
//! Linear Systems*, algorithm 9.1. Its stopping test is on the *true*
//! residual `b - A x`, recomputed from the matrix when the recurrence says
//! the tolerance is met, so drift in the recurred residual cannot report a
//! convergence that did not happen.
//!
//! [`Matrix`]: crate::Matrix

use crate::error::{all_finite, tolerance, NumericError, NumericResult, Status};

/// A sparse matrix in compressed sparse row form.
///
/// Rows are stored in order and the column indices within each row are
/// strictly increasing, so two matrices built from the same entries in any
/// order are equal and every product is computed in the same order.
#[derive(Debug, Clone, PartialEq)]
pub struct SparseMatrix {
    rows: usize,
    cols: usize,
    /// `row_start[i]..row_start[i + 1]` indexes row `i` in `columns` and `values`.
    row_start: Vec<usize>,
    columns: Vec<usize>,
    values: Vec<f64>,
}

impl SparseMatrix {
    /// Build from `(row, column, value)` triplets. Entries naming the same
    /// position are summed, in the order given; entries that sum to zero are
    /// kept as explicit zeros.
    ///
    /// # Errors
    ///
    /// [`NumericError::InvalidArgument`] for a triplet outside the
    /// `rows x cols` shape and [`NumericError::NonFiniteInput`] for a
    /// non-finite value.
    pub fn from_triplets(
        rows: usize,
        cols: usize,
        triplets: &[(usize, usize, f64)],
    ) -> NumericResult<Self> {
        let mut order: Vec<usize> = (0..triplets.len()).collect();
        for &(i, j, v) in triplets {
            if i >= rows || j >= cols {
                return Err(NumericError::InvalidArgument {
                    name: "triplet",
                    reason: "must lie inside the matrix shape",
                });
            }
            if !v.is_finite() {
                return Err(NumericError::NonFiniteInput {
                    name: "triplet value",
                });
            }
        }
        // Stable: equal positions keep their input order, so the summation
        // order (and the rounding) is the caller's.
        order.sort_by_key(|&k| (triplets[k].0, triplets[k].1));
        let mut row_start = vec![0usize; rows + 1];
        let mut columns: Vec<usize> = Vec::with_capacity(triplets.len());
        let mut values: Vec<f64> = Vec::with_capacity(triplets.len());
        let mut last: Option<(usize, usize)> = None;
        for k in order {
            let (i, j, v) = triplets[k];
            if last == Some((i, j)) {
                *values.last_mut().expect("an entry precedes its duplicate") += v;
            } else {
                columns.push(j);
                values.push(v);
                row_start[i + 1] += 1;
                last = Some((i, j));
            }
        }
        for i in 0..rows {
            row_start[i + 1] += row_start[i];
        }
        Ok(Self {
            rows,
            cols,
            row_start,
            columns,
            values,
        })
    }

    /// Number of rows.
    #[must_use]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    #[must_use]
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Number of stored entries.
    #[must_use]
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Entry `(i, j)`; zero when it is not stored.
    ///
    /// # Panics
    ///
    /// When `(i, j)` is outside the shape.
    #[must_use]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        assert!(i < self.rows && j < self.cols, "matrix index out of range");
        let (cols, vals) = self.row(i);
        cols.binary_search(&j).map_or(0.0, |k| vals[k])
    }

    /// Column indices (increasing) and values of row `i`.
    ///
    /// # Panics
    ///
    /// When `i` is not a row.
    #[must_use]
    pub fn row(&self, i: usize) -> (&[usize], &[f64]) {
        let range = self.row_start[i]..self.row_start[i + 1];
        (&self.columns[range.clone()], &self.values[range])
    }

    /// `A x`.
    ///
    /// # Errors
    ///
    /// Refuses `x` of the wrong length.
    pub fn mul_vec(&self, x: &[f64]) -> NumericResult<Vec<f64>> {
        if x.len() != self.cols {
            return Err(NumericError::DimensionMismatch {
                name: "vector",
                expected: self.cols,
                found: x.len(),
            });
        }
        let mut y = vec![0.0; self.rows];
        self.mul_into(x, &mut y);
        Ok(y)
    }

    fn mul_into(&self, x: &[f64], y: &mut [f64]) {
        for (i, out) in y.iter_mut().enumerate() {
            let (cols, vals) = self.row(i);
            *out = cols.iter().zip(vals).map(|(&j, &v)| v * x[j]).sum();
        }
    }

    /// Refuse a matrix that is not square or not symmetric to
    /// `16 eps max(|a_ij|, |a_ji|)`.
    fn require_symmetric(&self) -> NumericResult<()> {
        if self.rows != self.cols {
            return Err(NumericError::DimensionMismatch {
                name: "symmetric matrix",
                expected: self.rows,
                found: self.cols,
            });
        }
        for i in 0..self.rows {
            let (cols, vals) = self.row(i);
            for (&j, &v) in cols.iter().zip(vals) {
                if j <= i {
                    continue;
                }
                let w = self.get(j, i);
                if (v - w).abs() > 16.0 * f64::EPSILON * v.abs().max(w.abs()) {
                    return Err(NumericError::NotSymmetric { row: i, column: j });
                }
            }
        }
        Ok(())
    }
}

/// Stopping rule and budget for [`conjugate_gradient`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConjugateGradientOptions {
    /// Stop when `||b - A x||_2 <= relative_tolerance * ||b||_2`.
    pub relative_tolerance: f64,
    /// Iteration budget; `None` allows `10 n + 100` for an `n x n` system.
    pub max_iterations: Option<usize>,
}

impl Default for ConjugateGradientOptions {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-12,
            max_iterations: None,
        }
    }
}

/// The result of [`conjugate_gradient`].
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct ConjugateGradientSolution {
    /// The final iterate.
    pub x: Vec<f64>,
    /// `||b - A x||_2 / ||b||_2`, recomputed from the matrix at the final
    /// iterate (zero when `b` is zero). A computed value, not an estimate.
    pub relative_residual: f64,
    /// Iterations taken.
    pub iterations: usize,
    /// [`Status::Converged`] when the relative residual met the tolerance,
    /// [`Status::BudgetExhausted`] when the iteration budget ran out first.
    pub status: Status,
}

impl ConjugateGradientSolution {
    /// The solution when its tolerance was met.
    ///
    /// # Errors
    ///
    /// [`NumericError::NotConverged`] with the relative residual reached.
    pub fn converged(self) -> NumericResult<Self> {
        if self.status == Status::Converged {
            Ok(self)
        } else {
            Err(NumericError::NotConverged {
                name: "conjugate gradient",
                error_estimate: self.relative_residual,
            })
        }
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm2(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// Solve `A x = b` for a symmetric positive definite sparse `A` by the
/// Jacobi-preconditioned conjugate-gradient method, starting from
/// `initial` (zero when `None`).
///
/// A run that exhausts its budget still returns its last iterate and
/// residual with [`Status::BudgetExhausted`]; call
/// [`ConjugateGradientSolution::converged`] to refuse it. A zero `b` returns
/// the exact zero solution without iterating.
///
/// # Errors
///
/// - [`NumericError::DimensionMismatch`] for a non-square matrix or a
///   wrongly sized `b` or `initial`.
/// - [`NumericError::NonFiniteInput`] for a non-finite `b` or `initial`,
///   [`NumericError::InvalidArgument`] for a negative or non-finite
///   tolerance.
/// - [`NumericError::NotSymmetric`] for an asymmetric matrix.
/// - [`NumericError::NotPositiveDefinite`] when a diagonal entry is not
///   positive (`pivot` is its row) or when the iteration meets a direction
///   `p` with `p^T A p <= 0`, which proves `A` is not positive definite
///   (`pivot` is the iteration).
pub fn conjugate_gradient(
    a: &SparseMatrix,
    b: &[f64],
    initial: Option<&[f64]>,
    options: ConjugateGradientOptions,
) -> NumericResult<ConjugateGradientSolution> {
    a.require_symmetric()?;
    let n = a.rows;
    if b.len() != n {
        return Err(NumericError::DimensionMismatch {
            name: "right-hand side",
            expected: n,
            found: b.len(),
        });
    }
    all_finite(b, "right-hand side")?;
    let relative_tolerance = tolerance(options.relative_tolerance, "relative tolerance")?;
    let mut x = match initial {
        Some(x0) => {
            if x0.len() != n {
                return Err(NumericError::DimensionMismatch {
                    name: "initial guess",
                    expected: n,
                    found: x0.len(),
                });
            }
            all_finite(x0, "initial guess")?;
            x0.to_vec()
        }
        None => vec![0.0; n],
    };
    let mut inverse_diagonal = Vec::with_capacity(n);
    for i in 0..n {
        let d = a.get(i, i);
        if d <= 0.0 {
            return Err(NumericError::NotPositiveDefinite { pivot: i });
        }
        inverse_diagonal.push(1.0 / d);
    }
    let b_norm = norm2(b);
    if b_norm == 0.0 {
        return Ok(ConjugateGradientSolution {
            x: vec![0.0; n],
            relative_residual: 0.0,
            iterations: 0,
            status: Status::Converged,
        });
    }
    let target = relative_tolerance * b_norm;
    let budget = options.max_iterations.unwrap_or(10 * n + 100);

    let mut ax = vec![0.0; n];
    let true_residual = |x: &[f64], ax: &mut Vec<f64>, r: &mut Vec<f64>| {
        a.mul_into(x, ax);
        for i in 0..n {
            r[i] = b[i] - ax[i];
        }
        norm2(r)
    };
    let mut r = vec![0.0; n];
    let mut r_norm = true_residual(&x, &mut ax, &mut r);
    let mut z: Vec<f64> = r
        .iter()
        .zip(&inverse_diagonal)
        .map(|(r, d)| r * d)
        .collect();
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let mut ap = vec![0.0; n];
    let mut iterations = 0;
    loop {
        if r_norm <= target {
            // The recurred residual says done; confirm on the true one and
            // restart from it if rounding let the two drift apart.
            r_norm = true_residual(&x, &mut ax, &mut r);
            if r_norm <= target {
                return Ok(ConjugateGradientSolution {
                    x,
                    relative_residual: r_norm / b_norm,
                    iterations,
                    status: Status::Converged,
                });
            }
            for i in 0..n {
                z[i] = r[i] * inverse_diagonal[i];
            }
            p.clone_from(&z);
            rz = dot(&r, &z);
        }
        if iterations >= budget {
            r_norm = true_residual(&x, &mut ax, &mut r);
            let status = if r_norm <= target {
                Status::Converged
            } else {
                Status::BudgetExhausted
            };
            return Ok(ConjugateGradientSolution {
                x,
                relative_residual: r_norm / b_norm,
                iterations,
                status,
            });
        }
        a.mul_into(&p, &mut ap);
        let curvature = dot(&p, &ap);
        if curvature <= 0.0 {
            return Err(NumericError::NotPositiveDefinite { pivot: iterations });
        }
        let alpha = rz / curvature;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        r_norm = norm2(&r);
        for i in 0..n {
            z[i] = r[i] * inverse_diagonal[i];
        }
        let rz_next = dot(&r, &z);
        let beta = rz_next / rz;
        rz = rz_next;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        iterations += 1;
    }
}
