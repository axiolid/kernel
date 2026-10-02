//! Sparse storage and the conjugate-gradient solve: assembly, products,
//! solutions checked against known answers and a dense factorisation,
//! the residual it reports checked against the matrix, refusals.

use axiolid_numeric::{
    conjugate_gradient, Cholesky, ConjugateGradientOptions, Matrix, NumericError, SparseMatrix,
    Status,
};

/// The `n x n` second-difference matrix `tridiag(-1, 2, -1)`.
fn poisson_1d(n: usize) -> SparseMatrix {
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i, 2.0));
        if i + 1 < n {
            t.push((i, i + 1, -1.0));
            t.push((i + 1, i, -1.0));
        }
    }
    SparseMatrix::from_triplets(n, n, &t).expect("in shape")
}

/// The 5-point Laplacian on an `m x m` grid with Dirichlet boundary.
fn poisson_2d(m: usize) -> SparseMatrix {
    let n = m * m;
    let mut t = Vec::new();
    for i in 0..m {
        for j in 0..m {
            let k = i * m + j;
            t.push((k, k, 4.0));
            if i + 1 < m {
                t.push((k, k + m, -1.0));
                t.push((k + m, k, -1.0));
            }
            if j + 1 < m {
                t.push((k, k + 1, -1.0));
                t.push((k + 1, k, -1.0));
            }
        }
    }
    SparseMatrix::from_triplets(n, n, &t).expect("in shape")
}

fn residual(a: &SparseMatrix, x: &[f64], b: &[f64]) -> f64 {
    let ax = a.mul_vec(x).expect("sized");
    let r: f64 = ax.iter().zip(b).map(|(p, q)| (p - q) * (p - q)).sum();
    r.sqrt() / b.iter().map(|v| v * v).sum::<f64>().sqrt()
}

#[test]
fn triplets_are_sorted_and_duplicates_summed() {
    let a = SparseMatrix::from_triplets(
        3,
        4,
        &[
            (2, 3, 1.0),
            (0, 1, 2.0),
            (2, 0, -1.0),
            (0, 1, 0.5),
            (1, 2, 4.0),
        ],
    )
    .expect("in shape");
    assert_eq!((a.rows(), a.cols(), a.nnz()), (3, 4, 4));
    assert_eq!(a.get(0, 1), 2.5);
    assert_eq!(a.get(1, 2), 4.0);
    assert_eq!(a.get(2, 0), -1.0);
    assert_eq!(a.get(2, 3), 1.0);
    assert_eq!(a.get(1, 1), 0.0);
    assert_eq!(a.row(2), (&[0usize, 3][..], &[-1.0, 1.0][..]));
    assert_eq!(a.row(1), (&[2usize][..], &[4.0][..]));
    let b = SparseMatrix::from_triplets(
        3,
        4,
        &[
            (1, 2, 4.0),
            (0, 1, 2.0),
            (0, 1, 0.5),
            (2, 3, 1.0),
            (2, 0, -1.0),
        ],
    )
    .expect("in shape");
    assert_eq!(a, b, "storage does not depend on triplet order");
    assert_eq!(
        a.mul_vec(&[1.0, 2.0, 3.0, 4.0]).expect("sized"),
        vec![5.0, 12.0, 3.0]
    );
}

#[test]
fn assembly_refuses_bad_triplets() {
    assert!(matches!(
        SparseMatrix::from_triplets(2, 2, &[(2, 0, 1.0)]),
        Err(NumericError::InvalidArgument { .. })
    ));
    assert!(matches!(
        SparseMatrix::from_triplets(2, 2, &[(0, 2, 1.0)]),
        Err(NumericError::InvalidArgument { .. })
    ));
    assert!(matches!(
        SparseMatrix::from_triplets(2, 2, &[(0, 0, f64::NAN)]),
        Err(NumericError::NonFiniteInput { .. })
    ));
    let a = poisson_1d(3);
    assert!(matches!(
        a.mul_vec(&[1.0]),
        Err(NumericError::DimensionMismatch { .. })
    ));
}

#[test]
fn cg_solves_the_1d_poisson_problem() {
    // Exact solution of tridiag(-1,2,-1) x = e_1 + e_n is the all-ones vector.
    let n = 300;
    let a = poisson_1d(n);
    let mut b = vec![0.0; n];
    b[0] = 1.0;
    b[n - 1] = 1.0;
    let s = conjugate_gradient(&a, &b, None, ConjugateGradientOptions::default())
        .expect("SPD")
        .converged()
        .expect("met the tolerance");
    assert_eq!(s.status, Status::Converged);
    assert!(s.relative_residual <= 1e-12, "{}", s.relative_residual);
    assert!((residual(&a, &s.x, &b) - s.relative_residual).abs() <= 1e-15);
    let max_error = s.x.iter().map(|v| (v - 1.0).abs()).fold(0.0, f64::max);
    // Condition number ~ 4 n^2 / pi^2 ~ 3.6e4.
    assert!(max_error <= 1e-7, "max error {max_error}");
    // Exact arithmetic finishes in n iterations; allow some for rounding.
    assert!(s.iterations <= n + 10, "{} iterations", s.iterations);
}

#[test]
fn cg_matches_a_dense_cholesky_solve() {
    let m = 12;
    let a = poisson_2d(m);
    let n = m * m;
    let b: Vec<f64> = (0..n).map(|i| ((i * 37) % 11) as f64 - 5.0).collect();
    let dense = Matrix::from_fn(n, n, |i, j| a.get(i, j));
    let reference = Cholesky::new(&dense)
        .expect("SPD")
        .solve(&b)
        .expect("sized");
    let s = conjugate_gradient(&a, &b, None, ConjugateGradientOptions::default())
        .expect("SPD")
        .converged()
        .expect("met the tolerance");
    let scale = reference.x.iter().map(|v| v.abs()).fold(0.0, f64::max);
    for (p, q) in s.x.iter().zip(&reference.x) {
        assert!((p - q).abs() <= 1e-10 * scale, "{p} vs {q}");
    }
}

#[test]
fn cg_solves_a_squared_laplacian() {
    // The Gram form L^T L is how a bi-Laplacian fairing system arises; its
    // condition number is the square of the Laplacian's.
    let m = 10;
    let l = poisson_2d(m);
    let n = m * m;
    let mut t = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let v: f64 = (0..n).map(|k| l.get(k, i) * l.get(k, j)).sum();
            if v != 0.0 {
                t.push((i, j, v));
            }
        }
    }
    let a = SparseMatrix::from_triplets(n, n, &t).expect("in shape");
    let exact: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin()).collect();
    let b = a.mul_vec(&exact).expect("sized");
    let s = conjugate_gradient(&a, &b, None, ConjugateGradientOptions::default())
        .expect("SPD")
        .converged()
        .expect("met the tolerance");
    let max_error =
        s.x.iter()
            .zip(&exact)
            .map(|(p, q)| (p - q).abs())
            .fold(0.0, f64::max);
    assert!(max_error <= 1e-8, "max error {max_error}");
}

#[test]
fn zero_right_hand_side_returns_exact_zero() {
    let a = poisson_1d(5);
    let s = conjugate_gradient(
        &a,
        &[0.0; 5],
        Some(&[1.0, 2.0, 3.0, 4.0, 5.0]),
        ConjugateGradientOptions::default(),
    )
    .expect("SPD");
    assert_eq!(s.x, vec![0.0; 5]);
    assert_eq!((s.iterations, s.status), (0, Status::Converged));
}

#[test]
fn an_exact_initial_guess_needs_no_iteration() {
    let a = poisson_1d(4);
    let b = [1.0, 0.0, 0.0, 1.0];
    let s = conjugate_gradient(&a, &b, Some(&[1.0; 4]), ConjugateGradientOptions::default())
        .expect("SPD");
    assert_eq!(s.iterations, 0);
    assert_eq!(s.x, vec![1.0; 4]);
}

#[test]
fn an_exhausted_budget_is_reported_and_refusable() {
    let a = poisson_1d(100);
    let b = vec![1.0; 100];
    let s = conjugate_gradient(
        &a,
        &b,
        None,
        ConjugateGradientOptions {
            max_iterations: Some(3),
            ..ConjugateGradientOptions::default()
        },
    )
    .expect("SPD");
    assert_eq!((s.iterations, s.status), (3, Status::BudgetExhausted));
    assert!((residual(&a, &s.x, &b) - s.relative_residual).abs() <= 1e-15);
    assert!(s.relative_residual > 1e-3);
    assert!(matches!(
        s.converged(),
        Err(NumericError::NotConverged { .. })
    ));
}

#[test]
fn the_solve_is_bitwise_deterministic() {
    let a = poisson_2d(9);
    let b: Vec<f64> = (0..81).map(|i| (i as f64).cos()).collect();
    let first = conjugate_gradient(&a, &b, None, ConjugateGradientOptions::default()).expect("SPD");
    let second =
        conjugate_gradient(&a, &b, None, ConjugateGradientOptions::default()).expect("SPD");
    assert_eq!(first, second);
}

#[test]
fn cg_refuses_matrices_that_are_not_spd() {
    let asymmetric =
        SparseMatrix::from_triplets(2, 2, &[(0, 0, 2.0), (0, 1, 1.0), (1, 1, 2.0)]).expect("ok");
    assert_eq!(
        conjugate_gradient(&asymmetric, &[1.0, 1.0], None, Default::default()),
        Err(NumericError::NotSymmetric { row: 0, column: 1 })
    );
    let negative_diagonal =
        SparseMatrix::from_triplets(2, 2, &[(0, 0, 2.0), (1, 1, -1.0)]).expect("ok");
    assert_eq!(
        conjugate_gradient(&negative_diagonal, &[1.0, 1.0], None, Default::default()),
        Err(NumericError::NotPositiveDefinite { pivot: 1 })
    );
    // Positive diagonal, eigenvalues 3 and -1: the first direction already
    // has negative curvature.
    let indefinite =
        SparseMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (0, 1, 2.0), (1, 0, 2.0), (1, 1, 1.0)])
            .expect("ok");
    assert_eq!(
        conjugate_gradient(&indefinite, &[1.0, -1.0], None, Default::default()),
        Err(NumericError::NotPositiveDefinite { pivot: 0 })
    );
    let rectangular = SparseMatrix::from_triplets(2, 3, &[(0, 0, 1.0)]).expect("ok");
    assert!(matches!(
        conjugate_gradient(&rectangular, &[1.0, 1.0], None, Default::default()),
        Err(NumericError::DimensionMismatch { .. })
    ));
}

#[test]
fn cg_refuses_bad_vectors_and_options() {
    let a = poisson_1d(3);
    assert!(matches!(
        conjugate_gradient(&a, &[1.0, 1.0], None, Default::default()),
        Err(NumericError::DimensionMismatch { .. })
    ));
    assert!(matches!(
        conjugate_gradient(&a, &[1.0, f64::INFINITY, 1.0], None, Default::default()),
        Err(NumericError::NonFiniteInput { .. })
    ));
    assert!(matches!(
        conjugate_gradient(&a, &[1.0; 3], Some(&[0.0; 2]), Default::default()),
        Err(NumericError::DimensionMismatch { .. })
    ));
    assert!(matches!(
        conjugate_gradient(
            &a,
            &[1.0; 3],
            Some(&[0.0, f64::NAN, 0.0]),
            Default::default()
        ),
        Err(NumericError::NonFiniteInput { .. })
    ));
    assert!(matches!(
        conjugate_gradient(
            &a,
            &[1.0; 3],
            None,
            ConjugateGradientOptions {
                relative_tolerance: -1.0,
                max_iterations: None
            }
        ),
        Err(NumericError::InvalidArgument { .. })
    ));
}
