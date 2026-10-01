//! Dense solves and least squares: exact systems and fits, ill-conditioned
//! Hilbert matrices with the error estimate checked against the true
//! error, rank and condition reporting, refusals.

use axiolid_numeric::{
    constrained_least_squares, least_squares, Cholesky, Lu, Matrix, NumericError, Qr,
};

fn hilbert(n: usize) -> Matrix {
    Matrix::from_fn(n, n, |i, j| 1.0 / (i + j + 1) as f64)
}

fn relative_error(x: &[f64], exact: &[f64]) -> f64 {
    let diff: f64 = x.iter().zip(exact).map(|(a, b)| (a - b).abs()).sum();
    diff / exact.iter().map(|v| v.abs()).sum::<f64>()
}

#[test]
fn lu_solves_a_small_system_exactly() {
    let a = Matrix::from_rows(&[&[2.0, 1.0, -1.0], &[-3.0, -1.0, 2.0], &[-2.0, 1.0, 2.0]])
        .expect("rectangular");
    let s = Lu::new(&a)
        .expect("regular")
        .solve(&[8.0, -11.0, -3.0])
        .expect("sized");
    assert!(relative_error(&s.x, &[2.0, 3.0, -1.0]) <= 1e-15);
    assert!(s.residual_norm <= 1e-14);
    assert!(s.relative_error_estimate <= 1e-13);
}

#[test]
fn lu_pivots() {
    // Without row exchange the 1e-20 pivot destroys the answer.
    let a = Matrix::from_rows(&[&[1e-20, 1.0], &[1.0, 1.0]]).expect("rectangular");
    let s = Lu::new(&a)
        .expect("regular")
        .solve(&[1.0, 2.0])
        .expect("sized");
    assert!(
        (s.x[0] - 1.0).abs() <= 1e-15 && (s.x[1] - 1.0).abs() <= 1e-15,
        "{:?}",
        s.x
    );
}

#[test]
fn lu_estimates_condition_and_error_on_hilbert_matrices() {
    // True 1-norm condition numbers of the Hilbert matrices.
    for (n, kappa) in [(4usize, 2.837_5e4), (6, 2.907_028_7e7), (8, 3.387_279_1e10)] {
        let a = hilbert(n);
        let lu = Lu::new(&a).expect("regular");
        let estimate = lu.condition_estimate();
        assert!(
            estimate <= kappa * 1.0001 && estimate >= kappa / 3.0,
            "n = {n}: estimate {estimate} vs {kappa}"
        );
        let exact = vec![1.0; n];
        let b = a.mul_vec(&exact).expect("sized");
        let s = lu.solve(&b).expect("sized");
        let error = relative_error(&s.x, &exact);
        assert!(
            error <= s.relative_error_estimate,
            "n = {n}: error {error} exceeds estimate {}",
            s.relative_error_estimate
        );
        // The estimate is informative, not just huge.
        assert!(s.relative_error_estimate <= 1e4 * kappa * f64::EPSILON);
    }
}

#[test]
fn lu_refuses_singular_and_malformed_matrices() {
    let singular = Matrix::from_rows(&[&[1.0, 2.0], &[2.0, 4.0]]).expect("rectangular");
    assert_eq!(
        Lu::new(&singular),
        Err(NumericError::Singular {
            name: "system matrix"
        })
    );
    // Numerically singular: condition far above 1/eps.
    assert!(matches!(
        Lu::new(&hilbert(14)),
        Err(NumericError::Singular { .. })
    ));
    assert!(matches!(
        Lu::new(&Matrix::zeros(2, 3)),
        Err(NumericError::DimensionMismatch { .. })
    ));
    let mut bad = Matrix::identity(2);
    bad[(1, 0)] = f64::INFINITY;
    assert_eq!(
        Lu::new(&bad),
        Err(NumericError::NonFiniteInput {
            name: "system matrix"
        })
    );
    let lu = Lu::new(&Matrix::identity(2)).expect("regular");
    assert!(matches!(
        lu.solve(&[1.0]),
        Err(NumericError::DimensionMismatch { .. })
    ));
    assert_eq!(
        lu.solve(&[1.0, f64::NAN]),
        Err(NumericError::NonFiniteInput {
            name: "right-hand side"
        })
    );
}

#[test]
fn cholesky_solves_spd_systems_with_error_estimates() {
    let a = Matrix::from_rows(&[
        &[4.0, 12.0, -16.0],
        &[12.0, 37.0, -43.0],
        &[-16.0, -43.0, 98.0],
    ])
    .expect("rectangular");
    let c = Cholesky::new(&a).expect("SPD");
    let l = c.factor();
    let expected = [[2.0, 0.0, 0.0], [6.0, 1.0, 0.0], [-8.0, 5.0, 3.0]];
    for (i, row) in expected.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            assert_eq!(l[(i, j)], *v);
        }
    }
    let s = c.solve(&[0.0, 6.0, 39.0]).expect("sized");
    assert!(relative_error(&s.x, &[1.0, 1.0, 1.0]) <= 1e-14, "{:?}", s.x);

    for n in [4usize, 7, 9] {
        let a = hilbert(n);
        let exact: Vec<f64> = (0..n).map(|i| (i as f64 + 1.0).recip()).collect();
        let b = a.mul_vec(&exact).expect("sized");
        let s = Cholesky::new(&a).expect("SPD").solve(&b).expect("sized");
        let error = relative_error(&s.x, &exact);
        assert!(
            error <= s.relative_error_estimate,
            "n = {n}: {error} > {}",
            s.relative_error_estimate
        );
    }
}

#[test]
fn cholesky_refuses_asymmetric_and_indefinite_matrices() {
    let asymmetric = Matrix::from_rows(&[&[2.0, 1.0], &[0.0, 2.0]]).expect("rectangular");
    assert_eq!(
        Cholesky::new(&asymmetric),
        Err(NumericError::NotSymmetric { row: 1, column: 0 })
    );
    let indefinite = Matrix::from_rows(&[&[1.0, 2.0], &[2.0, 1.0]]).expect("rectangular");
    assert_eq!(
        Cholesky::new(&indefinite),
        Err(NumericError::NotPositiveDefinite { pivot: 1 })
    );
    let negative = Matrix::from_rows(&[&[-1.0, 0.0], &[0.0, 1.0]]).expect("rectangular");
    assert_eq!(
        Cholesky::new(&negative),
        Err(NumericError::NotPositiveDefinite { pivot: 0 })
    );
}

#[test]
fn least_squares_recovers_an_exact_fit() {
    // y = 1 + 2t - 0.5t^2 + 0.25t^3 sampled at 12 points.
    let coefficients = [1.0, 2.0, -0.5, 0.25];
    let ts: Vec<f64> = (0..12).map(|i| -1.0 + i as f64 * 0.25).collect();
    let a = Matrix::from_fn(ts.len(), 4, |i, j| ts[i].powi(j as i32));
    let b: Vec<f64> = ts
        .iter()
        .map(|t| coefficients.iter().rev().fold(0.0, |acc, c| acc * t + c))
        .collect();
    let s = least_squares(&a, &b).expect("full rank");
    let error = relative_error(&s.x, &coefficients);
    assert!(error <= 1e-14, "{:?}", s.x);
    assert!(error <= s.relative_error_estimate);
    assert!(s.residual_norm <= 1e-13);
}

#[test]
fn least_squares_matches_the_closed_form_line_fit() {
    // Points off a line: the normal-equation closed form for slope and
    // intercept is exact in rationals here.
    let points = [(0.0, 1.0), (1.0, 3.0), (2.0, 2.0), (3.0, 5.0), (4.0, 4.0)];
    let a = Matrix::from_fn(
        points.len(),
        2,
        |i, j| if j == 0 { 1.0 } else { points[i].0 },
    );
    let b: Vec<f64> = points.iter().map(|p| p.1).collect();
    let s = least_squares(&a, &b).expect("full rank");
    // n = 5, sum t = 10, sum t^2 = 30, sum y = 15, sum ty = 38.
    let slope = (5.0 * 38.0 - 10.0 * 15.0) / (5.0 * 30.0 - 100.0);
    let intercept = (15.0 - slope * 10.0) / 5.0;
    assert!(
        (s.x[0] - intercept).abs() <= 1e-13 && (s.x[1] - slope).abs() <= 1e-13,
        "{:?} vs {intercept}, {slope}",
        s.x
    );
    let residual: f64 = points
        .iter()
        .map(|(t, y)| (intercept + slope * t - y).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!((s.residual_norm - residual).abs() <= 1e-14);
}

#[test]
fn least_squares_error_estimate_bounds_ill_conditioned_fits() {
    // Monomial basis on [0, 1]: badly conditioned, increasingly with degree.
    for degree in [5usize, 8, 10] {
        let ts: Vec<f64> = (0..30).map(|i| i as f64 / 29.0).collect();
        let a = Matrix::from_fn(ts.len(), degree + 1, |i, j| ts[i].powi(j as i32));
        let exact: Vec<f64> = (0..=degree).map(|j| 1.0 / (j as f64 + 1.0)).collect();
        let b = a.mul_vec(&exact).expect("sized");
        let qr = Qr::new(&a).expect("finite");
        assert_eq!(qr.rank(), degree + 1);
        let s = qr.solve_least_squares(&b).expect("full rank");
        let error = relative_error(&s.x, &exact);
        assert!(
            error <= s.relative_error_estimate,
            "degree {degree}: {error} > {} (condition {})",
            s.relative_error_estimate,
            qr.condition_estimate()
        );
    }
}

#[test]
fn qr_detects_rank_and_refuses_rank_deficient_fits() {
    // Third column = first + second.
    let a = Matrix::from_fn(6, 3, |i, j| {
        let t = i as f64;
        match j {
            0 => 1.0,
            1 => t,
            _ => 1.0 + t,
        }
    });
    let qr = Qr::new(&a).expect("finite");
    assert_eq!(qr.rank(), 2);
    assert_eq!(
        qr.solve_least_squares(&[0.0; 6]),
        Err(NumericError::RankDeficient {
            name: "least-squares matrix",
            rank: 2,
            required: 3
        })
    );
    // Underdetermined.
    assert!(matches!(
        least_squares(
            &Matrix::from_fn(2, 3, |i, j| (i + j) as f64 + (i * j) as f64),
            &[1.0, 2.0]
        ),
        Err(NumericError::RankDeficient { .. })
    ));
    let mut bad = Matrix::identity(3);
    bad[(0, 2)] = f64::NAN;
    assert_eq!(
        Qr::new(&bad),
        Err(NumericError::NonFiniteInput {
            name: "least-squares matrix"
        })
    );
}

#[test]
fn constrained_fit_through_a_fixed_point() {
    // Fit y = c0 + c1 t through (0, 0): the constraint forces c0 = 0 and the
    // slope is then sum(t y) / sum(t^2).
    let points = [(1.0, 2.1), (2.0, 3.9), (3.0, 6.2), (4.0, 7.8)];
    let a = Matrix::from_fn(
        points.len(),
        2,
        |i, j| if j == 0 { 1.0 } else { points[i].0 },
    );
    let b: Vec<f64> = points.iter().map(|p| p.1).collect();
    let c = Matrix::from_rows(&[&[1.0, 0.0]]).expect("rectangular");
    let s = constrained_least_squares(&a, &b, &c, &[0.0]).expect("well posed");
    let sty: f64 = points.iter().map(|(t, y)| t * y).sum();
    let stt: f64 = points.iter().map(|(t, _)| t * t).sum();
    assert!(s.x[0].abs() <= 1e-15, "{:?}", s.x);
    assert!((s.x[1] - sty / stt).abs() <= 1e-14);
    assert!(s.constraint_residual_norm <= 1e-15);
    assert!(s.relative_error_estimate < 1e-12);
}

#[test]
fn constrained_fit_with_interpolated_ends() {
    // A quadratic through both end samples, best fit in between: the ends
    // are met to rounding while the unconstrained fit misses them.
    let ts: Vec<f64> = (0..9).map(|i| i as f64 / 8.0).collect();
    let ys: Vec<f64> = ts.iter().map(|t| (3.0 * t).sin()).collect();
    let a = Matrix::from_fn(ts.len(), 3, |i, j| ts[i].powi(j as i32));
    let c = Matrix::from_rows(&[&[1.0, 0.0, 0.0], &[1.0, 1.0, 1.0]]).expect("rectangular");
    let d = [ys[0], ys[8]];
    let s = constrained_least_squares(&a, &ys, &c, &d).expect("well posed");
    let at = |t: f64| s.x[0] + s.x[1] * t + s.x[2] * t * t;
    assert!((at(0.0) - ys[0]).abs() <= 1e-14 && (at(1.0) - ys[8]).abs() <= 1e-14);
    let free = least_squares(&a, &ys).expect("full rank");
    assert!(free.residual_norm < s.residual_norm);
    let free_end = free.x[0] + free.x[1] + free.x[2];
    assert!((free_end - ys[8]).abs() > 1e-3);
    // Optimality: any other feasible quadratic does worse. Perturb along the
    // constraints' null space, (0, t, -t).
    for step in [1e-3, -1e-3] {
        let other = [s.x[0], s.x[1] + step, s.x[2] - step];
        let r: f64 = ts
            .iter()
            .zip(&ys)
            .map(|(t, y)| (other[0] + other[1] * t + other[2] * t * t - y).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(r > s.residual_norm);
    }
}

#[test]
fn constrained_fit_fully_determined_by_constraints() {
    let a = Matrix::from_rows(&[&[1.0, 1.0]]).expect("rectangular");
    let c = Matrix::from_rows(&[&[1.0, 0.0], &[1.0, -1.0]]).expect("rectangular");
    let s = constrained_least_squares(&a, &[10.0], &c, &[2.0, 1.0]).expect("determined");
    assert!(
        (s.x[0] - 2.0).abs() <= 4e-15 && (s.x[1] - 1.0).abs() <= 4e-15,
        "{:?}",
        s.x
    );
    assert!((s.residual_norm - 7.0).abs() <= 1e-14);
}

#[test]
fn constrained_fit_refuses_redundant_constraints_and_free_directions() {
    let a = Matrix::identity(3);
    let redundant = Matrix::from_rows(&[&[1.0, 1.0, 0.0], &[2.0, 2.0, 0.0]]).expect("rectangular");
    assert_eq!(
        constrained_least_squares(&a, &[1.0, 1.0, 1.0], &redundant, &[1.0, 2.0]),
        Err(NumericError::RankDeficient {
            name: "constraint matrix",
            rank: 1,
            required: 2
        })
    );
    // The objective ignores x2, and the constraint does not fix it.
    let a = Matrix::from_rows(&[&[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0]]).expect("rectangular");
    let c = Matrix::from_rows(&[&[1.0, 1.0, 0.0]]).expect("rectangular");
    assert!(matches!(
        constrained_least_squares(&a, &[1.0, 1.0], &c, &[1.0]),
        Err(NumericError::RankDeficient {
            rank: 1,
            required: 2,
            ..
        })
    ));
    assert!(matches!(
        constrained_least_squares(&a, &[1.0, 1.0], &Matrix::identity(2), &[1.0, 1.0]),
        Err(NumericError::DimensionMismatch { .. })
    ));
}

struct XorShift(u64);
impl XorShift {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    }
}

/// True 1-norm condition number from the explicit inverse.
fn condition1(a: &Matrix, lu: &Lu) -> f64 {
    let n = a.rows();
    let mut inverse_norm = 0.0f64;
    for j in 0..n {
        let mut e = vec![0.0; n];
        e[j] = 1.0;
        let column = lu.solve(&e).expect("sized").x;
        inverse_norm = inverse_norm.max(column.iter().map(|v| v.abs()).sum());
    }
    a.norm1() * inverse_norm
}

#[test]
fn lu_condition_estimate_tracks_the_true_condition_on_pivoting_matrices() {
    let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
    for case in 0..200 {
        let n = 3 + case % 6;
        // Random nonsymmetric matrices, graded so row exchanges happen.
        let a = Matrix::from_fn(n, n, |i, j| {
            rng.next() * 10f64.powi(((i * 7 + j * 3) % 5) as i32 - 2)
        });
        let lu = Lu::new(&a).expect("regular");
        let truth = condition1(&a, &lu);
        let estimate = lu.condition_estimate();
        assert!(
            estimate <= truth * (1.0 + 1e-9) && estimate >= truth / 2.0,
            "case {case}: estimate {estimate} vs {truth}"
        );
    }
}

/// Reference solution by iterative refinement with the residual in
/// double-double arithmetic.
fn refined(a: &Matrix, b: &[f64], lu: &Lu) -> Vec<f64> {
    let n = b.len();
    let mut x = lu.solve(b).expect("sized").x;
    for _ in 0..6 {
        let r: Vec<f64> = (0..n)
            .map(|i| {
                let (mut hi, mut lo) = (b[i], 0.0f64);
                for j in 0..n {
                    let p = -a[(i, j)] * x[j];
                    let pe = (-a[(i, j)]).mul_add(x[j], -p);
                    let s = hi + p;
                    let bb = s - hi;
                    lo += (hi - (s - bb)) + (p - bb) + pe;
                    hi = s;
                }
                hi + lo
            })
            .collect();
        let d = lu.solve(&r).expect("sized").x;
        for (xi, di) in x.iter_mut().zip(&d) {
            *xi += di;
        }
    }
    x
}

#[test]
fn lu_error_estimate_bounds_the_error_on_random_ill_conditioned_systems() {
    let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
    for case in 0..300 {
        let n = 4 + case % 5;
        // Nearly dependent rows: row i = base + 10^-k * noise.
        let scale = 10f64.powi(-((case % 10) as i32));
        let base: Vec<f64> = (0..n).map(|_| rng.next()).collect();
        let a = Matrix::from_fn(n, n, |i, j| {
            if i == 0 {
                base[j]
            } else {
                base[j] + scale * rng.next() + if i == j { scale } else { 0.0 }
            }
        });
        let b: Vec<f64> = (0..n).map(|_| rng.next()).collect();
        let Ok(lu) = Lu::new(&a) else { continue };
        let s = lu.solve(&b).expect("sized");
        let reference = refined(&a, &b, &lu);
        let error = relative_error(&s.x, &reference);
        assert!(
            error <= s.relative_error_estimate,
            "case {case}: error {error} exceeds estimate {} (condition {})",
            s.relative_error_estimate,
            lu.condition_estimate()
        );
    }
}
