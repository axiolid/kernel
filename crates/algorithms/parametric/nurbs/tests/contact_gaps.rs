//! Roots next to a flat stretch of a second field (#427).
//!
//! Two plane curves, `v = 0.3 u^2` and `v = 0.3 u^2 - p(u)`, meet where
//! `p(u)` vanishes. Each `p` here is a high-order contact displaced by a
//! little, so that the field is flat and the enclosures straddle zero over a
//! stretch of the curve, while the true roots sit outside that stretch: a
//! chain that merged the stretch into one multiple root without searching
//! the gaps beside it reported one root of order 8 between them, or
//! swallowed a crossing. The flat stretch is below what the enclosures
//! resolve, so refusing it as `Undecided` is allowed; a decided answer is
//! exactly the crossings.

use axiolid_core::{Interval, Tolerance};
use axiolid_curve::implicit::{Basis, SeriesField2};
use axiolid_curve::{Axis, Curve2, Field2, ImplicitCell, ImplicitCurve2};
use axiolid_nurbs::{section_curve_curve_intersection2, ExactCurveRefusal};

/// The coefficients (constant first) of `prod (u - r)^m`.
fn product(factors: &[(f64, usize)]) -> Vec<f64> {
    let mut p = vec![1.0];
    for &(r, m) in factors {
        for _ in 0..m {
            let mut q = vec![0.0; p.len() + 1];
            for (i, c) in p.iter().enumerate() {
                q[i + 1] += c;
                q[i] -= c * r;
            }
            p = q;
        }
    }
    p
}

fn field(p: &[f64]) -> Field2 {
    let n = p.len().max(3);
    let mut c = vec![vec![0.0; 2]; n];
    c[0][1] = 1.0;
    c[2][0] = -0.3;
    for (i, &x) in p.iter().enumerate() {
        c[i][0] += x;
    }
    Field2::Series(SeriesField2 {
        u: Basis::Power,
        v: Basis::Power,
        coefficients: c,
    })
}

fn curve(p: &[f64], low: f64, high: f64) -> Curve2 {
    Curve2::Implicit(ImplicitCurve2 {
        field: field(p),
        cells: vec![ImplicitCell {
            axis: Axis::U,
            from: -1.0,
            to: 1.0,
            low,
            high,
            bridge: None,
        }],
    })
}

/// The `(u, multiplicity)` of every meeting of the curve of `p` with
/// `v = 0.3 u^2`.
fn meetings(p: &[f64]) -> Result<Vec<(f64, usize)>, ExactCurveRefusal> {
    let first = curve(&[], -0.1, 0.4);
    let second = curve(p, -40.0, 40.0);
    let tolerance = Tolerance::new(1e-5, 1e-9).unwrap();
    let mut hits = section_curve_curve_intersection2(
        &first,
        Interval::new(0.0, 1.0),
        &second,
        Interval::new(0.0, 1.0),
        tolerance,
    )?
    .into_iter()
    .map(|h| (h.point.x, h.multiplicity))
    .collect::<Vec<_>>();
    hits.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(hits)
}

fn displaced(factors: &[(f64, usize)], by: f64) -> Vec<f64> {
    let mut p = product(factors);
    p[0] -= by;
    p
}

fn assert_simple(found: Result<Vec<(f64, usize)>, ExactCurveRefusal>, expected: &[f64]) {
    let hits = match found {
        Ok(hits) => hits,
        Err(refusal) => {
            assert_eq!(refusal, ExactCurveRefusal::Undecided);
            return;
        }
    };
    assert_eq!(hits.len(), expected.len(), "{hits:?}");
    for ((u, m), x) in hits.iter().zip(expected) {
        assert_eq!(*m, 1, "{hits:?}");
        assert!((u - x).abs() < 1e-5, "{hits:?} vs {expected:?}");
    }
}

#[test]
fn a_flat_stretch_between_two_crossings_does_not_swallow_them() {
    let p = displaced(&[(29.0 / 64.0, 5), (33.0 / 64.0, 3)], 2f64.powi(-20));
    assert_simple(meetings(&p), &[0.297_355, 0.656_069]);
}

#[test]
fn the_same_on_the_other_side_of_the_origin() {
    let p = displaced(&[(-11.0 / 32.0, 5), (-23.0 / 64.0, 3)], 2f64.powi(-20));
    assert_simple(meetings(&p), &[-0.526_550, -0.172_673]);
}

#[test]
fn a_flat_stretch_before_a_crossing_is_not_a_root_of_order_seven() {
    let p = displaced(&[(35.0 / 64.0, 6), (38.0 / 64.0, 1)], 2f64.powi(-29));
    assert_simple(meetings(&p), &[0.614_040]);
}
