//! `Curve2::Intrinsic` resolves through the generic 2D curve API.
//!
//! ADR 0060 added the type and a direct evaluator, but `evaluate2`,
//! `derivative2` and `domain2` still refused it by name, so every generic
//! consumer -- flattening, profiles, sweeps -- inherited that refusal.
//! These pin the dispatch, not the quadrature: correctness of the
//! integrator itself is pinned against Fresnel in `arc_length.rs`.

use axiolid_core::{Frame2, Point2, Scalar, Vec2};
use axiolid_curve::{CurvatureLaw, Curve2, Intrinsic2};
use axiolid_evaluate::curve::domain2;
use axiolid_evaluate::{derivative2, evaluate2, flatten2};

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(0.0, 0.0),
        x: Vec2::X,
        y: Vec2::Y,
    }
}

/// A circular law of curvature `1/r` is an arc of radius `r`.
///
/// The closed form is elementary here, which is what makes the generic
/// path checkable without trusting the quadrature under test.
fn arc(radius: Scalar, length: Scalar) -> Curve2 {
    Curve2::Intrinsic(Intrinsic2::new(
        frame(),
        CurvatureLaw::circular(1.0 / radius),
        length,
    ))
}

#[test]
fn the_generic_api_dispatches_an_intrinsic_curve() {
    // Quarter arc of radius 4: the generic entry point must return the
    // arc's own point, not a refusal.
    let radius = 4.0;
    let quarter = radius * core::f64::consts::FRAC_PI_2;
    let curve = arc(radius, quarter);

    let end = evaluate2(&curve, quarter).expect("evaluate2 must dispatch Intrinsic");
    // Start frame has x along the tangent, so the centre is at +y and the
    // quarter arc ends at (r, r).
    let want = Point2::new(radius, radius);
    let error = (end - want).length();
    assert!(
        error < 1e-9,
        "quarter arc ended at {end:?}, want {want:?} (error {error:e})"
    );
}

#[test]
fn the_generic_derivative_is_a_unit_tangent() {
    let curve = arc(4.0, 6.0);
    // Arc length parameterisation means |dC/ds| == 1 everywhere, which a
    // parametric arm would not satisfy.
    for step in 0..=6 {
        let s = 6.0 * Scalar::from(step) / 6.0;
        let d = derivative2(&curve, s).expect("derivative2 must dispatch Intrinsic");
        assert!(
            (d.length() - 1.0).abs() < 1e-12,
            "|d| = {} at s={s}",
            d.length()
        );
    }
}

#[test]
fn the_domain_is_the_declared_arc_length() {
    let curve = arc(4.0, 7.5);
    let domain = domain2(&curve);
    assert!((domain.start - 0.0).abs() < 1e-12);
    assert!(
        (domain.end - 7.5).abs() < 1e-12,
        "domain end {}",
        domain.end
    );
}

#[test]
fn a_length_that_claims_no_domain_is_reported_as_empty() {
    // A non-positive length has no knowable span; claiming the unit
    // interval would invent one. Zero and NEGATIVE are both checked: a
    // zero-length curve yields an empty interval under either the guarded
    // or the unguarded rule, so only a negative length actually
    // discriminates a missing `length > 0.0` guard.
    for length in [0.0, -5.0, Scalar::NAN] {
        let curve = arc(4.0, length);
        let domain = domain2(&curve);
        assert!(
            (domain.end - domain.start).abs() < 1e-12,
            "length {length} must claim no domain, got [{}, {}]",
            domain.start,
            domain.end
        );
        assert!(
            domain.start == 0.0 && domain.end == 0.0,
            "an unknowable domain must be the empty interval at zero, got [{}, {}]",
            domain.start,
            domain.end
        );
    }
}

#[test]
fn a_clothoid_flattens_through_the_generic_path() {
    // The payoff of dispatch: a transition spiral becomes points for any
    // downstream consumer. `rate` is the sharpness, k(s) = rate * s.
    let rate = 1.0 / (300.0 * 120.0);
    let length = 120.0;
    let clothoid = Curve2::Intrinsic(Intrinsic2::new(
        frame(),
        CurvatureLaw::Polynomial {
            coefficients: vec![0.0, rate],
        },
        length,
    ));

    let points = flatten2(&clothoid, domain2(&clothoid), 1e-3, 24)
        .expect("a clothoid must flatten, not refuse");
    assert!(
        points.len() > 2,
        "a curved spiral needs interior points, got {}",
        points.len()
    );

    // Every flattened point must lie on the curve itself.
    let first = points.first().copied().expect("first");
    let last = points.last().copied().expect("last");
    let start = evaluate2(&clothoid, 0.0).expect("start");
    let end = evaluate2(&clothoid, length).expect("end");
    assert!(
        (first - start).length() < 1e-9,
        "flatten start {first:?} vs {start:?}"
    );
    assert!(
        (last - end).length() < 1e-9,
        "flatten end {last:?} vs {end:?}"
    );
}

#[test]
fn a_straight_law_is_exactly_the_line_through_the_generic_path() {
    // Zero curvature must reduce to the start tangent exactly: a
    // degenerate case the quadrature must not perturb.
    let curve = Curve2::Intrinsic(Intrinsic2::new(frame(), CurvatureLaw::straight(), 10.0));
    let point = evaluate2(&curve, 7.0).expect("straight law");
    assert!(
        (point - Point2::new(7.0, 0.0)).length() < 1e-12,
        "got {point:?}"
    );
}
