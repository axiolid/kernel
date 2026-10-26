//! Adversarial coverage for globally certified surface inversion (#7).
//!
//! Fixtures target the cases the issue names: poles (a whole parameter
//! family collapsing to one point), seams, and near-degenerate
//! parameterisations.

use axiolid_core::{Point3, Scalar, Tolerance};
use axiolid_curve::KnotSpec;
use axiolid_evaluate::surface::bspline_jet;
use axiolid_nurbs::{
    invert_periodic_surface_certified, invert_surface_certified, CertifiedSurfaceProjectionOptions,
    PeriodicBSplineSurface, SurfaceInversionRefusal,
};
use axiolid_surface::BSplineSurface;

fn options() -> CertifiedSurfaceProjectionOptions {
    CertifiedSurfaceProjectionOptions::new(Tolerance::new(1e-8, 1e-12).unwrap(), 1e-4, 100_000, 64)
        .unwrap()
}

fn plane() -> BSplineSurface {
    BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        control_points: vec![
            vec![Point3::new(-1.0, -1.0, 0.0), Point3::new(-1.0, 1.0, 0.0)],
            vec![Point3::new(1.0, -1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ],
        u_knots: vec![2.0, 5.0],
        u_multiplicities: vec![2, 2],
        v_knots: vec![-3.0, 1.0],
        v_multiplicities: vec![2, 2],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: Some(false),
    }
}

/// A cone-like patch whose entire v=0 row is the SAME apex point.
///
/// Every u names the apex, so inversion there has no unique answer. This is
/// the degenerate parameterisation #7 requires be refused, not guessed.
fn poled_patch() -> BSplineSurface {
    let apex = Point3::new(0.0, 0.0, 1.0);
    BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        control_points: vec![
            vec![apex, Point3::new(-1.0, 0.0, 0.0)],
            vec![apex, Point3::new(0.0, 1.0, 0.0)],
            vec![apex, Point3::new(1.0, 0.0, 0.0)],
        ],
        u_knots: vec![0.0, 1.0, 2.0],
        u_multiplicities: vec![2, 1, 2],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![2, 2],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: Some(false),
    }
}

fn periodic_ring() -> PeriodicBSplineSurface {
    let ring = [
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(0.0, -1.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    ];
    let control_points = ring
        .into_iter()
        .map(|point| vec![point, point + Point3::Z])
        .collect();
    PeriodicBSplineSurface::new(BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        control_points,
        u_knots: vec![-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        u_multiplicities: vec![1; 7],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![2, 2],
        weights: None,
        u_closed: true,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: Some(false),
    })
    .unwrap()
}

#[test]
fn an_interior_on_surface_point_inverts_uniquely() {
    let surface = plane();
    // Native domains are u in [2, 5] and v in [-3, 1]; pick a strict interior
    // point so no boundary tie can occur.
    let target = Point3::new(0.25, -0.5, 0.0);

    let certificate = invert_surface_certified(&surface, target, options())
        .expect("a valid query")
        .expect("an on-surface interior point has unique parameters");

    // The parameters must actually name the queried point.
    assert!(certificate.residual_upper_bound <= 1e-8);
    // Uniqueness means the enclosure really does own the answer.
    assert!(certificate.enclosure.u.start <= certificate.u);
    assert!(certificate.u <= certificate.enclosure.u.end);
    assert!(certificate.enclosure.v.start <= certificate.v);
    assert!(certificate.v <= certificate.enclosure.v.end);
    // The cover may be split across adjacent cells; what matters is that it
    // is ONE connected region, which the enclosure hull represents.
    assert!(!certificate.projection.possible_minimizer_boxes.is_empty());
}

#[test]
fn a_point_off_the_surface_is_refused_with_a_certified_lower_bound() {
    let surface = plane();
    // Well clear of the z=0 patch: this is a separation proof, not a
    // marginal residual.
    let target = Point3::new(0.0, 0.0, 2.0);

    let refusal = invert_surface_certified(&surface, target, options())
        .expect("a valid query")
        .expect_err("a point 2 units away is not ON the surface");

    let SurfaceInversionRefusal::OffSurface {
        distance_lower_bound,
    } = refusal
    else {
        panic!("distance, not ambiguity, is the reason here");
    };
    // The bound is a real proof of separation, not a token value.
    assert!(
        distance_lower_bound > 1.9,
        "lower bound was {distance_lower_bound}"
    );
}

#[test]
fn a_pole_is_reported_ambiguous_rather_than_resolved_arbitrarily() {
    let surface = poled_patch();
    // The apex: EVERY u maps here at v=0. Any single (u, v) would be an
    // arbitrary pick presented as a certified inverse.
    let apex = Point3::new(0.0, 0.0, 1.0);

    // A pole cannot be resolved to a tight parameter tolerance: the whole
    // u-family ties, so the search would subdivide forever. Use a coarse
    // parameter tolerance so the search TERMINATES and the ambiguity itself
    // is what gets reported.
    let coarse = CertifiedSurfaceProjectionOptions::new(
        Tolerance::new(1e-2, 1e-12).unwrap(),
        0.5,
        100_000,
        8,
    )
    .unwrap();

    let refusal = invert_surface_certified(&surface, apex, coarse)
        .expect("a valid query")
        .expect_err("an apex has no unique parameters");

    let SurfaceInversionRefusal::Ambiguous { candidates } = refusal else {
        panic!("the apex is ambiguous, and it is ON the surface");
    };
    // The ambiguity is handed back for inspection, not swallowed.
    assert!(candidates.len() > 1, "got {} candidates", candidates.len());
}

/// The bicubic bump of the edge-inversion report: a 5x5 net whose interior
/// 3x3 block is lifted, with four Bezier patches meeting at (0.5, 0.5).
fn bump() -> BSplineSurface {
    let control_points = (0..5)
        .map(|i| {
            (0..5)
                .map(|j| {
                    let (x, y) = (i as Scalar * 0.25, j as Scalar * 0.25);
                    let lift = if (1..4).contains(&i) && (1..4).contains(&j) {
                        0.12
                    } else {
                        0.0
                    };
                    Point3::new(x, y, lift + 0.05 * x * y)
                })
                .collect()
        })
        .collect();
    BSplineSurface {
        u_degree: 3,
        v_degree: 3,
        control_points,
        u_knots: vec![0.0, 0.5, 1.0],
        u_multiplicities: vec![4, 1, 4],
        v_knots: vec![0.0, 0.5, 1.0],
        v_multiplicities: vec![4, 1, 4],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    }
}

/// A point off a curved patch used to exhaust the budget refining every
/// retained box down to the parameter resolution, and so surfaced
/// `BudgetExceeded` instead of the separation it had already proven (#320).
#[test]
fn a_point_far_off_a_curved_patch_is_refused_off_surface_not_over_budget() {
    let surface = bump();
    let jet = bspline_jet(&surface, 0.3, 0.7).expect("an interior jet");
    let normal = jet.du.cross(jet.dv).normalize();
    for offset in [1e-3, 0.05, -0.05] {
        let target = jet.point + normal * offset;
        let refusal = invert_surface_certified(&surface, target, options())
            .expect("an off-surface point is a verdict, not an error")
            .expect_err("a point off the surface is not invertible");
        let SurfaceInversionRefusal::OffSurface {
            distance_lower_bound,
        } = refusal
        else {
            panic!("offset {offset}: expected OffSurface, got {refusal:?}");
        };
        assert!(
            distance_lower_bound > 1e-8 && distance_lower_bound <= offset.abs() + 1e-12,
            "offset {offset}: lower bound {distance_lower_bound}"
        );
    }
}

#[test]
fn a_point_off_a_periodic_surface_is_refused_off_surface_not_over_budget() {
    let surface = periodic_ring();
    // The ring's chord midpoint lies 1/sqrt(2) from the axis; (0.55, 0.55)
    // is 0.1/sqrt(2) outside it.
    let target = Point3::new(0.55, 0.55, 0.5);
    let refusal = invert_periodic_surface_certified(&surface, target, options())
        .expect("an off-surface point is a verdict, not an error")
        .expect_err("a point off the surface is not invertible");
    let SurfaceInversionRefusal::OffSurface {
        distance_lower_bound,
    } = refusal
    else {
        panic!("expected OffSurface, got {refusal:?}");
    };
    assert!(
        distance_lower_bound > 1e-8 && distance_lower_bound <= 0.1 / 2f64.sqrt() + 1e-12,
        "lower bound {distance_lower_bound}"
    );
}

#[test]
fn a_seam_point_inverts_on_the_quotient_domain() {
    let surface = periodic_ring();
    // Exactly on the seam of the closed U axis. On the quotient domain this
    // is ONE location, not two rival endpoint boxes.
    let seam = Point3::new(1.0, 0.0, 0.5);

    let outcome =
        invert_periodic_surface_certified(&surface, seam, options()).expect("a valid query");

    // Whatever the verdict, it must be SOUND: either unique parameters that
    // reproduce the seam point, or an explicit ambiguity naming the rival
    // regions. Never a silent pick.
    match outcome {
        Ok(certificate) => {
            assert!(certificate.residual_upper_bound <= 1e-8);
        }
        Err(SurfaceInversionRefusal::Ambiguous { candidates }) => {
            assert!(candidates.len() > 1);
        }
        other => panic!(
            "a seam point on its own surface must not be off-surface or unresolved: {other:?}"
        ),
    }
}
