//! Certified chord bounds of the pcurves of a pair section (#414).
//!
//! Two hand-built Bezier patches meet along a curve known in closed form:
//! `S1 = (2u - 1, v, (2u - 1)^2)`, a parabolic cylinder (quadratic in `u`),
//! and `S2 = (2u - 1, v, 1/4 + v/2)`, a plane (bilinear). They meet where
//! `x^2 = 1/4 + y/2`, `x = 2u - 1`, `y = v`, on both patches at the same
//! parameters. Its nodes are at `v = 0, 1/4, ..., 1`, and the pcurve at `t`
//! is where the plane across the chord of the node's chord crosses the curve,
//! found here by bisection, independently of the section's own Newton solve.

use axiolid_core::{Point2, Point3, Scalar, Vec3};
use axiolid_curve::{
    BSplineSurface, Carrier, Curve2, Curve3, KnotSpec, LiftedCurve2, PairNode, PairSection3,
};
use axiolid_evaluate::bound::{chord_bound2, continuity_breaks2};

const SAMPLES: usize = 4000;
const USEFUL: Scalar = 10.0;

fn patch(
    u_degree: u16,
    v_degree: u16,
    rows: Vec<Vec<Point3>>,
    u_domain: (Scalar, Scalar),
) -> BSplineSurface {
    let (p, q) = (usize::from(u_degree), usize::from(v_degree));
    assert_eq!((rows.len(), rows[0].len()), (p + 1, q + 1));
    BSplineSurface {
        u_degree,
        v_degree,
        control_points: rows,
        u_knots: vec![u_domain.0, u_domain.1],
        u_multiplicities: vec![u32::from(u_degree) + 1; 2],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![u32::from(v_degree) + 1; 2],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    }
}

/// `(2u - 1, v, (2u - 1)^2)`: control `x = -1, 0, 1`, `z = 1, -1, 1`.
fn cylinder(u_domain: (Scalar, Scalar)) -> BSplineSurface {
    let row = |x: Scalar, z: Scalar| vec![Point3::new(x, 0.0, z), Point3::new(x, 1.0, z)];
    patch(
        2,
        1,
        vec![row(-1.0, 1.0), row(0.0, -1.0), row(1.0, 1.0)],
        u_domain,
    )
}

/// `(2u - 1, v, 1/4 + v/2)`.
fn plane() -> BSplineSurface {
    let row = |x: Scalar| vec![Point3::new(x, 0.0, 0.25), Point3::new(x, 1.0, 0.75)];
    patch(1, 1, vec![row(-1.0), row(1.0)], (0.0, 1.0))
}

/// The curve's point at height `v`, and its parameters on both patches.
fn node(v: Scalar) -> PairNode {
    let x = (0.25 + 0.5 * v).sqrt();
    let u = (1.0 + x) / 2.0;
    PairNode {
        point: Point3::new(x, v, x * x),
        first: Point2::new(u, v),
        second: Point2::new(u, v),
    }
}

fn heights() -> [Scalar; 5] {
    [0.0, 0.25, 0.5, 0.75, 1.0]
}

fn section(u_domain: (Scalar, Scalar)) -> PairSection3 {
    PairSection3 {
        first: Carrier::Spline(Box::new(cylinder(u_domain))),
        second: Carrier::Spline(Box::new(plane())),
        nodes: heights().map(node).to_vec(),
    }
}

fn lifted(section: &PairSection3, first: bool) -> Curve2 {
    Curve2::Lifted(LiftedCurve2 {
        curve: Box::new(Curve3::PairSection(section.clone())),
        carrier: if first {
            section.first.clone()
        } else {
            section.second.clone()
        },
        start: 0.0,
        end: 4.0,
        guide: Vec::new(),
    })
}

/// The pcurve at `t`: the point of the curve whose projection on the chord
/// is the fraction `t - i` of it, by bisection on the height.
fn point(t: Scalar) -> Point2 {
    let nodes: Vec<PairNode> = heights().map(node).to_vec();
    let i = (t.floor() as usize).min(3);
    let s = t - i as Scalar;
    let (p0, p1) = (nodes[i].point, nodes[i + 1].point);
    let d: Vec3 = p1 - p0;
    let f = |v: Scalar| {
        let x = (0.25 + 0.5 * v).sqrt();
        d.dot(Point3::new(x, v, x * x) - p0) - s * d.dot(d)
    };
    let (mut lo, mut hi) = (nodes[i].first.y, nodes[i + 1].first.y);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let v = 0.5 * (lo + hi);
    let x = (0.25 + 0.5 * v).sqrt();
    Point2::new((1.0 + x) / 2.0, v)
}

fn sampled(lo: Scalar, hi: Scalar) -> Scalar {
    let (p, q) = (point(lo), point(hi));
    let chord = q - p;
    (1..SAMPLES)
        .map(|i| {
            let t = lo + (hi - lo) * i as Scalar / SAMPLES as Scalar;
            let r = point(t) - p;
            let s = (r.dot(chord) / chord.length_squared()).clamp(0.0, 1.0);
            (r - chord * s).length()
        })
        .fold(0.0, Scalar::max)
}

#[test]
fn the_section_locus_is_the_curve_on_both_patches() {
    // The independent pcurve really is the pair's: the section reads the
    // same parameters at every sample.
    let section = section((0.0, 1.0));
    for first in [true, false] {
        let curve = lifted(&section, first);
        for i in 0..=200 {
            let t = 4.0 * i as Scalar / 200.0;
            let read = axiolid_evaluate::curve::evaluate2(&curve, t).unwrap();
            let truth = point(t);
            assert!((read - truth).length() < 1e-9, "{t}: {read:?} {truth:?}");
        }
    }
}

#[test]
fn a_pair_section_bound_holds_the_sampled_pcurve_and_stays_close() {
    let section = section((0.0, 1.0));
    let spans = [
        (0.1, 0.9),
        (0.0, 1.0),
        (1.2, 1.6),
        (2.0, 3.0),
        (3.4, 3.5),
        (2.0, 2.0001),
    ];
    for first in [true, false] {
        let curve = lifted(&section, first);
        for (lo, hi) in spans {
            let sample = sampled(lo, hi);
            let bound = chord_bound2(&curve, lo, hi)
                .unwrap_or_else(|| panic!("first {first} [{lo}, {hi}]: no bound"));
            assert!(
                sample <= bound,
                "first {first} [{lo}, {hi}]: sampled {sample} above {bound}"
            );
            let reference = sample.max((hi - lo) * (hi - lo) / 8.0);
            assert!(
                bound <= USEFUL * reference + 1e-12,
                "first {first} [{lo}, {hi}]: {bound} over {USEFUL}x {reference}"
            );
        }
    }
}

#[test]
fn a_span_across_a_node_or_ending_past_the_patch_domain_is_refused() {
    let whole = section((0.0, 1.0));
    let curve = lifted(&whole, true);
    assert_eq!(continuity_breaks2(&curve, 1), vec![1.0, 2.0, 3.0]);
    assert!(chord_bound2(&curve, 0.5, 1.5).is_none());
    assert!(chord_bound2(&curve, 0.0, 4.0).is_none());
    assert!(chord_bound2(&curve, 3.5, 4.5).is_none());
    // An evaluated end outside a patch's knot domain is not certified, on
    // either side of the section.
    let mut out = section((0.0, 1.0));
    out.nodes[4].first = Point2::new(1.2, 1.0);
    for first in [true, false] {
        assert!(chord_bound2(&lifted(&out, first), 3.5, 4.0).is_none());
        assert!(chord_bound2(&lifted(&out, first), 3.0, 3.5).is_some());
    }
}
