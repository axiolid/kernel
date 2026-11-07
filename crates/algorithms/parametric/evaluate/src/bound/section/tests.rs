//! The plane-curve and space-curve enclosures behind the lifted bounds:
//! dense samples of the derivatives never exceed them.

use super::*;
use crate::curve::{derivative2, derivative3, evaluate2, second_derivative2, second_derivative3};
use axiolid_core::{Frame3, Point2, Point3, Vec2, Vec3};
use axiolid_curve::{
    AngleGraph2, Axis, Basis, Branch, Curve2, Field2, ImplicitCell, ImplicitCurve2,
    ImplicitSection3, QuadraticGraph2, RuledCarrier, RuledSection3, SeriesField2, TorusCarrier,
    TorusSection3, Trig2,
};

const SAMPLES: usize = 400;

fn trig(constant: Scalar, cos: Scalar, sin: Scalar, cos2: Scalar, sin2: Scalar) -> Trig2 {
    Trig2 {
        constant,
        cos,
        sin,
        cos2,
        sin2,
    }
}

fn quadratic(branch: Branch) -> QuadraticGraph2 {
    QuadraticGraph2 {
        a: trig(1.0, 0.2, 0.0, 0.0, 0.0),
        b: trig(0.5, 0.0, -1.5, 0.25, 0.0),
        c: trig(-2.0, 0.75, 0.0, 0.0, 0.5),
        branch,
    }
}

fn angle(branch: Branch) -> AngleGraph2 {
    AngleGraph2 {
        a: trig(1.0, 0.2, 0.0, 0.0, 0.0),
        b: trig(0.0, 0.0, 0.5, 0.1, 0.0),
        c: trig(0.3, 0.0, 0.0, 0.0, 0.2),
        branch,
    }
}

fn frame() -> Frame3 {
    Frame3 {
        origin: Point3::new(0.3, -0.2, 0.5),
        x: Vec3::new(1.0, 0.2, -0.1),
        y: Vec3::new(0.15, 0.9, 0.3),
        z: Vec3::new(0.1, -0.2, 1.1),
    }
}

fn at(lo: Scalar, hi: Scalar, i: usize) -> Scalar {
    lo + (hi - lo) * i as Scalar / SAMPLES as Scalar
}

/// The sampled derivatives of a plane curve lie within `bounds`.
fn assert_plane(curve: &Curve2, bounds: &PlaneBounds, lo: Scalar, hi: Scalar, what: &str) {
    for i in 0..=SAMPLES {
        // The second derivative jumps at a cell join: read the end of a
        // span from inside it.
        let t = at(lo, hi, i);
        let t = if t == hi && hi.fract() == 0.0 {
            t.next_down()
        } else {
            t
        };
        let (p, d, dd) = (
            evaluate2(curve, t).unwrap(),
            derivative2(curve, t).unwrap(),
            second_derivative2(curve, t).unwrap(),
        );
        assert!(
            bounds.lo.x <= p.x && p.x <= bounds.hi.x && bounds.lo.y <= p.y && p.y <= bounds.hi.y,
            "{what} at {t}: {p:?} outside {:?}..{:?}",
            bounds.lo,
            bounds.hi
        );
        for (value, bound, name) in [
            (d.x, bounds.first[0], "first x"),
            (d.y, bounds.first[1], "first y"),
            (dd.x, bounds.second[0], "second x"),
            (dd.y, bounds.second[1], "second y"),
        ] {
            assert!(
                value.abs() <= bound,
                "{what} at {t}: {name} {value} over {bound}"
            );
        }
    }
}

#[test]
fn graph_bounds_cover_the_sampled_derivatives() {
    for branch in [Branch::Plus, Branch::Minus] {
        for (lo, hi) in [(0.0, 0.05), (0.2, 0.3), (-0.5, -0.4), (0.6, 0.75)] {
            let g = quadratic(branch);
            let b = quadratic_bounds(&g, lo, hi).expect("certified");
            assert_plane(&Curve2::QuadraticGraph(g), &b, lo, hi, "quadratic");
            let g = angle(branch);
            let b = angle_bounds(&g, lo, hi).expect("certified");
            assert_plane(&Curve2::AngleGraph(g), &b, lo, hi, "angle");
        }
    }
}

#[test]
fn a_graph_bound_over_a_branch_end_or_an_empty_span_is_refused() {
    let g = QuadraticGraph2 {
        a: trig(1.0, 0.0, 0.0, 0.0, 0.0),
        b: Trig2::default(),
        c: trig(-0.5, 0.5, 0.0, 0.0, 0.0),
        branch: Branch::Plus,
    };
    assert!(quadratic_bounds(&g, 0.0, 0.2).is_none());
    assert!(quadratic_bounds(&quadratic(Branch::Plus), 0.3, 0.3).is_none());
    assert!(angle_bounds(&angle(Branch::Plus), Scalar::NAN, 0.3).is_none());
}

/// The unit circle `u^2 + v^2 - 1` as an implicit curve: a regular cell
/// along `u` (upper arc), then a bridge cell.
fn circle() -> ImplicitCurve2 {
    let field = Field2::Series(SeriesField2 {
        u: Basis::Power,
        v: Basis::Power,
        coefficients: vec![vec![-1.0, 0.0, 1.0], vec![0.0], vec![1.0]],
    });
    let regular = ImplicitCell {
        axis: Axis::U,
        from: -0.6,
        to: 0.6,
        low: 0.5,
        high: 1.5,
        bridge: None,
    };
    let bridge = ImplicitCell::bridge(
        Point2::new(0.6, 0.8),
        Point2::new(0.9, 0.3),
        Vec2::new(0.8, -0.6),
        Vec2::new(1.0, -1.5),
    );
    ImplicitCurve2 {
        field,
        cells: vec![regular, bridge],
    }
}

#[test]
fn implicit_bounds_cover_the_sampled_derivatives_in_both_kinds_of_cell() {
    let c = circle();
    let curve = Curve2::Implicit(c.clone());
    for (lo, hi) in [(0.0, 0.3), (0.5, 1.0), (1.0, 1.4), (1.5, 2.0), (0.25, 0.75)] {
        let b = derivative_bounds_checked(&c, lo, hi);
        assert_plane(&curve, &b, lo, hi, "implicit");
    }
    // Across the join the two cells are not one curve of bounded bend.
    assert!(implicit::derivative_bounds(&c, 0.5, 1.5).is_none());
    assert!(implicit::derivative_bounds(&c, 0.5, 2.5).is_none());
    assert!(implicit::derivative_bounds(&c, 1.0, 1.0).is_none());
}

fn derivative_bounds_checked(c: &ImplicitCurve2, lo: Scalar, hi: Scalar) -> PlaneBounds {
    implicit::derivative_bounds(c, lo, hi).unwrap_or_else(|| panic!("[{lo}, {hi}] refused"))
}

/// `v - u^2 = 0` over `u` in `[-0.5, 0.5]`: one regular cell.
fn parabola() -> ImplicitCurve2 {
    ImplicitCurve2 {
        field: Field2::Series(SeriesField2 {
            u: Basis::Power,
            v: Basis::Power,
            coefficients: vec![vec![0.0, 1.0], vec![0.0], vec![-1.0]],
        }),
        cells: vec![ImplicitCell {
            axis: Axis::U,
            from: -0.5,
            to: 0.5,
            low: -0.1,
            high: 0.4,
            bridge: None,
        }],
    }
}

/// The sampled space derivatives lie within `bounds`.
fn assert_space(curve: &Curve3, lo: Scalar, hi: Scalar, what: &str) -> SpaceBounds {
    let b = space_bounds(curve, lo, hi).unwrap_or_else(|| panic!("{what}: no bounds"));
    let second = b
        .second
        .unwrap_or_else(|| panic!("{what}: no second bound"));
    let (mut c1, mut c2): (Scalar, Scalar) = (0.0, 0.0);
    for i in 0..=SAMPLES {
        let t = at(lo, hi, i);
        c1 = c1.max(derivative3(curve, t).unwrap().length());
        c2 = c2.max(second_derivative3(curve, t).unwrap().length());
    }
    assert!(c1 <= b.first, "{what}: |c'| {c1} over {}", b.first);
    assert!(c2 <= second, "{what}: |c''| {c2} over {second}");
    // Not slack beyond what the carrier's frame and the interval
    // enclosures cost.
    assert!(
        b.first <= 8.0 * c1 + 1e-12,
        "{what}: |c'| bound {}",
        b.first
    );
    assert!(second <= 40.0 * c2 + 1e-9, "{what}: |c''| bound {second}");
    b
}

#[test]
fn section_bounds_cover_the_sampled_space_derivatives() {
    let ruled = |branch| {
        Curve3::RuledSection(RuledSection3 {
            carrier: RuledCarrier {
                frame: frame(),
                x_radius: 0.7,
                y_radius: 0.5,
                slope: 0.3,
            },
            graph: quadratic(branch),
        })
    };
    let torus = |branch| {
        Curve3::TorusSection(TorusSection3 {
            torus: TorusCarrier {
                frame: frame(),
                major_radius: 0.9,
                minor_radius: 0.3,
            },
            graph: angle(branch),
        })
    };
    for branch in [Branch::Plus, Branch::Minus] {
        for (lo, hi) in [(0.2, 0.3), (-0.5, -0.4), (0.6, 0.75)] {
            assert_space(&ruled(branch), lo, hi, "ruled section");
            assert_space(&torus(branch), lo, hi, "torus section");
        }
    }
    let parabola = Curve3::ImplicitSection(ImplicitSection3 {
        carrier: Carrier::Plane(frame()),
        curve: parabola(),
    });
    for (lo, hi) in [(0.0, 0.4), (0.2, 0.9)] {
        assert_space(&parabola, lo, hi, "parabola on a plane");
    }
    // A bridge cell of the circle on a sphere.
    let circle = Curve3::ImplicitSection(ImplicitSection3 {
        carrier: Carrier::Sphere {
            frame: frame(),
            radius: 0.8,
        },
        curve: circle(),
    });
    for (lo, hi) in [(0.1, 0.4), (1.0, 1.4), (1.5, 2.0)] {
        assert_space(&circle, lo, hi, "circle on a sphere");
    }
}

#[test]
fn a_section_across_a_cell_join_or_without_a_bound_has_none() {
    let section = Curve3::ImplicitSection(ImplicitSection3 {
        carrier: Carrier::Plane(frame()),
        curve: circle(),
    });
    assert!(space_bounds(&section, 0.5, 1.5).is_none());
    assert!(space_bounds(&section, 0.4, 0.2).is_none());
    // A line is covered by the curve derivative bounds.
    let line = Curve3::Line(axiolid_curve::Line3 {
        origin: Point3::ZERO,
        direction: Vec3::X,
    });
    let b = space_bounds(&line, 0.0, 1.0).unwrap();
    assert!(b.first >= 1.0 && b.second == Some(0.0));
}

/// `v - k u = 0` over `u` in `[from, to]`, the cell's bracket the whole
/// image: one regular cell of rate `to - from` in `t`.
fn line_field(k: Scalar, from: Scalar, to: Scalar, low: Scalar, high: Scalar) -> ImplicitCurve2 {
    ImplicitCurve2 {
        field: Field2::Series(SeriesField2 {
            u: Basis::Power,
            v: Basis::Power,
            coefficients: vec![vec![0.0, 1.0], vec![-k]],
        }),
        cells: vec![ImplicitCell {
            axis: Axis::U,
            from,
            to,
            low,
            high,
            bridge: None,
        }],
    }
}

#[test]
fn an_implicit_cell_bounds_scale_with_its_rate_and_the_mixed_partial_counts_twice() {
    // `w = u^2` over `u` in `[-1.5, 1.5]`: a rate of 3, so `q'' = 9 (0, 2)`
    // exactly, which the bound is within rounding of.
    let wide = ImplicitCurve2 {
        cells: vec![ImplicitCell {
            from: -1.5,
            to: 1.5,
            low: -0.1,
            high: 2.5,
            ..parabola().cells[0]
        }],
        ..parabola()
    };
    let b = implicit::derivative_bounds(&wide, 0.0, 1.0).unwrap();
    assert!(b.second[1] >= 18.0 && b.second[1] <= 18.5, "{:?}", b.second);
    assert!(b.second[0] == 0.0 && b.first[0] >= 3.0 && b.first[0] <= 3.0001);

    // A line `v = 1.5 u` on a cone near its apex: `c''` is
    // `r^2 (S_uu + 2 S_uv 1.5)`, its radial and tangential parts both
    // counted, which is more than the radial part and the mixed partial once.
    let (rate, slope, k) = (0.6, 1.0, 1.5);
    let section = Curve3::ImplicitSection(ImplicitSection3 {
        carrier: Carrier::Ruled(RuledCarrier {
            frame: Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            x_radius: 0.05,
            y_radius: 0.05,
            slope,
        }),
        curve: line_field(k, 0.2, 0.2 + rate, 0.2, 1.4),
    });
    let sampled = (0..=SAMPLES)
        .map(|i| second_derivative3(&section, i as Scalar / SAMPLES as Scalar).unwrap())
        .map(|v| v.length())
        .fold(0.0, Scalar::max);
    let got = space_bounds(&section, 0.0, 1.0).unwrap().second.unwrap();
    assert!(sampled <= got, "{sampled} over {got}");
}
