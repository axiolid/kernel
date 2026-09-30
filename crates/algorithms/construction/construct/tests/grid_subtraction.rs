//! Long chains of grid-aligned subtraction through `boolean_polyhedra_exact`
//! (#199).
//!
//! The depth-2 Menger sponge is 147 consecutive differences from a unit
//! cube. It used to refuse at subtraction 82: split points were built with
//! an f64 formula, so a cut through the plane `z = 1/3` could land at
//! `0.33333333333333326`, two splits of one exact point came back ULPs
//! apart, and a later split through the pair emitted a ring enclosing no
//! area that every probe ray met edge-on.
//!
//! Every oracle here is exact and shares no code with the boolean: the
//! volume is summed in dyadic arithmetic and compared for EQUALITY with a
//! cell decomposition of the input boxes, and closure is checked by
//! cancelling the oriented coverage of every edge line, which tolerates the
//! T-junctions a fragment boolean legitimately leaves.

use axiolid_construct::polyhedron::{boolean_polyhedra_exact, BooleanOp, Polyhedron};
use axiolid_core::Point3;
use axiolid_exact::{Arith, Dyadic};
use std::collections::{BTreeMap, BTreeSet};

/// An axis-aligned box, outward-wound.
fn box_solid(min: [f64; 3], max: [f64; 3]) -> Polyhedron {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    Polyhedron::new(vec![
        vec![p(x0, y0, z0), p(x0, y1, z0), p(x1, y1, z0), p(x1, y0, z0)],
        vec![p(x0, y0, z1), p(x1, y0, z1), p(x1, y1, z1), p(x0, y1, z1)],
        vec![p(x0, y0, z0), p(x1, y0, z0), p(x1, y0, z1), p(x0, y0, z1)],
        vec![p(x0, y1, z0), p(x0, y1, z1), p(x1, y1, z1), p(x1, y1, z0)],
        vec![p(x0, y0, z0), p(x0, y0, z1), p(x0, y1, z1), p(x0, y1, z0)],
        vec![p(x1, y0, z0), p(x1, y1, z0), p(x1, y1, z1), p(x1, y0, z1)],
    ])
    .expect("box is a valid solid")
}

type Void = ([f64; 3], [f64; 3]);

/// The void boxes of a unit Menger sponge, generated with the same f64
/// arithmetic as the `axiolid/benchmarks` harness that found #199.
///
/// That arithmetic is part of the test: it produces both
/// `0.6666666666666666` and `0.6666666666666667` for "two thirds", so the
/// input itself carries one-ULP slivers the boolean must keep exactly.
fn menger_voids(depth: u32) -> Vec<Void> {
    fn carve(out: &mut Vec<Void>, origin: [f64; 3], size: f64, depth: u32) {
        if depth == 0 {
            return;
        }
        let third = size / 3.0;
        for i in 0..3u32 {
            for j in 0..3u32 {
                for k in 0..3u32 {
                    let centred = u32::from(i == 1) + u32::from(j == 1) + u32::from(k == 1);
                    let corner = [
                        origin[0] + third * f64::from(i),
                        origin[1] + third * f64::from(j),
                        origin[2] + third * f64::from(k),
                    ];
                    if centred >= 2 {
                        out.push((
                            corner,
                            [corner[0] + third, corner[1] + third, corner[2] + third],
                        ));
                    } else {
                        carve(out, corner, third, depth - 1);
                    }
                }
            }
        }
    }
    let mut out = Vec::new();
    carve(&mut out, [0.0; 3], 1.0, depth);
    out
}

fn dy(v: f64) -> Dyadic {
    Dyadic::try_from_f64(v).expect("finite")
}

/// Six times the solid's volume, exactly: the divergence theorem over each
/// face, fanned from its first vertex, summed in dyadic arithmetic.
fn six_volume(solid: &Polyhedron) -> Dyadic {
    let mut total = Dyadic::zero();
    for face in solid.faces() {
        let p = |q: Point3| [dy(q.x), dy(q.y), dy(q.z)];
        let a = p(face[0]);
        for i in 1..face.len() - 1 {
            let (b, c) = (p(face[i]), p(face[i + 1]));
            let det = a[0]
                .mul(&b[1].mul(&c[2]).sub(&b[2].mul(&c[1])))
                .sub(&a[1].mul(&b[0].mul(&c[2]).sub(&b[2].mul(&c[0]))))
                .add(&a[2].mul(&b[0].mul(&c[1]).sub(&b[1].mul(&c[0]))));
            total = total.add(&det);
        }
    }
    total
}

/// Six times the exact volume of `host` minus every void, from a cell
/// decomposition on the input coordinates. Cell membership is a comparison
/// of stored doubles, so a one-ULP cell is counted exactly like any other.
fn six_volume_oracle(host: Void, voids: &[Void]) -> Dyadic {
    let mut axes: [Vec<f64>; 3] = Default::default();
    for (axis, coords) in axes.iter_mut().enumerate() {
        let mut set = BTreeSet::new();
        for (n, x) in std::iter::once(&host).chain(voids) {
            set.insert(n[axis].to_bits());
            set.insert(x[axis].to_bits());
        }
        *coords = set.into_iter().map(f64::from_bits).collect();
        coords.sort_by(f64::total_cmp);
    }
    let inside =
        |b: &Void, lo: [f64; 3], hi: [f64; 3]| (0..3).all(|k| b.0[k] <= lo[k] && hi[k] <= b.1[k]);
    let mut total = Dyadic::zero();
    for i in 0..axes[0].len() - 1 {
        for j in 0..axes[1].len() - 1 {
            for k in 0..axes[2].len() - 1 {
                let lo = [axes[0][i], axes[1][j], axes[2][k]];
                let hi = [axes[0][i + 1], axes[1][j + 1], axes[2][k + 1]];
                if !inside(&host, lo, hi) || voids.iter().any(|v| inside(v, lo, hi)) {
                    continue;
                }
                let cell = (0..3).fold(dy(1.0), |acc, a| acc.mul(&dy(hi[a]).sub(&dy(lo[a]))));
                total = total.add(&cell);
            }
        }
    }
    total.mul(&dy(6.0))
}

/// Every face encloses area: its exact vector area is not zero.
fn assert_no_collapsed_face(solid: &Polyhedron) {
    for face in solid.faces() {
        let p = |q: Point3| [dy(q.x), dy(q.y), dy(q.z)];
        let a = p(face[0]);
        let mut sum = [Dyadic::zero(), Dyadic::zero(), Dyadic::zero()];
        for i in 1..face.len() - 1 {
            let (b, c) = (p(face[i]), p(face[i + 1]));
            let u = [b[0].sub(&a[0]), b[1].sub(&a[1]), b[2].sub(&a[2])];
            let v = [c[0].sub(&a[0]), c[1].sub(&a[1]), c[2].sub(&a[2])];
            let cross = [
                u[1].mul(&v[2]).sub(&u[2].mul(&v[1])),
                u[2].mul(&v[0]).sub(&u[0].mul(&v[2])),
                u[0].mul(&v[1]).sub(&u[1].mul(&v[0])),
            ];
            for axis in 0..3 {
                sum[axis] = sum[axis].add(&cross[axis]);
            }
        }
        assert!(
            sum.iter()
                .any(|s| s.sign() != Some(axiolid_guarantees::Sign::Zero)),
            "face {face:?} encloses no area"
        );
    }
}

/// The shell is closed: on every line carrying an edge, the oriented
/// coverage of the edges cancels everywhere.
///
/// A closed oriented surface uses every stretch of every edge once in each
/// direction. Summing interval endpoints per line (rather than matching
/// whole edges) accepts a T-junction, where one long edge faces two short
/// ones, and still catches a gap of any length. Every edge here is
/// axis-aligned, which the check asserts.
fn assert_closed(solid: &Polyhedron) {
    let mut events: BTreeMap<(usize, u64, u64, u64), i64> = BTreeMap::new();
    for face in solid.faces() {
        for i in 0..face.len() {
            let (p, q) = (face[i], face[(i + 1) % face.len()]);
            let (p, q) = ([p.x, p.y, p.z], [q.x, q.y, q.z]);
            let moving: Vec<usize> = (0..3).filter(|&k| p[k] != q[k]).collect();
            assert_eq!(moving.len(), 1, "edge {p:?} -> {q:?} is not axis-aligned");
            let axis = moving[0];
            let (s, t) = ((axis + 1) % 3, (axis + 2) % 3);
            // `+ 0.0` folds -0.0 into 0.0 so equal coordinates share a key.
            let bits = |v: f64| (v + 0.0).to_bits();
            let line = (axis, bits(p[s]), bits(p[t]));
            let sign = if q[axis] > p[axis] { 1 } else { -1 };
            let (lo, hi) = if sign > 0 {
                (p[axis], q[axis])
            } else {
                (q[axis], p[axis])
            };
            let key = |v: f64| (line.0, line.1, line.2, bits(v));
            *events.entry(key(lo)).or_default() += sign;
            *events.entry(key(hi)).or_default() -= sign;
        }
    }
    let open: Vec<String> = events
        .iter()
        .filter(|(_, &d)| d != 0)
        .map(|((axis, s, t, v), d)| {
            let f = |b: &u64| f64::from_bits(*b);
            format!("axis {axis} line ({}, {}) at {} by {d}", f(s), f(t), f(v))
        })
        .collect();
    assert!(open.is_empty(), "shell is open: {open:#?}");
}

/// The #199 reproduction: all 147 depth-2 voids, one difference each.
#[test]
fn a_depth_two_menger_chain_completes_with_the_exact_volume() {
    let host: Void = ([0.0; 3], [1.0; 3]);
    let voids = menger_voids(2);
    assert_eq!(voids.len(), 147);

    let mut solid = box_solid(host.0, host.1);
    for (step, (min, max)) in voids.iter().enumerate() {
        solid = boolean_polyhedra_exact(&solid, &box_solid(*min, *max), BooleanOp::Difference)
            .unwrap_or_else(|e| panic!("subtraction {step} of 147 refused: {e}"));
    }

    assert_no_collapsed_face(&solid);
    assert_closed(&solid);
    // Equality, not closeness: both sides are exact, so any difference is
    // a wrong solid rather than rounding.
    assert_eq!(
        six_volume(&solid),
        six_volume_oracle(host, &voids),
        "depth-2 sponge volume is not the exact volume of its inputs"
    );
}

/// The root cause in one step: a cut of an axis-aligned edge by an
/// axis-aligned plane is a double, and the boolean must return that double.
///
/// `1/3` and `2/3` are not dyadic, so the naive `a + (b - a) * t` rounds
/// twice and misses: before #199 the first Menger void already produced
/// coordinates that none of the inputs carried.
#[test]
fn split_points_on_a_grid_are_input_coordinates() {
    let t = 1.0 / 3.0;
    let host = box_solid([0.0; 3], [1.0; 3]);
    let cutters = [
        ([t, t, -1.0], [2.0 * t, 2.0 * t, 2.0]),
        ([t, -1.0, t], [2.0 * t, 2.0, 2.0 * t]),
        ([-1.0, t, t], [2.0, 2.0 * t, 2.0 * t]),
    ];
    let grid: BTreeSet<u64> = [0.0, t, 2.0 * t, 1.0]
        .iter()
        .map(|v: &f64| v.to_bits())
        .collect();
    for (min, max) in cutters {
        let result = boolean_polyhedra_exact(&host, &box_solid(min, max), BooleanOp::Difference)
            .expect("one grid-aligned difference");
        for face in result.faces() {
            for p in face {
                for c in [p.x, p.y, p.z] {
                    assert!(
                        grid.contains(&c.to_bits()),
                        "split point {p:?} carries {c:e}, which is not a grid coordinate"
                    );
                }
            }
        }
    }
}
