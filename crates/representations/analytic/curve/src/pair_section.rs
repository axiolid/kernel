//! Where two parametric surfaces meet, carried by nodes on both (ADR 0077).
//!
//! Two B-spline surfaces have no implicit equation to read one in the
//! other's parameters, so their section is carried directly: a chain of
//! [`PairNode`]s, each a point on both surfaces with its parameters on each,
//! found to the last bits. Between two nodes the curve is defined, not
//! interpolated: at local parameter `s` it is the point where both surfaces
//! meet on the plane across the chord at `P0 + s (P1 - P0)`, the unique
//! solution of four equations in the four parameters near the chord --
//! `S1(a) = S2(b)` and the plane -- found by Newton from the nodes'
//! parameters. Derivatives follow from the same system.

use axiolid_core::{Point2, Point3, Scalar, Vec3};

use crate::implicit::{Carrier, SurfaceJet};

/// A point on both surfaces, with its parameters on each.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct PairNode {
    /// The point.
    pub point: Point3,
    /// Its parameters on the first surface.
    pub first: Point2,
    /// Its parameters on the second surface.
    pub second: Point2,
}

/// A stretch of the curve where two surfaces meet. The parameter runs over
/// `[0, nodes.len() - 1]`, one unit per chord.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct PairSection3 {
    /// The first surface.
    pub first: Carrier,
    /// The second surface.
    pub second: Carrier,
    /// Nodes on both, in order along the curve.
    pub nodes: Vec<PairNode>,
}

/// Solve the 4x4 system `m x = r` by Gaussian elimination with partial
/// pivoting; `None` when it is singular.
#[must_use]
#[allow(clippy::needless_range_loop)]
pub fn solve4(mut m: [[Scalar; 4]; 4], mut r: [Scalar; 4]) -> Option<[Scalar; 4]> {
    for col in 0..4 {
        let pivot = (col..4).max_by(|&a, &b| m[a][col].abs().total_cmp(&m[b][col].abs()))?;
        if m[pivot][col] == 0.0 || !m[pivot][col].is_finite() {
            return None;
        }
        m.swap(col, pivot);
        r.swap(col, pivot);
        for row in col + 1..4 {
            let f = m[row][col] / m[col][col];
            for k in col..4 {
                m[row][k] -= f * m[col][k];
            }
            r[row] -= f * r[col];
        }
    }
    let mut x = [0.0; 4];
    for row in (0..4).rev() {
        let mut acc = r[row];
        for k in row + 1..4 {
            acc -= m[row][k] * x[k];
        }
        x[row] = acc / m[row][row];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// The system's matrix at `(a, b)` for chord direction `d`.
fn matrix(j1: &SurfaceJet, j2: &SurfaceJet, d: Vec3) -> [[Scalar; 4]; 4] {
    [
        [j1.u.x, j1.v.x, -j2.u.x, -j2.v.x],
        [j1.u.y, j1.v.y, -j2.u.y, -j2.v.y],
        [j1.u.z, j1.v.z, -j2.u.z, -j2.v.z],
        [d.dot(j1.u), d.dot(j1.v), 0.0, 0.0],
    ]
}

impl PairSection3 {
    /// The parameter range's end, `nodes.len() - 1`.
    #[must_use]
    pub fn end(&self) -> Scalar {
        self.nodes.len().saturating_sub(1) as Scalar
    }

    fn locate(&self, t: Scalar) -> Option<(usize, Scalar)> {
        let n = self.nodes.len();
        if n < 2 || !t.is_finite() {
            return None;
        }
        let slack = 1e-12 * (1.0 + self.end());
        if t < -slack || t > self.end() + slack {
            return None;
        }
        let i = (t.floor().max(0.0) as usize).min(n - 2);
        Some((i, (t - i as Scalar).clamp(0.0, 1.0)))
    }

    /// The parameters on both surfaces at `t`, and the point.
    #[must_use]
    pub fn solve(&self, t: Scalar) -> Option<(Point2, Point2, Point3)> {
        let (i, s) = self.locate(t)?;
        let (n0, n1) = (self.nodes[i], self.nodes[i + 1]);
        if s == 0.0 {
            return Some((n0.first, n0.second, n0.point));
        }
        if s == 1.0 {
            return Some((n1.first, n1.second, n1.point));
        }
        let d = n1.point - n0.point;
        let target = n0.point + d * s;
        let mut a = n0.first + (n1.first - n0.first) * s;
        let mut b = n0.second + (n1.second - n0.second) * s;
        let scale = 1.0 + n0.point.length() + d.length();
        for _ in 0..50 {
            let j1 = self.first.jet(a.x, a.y);
            let j2 = self.second.jet(b.x, b.y);
            let gap = j1.point - j2.point;
            let plane = d.dot(j1.point - target);
            let x = solve4(matrix(&j1, &j2, d), [-gap.x, -gap.y, -gap.z, -plane])?;
            a += Point2::new(x[0], x[1]);
            b += Point2::new(x[2], x[3]);
            let step = x.iter().map(|v| v.abs()).fold(0.0, Scalar::max);
            if step <= 4.0 * Scalar::EPSILON * (1.0 + a.length() + b.length()) {
                break;
            }
        }
        let (j1, j2) = (self.first.jet(a.x, a.y), self.second.jet(b.x, b.y));
        ((j1.point - j2.point).length() <= 1e-9 * scale).then_some((a, b, j1.point))
    }

    /// The point at `t`.
    #[must_use]
    pub fn point(&self, t: Scalar) -> Option<Point3> {
        self.solve(t).map(|(_, _, p)| p)
    }

    /// The rates of both surfaces' parameters and of the point at `t`:
    /// from `S1_a a' - S2_b b' = 0` and `d . S1_a a' = |d|^2`.
    #[must_use]
    pub fn rates(&self, t: Scalar) -> Option<(Point2, Point2, Vec3)> {
        let (i, _) = self.locate(t)?;
        let d = self.nodes[i + 1].point - self.nodes[i].point;
        let (a, b, _) = self.solve(t)?;
        let (j1, j2) = (self.first.jet(a.x, a.y), self.second.jet(b.x, b.y));
        let x = solve4(matrix(&j1, &j2, d), [0.0, 0.0, 0.0, d.dot(d)])?;
        let (da, db) = (Point2::new(x[0], x[1]), Point2::new(x[2], x[3]));
        Some((da, db, j1.u * da.x + j1.v * da.y))
    }

    /// `dP/dt`.
    #[must_use]
    pub fn tangent(&self, t: Scalar) -> Option<Vec3> {
        self.rates(t).map(|(_, _, d)| d)
    }

    /// The second rates of both surfaces' parameters and of the point at
    /// `t`, by central differences of [`Self::rates`] inside the chord.
    #[must_use]
    pub fn second_rates(&self, t: Scalar) -> Option<(Point2, Point2, Vec3)> {
        let (i, s) = self.locate(t)?;
        let h = 1e-5;
        let (lo, hi) = ((s - h).max(0.0), (s + h).min(1.0));
        let base = i as Scalar;
        let (p, q) = (self.rates(base + lo)?, self.rates(base + hi)?);
        let w = hi - lo;
        Some(((q.0 - p.0) / w, (q.1 - p.1) / w, (q.2 - p.2) / w))
    }

    /// `d2P/dt2`.
    #[must_use]
    pub fn bend(&self, t: Scalar) -> Option<Vec3> {
        self.second_rates(t).map(|(_, _, d)| d)
    }

    /// Which of the two surfaces `carrier` is: `Some(true)` for the first,
    /// `Some(false)` for the second.
    #[must_use]
    pub fn side(&self, carrier: &Carrier) -> Option<bool> {
        if *carrier == self.first {
            Some(true)
        } else if *carrier == self.second {
            Some(false)
        } else {
            None
        }
    }

    /// The stretch from `t0` to `t1` (`t0 < t1`) as a curve of its own over
    /// `[0, nodes]`: the nodes between, and new end nodes solved at `t0` and
    /// `t1`.
    #[must_use]
    pub fn sub(&self, t0: Scalar, t1: Scalar) -> Option<Self> {
        let (a0, b0, p0) = self.solve(t0)?;
        let (a1, b1, p1) = self.solve(t1)?;
        let mut nodes = vec![PairNode {
            point: p0,
            first: a0,
            second: b0,
        }];
        let first = t0.floor() as usize + 1;
        let last = t1.ceil() as usize;
        for i in first..last.min(self.nodes.len()) {
            let t = i as Scalar;
            if t > t0 + 1e-9 && t < t1 - 1e-9 {
                nodes.push(self.nodes[i]);
            }
        }
        nodes.push(PairNode {
            point: p1,
            first: a1,
            second: b1,
        });
        Some(Self {
            first: self.first.clone(),
            second: self.second.clone(),
            nodes,
        })
    }

    /// The same curve run backwards.
    #[must_use]
    pub fn reversed(&self) -> Self {
        let mut nodes = self.nodes.clone();
        nodes.reverse();
        Self {
            first: self.first.clone(),
            second: self.second.clone(),
            nodes,
        }
    }

    /// The parameter of a point on the curve: the chord nearest it, the
    /// plane across that chord through the point, then Newton along the
    /// curve's tangent within the chord.
    #[must_use]
    pub fn parameter_of(&self, p: Point3) -> Option<Scalar> {
        let mut best: Option<(Scalar, usize, Scalar)> = None;
        for (i, w) in self.nodes.windows(2).enumerate() {
            let d = w[1].point - w[0].point;
            let l2 = d.length_squared();
            if l2 == 0.0 {
                continue;
            }
            let s = ((p - w[0].point).dot(d) / l2).clamp(0.0, 1.0);
            let miss = (w[0].point + d * s - p).length();
            if best.is_none_or(|(m, _, _)| miss < m) {
                best = Some((miss, i, s));
            }
        }
        let (_, i, mut s) = best?;
        let base = i as Scalar;
        for _ in 0..30 {
            let t = base + s;
            let (q, d) = (self.point(t)?, self.tangent(t)?);
            let l2 = d.length_squared();
            if l2 == 0.0 {
                break;
            }
            let next = (s - (q - p).dot(d) / l2).clamp(0.0, 1.0);
            let step = (next - s).abs();
            s = next;
            if step <= 4.0 * Scalar::EPSILON * (1.0 + base) {
                break;
            }
        }
        Some(base + s)
    }

    /// Whether every number is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.first.is_finite()
            && self.second.is_finite()
            && self
                .nodes
                .iter()
                .all(|n| n.point.is_finite() && n.first.is_finite() && n.second.is_finite())
    }
}
