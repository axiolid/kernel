//! Where two B-spline surfaces meet (#119, ADR 0077): every component,
//! found with a certificate, carried as a [`PairSection3`].
//!
//! Neither surface has an equation to read in the other's parameters, so
//! the search runs over pairs of rational Bezier sub-patches:
//!
//! 1. A pair whose control hulls' boxes are apart holds no section.
//! 2. A pair whose normal cones are apart (no normal of one parallel to any
//!    normal of the other) holds no closed loop of the section: every piece
//!    of it there leaves through an edge of one sub-patch (Sederberg and
//!    Meyers). The cones are certain: they span the Bernstein coefficient
//!    vectors of the normal's own polynomial.
//! 3. Anything else is split, down to a smallest size, where the pair is
//!    refused: the surfaces touch or come closer than the search resolves.
//!
//! The edges of every pair of step 2 are intersected with the other
//! sub-patch (curve against patch by the same hull pruning, then Newton on
//! three unknowns), which seeds every component. From each seed the curve
//! is followed with a step held to a few degrees of turning, each node
//! corrected onto both surfaces to the last bits, until it leaves a window
//! or closes on itself.

use axiolid_core::{Point2, Point3, Scalar, Vec3};
use axiolid_curve::{BSplineSurface, Carrier, PairNode, PairSection3};

use crate::spline_field::bezier_net;

/// A homogeneous control point.
type H = [Scalar; 4];

/// A rational Bezier sub-patch: its homogeneous net and parameter box.
#[derive(Debug, Clone)]
struct Patch {
    net: Vec<Vec<H>>,
    lo: Point2,
    hi: Point2,
}

impl Patch {
    fn p(&self) -> usize {
        self.net.len() - 1
    }

    fn q(&self) -> usize {
        self.net[0].len() - 1
    }

    fn points(&self) -> impl Iterator<Item = Vec3> + '_ {
        self.net
            .iter()
            .flatten()
            .map(|h| Vec3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]))
    }

    fn aabb(&self) -> (Vec3, Vec3) {
        let mut lo = Vec3::splat(Scalar::INFINITY);
        let mut hi = Vec3::splat(Scalar::NEG_INFINITY);
        for p in self.points() {
            lo = lo.min(p);
            hi = hi.max(p);
        }
        let pad = 1e-12 * (1.0 + lo.abs().max(hi.abs()).max_element());
        (lo - Vec3::splat(pad), hi + Vec3::splat(pad))
    }

    fn split_u(&self) -> (Patch, Patch) {
        let p = self.p();
        let q = self.q();
        let mut left = vec![vec![[0.0; 4]; q + 1]; p + 1];
        let mut right = vec![vec![[0.0; 4]; q + 1]; p + 1];
        for b in 0..=q {
            let column: Vec<H> = (0..=p).map(|a| self.net[a][b]).collect();
            let (l, r) = casteljau(&column, 0.5);
            for a in 0..=p {
                left[a][b] = l[a];
                right[a][b] = r[a];
            }
        }
        let m = 0.5 * (self.lo.x + self.hi.x);
        (
            Patch {
                net: left,
                lo: self.lo,
                hi: Point2::new(m, self.hi.y),
            },
            Patch {
                net: right,
                lo: Point2::new(m, self.lo.y),
                hi: self.hi,
            },
        )
    }

    fn split_v(&self) -> (Patch, Patch) {
        let mut left = Vec::with_capacity(self.net.len());
        let mut right = Vec::with_capacity(self.net.len());
        for row in &self.net {
            let (l, r) = casteljau(row, 0.5);
            left.push(l);
            right.push(r);
        }
        let m = 0.5 * (self.lo.y + self.hi.y);
        (
            Patch {
                net: left,
                lo: self.lo,
                hi: Point2::new(self.hi.x, m),
            },
            Patch {
                net: right,
                lo: Point2::new(self.lo.x, m),
                hi: self.hi,
            },
        )
    }

    fn split(&self) -> Vec<Patch> {
        let (a, b) = self.split_u();
        let (a1, a2) = a.split_v();
        let (b1, b2) = b.split_v();
        vec![a1, a2, b1, b2]
    }

    /// The patch's four edges as rational Bezier curves with their
    /// parameter lines: `(net, fixed parameter is u?, fixed value, range)`.
    fn edges(&self) -> Vec<Edge> {
        let p = self.p();
        let q = self.q();
        vec![
            Edge {
                net: self.net[0].clone(),
                along_u: false,
                fixed: self.lo.x,
                from: self.lo.y,
                to: self.hi.y,
            },
            Edge {
                net: self.net[p].clone(),
                along_u: false,
                fixed: self.hi.x,
                from: self.lo.y,
                to: self.hi.y,
            },
            Edge {
                net: (0..=p).map(|a| self.net[a][0]).collect(),
                along_u: true,
                fixed: self.lo.y,
                from: self.lo.x,
                to: self.hi.x,
            },
            Edge {
                net: (0..=p).map(|a| self.net[a][q]).collect(),
                along_u: true,
                fixed: self.hi.y,
                from: self.lo.x,
                to: self.hi.x,
            },
        ]
    }

    /// A cone `(axis, half-angle)` holding every normal of the patch, from
    /// the Bernstein coefficient vectors of `N = (W X_u - W_u X) x (W X_v -
    /// W_v X)`; `None` when they span half a sphere or more.
    fn normal_cone(&self) -> Option<(Vec3, Scalar)> {
        let (p, q) = (self.p(), self.q());
        if p == 0 || q == 0 {
            return None;
        }
        let comp = |k: usize| -> Bern {
            Bern {
                p,
                q,
                c: self.net.iter().flatten().map(|h| h[k]).collect(),
            }
        };
        let (x, y, z, w) = (comp(0), comp(1), comp(2), comp(3));
        let (wu, wv) = (w.du(), w.dv());
        let tangent = |f: &Bern, along_u: bool| -> Bern {
            // W f_u - W_u f, both of degree (2p - 1, 2q).
            if along_u {
                w.mul(&f.du()).sub(&wu.mul(f))
            } else {
                w.mul(&f.dv()).sub(&wv.mul(f))
            }
        };
        let (tux, tuy, tuz) = (tangent(&x, true), tangent(&y, true), tangent(&z, true));
        let (tvx, tvy, tvz) = (tangent(&x, false), tangent(&y, false), tangent(&z, false));
        let nx = tuy.mul(&tvz).sub(&tuz.mul(&tvy));
        let ny = tuz.mul(&tvx).sub(&tux.mul(&tvz));
        let nz = tux.mul(&tvy).sub(&tuy.mul(&tvx));
        let vectors: Vec<Vec3> = (0..nx.c.len())
            .map(|k| Vec3::new(nx.c[k], ny.c[k], nz.c[k]))
            .filter(|v| v.length() > 0.0)
            .collect();
        if vectors.is_empty() {
            return None;
        }
        let axis = vectors
            .iter()
            .fold(Vec3::ZERO, |acc, v| acc + v.normalize());
        if axis.length() == 0.0 {
            return None;
        }
        let axis = axis.normalize();
        let mut half = 0.0 as Scalar;
        for v in &vectors {
            let c = v.normalize().dot(axis).clamp(-1.0, 1.0);
            if c <= 0.0 {
                return None;
            }
            half = half.max(c.acos());
        }
        Some((axis, half + 1e-12))
    }
}

/// One edge of a sub-patch, a rational Bezier curve over `[from, to]` of
/// the free parameter, the other parameter `fixed`.
#[derive(Debug, Clone)]
struct Edge {
    net: Vec<H>,
    along_u: bool,
    fixed: Scalar,
    from: Scalar,
    to: Scalar,
}

impl Edge {
    fn uv(&self, s: Scalar) -> Point2 {
        let free = self.from + (self.to - self.from) * s;
        if self.along_u {
            Point2::new(free, self.fixed)
        } else {
            Point2::new(self.fixed, free)
        }
    }
}

/// de Casteljau split of a homogeneous control polygon at `s`.
fn casteljau(c: &[H], s: Scalar) -> (Vec<H>, Vec<H>) {
    let n = c.len();
    let mut w = c.to_vec();
    let mut left = vec![[0.0; 4]; n];
    let mut right = vec![[0.0; 4]; n];
    left[0] = w[0];
    right[n - 1] = w[n - 1];
    for k in 1..n {
        for i in 0..n - k {
            for d in 0..4 {
                w[i][d] = w[i][d] * (1.0 - s) + w[i + 1][d] * s;
            }
        }
        left[k] = w[0];
        right[n - 1 - k] = w[n - 1 - k];
    }
    (left, right)
}

/// A tensor Bernstein polynomial (direction only matters here, so the
/// local parameter scale is left out of derivatives).
#[derive(Debug, Clone)]
struct Bern {
    p: usize,
    q: usize,
    c: Vec<Scalar>,
}

fn binomial(n: usize, k: usize) -> Scalar {
    let mut r = 1.0;
    for i in 0..k {
        r = r * (n - i) as Scalar / (i + 1) as Scalar;
    }
    r
}

impl Bern {
    fn at(&self, a: usize, b: usize) -> Scalar {
        self.c[a * (self.q + 1) + b]
    }

    fn du(&self) -> Bern {
        let (p, q) = (self.p, self.q);
        let mut c = Vec::with_capacity(p * (q + 1));
        for a in 0..p {
            for b in 0..=q {
                c.push(p as Scalar * (self.at(a + 1, b) - self.at(a, b)));
            }
        }
        Bern { p: p - 1, q, c }
    }

    fn dv(&self) -> Bern {
        let (p, q) = (self.p, self.q);
        let mut c = Vec::with_capacity((p + 1) * q);
        for a in 0..=p {
            for b in 0..q {
                c.push(q as Scalar * (self.at(a, b + 1) - self.at(a, b)));
            }
        }
        Bern { p, q: q - 1, c }
    }

    fn mul(&self, o: &Bern) -> Bern {
        let (p, q) = (self.p + o.p, self.q + o.q);
        let mut c = vec![0.0; (p + 1) * (q + 1)];
        for a1 in 0..=self.p {
            for b1 in 0..=self.q {
                let x = self.at(a1, b1) * binomial(self.p, a1) * binomial(self.q, b1);
                if x == 0.0 {
                    continue;
                }
                for a2 in 0..=o.p {
                    for b2 in 0..=o.q {
                        let y = o.at(a2, b2) * binomial(o.p, a2) * binomial(o.q, b2);
                        c[(a1 + a2) * (q + 1) + b1 + b2] += x * y;
                    }
                }
            }
        }
        for a in 0..=p {
            for b in 0..=q {
                c[a * (q + 1) + b] /= binomial(p, a) * binomial(q, b);
            }
        }
        Bern { p, q, c }
    }

    fn sub(&self, o: &Bern) -> Bern {
        Bern {
            p: self.p,
            q: self.q,
            c: self.c.iter().zip(&o.c).map(|(a, b)| a - b).collect(),
        }
    }
}

/// Why a pair trace was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PairRefusal {
    /// The surfaces touch, or come closer than the search resolves.
    Unresolved,
    /// A surface is not a clamped B-spline.
    Unsupported,
    /// The work budget ran out.
    Budget,
}

fn boxes_meet(a: &(Vec3, Vec3), b: &(Vec3, Vec3)) -> bool {
    a.0.x <= b.1.x
        && b.0.x <= a.1.x
        && a.0.y <= b.1.y
        && b.0.y <= a.1.y
        && a.0.z <= b.1.z
        && b.0.z <= a.1.z
}

fn patches(b: &BSplineSurface) -> Option<Vec<Patch>> {
    let (_, _, u_breaks, v_breaks, net) = bezier_net(b)?;
    let mut out = Vec::new();
    for (iu, row) in net.iter().enumerate() {
        for (jv, cell) in row.iter().enumerate() {
            out.push(Patch {
                net: cell.clone(),
                lo: Point2::new(u_breaks[iu], v_breaks[jv]),
                hi: Point2::new(u_breaks[iu + 1], v_breaks[jv + 1]),
            });
        }
    }
    Some(out)
}

/// The sub-patch pairs where the section lies, each free of closed loops.
fn resolve(first: &[Patch], second: &[Patch]) -> Result<Vec<(Patch, Patch)>, PairRefusal> {
    let mut queue: Vec<(Patch, Patch, u32)> = Vec::new();
    for a in first {
        for b in second {
            queue.push((a.clone(), b.clone(), 0));
        }
    }
    let mut out = Vec::new();
    let mut work = 0usize;
    while let Some((a, b, depth)) = queue.pop() {
        work += 1;
        if work > 200_000 {
            return Err(PairRefusal::Budget);
        }
        let (ba, bb) = (a.aabb(), b.aabb());
        if !boxes_meet(&ba, &bb) {
            continue;
        }
        if let (Some((na, ha)), Some((nb, hb))) = (a.normal_cone(), b.normal_cone()) {
            let angle = na.dot(nb).clamp(-1.0, 1.0).acos();
            if angle - ha - hb > 1e-9 && (core::f64::consts::PI - angle) - ha - hb > 1e-9 {
                out.push((a, b));
                continue;
            }
        }
        if depth > 24 {
            return Err(PairRefusal::Unresolved);
        }
        // Split the larger.
        let size = |x: &(Vec3, Vec3)| (x.1 - x.0).length();
        if size(&ba) >= size(&bb) {
            for piece in a.split() {
                queue.push((piece, b.clone(), depth + 1));
            }
        } else {
            for piece in b.split() {
                queue.push((a.clone(), piece, depth + 1));
            }
        }
    }
    Ok(out)
}

/// Points where an edge of one sub-patch crosses the other sub-patch, as
/// nodes: hull pruning down to small pieces, then Newton on three unknowns.
fn edge_hits(
    edge: &Edge,
    patch: &Patch,
    s1: &Carrier,
    s2: &Carrier,
    edge_on_first: bool,
) -> Vec<PairNode> {
    let mut out = Vec::new();
    let mut queue: Vec<(Vec<H>, Scalar, Scalar, Patch, u32)> =
        vec![(edge.net.clone(), 0.0, 1.0, patch.clone(), 0)];
    let curve_box = |net: &[H]| {
        let mut lo = Vec3::splat(Scalar::INFINITY);
        let mut hi = Vec3::splat(Scalar::NEG_INFINITY);
        for h in net {
            let p = Vec3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]);
            lo = lo.min(p);
            hi = hi.max(p);
        }
        let pad = 1e-12 * (1.0 + lo.abs().max(hi.abs()).max_element());
        (lo - Vec3::splat(pad), hi + Vec3::splat(pad))
    };
    let mut work = 0;
    while let Some((net, s0, s1_, sub, depth)) = queue.pop() {
        work += 1;
        if work > 20_000 {
            break;
        }
        let (bc, bp) = (curve_box(&net), sub.aabb());
        if !boxes_meet(&bc, &bp) {
            continue;
        }
        let small = (bc.1 - bc.0).length() + (bp.1 - bp.0).length()
            <= 1e-3 * (1.0 + bc.0.abs().max(bc.1.abs()).max_element());
        if small || depth > 30 {
            // Newton on (s, u, v): C(s) = S(u, v).
            let mut s = 0.5 * (s0 + s1_);
            let mut uv = (sub.lo + sub.hi) * 0.5;
            let (edge_surface, other) = if edge_on_first { (s1, s2) } else { (s2, s1) };
            let mut ok = false;
            for _ in 0..60 {
                let e = edge.uv(s);
                let je = edge_surface.jet(e.x, e.y);
                let tangent = if edge.along_u { je.u } else { je.v } * (edge.to - edge.from);
                let jo = other.jet(uv.x, uv.y);
                let r = je.point - jo.point;
                // Solve [tangent, -S_u, -S_v] (ds, du, dv) = -r.
                let m = [
                    [tangent.x, -jo.u.x, -jo.v.x],
                    [tangent.y, -jo.u.y, -jo.v.y],
                    [tangent.z, -jo.u.z, -jo.v.z],
                ];
                let Some(x) = solve3(m, [-r.x, -r.y, -r.z]) else {
                    break;
                };
                s += x[0];
                uv += Point2::new(x[1], x[2]);
                if x.iter().map(|v| v.abs()).fold(0.0, Scalar::max)
                    <= 4.0 * Scalar::EPSILON * (1.0 + s.abs() + uv.length())
                {
                    ok = true;
                    break;
                }
            }
            // Rounding can keep the steps above the stopping test at a
            // root: accept a residual at rounding level.
            if !ok && s.is_finite() && uv.is_finite() {
                let e = edge.uv(s);
                let p = edge_surface.jet(e.x, e.y).point;
                ok = (p - other.jet(uv.x, uv.y).point).length() <= 1e-12 * (1.0 + p.length());
            }
            let inside =
                |x: Scalar, a: Scalar, b: Scalar| x >= a.min(b) - 1e-9 && x <= a.max(b) + 1e-9;
            if ok
                && inside(s, 0.0, 1.0)
                && inside(uv.x, patch.lo.x, patch.hi.x)
                && inside(uv.y, patch.lo.y, patch.hi.y)
            {
                let e = edge.uv(s.clamp(0.0, 1.0));
                let p = edge_surface.jet(e.x, e.y).point;
                let node = if edge_on_first {
                    PairNode {
                        point: p,
                        first: e,
                        second: uv,
                    }
                } else {
                    PairNode {
                        point: p,
                        first: uv,
                        second: e,
                    }
                };
                if !out.iter().any(|n: &PairNode| {
                    (n.point - node.point).length() <= 1e-9 * (1.0 + p.length())
                }) {
                    out.push(node);
                }
            }
            continue;
        }
        // Split whichever is larger.
        if (bc.1 - bc.0).length() >= (bp.1 - bp.0).length() {
            let (l, r) = casteljau(&net, 0.5);
            let m = 0.5 * (s0 + s1_);
            queue.push((l, s0, m, sub.clone(), depth + 1));
            queue.push((r, m, s1_, sub, depth + 1));
        } else {
            for piece in sub.split() {
                queue.push((net.clone(), s0, s1_, piece, depth + 1));
            }
        }
    }
    out
}

fn solve3(m: [[Scalar; 3]; 3], r: [Scalar; 3]) -> Option<[Scalar; 3]> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det == 0.0 || !det.is_finite() {
        return None;
    }
    let col = |k: usize| -> Scalar {
        let mut mm = m;
        for row in 0..3 {
            mm[row][k] = r[row];
        }
        mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1])
            - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0])
            + mm[0][2] * (mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0])
    };
    Some([col(0) / det, col(1) / det, col(2) / det])
}

/// Correct a guess onto both surfaces on the plane through `target` across
/// `direction`: Newton on `S1(a) = S2(b)`, `direction . (S1(a) - target) = 0`.
fn correct(
    s1: &Carrier,
    s2: &Carrier,
    mut a: Point2,
    mut b: Point2,
    target: Point3,
    direction: Vec3,
) -> Option<PairNode> {
    for _ in 0..40 {
        let (j1, j2) = (s1.jet(a.x, a.y), s2.jet(b.x, b.y));
        let gap = j1.point - j2.point;
        let plane = direction.dot(j1.point - target);
        let m = [
            [j1.u.x, j1.v.x, -j2.u.x, -j2.v.x],
            [j1.u.y, j1.v.y, -j2.u.y, -j2.v.y],
            [j1.u.z, j1.v.z, -j2.u.z, -j2.v.z],
            [direction.dot(j1.u), direction.dot(j1.v), 0.0, 0.0],
        ];
        let x = axiolid_curve::pair_section::solve4(m, [-gap.x, -gap.y, -gap.z, -plane])?;
        a += Point2::new(x[0], x[1]);
        b += Point2::new(x[2], x[3]);
        if x.iter().map(|v| v.abs()).fold(0.0, Scalar::max)
            <= 4.0 * Scalar::EPSILON * (1.0 + a.length() + b.length())
        {
            let p = s1.jet(a.x, a.y).point;
            let q = s2.jet(b.x, b.y).point;
            return ((p - q).length() <= 1e-10 * (1.0 + p.length())).then_some(PairNode {
                point: p,
                first: a,
                second: b,
            });
        }
    }
    None
}

/// The curve's unit tangent at a node: the two normals' cross product.
fn tangent_at(s1: &Carrier, s2: &Carrier, n: &PairNode) -> Option<Vec3> {
    let (j1, j2) = (s1.jet(n.first.x, n.first.y), s2.jet(n.second.x, n.second.y));
    let t = j1.u.cross(j1.v).cross(j2.u.cross(j2.v));
    let l = t.length();
    (l > 0.0 && l.is_finite()).then(|| t / l)
}

/// Whether a node's parameters lie in both windows.
fn inside(n: &PairNode, w1: (Point2, Point2), w2: (Point2, Point2)) -> bool {
    let within = |p: Point2, w: (Point2, Point2)| {
        p.x >= w.0.x && p.x <= w.1.x && p.y >= w.0.y && p.y <= w.1.y
    };
    within(n.first, w1) && within(n.second, w2)
}

/// A bound of a window: which surface (`true` for the first), which
/// parameter (`0` for `u`), and its value.
#[derive(Debug, Clone, Copy)]
struct Bound {
    first: bool,
    axis: usize,
    value: Scalar,
}

/// The first window bound the straight move from `(a0, b0)` to `(a1, b1)`
/// crosses, with the fraction of the move where it does.
fn crossing(
    (a0, b0): (Point2, Point2),
    (a1, b1): (Point2, Point2),
    w1: (Point2, Point2),
    w2: (Point2, Point2),
) -> Option<(Bound, Scalar)> {
    let mut best: Option<(Bound, Scalar)> = None;
    for (first, from, to, w) in [(true, a0, a1, w1), (false, b0, b1, w2)] {
        for axis in 0..2 {
            let (x0, x1) = (from[axis], to[axis]);
            for value in [w.0[axis], w.1[axis]] {
                let outside = if value == w.0[axis] {
                    x1 < value
                } else {
                    x1 > value
                };
                if !outside || x1 == x0 {
                    continue;
                }
                let f = ((value - x0) / (x1 - x0)).clamp(0.0, 1.0);
                if best.is_none_or(|(_, g)| f < g) {
                    best = Some((Bound { first, axis, value }, f));
                }
            }
        }
    }
    best
}

/// Correct a guess onto both surfaces with one parameter held on a window
/// bound: Newton on `S1(a) = S2(b)` and that parameter's value.
fn correct_on(
    s1: &Carrier,
    s2: &Carrier,
    mut a: Point2,
    mut b: Point2,
    bound: Bound,
) -> Option<PairNode> {
    let pin = |a: &mut Point2, b: &mut Point2| {
        let p = if bound.first { a } else { b };
        p[bound.axis] = bound.value;
    };
    pin(&mut a, &mut b);
    for _ in 0..40 {
        let (j1, j2) = (s1.jet(a.x, a.y), s2.jet(b.x, b.y));
        let gap = j1.point - j2.point;
        let mut row = [0.0; 4];
        row[usize::from(!bound.first) * 2 + bound.axis] = 1.0;
        let m = [
            [j1.u.x, j1.v.x, -j2.u.x, -j2.v.x],
            [j1.u.y, j1.v.y, -j2.u.y, -j2.v.y],
            [j1.u.z, j1.v.z, -j2.u.z, -j2.v.z],
            row,
        ];
        let x = axiolid_curve::pair_section::solve4(m, [-gap.x, -gap.y, -gap.z, 0.0])?;
        a += Point2::new(x[0], x[1]);
        b += Point2::new(x[2], x[3]);
        pin(&mut a, &mut b);
        if x.iter().map(|v| v.abs()).fold(0.0, Scalar::max)
            <= 4.0 * Scalar::EPSILON * (1.0 + a.length() + b.length())
        {
            let p = s1.jet(a.x, a.y).point;
            let q = s2.jet(b.x, b.y).point;
            return ((p - q).length() <= 1e-10 * (1.0 + p.length())).then_some(PairNode {
                point: p,
                first: a,
                second: b,
            });
        }
    }
    None
}

/// Whether a node's parameters lie in both windows, to rounding.
fn inside_loosely(n: &PairNode, w1: (Point2, Point2), w2: (Point2, Point2)) -> bool {
    let within = |p: Point2, w: (Point2, Point2)| {
        let e = 1e-12 * (1.0 + w.0.abs().max(w.1.abs()).max_element());
        p.x >= w.0.x - e && p.x <= w.1.x + e && p.y >= w.0.y - e && p.y <= w.1.y + e
    };
    within(n.first, w1) && within(n.second, w2)
}

/// Follow the curve from `seed` in direction `sign`, until it leaves a
/// window (its last node then lies on the window's edge) or comes back to
/// `seed`; returns the nodes after the seed and whether it closed.
fn march(
    s1: &Carrier,
    s2: &Carrier,
    seed: PairNode,
    sign: Scalar,
    w1: (Point2, Point2),
    w2: (Point2, Point2),
    scale: Scalar,
) -> Option<(Vec<PairNode>, bool)> {
    let mut out = Vec::new();
    let mut here = seed;
    let mut t = tangent_at(s1, s2, &here)? * sign;
    let mut h = 0.02 * scale;
    let min_h = 1e-9 * scale;
    for _ in 0..20_000 {
        if h < min_h {
            return None;
        }
        // Rates of both parameter pairs along the step, from the system.
        let (j1, j2) = (
            s1.jet(here.first.x, here.first.y),
            s2.jet(here.second.x, here.second.y),
        );
        let project = |j: &axiolid_curve::SurfaceJet| {
            let (a, b, c) = (j.u.dot(j.u), j.u.dot(j.v), j.v.dot(j.v));
            let det = a * c - b * b;
            let (g0, g1) = (j.u.dot(t), j.v.dot(t));
            Point2::new((c * g0 - b * g1) / det, (a * g1 - b * g0) / det)
        };
        let (da, db) = (project(&j1), project(&j2));
        let guess = (here.first + da * h, here.second + db * h);
        // A step past a window's edge ends on it, solved there exactly.
        let exit = |to: (Point2, Point2)| -> Option<Option<PairNode>> {
            let (bound, f) = crossing((here.first, here.second), to, w1, w2)?;
            let at = (
                here.first + (to.0 - here.first) * f,
                here.second + (to.1 - here.second) * f,
            );
            let n = correct_on(s1, s2, at.0, at.1, bound).filter(|n| {
                inside_loosely(n, w1, w2)
                    && (n.point - here.point).dot(t) >= 0.0
                    && (n.point - here.point).length() <= 2.0 * h
            });
            Some(n)
        };
        let finish = |n: PairNode, mut out: Vec<PairNode>| {
            if (n.point - here.point).length() > 1e-12 * scale {
                out.push(n);
            }
            Some((out, false))
        };
        match exit(guess) {
            Some(Some(n)) => return finish(n, out),
            Some(None) => {
                h *= 0.5;
                continue;
            }
            None => {}
        }
        let target = here.point + t * h;
        let Some(next) = correct(s1, s2, guess.0, guess.1, target, t) else {
            h *= 0.5;
            continue;
        };
        let nt = tangent_at(s1, s2, &next)?;
        let nt = if nt.dot(t) < 0.0 { -nt } else { nt };
        // Turning held to a few degrees per step, and the step a true
        // advance.
        if nt.dot(t) < (6.0f64).to_radians().cos() || (next.point - here.point).dot(t) <= 0.0 {
            h *= 0.5;
            continue;
        }
        // Back at the seed: a closed loop.
        if out.len() > 3
            && (next.point - seed.point).length() <= h * 1.01
            && (seed.point - here.point).dot(t) > 0.0
        {
            return Some((out, true));
        }
        if !inside(&next, w1, w2) {
            match exit((next.first, next.second)) {
                Some(Some(n)) => return finish(n, out),
                _ => {
                    h *= 0.5;
                    continue;
                }
            }
        }
        out.push(next);
        here = next;
        t = nt;
        h = (h * 1.5).min(0.05 * scale);
    }
    None
}

/// Every component of the section of two B-spline surfaces within windows
/// of their parameters (their whole domains when `None`).
///
/// # Errors
///
/// `Unsupported` for an unclamped or malformed surface; `Unresolved` where
/// the surfaces touch or come closer than the search resolves; `Budget`.
pub(crate) fn pair_trace(
    b1: &BSplineSurface,
    b2: &BSplineSurface,
    windows: Option<((Point2, Point2), (Point2, Point2))>,
) -> Result<Vec<PairSection3>, PairRefusal> {
    let (s1, s2) = (
        Carrier::Spline(Box::new(b1.clone())),
        Carrier::Spline(Box::new(b2.clone())),
    );
    // The windows, never past the surfaces' own domains.
    let d1 = b1.domain().ok_or(PairRefusal::Unsupported)?;
    let d2 = b2.domain().ok_or(PairRefusal::Unsupported)?;
    let whole1 = (Point2::new(d1.0 .0, d1.1 .0), Point2::new(d1.0 .1, d1.1 .1));
    let whole2 = (Point2::new(d2.0 .0, d2.1 .0), Point2::new(d2.0 .1, d2.1 .1));
    let clip = |w: (Point2, Point2), d: (Point2, Point2)| (w.0.max(d.0), w.1.min(d.1));
    let (w1, w2) = match windows {
        Some((a, b)) => (clip(a, whole1), clip(b, whole2)),
        None => (whole1, whole2),
    };
    let (p1, p2) = (
        patches(b1).ok_or(PairRefusal::Unsupported)?,
        patches(b2).ok_or(PairRefusal::Unsupported)?,
    );
    let pairs = resolve(&p1, &p2)?;
    // Seeds: every edge of every resolved pair against the other patch.
    let mut seeds: Vec<PairNode> = Vec::new();
    for (a, b) in &pairs {
        for e in a.edges() {
            seeds.extend(edge_hits(&e, b, &s1, &s2, true));
        }
        for e in b.edges() {
            seeds.extend(edge_hits(&e, a, &s1, &s2, false));
        }
    }
    seeds.retain(|n| inside(n, w1, w2));
    let scale = {
        let mut lo = Vec3::splat(Scalar::INFINITY);
        let mut hi = Vec3::splat(Scalar::NEG_INFINITY);
        for p in p1.iter().flat_map(|p| p.points().collect::<Vec<_>>()) {
            lo = lo.min(p);
            hi = hi.max(p);
        }
        (hi - lo).length().max(1e-9)
    };
    let mut out: Vec<PairSection3> = Vec::new();
    // A seed on a curve already traced: the curve itself passes through it
    // (chords sag, so the test is on the curve, not its polygon).
    let covered = |out: &[PairSection3], n: &PairNode| {
        out.iter().any(|c| {
            c.parameter_of(n.point)
                .and_then(|t| c.point(t))
                .is_some_and(|p| (p - n.point).length() <= 1e-7 * scale)
        })
    };
    for seed in seeds {
        if covered(&out, &seed) {
            continue;
        }
        let (forward, closed) =
            march(&s1, &s2, seed, 1.0, w1, w2, scale).ok_or(PairRefusal::Unresolved)?;
        let mut nodes = Vec::new();
        if closed {
            nodes.push(seed);
            nodes.extend(forward);
            nodes.push(seed);
        } else {
            let (backward, _) =
                march(&s1, &s2, seed, -1.0, w1, w2, scale).ok_or(PairRefusal::Unresolved)?;
            nodes.extend(backward.into_iter().rev());
            nodes.push(seed);
            nodes.extend(forward);
        }
        if nodes.len() >= 2 {
            out.push(PairSection3 {
                first: s1.clone(),
                second: s2.clone(),
                nodes,
            });
        }
    }
    Ok(out)
}

/// A rational Bezier curve's point and derivative at `s` in `[0, 1]`.
fn bezier_jet(net: &[H], s: Scalar) -> (Point3, Vec3) {
    let n = net.len() - 1;
    let (left, _) = casteljau(net, s);
    // The point is the last left control; the derivative from the last two
    // of the degree n - 1 level: n (P_n - P_{n-1}) in homogeneous terms.
    let mut w = net.to_vec();
    for k in 1..n {
        for i in 0..=n - k {
            for d in 0..4 {
                w[i][d] = w[i][d] * (1.0 - s) + w[i + 1][d] * s;
            }
        }
    }
    let (a, b) = (w[0], w[1]);
    let h = left[n];
    let point = Vec3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]);
    let dh = [
        n as Scalar * (b[0] - a[0]),
        n as Scalar * (b[1] - a[1]),
        n as Scalar * (b[2] - a[2]),
        n as Scalar * (b[3] - a[3]),
    ];
    let derivative = (Vec3::new(dh[0], dh[1], dh[2]) - point * dh[3]) / h[3];
    (point, derivative)
}

/// Where a B-spline curve meets a B-spline surface: the curve's parameter
/// and the point, for every crossing (hull pruning over Bezier pieces of
/// both, then Newton on three unknowns). `None` for unclamped operands.
pub(crate) fn spline_curve_surface_hits(
    curve: &axiolid_curve::BSplineCurve3,
    surface: &BSplineSurface,
) -> Option<Vec<(Scalar, Point3)>> {
    let control: Vec<H> = curve
        .control_points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let w = curve.weights.as_ref().map_or(1.0, |ws| ws[i]);
            [p.x * w, p.y * w, p.z * w, w]
        })
        .collect();
    let knots = {
        let mut out = Vec::new();
        for (&k, &m) in curve.knots.iter().zip(&curve.multiplicities) {
            out.extend(core::iter::repeat_n(k, m as usize));
        }
        out
    };
    let (breaks, segments) =
        crate::spline_field::decompose(&knots, usize::from(curve.degree), &control)?;
    let carrier = Carrier::Spline(Box::new(surface.clone()));
    let mut out: Vec<(Scalar, Point3)> = Vec::new();
    for (i, seg) in segments.iter().enumerate() {
        let (k0, k1) = (breaks[i], breaks[i + 1]);
        if k1 <= k0 {
            continue;
        }
        for patch in patches(surface)? {
            let mut queue = vec![(seg.clone(), 0.0, 1.0, patch.clone(), 0u32)];
            let mut work = 0;
            while let Some((net, s0, s1, sub, depth)) = queue.pop() {
                work += 1;
                if work > 20_000 {
                    break;
                }
                let cb = {
                    let mut lo = Vec3::splat(Scalar::INFINITY);
                    let mut hi = Vec3::splat(Scalar::NEG_INFINITY);
                    for h in &net {
                        let p = Vec3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]);
                        lo = lo.min(p);
                        hi = hi.max(p);
                    }
                    let pad = 1e-12 * (1.0 + lo.abs().max(hi.abs()).max_element());
                    (lo - Vec3::splat(pad), hi + Vec3::splat(pad))
                };
                let pb = sub.aabb();
                if !boxes_meet(&cb, &pb) {
                    continue;
                }
                let small = (cb.1 - cb.0).length() + (pb.1 - pb.0).length()
                    <= 1e-3 * (1.0 + cb.0.abs().max(cb.1.abs()).max_element());
                if small || depth > 30 {
                    let mut s = 0.5;
                    let mut uv = (sub.lo + sub.hi) * 0.5;
                    let mut ok = false;
                    for _ in 0..60 {
                        let (c, dc) = bezier_jet(seg, s0 + (s1 - s0) * s);
                        let dc = dc * (s1 - s0);
                        let j = carrier.jet(uv.x, uv.y);
                        let r = c - j.point;
                        let m = [
                            [dc.x, -j.u.x, -j.v.x],
                            [dc.y, -j.u.y, -j.v.y],
                            [dc.z, -j.u.z, -j.v.z],
                        ];
                        let Some(x) = solve3(m, [-r.x, -r.y, -r.z]) else {
                            break;
                        };
                        s += x[0];
                        uv += Point2::new(x[1], x[2]);
                        if x.iter().map(|v| v.abs()).fold(0.0, Scalar::max)
                            <= 4.0 * Scalar::EPSILON * (1.0 + s.abs() + uv.length())
                        {
                            ok = true;
                            break;
                        }
                    }
                    // Rounding can keep the steps above the stopping
                    // test at a root: accept a residual at rounding level.
                    if !ok {
                        let (c, _) = bezier_jet(seg, s0 + (s1 - s0) * s);
                        let miss = (c - carrier.jet(uv.x, uv.y).point).length();
                        ok = miss <= 1e-12 * (1.0 + c.length());
                    }
                    let local = s0 + (s1 - s0) * s;
                    if ok
                        && (-1e-9..=1.0 + 1e-9).contains(&local)
                        && uv.x >= patch.lo.x - 1e-9
                        && uv.x <= patch.hi.x + 1e-9
                        && uv.y >= patch.lo.y - 1e-9
                        && uv.y <= patch.hi.y + 1e-9
                    {
                        let t = k0 + (k1 - k0) * local.clamp(0.0, 1.0);
                        let (p, _) = bezier_jet(seg, local.clamp(0.0, 1.0));
                        if !out
                            .iter()
                            .any(|(tt, _)| (tt - t).abs() <= 1e-9 * (1.0 + t.abs()))
                        {
                            out.push((t, p));
                        }
                    }
                    continue;
                }
                if (cb.1 - cb.0).length() >= (pb.1 - pb.0).length() {
                    let (l, r) = casteljau(&net, 0.5);
                    let m = 0.5 * (s0 + s1);
                    queue.push((l, s0, m, sub.clone(), depth + 1));
                    queue.push((r, m, s1, sub, depth + 1));
                } else {
                    for piece in sub.split() {
                        queue.push((net.clone(), s0, s1, piece, depth + 1));
                    }
                }
            }
        }
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    Some(out)
}

/// Every component of the section of two B-spline surfaces within windows
/// of their parameters (their whole domains when `None`), each a
/// [`PairSection3`] whose first surface is `first` (ADR 0077).
///
/// Every component is found: the search splits Bezier sub-patch pairs until
/// their normal cones are apart, where no closed loop can hide, and seeds
/// from every crossing of a sub-patch edge.
///
/// # Errors
///
/// - `UnsupportedPair`: an unclamped or malformed surface, or the search
///   ran out of budget.
/// - `NotRegularCurve`: the surfaces touch, or come closer than the search
///   resolves.
/// - `Disjoint`: no section in the windows.
pub fn spline_pair_intersection(
    first: &BSplineSurface,
    second: &BSplineSurface,
    windows: Option<((Point2, Point2), (Point2, Point2))>,
) -> Result<Vec<PairSection3>, crate::ExactIntersectionRefusal> {
    use crate::ExactIntersectionRefusal as R;
    let curves = pair_trace(first, second, windows).map_err(|refusal| match refusal {
        PairRefusal::Unresolved => R::NotRegularCurve,
        PairRefusal::Unsupported | PairRefusal::Budget => R::UnsupportedPair,
    })?;
    if curves.is_empty() {
        return Err(R::Disjoint);
    }
    Ok(curves)
}
