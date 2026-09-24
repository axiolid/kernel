//! Boundary edges, exactly: segments and circular arcs from bulges.
//!
//! An arc edge from `P0` to `P1` with bulge `b` (`b = tan(theta / 4)`) lies
//! on the circle through both endpoints whose centre is
//! `C = M + (1 - b^2) / (4b) * perp(D)`, with `D = P1 - P0`, `M` the chord
//! midpoint and `perp(x, y) = (-y, x)`. Multiplying the circle equation by
//! `(4b)^2` makes every coefficient an exact dyadic:
//!
//! ```text
//! k = 4b,  kC = k*M + (1 - b^2) * perp(D)
//! k^2 |X|^2 - 2k X.(kC) + |kC|^2 - |D|^2 (1 + b^2)^2 = 0
//! ```
//!
//! Points on the arc are exactly those on the circle that lie on the side
//! `-sign(b)` of the chord line, plus the two endpoints. Along the arc,
//! `A` comes before `B` iff `sign(b) * orient(P0, A, B) > 0`: seen from
//! `P0`, the direction to a point sweeping the arc turns monotonically
//! (inscribed angle theorem), whatever the arc's size.
//!
//! Rational points on the arc come from the half-angle parametrization
//! `s = b (1 - u)`, `u in [0, 1]`:
//!
//! ```text
//! w = b (1 + s^2)^2,  g = (b - s)(1 + b s)
//! X(u) = P0 + (g / w) * ((1 - s^2) D - 2 s perp(D))
//! ```
//!
//! `u = 0` is `P0`, `u = 1` is `P1`, and a dyadic `u` gives an exact
//! dyadic homogeneous point.

use axiolid_core::Point2;
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;

use super::point::{dy, orient, same_point, sgn, sign, Circle, Pred, XPoint};

/// Which curve an edge follows.
// Edges live in short per-call vectors; boxing the arc would add an
// allocation per edge to save a few hundred bytes.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone)]
pub(crate) enum Carrier {
    /// A straight segment.
    Segment,
    /// A circular arc; `turn` is the sign of the bulge.
    Arc {
        circle: Circle,
        turn: Sign,
        bulge: Dyadic,
    },
}

/// One input edge.
#[derive(Debug, Clone)]
pub(crate) struct Edge {
    pub(crate) p0: XPoint,
    pub(crate) p1: XPoint,
    p0f: Point2,
    d: (Dyadic, Dyadic),
    pub(crate) carrier: Carrier,
}

fn half(value: &Dyadic) -> Dyadic {
    value.mul(&dy(0.5))
}

impl Edge {
    /// The edge from `from` to `to` with the given bulge.
    pub(crate) fn new(from: Point2, to: Point2, bulge: f64) -> Self {
        let (p0x, p0y) = (dy(from.x), dy(from.y));
        let d = (dy(to.x).sub(&p0x), dy(to.y).sub(&p0y));
        let carrier = if bulge == 0.0 {
            Carrier::Segment
        } else {
            let b = dy(bulge);
            let k = b.mul(&dy(4.0));
            let one_minus = dy(1.0).sub(&b.square());
            let mx = half(&p0x.add(&dy(to.x)));
            let my = half(&p0y.add(&dy(to.y)));
            // perp(D) = (-dy, dx)
            let kcx = k.mul(&mx).sub(&one_minus.mul(&d.1));
            let kcy = k.mul(&my).add(&one_minus.mul(&d.0));
            let chord2 = d.0.square().add(&d.1.square());
            let one_plus = dy(1.0).add(&b.square());
            let circle = Circle {
                alpha: k.square(),
                bx: dy(-2.0).mul(&k).mul(&kcx),
                by: dy(-2.0).mul(&k).mul(&kcy),
                gamma: kcx
                    .square()
                    .add(&kcy.square())
                    .sub(&chord2.mul(&one_plus.square())),
            };
            Carrier::Arc {
                circle,
                turn: sgn(&b),
                bulge: b,
            }
        };
        Self {
            p0: XPoint::from_f64(from),
            p1: XPoint::from_f64(to),
            p0f: from,
            d,
            carrier,
        }
    }

    pub(crate) fn is_arc(&self) -> bool {
        matches!(self.carrier, Carrier::Arc { .. })
    }

    /// The exact rational point at parameter `u` in `[0, 1]`.
    pub(crate) fn point_at(&self, u: &Dyadic) -> XPoint {
        let (p0x, p0y) = (dy(self.p0f.x), dy(self.p0f.y));
        match &self.carrier {
            Carrier::Segment => XPoint::rational(
                p0x.add(&u.mul(&self.d.0)),
                p0y.add(&u.mul(&self.d.1)),
                dy(1.0),
            ),
            Carrier::Arc { bulge: b, .. } => {
                let s = b.mul(&dy(1.0).sub(u));
                let s2 = s.square();
                let w = b.mul(&dy(1.0).add(&s2).square());
                let g = b.sub(&s).mul(&dy(1.0).add(&b.mul(&s)));
                let one_minus = dy(1.0).sub(&s2);
                let two_s = dy(2.0).mul(&s);
                // V = (1 - s^2) D - 2 s perp(D), perp(D) = (-Dy, Dx)
                let vx = one_minus.mul(&self.d.0).add(&two_s.mul(&self.d.1));
                let vy = one_minus.mul(&self.d.1).sub(&two_s.mul(&self.d.0));
                XPoint::rational(
                    p0x.mul(&w).add(&g.mul(&vx)),
                    p0y.mul(&w).add(&g.mul(&vy)),
                    w,
                )
            }
        }
    }

    /// [`Edge::point_at`] in plain `f64`, for seeding searches only.
    pub(crate) fn approx_at(&self, u: f64) -> Point2 {
        let p1 = self.p1.approx();
        let (dx, dy) = (p1.x - self.p0f.x, p1.y - self.p0f.y);
        match &self.carrier {
            Carrier::Segment => Point2::new(self.p0f.x + u * dx, self.p0f.y + u * dy),
            Carrier::Arc { bulge, .. } => {
                let b = bulge.to_f64();
                let s = b * (1.0 - u);
                let w = b * (1.0 + s * s) * (1.0 + s * s);
                let f = (b - s) * (1.0 + b * s) / w;
                let one_minus = 1.0 - s * s;
                Point2::new(
                    self.p0f.x + f * (one_minus * dx + 2.0 * s * dy),
                    self.p0f.y + f * (one_minus * dy - 2.0 * s * dx),
                )
            }
        }
    }

    /// Whether `x`, known to lie on this edge's line or circle, is on the
    /// edge itself (endpoints included).
    pub(crate) fn holds(&self, x: &XPoint) -> bool {
        match &self.carrier {
            Carrier::Segment => {
                sign(Pred::Dot(&self.p0, x, &self.d)) != Sign::Negative
                    && sign(Pred::Dot(x, &self.p1, &self.d)) != Sign::Negative
            }
            Carrier::Arc { turn, .. } => {
                same_point(x, &self.p0)
                    || same_point(x, &self.p1)
                    || orient(&self.p0, &self.p1, x) == turn.flip()
            }
        }
    }

    /// Whether `x` lies on this edge (any point, carrier not assumed).
    pub(crate) fn contains(&self, x: &XPoint) -> bool {
        let on_carrier = match &self.carrier {
            Carrier::Segment => orient(&self.p0, &self.p1, x) == Sign::Zero,
            Carrier::Arc { circle, .. } => sign(Pred::OnCircle(x, circle)) == Sign::Zero,
        };
        on_carrier && self.holds(x)
    }

    /// Sign of `position(a) - position(b)` along the edge, both on it.
    pub(crate) fn order(&self, a: &XPoint, b: &XPoint) -> Sign {
        if same_point(a, b) {
            return Sign::Zero;
        }
        match &self.carrier {
            // position(b) - position(a) has the sign of dot(b - a, D).
            Carrier::Segment => sign(Pred::Dot(a, b, &self.d)).flip(),
            Carrier::Arc { turn, .. } => {
                if same_point(a, &self.p0) || same_point(b, &self.p1) {
                    return Sign::Negative;
                }
                if same_point(b, &self.p0) || same_point(a, &self.p1) {
                    return Sign::Positive;
                }
                // a before b iff turn * orient(P0, a, b) > 0
                let o = orient(&self.p0, a, b);
                if o == *turn {
                    Sign::Negative
                } else {
                    Sign::Positive
                }
            }
        }
    }
}

/// Points where the line `q / qw + t * e` meets `circle`, in increasing
/// `t`: none, one (tangent, `true`), or two.
///
/// Scaled by `qw^2`, the quadratic is `A t^2 + 2 B t + C` with
/// `A = alpha |e|^2 qw^2`, `B = alpha (q.e) qw + (bx ex + by ey) qw^2 / 2`,
/// `C = alpha |q|^2 + (bx qx + by qy) qw + gamma qw^2`, and the roots
/// `t = (-B +- sqrt(B^2 - A C)) / A` give points
/// `x = (qx alpha |e|^2 qw - ex B +- ex sqrt(disc)) / A`.
pub(crate) fn line_meets_circle(
    q: (&Dyadic, &Dyadic, &Dyadic),
    e: &(Dyadic, Dyadic),
    circle: &Circle,
) -> (Vec<XPoint>, bool) {
    let (qx, qy, qw) = q;
    let e2 = e.0.square().add(&e.1.square());
    let qw2 = qw.square();
    let a = circle.alpha.mul(&e2).mul(&qw2);
    let qe = qx.mul(&e.0).add(&qy.mul(&e.1));
    let b = circle
        .alpha
        .mul(&qe)
        .mul(qw)
        .add(&half(&circle.bx.mul(&e.0).add(&circle.by.mul(&e.1))).mul(&qw2));
    let c = circle
        .alpha
        .mul(&qx.square().add(&qy.square()))
        .add(&circle.bx.mul(qx).add(&circle.by.mul(qy)).mul(qw))
        .add(&circle.gamma.mul(&qw2));
    let disc = b.square().sub(&a.mul(&c));
    let base = circle.alpha.mul(&e2).mul(qw);
    let xa = qx.mul(&base).sub(&e.0.mul(&b));
    let ya = qy.mul(&base).sub(&e.1.mul(&b));
    match sgn(&disc) {
        Sign::Negative => (Vec::new(), false),
        Sign::Zero => (vec![XPoint::rational(xa, ya, a)], true),
        _ => {
            let minus = XPoint::new(
                xa.clone(),
                e.0.neg(),
                ya.clone(),
                e.1.neg(),
                a.clone(),
                disc.clone(),
            );
            let plus = XPoint::new(xa, e.0.clone(), ya, e.1.clone(), a, disc);
            (vec![minus, plus], false)
        }
    }
}

/// Every point where two edges meet, endpoints included. Overlapping
/// pieces contribute the endpoints of each edge lying on the other.
pub(crate) fn crossings(first: &Edge, second: &Edge) -> Vec<XPoint> {
    let mut out: Vec<XPoint> = Vec::new();
    let mut push = |x: XPoint| {
        if first.holds(&x) && second.holds(&x) && !out.iter().any(|y| same_point(y, &x)) {
            out.push(x);
        }
    };
    let endpoints_on_each_other = |push: &mut dyn FnMut(XPoint)| {
        for p in [&first.p0, &first.p1] {
            if second.contains(p) {
                push(p.clone());
            }
        }
        for p in [&second.p0, &second.p1] {
            if first.contains(p) {
                push(p.clone());
            }
        }
    };
    match (&first.carrier, &second.carrier) {
        (Carrier::Segment, Carrier::Segment) => {
            // Lines P0 + t D and Q0 + s E: cross(D, E) = 0 means parallel.
            let (d, e) = (&first.d, &second.d);
            let den = d.0.mul(&e.1).sub(&d.1.mul(&e.0));
            if sgn(&den) == Sign::Zero {
                endpoints_on_each_other(&mut push);
            } else {
                // t = cross(Q0 - P0, E) / cross(D, E)
                let (p0x, p0y) = (dy(first.p0f.x), dy(first.p0f.y));
                let (q0x, q0y) = (dy(second.p0f.x), dy(second.p0f.y));
                let (rx, ry) = (q0x.sub(&p0x), q0y.sub(&p0y));
                let num = rx.mul(&e.1).sub(&ry.mul(&e.0));
                push(XPoint::rational(
                    p0x.mul(&den).add(&num.mul(&d.0)),
                    p0y.mul(&den).add(&num.mul(&d.1)),
                    den,
                ));
            }
        }
        (Carrier::Segment, Carrier::Arc { circle, .. })
        | (Carrier::Arc { circle, .. }, Carrier::Segment) => {
            let seg = if first.is_arc() { second } else { first };
            let (px, py) = (dy(seg.p0f.x), dy(seg.p0f.y));
            let (points, _) = line_meets_circle((&px, &py, &dy(1.0)), &seg.d, circle);
            for x in points {
                push(x);
            }
            // Endpoints landing exactly on the other edge are found by the
            // quadratic too; this keeps the rule uniform with the others.
            endpoints_on_each_other(&mut push);
        }
        (Carrier::Arc { circle: c1, .. }, Carrier::Arc { circle: c2, .. }) => {
            if c1.same_as(c2) {
                endpoints_on_each_other(&mut push);
            } else {
                // Radical line: alpha2 F1 - alpha1 F2 = u x + v y + g = 0.
                let u = c2.alpha.mul(&c1.bx).sub(&c1.alpha.mul(&c2.bx));
                let v = c2.alpha.mul(&c1.by).sub(&c1.alpha.mul(&c2.by));
                let g = c2.alpha.mul(&c1.gamma).sub(&c1.alpha.mul(&c2.gamma));
                if sgn(&u) != Sign::Zero || sgn(&v) != Sign::Zero {
                    // Q = (-u g, -v g) / n, n = u^2 + v^2; direction (-v, u).
                    let n = u.square().add(&v.square());
                    let (qx, qy) = (u.mul(&g).neg(), v.mul(&g).neg());
                    let (points, _) = line_meets_circle((&qx, &qy, &n), &(v.neg(), u), c1);
                    for x in points {
                        push(x);
                    }
                }
                // Concentric distinct circles never meet.
                endpoints_on_each_other(&mut push);
            }
        }
    }
    out
}
