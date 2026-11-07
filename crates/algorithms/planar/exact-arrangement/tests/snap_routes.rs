//! The route of an edge through the hot cells, against rounding the points of
//! the edge one by one (ADR 0100).
//!
//! The reference shares no code with `snap_round`: it cuts the edge at every
//! line halfway between neighbouring binary64 numbers of every hot centre,
//! rounds a sample of each stretch and at each cut with `ExactPoints::round`,
//! and reads the cells off in order. Ends sit on and beside cell corners and
//! ties (a quarter, half and whole step off a binary64 number), with exact
//! weights 1, 3, 5 and 7, near zero, subnormals, powers of two and the
//! largest numbers. The cases also pin an edge whose far end sits on the
//! corner of a hot cell it does not enter.
use axiolid_contracts::{ExecutionOptions, Sign};
use axiolid_core::{Point2, Tolerance};
use axiolid_exact::{Arith, Dyadic};
use axiolid_exact_arrangement::{snap_round, ExactPoints, PairBudget, PointId, Segment};
use num_bigint::BigInt;
use std::cmp::Ordering;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}
fn d(v: f64) -> Dyadic {
    Dyadic::try_from_f64(v).unwrap()
}
fn sgn(x: &Dyadic) -> Sign {
    x.sign().unwrap()
}

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 11
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn step(mut v: f64, j: i64) -> f64 {
    if j > 0 {
        for _ in 0..j {
            v = v.next_up();
        }
    } else {
        for _ in 0..(-j) {
            v = v.next_down();
        }
    }
    v
}

/// exact coordinate near grid value g(j) of base: (g*w + m*q)/w scaled into homogeneous with weight w
fn coord(base: f64, j: i64, m: i64, w: i64) -> Option<Dyadic> {
    let g = step(base, j);
    let n = step(base, j + 1);
    if !g.is_finite() || !n.is_finite() {
        return None;
    }
    let q = d(n).sub(&d(g)).mul(&d(0.25));
    Some(d(g).mul(&d(w as f64)).add(&q.mul(&d(m as f64))))
}

fn rnd(h: [Dyadic; 3]) -> Option<Point2> {
    let mut s = ExactPoints::new();
    let id = s.push(h).ok()?;
    let p = s.round(id).ok()?;
    Some(Point2::new(p.x + 0.0, p.y + 0.0))
}

#[derive(Clone)]
struct T {
    n: Dyadic,
    d: Dyadic,
}
fn cmp_t(a: &T, b: &T) -> Ordering {
    match sgn(&a.n.mul(&b.d).sub(&b.n.mul(&a.d))) {
        Sign::Positive => Ordering::Greater,
        Sign::Negative => Ordering::Less,
        _ => Ordering::Equal,
    }
}
fn at(from: &[Dyadic; 3], to: &[Dyadic; 3], t: &T) -> [Dyadic; 3] {
    let r = t.d.sub(&t.n);
    [0, 1, 2].map(|i| r.mul(&from[i]).add(&t.n.mul(&to[i])))
}
fn lines(c: f64) -> Vec<Dyadic> {
    let (p, n) = (c.next_down(), c.next_up());
    let mut v = Vec::new();
    if p.is_finite() {
        v.push(d(c).add(&d(p)).mul(&d(0.5)));
    }
    if n.is_finite() {
        v.push(d(c).add(&d(n)).mul(&d(0.5)));
    }
    v
}

fn reference(from: &[Dyadic; 3], to: &[Dyadic; 3], hot: &[Point2]) -> Vec<Point2> {
    let mut ts = vec![
        T {
            n: Dyadic::zero(),
            d: d(1.0),
        },
        T {
            n: d(1.0),
            d: d(1.0),
        },
    ];
    for c in hot {
        for axis in 0..2 {
            let cc = if axis == 0 { c.x } else { c.y };
            for b in lines(cc) {
                let h0 = from[axis].sub(&b.mul(&from[2]));
                let h1 = to[axis].sub(&b.mul(&to[2]));
                let mut den = h0.sub(&h1);
                let mut num = h0.clone();
                match sgn(&den) {
                    Sign::Zero => continue,
                    Sign::Negative => {
                        den = den.neg();
                        num = num.neg();
                    }
                    _ => {}
                }
                if sgn(&num) == Sign::Negative || sgn(&den.sub(&num)) == Sign::Negative {
                    continue;
                }
                ts.push(T { n: num, d: den });
            }
        }
    }
    ts.sort_by(cmp_t);
    ts.dedup_by(|a, b| cmp_t(a, b).is_eq());
    let mut samples = Vec::new();
    for i in 0..ts.len() {
        samples.push(ts[i].clone());
        if i + 1 < ts.len() {
            let (a, b) = (&ts[i], &ts[i + 1]);
            samples.push(T {
                n: a.n.mul(&b.d).add(&b.n.mul(&a.d)),
                d: a.d.mul(&b.d).mul(&d(2.0)),
            });
        }
    }
    let mut seq: Vec<Point2> = Vec::new();
    for t in &samples {
        let p = rnd(at(from, to, t)).unwrap();
        if hot.contains(&p) && seq.last() != Some(&p) {
            assert!(!seq.contains(&p), "reference: a convex cell re-entered");
            seq.push(p);
        }
    }
    seq
}

fn on_open(a: Point2, b: Point2, p: Point2) -> bool {
    let (ax, ay, bx, by, px, py) = (d(a.x), d(a.y), d(b.x), d(b.y), d(p.x), d(p.y));
    let (ux, uy, vx, vy) = (bx.sub(&ax), by.sub(&ay), px.sub(&ax), py.sub(&ay));
    let cross = ux.mul(&vy).sub(&uy.mul(&vx));
    let dot = ux.mul(&vx).add(&uy.mul(&vy));
    let len = ux.mul(&ux).add(&uy.mul(&uy));
    sgn(&cross) == Sign::Zero
        && sgn(&dot) == Sign::Positive
        && sgn(&len.sub(&dot)) == Sign::Positive
}

#[test]
fn fuzz_routes_against_rounding() {
    let bases = [
        0.0,
        1.0,
        2.0,
        0.5,
        -1.0,
        -2.0,
        f64::MIN_POSITIVE,
        -f64::MIN_POSITIVE,
        f64::from_bits(3),
        f64::MAX
            .next_down()
            .next_down()
            .next_down()
            .next_down()
            .next_down()
            .next_down()
            .next_down()
            .next_down(),
        3.0,
        1e-300,
    ];
    let mut rng = Lcg(0xfeed_beef);
    let (mut ok, mut bad, mut errs, mut rer, mut skipped) = (0, 0, 0, 0, 0);
    let mut reported = 0;
    for case in 0..2000u32 {
        let bx = bases[rng.below(bases.len() as u64) as usize];
        let by = bases[rng.below(bases.len() as u64) as usize];
        let mut pts = ExactPoints::new();
        let endp = |rng: &mut Lcg, pts: &mut ExactPoints| -> Option<PointId> {
            let w = [1i64, 1, 1, 3, 5, 7][rng.below(6) as usize];
            let m = if w == 1 {
                [0i64, 2, -2, 1, -1, 2][rng.below(6) as usize]
            } else {
                rng.below(9) as i64 - 4
            };
            let m2 = if w == 1 {
                [0i64, 2, -2, 1, -1, 2][rng.below(6) as usize]
            } else {
                rng.below(9) as i64 - 4
            };
            let x = coord(bx, rng.below(13) as i64 - 6, m, w)?;
            let y = coord(by, rng.below(13) as i64 - 6, m2, w)?;
            pts.push([x, y, d(w as f64)]).ok()
        };
        let (Some(a), Some(b)) = (endp(&mut rng, &mut pts), endp(&mut rng, &mut pts)) else {
            skipped += 1;
            continue;
        };
        let (Ok(ra), Ok(rb)) = (pts.round(a), pts.round(b)) else {
            skipped += 1;
            continue;
        };
        if pts.exact(a) == pts.exact(b) {
            skipped += 1;
            continue;
        }
        let _ = (ra, rb);
        let mut edges = vec![Segment { ends: [a, b] }];
        // dots: hot cells without fragments
        let ndots = rng.below(12);
        for _ in 0..ndots {
            let gx = step(bx, rng.below(13) as i64 - 6);
            let gy = step(by, rng.below(13) as i64 - 6);
            if !gx.is_finite() || !gy.is_finite() {
                continue;
            }
            let p = pts.push_flat(Point2::new(gx, gy)).unwrap();
            if p == a || p == b {
                continue;
            }
            // second end: a 3-weighted point off p by 2^-1100 in x: rounds to p
            let tiny = Dyadic::from_parts(BigInt::from(1), -1100);
            let q = pts
                .push([d(gx).mul(&d(3.0)).add(&tiny), d(gy).mul(&d(3.0)), d(3.0)])
                .unwrap();
            if pts.round(q).unwrap() != Point2::new(gx, gy) {
                continue;
            }
            edges.push(Segment { ends: [p, q] });
        }
        let mut hot: Vec<Point2> = Vec::new();
        for e in &edges {
            for &id in &e.ends {
                let p = pts.round(id).unwrap();
                let p = Point2::new(p.x + 0.0, p.y + 0.0);
                if !hot.contains(&p) {
                    hot.push(p);
                }
            }
        }
        let r = reference(pts.exact(a), pts.exact(b), &hot);
        let s = match snap_round(
            &pts,
            &edges,
            &mut PairBudget::new(10_000_000, "t"),
            0,
            &options(),
        ) {
            Ok(s) => s,
            Err(e) => {
                errs += 1;
                if reported < 6 {
                    reported += 1;
                    eprintln!("case {case}: error {e:?}; reference {r:?}");
                }
                continue;
            }
        };
        let route = &s.routes[0];
        let chain: Vec<Point2> = if route.is_empty() {
            vec![]
        } else {
            let piece = |p: usize| s.arrangement.pieces[p].ends;
            let mut c = Vec::new();
            let (p0, f0) = route[0];
            let e0 = piece(p0.0);
            c.push(s.point(if f0 { e0[0] } else { e0[1] }).unwrap());
            for &(p, f) in route {
                let e = piece(p.0);
                c.push(s.point(if f { e[1] } else { e[0] }).unwrap());
            }
            c
        };
        let good = if chain.is_empty() {
            r.len() == 1
        } else {
            let mut j = 0usize;
            let mut fine = chain[0] == r[0];
            for &p in &chain[1..] {
                if !fine {
                    break;
                }
                if j + 1 < r.len() && p == r[j + 1] {
                    j += 1;
                } else if j + 1 < r.len() && on_open(r[j], r[j + 1], p) {
                } else {
                    fine = false;
                }
            }
            fine && j + 1 == r.len()
        };
        let interior = r.len().saturating_sub(2);
        assert_eq!(
            s.rerouted,
            usize::from(interior > 0),
            "case {case}: rerouted against {interior} interior cells"
        );
        rer += usize::from(interior > 0);
        if good {
            ok += 1
        } else {
            bad += 1;
            if reported < 12 {
                reported += 1;
                eprintln!(
                    "case {case}: route {chain:?}\n  ref {r:?}\n  a={:?} b={:?}",
                    pts.exact(a),
                    pts.exact(b)
                );
            }
        }
    }
    eprintln!("ok={ok} bad={bad} errs={errs} rerouted_cases={rer} skipped={skipped}");
    assert_eq!(bad + errs, 0);
}

fn flat(points: &mut ExactPoints, x: f64, y: f64) -> PointId {
    points.push_flat(Point2::new(x, y)).unwrap()
}

/// `A` is binary64 and `Q = (1 + e/2, 1 + e/2)` is not; `Q` lies on the corner
/// of the cells of `(1, 1 + e)` and `(1 + e, 1)`, which the corner does not
/// belong to (it is a tie that the even centre `(1, 1)` owns). An open bound
/// reached exactly at the far end of the edge must exclude the cell.
#[test]
fn an_end_on_a_cell_corner_does_not_enter_the_cell() {
    let e = f64::EPSILON;
    let mut points = ExactPoints::new();
    let a = flat(&mut points, 1.0 + 4.0 * e, 1.0 + 4.0 * e);
    let q = points.push([d(2.0 + e), d(2.0 + e), d(2.0)]).unwrap();
    let left = (
        flat(&mut points, 1.0, 1.0 + e),
        flat(&mut points, 0.5, 1.0 + e),
    );
    let below = (
        flat(&mut points, 1.0 + e, 1.0),
        flat(&mut points, 1.0 + e, 0.5),
    );
    let mut edges = vec![
        Segment { ends: [a, q] },
        Segment {
            ends: [left.0, left.1],
        },
    ];
    let run = |edges: &[Segment]| {
        snap_round(
            &points,
            edges,
            &mut PairBudget::new(1000, "t"),
            0,
            &options(),
        )
    };
    let one = run(&edges).expect("one side cell is hot");
    assert_eq!(
        one.rerouted, 0,
        "the corner is not in the cell of (1, 1 + e)"
    );
    assert_eq!(one.routes[0].len(), 1);
    edges.push(Segment {
        ends: [below.0, below.1],
    });
    let both = run(&edges).expect("both side cells are hot");
    assert_eq!(both.rerouted, 0);
    assert_eq!(both.routes[0].len(), 1);
}
