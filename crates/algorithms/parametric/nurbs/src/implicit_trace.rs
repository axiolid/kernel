//! Every component of a field's zero set in a parameter box, as
//! [`ImplicitCurve2`]s (ADR 0077).
//!
//! The box is subdivided until each piece either certainly misses the
//! curve (the field's bound excludes zero) or is *regular*: one partial
//! derivative's bound excludes zero, so along that parameter the field is
//! strictly monotone and the curve crosses each line across the piece at
//! most once. A regular piece holds the curve as a graph over the other
//! parameter, existing exactly where the field has opposite signs on the
//! piece's two sides -- decided by isolating the field's roots along those
//! sides, each root certified simple. Those stretches are the curve's
//! cells. Cells of neighbouring pieces meet where the curve crosses their
//! common side, so chaining them by their end points gives the components,
//! closed loops or open chains leaving the box. Nothing is marched and no
//! step size is guessed: a component is found because some regular piece
//! must contain part of it, and pieces cover the box.
//!
//! A piece where neither partial can be bounded away from zero while the
//! field may vanish, down to the smallest size allowed, may hold a singular
//! point of the curve (where the surfaces touch): the trace is refused
//! there by name rather than guessed through.

use axiolid_core::{Point2, Scalar, Vec2};
use axiolid_curve::{Axis, Field2, ImplicitCell, ImplicitCurve2};
use core::f64::consts::TAU;

use axiolid_curve::implicit::{bound, bound_simple, partial, Cell, Range};
use axiolid_guarantees::Sign;

use crate::exact_field::Exact;

/// Why a trace was refused.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum TraceRefusal {
    /// The field and both partials may vanish together near this point: a
    /// singular point of the curve (the surfaces touch), or two branches
    /// closer than the trace resolves.
    Singular(Point2),
    /// The subdivision exceeded its work budget.
    Budget,
    /// The zero set is only isolated points: the surfaces touch there
    /// without crossing.
    Touching(Point2),
    /// More singular points than isolated ones come to (a whole curve of
    /// them, where the surfaces may be tangent along it).
    Many(Point2),
}

/// How many singular points a first trace takes before it hands over to
/// [`along_contact`].
const FEW: usize = 32;

/// Periodicity of the box's parameters: a component leaving through one
/// side of a periodic parameter comes back through the other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Periodic {
    pub(crate) u: bool,
    pub(crate) v: bool,
}

/// Every component of `field = 0` in `domain`.
pub(crate) fn trace(
    field: &Field2,
    domain: Cell,
    periodic: Periodic,
) -> Result<Vec<ImplicitCurve2>, TraceRefusal> {
    trace_within(field, domain, periodic, BUDGET)
}

/// The work budget of a trace, in pieces examined.
pub(crate) const BUDGET: usize = 400_000;

/// [`trace`] within a work budget of its own.
pub(crate) fn trace_within(
    field: &Field2,
    domain: Cell,
    periodic: Periodic,
    budget: usize,
) -> Result<Vec<ImplicitCurve2>, TraceRefusal> {
    let (curves, touches) =
        match trace_with_touches(field, domain, periodic, budget, &Tubes::none(), FEW) {
            Ok(found) => found,
            // Undecidable pieces all along a curve: the surfaces may be
            // tangent along it.
            Err(TraceRefusal::Many(_)) => match along_contact(field, domain, periodic, budget) {
                Some(found) => found?,
                // Not a line of contact after all (a pole, a curve of
                // touching points on the window's edge): the singular
                // points one by one.
                None => {
                    trace_with_touches(field, domain, periodic, budget, &Tubes::none(), usize::MAX)?
                }
            },
            Err(refusal @ (TraceRefusal::Singular(_) | TraceRefusal::Budget)) => {
                along_contact(field, domain, periodic, budget).ok_or(refusal)??
            }
            Err(refusal) => return Err(refusal),
        };
    // Only isolated touching points: the surfaces touch, they do not cross.
    if curves.is_empty() && !touches.is_empty() {
        return Err(TraceRefusal::Touching(touches[0]));
    }
    Ok(curves)
}

/// Thin tubes about curves along which the surfaces are tangent: pieces
/// inside one are left out of a trace, since the field vanishes there to
/// higher order and no piece can be certified. Sample points of the curves,
/// bucketed by the tube's radius.
#[derive(Debug, Clone, Default)]
pub(crate) struct Tubes {
    radius: Scalar,
    buckets: std::collections::HashMap<(i64, i64), Vec<Point2>>,
    /// Set by a trace that met a degenerate singular point (its branches
    /// there bridged by the field's signs, not by certified directions).
    degenerate: core::cell::Cell<bool>,
}

impl Tubes {
    fn none() -> Self {
        Self::default()
    }

    fn new(radius: Scalar, points: impl Iterator<Item = Point2>) -> Self {
        let mut buckets: std::collections::HashMap<(i64, i64), Vec<Point2>> =
            std::collections::HashMap::new();
        for p in points {
            let key = ((p.x / radius).floor() as i64, (p.y / radius).floor() as i64);
            buckets.entry(key).or_default().push(p);
        }
        Self {
            radius,
            buckets,
            degenerate: core::cell::Cell::new(false),
        }
    }

    /// Whether a sample lies within the radius of the box.
    fn near(&self, lo: Point2, hi: Point2) -> bool {
        if self.buckets.is_empty() {
            return false;
        }
        let r = self.radius;
        let (a, b) = (lo - Point2::splat(r), hi + Point2::splat(r));
        let key = |x: Scalar| (x / r).floor() as i64;
        if key(b.x) - key(a.x) > 8 || key(b.y) - key(a.y) > 8 {
            return false;
        }
        for i in key(a.x)..=key(b.x) {
            for j in key(a.y)..=key(b.y) {
                if let Some(points) = self.buckets.get(&(i, j)) {
                    if points
                        .iter()
                        .any(|p| p.x >= a.x && p.x <= b.x && p.y >= a.y && p.y <= b.y)
                    {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// Whether a piece lies in a tube: no wider than the radius, and near
    /// a sample.
    fn holds(&self, cell: &Cell) -> bool {
        let side = (cell.hi.x - cell.lo.x).max(cell.hi.y - cell.lo.y);
        side <= self.radius && self.near(cell.lo, cell.hi)
    }
}

/// How a curve of a second derivative's zeros sits in the field's own.
enum Contact {
    /// On it, the field changing sign across: the surfaces are tangent
    /// along the curve and cross there.
    Crossing,
    /// On it, the field keeping its sign: tangent without crossing.
    Touching,
    /// Not on it throughout.
    Not,
}

/// The field's zeros where the surfaces are tangent along whole curves.
///
/// Along such a curve the field vanishes with its gradient (and, where the
/// surfaces cross, its Hessian), as `c n^2` or `c n^3` in the distance `n`
/// across it, so no piece near it can be certified. But then every first
/// (or second) derivative vanishes exactly on the curve, and generically
/// one of them crosses zero there regularly (`F_uu ~ 6 c a^2 n`): its
/// trace, certified as usual, holds the curve.
/// Each of its curves is kept as a contact curve only where the field
/// vanishes along the whole of it (to rounding) -- as crossing where the
/// field changes sign across it, touching where not. The field is then
/// traced again with thin tubes about the contact curves left out; a
/// branch of it ending at a tube runs into a contact curve, and is bridged
/// to it there, the contact curve cut: a vertex, as at a crossing.
#[allow(clippy::type_complexity)]
fn along_contact(
    field: &Field2,
    domain: Cell,
    periodic: Periodic,
    budget: usize,
) -> Option<Result<(Vec<ImplicitCurve2>, Vec<Point2>), TraceRefusal>> {
    let extent = (domain.hi - domain.lo).abs();
    let size = extent.x.max(extent.y);
    // Where the field vanishes as `c n^m` across the curve, every
    // derivative of order `m - 1` vanishes exactly on it, and generically
    // one of them crosses zero there regularly. Orders are tried from the
    // first up: `level[i]` is the derivative taken `i` times along `u` and
    // the rest along `v`.
    let mut level = vec![field.clone()];
    // A B-spline field's derivatives run out at its degree, and its signs
    // are exact: every order is tried. A series field's never run out.
    let orders = match field {
        Field2::Patches(f) => f.u_degree + f.v_degree,
        Field2::Series(_) => MAX_ORDER,
    };
    for order in 1..=orders {
        let next: Vec<Field2> = (0..=order)
            .map(|i| {
                if i == 0 {
                    partial(&level[0], false)
                } else {
                    partial(&level[i - 1], true)
                }
            })
            .collect();
        if next.iter().all(|g| g.magnitude() == 0.0) {
            return None;
        }
        level = next;
        if let Some(found) =
            contact_of_order(field, &level, order + 1, domain, periodic, budget, size)
        {
            return found;
        }
    }
    None
}

/// The highest order of derivative [`along_contact`] tries on a series
/// field (a B-spline field's run out at its degree): contact up to
/// `c n^(MAX_ORDER + 1)`.
const MAX_ORDER: usize = 16;

/// [`along_contact`] with the derivatives of one order, `level`, for
/// contact of order `m`. `None` where none of them yields a line.
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn contact_of_order(
    field: &Field2,
    level: &[Field2],
    m: usize,
    domain: Cell,
    periodic: Periodic,
    budget: usize,
    size: Scalar,
) -> Option<Option<Result<(Vec<ImplicitCurve2>, Vec<Point2>), TraceRefusal>>> {
    for g in level {
        if g.magnitude() == 0.0 {
            continue;
        }
        // Only a derivative that crosses zero regularly holds a line: one
        // whose own trace met a degenerate point there is itself tangent
        // to the line, and its bridged branches are no better than the
        // field's.
        let watch = Tubes::none();
        let Ok((candidates, _)) = trace_with_touches(g, domain, periodic, budget, &watch, FEW)
        else {
            continue;
        };
        if watch.degenerate.get() {
            continue;
        }
        let (mut crossing, mut touching) = (Vec::new(), Vec::new());
        for curve in candidates {
            match contact(field, g, &curve, size) {
                Contact::Crossing => crossing.push(curve),
                Contact::Touching => touching.push(curve),
                Contact::Not => {}
            }
        }
        if crossing.is_empty() && touching.is_empty() {
            continue;
        }
        // The tubes: out to where the field, `c n^m`, is 64 times its
        // rounding, with `c` read off the derivative's slope across the
        // line (`m! c`); at least a thousandth of the window. Wider than a
        // twentieth, the lines are not followed: rounding hides them.
        let mut radius = 1e-3 * size;
        let factorial: Scalar = (1..=m).map(|k| k as Scalar).product();
        // With a certified tier (every field of finite coefficients has
        // one, B-spline or series) nothing is hidden: the least width
        // serves.
        let exact = Exact::of(field).is_some();
        for curve in crossing.iter().chain(&touching).filter(|_| !exact) {
            let p = curve.point(0.5 * curve.end())?;
            let slope = g.jet(p).gradient.length();
            let rounding = 64.0 * Scalar::EPSILON * field.scale_at(p).max(field.magnitude());
            if slope > 0.0 && slope.is_finite() {
                // Where `c n^m` reaches 64 times the rounding: decidable.
                let clear = (64.0 * rounding * factorial / slope).powf(1.0 / m as Scalar);
                radius = radius.max(clear);
            }
        }
        if radius > 0.05 * size {
            return Some(None);
        }
        // Samples of every contact curve, a quarter of the radius apart,
        // wrapped into the window along periodic parameters.
        let wrap = |mut p: Point2| {
            if periodic.u {
                p.x = domain.lo.x + (p.x - domain.lo.x).rem_euclid(TAU);
            }
            if periodic.v {
                p.y = domain.lo.y + (p.y - domain.lo.y).rem_euclid(TAU);
            }
            p
        };
        let mut samples = Vec::new();
        for curve in crossing.iter().chain(&touching) {
            for i in 0..curve.cells.len() {
                let (a, b) = (curve.point(i as Scalar)?, curve.point(i as Scalar + 1.0)?);
                let n = (((b - a).length() / (0.25 * radius)).ceil() as usize).clamp(1, 100_000);
                for k in 0..=n {
                    samples.push(wrap(curve.point(i as Scalar + k as Scalar / n as Scalar)?));
                }
            }
        }
        let tubes = Tubes::new(radius, samples.iter().copied());
        let (curves, mut touches) =
            match trace_with_touches(field, domain, periodic, budget, &tubes, usize::MAX) {
                Ok(found) => found,
                Err(_) => return Some(None),
            };
        // A branch running into a contact curve meets it at a vertex: the
        // branch is bridged to the curve's nearest point, and the curve cut
        // there.
        let mut curves = curves;
        let mut cuts: Vec<Vec<(Scalar, Point2)>> =
            vec![Vec::new(); crossing.len() + touching.len()];
        for c in &mut curves {
            for at_start in [true, false] {
                let t = if at_start { 0.0 } else { c.end() };
                let p = c.point(t)?;
                // An end on the window's own edge is where the window cuts
                // the branch, not a tube.
                let edge = 1e-9 * size;
                let on_edge = (!periodic.u
                    && ((p.x - domain.lo.x).abs() <= edge || (p.x - domain.hi.x).abs() <= edge))
                    || (!periodic.v
                        && ((p.y - domain.lo.y).abs() <= edge
                            || (p.y - domain.hi.y).abs() <= edge));
                if on_edge {
                    continue;
                }
                // Ends at a tube: pieces up to its radius were left out
                // about the curve, so the branch stops within a few radii.
                let (k, mut s, mut q) =
                    nearest(crossing.iter().chain(&touching), p, periodic, size)?;
                if (q - p).length() > 4.0 * radius {
                    continue;
                }
                // Where another branch already meets the line within
                // rounding, this one meets it there too: one vertex.
                if let Some(&(t, at)) = cuts[k]
                    .iter()
                    .find(|(_, at)| (*at - q).length() <= 1e-6 * size)
                {
                    (s, q) = (t, at);
                }
                let g = field.jet(p).gradient;
                let leaving = Vec2::new(-g.y, g.x);
                let chord = q - p;
                let cell = ImplicitCell::bridge(p, q, leaving, Vec2::new(chord.x, chord.y));
                if at_start {
                    c.cells.insert(0, cell.reversed());
                } else {
                    c.cells.push(cell);
                }
                cuts[k].push((s, q));
            }
        }
        // Touching curves cut by a branch are not sections themselves.
        let mut all = curves;
        for (k, curve) in crossing.into_iter().enumerate() {
            let mut at: Vec<Scalar> = cuts[k].iter().map(|(t, _)| *t).collect();
            all.extend(cut_at(curve, &mut at, periodic)?);
        }
        if all.is_empty() {
            if let Some(p) = touching.first().and_then(|c| c.point(0.0)) {
                touches.push(p);
            }
        }
        return Some(Some(Ok((all, touches))));
    }
    None
}

/// The curve, its parameter and its point nearest `p` (up to whole turns
/// along periodic parameters): the nearest sample, refined by golden
/// section on the distance.
fn nearest<'a>(
    curves: impl Iterator<Item = &'a ImplicitCurve2>,
    p: Point2,
    periodic: Periodic,
    size: Scalar,
) -> Option<(usize, Scalar, Point2)> {
    let turn = |q: Point2| {
        let mut d = q - p;
        if periodic.u {
            d.x -= (d.x / TAU).round() * TAU;
        }
        if periodic.v {
            d.y -= (d.y / TAU).round() * TAU;
        }
        d.length()
    };
    let mut best: Option<(Scalar, usize, Scalar)> = None;
    let curves: Vec<&ImplicitCurve2> = curves.collect();
    for (k, c) in curves.iter().enumerate() {
        let n = 64 * c.cells.len().max(1);
        for i in 0..=n {
            let t = c.end() * i as Scalar / n as Scalar;
            let Some(q) = c.point(t) else { continue };
            let d = turn(q);
            if best.is_none_or(|(b, _, _)| d < b) {
                best = Some((d, k, t));
            }
        }
    }
    let (_, k, t0) = best?;
    let c = curves[k];
    let step = c.end() / (64 * c.cells.len().max(1)) as Scalar;
    let (mut a, mut b) = ((t0 - step).max(0.0), (t0 + step).min(c.end()));
    let f = |t: Scalar| c.point(t).map_or(Scalar::INFINITY, turn);
    for _ in 0..100 {
        let (m1, m2) = (a + 0.382 * (b - a), b - 0.382 * (b - a));
        if f(m1) <= f(m2) {
            b = m2;
        } else {
            a = m1;
        }
    }
    let mut t = 0.5 * (a + b);
    // An end of the curve (a vertex it already has) within rounding of the
    // search is the point.
    for end in [0.0, c.end()] {
        if let (Some(e), Some(q)) = (c.point(end), c.point(t)) {
            if (e - q).length() <= 1e-6 * size {
                t = end;
            }
        }
    }
    let q = c.point(t)?;
    // Onto the turn nearest `p`.
    let mut shift = Point2::ZERO;
    if periodic.u {
        shift.x = ((p.x - q.x) / TAU).round() * TAU;
    }
    if periodic.v {
        shift.y = ((p.y - q.y) / TAU).round() * TAU;
    }
    Some((k, t, q + shift))
}

/// A contact curve cut at the parameters where branches meet it.
fn cut_at(
    curve: ImplicitCurve2,
    cuts: &mut [Scalar],
    periodic: Periodic,
) -> Option<Vec<ImplicitCurve2>> {
    let n = curve.end();
    let slack = 1e-9 * (1.0 + n);
    cuts.sort_by(Scalar::total_cmp);
    let inner: Vec<Scalar> = cuts
        .iter()
        .copied()
        .filter(|&t| t > slack && t < n - slack)
        .collect();
    if inner.is_empty() {
        return Some(vec![curve]);
    }
    let mut out = Vec::new();
    match curve.closure(periodic.u, periodic.v) {
        Some(closure) => {
            for w in inner.windows(2) {
                out.push(curve.sub(w[0], w[1], None)?);
            }
            out.push(curve.sub(inner[inner.len() - 1], inner[0], Some(closure))?);
        }
        None => {
            let mut ends = vec![0.0];
            ends.extend(inner);
            ends.push(n);
            for w in ends.windows(2) {
                if w[1] - w[0] > slack {
                    out.push(curve.sub(w[0], w[1], None)?);
                }
            }
        }
    }
    Some(out)
}

/// Whether the zeros of a second derivative `g` along `curve` are the
/// field's own, and whether the field changes sign across them.
fn contact(field: &Field2, g: &Field2, curve: &ImplicitCurve2, size: Scalar) -> Contact {
    let exact = Exact::of(field);
    let n = 16 * curve.cells.len().max(1);
    let (mut crossing, mut touching) = (0, 0);
    for k in 0..=n {
        let Some(p) = curve.point(curve.end() * k as Scalar / n as Scalar) else {
            return Contact::Not;
        };
        let rounding = 64.0 * Scalar::EPSILON * field.scale_at(p).max(field.magnitude());
        // On the field's zeros, to rounding (a thousandfold for the cube's
        // cancellation).
        if field.value(p).abs() > 1e3 * rounding {
            return Contact::Not;
        }
        // Across the curve, along `g`'s gradient, far enough for the field
        // to show its sign.
        // Where `g` is itself singular (a vertex of its own, on a bridge)
        // there is no direction across: the field's zero is all to check.
        let t = curve.end() * k as Scalar / n as Scalar;
        let on_bridge = curve
            .cells
            .get((t.floor() as usize).min(curve.cells.len().saturating_sub(1)))
            .is_some_and(|c| c.bridge.is_some());
        let grad = g.jet(p).gradient;
        let l = grad.length();
        let scale = g.jet(p).uu.abs() + g.jet(p).vv.abs() + g.magnitude() * Scalar::EPSILON;
        if on_bridge || !l.is_finite() || l <= 1e-6 * scale {
            continue;
        }
        let normal = grad / l;
        let mut decided = false;
        // Far enough for the field, vanishing as `c n^m`, to show its sign:
        // in `f64` once clear of the rounding, exactly at once where the
        // field has an exact tier.
        let sign = |q: Point2| -> Option<Sign> {
            let f = field.value(q);
            if f.abs() > 1e3 * rounding {
                return Some(if f < 0.0 {
                    Sign::Negative
                } else {
                    Sign::Positive
                });
            }
            exact.and_then(|e| e.sign_at(q, (0, 0)))
        };
        for step in [1e-4, 1e-3, 1e-2, 3e-2, 1e-1] {
            let d = normal * (step * size);
            let (qa, qb) = (p + Point2::new(d.x, d.y), p - Point2::new(d.x, d.y));
            let (Some(a), Some(b)) = (sign(qa), sign(qb)) else {
                continue;
            };
            if a == Sign::Zero || b == Sign::Zero {
                continue;
            }
            // Tangent there: the field's gradient vanishes on the curve
            // (as `n^(m-1)`), far below its size just beside it; on a
            // regular zero it is much the same. Below its own rounding the
            // gradient tells nothing, and the check is passed over.
            let beside = field
                .jet(qa)
                .gradient
                .length()
                .max(field.jet(qb).gradient.length());
            let floor = 64.0 * Scalar::EPSILON * field.magnitude() * size.recip();
            if beside > 1e3 * floor && field.jet(p).gradient.length() > 0.25 * beside {
                return Contact::Not;
            }
            if a != b {
                crossing += 1;
            } else {
                touching += 1;
            }
            decided = true;
            break;
        }
        if !decided {
            return Contact::Not;
        }
    }
    match (crossing, touching) {
        (_, 0) => Contact::Crossing,
        (0, _) => Contact::Touching,
        _ => Contact::Not,
    }
}

/// [`trace`], also returning the isolated points where the zero set is a
/// single point (the surfaces touch there without crossing).
pub(crate) fn trace_with_touches(
    field: &Field2,
    domain: Cell,
    periodic: Periodic,
    budget: usize,
    tubes: &Tubes,
    cap: usize,
) -> Result<(Vec<ImplicitCurve2>, Vec<Point2>), TraceRefusal> {
    let du = partial(field, true);
    let dv = partial(field, false);
    // Signs `f64` cannot decide, decided exactly (B-spline fields).
    let exact = Exact::of(field);
    let extent = (domain.hi - domain.lo).abs();
    let smallest = extent * 1e-9;
    // How far from a singular point a piece may still be its own: a bridge
    // into a crossing is at most this long.
    let limit = 2e-3 * extent.x.max(extent.y);
    // Start from a grid so that the bounds are tight from the outset; its
    // lines sit off the round fractions a symmetric input would put a
    // crossing on.
    let mut queue = split_grid(domain, 8);
    let mut all: Vec<ImplicitCell> = Vec::new();
    // Singular points of the curve found so far, each with the square about
    // it that rounding keeps from deciding: isolated touching points
    // (extrema of the field on its zero set) and crossings of two branches
    // (saddles). Pieces inside a square hold no cells; the branches into a
    // crossing are bridged to it afterwards.
    let mut singular: Vec<Singular> = Vec::new();
    let inside_one = |singular: &[Singular], cell: &Cell| {
        singular.iter().any(|z| {
            cell.lo.x >= z.at.x - z.reach
                && cell.hi.x <= z.at.x + z.reach
                && cell.lo.y >= z.at.y - z.reach
                && cell.hi.y <= z.at.y + z.reach
        })
    };
    // An undecidable piece: a singular point already known near it grows
    // its square to hold the piece; a new one is recorded; anything else is
    // refused.
    let undecided = |singular: &mut Vec<Singular>, cell: &Cell| -> Result<(), TraceRefusal> {
        let far = |z: &Singular| {
            let d = (cell.lo - z.at).abs().max((cell.hi - z.at).abs());
            d.x.max(d.y)
        };
        let Some((at, kind, reach, hessian)) = critical_point(field, cell, limit, limit) else {
            return Err(TraceRefusal::Singular(cell.centre()));
        };
        if let Some(z) = singular
            .iter_mut()
            .find(|z| (z.at - at).abs().max_element() <= z.reach.max(reach))
        {
            z.reach = z.reach.max(far(z)) * 1.000_001;
            // Beyond this the bridges would be too long to trust.
            if z.reach > 2.0 * limit {
                return Err(TraceRefusal::Singular(z.at));
            }
            return Ok(());
        }
        let mut z = Singular {
            at,
            reach,
            kind,
            hessian,
        };
        z.reach = z.reach.max(far(&z)) * 1.000_001;
        singular.push(z);
        many(singular, cap)
    };
    let mut work = 0usize;
    while let Some((cell, depth)) = queue.pop() {
        work += 1;
        if work > budget {
            return Err(TraceRefusal::Budget);
        }
        if !bound(field, &du, &dv, &cell).straddles_zero()
            || exact.is_some_and(|e| e.keeps_sign(cell.lo, cell.hi, (0, 0)).is_some())
        {
            continue;
        }
        if inside_one(&singular, &cell) || tubes.holds(&cell) {
            continue;
        }
        let fu = bound_simple(&du, &cell);
        let fv = bound_simple(&dv, &cell);
        // Monotone along a parameter: by the partial's bound, else exactly.
        let mono = |r: Range, d: (usize, usize)| {
            !r.straddles_zero()
                || exact.is_some_and(|e| e.keeps_sign(cell.lo, cell.hi, d).is_some())
        };
        // Monotone in v: a graph over u. Where both hold, the steeper.
        let axis = match (mono(fu, (1, 0)), mono(fv, (0, 1))) {
            (true, true) => {
                if steeper(fv, &cell, false) >= steeper(fu, &cell, true) {
                    Some(Axis::U)
                } else {
                    Some(Axis::V)
                }
            }
            (false, true) => Some(Axis::U),
            (true, false) => Some(Axis::V),
            (false, false) => None,
        };
        match axis {
            Some(axis) => match cells_in(field, &du, &dv, &cell, axis, exact) {
                Ok(found) => all.extend(found),
                Err(()) => {
                    // A root on a side could not be certified simple: move
                    // the side by splitting off-centre.
                    if depth > 60 {
                        // Where the sides' roots merge: a singular point.
                        undecided(&mut singular, &cell)?;
                        continue;
                    }
                    for piece in split(cell, 0.4129) {
                        queue.push((piece, depth + 1));
                    }
                }
            },
            None => {
                if (cell.hi.x - cell.lo.x) <= smallest.x && (cell.hi.y - cell.lo.y) <= smallest.y {
                    undecided(&mut singular, &cell)?;
                    continue;
                }
                // A small piece where both partials may vanish: a singular
                // point of the curve may lie in it. Found now, its square
                // spares the subdivision that would otherwise run down to
                // it along branches too close to tell apart cheaply (where
                // two branches touch, the cells shrink as the square of
                // the distance).
                let side = (cell.hi.x - cell.lo.x).max(cell.hi.y - cell.lo.y);
                if side <= 0.25 * limit {
                    if let Some((at, kind, reach, hessian)) =
                        critical_point(field, &cell, limit, 4.0 * side)
                    {
                        let off = (at - cell.centre()).abs().max_element();
                        let known = singular
                            .iter()
                            .any(|z| (z.at - at).abs().max_element() <= z.reach);
                        if off <= side && !known {
                            let reach = match kind {
                                Critical::Degenerate => limit,
                                _ => reach.min(limit),
                            };
                            singular.push(Singular {
                                at,
                                reach: reach.max(side),
                                kind,
                                hessian,
                            });
                            many(&singular, cap)?;
                            if inside_one(&singular, &cell) {
                                continue;
                            }
                        }
                    }
                }
                for piece in split(cell, 0.5) {
                    queue.push((piece, depth + 1));
                }
            }
        }
    }
    // Cells made before a square grew to hold them are its now: every
    // branch then stops at the square, where it is bridged in.
    all.retain(|cell| {
        let Some((a, b)) = ends(field, cell) else {
            return true;
        };
        let m = (a + b) * 0.5;
        !singular
            .iter()
            .any(|z| (m - z.at).abs().max_element() < z.reach)
    });
    let mut curves = chain(field, all, domain, periodic);
    let mut touches = Vec::new();
    for z in &singular {
        if bridge(&mut curves, z, periodic)? {
            tubes.degenerate.set(true);
        }
        // A point with no branch into it is where the surfaces only touch.
        let into = curves.iter().any(|c| {
            [c.point(0.0), c.point(c.end())]
                .into_iter()
                .flatten()
                .any(|e| (e - z.at).abs().max_element() <= 1e-12 * (1.0 + z.at.abs().max_element()))
        });
        if !into {
            touches.push(z.at);
        }
    }
    Ok((curves, touches))
}

/// Bridges into a point where the surfaces touch to higher order. The
/// Hessian gives no directions there, so the ends are checked against the
/// field itself: as many as its sign changes round a square about the
/// point, outside the part rounding hides. Each end is bridged along its
/// own direction to the point; the bridge lies where the field is below
/// its rounding, on both surfaces to that rounding.
fn bridge_degenerate(
    curves: &mut [ImplicitCurve2],
    z: &Singular,
    ends: &[(usize, bool, Point2)],
    near: impl Fn(Point2) -> Point2,
) -> Result<(), TraceRefusal> {
    let refuse = || TraceRefusal::Singular(z.at);
    let Some(field) = curves.first().map(|c| c.field.clone()) else {
        return if ends.is_empty() {
            Ok(())
        } else {
            Err(refuse())
        };
    };
    // Sign changes round the square three times the hidden part's size.
    let r = 3.0 * z.reach;
    let n = 4096;
    let at = |k: usize| {
        let s = 4.0 * k as Scalar / n as Scalar;
        let (side, f) = ((s.floor() as usize) % 4, s.fract());
        let w = -1.0 + 2.0 * f;
        let offset = match side {
            0 => Point2::new(w, -1.0),
            1 => Point2::new(1.0, w),
            2 => Point2::new(-w, 1.0),
            _ => Point2::new(-1.0, -w),
        };
        z.at + offset * r
    };
    let mut changes = 0;
    let mut last = field.value(at(0));
    for k in 1..=n {
        let now = field.value(at(k % n));
        if now != 0.0 && last != 0.0 && (now < 0.0) != (last < 0.0) {
            changes += 1;
        }
        if now != 0.0 {
            last = now;
        }
    }
    if changes != ends.len() || changes % 2 != 0 {
        return Err(refuse());
    }
    // Where the Hessian keeps one direction (a tacnode: branches tangent
    // to each other), every branch arrives along the other, its null
    // direction; with none (a higher crossing), along its own chord.
    let (uu, uv, vv) = z.hessian;
    let trace = uu + vv;
    let root = ((uu - vv) * (uu - vv) + 4.0 * uv * uv).sqrt();
    let (big, small) = (
        0.5 * (trace + root.copysign(trace)),
        0.5 * (trace - root.copysign(trace)),
    );
    let null = if big.abs() > 1e3 * small.abs() && big != 0.0 {
        // The eigenvector of the small eigenvalue.
        let v = if (uu - small).abs() >= (vv - small).abs() {
            Vec2::new(-uv, uu - small)
        } else {
            Vec2::new(vv - small, -uv)
        };
        (v.length() > 0.0).then(|| v.normalize())
    } else {
        None
    };
    for &(k, at_start, p) in ends {
        let c = near(p);
        let curve = &mut curves[k];
        let g = curve.field.jet(p).gradient;
        let leaving = Vec2::new(-g.y, g.x);
        let chord = c - p;
        let into = null.unwrap_or(Vec2::new(chord.x, chord.y));
        let cell = ImplicitCell::bridge(p, c, leaving, into);
        if at_start {
            curve.cells.insert(0, cell.reversed());
        } else {
            curve.cells.push(cell);
        }
    }
    Ok(())
}

/// More singular points than `cap`: a whole curve of them, where the
/// surfaces may be tangent along it. Refused at once, for
/// [`along_contact`] to take over.
fn many(singular: &[Singular], cap: usize) -> Result<(), TraceRefusal> {
    match singular.last() {
        Some(z) if singular.len() > cap => Err(TraceRefusal::Many(z.at)),
        _ => Ok(()),
    }
}

/// A singular point of the curve and the square about it that rounding
/// keeps from deciding.
#[derive(Debug, Clone, Copy)]
struct Singular {
    at: Point2,
    /// Half the square's side.
    reach: Scalar,
    kind: Critical,
    /// `(F_uu, F_uv, F_vv)` there.
    hessian: (Scalar, Scalar, Scalar),
}

/// Join the branches that stop at a singular point's square to the point
/// itself, each by a straight bridge ([`ImplicitCell::bridge`]); every
/// branch then ends at the point, a vertex of the section.
///
/// A crossing must have exactly four branch ends at its square, two along
/// each of the directions where the Hessian's form vanishes, on opposite
/// sides; a touching point must have none. Anything else is refused.
fn bridge(
    curves: &mut [ImplicitCurve2],
    z: &Singular,
    periodic: Periodic,
) -> Result<bool, TraceRefusal> {
    let refuse = || TraceRefusal::Singular(z.at);
    // The point at the turns nearest `p`.
    let near = |p: Point2| {
        let mut c = z.at;
        if periodic.u {
            c.x += ((p.x - c.x) / TAU).round() * TAU;
        }
        if periodic.v {
            c.y += ((p.y - c.y) / TAU).round() * TAU;
        }
        c
    };
    let slack = 2.0 * z.reach;
    let close = |p: Point2| (p - near(p)).abs().max_element() <= slack;
    // (curve, at its start?, end point)
    let mut ends: Vec<(usize, bool, Point2)> = Vec::new();
    for (k, curve) in curves.iter().enumerate() {
        let (Some(a), Some(b)) = (curve.point(0.0), curve.point(curve.end())) else {
            continue;
        };
        // A closed loop has no ends.
        if (a - b).length() <= 1e-9 * (1.0 + a.length()) {
            continue;
        }
        if close(a) {
            ends.push((k, true, a));
        }
        if close(b) {
            ends.push((k, false, b));
        }
    }
    match z.kind {
        Critical::Extremum if ends.is_empty() => return Ok(false),
        Critical::Extremum => return Err(refuse()),
        // A crossing whose ends do not follow the Hessian's directions -- a
        // near-degenerate one, whose point rounding moves off the true one
        // -- is checked against the field instead.
        Critical::Saddle if ends.len() != 4 => {
            return bridge_degenerate(curves, z, &ends, near).map(|()| true)
        }
        Critical::Saddle => {}
        Critical::Degenerate => return bridge_degenerate(curves, z, &ends, near).map(|()| true),
    }
    // The two directions where `uu x^2 + 2 uv x y + vv y^2 = 0`.
    let (uu, uv, vv) = z.hessian;
    let lines: [Point2; 2] = if uu.abs() >= vv.abs() {
        let root = (uv * uv - uu * vv).max(0.0).sqrt();
        [
            Point2::new(-uv + root, uu).normalize(),
            Point2::new(-uv - root, uu).normalize(),
        ]
    } else {
        let root = (uv * uv - uu * vv).max(0.0).sqrt();
        [
            Point2::new(vv, -uv + root).normalize(),
            Point2::new(vv, -uv - root).normalize(),
        ]
    };
    // Each end along one of them, two per direction, on opposite sides.
    let mut sides = [[0usize; 2]; 2];
    let mut along = Vec::with_capacity(4);
    for &(_, _, p) in &ends {
        let d = (p - near(p)).normalize();
        let (k, cos) = if d.dot(lines[0]).abs() >= d.dot(lines[1]).abs() {
            (0, d.dot(lines[0]))
        } else {
            (1, d.dot(lines[1]))
        };
        // Within about 25 degrees of the direction.
        if cos.abs() < 0.9 {
            return bridge_degenerate(curves, z, &ends, near).map(|()| true);
        }
        sides[k][usize::from(cos > 0.0)] += 1;
        along.push(lines[k]);
    }
    if sides != [[1, 1], [1, 1]] {
        return bridge_degenerate(curves, z, &ends, near).map(|()| true);
    }
    for ((k, at_start, p), into) in ends.into_iter().zip(along) {
        let c = near(p);
        let curve = &mut curves[k];
        // The branch's tangent at its certified end, from the field.
        let g = curve.field.jet(p).gradient;
        let leaving = Vec2::new(-g.y, g.x);
        let cell = ImplicitCell::bridge(p, c, leaving, into);
        if at_start {
            curve.cells.insert(0, cell.reversed());
        } else {
            curve.cells.push(cell);
        }
    }
    Ok(false)
}

/// What kind of critical point of the field lies on the curve.
#[derive(Debug, Clone, Copy)]
enum Critical {
    /// Indefinite Hessian: two branches cross.
    Saddle,
    /// Definite Hessian: an isolated point of the zero set.
    Extremum,
    /// Singular Hessian: the surfaces touch to higher order. Any even
    /// number of branches may end there, or none.
    Degenerate,
}

/// The field's critical point on its zero set near `cell`, found by Newton
/// on the gradient from the cell's centre, when it is non-degenerate: the
/// field vanishes there to its rounding, the gradient too, and the Hessian
/// is regular. `None` for anything else (a degenerate singularity, or two
/// branches closer than the trace resolves), which stays a refusal.
#[allow(clippy::type_complexity)]
fn critical_point(
    field: &Field2,
    cell: &Cell,
    limit: Scalar,
    leash: Scalar,
) -> Option<(Point2, Critical, Scalar, (Scalar, Scalar, Scalar))> {
    let mut p = cell.centre();
    let size = (cell.hi - cell.lo).length().max(1e-300);
    // Newton on the gradient, damped (Levenberg-Marquardt) so that it
    // still converges, if only linearly, where the Hessian is singular.
    for _ in 0..200 {
        let jet = field.jet(p);
        let g = jet.gradient;
        let (a, b, c) = (jet.uu, jet.uv, jet.vv);
        let norm = a.abs() + b.abs() + c.abs();
        if !norm.is_finite() {
            return None;
        }
        // (H^T H + mu I) step = H^T g, H symmetric.
        let mu = 1e-12 * norm * norm;
        let (m00, m01, m11) = (a * a + b * b + mu, a * b + b * c, b * b + c * c + mu);
        let (r0, r1) = (a * g.x + b * g.y, b * g.x + c * g.y);
        let det = m00 * m11 - m01 * m01;
        if det == 0.0 || !det.is_finite() {
            break;
        }
        let step = Point2::new((m11 * r0 - m01 * r1) / det, (m00 * r1 - m01 * r0) / det);
        p -= step;
        // Wandered off: no critical point of this piece's.
        if (p - cell.centre()).abs().max_element() > leash {
            return None;
        }
        if step.length() <= 1e-15 * (1.0 + p.length()) {
            break;
        }
    }
    // On the curve to its rounding. A regular Hessian hides the zero set
    // within `sqrt(2 m / lambda)` of the point (m the rounding, lambda the
    // smaller curvature); a singular one further, which the square's growth
    // below takes care of.
    let jet = field.jet(p);
    let rounding = 64.0 * Scalar::EPSILON * field.magnitude().max(field.scale_at(p)).max(1.0);
    if jet.value.abs() > 16.0 * rounding || !p.is_finite() {
        return None;
    }
    let det = jet.uu * jet.vv - jet.uv * jet.uv;
    let norm = jet.uu.abs() + jet.vv.abs() + jet.uv.abs();
    let kind = if det < -1e-9 * norm * norm {
        Critical::Saddle
    } else if det > 1e-9 * norm * norm {
        Critical::Extremum
    } else {
        Critical::Degenerate
    };
    let reach = match kind {
        Critical::Degenerate => 8.0 * size,
        _ => {
            let trace = jet.uu + jet.vv;
            let root = (trace * trace - 4.0 * det).max(0.0).sqrt();
            let lambda = (0.5 * (trace.abs() - root)).abs().max(1e-300);
            4.0 * (2.0 * rounding / lambda).sqrt() + 8.0 * size
        }
    };
    // Interval bounds can be too loose to certify cells some way further
    // out, where the branches are still well apart: the square then grows
    // to that distance, up to `limit`.
    let away = (p - cell.centre()).length();
    if away > limit {
        return None;
    }
    Some((p, kind, reach.max(away + size), (jet.uu, jet.uv, jet.vv)))
}

/// How strongly a partial's bound keeps the field monotone, scaled by the
/// cell's size along that parameter.
fn steeper(r: Range, cell: &Cell, along_u: bool) -> Scalar {
    let width = if along_u {
        cell.hi.x - cell.lo.x
    } else {
        cell.hi.y - cell.lo.y
    };
    r.lo.abs().min(r.hi.abs()) * width
}

fn split_grid(domain: Cell, n: usize) -> Vec<(Cell, u32)> {
    let line = |k: usize, lo: Scalar, hi: Scalar| -> Scalar {
        if k == 0 {
            lo
        } else if k == n {
            hi
        } else {
            lo + (hi - lo) * ((k as Scalar + 0.0371) / n as Scalar)
        }
    };
    let mut out = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let lo = Point2::new(
                line(i, domain.lo.x, domain.hi.x),
                line(j, domain.lo.y, domain.hi.y),
            );
            let hi = Point2::new(
                line(i + 1, domain.lo.x, domain.hi.x),
                line(j + 1, domain.lo.y, domain.hi.y),
            );
            out.push((Cell { lo, hi }, 0));
        }
    }
    out
}

/// Four pieces, split at `fraction` of each side.
fn split(cell: Cell, fraction: Scalar) -> [Cell; 4] {
    let m = cell.lo + (cell.hi - cell.lo) * fraction;
    [
        Cell { lo: cell.lo, hi: m },
        Cell {
            lo: Point2::new(m.x, cell.lo.y),
            hi: Point2::new(cell.hi.x, m.y),
        },
        Cell {
            lo: Point2::new(cell.lo.x, m.y),
            hi: Point2::new(m.x, cell.hi.y),
        },
        Cell { lo: m, hi: cell.hi },
    ]
}

/// The field along one side of a regular piece, as a function of the free
/// parameter: its value at `x` and a bound over `[a, b]`.
struct Side<'a> {
    field: &'a Field2,
    d_free: &'a Field2,
    axis: Axis,
    fixed: Scalar,
    /// The shortest stretch of the side worth halving: below it, roots
    /// too close to separate make the piece split instead, which sees
    /// them further apart relative to its size.
    resolution: Scalar,
    /// Stretches examined so far, against [`SIDE_BUDGET`]: where the
    /// field's bound straddles zero along much of the side (branches
    /// nearly touching), a smaller piece is cheaper than a deep search.
    work: core::cell::Cell<usize>,
    /// The field's exact tier, where it has one, and the free
    /// parameter's derivative order in it.
    exact: Option<Exact<'a>>,
    d_order: (usize, usize),
}

/// Stretches of one side examined before its piece is split instead.
const SIDE_BUDGET: usize = 256;

impl Side<'_> {
    fn at(&self, x: Scalar) -> Point2 {
        match self.axis {
            Axis::U => Point2::new(x, self.fixed),
            Axis::V => Point2::new(self.fixed, x),
        }
    }

    fn value(&self, x: Scalar) -> Scalar {
        self.field.value(self.at(x))
    }

    /// The field's sign at `x`: from its `f64` value where that clears the
    /// rounding, else exactly where the field has an exact tier, else
    /// `None`.
    fn sign(&self, x: Scalar) -> Option<Sign> {
        let p = self.at(x);
        let f = self.field.value(p);
        let rounding = 64.0 * Scalar::EPSILON * self.field.scale_at(p).max(self.field.magnitude());
        if f.abs() > rounding {
            return Some(if f < 0.0 {
                Sign::Negative
            } else {
                Sign::Positive
            });
        }
        self.exact?.sign_at(p, (0, 0))
    }

    /// Whether the field keeps one sign over `[a, b]`: its bound, else its
    /// exact tier.
    fn excludes_zero(&self, a: Scalar, b: Scalar) -> bool {
        if !self.range(a, b).straddles_zero() {
            return true;
        }
        let c = self.cell(a, b);
        self.exact
            .is_some_and(|e| e.keeps_sign(c.lo, c.hi, (0, 0)).is_some())
    }

    /// Whether the field is strictly monotone over `[a, b]`.
    fn monotone(&self, a: Scalar, b: Scalar) -> bool {
        if !self.slope(a, b).straddles_zero() {
            return true;
        }
        let c = self.cell(a, b);
        self.exact
            .is_some_and(|e| e.keeps_sign(c.lo, c.hi, self.d_order).is_some())
    }

    fn cell(&self, a: Scalar, b: Scalar) -> Cell {
        let (p, q) = (self.at(a), self.at(b));
        Cell {
            lo: p.min(q),
            hi: p.max(q),
        }
    }

    /// A bound of the field over `[a, b]`: the tighter of the direct
    /// bound and the mean-value one (the middle's value, give or take the
    /// slope's bound times the half-width). Where the field nearly cancels
    /// over a wide stretch, only the second is narrow.
    fn range(&self, a: Scalar, b: Scalar) -> Range {
        let direct = bound_simple(self.field, &self.cell(a, b));
        let m = 0.5 * (a + b);
        let slope = self.slope(a, b);
        let reach = slope.lo.abs().max(slope.hi.abs()) * 0.5 * (b - a).abs();
        let at = self.at(m);
        let margin = 64.0 * Scalar::EPSILON * self.field.scale_at(at).max(self.field.magnitude());
        let f = self.value(m);
        let spread = reach + margin;
        if !spread.is_finite() {
            return direct;
        }
        Range {
            lo: direct.lo.max(f - spread),
            hi: direct.hi.min(f + spread),
        }
    }

    fn slope(&self, a: Scalar, b: Scalar) -> Range {
        bound_simple(self.d_free, &self.cell(a, b))
    }

    /// The roots in `(a, b)`, each certified simple; `Err` where a root
    /// cannot be separated or certified.
    fn roots(&self, a: Scalar, b: Scalar, out: &mut Vec<Scalar>, depth: u32) -> Result<(), ()> {
        self.work.set(self.work.get() + 1);
        if self.work.get() > SIDE_BUDGET {
            return Err(());
        }
        if self.excludes_zero(a, b) {
            return Ok(());
        }
        if self.monotone(a, b) {
            // Monotone: one root exactly where the signs differ -- signs
            // taken from values clear of their rounding, or exactly.
            let (Some(sa), Some(sb)) = (self.sign(a), self.sign(b)) else {
                // Below rounding at an end, with no exact tier: which side
                // the root lies on cannot be told here; a smaller piece may.
                return Err(());
            };
            if sa == Sign::Zero || sb == Sign::Zero {
                // A root at an end: it belongs to the neighbouring piece's
                // side as much as to this one, and is found as an end there.
                return Err(());
            }
            if sa != sb {
                out.push(self.refine(a, b, sa));
            }
            return Ok(());
        }
        if depth > 60
            || (b - a).abs() <= self.resolution
            || (b - a).abs() <= 1e-13 * (1.0 + a.abs().max(b.abs()))
        {
            return Err(());
        }
        let m = 0.5 * (a + b);
        self.roots(a, m, out, depth + 1)?;
        self.roots(m, b, out, depth + 1)
    }

    /// The root in `[a, b]`, where the field is monotone and changes sign
    /// (`sa` its sign at `a`), bisected on certain signs where they can be
    /// had.
    fn refine(&self, mut a: Scalar, mut b: Scalar, sa: Sign) -> Scalar {
        for _ in 0..200 {
            let m = 0.5 * (a + b);
            if m <= a.min(b) || m >= a.max(b) {
                break;
            }
            // Below rounding with no exact tier, the `f64` sign is the best
            // there is: the root is placed to within that noise.
            let sign = self.sign(m).unwrap_or_else(|| {
                let f = self.value(m);
                if f == 0.0 {
                    Sign::Zero
                } else if f < 0.0 {
                    Sign::Negative
                } else {
                    Sign::Positive
                }
            });
            match sign {
                Sign::Zero => return m,
                s if s == sa => a = m,
                _ => b = m,
            }
        }
        0.5 * (a + b)
    }
}

/// The cells of the curve in a regular piece, where the field is monotone
/// along the solved parameter (the one that is not `axis`).
fn cells_in(
    field: &Field2,
    du: &Field2,
    dv: &Field2,
    cell: &Cell,
    axis: Axis,
    exact: Option<Exact<'_>>,
) -> Result<Vec<ImplicitCell>, ()> {
    let d_order = match axis {
        Axis::U => (1, 0),
        Axis::V => (0, 1),
    };
    let (d_free, lo_free, hi_free, lo_solved, hi_solved) = match axis {
        Axis::U => (du, cell.lo.x, cell.hi.x, cell.lo.y, cell.hi.y),
        Axis::V => (dv, cell.lo.y, cell.hi.y, cell.lo.x, cell.hi.x),
    };
    let resolution = 1e-7 * (hi_free - lo_free).abs();
    let low = Side {
        field,
        d_free,
        axis,
        fixed: lo_solved,
        resolution,
        work: core::cell::Cell::new(0),
        exact,
        d_order,
    };
    let high = Side {
        field,
        d_free,
        axis,
        fixed: hi_solved,
        resolution,
        work: core::cell::Cell::new(0),
        exact,
        d_order,
    };
    let mut breaks = vec![lo_free, hi_free];
    low.roots(lo_free, hi_free, &mut breaks, 0)?;
    high.roots(lo_free, hi_free, &mut breaks, 0)?;
    breaks.sort_by(Scalar::total_cmp);
    let mut out = Vec::new();
    for pair in breaks.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        if b - a <= 1e-14 * (1.0 + a.abs().max(b.abs())) {
            continue;
        }
        let m = 0.5 * (a + b);
        // Signs from values clear of their rounding, or exactly; neither
        // to be had, a smaller piece may tell.
        let (Some(s0), Some(s1)) = (low.sign(m), high.sign(m)) else {
            return Err(());
        };
        if s0 == Sign::Zero || s1 == Sign::Zero {
            return Err(());
        }
        if s0 != s1 {
            out.push(ImplicitCell {
                axis,
                from: a,
                to: b,
                low: lo_solved,
                high: hi_solved,
                bridge: None,
            });
        }
    }
    Ok(out)
}

/// A cell's end points in `(u, v)`, at its `from` and `to`.
fn ends(field: &Field2, cell: &ImplicitCell) -> Option<(Point2, Point2)> {
    let curve = ImplicitCurve2 {
        field: field.clone(),
        cells: vec![*cell],
    };
    Some((curve.point(0.0)?, curve.point(1.0)?))
}

/// Shift a cell by whole periods.
fn shifted(cell: &ImplicitCell, du: Scalar, dv: Scalar) -> ImplicitCell {
    let (df, ds) = match cell.axis {
        Axis::U => (du, dv),
        Axis::V => (dv, du),
    };
    ImplicitCell {
        from: cell.from + df,
        to: cell.to + df,
        low: cell.low + ds,
        high: cell.high + ds,
        ..*cell
    }
}

fn reversed(cell: &ImplicitCell) -> ImplicitCell {
    cell.reversed()
}

/// Chain cells end to end into components.
fn chain(
    field: &Field2,
    cells: Vec<ImplicitCell>,
    domain: Cell,
    periodic: Periodic,
) -> Vec<ImplicitCurve2> {
    let points: Vec<(Point2, Point2)> = cells
        .iter()
        .map(|c| ends(field, c).unwrap_or((Point2::splat(Scalar::NAN), Point2::splat(Scalar::NAN))))
        .collect();
    let extent = (domain.hi - domain.lo).abs();
    let eps = 1e-8 * (1.0 + extent.x.max(extent.y));
    // Wrap a point into the domain for matching across periods.
    let wrap = |p: Point2| -> Point2 {
        let mut q = p;
        if periodic.u {
            q.x = domain.lo.x + (q.x - domain.lo.x).rem_euclid(TAU);
            if (q.x - domain.lo.x - TAU).abs() <= eps {
                q.x = domain.lo.x;
            }
        }
        if periodic.v {
            q.y = domain.lo.y + (q.y - domain.lo.y).rem_euclid(TAU);
            if (q.y - domain.lo.y - TAU).abs() <= eps {
                q.y = domain.lo.y;
            }
        }
        q
    };
    let near = |p: Point2, q: Point2| {
        let (a, b) = (wrap(p), wrap(q));
        let mut d = a - b;
        if periodic.u && d.x.abs() > TAU - eps {
            d.x = d.x.abs() - TAU;
        }
        if periodic.v && d.y.abs() > TAU - eps {
            d.y = d.y.abs() - TAU;
        }
        d.length() <= eps
    };
    // Which other cell end meets this one.
    let find = |me: usize, p: Point2, used: &[bool]| -> Option<(usize, bool)> {
        let mut best: Option<(usize, bool, Scalar)> = None;
        for (k, (a, b)) in points.iter().enumerate() {
            if k == me || used[k] {
                continue;
            }
            for (at_start, q) in [(true, *a), (false, *b)] {
                if near(p, q) {
                    let d = (wrap(p) - wrap(q)).length();
                    if best.is_none_or(|(_, _, e)| d < e) {
                        best = Some((k, at_start, d));
                    }
                }
            }
        }
        best.map(|(k, s, _)| (k, s))
    };
    let mut used = vec![false; cells.len()];
    let mut out = Vec::new();
    for first in 0..cells.len() {
        if used[first] || points[first].0.x.is_nan() {
            continue;
        }
        used[first] = true;
        let mut chain = vec![cells[first]];
        let mut tail = points[first].1;
        let mut closed = false;
        // Forward.
        loop {
            if near(tail, points[first].0) && chain.len() > 1 {
                closed = true;
                break;
            }
            let Some((k, at_start)) = find(usize::MAX, tail, &used) else {
                break;
            };
            used[k] = true;
            let (next, (start, end)) = if at_start {
                (cells[k], points[k])
            } else {
                (reversed(&cells[k]), (points[k].1, points[k].0))
            };
            // Unwrap: shift the cell so it starts where the chain is.
            let du = if periodic.u {
                ((tail.x - start.x) / TAU).round() * TAU
            } else {
                0.0
            };
            let dv = if periodic.v {
                ((tail.y - start.y) / TAU).round() * TAU
            } else {
                0.0
            };
            chain.push(shifted(&next, du, dv));
            tail = Point2::new(end.x + du, end.y + dv);
        }
        if !closed {
            // Backward from the first cell's start.
            let mut head = points[first].0;
            loop {
                let Some((k, at_start)) = find(usize::MAX, head, &used) else {
                    break;
                };
                used[k] = true;
                // The cell must end at `head`.
                let (prev, (start, end)) = if at_start {
                    (reversed(&cells[k]), (points[k].1, points[k].0))
                } else {
                    (cells[k], points[k])
                };
                let du = if periodic.u {
                    ((head.x - end.x) / TAU).round() * TAU
                } else {
                    0.0
                };
                let dv = if periodic.v {
                    ((head.y - end.y) / TAU).round() * TAU
                } else {
                    0.0
                };
                chain.insert(0, shifted(&prev, du, dv));
                head = Point2::new(start.x + du, start.y + dv);
            }
        }
        out.push(ImplicitCurve2 {
            field: field.clone(),
            cells: chain,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_curve::implicit::{Basis, SeriesField2};

    /// `sin(u)^k` in harmonics along `u` (`k` odd:
    /// `2^(1-k) sum_j (-1)^((k-1)/2 - j) C(k, j) sin((k - 2j) u)`; `k` even:
    /// `2^(1-k) (C(k, k/2)/2 + sum_j (-1)^(k/2 - j) C(k, j) cos((k - 2j) u))`),
    /// times `along_v` (coefficients of `1, cos v, sin v`).
    fn sine_power(k: usize, along_v: &[Scalar]) -> Field2 {
        let mut u = vec![0.0; 2 * k + 1];
        let mut binomial = 1.0;
        let scale = 2f64.powi(1 - k as i32);
        for j in 0..=k / 2 {
            let w = k - 2 * j;
            if k % 2 == 1 {
                let sign = if ((k - 1) / 2 - j) % 2 == 0 {
                    1.0
                } else {
                    -1.0
                };
                u[2 * w] = sign * binomial * scale;
            } else if w == 0 {
                u[0] = binomial * scale / 2.0;
            } else {
                let sign = if (k / 2 - j) % 2 == 0 { 1.0 } else { -1.0 };
                u[2 * w - 1] = sign * binomial * scale;
            }
            binomial = binomial * (k - j) as Scalar / (j + 1) as Scalar;
        }
        Field2::Series(SeriesField2 {
            u: Basis::Fourier,
            v: if along_v.len() > 1 {
                Basis::Fourier
            } else {
                Basis::Power
            },
            coefficients: u
                .iter()
                .map(|c| along_v.iter().map(|b| c * b).collect())
                .collect(),
        })
    }

    fn square() -> Cell {
        Cell {
            lo: Point2::new(-1.0, -1.0),
            hi: Point2::new(1.0, 1.0),
        }
    }

    const OPEN: Periodic = Periodic { u: false, v: false };

    #[test]
    fn sine_powers_are_what_they_say() {
        for k in [9, 12, 13] {
            let f = sine_power(k, &[2.0, 1.0]);
            for (x, y) in [(0.3, 0.2), (-0.7, 0.9), (1.1, -0.4)] {
                let exact = f64::sin(x).powi(k as i32) * (2.0 + f64::cos(y));
                assert!((f.value(Point2::new(x, y)) - exact).abs() < 1e-14, "{k}");
            }
        }
    }

    #[test]
    fn a_line_of_contact_hidden_by_rounding_is_traced() {
        // sin(u)^13 crosses zero along u = 0, vanishing to 13th order: in
        // `f64` its value is below rounding out to |u| ~ 0.12, more than a
        // twentieth of the window each side. Refused before (#181); now
        // its signs there are certified, and the line is traced.
        let curves = trace(&sine_power(13, &[1.0]), square(), OPEN).expect("the line");
        assert_eq!(curves.len(), 1);
        let (mut lo, mut hi) = (Scalar::INFINITY, Scalar::NEG_INFINITY);
        for i in 0..=200 {
            let p = curves[0]
                .point(curves[0].end() * i as Scalar / 200.0)
                .expect("a point");
            assert!(p.x.abs() < 1e-12, "{p:?}");
            lo = lo.min(p.y);
            hi = hi.max(p.y);
        }
        assert!(lo < -1.0 + 1e-9 && hi > 1.0 - 1e-9, "{lo} .. {hi}");
    }

    #[test]
    fn a_hidden_line_of_contact_varying_along_it_is_traced() {
        // sin(u)^9 (2 + cos v): harmonics along both parameters, so the
        // certified boxes are two-dimensional.
        let curves = trace(&sine_power(9, &[2.0, 1.0]), square(), OPEN).expect("the line");
        assert_eq!(curves.len(), 1);
        for i in 0..=200 {
            let p = curves[0]
                .point(curves[0].end() * i as Scalar / 200.0)
                .expect("a point");
            assert!(p.x.abs() < 1e-12, "{p:?}");
        }
    }

    #[test]
    fn a_hidden_line_of_even_contact_only_touches() {
        // sin(u)^12 is zero along u = 0 and positive either side.
        assert!(matches!(
            trace(&sine_power(12, &[1.0]), square(), OPEN),
            Err(TraceRefusal::Touching(_))
        ));
    }
}
