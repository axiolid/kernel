//! Runs of sections whose planes fold (#361, ADR 0082 amendment).
//!
//! The exact surface of the [parent module](super) is `A_k(s) = O + a L +
//! b U + c T + x_k L_t + y_k U_t`. A run of sections that stands in planes
//! normal to the frame's tangent `N` is a solid, or a sheet, only while the
//! planes keep advancing past every vertex; where a vertex lies on the inner
//! side of a bend at or beyond its centre of curvature, the planes stop
//! advancing, the sections cross and the loft is not a solid. That is
//! invalid input, as a swept disk whose radius reaches the bend's
//! (`pipe::horn_torus`) is.
//!
//! # The condition
//!
//! A vertex at `(a, b)` in the unturned section frame (`a = lat + x_k`, `b =
//! vert + y_k`) with the longitudinal offset `c` moves along `N` at
//!
//! `F = sigma + c' - a beta_L - b beta_U`
//!
//! (`axiolid_reference::station::FrameLaw::bend`: `sigma = O'.N`, `beta_L =
//! L.N'`, `beta_U = U.N'`, in station distance). `F` is affine in `(a, b)`
//! at a fixed `s`, so the smallest `F` of a section stands at a vertex, and
//! each vertex's `(a, b)` is affine in `s` between two stations. The
//! planes fold where `F <= 0`. On a plan arc `F = 1 - a / R` (the curvature
//! radius `R` against the section's lateral extent on the inner side); a
//! vertical arc folds a section whose up extent reaches the vertical
//! radius, in the section frame only, since the plan frame's `U = +Z` never
//! turns toward the horizontal tangent; and a longitudinal offset that
//! recedes faster than the basis advances (`c' <= -sigma`) folds a section
//! on a line.
//!
//! # Exact decision
//!
//! An unturned run whose basis has a bend (a line, circle, polyline, an
//! elevated or banked curve over a line or a circle, an intrinsic curve) is
//! decided per seam-free span of a station interval by bisection. A span is
//! clear once the lower bound of `F` is positive at every vertex over the
//! span's certified coefficient ranges, less the profile's flattening bound
//! `d` (a flattened ring's vertex stands up to `d` off the profile it stands
//! for, and `F` moves by at most `d |beta|`) and the rounding of the
//! arithmetic. Otherwise the span's middle is a witness: the coefficients
//! there are known to their rounding, and a vertex whose `F` is not
//! certified positive within that rounding folds or touches, so the run is
//! refused. A touch is a non-manifold solid: refused alike. A span neither
//! cleared nor refused is halved, to [`MAX_FOLD_DEPTH`] halvings; a span
//! left undecided there lies within rounding or the profile's flattening bound
//! of a fold and is refused too.
//! The bisection visits at most [`MAX_FOLD_SPANS`] spans, then fails with a
//! budget error by name.
//!
//! # Sheets (#407)
//!
//! A sectioned surface or open-section run is decided by the same condition
//! at the vertices of its polylines with margin 0: an open section is its own
//! exact polyline, and the rate is affine in the section coordinates, so the
//! minimum stands at a polyline vertex. The rule, planes that advance past
//! every section point, is stricter than the sheet merely being regular (a
//! sheet wholly past the centre of a bend is regular yet refused; ADR 0082
//! amendment).
//!
//! # Conservative bound
//!
//! A run with a turned section, or whose frames are read by quadrature (an
//! ellipse, a B-spline, an offset) has no closed-form bend. Its turned
//! section plane has the normal `N_t`, with components `n` in the base
//! frame, and a vertex advances along it at
//!
//! `O'.N_t + lat' n_L + vert' n_U + c' n_T + lat L'.N_t + vert U'.N_t + c T'.N_t + x L_t'.N_t + y U_t'.N_t`,
//!
//! at least `sigma n_T - w |n_U| + lo(lat' n_L) + lo(vert' n_U) + lo(c'
//! n_T) - (|lat| |L'| + |vert| |U'| + |c| |T'|) - |x_k| |L_t'| - |y_k|
//! |U_t'| - d (|L_t'| + |U_t'|)`: `sigma` between the certified forward
//! floor and `sup |O'|`, `w` the vertical part of `O'` that the plan frame
//! leaves out (`sqrt(sup |O'|^2 - forward^2)` under the plan reading, where
//! `|O'_h| = O'.N` is at least the forward floor; zero otherwise), `n` from
//! the orientation interpolated at the span's middle and widened by what it
//! can turn within the span. This bound never refuses: where it is not
//! positive after [`CONSERVATIVE_FOLD_DEPTH`] halvings, the run is
//! unbounded by name ([`FOLD_NOT_EXCLUDED`]). Likewise a span without
//! certified frames ([`UNCERTIFIED_FRAMES`]) or with a turn it cannot bound
//! ([`UNBOUNDED_TURN`]).
//!
//! # Local only
//!
//! The condition is local: the section planes advance at every station
//! distance. A run that meets itself far away (a loop of the basis, two legs
//! of a spiral) is not detected here (#407). A turned section or a frame read
//! by quadrature is never refused, only checked conservatively (#406 asks for
//! an exact decision).

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_curve::{Curve3, Line3};
use axiolid_reference::station::{station_section3, BendBounds, FrameReading, SectionFrame};

use super::{
    lerp, lerp_offsets, oriented, seam_slack, turn, turned, SectionLaw, FOLD_NOT_EXCLUDED,
    UNBOUNDED_TURN, UNCERTIFIED_FRAMES,
};

/// The refusal's prefix.
const FOLD: &str = "station-placed sections fold";

/// The halvings of an exact span before it counts as undecided.
const MAX_FOLD_DEPTH: u32 = 48;

/// The spans the exact bisection of a run may visit.
const MAX_FOLD_SPANS: usize = 1 << 16;

/// The halvings of a span the conservative bound may use.
const CONSERVATIVE_FOLD_DEPTH: u32 = 8;

/// A rounding allowance factor, as in the frame bounds.
const ROUNDING: Scalar = 64.0 * Scalar::EPSILON;

/// What a span's decision is: `Ok` clear, `Err` unbounded by name.
type Verdict = Result<(), &'static str>;

/// The smaller and larger of two values.
fn ordered(p: Scalar, q: Scalar) -> [Scalar; 2] {
    [p.min(q), p.max(q)]
}

/// The largest product of a value in `a` and a value in `b`.
fn product_hi(a: [Scalar; 2], b: [Scalar; 2]) -> Scalar {
    let [p, q] = [a[0] * b[0], a[0] * b[1]];
    let [r, t] = [a[1] * b[0], a[1] * b[1]];
    p.max(q).max(r).max(t)
}

/// The smallest product of a value in `a` and a value in `b`.
fn product_lo(a: [Scalar; 2], b: [Scalar; 2]) -> Scalar {
    -product_hi([-a[1], -a[0]], b)
}

fn magnitude(range: [Scalar; 2]) -> Scalar {
    range[0].abs().max(range[1].abs())
}

/// The frame whose axes are the unit vectors of the base: the orientation
/// interpolated onto it presents the turned tangent's components in any
/// base frame.
fn identity() -> GeomResult<SectionFrame> {
    let frame = station_section3(
        &Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::X,
        }),
        0.0,
    )?;
    if frame.tangent == Vec3::X && frame.lateral == Vec3::Y && frame.up == Vec3::Z {
        Ok(frame)
    } else {
        Err(GeomError::InvalidInput(
            "station: the reference frame of +X is not the identity".into(),
        ))
    }
}

/// A span of one station interval, between two distances.
#[derive(Clone, Copy)]
struct Span {
    /// The interval's index.
    i: usize,
    lo: Scalar,
    hi: Scalar,
}

impl<'r, 'g> SectionLaw<'r, 'g> {
    /// Whether the section planes of the run keep advancing past every
    /// vertex, within the profile's flattening bound `margin`: `Ok(())`
    /// decided, `Err` names why it is not (never a fold: see the
    /// [module documentation](self)).
    ///
    /// # Errors
    ///
    /// A fold, or a touch within rounding, as invalid input; a bisection
    /// past its span budget, by name.
    pub(in crate::station) fn fold(&self, margin: Scalar) -> GeomResult<Verdict> {
        let Some(reading) = self.reading else {
            return Ok(Err(UNCERTIFIED_FRAMES));
        };
        let identity = identity()?;
        let mut visited = 0;
        let mut verdict: Verdict = Ok(());
        for i in 0..self.distances.len() - 1 {
            let (d0, d1) = (self.distances[i], self.distances[i + 1]);
            let mut inside: Vec<Scalar> = self
                .seams
                .iter()
                .copied()
                .filter(|&c| c > d0 + seam_slack(d0) && c < d1 - seam_slack(d1))
                .collect();
            inside.sort_by(Scalar::total_cmp);
            let mut cuts = vec![d0];
            cuts.extend(inside);
            cuts.push(d1);
            for pair in cuts.windows(2) {
                let span = Span {
                    i,
                    lo: pair[0],
                    hi: pair[1],
                };
                let one = self.fold_span(span, reading, margin, &identity, &mut visited)?;
                verdict = verdict.and(one);
            }
        }
        Ok(verdict)
    }

    fn fold_span(
        &self,
        span: Span,
        reading: FrameReading,
        margin: Scalar,
        identity: &SectionFrame,
        visited: &mut usize,
    ) -> GeomResult<Verdict> {
        let plain = self.orientations[span.i].is_base() && self.orientations[span.i + 1].is_base();
        // A frame with a closed-form bend is bisected even where its bend
        // is not bounded over the whole span (a span longer than twice the
        // bend's radius); only a frame without one is checked conservatively.
        let verdict = if plain && self.fold_frames.has_bend(span.lo, span.hi) {
            let first = self.fold_frames.bend(span.lo, span.hi, reading);
            self.bisect(span, first, reading, margin, visited)?
        } else {
            self.conservative(span, reading, margin, identity)?
        };
        // The fold is decided on the frames as read, which a placement
        // read inexactly does not change in their section planes (#408);
        // a fold found there is refused, while the run is certified against
        // the exact surface by no cell (`SectionLaw::cell` bounds the frames
        // strictly), so it is never cleared to a bound.
        Ok(verdict)
    }

    /// The vertices' section coordinates `(a, b)` (lateral, up in the
    /// unturned frame) at `s` in station interval `i`.
    fn coordinates(&self, i: usize, s: Scalar) -> Vec<(Scalar, Scalar)> {
        let (d0, d1) = (self.distances[i], self.distances[i + 1]);
        let u = ((s - d0) / (d1 - d0)).clamp(0.0, 1.0);
        let offsets = lerp_offsets(&self.offsets[i], &self.offsets[i + 1], u);
        self.locals[i]
            .iter()
            .zip(&self.locals[i + 1])
            .map(|(p, q)| {
                (
                    offsets.lateral + lerp(p.x, q.x, u),
                    offsets.vertical + lerp(p.y, q.y, u),
                )
            })
            .collect()
    }

    /// The rate of the longitudinal offset in station interval `i`.
    fn long_rate(&self, i: usize) -> Scalar {
        (self.offsets[i + 1].longitudinal - self.offsets[i].longitudinal)
            / (self.distances[i + 1] - self.distances[i])
    }

    /// Bisects `span` over its certified bend `first`; see the [module
    /// documentation](self#exact-decision).
    fn bisect(
        &self,
        span: Span,
        first: Option<BendBounds>,
        reading: FrameReading,
        margin: Scalar,
        visited: &mut usize,
    ) -> GeomResult<Verdict> {
        let long_rate = self.long_rate(span.i);
        let mut verdict: Verdict = Ok(());
        let mut stack = vec![(span.lo, span.hi, 0_u32, first)];
        while let Some((lo, hi, depth, known)) = stack.pop() {
            *visited += 1;
            if *visited > MAX_FOLD_SPANS {
                return Err(GeomError::BudgetExceeded {
                    resource: "station fold check",
                });
            }
            let Some(bend) = known.or_else(|| self.fold_frames.bend(lo, hi, reading)) else {
                // Not bounded over this width: halve it, and where it stays
                // so (a vertical tangent) leave it undecided, the other
                // spans still decided.
                let mid = 0.5 * (lo + hi);
                if depth >= CONSERVATIVE_FOLD_DEPTH || !(mid > lo && mid < hi) {
                    verdict = Err(FOLD_NOT_EXCLUDED);
                } else {
                    stack.push((mid, hi, depth + 1, None));
                    stack.push((lo, mid, depth + 1, None));
                }
                continue;
            };
            let ends = [self.coordinates(span.i, lo), self.coordinates(span.i, hi)];
            if clear(&bend, &ends, long_rate, margin) {
                continue;
            }
            let mid = 0.5 * (lo + hi);
            witness(&bend, &self.coordinates(span.i, mid), long_rate, mid)?;
            if depth >= MAX_FOLD_DEPTH || !(mid > lo && mid < hi) {
                return Err(GeomError::InvalidInput(format!(
                    "{FOLD}, or touch within rounding or the profile's flattening bound of a \
                     fold, between distances {lo} and {hi}: the forward rate of a section \
                     point is not certified positive there, and halving the span {depth} \
                     times cannot decide it; a touch within rounding is refused alike"
                )));
            }
            stack.push((mid, hi, depth + 1, None));
            stack.push((lo, mid, depth + 1, None));
        }
        Ok(verdict)
    }

    /// The conservative decision of `span`; see the [module
    /// documentation](self#conservative-bound).
    fn conservative(
        &self,
        span: Span,
        reading: FrameReading,
        margin: Scalar,
        identity: &SectionFrame,
    ) -> GeomResult<Verdict> {
        let mut stack = vec![(span.lo, span.hi, 0_u32)];
        while let Some((lo, hi, depth)) = stack.pop() {
            let part = Span { lo, hi, ..span };
            // A span the frames cannot be bounded over, or turned through,
            // may be halved into spans they can, like one not cleared.
            let reason = match self.advancing(part, reading, margin, identity)? {
                Ok(true) => continue,
                Ok(false) => FOLD_NOT_EXCLUDED,
                Err(reason) => reason,
            };
            if depth >= CONSERVATIVE_FOLD_DEPTH {
                return Ok(Err(reason));
            }
            let mid = 0.5 * (lo + hi);
            stack.push((mid, hi, depth + 1));
            stack.push((lo, mid, depth + 1));
        }
        Ok(Ok(()))
    }

    /// Whether a lower bound of every vertex's rate along the turned
    /// section's normal is positive over `span`; `Err` names why it cannot
    /// be told.
    fn advancing(
        &self,
        span: Span,
        reading: FrameReading,
        margin: Scalar,
        identity: &SectionFrame,
    ) -> GeomResult<Result<bool, &'static str>> {
        let Span { i, lo, hi } = span;
        let Some(bounds) = self.fold_frames.bounds(lo, hi, reading) else {
            return Ok(Err(UNCERTIFIED_FRAMES));
        };
        let (d0, d1) = (self.distances[i], self.distances[i + 1]);
        let delta = d1 - d0;
        let (u0, u1) = (
            ((lo - d0) / delta).clamp(0.0, 1.0),
            ((hi - d0) / delta).clamp(0.0, 1.0),
        );
        let unit = |o: &axiolid_model::StationOrientation| o.unit_axes().ok();
        let (Some(ta), Some(tb)) = (unit(&self.orientations[i]), unit(&self.orientations[i + 1]))
        else {
            return Ok(Err(UNBOUNDED_TURN));
        };
        let Some(turn) = turn(ta, tb, u0, u1, delta) else {
            return Ok(Err(UNBOUNDED_TURN));
        };
        let plain = self.orientations[i].is_base() && self.orientations[i + 1].is_base();
        // The turned tangent's components in the base frame, widened by
        // what they can move within the span.
        let normal = oriented(
            identity,
            &self.orientations[i],
            &self.orientations[i + 1],
            0.5 * (u0 + u1),
        )?
        .tangent;
        let widen = 0.5 * (hi - lo) * (turn.lateral[0] + turn.up[0]);
        let range = |n: Scalar| [n - widen, n + widen];
        let (n_t, n_l, n_u) = (range(normal.x), range(normal.y), range(normal.z));
        // The vertical part of `O'` the plan frame leaves out.
        let vertical = if reading == FrameReading::Plan {
            // `|O'_h| = O'.N >= forward`, so what is left is vertical.
            let speed = bounds.point[0];
            (speed * speed - bounds.forward * bounds.forward)
                .max(0.0)
                .sqrt()
                * (1.0 + ROUNDING)
                + ROUNDING * speed
        } else {
            0.0
        };
        let forward = [bounds.forward, bounds.point[0]];
        let (a, b) = (&self.offsets[i], &self.offsets[i + 1]);
        let rates = [
            (b.lateral - a.lateral) / delta,
            (b.vertical - a.vertical) / delta,
            (b.longitudinal - a.longitudinal) / delta,
        ];
        let (o0, o1) = (lerp_offsets(a, b, u0), lerp_offsets(a, b, u1));
        let big = |p: Scalar, q: Scalar| p.abs().max(q.abs());
        let base = [
            big(o0.lateral, o1.lateral),
            big(o0.vertical, o1.vertical),
            big(o0.longitudinal, o1.longitudinal),
        ];
        let swing =
            base[0] * bounds.lateral[0] + base[1] * bounds.up[0] + base[2] * bounds.tangent[0];
        let (f1, f2) = (super::rss(&bounds, 0), super::rss(&bounds, 1));
        let (lt, ut) = if plain {
            (bounds.lateral, bounds.up)
        } else {
            (turned(turn.lateral, [f1, f2]), turned(turn.up, [f1, f2]))
        };
        let common = product_lo(forward, n_t) - vertical * magnitude(n_u)
            + (rates[0] * n_l[0]).min(rates[0] * n_l[1])
            + (rates[1] * n_u[0]).min(rates[1] * n_u[1])
            + (rates[2] * n_t[0]).min(rates[2] * n_t[1])
            - swing
            - margin * (lt[0] + ut[0]);
        let allowance = ROUNDING
            * (forward[1] * magnitude(n_t).max(1.0)
                + vertical
                + rates.iter().map(|r| r.abs()).sum::<Scalar>()
                + swing
                + margin * (lt[0] + ut[0]));
        for (p, q) in self.locals[i].iter().zip(&self.locals[i + 1]) {
            let (x, y) = (
                big(lerp(p.x, q.x, u0), lerp(p.x, q.x, u1)),
                big(lerp(p.y, q.y, u0), lerp(p.y, q.y, u1)),
            );
            let floor =
                common - x * lt[0] - y * ut[0] - allowance - ROUNDING * (x * lt[0] + y * ut[0]);
            if floor.is_nan() || floor <= 0.0 {
                return Ok(Ok(false));
            }
        }
        Ok(Ok(true))
    }
}

/// Whether the forward rate of every vertex is positive over the whole span
/// of `bend`, the vertices moving between `ends` (their coordinates at the
/// span's ends), less the profile's flattening bound and the rounding.
fn clear(
    bend: &BendBounds,
    ends: &[Vec<(Scalar, Scalar)>; 2],
    long_rate: Scalar,
    margin: Scalar,
) -> bool {
    let reach = margin * magnitude(bend.lateral).hypot(magnitude(bend.up));
    ends[0].iter().zip(&ends[1]).all(|(p, q)| {
        let (a, b) = (ordered(p.0, q.0), ordered(p.1, q.1));
        let magnitudes = magnitude(bend.forward)
            + long_rate.abs()
            + magnitude(a) * magnitude(bend.lateral)
            + magnitude(b) * magnitude(bend.up)
            + reach;
        let floor = bend.forward[0] + long_rate
            - product_hi(a, bend.lateral)
            - product_hi(b, bend.up)
            - reach
            - ROUNDING * magnitudes;
        floor > 0.0
    })
}

/// Refuses a vertex of the section at `s` (coordinates `at`) whose forward
/// rate is not certified positive within the rounding of the bend's middle
/// values.
fn witness(
    bend: &BendBounds,
    at: &[(Scalar, Scalar)],
    long_rate: Scalar,
    s: Scalar,
) -> GeomResult<()> {
    let around = |j: usize| {
        [
            bend.middle[j] - bend.middle_error[j],
            bend.middle[j] + bend.middle_error[j],
        ]
    };
    let (forward, lateral, up) = (around(0), around(1), around(2));
    for &(a, b) in at {
        let magnitudes = magnitude(forward)
            + long_rate.abs()
            + a.abs() * magnitude(lateral)
            + b.abs() * magnitude(up);
        let floor = forward[0] + long_rate
            - product_hi([a, a], lateral)
            - product_hi([b, b], up)
            - ROUNDING * magnitudes;
        if floor.is_nan() || floor <= 0.0 {
            let rate = bend.middle[0] + long_rate - a * bend.middle[1] - b * bend.middle[2];
            return Err(GeomError::InvalidInput(format!(
                "{FOLD}: at distance {s} the section point ({a}, {b}) (lateral, up in the \
                 section frame): the section planes stop advancing there (forward rate \
                 {rate:e}, rounding {:e}); a touch within rounding is refused alike",
                rate - floor
            )));
        }
    }
    Ok(())
}
