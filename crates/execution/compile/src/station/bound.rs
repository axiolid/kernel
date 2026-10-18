//! A certified bound on how far a run of sections at stations lies from
//! its mesh (#247, ADR 0082 amendment).
//!
//! # The exact surface
//!
//! Between two stations every section point is interpolated linearly in
//! distance (the parent module's [meshing](super#meshing-between-stations)):
//! vertex `k` of the run stands at
//!
//! `A_k(s) = O(s) + lat L + vert U + long T + x_k L_t + y_k U_t`,
//!
//! `O, T, L, U` the basis point and section frame at `s` (or the plan frame),
//! `L_t, U_t` the frame turned by the interpolated orientation, and the
//! offsets and the profile coordinates `x_k, y_k` linear in `s`. At a fixed
//! `s` the placement is rigid, so the exact section's edge between two
//! vertices is the segment between them: the exact surface is
//! `P(s, v) = (1 - v) A_k(s) + v A_l(s)` over every edge `(k, l)`.
//!
//! # The bound per cell
//!
//! A mesh cell joins the placed sections at two meshed distances `s0 < s1`,
//! each wall quad cut into two triangles (`axiolid_construct::loft`). Let
//! `l_k` be the chord of `A_k` between the mesh's vertices and `B` the
//! bilinear patch over the quad. Then
//!
//! - `|P - B| <= max(e_k, e_l)` with `e_k >= sup |A_k - l_k|`, since `P`
//!   and `B` are the same convex combination of `A` and `l`;
//! - `B` lies within `axiolid_construct::loft::quad_deviation` of the two
//!   triangles.
//!
//! `e_k` is split where a seam of the basis lies inside the cell (a joint
//! whose tangents agree, a change of curvature, a jump of a banked
//! curve's roll), and at the cell's ends where the mesh's vertex was read
//! from another side than the cell's (#301): `A_k` is chorded over each
//! piece, within `h^2 / 8 sup |A_k''|` (`axiolid_reference::bound`), and
//! the polyline of those chords lies within the largest distance of a
//! breakpoint value `A_k(c+-)` from `l_k(c)` of `l_k`, since their
//! difference is piecewise linear. Over a piece,
//!
//! `|A''| <= |O''| + |lat| |L''| + |vert| |U''| + |long| |T''|
//!           + 2 (|lat'| |L'| + |vert'| |U'| + |long'| |T'|)
//!           + |x| |L_t''| + |y| |U_t''| + 2 (|x'| |L_t'| + |y'| |U_t'|)`,
//!
//! the local coordinates at their larger end value and their rates exact,
//! the frame's derivative bounds from
//! `axiolid_reference::station::piece_frame_bounds`. A turned axis `W = sum
//! w_i F_i` (components `w` in the orthonormal base frame `F`) has `|W'| <=
//! |w'| + F_1` and `|W''| <= |w''| + 2 |w'| F_1 + F_2`, `F_j` the
//! root-sum-square of the base axes' bounds; between two orientations
//! `w` comes from unit vectors interpolated linearly and turned by
//! Gram-Schmidt, bounded by the unit-vector rule of
//! `axiolid_reference::station::frame_bound` (a constant orientation
//! adds nothing).
//!
//! Where the basis reads its distance by the arc-length inverse (a circle,
//! a polyline), the frame stands at a distance within `delta` of the one
//! asked for while the local coordinates are read at the latter: each
//! vertex is within `delta Q` of the exact surface (`Q` the rate of its
//! local coordinates), the cells cover the surface up to slivers `delta`
//! long at the run's ends, and each is up to `2 delta` longer: `delta (S +
//! 2 Q)` per vertex, `S` its speed, covers all three.
//!
//! Where it reads its frame by quadrature (an ellipse, a B-spline, an
//! offset beside a line, a circle, an ellipse or a B-spline, #316), the
//! exact surface is the one on the exact frame law
//! (`axiolid_reference::station::FrameLaw`): `delta` bounds how far a
//! read's true distance lies from the one asked for, its certified length
//! measured to a thousandth of the run's budget, and each read point and
//! axis lies within `e_p` and `e_a` of the exact frame's there. A vertex is
//! then within `e_p + e_a (|lat| + |vert| + |long| + c (|x| + |y|))` of
//! its exact position besides `delta (S + 2 Q)`, `c` one for the base
//! frame and `sqrt 3` for a turned one (`|sum w_i dF_i| <= sqrt 3 e_a` for
//! unit components `w`).
//!
//! # Mitres
//!
//! A cell beside a mitre ends on the mitre's plane. Where the piece it
//! lies on is straight and the section constant across the cell (a
//! constant section along a polyline's legs), every vertex runs on a
//! line through the plane and the wall is a planar strip cut by it: the
//! mesh quad is the exact wall, and only its `quad_deviation` counts.
//! Any other cell beside a mitre is not bounded:
//! [`MITRED_SECTIONS`].

use axiolid_contracts::GeomResult;
use axiolid_core::{Point2, Point3, Scalar, Vec3};
use axiolid_model::{SeamSide, StationFrame, StationOffsets, StationOrientation};
use axiolid_reference::station::{FrameLaw, FrameLawBounds, FrameReading};

use super::{interval, lerp, lerp_offsets, oriented, place_turned, seam_slack, Run};

/// Why a run crossing a mitre is not bounded.
pub(crate) const MITRED_SECTIONS: &str =
    "stationed sections: a cell beside a mitred seam on a curved piece, or with a section that \
     changes across it: the mitre's cut of the surface is not bounded";

/// Why a run along a basis without a certified frame law is not bounded.
pub(crate) const UNCERTIFIED_FRAMES: &str =
    "stationed sections: the frame law along this basis has no certified derivative bounds (a \
     frame read on an intrinsic curve or a chain, an offset beside one of them or beside an \
     elevated or banked curve, a skewed circle, a steep tangent, a curve read by quadrature \
     placed in a frame that tilts +Z, or a piece after joints whose distances cannot be \
     certified)";

/// Why an orientation that turns too far between stations is not bounded.
pub(crate) const UNBOUNDED_TURN: &str =
    "stationed sections: the orientation turns between stations through a direction its \
     interpolation cannot bound";

/// The frame reading a run's frame names, `None` for one this compiler
/// does not know.
pub(super) fn reading(frame: StationFrame) -> Option<FrameReading> {
    match frame {
        StationFrame::Section => Some(FrameReading::Section),
        StationFrame::Plan => Some(FrameReading::Plan),
        _ => None,
    }
}

/// A run of sections, as its bound reads it: the run, its stations'
/// distances, offsets and orientations, every vertex's profile point per
/// station, and the wall edges.
pub(super) struct SectionLaw<'r, 'g> {
    pub(super) run: &'r Run<'g>,
    pub(super) distances: &'r [Scalar],
    pub(super) offsets: Vec<StationOffsets>,
    pub(super) orientations: Vec<StationOrientation>,
    pub(super) locals: Vec<Vec<Point2>>,
    pub(super) edges: &'r [(usize, usize)],
    reading: Option<FrameReading>,
    frames: FrameLaw<'r, 'g>,
    seams: Vec<Scalar>,
}

/// The share of a run's budget its frames' quadratures are certified to
/// (#316): their lengths, so a read's measure, to a thousandth of it.
const CERTIFICATION_SHARE: Scalar = 1e-3;

/// `sup |W'|` and `sup |W''|` of a turned axis.
type Pair = [Scalar; 2];

/// The bounds of an orientation interpolated between two stations: the
/// turned lateral's and up's components (`[y1, y2]`, `[n1, n2]`).
struct Turn {
    lateral: Pair,
    up: Pair,
}

/// The smallest length of `p + t (q - p)` over `t` in `[0, 1]`.
fn segment_floor(p: Vec3, q: Vec3) -> Scalar {
    let d = q - p;
    let t = if d.length_squared() > 0.0 {
        (-p.dot(d) / d.length_squared()).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p + d * t).length() * (1.0 - 1e-12)
}

/// The turned frame's component bounds over `[u0, u1]` of an interval
/// `delta` long, between orientations `a` and `b` (unit axis and
/// reference direction each); see the module documentation.
fn turn(a: (Vec3, Vec3), b: (Vec3, Vec3), u0: Scalar, u1: Scalar, delta: Scalar) -> Option<Turn> {
    let lerp3 = |p: Vec3, q: Vec3, u: Scalar| p + (q - p) * u;
    let rate_a = (b.0 - a.0).length() / delta;
    let rate_r = (b.1 - a.1).length() / delta;
    if rate_a == 0.0 && rate_r == 0.0 {
        return Some(Turn {
            lateral: [0.0; 2],
            up: [0.0; 2],
        });
    }
    let rho_a = segment_floor(lerp3(a.0, b.0, u0), lerp3(a.0, b.0, u1));
    let rho_r = segment_floor(lerp3(a.1, b.1, u0), lerp3(a.1, b.1, u1));
    if !(rho_a > 0.0 && rho_r > 0.0) {
        return None;
    }
    let (n1, n2) = (rate_a / rho_a, 6.0 * rate_a * rate_a / (rho_a * rho_a));
    let (r1, r2) = (rate_r / rho_r, 6.0 * rate_r * rate_r / (rho_r * rho_r));
    let p1 = 2.0 * r1 + 2.0 * n1;
    let p2 = 2.0 * r2 + 2.0 * n2 + 4.0 * r1 * n1 + 2.0 * n1 * n1;
    let middle = 0.5 * (u0 + u1);
    let n = lerp3(a.0, b.0, middle).normalize();
    let r = lerp3(a.1, b.1, middle).normalize();
    let rho_p = (r - n * r.dot(n)).length() - 0.5 * (u1 - u0) * delta * p1;
    if !(rho_p.is_finite() && rho_p > 0.0) {
        return None;
    }
    let x1 = p1 / rho_p;
    let x2 = 2.0 * p2 / rho_p + 6.0 * p1 * p1 / (rho_p * rho_p);
    Some(Turn {
        lateral: [n1 + x1, n2 + 2.0 * n1 * x1 + x2],
        up: [n1, n2],
    })
}

/// A turned axis's bounds from its components' and the base frame's
/// root-sum-square `f`.
fn turned(w: Pair, f: Pair) -> Pair {
    [w[0] + f[0], w[1] + 2.0 * w[0] * f[0] + f[1]]
}

fn rss(bounds: &FrameLawBounds, j: usize) -> Scalar {
    (bounds.tangent[j].powi(2) + bounds.lateral[j].powi(2) + bounds.up[j].powi(2)).sqrt()
}

impl<'r, 'g> SectionLaw<'r, 'g> {
    pub(super) fn new(
        run: &'r Run<'g>,
        distances: &'r [Scalar],
        offsets: Vec<StationOffsets>,
        orientations: Vec<StationOrientation>,
        locals: Vec<Vec<Point2>>,
        edges: &'r [(usize, usize)],
        budget: Scalar,
    ) -> GeomResult<Self> {
        let seams = run
            .basis
            .seams()?
            .iter()
            .map(|seam| seam.distance)
            .collect();
        Ok(Self {
            reading: reading(run.frame),
            frames: run.basis.frame_law(CERTIFICATION_SHARE * budget)?,
            run,
            distances,
            offsets,
            orientations,
            locals,
            edges,
            seams,
        })
    }

    /// The section's points at `s`, read from `side` there: the exact
    /// surface's vertices on that side.
    fn exact_at(&self, s: Scalar, side: SeamSide) -> GeomResult<Vec<Point3>> {
        let (i, u) = interval(self.distances, s);
        let base = self.run.basis.section_on(s, self.run.frame, side)?;
        let turned = oriented(&base, &self.orientations[i], &self.orientations[i + 1], u)?;
        let offsets = lerp_offsets(&self.offsets[i], &self.offsets[i + 1], u);
        Ok(self.locals[i]
            .iter()
            .zip(&self.locals[i + 1])
            .map(|(p, q)| {
                let local = Point2::new(lerp(p.x, q.x, u), lerp(p.y, q.y, u));
                place_turned(&base, &turned, &offsets, local)
            })
            .collect())
    }

    fn on_seam(&self, s: Scalar) -> bool {
        self.seams
            .iter()
            .any(|seam| (seam - s).abs() <= seam_slack(s))
    }

    fn is_mitre(&self, s: Scalar) -> bool {
        self.run.mitres.iter().any(|(knot, _)| *knot == s)
    }

    /// A certified bound on the distance from the exact surface between
    /// the meshed distances `lo < hi` to the cell's triangles, whose
    /// vertices are `at_lo` and `at_hi`; `Err` names why it is unbounded.
    pub(super) fn cell(
        &self,
        lo: Scalar,
        hi: Scalar,
        at_lo: &[Point3],
        at_hi: &[Point3],
    ) -> GeomResult<Result<Scalar, &'static str>> {
        let Some(reading) = self.reading else {
            return Ok(Err(UNCERTIFIED_FRAMES));
        };
        let (i, _) = interval(self.distances, 0.5 * (lo + hi));
        let (d0, d1) = (self.distances[i], self.distances[i + 1]);
        let delta = d1 - d0;
        let (u0, u1) = (
            ((lo - d0) / delta).clamp(0.0, 1.0),
            ((hi - d0) / delta).clamp(0.0, 1.0),
        );
        let mitred = self.is_mitre(lo) || self.is_mitre(hi);
        // The breakpoints strictly inside the cell.
        let inside: Vec<Scalar> = self
            .seams
            .iter()
            .copied()
            .filter(|&c| c > lo + seam_slack(lo) && c < hi - seam_slack(hi))
            .collect();
        if mitred && !inside.is_empty() {
            return Ok(Err(MITRED_SECTIONS));
        }
        let mut cuts = Vec::with_capacity(inside.len() + 2);
        cuts.push(lo);
        cuts.extend(&inside);
        cuts.push(hi);
        let mut pieces = Vec::with_capacity(cuts.len() - 1);
        for pair in cuts.windows(2) {
            let Some(bounds) = self.frames.bounds(pair[0], pair[1], reading) else {
                return Ok(Err(UNCERTIFIED_FRAMES));
            };
            if mitred && !bounds.is_straight() {
                return Ok(Err(MITRED_SECTIONS));
            }
            pieces.push((pair[1] - pair[0], bounds));
        }
        // The local coordinates: offsets in the base frame, profile points
        // in the turned one, each linear in distance.
        let (a, b) = (&self.offsets[i], &self.offsets[i + 1]);
        let (o0, o1) = (lerp_offsets(a, b, u0), lerp_offsets(a, b, u1));
        let big = |p: Scalar, q: Scalar| p.abs().max(q.abs());
        let base = [
            big(o0.lateral, o1.lateral),
            big(o0.vertical, o1.vertical),
            big(o0.longitudinal, o1.longitudinal),
        ];
        let base_rate = [
            (b.lateral - a.lateral).abs() / delta,
            (b.vertical - a.vertical).abs() / delta,
            (b.longitudinal - a.longitudinal).abs() / delta,
        ];
        let offset_rate = Vec3::new(
            b.lateral - a.lateral,
            b.vertical - a.vertical,
            b.longitudinal - a.longitudinal,
        )
        .length()
            / delta;
        let unit = |o: &StationOrientation| o.unit_axes().ok();
        let (Some(ta), Some(tb)) = (unit(&self.orientations[i]), unit(&self.orientations[i + 1]))
        else {
            return Ok(Err(UNBOUNDED_TURN));
        };
        let plain = self.orientations[i].is_base() && self.orientations[i + 1].is_base();
        let Some(turn) = turn(ta, tb, u0, u1, delta) else {
            return Ok(Err(UNBOUNDED_TURN));
        };
        let changes = offset_rate > 0.0
            || turn.up[0] > 0.0
            || self.locals[i]
                .iter()
                .zip(&self.locals[i + 1])
                .any(|(p, q)| p != q);
        if mitred && changes {
            return Ok(Err(MITRED_SECTIONS));
        }
        let count = at_lo.len();
        let mut chord = vec![0.0 as Scalar; count];
        let mut slide = vec![0.0 as Scalar; count];
        for (length, bounds) in &pieces {
            let (f1, f2) = (rss(bounds, 0), rss(bounds, 1));
            let (lt, ut) = if plain {
                (bounds.lateral, bounds.up)
            } else {
                (turned(turn.lateral, [f1, f2]), turned(turn.up, [f1, f2]))
            };
            let frame = |j: usize| {
                base[0] * bounds.lateral[j] + base[1] * bounds.up[j] + base[2] * bounds.tangent[j]
            };
            let swing = base_rate[0] * bounds.lateral[0]
                + base_rate[1] * bounds.up[0]
                + base_rate[2] * bounds.tangent[0];
            let span = length + 2.0 * bounds.measure;
            for k in 0..count {
                let (p, q) = (self.locals[i][k], self.locals[i + 1][k]);
                let (x0, x1) = (lerp(p.x, q.x, u0), lerp(p.x, q.x, u1));
                let (y0, y1) = (lerp(p.y, q.y, u0), lerp(p.y, q.y, u1));
                let (x, y) = (big(x0, x1), big(y0, y1));
                let (xr, yr) = ((q.x - p.x).abs() / delta, (q.y - p.y).abs() / delta);
                let second = bounds.point[1]
                    + frame(1)
                    + 2.0 * swing
                    + x * lt[1]
                    + y * ut[1]
                    + 2.0 * (xr * lt[0] + yr * ut[0]);
                chord[k] = chord[k].max(0.125 * span * span * second);
                // A read point and axes off the exact law's (#316).
                let lever = if plain { 1.0 } else { Scalar::sqrt(3.0) };
                let read = bounds.point_error
                    + bounds.axis_error * (base.iter().sum::<Scalar>() + lever * (x + y));
                slide[k] = slide[k].max(read);
                if bounds.measure > 0.0 {
                    let speed = bounds.point[0]
                        + frame(0)
                        + base_rate.iter().sum::<Scalar>()
                        + x * lt[0]
                        + y * ut[0]
                        + xr
                        + yr;
                    let reach = x.hypot(y);
                    let rate = offset_rate
                        + (q - p).length() / delta
                        + reach * (turn.lateral[0] + turn.up[0]);
                    slide[k] = slide[k].max(bounds.measure * (speed + 2.0 * rate) + read);
                }
            }
        }
        // Breakpoints: seams inside the cell from both sides, and an end
        // on a seam from inside the cell, against the mesh's chord.
        let mut ends: Vec<(Scalar, SeamSide)> = Vec::new();
        if !mitred {
            if self.on_seam(lo) {
                ends.push((lo, SeamSide::Outgoing));
            }
            for &c in &inside {
                ends.push((c, SeamSide::Incoming));
                ends.push((c, SeamSide::Outgoing));
            }
            if self.on_seam(hi) {
                ends.push((hi, SeamSide::Incoming));
            }
        }
        let mut corner = vec![0.0 as Scalar; count];
        for (c, side) in ends {
            let exact = self.exact_at(c, side)?;
            let t = (c - lo) / (hi - lo);
            for k in 0..count {
                let on_chord = at_lo[k] + (at_hi[k] - at_lo[k]) * t;
                corner[k] = corner[k].max((exact[k] - on_chord).length());
            }
        }
        let mut worst: Scalar = 0.0;
        for &(k, l) in self.edges {
            let along = (chord[k] + corner[k]).max(chord[l] + corner[l]);
            let wall =
                axiolid_construct::loft::quad_deviation(at_lo[k], at_lo[l], at_hi[k], at_hi[l]);
            worst = worst.max(along + wall + slide[k].max(slide[l]));
        }
        let worst = worst * (1.0 + 1e-12);
        Ok(if worst.is_finite() {
            Ok(worst)
        } else {
            Err(UNCERTIFIED_FRAMES)
        })
    }
}
