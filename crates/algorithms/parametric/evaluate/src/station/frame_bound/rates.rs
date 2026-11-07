//! The rates of a section frame law (#406, ADR 0082 amendment).
//!
//! A section placed at the distance `s` stands in the frame `E = [T, L, U]`
//! (tangent, lateral, up) at the basis point `O`. In station distance the
//! frame turns by the Darboux form `E' = E Omega`:
//!
//! - `T' = beta_L L + beta_U U`
//! - `L' = -beta_L T + tau U`
//! - `U' = -beta_U T - tau L`
//! - `O' = sigma T + rise U`
//!
//! ([`FrameLaw::rates`]). `sigma = O'.T` is the forward rate, `beta_L = L.T'`
//! and `beta_U = U.T'` the bend toward the lateral and the up axis, `tau =
//! U.L'` the twist and `rise = O'.U` the basis point's rate along the up axis.
//! The Darboux vector in components is `(tau, -beta_U, beta_L)`.
//!
//! # Closed forms
//!
//! With `v = O'(s)` and `w = O''(s)` in the world, the Section reading has
//! `T = v / |v|`, so `rise = 0`, and every section frame is a rolled
//! reference-up frame: `L = cos(phi) n + sin(phi) u` and `U = -sin(phi) n +
//! cos(phi) u` for the horizontal left normal `n = Z x T / |T_xy|`, `u = T x
//! n` and the roll `phi` (zero off a banked curve). Then
//!
//! `tau = phi' + T'.(u x Z) / |T_xy|`
//!
//! with `T' = (w - (w.T) T) / |v|`: the unrolled frame's twist, which a level
//! curve has not (`u = Z`), plus the roll's rate. A reversed piece and a
//! placement keep the formula (the presented jet is used); a placement that
//! tilts `+Z` reads the reference-up frame of the moved tangent, so `phi' =
//! 0` there.
//!
//! The Plan reading has `T` the horizontal part of the tangent, `U = +Z`, so
//! `tau = beta_U = 0`, and `rise = O'.Z`.
//!
//! The roll's rate comes from the cant law (`crate::banked::roll_rate`).
//!
//! # Bounds over a span
//!
//! At the span's middle each rate is computed from the jet, within its
//! rounding, the frames' certified axis error and the jet's own error `ev`
//! of `v`, `ew` of `w` (zero in closed form; a frame read by quadrature
//! carries the error of its certified second derivative, see
//! `PieceReader::jet`). Over the span it moves by at most half the span
//! (plus the distance a read may lie off, the measure) times the rate's
//! Lipschitz constant, from the span's [`FrameLawBounds`]:
//!
//! - `sigma`: `sup |O''| + sup |O'| sup |T'|`;
//! - `beta_L`: `sup |L'| sup |T'| + sup |T''|`, and `beta_U` likewise;
//! - `rise`: `sup |O''|`;
//! - `tau`: `sup |U'| sup |L'| + sup |L''|`.

use axiolid_core::{Scalar, Vec3};
use axiolid_curve::{Curve3, SeamSide};

use super::bend::{closed_raw_jet, lipschitz, presented};
use super::{read_frame_bounds, FrameLaw, FrameReading, Slot, KEEPS_UP_TOLERANCE};
use crate::station::StationCurve;

/// One rate of a section frame law over a span (#406).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RateBounds {
    /// Certified `[lo, hi]` over the span.
    pub range: [Scalar; 2],
    /// The value at the span's middle read.
    pub middle: Scalar,
    /// How far `middle` may lie from the exact law's rate at the true
    /// distance that read stands at.
    pub error: Scalar,
}

impl RateBounds {
    /// A rate that is zero by the frame's construction.
    const ZERO: Self = Self {
        range: [0.0, 0.0],
        middle: 0.0,
        error: 0.0,
    };

    /// `middle +- (error + reach * lipschitz)`.
    fn around(middle: Scalar, error: Scalar, reach: Scalar, lipschitz: Scalar) -> Self {
        let pad = error + reach * lipschitz;
        Self {
            range: [middle - pad, middle + pad],
            middle,
            error,
        }
    }

    fn is_finite(&self) -> bool {
        self.range
            .iter()
            .chain([&self.middle, &self.error])
            .all(|x| x.is_finite())
    }
}

/// Certified bounds of a section frame law's rates over a span (#406):
/// `O' = sigma T + rise U`, `T' = beta_L L + beta_U U`, `L' = -beta_L T +
/// tau U`, in station distance; see the [module documentation](self).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameRates {
    /// `sigma = O'.T`, the forward rate.
    pub forward: RateBounds,
    /// `rise = O'.U`: exactly zero under [`FrameReading::Section`].
    pub rise: RateBounds,
    /// `beta_L = L.T'`, the bend toward the lateral axis.
    pub lateral: RateBounds,
    /// `beta_U = U.T'`, the bend toward the up axis: exactly zero under
    /// [`FrameReading::Plan`].
    pub up: RateBounds,
    /// `tau = U.L'`, the twist: exactly zero under [`FrameReading::Plan`].
    pub twist: RateBounds,
    /// The span's [`FrameLawBounds::measure`](super::FrameLawBounds::measure):
    /// a middle read's true distance lies within it of the distance asked.
    pub measure: Scalar,
}

/// The jet of the exact law at the span's middle in station distance, with
/// the certified error of each part.
pub(super) struct Jet {
    pub(super) v: Vec3,
    pub(super) w: Vec3,
    pub(super) velocity_error: Scalar,
    pub(super) acceleration_error: Scalar,
}

impl FrameLaw<'_, '_> {
    /// Certified bounds of the rates of the frame law over the distances
    /// `[a, b]` along the basis, a span that crosses no seam, read as
    /// `reading`; see the [module documentation](self).
    ///
    /// `None` wherever [`Self::bounds`] is `None`, for a chain, for a span
    /// the jet is not bounded over at this width, and for a frame read by
    /// quadrature whose span crosses a point where its curvature may jump.
    #[must_use]
    pub fn rates(&self, a: Scalar, b: Scalar, reading: FrameReading) -> Option<FrameRates> {
        let bounds = self.bounds(a, b, reading)?;
        let at = self.locate(a, b)?;
        let mid = 0.5 * (at.from + at.to);
        let (m, _) = at.piece.read_at(mid, SeamSide::Outgoing);
        let (jet, roll, measure, spanned) = match self.slots[at.slot].get()? {
            Slot::Closed => {
                let raw = closed_raw_jet(&at.piece, mid)?;
                let banked = match at.piece.curve {
                    StationCurve::Three(Curve3::Banked(banked)) => Some(banked),
                    _ => None,
                };
                (
                    Jet::exact(presented(&at.piece, raw)?),
                    roll_of(&at.piece, banked, m, raw),
                    bounds.measure,
                    bounds,
                )
            }
            Slot::Intrinsic(reader) => {
                let raw = reader.jet(m)?;
                (
                    Jet::exact(presented(&at.piece, raw)?),
                    roll_of(&at.piece, reader.banked(), m, raw),
                    bounds.measure,
                    bounds,
                )
            }
            Slot::Read(reader) => {
                // The exact law's true distance may lie `measure` off the
                // one read: its Lipschitz constants are taken over the span
                // widened by it, which must not cross a point where the
                // curvature may jump.
                let widening = bounds.measure;
                let from = (at.from - widening).max(0.0);
                let to = (at.to + widening).min(at.piece.length());
                let wide = read_frame_bounds(reader, &at.piece, from, to, at.table, reading)?;
                let (ta, tb) = reader
                    .span_parameters(&at.piece, from, to, at.table)
                    .ok()??;
                let (lo, hi) = (ta.min(tb), ta.max(tb));
                if reader.bend_breaks().iter().any(|&k| lo < k && k < hi) {
                    return None;
                }
                let jet = reader.jet(&at.piece, mid, at.table).ok()??;
                (
                    Jet {
                        velocity_error: jet.tangent_error,
                        acceleration_error: jet.curvature_error,
                        ..Jet::exact(presented(&at.piece, (jet.tangent, jet.curvature))?)
                    },
                    Some((0.0, 0.0)),
                    widening.max(wide.measure),
                    wide,
                )
            }
            Slot::Unbounded => return None,
        };
        let Jet {
            v,
            w,
            velocity_error: ev,
            acceleration_error: ew,
        } = jet;
        let mut frame = at.piece.section_on(mid, SeamSide::Outgoing).ok()?;
        let plan = reading == FrameReading::Plan;
        if plan {
            frame = frame.plan().ok()?;
        }
        let length = v.length();
        let speed = if plan { v.x.hypot(v.y) } else { length };
        if !(speed.is_finite() && speed > ev && ev >= 0.0 && ew >= 0.0) {
            return None;
        }
        let (rounding, axis) = (64.0 * Scalar::EPSILON, bounds.axis_error);
        let reach = 0.5 * (b - a) + measure;
        let lip = lipschitz(&spanned, plan);

        let forward_error = (rounding + 2.0 * axis) * length + ev;
        let floor = speed - ev;
        let turn_in = (ew + 2.0 * w.length() * ev / floor) / floor;
        let turn = (rounding + 3.0 * axis) * w.length() / floor + turn_in;
        let mut forward = RateBounds::around(v.dot(frame.tangent), forward_error, reach, lip[0]);
        forward.range[0] = forward.range[0].max(bounds.forward);
        let lateral = RateBounds::around(frame.lateral.dot(w) / speed, turn, reach, lip[1]);
        let up = if plan {
            RateBounds::ZERO
        } else {
            RateBounds::around(frame.up.dot(w) / speed, turn, reach, lip[2])
        };
        let rise = if plan {
            RateBounds::around(
                v.dot(frame.up),
                rounding * length + ev,
                reach,
                spanned.point[1],
            )
        } else {
            RateBounds::ZERO
        };
        let twist = if plan {
            RateBounds::ZERO
        } else {
            let (phi, phi_error) = roll?;
            let tn = v / length;
            let rho = tn.x.hypot(tn.y);
            let normal = Vec3::Z.cross(tn) / rho;
            let up0 = tn.cross(normal);
            let dt = (w - tn * w.dot(tn)) / length;
            let tangent_error = 2.0 * ev / (length - ev);
            let rho_floor = rho - tangent_error - rounding;
            if !(rho_floor.is_finite() && rho_floor > 0.0) {
                return None;
            }
            let middle = phi + dt.dot(up0.cross(Vec3::Z)) / rho;
            let error = rounding * (phi.abs() + 4.0 * dt.length() / rho_floor)
                + (turn_in + 4.0 * dt.length() * tangent_error) / rho_floor
                + phi_error;
            let slope = spanned.up[0] * spanned.lateral[0] + spanned.lateral[1];
            RateBounds::around(middle, error, reach, slope)
        };
        let out = FrameRates {
            forward,
            rise,
            lateral,
            up,
            twist,
            measure,
        };
        [out.forward, out.rise, out.lateral, out.up, out.twist]
            .iter()
            .all(RateBounds::is_finite)
            .then_some(out)
    }
}

impl Jet {
    /// A jet known in closed form.
    fn exact((v, w): (Vec3, Vec3)) -> Self {
        Self {
            v,
            w,
            velocity_error: 0.0,
            acceleration_error: 0.0,
        }
    }
}

/// The roll's rate and its error at the piece's measure `m` for the raw jet
/// `raw` of the curve: zero off a banked curve and under a placement that
/// tilts `+Z` (which reads the reference-up frame), `None` where a banked
/// curve's cant has no rate (a held-rail piece).
fn roll_of(
    piece: &crate::station::StationPiece<'_>,
    banked: Option<&axiolid_curve::Banked3>,
    m: Scalar,
    raw: (Vec3, Vec3),
) -> Option<(Scalar, Scalar)> {
    let Some(banked) = banked else {
        return Some((0.0, 0.0));
    };
    if let Some(rigid) = piece.placement {
        if (rigid.transform_vector3(Vec3::Z) - Vec3::Z).length() > KEEPS_UP_TOLERANCE {
            return Some((0.0, 0.0));
        }
    }
    crate::banked::roll_rate(banked, m, raw)
}
