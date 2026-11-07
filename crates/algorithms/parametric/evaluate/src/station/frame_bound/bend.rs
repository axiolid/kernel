//! The signed bend of a section frame law (#361, ADR 0082 amendment).
//!
//! A section placed at the distance `s` lies in the plane through its
//! point normal to the frame's tangent `N`. A point `A = O + a L + b U + c
//! T` of the section, with the offsets `a`, `b`, `c` and the frame `(T, L,
//! U)` read at `s`, moves along `N` at
//!
//! `F = O'.N + c' - a L.N' - b U.N'`
//!
//! (`T.T' = 0`, and `L.N = U.N = 0`, so `L' . N = -L . N'`): the section
//! planes advance past the point while `F` is positive, and stop where it
//! is not, the point on the inner side of the bend at or beyond its centre
//! of curvature. [`FrameLaw::bend`] bounds the three coefficients
//!
//! - `sigma = O'.N`, the forward rate of the basis point;
//! - `beta_L = L.N'` and `beta_U = U.N'`, how fast the frame's tangent
//!   turns toward the lateral and the up axis,
//!
//! over a span, and a caller decides `F` for the points of its sections.
//! They are affine in the point's coordinates, so over a section the
//! smallest `F` stands at a vertex. A section turned about the frame needs
//! the frame's twist and rise too: [`FrameLaw::rates`] (#406), of which
//! [`FrameLaw::bend`] is the projection.
//!
//! # Closed forms
//!
//! With `v = O'(s)` and `w = O''(s)` in the world, per unit of station
//! distance, the Section reading has `N = v / |v|`, so `sigma = v . N`,
//! `beta_L = L . w / |v|` and `beta_U = U . w / |v|`. The plan reading has
//! `N` the horizontal part of the tangent, `U = +Z` and `N' `horizontal, so
//! `beta_U = 0` and `beta_L = L_p . w / |v_h|`. On a plan circle `beta_L`
//! is the curvature, on a gradient curve `beta_U` is `z'' / |v|^2`.
//!
//! The jet `(v, w)` is closed form where the point is: a line, a circle
//! with an orthonormal frame, a polyline's segment, an elevated or a banked
//! curve over a line or a circle, and a piece of an intrinsic curve (a
//! banked curve over an intrinsic plan included, #396), whose heading and
//! curvature are closed form (`IntrinsicReader::jet`).
//! A reversed piece negates `v` and keeps `w`; a placement maps both.
//!
//! # Bounds over a span
//!
//! At the span's middle each coefficient is computed from the jet and the
//! frame the kernel reads, within `(64 eps + 2 e) |v|` for `sigma` and
//! `(64 eps + 3 e) |w| / |v|` for the other two, `e` the frames' certified
//! axis error. Over the span it moves by at most half the span (plus the
//! distance a read may lie off) times the coefficient's Lipschitz
//! constant, from the span's [`FrameLawBounds`]:
//!
//! - `sigma`: `sup |O''| + sup |O'| sup |N'|`;
//! - `beta_L`: `sup |L'| sup |N'| + sup |N''|`, and `beta_U` likewise.
//!
//! Frames read by quadrature (an ellipse, a B-spline, an offset) and
//! every frame without certified bounds have no bend.

use axiolid_core::{Scalar, Vec3};
use axiolid_curve::{Curve2, Curve3, SeamSide};

use super::{
    closed_form_elevated, lift, orthonormal, polyline_segment, unit, FrameLaw, FrameLawBounds,
    FrameReading, Slot,
};
use crate::station::{StationCurve, StationPiece};

/// Certified bounds, over a span of station distance, of the coefficients
/// of the signed bend of a section frame law; see the
/// [module documentation](super#the-signed-bend-361).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BendBounds {
    /// Certified `[lo, hi]` over the span of `O'.N`.
    pub forward: [Scalar; 2],
    /// Certified `[lo, hi]` over the span of `L.N'`.
    pub lateral: [Scalar; 2],
    /// Certified `[lo, hi]` over the span of `U.N'`.
    pub up: [Scalar; 2],
    /// `[forward, lateral, up]` at the span's middle, each within
    /// `middle_error`.
    pub middle: [Scalar; 3],
    /// How far each `middle` value may lie from the coefficient there.
    pub middle_error: [Scalar; 3],
}

/// The jet `(O', O'')` of a closed-form piece at distance `u` into it, in
/// station distance, as the piece presents it (reversed, placed); `None`
/// for a curve whose point is not closed form.
pub(crate) fn closed_jet(piece: &StationPiece<'_>, u: Scalar) -> Option<(Vec3, Vec3)> {
    presented(piece, closed_raw_jet(piece, u)?)
}

/// [`closed_jet`] in the curve's own measure, before the piece reverses or
/// places it.
pub(super) fn closed_raw_jet(piece: &StationPiece<'_>, u: Scalar) -> Option<(Vec3, Vec3)> {
    let (m, _) = piece.read_at(u, SeamSide::Outgoing);
    match piece.curve {
        StationCurve::Two(curve) => jet2(curve, m),
        StationCurve::Three(curve) => jet3(curve, m),
        StationCurve::Offset(_) => None,
    }
}

/// A curve's jet in its own measure as `piece` presents it: negated
/// velocity when reversed, then placed.
pub(crate) fn presented(piece: &StationPiece<'_>, (v, w): (Vec3, Vec3)) -> Option<(Vec3, Vec3)> {
    let (v, w) = if piece.reversed { (-v, w) } else { (v, w) };
    let (v, w) = match piece.placement {
        Some(rigid) => (rigid.transform_vector3(v), rigid.transform_vector3(w)),
        None => (v, w),
    };
    (v.is_finite() && w.is_finite()).then_some((v, w))
}

fn jet2(curve: &Curve2, m: Scalar) -> Option<(Vec3, Vec3)> {
    match curve {
        Curve2::Line(line) => Some((unit(lift(line.direction))?, Vec3::ZERO)),
        Curve2::Circle(c) => {
            if !orthonormal(lift(c.frame.x), lift(c.frame.y)) {
                return None;
            }
            let t = m / c.radius;
            let v = lift(crate::curve::derivative2(curve, t).ok()?) / c.radius;
            let w = lift(crate::curve::second_derivative2(curve, t).ok()?) / (c.radius * c.radius);
            Some((v, w))
        }
        Curve2::Polyline(polyline) => {
            let step =
                polyline_segment(&polyline.points, polyline.closed, m, m, |a, b| lift(b - a))?;
            Some((unit(step)?, Vec3::ZERO))
        }
        _ => None,
    }
}

fn jet3(curve: &Curve3, m: Scalar) -> Option<(Vec3, Vec3)> {
    match curve {
        Curve3::Line(line) => Some((unit(line.direction)?, Vec3::ZERO)),
        Curve3::Circle(c) => {
            if !orthonormal(c.frame.x, c.frame.y) {
                return None;
            }
            let t = m / c.radius;
            let v = crate::curve::derivative3(curve, t).ok()? / c.radius;
            let w = crate::curve::second_derivative3(curve, t).ok()? / (c.radius * c.radius);
            Some((v, w))
        }
        Curve3::Polyline(polyline) => {
            let step = polyline_segment(&polyline.points, polyline.closed, m, m, |a, b| b - a)?;
            Some((unit(step)?, Vec3::ZERO))
        }
        Curve3::Elevated(elevated) => {
            if !closed_form_elevated(elevated) {
                return None;
            }
            Some((
                crate::elevated::elevated_derivative(elevated, m).ok()?,
                crate::elevated::elevated_second_derivative(elevated, m).ok()?,
            ))
        }
        Curve3::Banked(banked) => {
            if !closed_form_elevated(&banked.base) {
                return None;
            }
            Some((
                crate::banked::banked_derivative(banked, m).ok()?,
                crate::banked::banked_second_derivative(banked, m).ok()?,
            ))
        }
        _ => None,
    }
}

/// The Lipschitz constants, in station distance, of `sigma`, `beta_L` and
/// `beta_U` from the span's frame bounds (`beta_U` is zero in the plan
/// reading, where `U = +Z` and `N'` is horizontal).
pub(super) fn lipschitz(bounds: &FrameLawBounds, plan: bool) -> [Scalar; 3] {
    let beta = |axis: [Scalar; 2]| axis[0] * bounds.tangent[0] + bounds.tangent[1];
    [
        bounds.point[1] + bounds.point[0] * bounds.tangent[0],
        beta(bounds.lateral),
        if plan { 0.0 } else { beta(bounds.up) },
    ]
}

impl FrameLaw<'_, '_> {
    /// Whether the frame over the span `[a, b]` has a closed-form jet, so
    /// that [`Self::bend`] may bound it over narrower spans: `None` from
    /// `bend` can then be "not bounded at this width" (a span longer than
    /// twice a bend's radius), which halving the span usually decides. It
    /// can stay `None` however narrow the span (a vertical tangent under the
    /// plan reading, joints after a piece read by quadrature that cannot be
    /// certified, a placement whose columns are not orthogonal). `false` for
    /// frames read by quadrature, every frame without a closed form and a
    /// piece placed at a station whose frame is not exact (unless the law
    /// reads placements as exact, [`FrameLaw::with_placements_read_as_exact`]).
    #[must_use]
    pub fn has_bend(&self, a: Scalar, b: Scalar) -> bool {
        let Some(at) = self.locate(a, b) else {
            return false;
        };
        // A piece placed inexactly bounds nothing (#408), whatever its jet.
        if !at.piece.placement_exact && !self.placements_read_as_exact {
            return false;
        }
        // Builds the piece's slot, whatever the span's bounds are.
        let _ = self.read(
            at.slot,
            &at.piece,
            at.from,
            at.to,
            at.table,
            FrameReading::Section,
        );
        let mid = 0.5 * (at.from + at.to);
        match self.slots[at.slot].get() {
            Some(Slot::Closed) => closed_jet(&at.piece, mid).is_some(),
            Some(Slot::Intrinsic(reader)) => {
                let (m, _) = at.piece.read_at(mid, SeamSide::Outgoing);
                reader.jet(m).is_some()
            }
            _ => false,
        }
    }

    /// Certified bounds of the coefficients of the signed bend over the
    /// distances `[a, b]` along the basis, a span that crosses no seam,
    /// read as `reading`; see the [module documentation](super#the-signed-bend-361).
    /// A projection of [`Self::rates`].
    ///
    /// `None` wherever [`Self::rates`] is: a chain, a span the jet is not
    /// bounded over, and a frame read by quadrature whose span crosses a
    /// point where its curvature may jump.
    #[must_use]
    pub fn bend(&self, a: Scalar, b: Scalar, reading: FrameReading) -> Option<BendBounds> {
        let rates = self.rates(a, b, reading)?;
        let [forward, lateral, up] = [rates.forward, rates.lateral, rates.up];
        Some(BendBounds {
            forward: forward.range,
            lateral: lateral.range,
            up: up.range,
            middle: [forward.middle, lateral.middle, up.middle],
            middle_error: [forward.error, lateral.error, up.error],
        })
    }
}
