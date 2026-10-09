//! Portable curve-evaluation provider contract.

use axiolid_contracts::{Backend, Determinism, GeomError, GeomResult, Operation};
use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_curve::{Curve3, SeamSide};

use crate::{CurveMeasure, DistanceConvention};

/// Curve evaluation provider.
///
/// Implementing this trait is the capability declaration: it lets a
/// consumer name "evaluate a curve at a distance" without depending on a
/// particular engine. A provider that cannot evaluate curves must not
/// implement it.
///
/// # Why a frame and not just a tangent
///
/// A tangent fixes direction but not roll, so placing an object needs a
/// full oriented frame. The convention has to be decided ONCE, here,
/// because a wrong one tilts the placed object rather than failing
/// loudly -- and every consumer re-deriving it would fragment slightly.
///
/// # Why the frame is not Frenet
///
/// The Frenet normal is the wrong tool for placement even though it is
/// the classical one. On an alignment with a crest followed by a sag it
/// points DOWN on the crest and UP in the sag, flipping at the
/// inflection, and on any straight run the curvature is zero and the
/// normal is undefined entirely. An object placed by it would be upright,
/// then inverted, then unplaceable.
///
/// [`frame_at`](Self::frame_at) therefore uses a REFERENCE-UP frame:
/// `right = normalise(tangent x up)` and `up' = right x tangent`, which
/// is continuous through inflections and defined on straights. It is
/// undefined only when the tangent is parallel to the reference
/// direction -- a truly vertical curve -- where a provider must refuse
/// rather than return an arbitrary roll.
///
/// # Seams (#286)
///
/// A curve made of pieces -- a polyline's vertex, an elevated curve's
/// grade break, a banked curve's cant jump -- has two tangents or two
/// rolls where its pieces meet. [`point_at`](Self::point_at),
/// [`tangent_at`](Self::tangent_at) and [`frame_at`](Self::frame_at) read
/// the piece the measure falls in, which on a seam is the piece that
/// STARTS there. [`point_at_on`](Self::point_at_on),
/// [`tangent_at_on`](Self::tangent_at_on) and
/// [`frame_at_on`](Self::frame_at_on) name the side instead, by the
/// station seam rule of ADR 0082 (#263): a measure within the arc-length
/// tolerance of a seam whose frame may jump is ON it, read at the seam
/// from the piece [`SeamSide`] names -- [`SeamSide::Outgoing`], the one
/// starting there, or [`SeamSide::Incoming`], the one ending there. Off a
/// seam, and at the curve's own start and end, both sides read the same
/// piece and agree.
///
/// The sided methods have defaults so that adding them broke no provider:
/// `Outgoing` delegates to the side-less method, and every other side is
/// refused by a typed [`GeomError::UnsupportedInput`] naming
/// [`SEAM_SIDE_UNSUPPORTED`]. A provider that does not implement sides
/// never answers `Incoming` with the outgoing frame.
pub trait CurveEvaluator: Backend {
    /// Which distance this provider measures for `curve`.
    ///
    /// Callers MUST consult this before trusting a distance. Defaults to
    /// [`DistanceConvention::Unsupported`] so a provider that has not
    /// declared a convention is treated as unable rather than assumed to
    /// mean 3D arc length.
    fn distance_convention(&self, curve: &Curve3) -> DistanceConvention {
        let _ = curve;
        DistanceConvention::Unsupported
    }

    /// Reproducibility this provider guarantees.
    ///
    /// Defaults to [`Determinism::BestEffort`], the weakest level, so an
    /// unaudited provider cannot silently satisfy a stronger request.
    fn determinism(&self) -> Determinism {
        Determinism::BestEffort
    }

    /// Position at `at` along `curve`.
    ///
    /// [`CurveMeasure`] carries WHICH method of measurement the caller
    /// means, so a native parameter cannot be mistaken for a length. A
    /// [`Distance`](CurveMeasure::Distance) is interpreted in this
    /// provider's [`distance_convention`](Self::distance_convention) for
    /// this curve; a [`Parameter`](CurveMeasure::Parameter) is the curve's
    /// own parameter and is independent of that convention.
    ///
    /// Refuses a non-finite value, a distance on a curve whose convention
    /// is [`Unsupported`](DistanceConvention::Unsupported), and a value
    /// outside the curve.
    fn point_at(&self, curve: &Curve3, at: CurveMeasure) -> GeomResult<Point3>;

    /// Unit tangent at `at` along `curve`.
    ///
    /// Unit length is part of the contract: a caller composing a rotation
    /// from this must not have to renormalise, and a non-unit result
    /// would silently scale whatever it is applied to.
    fn tangent_at(&self, curve: &Curve3, at: CurveMeasure) -> GeomResult<Vec3>;

    /// Oriented frame at `at` along `curve`.
    ///
    /// `x` is the unit tangent, `y` is the reference up made perpendicular
    /// to it (`y = z x x`, leaning back with the grade), and
    /// `z = normalise(x x up)` points to the RIGHT of the tangent: a
    /// right-handed orthonormal triad whose second axis is up (ADR 0063;
    /// #242 corrected this text, which used to call `z` up). On a level
    /// line along `+X` with `up = +Z` the frame is `x = +X`, `y = +Z`,
    /// `z = -Y`. See the trait docs for why this is not the Frenet frame.
    ///
    /// Refuses when the tangent is parallel to `up`, where roll is
    /// genuinely undetermined.
    ///
    /// A curve that carries its own roll ([`Curve3::Banked`]) is framed by
    /// it: the section frame rolled about the tangent by the curve's bank
    /// convention -- `x` the tangent, `y` the section up, `z` minus its
    /// left lateral -- which at zero cant is the reference-up frame with
    /// `up = +Z` (ADR 0081).
    fn frame_at(&self, curve: &Curve3, at: CurveMeasure) -> GeomResult<Frame3>;

    /// [`point_at`](Self::point_at) reading `side` of a seam `at` lies on
    /// (see the [trait documentation](Self#seams-286)).
    ///
    /// On a seam both sides share the point up to the curve's own
    /// continuity; a provider that applies the seam tolerance reads the
    /// point at the seam, not at `at`.
    ///
    /// # Errors
    ///
    /// As [`point_at`](Self::point_at). By default, every side but
    /// [`SeamSide::Outgoing`] by [`GeomError::UnsupportedInput`] naming
    /// [`SEAM_SIDE_UNSUPPORTED`].
    fn point_at_on(&self, curve: &Curve3, at: CurveMeasure, side: SeamSide) -> GeomResult<Point3> {
        match side {
            SeamSide::Outgoing => self.point_at(curve, at),
            _ => Err(side_unsupported(self)),
        }
    }

    /// [`tangent_at`](Self::tangent_at) reading `side` of a seam `at` lies
    /// on (see the [trait documentation](Self#seams-286)): at a polyline's
    /// vertex, the direction of the segment ending there (`Incoming`) or
    /// starting there (`Outgoing`).
    ///
    /// # Errors
    ///
    /// As [`tangent_at`](Self::tangent_at). By default, every side but
    /// [`SeamSide::Outgoing`] by [`GeomError::UnsupportedInput`] naming
    /// [`SEAM_SIDE_UNSUPPORTED`].
    fn tangent_at_on(&self, curve: &Curve3, at: CurveMeasure, side: SeamSide) -> GeomResult<Vec3> {
        match side {
            SeamSide::Outgoing => self.tangent_at(curve, at),
            _ => Err(side_unsupported(self)),
        }
    }

    /// [`frame_at`](Self::frame_at) reading `side` of a seam `at` lies on
    /// (see the [trait documentation](Self#seams-286)), in the same
    /// layout: `x` the side's tangent, `y` up, `z` to the right; on a
    /// banked curve's cant jump, the side's roll.
    ///
    /// # Errors
    ///
    /// As [`frame_at`](Self::frame_at). By default, every side but
    /// [`SeamSide::Outgoing`] by [`GeomError::UnsupportedInput`] naming
    /// [`SEAM_SIDE_UNSUPPORTED`].
    fn frame_at_on(&self, curve: &Curve3, at: CurveMeasure, side: SeamSide) -> GeomResult<Frame3> {
        match side {
            SeamSide::Outgoing => self.frame_at(curve, at),
            _ => Err(side_unsupported(self)),
        }
    }
}

/// The input a provider that does not implement seam sides refuses, in
/// the [`GeomError::UnsupportedInput`] its default sided methods return.
pub const SEAM_SIDE_UNSUPPORTED: &str =
    "a seam side other than Outgoing: this provider reads only the piece that starts at a seam";

/// The default refusal of a seam side, naming `provider`.
fn side_unsupported<E: Backend + ?Sized>(provider: &E) -> GeomError {
    GeomError::UnsupportedInput {
        backend: provider.descriptor().id,
        operation: Operation::CurveEvaluation,
        input: SEAM_SIDE_UNSUPPORTED,
    }
}
