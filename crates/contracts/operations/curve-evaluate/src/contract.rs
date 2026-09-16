//! Portable curve-evaluation provider contract.

use axiolid_contracts::{Backend, Determinism, GeomResult};
use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::Curve3;

use crate::DistanceConvention;

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

    /// Position at `distance` along `curve`.
    ///
    /// `distance` is measured in this provider's
    /// [`distance_convention`](Self::distance_convention) for this curve.
    /// Refuses when the convention is
    /// [`Unsupported`](DistanceConvention::Unsupported), when the distance
    /// is not finite, or when it falls outside the curve.
    fn point_at(&self, curve: &Curve3, distance: Scalar) -> GeomResult<Point3>;

    /// Unit tangent at `distance` along `curve`.
    ///
    /// Unit length is part of the contract: a caller composing a rotation
    /// from this must not have to renormalise, and a non-unit result
    /// would silently scale whatever it is applied to.
    fn tangent_at(&self, curve: &Curve3, distance: Scalar) -> GeomResult<Vec3>;

    /// Oriented frame at `distance` along `curve`.
    ///
    /// `x` is the unit tangent, `z` is the reference-up direction, and
    /// `y = z x x` completes a right-handed orthonormal triad. See the
    /// trait docs for why this is not the Frenet frame.
    ///
    /// Refuses when the tangent is parallel to `up`, where roll is
    /// genuinely undetermined.
    fn frame_at(&self, curve: &Curve3, distance: Scalar) -> GeomResult<Frame3>;
}
