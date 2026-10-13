//! Portable curve-evaluation provider contract.

use axiolid_contracts::{Backend, Determinism, GeomError, GeomResult, Operation};
use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::{Curve3, CurvePath, SeamSide};

use crate::{CurveMeasure, DistanceConvention, DistanceInterval};

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
///
/// # Curve paths (#290)
///
/// A curve relation a distance runs along -- a composite curve, a trim of
/// one, segments placed at stations of another curve -- reaches a provider
/// as a [`CurvePath`]: spans of atomic curves laid end to end, each read
/// forwards or backwards and carried by a rigid placement
/// (`axiolid-curve`'s neutral value; a graph compiler flattens a relation
/// into one). The `path_*` queries read it as ADR 0082 reads a composite
/// station basis (#285):
///
/// - The distance runs end to end, each piece in its own station measure
///   (plan distance on an elevated or banked curve, arc length on every
///   other); every piece must measure alike
///   ([`path_distance_convention`](Self::path_distance_convention)).
///   Consecutive pieces must meet; a path that does not is refused by
///   name.
/// - Every interior joint is a seam, never smooth: a distance within the
///   arc-length tolerance `1e-12 * max(1, s)` of a joint is ON it and is
///   read at the joint from the piece [`SeamSide`] names, the incoming one
///   at its end, the outgoing one at its start, each at its own point. A
///   piece's own seams strictly inside it are read by its curve's rule; at
///   a piece's ends it is read from inside. The plain queries read
///   [`SeamSide::Outgoing`].
/// - A reversed piece reads its curve at `end - u`, seam sides swapped,
///   tangent and left lateral negated, up kept; a placed piece's frame is
///   its curve's carried by the placement.
/// - The frame is the station section in this trait's layout (`x`
///   tangent, `y` up, `z` right); a frame is exact, rounding aside, only
///   where every piece up to and including the one read is a line placed,
///   if at all, exactly
///   ([`path_frame_is_exact_at`](Self::path_frame_is_exact_at)).
///
/// A piece may leave its end open ([`PathPiece::to_end`](axiolid_curve::PathPiece::to_end),
/// [`PathPiece::whole`](axiolid_curve::PathPiece::whole), #298): a caller
/// with no provider cannot state the length of an offset beside a gradient
/// curve, a clothoid or a B-spline. A provider that reads paths resolves
/// the end to where the curve ends by its own length of it (exact where
/// that is read from stored data, its numerical estimate otherwise), so
/// the distance runs through the pieces as it would through the stated
/// path holding those lengths, and
/// [`path_joint_distances`](Self::path_joint_distances) certifies where
/// the joints after it truly fall. A provider that reads stated paths but
/// not open ends refuses a path with one by a typed
/// [`GeomError::UnsupportedInput`] naming [`OPEN_PATH_END_UNSUPPORTED`],
/// never reading the `end` the piece leaves unset.
///
/// A piece may also start or end at a [`PathBound`](axiolid_curve::PathBound)
/// (#304): a native parameter of its curve, or a point on it in the
/// curve's own coordinates ([`PathPiece::between_parameters`](axiolid_curve::PathPiece::between_parameters),
/// [`PathPiece::bounded`](axiolid_curve::PathPiece::bounded)), since the
/// distance to an interior parameter of an ellipse or a B-spline is a
/// quadrature too. A provider that reads paths resolves a point to the
/// one parameter naming it, refusing by a typed
/// [`GeomError::UnsupportedInput`] a point off the curve
/// ([`PATH_POINT_OFF_CURVE`]) or one several parameters name
/// ([`PATH_POINT_AMBIGUOUS`]), and a parameter to its distance in the
/// curve's station measure: exact where that is read from stored data,
/// otherwise the same estimate it measures a span between parameters by,
/// so the path reads as the stated path holding those distances and
/// [`path_joint_distances`](Self::path_joint_distances) and
/// [`path_joint_distances_within`](Self::path_joint_distances_within)
/// certify the joints after it. A provider that does not read bounds refuses a path
/// with one naming [`PATH_PIECE_BOUND_UNSUPPORTED`], never reading the
/// `start` or `end` the piece leaves unset.
///
/// They have defaults so that adding them broke no provider: a provider
/// that does not implement paths refuses every path query by a typed
/// [`GeomError::UnsupportedInput`] naming [`CURVE_PATH_UNSUPPORTED`],
/// reports [`DistanceConvention::Unsupported`] and claims no frame exact.
/// It never reads a path as its first piece, or a joint from the wrong
/// side.
///
/// # Joint distances (#294)
///
/// A caller that holds a distance authored against the true geometry --
/// a station a model's precision puts on a joint -- needs to know where
/// the joints lie, and a provider's own joint distances are partly
/// numerical: the length of an ellipse or a B-spline piece, and of an
/// offset other than a line beside a line or a circle beside a circle, is
/// a quadrature. [`path_joint_distances`](Self::path_joint_distances)
/// answers each interior joint with a [`DistanceInterval`] that CONTAINS
/// its true distance and is no wider than the requested tolerance: the
/// sum of the true lengths of the pieces before it, each in the path's
/// convention (its own plan length for an offset beside a plan-measured
/// curve), between the points the provider reads as the piece's ends. A
/// piece whose length is read from stored data (a line, a circle, a
/// polyline, an intrinsic curve or a chain, an elevated or banked curve,
/// a constant offset of a line, or of a circle in planes normal to its
/// axis) contributes exactly, so a path of them has point intervals.
///
/// The default refuses by a typed [`GeomError::UnsupportedInput`] naming
/// [`PATH_JOINT_DISTANCES_UNSUPPORTED`]: a provider that does not bound its
/// pieces' lengths never hands out an estimate as a certificate.
///
/// # Joints near a window (#305)
///
/// A caller placing something at a distance `around` needs only the joints
/// that can fall within its precision `window` of it, and one piece whose
/// length cannot be certified anywhere on the path (an offset by distances
/// along a B-spline) refuses
/// [`path_joint_distances`](Self::path_joint_distances) for all of them.
/// [`path_joint_distances_within`](Self::path_joint_distances_within)
/// answers only the joints whose certified interval meets `[around -
/// window, around + window]`, each with the interval
/// [`path_joint_distances`](Self::path_joint_distances) would give it, so
/// it certifies only the pieces those joints lie after:
///
/// - **Before the window,** every piece is certified, since its length
///   moves every later joint; the intervals of the sum may widen with
///   each numerical piece, but stay within the tolerance. A piece before
///   the window whose length cannot be certified refuses the query by
///   name: the window's position against the joints after it is unknown.
///   No coarser bound is substituted; one loose enough to reach every
///   such piece would be far wider than any tolerance a placement asks.
/// - **Within the window,** each joint that can fall in it is certified
///   or the query refused by name: the uncertainty reaches the window.
/// - **After the window,** no piece is read: a joint is known to lie past
///   it once the joint before it does, or once a certified lower bound of
///   its piece's length (its chord) carries it past.
///
/// A joint left out therefore truly lies outside the window, and a joint
/// whose certified interval lies inside it is always returned. The
/// default delegates to [`path_joint_distances`](Self::path_joint_distances)
/// and keeps the joints whose intervals meet the window, which cannot be
/// wrong where that query is right, so a provider that does not certify
/// joints refuses it by the same name.
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

    /// Which distance this provider measures along `path` (see the
    /// [trait documentation](Self#curve-paths-290)): the convention every
    /// piece measures in, or [`DistanceConvention::Unsupported`] where the
    /// provider cannot measure the path (pieces measured differently, a
    /// joint whose pieces do not meet, a piece it does not know).
    ///
    /// Defaults to [`DistanceConvention::Unsupported`].
    fn path_distance_convention(&self, path: &CurvePath) -> DistanceConvention {
        let _ = path;
        DistanceConvention::Unsupported
    }

    /// Position at `at` along `path`, reading a joint from its outgoing
    /// piece: [`path_point_at_on`](Self::path_point_at_on) with
    /// [`SeamSide::Outgoing`].
    ///
    /// # Errors
    ///
    /// As [`path_point_at_on`](Self::path_point_at_on).
    fn path_point_at(&self, path: &CurvePath, at: CurveMeasure) -> GeomResult<Point3> {
        self.path_point_at_on(path, at, SeamSide::Outgoing)
    }

    /// Unit tangent at `at` along `path`, reading a joint from its outgoing
    /// piece: [`path_tangent_at_on`](Self::path_tangent_at_on) with
    /// [`SeamSide::Outgoing`].
    ///
    /// # Errors
    ///
    /// As [`path_tangent_at_on`](Self::path_tangent_at_on).
    fn path_tangent_at(&self, path: &CurvePath, at: CurveMeasure) -> GeomResult<Vec3> {
        self.path_tangent_at_on(path, at, SeamSide::Outgoing)
    }

    /// Oriented frame at `at` along `path`, reading a joint from its
    /// outgoing piece: [`path_frame_at_on`](Self::path_frame_at_on) with
    /// [`SeamSide::Outgoing`].
    ///
    /// # Errors
    ///
    /// As [`path_frame_at_on`](Self::path_frame_at_on).
    fn path_frame_at(&self, path: &CurvePath, at: CurveMeasure) -> GeomResult<Frame3> {
        self.path_frame_at_on(path, at, SeamSide::Outgoing)
    }

    /// Position at `at` along `path`, reading `side` of a joint or a seam
    /// it lies on (see the [trait documentation](Self#curve-paths-290)).
    ///
    /// # Errors
    ///
    /// A measure that is not a distance in the path's convention, not
    /// finite, or outside `[0, L]`; a path the provider cannot measure; the
    /// pieces' own refusals. By default, everything by
    /// [`GeomError::UnsupportedInput`] naming [`CURVE_PATH_UNSUPPORTED`].
    fn path_point_at_on(
        &self,
        path: &CurvePath,
        at: CurveMeasure,
        side: SeamSide,
    ) -> GeomResult<Point3> {
        let _ = (path, at, side);
        Err(path_unsupported(self))
    }

    /// Unit tangent at `at` along `path`, reading `side` of a joint or a
    /// seam it lies on (see the [trait documentation](Self#curve-paths-290)).
    ///
    /// # Errors
    ///
    /// As [`path_point_at_on`](Self::path_point_at_on).
    fn path_tangent_at_on(
        &self,
        path: &CurvePath,
        at: CurveMeasure,
        side: SeamSide,
    ) -> GeomResult<Vec3> {
        let _ = (path, at, side);
        Err(path_unsupported(self))
    }

    /// Oriented frame at `at` along `path`, reading `side` of a joint or a
    /// seam it lies on (see the [trait documentation](Self#curve-paths-290)),
    /// in [`frame_at`](Self::frame_at)'s layout: `x` the side's tangent,
    /// `y` up, `z` to the right; on a banked piece, its rolled section.
    ///
    /// # Errors
    ///
    /// As [`path_point_at_on`](Self::path_point_at_on), and a tangent
    /// parallel to the reference up.
    fn path_frame_at_on(
        &self,
        path: &CurvePath,
        at: CurveMeasure,
        side: SeamSide,
    ) -> GeomResult<Frame3> {
        let _ = (path, at, side);
        Err(path_unsupported(self))
    }

    /// Whether the frame at `at` along `path`, read from `side`, is exact,
    /// rounding aside: only where every piece up to and including the one
    /// read is a line placed, if at all, exactly (see the
    /// [trait documentation](Self#curve-paths-290)). `false` wherever the
    /// provider would refuse the query.
    ///
    /// Defaults to `false`: a provider that does not implement paths
    /// claims nothing exact.
    fn path_frame_is_exact_at(&self, path: &CurvePath, at: CurveMeasure, side: SeamSide) -> bool {
        let _ = (path, at, side);
        false
    }

    /// Where the interior joints of `path` lie, each certified to a
    /// [`DistanceInterval`] no wider than `tolerance` (see the
    /// [trait documentation](Self#joint-distances-294)): one interval per
    /// joint, in order, `pieces - 1` in all.
    ///
    /// # Errors
    ///
    /// A `tolerance` that is not finite and positive, a path the provider
    /// cannot measure (as [`path_point_at_on`](Self::path_point_at_on)),
    /// and a piece whose length it cannot certify to `tolerance`, by name.
    /// By default, everything by [`GeomError::UnsupportedInput`] naming
    /// [`PATH_JOINT_DISTANCES_UNSUPPORTED`].
    fn path_joint_distances(
        &self,
        path: &CurvePath,
        tolerance: Scalar,
    ) -> GeomResult<Vec<DistanceInterval>> {
        let _ = (path, tolerance);
        Err(GeomError::UnsupportedInput {
            backend: self.descriptor().id,
            operation: Operation::CurveEvaluation,
            input: PATH_JOINT_DISTANCES_UNSUPPORTED,
        })
    }

    /// The interior joints of `path` that can fall within `window` of
    /// `around` (see the [trait documentation](Self#joints-near-a-window-305)):
    /// each as `(joint, interval)`, `joint` its index among the path's
    /// joints (the joint where piece `joint` ends and piece `joint + 1`
    /// starts, its position in
    /// [`path_joint_distances`](Self::path_joint_distances)' answer), in
    /// order, `interval` certified as that query certifies it, no wider
    /// than `tolerance`. Returned: every joint whose interval meets
    /// `[around - window, around + window]` (the edges included). Left
    /// out: only joints certified to lie outside it. Empty when no joint
    /// can fall within it.
    ///
    /// # Errors
    ///
    /// An `around` that is not finite or a `window` that is negative or not
    /// finite ([`GeomError::InvalidInput`]); as
    /// [`path_joint_distances`](Self::path_joint_distances) otherwise, but
    /// only for the pieces before and within the window. By default,
    /// whatever [`path_joint_distances`](Self::path_joint_distances)
    /// refuses, so [`PATH_JOINT_DISTANCES_UNSUPPORTED`] where that is not
    /// implemented.
    fn path_joint_distances_within(
        &self,
        path: &CurvePath,
        around: Scalar,
        window: Scalar,
        tolerance: Scalar,
    ) -> GeomResult<Vec<(usize, DistanceInterval)>> {
        let span = window_span(around, window)?;
        Ok(self
            .path_joint_distances(path, tolerance)?
            .into_iter()
            .enumerate()
            .filter(|(_, interval)| interval.meets(&span))
            .collect())
    }
}

/// The window `[around - window, around + window]` of
/// [`path_joint_distances_within`](CurveEvaluator::path_joint_distances_within),
/// checked: `around` finite, `window` finite and not negative.
///
/// # Errors
///
/// Any other, by [`GeomError::InvalidInput`].
pub fn window_span(around: Scalar, window: Scalar) -> GeomResult<DistanceInterval> {
    if !around.is_finite() {
        return Err(GeomError::InvalidInput(format!(
            "joints near a window: the window's centre {around} is not finite"
        )));
    }
    if !(window.is_finite() && window >= 0.0) {
        return Err(GeomError::InvalidInput(format!(
            "joints near a window: the window {window} is negative or not finite"
        )));
    }
    Ok(DistanceInterval::new(around - window, around + window))
}

/// The input a provider that does not certify joint distances refuses, in
/// the [`GeomError::UnsupportedInput`] its default
/// [`path_joint_distances`](CurveEvaluator::path_joint_distances) returns.
pub const PATH_JOINT_DISTANCES_UNSUPPORTED: &str =
    "certified joint distances of a curve path: this provider does not bound the lengths of a \
     path's pieces";

/// The input a provider that reads curve paths but not open piece ends
/// (#298) refuses a path with one for.
pub const OPEN_PATH_END_UNSUPPORTED: &str =
    "a curve path piece whose end is open (runs to its curve's end): this provider reads stated \
     spans only";

/// The input a provider that reads curve paths but not piece bounds given
/// as a native parameter or a point (#304) refuses a path with one for.
pub const PATH_PIECE_BOUND_UNSUPPORTED: &str =
    "a curve path piece bounded by a native parameter or a point: this provider reads spans \
     stated as distances only";

/// The input a point bound of a curve path piece (#304) is refused for
/// when it does not lie on the piece's curve within the joint tolerance
/// `1e-9 * max(1, |p|)`: a provider never projects it onto the curve.
pub const PATH_POINT_OFF_CURVE: &str =
    "a curve path piece bounded by a point that does not lie on its curve within the joint \
     tolerance";

/// The input a point bound of a curve path piece (#304) is refused for
/// when several parameters of the curve name it: a curve that passes
/// through it more than once, or the seam of a closed curve other than a
/// circle or an ellipse (whose angle runs round).
pub const PATH_POINT_AMBIGUOUS: &str =
    "a curve path piece bounded by a point that several parameters of its curve name (the curve \
     passes through it more than once, or it is the seam of a closed curve other than a circle \
     or an ellipse)";

/// The input a provider that does not implement curve paths refuses, in
/// the [`GeomError::UnsupportedInput`] its default `path_*` queries return.
pub const CURVE_PATH_UNSUPPORTED: &str =
    "a curve path (pieces of curves laid end to end): this provider evaluates single curves only";

/// The default refusal of a curve path, naming `provider`.
fn path_unsupported<E: Backend + ?Sized>(provider: &E) -> GeomError {
    GeomError::UnsupportedInput {
        backend: provider.descriptor().id,
        operation: Operation::CurveEvaluation,
        input: CURVE_PATH_UNSUPPORTED,
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
