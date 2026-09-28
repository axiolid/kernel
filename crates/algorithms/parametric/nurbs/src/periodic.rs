//! Verified periodic parameter semantics for closed NURBS curves.
//!
//! A `closed` flag is metadata, not proof of periodicity. Wrapping a
//! parameter goes through a view whose constructor has checked the seam
//! (at least positional continuity), never through the flag alone.

use crate::{
    axis::active_spans,
    transform::{insert_knot2, insert_knot3, split2, split3},
};
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{BSplineCurve2, BSplineCurve3};
use axiolid_evaluate::curve::{bspline_jet2, bspline_jet3, CurveJet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
/// Highest endpoint continuity verified in the curve's native parameter.
pub enum SeamContinuity {
    /// Endpoint positions differ beyond tolerance.
    Discontinuous,
    /// Endpoint positions agree, but first derivatives do not.
    Position,
    /// Positions and first derivatives agree, but second derivatives do not.
    FirstDerivative,
    /// Positions and first two derivatives agree.
    SecondDerivative,
}

/// Verified closed-seam extension of a planar B-spline curve.
///
/// This view does not reinterpret the curve as an algebraically periodic knot
/// vector. It validates the declared geometric seam and wraps only through the
/// explicit methods below. The borrowed neutral curve remains unchanged.
#[derive(Debug, Clone, Copy)]
pub struct PeriodicCurve2<'a> {
    curve: &'a BSplineCurve2,
    tolerance: Tolerance,
    domain: (Scalar, Scalar),
    continuity: SeamContinuity,
}

impl<'a> PeriodicCurve2<'a> {
    /// Validate `closed` metadata and at least positional seam continuity.
    pub fn try_new(curve: &'a BSplineCurve2, tolerance: Tolerance) -> GeomResult<Self> {
        let domain = domain(
            curve.degree,
            curve.control_points.len(),
            &curve.knots,
            &curve.multiplicities,
        )?;
        let continuity = curve2_seam_continuity(curve, tolerance)?;
        require_verified_seam(curve.closed, continuity)?;
        Ok(Self {
            curve,
            tolerance,
            domain,
            continuity,
        })
    }

    /// Native active parameter domain retained by this view.
    pub const fn domain(self) -> (Scalar, Scalar) {
        self.domain
    }

    /// Highest endpoint continuity verified during construction.
    pub const fn seam_continuity(self) -> SeamContinuity {
        self.continuity
    }

    /// Wrap a finite parameter into the native active domain.
    ///
    /// An in-domain upper endpoint remains the upper endpoint so one-sided
    /// endpoint jets retain the neutral evaluator's existing semantics.
    pub fn wrap_parameter(self, parameter: Scalar) -> GeomResult<Scalar> {
        wrap(parameter, self.domain)
    }

    /// Evaluate the closed-seam extension at any finite parameter.
    pub fn evaluate(self, parameter: Scalar) -> GeomResult<Point2> {
        Ok(self.jet(parameter)?.point)
    }

    /// Evaluate point, first derivative, and second derivative after wrapping.
    pub fn jet(self, parameter: Scalar) -> GeomResult<CurveJet<Point2, Vec2>> {
        bspline_jet2(self.curve, self.wrap_parameter(parameter)?)
    }

    /// Insert a shape-preserving knot at a periodic-equivalent interior parameter.
    ///
    /// Exterior parameters are canonicalized. Parameters equivalent to either
    /// seam endpoint are rejected because clamped endpoint multiplicity cannot
    /// be increased as an interior edit.
    pub fn insert_knot(self, parameter: Scalar) -> GeomResult<BSplineCurve2> {
        let native = interior_periodic_parameter(parameter, self.domain)?;
        let edited = insert_knot2(self.curve, native)?;
        PeriodicCurve2::try_new(&edited, self.tolerance)?;
        Ok(edited)
    }

    /// Split at a periodic-equivalent interior parameter into two open curves.
    ///
    /// The operation cuts the cycle; neither output retains `closed` metadata.
    pub fn split_at(self, parameter: Scalar) -> GeomResult<(BSplineCurve2, BSplineCurve2)> {
        split2(
            self.curve,
            interior_periodic_parameter(parameter, self.domain)?,
        )
    }
}

/// Verified closed-seam extension of a spatial B-spline curve.
///
/// As with [`PeriodicCurve2`], this is an explicit evaluator/editing view over
/// the existing clamped representation, not an inferred periodic knot schema.
#[derive(Debug, Clone, Copy)]
pub struct PeriodicCurve3<'a> {
    curve: &'a BSplineCurve3,
    tolerance: Tolerance,
    domain: (Scalar, Scalar),
    continuity: SeamContinuity,
}

impl<'a> PeriodicCurve3<'a> {
    /// Validate `closed` metadata and at least positional seam continuity.
    pub fn try_new(curve: &'a BSplineCurve3, tolerance: Tolerance) -> GeomResult<Self> {
        let domain = domain(
            curve.degree,
            curve.control_points.len(),
            &curve.knots,
            &curve.multiplicities,
        )?;
        let continuity = curve3_seam_continuity(curve, tolerance)?;
        require_verified_seam(curve.closed, continuity)?;
        Ok(Self {
            curve,
            tolerance,
            domain,
            continuity,
        })
    }

    /// Native active parameter domain retained by this view.
    pub const fn domain(self) -> (Scalar, Scalar) {
        self.domain
    }

    /// Highest endpoint continuity verified during construction.
    pub const fn seam_continuity(self) -> SeamContinuity {
        self.continuity
    }

    /// Wrap a finite parameter into the native active domain.
    pub fn wrap_parameter(self, parameter: Scalar) -> GeomResult<Scalar> {
        wrap(parameter, self.domain)
    }

    /// Evaluate the closed-seam extension at any finite parameter.
    pub fn evaluate(self, parameter: Scalar) -> GeomResult<Point3> {
        Ok(self.jet(parameter)?.point)
    }

    /// Evaluate point, first derivative, and second derivative after wrapping.
    pub fn jet(self, parameter: Scalar) -> GeomResult<CurveJet<Point3, Vec3>> {
        bspline_jet3(self.curve, self.wrap_parameter(parameter)?)
    }

    /// Insert a shape-preserving knot at a periodic-equivalent interior parameter.
    pub fn insert_knot(self, parameter: Scalar) -> GeomResult<BSplineCurve3> {
        let native = interior_periodic_parameter(parameter, self.domain)?;
        let edited = insert_knot3(self.curve, native)?;
        PeriodicCurve3::try_new(&edited, self.tolerance)?;
        Ok(edited)
    }

    /// Split at a periodic-equivalent interior parameter into two open curves.
    pub fn split_at(self, parameter: Scalar) -> GeomResult<(BSplineCurve3, BSplineCurve3)> {
        split3(
            self.curve,
            interior_periodic_parameter(parameter, self.domain)?,
        )
    }
}

fn require_verified_seam(closed: bool, continuity: SeamContinuity) -> GeomResult<()> {
    if !closed || continuity < SeamContinuity::Position {
        return Err(GeomError::InvalidInput(
            "curve is not a verified closed seam".to_owned(),
        ));
    }
    Ok(())
}

fn interior_periodic_parameter(parameter: Scalar, domain: (Scalar, Scalar)) -> GeomResult<Scalar> {
    let native = wrap(parameter, domain)?;
    if native <= domain.0 || native >= domain.1 {
        return Err(GeomError::InvalidInput(
            "periodic edit parameter must not be seam-equivalent".to_owned(),
        ));
    }
    Ok(native)
}

/// Verify native-parameter endpoint continuity of a planar B-spline curve.
pub fn curve2_seam_continuity(
    curve: &BSplineCurve2,
    tolerance: Tolerance,
) -> GeomResult<SeamContinuity> {
    let (lo, hi) = domain(
        curve.degree,
        curve.control_points.len(),
        &curve.knots,
        &curve.multiplicities,
    )?;
    let a = bspline_jet2(curve, lo)?;
    let b = bspline_jet2(curve, hi)?;
    Ok(classify(a, b, tolerance))
}

/// Verify native-parameter endpoint continuity of a spatial B-spline curve.
pub fn curve3_seam_continuity(
    curve: &BSplineCurve3,
    tolerance: Tolerance,
) -> GeomResult<SeamContinuity> {
    let (lo, hi) = domain(
        curve.degree,
        curve.control_points.len(),
        &curve.knots,
        &curve.multiplicities,
    )?;
    let a = bspline_jet3(curve, lo)?;
    let b = bspline_jet3(curve, hi)?;
    Ok(classify(a, b, tolerance))
}

/// Wrap a finite planar-curve parameter into its active domain.
///
/// Wrapping is rejected unless `closed` is set and position continuity is
/// independently verified. An in-domain upper endpoint remains the upper
/// endpoint rather than being remapped to the lower endpoint.
pub fn wrap_curve2_parameter(
    curve: &BSplineCurve2,
    parameter: Scalar,
    tolerance: Tolerance,
) -> GeomResult<Scalar> {
    if !curve.closed || curve2_seam_continuity(curve, tolerance)? < SeamContinuity::Position {
        return Err(GeomError::InvalidInput(
            "curve is not a verified closed seam".to_owned(),
        ));
    }
    wrap(
        parameter,
        domain(
            curve.degree,
            curve.control_points.len(),
            &curve.knots,
            &curve.multiplicities,
        )?,
    )
}

/// Wrap a finite spatial-curve parameter into its active domain.
///
/// The same verified-closed precondition as [`wrap_curve2_parameter`] applies.
pub fn wrap_curve3_parameter(
    curve: &BSplineCurve3,
    parameter: Scalar,
    tolerance: Tolerance,
) -> GeomResult<Scalar> {
    if !curve.closed || curve3_seam_continuity(curve, tolerance)? < SeamContinuity::Position {
        return Err(GeomError::InvalidInput(
            "curve is not a verified closed seam".to_owned(),
        ));
    }
    wrap(
        parameter,
        domain(
            curve.degree,
            curve.control_points.len(),
            &curve.knots,
            &curve.multiplicities,
        )?,
    )
}

fn classify<P, D>(a: CurveJet<P, D>, b: CurveJet<P, D>, tolerance: Tolerance) -> SeamContinuity
where
    P: Copy + core::ops::Sub<P, Output = D>,
    D: Copy + core::ops::Sub<D, Output = D> + Length,
{
    if (a.point - b.point).length() > tolerance.linear() {
        return SeamContinuity::Discontinuous;
    }
    if !close_vector(a.first, b.first, tolerance) {
        return SeamContinuity::Position;
    }
    if !close_vector(a.second, b.second, tolerance) {
        return SeamContinuity::FirstDerivative;
    }
    SeamContinuity::SecondDerivative
}

fn close_vector<D: Copy + core::ops::Sub<D, Output = D> + Length>(
    a: D,
    b: D,
    tolerance: Tolerance,
) -> bool {
    let scale = a.length().max(b.length()).max(1.0);
    (a - b).length() <= tolerance.linear() + tolerance.angular() * scale
}

trait Length {
    fn length(self) -> Scalar;
}
impl Length for axiolid_core::Vec2 {
    fn length(self) -> Scalar {
        self.length()
    }
}
impl Length for axiolid_core::Vec3 {
    fn length(self) -> Scalar {
        self.length()
    }
}

fn domain(
    degree: u16,
    count: usize,
    knots: &[Scalar],
    multiplicities: &[u32],
) -> GeomResult<(Scalar, Scalar)> {
    let spans = active_spans(knots, multiplicities, degree, count)?;
    Ok((spans[0].0, spans[spans.len() - 1].1))
}

fn wrap(parameter: Scalar, (lo, hi): (Scalar, Scalar)) -> GeomResult<Scalar> {
    if !parameter.is_finite() {
        return Err(GeomError::InvalidInput(
            "parameter must be finite".to_owned(),
        ));
    }
    let period = hi - lo;
    if !period.is_finite() || period <= 0.0 {
        return Err(GeomError::InvalidInput(
            "periodic domain must have finite positive length".to_owned(),
        ));
    }
    if parameter >= lo && parameter <= hi {
        return Ok(parameter);
    }
    let offset = parameter - lo;
    if !offset.is_finite() {
        return Err(GeomError::InvalidInput(
            "parameter offset exceeds finite periodic arithmetic".to_owned(),
        ));
    }
    let wrapped = lo + offset.rem_euclid(period);
    if !wrapped.is_finite() || wrapped < lo || wrapped >= hi {
        return Err(GeomError::InvalidInput(
            "parameter could not be wrapped into the periodic domain".to_owned(),
        ));
    }
    Ok(wrapped)
}
