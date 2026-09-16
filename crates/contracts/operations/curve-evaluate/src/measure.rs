//! How a caller's number locates a place on a curve.

use axiolid_core::Scalar;

/// Where along a curve to evaluate, and by WHICH method of measurement.
///
/// Two different quantities are routinely written into exchange files for
/// the same purpose, and they are not interchangeable:
///
/// - a LENGTH from the start of the curve, and
/// - a dimensionless native PARAMETER, which is not a distance at all.
///
/// IFC4x3 makes the choice explicit in the file
/// (`IfcCurveMeasureSelect`, either `IfcNonNegativeLengthMeasure` or
/// `IfcParameterValue`), and STEP carries the same distinction. A
/// consumer maps the authored value to this enum ONCE and passes it
/// through; it never has to decide which evaluator to call.
///
/// Why this is a value and not two method names: with separate
/// `*_at_parameter` methods the caller still branches, and a caller that
/// branches can branch wrongly. Passing `1.5` as a distance to a circle
/// authored in radians places an object 1.5 m along a curve where the
/// author meant about 86 degrees around it -- wrong, finite, and
/// plausible, which is the worst combination. Carrying the method of
/// measurement in the value makes that mistake unrepresentable.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CurveMeasure {
    /// A length from the start of the curve.
    ///
    /// Interpreted in the provider's
    /// [`DistanceConvention`](crate::DistanceConvention) for that curve,
    /// which is why the convention must be consulted before trusting it.
    Distance(Scalar),
    /// The curve's own parameter: an angle for a circle, a knot value for
    /// a spline, arc length for an intrinsic curve.
    ///
    /// Meaningful for every curve family, including those whose arc length
    /// has no closed form, so this route stays available where
    /// [`Distance`](Self::Distance) is refused.
    Parameter(Scalar),
}

impl CurveMeasure {
    /// The carried number, whichever method of measurement it uses.
    ///
    /// For validity checks such as finiteness that apply to both.
    #[must_use]
    pub fn value(self) -> Scalar {
        match self {
            Self::Distance(v) | Self::Parameter(v) => v,
        }
    }

    /// Whether this is a length rather than a native parameter.
    #[must_use]
    pub fn is_distance(self) -> bool {
        matches!(self, Self::Distance(_))
    }
}
