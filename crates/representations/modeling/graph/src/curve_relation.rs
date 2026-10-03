//! Curve relationships that require graph references.

use axiolid_core::{Point2, Point3, Scalar, Vec3};

use crate::NodeId;

/// One trim selector preserved from a source representation.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TrimSelector {
    /// Curve parameter.
    Parameter(Scalar),
    /// Two-dimensional point.
    Point2(Point2),
    /// Three-dimensional point.
    Point3(Point3),
    /// Arc length along the basis curve.
    ///
    /// Selects the basis parameter at which the arc length along the basis,
    /// measured from the basis curve's parameter `0` in the sense of the
    /// trim, equals the value: in increasing parameter when the trim agrees
    /// with the basis, in decreasing parameter when it does not. A negative
    /// value measures the other way.
    ///
    /// The value is exact stored data, and it is a parameter-kind selector:
    /// it names a position by a measure along the curve, not by a point, so
    /// it satisfies [`TrimmingPreference::Parameter`]. Most curves have no
    /// closed-form arc length (a cubic parabola's is elliptic), so an
    /// evaluator resolves the parameter by quadrature and a root find to a
    /// tolerance it states, and a path that cannot refuses by name.
    ArcLength(Scalar),
}

/// Preference when both parameter and Cartesian trim selectors exist.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TrimmingPreference {
    /// Prefer parameter values.
    Parameter,
    /// Prefer Cartesian points.
    Cartesian,
    /// Use source order when no preference was stated.
    Unspecified,
}

/// Continuity declared between consecutive composite segments.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Transition {
    /// Discontinuous.
    Discontinuous,
    /// Position continuous.
    Continuous,
    /// Position and tangent continuous.
    ContinuousSameGradient,
    /// Position, tangent, and curvature continuous.
    ContinuousSameGradientSameCurvature,
}

/// One oriented curve in a composite.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CurveSegment {
    /// Child curve.
    pub curve: NodeId,
    /// Whether child parameterization agrees with composite orientation.
    pub same_sense: bool,
    /// Transition from the preceding segment.
    pub transition: Transition,
}

/// Relationship between curve nodes.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum CurveRelation {
    /// Ordered composite curve.
    Composite { segments: Vec<CurveSegment> },
    /// Trimmed view of a basis curve.
    Trimmed {
        basis: NodeId,
        start: Vec<TrimSelector>,
        end: Vec<TrimSelector>,
        sense_agreement: bool,
        preference: TrimmingPreference,
    },
    /// Constant-distance offset.
    Offset {
        basis: NodeId,
        distance: Scalar,
        reference_direction: Option<Vec3>,
    },
    /// Three-dimensional curve associated with one or more surfaces/pcurves.
    SurfaceCurve {
        curve_3d: NodeId,
        /// The parametric sides, each pairing a surface with its own p-curve.
        sides: SurfaceSides,
        master: MasterRepresentation,
    },
    /// Two-dimensional parameter curve on a surface.
    ParameterCurve {
        basis_surface: NodeId,
        reference_curve: NodeId,
    },
}

/// Which representation governs a redundant surface-curve definition.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MasterRepresentation {
    /// Three-dimensional curve.
    Curve3d,
    /// The p-curve on the FIRST parametric side governs.
    ParameterCurveS1,
    /// The p-curve on the SECOND parametric side governs.
    ///
    /// Only meaningful on a two-sided curve; naming it on a single-sided one
    /// is contradictory rather than merely unusual, and is refused.
    ParameterCurveS2,
    /// Both are authoritative and must agree.
    Both,
    /// Unspecified.
    Unspecified,
}

/// Which parametric side of a surface curve a p-curve belongs to.
///
/// A surface curve is the intersection of two surfaces, so each side owns one
/// surface and the p-curve that is that curve's image in the surface's own
/// parameter domain. Pairing them here is what stops a consumer from having to
/// guess which p-curve to trim with -- a guess that would otherwise require
/// re-inverting the surface, the exact operation the p-curve exists to avoid.
///
/// A single side is legitimate: not every edge has two parametric images.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SurfaceSides {
    first: (NodeId, NodeId),
    second: Option<(NodeId, NodeId)>,
}

impl SurfaceSides {
    /// One parametric side: a surface and this curve's image in it.
    #[must_use]
    pub const fn one(surface: NodeId, pcurve: NodeId) -> Self {
        Self {
            first: (surface, pcurve),
            second: None,
        }
    }

    /// Both parametric sides, in the order the authoring format states them.
    ///
    /// The order is load-bearing: it is what `ParameterCurveS1` and
    /// `ParameterCurveS2` name.
    #[must_use]
    pub const fn two(
        first_surface: NodeId,
        first_pcurve: NodeId,
        second_surface: NodeId,
        second_pcurve: NodeId,
    ) -> Self {
        Self {
            first: (first_surface, first_pcurve),
            second: Some((second_surface, second_pcurve)),
        }
    }

    /// The first side as `(surface, pcurve)`.
    #[must_use]
    pub const fn first(&self) -> (NodeId, NodeId) {
        self.first
    }

    /// The second side as `(surface, pcurve)`, if this curve has one.
    #[must_use]
    pub const fn second(&self) -> Option<(NodeId, NodeId)> {
        self.second
    }

    /// Whether both parametric sides are present.
    #[must_use]
    pub const fn is_two_sided(&self) -> bool {
        self.second.is_some()
    }

    /// Every node this pairing references, surfaces and p-curves alike.
    #[must_use]
    pub fn references(&self) -> Vec<NodeId> {
        let mut out = vec![self.first.0, self.first.1];
        if let Some((surface, pcurve)) = self.second {
            out.extend([surface, pcurve]);
        }
        out
    }
}

impl CurveRelation {
    pub(crate) fn references(&self, out: &mut Vec<NodeId>) {
        match self {
            Self::Composite { segments } => out.extend(segments.iter().map(|item| item.curve)),
            Self::Trimmed { basis, .. } | Self::Offset { basis, .. } => out.push(*basis),
            Self::SurfaceCurve {
                curve_3d, sides, ..
            } => {
                out.push(*curve_3d);
                out.extend(sides.references());
            }
            Self::ParameterCurve {
                basis_surface,
                reference_curve,
            } => out.extend([*basis_surface, *reference_curve]),
        }
    }
}
