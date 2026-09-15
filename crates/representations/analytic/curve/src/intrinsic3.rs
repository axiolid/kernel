//! Space curves given by their natural equations: curvature AND torsion as
//! functions of arc length, anchored to a start frame.
//!
//! # Why this is a separate type from `Intrinsic2`
//!
//! In the plane a curve's shape is fixed by curvature alone, and the frame is
//! a single angle: `theta(s) = theta_0 + int k`. Angles commute, so the
//! heading has an elementary closed form and only POSITION needs quadrature.
//!
//! In space the frame is a rotation, and rotations do NOT commute. The
//! Frenet-Serret system
//!
//! ```text
//! T' =        k N
//! N' = -k T       + tau B
//! B' =     -tau N
//! ```
//!
//! is a matrix ODE `R'(s) = R(s) Omega(s)` on SO(3). Its solution is a
//! product integral, not `exp(int Omega)` -- those agree only when the
//! generators at different arc lengths commute, which happens exactly when
//! the ratio `tau/k` is constant (a helix, or a plane curve). So a space
//! curve cannot reuse the 2D trick: the FRAME needs integration too, not
//! just the position.
//!
//! This is why `Intrinsic3` is its own type rather than an `Intrinsic2` with
//! a torsion field bolted on -- the mathematics of recovering a frame is
//! different in kind, not in degree.

use crate::intrinsic::CurvatureLaw;
use axiolid_core::{Frame3, Scalar};

/// A space curve given by its natural equations.
///
/// `curvature` is `k(s) >= 0` by the Frenet convention; a signed law is
/// storable, because dirty imported data stays representable here as
/// everywhere else in this crate, and naming it is a validator's job.
///
/// `torsion` is `tau(s)`, signed: positive is a right-handed screw. Torsion
/// identically zero is a plane curve, and the plane is the one spanned by
/// the start frame's `x` and `y` axes.
///
/// Both laws reuse `CurvatureLaw`, which is the general "scalar function of
/// arc length" in this crate: constant, polynomial, sinusoid, their sum, and
/// piecewise combinations. Nothing about it is curvature-specific, and a
/// second near-identical enum for torsion would be duplication with a
/// different name.
#[derive(Debug, Clone, PartialEq)]
pub struct Intrinsic3 {
    /// Start frame: origin at the curve start, `x` along the start tangent,
    /// `y` along the start normal, `z` along the start binormal.
    pub start: Frame3,
    /// Curvature as a function of arc length from `start`.
    pub curvature: CurvatureLaw,
    /// Torsion as a function of arc length from `start`.
    pub torsion: CurvatureLaw,
    /// Arc length the laws are defined over.
    pub length: Scalar,
}

impl Intrinsic3 {
    /// Anchor a curvature and a torsion law to a start frame over an arc length.
    #[must_use]
    pub const fn new(
        start: Frame3,
        curvature: CurvatureLaw,
        torsion: CurvatureLaw,
        length: Scalar,
    ) -> Self {
        Self {
            start,
            curvature,
            torsion,
            length,
        }
    }

    /// Whether the curve is planar: torsion identically zero.
    ///
    /// A planar space curve is exactly a 2D curve embedded in the start
    /// frame's plane, so this is the predicate that decides whether the
    /// cheaper 2D path applies.
    #[must_use]
    pub fn is_planar(&self) -> bool {
        self.torsion.is_straight()
    }

    /// Whether the curve is a helix: curvature and torsion both constant.
    ///
    /// This is the case where the Frenet generators commute at every arc
    /// length, so the product integral collapses to a single matrix
    /// exponential and the curve has an elementary closed form. Worth naming
    /// because it is both the classical example and the exactly-solvable case.
    #[must_use]
    pub fn is_helical(&self) -> bool {
        self.curvature.is_constant() && self.torsion.is_constant()
    }

    /// Total tangent turning over the curve, in radians, in closed form.
    ///
    /// The integral of `k` -- exact, as in 2D. Note this is NOT enough to
    /// recover the frame in space, only the amount the tangent has swung.
    #[must_use]
    pub fn total_turning(&self) -> Option<Scalar> {
        if !self.length.is_finite() {
            return None;
        }
        crate::Intrinsic2::new(
            axiolid_core::Frame2 {
                origin: axiolid_core::Point2::new(0.0, 0.0),
                x: axiolid_core::Vec2::X,
                y: axiolid_core::Vec2::Y,
            },
            self.curvature.clone(),
            self.length,
        )
        .total_turning()
    }

    /// Total torsion over the curve, in radians, in closed form.
    ///
    /// The integral of `tau`. For a helix this is the angle the binormal has
    /// swung about the axis.
    #[must_use]
    pub fn total_torsion(&self) -> Option<Scalar> {
        if !self.length.is_finite() {
            return None;
        }
        crate::Intrinsic2::new(
            axiolid_core::Frame2 {
                origin: axiolid_core::Point2::new(0.0, 0.0),
                x: axiolid_core::Vec2::X,
                y: axiolid_core::Vec2::Y,
            },
            self.torsion.clone(),
            self.length,
        )
        .total_turning()
    }
}
