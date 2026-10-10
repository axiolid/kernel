#![forbid(unsafe_code)]

//! Exact two-dimensional section profiles.
//!
//! The crate stores profile intent. Boolean cleanup and triangulation are
//! algorithms in higher tiers so a consumer can use profile data without them.

pub mod center_line;
pub mod contour;
pub mod parameterized;
pub mod validate;

use axiolid_core::Transform2;

pub use center_line::CenterLineProfile;
pub use contour::{Contour, ContourProfile, ProfileSegment};
pub use parameterized::{CircleProfile, EllipseProfile, RectangleProfile, SectionProfile};
pub use validate::ValidateProfile;

/// Format-neutral profile representation.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub enum Profile {
    /// Rectangle or rounded rectangle.
    Rectangle(RectangleProfile),
    /// Circle or annulus.
    Circle(CircleProfile),
    /// Ellipse.
    Ellipse(EllipseProfile),
    /// Structural parameterized section.
    Section(SectionProfile),
    /// Arbitrary exact contour with holes.
    Contour(ContourProfile),
    /// Profile transformed from another profile.
    ///
    /// Placement is recorded as a transform over the exact basis, not by
    /// transforming tessellated points, so a placed circle stays a circle.
    Derived {
        /// Base profile.
        basis: Box<Profile>,
        /// Two-dimensional transform.
        transform: Transform2,
    },
    /// Ordered collection of profiles used as one section.
    Composite(Vec<Profile>),
    /// Open path plus a constant width, centred on the path.
    ///
    /// Distinct from [`Profile::Contour`] because the area is *implied* by an
    /// offset rather than bounded explicitly. A consumer that cannot offset
    /// must refuse this rather than treat the path as a boundary: the path is
    /// open, so reading it as a contour yields a degenerate zero-area region.
    CenterLine(CenterLineProfile),
}
