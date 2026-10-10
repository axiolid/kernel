//! Compact parameterized profile families.

use axiolid_core::Scalar;

/// Rectangle, optionally rounded at the corners.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct RectangleProfile {
    /// Extent along local x.
    pub x: Scalar,
    /// Extent along local y.
    pub y: Scalar,
    /// Optional wall thickness. `None` denotes a filled section.
    pub thickness: Option<Scalar>,
    /// Optional outer corner radius.
    pub outer_radius: Option<Scalar>,
    /// Optional inner corner radius for hollow sections.
    pub inner_radius: Option<Scalar>,
}

/// Circle or annulus.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct CircleProfile {
    /// Outer radius.
    pub radius: Scalar,
    /// Optional wall thickness. `None` denotes a filled disk.
    pub thickness: Option<Scalar>,
}

/// Ellipse.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct EllipseProfile {
    /// Semi-axis along local x.
    pub semi_axis_x: Scalar,
    /// Semi-axis along local y.
    pub semi_axis_y: Scalar,
}

/// Generic structural section dimensions.
///
/// Each variant carries the dimensions the corresponding source entity
/// declares, including the optional fillet radii, edge radii and flange
/// slopes. Those optionals are not decoration: for a rolled steel section the
/// root fillet is a large share of the cross-sectional area, so discarding it
/// yields a section whose area, second moment and mass are all wrong while
/// still looking like the right shape.
///
/// Angles are radians, lengths are model units. `None` means the source did
/// not state the value, which is distinct from stating zero.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub enum SectionProfile {
    /// Symmetric I section: equal top and bottom flanges.
    I {
        depth: Scalar,
        width: Scalar,
        web_thickness: Scalar,
        flange_thickness: Scalar,
        /// Web-to-flange root fillet.
        fillet_radius: Option<Scalar>,
        /// Radius at the flange toe.
        flange_edge_radius: Option<Scalar>,
        /// Flange taper, radians from the horizontal.
        flange_slope: Option<Scalar>,
    },
    /// I section with independent top and bottom flanges.
    ///
    /// Separate from `I` because collapsing the two flange widths into one
    /// silently symmetrises the section: the shape stays plausible while the
    /// neutral axis, and therefore every bending result, moves.
    AsymmetricI {
        depth: Scalar,
        web_thickness: Scalar,
        bottom_flange_width: Scalar,
        bottom_flange_thickness: Scalar,
        bottom_fillet_radius: Option<Scalar>,
        bottom_flange_edge_radius: Option<Scalar>,
        bottom_flange_slope: Option<Scalar>,
        top_flange_width: Scalar,
        /// Defaults to the bottom thickness when absent.
        top_flange_thickness: Option<Scalar>,
        top_fillet_radius: Option<Scalar>,
        top_flange_edge_radius: Option<Scalar>,
        top_flange_slope: Option<Scalar>,
    },
    /// L or angle section.
    L {
        depth: Scalar,
        /// Defaults to the depth when absent, giving an equal angle.
        width: Option<Scalar>,
        thickness: Scalar,
        fillet_radius: Option<Scalar>,
        edge_radius: Option<Scalar>,
        /// Leg taper, radians.
        leg_slope: Option<Scalar>,
    },
    /// T section.
    T {
        depth: Scalar,
        flange_width: Scalar,
        web_thickness: Scalar,
        flange_thickness: Scalar,
        fillet_radius: Option<Scalar>,
        flange_edge_radius: Option<Scalar>,
        web_edge_radius: Option<Scalar>,
        web_slope: Option<Scalar>,
        flange_slope: Option<Scalar>,
    },
    /// U or channel section.
    U {
        depth: Scalar,
        flange_width: Scalar,
        web_thickness: Scalar,
        flange_thickness: Scalar,
        fillet_radius: Option<Scalar>,
        edge_radius: Option<Scalar>,
        flange_slope: Option<Scalar>,
    },
    /// C or lipped-channel section.
    C {
        depth: Scalar,
        width: Scalar,
        wall_thickness: Scalar,
        /// Length of the returned lip.
        girth: Scalar,
        internal_fillet_radius: Option<Scalar>,
    },
    /// Z section.
    Z {
        depth: Scalar,
        flange_width: Scalar,
        web_thickness: Scalar,
        flange_thickness: Scalar,
        fillet_radius: Option<Scalar>,
        edge_radius: Option<Scalar>,
    },
    /// Trapezium.
    Trapezium {
        bottom_x: Scalar,
        top_x: Scalar,
        y: Scalar,
        /// Offset of the top edge relative to the bottom; may be negative.
        top_offset: Scalar,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An asymmetric I is not representable as a symmetric one.
    ///
    /// This is the whole reason `AsymmetricI` is a separate variant. Before it
    /// existed an adapter had to pick one of the two flange widths, producing
    /// a section that looks right and whose neutral axis is in the wrong
    /// place. Distinct variants make that substitution impossible to express.
    #[test]
    fn an_asymmetric_i_is_a_distinct_variant_from_a_symmetric_one() {
        let symmetric = SectionProfile::I {
            depth: 0.4,
            width: 0.3,
            web_thickness: 0.011,
            flange_thickness: 0.019,
            fillet_radius: Some(0.021),
            flange_edge_radius: None,
            flange_slope: None,
        };
        let asymmetric = SectionProfile::AsymmetricI {
            depth: 0.4,
            web_thickness: 0.011,
            bottom_flange_width: 0.3,
            bottom_flange_thickness: 0.019,
            bottom_fillet_radius: Some(0.021),
            bottom_flange_edge_radius: None,
            bottom_flange_slope: None,
            // The distinguishing dimension: a wider top flange.
            top_flange_width: 0.2,
            top_flange_thickness: Some(0.019),
            top_fillet_radius: Some(0.021),
            top_flange_edge_radius: None,
            top_flange_slope: None,
        };
        assert_ne!(symmetric, asymmetric);
    }

    /// A stated fillet radius is not the same as an absent one.
    ///
    /// `None` means the source did not declare the value; `Some(0.0)` means it
    /// declared a sharp corner. For a rolled section the root fillet carries
    /// real area, so the two must not compare equal.
    #[test]
    fn an_absent_fillet_differs_from_a_declared_zero() {
        let make = |fillet| SectionProfile::U {
            depth: 0.3,
            flange_width: 0.1,
            web_thickness: 0.0075,
            flange_thickness: 0.0125,
            fillet_radius: fillet,
            edge_radius: None,
            flange_slope: None,
        };
        assert_ne!(make(None), make(Some(0.0)));
        assert_eq!(make(Some(0.015)), make(Some(0.015)));
    }

    /// Every section variant survives a round trip through `Profile`.
    #[test]
    fn a_section_wraps_into_a_profile_without_loss() {
        let section = SectionProfile::L {
            depth: 0.1,
            width: None,
            thickness: 0.008,
            fillet_radius: Some(0.012),
            edge_radius: Some(0.006),
            leg_slope: None,
        };
        let profile = crate::Profile::Section(section.clone());
        match profile {
            crate::Profile::Section(round_tripped) => assert_eq!(section, round_tripped),
            other => panic!("expected a section, got {other:?}"),
        }
    }
}
