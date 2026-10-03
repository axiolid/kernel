//! Solid construction relationships and CSG instructions.

use axiolid_core::{BooleanOperator, Point3, Scalar, Transform3, Vec3};

use crate::{NodeId, StationFrame, StationedSection};

/// Position of one section along a sectioned sweep.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Section {
    /// Profile node.
    pub profile: NodeId,
    /// Local placement of the profile.
    pub placement: Transform3,
}

/// Relationship that constructs a solid from lower-level geometry.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum SolidOperation {
    /// Linear extrusion of a profile.
    Extrusion {
        profile: NodeId,
        direction: Vec3,
        depth: Scalar,
    },
    /// Tapered linear extrusion between two profiles.
    TaperedExtrusion {
        start_profile: NodeId,
        end_profile: NodeId,
        direction: Vec3,
        depth: Scalar,
    },
    /// Revolution of a profile.
    Revolution {
        profile: NodeId,
        axis_origin: Point3,
        axis_direction: Vec3,
        angle: Scalar,
    },
    /// Tapered revolution between two profiles.
    TaperedRevolution {
        start_profile: NodeId,
        end_profile: NodeId,
        axis_origin: Point3,
        axis_direction: Vec3,
        angle: Scalar,
    },
    /// Disk swept along a directrix curve.
    ///
    /// `fillet_radius` rounds the corners where consecutive directrix segments
    /// meet, and is meaningful only on a piecewise-linear directrix: a smooth
    /// curve has no corners to round. `None` means sharp corners, which is
    /// also the correct reading for any directrix that is already smooth.
    ///
    /// It is a property of the SWEEP, not of the disk: the disk stays circular
    /// and it is the swept path whose corners are filleted. A consumer that
    /// cannot round corners must refuse a `Some` rather than drop it, because
    /// silently sharpening a pipe run produces geometry that builds, renders,
    /// and is wrong.
    SweptDisk {
        directrix: NodeId,
        radius: Scalar,
        inner_radius: Option<Scalar>,
        parameter_range: Option<(Scalar, Scalar)>,
        /// Corner rounding radius; `None` means sharp corners.
        fillet_radius: Option<Scalar>,
    },
    /// Profile swept along a directrix using a fixed reference direction.
    FixedReferenceSweep {
        profile: NodeId,
        directrix: NodeId,
        reference_direction: Vec3,
        parameter_range: Option<(Scalar, Scalar)>,
    },
    /// Profile swept along a directrix constrained by a reference surface.
    SurfaceCurveSweep {
        profile: NodeId,
        directrix: NodeId,
        reference_surface: NodeId,
        parameter_range: Option<(Scalar, Scalar)>,
    },
    /// Sections interpolated along a spine.
    SectionedSpine {
        spine: NodeId,
        sections: Vec<Section>,
    },
    /// General CSG binary operation.
    Boolean {
        left: NodeId,
        right: NodeId,
        operator: BooleanOperator,
    },
    /// Unbounded half-space clipped by a finite boundary geometry.
    BoundedHalfSpace {
        half_space: NodeId,
        /// A closed 2D curve in `placement`'s XY plane. A 3D curve is
        /// refused when the graph is built.
        boundary: NodeId,
        /// The boundary's own frame, independent of the clip plane.
        ///
        /// The boundary profile is authored in this frame, so its rotation
        /// orients the profile itself rather than the finished solid. A
        /// source format may place the boundary independently of the base
        /// surface, and only the rotation about the clip normal is
        /// meaningful: the compiler projects this frame's axes into the
        /// plane, so a component along the normal is dropped rather than
        /// tilting the profile out of its own plane.
        placement: Transform3,
    },
    /// Closed profiles standing at stations along a directrix (#241).
    ///
    /// Each profile is placed in the directrix's `frame` at its station
    /// (profile `x` to the lateral axis, `y` to up, normal along the
    /// tangent; see [`crate::station`]). Between two stations the profiles
    /// are matched vertex for vertex, by ring and index, and interpolated
    /// linearly in distance, as are the offsets, so every section must
    /// flatten to the same ring structure. At least two sections, strictly
    /// increasing in distance, are required.
    StationedSpine {
        directrix: NodeId,
        sections: Vec<StationedSection>,
        frame: StationFrame,
    },
}

impl SolidOperation {
    pub(crate) fn references(&self, out: &mut Vec<NodeId>) {
        match self {
            Self::Extrusion { profile, .. } | Self::Revolution { profile, .. } => {
                out.push(*profile)
            }
            Self::TaperedExtrusion {
                start_profile,
                end_profile,
                ..
            }
            | Self::TaperedRevolution {
                start_profile,
                end_profile,
                ..
            } => out.extend([*start_profile, *end_profile]),
            Self::SweptDisk { directrix, .. } => out.push(*directrix),
            Self::FixedReferenceSweep {
                profile, directrix, ..
            } => out.extend([*profile, *directrix]),
            Self::SurfaceCurveSweep {
                profile,
                directrix,
                reference_surface,
                ..
            } => out.extend([*profile, *directrix, *reference_surface]),
            Self::SectionedSpine { spine, sections } => {
                out.push(*spine);
                out.extend(sections.iter().map(|section| section.profile));
            }
            Self::StationedSpine {
                directrix,
                sections,
                ..
            } => {
                out.push(*directrix);
                out.extend(sections.iter().map(|section| section.profile));
            }
            Self::Boolean { left, right, .. } => out.extend([*left, *right]),
            Self::BoundedHalfSpace {
                half_space,
                boundary,
                ..
            } => out.extend([*half_space, *boundary]),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::GeometryGraphBuilder;
    use crate::node::GeometryNode;
    use axiolid_core::Vec3;

    fn directrix(builder: &mut GeometryGraphBuilder) -> NodeId {
        builder.push(GeometryNode::Point3(Vec3::ZERO)).unwrap()
    }

    /// A fillet radius does not add a node reference.
    ///
    /// `references` drives graph traversal and validation, so a scalar that
    /// leaked into it would be read as a `NodeId` and either dangle or alias
    /// an unrelated node. The fillet is geometry data, not a reference.
    #[test]
    fn a_fillet_radius_is_not_a_node_reference() {
        let mut builder = GeometryGraphBuilder::default();
        let curve = directrix(&mut builder);

        let sharp = SolidOperation::SweptDisk {
            directrix: curve,
            radius: 0.05,
            inner_radius: None,
            parameter_range: None,
            fillet_radius: None,
        };
        let rounded = SolidOperation::SweptDisk {
            directrix: curve,
            radius: 0.05,
            inner_radius: None,
            parameter_range: None,
            fillet_radius: Some(0.09),
        };

        let mut sharp_refs = Vec::new();
        sharp.references(&mut sharp_refs);
        let mut rounded_refs = Vec::new();
        rounded.references(&mut rounded_refs);

        assert_eq!(sharp_refs, vec![curve]);
        assert_eq!(
            sharp_refs, rounded_refs,
            "a fillet changes geometry, not the reference graph"
        );
    }

    /// Sharp and rounded sweeps are distinguishable.
    ///
    /// This is the whole reason the field exists. If they compared equal, a
    /// consumer could not tell a filleted pipe run from a mitred one, and
    /// dropping the fillet would be undetectable downstream.
    #[test]
    fn a_fillet_radius_distinguishes_two_otherwise_identical_sweeps() {
        let mut builder = GeometryGraphBuilder::default();
        let curve = directrix(&mut builder);

        let common = |fillet| SolidOperation::SweptDisk {
            directrix: curve,
            radius: 0.05,
            inner_radius: Some(0.04),
            parameter_range: Some((0.0, 2.0)),
            fillet_radius: fillet,
        };
        assert_ne!(common(None), common(Some(0.09)));
        assert_eq!(common(Some(0.09)), common(Some(0.09)));
    }
}
