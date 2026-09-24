//! Semantic validation for graph references.

use std::collections::HashSet;

use axiolid_curve::{BSplineCurve2, Curve2};

use crate::{
    CurveRelation, GeometryNode, GraphError, MasterRepresentation, NodeId, SolidOperation,
    SurfaceRelation, TrimSelector, TrimmingPreference,
};

#[derive(Debug, Clone, Copy)]
enum ExpectedReference {
    Curve,
    Curve2,
    BoundedOpenCurve2,
    Curve3,
    Surface,
    Profile,
    Solid,
    HalfSpace,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CurveDimension {
    Two,
    Three,
}

fn curve_has_dimension(root: NodeId, nodes: &[GeometryNode], dimension: CurveDimension) -> bool {
    let mut pending = vec![root];
    let mut visited = HashSet::new();

    while let Some(node_id) = pending.pop() {
        if !visited.insert(node_id.index()) {
            continue;
        }
        match &nodes[node_id.index()] {
            GeometryNode::Instance(instance) => pending.push(instance.source),
            GeometryNode::Curve2(_) => {
                if dimension != CurveDimension::Two {
                    return false;
                }
            }
            GeometryNode::Curve3(_) => {
                if dimension != CurveDimension::Three {
                    return false;
                }
            }
            GeometryNode::CurveRelation(CurveRelation::Trimmed { basis, .. })
            | GeometryNode::CurveRelation(CurveRelation::Offset { basis, .. }) => {
                pending.push(*basis);
            }
            GeometryNode::CurveRelation(CurveRelation::Composite { segments }) => {
                pending.extend(segments.iter().map(|segment| segment.curve));
            }
            GeometryNode::CurveRelation(CurveRelation::SurfaceCurve { .. }) => {
                if dimension != CurveDimension::Three {
                    return false;
                }
            }
            GeometryNode::CurveRelation(CurveRelation::ParameterCurve { .. }) => {
                if dimension != CurveDimension::Two {
                    return false;
                }
            }
            _ => return false,
        }
    }
    true
}

fn bspline_is_structurally_valid_2d(curve: &BSplineCurve2) -> bool {
    let degree = usize::from(curve.degree);
    let expected_sum = curve
        .control_points
        .len()
        .checked_add(degree)
        .and_then(|value| value.checked_add(1));
    let actual_sum = curve
        .multiplicities
        .iter()
        .try_fold(0usize, |sum, value| sum.checked_add(*value as usize));
    let weights_are_valid = curve.weights.as_ref().is_none_or(|weights| {
        weights.len() == curve.control_points.len()
            && weights
                .iter()
                .all(|weight| weight.is_finite() && *weight > 0.0)
    });

    degree > 0
        && curve.control_points.len() > degree
        && curve.control_points.iter().all(|point| point.is_finite())
        && !curve.knots.is_empty()
        && curve.knots.iter().all(|knot| knot.is_finite())
        && curve.knots.windows(2).all(|pair| pair[0] < pair[1])
        && curve.multiplicities.len() == curve.knots.len()
        && curve
            .multiplicities
            .iter()
            .all(|value| *value > 0 && *value <= u32::from(curve.degree) + 1)
        && actual_sum == expected_sum
        && weights_are_valid
}

fn curve2_is_structurally_valid_trim_basis(curve: &Curve2) -> bool {
    match curve {
        Curve2::Line(line) => {
            line.origin.is_finite()
                && line.direction.is_finite()
                && line.direction.length_squared() > 0.0
        }
        Curve2::Circle(circle) => {
            circle.frame.origin.is_finite()
                && circle.frame.x.is_finite()
                && circle.frame.y.is_finite()
                && circle.frame.x.perp_dot(circle.frame.y) != 0.0
                && circle.radius.is_finite()
                && circle.radius > 0.0
        }
        Curve2::Ellipse(ellipse) => {
            ellipse.frame.origin.is_finite()
                && ellipse.frame.x.is_finite()
                && ellipse.frame.y.is_finite()
                && ellipse.frame.x.perp_dot(ellipse.frame.y) != 0.0
                && ellipse.semi_axis_x.is_finite()
                && ellipse.semi_axis_x > 0.0
                && ellipse.semi_axis_y.is_finite()
                && ellipse.semi_axis_y > 0.0
        }
        Curve2::Polyline(polyline) => {
            polyline.points.len() >= 2 && polyline.points.iter().all(|point| point.is_finite())
        }
        Curve2::BSpline(spline) => bspline_is_structurally_valid_2d(spline),
        Curve2::Sinusoid(wave) => wave.is_finite(),
        _ => false,
    }
}

fn trim_selector_is_finite_2d(selector: &TrimSelector) -> bool {
    match selector {
        TrimSelector::Parameter(value) => value.is_finite(),
        TrimSelector::Point2(point) => point.is_finite(),
        TrimSelector::Point3(_) => false,
    }
}

fn trim_end_supports_preference(
    selectors: &[TrimSelector],
    preference: TrimmingPreference,
) -> bool {
    match preference {
        TrimmingPreference::Parameter => selectors
            .iter()
            .any(|selector| matches!(selector, TrimSelector::Parameter(_))),
        TrimmingPreference::Cartesian => selectors
            .iter()
            .any(|selector| matches!(selector, TrimSelector::Point2(_))),
        TrimmingPreference::Unspecified => true,
    }
}

fn trim_selectors_definitely_equal(start: &[TrimSelector], end: &[TrimSelector]) -> bool {
    start.iter().any(|left| {
        end.iter().any(|right| match (left, right) {
            (TrimSelector::Parameter(a), TrimSelector::Parameter(b)) => a == b,
            (TrimSelector::Point2(a), TrimSelector::Point2(b)) => a == b,
            _ => false,
        })
    })
}

fn trim_declaration_is_structurally_open_2d(
    start: &[TrimSelector],
    end: &[TrimSelector],
    preference: TrimmingPreference,
) -> bool {
    !start.is_empty()
        && !end.is_empty()
        && start.iter().all(trim_selector_is_finite_2d)
        && end.iter().all(trim_selector_is_finite_2d)
        && trim_end_supports_preference(start, preference)
        && trim_end_supports_preference(end, preference)
        && !trim_selectors_definitely_equal(start, end)
}

fn curve_is_valid_2d_trim_basis(root: NodeId, nodes: &[GeometryNode]) -> bool {
    let mut pending = vec![root];
    let mut visited = HashSet::new();

    'pending: while let Some(node_id) = pending.pop() {
        if !visited.insert(node_id.index()) {
            continue;
        }
        let mut node = &nodes[node_id.index()];
        while let GeometryNode::Instance(instance) = node {
            if !instance.transform.is_finite() {
                return false;
            }
            if !visited.insert(instance.source.index()) {
                continue 'pending;
            }
            node = &nodes[instance.source.index()];
        }
        match node {
            GeometryNode::Curve2(curve) => {
                if !curve2_is_structurally_valid_trim_basis(curve) {
                    return false;
                }
            }
            GeometryNode::CurveRelation(CurveRelation::Trimmed {
                basis,
                start,
                end,
                preference,
                ..
            }) => {
                if !trim_declaration_is_structurally_open_2d(start, end, *preference) {
                    return false;
                }
                pending.push(*basis);
            }
            GeometryNode::CurveRelation(CurveRelation::Composite { segments }) => {
                if segments.is_empty() {
                    return false;
                }
                pending.extend(segments.iter().map(|segment| segment.curve));
            }
            GeometryNode::CurveRelation(CurveRelation::Offset {
                basis,
                distance,
                reference_direction,
            }) => {
                if !distance.is_finite() || reference_direction.is_some() {
                    return false;
                }
                pending.push(*basis);
            }
            GeometryNode::CurveRelation(CurveRelation::ParameterCurve {
                reference_curve, ..
            }) => pending.push(*reference_curve),
            GeometryNode::CurveRelation(CurveRelation::SurfaceCurve { .. }) | _ => return false,
        }
    }
    true
}

fn trimmed_curve_is_structurally_open_2d(
    basis: NodeId,
    start: &[TrimSelector],
    end: &[TrimSelector],
    preference: TrimmingPreference,
    nodes: &[GeometryNode],
) -> bool {
    trim_declaration_is_structurally_open_2d(start, end, preference)
        && curve_is_valid_2d_trim_basis(basis, nodes)
}

fn curve_is_bounded_open_2d(root: NodeId, nodes: &[GeometryNode]) -> bool {
    let mut pending = vec![root];
    let mut visited = HashSet::new();

    'pending: while let Some(node_id) = pending.pop() {
        if !visited.insert(node_id.index()) {
            continue;
        }
        let mut node = &nodes[node_id.index()];
        while let GeometryNode::Instance(instance) = node {
            if !instance.transform.is_finite() {
                return false;
            }
            if !visited.insert(instance.source.index()) {
                continue 'pending;
            }
            node = &nodes[instance.source.index()];
        }
        match node {
            GeometryNode::Curve2(Curve2::Polyline(curve)) => {
                if curve.closed
                    || curve.points.len() < 2
                    || !curve.points.iter().all(|point| point.is_finite())
                    || curve.points.first() == curve.points.last()
                {
                    return false;
                }
            }
            GeometryNode::Curve2(Curve2::BSpline(curve)) => {
                if curve.closed || !bspline_is_structurally_valid_2d(curve) {
                    return false;
                }
            }
            GeometryNode::Curve2(_) => return false,
            GeometryNode::CurveRelation(CurveRelation::Trimmed {
                basis,
                start,
                end,
                preference,
                ..
            }) => {
                if !trimmed_curve_is_structurally_open_2d(*basis, start, end, *preference, nodes) {
                    return false;
                }
            }
            GeometryNode::CurveRelation(CurveRelation::Composite { segments }) => {
                if segments.is_empty() {
                    return false;
                }
                pending.extend(segments.iter().map(|segment| segment.curve));
            }
            GeometryNode::CurveRelation(CurveRelation::Offset {
                basis,
                distance,
                reference_direction,
            }) => {
                if !distance.is_finite() || reference_direction.is_some() {
                    return false;
                }
                pending.push(*basis);
            }
            GeometryNode::CurveRelation(CurveRelation::ParameterCurve {
                reference_curve, ..
            }) => pending.push(*reference_curve),
            GeometryNode::CurveRelation(CurveRelation::SurfaceCurve { .. }) | _ => return false,
        }
    }
    true
}

impl ExpectedReference {
    const fn description(self) -> &'static str {
        match self {
            Self::Curve => "curve",
            Self::Curve2 => "curve2",
            Self::BoundedOpenCurve2 => "bounded open curve2",
            Self::Curve3 => "curve3",
            Self::Surface => "surface",
            Self::Profile => "profile",
            Self::Solid => "solid",
            Self::HalfSpace => "half-space",
        }
    }

    fn accepts<'a>(
        self,
        reference: NodeId,
        mut node: &'a GeometryNode,
        nodes: &'a [GeometryNode],
    ) -> bool {
        match self {
            Self::Curve => {
                return curve_has_dimension(reference, nodes, CurveDimension::Two)
                    || curve_has_dimension(reference, nodes, CurveDimension::Three);
            }
            Self::Curve2 => {
                return curve_has_dimension(reference, nodes, CurveDimension::Two);
            }
            Self::BoundedOpenCurve2 => return curve_is_bounded_open_2d(reference, nodes),
            Self::Curve3 => {
                return curve_has_dimension(reference, nodes, CurveDimension::Three);
            }
            _ => {}
        }

        // An instance preserves the dimensional/reference family of its source.
        while let GeometryNode::Instance(instance) = node {
            node = &nodes[instance.source.index()];
        }

        let surface = matches!(
            node,
            GeometryNode::Surface(_) | GeometryNode::SurfaceRelation(_)
        );
        match self {
            Self::Curve | Self::Curve2 | Self::BoundedOpenCurve2 | Self::Curve3 => false,
            Self::Surface => surface,
            Self::Profile => matches!(node, GeometryNode::Profile(_)),
            Self::Solid => matches!(
                node,
                GeometryNode::Primitive(_)
                    | GeometryNode::HalfSpace(_)
                    | GeometryNode::SolidOperation(_)
                    | GeometryNode::BRep(_)
                    | GeometryNode::PolygonMesh(_)
                    | GeometryNode::TriMesh(_)
            ),
            Self::HalfSpace => matches!(node, GeometryNode::HalfSpace(_)),
        }
    }
}

pub(crate) fn validate_reference_types(
    node: &GeometryNode,
    nodes: &[GeometryNode],
) -> Result<(), GraphError> {
    // Keep this match exhaustive so adding a node variant cannot silently bypass
    // semantic reference validation.
    match node {
        GeometryNode::CurveRelation(value) => validate_curve_relation(value, nodes),
        GeometryNode::PointOnCurve(value) => {
            expect_reference(nodes, value.curve, ExpectedReference::Curve)
        }
        GeometryNode::SurfaceRelation(value) => validate_surface_relation(value, nodes),
        GeometryNode::PointOnSurface(value) => {
            expect_reference(nodes, value.surface, ExpectedReference::Surface)
        }
        GeometryNode::OpenProfile(value) => {
            expect_reference(nodes, value.path, ExpectedReference::BoundedOpenCurve2)
        }
        GeometryNode::SolidOperation(value) => validate_solid_operation(value, nodes),
        GeometryNode::BRep(value) => {
            for edge in value.edges() {
                if let Some(curve) = edge.curve {
                    expect_reference(nodes, curve, ExpectedReference::Curve)?;
                }
            }
            for face in value.faces() {
                if let Some(surface) = face.surface {
                    expect_reference(nodes, surface, ExpectedReference::Surface)?;
                }
            }
            // A pcurve is a 2D curve in a surface's parameter domain, so it
            // must resolve to a Curve2 node. Accepting a Curve3 here would
            // let a model claim a trim it cannot supply.
            for wire in value.loops() {
                for use_ in &wire.edges {
                    if let Some(pcurve) = use_.pcurve {
                        expect_reference(nodes, pcurve, ExpectedReference::Curve)?;
                    }
                }
            }
            Ok(())
        }
        GeometryNode::Point2(_)
        | GeometryNode::Point3(_)
        | GeometryNode::Vector2(_)
        | GeometryNode::Vector3(_)
        | GeometryNode::Frame2(_)
        | GeometryNode::Frame3(_)
        | GeometryNode::Transform(_)
        | GeometryNode::PointList2(_)
        | GeometryNode::PointList3(_)
        | GeometryNode::Curve2(_)
        | GeometryNode::Curve3(_)
        | GeometryNode::Surface(_)
        | GeometryNode::Profile(_) => Ok(()),
        GeometryNode::Primitive(_) | GeometryNode::HalfSpace(_) | GeometryNode::PolygonMesh(_) => {
            Ok(())
        }
        GeometryNode::TriMesh(_)
        | GeometryNode::BoundingBox(_)
        | GeometryNode::Instance(_)
        | GeometryNode::Collection(_) => Ok(()),
    }
}

fn expect_reference(
    nodes: &[GeometryNode],
    reference: NodeId,
    expected: ExpectedReference,
) -> Result<(), GraphError> {
    let actual = &nodes[reference.index()];
    if expected.accepts(reference, actual, nodes) {
        return Ok(());
    }
    Err(GraphError::InvalidReferenceType {
        reference,
        expected: expected.description(),
        actual: node_kind(actual),
    })
}

fn node_kind(node: &GeometryNode) -> &'static str {
    match node {
        GeometryNode::Point2(_) => "point2",
        GeometryNode::Point3(_) => "point3",
        GeometryNode::Vector2(_) => "vector2",
        GeometryNode::Vector3(_) => "vector3",
        GeometryNode::Frame2(_) => "frame2",
        GeometryNode::Frame3(_) => "frame3",
        GeometryNode::Transform(_) => "transform",
        GeometryNode::PointList2(_) => "point-list2",
        GeometryNode::PointList3(_) => "point-list3",
        GeometryNode::Curve2(_) => "curve2",
        GeometryNode::Curve3(_) => "curve3",
        GeometryNode::CurveRelation(_) => "curve-relation",
        GeometryNode::PointOnCurve(_) => "point-on-curve",
        GeometryNode::Surface(_) => "surface",
        GeometryNode::SurfaceRelation(_) => "surface-relation",
        GeometryNode::PointOnSurface(_) => "point-on-surface",
        GeometryNode::Profile(_) => "profile",
        GeometryNode::OpenProfile(_) => "open-profile",
        GeometryNode::Primitive(_) => "primitive",
        GeometryNode::HalfSpace(_) => "half-space",
        GeometryNode::SolidOperation(_) => "solid-operation",
        GeometryNode::BRep(_) => "brep",
        GeometryNode::PolygonMesh(_) => "polygon-mesh",
        GeometryNode::TriMesh(_) => "triangle-mesh",
        GeometryNode::BoundingBox(_) => "bounding-box",
        GeometryNode::Instance(_) => "instance",
        GeometryNode::Collection(_) => "collection",
    }
}

fn validate_curve_relation(
    relation: &CurveRelation,
    nodes: &[GeometryNode],
) -> Result<(), GraphError> {
    match relation {
        CurveRelation::Composite { segments } => {
            for segment in segments {
                expect_reference(nodes, segment.curve, ExpectedReference::Curve)?;
            }
            Ok(())
        }
        CurveRelation::Trimmed { basis, .. } | CurveRelation::Offset { basis, .. } => {
            expect_reference(nodes, *basis, ExpectedReference::Curve)
        }
        CurveRelation::SurfaceCurve {
            curve_3d,
            sides,
            master,
        } => {
            expect_reference(nodes, *curve_3d, ExpectedReference::Curve3)?;
            // Each side is a (surface, pcurve) pair, so the two references are
            // checked against their own kinds rather than a permissive
            // "curve or surface" that accepted any mix.
            let (first_surface, first_pcurve) = sides.first();
            expect_reference(nodes, first_surface, ExpectedReference::Surface)?;
            expect_reference(nodes, first_pcurve, ExpectedReference::Curve2)?;
            if let Some((second_surface, second_pcurve)) = sides.second() {
                expect_reference(nodes, second_surface, ExpectedReference::Surface)?;
                expect_reference(nodes, second_pcurve, ExpectedReference::Curve2)?;
            }
            // A master naming a side the curve does not have is contradictory:
            // it says "believe the second image" when there is no second image.
            if *master == MasterRepresentation::ParameterCurveS2 && !sides.is_two_sided() {
                return Err(GraphError::ContradictoryMaster {
                    detail: "master names the second parametric side, but the \
                             surface curve has only one",
                });
            }
            Ok(())
        }
        CurveRelation::ParameterCurve {
            basis_surface,
            reference_curve,
        } => {
            expect_reference(nodes, *basis_surface, ExpectedReference::Surface)?;
            expect_reference(nodes, *reference_curve, ExpectedReference::Curve2)
        }
    }
}

fn validate_surface_relation(
    relation: &SurfaceRelation,
    nodes: &[GeometryNode],
) -> Result<(), GraphError> {
    match relation {
        SurfaceRelation::CurveBounded {
            basis, boundaries, ..
        } => {
            expect_reference(nodes, *basis, ExpectedReference::Surface)?;
            for boundary in boundaries {
                expect_reference(nodes, *boundary, ExpectedReference::Curve)?;
            }
            Ok(())
        }
        SurfaceRelation::RectangularTrimmed { basis, .. }
        | SurfaceRelation::Offset { basis, .. } => {
            expect_reference(nodes, *basis, ExpectedReference::Surface)
        }
        SurfaceRelation::LinearExtrusion { swept_curve, .. }
        | SurfaceRelation::Revolution { swept_curve, .. } => {
            expect_reference(nodes, *swept_curve, ExpectedReference::Curve)
        }
    }
}

fn validate_solid_operation(
    operation: &SolidOperation,
    nodes: &[GeometryNode],
) -> Result<(), GraphError> {
    match operation {
        SolidOperation::Extrusion { profile, .. } | SolidOperation::Revolution { profile, .. } => {
            expect_reference(nodes, *profile, ExpectedReference::Profile)
        }
        SolidOperation::TaperedExtrusion {
            start_profile,
            end_profile,
            ..
        }
        | SolidOperation::TaperedRevolution {
            start_profile,
            end_profile,
            ..
        } => {
            expect_reference(nodes, *start_profile, ExpectedReference::Profile)?;
            expect_reference(nodes, *end_profile, ExpectedReference::Profile)
        }
        SolidOperation::SweptDisk { directrix, .. } => {
            expect_reference(nodes, *directrix, ExpectedReference::Curve)
        }
        SolidOperation::FixedReferenceSweep {
            profile, directrix, ..
        } => {
            expect_reference(nodes, *profile, ExpectedReference::Profile)?;
            expect_reference(nodes, *directrix, ExpectedReference::Curve)
        }
        SolidOperation::SurfaceCurveSweep {
            profile,
            directrix,
            reference_surface,
            ..
        } => {
            expect_reference(nodes, *profile, ExpectedReference::Profile)?;
            expect_reference(nodes, *directrix, ExpectedReference::Curve)?;
            expect_reference(nodes, *reference_surface, ExpectedReference::Surface)
        }
        SolidOperation::SectionedSpine { spine, sections } => {
            expect_reference(nodes, *spine, ExpectedReference::Curve)?;
            for section in sections {
                expect_reference(nodes, section.profile, ExpectedReference::Profile)?;
            }
            Ok(())
        }
        SolidOperation::Boolean { left, right, .. } => {
            expect_reference(nodes, *left, ExpectedReference::Solid)?;
            expect_reference(nodes, *right, ExpectedReference::Solid)
        }
        SolidOperation::BoundedHalfSpace {
            half_space,
            boundary,
            ..
        } => {
            expect_reference(nodes, *half_space, ExpectedReference::HalfSpace)?;
            expect_reference(nodes, *boundary, ExpectedReference::Curve)
        }
    }
}
