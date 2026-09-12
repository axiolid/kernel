use axiolid_contracts::GeomResult;
use axiolid_core::Interval;
use axiolid_nurbs::TransverseSurfaceSurfaceTrace3;
use axiolid_surface::{BSplineSurface, Surface};

use crate::trimmed_intersection_builder::ArrangementBuilder;
use crate::trimmed_intersection_classify::{DualClassification, SplitClassification};
use crate::trimmed_intersection_clone_surface::clone_surface;
use crate::trimmed_intersection_rectangle::{add_split_faces, add_unsplit_face};
use crate::trimmed_intersection_types::{
    CertifiedDualTrimmedSurfacePair3, CertifiedTrimmedSurfacePair3, EmbeddedFaceCurve,
    SurfacePairMember,
};

pub(super) fn assemble(
    first: &BSplineSurface,
    second: &BSplineSurface,
    trace: TransverseSurfaceSurfaceTrace3,
    classification: SplitClassification,
    residual_upper_bound: f64,
    visited_patch_pairs: u32,
    boundary_queries: u8,
) -> GeomResult<CertifiedTrimmedSurfacePair3> {
    let mut builder = ArrangementBuilder::new()?;
    let first_support = builder.add_surface(Surface::BSpline(clone_surface(first)?));
    let second_support = builder.add_surface(Surface::BSpline(clone_surface(second)?));

    let start_vertex = builder.add_vertex(trace.start.point);
    let end_vertex = builder.add_vertex(trace.end.point);
    let intersection_edge =
        builder.add_line_edge(start_vertex, end_vertex, trace.start.point, trace.end.point)?;

    let (owner_surface, owner_support, embedded_surface, embedded_support) =
        match classification.member {
            SurfacePairMember::First => (first, first_support, second, second_support),
            SurfacePairMember::Second => (second, second_support, first, first_support),
        };

    let split = add_split_faces(
        &mut builder,
        owner_surface,
        owner_support,
        classification.owner_domain,
        classification.owner_start,
        classification.owner_end,
        trace.start.point,
        trace.end.point,
        start_vertex,
        end_vertex,
        intersection_edge,
    )?;
    let embedded_face = add_unsplit_face(
        &mut builder,
        embedded_surface,
        embedded_support,
        classification.embedded_domain,
    )?;
    let embedded_pcurve = builder.add_pcurve(
        classification.embedded_start.uv,
        classification.embedded_end.uv,
    )?;
    let embedded_curve = EmbeddedFaceCurve {
        face: embedded_face,
        edge: intersection_edge,
        pcurve: embedded_pcurve,
        interval: Interval::UNIT,
    };
    let brep = builder.finish()?;

    Ok(CertifiedTrimmedSurfacePair3 {
        brep,
        trace,
        intersection_edge,
        split_surface: classification.member,
        split_faces: split.faces,
        unsplit_face: embedded_face,
        embedded_curve,
        max_surface_residual_upper_bound: residual_upper_bound,
        visited_patch_pairs,
        boundary_queries: u32::from(boundary_queries),
    })
}

/// Assemble four closed faces for a chord that partitions both patches.
///
/// The same model-space edge and its two vertices are shared by all four
/// loops: creating them once is what makes the result a single connected
/// arrangement rather than two arrangements that merely coincide.
pub(super) fn assemble_dual(
    first: &BSplineSurface,
    second: &BSplineSurface,
    trace: TransverseSurfaceSurfaceTrace3,
    classification: DualClassification,
    residual_upper_bound: f64,
    visited_patch_pairs: u32,
    boundary_queries: u8,
) -> GeomResult<CertifiedDualTrimmedSurfacePair3> {
    let mut builder = ArrangementBuilder::new()?;
    let first_support = builder.add_surface(Surface::BSpline(clone_surface(first)?));
    let second_support = builder.add_surface(Surface::BSpline(clone_surface(second)?));

    let start_vertex = builder.add_vertex(trace.start.point);
    let end_vertex = builder.add_vertex(trace.end.point);
    let intersection_edge =
        builder.add_line_edge(start_vertex, end_vertex, trace.start.point, trace.end.point)?;

    let first_split = add_split_faces(
        &mut builder,
        first,
        first_support,
        classification.first_domain,
        classification.first_start,
        classification.first_end,
        trace.start.point,
        trace.end.point,
        start_vertex,
        end_vertex,
        intersection_edge,
    )?;
    let second_split = add_split_faces(
        &mut builder,
        second,
        second_support,
        classification.second_domain,
        classification.second_start,
        classification.second_end,
        trace.start.point,
        trace.end.point,
        start_vertex,
        end_vertex,
        intersection_edge,
    )?;
    let brep = builder.finish()?;

    Ok(CertifiedDualTrimmedSurfacePair3 {
        brep,
        intersection_edge,
        first_faces: first_split.faces,
        second_faces: second_split.faces,
        trace,
        max_surface_residual_upper_bound: residual_upper_bound,
        visited_patch_pairs,
        boundary_queries: u32::from(boundary_queries),
    })
}
