#![forbid(unsafe_code)]
#![warn(missing_docs)]

//! General NURBS algorithms over Axiolid's format-neutral B-spline values.
//!
//! This crate builds on the portable scalar oracle. It owns differential
//! geometry and exact shape-preserving transformations; importers and
//! tessellators are consumers, not the capability boundary.

mod axis;
mod certified_bezier;
mod certified_curve_distance;
mod certified_curve_intersection;
mod certified_curve_projection;
mod certified_curve_surface_intersection;
mod certified_projection;
mod certified_refinement;
mod certified_surface_arcs;
mod certified_surface_bezier;
mod certified_surface_inversion;
mod certified_surface_projection;
mod certified_surface_surface_intersection;
mod curve_analysis;
mod curve_projection;
mod degree;
mod fit;
mod intersection_curve;
mod periodic;
mod periodic_surface;
mod projection;
mod surface_analysis;
mod surface_projection;
mod surface_transform;
mod transform;

pub use certified_curve_distance::{distance_curve2_certified, distance_curve3_certified};
pub use certified_curve_intersection::{
    intersect_curve2_certified, CertifiedCurveIntersection2, CertifiedCurveIntersectionOptions,
    ClassifiedCurveContact2, CurveIntersectionDegeneracy, TransverseCurveIntersection2,
};
pub use certified_curve_projection::{project_curve2_certified, project_curve3_certified};
pub use certified_curve_surface_intersection::{
    intersect_curve_surface_certified, CertifiedCurveSurfaceIntersection3,
    CertifiedCurveSurfaceIntersectionOptions, CurveSurfaceParameterBox,
    TransverseCurveSurfaceIntersection3,
};
pub use certified_projection::{
    CertifiedProjectionOptions, CertifiedSurfaceProjection3, CertifiedSurfaceProjectionOptions,
    CurveDistanceCertificate2, CurveDistanceCertificate3, CurvePairParameterBox,
    CurveProjectionCertificate2, CurveProjectionCertificate3, ParameterInterval,
    SurfaceParameterBox, SurfaceProjectionCertificate3, SurfaceProjectionUnresolvedReason,
    MAX_CERTIFIED_SURFACE_PROJECTION_DEPTH, MAX_CERTIFIED_SURFACE_PROJECTION_WORK,
};
pub use certified_surface_arcs::{
    audit_coverage, certify_surface_arcs, CertifiedRegion, CertifiedSurfaceArcs3,
    CertifiedSurfaceArcsOptions, CoverageFault, RegionKind, MAX_AUDIT_DEPTH,
};
pub use certified_surface_inversion::{
    invert_periodic_surface_certified, invert_surface_certified, SurfaceInversionCertificate3,
    SurfaceInversionRefusal,
};
pub use certified_surface_projection::{
    project_periodic_surface_certified, project_surface_certified,
};
pub use certified_surface_surface_intersection::{
    intersect_surface_surface_certified, CertifiedSurfaceSurfaceIntersection3,
    CertifiedSurfaceSurfaceIntersectionOptions, SurfaceSurfaceParameterBox,
    SurfaceSurfaceTraceEndpoint3, TransverseSurfaceSurfaceTrace3,
};
pub use curve_analysis::{analyze_curve2, analyze_curve3, CurveDifferential2, CurveDifferential3};
pub use curve_projection::{project_curve2, project_curve3};
pub use degree::{
    elevate_degree2, elevate_degree3, reduce_degree2, reduce_degree3, remove_knot2, remove_knot3,
    BoundedResult,
};
pub use fit::{interpolate_curve3, loft_surface};
pub use intersection_curve::{
    construct_curve_surface_points, construct_surface_surface_curves,
    ConstructedCurveSurfacePoint3, ConstructedIntersectionCurve3, IntersectionCurveRefusal,
};
pub use periodic::{
    curve2_seam_continuity, curve3_seam_continuity, wrap_curve2_parameter, wrap_curve3_parameter,
    PeriodicCurve2, PeriodicCurve3, SeamContinuity,
};
pub use periodic_surface::PeriodicBSplineSurface;
pub use projection::{
    CurveProjection2, CurveProjection3, ProjectionOptions, ProjectionStatus, SurfaceProjection,
};
pub use surface_analysis::{analyze_surface, FundamentalForm, SurfaceDifferential};
pub use surface_projection::project_surface;
pub use surface_transform::{
    insert_surface_knot_u, insert_surface_knot_v, reverse_surface_u, reverse_surface_v,
};
pub use transform::{
    bezier_segments2, bezier_segments3, insert_knot2, insert_knot3, reverse2, reverse3, split2,
    split3,
};
