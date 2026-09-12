#![forbid(unsafe_code)]

//! Geometry generation from exact inputs into explicit mesh or focused analytic results.
//!
//! # Why this is its own crate
//!
//! Everything here answers one question: given an exact profile and a path,
//! what solid does that denote? Extrusion, revolution, the sweep families,
//! lofting and half-space clipping are all the same problem with different
//! path kinds, and they share one stitching implementation so that winding,
//! hole orientation and cap pairing are decided once.
//!
//! None of that needs an operation graph. These functions take geometry and
//! return geometry; they do not walk a DAG, cache results, resolve node
//! references, or dispatch to a backend. That is why they live here and not
//! in `axiolid-mesh-compile`, which does all four.
//!
//! The split matters beyond tidiness. A caller that wants a swept solid --
//! a CAD front end, a test, a future exact B-rep generator -- should not have
//! to construct a `SolidOperation` graph and run a compiler to get one. Under
//! the old layout that was the only way to reach this code.
//!
//! # What this crate does not do
//!
//! Broad sweep/profile families still produce meshes by default. Exact generation is
//! intentionally narrow: sharp rectangular (including through-hole) and axial circular
//! extrusions populate all supports, pcurves, and spans, while unsupported families
//! refuse. Certified affine surface-pair arrangement remains the other focused analytic
//! constructor. Neither slice implies exact booleans or general sweeps.

use axiolid_contracts::BackendId;

/// Identity these generators report in diagnostics.
///
/// Distinct from the compiler's: a failure raised while building a swept
/// solid comes from this crate, and attributing it to `axiolid-mesh-compile`
/// would send a reader to the wrong place. Sweeps already reported a
/// separate identity before the split; this makes every generator
/// consistent with that.
pub const BACKEND_ID: BackendId = BackendId::new("scalar-generate");

pub mod boolean_exact;
mod boolean_provenance;
pub mod center_line;
pub mod extrude;
mod extrude_exact;
pub mod feature;
pub mod half_space;
pub mod hull;
pub mod loft;
pub mod offset;
pub mod polyhedron;
pub mod profile;
pub mod result;
pub mod revolve;
pub mod revolve_exact;
pub mod sweep;
pub mod trimmed_intersection;
mod trimmed_intersection_assembly;
mod trimmed_intersection_builder;
mod trimmed_intersection_classify;
mod trimmed_intersection_clone_surface;
mod trimmed_intersection_rectangle;
mod trimmed_intersection_types;

pub use axiolid_brep::{
    Curve2Id, Curve3Id, ExactBRep, ExactBRepBuilder, ExactBRepError, ExactTopology, SurfaceId,
};
pub use result::{GeneratedGeometry, GenerationOutput, GenerationRequest, TessellationRequest};
pub use trimmed_intersection::{
    split_surface_pair_certified, CertifiedDualTrimmedSurfacePair3, CertifiedSurfacePairSplit3,
    CertifiedSurfacePairSplitOptions, CertifiedTrimmedSurfacePair3, EmbeddedFaceCurve,
    SurfacePairMember, SurfacePairSplitUnresolvedReason,
};
