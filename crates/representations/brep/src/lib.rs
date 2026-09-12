#![forbid(unsafe_code)]

//! Strict exact B-rep results over neutral analytic curve and surface values.
//!
//! [`axiolid_topology::BRep`] remains generic so graph/import clients can link
//! their own handles. [`ExactBRep`] binds that graph to owned `Curve3`, `Curve2`,
//! and `Surface` catalogs, then requires every support and trim span explicitly.
//! It does not evaluate, intersect, or tessellate geometry.

use std::collections::HashMap;
use std::fmt;

use axiolid_core::Interval;
use axiolid_curve::{Curve2, Curve3};
use axiolid_surface::Surface;
use axiolid_topology::{audit_brep, BRep, BRepHealth, EdgeId, FaceId, LoopId};

/// Persistent structural names for faces and edges.
pub mod name;

pub use name::{EdgeName, FaceName, Operand, SweptFace};

macro_rules! geometry_id {
    ($name:ident, $label:literal) => {
        #[doc = concat!("Typed handle into the exact B-rep ", $label, " catalog.")]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(u32);

        impl $name {
            fn from_index(index: usize) -> Self {
                Self(u32::try_from(index).expect("exact B-rep catalog exceeds u32 capacity"))
            }

            /// Zero-based catalog index.
            pub const fn index(self) -> usize {
                self.0 as usize
            }
        }
    };
}

geometry_id!(Curve3Id, "3D curve");
geometry_id!(Curve2Id, "2D pcurve");
geometry_id!(SurfaceId, "surface");

/// Topology whose supports are typed references into [`ExactBRep`] catalogs.
pub type ExactTopology = BRep<Curve3Id, Curve2Id, SurfaceId>;

/// Owned, analytic boundary representation.
///
/// Every edge and pcurve use has a finite, non-zero native parameter span. An
/// interval is deliberately owned here rather than inferred from endpoints:
/// periodic supports, reversed pcurves, and parameter re-mapping are not
/// recoverable from coordinates. The result may be an open sheet; an eventual
/// exact-solid contract will additionally require a selected closed shell.
#[derive(Debug, Clone, PartialEq)]
pub struct ExactBRep {
    topology: ExactTopology,
    curves3: Vec<Curve3>,
    curves2: Vec<Curve2>,
    surfaces: Vec<Surface>,
    edge_intervals: HashMap<EdgeId, Interval>,
    pcurve_intervals: HashMap<(LoopId, usize), Interval>,
    face_names: HashMap<FaceId, FaceName>,
    edge_names: HashMap<EdgeId, EdgeName>,
}

impl ExactBRep {
    /// Persistent structural name for a face, when the producer supplied one.
    ///
    /// Absence means the producing operation did not track provenance, not
    /// that the face is invalid. Callers that require a name should treat
    /// `None` as a refusal rather than substituting an index.
    pub fn face_name(&self, face: FaceId) -> Option<&FaceName> {
        self.face_names.get(&face)
    }

    /// Persistent structural name for an edge, when the producer supplied one.
    pub fn edge_name(&self, edge: EdgeId) -> Option<&EdgeName> {
        self.edge_names.get(&edge)
    }

    /// Find the face carrying `name`.
    ///
    /// This is the lookup that makes names useful across operations: a caller
    /// holding a name from before an edit resolves it against the new result
    /// instead of hoping an arena index still points at the same face.
    pub fn face_by_name(&self, name: &FaceName) -> Option<FaceId> {
        self.face_names
            .iter()
            .find(|(_, candidate)| *candidate == name)
            .map(|(id, _)| *id)
    }

    /// Find the edge carrying `name`.
    pub fn edge_by_name(&self, name: &EdgeName) -> Option<EdgeId> {
        self.edge_names
            .iter()
            .find(|(_, candidate)| *candidate == name)
            .map(|(id, _)| *id)
    }
    /// Structural topology with typed support handles.
    pub fn topology(&self) -> &ExactTopology {
        &self.topology
    }

    /// Owned exact 3D curve supports.
    pub fn curves3(&self) -> &[Curve3] {
        &self.curves3
    }

    /// Owned exact 2D trim-curve supports.
    pub fn curves2(&self) -> &[Curve2] {
        &self.curves2
    }

    /// Owned exact support surfaces.
    pub fn surfaces(&self) -> &[Surface] {
        &self.surfaces
    }

    /// Native parameter span oriented from an edge's start vertex to its end.
    pub fn edge_interval(&self, edge: EdgeId) -> Option<Interval> {
        self.edge_intervals.get(&edge).copied()
    }

    /// Native parameter span for an edge use's pcurve in loop traversal order.
    pub fn pcurve_interval(&self, loop_id: LoopId, use_index: usize) -> Option<Interval> {
        self.pcurve_intervals.get(&(loop_id, use_index)).copied()
    }
}

/// Mutable assembly state that can only yield an [`ExactBRep`] after validation.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ExactBRepBuilder {
    topology: ExactTopology,
    curves3: Vec<Curve3>,
    curves2: Vec<Curve2>,
    surfaces: Vec<Surface>,
    edge_intervals: HashMap<EdgeId, Interval>,
    pcurve_intervals: HashMap<(LoopId, usize), Interval>,
    face_names: HashMap<FaceId, FaceName>,
    edge_names: HashMap<EdgeId, EdgeName>,
}

impl ExactBRepBuilder {
    /// Record the persistent structural name of a face.
    ///
    /// Naming is opt-in per producer: an operation that cannot honestly say
    /// where a face came from simply does not call this, and the face reports
    /// no name rather than a misleading one.
    pub fn set_face_name(&mut self, face: FaceId, name: FaceName) {
        self.face_names.insert(face, name);
    }

    /// Record the persistent structural name of an edge.
    pub fn set_edge_name(&mut self, edge: EdgeId, name: EdgeName) {
        self.edge_names.insert(edge, name);
    }
    /// Fallibly reserve owned support and interval catalogs for bounded assembly.
    pub fn try_reserve(
        &mut self,
        curves3: usize,
        curves2: usize,
        surfaces: usize,
        edge_intervals: usize,
        pcurve_intervals: usize,
    ) -> Result<(), std::collections::TryReserveError> {
        self.curves3.try_reserve_exact(curves3)?;
        self.curves2.try_reserve_exact(curves2)?;
        self.surfaces.try_reserve_exact(surfaces)?;
        self.edge_intervals.try_reserve(edge_intervals)?;
        self.pcurve_intervals.try_reserve(pcurve_intervals)?;
        Ok(())
    }

    /// Store an exact 3D curve support.
    pub fn add_curve3(&mut self, curve: Curve3) -> Curve3Id {
        let id = Curve3Id::from_index(self.curves3.len());
        self.curves3.push(curve);
        id
    }

    /// Store an exact 2D pcurve support.
    pub fn add_curve2(&mut self, curve: Curve2) -> Curve2Id {
        let id = Curve2Id::from_index(self.curves2.len());
        self.curves2.push(curve);
        id
    }

    /// Store an exact surface support.
    pub fn add_surface(&mut self, surface: Surface) -> SurfaceId {
        let id = SurfaceId::from_index(self.surfaces.len());
        self.surfaces.push(surface);
        id
    }

    /// Mutable typed topology under construction. [`Self::finish`] validates it.
    pub fn topology_mut(&mut self) -> &mut ExactTopology {
        &mut self.topology
    }

    /// State an edge's finite native curve span.
    pub fn set_edge_interval(&mut self, edge: EdgeId, interval: Interval) {
        self.edge_intervals.insert(edge, interval);
    }

    /// State an edge use's finite native pcurve span in its loop traversal order.
    pub fn set_pcurve_interval(&mut self, loop_id: LoopId, use_index: usize, interval: Interval) {
        self.pcurve_intervals.insert((loop_id, use_index), interval);
    }

    /// Validate and freeze the exact B-rep result.
    pub fn finish(self) -> Result<ExactBRep, ExactBRepError> {
        validate(&self)?;
        Ok(ExactBRep {
            topology: self.topology,
            curves3: self.curves3,
            curves2: self.curves2,
            surfaces: self.surfaces,
            edge_intervals: self.edge_intervals,
            pcurve_intervals: self.pcurve_intervals,
            face_names: self.face_names,
            edge_names: self.edge_names,
        })
    }
}

/// Why exact B-rep assembly was refused.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExactBRepError {
    /// A result without faces is not a boundary representation.
    Empty,
    /// Generic topology has unresolved handles or invalid loop/face structure.
    Topology(BRepHealth),
    /// An edge did not state its exact three-dimensional support.
    MissingEdgeCurve { edge_index: usize },
    /// An edge support handle does not resolve in the 3D curve catalog.
    UnknownCurve3 { edge_index: usize },
    /// An edge did not state its native support-curve interval.
    MissingEdgeInterval { edge_index: usize },
    /// An edge support interval was non-finite or zero-length.
    InvalidEdgeInterval { edge_index: usize },
    /// An edge use did not state its trim curve in the owning face's parameters.
    MissingPcurve { loop_index: usize, use_index: usize },
    /// A pcurve support handle does not resolve in the 2D curve catalog.
    UnknownCurve2 { loop_index: usize, use_index: usize },
    /// A pcurve did not state its native interval.
    MissingPcurveInterval { loop_index: usize, use_index: usize },
    /// A pcurve interval was non-finite or zero-length.
    InvalidPcurveInterval { loop_index: usize, use_index: usize },
    /// A face did not state its exact support surface.
    MissingFaceSurface { face_index: usize },
    /// A face support handle does not resolve in the surface catalog.
    UnknownSurface { face_index: usize },
}

impl fmt::Display for ExactBRepError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("exact B-rep requires at least one face"),
            Self::Topology(_) => f.write_str("exact B-rep topology is structurally invalid"),
            Self::MissingEdgeCurve { edge_index } => {
                write!(f, "edge {edge_index} lacks a 3D curve")
            }
            Self::UnknownCurve3 { edge_index } => {
                write!(f, "edge {edge_index} has an unknown 3D curve")
            }
            Self::MissingEdgeInterval { edge_index } => {
                write!(f, "edge {edge_index} lacks a curve interval")
            }
            Self::InvalidEdgeInterval { edge_index } => {
                write!(f, "edge {edge_index} has an invalid curve interval")
            }
            Self::MissingPcurve {
                loop_index,
                use_index,
            } => write!(f, "loop {loop_index} use {use_index} lacks a pcurve"),
            Self::UnknownCurve2 {
                loop_index,
                use_index,
            } => write!(f, "loop {loop_index} use {use_index} has an unknown pcurve"),
            Self::MissingPcurveInterval {
                loop_index,
                use_index,
            } => write!(
                f,
                "loop {loop_index} use {use_index} lacks a pcurve interval"
            ),
            Self::InvalidPcurveInterval {
                loop_index,
                use_index,
            } => write!(
                f,
                "loop {loop_index} use {use_index} has an invalid pcurve interval"
            ),
            Self::MissingFaceSurface { face_index } => {
                write!(f, "face {face_index} lacks a support surface")
            }
            Self::UnknownSurface { face_index } => {
                write!(f, "face {face_index} has an unknown support surface")
            }
        }
    }
}

impl std::error::Error for ExactBRepError {}

fn validate(value: &ExactBRepBuilder) -> Result<(), ExactBRepError> {
    if value.topology.faces().is_empty() {
        return Err(ExactBRepError::Empty);
    }
    let health = audit_brep(&value.topology);
    if !health.is_tessellable() {
        return Err(ExactBRepError::Topology(health));
    }
    for (edge_index, edge) in value.topology.edges().iter().enumerate() {
        let Some(curve) = edge.curve else {
            return Err(ExactBRepError::MissingEdgeCurve { edge_index });
        };
        if curve.index() >= value.curves3.len() {
            return Err(ExactBRepError::UnknownCurve3 { edge_index });
        }
        let Some(edge_id) = value.topology.edge_id_at(edge_index) else {
            return Err(ExactBRepError::MissingEdgeInterval { edge_index });
        };
        let Some(interval) = value.edge_intervals.get(&edge_id) else {
            return Err(ExactBRepError::MissingEdgeInterval { edge_index });
        };
        if !valid_interval(*interval) {
            return Err(ExactBRepError::InvalidEdgeInterval { edge_index });
        }
    }
    for (loop_index, loop_) in value.topology.loops().iter().enumerate() {
        let Some(loop_id) = value.topology.loop_id_at(loop_index) else {
            return Err(ExactBRepError::Topology(health));
        };
        for (use_index, use_) in loop_.edges.iter().enumerate() {
            let Some(curve) = use_.pcurve else {
                return Err(ExactBRepError::MissingPcurve {
                    loop_index,
                    use_index,
                });
            };
            if curve.index() >= value.curves2.len() {
                return Err(ExactBRepError::UnknownCurve2 {
                    loop_index,
                    use_index,
                });
            }
            let Some(interval) = value.pcurve_intervals.get(&(loop_id, use_index)) else {
                return Err(ExactBRepError::MissingPcurveInterval {
                    loop_index,
                    use_index,
                });
            };
            if !valid_interval(*interval) {
                return Err(ExactBRepError::InvalidPcurveInterval {
                    loop_index,
                    use_index,
                });
            }
        }
    }
    for (face_index, face) in value.topology.faces().iter().enumerate() {
        let Some(surface) = face.surface else {
            return Err(ExactBRepError::MissingFaceSurface { face_index });
        };
        if surface.index() >= value.surfaces.len() {
            return Err(ExactBRepError::UnknownSurface { face_index });
        }
    }
    Ok(())
}

fn valid_interval(interval: Interval) -> bool {
    interval.start.is_finite() && interval.end.is_finite() && interval.length() > 0.0
}
