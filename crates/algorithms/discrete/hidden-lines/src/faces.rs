//! Exact visible face regions of meshes in an orthographic view (ADR 0094).
//!
//! The image is partitioned by an exact arrangement of segments, each piece
//! is labelled on both sides with the sources frontmost beside it, and the
//! pieces that separate different labels bound the regions.
//!
//! 1. `prepare` (shared with hidden lines) audits the meshes and builds the
//!    exact view-space vertices, triangles with their image orientation, and
//!    the cut solids. A triangle is *removed* when, decided exactly, it lies
//!    wholly beyond `near`, wholly beyond `far`, or its image misses the crop
//!    box; a removed triangle stays in the broad phase because cap parity
//!    counts the triangles beyond `near`.
//! 2. `segments` and `intersect` give the arrangement's input: the visible
//!    pieces of mesh edges (the hidden-line resolver with edge-on triangles
//!    not occluding, since they cover no area), the near and far clip
//!    segments, the cross-section segments of interpenetrating triangles, and
//!    the crop box. Hidden pieces of edges are dropped: a hidden piece lies
//!    under a frontmost surface that either covers both of its sides, so both
//!    carry its key, or covers one side only, so the piece lies on that
//!    surface's own visible edge, which is kept.
//! 3. `axiolid-exact-arrangement` cuts the segments into pieces.
//! 4. `label` decides, for each side of each piece, the frontmost sources at
//!    a point just beside it (`side`): an infinitesimal point
//!    `A + ε·d + ε²·n` whose affine signs are lexicographic and exact. Equal
//!    depth keeps every tied source; a triangle seen edge-on covers nothing;
//!    the cap of a cut solid is a surface at depth `near` found by ray
//!    parity.
//! 5. `regions` closes the boundary pieces of each source into exact rings,
//!    nests holes and attributes the frontmost triangles; then snap-rounds
//!    the edges of all rings together (`axiolid_exact_arrangement::snap_round`,
//!    ADR 0100), nets each source's boundary on the snapped fragments and
//!    traces it again into rings, which are certified.
//!
//! No tolerance enters any decision. Where no edge passes through the
//! rounding cell of another vertex, every vertex is an exact point rounded
//! once to nearest, so a vertex shared by two regions is bit-identical in
//! both and equals the hidden-line endpoint of the same point. Otherwise the
//! edge is routed through that cell's centre, a region that collapses is
//! dropped, and the evidence reports how far anything moved and how much
//! collapsed.

mod intersect;
mod label;
mod points;
mod provider;
mod regions;
mod segments;
mod side;
#[cfg(test)]
mod testing;

pub use provider::ScalarVisibleFaces;
#[doc(hidden)]
pub use provider::{ScanVisibleFaces, UnprunedVisibleFaces};

use std::cmp::Ordering;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::Aabb2;
use axiolid_exact_arrangement::PairBudget;
use axiolid_mesh_hidden_lines_contract::{FaceGroup, HiddenLineMesh, MeshId};

use crate::broad::Broad;
use crate::caps::Cap;
use crate::decide::decide;
use crate::edges::Tri;
use crate::groups::resolve_groups;
use crate::occlusion::{kept_constraints, Scene};
use crate::prepare::Prepared;
use crate::recipe::{Constraint, ValueSign};
use crate::view::Vert;
use side::{Form, Plane};

/// Counters shared by all steps of one request.
pub(crate) struct Work {
    pub budget: PairBudget,
    pub escalations: usize,
}

/// The prepared scene with what the face operation adds to it.
pub(crate) struct Faces<'a> {
    pub verts: &'a [Vert],
    pub tris: &'a [Tri],
    pub caps: &'a [Cap],
    pub broad: &'a Broad,
    pub near: Option<f64>,
    pub far: Option<f64>,
    pub region: Option<Aabb2>,
    /// Per triangle: wholly beyond a depth bound or outside the crop.
    pub removed: Vec<bool>,
    /// Per triangle: its plane, `None` when seen edge-on.
    pub planes: Vec<Option<Plane>>,
    /// Per triangle: the position of its mesh in the request.
    pub mesh_of: Vec<usize>,
    /// Per triangle: its face group.
    pub group_of: Vec<FaceGroup>,
    /// Per mesh: the scene id of its first vertex.
    pub bases: Vec<u32>,
    /// Per mesh: its id.
    pub ids: Vec<MeshId>,
    /// Per mesh: its first triangle in the scene's triangle list.
    pub firsts: Vec<usize>,
    /// The forms that are negative outside the crop box.
    pub crop: Vec<Form>,
}

impl<'a> Faces<'a> {
    pub(crate) fn new(
        prepared: &'a Prepared,
        broad: &'a Broad,
        meshes: &[HiddenLineMesh<'_>],
        (near, far, region): (Option<f64>, Option<f64>, Option<Aabb2>),
        work: &mut Work,
        options: &ExecutionOptions,
    ) -> Result<Self, Error> {
        let (verts, tris) = (&prepared.verts[..], &prepared.tris[..]);
        let mut mesh_of = Vec::with_capacity(tris.len());
        let mut group_of = Vec::with_capacity(tris.len());
        for (index, (entry, range)) in meshes.iter().zip(&prepared.ranges).enumerate() {
            mesh_of.resize(mesh_of.len() + range.triangles.len(), index);
            group_of.extend(resolve_groups(entry, options, &mut work.escalations)?);
        }
        let mut faces = Self {
            verts,
            tris,
            caps: &prepared.caps,
            broad,
            near,
            far,
            region,
            removed: Vec::with_capacity(tris.len()),
            planes: Vec::with_capacity(tris.len()),
            mesh_of,
            group_of,
            bases: prepared.ranges.iter().map(|r| r.base).collect(),
            ids: meshes.iter().map(|m| m.id).collect(),
            firsts: prepared.ranges.iter().map(|r| r.triangles.start).collect(),
            crop: region.map_or_else(Vec::new, |r| {
                [
                    (0, r.min.x, false),
                    (0, r.max.x, true),
                    (1, r.min.y, false),
                    (1, r.max.y, true),
                ]
                .map(|(axis, bound, upper)| Form::bound(axis, upper, bound))
                .to_vec()
            }),
        };
        let bounds = kept_constraints(&faces.scene());
        for (t, tri) in tris.iter().enumerate() {
            if t % 1024 == 0 {
                options.check_cancelled()?;
            }
            let id = faces.ids[faces.mesh_of[t]];
            let removed = Self::is_removed(tri, verts, &bounds, &mut work.escalations)
                .map_err(|error| Error::mesh(id, error))?;
            faces.removed.push(removed);
            faces.planes.push(match tri.sigma {
                Ordering::Equal => None,
                sigma => Some(Plane::new(verts, tri.v, sigma == Ordering::Less)),
            });
        }
        Ok(faces)
    }

    /// The occlusion scene with edge-on triangles not occluding.
    pub(crate) fn scene(&self) -> Scene<'_> {
        Scene {
            verts: self.verts,
            tris: self.tris,
            near: self.near,
            far: self.far,
            region: self.region,
            caps: self.caps,
            broad: self.broad,
            edge_on_occludes: false,
            count_see_through: false,
            perspective: false,
        }
    }

    /// Whether triangle `t` covers area that can be kept: not edge-on and not
    /// removed.
    pub(crate) fn covers(&self, t: usize) -> bool {
        self.planes[t].is_some() && !self.removed[t]
    }

    /// The error a failure in the geometry of triangle `t`'s mesh becomes.
    pub(crate) fn blame(&self, t: usize, error: GeomError) -> Error {
        match error {
            GeomError::Degenerate(_) => Error::mesh(self.ids[self.mesh_of[t]], error),
            other => Error::from(other),
        }
    }

    /// Decided exactly: all three vertices strictly outside one kept bound.
    fn is_removed(
        tri: &Tri,
        verts: &[Vert],
        bounds: &[Constraint],
        escalations: &mut usize,
    ) -> GeomResult<bool> {
        for constraint in bounds {
            let mut outside = true;
            for at in tri.v {
                let value = ValueSign {
                    verts,
                    constraint,
                    at,
                };
                if decide(&value, escalations)? != Ordering::Less {
                    outside = false;
                    break;
                }
            }
            if outside {
                return Ok(true);
            }
        }
        Ok(false)
    }
}

pub(crate) type Error = axiolid_mesh_hidden_lines_contract::HiddenLineError;
