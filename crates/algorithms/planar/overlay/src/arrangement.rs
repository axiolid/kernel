//! Planar subdivision by several arc rings at once (#120).
//!
//! [`arc_overlay`](crate::arc_overlay) answers one boolean of two rings.
//! A stepped or stacked solid needs more: the section changes with height,
//! and every face of the solid (each wall, ledge and cap) is bounded by
//! pieces of the same boundaries. [`ArcArrangement`] cuts the plane by all
//! the rings at once, once, and records for each piece which rings contain
//! the region on either side. Any region built from those pieces then
//! shares vertices with every other region by index, not by two roundings
//! happening to agree.
//!
//! # Exact, then rounded once
//!
//! Where boundaries cross, which pieces coincide, and which ring contains
//! what are exact decisions (ADR 0070). Crossing points are rounded to
//! `f64` once, into [`ArcArrangement::vertices`], and every piece refers to
//! them by index.

use axiolid_core::{Point2, Tolerance};

use crate::arc::{arc_ring_area, reverse_arc_ring, validate_arc_ring, ArcRing};
use crate::exact_arc::arrangement::{self, Raw};
use crate::OverlayError;

/// One piece of the subdivision: part of one or more input edges, running
/// between two vertices with no other ring's boundary crossing it.
#[derive(Debug, Clone, PartialEq)]
pub struct ArrangementEdge {
    /// Start vertex, an index into [`ArcArrangement::vertices`].
    pub from: usize,
    /// End vertex.
    pub to: usize,
    /// Bulge from `from` to `to`; `0` for a straight piece.
    pub bulge: f64,
    /// Rings whose boundary carries this piece.
    pub sources: Vec<EdgeSource>,
    left: Vec<bool>,
    right: Vec<bool>,
}

impl ArrangementEdge {
    /// Whether ring `ring` contains the region on the left of the piece.
    ///
    /// Left and right are taken along `from -> to`. A ring whose boundary
    /// carries the piece contains exactly one side.
    pub fn inside_left(&self, ring: usize) -> bool {
        self.left.get(ring).copied().unwrap_or(false)
    }

    /// Whether ring `ring` contains the region on the right of the piece.
    pub fn inside_right(&self, ring: usize) -> bool {
        self.right.get(ring).copied().unwrap_or(false)
    }
}

/// Where a piece of the subdivision came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdgeSource {
    /// Index of the input ring.
    pub ring: usize,
    /// Index of the edge in that ring, as the caller passed it: the edge
    /// leaving `vertices[edge]`.
    pub edge: usize,
    /// Whether the piece runs the same way as that input edge.
    pub forward: bool,
}

/// One use of a piece in a region boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdgeUse {
    /// Index into [`ArcArrangement::edges`].
    pub edge: usize,
    /// Whether the boundary traverses the piece from `to` to `from`.
    pub reversed: bool,
}

/// A region of the subdivision: an outer boundary and its holes.
///
/// The outer boundary runs counter-clockwise and each hole clockwise, so
/// the region always lies to the left of its boundary.
#[derive(Debug, Clone, PartialEq)]
pub struct ArrangementRegion {
    /// Outer boundary, in travel order.
    pub outer: Vec<EdgeUse>,
    /// Hole boundaries, in travel order.
    pub holes: Vec<Vec<EdgeUse>>,
}

/// The plane cut by several simple arc rings.
///
/// Built once, queried many times: [`Self::regions`] selects any set of
/// faces by a predicate over ring membership and links them into regions,
/// all over the same vertices.
pub struct ArcArrangement {
    raw: Raw,
    vertices: Vec<Point2>,
    edges: Vec<ArrangementEdge>,
    counts: Vec<usize>,
}

impl std::fmt::Debug for ArcArrangement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ArcArrangement")
            .field("vertices", &self.vertices)
            .field("edges", &self.edges)
            .finish_non_exhaustive()
    }
}

impl ArcArrangement {
    /// Subdivide the plane by `rings`.
    ///
    /// Each ring must pass [`validate_arc_ring`] and be simple; its winding
    /// does not matter (each is read as the region it encloses). Rings may
    /// cross, touch and share boundary pieces with each other.
    ///
    /// # Errors
    ///
    /// Any [`validate_arc_ring`] refusal, with the ring's own reason.
    pub fn new(rings: &[ArcRing], tolerance: Tolerance) -> Result<Self, OverlayError> {
        let mut oriented = Vec::with_capacity(rings.len());
        let mut reversed = Vec::with_capacity(rings.len());
        for ring in rings {
            validate_arc_ring(ring, tolerance)?;
            let flip = arc_ring_area(ring) < 0.0;
            oriented.push(if flip {
                reverse_arc_ring(ring)
            } else {
                ring.clone()
            });
            reversed.push(flip);
        }
        let raw = arrangement::build(&oriented);
        let counts: Vec<usize> = rings.iter().map(|ring| ring.vertices.len()).collect();
        let edges = raw
            .edges
            .iter()
            .map(|edge| ArrangementEdge {
                from: edge.from,
                to: edge.to,
                bulge: edge.bulge,
                sources: edge
                    .sources
                    .iter()
                    .map(|&(ring, index, forward)| {
                        // Reversal re-indexes edges: reversed edge `k`
                        // is original edge `n - 2 - k` (mod n), run the
                        // other way.
                        let n = counts[ring];
                        if reversed[ring] {
                            EdgeSource {
                                ring,
                                edge: (2 * n - 2 - index) % n,
                                forward: !forward,
                            }
                        } else {
                            EdgeSource {
                                ring,
                                edge: index,
                                forward,
                            }
                        }
                    })
                    .collect(),
                left: edge.left.clone(),
                right: edge.right.clone(),
            })
            .collect();
        Ok(Self {
            vertices: raw.vertex_positions(),
            raw,
            edges,
            counts,
        })
    }

    /// Distinct vertices, each an exact point rounded once.
    pub fn vertices(&self) -> &[Point2] {
        &self.vertices
    }

    /// Pieces of the subdivision.
    pub fn edges(&self) -> &[ArrangementEdge] {
        &self.edges
    }

    /// Number of input rings.
    pub fn ring_count(&self) -> usize {
        self.counts.len()
    }

    /// The regions where `inside` holds, as linked boundaries.
    ///
    /// `inside` receives one flag per input ring (whether a point lies in
    /// that ring) and says whether the point belongs to the wanted set. A
    /// piece bounds the set exactly when `inside` differs across it; it is
    /// traversed so that the set lies on its left.
    ///
    /// # Errors
    ///
    /// [`OverlayError::SelfIntersection`] if the boundary cannot be linked,
    /// which simple input rings cannot produce.
    pub fn regions(
        &self,
        inside: impl Fn(&[bool]) -> bool,
    ) -> Result<Vec<ArrangementRegion>, OverlayError> {
        let keep: Vec<Option<bool>> = self
            .edges
            .iter()
            .map(|edge| match (inside(&edge.left), inside(&edge.right)) {
                (true, false) => Some(false),
                (false, true) => Some(true),
                _ => None,
            })
            .collect();
        let uses = |ring: Vec<(usize, bool)>| -> Vec<EdgeUse> {
            ring.into_iter()
                .map(|(edge, reversed)| EdgeUse { edge, reversed })
                .collect()
        };
        Ok(self
            .raw
            .regions(&keep)?
            .into_iter()
            .map(|(outer, holes)| ArrangementRegion {
                outer: uses(outer),
                holes: holes.into_iter().map(uses).collect(),
            })
            .collect())
    }

    /// A region boundary as an [`ArcRing`] over the rounded vertices.
    pub fn ring(&self, uses: &[EdgeUse]) -> ArcRing {
        ArcRing::new(
            uses.iter()
                .map(|u| {
                    let edge = &self.edges[u.edge];
                    let (from, bulge) = if u.reversed {
                        (edge.to, -edge.bulge)
                    } else {
                        (edge.from, edge.bulge)
                    };
                    crate::arc::ArcVertex::bulged(self.vertices[from], bulge)
                })
                .collect(),
        )
    }
}
