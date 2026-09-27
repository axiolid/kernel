//! A persistent planar region with set algebra, morphology and components.
//!
//! # Why a type rather than more free functions
//!
//! [`crate::overlay`] is stateless: it validates both operands on every
//! call and returns a fresh result. That is the right primitive, but a
//! consumer that builds free space, subtracts obstacles, erodes by a body
//! radius and then counts what is left has to thread raw polygon vectors
//! between calls and pay revalidation each time. Area and component count
//! are not reachable at all.
//!
//! A [`Region`] is validated once at construction. Every operation on it
//! consumes already-valid geometry and produces already-valid geometry, so
//! the invariant is established at the boundary rather than re-proven at
//! each step.
//!
//! # Emptiness is a reported fact, not an absence
//!
//! Every operation that can annihilate a region reports whether it did.
//! An erosion that removes the last of a corridor and an intersection of
//! two disjoint regions both yield no polygons, but they are different
//! facts, and a bare empty vector cannot tell a caller which happened.
//! [`RegionEvidence::emptied`] records that the operation consumed a
//! non-empty input, so "nothing was there" stays distinguishable from
//! "this operation destroyed it".
//!
//! # Boundary
//!
//! The kernel owns the region, the operations and the reported measures.
//! It does not own what the region represents, which radius is correct, or
//! whether a resulting area or component count is acceptable.

use axiolid_core::{Point2, Tolerance, Vec2};

use crate::minkowski::MorphologyBound;
use crate::offset::{offset_polygons, total_area, JoinStyle};
use crate::{canonical, validate_ring, OverlayError, Polygon, Ring};

/// What an operation did to produce a region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct RegionEvidence {
    /// Number of disjoint polygons in the result.
    pub polygons: usize,
    /// Number of inner boundary components across all polygons.
    pub holes: usize,
    /// The operation consumed a non-empty input and produced nothing.
    ///
    /// Distinguishes "eroded out of existence" or "disjoint operands" from
    /// "the input was already empty", which a bare empty result cannot.
    pub emptied: bool,
}

/// A validated planar region: a normalised set of polygons with holes.
///
/// Construction validates; operations preserve validity. Ring ordering is
/// canonical (outer counter-clockwise, holes clockwise, each rotated to its
/// lexicographically smallest vertex) and polygons are sorted, so equal
/// regions compare equal regardless of how they were built.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Region {
    polygons: Vec<Polygon>,
    evidence: RegionEvidence,
    /// For a disc morphology with a stated side, that side and the
    /// deviation from the exact result (#163).
    bound: Option<MorphologyBound>,
}

/// Normalise a polygon set: canonical ring order, deterministic sorting.
///
/// Applied to every operation result so that two regions holding the same
/// point set are structurally equal, which is what makes idempotence and
/// round-trip assertions meaningful rather than incidental.
fn normalise(mut polygons: Vec<Polygon>) -> Vec<Polygon> {
    polygons = polygons
        .into_iter()
        .map(|polygon| Polygon {
            outer: canonical(polygon.outer, true),
            holes: polygon
                .holes
                .into_iter()
                .map(|hole| canonical(hole, false))
                .collect(),
        })
        .collect();
    for polygon in &mut polygons {
        polygon.holes.sort_by(|a, b| {
            a.points[0]
                .x
                .total_cmp(&b.points[0].x)
                .then(a.points[0].y.total_cmp(&b.points[0].y))
        });
    }
    polygons.sort_by(|a, b| {
        a.outer.points[0]
            .x
            .total_cmp(&b.outer.points[0].x)
            .then(a.outer.points[0].y.total_cmp(&b.outer.points[0].y))
    });
    polygons
}

impl Region {
    /// The empty region. Not a failure: the identity for union.
    #[must_use]
    pub fn empty() -> Self {
        Self::default()
    }

    /// Build a region from polygons, validating every ring.
    ///
    /// Validation is the same used by [`crate::overlay`], so the two cannot
    /// disagree about what a well-formed polygon is.
    pub fn new(polygons: Vec<Polygon>, tolerance: Tolerance) -> Result<Self, OverlayError> {
        for polygon in &polygons {
            validate_ring(&polygon.outer, tolerance)?;
            for hole in &polygon.holes {
                validate_ring(hole, tolerance)?;
            }
        }
        Ok(Self::from_valid(normalise(polygons), false))
    }

    /// Wrap already-valid geometry, recording whether the producing operation
    /// annihilated a non-empty input.
    fn from_valid(polygons: Vec<Polygon>, had_input: bool) -> Self {
        let evidence = RegionEvidence {
            polygons: polygons.len(),
            holes: polygons.iter().map(|p| p.holes.len()).sum(),
            emptied: had_input && polygons.is_empty(),
        };
        Self {
            polygons: normalise(polygons),
            evidence,
            bound: None,
        }
    }

    /// [`Self::from_valid`] for the other modules of the crate.
    pub(crate) fn from_valid_polygons(polygons: Vec<Polygon>, had_input: bool) -> Self {
        Self::from_valid(polygons, had_input)
    }

    /// Record the side and deviation of a disc morphology.
    pub(crate) fn set_bound(&mut self, bound: MorphologyBound) {
        self.bound = Some(bound);
    }

    /// For a result of [`Self::dilate_inner`], [`Self::dilate_outer`],
    /// [`Self::erode_inner`] or [`Self::erode_outer`]: which side of the
    /// exact disc morphology it lies on and how far from it it can be.
    /// `None` for every other region, including [`Self::dilate`] and
    /// [`Self::erode`], whose side is not stated.
    #[must_use]
    pub const fn bound(&self) -> Option<MorphologyBound> {
        self.bound
    }

    /// The region's polygons, in canonical order.
    #[must_use]
    pub fn polygons(&self) -> &[Polygon] {
        &self.polygons
    }

    /// What the producing operation did.
    #[must_use]
    pub const fn evidence(&self) -> RegionEvidence {
        self.evidence
    }

    /// True when the region covers no area.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.polygons.is_empty()
    }

    /// Total covered area. Holes subtract.
    #[must_use]
    pub fn area(&self) -> f64 {
        total_area(&self.polygons)
    }

    /// Number of disjoint connected components.
    ///
    /// One polygon is one component: the overlay backend already resolves
    /// touching and overlapping input into disjoint output polygons, so
    /// counting them is the component count rather than an approximation
    /// of it.
    #[must_use]
    pub fn component_count(&self) -> usize {
        self.polygons.len()
    }

    /// Every boundary ring, outer boundaries first then holes.
    #[must_use]
    pub fn boundary_rings(&self) -> Vec<Ring> {
        let mut rings: Vec<Ring> = self.polygons.iter().map(|p| p.outer.clone()).collect();
        for polygon in &self.polygons {
            rings.extend(polygon.holes.iter().cloned());
        }
        rings
    }

    /// Set union.
    pub fn union(&self, other: &Self, tolerance: Tolerance) -> Result<Self, OverlayError> {
        self.combine(other, crate::OverlayOperation::Union, tolerance)
    }

    /// Set intersection.
    pub fn intersection(&self, other: &Self, tolerance: Tolerance) -> Result<Self, OverlayError> {
        self.combine(other, crate::OverlayOperation::Intersection, tolerance)
    }

    /// Set difference: this region minus `other`.
    pub fn difference(&self, other: &Self, tolerance: Tolerance) -> Result<Self, OverlayError> {
        self.combine(other, crate::OverlayOperation::Difference, tolerance)
    }

    /// Apply a boolean operation through the validated overlay primitive.
    ///
    /// Empty operands are handled here rather than pushed into the backend:
    /// the identities are exact, and routing them through a general overlay
    /// would risk an incidental simplification pass changing geometry that
    /// set algebra says must be returned unchanged.
    fn combine(
        &self,
        other: &Self,
        operation: crate::OverlayOperation,
        tolerance: Tolerance,
    ) -> Result<Self, OverlayError> {
        let had_input = !self.is_empty() || !other.is_empty();
        match operation {
            crate::OverlayOperation::Union if other.is_empty() => {
                return Ok(Self::from_valid(self.polygons.clone(), had_input));
            }
            crate::OverlayOperation::Union if self.is_empty() => {
                return Ok(Self::from_valid(other.polygons.clone(), had_input));
            }
            crate::OverlayOperation::Intersection if self.is_empty() || other.is_empty() => {
                return Ok(Self::from_valid(Vec::new(), had_input));
            }
            crate::OverlayOperation::Difference if other.is_empty() => {
                return Ok(Self::from_valid(self.polygons.clone(), had_input));
            }
            crate::OverlayOperation::Difference if self.is_empty() => {
                return Ok(Self::from_valid(Vec::new(), had_input));
            }
            _ => {}
        }
        // The identity frame. A region carries no frame of its own: it is a
        // point set in whatever frame the caller established, and overlay
        // only requires both operands to agree.
        let frame = axiolid_core::Frame2 {
            origin: Point2::new(0.0, 0.0),
            x: Vec2::new(1.0, 0.0),
            y: Vec2::new(0.0, 1.0),
        };
        let subject = crate::OverlayInput {
            frame,
            polygons: self.polygons.clone(),
        };
        let clip = crate::OverlayInput {
            frame,
            polygons: other.polygons.clone(),
        };
        let result = crate::overlay(
            &subject,
            &clip,
            operation,
            crate::FillRule::NonZero,
            tolerance,
        )?;
        Ok(Self::from_valid(result.polygons, had_input))
    }

    /// Dilate by `radius`: the Minkowski sum with a disc, approximated by
    /// round joins on an unstated side of the exact result. Where a verdict
    /// must be proven, use [`Self::dilate_inner`] or [`Self::dilate_outer`].
    ///
    /// This is the disc-expansion form used for clearance envelopes. A zero
    /// radius is the identity.
    pub fn dilate(&self, radius: f64, tolerance: Tolerance) -> Result<Self, OverlayError> {
        self.morphology(radius, tolerance)
    }

    /// Erode by `radius`: the Minkowski erosion by a disc, approximated on
    /// an unstated side of the exact result. Where a verdict must be
    /// proven, use [`Self::erode_inner`] (a route found proves
    /// reachability) or [`Self::erode_outer`] (no route proves
    /// unreachability).
    ///
    /// A region thinner than `2 * radius` anywhere is cut there, which is
    /// how a corridor narrower than a body radius becomes impassable. If
    /// that removes everything the result is empty and
    /// [`RegionEvidence::emptied`] is set.
    pub fn erode(&self, radius: f64, tolerance: Tolerance) -> Result<Self, OverlayError> {
        if radius < 0.0 {
            return Err(OverlayError::InvalidOffsetDistance);
        }
        self.morphology(-radius, tolerance)
    }

    /// Shared offset path. Round joins approximate the disc, which is what
    /// makes this morphology rather than a polygonal offset.
    fn morphology(&self, distance: f64, tolerance: Tolerance) -> Result<Self, OverlayError> {
        if !distance.is_finite() {
            return Err(OverlayError::InvalidOffsetDistance);
        }
        if self.is_empty() || distance == 0.0 {
            return Ok(Self::from_valid(self.polygons.clone(), false));
        }
        // A disc is approximated by round joins; the parameter is the maximum
        // segment-length-to-radius ratio, so smaller is closer to a true disc.
        let join = JoinStyle::Round {
            max_segment_ratio: 0.1,
        };
        let result = offset_polygons(&self.polygons, distance, join, tolerance)?;
        Ok(Self::from_valid(result.polygons, true))
    }

    /// Translate by `offset`. Rigid: area and component count are preserved.
    pub fn translate(&self, offset: Vec2) -> Result<Self, OverlayError> {
        if !offset.is_finite() {
            return Err(OverlayError::NonFinitePoint);
        }
        let shift = |ring: &Ring| Ring {
            points: ring
                .points
                .iter()
                .map(|p| Point2::new(p.x + offset.x, p.y + offset.y))
                .collect(),
        };
        let polygons = self
            .polygons
            .iter()
            .map(|polygon| Polygon {
                outer: shift(&polygon.outer),
                holes: polygon.holes.iter().map(shift).collect(),
            })
            .collect();
        Ok(Self::from_valid(polygons, !self.is_empty()))
    }

    /// Sweep along `direction`: the union of the region with every
    /// translate of itself along the vector.
    ///
    /// Equivalent to the Minkowski sum with the segment `[0, direction]`.
    /// Built as the union of the region and its translate plus the hull
    /// swept between them, which for a polygon set is exactly the union of
    /// the two end positions with the stroke of each boundary edge; using
    /// the union of endpoints alone would miss the swept middle whenever
    /// the translation exceeds the region's own extent.
    pub fn sweep(&self, direction: Vec2, tolerance: Tolerance) -> Result<Self, OverlayError> {
        if !direction.is_finite() {
            return Err(OverlayError::NonFinitePoint);
        }
        if self.is_empty() || direction.length() == 0.0 {
            return Ok(Self::from_valid(self.polygons.clone(), false));
        }
        let moved = self.translate(direction)?;
        let mut swept = self.union(&moved, tolerance)?;
        // Fill the corridor between the two positions: each boundary edge
        // sweeps a parallelogram, and their union with the endpoints is the
        // Minkowski sum. Built from quads so no offset approximation enters.
        for ring in self.boundary_rings() {
            let count = ring.points.len();
            for index in 0..count {
                let a = ring.points[index];
                let b = ring.points[(index + 1) % count];
                let quad = Polygon {
                    outer: Ring {
                        points: vec![
                            a,
                            b,
                            Point2::new(b.x + direction.x, b.y + direction.y),
                            Point2::new(a.x + direction.x, a.y + direction.y),
                        ],
                    },
                    holes: Vec::new(),
                };
                // A degenerate quad (edge parallel to the sweep) contributes
                // nothing and is skipped rather than rejected.
                if let Ok(band) = Self::new(vec![quad], tolerance) {
                    swept = swept.union(&band, tolerance)?;
                }
            }
        }
        Ok(swept)
    }
}
