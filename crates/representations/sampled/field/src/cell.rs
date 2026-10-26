//! One `(x, y)` column of the field: surface crossings and occupancy spans.
//!
//! A triangle is a zero-thickness surface, so triangle coverage produces
//! [`SurfaceHit`]s. Positive-length [`Interval`]s are only produced by an
//! explicit occupancy construction over a closed shell. Keeping the two
//! channels separate stops a single facet from being mistaken for filled space.

use axiolid_core::{Interval, Scalar, Tolerance};

use crate::LayeredFieldError;

/// Orientation of a surface crossing relative to the field's local `z` axis.
///
/// This is a geometric fact about winding, not a claim that a facet is a floor,
/// a ceiling, or anything walkable. NOT-A-VERDICT
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SurfaceFacing {
    /// The triangle normal opposes the sampling direction (an entering crossing).
    AgainstNormal,
    /// The triangle normal agrees with the sampling direction (an exiting crossing).
    WithNormal,
}

/// A zero-thickness crossing of a sampling line at local coordinate `w`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceHit {
    w: Scalar,
    facing: SurfaceFacing,
}

impl SurfaceHit {
    /// Construct a crossing. `w` is validated when the cell is built.
    pub const fn new(w: Scalar, facing: SurfaceFacing) -> Self {
        Self { w, facing }
    }

    /// Layer coordinate along the field's local `z` axis.
    pub const fn w(self) -> Scalar {
        self.w
    }

    /// Winding orientation of the crossing.
    pub const fn facing(self) -> SurfaceFacing {
        self.facing
    }
}

/// Ordered surface crossings plus sorted, strictly disjoint occupancy spans.
///
/// Both channels may be empty, hold one entry, or hold many; a cell never
/// selects a "primary" layer.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct LayeredCell {
    surfaces: Vec<SurfaceHit>,
    occupancy: Vec<Interval>,
}

impl LayeredCell {
    /// An empty column.
    pub const fn empty() -> Self {
        Self {
            surfaces: Vec::new(),
            occupancy: Vec::new(),
        }
    }

    /// Construct an occupancy-only column.
    pub fn new(occupancy: Vec<Interval>) -> Result<Self, LayeredFieldError> {
        Self::with_layers(Vec::new(), occupancy)
    }

    /// Construct a validated column with both channels.
    ///
    /// Surfaces are sorted by `w`, then by facing so equal-coordinate ties are
    /// deterministic. Occupancy is sorted and must be strictly disjoint;
    /// touching or overlapping spans are rejected rather than merged, because a
    /// merge would silently erase a topology fact the caller may need.
    pub fn with_layers(
        mut surfaces: Vec<SurfaceHit>,
        mut occupancy: Vec<Interval>,
    ) -> Result<Self, LayeredFieldError> {
        if surfaces.iter().any(|hit| !hit.w.is_finite()) {
            return Err(LayeredFieldError::InvalidInterval);
        }
        if occupancy
            .iter()
            .any(|span| !span.start.is_finite() || !span.end.is_finite() || span.start >= span.end)
        {
            return Err(LayeredFieldError::InvalidInterval);
        }
        surfaces.sort_by(|left, right| {
            left.w
                .total_cmp(&right.w)
                .then_with(|| left.facing.cmp(&right.facing))
        });
        occupancy.sort_by(|left, right| {
            left.start
                .total_cmp(&right.start)
                .then_with(|| left.end.total_cmp(&right.end))
        });
        if occupancy
            .windows(2)
            .any(|pair| pair[0].end >= pair[1].start)
        {
            return Err(LayeredFieldError::NonDisjointIntervals);
        }
        Ok(Self {
            surfaces,
            occupancy,
        })
    }

    /// Crossings in increasing layer order.
    pub fn surfaces(&self) -> &[SurfaceHit] {
        &self.surfaces
    }

    /// Occupied spans in increasing layer order.
    pub fn occupancy(&self) -> &[Interval] {
        &self.occupancy
    }

    /// Total stored layers in both channels.
    pub fn layer_count(&self) -> usize {
        self.surfaces.len() + self.occupancy.len()
    }

    /// Whether both channels are empty.
    pub fn is_empty(&self) -> bool {
        self.surfaces.is_empty() && self.occupancy.is_empty()
    }

    /// Pair alternating crossings into occupancy spans for a closed shell.
    ///
    /// The sequence must alternate `AgainstNormal` (enter) then `WithNormal`
    /// (exit); anything else is reported as [`LayeredFieldError::UnbalancedCrossings`]
    /// rather than guessed. A pair whose span is within the linear tolerance is
    /// reported as [`LayeredFieldError::DegenerateOccupancy`].
    pub fn derive_occupancy(&self, tolerance: Tolerance) -> Result<Self, LayeredFieldError> {
        if !self.surfaces.len().is_multiple_of(2) {
            return Err(LayeredFieldError::UnbalancedCrossings);
        }
        let mut spans = Vec::with_capacity(self.surfaces.len() / 2);
        for pair in self.surfaces.as_chunks::<2>().0 {
            if pair[0].facing != SurfaceFacing::AgainstNormal
                || pair[1].facing != SurfaceFacing::WithNormal
            {
                return Err(LayeredFieldError::UnbalancedCrossings);
            }
            if (pair[1].w - pair[0].w).abs() <= tolerance.linear() {
                return Err(LayeredFieldError::DegenerateOccupancy);
            }
            spans.push(Interval::new(pair[0].w, pair[1].w));
        }
        Self::with_layers(self.surfaces.clone(), spans)
    }

    /// Largest unoccupied span strictly inside `search`, or `None` when the
    /// window is fully occupied.
    ///
    /// This reports a distance. It does not decide whether that distance is
    /// sufficient for any purpose.
    pub fn largest_free_span(&self, search: Interval) -> Option<Interval> {
        let (low, high) = if search.start <= search.end {
            (search.start, search.end)
        } else {
            (search.end, search.start)
        };
        let mut cursor = low;
        let mut best: Option<Interval> = None;
        for span in &self.occupancy {
            let start = span.start.max(low);
            let end = span.end.min(high);
            if start >= end {
                continue;
            }
            if start > cursor {
                best = keep_longer(best, Interval::new(cursor, start));
            }
            cursor = cursor.max(end);
        }
        if cursor < high {
            best = keep_longer(best, Interval::new(cursor, high));
        }
        best
    }
}

fn keep_longer(best: Option<Interval>, candidate: Interval) -> Option<Interval> {
    match best {
        Some(current) if current.length() >= candidate.length() => Some(current),
        _ => Some(candidate),
    }
}
