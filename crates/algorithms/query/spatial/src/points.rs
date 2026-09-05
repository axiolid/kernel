//! Nearest-neighbour and radius queries over point sets.
//!
//! # Exact, not broad phase
//!
//! The BVH in this crate answers with AABB *lower bounds*: a candidate may
//! be further away than its bound suggests, so its results are broad-phase
//! candidates that a narrow phase must confirm.
//!
//! A point has no extent, so the distance to it is exact. These queries
//! therefore return real distances and are complete: a radius query returns
//! every point inside the radius and nothing else. That difference is why
//! the result types here are distinct from [`NearestCandidate`] — a caller
//! must never mistake an exact hit for a candidate needing confirmation, or
//! the reverse.
//!
//! [`NearestCandidate`]: crate::NearestCandidate
//!
//! # Allocation
//!
//! Queries take a callback and allocate nothing per hit. The index itself
//! is built once and reused. `*_into` helpers exist for callers that do
//! want a vector, but they are a convenience over the callback form rather
//! than the primitive.
//!
//! # Determinism
//!
//! Ties are broken by point index, so equal distances always resolve the
//! same way. Sorted variants order by distance then index. Without this a
//! caller could get different neighbours across runs of identical input,
//! which would make any downstream reconstruction irreproducible.

use axiolid_core::{Point3, Scalar};

/// A point found by a query, with its exact distance.
///
/// Distinct from a broad-phase candidate: `distance` is measured, not
/// bounded, and needs no narrow-phase confirmation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointHit {
    /// Index of the point in the source cloud.
    pub index: usize,
    /// Exact distance from the query position.
    pub distance: Scalar,
}

/// Why a query could not be answered.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum PointQueryError {
    /// The query position is not finite.
    NonFiniteQuery,
    /// A radius is negative or not finite.
    InvalidRadius,
}

impl core::fmt::Display for PointQueryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFiniteQuery => f.write_str("query position is not finite"),
            Self::InvalidRadius => f.write_str("radius must be finite and non-negative"),
        }
    }
}

impl core::error::Error for PointQueryError {}

/// A uniform-grid index over a point set.
///
/// A grid rather than a k-d tree because scan data is dense and roughly
/// uniform: bucketing is O(n) to build with no comparisons, and a radius
/// query touches only the cells the sphere overlaps. A k-d tree wins on
/// wildly non-uniform data, which is why this type is deliberately not the
/// only shape the API could take — the query methods are what callers use,
/// so the structure can change without moving them.
#[derive(Debug, Clone)]
pub struct PointIndex {
    points: Vec<Point3>,
    /// Cell edge length. Always positive and finite once constructed.
    cell: Scalar,
    origin: Point3,
    /// Grid dimensions in cells.
    dims: [usize; 3],
    /// Start offset per cell into `ordered`, with a trailing total.
    starts: Vec<u32>,
    /// Point indices grouped by cell.
    ordered: Vec<u32>,
}

impl PointIndex {
    /// Build an index over the given positions.
    ///
    /// Non-finite positions are dropped from the index and reported by
    /// [`rejected`](Self::rejected); they are not silently treated as
    /// present at the origin, which would corrupt every query near it.
    pub fn build(points: &[Point3]) -> Self {
        let finite: Vec<u32> = points
            .iter()
            .enumerate()
            .filter(|(_, p)| p.is_finite())
            .map(|(index, _)| index as u32)
            .collect();

        if finite.is_empty() {
            return Self {
                points: points.to_vec(),
                cell: 1.0,
                origin: Point3::ZERO,
                dims: [1, 1, 1],
                starts: vec![0, 0],
                ordered: Vec::new(),
            };
        }

        let mut min = points[finite[0] as usize];
        let mut max = min;
        for &i in &finite {
            let p = points[i as usize];
            min = Point3::new(min.x.min(p.x), min.y.min(p.y), min.z.min(p.z));
            max = Point3::new(max.x.max(p.x), max.y.max(p.y), max.z.max(p.z));
        }
        let span = max - min;

        // Target roughly one point per cell: cube-root the count and divide
        // the extent by it. A degenerate extent (all points coincident, or
        // a planar sheet) collapses to a single cell rather than dividing
        // by zero.
        let target = (finite.len() as Scalar).cbrt().max(1.0);
        let longest = span.x.max(span.y).max(span.z);
        let cell = if longest > 0.0 {
            (longest / target).max(Scalar::MIN_POSITIVE)
        } else {
            1.0
        };

        let dims = [
            ((span.x / cell).ceil() as usize + 1).max(1),
            ((span.y / cell).ceil() as usize + 1).max(1),
            ((span.z / cell).ceil() as usize + 1).max(1),
        ];
        let cell_count = dims[0] * dims[1] * dims[2];

        // Counting sort into cells: one pass to count, one to place. No
        // per-cell Vec, so building costs two linear passes and one
        // allocation for each of the two arrays.
        let mut counts = vec![0u32; cell_count + 1];
        let locate = |p: Point3| -> usize {
            let ix = (((p.x - min.x) / cell) as usize).min(dims[0] - 1);
            let iy = (((p.y - min.y) / cell) as usize).min(dims[1] - 1);
            let iz = (((p.z - min.z) / cell) as usize).min(dims[2] - 1);
            (iz * dims[1] + iy) * dims[0] + ix
        };
        for &i in &finite {
            counts[locate(points[i as usize]) + 1] += 1;
        }
        for k in 1..counts.len() {
            counts[k] += counts[k - 1];
        }
        let starts = counts.clone();

        let mut cursor = counts;
        let mut ordered = vec![0u32; finite.len()];
        for &i in &finite {
            let cell_index = locate(points[i as usize]);
            ordered[cursor[cell_index] as usize] = i;
            cursor[cell_index] += 1;
        }

        Self {
            points: points.to_vec(),
            cell,
            origin: min,
            dims,
            starts,
            ordered,
        }
    }

    /// Number of positions the index was built over, including rejected.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether the index holds no positions.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Positions dropped because they were not finite.
    pub fn rejected(&self) -> usize {
        self.points.len() - self.ordered.len()
    }

    /// Visit every point within `radius` of `query`.
    ///
    /// Results are exact and complete. Visit order is unspecified; use
    /// [`radius_into`](Self::radius_into) when order matters.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite query position and a negative or non-finite
    /// radius, rather than returning an empty result that a caller could
    /// mistake for "nothing nearby".
    pub fn for_each_within(
        &self,
        query: Point3,
        radius: Scalar,
        mut visit: impl FnMut(PointHit),
    ) -> Result<(), PointQueryError> {
        if !query.is_finite() {
            return Err(PointQueryError::NonFiniteQuery);
        }
        if !radius.is_finite() || radius < 0.0 {
            return Err(PointQueryError::InvalidRadius);
        }
        if self.ordered.is_empty() {
            return Ok(());
        }

        let radius_squared = radius * radius;
        let lo = self.cell_of(query - Point3::splat(radius));
        let hi = self.cell_of(query + Point3::splat(radius));

        for iz in lo[2]..=hi[2] {
            for iy in lo[1]..=hi[1] {
                for ix in lo[0]..=hi[0] {
                    let cell_index = (iz * self.dims[1] + iy) * self.dims[0] + ix;
                    let from = self.starts[cell_index] as usize;
                    let to = self.starts[cell_index + 1] as usize;
                    for &point_index in &self.ordered[from..to] {
                        let point = self.points[point_index as usize];
                        // Compare squared distances: the square root is
                        // only paid for points that actually qualify.
                        let squared = (point - query).length_squared();
                        if squared <= radius_squared {
                            visit(PointHit {
                                index: point_index as usize,
                                distance: squared.sqrt(),
                            });
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Collect every point within `radius`, ordered by distance then index.
    ///
    /// The vector is cleared first, so a caller may reuse one buffer across
    /// many queries and allocate once.
    ///
    /// # Errors
    ///
    /// As [`for_each_within`](Self::for_each_within).
    pub fn radius_into(
        &self,
        query: Point3,
        radius: Scalar,
        out: &mut Vec<PointHit>,
    ) -> Result<(), PointQueryError> {
        out.clear();
        self.for_each_within(query, radius, |hit| out.push(hit))?;
        sort_hits(out);
        Ok(())
    }

    /// Collect the `k` nearest points, ordered by distance then index.
    ///
    /// Fewer than `k` are returned when the cloud holds fewer usable
    /// points; that is a complete answer, not a truncated one.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite query position.
    pub fn nearest_into(
        &self,
        query: Point3,
        k: usize,
        out: &mut Vec<PointHit>,
    ) -> Result<(), PointQueryError> {
        out.clear();
        if !query.is_finite() {
            return Err(PointQueryError::NonFiniteQuery);
        }
        if k == 0 || self.ordered.is_empty() {
            return Ok(());
        }

        // Grow a search sphere until it holds k points, then confirm with
        // one exact radius query. Doubling keeps the number of rounds
        // logarithmic in how badly the first guess was scaled, and the
        // final query guarantees completeness -- expanding rings alone can
        // miss a nearer point sitting just outside the last ring visited.
        let mut radius = self.cell * (k as Scalar).cbrt().max(1.0);
        let ceiling = self.diagonal();
        loop {
            self.radius_into(query, radius.min(ceiling), out)?;
            if out.len() >= k || radius >= ceiling {
                break;
            }
            radius *= 2.0;
        }
        out.truncate(k);
        Ok(())
    }

    /// The single nearest point, if the cloud has a usable one.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite query position.
    pub fn nearest(&self, query: Point3) -> Result<Option<PointHit>, PointQueryError> {
        let mut out = Vec::new();
        self.nearest_into(query, 1, &mut out)?;
        Ok(out.into_iter().next())
    }

    /// Clamp a position to a cell coordinate.
    fn cell_of(&self, p: Point3) -> [usize; 3] {
        let axis = |value: Scalar, origin: Scalar, dim: usize| -> usize {
            let raw = (value - origin) / self.cell;
            if raw < 0.0 {
                0
            } else {
                (raw as usize).min(dim - 1)
            }
        };
        [
            axis(p.x, self.origin.x, self.dims[0]),
            axis(p.y, self.origin.y, self.dims[1]),
            axis(p.z, self.origin.z, self.dims[2]),
        ]
    }

    /// Diagonal of the indexed extent, used as a search ceiling.
    fn diagonal(&self) -> Scalar {
        let span = Point3::new(
            self.dims[0] as Scalar * self.cell,
            self.dims[1] as Scalar * self.cell,
            self.dims[2] as Scalar * self.cell,
        );
        (span.x * span.x + span.y * span.y + span.z * span.z).sqrt()
    }
}

/// Order hits by distance, breaking ties by index.
fn sort_hits(hits: &mut [PointHit]) {
    hits.sort_by(|a, b| {
        a.distance
            .partial_cmp(&b.distance)
            .unwrap_or(core::cmp::Ordering::Equal)
            .then(a.index.cmp(&b.index))
    });
}
