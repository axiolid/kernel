// SPDX-License-Identifier: MPL-2.0

//! Bounded Ruppert/Chew quality refinement.
//!
//! # What refinement can and cannot promise
//!
//! Ruppert's algorithm inserts circumcentres of skinny triangles until every
//! angle meets a bound. It provably terminates for a minimum angle up to
//! about 20.7 degrees when the input has no small input angles. Neither
//! condition holds for arbitrary building geometry: two walls meeting at 10
//! degrees is a small *input* angle, and no amount of interior insertion can
//! open it, because the offending angle is pinned by two constraint edges the
//! algorithm is forbidden to move.
//!
//! So this implementation is explicitly *bounded*. It carries a Steiner point
//! budget and reports which of the two outcomes occurred:
//!
//! - [`RefineOutcome::Achieved`] -- every non-constrained angle meets the
//!   bound,
//! - [`RefineOutcome::Capped`] -- the budget ran out first; the mesh is still
//!   a valid constrained Delaunay triangulation, just not as good as asked.
//!
//! Returning `Capped` rather than looping is the difference between a slow
//! call and a hung one, and the difference between a caller who knows the
//! mesh is rough and one who assumes a guarantee that was never met.

use axiolid_core::Point2;

use crate::build::triangulate;
use crate::mesh::{Triangulation, TriangulationError};
use crate::{collinear, Constraint};

/// How much the global worst angle may drop for a single insertion to still
/// be accepted.
///
/// Zero would be ideal but rejects almost everything: retriangulating after
/// an insertion perturbs unrelated triangles by rounding, so a strictly
/// non-decreasing test stalls on noise. This tolerance is small enough that
/// the mesh cannot drift meaningfully worse over a bounded number of
/// insertions, and the end-to-end monotonicity is asserted by test.
const DEGRADATION_TOLERANCE_DEGREES: f64 = 1e-6;

/// A quality target for refinement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quality {
    /// Smallest acceptable interior angle, in degrees.
    ///
    /// Values above ~20.7 are not guaranteed to terminate for every input;
    /// the budget is what keeps the call bounded.
    pub min_angle_degrees: f64,
    /// Largest number of Steiner points to insert.
    ///
    /// Expressed as an absolute count rather than a multiple of the input so
    /// a caller can bound memory directly.
    pub max_steiner_points: usize,
}

impl Default for Quality {
    /// 20 degrees and a budget proportional to a typical wall face.
    ///
    /// 20 sits just under Ruppert's proven bound, so the common case
    /// terminates by the theorem rather than by the budget.
    fn default() -> Self {
        Self {
            min_angle_degrees: 20.0,
            max_steiner_points: 4096,
        }
    }
}

/// What refinement actually achieved.
///
/// Not `Eq`: the capped variant carries a measured angle, and comparing
/// floating-point measurements for exact equality is not a meaningful
/// operation on a quality report.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum RefineOutcome {
    /// Every unconstrained angle meets the requested bound.
    Achieved {
        /// Steiner points inserted.
        inserted: usize,
    },
    /// The budget was exhausted before the bound was met.
    ///
    /// The triangulation is valid and constrained-Delaunay; only the angle
    /// target is unmet.
    Capped {
        /// Steiner points inserted.
        inserted: usize,
        /// Worst angle still present, in degrees.
        worst_angle_degrees: f64,
    },
}

impl RefineOutcome {
    /// Whether the requested quality bound was met.
    #[must_use]
    pub const fn achieved(self) -> bool {
        matches!(self, Self::Achieved { .. })
    }
}

/// Build a constrained Delaunay triangulation and refine it toward `quality`.
///
/// # Errors
///
/// Propagates [`TriangulationError`] from the underlying triangulation.
pub fn triangulate_refined(
    points: &[Point2],
    constraints: &[Constraint],
    quality: Quality,
) -> Result<(Triangulation, RefineOutcome), TriangulationError> {
    let mut tri = triangulate(points, constraints)?;
    let outcome = refine(&mut tri, quality)?;
    Ok((tri, outcome))
}

/// Refine an existing triangulation in place.
///
/// # Errors
///
/// Propagates [`TriangulationError`] if a re-triangulation step fails.
pub fn refine(
    tri: &mut Triangulation,
    quality: Quality,
) -> Result<RefineOutcome, TriangulationError> {
    let threshold = quality.min_angle_degrees.to_radians().cos();
    let mut inserted = 0usize;
    // Triangles whose circumcentre cannot be inserted (it falls outside the
    // hull, or the triangle is degenerate). Skipping them individually is the
    // difference between "this one sliver is unfixable" and "stop refining":
    // abandoning the whole loop on the first such triangle left every other
    // fixable sliver in the mesh untouched.
    let mut skipped: Vec<usize> = Vec::new();

    while inserted < quality.max_steiner_points {
        let Some(bad) = worst_triangle(tri, threshold, &skipped) else {
            break;
        };
        let Some(centre) = circumcentre(tri, bad) else {
            skipped.push(bad);
            continue;
        };
        // Only insert strictly inside the existing hull. A circumcentre
        // outside it would need boundary splitting, which moves a constraint
        // edge -- forbidden here, and the reason this stays "interior-only".
        if !inside_hull(tri, centre) {
            skipped.push(bad);
            continue;
        }
        let mut points = tri.points.clone();
        points.push(centre);
        let constraints = tri.constraints.clone();
        let candidate = triangulate(&points, &constraints)?;
        // Refinement must never hand back a worse mesh than it was given.
        // Inserting the circumcentre of a sliver whose small angle is pinned
        // by input vertices can create an even thinner triangle beside it --
        // measured on a real fixture, 0.51 degrees became 0.22.
        //
        // The test is "does not degrade", not "improves the global worst
        // angle". Demanding global improvement rejects every individual
        // insertion, because fixing one bad triangle usually leaves the
        // single worst one elsewhere untouched -- which silently turned
        // refinement into a no-op.
        let before_worst = worst_angle_degrees(tri);
        let after_worst = worst_angle_degrees(&candidate);
        if after_worst < before_worst - DEGRADATION_TOLERANCE_DEGREES {
            skipped.push(bad);
            continue;
        }
        *tri = candidate;
        inserted += 1;
        // Triangle indices are meaningless after a rebuild.
        skipped.clear();
    }

    let worst = worst_angle_degrees(tri);
    if worst >= quality.min_angle_degrees {
        Ok(RefineOutcome::Achieved { inserted })
    } else {
        Ok(RefineOutcome::Capped {
            inserted,
            worst_angle_degrees: worst,
        })
    }
}

/// The first triangle whose smallest angle is below the bound, ignoring any
/// already known to be un-insertable.
fn worst_triangle(tri: &Triangulation, cos_threshold: f64, skipped: &[usize]) -> Option<usize> {
    (0..tri.triangle_count()).find(|t| {
        if skipped.contains(t) {
            return false;
        }
        let (a, b, c) = corners(tri, *t);
        max_cos(a, b, c) > cos_threshold
    })
}

/// Largest cosine among a triangle's three angles.
///
/// The largest cosine corresponds to the smallest angle, so one comparison
/// against `cos(min_angle)` decides the whole triangle.
fn max_cos(a: Point2, b: Point2, c: Point2) -> f64 {
    let ab = (b.x - a.x, b.y - a.y);
    let bc = (c.x - b.x, c.y - b.y);
    let ca = (a.x - c.x, a.y - c.y);
    let at = angle_cos((-ca.0, -ca.1), ab);
    let bt = angle_cos((-ab.0, -ab.1), bc);
    let ct = angle_cos((-bc.0, -bc.1), ca);
    at.max(bt).max(ct)
}

fn angle_cos(u: (f64, f64), v: (f64, f64)) -> f64 {
    let dot = u.0 * v.0 + u.1 * v.1;
    let nu = (u.0 * u.0 + u.1 * u.1).sqrt();
    let nv = (v.0 * v.0 + v.1 * v.1).sqrt();
    if nu == 0.0 || nv == 0.0 {
        return 1.0;
    }
    (dot / (nu * nv)).clamp(-1.0, 1.0)
}

/// Smallest interior angle anywhere in the mesh, in degrees.
fn worst_angle_degrees(tri: &Triangulation) -> f64 {
    let mut worst: f64 = 180.0;
    for t in 0..tri.triangle_count() {
        let (a, b, c) = corners(tri, t);
        let angle = max_cos(a, b, c).clamp(-1.0, 1.0).acos().to_degrees();
        worst = worst.min(angle);
    }
    worst
}

fn corners(tri: &Triangulation, t: usize) -> (Point2, Point2, Point2) {
    (
        tri.points[tri.triangles[3 * t] as usize],
        tri.points[tri.triangles[3 * t + 1] as usize],
        tri.points[tri.triangles[3 * t + 2] as usize],
    )
}

/// Circumcentre of triangle `t`, or `None` when it is degenerate.
fn circumcentre(tri: &Triangulation, t: usize) -> Option<Point2> {
    let (a, b, c) = corners(tri, t);
    if collinear(a, b, c) {
        return None;
    }
    let d = 2.0 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    if d == 0.0 || !d.is_finite() {
        return None;
    }
    let a2 = a.x * a.x + a.y * a.y;
    let b2 = b.x * b.x + b.y * b.y;
    let c2 = c.x * c.x + c.y * c.y;
    let ux = (a2 * (b.y - c.y) + b2 * (c.y - a.y) + c2 * (a.y - b.y)) / d;
    let uy = (a2 * (c.x - b.x) + b2 * (a.x - c.x) + c2 * (b.x - a.x)) / d;
    if ux.is_finite() && uy.is_finite() {
        Some(Point2::new(ux, uy))
    } else {
        None
    }
}

/// Whether `p` lies inside some existing triangle.
///
/// Used as the interior test for candidate Steiner points: a point inside a
/// triangle is inside the hull by construction, and the mesh is small enough
/// at refinement scale that the linear scan is not the bottleneck.
fn inside_hull(tri: &Triangulation, p: Point2) -> bool {
    (0..tri.triangle_count()).any(|t| {
        let (a, b, c) = corners(tri, t);
        let ab = crate::turns_left(a, b, p);
        let bc = crate::turns_left(b, c, p);
        let ca = crate::turns_left(c, a, p);
        ab == bc && bc == ca
    })
}
