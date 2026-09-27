//! Minkowski sums and erosions of a region by a convex polygon (#145), and
//! disc morphology with a known side of error (#163).
//!
//! # One arrangement per operation
//!
//! For a region `R` and a convex polygon `K`, and any vertex `k0` of `K`,
//!
//! ```text
//! R + K = (R + k0)  union  (union over boundary edges e of R: e + K)
//! ```
//!
//! since a point `x` has `x - K` meeting `R` either through the boundary
//! (then `x` lies in some `e + K`) or wholly inside (then `x - k0` is in
//! `R`). Each `e + K` is the convex hull of `K + a` and `K + b`, `e` running
//! from `a` to `b`. The erosion `R - K` (the points `x` with `x + K` inside
//! `R`) is `R` minus the sum of its complement with `-K`, the complement
//! taken within a box large enough that nothing beyond it matters.
//!
//! All the rings -- the region's own, the translated ones and the pieces --
//! go into one [`ArcArrangement`] (ADR 0070): where boundaries cross, which
//! pieces coincide and which ring holds what are decided exactly, and the
//! result's faces are chosen by ring membership. Vertices are rounded once.
//! The one rounding before that is each vertex sum `a + k`, a sum of two
//! `f64`s, within half an ulp.
//!
//! # Discs: inner and outer
//!
//! No polygon is a disc, so [`Region::dilate`] and [`Region::erode`] are
//! approximations of unstated side. The inscribed polygon `P_in` of a disc
//! `D` lies inside it and the circumscribed `P_out` around it, so
//!
//! ```text
//! R + P_in  <=  R + D  <=  R + P_out        R - P_out  <=  R - D  <=  R - P_in
//! ```
//!
//! Each polygon is moved a margin further to its side, so that the
//! roundings above -- and the dropping of edges shorter than the tolerance
//! when the result is presented -- cannot carry a point across. What each
//! side proves, and how far it can be from the true disc morphology, is in
//! [`MorphologyBound`].

use axiolid_core::{Point2, Tolerance};
use axiolid_exact::{certify, Arith, SignExpr};
use axiolid_guarantees::Sign;

use crate::arc::ArcRing;
use crate::arrangement::ArcArrangement;
use crate::region::Region;
use crate::{validate_ring, OverlayError, Polygon, Ring};

/// Why a Minkowski operation was refused.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum MinkowskiError {
    /// The structuring polygon is not convex. Sums with a non-convex
    /// polygon are not built yet (#145).
    NotConvex,
    /// An operand or the result failed the overlay's own checks.
    Overlay(OverlayError),
}

impl From<OverlayError> for MinkowskiError {
    fn from(error: OverlayError) -> Self {
        Self::Overlay(error)
    }
}

/// Which side of the exact disc morphology a region lies on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum BoundSide {
    /// Contained in the exact result. A route found in an inner erosion
    /// proves reachability; "unreachable" proves nothing.
    Inner,
    /// Containing the exact result. "Unreachable" in an outer erosion is a
    /// proof; a route found in it may be invalid.
    Outer,
}

/// How a disc morphology's result relates to the exact one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MorphologyBound {
    /// Which side of the exact result it lies on.
    pub side: BoundSide,
    /// The greatest distance between the two boundaries: the polygonal
    /// disc's deviation from the true disc, plus the margin kept for
    /// rounding and presentation.
    pub deviation: f64,
}

/// Orientation of `c` against the directed line `a -> b`, exactly.
struct Orient {
    a: Point2,
    b: Point2,
    c: Point2,
}

impl SignExpr for Orient {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let f = T::from_f64;
        let (ux, uy) = (f(self.b.x).sub(&f(self.a.x)), f(self.b.y).sub(&f(self.a.y)));
        let (vx, vy) = (f(self.c.x).sub(&f(self.a.x)), f(self.c.y).sub(&f(self.a.y)));
        ux.mul(&vy).sub(&uy.mul(&vx)).sign()
    }
}

fn orient(a: Point2, b: Point2, c: Point2) -> Sign {
    certify(&Orient { a, b, c }).unwrap_or(Sign::Zero)
}

/// A convex polygon's vertices, counter-clockwise, collinear ones dropped.
fn convex(ring: &Ring, tolerance: Tolerance) -> Result<Vec<Point2>, MinkowskiError> {
    validate_ring(ring, tolerance)?;
    let hull = hull(ring.points.clone());
    // Convex exactly when every vertex is on its hull.
    let on_hull = ring.points.iter().all(|p| {
        hull.contains(p) || {
            // Or on a hull edge (a collinear vertex).
            (0..hull.len()).any(|i| {
                let (a, b) = (hull[i], hull[(i + 1) % hull.len()]);
                orient(a, b, *p) == Sign::Zero
                    && p.x >= a.x.min(b.x)
                    && p.x <= a.x.max(b.x)
                    && p.y >= a.y.min(b.y)
                    && p.y <= a.y.max(b.y)
            })
        }
    });
    if !on_hull || hull.len() < 3 {
        return Err(MinkowskiError::NotConvex);
    }
    Ok(hull)
}

/// The convex hull, counter-clockwise, without collinear vertices
/// (Andrew's monotone chain on exact orientations).
fn hull(mut points: Vec<Point2>) -> Vec<Point2> {
    points.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));
    points.dedup();
    if points.len() < 3 {
        return points;
    }
    let mut lower: Vec<Point2> = Vec::new();
    for &p in &points {
        while lower.len() >= 2
            && orient(lower[lower.len() - 2], lower[lower.len() - 1], p) != Sign::Positive
        {
            lower.pop();
        }
        lower.push(p);
    }
    let mut upper: Vec<Point2> = Vec::new();
    for &p in points.iter().rev() {
        while upper.len() >= 2
            && orient(upper[upper.len() - 2], upper[upper.len() - 1], p) != Sign::Positive
        {
            upper.pop();
        }
        upper.push(p);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

/// The region's rings in order (per polygon: outer, then holes), and
/// whether a point with the given ring flags lies in the region.
fn rings_of(polygons: &[Polygon]) -> Vec<Vec<Point2>> {
    polygons
        .iter()
        .flat_map(|p| std::iter::once(&p.outer).chain(&p.holes))
        .map(|r| r.points.clone())
        .collect()
}

/// Membership in a polygon set, from per-ring flags laid out as
/// [`rings_of`] lays them out.
fn member(polygons: &[usize], flags: &[bool]) -> bool {
    let mut at = 0;
    for &holes in polygons {
        let inside = flags[at] && !flags[at + 1..at + 1 + holes].iter().any(|&f| f);
        if inside {
            return true;
        }
        at += 1 + holes;
    }
    false
}

/// A ring's points, with neighbours closer than the tolerance merged (a
/// hull of two nearly equal translates has such pairs).
fn tidy(points: Vec<Point2>, tolerance: Tolerance) -> Option<ArcRing> {
    let mut out: Vec<Point2> = Vec::with_capacity(points.len());
    for p in points {
        if out
            .last()
            .is_none_or(|q: &Point2| (p - *q).length() > tolerance.linear())
        {
            out.push(p);
        }
    }
    while out.len() > 1 && (out[0] - out[out.len() - 1]).length() <= tolerance.linear() {
        out.pop();
    }
    (out.len() >= 3).then(|| ArcRing::from_points(&out))
}

/// The rings whose union with the translated set is the sum of the set
/// bounded by `rings` with `k`: the rings translated by `k[0]`, and each
/// boundary edge's piece.
fn sum_rings(
    rings: &[Vec<Point2>],
    k: &[Point2],
    tolerance: Tolerance,
) -> (Vec<ArcRing>, Vec<ArcRing>) {
    let shift = |p: Point2, by: Point2| Point2::new(p.x + by.x, p.y + by.y);
    let translated = rings
        .iter()
        .map(|r| ArcRing::from_points(&r.iter().map(|&p| shift(p, k[0])).collect::<Vec<_>>()))
        .collect();
    let mut pieces = Vec::new();
    for ring in rings {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            let points: Vec<Point2> = k.iter().flat_map(|&q| [shift(a, q), shift(b, q)]).collect();
            if let Some(piece) = tidy(hull(points), tolerance) {
                pieces.push(piece);
            }
        }
    }
    (translated, pieces)
}

/// The faces of `arrangement` where `inside` holds, as a region.
fn region_of(
    arrangement: &ArcArrangement,
    inside: impl Fn(&[bool]) -> bool,
    had_input: bool,
    tolerance: Tolerance,
) -> Result<Region, OverlayError> {
    let straight = |ring: ArcRing| Ring {
        points: ring.vertices.iter().map(|v| v.point).collect(),
    };
    let mut polygons = Vec::new();
    for region in arrangement.regions(inside)? {
        let Some(outer) = crate::arc_overlay::presented(arrangement.ring(&region.outer), tolerance)
        else {
            continue;
        };
        let holes = region
            .holes
            .iter()
            .filter_map(|h| crate::arc_overlay::presented(arrangement.ring(h), tolerance))
            .map(straight)
            .collect();
        polygons.push(Polygon {
            outer: straight(outer),
            holes,
        });
    }
    Ok(Region::from_valid_polygons(polygons, had_input))
}

impl Region {
    /// The Minkowski sum with a convex polygon `convex`: every point of the
    /// region moved by every point of the polygon.
    ///
    /// Decided exactly (ADR 0070), with each vertex sum `a + k` rounded
    /// once to `f64` and output vertices rounded once.
    ///
    /// # Errors
    ///
    /// [`MinkowskiError::NotConvex`] for a polygon that is not convex; the
    /// sum with a non-convex one is not built yet (#145).
    pub fn minkowski_sum(
        &self,
        convex_ring: &Ring,
        tolerance: Tolerance,
    ) -> Result<Self, MinkowskiError> {
        let k = convex(convex_ring, tolerance)?;
        Ok(self.sum_with(&k, tolerance)?)
    }

    /// The Minkowski erosion by a convex polygon: the points `x` for which
    /// the polygon moved by `x` lies within the region.
    ///
    /// # Errors
    ///
    /// As [`Self::minkowski_sum`].
    pub fn minkowski_erosion(
        &self,
        convex_ring: &Ring,
        tolerance: Tolerance,
    ) -> Result<Self, MinkowskiError> {
        let k = convex(convex_ring, tolerance)?;
        Ok(self.erode_with(&k, tolerance)?)
    }

    fn sum_with(&self, k: &[Point2], tolerance: Tolerance) -> Result<Self, OverlayError> {
        if self.is_empty() {
            return Ok(Self::empty());
        }
        let own = rings_of(self.polygons());
        let shape: Vec<usize> = self.polygons().iter().map(|p| p.holes.len()).collect();
        let (translated, pieces) = sum_rings(&own, k, tolerance);
        let n = translated.len();
        let mut rings = translated;
        rings.extend(pieces);
        let arrangement = ArcArrangement::new(&rings, tolerance)?;
        region_of(
            &arrangement,
            |flags| member(&shape, &flags[..n]) || flags[n..].iter().any(|&f| f),
            true,
            tolerance,
        )
    }

    fn erode_with(&self, k: &[Point2], tolerance: Tolerance) -> Result<Self, OverlayError> {
        if self.is_empty() {
            return Ok(Self::empty());
        }
        let own = rings_of(self.polygons());
        let shape: Vec<usize> = self.polygons().iter().map(|p| p.holes.len()).collect();
        // A box holding the region: its complement within the box, summed
        // with -K, reaches every point of the region that K cannot sit
        // around (a translate of K leaving the box crosses its boundary,
        // whose own edge pieces cover that). The room about the region only
        // keeps the box's edges off the region's.
        let (mut lo, mut hi) = (own[0][0], own[0][0]);
        for p in own.iter().flatten() {
            lo = lo.min(*p);
            hi = hi.max(*p);
        }
        let reach = k
            .iter()
            .fold(0.0_f64, |m, q| m.max(q.x.abs()).max(q.y.abs()));
        let pad = 2.0 * reach + (hi - lo).max_element() + 1.0;
        let frame = vec![
            Point2::new(lo.x - pad, lo.y - pad),
            Point2::new(hi.x + pad, lo.y - pad),
            Point2::new(hi.x + pad, hi.y + pad),
            Point2::new(lo.x - pad, hi.y + pad),
        ];
        let mut outside = vec![frame];
        outside.extend(own.iter().cloned());
        let flipped: Vec<Point2> = k.iter().map(|q| Point2::new(-q.x, -q.y)).collect();
        let (translated, pieces) = sum_rings(&outside, &flipped, tolerance);
        let m = own.len();
        let t = translated.len();
        let mut rings: Vec<ArcRing> = own.iter().map(|r| ArcRing::from_points(r)).collect();
        rings.extend(translated);
        rings.extend(pieces);
        let arrangement = ArcArrangement::new(&rings, tolerance)?;
        let in_outside = |flags: &[bool]| flags[0] && !member(&shape, &flags[1..]);
        region_of(
            &arrangement,
            |flags| {
                member(&shape, &flags[..m])
                    && !(in_outside(&flags[m..m + t]) || flags[m + t..].iter().any(|&f| f))
            },
            true,
            tolerance,
        )
    }

    /// Dilation by a disc of `radius`, contained in the exact one: the sum
    /// with an inscribed polygon, drawn in by a margin.
    ///
    /// # Errors
    ///
    /// [`OverlayError::InvalidOffsetDistance`] for a radius that is not
    /// finite and nonnegative, and overlay refusals.
    pub fn dilate_inner(&self, radius: f64, tolerance: Tolerance) -> Result<Self, OverlayError> {
        self.disc(radius, BoundSide::Inner, false, tolerance)
    }

    /// Dilation by a disc of `radius`, containing the exact one: the sum
    /// with a circumscribed polygon, pushed out by a margin.
    ///
    /// # Errors
    ///
    /// As [`Self::dilate_inner`].
    pub fn dilate_outer(&self, radius: f64, tolerance: Tolerance) -> Result<Self, OverlayError> {
        self.disc(radius, BoundSide::Outer, false, tolerance)
    }

    /// Erosion by a disc of `radius`, contained in the exact one: eroded by
    /// a circumscribed polygon, pushed out by a margin. A route found in it
    /// proves reachability.
    ///
    /// # Errors
    ///
    /// As [`Self::dilate_inner`].
    pub fn erode_inner(&self, radius: f64, tolerance: Tolerance) -> Result<Self, OverlayError> {
        self.disc(radius, BoundSide::Inner, true, tolerance)
    }

    /// Erosion by a disc of `radius`, containing the exact one: eroded by
    /// an inscribed polygon, drawn in by a margin. No route in it proves
    /// unreachability.
    ///
    /// # Errors
    ///
    /// As [`Self::dilate_inner`].
    pub fn erode_outer(&self, radius: f64, tolerance: Tolerance) -> Result<Self, OverlayError> {
        self.disc(radius, BoundSide::Outer, true, tolerance)
    }

    fn disc(
        &self,
        radius: f64,
        side: BoundSide,
        erode: bool,
        tolerance: Tolerance,
    ) -> Result<Self, OverlayError> {
        if !radius.is_finite() || radius < 0.0 {
            return Err(OverlayError::InvalidOffsetDistance);
        }
        if self.is_empty() || radius == 0.0 {
            let mut same = Self::from_valid_polygons(self.polygons().to_vec(), false);
            same.set_bound(MorphologyBound {
                side,
                deviation: 0.0,
            });
            return Ok(same);
        }
        // Sixty-four sides: the inscribed polygon's sagitta is
        // r (1 - cos(pi / 64)), about 0.12% of the radius.
        let n = 64;
        let half = core::f64::consts::PI / n as f64;
        // The margin: past the vertex sums' and output vertices' rounding
        // (a few ulps of the largest coordinate), the sin/cos of the
        // polygon's own vertices (relative), and the presentation's dropping
        // of edges shorter than the tolerance.
        let size = self
            .polygons()
            .iter()
            .flat_map(|p| std::iter::once(&p.outer).chain(&p.holes))
            .flat_map(|r| r.points.iter())
            .fold(0.0_f64, |m, p| m.max(p.x.abs()).max(p.y.abs()));
        let margin =
            2.0 * tolerance.linear() + 1e-12 * radius + 64.0 * f64::EPSILON * (size + radius);
        // Inner dilation and outer erosion use the inscribed polygon, the
        // other two the circumscribed one.
        let inscribed = (side == BoundSide::Inner) != erode;
        let (reach, deviation) = if inscribed {
            let reach = radius - margin;
            (reach, radius - reach * half.cos())
        } else {
            let reach = (radius + margin) / half.cos();
            (reach, reach - radius)
        };
        if reach <= 0.0 {
            // The margin swallows the disc: the polygon is a point for the
            // inner dilation (the region itself) and the region itself is an
            // outer erosion.
            let mut same = Self::from_valid_polygons(self.polygons().to_vec(), true);
            same.set_bound(MorphologyBound {
                side,
                deviation: radius + margin,
            });
            return Ok(same);
        }
        let k: Vec<Point2> = (0..n)
            .map(|i| {
                let t = 2.0 * half * i as f64;
                Point2::new(reach * t.cos(), reach * t.sin())
            })
            .collect();
        let mut out = if erode {
            self.erode_with(&k, tolerance)?
        } else {
            self.sum_with(&k, tolerance)?
        };
        out.set_bound(MorphologyBound { side, deviation });
        Ok(out)
    }
}
