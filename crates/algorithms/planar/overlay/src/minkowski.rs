//! Minkowski sums and erosions of a region by a polygon or a region,
//! convex or not (#145), and disc morphology with a known side of error
//! (#163).
//!
//! # One arrangement per operation
//!
//! For a region `R`, a connected shape `K` and any point `k0` of `K`,
//!
//! ```text
//! R + K = (R + k0)  union  (union over boundary edges e of R: e + K)
//! ```
//!
//! since a point `x` has `x - K` meeting `R` either through the boundary
//! (then `x` lies in some `e + K`) or wholly inside (then `x - k0` is in
//! `R`). For a convex `K`, `e + K` is the convex hull of `K + a` and
//! `K + b`, `e` running from `a` to `b`. A non-convex `K` -- with holes, if
//! it is a region's polygon -- is cut into convex pieces `P` with its own
//! vertices (`convex_parts`: ear clipping and Hertel-Mehlhorn
//! merging, the cut certified to tile `K`), and `e + K` is the union of the
//! hulls `e + P`. A shape of several polygons is the union of their sums.
//! Should no certified cut be found, `R + B` is taken from the
//! parallelograms `f + g` of boundary edge pairs, `R` moved by a vertex of
//! `B`, and `B` moved by a vertex of every ring of `R`: if neither
//! boundary meets the other yet the two sets do, an outer ring of one lies
//! inside the other.
//!
//! The erosion `R - K` (the points `x` with `x + K` inside `R`) is `R - k0`
//! minus the sum of the complement of `R` with `-K`, the complement taken
//! within a box large enough that nothing beyond it matters; when `K`
//! holds the origin, `R` itself stands in for `R - k0`.
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
    /// The structuring polygon is not convex. No longer returned: sums
    /// and erosions with a non-convex polygon are built since 0.3.10
    /// (#145); kept so that matches on it still compile.
    NotConvex,
    /// An operand or the result failed the overlay's own checks.
    Overlay(OverlayError),
    /// An erosion by an empty region, under which every point qualifies.
    EmptyStructuring,
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

pub(crate) fn orient(a: Point2, b: Point2, c: Point2) -> Sign {
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

/// One polygon of a structuring shape, as the sum takes it.
#[derive(Debug, Clone)]
pub(crate) enum Part {
    /// Convex pieces, counter-clockwise, whose union is the polygon.
    Pieces(Vec<Vec<Point2>>),
    /// The polygon's rings, outer first: when no certified convex
    /// decomposition was found, the sum is taken edge by edge instead.
    Rings(Vec<Vec<Point2>>),
}

impl Part {
    /// The part reflected through the origin.
    fn flipped(&self) -> Self {
        let flip = |v: &Vec<Point2>| v.iter().map(|q| Point2::new(-q.x, -q.y)).collect();
        match self {
            Self::Pieces(p) => Self::Pieces(p.iter().map(flip).collect()),
            Self::Rings(r) => Self::Rings(r.iter().map(flip).collect()),
        }
    }

    /// A point of the part, by which the other operand is translated.
    fn anchor(&self) -> Point2 {
        match self {
            Self::Pieces(p) => p[0][0],
            Self::Rings(r) => r[0][0],
        }
    }

    fn points(&self) -> impl Iterator<Item = &Point2> {
        match self {
            Self::Pieces(p) | Self::Rings(p) => p.iter().flatten(),
        }
    }

    /// Whether the origin lies in one of the convex pieces (closed). A
    /// part taken by its rings answers `false`, which is always safe.
    fn holds_origin(&self) -> bool {
        let o = Point2::new(0.0, 0.0);
        match self {
            Self::Pieces(pieces) => pieces.iter().any(|piece| {
                let h = hull(piece.clone());
                h.len() >= 3
                    && (0..h.len()).all(|i| orient(h[i], h[(i + 1) % h.len()], o) != Sign::Negative)
            }),
            Self::Rings(_) => false,
        }
    }
}

/// The part for one polygon: its convex hull when it is convex (the
/// convex path, unchanged), else certified convex pieces, else its rings.
fn polygon_part(
    outer: &Ring,
    holes: &[Ring],
    tolerance: Tolerance,
) -> Result<Part, MinkowskiError> {
    if holes.is_empty() {
        match convex(outer, tolerance) {
            Ok(k) => return Ok(Part::Pieces(vec![k])),
            Err(MinkowskiError::NotConvex) => {}
            Err(other) => return Err(other),
        }
    }
    let holes: Vec<Vec<Point2>> = holes.iter().map(|h| h.points.clone()).collect();
    Ok(
        match crate::convex_parts::convex_parts(&outer.points, &holes) {
            Some(pieces) => Part::Pieces(pieces),
            None => Part::Rings(std::iter::once(outer.points.clone()).chain(holes).collect()),
        },
    )
}

/// A set given by rings in the arrangement, read from the flags starting
/// at its first ring.
#[derive(Debug, Clone)]
enum Set {
    /// Polygons with holes laid out as [`rings_of`] lays them out (the
    /// hole count per polygon).
    Region(Vec<usize>),
    /// A frame (the first ring) minus such polygons (the rings after it).
    Outside(Vec<usize>),
}

impl Set {
    fn holds(&self, flags: &[bool]) -> bool {
        match self {
            Self::Region(shape) => member(shape, flags),
            Self::Outside(shape) => flags[0] && !member(shape, &flags[1..]),
        }
    }
}

/// The rings whose union is the sum of a set with a structuring shape,
/// and how to read it from an arrangement's ring flags.
struct Terms {
    rings: Vec<ArcRing>,
    /// Translated copies of a set: the first ring's index, and the set.
    sets: Vec<(usize, Set)>,
    /// Index of the first convex piece; every ring from there on is one.
    pieces: usize,
}

impl Terms {
    fn holds(&self, flags: &[bool]) -> bool {
        self.sets.iter().any(|(at, set)| set.holds(&flags[*at..]))
            || flags[self.pieces..].iter().any(|&f| f)
    }
}

/// The sum of the set `set`, bounded by `rings`, with the union of `parts`.
///
/// For a connected part `K` and a point `k0` of it,
/// `S + K = (S + k0) union (union over boundary edges e of S: e + K)`, and
/// `e + K` is the union of `e + P` over convex pieces `P` of `K`, each the
/// convex hull of `P + a` and `P + b`. A part taken by its rings `B` uses
/// `S + B = (dS + dB) union (S + b0) union (a + B for a vertex a of every
/// ring of S)`, `dS + dB` the parallelograms of edge pairs: if the
/// boundaries of `x - B` and `S` do not meet yet the two do, an outer ring
/// of one lies inside the other.
fn sum_terms(rings: &[Vec<Point2>], set: &Set, parts: &[Part], tolerance: Tolerance) -> Terms {
    let shift = |p: Point2, by: Point2| Point2::new(p.x + by.x, p.y + by.y);
    let moved = |r: &[Point2], by: Point2| {
        ArcRing::from_points(&r.iter().map(|&p| shift(p, by)).collect::<Vec<_>>())
    };
    let mut out = Vec::new();
    let mut sets = Vec::new();
    for part in parts {
        sets.push((out.len(), set.clone()));
        out.extend(rings.iter().map(|r| moved(r, part.anchor())));
        if let Part::Rings(own) = part {
            let shape = vec![own.len() - 1];
            for ring in rings {
                sets.push((out.len(), Set::Region(shape.clone())));
                out.extend(own.iter().map(|r| moved(r, ring[0])));
            }
        }
    }
    let pieces = out.len();
    for ring in rings {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            for part in parts {
                match part {
                    Part::Pieces(convex) => {
                        for k in convex {
                            let points: Vec<Point2> =
                                k.iter().flat_map(|&q| [shift(a, q), shift(b, q)]).collect();
                            out.extend(tidy(hull(points), tolerance));
                        }
                    }
                    Part::Rings(own) => {
                        for r in own {
                            for j in 0..r.len() {
                                let (p, q) = (r[j], r[(j + 1) % r.len()]);
                                let points =
                                    vec![shift(a, p), shift(a, q), shift(b, p), shift(b, q)];
                                out.extend(tidy(hull(points), tolerance));
                            }
                        }
                    }
                }
            }
        }
    }
    Terms {
        rings: out,
        sets,
        pieces,
    }
}

/// One part per polygon of a region.
fn parts_of(region: &Region, tolerance: Tolerance) -> Result<Vec<Part>, MinkowskiError> {
    region
        .polygons()
        .iter()
        .map(|p| polygon_part(&p.outer, &p.holes, tolerance))
        .collect()
}

/// The faces of `arrangement` where `inside` holds, as a region.
pub(crate) fn region_of(
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
    Ok(Region::from_valid_polygons(
        crate::settle::settle(polygons, tolerance),
        had_input,
    ))
}

impl Region {
    /// The Minkowski sum with a polygon `ring`: every point of the region
    /// moved by every point of the polygon. The polygon may be non-convex.
    ///
    /// Decided exactly (ADR 0070), with each vertex sum `a + k` rounded
    /// once to `f64` and output vertices rounded once. A convex polygon is
    /// summed whole; a non-convex one is cut into convex pieces (its own
    /// vertices, the cut certified to tile it) whose sums are united in
    /// the same arrangement.
    ///
    /// # Errors
    ///
    /// [`MinkowskiError::Overlay`] for a polygon the overlay refuses or a
    /// refusal of the arrangement.
    pub fn minkowski_sum(&self, ring: &Ring, tolerance: Tolerance) -> Result<Self, MinkowskiError> {
        let part = polygon_part(ring, &[], tolerance)?;
        Ok(self.sum_parts(&[part], tolerance)?)
    }

    /// The Minkowski erosion by a polygon `ring`, convex or not: the points
    /// `x` for which the polygon moved by `x` lies within the region.
    ///
    /// # Errors
    ///
    /// As [`Self::minkowski_sum`].
    pub fn minkowski_erosion(
        &self,
        ring: &Ring,
        tolerance: Tolerance,
    ) -> Result<Self, MinkowskiError> {
        let part = polygon_part(ring, &[], tolerance)?;
        Ok(self.erode_parts(&[part], tolerance)?)
    }

    /// The Minkowski sum with another region -- both may be non-convex,
    /// have holes and several components: every point of one moved by
    /// every point of the other. Symmetric in the two regions.
    ///
    /// The region with fewer vertices is cut into convex pieces and the
    /// other summed with each, all in one exact arrangement, rounded as
    /// [`Self::minkowski_sum`] rounds.
    ///
    /// # Errors
    ///
    /// [`MinkowskiError::Overlay`] for a refusal of the arrangement.
    pub fn minkowski_sum_region(
        &self,
        other: &Self,
        tolerance: Tolerance,
    ) -> Result<Self, MinkowskiError> {
        if self.is_empty() || other.is_empty() {
            return Ok(Self::empty());
        }
        let size = |r: &Self| {
            r.polygons()
                .iter()
                .flat_map(|p| std::iter::once(&p.outer).chain(&p.holes))
                .map(|r| r.points.len())
                .sum::<usize>()
        };
        let (base, cut) = if size(self) < size(other) {
            (other, self)
        } else {
            (self, other)
        };
        let parts = parts_of(cut, tolerance)?;
        Ok(base.sum_parts(&parts, tolerance)?)
    }

    /// The Minkowski erosion by another region: the points `x` for which
    /// the other region moved by `x` lies within this one.
    ///
    /// # Errors
    ///
    /// [`MinkowskiError::EmptyStructuring`] for an empty `other` (every
    /// point would qualify), and as [`Self::minkowski_sum_region`].
    pub fn minkowski_erosion_region(
        &self,
        other: &Self,
        tolerance: Tolerance,
    ) -> Result<Self, MinkowskiError> {
        if other.is_empty() {
            return Err(MinkowskiError::EmptyStructuring);
        }
        let parts = parts_of(other, tolerance)?;
        Ok(self.erode_parts(&parts, tolerance)?)
    }

    fn sum_with(&self, k: &[Point2], tolerance: Tolerance) -> Result<Self, OverlayError> {
        self.sum_parts(&[Part::Pieces(vec![k.to_vec()])], tolerance)
    }

    fn erode_with(&self, k: &[Point2], tolerance: Tolerance) -> Result<Self, OverlayError> {
        self.erode_parts(&[Part::Pieces(vec![k.to_vec()])], tolerance)
    }

    pub(crate) fn sum_parts(
        &self,
        parts: &[Part],
        tolerance: Tolerance,
    ) -> Result<Self, OverlayError> {
        if self.is_empty() {
            return Ok(Self::empty());
        }
        let own = rings_of(self.polygons());
        let shape: Vec<usize> = self.polygons().iter().map(|p| p.holes.len()).collect();
        let terms = sum_terms(&own, &Set::Region(shape), parts, tolerance);
        let arrangement = ArcArrangement::new(&terms.rings, tolerance)?;
        region_of(&arrangement, |flags| terms.holds(flags), true, tolerance)
    }

    pub(crate) fn erode_parts(
        &self,
        parts: &[Part],
        tolerance: Tolerance,
    ) -> Result<Self, OverlayError> {
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
        let reach = parts
            .iter()
            .flat_map(Part::points)
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
        let flipped: Vec<Part> = parts.iter().map(Part::flipped).collect();
        let terms = sum_terms(&outside, &Set::Outside(shape.clone()), &flipped, tolerance);
        let m = own.len();
        let mut rings: Vec<ArcRing> = own.iter().map(|r| ArcRing::from_points(r)).collect();
        rings.extend(terms.rings.iter().cloned());
        let arrangement = ArcArrangement::new(&rings, tolerance)?;
        // A point of the erosion moved by the anchor `k0` of the first
        // part lies in the region: when the origin is in K that is implied
        // by the region itself; otherwise it is read from the region moved
        // by `-k0`, whose rings follow the frame's in the first translated
        // copy (frame minus region).
        let origin = parts.iter().any(Part::holds_origin);
        let first = m + terms.sets[0].0 + 1;
        region_of(
            &arrangement,
            |flags| {
                let anchored = if origin {
                    member(&shape, &flags[..m])
                } else {
                    member(&shape, &flags[first..])
                };
                anchored && !terms.holds(&flags[m..])
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Polygon;

    fn ring(points: &[(f64, f64)]) -> Ring {
        Ring {
            points: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
        }
    }

    fn u_with_hole() -> Region {
        Region::new(
            vec![Polygon {
                outer: ring(&[
                    (0.0, 0.0),
                    (6.0, 0.0),
                    (6.0, 5.0),
                    (4.0, 5.0),
                    (4.0, 2.0),
                    (2.0, 2.0),
                    (2.0, 5.0),
                    (0.0, 5.0),
                ]),
                holes: vec![ring(&[(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)])],
            }],
            Tolerance::METRE,
        )
        .unwrap()
    }

    fn shapes() -> Vec<(Ring, Vec<Ring>)> {
        vec![
            // An L-shape.
            (
                ring(&[
                    (0., 0.),
                    (0.8, 0.),
                    (0.8, 0.4),
                    (0.4, 0.4),
                    (0.4, 0.8),
                    (0., 0.8),
                ]),
                vec![],
            ),
            // A dart off the origin, clockwise.
            (
                ring(&[(0.6, 0.4), (0.3, 1.0), (1.5, 0.4), (0.3, -0.2)]),
                vec![],
            ),
            // A square ring.
            (
                ring(&[(-0.6, -0.6), (0.6, -0.6), (0.6, 0.6), (-0.6, 0.6)]),
                vec![ring(&[(-0.3, -0.3), (0.3, -0.3), (0.3, 0.3), (-0.3, 0.3)])],
            ),
        ]
    }

    fn apart(a: &Region, b: &Region) -> f64 {
        let t = Tolerance::METRE;
        a.difference(b, t).unwrap().area() + b.difference(a, t).unwrap().area()
    }

    #[test]
    fn non_convex_polygons_are_cut_into_few_convex_pieces() {
        let pieces: Vec<usize> = shapes()
            .iter()
            .map(
                |(outer, holes)| match polygon_part(outer, holes, Tolerance::METRE).unwrap() {
                    Part::Pieces(p) => p.len(),
                    Part::Rings(_) => panic!("not decomposed"),
                },
            )
            .collect();
        assert_eq!(pieces[0], 2);
        assert_eq!(pieces[1], 2);
        assert!(pieces[2] <= 6, "{pieces:?}");
    }

    #[test]
    fn the_edge_pair_route_agrees_with_the_convex_pieces() {
        // The route taken when no certified decomposition is found. The
        // small square fits wholly inside the shapes, where only the shape
        // moved by its vertex finds it.
        let t = Tolerance::METRE;
        let small = Region::new(
            vec![Polygon {
                outer: ring(&[(3.0, 3.0), (3.1, 3.0), (3.1, 3.1), (3.0, 3.1)]),
                holes: vec![],
            }],
            t,
        )
        .unwrap();
        let r = u_with_hole().union(&small, t).unwrap();
        for (outer, holes) in shapes() {
            let cut = polygon_part(&outer, &holes, t).unwrap();
            let whole = Part::Rings(
                std::iter::once(outer.points.clone())
                    .chain(holes.iter().map(|h| h.points.clone()))
                    .collect(),
            );
            let (a, b) = (
                r.sum_parts(std::slice::from_ref(&cut), t).unwrap(),
                r.sum_parts(std::slice::from_ref(&whole), t).unwrap(),
            );
            assert!(!a.is_empty() && apart(&a, &b) < 1e-9, "{}", apart(&a, &b));
            let (a, b) = (
                r.erode_parts(std::slice::from_ref(&cut), t).unwrap(),
                r.erode_parts(std::slice::from_ref(&whole), t).unwrap(),
            );
            assert!(!a.is_empty() && apart(&a, &b) < 1e-9, "{}", apart(&a, &b));
        }
    }
}
