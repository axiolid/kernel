//! The bound of a cell beside a mitre (#335, ADR 0082 amendment).
//!
//! # The exact solid at a mitre
//!
//! Where a run crosses a seam whose tangents differ, the parent module's
//! [seams](super::super#seams-263) mitre it: the surface `P_in` of the piece
//! ending at the seam and `P_out` of the piece starting there (the surfaces of
//! the [parent module](super)) meet in the plane `H` through the seam point
//! moved along its unit normal `n` by the longitudinal offset there. Write
//! `h(x) = (x - q0) . n` for the height above `H`. The exact solid is
//!
//! - `P_in` continued past the seam by extruding its end section along its end
//!   tangent `t_in`, cut by `H` and kept where `h <= 0`; and `P_out`
//!   continued backwards along `t_out`, kept where `h >= 0` (the projection
//!   of the amendment on seams, along each side's own tangent);
//! - the planar faces of `H` inside one side's cut and not the other's,
//!   where the two cuts differ: on the inner side of a curved piece, in a
//!   turn on a grade in the section frame, and with a longitudinal offset.
//!
//! The mesh keeps one mitre row `M`, the midpoint of the two sides' projections
//! moved by the offset (`axiolid_reference::station::Mitre::place`), so no
//! topology changes; the bound measures `M` against the solid above.
//!
//! # The cell
//!
//! Take a cell `[lo, hi]` with a mitre at `hi`; at `lo` everything is
//! mirrored (the outgoing side, `toward = -t_out`, `normal = -n`). Let `a`
//! be the mesh's far row, `c` the exact section of the cell's own side at
//! the knot, `P = c - toward h(c) / cos` its projection onto `H` along
//! `toward = t_in` (`cos = toward . n`), `M` the mesh's mitre row and
//! `w = c - a` the chord of each vertex. Each term below holds along a wall
//! edge `(k, l)` because everything in it is affine in `v`.
//!
//! 1. The per-vertex bound `eta` of the parent (chord, breakpoints at seams
//!    inside the cell, read errors) holds against the quad `V = (a, c)`,
//!    measured to the exact section `c` and not to the mesh's mitre row: the
//!    exact surface up to the knot is `V` up to `eta`. The extrusion is
//!    exactly the strip `(c, P)`.
//! 2. A point of `V` with `h <= 0` and the point of `(a, P)` at the same `v`
//!    and the same height differ only in `H`'s direction, by at most the
//!    **kink** `K = |h(c)| |w_H| / (w . n)`, `w_H = w - (w . n / cos) toward`
//!    ([`kink`], over eight stretches of `v`). It is second order on a
//!    curved piece and first order in a changing section's rate. Pairing the
//!    two vertices of an edge at equal distance instead would charge an
//!    in-plane slide of about the section's width times the half turn, since
//!    they cross `H` at different distances.
//! 3. A point of `V` up to `eta'` past `H` (the read errors) is clamped back
//!    onto it for at most `eta' Gamma`, `Gamma = max(|w_k| / w_k . n, |w_l| /
//!    w_l . n, 1 / cos)`.
//! 4. The quads `(a, P)` and `(a, M)` differ by `D = |M - P|` per vertex; the
//!    triangles lie within `quad_deviation` of their quad.
//! 5. A face of `H` is a point inside one cut ring and not the other. By the
//!    winding number of the two rings it lies within the crossing error of
//!    the mesh's mitre row, whose edges are mesh edges, so within `D`, the
//!    row's distance from `H` and the reach of the extrusion past the plane
//!    ([`reach_past`]).
//!
//! The read errors are the plane's: `e_p` of the point and `e_t`, `e_t'` of
//! the tangents on the cell's side and across the knot (the frame laws'
//! `point_error`, `axis_error` and measure), `e_n = (e_t + e_t') / cos` of
//! the normal, so `eta_H = e_p + rho e_n` at the reach `rho` of the points
//! involved.
//!
//! The wall and the cut are
//!
//! `wall = eta_s + K + D + quad_deviation + eta' Gamma`,
//! `cut = eta_s + K + D + max(eta' Gamma, reach) + hm + eta_H`,
//!
//! `eta_s = eta + lambda e_t` (the vertex's bound plus the read error of the
//! extrusion's length, `lambda` the largest extension), `eta' = eta_s +
//! eta_H`, `hm` the mesh row's largest distance from `H`; the cell's bound is
//! the largest over its edges, times `1 + 1e-12`. The cell's **floor** is
//! the largest `D`: splitting a cell never removes it, so refinement stops at
//! a cell whose floor exceeds the budget.
//!
//! A cell on a straight piece with a constant section (every vertex runs on a
//! line through `H`: `K = 0`, `P` is exact and the quad is the wall) keeps
//! only `D`, `hm`, `eta_H` and the quad: [`SectionLaw::straight_mitre_cell`].
//! It is exact for `D = 0`, which holds when the two sides' cuts coincide (a
//! flat turn without an offset), and reports the real gap otherwise.
//!
//! # Refusals
//!
//! The argument needs each vertex to cross `H` once over the cell
//! (`mu_j = (w_j . n - 2 slide_j - |w_j| e_n) / span - span second_j > 0`,
//! `span = hi - lo + 2 measure`: the advance across `H` at the chord's rate
//! less the curvature's, so the crossing is unique), the plane to advance
//! (`cos - e_t - e_n > 0`) and the far row to lie strictly on the cell's side
//! (`h(a_j) + slide_j + eta_H < 0`). A cell where it does not hold is
//! unbounded by name, `MITRE_NOT_TRANSVERSAL`, and a curved or changing cell
//! between two mitres, `MITRE_ON_BOTH_ENDS`; both are refinable
//! ([`refinable`](super::refinable)): a shorter cell crosses `H` more
//! transversally, and a half of a cell between two mitres has one.

use axiolid_contracts::GeomResult;
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_model::SeamSide;
use axiolid_reference::station::{FrameLawBounds, FrameReading, Mitre};

use super::super::seam_slack;
use super::{CellBound, SectionLaw, MITRE_NOT_TRANSVERSAL, UNCERTIFIED_FRAMES};

/// What `SectionLaw::cell` knows of a cell beside a mitre.
pub(super) struct MitreInputs<'a> {
    pub(super) lo: Scalar,
    pub(super) hi: Scalar,
    pub(super) at_lo: &'a [Point3],
    pub(super) at_hi: &'a [Point3],
    /// Per vertex, the chord, breakpoint and read terms.
    pub(super) eta: &'a [Scalar],
    pub(super) slide: &'a [Scalar],
    /// Per vertex, the largest second derivative over the cell.
    pub(super) second: &'a [Scalar],
    /// The largest distance a frame is read away from the one asked for.
    pub(super) measure: Scalar,
    /// The frame bounds of the pieces next to the cell's ends.
    pub(super) first: &'a FrameLawBounds,
    pub(super) last: &'a FrameLawBounds,
    pub(super) reading: FrameReading,
    /// The longitudinal offset at the cell's ends.
    pub(super) long: [Scalar; 2],
}

/// One side of a mitre, as the cell beside it reads it.
struct SideCut {
    /// The point of the mitre plane the exact solid is cut at: the seam
    /// point moved by the longitudinal offset along the normal.
    q0: Point3,
    /// The plane's normal, pointing out of the cell.
    normal: Vec3,
    /// The cell's end tangent, pointing at the plane.
    toward: Vec3,
    cos: Scalar,
    /// The cell's own exact section at the knot.
    exact: Vec<Point3>,
    /// That section extruded along `toward` onto the plane.
    projected: Vec<Point3>,
    /// Per vertex, the mesh's mitre row's distance from `projected`.
    gap: Vec<Scalar>,
    /// The mesh's mitre row's largest distance from the plane.
    off_plane: Scalar,
    /// How far the plane's point may lie from the exact one.
    point_error: Scalar,
    /// How far the end tangent may lie from the exact one.
    tangent_error: Scalar,
    /// How far the plane's normal may lie from the exact one.
    normal_error: Scalar,
}

/// Whether `x` is positive: a NaN is not, so it refuses.
fn positive(x: Scalar) -> bool {
    x.partial_cmp(&0.0) == Some(core::cmp::Ordering::Greater)
}

/// How many equal stretches of an edge [`kink`] splits it into.
const STRETCHES: usize = 8;

/// A bound on `sup |h(v)| |oblique(w(v))| / (w(v) . n)` over an edge: `h`
/// the height above the plane and `w` the vertex chord's, both affine in
/// `v`, `oblique(x) = x - toward (x . n) / cos` the part of `x` the
/// extrusion along `toward` does not remove. Over a stretch of `v` the
/// factor `|h|` is at most its larger end value, and the ratio of an affine
/// vector's length to a positive affine scalar is quasi-convex, so at most
/// its larger end value: the product of the two bounds the stretch.
fn kink(h: [Scalar; 2], w: [Vec3; 2], normal: Vec3, toward: Vec3, cos: Scalar) -> Scalar {
    let ratio = |v: Scalar| {
        let x = w[0] + (w[1] - w[0]) * v;
        (x - toward * (x.dot(normal) / cos)).length() / x.dot(normal)
    };
    let height = |v: Scalar| (h[0] + (h[1] - h[0]) * v).abs();
    let mut worst: Scalar = 0.0;
    for i in 0..STRETCHES {
        let (v0, v1) = (
            i as Scalar / STRETCHES as Scalar,
            (i + 1) as Scalar / STRETCHES as Scalar,
        );
        worst = worst.max(height(v0).max(height(v1)) * ratio(v0).max(ratio(v1)));
    }
    worst
}

/// A bound on `sup |e(v)| min(1, eta / en(v))` over an edge, `e` an affine
/// vector and `en > 0` an affine scalar: the part of the extrusion `e`
/// that a point within `eta` of the plane can still reach. `|e|` is
/// convex and `|e| / en` quasi-convex, so the supremum sits at an end of
/// the edge or where `en` crosses `eta`.
fn reach_past(e: [Vec3; 2], en: [Scalar; 2], eta: Scalar) -> Scalar {
    let end = |j: usize| e[j].length() * (eta / en[j]).min(1.0);
    let mut worst = end(0).max(end(1));
    if (en[0] - eta) * (en[1] - eta) < 0.0 {
        let v = (eta - en[0]) / (en[1] - en[0]);
        worst = worst.max((e[0] + (e[1] - e[0]) * v).length());
    }
    worst
}

impl SideCut {
    fn height(&self, x: Point3) -> Scalar {
        (x - self.q0).dot(self.normal)
    }

    /// The bound on how far the plane read from the frames lies from the
    /// exact one at the points involved: `far` is the cell's other row.
    fn plane_error(&self, far: &[Point3], eta_max: Scalar) -> Scalar {
        let reach = far
            .iter()
            .chain(&self.exact)
            .chain(&self.projected)
            .map(|x| (*x - self.q0).length())
            .fold(0.0, Scalar::max)
            + eta_max;
        self.point_error + reach * self.normal_error
    }

    fn lift(&self) -> Scalar {
        self.exact
            .iter()
            .map(|c| self.height(*c).abs())
            .fold(0.0, Scalar::max)
    }
}

impl SectionLaw<'_, '_> {
    /// The frame bounds over a span of `length` on the other side of the
    /// knot `at`, clipped to the next seam, station or run end.
    fn across(
        &self,
        at: Scalar,
        side: SeamSide,
        length: Scalar,
        reading: FrameReading,
    ) -> Option<FrameLawBounds> {
        let mut marks: Vec<Scalar> = self
            .seams
            .iter()
            .chain(self.distances.iter())
            .copied()
            .collect();
        marks.retain(|m| (m - at).abs() > seam_slack(at));
        let (first, last) = (self.distances[0], self.distances[self.distances.len() - 1]);
        match side {
            // The cell is incoming: the other side starts at the knot.
            SeamSide::Incoming => {
                let end = marks
                    .iter()
                    .copied()
                    .filter(|m| *m > at)
                    .fold((at + length).min(last), Scalar::min);
                self.frames.bounds(at, end, reading)
            }
            _ => {
                let start = marks
                    .iter()
                    .copied()
                    .filter(|m| *m < at)
                    .fold((at - length).max(first), Scalar::max);
                self.frames.bounds(start, at, reading)
            }
        }
    }

    /// The cell's side of `mitre`, `side` the cell's own: `Err` names why
    /// the frames beside it are not bounded.
    fn side_cut(
        &self,
        inputs: &MitreInputs<'_>,
        mitre: &Mitre,
        side: SeamSide,
    ) -> GeomResult<Result<SideCut, &'static str>> {
        let incoming = side == SeamSide::Incoming;
        let (at, row, near, long) = if incoming {
            (inputs.hi, inputs.at_hi, inputs.last, inputs.long[1])
        } else {
            (inputs.lo, inputs.at_lo, inputs.first, inputs.long[0])
        };
        let (toward, normal, base) = if incoming {
            (mitre.incoming.tangent, mitre.normal, mitre.incoming.point)
        } else {
            (-mitre.outgoing.tangent, -mitre.normal, mitre.outgoing.point)
        };
        let Some(other) = self.across(at, side, inputs.hi - inputs.lo, inputs.reading) else {
            return Ok(Err(UNCERTIFIED_FRAMES));
        };
        let cos = toward.dot(normal);
        let tangent_error = near.axis_error + near.measure * near.tangent[0];
        let other_error = other.axis_error + other.measure * other.tangent[0];
        let mut cut = SideCut {
            // The shift is along the mitre's own normal, whichever side.
            q0: base + mitre.normal * long,
            normal,
            toward,
            cos,
            exact: self.exact_at(at, side)?,
            projected: Vec::new(),
            gap: Vec::new(),
            off_plane: 0.0,
            point_error: near.point_error + near.measure * near.point[0],
            tangent_error,
            normal_error: (tangent_error + other_error) / cos,
        };
        cut.projected = cut
            .exact
            .iter()
            .map(|&c| c - cut.toward * (cut.height(c) / cut.cos))
            .collect();
        cut.gap = row
            .iter()
            .zip(&cut.projected)
            .map(|(m, p)| (*m - *p).length())
            .collect();
        cut.off_plane = row
            .iter()
            .map(|m| cut.height(*m).abs())
            .fold(0.0, Scalar::max);
        Ok(Ok(cut))
    }

    /// A cell on a straight piece with a constant section, beside one or
    /// both of its mitres: every vertex runs on a line through the plane,
    /// so the exact wall is the quad up to the plane and the mesh's mitre
    /// row differs from it by the rows' gap alone.
    pub(super) fn straight_mitre_cell(
        &self,
        inputs: &MitreInputs<'_>,
        lo: Option<&Mitre>,
        hi: Option<&Mitre>,
    ) -> GeomResult<Result<CellBound, &'static str>> {
        let count = inputs.at_lo.len();
        let mut gap = vec![0.0 as Scalar; count];
        let (mut off_plane, mut plane_error, mut extension): (Scalar, Scalar, Scalar) =
            (0.0, 0.0, 0.0);
        let eta_max = inputs.eta.iter().copied().fold(0.0, Scalar::max);
        for (mitre, side, far) in [
            (lo, SeamSide::Outgoing, inputs.at_hi),
            (hi, SeamSide::Incoming, inputs.at_lo),
        ] {
            let Some(mitre) = mitre else { continue };
            let cut = match self.side_cut(inputs, mitre, side)? {
                Ok(cut) => cut,
                Err(reason) => return Ok(Err(reason)),
            };
            for (g, d) in gap.iter_mut().zip(&cut.gap) {
                *g = g.max(*d);
            }
            off_plane = off_plane.max(cut.off_plane);
            plane_error = plane_error.max(cut.plane_error(far, eta_max));
            extension = extension.max(cut.lift() / cut.cos * cut.tangent_error);
        }
        let mut worst: Scalar = 0.0;
        for &(k, l) in self.edges {
            let quad = axiolid_construct::loft::quad_deviation(
                inputs.at_lo[k],
                inputs.at_lo[l],
                inputs.at_hi[k],
                inputs.at_hi[l],
            );
            let gap = gap[k].max(gap[l]);
            let slide = inputs.slide[k].max(inputs.slide[l]) + extension;
            worst = worst.max(slide + (quad + gap).max(gap + off_plane + plane_error));
        }
        let floor = gap.iter().copied().fold(0.0, Scalar::max);
        let worst = worst * (1.0 + 1e-12);
        Ok(if worst.is_finite() {
            Ok(CellBound {
                value: worst,
                floor,
            })
        } else {
            Err(UNCERTIFIED_FRAMES)
        })
    }

    /// A cell beside one mitre on a curved piece, a changing section or a
    /// seam inside it: the bound through the plane's cut of the surface.
    pub(super) fn mitre_cell(
        &self,
        inputs: &MitreInputs<'_>,
        mitre: &Mitre,
        side: SeamSide,
    ) -> GeomResult<Result<CellBound, &'static str>> {
        let cut = match self.side_cut(inputs, mitre, side)? {
            Ok(cut) => cut,
            Err(reason) => return Ok(Err(reason)),
        };
        let far = if side == SeamSide::Incoming {
            inputs.at_lo
        } else {
            inputs.at_hi
        };
        let count = far.len();
        let eta_max = inputs.eta.iter().copied().fold(0.0, Scalar::max);
        let plane_error = cut.plane_error(far, eta_max);
        let share = cut.cos - cut.tangent_error - cut.normal_error;
        if !positive(share) {
            return Ok(Err(MITRE_NOT_TRANSVERSAL));
        }
        let reach = (cut.lift() + eta_max + plane_error) / share;
        // Each vertex must cross the plane once over the cell, and the far
        // row lie strictly on the cell's side.
        let span = inputs.hi - inputs.lo + 2.0 * inputs.measure;
        let chord: Vec<Vec3> = (0..count).map(|j| cut.exact[j] - far[j]).collect();
        let advance: Vec<Scalar> = chord.iter().map(|w| w.dot(cut.normal)).collect();
        for j in 0..count {
            let rate = (advance[j] - 2.0 * inputs.slide[j] - chord[j].length() * cut.normal_error)
                / span
                - span * inputs.second[j];
            let clear = cut.height(far[j]) + inputs.slide[j] + plane_error < 0.0;
            if !positive(rate) || !clear {
                return Ok(Err(MITRE_NOT_TRANSVERSAL));
            }
        }
        let extruded: Vec<Vec3> = (0..count).map(|j| cut.projected[j] - far[j]).collect();
        let depth: Vec<Scalar> = (0..count).map(|j| -cut.height(far[j])).collect();
        let mut worst: Scalar = 0.0;
        for &(k, l) in self.edges {
            let along = inputs.eta[k].max(inputs.eta[l]) + reach * cut.tangent_error;
            let near = along + plane_error;
            let slope = (chord[k].length() / advance[k])
                .max(chord[l].length() / advance[l])
                .max(1.0 / cut.cos);
            let kink = kink(
                [cut.height(cut.exact[k]), cut.height(cut.exact[l])],
                [chord[k], chord[l]],
                cut.normal,
                cut.toward,
                cut.cos,
            );
            let gap = cut.gap[k].max(cut.gap[l]);
            let quad = axiolid_construct::loft::quad_deviation(
                inputs.at_lo[k],
                inputs.at_lo[l],
                inputs.at_hi[k],
                inputs.at_hi[l],
            );
            let clamp = near * slope;
            let past = reach_past([extruded[k], extruded[l]], [depth[k], depth[l]], near);
            let wall = along + kink + gap + quad + clamp;
            let face = along + kink + gap + clamp.max(past) + cut.off_plane + plane_error;
            worst = worst.max(wall.max(face));
        }
        let floor = cut.gap.iter().copied().fold(0.0, Scalar::max);
        let worst = worst * (1.0 + 1e-12);
        Ok(if worst.is_finite() {
            Ok(CellBound {
                value: worst,
                floor,
            })
        } else {
            Err(UNCERTIFIED_FRAMES)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DENSE: usize = 10_000;

    #[test]
    fn the_kink_bounds_its_dense_supremum_closely() {
        let normal = Vec3::new(0.6, 0.8, 0.0);
        let toward = Vec3::new(0.8, 0.6, 0.0);
        let cos = toward.dot(normal);
        // The height changes sign along the first edge; the chords advance
        // across the plane at every point of both.
        let cases = [
            (
                [0.4, -0.7],
                [Vec3::new(0.5, 0.3, 0.2), Vec3::new(0.1, 0.6, -0.4)],
            ),
            (
                [-0.9, -0.1],
                [Vec3::new(1.0, 0.5, -0.3), Vec3::new(0.2, 0.2, 0.9)],
            ),
        ];
        for (h, w) in cases {
            let dense = (0..=DENSE)
                .map(|i| {
                    let v = i as Scalar / DENSE as Scalar;
                    let x = w[0] + (w[1] - w[0]) * v;
                    assert!(x.dot(normal) > 0.0);
                    let height = (h[0] + (h[1] - h[0]) * v).abs();
                    height * (x - toward * (x.dot(normal) / cos)).length() / x.dot(normal)
                })
                .fold(0.0, Scalar::max);
            let bound = kink(h, w, normal, toward, cos);
            assert!(bound >= dense, "{bound} < {dense}");
            // Eight stretches: a product of end maxima, loose by the stretch.
            assert!(bound <= 2.0 * dense, "{bound} is loose against {dense}");
        }
    }

    #[test]
    fn the_reach_past_the_plane_bounds_its_dense_supremum() {
        let cases = [
            // The band within `eta` of the plane is crossed inside the
            // edge, and the crossing dominates both ends.
            (
                [Vec3::new(-0.6, 0.3, 0.2), Vec3::new(-0.1, 0.0, 0.06)],
                [0.7, 0.0025],
                0.125,
                true,
            ),
            // The band lies past both ends, or holds them both.
            (
                [Vec3::new(0.3, 0.1, 0.0), Vec3::new(-0.2, 0.4, 0.1)],
                [0.5, 0.4],
                0.2,
                false,
            ),
            (
                [Vec3::new(0.3, 0.1, 0.0), Vec3::new(0.2, 0.4, 0.1)],
                [0.1, 0.05],
                0.2,
                false,
            ),
        ];
        for (e, en, eta, crossing) in cases {
            let at = |v: Scalar| (e[0] + (e[1] - e[0]) * v, en[0] + (en[1] - en[0]) * v);
            let dense = (0..=DENSE)
                .map(|i| {
                    let (x, depth) = at(i as Scalar / DENSE as Scalar);
                    x.length() * (eta / depth).min(1.0)
                })
                .fold(0.0, Scalar::max);
            let bound = reach_past(e, en, eta);
            assert!(bound >= dense, "{bound} < {dense}");
            assert!(bound <= dense + 1e-3, "{bound} is loose against {dense}");
            if crossing {
                // The crossing is what the ends alone miss.
                let ends = (0..2)
                    .map(|j| e[j].length() * (eta / en[j]).min(1.0))
                    .fold(0.0, Scalar::max);
                assert!(ends + 0.05 < dense, "{ends} vs {dense}");
            }
        }
    }
}
