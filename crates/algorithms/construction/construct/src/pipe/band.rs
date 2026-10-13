//! The band of a bend's tube next to a corner mitred beside an arc (#248).
//!
//! A corner between two straight segments shares one mitre ring that lies
//! on both exact cylinders (#245). Beside an arc it cannot: the mitre plane
//! cuts the bend's torus in a different curve than the tube on the other
//! side, so the exact solid, each tube cut by the plane, has a planar face
//! between the two section curves where one tube ends and the other does
//! not begin. The pipe therefore meshes the end of the bend as a band of
//! stations whose vertices walk along the torus itself:
//!
//! - The band starts at a square section `back` short of the corner, and
//!   each profile vertex of it turns about the bend's axis, along its own
//!   circle of the torus, toward the mitre plane. The torus continued past
//!   the arc's end is the exact tube there (the format's mitre at half
//!   angle extends the tube to the plane), so every band vertex lies on the
//!   exact surface.
//! - `back` is the least turn at which the section short of the corner
//!   lies wholly behind the plane, and `ahead` the least past it at which
//!   the section lies wholly beyond, each by at least `MARGIN r`: each
//!   vertex's circle then crosses the plane exactly once between them (its
//!   height above the plane is `A0 + A1 cos psi + A2 sin psi`, with
//!   opposite signs at the band's two ends less than a half turn apart),
//!   and [`reach`] finds both from the closed form of the height's extremes
//!   over the section circle. A corner whose plane does not clear the tube
//!   within a quarter turn either way is refused.
//! - Where the other side's ring is not this band's own cut, the band stops
//!   a little short of the plane (`offset`) and its last span runs, nearly
//!   in the plane, to the shared ring: that strip is the mesh of the planar
//!   face between the two section curves.
//!
//! No proof by construction covers the band and the planar face; the
//! compiler certifies them against the mesh by branch and bound
//! (`axiolid_mesh_compile`'s deviation report), which is why the band's
//! geometry is public.

use core::f64::consts::FRAC_PI_2;

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Vec3};

use crate::loft::Station;
use crate::revolve::rotate;

/// The tube of a bend next to a corner mitred beside it (#248): the torus
/// about the circle through `corner` round `centre` and `axis`, from the
/// square section `back` short of the corner to the mitre plane.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct MitreBand {
    /// Centre of the bend's circle.
    pub centre: Point3,
    /// Unit axis of the bend, oriented so that a positive turn about it
    /// carries the arc toward the corner and on past it.
    pub axis: Vec3,
    /// The arc's end at the corner. The mitre plane passes through it.
    pub corner: Point3,
    /// Unit normal of the mitre plane, pointing away from the band: the
    /// band's part of the exact solid has non-positive height along it.
    pub normal: Vec3,
    /// How far the band reaches back from the corner, in radians: the
    /// section turned by `-back` about `axis` from the corner's lies wholly
    /// behind the mitre plane.
    pub back: Scalar,
    /// How far past the corner the tube is continued, in radians: the
    /// section turned by `ahead` lies wholly beyond the plane, and the
    /// plane cuts every circle of the torus strictly between the two.
    pub ahead: Scalar,
}

/// How far, as a fraction of the disk radius, the band's two end sections
/// must clear the mitre plane.
pub(crate) const MARGIN: Scalar = 0.01;

/// Angles tried for [`reach`], up to a quarter turn.
const REACH_STEPS: usize = 1024;

impl MitreBand {
    /// The radius of the bend's circle.
    #[must_use]
    pub fn bend_radius(&self) -> Scalar {
        (self.corner - self.centre).length()
    }

    /// The signed height of `p` above the mitre plane, along `normal`.
    #[must_use]
    pub fn height(&self, p: Point3) -> Scalar {
        (p - self.corner).dot(self.normal)
    }

    /// `p` turned by `angle` about the bend's axis.
    #[must_use]
    pub fn turned(&self, p: Point3, angle: Scalar) -> Point3 {
        rotate(p, self.centre, self.axis, angle)
    }

    /// The turn in `[lo, hi]` that carries `p` along its circle to height
    /// `level`, found by bisection; `None` unless the height is below
    /// `level` at `lo` and above it at `hi`.
    pub(crate) fn turn_to(
        &self,
        p: Point3,
        level: Scalar,
        lo: Scalar,
        hi: Scalar,
    ) -> Option<Scalar> {
        let f = |t: Scalar| self.height(self.turned(p, t)) - level;
        let (mut a, mut b) = (lo, hi);
        if !(f(a) < 0.0 && f(b) > 0.0) {
            return None;
        }
        for _ in 0..200 {
            let m = 0.5 * (a + b);
            if m <= a || m >= b {
                break;
            }
            if f(m) < 0.0 {
                a = m;
            } else {
                b = m;
            }
        }
        Some(0.5 * (a + b))
    }

    /// The band's stations from `from` (excluded), the square section at
    /// its far end, toward the corner, ending with `target`, the ring the
    /// pipe shares across the mitre.
    ///
    /// Each vertex of `from` turns along its circle by a fraction of the
    /// turn that takes it to the plane (`offset` = `None`) or to `offset`
    /// short of it (at most half its own height below the plane); the turns
    /// are split into equal steps whose sagitta at `reach_radius` is within
    /// `tol`. With `offset` the last band station stands that close to the
    /// plane and `target` follows it; without, `target` takes the place of
    /// the last one, which is the cut itself.
    pub(crate) fn stations_toward(
        &self,
        from: &Station,
        target: &Station,
        offset: Option<Scalar>,
        reach_radius: Scalar,
        tol: Scalar,
    ) -> GeomResult<Vec<Station>> {
        let mut turns: Vec<Vec<Scalar>> = Vec::with_capacity(from.loops.len());
        let mut widest: Scalar = 0.0;
        for ring in &from.loops {
            let mut row = Vec::with_capacity(ring.len());
            for &p in ring {
                let below = -self.height(p);
                let level = offset.map_or(0.0, |cap| -cap.min(0.5 * below));
                let t = self
                    .turn_to(p, level, 0.0, self.back + self.ahead)
                    .ok_or_else(|| {
                        GeomError::Degenerate(
                            "swept disk mitre band: a section vertex does not cross the mitre \
                             plane within the band"
                                .to_owned(),
                        )
                    })?;
                widest = widest.max(t);
                row.push(t);
            }
            turns.push(row);
        }
        let steps = crate::revolve::steps(reach_radius, widest, tol).map_err(|_| {
            GeomError::BudgetExceeded {
                resource: "swept disk mitre band steps",
            }
        })?;
        let last = if offset.is_some() { steps } else { steps - 1 };
        let mut out = Vec::with_capacity(last + 1);
        for j in 1..=last {
            let s = j as Scalar / steps as Scalar;
            out.push(Station {
                loops: from
                    .loops
                    .iter()
                    .zip(&turns)
                    .map(|(ring, row)| {
                        ring.iter()
                            .zip(row)
                            .map(|(&p, &t)| self.turned(p, s * t))
                            .collect()
                    })
                    .collect(),
            });
        }
        out.push(target.clone());
        Ok(out)
    }

    /// The section of the band's tube by the mitre plane, at the vertices
    /// of `square`, a section standing square to the arc at the corner:
    /// each vertex turned along its circle onto the plane.
    pub(crate) fn cut(&self, square: &Station) -> GeomResult<Station> {
        let loops = square
            .loops
            .iter()
            .map(|ring| {
                ring.iter()
                    .map(|&p| {
                        self.turn_to(p, 0.0, -self.back, self.ahead)
                            .map(|t| self.turned(p, t))
                            .ok_or_else(|| {
                                GeomError::Degenerate(
                                    "swept disk mitre band: a section vertex does not cross \
                                     the mitre plane within the band"
                                        .to_owned(),
                                )
                            })
                    })
                    .collect::<GeomResult<Vec<_>>>()
            })
            .collect::<GeomResult<Vec<_>>>()?;
        Ok(Station { loops })
    }
}

/// The least turns `(back, ahead)` (multiples of a quarter turn over
/// [`REACH_STEPS`]) at which the section `back` short of the corner lies
/// wholly at least `MARGIN radius` behind the plane through `corner` with
/// unit normal `normal`, and the section `ahead` past it wholly that far
/// beyond; `None` when either needs more than a quarter turn.
///
/// With `e` the unit vector from `centre` to `corner` (bend radius `R`),
/// `t = axis x e` and a point `corner + r (cos b e + sin b axis)` of the
/// corner's section turned by `psi`, the height above the plane is
/// `R ((cos psi - 1) e.n + sin psi t.n) + r (cos b (cos psi e.n + sin psi
/// t.n) + sin b axis.n)`; its extremes over `b` are the first term plus or
/// minus `r` times the length of the second's coefficients. Any smaller
/// radius (a bore, the centre) lies between them.
pub(crate) fn reach(
    centre: Point3,
    axis: Vec3,
    corner: Point3,
    normal: Vec3,
    radius: Scalar,
) -> Option<(Scalar, Scalar)> {
    let big = (corner - centre).length();
    let e = (corner - centre) / big;
    let t = axis.cross(e);
    let (en, tn, an) = (e.dot(normal), t.dot(normal), axis.dot(normal));
    let margin = MARGIN * radius;
    let turns = || (1..=REACH_STEPS).map(|k| FRAC_PI_2 * k as Scalar / REACH_STEPS as Scalar);
    let back = turns().find(|&psi| {
        let (s, c) = psi.sin_cos();
        big * ((c - 1.0) * en - s * tn) + radius * ((c * en - s * tn).powi(2) + an * an).sqrt()
            <= -margin
    })?;
    let ahead = turns().find(|&psi| {
        let (s, c) = psi.sin_cos();
        big * ((c - 1.0) * en + s * tn) - radius * ((c * en + s * tn).powi(2) + an * an).sqrt()
            >= margin
    })?;
    Some((back, ahead))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn band(turn: Scalar, tilt: Vec3, big: Scalar, r: Scalar) -> MitreBand {
        // A bend in the xy plane about +z ending at (big, 0, 0) heading +y;
        // the mitre normal bisects +y and the outgoing direction.
        let centre = Point3::ZERO;
        let axis = Vec3::Z;
        let corner = Point3::new(big, 0.0, 0.0);
        let u = Vec3::Y;
        let v = (u * turn.cos() + tilt.normalize() * turn.sin()).normalize();
        let normal = (u + v).normalize();
        let (back, ahead) = reach(centre, axis, corner, normal, r).expect("a band");
        MitreBand {
            centre,
            axis,
            corner,
            normal,
            back,
            ahead,
        }
    }

    #[test]
    fn every_circle_of_the_band_crosses_the_plane_once_inside_it() {
        let r = 0.05;
        for (turn, tilt) in [
            (0.3, Vec3::X),
            (1.2, -Vec3::X),
            (FRAC_PI_2, Vec3::X),
            (FRAC_PI_2, -Vec3::X),
            (1.0, Vec3::Z),
            (0.8, Vec3::new(1.0, 0.0, 1.0)),
        ] {
            let b = band(turn, tilt, 0.2, r);
            let e = Vec3::X;
            for k in 0..360 {
                let beta = core::f64::consts::TAU * k as Scalar / 360.0;
                for rho in [r, 0.5 * r, 0.0] {
                    let p = b.corner + (e * beta.cos() + Vec3::Z * beta.sin()) * rho;
                    let lo = b.height(b.turned(p, -b.back));
                    let hi = b.height(b.turned(p, b.ahead));
                    assert!(lo <= -MARGIN * r * (1.0 - 1e-9), "{turn} {tilt:?}: {lo}");
                    assert!(hi >= MARGIN * r * (1.0 - 1e-9), "{turn} {tilt:?}: {hi}");
                    // Dense: exactly one sign change.
                    let n = 400;
                    let changes = (0..n)
                        .filter(|&j| {
                            let at = |j: usize| {
                                b.height(b.turned(
                                    p,
                                    -b.back + (b.back + b.ahead) * j as Scalar / n as Scalar,
                                ))
                            };
                            (at(j) < 0.0) != (at(j + 1) < 0.0)
                        })
                        .count();
                    assert_eq!(changes, 1, "{turn} {tilt:?} {beta} {rho}");
                }
            }
        }
    }

    #[test]
    fn a_plane_that_wraps_round_the_bend_has_no_band() {
        // A bend of radius 1.1 r turning a corner of nearly half a turn.
        let r = 0.05;
        let corner = Point3::new(1.1 * r, 0.0, 0.0);
        let u = Vec3::Y;
        let v = (u * 3.0f64.cos() + Vec3::X * 3.0f64.sin()).normalize();
        let normal = (u + v).normalize();
        assert!(reach(Point3::ZERO, Vec3::Z, corner, normal, r).is_none());
    }
}
