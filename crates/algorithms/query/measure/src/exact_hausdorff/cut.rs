//! Matched bounds for faces cut down to their free regions (#229).
//!
//! A body whose items share patches of face is measured on the boundary of
//! the union: each face in contact replaced by its free region (see
//! `exact_bodies`). The free regions of a translated copy are cut from
//! translated faces, but an exact cut of nearly coincident outlines can
//! fall out differently in the last bit -- a sliver here, an extra vertex
//! there -- and then no free region of one is a translate of one of the
//! other. This bounds them through the faces they were cut from instead.
//!
//! # The bound
//!
//! Let `x` be a point of the free region of face `F` of body `A`: on `F`,
//! in no face `G` of another item lying against `F`. Let `g` be a face of
//! `B`'s items matched to `F` (the same trims, or a translate's, as for
//! any pair), and `t` the displacement of the match at the chart origin.
//! The matched bound gives a point `y` of `g` with `|x - y| <= m`, and
//! the same bound for `S_F + t` against `g` gives `|y - (x + t)| <= e`.
//! Either `y` is in no face lying against `g` -- then it is on `B`'s union
//! boundary -- or it is in one, `G'`. Every such `G'` must be matched to a
//! face `G` lying against `F`, and its boundary loops lie within `zeta` of
//! those of `G + t` (the matched bound of `G + t` against `G'` over a box
//! holding `G`'s trims, with the trim residue of a translate). Then, by
//! the straight homotopy between the two loops (as in `translate`), a
//! point farther than `zeta` from `G`'s loops is in `G` exactly when its
//! translate is in `G'`. `y - t` is in `G' - t` and within `e` of `x`,
//! which is not in `G`, so `y - t` is within `max(zeta, e)` of `G`'s
//! loops and `y` within `max(zeta, e) + zeta` of `G'`'s. Those loops are
//! on `B`'s union boundary when each of their edges is also an edge of a
//! face of its item with nothing lying against it. So
//!
//! `d(x, dB) <= m + max(zeta, e) + zeta`.
//!
//! For a translated copy `m` is `|t|` and `e` and `zeta` are rounding, so
//! every free region closes at once whatever its cut looks like.

use std::collections::HashMap;

use axiolid_brep::ExactBRep;
use axiolid_core::{Point2, Scalar, Vec3};
use axiolid_curve::Curve2;
use axiolid_surface::{Plane, Surface};

use super::{matched_bound, matches, Match};
use crate::exact::ExactMeasureError;
use crate::exact_distance::surface_of;

/// How a body's boundary was cut from its items.
pub(crate) struct Cut<'a> {
    /// All the body's items, uncut, as one B-rep.
    pub(crate) items: &'a ExactBRep,
    /// For each face of the measured boundary cut from a face of `items`,
    /// that face.
    pub(crate) origin: HashMap<usize, usize>,
    /// The faces of `items` lying against each face in contact.
    pub(crate) partners: HashMap<usize, Vec<usize>>,
}

/// The bound of the module docs between two cut bodies, with what it needs
/// prepared once.
pub(crate) struct CutBounds<'a> {
    from: &'a Cut<'a>,
    to: &'a Cut<'a>,
    matched: Vec<Vec<Match>>,
    /// Per pair `(F, index into matched[F])`: the displacement and `zeta`,
    /// or `None` where the bound does not apply.
    pairs: HashMap<(usize, usize), Option<(Vec3, Scalar)>>,
    /// Whether every edge of a face of `to.items` is an edge of another
    /// face of its item with nothing lying against it.
    rim_free: HashMap<usize, bool>,
}

impl<'a> CutBounds<'a> {
    pub(crate) fn new(from: &'a Cut<'a>, to: &'a Cut<'a>) -> Self {
        Self {
            from,
            to,
            matched: matches(from.items, to.items),
            pairs: HashMap::new(),
            rim_free: HashMap::new(),
        }
    }

    /// A bound on `d(x, dB)` over the patch `[lo, hi]` of face `face` of the
    /// measured boundary `model`, or `None`.
    pub(crate) fn bound(
        &mut self,
        model: &ExactBRep,
        face: usize,
        lo: Point2,
        hi: Point2,
    ) -> Result<Option<Scalar>, ExactMeasureError> {
        let Some(&origin) = self.from.origin.get(&face) else {
            return Ok(None);
        };
        let Surface::Plane(cut) = surface_of(model, model.topology().faces()[face].surface)? else {
            return Ok(None);
        };
        let items = self.from.items;
        let Surface::Plane(whole) =
            surface_of(items, items.topology().faces()[origin].surface)?.clone()
        else {
            return Ok(None);
        };
        // The patch in the parameters of the face it was cut from.
        let Some((lo, hi)) = recharted(cut, &whole, lo, hi) else {
            return Ok(None);
        };
        let mut best: Option<Scalar> = None;
        for index in 0..self.matched[origin].len() {
            let Some((t, zeta)) = self.pair(origin, index)? else {
                continue;
            };
            let found = &self.matched[origin][index];
            let target = target_surface(self.to.items, found)?;
            let moved = Surface::Plane(Plane {
                frame: axiolid_core::Frame3 {
                    origin: whole.frame.origin + t,
                    ..whole.frame
                },
            });
            let surface = Surface::Plane(whole);
            let (Some(m), Some(e)) = (
                finish(found, matched_bound(&surface, &target, lo, hi), lo, hi),
                finish(found, matched_bound(&moved, &target, lo, hi), lo, hi),
            ) else {
                continue;
            };
            let bound = m + e.max(zeta) + zeta;
            if bound.is_finite() {
                best = Some(best.map_or(bound, |b: Scalar| b.min(bound)));
            }
        }
        Ok(best)
    }

    /// The displacement and `zeta` of the `index`th match of face `face`.
    fn pair(
        &mut self,
        face: usize,
        index: usize,
    ) -> Result<Option<(Vec3, Scalar)>, ExactMeasureError> {
        if let Some(known) = self.pairs.get(&(face, index)) {
            return Ok(*known);
        }
        let found = self.pair_uncached(face, index)?;
        self.pairs.insert((face, index), found);
        Ok(found)
    }

    fn pair_uncached(
        &mut self,
        face: usize,
        index: usize,
    ) -> Result<Option<(Vec3, Scalar)>, ExactMeasureError> {
        let (from, to) = (self.from.items, self.to.items);
        let found = self.matched[face][index].clone();
        let Surface::Plane(whole) = surface_of(from, from.topology().faces()[face].surface)? else {
            return Ok(None);
        };
        let Surface::Plane(target) = target_surface(to, &found)? else {
            return Ok(None);
        };
        let t = target.frame.origin - whole.frame.origin;
        if !t.is_finite() {
            return Ok(None);
        }
        let mut zeta: Scalar = 0.0;
        let theirs = self
            .to
            .partners
            .get(&found.face)
            .cloned()
            .unwrap_or_default();
        let mine = self.from.partners.get(&face).cloned().unwrap_or_default();
        for against in theirs {
            if !self.rim_free(against) {
                return Ok(None);
            }
            // A face lying against `F` that `against` is a translate of.
            let mut best: Option<Scalar> = None;
            for &candidate in &mine {
                for other in &self.matched[candidate] {
                    if other.face != against {
                        continue;
                    }
                    if let Some(residue) = residue(from, candidate, to, other, t)? {
                        best = Some(best.map_or(residue, |b: Scalar| b.min(residue)));
                    }
                }
            }
            let Some(residue) = best else {
                return Ok(None);
            };
            zeta = zeta.max(residue);
        }
        Ok(Some((t, zeta)))
    }

    /// Whether every edge of face `face` of `to.items` is also an edge of a
    /// face of its item with nothing lying against it.
    fn rim_free(&mut self, face: usize) -> bool {
        if let Some(known) = self.rim_free.get(&face) {
            return *known;
        }
        let topology = self.to.items.topology();
        let mut users: HashMap<usize, Vec<usize>> = HashMap::new();
        for (index, f) in topology.faces().iter().enumerate() {
            for bound in &f.bounds {
                if let Some(wire) = topology.loops().get(bound.loop_id.index()) {
                    for use_ in &wire.edges {
                        users.entry(use_.edge.index()).or_default().push(index);
                    }
                }
            }
        }
        let free = topology.faces().get(face).is_some_and(|f| {
            f.bounds.iter().all(|bound| {
                topology
                    .loops()
                    .get(bound.loop_id.index())
                    .is_some_and(|wire| {
                        wire.edges.iter().all(|use_| {
                            users.get(&use_.edge.index()).is_some_and(|faces| {
                                faces.iter().any(|&other| {
                                    other != face && !self.to.partners.contains_key(&other)
                                })
                            })
                        })
                    })
            })
        });
        self.rim_free.insert(face, free);
        free
    }
}

/// The surface of the matched face, re-charted onto `F`'s parameters for a
/// translate trimmed in another chart.
fn target_surface(to: &ExactBRep, found: &Match) -> Result<Surface, ExactMeasureError> {
    Ok(match &found.shifted {
        Some(shifted) => shifted.surface.clone(),
        None => surface_of(to, to.topology().faces()[found.face].surface)?.clone(),
    })
}

/// A matched bound with a translate's trim residue added.
fn finish(found: &Match, matched: Option<Scalar>, lo: Point2, hi: Point2) -> Option<Scalar> {
    let matched = matched?;
    match &found.shifted {
        Some(shifted) => shifted.bound(matched, lo, hi),
        None => Some(matched),
    }
}

/// How far the loops of `B`'s face `other` are from those of `A`'s face
/// `face` moved by `t`: the matched bound of the moved face against it over
/// a box holding the face's trims, with a translate's trim residue.
fn residue(
    from: &ExactBRep,
    face: usize,
    to: &ExactBRep,
    other: &Match,
    t: Vec3,
) -> Result<Option<Scalar>, ExactMeasureError> {
    let Surface::Plane(plane) = surface_of(from, from.topology().faces()[face].surface)? else {
        return Ok(None);
    };
    let Some((lo, hi)) = trim_box(from, face) else {
        return Ok(None);
    };
    let moved = Surface::Plane(Plane {
        frame: axiolid_core::Frame3 {
            origin: plane.frame.origin + t,
            ..plane.frame
        },
    });
    let target = target_surface(to, other)?;
    Ok(finish(
        other,
        matched_bound(&moved, &target, lo, hi),
        lo,
        hi,
    ))
}

/// A parameter box holding every pcurve of a face's loops, for lines and
/// circles.
fn trim_box(brep: &ExactBRep, face: usize) -> Option<(Point2, Point2)> {
    let topology = brep.topology();
    let mut lo = Point2::new(Scalar::INFINITY, Scalar::INFINITY);
    let mut hi = Point2::new(Scalar::NEG_INFINITY, Scalar::NEG_INFINITY);
    let mut grow = |p: Point2, r: Scalar| {
        lo = Point2::new(lo.x.min(p.x - r), lo.y.min(p.y - r));
        hi = Point2::new(hi.x.max(p.x + r), hi.y.max(p.y + r));
    };
    for bound in &topology.faces().get(face)?.bounds {
        let wire = topology.loops().get(bound.loop_id.index())?;
        for (index, use_) in wire.edges.iter().enumerate() {
            match brep.curves2().get(use_.pcurve?.index())? {
                Curve2::Line(line) => {
                    let span = brep.pcurve_interval(bound.loop_id, index)?;
                    for s in [span.start, span.end] {
                        grow(
                            Point2::new(
                                line.origin.x + line.direction.x * s,
                                line.origin.y + line.direction.y * s,
                            ),
                            0.0,
                        );
                    }
                }
                Curve2::Circle(circle) => {
                    let scale = circle.frame.x.length().max(circle.frame.y.length());
                    grow(circle.frame.origin, circle.radius.abs() * scale);
                }
                _ => return None,
            }
        }
    }
    let pad = 1e-12 * (1.0 + lo.x.abs().max(lo.y.abs()).max(hi.x.abs()).max(hi.y.abs()));
    let (lo, hi) = (
        Point2::new(lo.x - pad, lo.y - pad),
        Point2::new(hi.x + pad, hi.y + pad),
    );
    (lo.x.is_finite() && hi.x.is_finite() && lo.y.is_finite() && hi.y.is_finite())
        .then_some((lo, hi))
}

/// A box in the parameters of plane `whole` holding the patch `[lo, hi]` of
/// plane `cut` (the same plane, charted otherwise).
fn recharted(cut: &Plane, whole: &Plane, lo: Point2, hi: Point2) -> Option<(Point2, Point2)> {
    let f = whole.frame;
    let (xx, yy, xy) = (f.x.dot(f.x), f.y.dot(f.y), f.x.dot(f.y));
    let det = xx * yy - xy * xy;
    if !det.is_finite() || det <= 0.0 {
        return None;
    }
    let chart = |p: axiolid_core::Point3| {
        let d = p - f.origin;
        let (a, b) = (d.dot(f.x), d.dot(f.y));
        Point2::new((a * yy - b * xy) / det, (b * xx - a * xy) / det)
    };
    let corners = [
        Point2::new(lo.x, lo.y),
        Point2::new(hi.x, lo.y),
        Point2::new(lo.x, hi.y),
        Point2::new(hi.x, hi.y),
    ]
    .map(|q| chart(cut.frame.origin + cut.frame.x * q.x + cut.frame.y * q.y));
    let mut a = corners[0];
    let mut b = corners[0];
    for c in &corners[1..] {
        a = Point2::new(a.x.min(c.x), a.y.min(c.y));
        b = Point2::new(b.x.max(c.x), b.y.max(c.y));
    }
    // Rounding in the chart change: it holds the patch's off-plane
    // residue only through `pad`, which the bounds' own pads exceed.
    let pad = 1e-12 * (1.0 + a.x.abs().max(a.y.abs()).max(b.x.abs()).max(b.y.abs()));
    let (a, b) = (
        Point2::new(a.x - pad, a.y - pad),
        Point2::new(b.x + pad, b.y + pad),
    );
    (a.x.is_finite() && b.x.is_finite() && a.y.is_finite() && b.y.is_finite()).then_some((a, b))
}
