//! A region's rings are apart and nested (#448).
//!
//! [`crate::validate_ring`] judges each ring alone, so `Region::new`
//! accepted a hole crossing or leaving its outer ring, holes overlapping
//! and polygons overlapping. The operations read such a set three ways: the
//! boolean by winding numbers, the morphology by ring membership, and the
//! fast morphology's notch fill (`crate::fill`) assumes rings apart and
//! nested, so it took a hole's needles beyond its outer ring for the set's
//! boundary and left out a point 0.9 m inside the result (fuzz seed 21351).
//! [`check`] refuses such a set by name.
//!
//! # What is refused
//!
//! Each ring is read as the region it encloses, whatever its winding, and a
//! polygon is its outer ring's inside less its holes'. Refused:
//!
//! - [`OverlayError::HoleOutsideOuter`]: part of a hole's inside lies
//!   outside its own outer ring.
//! - [`OverlayError::HolesOverlap`]: part of a hole's inside lies in
//!   another hole of its polygon (crossing, nested or repeated).
//! - [`OverlayError::PolygonsOverlap`]: a point lies in two polygons.
//!
//! Rings may meet: at points, and along a stretch two rings run with their
//! insides on opposite sides. Settled outputs do both: a hole touching its
//! outer ring at a vertex (#191), and a union of cells returning two
//! polygons that share an edge (`tests/roundtrip.rs`), which every
//! operation reads as the point set they bound. A ring whose vertices all
//! lie on one line encloses nothing and is left out (#219); a hole of an
//! outer ring on one line lies outside it. When a set has several defects,
//! one of them is named, the same one for the same input. Polygons whose
//! contacts cut a polygon's inside apart are not refused.
//!
//! # How
//!
//! Every decision is an exact sign and no tolerance enters, as in
//! `check_arc_polygon` (ADR 0095). First the edges of different rings are
//! paired through the crate's box tree, and a proper crossing whose defect
//! the two rings' roles alone name -- a hole with its own outer ring or with
//! another hole of its polygon, or the outer rings of two polygons without
//! holes -- refuses at once, so a soup given as a region is refused before
//! its crossings are cut. Then every ring goes into one exact arrangement
//! (`exact_arc::arrangement`), and each side of each piece is read from the
//! rings that hold it, as the morphology reads a side (`minkowski`'s
//! `Layout`). A crossing at a vertex, a hole leaving through two vertices
//! on its outer ring's edge, and a stretch two rings run the same way all
//! show as a side held wrongly. A region of one ring is not checked.
//!
//! The arrangement of a valid region has about one piece per edge, so the
//! check costs about what a boolean on the region costs to cut it.
//!
//! `scripts/probe_region_nesting_mutants.py` lists the faults the tests in
//! `tests/region_nesting.rs` and here must catch.

use axiolid_core::Point2;
use axiolid_guarantees::Sign;

use crate::arc::ArcRing;
use crate::contact::proper_crossing;
use crate::exact_arc::boxes::BoxTree;
use crate::exact_arc::{arrangement, Bounds};
use crate::exact_overlay::on_one_line;
use crate::{orientation, OverlayError, Polygon};

/// Which polygon a ring belongs to, and whether it is one of its holes.
#[derive(Debug, Clone, Copy)]
struct Role {
    polygon: usize,
    hole: bool,
}

/// Refuse `polygons` unless every hole lies inside its own outer ring and
/// in no other of its holes, and no point lies in two polygons (module
/// docs). Every ring must already pass [`crate::validate_ring`].
pub(crate) fn check(polygons: &[Polygon]) -> Result<(), OverlayError> {
    if polygons.iter().map(|p| 1 + p.holes.len()).sum::<usize>() < 2 {
        return Ok(());
    }
    let mut rings: Vec<Vec<Point2>> = Vec::new();
    let mut roles: Vec<Role> = Vec::new();
    for (polygon, p) in polygons.iter().enumerate() {
        let all = std::iter::once((&p.outer, false)).chain(p.holes.iter().map(|h| (h, true)));
        for (ring, hole) in all {
            // A ring on one line encloses nothing (#219).
            if on_one_line(&ring.points) {
                continue;
            }
            let mut points = ring.points.clone();
            // Counter-clockwise: a ring holds the left side of a piece it
            // runs the same way.
            if orientation(&points) == Sign::Negative {
                points.reverse();
            }
            rings.push(points);
            roles.push(Role { polygon, hole });
        }
    }
    let mut outer_of: Vec<Option<usize>> = vec![None; polygons.len()];
    for (index, role) in roles.iter().enumerate() {
        if !role.hole {
            outer_of[role.polygon] = Some(index);
        }
    }
    // An outer ring on one line encloses no hole.
    if roles
        .iter()
        .any(|r| r.hole && outer_of[r.polygon].is_none())
    {
        return Err(OverlayError::HoleOutsideOuter);
    }
    if rings.len() < 2 {
        return Ok(());
    }
    // A crossing names its defect by the rings' roles alone when it is a
    // hole's with its own outer ring or another of its polygon's holes, or
    // two outer rings' of polygons without holes: refused before the
    // arrangement cuts what may be a soup's many crossings.
    let mut holed = vec![false; polygons.len()];
    for role in roles.iter().filter(|r| r.hole) {
        holed[role.polygon] = true;
    }
    let named = |a: Role, b: Role| {
        if a.polygon == b.polygon {
            true
        } else {
            !a.hole && !b.hole && !holed[a.polygon] && !holed[b.polygon]
        }
    };
    if let Some((a, b)) = first_crossing(&rings, |a, b| named(roles[a], roles[b])) {
        return Err(crossed(roles[a], roles[b]));
    }
    sides(&rings, &roles, &outer_of)
}

/// The refusal for two rings that cross, by their roles.
fn crossed(a: Role, b: Role) -> OverlayError {
    if a.polygon != b.polygon {
        OverlayError::PolygonsOverlap
    } else if a.hole && b.hole {
        OverlayError::HolesOverlap
    } else {
        OverlayError::HoleOutsideOuter
    }
}

/// The rings of the first pair of edges of different rings that cross
/// properly, by edge order.
fn first_crossing(
    rings: &[Vec<Point2>],
    named: impl Fn(usize, usize) -> bool,
) -> Option<(usize, usize)> {
    let edges: Vec<(usize, usize)> = rings
        .iter()
        .enumerate()
        .flat_map(|(ring, points)| (0..points.len()).map(move |edge| (ring, edge)))
        .collect();
    let ends = |(ring, edge): (usize, usize)| {
        let points = &rings[ring];
        (points[edge], points[(edge + 1) % points.len()])
    };
    let boxes: Vec<Bounds> = edges
        .iter()
        .map(|&e| {
            let (a, b) = ends(e);
            Bounds::around(&[a, b], 0.0)
        })
        .collect();
    let tree = BoxTree::new(boxes.clone());
    let (mut found, mut stack) = (Vec::new(), Vec::new());
    for (i, own) in boxes.iter().enumerate() {
        tree.query_into(|b| b.overlaps(own), &mut found, &mut stack);
        let (a, b) = ends(edges[i]);
        for &j in found
            .iter()
            .filter(|&&j| j > i && edges[j].0 != edges[i].0 && named(edges[i].0, edges[j].0))
        {
            let (c, d) = ends(edges[j]);
            if proper_crossing(a, b, c, d) {
                return Some((edges[i].0, edges[j].0));
            }
        }
    }
    None
}

/// The number of polygons holding a side held by the rings `held`: a
/// polygon holds it when its outer ring does and none of its holes does.
fn polygons_holding(held: &[usize], roles: &[Role]) -> usize {
    held.iter()
        .filter(|&&r| {
            let role = roles[r];
            !role.hole
                && !held
                    .iter()
                    .any(|&h| roles[h].hole && roles[h].polygon == role.polygon)
        })
        .count()
}

/// Every side of every piece of the rings' exact arrangement, read from the
/// rings that hold it.
fn sides(
    rings: &[Vec<Point2>],
    roles: &[Role],
    outer_of: &[Option<usize>],
) -> Result<(), OverlayError> {
    #[cfg(test)]
    tests::ARRANGED.with(|count| count.set(count.get() + 1));
    let arcs: Vec<ArcRing> = rings.iter().map(|r| ArcRing::from_points(r)).collect();
    let raw = arrangement::build(&arcs)?;
    let (mut outside, mut crowded, mut overlap) = (false, false, false);
    let mut held: Vec<usize> = Vec::new();
    for edge in &raw.edges {
        for left in [true, false] {
            held.clear();
            held.extend_from_slice(&edge.inside);
            held.extend(edge.sources.iter().filter(|s| s.2 == left).map(|s| s.0));
            held.sort_unstable();
            held.dedup();
            overlap |= polygons_holding(&held, roles) > 1;
            // Each hole whose inside this side is: inside its outer ring,
            // and in no other of its holes.
            for source in edge.sources.iter().filter(|s| s.2 == left) {
                let role = roles[source.0];
                if !role.hole {
                    continue;
                }
                let outer = outer_of[role.polygon].expect("every hole's polygon has an outer ring");
                outside |= held.binary_search(&outer).is_err();
                crowded |= held
                    .iter()
                    .any(|&r| r != source.0 && roles[r].hole && roles[r].polygon == role.polygon);
            }
        }
    }
    if outside {
        Err(OverlayError::HoleOutsideOuter)
    } else if crowded {
        Err(OverlayError::HolesOverlap)
    } else if overlap {
        Err(OverlayError::PolygonsOverlap)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Ring;
    use std::cell::Cell;
    thread_local! {
        /// Arrangements built on this thread.
        pub(super) static ARRANGED: Cell<usize> = const { Cell::new(0) };
    }
    fn arranged() -> usize {
        ARRANGED.with(Cell::get)
    }
    fn polygon(points: &[(f64, f64)], holes: &[&[(f64, f64)]]) -> Polygon {
        let ring = |p: &[(f64, f64)]| Ring {
            points: p.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
        };
        Polygon {
            outer: ring(points),
            holes: holes.iter().map(|h| ring(h)).collect(),
        }
    }
    #[test]
    fn a_soup_of_crossing_slivers_is_refused_before_it_is_cut() {
        let slivers: Vec<Polygon> = (0..200)
            .map(|k| {
                let y = f64::from(k) * 1e-4;
                polygon(&[(0.0, y), (10.0, y + 0.5), (0.0, y + 1.0)], &[])
            })
            .collect();
        let before = arranged();
        assert_eq!(check(&slivers), Err(OverlayError::PolygonsOverlap));
        assert_eq!(arranged(), before);
        let ring = polygon(
            &[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
            &[&[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)]],
        );
        assert_eq!(check(&[ring]), Ok(()));
        assert_eq!(arranged(), before + 1);
    }
}
