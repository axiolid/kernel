//! A profile swept with a fixed reference along a tangent-continuous chain
//! of segments and arcs (#122).
//!
//! A fixed-reference sweep places the profile's x axis along the reference
//! direction projected square to the directrix tangent, z along the tangent
//! and y = z x x. Along a chain of segments and circular arcs that frame is
//! the rotation-minimising frame [`super`] carries for a swept disk
//! whenever, on every arc, the reference is normal to the arc's plane or
//! lies in it:
//!
//! - Normal to the plane, the projected reference is the plane's normal
//!   itself, which stays put while the tangent turns about it.
//! - In the plane, the projected reference is a direction in the plane
//!   that turns with the tangent, and the plane's normal is the profile's y.
//!
//! In both cases consecutive stations of an arc are rotations of one ring
//! about the arc's axis, so each arc is an exact revolution of the placed
//! profile, and each segment is an exact prism. [`super`]'s per-piece proof
//! (the span bound with the profile's reach, joint stations shared between
//! pieces, the measured joint shift) then carries over with the profile's
//! rings in place of the disk's: every point of the exact swept surface lies
//! within `chord / 2` of the triangles, a lone segment within nothing. The
//! caller adds the rings' own distance from a curved profile.
//!
//! # The twist term
//!
//! A reference that only nearly is normal to the plane or in it is no
//! exception to be waved through: the exact frame then turns about the
//! tangent off the rotation-minimising stations by an angle `delta` that
//! grows along the arc. Both frames lie in the plane of the arc's axis `n`
//! and the tangent's in-plane normal `nu`; the rotation-minimising one at
//! fixed coordinates `(a, q(0))`, the exact one at `(a, q(phi))`, with `a`
//! the reference's component along `n` and `q` its component along `nu`.
//! The angle between them is monotone in `q`, so its largest value is read
//! at the arc's ends and where `q` has an extremum (`arc_twist`). A profile
//! vertex at `reach` from the origin moves by at most `2 reach sin(delta/2)`
//! between the two, and that is added to the walls' share of the budget. A
//! term above an eighth of the budget leaves the arc to the sampled sweep
//! ([`FixedReferencePlan::Twisting`], reported unbounded).
//!
//! Each piece starts from the exact frame at its own start, not carried over
//! a joint. At a kink `theta`, within the joint tolerance, the re-seeded
//! frame turns about the tangent by about `theta cot(beta)` against the
//! previous piece's end, `beta` being the reference's angle to the tangent.
//! That and the previous arc's twist at its end are paid by the first span
//! of the next piece; above a quarter of the budget the joint is refused as
//! a corner for this reference.
//!
//! # Refusals
//!
//! The IFC rule is that the directrix of a fixed-reference sweep shall be
//! tangent continuous, so a corner (a joint turning by more than
//! [`joint_tolerance`] at the reach of the profile) is refused by name, as
//! is a gap between pieces above a quarter of the budget. The reference may
//! never come within the angular tolerance of the tangent on a piece, or
//! the frame degenerates there: a segment is checked whole, an arc over its
//! turn ([`reference_clearance`]). A profile that reaches the axis of a bend
//! is refused: the solid would fold through it, and the revolution would
//! pinch.
//!
//! An arc oblique to the reference twists the frame about the tangent by
//! more than the term allows: [`fixed_reference_plan`] reports it as
//! [`FixedReferencePlan::Twisting`] and the caller sweeps it as sampled and
//! reports it unbounded; nothing here vouches for it.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Scalar, Tolerance, Vec3};
use axiolid_mesh::TriMesh;

use super::{bounded_piece, joint_tolerance, turn, Close, PathPiece, Sampled, Station};
use crate::loft::{self, Frame};
use crate::profile::Rings;

/// How [`fixed_reference_plan`] reads a chain.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum FixedReferencePlan {
    /// No arc twists the frame: swept by [`fixed_reference_along_pieces`],
    /// the bound proven.
    Untwisted(FixedReferenceChain),
    /// Piece `piece` is an arc oblique to the reference: the frame twists
    /// there; not bounded here.
    Twisting {
        /// Index of the first such arc.
        piece: usize,
    },
}

/// A chain checked for a fixed-reference sweep; only [`fixed_reference_plan`]
/// makes one, so a chain in hand has passed its checks.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct FixedReferenceChain {
    pieces: Vec<PathPiece>,
    reference: Vec3,
    closed: bool,
    chord: Scalar,
    reach: Scalar,
    walls: Scalar,
}

impl FixedReferenceChain {
    /// The pieces, in the direction of travel.
    #[must_use]
    pub fn pieces(&self) -> &[PathPiece] {
        &self.pieces
    }

    /// The unit reference direction.
    #[must_use]
    pub const fn reference(&self) -> Vec3 {
        self.reference
    }

    /// Whether the last piece joins the first as a loop (see
    /// [`fixed_reference_plan`]).
    #[must_use]
    pub const fn closed(&self) -> bool {
        self.closed
    }

    /// The chord budget of the walls.
    #[must_use]
    pub const fn chord(&self) -> Scalar {
        self.chord
    }

    /// The profile's furthest ring vertex from its origin.
    #[must_use]
    pub const fn reach(&self) -> Scalar {
        self.reach
    }

    /// The proven distance of the walls from the exact surface: zero for
    /// exactly one segment, else half the chord budget plus the largest
    /// twist term of an arc (see the module notes).
    #[must_use]
    pub const fn walls(&self) -> Scalar {
        self.walls
    }
}

/// The furthest ring vertex from the profile's origin.
#[must_use]
pub fn ring_reach(rings: &Rings) -> Scalar {
    rings
        .outer
        .iter()
        .chain(rings.holes.iter().flatten())
        .map(|p| p.x.hypot(p.y))
        .fold(0.0, Scalar::max)
}

fn refuse(what: String) -> GeomError {
    GeomError::InvalidInput(format!("fixed-reference sweep: {what}"))
}

/// The largest share of the chord budget the exact frame's twist off the
/// rotation-minimising one may take before an arc is left to the sampled
/// sweep, and the share of it a joint's re-seeding may take together with
/// its kink.
const TWIST_SHARE: Scalar = 0.125;
const JOINT_SHARE: Scalar = 0.25;

/// The turn about the tangent between the exact fixed-reference frame and
/// the rotation-minimising frame carried from the arc's start, largest over
/// the arc and at its end.
///
/// With `n` the axis, `a` the reference's component along it and `q(phi)`
/// the reference's component along the in-plane normal `nu(phi)` of the
/// tangent, both frames lie in the plane of `n` and `nu`: the rotation-
/// minimising one at fixed coordinates `(a, q(0))`, the exact one at
/// `(a, q(phi))`, `q(phi) = rho sin(alpha - phi)` for the in-plane part
/// `rho`. The angle is monotone in `q`, so its largest value is at an end of
/// the arc or where `q` has an extremum.
fn arc_twist(piece: &PathPiece, r: Vec3) -> (Scalar, Scalar) {
    let PathPiece::Arc { axis, angle, .. } = *piece else {
        return (0.0, 0.0);
    };
    let t0 = piece.start_tangent();
    let e2 = (axis * angle.signum()).cross(t0);
    let p = r - axis * r.dot(axis);
    let a = r.dot(axis).abs();
    let (qc, qs) = (p.dot(e2), p.dot(t0));
    let q = |phi: Scalar| qc * phi.cos() - qs * phi.sin();
    let theta = |phi: Scalar| q(phi).atan2(a);
    let turned = angle.abs();
    let alpha = qc.atan2(qs);
    let first = (alpha - core::f64::consts::FRAC_PI_2).rem_euclid(core::f64::consts::PI);
    let start = theta(0.0);
    let end = (theta(turned) - start).abs();
    let largest = (0..3)
        .map(|j| first + j as Scalar * core::f64::consts::PI)
        .filter(|phi| *phi <= turned)
        .map(|phi| (theta(phi) - start).abs())
        .fold(end, Scalar::max);
    let safe = 1.0 + 1e-9;
    (largest * safe, end * safe)
}

/// How far a profile vertex at `reach` moves when turned by `delta` about
/// the tangent.
fn twist_term(reach: Scalar, delta: Scalar) -> Scalar {
    2.0 * reach * (0.5 * delta).min(core::f64::consts::FRAC_PI_2).sin()
}

/// The smallest angle between `reference` and the line of the tangent over
/// the piece, in `[0, pi/2]`. Zero where the frame degenerates.
#[must_use]
pub fn reference_clearance(piece: &PathPiece, reference: Vec3) -> Scalar {
    let r = reference.normalize_or_zero();
    let between = |t: Vec3| r.cross(t).length().atan2(r.dot(t).abs());
    match *piece {
        PathPiece::Segment { .. } => between(piece.start_tangent()),
        PathPiece::Arc { axis, angle, .. } => {
            let t0 = piece.start_tangent();
            let k = axis * angle.signum();
            let e2 = k.cross(t0);
            // The tangent at turn `phi` is `t0 cos(phi) + e2 sin(phi)`: it
            // is parallel to the reference's in-plane part `p` where `phi`
            // is `alpha` plus a multiple of pi.
            let p = r - axis * r.dot(axis);
            let alpha = p.dot(e2).atan2(p.dot(t0));
            let turned = angle.abs();
            if alpha.rem_euclid(core::f64::consts::PI) <= turned {
                return r.dot(axis).abs().atan2(p.length());
            }
            let at = |phi: Scalar| t0 * phi.cos() + e2 * phi.sin();
            let end = if (turned - alpha).cos().abs() > alpha.cos().abs() {
                turned
            } else {
                0.0
            };
            between(at(end))
        }
    }
}

/// Check a chain of pieces for a fixed-reference sweep of a profile of
/// `reach` at wall budget `chord` and decide whether the bound is proven.
///
/// `closed` is whether the path declares itself a loop (a closed polyline):
/// its closing joint is then checked like any other. An undeclared chain is
/// a loop all the same when its last piece ends where its first starts, to
/// rounding, and the joint there is tangent continuous (a full circle, a
/// composite closing on itself); where it is not, the ends just meet and the
/// chain is open, nothing refused.
///
/// Refused by name, in this order: no pieces or a reference, chord or reach
/// that is zero or not finite; a malformed piece; a joint with a gap above a
/// quarter of `chord` or a corner; a piece whose tangent comes within the
/// angular tolerance of the reference's direction; a joint whose kink
/// turns the re-seeded profile (by about the kink times the cotangent of the
/// reference's angle to the tangent) further than a quarter of `chord` at
/// `reach`. An arc whose exact frame twists off the rotation-minimising one
/// by more than an eighth of `chord` at `reach` is reported
/// [`FixedReferencePlan::Twisting`]; a smaller twist is added to the walls.
pub fn fixed_reference_plan(
    pieces: &[PathPiece],
    closed: bool,
    reference: Vec3,
    reach: Scalar,
    chord: Scalar,
    tolerance: Tolerance,
) -> GeomResult<FixedReferencePlan> {
    if pieces.is_empty() {
        return Err(refuse("the directrix needs at least one piece".to_owned()));
    }
    let r = reference.normalize_or_zero();
    if !reference.is_finite() || r == Vec3::ZERO {
        return Err(refuse(format!(
            "the reference direction must be finite and non-zero, got {reference:?}"
        )));
    }
    for (what, value) in [("chord budget", chord), ("profile reach", reach)] {
        if !(value.is_finite() && value > 0.0) {
            return Err(refuse(format!(
                "the {what} must be positive and finite, got {value}"
            )));
        }
    }
    for (i, piece) in pieces.iter().enumerate() {
        check_shape(i, piece)?;
    }
    let angular = tolerance.angular();
    let m = pieces.len();
    let allowed = joint_tolerance(chord, reach);
    let rounding = 1e-12 * pieces.iter().map(PathPiece::length).sum::<Scalar>();
    let (last, first) = (&pieces[m - 1], &pieces[0]);
    let loops = closed
        || (last.end_point().distance(first.start_point()) <= rounding
            && turn(last.end_tangent(), first.start_tangent()) <= allowed);
    let joints = if loops { m } else { m - 1 };
    for k in 0..joints {
        let (a, b) = (&pieces[k], &pieces[(k + 1) % m]);
        let gap = a.end_point().distance(b.start_point());
        if gap > 0.25 * chord {
            return Err(refuse(format!(
                "joint {k} leaves a gap of {gap}, above a quarter of the chord budget {chord}"
            )));
        }
        let theta = turn(a.end_tangent(), b.start_tangent());
        if theta > allowed {
            return Err(refuse(format!(
                "the directrix turns a corner of {theta} rad at joint {k}, above {allowed}: the \
                 directrix must be tangent continuous"
            )));
        }
    }
    for (i, piece) in pieces.iter().enumerate() {
        let clearance = reference_clearance(piece, r);
        if clearance <= angular {
            return Err(refuse(format!(
                "the reference direction comes within {clearance} rad of the directrix tangent \
                 on piece {i}"
            )));
        }
    }
    let mut twist: Scalar = 0.0;
    let mut ends = vec![0.0; m];
    for (i, piece) in pieces.iter().enumerate() {
        let (largest, end) = arc_twist(piece, r);
        let term = twist_term(reach, largest);
        if term > TWIST_SHARE * chord {
            return Ok(FixedReferencePlan::Twisting { piece: i });
        }
        twist = twist.max(term);
        ends[i] = twist_term(reach, end);
    }
    for k in 0..joints {
        let (a, b) = (&pieces[k], &pieces[(k + 1) % m]);
        let (ta, tb) = (a.end_tangent(), b.start_tangent());
        let xa = Frame::from_reference(a.end_point(), ta, r)?.x;
        let xb = Frame::from_reference(b.start_point(), tb, r)?.x;
        let kink = turn(ta, tb);
        let reseed = turn(super::least_rotation(xa, ta, tb), xb);
        let shift = reach * (kink + reseed) + ends[k];
        if shift > JOINT_SHARE * chord {
            return Err(refuse(format!(
                "joint {k} turns the profile by {reseed} rad for this reference (a kink of \
                 {kink} rad), moving it {shift}: a corner for this reference, above a quarter \
                 of the chord budget {chord}"
            )));
        }
    }
    let walls = if matches!(pieces, [PathPiece::Segment { .. }]) {
        0.0
    } else {
        0.5 * chord + twist
    };
    Ok(FixedReferencePlan::Untwisted(FixedReferenceChain {
        pieces: pieces.to_vec(),
        reference: r,
        closed: loops,
        chord,
        reach,
        walls,
    }))
}

/// Refuse a piece with no direction or a malformed arc.
fn check_shape(i: usize, piece: &PathPiece) -> GeomResult<()> {
    match *piece {
        PathPiece::Segment { start, end } => {
            if !(start.is_finite() && end.is_finite()) {
                return Err(refuse(format!("piece {i} must be finite")));
            }
            if start == end {
                return Err(GeomError::Degenerate(format!(
                    "fixed-reference sweep: segment {i} has no length"
                )));
            }
        }
        PathPiece::Arc {
            centre,
            axis,
            start,
            angle,
        } => {
            if !(centre.is_finite() && axis.is_finite() && start.is_finite() && angle.is_finite()) {
                return Err(refuse(format!("piece {i} must be finite")));
            }
            let bend = (start - centre).length();
            if (axis.length() - 1.0).abs() > 1e-9
                || bend == 0.0
                || (start - centre).dot(axis).abs() > 1e-9 * bend
            {
                return Err(refuse(format!(
                    "arc {i} needs a unit axis normal to a non-zero radius"
                )));
            }
            if angle == 0.0 || angle.abs() > core::f64::consts::TAU * (1.0 + 1e-12) {
                return Err(refuse(format!(
                    "arc {i} must turn by a non-zero angle of at most a full turn, got {angle}"
                )));
            }
        }
    }
    Ok(())
}

/// The budget errors a piece's refinement raises, worded for this sweep.
fn named(error: GeomError) -> GeomError {
    match error {
        GeomError::BudgetExceeded {
            resource: "swept disk directrix refinement",
        } => GeomError::BudgetExceeded {
            resource: "fixed-reference sweep directrix refinement",
        },
        GeomError::BudgetExceeded {
            resource: "swept disk bend steps",
        } => GeomError::BudgetExceeded {
            resource: "fixed-reference sweep bend steps",
        },
        other => other,
    }
}

/// Sweep `rings` along a planned chain so that every point of the exact
/// surface lies within `chain.walls()` of the mesh (see the module notes).
///
/// Each piece starts from the fixed-reference frame at its own start, never
/// carried across a joint; consecutive pieces share one station, and a loop
/// closes on its first station without caps.
pub fn fixed_reference_along_pieces(
    rings: &Rings,
    chain: &FixedReferenceChain,
) -> GeomResult<TriMesh> {
    let reference = chain.reference();
    let share = 0.5 * chain.chord();
    for (i, piece) in chain.pieces().iter().enumerate() {
        if let PathPiece::Arc { centre, start, .. } = *piece {
            let frame = Frame::from_reference(start, piece.start_tangent(), reference)?;
            let bend = (start - centre).length();
            let inward = (centre - start).normalize_or_zero();
            let depth = rings
                .outer
                .iter()
                .chain(rings.holes.iter().flatten())
                .map(|p| (frame.x * p.x + frame.y * p.y).dot(inward))
                .fold(Scalar::NEG_INFINITY, Scalar::max)
                + share;
            if depth >= bend {
                return Err(refuse(format!(
                    "the profile reaches {depth} toward the axis of arc {i}, not below its bend \
                     radius {bend}: the solid would fold through the bend's axis"
                )));
            }
        }
    }
    let m = chain.pieces().len();
    let mut stations: Vec<Station> = Vec::new();
    for (i, piece) in chain.pieces().iter().enumerate() {
        let seed = Frame::from_reference(piece.start_point(), piece.start_tangent(), reference)?.x;
        let closing = chain.closed() && i + 1 == m;
        let close = match (closing, stations.first()) {
            (false, _) => Close::Open,
            (true, Some(first)) => Close::To(first),
            (true, None) => Close::Own,
        };
        let sample = |budget: Scalar| -> GeomResult<Sampled> {
            let n = piece.steps(budget)?;
            let points: Vec<_> = (0..=n)
                .map(|k| piece.point(k as Scalar / n as Scalar))
                .collect();
            let tangents = points.iter().map(|p| piece.tangent_at(*p)).collect();
            Ok(Sampled {
                twists: vec![0.0; points.len()],
                points,
                tangents,
            })
        };
        let (_, placed) = bounded_piece(
            sample,
            rings,
            seed,
            chain.reach(),
            share,
            stations.last(),
            [None, None],
            close,
        )
        .map_err(named)?;
        let skip = usize::from(!stations.is_empty());
        let keep = placed.len() - usize::from(closing);
        stations.extend(placed.into_iter().take(keep).skip(skip));
    }
    let mut mesh = loft::loft(rings, &stations, chain.closed())?;
    if matches!(crate::extrude::outward_orientation(&mesh), Some(false)) {
        for triangle in mesh.indices.as_chunks_mut::<3>().0 {
            triangle.swap(1, 2);
        }
    }
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

    use axiolid_core::{Point2, Point3};
    use axiolid_measure::volume_properties;

    use super::*;

    const CHORD: Scalar = 1e-3;

    fn tol() -> Tolerance {
        Tolerance::new(1e-6, 1e-2).unwrap()
    }

    fn square(lo: Scalar, hi: Scalar) -> Rings {
        Rings {
            outer: vec![
                Point2::new(lo, lo),
                Point2::new(hi, lo),
                Point2::new(hi, hi),
                Point2::new(lo, hi),
            ],
            holes: Vec::new(),
        }
    }

    fn segment(from: [Scalar; 3], to: [Scalar; 3]) -> PathPiece {
        PathPiece::Segment {
            start: Point3::new(from[0], from[1], from[2]),
            end: Point3::new(to[0], to[1], to[2]),
        }
    }

    /// A unit arc about Z from (1, 0, 0), counter-clockwise for a positive
    /// `angle`.
    fn arc(angle: Scalar) -> PathPiece {
        PathPiece::Arc {
            centre: Point3::ZERO,
            axis: Vec3::Z,
            start: Point3::new(1.0, 0.0, 0.0),
            angle,
        }
    }

    fn plan(pieces: &[PathPiece], reference: Vec3) -> GeomResult<FixedReferencePlan> {
        fixed_reference_plan(pieces, false, reference, 0.07, CHORD, tol())
    }

    fn refusal(result: GeomResult<impl core::fmt::Debug>) -> String {
        let message = result.unwrap_err().to_string();
        assert!(message.contains("fixed-reference sweep"), "{message}");
        message
    }

    #[test]
    fn a_segment_is_cleared_by_its_angle_to_the_reference_line() {
        let s = segment([0.0; 3], [2.0, 0.0, 0.0]);
        let r = Vec3::new(1.0, 1.0, 0.0);
        assert!((reference_clearance(&s, r) - FRAC_PI_4).abs() < 1e-12);
        assert!((reference_clearance(&s.reversed(), r) - FRAC_PI_4).abs() < 1e-12);
        assert!((reference_clearance(&s, -r) - FRAC_PI_4).abs() < 1e-12);
        assert!(reference_clearance(&s, Vec3::X) < 1e-12);
    }

    #[test]
    fn an_in_plane_arc_is_cleared_by_where_its_tangent_comes_closest() {
        let r = Vec3::new(1.0, 1.0, 0.0);
        // Crossing: the tangent is parallel to the reference at 3 pi / 4.
        assert!(reference_clearance(&arc(PI), r) < 1e-12);
        // Not crossing: the start is nearest.
        assert!((reference_clearance(&arc(FRAC_PI_4), r) - FRAC_PI_4).abs() < 1e-12);
        // Ending at it: the tangent at pi / 2 is -x.
        assert!(reference_clearance(&arc(FRAC_PI_2), Vec3::X) < 1e-12);
        // Walked the other way, the same set of tangents.
        assert!(reference_clearance(&arc(-PI), Vec3::new(1.0, -1.0, 0.0)) < 1e-12);
        // Leaning out of the plane: the angle to the in-plane part.
        let lean = Vec3::new(0.0, 1.0, 1.0);
        assert!((reference_clearance(&arc(1.0), lean) - FRAC_PI_4).abs() < 1e-12);
        // Normal to the plane: square to every tangent.
        assert!((reference_clearance(&arc(1.0), Vec3::Z) - FRAC_PI_2).abs() < 1e-12);
    }

    #[test]
    fn an_arc_pays_for_the_twist_of_its_lean() {
        let (reach, exact) = (0.07, 0.5 * CHORD);
        let walls = |pieces: &[PathPiece], reference: Vec3| match fixed_reference_plan(
            pieces,
            false,
            reference,
            reach,
            CHORD,
            tol(),
        )
        .unwrap()
        {
            FixedReferencePlan::Untwisted(chain) => Some(chain.walls()),
            FixedReferencePlan::Twisting { .. } => None,
        };
        let pieces = [arc(1.0)];
        assert_eq!(walls(&pieces, Vec3::Z), Some(exact));
        for lean in [
            |l: Scalar| Vec3::new(0.0, l.sin(), l.cos()),
            |l: Scalar| Vec3::new(l.sin(), 0.0, l.cos()),
        ] {
            // A small lean is paid for: at most the rim's chord at twice
            // the lean, and more than nothing.
            let small = 1e-4;
            let w = walls(&pieces, lean(small)).expect("untwisted");
            assert!(w > exact && w <= exact + reach * 2.0 * small * 1.01, "{w}");
            // A large one is left to the sampled sweep.
            assert_eq!(walls(&pieces, lean(1e-2)), None);
        }
        let twisting = plan(&pieces, Vec3::new(1.0, 1.0, 1.0)).unwrap();
        assert_eq!(twisting, FixedReferencePlan::Twisting { piece: 0 });
    }

    #[test]
    fn the_twist_is_the_measured_turn_of_the_exact_frame() {
        // A reference leaning from the normal by `lean`: over a quarter turn
        // the exact frame turns off the rotation-minimising one by the lean.
        let lean: Scalar = 9e-4;
        let r = Vec3::new(lean.sin(), 0.0, lean.cos());
        let (largest, end) = arc_twist(&arc(FRAC_PI_2), r);
        assert!((largest - lean).abs() < 1e-9 && (end - lean).abs() < 1e-9);
        // Over a full turn it peaks at twice the lean, where q changes sign.
        let (largest, _) = arc_twist(&arc(TAU), r);
        assert!((largest - 2.0 * lean).abs() < 1e-8, "{largest}");
        // Exactly normal: none.
        assert_eq!(arc_twist(&arc(1.0), Vec3::Z), (0.0, 0.0));
    }

    #[test]
    fn a_plan_refuses_by_name() {
        let bad_inputs = [
            (plan(&[], Vec3::Z), "at least one piece"),
            (plan(&[arc(1.0)], Vec3::ZERO), "non-zero"),
            (
                plan(&[arc(1.0)], Vec3::new(Scalar::NAN, 0.0, 1.0)),
                "finite",
            ),
            (
                fixed_reference_plan(&[arc(1.0)], false, Vec3::Z, 0.07, 0.0, tol()),
                "chord budget",
            ),
            (
                fixed_reference_plan(&[arc(1.0)], false, Vec3::Z, Scalar::NAN, CHORD, tol()),
                "profile reach",
            ),
            (plan(&[segment([0.0; 3], [0.0; 3])], Vec3::Z), "no length"),
            (plan(&[arc(0.0)], Vec3::Z), "arc 0"),
            (plan(&[arc(7.0)], Vec3::Z), "full turn"),
        ];
        for (result, needle) in bad_inputs {
            let message = refusal(result);
            assert!(message.contains(needle), "{needle}: {message}");
        }
        // A corner, by joint.
        let corner = [
            segment([0.0; 3], [1.0, 0.0, 0.0]),
            segment([1.0, 0.0, 0.0], [1.0, 1.0, 0.0]),
        ];
        let message = refusal(plan(&corner, Vec3::Z));
        assert!(
            message.contains("joint 0") && message.contains("tangent continuous"),
            "{message}"
        );
        // The corner a loop closes with.
        let triangle = [
            segment([0.0; 3], [1.0, 0.0, 0.0]),
            segment([1.0, 0.0, 0.0], [1.0, 1.0, 0.0]),
            segment([1.0, 1.0, 0.0], [0.0; 3]),
        ];
        let message = refusal(fixed_reference_plan(
            &triangle,
            true,
            Vec3::Z,
            0.07,
            CHORD,
            tol(),
        ));
        assert!(message.contains("joint 0"), "{message}");
        // A gap.
        let gap = [
            segment([0.0; 3], [1.0, 0.0, 0.0]),
            segment([1.5, 0.0, 0.0], [2.0, 0.0, 0.0]),
        ];
        let message = refusal(plan(&gap, Vec3::Z));
        assert!(
            message.contains("gap") && message.contains("joint 0"),
            "{message}"
        );
        // A reference parallel to a segment.
        let line = [
            segment([0.0; 3], [1.0, 0.0, 0.0]),
            segment([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]),
        ];
        let message = refusal(plan(&line, Vec3::X));
        assert!(
            message.contains("piece 0") && message.contains("tangent"),
            "{message}"
        );
        // An in-plane arc whose tangent reaches the reference.
        let message = refusal(plan(&[arc(PI + 0.5)], Vec3::new(1.0, 1.0, 0.0)));
        assert!(message.contains("piece 0"), "{message}");
        // A kink within the joint tolerance that turns the re-seeded profile
        // far: the reference nearly along the tangent.
        let kinked = [
            segment([0.0; 3], [1.0, 0.0, 0.0]),
            segment([1.0, 0.0, 0.0], [2.0, 0.003, 0.0]),
        ];
        let flat = Vec3::new(0.05f64.cos(), 0.0, 0.05f64.sin());
        let message = refusal(plan(&kinked, flat));
        assert!(
            message.contains("joint 0") && message.contains("corner for this reference"),
            "{message}"
        );
        assert!(plan(&kinked, Vec3::Z).is_ok());
    }

    #[test]
    fn walls_are_free_only_along_one_segment() {
        let one = segment([0.0; 3], [1.0, 0.0, 0.0]);
        let two = segment([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        let walls = |pieces: &[PathPiece]| match plan(pieces, Vec3::Z).unwrap() {
            FixedReferencePlan::Untwisted(chain) => chain.walls(),
            other => panic!("{other:?}"),
        };
        assert_eq!(walls(&[one]), 0.0);
        assert_eq!(walls(&[one, two]), 0.5 * CHORD);
        assert_eq!(walls(&[arc(1.0)]), 0.5 * CHORD);
    }

    fn untwisted(pieces: &[PathPiece], closed: bool, reference: Vec3) -> FixedReferenceChain {
        let reach = ring_reach(&square(-0.05, 0.05));
        match fixed_reference_plan(pieces, closed, reference, reach, CHORD, tol()).unwrap() {
            FixedReferencePlan::Untwisted(chain) => chain,
            other => panic!("{other:?}"),
        }
    }

    fn volume(mesh: &TriMesh) -> Scalar {
        volume_properties(mesh, Tolerance::new(1e-9, 1e-9).unwrap())
            .expect("a closed two-manifold")
            .signed_volume
    }

    #[test]
    fn a_full_circle_with_a_normal_reference_is_a_closed_torus() {
        let chain = untwisted(&[arc(TAU)], false, Vec3::Z);
        assert!(chain.closed());
        let mesh = fixed_reference_along_pieces(&square(-0.05, 0.05), &chain).unwrap();
        // Area times the centreline length (a centred profile).
        let got = volume(&mesh);
        assert!((got - 0.01 * TAU).abs() < 0.01 * TAU * 1e-2, "{got}");
    }

    #[test]
    fn a_stadium_loop_closes_without_caps() {
        let stadium = [
            segment([0.0, -1.0, 0.0], [2.0, -1.0, 0.0]),
            PathPiece::Arc {
                centre: Point3::new(2.0, 0.0, 0.0),
                axis: Vec3::Z,
                start: Point3::new(2.0, -1.0, 0.0),
                angle: PI,
            },
            segment([2.0, 1.0, 0.0], [0.0, 1.0, 0.0]),
            PathPiece::Arc {
                centre: Point3::ZERO,
                axis: Vec3::Z,
                start: Point3::new(0.0, 1.0, 0.0),
                angle: PI,
            },
        ];
        let chain = untwisted(&stadium, true, Vec3::Z);
        let mesh = fixed_reference_along_pieces(&square(-0.05, 0.05), &chain).unwrap();
        let length = 4.0 + 2.0 * PI;
        let got = volume(&mesh);
        assert!((got - 0.01 * length).abs() < 0.01 * length * 1e-2, "{got}");
    }

    #[test]
    fn a_profile_reaching_a_bends_axis_is_refused() {
        let at = |lowest: Scalar| {
            let rings = Rings {
                outer: vec![
                    Point2::new(-0.05, lowest),
                    Point2::new(0.05, lowest),
                    Point2::new(0.05, 0.05),
                    Point2::new(-0.05, 0.05),
                ],
                holes: Vec::new(),
            };
            let reach = ring_reach(&rings);
            let plan = fixed_reference_plan(&[arc(1.0)], false, Vec3::Z, reach, CHORD, tol());
            let Ok(FixedReferencePlan::Untwisted(chain)) = plan else {
                panic!("{plan:?}")
            };
            fixed_reference_along_pieces(&rings, &chain)
        };
        // In the normal-reference frame the profile's -y points at the axis.
        let message = refusal(at(-1.0));
        assert!(
            message.contains("fold through the bend's axis"),
            "{message}"
        );
        // Within the half-budget slack of the axis.
        refusal(at(-1.0 + 0.25 * CHORD));
        assert!(at(-0.9).is_ok());
    }
}
