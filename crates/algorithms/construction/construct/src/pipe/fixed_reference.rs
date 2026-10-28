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
//! An arc oblique to the reference (neither normal to it nor containing it)
//! twists the frame about the tangent, and the fixed-reference frame is no
//! longer the rotation-minimising one: [`fixed_reference_plan`] reports it
//! as [`FixedReferencePlan::Twisting`] and the caller sweeps it as sampled
//! and reports it unbounded; nothing here vouches for it.

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

/// A chain checked for a fixed-reference sweep.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct FixedReferenceChain {
    /// The pieces, in the direction of travel.
    pub pieces: Vec<PathPiece>,
    /// The unit reference direction.
    pub reference: Vec3,
    /// Whether the last piece ends where the first starts (the path's own
    /// flag, or its end within the linear tolerance of its start).
    pub closed: bool,
    /// The chord budget of the walls.
    pub chord: Scalar,
    /// The profile's furthest ring vertex from its origin.
    pub reach: Scalar,
    /// The proven distance of the walls from the exact surface: zero for
    /// exactly one segment, else half the chord budget.
    pub walls: Scalar,
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

/// How an arc lies to the reference.
enum Lie {
    /// The reference is normal to the arc's plane.
    Normal,
    /// The reference lies in the arc's plane.
    InPlane,
    /// Neither: the frame twists.
    Oblique,
}

fn lie(axis: Vec3, reference: Vec3, angular: Scalar) -> Lie {
    let (across, along) = (reference.cross(axis).length(), reference.dot(axis).abs());
    if across.atan2(along) <= angular {
        Lie::Normal
    } else if along.atan2(across) <= angular {
        Lie::InPlane
    } else {
        Lie::Oblique
    }
}

/// The first arc oblique to `reference` (a direction, normalised here) by
/// more than `angular` from both its normal and its plane: where the frame
/// twists.
#[must_use]
pub fn fixed_reference_twist(
    pieces: &[PathPiece],
    reference: Vec3,
    angular: Scalar,
) -> Option<usize> {
    let r = reference.normalize_or_zero();
    pieces.iter().position(|piece| match *piece {
        PathPiece::Arc { axis, .. } => matches!(lie(axis, r, angular), Lie::Oblique),
        _ => false,
    })
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
/// Refused by name, in this order: no pieces or a reference, chord or reach
/// that is zero or not finite; a malformed piece; a joint with a gap above a
/// quarter of `chord` or a corner (the closing joint of a loop included);
/// a segment, or an arc normal to or containing the reference, that comes
/// within the angular tolerance of the reference's direction.
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
    let loops =
        closed || pieces[m - 1].end_point().distance(pieces[0].start_point()) <= tolerance.linear();
    let allowed = joint_tolerance(chord, reach);
    for k in 0..if loops { m } else { m - 1 } {
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
        let straight = match *piece {
            PathPiece::Arc { axis, .. } => !matches!(lie(axis, r, angular), Lie::Oblique),
            PathPiece::Segment { .. } => true,
        };
        if !straight {
            continue;
        }
        let clearance = reference_clearance(piece, r);
        if clearance <= angular {
            return Err(refuse(format!(
                "the reference direction comes within {clearance} rad of the directrix tangent \
                 on piece {i}"
            )));
        }
    }
    if let Some(piece) = fixed_reference_twist(pieces, r, angular) {
        return Ok(FixedReferencePlan::Twisting { piece });
    }
    let walls = if matches!(pieces, [PathPiece::Segment { .. }]) {
        0.0
    } else {
        0.5 * chord
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

/// Sweep `rings` along a planned chain so that every point of the exact
/// surface lies within `chain.walls` of the mesh (see the module notes).
///
/// Each piece starts from the fixed-reference frame at its own start, never
/// carried across a joint; consecutive pieces share one station, and a loop
/// closes on its first station without caps.
pub fn fixed_reference_along_pieces(
    rings: &Rings,
    chain: &FixedReferenceChain,
) -> GeomResult<TriMesh> {
    let reference = chain.reference;
    let share = 0.5 * chain.chord;
    for (i, piece) in chain.pieces.iter().enumerate() {
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
    let m = chain.pieces.len();
    let mut stations: Vec<Station> = Vec::new();
    for (i, piece) in chain.pieces.iter().enumerate() {
        let seed = Frame::from_reference(piece.start_point(), piece.start_tangent(), reference)?.x;
        let closing = chain.closed && i + 1 == m;
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
            chain.reach,
            share,
            stations.last(),
            [None, None],
            close,
        )?;
        let skip = usize::from(!stations.is_empty());
        let keep = placed.len() - usize::from(closing);
        stations.extend(placed.into_iter().take(keep).skip(skip));
    }
    let mut mesh = loft::loft(rings, &stations, chain.closed)?;
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
    fn an_arc_is_oblique_just_beyond_the_angular_tolerance() {
        let a = tol().angular();
        let lean = |theta: Scalar| Vec3::new(0.0, theta.sin(), theta.cos());
        let tilt = |theta: Scalar| Vec3::new(theta.cos(), 0.0, theta.sin());
        let pieces = [arc(1.0)];
        for (reference, twists) in [
            (lean(a * (1.0 - 1e-6)), false),
            (lean(a * (1.0 + 1e-6)), true),
            (tilt(a * (1.0 - 1e-6)), false),
            (tilt(a * (1.0 + 1e-6)), true),
        ] {
            let found = fixed_reference_twist(&pieces, reference, a);
            assert_eq!(found, twists.then_some(0), "{reference:?}");
        }
        let twisting = plan(&pieces, Vec3::new(1.0, 1.0, 1.0)).unwrap();
        assert_eq!(twisting, FixedReferencePlan::Twisting { piece: 0 });
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
    }

    #[test]
    fn walls_are_free_only_along_one_segment() {
        let one = segment([0.0; 3], [1.0, 0.0, 0.0]);
        let two = segment([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        let walls = |pieces: &[PathPiece]| match plan(pieces, Vec3::Z).unwrap() {
            FixedReferencePlan::Untwisted(chain) => chain.walls,
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
        assert!(chain.closed);
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
