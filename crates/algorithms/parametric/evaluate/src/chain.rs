//! Evaluation of arc-length chains ([`Chain2`]).
//!
//! A chain's parameter is cumulative arc length. Each piece is evaluated in
//! its own local frame and placed rigidly at the end of the one before it,
//! so reading a point resolves every earlier piece's end first:
//!
//! - an [`ChainPiece2::Intrinsic`] piece through [`intrinsic_point`] and
//!   [`intrinsic_tangent`]: heading exact, position Gauss-Legendre
//!   quadrature, as for [`Curve2::Intrinsic`];
//! - a [`ChainPiece2::Parametric`] piece through the arc-length inverse of
//!   [`crate::arc_parameter`]: the curve parameter whose arc length from
//!   `start` is the distance asked for, to [`ARC_LENGTH_TOLERANCE`]
//!   relative.
//!
//! # Accuracy contract
//!
//! A point at distance `s` is the exact chain's point at a distance within
//! `ARC_LENGTH_TOLERANCE * max(1, piece length)` of `s` on each parametric
//! piece it crosses, and within the intrinsic quadrature's accuracy on each
//! intrinsic piece (machine precision on the closed forms it is pinned
//! against, ADR 0060). The tangent is the exact unit tangent at that point:
//! a parametric curve's derivative or an intrinsic piece's closed-form
//! heading. Positions after a parametric piece carry its end error, rotated
//! by the remaining placement; nothing is re-fitted.
//!
//! # Refused by name
//!
//! A malformed chain (no pieces, a non-finite or degenerate start frame, a
//! non-finite or non-positive length, a malformed curvature law), a
//! parametric piece that does not start at its local origin with tangent
//! `+x` within [`PIECE_FRAME_TOLERANCE`], a piece longer than its curve
//! (named with the arc length available), and a distance outside
//! `[0, length]`.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Point2, Scalar, Vec2};
use axiolid_curve::{Chain2, ChainPiece2, CurvatureLaw, Curve2, Intrinsic2};

use crate::arc_length::{intrinsic_point, intrinsic_tangent};
use crate::arc_parameter::{parameter_at_arc_length2, ARC_LENGTH_TOLERANCE};
use crate::curve::{derivative2, evaluate2};

/// How far a parametric piece's curve may sit from its local frame at its
/// start parameter: `PIECE_FRAME_TOLERANCE * max(1, length)` from the
/// origin, and `PIECE_FRAME_TOLERANCE` radians from `+x`.
pub const PIECE_FRAME_TOLERANCE: Scalar = 1e-9;

fn invalid(detail: String) -> GeomError {
    GeomError::InvalidInput(detail)
}

/// Where a piece starts, in the chain's start-frame coordinates: its origin
/// and its unit `+x` direction.
#[derive(Debug, Clone, Copy)]
struct Placement {
    origin: Vec2,
    x: Vec2,
}

impl Placement {
    const IDENTITY: Self = Self {
        origin: Vec2::ZERO,
        x: Vec2::X,
    };

    fn turn(&self, v: Vec2) -> Vec2 {
        Vec2::new(
            self.x.x * v.x - self.x.y * v.y,
            self.x.y * v.x + self.x.x * v.y,
        )
    }

    fn place(&self, p: Vec2) -> Vec2 {
        self.origin + self.turn(p)
    }
}

/// The identity frame an intrinsic piece is integrated in.
fn local_frame() -> Frame2 {
    Frame2 {
        origin: Point2::ZERO,
        x: Vec2::X,
        y: Vec2::Y,
    }
}

fn unit(v: Vec2, what: &str) -> GeomResult<Vec2> {
    let length = v.length();
    if !(length.is_finite() && length > 0.0) {
        return Err(GeomError::Degenerate(format!(
            "{what} has no tangent direction"
        )));
    }
    Ok(v / length)
}

/// Refuse a structurally malformed chain, naming the piece at fault.
fn check_well_formed(chain: &Chain2) -> GeomResult<()> {
    if chain.pieces.is_empty() {
        return Err(invalid("arc-length chain has no pieces".into()));
    }
    let frame = &chain.start;
    if !(frame.origin.is_finite() && frame.x.is_finite() && frame.y.is_finite())
        || frame.x.perp_dot(frame.y) == 0.0
    {
        return Err(invalid(
            "arc-length chain start frame must be finite with independent axes".into(),
        ));
    }
    for (index, piece) in chain.pieces.iter().enumerate() {
        let length = piece.length();
        if !(length.is_finite() && length > 0.0) {
            return Err(invalid(format!(
                "arc-length chain piece {index} has length {length}; it must be finite and positive"
            )));
        }
        match piece {
            ChainPiece2::Intrinsic { curvature, .. } if !curvature.is_well_formed() => {
                return Err(invalid(format!(
                    "arc-length chain piece {index} has a malformed curvature law"
                )));
            }
            ChainPiece2::Parametric { start, .. } if !start.is_finite() => {
                return Err(invalid(format!(
                    "arc-length chain piece {index} has a non-finite start parameter"
                )));
            }
            _ => {}
        }
    }
    if !chain.is_well_formed() {
        return Err(invalid("arc-length chain is malformed".into()));
    }
    Ok(())
}

/// A parametric piece's curve must start at its local origin with tangent
/// `+x`: the chain places it rigidly and never moves it there.
fn check_local_start(
    index: usize,
    curve: &Curve2,
    start: Scalar,
    length: Scalar,
) -> GeomResult<()> {
    let origin = evaluate2(curve, start)?;
    let tangent = unit(derivative2(curve, start)?, "chain piece curve")?;
    let offset = Vec2::new(origin.x, origin.y).length();
    let angle = tangent.y.atan2(tangent.x).abs();
    if offset > PIECE_FRAME_TOLERANCE * length.max(1.0) || angle > PIECE_FRAME_TOLERANCE {
        return Err(invalid(format!(
            "arc-length chain piece {index} does not start at its local origin with tangent +x \
             (offset {offset}, tangent angle {angle} rad)"
        )));
    }
    Ok(())
}

/// Point and unit tangent of one piece at local arc length `s`, in the
/// piece's local frame.
fn piece_local(index: usize, piece: &ChainPiece2, s: Scalar) -> GeomResult<(Vec2, Vec2)> {
    match piece {
        ChainPiece2::Intrinsic { curvature, length } => {
            let curve = Intrinsic2::new(local_frame(), curvature.clone(), *length);
            let p = intrinsic_point(&curve, s)?;
            Ok((Vec2::new(p.x, p.y), intrinsic_tangent(&curve, s)?))
        }
        ChainPiece2::Parametric {
            curve,
            start,
            length,
        } => {
            check_local_start(index, curve, *start, *length)?;
            // The end is checked on every read, so a piece longer than its
            // curve is refused wherever the chain is read, not only past it.
            let end = parameter_at_arc_length2(curve, *start, *length).map_err(|error| {
                invalid(format!(
                    "arc-length chain piece {index} asks for arc length {length} from parameter \
                     {start}: {error}"
                ))
            })?;
            let u = if s == *length {
                end
            } else {
                parameter_at_arc_length2(curve, *start, s)?
            };
            let p = evaluate2(curve, u)?;
            Ok((
                Vec2::new(p.x, p.y),
                unit(derivative2(curve, u)?, "chain piece curve")?,
            ))
        }
        _ => Err(GeomError::Unsupported {
            backend: axiolid_contracts::BackendId::new("axiolid-evaluate"),
            operation: axiolid_contracts::Operation::CurveEvaluation,
        }),
    }
}

/// The piece at distance `s`, its local distance, and its placement.
fn locate(chain: &Chain2, s: Scalar) -> GeomResult<(usize, Scalar, Placement)> {
    check_well_formed(chain)?;
    if !s.is_finite() {
        return Err(invalid(
            "distance along an arc-length chain must be finite".into(),
        ));
    }
    let total = chain
        .length()
        .ok_or_else(|| invalid("arc-length chain has no finite length".into()))?;
    // A distance a rounding step past the end is the end.
    let s = if s > total && s <= total + ARC_LENGTH_TOLERANCE * total.max(1.0) {
        total
    } else {
        s
    };
    let Some((index, _, local)) = chain.piece_at(s) else {
        return Err(invalid(format!(
            "distance {s} is outside the arc-length chain [0, {total}]"
        )));
    };
    let mut placement = Placement::IDENTITY;
    for (before, piece) in chain.pieces[..index].iter().enumerate() {
        let (end, tangent) = piece_local(before, piece, piece.length())?;
        placement = Placement {
            origin: placement.place(end),
            x: unit(placement.turn(tangent), "chain join")?,
        };
    }
    Ok((index, local, placement))
}

fn to_world(frame: &Frame2, v: Vec2) -> Vec2 {
    Vec2::new(
        frame.x.x * v.x + frame.y.x * v.y,
        frame.x.y * v.x + frame.y.y * v.y,
    )
}

/// Point on a chain at arc length `s` from its start. See the
/// [module documentation](self) for the accuracy contract and refusals.
pub fn chain_point(chain: &Chain2, s: Scalar) -> GeomResult<Point2> {
    let (index, local, placement) = locate(chain, s)?;
    let (p, _) = piece_local(index, &chain.pieces[index], local)?;
    let v = to_world(&chain.start, placement.place(p));
    Ok(Point2::new(
        chain.start.origin.x + v.x,
        chain.start.origin.y + v.y,
    ))
}

/// Unit tangent of a chain at arc length `s` (unit when the start frame is
/// orthonormal; a skewed frame maps it as it maps an intrinsic curve's).
pub fn chain_tangent(chain: &Chain2, s: Scalar) -> GeomResult<Vec2> {
    let (index, local, placement) = locate(chain, s)?;
    let (_, t) = piece_local(index, &chain.pieces[index], local)?;
    Ok(to_world(&chain.start, placement.turn(t)))
}

/// The parameter a parametric piece is read at, for a chain span inside
/// that piece; `None` for an intrinsic piece.
fn piece_parameter(piece: &ChainPiece2, s: Scalar) -> Option<GeomResult<Scalar>> {
    match piece {
        ChainPiece2::Parametric {
            curve,
            start,
            length,
        } if s == *length => Some(parameter_at_arc_length2(curve, *start, *length)),
        ChainPiece2::Parametric { curve, start, .. } => {
            Some(parameter_at_arc_length2(curve, *start, s))
        }
        _ => None,
    }
}

/// Certified chord bound of a chain over arc lengths `[a, b]` inside one
/// piece (#232); `None` across a join, for a span outside the chain, a
/// piece whose curve has no certified bound, or a chain that does not
/// evaluate.
///
/// A chain is parameterised by arc length, so `|c''|` is its curvature.
/// An intrinsic piece bounds that from its law; a parametric piece is the
/// same point set as its curve between the two parameters the span reads,
/// so its chord bound is the curve's. Either is then mapped through the
/// rigid placement (a rotation) and the chain's start frame. The bound is
/// about the chain as evaluated: vertices carry the accuracy contract of
/// the [module documentation](self).
#[must_use]
pub fn chain_chord_bound(chain: &Chain2, a: Scalar, b: Scalar) -> Option<Scalar> {
    let (lo, hi) = (a.min(b), a.max(b));
    let (index, local_lo, _) = locate(chain, lo).ok()?;
    let piece = &chain.pieces[index];
    // The piece start summed exactly as `Chain2::piece_at` sums it, so both
    // ends are read at the parameters evaluation reads them at.
    let begin = chain.pieces[..index]
        .iter()
        .fold(0.0, |begin, piece| begin + piece.length());
    let length = piece.length();
    let local_hi = hi - begin;
    if local_hi > length + ARC_LENGTH_TOLERANCE * length.max(1.0) {
        return None;
    }
    // At (or a rounding step from) the piece end, read the end itself: the
    // next piece starts exactly there.
    let local_hi = if length - local_hi <= 4.0 * Scalar::EPSILON * length.max(1.0) {
        length
    } else {
        local_hi
    };
    // Validate the piece itself (its frame, its available length).
    piece_local(index, piece, local_hi).ok()?;
    let stretch = crate::bound::frame_stretch2(&chain.start);
    let local = match piece {
        ChainPiece2::Intrinsic { curvature, .. } => {
            let k = curvature_bound(curvature, local_lo, local_hi)?;
            (local_hi - local_lo).powi(2) * 0.125 * k
        }
        ChainPiece2::Parametric { curve, .. } => {
            let u0 = piece_parameter(piece, local_lo)?.ok()?;
            let u1 = piece_parameter(piece, local_hi)?.ok()?;
            crate::bound::chord_bound2(curve, u0, u1)?
        }
        _ => return None,
    };
    let bound = local * stretch * (1.0 + crate::bound::ROUNDING);
    bound.is_finite().then_some(bound)
}

/// Whether [`crate::curve::flatten2`] certifies a chain's chords: every
/// intrinsic piece's curvature is bounded and every parametric piece's curve
/// is a certified family with no corner anywhere in its domain (#232).
///
/// Conservative: a corner outside the stretch a piece reads still answers
/// `false`, since locating the stretch needs the arc-length inverse.
#[must_use]
pub fn chain_certifies(chain: &Chain2) -> bool {
    chain.is_well_formed()
        && chain.pieces.iter().all(|piece| match piece {
            ChainPiece2::Intrinsic { curvature, length } => {
                curvature_bound(curvature, 0.0, *length).is_some()
            }
            ChainPiece2::Parametric { curve, .. } => {
                crate::bound::certifies_flattening2(curve)
                    && crate::bound::continuity_breaks2(curve, 1).is_empty()
            }
            _ => false,
        })
}

/// Upper bound on `|k(s)|` over `[lo, hi]` (`0 <= lo <= hi`), from the
/// triangle inequality term by term; `None` for a law it cannot bound.
fn curvature_bound(law: &CurvatureLaw, lo: Scalar, hi: Scalar) -> Option<Scalar> {
    let reach = lo.abs().max(hi.abs());
    let polynomial = |coefficients: &[Scalar]| -> Scalar {
        coefficients
            .iter()
            .enumerate()
            .map(|(power, c)| c.abs() * reach.powi(power as i32))
            .sum()
    };
    let bound = match law {
        CurvatureLaw::Constant { curvature } => curvature.abs(),
        CurvatureLaw::Polynomial { coefficients } => polynomial(coefficients),
        CurvatureLaw::Sinusoid {
            mean, amplitude, ..
        } => mean.abs() + amplitude.abs(),
        CurvatureLaw::Composite {
            polynomial: coefficients,
            harmonics,
        } => polynomial(coefficients) + harmonics.iter().map(|h| h.amplitude.abs()).sum::<Scalar>(),
        CurvatureLaw::Piecewise { breaks, laws } => {
            if laws.len() != breaks.len() + 1 {
                return None;
            }
            // Each piece is written in its own arc length from its seam.
            let mut worst: Scalar = 0.0;
            let mut begin = 0.0;
            for (index, piece) in laws.iter().enumerate() {
                let end = breaks.get(index).copied().unwrap_or(Scalar::INFINITY);
                if end >= lo && begin <= hi {
                    let from = (lo - begin).max(0.0);
                    let to = (hi.min(end) - begin).max(from);
                    worst = worst.max(curvature_bound(piece, from, to)?);
                }
                begin = end;
            }
            worst
        }
        _ => return None,
    };
    let bound = bound * (1.0 + crate::bound::ROUNDING);
    bound.is_finite().then_some(bound)
}
