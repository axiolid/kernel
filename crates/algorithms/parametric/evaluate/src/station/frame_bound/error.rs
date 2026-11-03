//! How far one station's frame lies from the exact frame law (#405); the
//! derivation is the parent module's section "One station".

use axiolid_core::Scalar;
use axiolid_curve::SeamSide;

use super::{FrameLaw, FrameReading, LawBasis, Slot};
use crate::arc_parameter::ARC_LENGTH_TOLERANCE;
use crate::bound::ROUNDING;
use crate::station::certify::{AXIS_ROUNDING, POINT_ROUNDING};
use crate::station::offset::OffsetTable;
use crate::station::{admitted, StationCurve, StationPiece};

/// How far a frame read at a station lies from the exact frame law's, see
/// the [module documentation](super#one-station-405).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StationFrameError {
    /// The largest distance of the read frame's point from the exact
    /// frame's.
    pub point: Scalar,
    /// The largest distance of one of the read frame's axes from the exact
    /// frame's.
    pub axis: Scalar,
}

/// The arc-length tolerance around a distance: what counts as ON a seam.
fn slack(distance: Scalar) -> Scalar {
    ARC_LENGTH_TOLERANCE * distance.abs().max(1.0)
}

/// A station read: its distance along the basis and into the piece it is
/// read on, and how it is read.
#[derive(Clone, Copy)]
struct Stand {
    at: Scalar,
    into: Scalar,
    side: SeamSide,
    reading: FrameReading,
}

/// The largest magnitude the closed-form read of `piece` works with: its
/// curve's stored data (a circle's `|origin|_inf + r`, a polyline's points,
/// a line's origin and direction, a plan's data under an elevation) and the
/// placement's translation. The read's cancellations (`r - r cos a`) round
/// at this scale, not at the point's.
fn data_scale(piece: &StationPiece<'_>) -> Scalar {
    use axiolid_curve::{Curve2, Curve3};
    fn plan(curve: &Curve2) -> Scalar {
        match curve {
            Curve2::Line(l) => l.origin.abs().max_element() + l.direction.abs().max_element(),
            Curve2::Circle(c) => c.frame.origin.abs().max_element() + c.radius.abs(),
            Curve2::Polyline(p) => p
                .points
                .iter()
                .map(|q| q.abs().max_element())
                .fold(0.0, Scalar::max),
            _ => Scalar::INFINITY,
        }
    }
    let curve = match piece.curve {
        StationCurve::Two(c) => plan(c),
        StationCurve::Three(c) => match c {
            Curve3::Line(l) => l.origin.abs().max_element() + l.direction.abs().max_element(),
            Curve3::Circle(c) => c.frame.origin.abs().max_element() + c.radius.abs(),
            Curve3::Polyline(p) => p
                .points
                .iter()
                .map(|q| q.abs().max_element())
                .fold(0.0, Scalar::max),
            Curve3::Elevated(e) => plan(&e.plan),
            Curve3::Banked(b) => plan(&b.base.plan),
            _ => Scalar::INFINITY,
        },
        StationCurve::Offset(_) => Scalar::INFINITY,
    };
    let moved = piece
        .placement
        .map_or(0.0, |rigid| rigid.translation.abs().max_element());
    curve.max(moved)
}

impl<'c> FrameLaw<'_, 'c> {
    /// How far the frame the kernel reads at `distance` along the basis,
    /// from `side` of a seam or joint the station is on, lies from the
    /// exact frame law's there, read as `reading`; `None` where it is not
    /// certified (see [one station](super#one-station-405)).
    #[must_use]
    pub fn error_at(
        &self,
        distance: Scalar,
        side: SeamSide,
        reading: FrameReading,
    ) -> Option<StationFrameError> {
        match &self.basis {
            LawBasis::Atomic {
                curve,
                placement,
                whole,
            } => {
                let owned;
                let (piece, length) = match whole {
                    Some(whole) => (whole, Some(whole.length())),
                    None => {
                        // A curve with no bounded length (a line plan) is
                        // read on a span about the station.
                        let (mut piece, length) = match StationPiece::whole(*curve) {
                            Ok(whole) => (whole, Some(whole.length())),
                            Err(_) => {
                                let lo = (distance - 1.0).max(0.0);
                                (
                                    StationPiece::between(*curve, lo, distance + 1.0).ok()?,
                                    None,
                                )
                            }
                        };
                        if let Some(rigid) = placement {
                            piece = piece.placed(*rigid, true);
                        }
                        owned = piece;
                        (&owned, length)
                    }
                };
                let at = admitted(distance, length).ok()?;
                let stand = Stand {
                    at,
                    into: at - piece.start,
                    side,
                    reading,
                };
                self.error_in(0, piece, None, stand, (0.0, [false; 2]))
            }
            LawBasis::Composite(composite) => {
                let at = admitted(distance, Some(composite.length())).ok()?;
                let (k, into) = composite.locate_admitted(distance, side).ok()?;
                let piece = &composite.pieces()[k];
                if !piece.placement_exact {
                    return None;
                }
                let joint = self.joint_error(composite, k)?;
                let ends = [k > 0, k + 1 < composite.pieces().len()];
                let table = composite.tables()[k].as_ref();
                let stand = Stand {
                    at,
                    into,
                    side,
                    reading,
                };
                self.error_in(k, piece, table, stand, (joint, ends))
            }
        }
    }

    /// The error of the frame read at `stand` on piece `k`; `joint` the
    /// distance of the piece's start from the true joint and `ends`
    /// whether its start and end are joints of a composite.
    fn error_in(
        &self,
        k: usize,
        piece: &StationPiece<'c>,
        table: Option<&OffsetTable>,
        stand: Stand,
        (joint, ends): (Scalar, [bool; 2]),
    ) -> Option<StationFrameError> {
        let Stand {
            at, side, reading, ..
        } = stand;
        let length = piece.length();
        let u = stand.into.clamp(0.0, length);
        let here = self.read(k, piece, u, u, table, reading)?;
        let (e_p, e_a) = match self.slot(k, piece, table) {
            Slot::Closed => {
                let point = piece.section_with(u, side, table).ok()?.point;
                // `POINT_ROUNDING` of the largest magnitude in the read:
                // the point, the curve's data and the placement.
                let scale = point.abs().max_element().max(data_scale(piece)).max(1.0);
                (POINT_ROUNDING * scale, AXIS_ROUNDING)
            }
            Slot::Intrinsic(_) => (here.point_error, here.axis_error),
            Slot::Read(_) | Slot::Unbounded => return None,
        };
        let delta = here.measure + joint;
        let (mut s_point, mut s_axis) = (0.0, 0.0);
        if delta > 0.0 {
            let (mut lo, mut hi) = (u - 2.0 * delta, u + 2.0 * delta);
            let (measure, _) = piece.read_at(u, side);
            let on = slack(at).min(slack(measure));
            let mut seams: Vec<Scalar> = piece.seams().ok()?.iter().map(|s| s.distance).collect();
            if ends[0] {
                seams.push(0.0);
            }
            if ends[1] {
                seams.push(length);
            }
            for c in seams {
                if (c - u).abs() <= on {
                    // On the seam: the side named.
                    match side {
                        SeamSide::Incoming => hi = hi.min(c),
                        _ => lo = lo.max(c),
                    }
                } else if (c - u).abs() <= 2.0 * delta {
                    return None;
                }
            }
            let span = self.read(k, piece, lo.max(0.0), hi.min(length), table, reading)?;
            s_point = span.point[0];
            s_axis = span.tangent[0].max(span.lateral[0]).max(span.up[0]);
        }
        let grow = 1.0 + ROUNDING;
        let out = StationFrameError {
            point: (e_p + delta * s_point) * grow,
            axis: (e_a + delta * s_axis) * grow,
        };
        (out.point.is_finite() && out.axis.is_finite()).then_some(out)
    }
}
