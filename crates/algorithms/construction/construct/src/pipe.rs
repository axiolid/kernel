//! Swept disks along chains of straight segments and circular arcs (#232).
//!
//! A pipe run, a rebar or a cable tray edge is a disk swept along legs
//! joined by bends: authoring formats write it as a disk swept along a
//! composite of lines and arcs, or along a polyline whose corners are
//! rounded to a fillet radius.
//! [`swept_disk_along_pieces`] meshes both so that every point of the
//! exact tube lies within the chord budget of the triangles.
//!
//! # The bound
//!
//! Each piece is bounded on its own, with the split [`crate::sweep`] uses
//! for one arc (#231): half the budget `c` chords the disk, so a ring of
//! `n` points inscribed in the circle of radius `r` lies within
//! `r (1 - cos(pi/n)) <= c/2` of it, and every point of the exact tube is
//! within `c/2` of the tube the ring polygon traces when it is carried
//! rigidly along the exact centreline by the piece's own frames (each
//! exact section circle and its polygon lie in the same plane). The
//! other half bounds that polygon tube against the triangles:
//!
//! - Along a segment the frames do not turn, so the polygon tube is a
//!   prism whose walls are the mesh's own planar rectangles: zero.
//! - Along an arc of radius `R` the frames are rotation-minimising, which
//!   along a circle means they turn with it: consecutive stations are
//!   rotations of one ring about the arc's axis, and the ring's plane
//!   contains that axis. The polygon tube between two stations is then the
//!   revolution [`crate::revolve`] bounds, and its triangles are within
//!   `(R + r)(1 - cos(phi/2))` of it for a step `phi`: the sagitta of the
//!   outer side of the bend, not of the centreline. Its wall quads are
//!   planar trapezoids, so their twist adds nothing; the span bound
//!   [`crate::sweep`] computes, directrix sagitta plus the section's reach
//!   plus [`crate::loft`]'s quad bound, is exactly that sum here. The step
//!   halves until every span fits `c/2`; more than 4096 steps on a bend
//!   is refused with `BudgetExceeded`.
//!
//! # Joints
//!
//! Consecutive pieces share one station, so the walls meet watertight and
//! wound one way: the incoming piece's end section. The outgoing piece's
//! own start frame is that frame turned by the least rotation that takes
//! the incoming tangent to the outgoing one, moved to its start point, and
//! its stations follow from it. Where the pieces are tangent continuous
//! and meet, the two frames coincide. Where they turn by `theta` or leave
//! a gap `g`, every vertex of the shared station lies at most
//! `delta <= g + 2 r sin(theta/2)` from the station the outgoing piece's
//! bound was proved for; moving a triangle's vertices by at most `delta`
//! moves each of its points by at most `delta` (same barycentric weights),
//! so the first span of the outgoing piece adds the measured `delta` to its
//! bound and is refined until the sum fits.
//!
//! A swept disk along a directrix that is not tangent continuous is
//! undefined (the section has no single plane at a corner), so a joint
//! turning by more than
//!
//! `theta_max = 2 asin(c / (8 r))` (pi once `c >= 8 r`)
//!
//! is a corner and is refused by name: at that turn the shared section
//! leans off the outgoing piece's own by a quarter of the budget, and any
//! larger turn would be a kink no chord budget can hide. A turn at or below
//! it is tangent continuous at the requested accuracy and costs at most
//! `c/4` of the directrix's half. A gap above `c/4` is refused too.
//!
//! # Fillets
//!
//! A fillet radius `f` rounds every corner between two straight segments
//! (the swept disk's `fillet_radius`): a corner turning by
//! `theta` is replaced by the arc of radius `f` tangent to both segments,
//! which cuts `f tan(theta/2)` off each. Segments shortened to nothing are
//! dropped (two fillets then meet tangentially). Refused by name: a fillet
//! that needs more of a segment than it has, a corner beside an arc (only
//! corners between straight segments are rounded), a disk radius at or
//! above the fillet radius, or any bend radius, whose tube folds through
//! the bend's axis on its inside (equality is a horn torus, which some
//! formats permit but which bounds no two-manifold solid), and a reversal.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_mesh::TriMesh;

use crate::loft::{self, Frame, Station};
use crate::revolve::rotate;
use crate::sweep;

/// One smooth piece of a swept disk's directrix.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum PathPiece {
    /// A straight segment from `start` to `end`.
    Segment {
        /// First point, in the direction of travel.
        start: Point3,
        /// Last point.
        end: Point3,
    },
    /// A circular arc: `start` turned by `angle` radians, right-handed,
    /// about the line through `centre` along the unit `axis`. `start -
    /// centre` must be perpendicular to `axis`; a negative `angle` turns
    /// the other way.
    Arc {
        /// Centre of the arc's circle.
        centre: Point3,
        /// Unit normal of the arc's plane.
        axis: Vec3,
        /// First point, in the direction of travel.
        start: Point3,
        /// Signed turn, at most a full turn either way.
        angle: Scalar,
    },
}

impl PathPiece {
    /// The first point, in the direction of travel.
    #[must_use]
    pub fn start_point(&self) -> Point3 {
        match *self {
            Self::Segment { start, .. } | Self::Arc { start, .. } => start,
        }
    }

    /// The last point, in the direction of travel.
    #[must_use]
    pub fn end_point(&self) -> Point3 {
        self.point(1.0)
    }

    /// The point a fraction `s` of the way along, by length.
    #[must_use]
    pub fn point(&self, s: Scalar) -> Point3 {
        match *self {
            Self::Segment { start, end } => start + (end - start) * s,
            Self::Arc {
                centre,
                axis,
                start,
                angle,
            } => rotate(start, centre, axis, angle * s),
        }
    }

    /// The unit tangent, in the direction of travel, at `point` on the
    /// piece; zero for a degenerate piece.
    fn tangent_at(&self, point: Point3) -> Vec3 {
        match *self {
            Self::Segment { start, end } => (end - start).normalize_or_zero(),
            Self::Arc {
                centre,
                axis,
                angle,
                ..
            } => (axis.cross(point - centre) * angle.signum()).normalize_or_zero(),
        }
    }

    /// The unit tangent at the start, in the direction of travel.
    #[must_use]
    pub fn start_tangent(&self) -> Vec3 {
        self.tangent_at(self.start_point())
    }

    /// The unit tangent at the end, in the direction of travel.
    #[must_use]
    pub fn end_tangent(&self) -> Vec3 {
        self.tangent_at(self.end_point())
    }

    /// The length along the piece.
    #[must_use]
    pub fn length(&self) -> Scalar {
        match *self {
            Self::Segment { start, end } => (end - start).length(),
            Self::Arc {
                centre,
                start,
                angle,
                ..
            } => (start - centre).length() * angle.abs(),
        }
    }

    /// The same piece walked the other way.
    #[must_use]
    pub fn reversed(&self) -> Self {
        match *self {
            Self::Segment { start, end } => Self::Segment {
                start: end,
                end: start,
            },
            Self::Arc {
                centre,
                axis,
                angle,
                ..
            } => Self::Arc {
                centre,
                axis,
                start: self.end_point(),
                angle: -angle,
            },
        }
    }

    /// The part between the fractions `from` and `to` of its length, in
    /// that order (`to < from` walks it backwards).
    #[must_use]
    pub fn between(&self, from: Scalar, to: Scalar) -> Self {
        match *self {
            Self::Segment { .. } => Self::Segment {
                start: self.point(from),
                end: self.point(to),
            },
            Self::Arc {
                centre,
                axis,
                angle,
                ..
            } => Self::Arc {
                centre,
                axis,
                start: self.point(from),
                angle: angle * (to - from),
            },
        }
    }

    /// Refuse a piece with no direction or a malformed arc.
    fn check(&self, radius: Scalar) -> GeomResult<()> {
        let finite = |p: Point3| p.is_finite();
        match *self {
            Self::Segment { start, end } => {
                if !(finite(start) && finite(end)) {
                    return Err(GeomError::InvalidInput(
                        "swept disk directrix segment must be finite".to_owned(),
                    ));
                }
                if start == end {
                    return Err(GeomError::Degenerate(
                        "swept disk directrix segment has no length".to_owned(),
                    ));
                }
            }
            Self::Arc {
                centre,
                axis,
                start,
                angle,
            } => {
                if !(finite(centre) && finite(axis) && finite(start) && angle.is_finite()) {
                    return Err(GeomError::InvalidInput(
                        "swept disk directrix arc must be finite".to_owned(),
                    ));
                }
                let bend = (start - centre).length();
                if (axis.length() - 1.0).abs() > 1e-9
                    || bend == 0.0
                    || (start - centre).dot(axis).abs() > 1e-9 * bend
                {
                    return Err(GeomError::InvalidInput(
                        "swept disk directrix arc needs a unit axis normal to a non-zero \
                         radius"
                            .to_owned(),
                    ));
                }
                if angle == 0.0 || angle.abs() > core::f64::consts::TAU * (1.0 + 1e-12) {
                    return Err(GeomError::InvalidInput(format!(
                        "swept disk directrix arc must turn by a non-zero angle of at most \
                         a full turn, got {angle}"
                    )));
                }
                if bend <= radius {
                    return Err(GeomError::InvalidInput(format!(
                        "swept disk radius {radius} is not below the bend radius {bend}: \
                         the tube would fold through the bend's axis"
                    )));
                }
            }
        }
        Ok(())
    }

    /// Samples to chord budget `budget`, and the exact unit tangent at
    /// each: a segment's two ends, or an arc in equal steps whose sagitta
    /// is within the budget.
    fn sample(&self, budget: Scalar) -> GeomResult<(Vec<Point3>, Vec<Vec3>)> {
        let n = match *self {
            Self::Segment { .. } => 1,
            Self::Arc {
                centre,
                start,
                angle,
                ..
            } => crate::revolve::steps((start - centre).length(), angle, budget).map_err(|_| {
                GeomError::BudgetExceeded {
                    resource: "swept disk bend steps",
                }
            })?,
        };
        let points: Vec<Point3> = (0..=n)
            .map(|k| self.point(k as Scalar / n as Scalar))
            .collect();
        let tangents = points.iter().map(|p| self.tangent_at(*p)).collect();
        Ok((points, tangents))
    }
}

/// The largest turn a joint may make and still be tangent continuous at
/// chord budget `chord` for a disk of radius `radius`: `2 asin(c / (8 r))`
/// (see the module notes).
#[must_use]
pub fn joint_tolerance(chord: Scalar, radius: Scalar) -> Scalar {
    2.0 * (chord / (8.0 * radius)).min(1.0).asin()
}

/// The angle between two unit vectors, accurate when it is small.
fn turn(a: Vec3, b: Vec3) -> Scalar {
    a.cross(b).length().atan2(a.dot(b))
}

/// Sweep a disk, optionally hollow, along a chain of segments and arcs so
/// that every point of the exact tube lies within `chord` of the mesh.
///
/// `fillet_radius` rounds each corner between two segments; see the
/// module notes for the bound, the joint rule, and what is refused. The
/// pieces run in the direction of travel, each starting where the last
/// ended.
pub fn swept_disk_along_pieces(
    pieces: &[PathPiece],
    radius: Scalar,
    inner_radius: Option<Scalar>,
    fillet_radius: Option<Scalar>,
    chord: Scalar,
) -> GeomResult<TriMesh> {
    sweep::check_radii(radius, inner_radius)?;
    if !(chord.is_finite() && chord > 0.0) {
        return Err(GeomError::InvalidInput(format!(
            "swept disk chord budget must be positive and finite, got {chord}"
        )));
    }
    if let Some(f) = fillet_radius {
        if !(f.is_finite() && f > 0.0) {
            return Err(GeomError::InvalidInput(format!(
                "swept disk fillet radius must be positive and finite, got {f}"
            )));
        }
    }
    if pieces.is_empty() {
        return Err(GeomError::InvalidInput(
            "a swept disk directrix needs at least one piece".to_owned(),
        ));
    }
    let tolerance = joint_tolerance(chord, radius);
    let path = match fillet_radius {
        Some(f) => fillet(pieces, f, radius, tolerance)?,
        None => pieces.to_vec(),
    };
    for piece in &path {
        piece.check(radius)?;
    }
    for (k, pair) in path.windows(2).enumerate() {
        let theta = turn(pair[0].end_tangent(), pair[1].start_tangent());
        if theta > tolerance {
            let advice = if fillet_radius.is_some() {
                "only corners between two straight segments are filleted"
            } else {
                "a swept disk along it is undefined; give a fillet radius or a \
                 tangent-continuous directrix"
            };
            return Err(GeomError::InvalidInput(format!(
                "swept disk directrix turns a corner of {theta} rad at joint {k}, above \
                 the {tolerance} rad a tangent-continuous joint may turn at this chord \
                 budget: {advice}"
            )));
        }
        let gap = (pair[1].start_point() - pair[0].end_point()).length();
        if gap > 0.25 * chord {
            return Err(GeomError::InvalidInput(format!(
                "swept disk directrix pieces leave a {gap} gap at joint {k}, more than a \
                 quarter of the {chord} chord budget"
            )));
        }
    }
    let rings = sweep::disk_rings(radius, inner_radius, 0.5 * chord)?;
    let share = 0.5 * chord;
    let mut stations: Vec<Station> = Vec::new();
    let mut last: Option<Frame> = None;
    for piece in &path {
        let tangent = piece.start_tangent();
        // The outgoing piece's own start frame: the incoming end frame
        // turned by the least rotation onto this tangent.
        let seed = match &last {
            None => seed_reference(tangent),
            Some(frame) => least_rotation(frame.x, frame.x.cross(frame.y), tangent),
        };
        let joint = stations.last();
        let (frames, placed) = bounded_piece(piece, &rings, seed, radius, share, joint)?;
        let skip = usize::from(joint.is_some());
        stations.extend(placed.into_iter().skip(skip));
        last = frames.into_iter().last();
    }
    let mut mesh = loft::loft(&rings, &stations, false)?;
    // As `revolve::revolve` (#221): settle the winding from the built
    // mesh rather than from the frames' handedness.
    if matches!(crate::extrude::outward_orientation(&mesh), Some(false)) {
        for triangle in mesh.indices.chunks_exact_mut(3) {
            triangle.swap(1, 2);
        }
    }
    Ok(mesh)
}

/// Rounds of halving a piece's directrix budget before it is refused.
const MAX_REFINEMENTS: usize = 16;

/// A piece's frames and stations, its directrix sampled until each span's
/// bound fits `share`; the first span also carries the measured distance
/// from the shared `joint` station to the piece's own first station.
fn bounded_piece(
    piece: &PathPiece,
    rings: &crate::profile::Rings,
    seed: Vec3,
    reach: Scalar,
    share: Scalar,
    joint: Option<&Station>,
) -> GeomResult<(Vec<Frame>, Vec<Station>)> {
    let mut b = share;
    for _ in 0..MAX_REFINEMENTS {
        let (points, tangents) = piece.sample(b)?;
        let frames = sweep::frames_carried(&points, &tangents, seed)?;
        let placed: Vec<Station> = frames
            .iter()
            .map(|f| loft::place(rings, |p| loft::at(f, p)))
            .collect();
        let shift = joint.map_or(0.0, |j| displacement(j, &placed[0]));
        let worst = frames
            .windows(2)
            .zip(placed.windows(2))
            .enumerate()
            .map(|(k, (f, s))| {
                let moved = if k == 0 { shift } else { 0.0 };
                sweep::span_bound(&f[0], &f[1], b, reach)
                    + loft::span_deviation(&s[0], &s[1])
                    + moved
            })
            .fold(0.0, Scalar::max);
        if worst <= share {
            return Ok((frames, placed));
        }
        b *= 0.5;
    }
    Err(GeomError::BudgetExceeded {
        resource: "swept disk directrix refinement",
    })
}

/// The largest distance between corresponding vertices of two stations.
fn displacement(a: &Station, b: &Station) -> Scalar {
    a.loops
        .iter()
        .zip(&b.loops)
        .flat_map(|(p, q)| p.iter().zip(q).map(|(u, v)| (*u - *v).length()))
        .fold(0.0, Scalar::max)
}

/// `v` turned by the least rotation that takes unit `from` to unit `to`.
fn least_rotation(v: Vec3, from: Vec3, to: Vec3) -> Vec3 {
    let axis = from.cross(to);
    let sin = axis.length();
    if sin == 0.0 {
        return v;
    }
    rotate(v, Point3::ZERO, axis / sin, turn(from, to))
}

/// Any direction perpendicular to unit `t`.
fn seed_reference(t: Vec3) -> Vec3 {
    let candidate = if t.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    candidate - t * t.dot(candidate)
}

/// Round every corner between two segments that turns by more than
/// `tolerance` with an arc of radius `fillet` (see the module notes).
fn fillet(
    pieces: &[PathPiece],
    fillet: Scalar,
    radius: Scalar,
    tolerance: Scalar,
) -> GeomResult<Vec<PathPiece>> {
    // Each segment keeps its direction even when the fillets either side
    // use it up, so the corner after it is still measured against it.
    struct Leg {
        piece: PathPiece,
        dir: Vec3,
        trim_start: Scalar,
        trim_end: Scalar,
    }
    let mut legs: Vec<Leg> = Vec::with_capacity(pieces.len());
    let mut arcs: Vec<Option<PathPiece>> = Vec::with_capacity(pieces.len());
    for piece in pieces {
        piece.check(radius)?;
        legs.push(Leg {
            piece: *piece,
            dir: piece.start_tangent(),
            trim_start: 0.0,
            trim_end: 0.0,
        });
    }
    for k in 1..legs.len() {
        let (PathPiece::Segment { .. }, PathPiece::Segment { start: corner, .. }) =
            (legs[k - 1].piece, legs[k].piece)
        else {
            arcs.push(None);
            continue;
        };
        let (u, v) = (legs[k - 1].dir, legs[k].dir);
        let theta = turn(u, v);
        if theta <= tolerance {
            arcs.push(None);
            continue;
        }
        if radius >= fillet {
            return Err(GeomError::InvalidInput(format!(
                "swept disk radius {radius} is not below the fillet radius {fillet}: the \
                 tube would fold through the bend's axis"
            )));
        }
        let normal = u.cross(v);
        if normal.length() <= 1e-12 {
            return Err(GeomError::InvalidInput(format!(
                "swept disk directrix reverses at corner {k}; no fillet rounds it"
            )));
        }
        let cut = fillet * (0.5 * theta).tan();
        legs[k - 1].trim_end = cut;
        legs[k].trim_start = cut;
        let tangent_point = corner - u * cut;
        let inward = (v - u * u.dot(v)).normalize();
        arcs.push(Some(PathPiece::Arc {
            centre: tangent_point + inward * fillet,
            axis: normal.normalize(),
            start: tangent_point,
            angle: theta,
        }));
    }
    let mut out = Vec::with_capacity(2 * legs.len());
    for (k, leg) in legs.iter().enumerate() {
        if k > 0 {
            if let Some(arc) = arcs[k - 1] {
                out.push(arc);
            }
        }
        let length = leg.piece.length();
        let needed = leg.trim_start + leg.trim_end;
        if needed == 0.0 {
            out.push(leg.piece);
            continue;
        }
        // A fillet using up exactly its segment meets the next one
        // tangentially; rounding in `tan` may overshoot by an ulp or two.
        let slack = 1e-9 * length;
        if needed > length + slack {
            return Err(GeomError::InvalidInput(format!(
                "swept disk fillet radius {fillet} does not fit: segment {k} is {length} \
                 long but its fillets need {needed} of it"
            )));
        }
        if needed < length - slack {
            out.push(
                leg.piece
                    .between(leg.trim_start / length, 1.0 - leg.trim_end / length),
            );
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_fillet_is_tangent_to_both_segments() {
        let corner = Point3::new(1.0, 0.0, 0.0);
        let pieces = [
            PathPiece::Segment {
                start: Point3::ZERO,
                end: corner,
            },
            PathPiece::Segment {
                start: corner,
                end: Point3::new(1.0, 1.0, 1.0),
            },
        ];
        let path = fillet(&pieces, 0.2, 0.05, 1e-6).unwrap();
        assert_eq!(path.len(), 3);
        for pair in path.windows(2) {
            assert!((pair[0].end_point() - pair[1].start_point()).length() < 1e-12);
            assert!(turn(pair[0].end_tangent(), pair[1].start_tangent()) < 1e-12);
        }
    }

    #[test]
    fn an_arc_walked_backwards_retraces_itself() {
        let arc = PathPiece::Arc {
            centre: Point3::new(0.0, 1.0, 0.0),
            axis: Vec3::Z,
            start: Point3::ZERO,
            angle: 1.0,
        };
        let back = arc.reversed();
        assert!((back.end_point() - arc.start_point()).length() < 1e-12);
        assert!((back.start_tangent() + arc.end_tangent()).length() < 1e-12);
        let part = arc.between(0.25, 0.75);
        assert!((part.length() - 0.5).abs() < 1e-12);
        assert!((part.start_point() - arc.point(0.25)).length() < 1e-12);
    }
}
