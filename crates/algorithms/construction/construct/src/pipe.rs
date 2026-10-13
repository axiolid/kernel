//! Swept disks along chains of straight segments and circular arcs (#232,
//! #245, #248).
//!
//! A pipe run, a rebar or a cable tray edge is a disk swept along legs
//! joined by bends or corners: authoring formats write it as a disk swept
//! along a composite of lines and arcs, or along a polyline whose corners
//! are mitred or rounded to a fillet radius.
//! [`swept_disk_along_pieces`] meshes all of them so that every point of
//! the exact tube lies within the chord budget of the triangles.
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
//! A joint turning by at most
//!
//! `theta_max = 2 asin(c / (8 r))` (pi once `c >= 8 r`)
//!
//! is tangent continuous at the requested accuracy: its shared section
//! leans off the outgoing piece's own by at most a quarter of the budget,
//! which that piece's first span pays as above. A gap above `c/4` is
//! refused.
//!
//! # Mitres
//!
//! A sharper joint between two straight segments is a corner, and a disk
//! swept round a corner is mitred at half angle (#245): both tubes are cut
//! by the plane through the corner whose normal `n = (u + v) / |u + v|`
//! bisects the incoming and outgoing tangents `u` and `v`. Reflection in
//! that plane maps `u` to `-v`, so it swaps the two axis lines and with
//! them the two cylinders of radius `r` (and the two bores); it fixes the
//! plane, so both cylinders cut it in the same ellipse. The outgoing frame
//! is the incoming one turned by the least rotation about `u x v`, which
//! agrees with that reflection on the plane square to `u`. So profile point
//! `k` of the incoming end section and of the outgoing start section, each
//! projected along its own axis onto the mitre plane, land on one point:
//! the two pieces share that mitre ring, and every vertex of it lies on
//! both exact cylinders. Each wall quad of a segment joins two generator
//! lines of the prism over the ring polygon, so it is a planar trapezoid,
//! exactly a face of that prism cut by the end planes: the walls add
//! nothing, and the span bound measures them as zero.
//!
//! The section pays instead. A point `w + s u` of the exact cylinder (`w`
//! square to `u`, `|w| = r`; the mitre plane cuts its generator at
//! `s_cut(w) = -w.n / u.n`, at most `r tan(theta/2)`) lies within
//! `(1 - l) r <= e = r (1 - cos(pi/m))` of the prism point `w' + s u` at
//! its azimuth, `w' = l w` on the ring polygon of `m` points. When that
//! point lies past the plane, the plane crosses the prism's generator at
//! `l s_cut(w)`, which moves the point by at most
//! `(1 - l) s_cut(w) <= (1 - l) r tan(theta/2)` along `u`, square to the
//! first move: the distance is at most `e / cos(theta/2)`, attained where
//! the ellipse's major axis ends. The rings are therefore chorded to
//! `(c/2) cos(theta_m/2)` for the sharpest mitre `theta_m` (only finer
//! elsewhere), and a mitred pipe keeps every point of the exact tube
//! within `c` of its triangles, as a smooth one does. A hollow disk's bore
//! is mitred by the same plane and bounded the same way.
//!
//! Along a segment of length `L` the tube runs, on the generator at `w`,
//! from its start plane to its end plane; both cuts are linear in `w`, so
//! the shortest generator is `L - r |g_perp|` with
//! `g = n1 / u.n1 - n0 / u.n0` (`n0 = u` at a square end), and it must be
//! positive (one mitre: `L > r tan(theta/2)`). Refused by name: a segment
//! whose mitres meet inside the tube, which would cut through itself, and a
//! reversal (`theta = pi`, where the mitre plane holds both legs).
//!
//! # Mitres beside an arc
//!
//! A corner beside an arc is mitred at half angle too (#248), by the same
//! bisector plane, but the plane cuts the bend's torus in a different curve
//! than the tube on its other side: no ring lies on both exact surfaces,
//! and the exact solid (each tube continued to the plane and cut by it) has
//! a planar face between the two section curves. The end of the bend is
//! meshed as a [`MitreBand`] (see its module): from a square section
//! `back` short of the corner, each profile vertex walks along its own
//! circle of the torus to the plane, so the band's vertices lie on the
//! exact torus, and where the shared ring is not the band's own cut, the
//! band stops just short of the plane and a strip, nearly in the plane,
//! joins it to that ring: the mesh of the planar face. The shared ring is
//! the outgoing piece's: a segment's mitre ring as above (so a segment
//! beside an arc keeps the proof of the previous section, its ring on its
//! own cylinder), or, between two arcs, the outgoing band's cut. Nothing
//! proves the band and the strip within the budget by construction; the
//! compiler certifies them against the mesh and reports that bound beside
//! the proven one ([`pipe_plan`] gives it the bands). Refused by name: a
//! corner whose plane does not clear the bend's tube within a quarter turn
//! of it, and bands that need more of an arc than it turns.
//!
//! # Fillets
//!
//! A fillet radius `f` rounds every corner between two straight segments
//! (the swept disk's `fillet_radius`): a corner turning by
//! `theta` is replaced by the arc of radius `f` tangent to both segments,
//! which cuts `f tan(theta/2)` off each. Segments shortened to nothing are
//! dropped (two fillets then meet tangentially). Refused by name: a fillet
//! that needs more of a segment than it has, a corner beside an arc (only
//! corners between straight segments are rounded), a disk radius above the
//! fillet radius, or any bend radius, whose tube folds through the bend's
//! axis on its inside, and a reversal.
//!
//! A disk radius equal to the fillet (or bend) radius stays refused by name
//! (#248), although the format rule permits it (fillet radius at least the
//! disk radius): the bend is then a horn torus whose inner wall pinches to
//! one point on the bend's axis, where every section touches. That boundary
//! is not a two-manifold; a closed mesh within the budget would have to
//! either meet itself there (a non-manifold vertex, degenerate triangles
//! along the pinch) or stand off the pinch and misstate it.
//!
//! # Closed directrices
//!
//! [`swept_disk_along_loop`] sweeps a loop: the last piece ends where the
//! first starts, and that joint is mitred, filleted or tangent like any
//! other (#248). The walls close on the first station, so there are no
//! caps. A frame carried round a loop that leaves its plane does not come
//! back to itself: least rotations at the joints and rotation-minimising
//! frames along the bends return it turned about the first tangent by the
//! loop's holonomy `tau`, in `(-pi, pi]` ([`PipePlan::twist`], the
//! measure reported). The pipe turns its frames back by `tau`, distributed
//! uniformly over the loop's clear length: every segment's length outside
//! the mitre zones at its ends (`r tan(theta/2)` from each mitred corner,
//! where the plane reaches along it), and every bend's outside its mitre
//! bands. In the zones and bands the frame is held, so a mitre still shares
//! one ring on both cylinders and its walls stay planar. The twist is part
//! of the proven bound: [`crate::sweep`]'s span bound measures the turn
//! between two frames, twist included, and a twisted segment gets stations
//! until each span's twist fits (more than 4096 is refused with
//! `BudgetExceeded`). The station closing the loop is the first one; the
//! last span pays its measured distance from the station it was proved for.
//! Refused by name: a loop whose frame comes back turned when no part of it
//! is clear to take the twist.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_mesh::TriMesh;

use crate::loft::{self, Frame, Station};
use crate::profile::Rings;
use crate::revolve::rotate;
use crate::sweep;

mod band;
pub use band::MitreBand;

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
                if bend < radius {
                    return Err(GeomError::InvalidInput(format!(
                        "swept disk radius {radius} is not below the bend radius {bend}: \
                         the tube would fold through the bend's axis"
                    )));
                }
                if bend == radius {
                    return Err(horn_torus("bend", radius));
                }
            }
        }
        Ok(())
    }

    /// Equal steps to chord budget `budget`: one for a segment, for an arc
    /// as many as keep each step's sagitta within the budget.
    fn steps(&self, budget: Scalar) -> GeomResult<usize> {
        match *self {
            Self::Segment { .. } => Ok(1),
            Self::Arc {
                centre,
                start,
                angle,
                ..
            } => crate::revolve::steps((start - centre).length(), angle, budget).map_err(|_| {
                GeomError::BudgetExceeded {
                    resource: "swept disk bend steps",
                }
            }),
        }
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

/// One side of a corner mitred beside an arc (#248).
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum MitreSide {
    /// A straight segment: near the corner its tube is the cylinder of the
    /// disk radius about the line through `corner` along `direction`.
    Segment {
        /// The segment's end at the corner.
        corner: Point3,
        /// Its unit direction of travel.
        direction: Vec3,
    },
    /// An arc: near the corner its tube is the band.
    Arc(MitreBand),
}

/// A corner mitred beside an arc (#248): the mitre plane and the tubes on
/// its two sides.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct ArcMitre {
    /// The joint, between piece `joint` of [`PipePlan::pieces`] and the
    /// next one (the first, for the joint that closes a loop).
    pub joint: usize,
    /// Unit normal of the mitre plane: the bisector of the incoming and the
    /// outgoing tangent.
    pub normal: Vec3,
    /// The incoming side, then the outgoing one.
    pub sides: [MitreSide; 2],
}

/// How [`swept_disk_along_pieces`] and [`swept_disk_along_loop`] sweep a
/// directrix (#248): what a caller needs to bound the parts of the mesh no
/// construction proof covers.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct PipePlan {
    /// The pieces swept, fillets included.
    pub pieces: Vec<PathPiece>,
    /// Whether the last piece joins the first.
    pub closed: bool,
    /// Every corner mitred beside an arc, in joint order.
    pub arc_mitres: Vec<ArcMitre>,
    /// For a loop, the angle in `(-pi, pi]` about the first piece's start
    /// tangent by which a frame carried once round the loop comes back
    /// turned (its holonomy); the pipe turns its frames back by it,
    /// uniformly over the loop's clear length (see the module notes). Zero
    /// for an open directrix.
    pub twist: Scalar,
}

/// Sweep a disk, optionally hollow, along a chain of segments and arcs so
/// that every point of the exact tube lies within `chord` of the mesh.
///
/// `fillet_radius` rounds each corner between two segments; without it,
/// each corner is mitred at half angle (#245, beside an arc #248). See the
/// module notes for the bound, the joint rule, the mitres, and what is
/// refused; a corner beside an arc is not proven but certified by the
/// caller ([`pipe_plan`]). The pieces run in the direction of travel, each
/// starting where the last ended.
pub fn swept_disk_along_pieces(
    pieces: &[PathPiece],
    radius: Scalar,
    inner_radius: Option<Scalar>,
    fillet_radius: Option<Scalar>,
    chord: Scalar,
) -> GeomResult<TriMesh> {
    Pipe::new(pieces, false, radius, inner_radius, fillet_radius, chord)?.mesh()
}

/// [`swept_disk_along_pieces`] round a loop (#248): the last piece ends
/// where the first starts, that joint is mitred, filleted or tangent like
/// any other, and the walls close on the first station without caps. A
/// frame that comes back turned is twisted back uniformly over the loop's
/// clear length (see the module notes and [`PipePlan::twist`]).
pub fn swept_disk_along_loop(
    pieces: &[PathPiece],
    radius: Scalar,
    inner_radius: Option<Scalar>,
    fillet_radius: Option<Scalar>,
    chord: Scalar,
) -> GeomResult<TriMesh> {
    Pipe::new(pieces, true, radius, inner_radius, fillet_radius, chord)?.mesh()
}

/// How [`swept_disk_along_pieces`] (`closed` false) or
/// [`swept_disk_along_loop`] (`closed` true) sweeps these pieces at this
/// budget, refused exactly as they refuse it.
pub fn pipe_plan(
    pieces: &[PathPiece],
    closed: bool,
    radius: Scalar,
    inner_radius: Option<Scalar>,
    fillet_radius: Option<Scalar>,
    chord: Scalar,
) -> GeomResult<PipePlan> {
    Ok(Pipe::new(pieces, closed, radius, inner_radius, fillet_radius, chord)?.plan())
}

/// Rounds of halving a piece's directrix budget before it is refused.
const MAX_REFINEMENTS: usize = 16;

/// Most stations the twist of one piece may need.
const MAX_TWIST_STEPS: usize = 4096;

/// A closing twist whose turn at the disk's rim is below this fraction of
/// the chord budget is rounding: the last span pays it as a measured
/// shift instead of the loop being twisted.
const ROUNDING_TWIST: Scalar = 1e-3;

/// How far short of the mitre plane a band stops before its strip, as a
/// fraction of the chord budget.
const STRIP_OFFSET: Scalar = 1.0 / 16.0;

/// A directrix checked and planned for sweeping.
struct Pipe {
    path: Vec<PathPiece>,
    closed: bool,
    radius: Scalar,
    inner: Option<Scalar>,
    chord: Scalar,
    /// The mitre plane's unit normal at each joint, `None` where the joint
    /// is tangent continuous at this budget. Joint `k` is between pieces
    /// `k` and `k + 1` (the first, closing a loop).
    mitres: Vec<Option<Vec3>>,
    /// The sharpest mitre's turn.
    sharpest: Scalar,
    /// Each piece's mitre bands at its start and its end.
    bands: Vec<[Option<MitreBand>; 2]>,
    /// The closing twist, and the rate it is turned back at per unit of
    /// clear length (zero where it is rounding).
    twist: Scalar,
    rate: Scalar,
}

/// A piece's stations, the untwisted frame at its end, and the twist
/// reached there.
struct Built {
    stations: Vec<Station>,
    end: Frame,
    turned: Scalar,
}

/// Where a piece's last station must land.
#[derive(Clone, Copy)]
enum Close<'a> {
    /// Anywhere: an open end, or a joint the next piece pays for.
    Open,
    /// On this station, the first of the loop.
    To(&'a Station),
    /// On the piece's own first station: a loop of one piece.
    Own,
}

/// A piece sampled for a budget: points, the exact unit tangent at each,
/// and the twist about it.
struct Sampled {
    points: Vec<Point3>,
    tangents: Vec<Vec3>,
    twists: Vec<Scalar>,
}

impl Pipe {
    fn new(
        pieces: &[PathPiece],
        closed: bool,
        radius: Scalar,
        inner: Option<Scalar>,
        fillet_radius: Option<Scalar>,
        chord: Scalar,
    ) -> GeomResult<Self> {
        sweep::check_radii(radius, inner)?;
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
            Some(f) => fillet(pieces, f, radius, tolerance, closed)?,
            None => pieces.to_vec(),
        };
        for piece in &path {
            piece.check(radius)?;
        }
        let m = path.len();
        let joints = if closed { m } else { m - 1 };
        let mut mitres: Vec<Option<Vec3>> = Vec::with_capacity(joints);
        let mut sharpest: Scalar = 0.0;
        for k in 0..joints {
            let (a, b) = (path[k], path[(k + 1) % m]);
            let gap = (b.start_point() - a.end_point()).length();
            if gap > 0.25 * chord {
                return Err(GeomError::InvalidInput(format!(
                    "swept disk directrix pieces leave a {gap} gap at joint {k}, more than a \
                     quarter of the {chord} chord budget"
                )));
            }
            let (u, v) = (a.end_tangent(), b.start_tangent());
            let theta = turn(u, v);
            if theta <= tolerance {
                mitres.push(None);
                continue;
            }
            let straight = matches!(
                (a, b),
                (PathPiece::Segment { .. }, PathPiece::Segment { .. })
            );
            if fillet_radius.is_some() && !straight {
                return Err(GeomError::InvalidInput(format!(
                    "swept disk directrix turns a corner of {theta} rad beside an arc at joint \
                     {k}, above the {tolerance} rad a tangent-continuous joint may turn at this \
                     chord budget: only corners between two straight segments are filleted"
                )));
            }
            let bisector = u + v;
            if theta >= core::f64::consts::PI - 1e-9 || bisector.length() <= 1e-9 {
                return Err(GeomError::InvalidInput(format!(
                    "swept disk directrix reverses at joint {k}: a mitre at half angle would \
                     hold both legs"
                )));
            }
            sharpest = sharpest.max(theta);
            mitres.push(Some(bisector.normalize()));
        }
        let mut bands: Vec<[Option<MitreBand>; 2]> = vec![[None, None]; m];
        for (k, normal) in mitres.iter().enumerate() {
            let Some(n) = *normal else { continue };
            for (i, end) in [(k, 1), ((k + 1) % m, 0)] {
                let PathPiece::Arc {
                    centre,
                    axis,
                    angle,
                    ..
                } = path[i]
                else {
                    continue;
                };
                // Toward the corner: along the travel at the arc's end,
                // against it at its start.
                let (toward, corner, normal) = if end == 1 {
                    (axis * angle.signum(), path[i].end_point(), n)
                } else {
                    (-axis * angle.signum(), path[i].start_point(), -n)
                };
                let (back, ahead) = band::reach(centre, toward, corner, normal, radius)
                    .ok_or_else(|| {
                        GeomError::InvalidInput(format!(
                            "swept disk directrix mitres a corner beside an arc at joint {k}, but \
                         the mitre plane does not clear the bend's tube within a quarter turn \
                         of the corner: the bend is too tight for a corner this sharp"
                        ))
                    })?;
                bands[i][end] = Some(MitreBand {
                    centre,
                    axis: toward,
                    corner,
                    normal,
                    back,
                    ahead,
                });
            }
        }
        for (i, piece) in path.iter().enumerate() {
            let PathPiece::Arc { angle, .. } = *piece else {
                continue;
            };
            let need: Scalar = bands[i].iter().flatten().map(|b| b.back).sum();
            if need > angle.abs() {
                return Err(GeomError::InvalidInput(format!(
                    "swept disk directrix arc {i} turns {} rad, but the mitre bands beside it \
                     need {need} rad of it: the mitres reach past the arc",
                    angle.abs()
                )));
            }
        }
        let mut pipe = Self {
            path,
            closed,
            radius,
            inner,
            chord,
            mitres,
            sharpest,
            bands,
            twist: 0.0,
            rate: 0.0,
        };
        pipe.check_segment_reach()?;
        if closed {
            pipe.twist = pipe.closing_twist();
            if pipe.twist.abs() * radius > ROUNDING_TWIST * chord {
                let clear: Scalar = (0..m).map(|i| pipe.clear(i).1).sum();
                if clear <= 0.0 {
                    return Err(GeomError::InvalidInput(format!(
                        "swept disk along a closed directrix: a frame carried round the loop \
                         comes back turned by {} rad about the tangent, and no part of the \
                         loop is clear of its mitres to take that twist back",
                        pipe.twist
                    )));
                }
                pipe.rate = pipe.twist / clear;
            }
        }
        Ok(pipe)
    }

    /// The joint before piece `i`, if any.
    fn before(&self, i: usize) -> Option<usize> {
        match i.checked_sub(1) {
            Some(k) => Some(k),
            None => self.closed.then(|| self.path.len() - 1),
        }
    }

    /// The joint after piece `i`, if any.
    fn after(&self, i: usize) -> Option<usize> {
        (self.closed || i + 1 < self.path.len()).then_some(i)
    }

    /// The mitre normals at the start and the end of piece `i`.
    fn mitre_ends(&self, i: usize) -> [Option<Vec3>; 2] {
        [
            self.before(i).and_then(|k| self.mitres[k]),
            self.after(i).and_then(|k| self.mitres[k]),
        ]
    }

    /// A segment's two end planes must not meet inside the tube.
    fn check_segment_reach(&self) -> GeomResult<()> {
        let radius = self.radius;
        for (i, piece) in self.path.iter().enumerate() {
            let PathPiece::Segment { start, end } = *piece else {
                continue;
            };
            let ends = self.mitre_ends(i);
            if ends == [None, None] {
                continue;
            }
            let u = (end - start).normalize();
            let lean = |n: Option<Vec3>| n.map_or(u, |n| n / u.dot(n));
            let g = lean(ends[1]) - lean(ends[0]);
            let reach = radius * (g - u * u.dot(g)).length();
            let length = (end - start).length();
            if reach >= length {
                return Err(GeomError::InvalidInput(format!(
                    "swept disk directrix segment {i} is {length} long, but the mitres at its \
                     ends reach {reach} along it at the disk radius {radius}: the tube would cut \
                     through itself"
                )));
            }
        }
        Ok(())
    }

    /// How far a mitre plane with unit normal `n` through a segment's end
    /// reaches along its direction `u` at the disk radius:
    /// `r tan(theta/2)`.
    fn hold(&self, u: Vec3, n: Vec3) -> Scalar {
        let c = u.dot(n).abs();
        self.radius * (1.0 - c * c).max(0.0).sqrt() / c
    }

    /// Piece `i`'s clear part, where a loop's twist is turned back: a
    /// segment's length outside the mitre zones at its ends, as the length
    /// held before it and its own length (both zero when the zones meet);
    /// a bend's length outside its bands, with nothing held before it.
    fn clear(&self, i: usize) -> (Scalar, Scalar) {
        match self.path[i] {
            PathPiece::Segment { start, end } => {
                let u = (end - start).normalize();
                let [z0, z1] = self
                    .mitre_ends(i)
                    .map(|n| n.map_or(0.0, |n| self.hold(u, n)));
                let clear = (end - start).length() - z0 - z1;
                if clear > 0.0 {
                    (z0, clear)
                } else {
                    (0.0, 0.0)
                }
            }
            PathPiece::Arc {
                centre,
                start,
                angle,
                ..
            } => {
                let banded: Scalar = self.bands[i].iter().flatten().map(|b| b.back).sum();
                (
                    0.0,
                    (start - centre).length() * (angle.abs() - banded).max(0.0),
                )
            }
        }
    }

    /// The holonomy of the loop: a frame carried round it (turning with
    /// each bend, by the least rotation at each joint) comes back to the
    /// first piece's start turned by this angle about its tangent, in
    /// `(-pi, pi]`; turning the returned frame on by it gives the first.
    fn closing_twist(&self) -> Scalar {
        let m = self.path.len();
        let t0 = self.path[0].start_tangent();
        let x0 = seed_reference(t0).normalize();
        let mut x = x0;
        for (i, piece) in self.path.iter().enumerate() {
            if let PathPiece::Arc { axis, angle, .. } = *piece {
                x = rotate(x, Point3::ZERO, axis, angle);
            }
            x = least_rotation(
                x,
                piece.end_tangent(),
                self.path[(i + 1) % m].start_tangent(),
            );
        }
        let x = (x - t0 * t0.dot(x)).normalize();
        x.cross(x0).dot(t0).atan2(x.dot(x0))
    }

    fn plan(&self) -> PipePlan {
        let m = self.path.len();
        let side = |i: usize, end: usize| -> Option<MitreSide> {
            let piece = &self.path[i];
            match piece {
                PathPiece::Segment { .. } => Some(MitreSide::Segment {
                    corner: if end == 1 {
                        piece.end_point()
                    } else {
                        piece.start_point()
                    },
                    direction: piece.start_tangent(),
                }),
                PathPiece::Arc { .. } => self.bands[i][end].map(MitreSide::Arc),
            }
        };
        let arc_mitres = self
            .mitres
            .iter()
            .enumerate()
            .filter_map(|(k, normal)| {
                let normal = (*normal)?;
                let sides = [side(k, 1)?, side((k + 1) % m, 0)?];
                matches!(sides, [MitreSide::Arc(_), _] | [_, MitreSide::Arc(_)]).then_some(
                    ArcMitre {
                        joint: k,
                        normal,
                        sides,
                    },
                )
            })
            .collect();
        PipePlan {
            pieces: self.path.clone(),
            closed: self.closed,
            arc_mitres,
            twist: self.twist,
        }
    }

    fn mesh(&self) -> GeomResult<TriMesh> {
        let rings = sweep::disk_rings(
            self.radius,
            self.inner,
            0.5 * self.chord * (0.5 * self.sharpest).cos(),
        )?;
        let m = self.path.len();
        let mut stations: Vec<Station> = Vec::new();
        let mut carried: Option<Frame> = None;
        let mut turned: Scalar = 0.0;
        for i in 0..m {
            let start = self.start_frame(i, carried.as_ref())?;
            let closing = self.closed && i + 1 == m;
            let close = match (closing, stations.first()) {
                (false, _) => Close::Open,
                (true, Some(first)) => Close::To(first),
                (true, None) => Close::Own,
            };
            let built = self.piece_stations(i, &rings, &start, stations.last(), close, turned)?;
            // Consecutive pieces share one station; a loop's last one is
            // its first.
            let skip = usize::from(!stations.is_empty());
            let keep = built.stations.len() - usize::from(closing);
            stations.extend(built.stations.into_iter().take(keep).skip(skip));
            carried = Some(built.end);
            turned = built.turned;
        }
        let mut mesh = loft::loft(&rings, &stations, self.closed)?;
        // As `revolve::revolve` (#221): settle the winding from the built
        // mesh rather than from the frames' handedness.
        if matches!(crate::extrude::outward_orientation(&mesh), Some(false)) {
            for triangle in mesh.indices.chunks_exact_mut(3) {
                triangle.swap(1, 2);
            }
        }
        Ok(mesh)
    }

    /// Piece `i`'s own start frame: the incoming end frame turned by the
    /// least rotation onto its tangent.
    fn start_frame(&self, i: usize, carried: Option<&Frame>) -> GeomResult<Frame> {
        let piece = &self.path[i];
        let tangent = piece.start_tangent();
        let seed = carried.map_or_else(
            || seed_reference(tangent),
            |f| least_rotation(f.x, f.x.cross(f.y), tangent),
        );
        Frame::from_reference(piece.start_point(), tangent, seed)
    }

    fn piece_stations(
        &self,
        i: usize,
        rings: &Rings,
        start: &Frame,
        joint: Option<&Station>,
        close: Close<'_>,
        turned: Scalar,
    ) -> GeomResult<Built> {
        match self.path[i] {
            PathPiece::Segment { .. } => {
                self.segment_stations(i, rings, start, joint, close, turned)
            }
            PathPiece::Arc { .. } => self.arc_stations(i, rings, start, joint, close, turned),
        }
    }

    /// A segment's stations: its two ends (mitred where a corner is), and
    /// in a twisted loop the ends of its mitre zones and enough stations
    /// between them for the twist.
    fn segment_stations(
        &self,
        i: usize,
        rings: &Rings,
        start: &Frame,
        joint: Option<&Station>,
        close: Close<'_>,
        turned: Scalar,
    ) -> GeomResult<Built> {
        let piece = self.path[i];
        let length = piece.length();
        let tangent = piece.start_tangent();
        let (held, clear) = self.clear(i);
        let tail = if clear > 0.0 {
            length - held - clear
        } else {
            0.0
        };
        let twist = self.rate * clear;
        let sample = |budget: Scalar| -> GeomResult<Sampled> {
            let mut at = vec![0.0];
            if twist != 0.0 {
                let n = twist_steps(twist, budget, self.radius)?;
                if held > 0.0 {
                    at.push(held);
                }
                at.extend((1..n).map(|k| held + clear * k as Scalar / n as Scalar));
                if tail > 0.0 {
                    at.push(length - tail);
                }
            }
            at.push(length);
            Ok(Sampled {
                points: at.iter().map(|s| piece.point(s / length)).collect(),
                tangents: vec![tangent; at.len()],
                twists: at
                    .iter()
                    .map(|s| turned + self.rate * (s - held).clamp(0.0, clear))
                    .collect(),
            })
        };
        let (frames, placed) = bounded_piece(
            sample,
            rings,
            start.x,
            self.radius,
            0.5 * self.chord,
            joint,
            self.mitre_ends(i),
            close,
        )?;
        Ok(Built {
            stations: placed,
            end: last_frame(frames)?,
            turned: turned + twist,
        })
    }

    /// An arc's stations: a mitre band at each end beside a mitred corner,
    /// and the part between proven like any bend.
    fn arc_stations(
        &self,
        i: usize,
        rings: &Rings,
        start: &Frame,
        joint: Option<&Station>,
        close: Close<'_>,
        turned: Scalar,
    ) -> GeomResult<Built> {
        let piece = self.path[i];
        let PathPiece::Arc {
            centre,
            start: from,
            angle,
            axis,
        } = piece
        else {
            return Err(not_an_arc());
        };
        let total = angle.abs();
        let along = axis * angle.signum();
        let reach_radius = (from - centre).length() + self.radius;
        let share = 0.5 * self.chord;
        let band_tol = 0.5 * share;
        let offset = STRIP_OFFSET * self.chord;
        let [first, last] = self.bands[i];
        let (head, tail) = (first.map_or(0.0, |b| b.back), last.map_or(0.0, |b| b.back));
        let place = |frame: &Frame| loft::place(rings, |p| loft::at(&twisted(frame, turned), p));
        let carried = |frame: &Frame, by: Scalar| Frame {
            origin: rotate(frame.origin, centre, along, by),
            x: rotate(frame.x, Point3::ZERO, along, by),
            y: rotate(frame.y, Point3::ZERO, along, by),
        };
        let mut out: Vec<Station> = Vec::new();
        let mut frame = start.clone();
        if let Some(band) = first {
            let k = self.before(i).ok_or_else(not_an_arc)?;
            let shared = self.shared_ring(k, &twisted(start, turned), rings)?;
            frame = carried(start, head);
            let far = place(&frame);
            // The shared ring is this band's own cut unless a segment
            // comes in.
            let strip = matches!(self.path[k], PathPiece::Segment { .. });
            let mut toward = band.stations_toward(
                &far,
                &shared,
                strip.then_some(offset),
                reach_radius,
                band_tol,
            )?;
            toward.reverse();
            out.extend(toward);
            out.push(far);
        }
        let middle = total - head - tail;
        let mut twist = 0.0;
        if middle > 0.0 {
            let part = piece.between(head / total, 1.0 - tail / total);
            twist = self.rate * part.length();
            let sample = |budget: Scalar| -> GeomResult<Sampled> {
                let n = part
                    .steps(budget)?
                    .max(twist_steps(twist, budget, self.radius)?);
                let points: Vec<Point3> = (0..=n)
                    .map(|k| part.point(k as Scalar / n as Scalar))
                    .collect();
                let tangents = points.iter().map(|p| part.tangent_at(*p)).collect();
                let twists = (0..=n)
                    .map(|k| turned + twist * k as Scalar / n as Scalar)
                    .collect();
                Ok(Sampled {
                    points,
                    tangents,
                    twists,
                })
            };
            let (joint, close) = (
                if first.is_some() { out.last() } else { joint },
                if last.is_some() { Close::Open } else { close },
            );
            let (frames, placed) = bounded_piece(
                sample,
                rings,
                frame.x,
                self.radius,
                share,
                joint,
                [None, None],
                close,
            )?;
            let skip = usize::from(first.is_some());
            out.extend(placed.into_iter().skip(skip));
            frame = last_frame(frames)?;
        } else if out.is_empty() {
            out.push(place(start));
        }
        let turned_end = turned + twist;
        if let Some(band) = last {
            let k = self.after(i).ok_or_else(not_an_arc)?;
            let end = carried(&frame, tail);
            let target = match close {
                Close::To(first) => first.clone(),
                Close::Own => out[0].clone(),
                Close::Open => {
                    let next = self.start_frame((k + 1) % self.path.len(), Some(&end))?;
                    self.shared_ring(k, &twisted(&next, turned_end), rings)?
                }
            };
            let far = out.last().cloned().ok_or_else(not_an_arc)?;
            out.extend(band.stations_toward(
                &far,
                &target,
                Some(offset),
                reach_radius,
                band_tol,
            )?);
            frame = end;
        }
        Ok(Built {
            stations: out,
            end: frame,
            turned: turned_end,
        })
    }

    /// The ring shared across mitred joint `k`, from the outgoing piece's
    /// start frame `frame` (twist included): a segment's mitre ring, or
    /// after a segment an arc's section projected the same way (the
    /// incoming segment's own ring, #245); between two arcs the outgoing
    /// band's cut.
    fn shared_ring(&self, k: usize, frame: &Frame, rings: &Rings) -> GeomResult<Station> {
        let next = (k + 1) % self.path.len();
        let normal = self.mitres[k].ok_or_else(not_an_arc)?;
        let square = loft::place(rings, |p| loft::at(frame, p));
        match (self.path[k], self.path[next]) {
            (PathPiece::Arc { .. }, PathPiece::Arc { .. }) => {
                self.bands[next][0].ok_or_else(not_an_arc)?.cut(&square)
            }
            _ => Ok(mitred(&square, frame, normal)),
        }
    }
}

/// The refusal for a plan that does not hold together; never reached.
fn not_an_arc() -> GeomError {
    GeomError::Degenerate("swept disk: a mitre band without its arc or corner".to_owned())
}

fn last_frame(frames: Vec<Frame>) -> GeomResult<Frame> {
    frames
        .into_iter()
        .last()
        .ok_or_else(|| GeomError::Degenerate("swept disk: a piece without a station".to_owned()))
}

/// `frame` turned by `angle` about its own tangent.
fn twisted(frame: &Frame, angle: Scalar) -> Frame {
    if angle == 0.0 {
        return frame.clone();
    }
    let (s, c) = angle.sin_cos();
    Frame {
        origin: frame.origin,
        x: frame.x * c + frame.y * s,
        y: frame.y * c - frame.x * s,
    }
}

/// Spans a twist of `twist` needs for the turn between two frames to move
/// the disk's rim by at most `budget`: `r (1 - cos(phi/2)) <= budget`.
fn twist_steps(twist: Scalar, budget: Scalar, radius: Scalar) -> GeomResult<usize> {
    if twist == 0.0 {
        return Ok(1);
    }
    let per = 2.0 * (1.0 - (budget / radius).min(1.0)).acos();
    let n = (twist.abs() / per).ceil();
    if n.is_nan() || n > MAX_TWIST_STEPS as Scalar {
        return Err(GeomError::BudgetExceeded {
            resource: "swept disk twist steps",
        });
    }
    Ok((n as usize).max(1))
}

/// A piece's frames (untwisted) and stations, its directrix sampled until
/// each span's bound fits `share`; the first span also carries the
/// measured distance from the shared `joint` station to the piece's own
/// first station, and the last the distance from its last station to the
/// one it `close`s on. An end with a mitre normal in `ends` (start, end)
/// has its station projected along the piece's tangent onto that mitre
/// plane.
#[allow(clippy::too_many_arguments)]
fn bounded_piece(
    sample: impl Fn(Scalar) -> GeomResult<Sampled>,
    rings: &Rings,
    seed: Vec3,
    reach: Scalar,
    share: Scalar,
    joint: Option<&Station>,
    ends: [Option<Vec3>; 2],
    close: Close<'_>,
) -> GeomResult<(Vec<Frame>, Vec<Station>)> {
    let mut b = share;
    for _ in 0..MAX_REFINEMENTS {
        let Sampled {
            points,
            tangents,
            twists,
        } = sample(b)?;
        let frames = sweep::frames_carried(&points, &tangents, seed)?;
        let placing: Vec<Frame> = frames
            .iter()
            .zip(&twists)
            .map(|(f, a)| twisted(f, *a))
            .collect();
        let mut placed: Vec<Station> = placing
            .iter()
            .map(|f| loft::place(rings, |p| loft::at(f, p)))
            .collect();
        let last = placed.len() - 1;
        for (k, normal) in [(0, ends[0]), (last, ends[1])] {
            if let Some(n) = normal {
                placed[k] = mitred(&placed[k], &placing[k], n);
            }
        }
        let shift = joint.map_or(0.0, |j| displacement(j, &placed[0]));
        let closing = match close {
            Close::Open => 0.0,
            Close::To(station) => displacement(station, &placed[last]),
            Close::Own => displacement(&placed[0], &placed[last]),
        };
        let worst = placing
            .windows(2)
            .zip(placed.windows(2))
            .enumerate()
            .map(|(k, (f, s))| {
                let moved = if k == 0 { shift } else { 0.0 };
                let closed = if k + 1 == last { closing } else { 0.0 };
                sweep::span_bound(&f[0], &f[1], b, reach)
                    + loft::span_deviation(&s[0], &s[1])
                    + moved
                    + closed
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

/// `station`, square to `frame`'s tangent at its origin, projected along
/// that tangent onto the mitre plane through the origin with unit normal
/// `normal` (see the module notes).
fn mitred(station: &Station, frame: &Frame, normal: Vec3) -> Station {
    let tangent = frame.x.cross(frame.y);
    let along = tangent.dot(normal);
    Station {
        loops: station
            .loops
            .iter()
            .map(|ring| {
                ring.iter()
                    .map(|q| *q - tangent * ((*q - frame.origin).dot(normal) / along))
                    .collect()
            })
            .collect(),
    }
}

/// The refusal of a disk radius equal to a bend's or fillet's radius
/// (see the module notes).
fn horn_torus(what: &str, radius: Scalar) -> GeomError {
    GeomError::InvalidInput(format!(
        "swept disk radius {radius} equals the {what} radius: the format rule permits a \
         fillet radius equal to the disk radius, but the bend is then a horn torus whose \
         inner wall degenerates to one point on the bend's axis, a boundary that is not a \
         two-manifold and that no closed mesh bounds without meeting itself there"
    ))
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
/// `tolerance` with an arc of radius `fillet` (see the module notes); in a
/// loop, the corner where it closes too.
fn fillet(
    pieces: &[PathPiece],
    fillet: Scalar,
    radius: Scalar,
    tolerance: Scalar,
    closed: bool,
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
    for piece in pieces {
        piece.check(radius)?;
        legs.push(Leg {
            piece: *piece,
            dir: piece.start_tangent(),
            trim_start: 0.0,
            trim_end: 0.0,
        });
    }
    let m = legs.len();
    // The fillet at the corner where leg `k` starts.
    let mut arcs: Vec<Option<PathPiece>> = vec![None; m];
    let first = if closed { 0 } else { 1 };
    for k in first..m {
        let prev = (k + m - 1) % m;
        if prev == k {
            continue;
        }
        let (PathPiece::Segment { .. }, PathPiece::Segment { start: corner, .. }) =
            (legs[prev].piece, legs[k].piece)
        else {
            continue;
        };
        let (u, v) = (legs[prev].dir, legs[k].dir);
        let theta = turn(u, v);
        if theta <= tolerance {
            continue;
        }
        if radius > fillet {
            return Err(GeomError::InvalidInput(format!(
                "swept disk radius {radius} is above the fillet radius {fillet}: the \
                 tube would fold through the bend's axis"
            )));
        }
        if radius == fillet {
            return Err(horn_torus("fillet", radius));
        }
        let normal = u.cross(v);
        if normal.length() <= 1e-12 {
            return Err(GeomError::InvalidInput(format!(
                "swept disk directrix reverses at corner {k}; no fillet rounds it"
            )));
        }
        let cut = fillet * (0.5 * theta).tan();
        legs[prev].trim_end = cut;
        legs[k].trim_start = cut;
        let tangent_point = corner - u * cut;
        let inward = (v - u * u.dot(v)).normalize();
        arcs[k] = Some(PathPiece::Arc {
            centre: tangent_point + inward * fillet,
            axis: normal.normalize(),
            start: tangent_point,
            angle: theta,
        });
    }
    let mut out = Vec::with_capacity(2 * m);
    for (k, leg) in legs.iter().enumerate() {
        if k > 0 {
            if let Some(arc) = arcs[k] {
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
    if closed {
        if let Some(arc) = arcs[0] {
            out.push(arc);
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
        let path = fillet(&pieces, 0.2, 0.05, 1e-6, false).unwrap();
        assert_eq!(path.len(), 3);
        for pair in path.windows(2) {
            assert!((pair[0].end_point() - pair[1].start_point()).length() < 1e-12);
            assert!(turn(pair[0].end_tangent(), pair[1].start_tangent()) < 1e-12);
        }
    }

    #[test]
    fn a_mitre_ring_lies_on_both_cylinders_and_on_the_bisector_plane() {
        let corner = Point3::new(1.0, 0.2, -0.3);
        let (u, v) = (
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.2, 0.7, -0.4).normalize(),
        );
        let pieces = [
            PathPiece::Segment {
                start: corner - u,
                end: corner,
            },
            PathPiece::Segment {
                start: corner,
                end: corner + v,
            },
        ];
        let (r, inner) = (0.05, 0.03);
        let mesh = swept_disk_along_pieces(&pieces, r, Some(inner), None, 1e-3).unwrap();
        // One station per vertex: start, mitre, end.
        let per = mesh.positions.len() / 3;
        let n = (u + v).normalize();
        let off = |p: Point3, axis: Vec3| {
            let d = p - corner;
            (d - axis * d.dot(axis)).length()
        };
        for p in &mesh.positions[per..2 * per] {
            assert!((*p - corner).dot(n).abs() < 1e-12);
            let (a, b) = (off(*p, u), off(*p, v));
            assert!((a - b).abs() < 1e-12, "{a} {b}");
            assert!((a - r).abs() < 1e-12 || (a - inner).abs() < 1e-12, "{a}");
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
