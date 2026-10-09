//! A bounded half-space's profile boundary, checked once for both
//! compilers (#277).
//!
//! A boundary of line and circular-arc segments is a profile node (see
//! `SolidOperation::BoundedHalfSpace`). Neither builder checks all of what
//! makes such a contour bound a region: the mesh flattener bridges a gap
//! between segments with a chord, and the exact extruder sweeps a contour
//! that crosses itself into walls that pass through each other, which the
//! general boolean then refuses under an unrelated name, or not at all
//! when the crossing lies outside the subject. So both compilers run
//! [`check`] first, and refuse the same contours under the same names:
//!
//! - **Open:** every contour is lowered to its ring of joints by
//!   `contour_to_arc_ring`, which refuses a gap wider than the tolerance
//!   (or the rounding of the two endpoint evaluations), a segment that is
//!   neither a line nor a circular arc, and an arc of no or more than one
//!   turn, each by name.
//! - **Zero radius:** an arc whose radius is not positive and finite
//!   bounds nothing; it is refused before lowering.
//! - **Self-crossing:** every two edges of the rings, lines and arcs below
//!   half a turn each, are intersected in closed form. Two edges that meet
//!   anywhere within the linear tolerance (or the rounding of those closed
//!   forms, when larger), other than at the joint two neighbours share,
//!   cross or touch: the contour is refused, naming them. Holes are checked
//!   against the outer ring and each other alike. This is a decision within
//!   the tolerance, not an exact predicate: a contour whose edges only come
//!   within the tolerance of each other is refused as touching, and a
//!   crossing is found wherever it lies further than the rounding from a
//!   tangency.
//!
//! The other profile families are well-formed by construction or checked
//! where they are built; a derived profile is checked through its basis (an
//! invertible affine map neither opens nor crosses a contour) and a
//! composite through its members, whose union the exact extruder forms.

use axiolid_construct::contour_lower::contour_to_arc_ring;
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Scalar, Tolerance, Vec2};
use axiolid_curve::Curve2;
use axiolid_profile::{Contour, Profile};

/// The rounding of the closed-form meetings, per unit of the magnitudes
/// they are computed from.
const ROUNDING: Scalar = 256.0 * Scalar::EPSILON;

/// Refuse a profile boundary that does not bound a region, by name.
pub(crate) fn check(profile: &Profile, tolerance: Tolerance) -> GeomResult<()> {
    match profile {
        Profile::Contour(contour) => {
            let mut edges = Vec::new();
            for (ring, contour) in core::iter::once(&contour.outer)
                .chain(&contour.holes)
                .enumerate()
            {
                ring_edges(contour, ring, tolerance, &mut edges)?;
            }
            uncrossed(&edges, tolerance)
        }
        Profile::Derived { basis, .. } => check(basis, tolerance),
        Profile::Composite(members) => members
            .iter()
            .try_for_each(|member| check(member, tolerance)),
        _ => Ok(()),
    }
}

/// One edge of a ring: a line, or an arc below half a turn.
#[derive(Debug, Clone, Copy)]
struct Edge {
    /// Which ring, and the edge's place in it.
    ring: usize,
    index: usize,
    /// Of edges in the ring.
    count: usize,
    from: Point2,
    to: Point2,
    /// `tan(sweep / 4)`, zero for a line.
    bulge: Scalar,
}

impl Edge {
    /// The arc's centre and radius, `None` for a line.
    fn circle(&self) -> Option<(Point2, Scalar)> {
        if self.bulge == 0.0 {
            return None;
        }
        let chord = self.to - self.from;
        let d = chord.length();
        let b = self.bulge;
        let left = Vec2::new(-chord.y, chord.x) / d;
        let centre = (self.from + self.to) * 0.5 + left * (d * (1.0 - b * b) / (4.0 * b));
        Some((centre, d * (1.0 + b * b) / (4.0 * b.abs())))
    }

    /// Whether `p`, on this arc's circle, lies on the arc: on the far side
    /// of the chord from the centre, which for less than half a turn is
    /// the arc, or within `slack` of an end.
    fn holds(&self, p: Point2, slack: Scalar) -> bool {
        if p.distance(self.from) <= slack || p.distance(self.to) <= slack {
            return true;
        }
        let chord = self.to - self.from;
        chord.perp_dot(p - self.from) * self.bulge.signum() <= 0.0
    }

    /// Whether `other` follows this edge in their ring, or this one it, and
    /// the joint they share.
    fn joint(&self, other: &Self) -> Option<Point2> {
        if self.ring != other.ring {
            return None;
        }
        if (self.index + 1) % self.count == other.index {
            Some(self.to)
        } else if (other.index + 1) % other.count == self.index {
            Some(self.from)
        } else {
            None
        }
    }
}

/// The edges of one contour, refused when it does not close or an arc has
/// no radius.
fn ring_edges(
    contour: &Contour,
    ring: usize,
    tolerance: Tolerance,
    edges: &mut Vec<Edge>,
) -> GeomResult<()> {
    for segment in &contour.segments {
        if let Curve2::Circle(circle) = &segment.curve {
            if !(circle.radius.is_finite() && circle.radius > 0.0) {
                return Err(GeomError::InvalidInput(format!(
                    "half-space boundary arc has radius {}, not a positive one",
                    circle.radius
                )));
            }
        }
    }
    let lowered = contour_to_arc_ring(contour, tolerance)?;
    let vertices = &lowered.vertices;
    let count = vertices.len();
    for (index, vertex) in vertices.iter().enumerate() {
        let to = vertices[(index + 1) % count].point;
        if vertex.point == to {
            return Err(GeomError::InvalidInput(format!(
                "half-space boundary has a segment of no length at {to:?}"
            )));
        }
        edges.push(Edge {
            ring,
            index,
            count,
            from: vertex.point,
            to,
            bulge: vertex.bulge,
        });
    }
    Ok(())
}

/// Refuse any two edges that meet, other than neighbours at their joint.
fn uncrossed(edges: &[Edge], tolerance: Tolerance) -> GeomResult<()> {
    // The tolerance, or the rounding of the closed forms when that is
    // larger (at `Tolerance::ZERO`): a few hundred ulps of the largest
    // coordinate or radius they are computed from.
    let scale = edges.iter().fold(0.0, |scale: Scalar, edge| {
        let radius = edge.circle().map_or(0.0, |(centre, r)| centre.length() + r);
        scale
            .max(edge.from.length())
            .max(edge.to.length())
            .max(radius)
    });
    let slack = tolerance.linear().max(ROUNDING * scale);
    for (i, a) in edges.iter().enumerate() {
        for b in &edges[i + 1..] {
            let joint = a.joint(b);
            // Neighbours meet at their joint, and a ring of two edges at
            // both ends; anywhere else they cross, touch or fold back.
            let elsewhere = |p: Point2| {
                let at_joint = joint.is_some_and(|q| p.distance(q) <= slack)
                    || (a.count == 2 && a.ring == b.ring && p.distance(a.from) <= slack);
                !at_joint
            };
            if let Some(p) = meeting(a, b, slack).into_iter().find(|&p| elsewhere(p)) {
                return Err(GeomError::InvalidInput(format!(
                    "half-space boundary crosses or touches itself: edge {} of ring {} and \
                     edge {} of ring {} meet at {p:?}",
                    a.index, a.ring, b.index, b.ring
                )));
            }
        }
    }
    Ok(())
}

/// Points where two edges meet within `slack`: their intersections, and
/// where they overlap, the ends of each that lie on the other.
fn meeting(a: &Edge, b: &Edge, slack: Scalar) -> Vec<Point2> {
    match (a.circle(), b.circle()) {
        (None, None) => line_line(a, b, slack),
        (Some(circle), None) => line_arc(b, a, circle, slack),
        (None, Some(circle)) => line_arc(a, b, circle, slack),
        (Some(ca), Some(cb)) => arc_arc(a, ca, b, cb, slack),
    }
}

/// The distance from `p` to the segment `from -> to`.
fn to_segment(p: Point2, from: Point2, to: Point2) -> Scalar {
    let d = to - from;
    let t = ((p - from).dot(d) / d.length_squared()).clamp(0.0, 1.0);
    p.distance(from + d * t)
}

fn line_line(a: &Edge, b: &Edge, slack: Scalar) -> Vec<Point2> {
    let mut out = Vec::new();
    let (u, v) = (a.to - a.from, b.to - b.from);
    let denominator = u.perp_dot(v);
    let w = b.from - a.from;
    if denominator != 0.0 {
        let s = w.perp_dot(v) / denominator;
        let t = w.perp_dot(u) / denominator;
        if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t) {
            out.push(a.from + u * s);
        }
    }
    // Near-parallel or touching within the slack: each end that lies on
    // the other edge.
    for (p, edge) in [(a.from, b), (a.to, b), (b.from, a), (b.to, a)] {
        if to_segment(p, edge.from, edge.to) <= slack {
            out.push(p);
        }
    }
    out
}

fn line_arc(
    line: &Edge,
    arc: &Edge,
    (centre, radius): (Point2, Scalar),
    slack: Scalar,
) -> Vec<Point2> {
    let mut out = Vec::new();
    let u = line.to - line.from;
    let w = line.from - centre;
    // |w + u s|^2 = r^2. A line passing within the slack outside the
    // circle has no root; the clamp puts both at its closest approach to
    // the centre, where it touches.
    let (qa, qb) = (u.length_squared(), w.dot(u));
    let qc = w.length_squared() - radius * radius;
    let discriminant = qb * qb - qa * qc;
    let root = discriminant.max(0.0).sqrt();
    for s in [(-qb - root) / qa, (-qb + root) / qa] {
        if (0.0..=1.0).contains(&s) {
            let p = line.from + u * s;
            if (p.distance(centre) - radius).abs() <= slack && arc.holds(p, slack) {
                out.push(p);
            }
        }
    }
    for p in [line.from, line.to] {
        if (p.distance(centre) - radius).abs() <= slack && arc.holds(p, slack) {
            out.push(p);
        }
    }
    for p in [arc.from, arc.to] {
        if to_segment(p, line.from, line.to) <= slack {
            out.push(p);
        }
    }
    out
}

fn arc_arc(
    a: &Edge,
    (ca, ra): (Point2, Scalar),
    b: &Edge,
    (cb, rb): (Point2, Scalar),
    slack: Scalar,
) -> Vec<Point2> {
    let mut out = Vec::new();
    let on_a = |p: Point2| (p.distance(ca) - ra).abs() <= slack && a.holds(p, slack);
    let on_b = |p: Point2| (p.distance(cb) - rb).abs() <= slack && b.holds(p, slack);
    let between = cb - ca;
    let d = between.length();
    if d > 0.0 && d <= ra + rb + slack && d + ra.min(rb) + slack >= ra.max(rb) {
        // Along the centre line to the chord of the two circles, clamped so
        // that circles within the slack of touching meet at one point.
        let x = (d * d + ra * ra - rb * rb) / (2.0 * d);
        let h = (ra * ra - x * x).max(0.0).sqrt();
        let along = between / d;
        let across = Vec2::new(-along.y, along.x);
        for p in [ca + along * x + across * h, ca + along * x - across * h] {
            if on_a(p) && on_b(p) {
                out.push(p);
            }
        }
    }
    // One circle within the slack of the other: the ends of each on the
    // other.
    for p in [b.from, b.to] {
        if on_a(p) {
            out.push(p);
        }
    }
    for p in [a.from, a.to] {
        if on_b(p) {
            out.push(p);
        }
    }
    out
}
