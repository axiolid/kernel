//! Profiles that change between stations, flattened together (#361).
//!
//! [`super::profile_rings`] flattens each profile on its own, so two
//! different curved profiles get unrelated vertex sets and no
//! correspondence between them: a run lofted from one to the other has no
//! bound on its profile share. Here the profiles whose rings pair segment
//! for segment are flattened *together*: segment `j` of ring `r` is cut
//! into the same number of equal steps of its own parameter in every
//! profile, and vertex `m` of one is matched to vertex `m` of the others.
//! The loft between two sections interpolates those matched vertices, which
//! is the chord polygon of the interpolated exact section sampled at the
//! same fractions, so its distance from that section is a bound that holds
//! for every interpolated section, not only the stations.
//!
//! # Lowering
//!
//! A profile lowers to segments: a rectangle, a section or a contour
//! through its exact contour; a circle to one full-turn circle (its bore,
//! if any, another); an ellipse to one full-turn ellipse; a derived profile
//! to its basis plus the transforms applied to the flattened points. Each
//! ring is oriented by the sign of its exact area (outer counter-clockwise,
//! holes clockwise) by reversing its segments. A ring reversed in one
//! profile and not in another (a mirrored profile against an unmirrored
//! one) does not pair: its start point stays mirrored, so the vertices would
//! match points of the opposite side. Rings and segments are then matched
//! by index. Anything else (a polyline or B-spline segment, a centre line, a composite,
//! another ring or segment count) does not pair, and the caller keeps the
//! unbounded result it names.
//!
//! # The bound
//!
//! With `m2` the supremum of `|C''|` of a segment over its parameter (zero
//! for a line, `r * hypot(|x|, |y|)` for a circle, `hypot(a |x|, b |y|)`
//! for an ellipse, times the largest stretch of a derived profile's
//! transforms) and `dt` its parameter span, a segment cut into `n` steps
//! lies within `(dt / n)^2 * m2 / 8` of its chords, and so does every
//! convex combination of two such segments. `n` is the smallest count that
//! keeps that within the chord error, taken over all profiles. Joints are
//! snapped as [`super::profile_rings`] does, and the distance a joint moved
//! (times the stretch) is added.
//!
//! A run whose every profile is a circle starts each ring half a step off
//! the axes, as the circle flattener does (#194).

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Interval, Point2, Scalar, Tolerance, Transform2, Vec2};
use axiolid_curve::{Circle2, Curve2, Ellipse2};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};
use axiolid_reference::curve::evaluate2;

use super::{near2, snap_joints, stretch2, Rings};

/// The most vertices the pairing produces over all rings of one profile.
const MAX_VERTICES: Scalar = (1u32 << 20) as Scalar;

/// Rings flattened together, and how far the interpolated sections lie from
/// them.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct PairedRings {
    /// One [`Rings`] per profile, in the order given, matching vertex for
    /// vertex.
    pub rings: Vec<Rings>,
    /// Every point of the exact section interpolated between two of the
    /// profiles, at any fraction, is within this of the polygon the
    /// interpolated rings form, and the reverse.
    pub deviation: Scalar,
}

/// A profile lowered to segments.
struct Lowered {
    /// Outer ring first, then the holes.
    rings: Vec<Vec<ProfileSegment>>,
    /// Applied to the flattened points, innermost first.
    transforms: Vec<Transform2>,
    /// Product of the transforms' largest stretches.
    stretch: Scalar,
    /// Whether [`orient`] reversed each ring.
    reversed: Vec<bool>,
}

/// Flatten `profiles` together under `chord_error`, or `None` when they do
/// not pair segment for segment (see the module docs).
///
/// # Errors
///
/// `InvalidInput` for a chord error that is not positive and finite or a
/// profile that does not validate; as the exact lowering of the profile
/// otherwise.
pub fn paired_profile_rings(
    profiles: &[&Profile],
    chord_error: Scalar,
    tolerance: Tolerance,
) -> GeomResult<Option<PairedRings>> {
    if !chord_error.is_finite() || chord_error <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "paired profile chord error must be positive and finite, got {chord_error}"
        )));
    }
    if profiles.is_empty() {
        return Ok(None);
    }
    let mut lowered = Vec::with_capacity(profiles.len());
    for profile in profiles {
        let Some(mut one) = lower(profile)? else {
            return Ok(None);
        };
        if !orient(&mut one)? {
            return Ok(None);
        }
        lowered.push(one);
    }
    // A ring reversed in one profile and not in another would pair each
    // vertex with a point of the opposite side (a mirrored start point
    // stays mirrored when the ring is reversed).
    if lowered.iter().any(|l| l.reversed != lowered[0].reversed) {
        return Ok(None);
    }
    let first = &lowered[0];
    if lowered.iter().any(|l| {
        l.rings.len() != first.rings.len()
            || l.rings
                .iter()
                .zip(&first.rings)
                .any(|(a, b)| a.len() != b.len())
    }) {
        return Ok(None);
    }
    // The common count of every segment, and the vertices it makes.
    let mut counts: Vec<Vec<usize>> = Vec::with_capacity(first.rings.len());
    let mut total: Scalar = 0.0;
    for (r, ring) in first.rings.iter().enumerate() {
        let mut row = Vec::with_capacity(ring.len());
        for j in 0..ring.len() {
            let mut n: Scalar = 1.0;
            for l in &lowered {
                let segment = &l.rings[r][j];
                let Some(m2) = second_bound(&segment.curve) else {
                    return Ok(None);
                };
                let steps = span(segment)? * (m2 * l.stretch / (8.0 * chord_error)).sqrt();
                n = n.max(steps.ceil());
            }
            total += n;
            if total.is_nan() || total > MAX_VERTICES {
                return Ok(None);
            }
            row.push(n as usize);
        }
        counts.push(row);
    }
    // Every profile a circle: a half step off the axes (#194).
    let circles = profiles.iter().all(|p| matches!(p, Profile::Circle(_)));
    let mut rings = Vec::with_capacity(lowered.len());
    let mut deviation: Scalar = 0.0;
    for l in &lowered {
        let mut moved: Scalar = 0.0;
        let mut chord: Scalar = 0.0;
        let mut flat: Vec<Vec<Point2>> = Vec::with_capacity(l.rings.len());
        for (ring, row) in l.rings.iter().zip(&counts) {
            let mut pieces: Vec<Vec<Point2>> = Vec::with_capacity(ring.len());
            for (segment, &n) in ring.iter().zip(row) {
                let step = span(segment)? / n as Scalar;
                let m2 = second_bound(&segment.curve).unwrap_or(0.0) * l.stretch;
                chord = chord.max(step * step * m2 / 8.0);
                let shift = if circles {
                    core::f64::consts::PI / n as Scalar
                } else {
                    0.0
                };
                pieces.push(sample(segment, n, shift)?);
            }
            let contour = Contour::new(ring.clone());
            snap_joints(&contour, &mut pieces, tolerance, &mut moved, true)?;
            let mut points = Vec::new();
            for mut piece in pieces {
                piece.pop();
                points.extend(piece);
            }
            transform(&mut points, &l.transforms);
            flat.push(points);
        }
        deviation = deviation.max(chord + moved * l.stretch);
        if flat.iter().any(|ring| !distinct(ring, tolerance)) {
            return Ok(None);
        }
        let mut it = flat.into_iter();
        let outer = it.next().unwrap_or_default();
        rings.push(Rings {
            outer,
            holes: it.collect(),
        });
    }
    Ok(Some(PairedRings {
        rings,
        deviation: deviation * (1.0 + 1e-12),
    }))
}

fn transform(points: &mut [Point2], transforms: &[Transform2]) {
    for p in points {
        for t in transforms {
            *p = t.transform_point2(*p);
        }
    }
}

/// Whether a closed ring has at least three vertices, none repeating the
/// one before it (nor the last the first) within the linear tolerance.
fn distinct(ring: &[Point2], tolerance: Tolerance) -> bool {
    ring.len() >= 3
        && (0..ring.len()).all(|i| !near2(ring[i], ring[(i + 1) % ring.len()], tolerance.linear()))
}

/// A segment's parameter span.
fn span(segment: &ProfileSegment) -> GeomResult<Scalar> {
    let span = (segment.domain.end - segment.domain.start).abs();
    if span.is_finite() {
        Ok(span)
    } else {
        Err(GeomError::InvalidInput(
            "a profile segment's parameter domain must be finite".into(),
        ))
    }
}

/// `sup |C''|` over a segment's parameter: a line has none, a circle
/// `r * hypot(|x|, |y|)`, an ellipse `hypot(a |x|, b |y|)`; other families
/// do not pair.
fn second_bound(curve: &Curve2) -> Option<Scalar> {
    match curve {
        Curve2::Line(_) => Some(0.0),
        Curve2::Circle(c) => Some(c.radius * c.frame.x.length().hypot(c.frame.y.length())),
        Curve2::Ellipse(e) => {
            Some((e.semi_axis_x * e.frame.x.length()).hypot(e.semi_axis_y * e.frame.y.length()))
        }
        _ => None,
    }
}

/// `count + 1` points of a segment in contour order, at equal fractions of
/// its parameter span from its start in that order, the domain moved by
/// `shift` first.
fn sample(segment: &ProfileSegment, count: usize, shift: Scalar) -> GeomResult<Vec<Point2>> {
    let (lo, hi) = (segment.domain.start + shift, segment.domain.end + shift);
    let (from, to) = if segment.same_sense {
        (lo, hi)
    } else {
        (hi, lo)
    };
    (0..=count)
        .map(|m| {
            let t = if m == 0 {
                from
            } else if m == count {
                to
            } else {
                from + (to - from) * (m as Scalar / count as Scalar)
            };
            evaluate2(&segment.curve, t)
        })
        .collect()
}

/// Reverse every ring whose area has the wrong sign: the outer
/// counter-clockwise, holes clockwise. `false` when a ring's area is zero or
/// not finite, so its orientation is not known.
fn orient(lowered: &mut Lowered) -> GeomResult<bool> {
    let det: Scalar = lowered
        .transforms
        .iter()
        .map(|t| t.matrix2.determinant())
        .product();
    for (r, ring) in lowered.rings.iter_mut().enumerate() {
        let mut area = 0.0;
        for segment in ring.iter() {
            area += segment_area(segment)?;
        }
        area *= det;
        if area == 0.0 || !area.is_finite() {
            return Ok(false);
        }
        let reverse = (r == 0 && area < 0.0) || (r > 0 && area > 0.0);
        if reverse {
            ring.reverse();
            for segment in ring.iter_mut() {
                segment.same_sense = !segment.same_sense;
            }
        }
        lowered.reversed.push(reverse);
    }
    Ok(true)
}

/// `1/2 * integral of (x dy - y dx)` along a segment in contour order: the
/// exact signed area it contributes to a closed ring. A line is the cross
/// product of its end points; a circle or ellipse `O + A cos t + B sin t`
/// gives `(A x B)(t1 - t0) + O x B (sin t1 - sin t0) + O x A (cos t1 - cos t0)`.
fn segment_area(segment: &ProfileSegment) -> GeomResult<Scalar> {
    let (t0, t1) = (segment.domain.start, segment.domain.end);
    let cross = |a: Vec2, b: Vec2| a.x * b.y - a.y * b.x;
    let conic = |o: Point2, a: Vec2, b: Vec2| {
        let o = Vec2::new(o.x, o.y);
        0.5 * (cross(a, b) * (t1 - t0)
            + cross(o, b) * (t1.sin() - t0.sin())
            + cross(o, a) * (t1.cos() - t0.cos()))
    };
    let forward = match &segment.curve {
        Curve2::Line(_) => {
            let p = evaluate2(&segment.curve, t0)?;
            let q = evaluate2(&segment.curve, t1)?;
            0.5 * (p.x * q.y - q.x * p.y)
        }
        Curve2::Circle(c) => conic(c.frame.origin, c.frame.x * c.radius, c.frame.y * c.radius),
        Curve2::Ellipse(e) => conic(
            e.frame.origin,
            e.frame.x * e.semi_axis_x,
            e.frame.y * e.semi_axis_y,
        ),
        _ => return Ok(Scalar::NAN),
    };
    Ok(if segment.same_sense {
        forward
    } else {
        -forward
    })
}

/// A profile as segments, or `None` for a family that does not pair.
fn lower(profile: &Profile) -> GeomResult<Option<Lowered>> {
    let plain = |rings: Vec<Vec<ProfileSegment>>| Lowered {
        rings,
        transforms: Vec::new(),
        stretch: 1.0,
        reversed: Vec::new(),
    };
    let contour = |c: &ContourProfile| {
        plain(
            core::iter::once(&c.outer)
                .chain(&c.holes)
                .map(|contour| contour.segments.clone())
                .collect(),
        )
    };
    Ok(match profile {
        Profile::Rectangle(r) => Some(contour(&crate::section_lower::rectangle_contour(r)?)),
        Profile::Section(s) => Some(contour(&crate::section_lower::section_contour(s)?)),
        Profile::Contour(c) => Some(contour(c)),
        Profile::Circle(c) => {
            if !c.radius.is_finite() || c.radius <= 0.0 {
                return Err(GeomError::InvalidInput(format!(
                    "circle profile radius must be positive and finite, got {}",
                    c.radius
                )));
            }
            let circle = |radius| {
                full_turn(Curve2::Circle(Circle2 {
                    frame: unit_frame(),
                    radius,
                }))
            };
            let mut rings = vec![vec![circle(c.radius)]];
            if let Some(t) = c.thickness {
                if !(t > 0.0 && t < c.radius) {
                    return Err(GeomError::InvalidInput(format!(
                        "annulus wall thickness {t} does not fit inside radius {}",
                        c.radius
                    )));
                }
                rings.push(vec![circle(c.radius - t)]);
            }
            Some(plain(rings))
        }
        Profile::Ellipse(e) => {
            if !(e.semi_axis_x.is_finite() && e.semi_axis_x > 0.0)
                || !(e.semi_axis_y.is_finite() && e.semi_axis_y > 0.0)
            {
                return Err(GeomError::InvalidInput(format!(
                    "ellipse semi-axes must be positive and finite, got {} and {}",
                    e.semi_axis_x, e.semi_axis_y
                )));
            }
            Some(plain(vec![vec![full_turn(Curve2::Ellipse(Ellipse2 {
                frame: unit_frame(),
                semi_axis_x: e.semi_axis_x,
                semi_axis_y: e.semi_axis_y,
            }))]]))
        }
        Profile::Derived { basis, transform } => lower(basis)?.map(|mut inner| {
            inner.transforms.push(*transform);
            inner.stretch *= stretch2(transform);
            inner
        }),
        _ => None,
    })
}

fn unit_frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(0.0, 0.0),
        x: Vec2::new(1.0, 0.0),
        y: Vec2::new(0.0, 1.0),
    }
}

fn full_turn(curve: Curve2) -> ProfileSegment {
    ProfileSegment {
        curve,
        domain: Interval {
            start: 0.0,
            end: core::f64::consts::TAU,
        },
        same_sense: true,
    }
}
