//! Lower `Derived` and `Composite` profiles onto the concrete builders
//! (ADR 0054).
//!
//! # Derived
//!
//! A derived profile is a basis plus an affine transform. Lowering it means
//! pushing the transform down onto the basis geometry, which is only sound
//! when the transform preserves the SHAPE KIND of what it moves:
//!
//! - A circle stays a circle only under a CONFORMAL (similarity) linear part.
//!   A non-uniform scale or a shear turns it into an ellipse, which this
//!   crate cannot extrude, so those are refused rather than silently
//!   relabelled. Note a shear has determinant 1, so a determinant check alone
//!   would let it through.
//! - A negative determinant MIRRORS, reversing ring orientation. The
//!   extruders assume counter-clockwise outer rings, so mirrored rings are
//!   reversed on lowering rather than handed over inside-out.
//!
//! # Composite
//!
//! `Profile::Composite` is the `IfcCompositeProfileDef` case: several
//! profiles acting as ONE section. Members may touch or overlap, so they are
//! unioned rather than concatenated -- concatenating overlapping members
//! would double-count the shared area and produce self-intersecting walls.
//!
//! A composite whose members are mutually DISJOINT is refused: the result is
//! two separate solids, and `Solid` holds one outer shell plus voids, so
//! there is nowhere honest to put the second body.

use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Interval, Point2, Scalar, Tolerance, Transform2, Vec2};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_overlay::{union_soup, Ring};
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile,
};

use crate::BACKEND_ID;

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BACKEND_ID,
        operation: Operation::Sweep,
        input,
    }
}

/// Whether the linear part of `transform` is a similarity.
///
/// `M^T M` is a positive multiple of the identity exactly when `M` scales all
/// directions equally and preserves angles. Checking the determinant is NOT
/// enough: a shear has determinant 1 and still distorts circles into
/// ellipses.
fn conformal_scale(transform: &Transform2, tolerance: Tolerance) -> Option<Scalar> {
    let x = transform.matrix2.x_axis;
    let y = transform.matrix2.y_axis;
    let xx = x.dot(x);
    let yy = y.dot(y);
    let xy = x.dot(y);
    if xx <= 0.0 || yy <= 0.0 {
        return None;
    }
    // Scale-relative comparison: the entries are squared lengths, so the
    // slack must be relative to them rather than an absolute epsilon.
    let slack = tolerance.linear() * xx.max(yy).max(1.0);
    if xy.abs() > slack || (xx - yy).abs() > slack {
        return None;
    }
    Some(xx.sqrt())
}

/// Apply an affine transform to a planar point.
fn apply(transform: &Transform2, point: Point2) -> Point2 {
    transform.transform_point2(point)
}

/// Push an affine transform down onto a basis profile.
///
/// Returns a profile that needs no further transformation. Nesting is handled
/// by composing transforms rather than recursing on already-lowered output,
/// so a deeply derived profile costs one pass.
pub fn lower_derived(
    basis: &Profile,
    transform: &Transform2,
    tolerance: Tolerance,
) -> GeomResult<Profile> {
    match basis {
        // Composing first keeps the recursion one level deep regardless of
        // how many times a profile was re-placed.
        Profile::Derived {
            basis: inner,
            transform: inner_transform,
        } => lower_derived(inner, &(*transform * *inner_transform), tolerance),

        Profile::Rectangle(rectangle) => lower_rectangle(rectangle, transform, tolerance),
        Profile::Circle(circle) => lower_circle(circle, transform, tolerance),
        Profile::Contour(contour) => Ok(Profile::Contour(ContourProfile {
            outer: lower_contour(&contour.outer, transform, tolerance)?,
            holes: contour
                .holes
                .iter()
                .map(|hole| lower_contour(hole, transform, tolerance))
                .collect::<GeomResult<Vec<_>>>()?,
        })),
        Profile::Composite(members) => Ok(Profile::Composite(
            members
                .iter()
                .map(|member| lower_derived(member, transform, tolerance))
                .collect::<GeomResult<Vec<_>>>()?,
        )),
        _ => Err(unsupported("derived profile over an unsupported basis")),
    }
}

/// A transformed rectangle is a rectangle only under a similarity.
///
/// Under a general affine map it becomes a parallelogram, which the rectangle
/// profile cannot express, so it is lowered to an explicit contour instead of
/// being refused: the contour path handles arbitrary polygons exactly.
fn lower_rectangle(
    rectangle: &RectangleProfile,
    transform: &Transform2,
    tolerance: Tolerance,
) -> GeomResult<Profile> {
    if rectangle.thickness.is_some()
        || rectangle.outer_radius.is_some()
        || rectangle.inner_radius.is_some()
    {
        // Rounded corners and a hollow core both lower to an exact contour.
        // Transforming THAT keeps each arc an arc under a similarity and
        // refuses the ellipse a shear would make, through the same check
        // every other contour goes through.
        let contour = crate::section_lower::rectangle_contour(rectangle)?;
        return lower_derived(&Profile::Contour(contour), transform, tolerance);
    }
    if !rectangle.x.is_finite()
        || !rectangle.y.is_finite()
        || rectangle.x <= 0.0
        || rectangle.y <= 0.0
    {
        return Err(GeomError::InvalidInput(format!(
            "derived rectangle extents must be positive and finite, got {} x {}",
            rectangle.x, rectangle.y
        )));
    }

    let (half_x, half_y) = (rectangle.x / 2.0, rectangle.y / 2.0);
    let corners = [
        Point2::new(-half_x, -half_y),
        Point2::new(half_x, -half_y),
        Point2::new(half_x, half_y),
        Point2::new(-half_x, half_y),
    ];
    let mut moved: Vec<Point2> = corners.iter().map(|p| apply(transform, *p)).collect();
    orient_counter_clockwise(&mut moved);
    let _ = tolerance;
    Ok(Profile::Contour(ContourProfile {
        outer: polygon_contour(&moved),
        holes: Vec::new(),
    }))
}

/// A transformed circle is a circle only under a similarity.
fn lower_circle(
    circle: &CircleProfile,
    transform: &Transform2,
    tolerance: Tolerance,
) -> GeomResult<Profile> {
    let Some(scale) = conformal_scale(transform, tolerance) else {
        // A shear or non-uniform scale makes this an ellipse. Relabelling it
        // a circle would report a wrong radius; building it as a polygon
        // would silently discard exactness.
        return Err(unsupported(
            "derived circle under a non-conformal transform is an ellipse",
        ));
    };
    if transform.translation != Vec2::ZERO {
        // The extruder places circles at the origin, so an off-origin circle
        // has nowhere to record its centre.
        return Err(unsupported("derived circle translated off the origin"));
    }
    Ok(Profile::Circle(CircleProfile {
        radius: circle.radius * scale,
        thickness: circle.thickness.map(|value| value * scale),
    }))
}

/// Transform every segment of a contour, preserving each segment's kind.
///
/// Segments are transformed as GEOMETRY, not sampled: a line stays a line and
/// a circular arc stays a circular arc. Taking only the segment endpoints
/// would quietly turn every arc into a chord, which is the same silent
/// approximation the contour path already refuses.
fn lower_contour(
    contour: &Contour,
    transform: &Transform2,
    tolerance: Tolerance,
) -> GeomResult<Contour> {
    let mirrors = transform.matrix2.determinant() < 0.0;
    let mut segments = Vec::with_capacity(contour.segments.len());
    for segment in &contour.segments {
        segments.push(lower_segment(segment, transform, tolerance, mirrors)?);
    }
    if mirrors {
        // A mirror reverses the boundary's sense, so the ring would come out
        // clockwise and the solid inside-out. Reversing the segment order and
        // each segment's own sense restores the original orientation.
        segments.reverse();
    }
    Ok(Contour::new(segments))
}

fn lower_segment(
    segment: &ProfileSegment,
    transform: &Transform2,
    tolerance: Tolerance,
    mirrors: bool,
) -> GeomResult<ProfileSegment> {
    let same_sense = segment.same_sense != mirrors;
    match &segment.curve {
        Curve2::Line(line) => Ok(ProfileSegment {
            curve: Curve2::Line(Line2 {
                origin: apply(transform, line.origin),
                // A direction is a vector, so it takes the linear part only;
                // translating it would move the line off its own points.
                direction: transform.matrix2 * line.direction,
            }),
            domain: segment.domain,
            same_sense,
        }),
        Curve2::Circle(circle) => {
            let Some(scale) = conformal_scale(transform, tolerance) else {
                return Err(unsupported(
                    "derived contour arc under a non-conformal transform is elliptical",
                ));
            };
            Ok(ProfileSegment {
                curve: Curve2::Circle(Circle2 {
                    frame: Frame2 {
                        origin: apply(transform, circle.frame.origin),
                        x: transform.matrix2 * circle.frame.x,
                        y: transform.matrix2 * circle.frame.y,
                    },
                    radius: circle.radius * scale,
                }),
                domain: segment.domain,
                same_sense,
            })
        }
        _ => Err(unsupported(
            "derived contour over a segment kind that cannot be transformed exactly",
        )),
    }
}

/// Build a closed contour of straight segments through `points`.
fn polygon_contour(points: &[Point2]) -> Contour {
    let count = points.len();
    let segments = (0..count)
        .map(|index| {
            let from = points[index];
            let to = points[(index + 1) % count];
            ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: from,
                    direction: to - from,
                }),
                domain: Interval::UNIT,
                same_sense: true,
            }
        })
        .collect();
    Contour::new(segments)
}

/// Signed area; positive when the ring runs counter-clockwise.
fn signed_area(points: &[Point2]) -> Scalar {
    let count = points.len();
    (0..count)
        .map(|index| {
            let a = points[index];
            let b = points[(index + 1) % count];
            a.perp_dot(b)
        })
        .sum::<Scalar>()
        / 2.0
}

/// Reverse a ring in place if it runs clockwise.
///
/// The extruders assume counter-clockwise outer rings; a mirrored transform
/// produces clockwise ones, which would build the solid inside-out.
fn orient_counter_clockwise(points: &mut [Point2]) {
    if signed_area(points) < 0.0 {
        points.reverse();
    }
}

/// Union the members of a composite profile into one section.
///
/// Members of an `IfcCompositeProfileDef` act as a single section and may
/// touch or overlap, so they are UNIONED. Concatenating them as rings would
/// double-count shared area and emit self-intersecting walls.
///
/// Measured behaviour of the union: overlapping members merge to one polygon,
/// edge-touching members merge to one, four bars arranged as a picture frame
/// merge to one polygon carrying one hole, and disjoint members stay two
/// polygons -- which is the case this refuses.
pub fn lower_composite(
    members: &[Profile],
    tolerance: Tolerance,
) -> GeomResult<(Vec<Point2>, Vec<Vec<Point2>>)> {
    if members.is_empty() {
        return Err(GeomError::InvalidInput(
            "a composite profile needs at least one member".to_owned(),
        ));
    }

    let mut rings = Vec::with_capacity(members.len());
    for member in members {
        rings.push(Ring {
            points: member_ring(member, tolerance)?,
        });
    }

    let polygons = union_soup(&rings, tolerance).map_err(|error| {
        GeomError::InvalidInput(format!("composite member union failed: {error:?}"))
    })?;

    match polygons.len() {
        0 => Err(GeomError::Degenerate(
            "composite profile members union to nothing".to_owned(),
        )),
        1 => {
            let polygon = &polygons[0];
            Ok((
                polygon.outer.points.clone(),
                polygon.holes.iter().map(|h| h.points.clone()).collect(),
            ))
        }
        // Disjoint members are two separate bodies. `Solid` holds one outer
        // shell plus voids, so there is nowhere honest to put the second.
        _ => Err(unsupported("composite profile whose members are disjoint")),
    }
}

/// One composite member as a closed polygon ring.
///
/// Members are reduced to rings because the union operates on polygons. A
/// member carrying arcs is refused rather than sampled: chord-sampling here
/// would defeat the exactness the contour path was built to preserve.
fn member_ring(member: &Profile, tolerance: Tolerance) -> GeomResult<Vec<Point2>> {
    let lowered;
    let resolved = match member {
        Profile::Derived { basis, transform } => {
            lowered = lower_derived(basis, transform, tolerance)?;
            &lowered
        }
        other => other,
    };

    match resolved {
        Profile::Rectangle(rectangle) => {
            if rectangle.thickness.is_some()
                || rectangle.outer_radius.is_some()
                || rectangle.inner_radius.is_some()
            {
                return Err(unsupported(
                    "composite member with a hollow or rounded rectangle",
                ));
            }
            let (half_x, half_y) = (rectangle.x / 2.0, rectangle.y / 2.0);
            Ok(vec![
                Point2::new(-half_x, -half_y),
                Point2::new(half_x, -half_y),
                Point2::new(half_x, half_y),
                Point2::new(-half_x, half_y),
            ])
        }
        Profile::Contour(contour) => {
            if !contour.holes.is_empty() {
                return Err(unsupported("composite member carrying its own holes"));
            }
            let mut points: Vec<Point2> = Vec::with_capacity(contour.outer.segments.len());
            for segment in &contour.outer.segments {
                match &segment.curve {
                    Curve2::Line(line) => {
                        let t = if segment.same_sense {
                            segment.domain.start
                        } else {
                            segment.domain.end
                        };
                        points.push(line.origin + line.direction * t);
                    }
                    _ => return Err(unsupported("composite member with a curved segment")),
                }
            }
            orient_counter_clockwise(&mut points);
            Ok(points)
        }
        _ => Err(unsupported(
            "composite member of an unsupported profile kind",
        )),
    }
}
