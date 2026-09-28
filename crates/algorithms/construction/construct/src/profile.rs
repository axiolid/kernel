//! Profile -> 2D polygon rings, then triangles.
//!
//! Curved boundaries are flattened to chords under an explicit
//! `TessellationOptions`-style budget; nothing here invents a default
//! tolerance.
//!
//! Triangulation of rings-with-holes is delegated to `earcut` (MIT/Apache-2.0,
//! pure Rust). ADR 0015 records why: hole bridging is a solved problem and a
//! hand-rolled version failed its own area gate on the two-hole case.
//! `axiolid_reference::triangulate_simple` is retained as the differential oracle for
//! the hole-free case, so the adopted implementation is audited, not trusted.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Scalar, Tolerance};
use axiolid_profile::{CircleProfile, EllipseProfile, Profile, RectangleProfile};

/// Outer ring plus holes, all CCW/CW normalised by the caller's contract:
/// outer counter-clockwise, holes clockwise.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Rings {
    /// Outer boundary, counter-clockwise.
    pub outer: Vec<Point2>,
    /// Inner boundaries, clockwise.
    pub holes: Vec<Vec<Point2>>,
}

/// Flatten a profile into rings under an explicit chord budget.
///
/// Only the families a format adapter currently emits are handled. Everything
/// else returns `Unsupported` rather than a silently wrong approximation --
/// a wrong wall is far more expensive than a missing one.
pub fn profile_rings(
    profile: &Profile,
    chord_error: Scalar,
    tolerance: Tolerance,
) -> GeomResult<Rings> {
    match profile {
        // Rounded corners are exact arcs in the contour, chorded like any
        // other arc; a sharp rectangle keeps its four exact corners.
        Profile::Rectangle(r) if r.outer_radius.is_some() || r.inner_radius.is_some() => {
            let contour = crate::section_lower::rectangle_contour(r)?;
            contour_rings(&contour, chord_error, tolerance)
        }
        Profile::Rectangle(r) => rectangle_rings(r, chord_error, tolerance),
        // Structural sections (#193): the same exact contour the exact
        // extruder builds, fillets and toe radii as arcs, chorded here.
        Profile::Section(section) => {
            let contour = crate::section_lower::section_contour(section)?;
            contour_rings(&contour, chord_error, tolerance)
        }
        Profile::Circle(c) => circle_rings(c, chord_error),
        // Reachable only since curve evaluation moved into `axiolid-reference`:
        // an ellipse has no closed-form segment count, so the old fixed-count
        // flattener could not express it at all.
        Profile::Ellipse(e) => ellipse_rings(e, chord_error),
        Profile::Contour(c) => contour_rings(c, chord_error, tolerance),
        Profile::Derived { basis, transform } => {
            let mut rings = profile_rings(basis, chord_error, tolerance)?;
            apply2(&mut rings.outer, transform);
            for hole in &mut rings.holes {
                apply2(hole, transform);
            }
            // A mirroring placement reverses ring orientation; earcut and the
            // extruder both rely on outer CCW / holes CW, so restore it here
            // rather than letting a silently inside-out solid reach them.
            if transform.matrix2.determinant() < 0.0 {
                rings.outer.reverse();
                for hole in &mut rings.holes {
                    hole.reverse();
                }
            }
            Ok(rings)
        }
        Profile::CenterLine(cl) => {
            crate::center_line::center_line_rings(cl, chord_error, tolerance, contour_points)
        }
        other => Err(GeomError::Unsupported {
            backend: crate::BACKEND_ID,
            operation: axiolid_contracts::Operation::ProfileTriangulation,
        })
        .inspect_err(|_| {
            let _ = other;
        }),
    }
}

/// An exact contour with holes, flattened under the chord budget.
fn contour_rings(
    c: &axiolid_profile::ContourProfile,
    chord_error: Scalar,
    tolerance: Tolerance,
) -> GeomResult<Rings> {
    let outer = contour_points(&c.outer, chord_error, tolerance)?;
    let mut holes = Vec::with_capacity(c.holes.len());
    for hole in &c.holes {
        holes.push(contour_points(hole, chord_error, tolerance)?);
    }
    Ok(orient_rings(outer, holes))
}

/// Rectangle, optionally hollow, with sharp corners (rounded ones go
/// through the exact contour).
fn rectangle_rings(
    r: &RectangleProfile,
    _chord_error: Scalar,
    tolerance: Tolerance,
) -> GeomResult<Rings> {
    if !(r.x.is_finite() && r.y.is_finite()) || r.x <= 0.0 || r.y <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "rectangle profile must have positive finite extents, got {} x {}",
            r.x, r.y
        )));
    }
    let (hx, hy) = (r.x / 2.0, r.y / 2.0);
    let outer = vec![
        Point2::new(-hx, -hy),
        Point2::new(hx, -hy),
        Point2::new(hx, hy),
        Point2::new(-hx, hy),
    ];
    let mut holes = Vec::new();
    if let Some(t) = r.thickness {
        if t <= 0.0 || 2.0 * t >= r.x || 2.0 * t >= r.y {
            return Err(GeomError::InvalidInput(format!(
                "hollow rectangle wall thickness {t} does not fit inside {} x {}",
                r.x, r.y
            )));
        }
        let (ix, iy) = (hx - t, hy - t);
        if !tolerance.eq(ix, 0.0) && !tolerance.eq(iy, 0.0) {
            // Clockwise: opposite winding to the outer ring marks it a hole.
            holes.push(vec![
                Point2::new(-ix, -iy),
                Point2::new(-ix, iy),
                Point2::new(ix, iy),
                Point2::new(ix, -iy),
            ]);
        }
    }
    Ok(Rings { outer, holes })
}

/// Circle, optionally annular.
fn circle_rings(c: &CircleProfile, chord_error: Scalar) -> GeomResult<Rings> {
    if !c.radius.is_finite() || c.radius <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "circle profile radius must be positive and finite, got {}",
            c.radius
        )));
    }
    // Flattened by the scalar evaluator so a parameterized circle and a
    // contour-declared circular arc obey exactly the same chord budget.
    let outer = flatten_circle(c.radius, chord_error)?;
    let mut holes = Vec::new();
    if let Some(t) = c.thickness {
        if t <= 0.0 || t >= c.radius {
            return Err(GeomError::InvalidInput(format!(
                "annulus wall thickness {t} does not fit inside radius {}",
                c.radius
            )));
        }
        let inner = c.radius - t;
        let mut ring = flatten_circle(inner, chord_error)?;
        ring.reverse(); // clockwise
        holes.push(ring);
    }
    Ok(Rings { outer, holes })
}

/// Ellipse profile, flattened under the chord budget.
fn ellipse_rings(e: &EllipseProfile, chord_error: Scalar) -> GeomResult<Rings> {
    if !e.semi_axis_x.is_finite() || e.semi_axis_x <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "ellipse semi-axis x must be positive and finite, got {}",
            e.semi_axis_x
        )));
    }
    if !e.semi_axis_y.is_finite() || e.semi_axis_y <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "ellipse semi-axis y must be positive and finite, got {}",
            e.semi_axis_y
        )));
    }
    use axiolid_core::{Frame2, Interval, Vec2};
    use axiolid_curve::{Curve2, Ellipse2};

    let curve = Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::new(0.0, 0.0),
            x: Vec2::new(1.0, 0.0),
            y: Vec2::new(0.0, 1.0),
        },
        semi_axis_x: e.semi_axis_x,
        semi_axis_y: e.semi_axis_y,
    });
    let mut ring = axiolid_reference::curve::flatten2(
        &curve,
        Interval {
            start: 0.0,
            end: core::f64::consts::TAU,
        },
        chord_error,
        MAX_SUBDIVISION_DEPTH,
    )?;
    // Drop the duplicate closing vertex; a ring is implicitly closed.
    ring.pop();
    Ok(Rings {
        outer: ring,
        holes: Vec::new(),
    })
}

/// Flatten a full circle of `radius` under the chord budget.
///
/// The closing duplicate of the start point is dropped: a ring is implicitly
/// closed, and a repeated vertex would create a zero-length edge that the
/// extruder would turn into a degenerate side quad.
fn flatten_circle(radius: Scalar, chord_error: Scalar) -> GeomResult<Vec<Point2>> {
    use axiolid_core::{Frame2, Interval, Vec2};
    use axiolid_curve::{Circle2, Curve2};

    let curve = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(0.0, 0.0),
            x: Vec2::new(1.0, 0.0),
            y: Vec2::new(0.0, 1.0),
        },
        radius,
    });
    let flatten = |start: Scalar| {
        axiolid_reference::curve::flatten2(
            &curve,
            Interval {
                start,
                end: start + core::f64::consts::TAU,
            },
            chord_error,
            MAX_SUBDIVISION_DEPTH,
        )
    };
    // Bisection puts a point at every quarter turn. Started half a step
    // later, the same chords put a chord's middle there instead, strictly
    // inside the circle: a void tangent to a face along an axis direction,
    // as openings are, then leaves a sliver of material rather than a
    // chord point on the face, which would pinch the solid there (#194).
    let steps = flatten(0.0)?.len().saturating_sub(1).max(1);
    let mut ring = flatten(core::f64::consts::PI / steps as Scalar)?;
    ring.pop();
    Ok(ring)
}

/// Triangulate rings into a flat index buffer over a single vertex list.
///
/// Delegates to `earcut`. The returned vertices are the concatenation
/// `outer ++ holes`, matching earcut's hole-index convention.
pub fn triangulate(rings: &Rings) -> GeomResult<(Vec<Point2>, Vec<[u32; 3]>)> {
    if rings.outer.len() < 3 {
        return Err(GeomError::InvalidInput(format!(
            "profile outer ring needs at least 3 vertices, got {}",
            rings.outer.len()
        )));
    }
    let mut verts: Vec<[Scalar; 2]> = rings.outer.iter().map(|p| [p.x, p.y]).collect();
    let mut hole_starts = Vec::with_capacity(rings.holes.len());
    for hole in &rings.holes {
        if hole.len() < 3 {
            return Err(GeomError::InvalidInput(format!(
                "profile hole needs at least 3 vertices, got {}",
                hole.len()
            )));
        }
        hole_starts.push(verts.len());
        verts.extend(hole.iter().map(|p| [p.x, p.y]));
    }

    let mut earcutter = earcut::Earcut::new();
    let mut flat: Vec<usize> = Vec::new();
    earcutter.earcut(verts.iter().copied(), &hole_starts, &mut flat);

    if flat.is_empty() || flat.len() % 3 != 0 {
        return Err(GeomError::Degenerate(format!(
            "triangulation produced {} indices for a {}-vertex profile",
            flat.len(),
            verts.len()
        )));
    }
    let tris = flat
        .chunks_exact(3)
        .map(|c| [c[0] as u32, c[1] as u32, c[2] as u32])
        .collect();
    let points = verts.into_iter().map(|v| Point2::new(v[0], v[1])).collect();
    Ok((points, tris))
}

/// Apply a 2D affine transform to a ring in place.
fn apply2(ring: &mut [Point2], t: &axiolid_core::Transform2) {
    for p in ring.iter_mut() {
        *p = t.transform_point2(*p);
    }
}

/// Flatten one closed contour into a point ring.
///
/// Consecutive duplicate points are dropped: adjoining segments share an
/// endpoint by construction, and earcut treats a repeated vertex as a
/// zero-length edge.
fn contour_points(
    contour: &axiolid_profile::Contour,
    chord_error: Scalar,
    tolerance: Tolerance,
) -> GeomResult<Vec<Point2>> {
    let mut out: Vec<Point2> = Vec::new();
    for segment in &contour.segments {
        let mut pts = segment_points(segment, chord_error)?;
        if !segment.same_sense {
            pts.reverse();
        }
        for p in pts {
            if out
                .last()
                .is_none_or(|last| !near2(*last, p, tolerance.linear()))
            {
                out.push(p);
            }
        }
    }
    // A closed ring must not repeat its first point as its last.
    while out.len() > 1 && near2(out[0], *out.last().expect("non-empty"), tolerance.linear()) {
        out.pop();
    }
    if out.len() < 3 {
        return Err(GeomError::Degenerate(format!(
            "contour flattened to {} points, need at least 3",
            out.len()
        )));
    }
    Ok(out)
}

/// Whether two points coincide within a linear tolerance.
fn near2(a: Point2, b: Point2, linear: Scalar) -> bool {
    (a.x - b.x).abs() <= linear && (a.y - b.y).abs() <= linear
}

/// Sample one bounded segment, flattening curves under the chord budget.
///
/// Delegates to `axiolid-reference`'s curve evaluator (ADR 0012). This crate used
/// to carry a private circle flattener with a closed-form segment count, and
/// refused ellipses and B-splines outright. Both limits are gone: the scalar
/// evaluator subdivides adaptively on measured sagitta, so every declared
/// `Curve2` family flattens under the same tolerance contract.
///
/// `MAX_SUBDIVISION_DEPTH` bounds the work. 24 levels is 16M potential
/// segments -- far beyond any real tolerance -- so it is a runaway guard, not
/// a quality knob.
const MAX_SUBDIVISION_DEPTH: u32 = 24;

fn segment_points(
    segment: &axiolid_profile::ProfileSegment,
    chord_error: Scalar,
) -> GeomResult<Vec<Point2>> {
    axiolid_reference::curve::flatten2(
        &segment.curve,
        segment.domain,
        chord_error,
        MAX_SUBDIVISION_DEPTH,
    )
}

/// Force the ring-orientation convention earcut and the extruder expect:
/// outer counter-clockwise, holes clockwise.
///
/// Source contours carry whatever orientation the authoring tool wrote, so
/// normalising here is cheaper than rejecting otherwise-valid geometry.
use axiolid_reference::signed_area2;

fn orient_rings(mut outer: Vec<Point2>, mut holes: Vec<Vec<Point2>>) -> Rings {
    if signed_area2(&outer) < 0.0 {
        outer.reverse();
    }
    for hole in &mut holes {
        if signed_area2(hole) > 0.0 {
            hole.reverse();
        }
    }
    Rings { outer, holes }
}
