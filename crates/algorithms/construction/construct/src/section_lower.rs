//! Lower parameterised structural sections into exact contours (ADR 0057).
//!
//! # Why the fillets are not optional
//!
//! For a rolled steel section the web-to-flange root fillet is real material.
//! Measured on a 0.4 x 0.3 I-section with an 0.021 root radius, the four
//! fillets carry **2.40%** of the cross-sectional area. Dropping them yields a
//! section that looks correct and whose area, second moment and mass are all
//! wrong, so they are built as exact arcs rather than ignored or approximated.
//!
//! # Structure
//!
//! Every variant reduces to a closed counter-clockwise ring of corners, each
//! carrying an optional radius. One shared router turns that into a contour,
//! inserting a tangent arc at each rounded corner. Concave root fillets and
//! convex toe radii are the SAME operation -- the sign of the turn decides
//! which way the arc bends -- so neither gets a special case that could drift
//! from the other.

use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Interval, Point2, Scalar, Vec2};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, ProfileSegment, RectangleProfile, SectionProfile,
};

use crate::BACKEND_ID;

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BACKEND_ID,
        operation: Operation::Sweep,
        input,
    }
}

/// One boundary corner: a position and the radius rounding it.
#[derive(Debug, Clone, Copy)]
struct Corner {
    point: Point2,
    /// Zero means a sharp corner.
    radius: Scalar,
}

fn sharp(x: Scalar, y: Scalar) -> Corner {
    Corner {
        point: Point2::new(x, y),
        radius: 0.0,
    }
}

fn rounded(x: Scalar, y: Scalar, radius: Option<Scalar>) -> Corner {
    Corner {
        point: Point2::new(x, y),
        // `None` means the source did not state a radius, which is a sharp
        // corner here; `Some(0.0)` states one explicitly and agrees.
        radius: radius.unwrap_or(0.0).max(0.0),
    }
}

/// Turn a closed counter-clockwise corner ring into an exact contour.
///
/// At each corner with a positive radius the boundary is cut back along both
/// adjacent edges by `r * tan(alpha / 2)`, where `alpha` is the turn angle,
/// and joined by an arc tangent to both. That setback is what makes the arc
/// tangent rather than merely near the corner.
fn route(corners: &[Corner]) -> GeomResult<Contour> {
    let count = corners.len();
    if count < 3 {
        return Err(GeomError::InvalidInput(format!(
            "a section outline needs at least three corners, got {count}"
        )));
    }

    // Setback and arc geometry per corner, or `None` when sharp.
    let mut cut = vec![0.0; count];
    let mut arcs: Vec<Option<(Point2, Point2, Point2, Scalar)>> = vec![None; count];

    for index in 0..count {
        let here = corners[index];
        if here.radius <= 0.0 {
            continue;
        }
        let previous = corners[(index + count - 1) % count].point;
        let next = corners[(index + 1) % count].point;
        let incoming = (here.point - previous).normalize_or_zero();
        let outgoing = (next - here.point).normalize_or_zero();
        if incoming == Vec2::ZERO || outgoing == Vec2::ZERO {
            return Err(GeomError::Degenerate(
                "section outline has a zero-length edge".to_owned(),
            ));
        }
        let cross = incoming.perp_dot(outgoing);
        if cross == 0.0 {
            // Collinear: there is no corner to round.
            continue;
        }
        let turn = incoming.dot(outgoing).clamp(-1.0, 1.0).acos();
        let setback = here.radius * (turn / 2.0).tan();
        // Distance from the corner to the arc centre along the interior
        // bisector.
        //
        // `turn` is the EXTERIOR deflection, so the interior angle is
        // `pi - turn` and the centre sits `r / sin(interior / 2)` away, which
        // is `r / cos(turn / 2)`. Using `sin(turn / 2)` here is wrong
        // everywhere EXCEPT at a right angle, where the two coincide -- and
        // every rounded corner in every section variant is a right angle, so
        // the error is invisible to them.
        let bisector = (outgoing - incoming).normalize_or_zero();
        if bisector == Vec2::ZERO {
            return Err(GeomError::Degenerate(
                "section outline reverses on itself".to_owned(),
            ));
        }
        let centre = here.point + bisector * (here.radius / (turn / 2.0).cos());
        let start = here.point - incoming * setback;
        let end = here.point + outgoing * setback;
        cut[index] = setback;
        arcs[index] = Some((centre, start, end, cross.signum() * turn));
    }

    // Both ends of an edge draw from the same edge length.
    for start in 0..count {
        let end = (start + 1) % count;
        let length = (corners[end].point - corners[start].point).length();
        if cut[start] + cut[end] > length + 1e-12 {
            return Err(unsupported(
                "section radii too large for the edge between two corners",
            ));
        }
    }

    Ok(Contour::new(emit(corners, &arcs)))
}
/// Emit segments: a straight run between consecutive corners, plus an arc at
/// each rounded corner.
fn emit(
    corners: &[Corner],
    arcs: &[Option<(Point2, Point2, Point2, Scalar)>],
) -> Vec<ProfileSegment> {
    let count = corners.len();
    let mut segments = Vec::with_capacity(count * 2);
    for index in 0..count {
        // Where this corner hands over to the straight run that follows.
        let leave = match arcs[index] {
            Some((centre, start, end, sweep)) => {
                segments.push(arc_segment(centre, start, sweep));
                let _ = end;
                end
            }
            None => corners[index].point,
        };
        let next = (index + 1) % count;
        let arrive = match arcs[next] {
            Some((_, start, _, _)) => start,
            None => corners[next].point,
        };
        // A fully consumed edge leaves the two arcs touching; emitting a
        // zero-length line there would be a degenerate segment.
        if (arrive - leave).length() > 1e-15 {
            segments.push(ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: leave,
                    direction: arrive - leave,
                }),
                domain: Interval::UNIT,
                same_sense: true,
            });
        }
    }
    segments
}

/// A tangent arc from `start`, about `centre`, turning by `sweep`.
///
/// The frame's x-axis points at the arc start, so the segment domain begins
/// at zero. A negative sweep is carried by a LEFT-handed frame rather than a
/// negative domain, matching the convention the contour lowering already
/// uses: the parameter always increases, and handedness says which way the
/// world turn goes.
fn arc_segment(centre: Point2, start: Point2, sweep: Scalar) -> ProfileSegment {
    let x = (start - centre).normalize_or_zero();
    let perpendicular = Vec2::new(-x.y, x.x);
    let y = if sweep >= 0.0 {
        perpendicular
    } else {
        -perpendicular
    };
    ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: centre,
                x,
                y,
            },
            radius: (start - centre).length(),
        }),
        domain: Interval::new(0.0, sweep.abs()),
        same_sense: true,
    }
}

/// Reject a stated taper.
///
/// Validate a declared taper and return it as an angle.
///
/// A slope is an angle from the horizontal, so it must stay well inside a
/// quarter turn: at a right angle the inner face would be parallel to the web
/// and the section would have no flange at all. The bound is deliberately
/// generous -- rolled sections taper by 5 to 14 degrees -- because the job
/// here is to exclude nonsense, not to second-guess a source that states an
/// unusual but buildable value.
///
/// `None` means the source did not state a taper, which is a parallel flange.
/// `Some(0.0)` states one explicitly and gives the same geometry.
fn checked_slope(slope: Option<Scalar>, what: &'static str) -> GeomResult<Scalar> {
    let value = slope.unwrap_or(0.0);
    if !value.is_finite() {
        return Err(GeomError::InvalidInput(format!(
            "{what} slope must be finite, got {value}"
        )));
    }
    // A quarter turn is the hard limit; stop short of it so the tangent stays
    // usable rather than exploding.
    let limit = core::f64::consts::FRAC_PI_2 * 0.9;
    if value.abs() >= limit {
        return Err(GeomError::Degenerate(format!(
            "{what} slope {value} rad is too steep to leave a flange"
        )));
    }
    Ok(value)
}

/// A circle profile, and its bore when it has a wall thickness, as exact
/// contours of four quarter arcs each (#111).
///
/// Four quarters, not one full turn: contour lowering refuses a segment of
/// half a turn or more (ADR 0053), and a quarter keeps every arc's endpoints
/// on the axes, where the seams of revolved and extruded walls sit.
pub fn circle_contour(circle: &CircleProfile) -> GeomResult<ContourProfile> {
    positive(circle.radius, "circle radius")?;
    let ring = |radius: Scalar| {
        let frame = Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        };
        let quarter = core::f64::consts::FRAC_PI_2;
        Contour::new(
            (0..4)
                .map(|index| ProfileSegment {
                    curve: Curve2::Circle(Circle2 { frame, radius }),
                    domain: Interval::new(
                        quarter * index as Scalar,
                        quarter * (index + 1) as Scalar,
                    ),
                    same_sense: true,
                })
                .collect(),
        )
    };
    let holes = match circle.thickness {
        None => Vec::new(),
        Some(thickness) => {
            positive(thickness, "circle wall thickness")?;
            if thickness >= circle.radius {
                return Err(GeomError::InvalidInput(format!(
                    "circle wall thickness {thickness} leaves no bore in radius {}",
                    circle.radius
                )));
            }
            vec![ring(circle.radius - thickness)]
        }
    };
    Ok(ContourProfile {
        outer: ring(circle.radius),
        holes,
    })
}

/// Lower a rectangle -- rounded, hollow, or both -- into an exact contour.
///
/// Corner radii become exact quarter arcs through the same router the
/// structural sections use, so a rounded corner extrudes to a cylinder wall
/// rather than a fan of chords. A hollow rectangle's inner boundary comes
/// back as a hole, built counter-clockwise like the outer ring; the extruder
/// re-orients holes itself.
///
/// Refuses, rather than clamps, a radius that is negative, non-finite, or
/// wider than the half-extent it rounds, and an inner radius on a filled
/// rectangle.
///
/// Also refuses a hollow section whose corners leave no wall. Two rounded
/// rectangles nest exactly when their support functions do; along the corner
/// diagonal that requires `outer - inner < (2 + sqrt 2) * thickness`. Past
/// that the inner corner reaches the outer arc and the profile crosses itself.
pub fn rectangle_contour(rectangle: &RectangleProfile) -> GeomResult<ContourProfile> {
    positive(rectangle.x, "rectangle x extent")?;
    positive(rectangle.y, "rectangle y extent")?;
    let (hx, hy) = (rectangle.x / 2.0, rectangle.y / 2.0);
    let outer_radius = corner_radius(rectangle.outer_radius, hx.min(hy), "outer")?;
    let outer = route(&box_corners(hx, hy, outer_radius))?;
    let Some(thickness) = rectangle.thickness else {
        if rectangle.inner_radius.is_some() {
            return Err(GeomError::InvalidInput(
                "an inner corner radius needs a hollow rectangle".to_owned(),
            ));
        }
        return Ok(ContourProfile {
            outer,
            holes: Vec::new(),
        });
    };
    positive(thickness, "rectangle wall thickness")?;
    if 2.0 * thickness >= rectangle.x.min(rectangle.y) {
        return Err(GeomError::Degenerate(format!(
            "wall thickness {thickness} leaves no opening in {} x {}",
            rectangle.x, rectangle.y
        )));
    }
    let (ix, iy) = (hx - thickness, hy - thickness);
    let inner_radius = corner_radius(rectangle.inner_radius, ix.min(iy), "inner")?;
    let wall_limit = (2.0 + core::f64::consts::SQRT_2) * thickness;
    if outer_radius - inner_radius >= wall_limit {
        return Err(GeomError::Degenerate(format!(
            "outer radius {outer_radius} minus inner radius {inner_radius} must stay \
             below {wall_limit} or the {thickness}-thick wall vanishes at the corners"
        )));
    }
    let hole = route(&box_corners(ix, iy, inner_radius))?;
    Ok(ContourProfile {
        outer,
        holes: vec![hole],
    })
}

/// Counter-clockwise corners of an origin-centred box, all rounded alike.
fn box_corners(hx: Scalar, hy: Scalar, radius: Scalar) -> [Corner; 4] {
    let radius = Some(radius);
    [
        rounded(-hx, -hy, radius),
        rounded(hx, -hy, radius),
        rounded(hx, hy, radius),
        rounded(-hx, hy, radius),
    ]
}

/// Validate a stated corner radius; `None` is a sharp corner.
fn corner_radius(radius: Option<Scalar>, limit: Scalar, which: &str) -> GeomResult<Scalar> {
    let value = radius.unwrap_or(0.0);
    if !value.is_finite() || value < 0.0 || value > limit {
        return Err(GeomError::InvalidInput(format!(
            "{which} corner radius must lie in [0, {limit}], got {value}"
        )));
    }
    Ok(value)
}

fn positive(value: Scalar, what: &str) -> GeomResult<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "section {what} must be positive and finite, got {value}"
        )));
    }
    Ok(())
}
/// Lower a parameterised section into an exact contour.
///
/// The outline is built counter-clockwise about the section centroid-ish
/// origin each source entity declares, so the result needs no re-orientation.
pub fn section_contour(section: &SectionProfile) -> GeomResult<ContourProfile> {
    let corners = match section {
        SectionProfile::I {
            depth,
            width,
            web_thickness,
            flange_thickness,
            fillet_radius,
            flange_edge_radius,
            flange_slope,
        } => {
            let slope = checked_slope(*flange_slope, "I section")?;
            i_corners(
                *depth,
                *width,
                *width,
                *web_thickness,
                *flange_thickness,
                *flange_thickness,
                *fillet_radius,
                *fillet_radius,
                *flange_edge_radius,
                *flange_edge_radius,
                slope,
                slope,
            )?
        }
        SectionProfile::AsymmetricI {
            depth,
            web_thickness,
            bottom_flange_width,
            bottom_flange_thickness,
            bottom_fillet_radius,
            bottom_flange_edge_radius,
            bottom_flange_slope,
            top_flange_width,
            top_flange_thickness,
            top_fillet_radius,
            top_flange_edge_radius,
            top_flange_slope,
        } => {
            let bottom_slope = checked_slope(*bottom_flange_slope, "I section")?;
            let top_slope = checked_slope(*top_flange_slope, "I section")?;
            i_corners(
                *depth,
                *bottom_flange_width,
                *top_flange_width,
                *web_thickness,
                *bottom_flange_thickness,
                // The source may omit the top thickness, meaning "same as
                // the bottom" rather than "zero".
                top_flange_thickness.unwrap_or(*bottom_flange_thickness),
                *bottom_fillet_radius,
                *top_fillet_radius,
                *bottom_flange_edge_radius,
                *top_flange_edge_radius,
                bottom_slope,
                top_slope,
            )?
        }
        SectionProfile::T {
            depth,
            flange_width,
            web_thickness,
            flange_thickness,
            fillet_radius,
            flange_edge_radius,
            web_edge_radius,
            web_slope,
            flange_slope,
        } => {
            // A T's web taper and flange taper are independent faces.
            let web = checked_slope(*web_slope, "T section web")?;
            let flange = checked_slope(*flange_slope, "T section flange")?;
            t_corners(
                *depth,
                *flange_width,
                *web_thickness,
                *flange_thickness,
                &RadiiT {
                    fillet: *fillet_radius,
                    flange_edge: *flange_edge_radius,
                    web_edge: *web_edge_radius,
                },
                &TaperT { web, flange },
            )?
        }
        SectionProfile::U {
            depth,
            flange_width,
            web_thickness,
            flange_thickness,
            fillet_radius,
            edge_radius,
            flange_slope,
        } => {
            let slope = checked_slope(*flange_slope, "U section")?;
            u_corners(
                *depth,
                *flange_width,
                *web_thickness,
                *flange_thickness,
                *fillet_radius,
                *edge_radius,
                slope,
            )?
        }
        SectionProfile::L {
            depth,
            width,
            thickness,
            fillet_radius,
            edge_radius,
            leg_slope,
        } => {
            let slope = checked_slope(*leg_slope, "L section")?;
            l_corners(
                *depth,
                // An absent width means an EQUAL angle, not a zero one.
                width.unwrap_or(*depth),
                *thickness,
                *fillet_radius,
                *edge_radius,
                slope,
            )?
        }
        SectionProfile::Z {
            depth,
            flange_width,
            web_thickness,
            flange_thickness,
            fillet_radius,
            edge_radius,
        } => z_corners(
            *depth,
            *flange_width,
            *web_thickness,
            *flange_thickness,
            *fillet_radius,
            *edge_radius,
        )?,
        SectionProfile::C {
            depth,
            width,
            wall_thickness,
            girth,
            internal_fillet_radius,
        } => c_corners(
            *depth,
            *width,
            *wall_thickness,
            *girth,
            *internal_fillet_radius,
        )?,
        SectionProfile::Trapezium {
            bottom_x,
            top_x,
            y,
            top_offset,
        } => trapezium_corners(*bottom_x, *top_x, *y, *top_offset)?,
        _ => return Err(unsupported("section profile of an unsupported kind")),
    };

    Ok(ContourProfile {
        outer: route(&corners)?,
        holes: Vec::new(),
    })
}
/// Twelve-corner I outline, walked counter-clockwise from the bottom-right.
///
/// Handles the asymmetric case directly; the symmetric variant passes equal
/// top and bottom dimensions rather than going through a separate routine
/// that could drift from this one.
#[allow(clippy::too_many_arguments)]
fn i_corners(
    depth: Scalar,
    bottom_width: Scalar,
    top_width: Scalar,
    web_thickness: Scalar,
    bottom_flange: Scalar,
    top_flange: Scalar,
    bottom_fillet: Option<Scalar>,
    top_fillet: Option<Scalar>,
    bottom_edge: Option<Scalar>,
    top_edge: Option<Scalar>,
    bottom_slope: Scalar,
    top_slope: Scalar,
) -> GeomResult<Vec<Corner>> {
    positive(depth, "depth")?;
    positive(bottom_width, "flange width")?;
    positive(top_width, "flange width")?;
    positive(web_thickness, "web thickness")?;
    positive(bottom_flange, "flange thickness")?;
    positive(top_flange, "flange thickness")?;
    if bottom_flange + top_flange >= depth {
        return Err(GeomError::Degenerate(format!(
            "flanges {bottom_flange} + {top_flange} leave no web in depth {depth}"
        )));
    }
    if web_thickness >= bottom_width.min(top_width) {
        return Err(GeomError::Degenerate(format!(
            "web thickness {web_thickness} is not narrower than the flange"
        )));
    }

    let (hd, hw) = (depth / 2.0, web_thickness / 2.0);
    let (hb, ht) = (bottom_width / 2.0, top_width / 2.0);
    let bottom_top = -hd + bottom_flange;
    let top_bottom = hd - top_flange;

    // A tapered flange's inner face is inclined, so its height depends on x.
    // The face pivots about the MID-POINT between the web face and the flange
    // tip, which keeps `flange_thickness` the MEAN thickness -- the value
    // section tables state. Pivoting about the tip or the web instead would
    // silently change the declared thickness and with it the area.
    let bottom_rise = |x: Scalar| (x - (hw + hb) / 2.0) * bottom_slope.tan();
    let top_rise = |x: Scalar| (x - (hw + ht) / 2.0) * top_slope.tan();

    Ok(vec![
        sharp(hb, -hd),
        rounded(hb, bottom_top - bottom_rise(hb), bottom_edge),
        rounded(hw, bottom_top - bottom_rise(hw), bottom_fillet),
        rounded(hw, top_bottom + top_rise(hw), top_fillet),
        rounded(ht, top_bottom + top_rise(ht), top_edge),
        sharp(ht, hd),
        sharp(-ht, hd),
        rounded(-ht, top_bottom + top_rise(ht), top_edge),
        rounded(-hw, top_bottom + top_rise(hw), top_fillet),
        rounded(-hw, bottom_top - bottom_rise(hw), bottom_fillet),
        rounded(-hb, bottom_top - bottom_rise(hb), bottom_edge),
        sharp(-hb, -hd),
    ])
}

/// Eight-corner T outline: flange on top, web hanging below.
/// The two independent tapers a T section can declare.
///
/// Grouped rather than passed loose so the web angle cannot be handed to the
/// flange by accident: the two are the same type and adjacent in the argument
/// list, which is exactly the shape of a silent swap.
/// The three optional radii a T section can declare.
struct RadiiT {
    fillet: Option<Scalar>,
    flange_edge: Option<Scalar>,
    web_edge: Option<Scalar>,
}

struct TaperT {
    web: Scalar,
    flange: Scalar,
}

fn t_corners(
    depth: Scalar,
    flange_width: Scalar,
    web_thickness: Scalar,
    flange_thickness: Scalar,
    radii: &RadiiT,
    taper: &TaperT,
) -> GeomResult<Vec<Corner>> {
    let (web_slope, flange_slope) = (taper.web, taper.flange);
    let (fillet, flange_edge, web_edge) = (radii.fillet, radii.flange_edge, radii.web_edge);
    positive(depth, "depth")?;
    positive(flange_width, "flange width")?;
    positive(web_thickness, "web thickness")?;
    positive(flange_thickness, "flange thickness")?;
    if flange_thickness >= depth {
        return Err(GeomError::Degenerate(format!(
            "flange {flange_thickness} leaves no web in depth {depth}"
        )));
    }
    if web_thickness >= flange_width {
        return Err(GeomError::Degenerate(format!(
            "web thickness {web_thickness} is not narrower than the flange"
        )));
    }

    let (hd, hw, hf) = (depth / 2.0, web_thickness / 2.0, flange_width / 2.0);
    let flange_bottom = hd - flange_thickness;

    // Flange underside inclines about the mid-point between web face and tip,
    // keeping `flange_thickness` the mean. The web's side faces incline about
    // the mid-height of the exposed web run for the same reason.
    let flange_rise = |x: Scalar| (x - (hw + hf) / 2.0) * flange_slope.tan();
    let web_mid = (-hd + flange_bottom) / 2.0;
    let web_out = |y: Scalar| (y - web_mid) * web_slope.tan();

    Ok(vec![
        rounded(hw + web_out(-hd), -hd, web_edge),
        rounded(hw + web_out(flange_bottom), flange_bottom, fillet),
        rounded(hf, flange_bottom - flange_rise(hf), flange_edge),
        sharp(hf, hd),
        sharp(-hf, hd),
        rounded(-hf, flange_bottom - flange_rise(hf), flange_edge),
        rounded(-hw - web_out(flange_bottom), flange_bottom, fillet),
        rounded(-hw - web_out(-hd), -hd, web_edge),
    ])
}

/// Eight-corner U (channel) outline: web on the left, flanges to the right.
fn u_corners(
    depth: Scalar,
    flange_width: Scalar,
    web_thickness: Scalar,
    flange_thickness: Scalar,
    fillet: Option<Scalar>,
    edge: Option<Scalar>,
    flange_slope: Scalar,
) -> GeomResult<Vec<Corner>> {
    positive(depth, "depth")?;
    positive(flange_width, "flange width")?;
    positive(web_thickness, "web thickness")?;
    positive(flange_thickness, "flange thickness")?;
    if 2.0 * flange_thickness >= depth {
        return Err(GeomError::Degenerate(format!(
            "flanges {flange_thickness} leave no web in depth {depth}"
        )));
    }
    if web_thickness >= flange_width {
        return Err(GeomError::Degenerate(format!(
            "web thickness {web_thickness} is not narrower than the flange"
        )));
    }

    let hd = depth / 2.0;
    let inner = web_thickness;
    // Same convention as the I: pivot about the mid-point of the inner face so
    // the stated flange thickness stays the mean.
    let rise = |x: Scalar| (x - (inner + flange_width) / 2.0) * flange_slope.tan();

    Ok(vec![
        sharp(flange_width, -hd),
        rounded(
            flange_width,
            -hd + flange_thickness - rise(flange_width),
            edge,
        ),
        rounded(inner, -hd + flange_thickness - rise(inner), fillet),
        rounded(inner, hd - flange_thickness + rise(inner), fillet),
        rounded(
            flange_width,
            hd - flange_thickness + rise(flange_width),
            edge,
        ),
        sharp(flange_width, hd),
        sharp(0.0, hd),
        sharp(0.0, -hd),
    ])
}

/// Six-corner L (angle) outline with the heel at the origin.
fn l_corners(
    depth: Scalar,
    width: Scalar,
    thickness: Scalar,
    fillet: Option<Scalar>,
    edge: Option<Scalar>,
    leg_slope: Scalar,
) -> GeomResult<Vec<Corner>> {
    positive(depth, "depth")?;
    positive(width, "width")?;
    positive(thickness, "thickness")?;
    if thickness >= depth.min(width) {
        return Err(GeomError::Degenerate(format!(
            "thickness {thickness} is not thinner than the legs"
        )));
    }

    // Each leg's inner face inclines about the mid-point of its run, so the
    // stated thickness remains the mean thickness of the leg.
    let tan = leg_slope.tan();
    let horizontal = |x: Scalar| (x - (thickness + width) / 2.0) * tan;
    let vertical = |y: Scalar| (y - (thickness + depth) / 2.0) * tan;

    Ok(vec![
        sharp(0.0, 0.0),
        sharp(width, 0.0),
        rounded(width, thickness - horizontal(width), edge),
        rounded(thickness, thickness, fillet),
        rounded(thickness - vertical(depth), depth, edge),
        sharp(0.0, depth),
    ])
}
/// Eight-corner Z outline: flanges point in opposite directions.
fn z_corners(
    depth: Scalar,
    flange_width: Scalar,
    web_thickness: Scalar,
    flange_thickness: Scalar,
    fillet: Option<Scalar>,
    edge: Option<Scalar>,
) -> GeomResult<Vec<Corner>> {
    positive(depth, "depth")?;
    positive(flange_width, "flange width")?;
    positive(web_thickness, "web thickness")?;
    positive(flange_thickness, "flange thickness")?;
    if 2.0 * flange_thickness >= depth {
        return Err(GeomError::Degenerate(format!(
            "flanges {flange_thickness} leave no web in depth {depth}"
        )));
    }

    let (hd, hw) = (depth / 2.0, web_thickness / 2.0);
    let bottom_top = -hd + flange_thickness;
    let top_bottom = hd - flange_thickness;

    // Bottom flange runs right, top flange runs left.
    Ok(vec![
        sharp(hw + flange_width, -hd),
        rounded(hw + flange_width, bottom_top, edge),
        rounded(hw, bottom_top, fillet),
        sharp(hw, hd),
        sharp(-hw - flange_width, hd),
        rounded(-hw - flange_width, top_bottom, edge),
        rounded(-hw, top_bottom, fillet),
        sharp(-hw, -hd),
    ])
}

/// Twelve-corner C (lipped channel) outline.
///
/// Unlike `U`, this is a THIN-WALLED section: the boundary follows the wall
/// all the way round, including the returned lips, so the enclosed area is
/// the wall material rather than the full channel envelope.
fn c_corners(
    depth: Scalar,
    width: Scalar,
    wall_thickness: Scalar,
    girth: Scalar,
    fillet: Option<Scalar>,
) -> GeomResult<Vec<Corner>> {
    positive(depth, "depth")?;
    positive(width, "width")?;
    positive(wall_thickness, "wall thickness")?;
    positive(girth, "girth")?;
    if 2.0 * wall_thickness >= depth || 2.0 * wall_thickness >= width {
        return Err(GeomError::Degenerate(format!(
            "wall thickness {wall_thickness} leaves no opening"
        )));
    }
    if girth <= wall_thickness {
        return Err(GeomError::Degenerate(format!(
            "lip girth {girth} is not longer than the wall thickness"
        )));
    }

    let hd = depth / 2.0;
    let t = wall_thickness;

    Ok(vec![
        // Outer boundary: up the web, out along each flange, down each lip.
        sharp(width, -hd),
        sharp(width, -hd + girth),
        sharp(width - t, -hd + girth),
        rounded(width - t, -hd + t, fillet),
        rounded(t, -hd + t, fillet),
        rounded(t, hd - t, fillet),
        rounded(width - t, hd - t, fillet),
        sharp(width - t, hd - girth),
        sharp(width, hd - girth),
        sharp(width, hd),
        sharp(0.0, hd),
        sharp(0.0, -hd),
    ])
}

/// Four-corner trapezium.
fn trapezium_corners(
    bottom_x: Scalar,
    top_x: Scalar,
    y: Scalar,
    top_offset: Scalar,
) -> GeomResult<Vec<Corner>> {
    positive(bottom_x, "bottom width")?;
    positive(top_x, "top width")?;
    positive(y, "height")?;
    if !top_offset.is_finite() {
        return Err(GeomError::InvalidInput(format!(
            "trapezium top offset must be finite, got {top_offset}"
        )));
    }

    Ok(vec![
        sharp(0.0, 0.0),
        sharp(bottom_x, 0.0),
        sharp(top_offset + top_x, y),
        sharp(top_offset, y),
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The corner router must round a NON-right corner correctly.
    ///
    /// Every rounded corner reachable through `SectionProfile` today is a
    /// right angle, and at 90 degrees `r * tan(45) == r`, so a setback that
    /// ignored the turn angle would be indistinguishable there. The router is
    /// general, so its generality is tested directly rather than left to a
    /// variant that happens not to exercise it.
    ///
    /// The invariant checked is TANGENCY: the arc must meet both adjacent
    /// edges at a point whose distance to the arc centre equals the radius,
    /// and the centre must sit exactly `radius` from each edge line.
    #[test]
    fn a_non_right_corner_is_rounded_tangentially() {
        // A 116.565-degree turn: setback is r*tan(58.28) = 1.618 r, not r.
        let radius = 0.02;
        let corners = vec![
            sharp(0.0, 0.0),
            Corner {
                point: Point2::new(0.4, 0.0),
                radius,
            },
            sharp(0.3, 0.2),
            sharp(0.05, 0.2),
        ];
        let contour = route(&corners).expect("a trapezoidal ring routes");

        let arc = contour
            .segments
            .iter()
            .find_map(|segment| match &segment.curve {
                Curve2::Circle(circle) => Some(*circle),
                _ => None,
            })
            .expect("the rounded corner produced an arc");
        assert!(
            (arc.radius - radius).abs() < 1e-12,
            "arc must carry the stated radius, got {}",
            arc.radius
        );

        // Distance from the arc centre to each adjacent edge LINE must equal
        // the radius. That is tangency, and it holds only for the correct
        // setback.
        let distance_to_line = |a: Point2, b: Point2| {
            let along = (b - a).normalize();
            let normal = Vec2::new(-along.y, along.x);
            (arc.frame.origin - a).dot(normal).abs()
        };
        let incoming = distance_to_line(Point2::new(0.0, 0.0), Point2::new(0.4, 0.0));
        let outgoing = distance_to_line(Point2::new(0.4, 0.0), Point2::new(0.3, 0.2));
        assert!(
            (incoming - radius).abs() < 1e-12,
            "arc must be tangent to the incoming edge, distance {incoming}"
        );
        assert!(
            (outgoing - radius).abs() < 1e-12,
            "arc must be tangent to the outgoing edge, distance {outgoing}"
        );
    }
}
