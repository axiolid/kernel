//! Exact B-reps of runs of sections at stations (#247, ADR 0082
//! amendment).
//!
//! Three runs have a closed form, and `crate::ReferenceExactCompiler`
//! builds them exactly:
//!
//! - **A constant section along a straight basis** (a line, or an open
//!   polyline whose vertices lie on one line and run forwards): every
//!   section is the first moved along the tangent, so the run is the
//!   profile extruded by the run's length along the tangent, in the
//!   section's plane (`axiolid_construct::extrude::extrude_profile_exact`,
//!   oblique where the orientation tilts the section). Its volume is `A
//!   L cos(tilt)`, its area `2 A + P L`-like.
//! - **A constant section along a horizontal circle** (a 2D circle, or a
//!   3D one in a horizontal plane), unturned and with no longitudinal
//!   offset: every section is the first turned about the circle's axis,
//!   which lies in the section's plane along its up, so the run is the
//!   profile revolved by `(s1 - s0) / r`
//!   (`axiolid_construct::revolve_exact::revolve_profile_exact`), with
//!   Pappus's volume `A theta R_c` and area.
//! - **A linear taper along a straight basis**: straight-edged sections
//!   (their rings are the exact profiles) in one constant orientation,
//!   joined vertex for vertex; every wall must be planar, and the run is
//!   the polyhedron of [`axiolid_construct::loft_exact`].
//!
//! The frames are read in closed form, never through the arc-length
//! inverse: on a line or a straight polyline by the line itself, on a
//! circle at the angle `s / r` from its frame's `x` (the circle's own
//! arc length, exact for an orthonormal frame). Everything else is
//! refused by name: another basis, a section that changes along a curve,
//! a curved section that changes, an orientation that turns between
//! stations, a warped wall, a turned or longitudinally offset section on
//! a circle.

use axiolid_brep::{ExactBRep, TransformError};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::loft_exact::loft_polygon_stations_exact;
use axiolid_construct::profile::{profile_deviation, ProfileDeviation};
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Frame3, Point2, Point3, Scalar, Transform3, Vec3};
use axiolid_curve::{Circle3, Curve2, Curve3, Line2, Line3};
use axiolid_model::{
    GeometryGraph, GeometryNode, NodeId, SeamSide, SectionAtStation, StationFrame,
};
use axiolid_reference::station::{station_section2, station_section3, SectionFrame};

use super::{oriented, spine_distances, spine_rings, Basis};

/// The closed-form bases.
enum Exact {
    /// A line through `origin` along the unit `direction`, measured from
    /// `origin`, in 2D or 3D.
    Straight(StraightLine),
    /// A circle with an orthonormal frame in a horizontal plane.
    Circle(Circle3),
}

enum StraightLine {
    Two(Curve2),
    Three(Curve3),
}

/// Why a basis has no exact run.
const OTHER_BASIS: &str = "an exact station-placed spine along a basis other than a line, a \
     straight polyline or a horizontal circle";

fn orthonormal(x: Vec3, y: Vec3) -> bool {
    let near = 1e-12;
    (x.length() - 1.0).abs() <= near && (y.length() - 1.0).abs() <= near && x.dot(y).abs() <= near
}

/// A polyline's vertices, if they lie on one line within `linear` and run
/// forwards along it: the line from the first, with a unit direction.
fn straight_polyline(points: &[Point3], linear: Scalar) -> Option<(Point3, Vec3)> {
    let (first, last) = (*points.first()?, *points.last()?);
    let direction = (last - first).normalize_or_zero();
    if direction == Vec3::ZERO {
        return None;
    }
    let mut along = 0.0;
    for p in points {
        let offset = *p - first;
        let t = offset.dot(direction);
        if (offset - direction * t).length() > linear || t < along - linear {
            return None;
        }
        along = t;
    }
    Some((first, direction))
}

fn exact_basis(graph: &GeometryGraph, directrix: NodeId, linear: Scalar) -> Option<Exact> {
    let lift = |p: Point2| Point3::new(p.x, p.y, 0.0);
    let line2 = |origin: Point3, direction: Vec3| {
        Exact::Straight(StraightLine::Two(Curve2::Line(Line2 {
            origin: Point2::new(origin.x, origin.y),
            direction: Point2::new(direction.x, direction.y),
        })))
    };
    let line3 = |origin: Point3, direction: Vec3| {
        Exact::Straight(StraightLine::Three(Curve3::Line(Line3 {
            origin,
            direction,
        })))
    };
    match graph.get(directrix)? {
        GeometryNode::Curve2(Curve2::Line(line)) => {
            let direction = Vec3::new(line.direction.x, line.direction.y, 0.0).normalize_or_zero();
            (direction != Vec3::ZERO).then(|| line2(lift(line.origin), direction))
        }
        GeometryNode::Curve3(Curve3::Line(line)) => {
            let direction = line.direction.normalize_or_zero();
            (direction != Vec3::ZERO).then(|| line3(line.origin, direction))
        }
        GeometryNode::Curve2(Curve2::Polyline(polyline)) if !polyline.closed => {
            let points: Vec<Point3> = polyline.points.iter().map(|p| lift(*p)).collect();
            straight_polyline(&points, linear).map(|(o, d)| line2(o, d))
        }
        GeometryNode::Curve3(Curve3::Polyline(polyline)) if !polyline.closed => {
            straight_polyline(&polyline.points, linear).map(|(o, d)| line3(o, d))
        }
        GeometryNode::Curve2(Curve2::Circle(c)) => {
            let (x, y) = (lift(c.frame.x), lift(c.frame.y));
            orthonormal(x, y).then(|| {
                Exact::Circle(Circle3 {
                    frame: Frame3 {
                        origin: lift(c.frame.origin),
                        x,
                        y,
                        z: x.cross(y),
                    },
                    radius: c.radius,
                })
            })
        }
        GeometryNode::Curve3(Curve3::Circle(c)) => {
            let normal = c.frame.x.cross(c.frame.y);
            (orthonormal(c.frame.x, c.frame.y) && normal.z.abs() >= 1.0 - 1e-12)
                .then_some(Exact::Circle(*c))
        }
        _ => None,
    }
}

impl Exact {
    /// The section frame at distance `s`, in closed form, in `frame`.
    fn section(&self, s: Scalar, frame: StationFrame) -> GeomResult<SectionFrame> {
        let section = match self {
            Self::Straight(StraightLine::Two(curve)) => station_section2(curve, s)?,
            Self::Straight(StraightLine::Three(curve)) => station_section3(curve, s)?,
            // A horizontal circle's reference-up frame is its plan frame:
            // up `+Z`, lateral the left normal.
            Self::Circle(circle) => {
                let (sin, cos) = (s / circle.radius).sin_cos();
                let point = circle.frame.origin
                    + circle.radius * (cos * circle.frame.x + sin * circle.frame.y);
                let tangent = (cos * circle.frame.y - sin * circle.frame.x).normalize();
                plan_frame(point, tangent)?
            }
        };
        match frame {
            StationFrame::Plan => section.plan(),
            _ => Ok(section),
        }
    }
}

/// The frame of a horizontal tangent at `point`: the reference-up frame
/// of the line along it, read where it starts (exact on a line).
fn plan_frame(point: Point3, tangent: Vec3) -> GeomResult<SectionFrame> {
    station_section3(
        &Curve3::Line(Line3 {
            origin: point,
            direction: tangent,
        }),
        0.0,
    )
}

/// The rigid motion taking local `x`, `y`, `z` onto `x`, `y`, `x x y` at
/// `origin`.
fn local(x: Vec3, y: Vec3, origin: Point3) -> Transform3 {
    Transform3::from_cols(x, y, x.cross(y), origin)
}

/// A run of sections at stations compiled exactly; see the
/// [module documentation](self). `unsupported` builds the exact
/// compiler's refusal.
pub(crate) fn stationed_spine_exact(
    graph: &GeometryGraph,
    directrix: NodeId,
    sections: &[SectionAtStation],
    frame: StationFrame,
    options: &ExecutionOptions,
    unsupported: fn(&'static str) -> GeomError,
) -> GeomResult<ExactBRep> {
    if !matches!(frame, StationFrame::Section | StationFrame::Plan) {
        return Err(unsupported("a station frame this compiler does not know"));
    }
    let tolerance = options.tolerance();
    let half = 0.5 * crate::compiler::chord_error(options);
    let (rings, _) = spine_rings(graph, sections, half, tolerance)?;
    let basis = Basis::of(graph, directrix)?;
    let (distances, _) = spine_distances(&basis, sections)?;
    // Every station must stand on the basis, as the mesh compiler reads it.
    for &s in &distances {
        basis.section_on(s, frame, SeamSide::Outgoing)?;
    }
    let exact = exact_basis(graph, directrix, tolerance.linear())
        .ok_or_else(|| unsupported(OTHER_BASIS))?;
    let first = &sections[0];
    let unit = |section: &SectionAtStation| {
        section
            .orientation
            .unit_axes()
            .map_err(|detail| GeomError::InvalidInput(format!("station: {detail}")))
    };
    let turn = unit(first)?;
    let mut constant_turn = true;
    for section in sections {
        constant_turn &= unit(section)? == turn;
    }
    if !constant_turn {
        return Err(unsupported(
            "an exact station-placed spine whose orientation turns between stations",
        ));
    }
    let profile = |section: &SectionAtStation| match graph.get(section.profile) {
        Some(GeometryNode::Profile(shape)) => Ok(shape),
        _ => Err(GeomError::InvalidInput(format!(
            "spine section profile {:?} is not a Profile node",
            section.profile
        ))),
    };
    let constant = sections.iter().all(|section| {
        (section.profile == first.profile || graph.get(section.profile) == graph.get(first.profile))
            && section.tags == first.tags
            && section.station.offsets == first.station.offsets
    });
    let (s0, s1) = (distances[0], distances[distances.len() - 1]);
    let base = exact.section(s0, frame)?;
    let turned = oriented(&base, &first.orientation, &first.orientation, 0.0)?;
    let offsets = first.station.offsets;
    let place_error = |error: TransformError| match error {
        TransformError::Unsupported(what) => unsupported(what),
        other => GeomError::InvalidInput(format!("exact station-placed spine: {other}")),
    };
    match (&exact, constant) {
        (Exact::Straight(_), true) => {
            let origin = base.place(offsets.lateral, offsets.vertical, offsets.longitudinal);
            let normal = turned.lateral.cross(turned.up);
            let direction = Vec3::new(
                base.tangent.dot(turned.lateral),
                base.tangent.dot(turned.up),
                base.tangent.dot(normal),
            );
            let prism = extrude_profile_exact(profile(first)?, direction, s1 - s0, tolerance)
                .map_err(remap(unsupported))?;
            prism
                .transformed(&local(turned.lateral, turned.up, origin))
                .map_err(place_error)
        }
        (Exact::Circle(circle), true) => {
            if !first.orientation.is_base() || offsets.longitudinal != 0.0 {
                return Err(unsupported(
                    "an exact station-placed spine along a circle with a turned section or a \
                     longitudinal offset: its section does not hold the circle's axis",
                ));
            }
            let origin = base.place(offsets.lateral, offsets.vertical, 0.0);
            let centre = circle.frame.origin - origin;
            let axis = circle.frame.x.cross(circle.frame.y);
            let along = axis.dot(base.up).signum();
            let prism = revolve_profile_exact(
                profile(first)?,
                Point3::new(centre.dot(base.lateral), 0.0, 0.0),
                Vec3::new(0.0, along, 0.0),
                (s1 - s0) / circle.radius,
                tolerance,
            )
            .map_err(remap(unsupported))?;
            prism
                .transformed(&local(base.lateral, base.up, origin))
                .map_err(place_error)
        }
        (Exact::Circle(_), false) => Err(unsupported(
            "an exact station-placed spine whose section changes along a circle: its walls are \
             not ruled by straight lines",
        )),
        (Exact::Straight(_), false) => {
            for section in sections {
                if !matches!(
                    profile_deviation(profile(section)?, half, tolerance)?,
                    ProfileDeviation::Bounded(d) if d == 0.0
                ) {
                    return Err(unsupported(
                        "an exact station-placed spine whose curved section changes between \
                         stations: no exact vertex correspondence",
                    ));
                }
            }
            let placed = sections
                .iter()
                .zip(&distances)
                .zip(&rings)
                .map(|((section, &s), rings)| {
                    let base = exact.section(s, frame)?;
                    let turned = oriented(&base, &section.orientation, &section.orientation, 0.0)?;
                    let o = section.station.offsets;
                    let frame = Frame3 {
                        origin: base.place(o.lateral, o.vertical, o.longitudinal),
                        x: turned.lateral,
                        y: turned.up,
                        z: turned.lateral.cross(turned.up),
                    };
                    let mut all = Vec::with_capacity(1 + rings.holes.len());
                    all.push(rings.outer.clone());
                    all.extend(rings.holes.iter().cloned());
                    Ok((frame, all))
                })
                .collect::<GeomResult<Vec<_>>>()?;
            if placed[0].0.z.dot(base.tangent) <= 0.0 {
                return Err(unsupported(
                    "an exact station-placed spine whose turned section faces away from its \
                     basis",
                ));
            }
            loft_polygon_stations_exact(&placed, tolerance).map_err(remap(unsupported))
        }
    }
}

/// A construction's typed refusal, renamed as the exact compiler's.
fn remap(unsupported: fn(&'static str) -> GeomError) -> impl Fn(GeomError) -> GeomError {
    move |error| match error {
        GeomError::UnsupportedInput { input, .. } => unsupported(input),
        error => error,
    }
}
