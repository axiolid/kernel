//! Conformance suite every curve-evaluation provider must pass.
//!
//! These are the properties a CALLER is entitled to assume. A provider
//! that passes them can be swapped for another without a consumer
//! noticing, which is the whole point of naming the capability.

use axiolid_contracts::GeomError;
use axiolid_core::{Frame3, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    Circle3, CurvatureLaw, Curve2, Curve3, CurvePath, Elevated3, ElevationLaw, Intrinsic3, Line2,
    Line3, PathCurve, PathPiece, Polyline3, SeamSide,
};

use crate::{CurveEvaluator, CurveMeasure, DistanceConvention};

/// One failed conformance expectation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConformanceFailure {
    /// Which expectation failed.
    pub check: &'static str,
    /// What went wrong.
    pub detail: String,
}

fn fail(check: &'static str, detail: impl Into<String>) -> ConformanceFailure {
    ConformanceFailure {
        check,
        detail: detail.into(),
    }
}

/// Run every conformance check against `provider`.
///
/// Returns the failures; an empty vector means conformant.
#[must_use]
pub fn check<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    out.extend(check_line(provider));
    out.extend(check_circle(provider));
    out.extend(check_tangent_is_unit(provider));
    out.extend(check_frame_is_orthonormal(provider));
    out.extend(check_refusals(provider));
    out.extend(check_measure_routes_differ(provider));
    out.extend(check_sides_agree_off_a_seam(provider));
    out.extend(check_sides_at_seams(provider));
    out.extend(check_paths(provider));
    out
}

/// A line with a NON-UNIT direction must still advance true distance.
///
/// This is the check that catches a provider handing back the native
/// parameter as though it were a distance: with `|direction| = 2` the two
/// differ by a factor of two, while a unit direction would hide the bug.
fn check_line<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    let line = Curve3::Line(Line3 {
        origin: Point3::new(1.0, 2.0, 3.0),
        direction: Vec3::new(2.0, 0.0, 0.0),
    });
    if !provider.distance_convention(&line).is_supported() {
        return out;
    }
    match provider.point_at(&line, CurveMeasure::Distance(5.0)) {
        Ok(p) => {
            let moved = (p - Point3::new(1.0, 2.0, 3.0)).length();
            if (moved - 5.0).abs() > 1e-9 {
                out.push(fail(
                    "line advances true distance",
                    format!("distance 5 moved {moved}, not 5"),
                ));
            }
        }
        Err(e) => out.push(fail("line advances true distance", format!("{e:?}"))),
    }
    out
}

/// A circle of radius r must reach angle d/r at distance d.
fn check_circle<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    let radius = 4.0;
    let circle = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius,
    });
    if !provider.distance_convention(&circle).is_supported() {
        return out;
    }
    // A quarter of the circumference must land on the +y axis.
    let quarter = core::f64::consts::FRAC_PI_2 * radius;
    match provider.point_at(&circle, CurveMeasure::Distance(quarter)) {
        Ok(p) => {
            let want = Point3::new(0.0, radius, 0.0);
            let error = (p - want).length();
            if error > 1e-9 {
                out.push(fail(
                    "circle reaches angle d/r",
                    format!("quarter circumference landed {error:e} away from the +y axis"),
                ));
            }
        }
        Err(e) => out.push(fail("circle reaches angle d/r", format!("{e:?}"))),
    }
    out
}

fn sample_curves() -> Vec<Curve3> {
    vec![
        Curve3::Line(Line3 {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(3.0, 4.0, 0.0),
        }),
        Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: Point3::new(2.0, 0.0, 1.0),
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius: 3.0,
        }),
        // An arc-length-native family. Included because a provider may route
        // it past the distance-to-parameter conversion entirely, so its
        // guards would otherwise go unexercised by this suite.
        Curve3::Intrinsic(Intrinsic3::new(
            Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            CurvatureLaw::circular(0.1),
            CurvatureLaw::circular(0.04),
            10.0,
        )),
    ]
}

/// The tangent must be unit length wherever it is defined.
fn check_tangent_is_unit<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    for curve in sample_curves() {
        if !provider.distance_convention(&curve).is_supported() {
            continue;
        }
        for distance in [0.0, 1.5, 4.0] {
            if let Ok(t) = provider.tangent_at(&curve, CurveMeasure::Distance(distance)) {
                if (t.length() - 1.0).abs() > 1e-9 {
                    out.push(fail(
                        "tangent is unit length",
                        format!("|t| = {} at distance {distance}", t.length()),
                    ));
                }
            }
        }
    }
    out
}

/// The frame must be right-handed orthonormal with `x` on the tangent.
///
/// A caller builds a rotation from this directly, so a frame that is
/// merely close to orthonormal shears whatever it places.
fn check_frame_is_orthonormal<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    for curve in sample_curves() {
        if !provider.distance_convention(&curve).is_supported() {
            continue;
        }
        for distance in [0.0, 1.5, 4.0] {
            let Ok(frame) = provider.frame_at(&curve, CurveMeasure::Distance(distance)) else {
                continue;
            };
            let checks = [
                ("x unit", frame.x.length() - 1.0),
                ("y unit", frame.y.length() - 1.0),
                ("z unit", frame.z.length() - 1.0),
                ("x.y", frame.x.dot(frame.y)),
                ("x.z", frame.x.dot(frame.z)),
                ("y.z", frame.y.dot(frame.z)),
            ];
            for (name, value) in checks {
                if value.abs() > 1e-9 {
                    out.push(fail(
                        "frame is orthonormal",
                        format!("{name} off by {value:e} at distance {distance}"),
                    ));
                }
            }
            // Right-handed: x cross y must be z, not -z.
            let handed = frame.x.cross(frame.y).dot(frame.z);
            if (handed - 1.0).abs() > 1e-9 {
                out.push(fail(
                    "frame is right-handed",
                    format!("x cross y . z = {handed} at distance {distance}"),
                ));
            }
            // The frame must sit ON the curve.
            if let Ok(point) = provider.point_at(&curve, CurveMeasure::Distance(distance)) {
                let off = (frame.origin - point).length();
                if off > 1e-9 {
                    out.push(fail(
                        "frame origin is the curve point",
                        format!("origin {off:e} away at distance {distance}"),
                    ));
                }
            }
            // And its x axis must be the tangent.
            if let Ok(tangent) = provider.tangent_at(&curve, CurveMeasure::Distance(distance)) {
                let off = (frame.x - tangent).length();
                if off > 1e-9 {
                    out.push(fail(
                        "frame x is the tangent",
                        format!("x differs from tangent by {off:e}"),
                    ));
                }
            }
        }
    }
    out
}

/// Bad input must be refused, not answered with a guess.
///
/// A NaN distance that returns a NaN point is the dangerous case: it
/// propagates silently into a placement instead of failing at the call.
fn check_refusals<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    for curve in sample_curves() {
        if !provider.distance_convention(&curve).is_supported() {
            continue;
        }
        for bad in [Scalar::NAN, Scalar::INFINITY, Scalar::NEG_INFINITY] {
            if provider
                .point_at(&curve, CurveMeasure::Distance(bad))
                .is_ok()
            {
                out.push(fail(
                    "non-finite distance is refused",
                    format!("point_at accepted {bad}"),
                ));
            }
            if provider
                .tangent_at(&curve, CurveMeasure::Distance(bad))
                .is_ok()
            {
                out.push(fail(
                    "non-finite distance is refused",
                    format!("tangent_at accepted {bad}"),
                ));
            }
            if provider
                .frame_at(&curve, CurveMeasure::Distance(bad))
                .is_ok()
            {
                out.push(fail(
                    "non-finite distance is refused",
                    format!("frame_at accepted {bad}"),
                ));
            }
        }
    }
    // A provider must not claim a convention it cannot honour: if it
    // reports Unsupported it must actually refuse.
    let vertical = Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::ZERO,
    });
    if provider.distance_convention(&vertical).is_supported()
        && provider
            .point_at(&vertical, CurveMeasure::Distance(1.0))
            .is_ok()
    {
        out.push(fail(
            "degenerate curve is refused",
            "a zero-direction line was evaluated instead of refused",
        ));
    }
    out
}

/// A parameter and a distance must not be silently interchangeable.
///
/// On a circle of radius 4 the same number means two different places:
/// `Parameter(1.5)` is 1.5 radians round, `Distance(1.5)` is 1.5 m along,
/// i.e. 0.375 rad. A provider that ignored the method of measurement would
/// return the same point for both -- wrong, finite, and plausible.
///
/// Also pins that the parameter route stays OPEN where the distance route
/// is refused: an ellipse has no closed-form arc length, but its native
/// parameter is perfectly meaningful, and a consumer holding an authored
/// parameter must still be able to evaluate it.
fn check_measure_routes_differ<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    let radius = 4.0;
    let circle = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius,
    });
    let value = 1.5;
    let by_parameter = provider.point_at(&circle, CurveMeasure::Parameter(value));
    let by_distance = provider.point_at(&circle, CurveMeasure::Distance(value));
    if let (Ok(p), Ok(d)) = (by_parameter, by_distance) {
        if (p - d).length() < 1e-9 {
            out.push(fail(
                "parameter and distance are distinct",
                format!("{value} gave the same point as a parameter and as a distance"),
            ));
        }
    }
    out
}

// --- seam sides (#286) --------------------------------------------------------

/// Whether a sided query was refused because the provider does not read
/// that side, as the contract's defaults refuse it: unsupported, typed.
fn side_refused(error: &GeomError) -> bool {
    matches!(
        error,
        GeomError::Unsupported { .. } | GeomError::UnsupportedInput { .. }
    )
}

fn frame_off(a: &Frame3, b: &Frame3) -> Scalar {
    [
        (a.origin - b.origin).length(),
        (a.x - b.x).length(),
        (a.y - b.y).length(),
        (a.z - b.z).length(),
    ]
    .into_iter()
    .fold(0.0, Scalar::max)
}

/// The three sided answers at `at`, or `None` when the provider refused
/// `Incoming` as unsupported on all three (it does not read sides).
type Sided = Option<(Point3, Vec3, Frame3)>;

fn sided<E: CurveEvaluator>(
    provider: &E,
    curve: &Curve3,
    at: CurveMeasure,
    side: SeamSide,
) -> Result<Sided, String> {
    match (
        provider.point_at_on(curve, at, side),
        provider.tangent_at_on(curve, at, side),
        provider.frame_at_on(curve, at, side),
    ) {
        (Ok(p), Ok(t), Ok(f)) => Ok(Some((p, t, f))),
        (Err(a), Err(b), Err(c))
            if side == SeamSide::Incoming
                && side_refused(&a)
                && side_refused(&b)
                && side_refused(&c) =>
        {
            Ok(None)
        }
        (p, t, f) => Err(format!(
            "point {:?}, tangent {:?}, frame {:?}",
            p.err(),
            t.err(),
            f.err()
        )),
    }
}

/// Off a seam both sides read the same piece: `Outgoing` must be the
/// side-less answer, and `Incoming`, where a provider reads it, too. A
/// provider that does not read `Incoming` must refuse it as UNSUPPORTED,
/// never answer it with something else or refuse it as bad input. A
/// sided query keeps the side-less refusals of a non-finite measure.
fn check_sides_agree_off_a_seam<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    for curve in sample_curves() {
        if !provider.distance_convention(&curve).is_supported() {
            continue;
        }
        for distance in [0.0, 1.5, 4.0] {
            let at = CurveMeasure::Distance(distance);
            let (Ok(point), Ok(tangent), Ok(frame)) = (
                provider.point_at(&curve, at),
                provider.tangent_at(&curve, at),
                provider.frame_at(&curve, at),
            ) else {
                continue;
            };
            for side in [SeamSide::Outgoing, SeamSide::Incoming] {
                match sided(provider, &curve, at, side) {
                    Ok(Some((p, t, f))) => {
                        let off = (p - point)
                            .length()
                            .max((t - tangent).length())
                            .max(frame_off(&f, &frame));
                        if off > 1e-9 {
                            out.push(fail(
                                "both sides agree off a seam",
                                format!(
                                    "{side:?} is {off:e} off the side-less answer at {distance}"
                                ),
                            ));
                        }
                    }
                    Ok(None) => {}
                    Err(detail) => out.push(fail(
                        "a side is answered or refused as unsupported",
                        format!("{side:?} at {distance}: {detail}"),
                    )),
                }
            }
        }
        for bad in [Scalar::NAN, Scalar::INFINITY] {
            let at = CurveMeasure::Distance(bad);
            for side in [SeamSide::Outgoing, SeamSide::Incoming] {
                if provider.point_at_on(&curve, at, side).is_ok()
                    || provider.tangent_at_on(&curve, at, side).is_ok()
                    || provider.frame_at_on(&curve, at, side).is_ok()
                {
                    out.push(fail(
                        "non-finite distance is refused",
                        format!("a {side:?} query accepted {bad}"),
                    ));
                }
            }
        }
    }
    out
}

/// A seam, its distance, the point there and the unit tangents of the
/// pieces ending and starting there, from closed forms.
struct Seam {
    what: &'static str,
    curve: Curve3,
    distance: Scalar,
    point: Point3,
    incoming: Vec3,
    outgoing: Vec3,
}

fn seams() -> Vec<Seam> {
    let grade = |g: Scalar| Vec3::new(1.0, 0.0, g) / g.hypot(1.0);
    vec![
        // A polyline's corner: 4 m along +x, then 3 m along +y.
        Seam {
            what: "polyline corner",
            curve: Curve3::Polyline(Polyline3 {
                points: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(4.0, 0.0, 0.0),
                    Point3::new(4.0, 3.0, 0.0),
                ],
                closed: false,
            }),
            distance: 4.0,
            point: Point3::new(4.0, 0.0, 0.0),
            incoming: Vec3::X,
            outgoing: Vec3::Y,
        },
        // A grade break at plan distance 50 on a straight plan: +2% up to
        // it, -1% after.
        Seam {
            what: "grade break",
            curve: Curve3::Elevated(Elevated3::new(
                Curve2::Line(Line2 {
                    origin: Point2::ZERO,
                    direction: Vec2::X,
                }),
                ElevationLaw::Piecewise {
                    breaks: vec![50.0],
                    laws: vec![
                        ElevationLaw::constant_grade(0.0, 0.02),
                        ElevationLaw::constant_grade(1.0, -0.01),
                    ],
                },
            )),
            distance: 50.0,
            point: Point3::new(50.0, 0.0, 1.0),
            incoming: grade(0.02),
            outgoing: grade(-0.01),
        },
    ]
}

/// On a seam each side reads its own piece -- for a provider that reads
/// sides also a hair either side of the seam, within the ADR 0082
/// tolerance `1e-12 * max(1, s)` -- and `Outgoing` is what the side-less
/// query reads exactly on it. A provider that does not read `Incoming`
/// must refuse it as unsupported rather than answer with the outgoing
/// piece.
fn check_sides_at_seams<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    for seam in seams() {
        if !provider.distance_convention(&seam.curve).is_supported() {
            continue;
        }
        let (what, curve) = (seam.what, &seam.curve);
        let on = CurveMeasure::Distance(seam.distance);
        if let (Ok(t), Ok(t_on)) = (
            provider.tangent_at(curve, on),
            provider.tangent_at_on(curve, on, SeamSide::Outgoing),
        ) {
            if (t - t_on).length() > 1e-9 {
                out.push(fail(
                    "outgoing is the side-less reading on a seam",
                    format!("{what}: {t:?} side-less, {t_on:?} outgoing"),
                ));
            }
        }
        // The tolerance is part of the seam rule a provider reading sides
        // implements; one that does not answers `Outgoing` side-lessly,
        // which reads the piece the measure falls in.
        let reads_sides = matches!(sided(provider, curve, on, SeamSide::Incoming), Ok(Some(_)));
        let hair = 0.25e-12 * seam.distance;
        let distances = if reads_sides {
            vec![seam.distance, seam.distance - hair, seam.distance + hair]
        } else {
            vec![seam.distance]
        };
        for distance in distances {
            let at = CurveMeasure::Distance(distance);
            for (side, want) in [
                (SeamSide::Outgoing, seam.outgoing),
                (SeamSide::Incoming, seam.incoming),
            ] {
                let (point, tangent, frame) = match sided(provider, curve, at, side) {
                    Ok(Some(answers)) => answers,
                    Ok(None) => continue,
                    Err(detail) => {
                        out.push(fail(
                            "a seam side reads its own piece",
                            format!("{what}, {side:?} at {distance}: {detail}"),
                        ));
                        continue;
                    }
                };
                let checks = [
                    ("point", (point - seam.point).length()),
                    ("tangent", (tangent - want).length()),
                    ("frame x", (frame.x - want).length()),
                    ("frame origin", (frame.origin - seam.point).length()),
                    ("frame x.y", frame.x.dot(frame.y)),
                    (
                        "frame handedness",
                        frame.x.cross(frame.y).dot(frame.z) - 1.0,
                    ),
                ];
                for (name, value) in checks {
                    if value.abs() > 1e-9 {
                        out.push(fail(
                            "a seam side reads its own piece",
                            format!("{what}, {side:?} at {distance}: {name} off by {value:e}"),
                        ));
                    }
                }
            }
        }
    }
    out
}

// --- curve paths (#290) -------------------------------------------------------

/// A place along a path and what each side reads there, from closed forms.
struct PathPlace {
    distance: Scalar,
    point: Point3,
    incoming: Vec3,
    outgoing: Vec3,
    /// Whether each side's frame is read on or after the arc, so never
    /// exact: `(incoming, outgoing)`.
    past_arc: (bool, bool),
}

/// Every path answer at `at` from `side`, or `None` when the provider
/// refused all of them as unsupported (it does not read paths).
type PathAnswers = Option<(Point3, Vec3, Frame3)>;

fn path_answers<E: CurveEvaluator>(
    provider: &E,
    path: &CurvePath,
    at: CurveMeasure,
    side: SeamSide,
) -> Result<PathAnswers, String> {
    match (
        provider.path_point_at_on(path, at, side),
        provider.path_tangent_at_on(path, at, side),
        provider.path_frame_at_on(path, at, side),
    ) {
        (Ok(p), Ok(t), Ok(f)) => Ok(Some((p, t, f))),
        (Err(a), Err(b), Err(c)) if side_refused(&a) && side_refused(&b) && side_refused(&c) => {
            Ok(None)
        }
        (p, t, f) => Err(format!(
            "point {:?}, tangent {:?}, frame {:?}",
            p.err(),
            t.err(),
            f.err()
        )),
    }
}

/// The places [`check_paths`] reads its path at, forwards and reversed.
fn path_places(radius: Scalar) -> [[PathPlace; 3]; 2] {
    let quarter = core::f64::consts::FRAC_PI_2 * radius;
    let length = 10.0 + quarter;
    let (sin, cos) = core::f64::consts::FRAC_PI_4.sin_cos();
    let mid = Point3::new(15.0 - radius * cos, radius * sin, 0.0);
    let mid_tangent = Vec3::new(sin, cos, 0.0);
    let joint = Point3::new(10.0, 0.0, 0.0);
    let four = Point3::new(4.0, 0.0, 0.0);
    [
        [
            PathPlace {
                distance: 4.0,
                point: four,
                incoming: Vec3::X,
                outgoing: Vec3::X,
                past_arc: (false, false),
            },
            PathPlace {
                distance: 10.0,
                point: joint,
                incoming: Vec3::X,
                outgoing: Vec3::Y,
                past_arc: (false, true),
            },
            PathPlace {
                distance: 10.0 + quarter / 2.0,
                point: mid,
                incoming: mid_tangent,
                outgoing: mid_tangent,
                past_arc: (true, true),
            },
        ],
        // Reversed: the arc from its end, then the line back to the origin.
        [
            PathPlace {
                distance: quarter / 2.0,
                point: mid,
                incoming: -mid_tangent,
                outgoing: -mid_tangent,
                past_arc: (true, true),
            },
            PathPlace {
                distance: quarter,
                point: joint,
                incoming: -Vec3::Y,
                outgoing: -Vec3::X,
                past_arc: (true, true),
            },
            PathPlace {
                distance: length - 4.0,
                point: four,
                incoming: -Vec3::X,
                outgoing: -Vec3::X,
                past_arc: (true, true),
            },
        ],
    ]
}

/// The failures of one path answer against its closed form.
fn path_answer_failures(
    what: &str,
    side: SeamSide,
    distance: Scalar,
    (point, tangent, frame): (Point3, Vec3, Frame3),
    place: &PathPlace,
    want: Vec3,
) -> Vec<ConformanceFailure> {
    let checks = [
        ("point", (point - place.point).length()),
        ("tangent", (tangent - want).length()),
        ("frame x", (frame.x - want).length()),
        ("frame origin", (frame.origin - place.point).length()),
        ("frame x.y", frame.x.dot(frame.y)),
        ("frame x.z", frame.x.dot(frame.z)),
        ("frame y.z", frame.y.dot(frame.z)),
        ("frame y unit", frame.y.length() - 1.0),
        (
            "frame handedness",
            frame.x.cross(frame.y).dot(frame.z) - 1.0,
        ),
    ];
    checks
        .into_iter()
        .filter(|(_, value)| value.abs() > 1e-9)
        .map(|(name, value)| {
            fail(
                "a path reads each side of a joint from its own piece",
                format!("{what} path, {side:?} at {distance}: {name} off by {value:e}"),
            )
        })
        .collect()
}

/// A path is read end to end as ADR 0082 reads a composite station basis:
/// a line along `+x` over `[0, 10]` (a non-unit direction, so its measure
/// is not its parameter), then a quarter of a circle of radius 5 about
/// `(15, 0, 0)` turning left to `(15, 5, 0)` -- a right-angle corner at
/// distance 10, incoming `+x`, outgoing `+y` -- read forwards and
/// reversed, on the joint and within the seam tolerance either side of
/// it. A provider that does not read paths must refuse every path query
/// as unsupported and report no convention; one that reads them must read
/// each side of the joint from its own piece, answer the plain queries as
/// the outgoing side, refuse a distance off the path and a path whose
/// pieces do not meet, measure it by arc length, and claim no frame on or
/// after the arc exact.
fn check_paths<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    let line = Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::new(2.0, 0.0, 0.0),
    });
    let radius = 5.0;
    let circle = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::new(15.0, 0.0, 0.0),
            x: -Vec3::X,
            y: Vec3::Y,
            z: -Vec3::Z,
        },
        radius,
    });
    let quarter = core::f64::consts::FRAC_PI_2 * radius;
    let length = 10.0 + quarter;
    let forwards = CurvePath::new(vec![
        PathPiece::new(PathCurve::Three(line.clone()), 0.0, 10.0),
        PathPiece::new(PathCurve::Three(circle.clone()), 0.0, quarter),
    ]);
    let [forwards_places, reversed_places] = path_places(radius);
    let paths = [
        ("forwards", forwards.clone(), forwards_places),
        ("reversed", forwards.clone().reversed(), reversed_places),
    ];
    let mut reads_paths = false;
    for (what, path, places) in &paths {
        for place in places {
            let hair = 0.25e-12 * place.distance;
            for distance in [place.distance - hair, place.distance, place.distance + hair] {
                let at = CurveMeasure::Distance(distance);
                for (side, want, past_arc) in [
                    (SeamSide::Outgoing, place.outgoing, place.past_arc.1),
                    (SeamSide::Incoming, place.incoming, place.past_arc.0),
                ] {
                    let answers = match path_answers(provider, path, at, side) {
                        Ok(Some(answers)) => answers,
                        Ok(None) => continue,
                        Err(detail) => {
                            out.push(fail(
                                "a path is answered or refused as unsupported",
                                format!("{what} path, {side:?} at {distance}: {detail}"),
                            ));
                            continue;
                        }
                    };
                    reads_paths = true;
                    out.extend(path_answer_failures(
                        what, side, distance, answers, place, want,
                    ));
                    if past_arc && provider.path_frame_is_exact_at(path, at, side) {
                        out.push(fail(
                            "a path frame is exact only on exactly placed lines",
                            format!("{what} path, {side:?} at {distance}: claimed exact"),
                        ));
                    }
                    if side != SeamSide::Outgoing {
                        continue;
                    }
                    // The plain queries are the outgoing reading.
                    let (point, tangent, frame) = answers;
                    let off = match (
                        provider.path_point_at(path, at),
                        provider.path_tangent_at(path, at),
                        provider.path_frame_at(path, at),
                    ) {
                        (Ok(p), Ok(t), Ok(f)) => (p - point)
                            .length()
                            .max((t - tangent).length())
                            .max(frame_off(&f, &frame)),
                        _ => Scalar::INFINITY,
                    };
                    if off > 1e-9 {
                        out.push(fail(
                            "a plain path query is the outgoing reading",
                            format!("{what} path at {distance}: {off:e} off"),
                        ));
                    }
                }
            }
        }
    }
    let convention = provider.path_distance_convention(&forwards);
    if !reads_paths {
        if convention.is_supported() {
            out.push(fail(
                "a provider that refuses paths reports them unsupported",
                format!("every path query was refused, but {convention:?} was reported"),
            ));
        }
        return out;
    }
    if convention != DistanceConvention::ArcLength3d {
        out.push(fail(
            "a path of lines and arcs is measured by arc length",
            format!("{convention:?}"),
        ));
    }
    for bad in [Scalar::NAN, Scalar::INFINITY, -1.0, length + 1e-3] {
        let at = CurveMeasure::Distance(bad);
        for side in [SeamSide::Outgoing, SeamSide::Incoming] {
            if provider.path_point_at_on(&forwards, at, side).is_ok()
                || provider.path_tangent_at_on(&forwards, at, side).is_ok()
                || provider.path_frame_at_on(&forwards, at, side).is_ok()
            {
                out.push(fail(
                    "a distance off the path is refused",
                    format!("a {side:?} path query accepted {bad}"),
                ));
            }
        }
    }
    // Pieces that do not meet: no distance runs across the gap.
    let gapped = CurvePath::new(vec![
        PathPiece::new(PathCurve::Three(line.clone()), 0.0, 9.0),
        PathPiece::new(PathCurve::Three(circle.clone()), 0.0, quarter),
    ]);
    if provider
        .path_point_at(&gapped, CurveMeasure::Distance(1.0))
        .is_ok()
    {
        out.push(fail(
            "a path whose pieces do not meet is refused",
            "a path with a 1 m gap was evaluated",
        ));
    }
    out
}
