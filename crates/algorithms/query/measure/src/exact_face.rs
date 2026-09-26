//! Mass properties of one curved face, by Green's theorem in its parameters.
//!
//! # Method
//!
//! The volume a closed boundary encloses, and its first and second moments,
//! are sums over the faces of `int F . n dA` for a field `F` whose divergence
//! is the wanted density (1, `x_i`, `x_i^2`). [`crate::exact`] uses the cone
//! fields `x / 3`, `x_i x / 4` and `x_i^2 x / 5` so each face contributes the
//! solid cone it subtends from the origin. A planar polygon's cone is a fan of
//! tetrahedra; a curved face's cone is the surface integral
//!
//! ```text
//! int int_D g(u, v) du dv,   g = (S . N) * [1/3, S/4, S^2/5],  N = S_u x S_v
//! ```
//!
//! over the face's domain `D` in the support surface's own `(u, v)`
//! parameters. Area is the same integral of `|N|`.
//!
//! # Orientation
//!
//! A loop anticlockwise in `(u, v)` runs anticlockwise about `N` in space,
//! whatever the handedness of the surface's frame. So the loops' winding
//! about `N` -- the sense the planar fan, the tessellator and edge-use
//! pairing all read -- is exactly the sign Green's theorem gives the domain,
//! and no separate normal is consulted.
//!
//! `D` is bounded by the pcurves of the face's loops, so Green's theorem turns
//! the double integral into a line integral round those pcurves:
//!
//! ```text
//! int int_D g du dv = - oint H(u, v) du,   H(u, v) = int_{v_ref}^{v} g(u, s) ds
//! ```
//!
//! Nothing is sampled into a polygon. Both integrals are adaptive
//! Gauss-Kronrod (G7/K15) over the exact surface and the exact pcurves, and a
//! piece is only accepted once its error estimate is below a relative bound
//! near machine precision; otherwise the face is refused with
//! [`ExactMeasureError::NotConverged`]. For a cylinder or cone the inner
//! integrand is a polynomial in `v` of degree at most four, which K15
//! integrates exactly; the angular direction converges geometrically.
//!
//! # Seams and poles
//!
//! A periodic surface's loops need not close in the parameter plane. A
//! cylinder wall bounded by two full circles has one loop running `u = 0 ..
//! 2 pi` at the bottom and one running back at the top, joined only through
//! the seam. Along a seam `u` is constant, so `- H du` vanishes there and the
//! seam contributes nothing: the two open loops already carry the whole
//! integral, and a missing seam edge is harmless. `H` is periodic in `u`, so
//! a pcurve stated one turn further round measures the same.
//!
//! When the loops' net winding in `u` is not zero, the domain reaches a pole
//! (a spherical cap, a cone's apex). The missing boundary is then the pole
//! itself, where `S` is one point; integrating `H` from `v_ref = v_pole`
//! makes that boundary contribute exactly zero. A net winding with no pole on
//! the domain's side (a lone circle on a cylinder) does not bound anything and
//! is refused.
//!
//! A surface periodic in `v` (a torus) whose loops wind in `v` instead is
//! integrated the other way round, `oint G dv` with `G` an integral in `u`.

use axiolid_brep::ExactBRep;
use axiolid_core::{Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::Curve2;
use axiolid_evaluate::surface::{evaluate, partials};
use axiolid_evaluate::{derivative2, evaluate2};
use axiolid_surface::Surface;
use axiolid_topology::{Face, Orientation};
use core::f64::consts::{FRAC_PI_2, TAU};

use crate::exact::{ExactMeasureError, Sums, COMPONENTS};

/// Relative error accepted on each Gauss-Kronrod piece.
const RELATIVE: Scalar = 1e-13;

/// Pieces one adaptive integral may split into before it is refused.
const MAX_PIECES: usize = 2048;

/// Length dimension of each component, for the absolute noise floor:
/// area, volume, three first moments, three second moments.
const DIMENSION: [i32; COMPONENTS] = [2, 3, 4, 4, 4, 5, 5, 5];

/// How a face's surface is parameterised: its periods and its poles.
struct Chart {
    u_period: Option<Scalar>,
    v_period: Option<Scalar>,
    /// `v` values where the surface collapses to one point for every `u`.
    poles: Vec<Scalar>,
}

fn chart(surface: &Surface) -> Result<Chart, ExactMeasureError> {
    let chart = match surface {
        Surface::Plane(_) | Surface::BSpline(_) => Chart {
            u_period: None,
            v_period: None,
            poles: Vec::new(),
        },
        Surface::Cylinder(_) | Surface::EllipticalCylinder(_) => Chart {
            u_period: Some(TAU),
            v_period: None,
            poles: Vec::new(),
        },
        Surface::Cone(cone) => {
            let slope = cone.semi_angle.tan();
            let poles = if slope.is_finite() && slope != 0.0 {
                vec![-cone.radius / slope]
            } else {
                Vec::new()
            };
            Chart {
                u_period: Some(TAU),
                v_period: None,
                poles,
            }
        }
        Surface::Sphere(_) => Chart {
            u_period: Some(TAU),
            v_period: None,
            poles: vec![-FRAC_PI_2, FRAC_PI_2],
        },
        Surface::Torus(_) => Chart {
            u_period: Some(TAU),
            v_period: Some(TAU),
            poles: Vec::new(),
        },
        _ => {
            return Err(ExactMeasureError::NonPlanarFace(crate::exact::family(
                surface,
            )))
        }
    };
    Ok(chart)
}

/// One stretch of a loop in the parameter plane.
enum Piece<'a> {
    /// A pcurve over `[start, end]`, shifted by whole periods.
    Curve {
        curve: &'a Curve2,
        start: Scalar,
        end: Scalar,
        offset: Vec2,
    },
    /// A straight closing stretch: across a pole, or a gap within tolerance.
    Segment { from: Point2, to: Point2 },
}

impl Piece<'_> {
    fn span(&self) -> (Scalar, Scalar) {
        match self {
            Self::Curve { start, end, .. } => (*start, *end),
            Self::Segment { .. } => (0.0, 1.0),
        }
    }

    fn at(&self, t: Scalar) -> Result<(Point2, Vec2), ExactMeasureError> {
        match self {
            Self::Curve { curve, offset, .. } => {
                let point = evaluate2(curve, t).map_err(|_| ExactMeasureError::Evaluation)?;
                let tangent = derivative2(curve, t).map_err(|_| ExactMeasureError::Evaluation)?;
                Ok((point + *offset, tangent))
            }
            Self::Segment { from, to } => Ok((*from + (*to - *from) * t, *to - *from)),
        }
    }

    fn end_point(&self) -> Result<Point2, ExactMeasureError> {
        Ok(self.at(self.span().1)?.0)
    }
}

/// A face's boundary, assembled in the parameter plane.
struct Boundary<'a> {
    pieces: Vec<Piece<'a>>,
    /// Net whole periods each loop winds, summed over the face.
    winding: [i64; 2],
    /// Whether any single loop winds in `u` / in `v`.
    wraps: [bool; 2],
    /// A point on the boundary, for choosing the reference line.
    anchor: Point2,
}

/// Integrate one face's contribution to every component.
///
/// Returns the domain integral with the loops read as the face uses them:
/// positive for a face whose outer loop runs anticlockwise in `(u, v)`.
/// The caller applies the face's own orientation to the volume terms.
pub(crate) fn face_sums(
    brep: &ExactBRep,
    face: &Face<axiolid_brep::SurfaceId>,
    surface: &Surface,
    linear: Scalar,
    scale: Scalar,
) -> Result<Sums, ExactMeasureError> {
    let chart = chart(surface)?;
    let boundary = assemble(brep, face, surface, &chart, linear)?;

    let floor = noise_floor(scale);
    let along_v = boundary.wraps[1];
    if along_v {
        // Loops wind round the tube: integrate in `u` first, `oint G dv`.
        if boundary.wraps[0] || boundary.winding[1] != 0 {
            return Err(ExactMeasureError::ParameterDomain(
                "face boundary winds around the surface in both directions",
            ));
        }
    }

    let reference = if along_v {
        boundary.anchor.x
    } else if boundary.winding[0] != 0 {
        pole_on_domain_side(&chart, &boundary)?
    } else {
        boundary.anchor.y
    };

    let mut total = [0.0; COMPONENTS];
    for piece in &boundary.pieces {
        let (a, b) = piece.span();
        let sums = adaptive(a, b, &floor, &mut |t| {
            let (point, tangent) = piece.at(t)?;
            // `- H du` in the default direction, `+ G dv` along the tube.
            let (weight, inner) = if along_v {
                (tangent.y, inner_along_u(surface, reference, point, &floor)?)
            } else {
                (
                    -tangent.x,
                    inner_along_v(surface, reference, point, &floor)?,
                )
            };
            let mut value = inner;
            for slot in &mut value {
                *slot *= weight;
            }
            Ok(value)
        })?;
        for (slot, value) in total.iter_mut().zip(sums) {
            *slot += value;
        }
    }
    Ok(total)
}

/// The pole the face's domain reaches, used as the reference line so the
/// pole itself -- the boundary no loop states -- contributes zero.
fn pole_on_domain_side(
    chart: &Chart,
    boundary: &Boundary<'_>,
) -> Result<Scalar, ExactMeasureError> {
    // A loop running +u keeps its domain on its left, towards +v.
    let upward = match boundary.winding[0] {
        1 => true,
        -1 => false,
        _ => {
            return Err(ExactMeasureError::ParameterDomain(
                "face boundary winds around the surface more than once",
            ))
        }
    };
    chart
        .poles
        .iter()
        .copied()
        .filter(|pole| (*pole > boundary.anchor.y) == upward)
        .min_by(|a, b| {
            (a - boundary.anchor.y)
                .abs()
                .total_cmp(&(b - boundary.anchor.y).abs())
        })
        .ok_or(ExactMeasureError::ParameterDomain(
            "face boundary winds around a surface with no pole on the domain side",
        ))
}

/// Assemble every bound's pcurves into one boundary, joining each use to the
/// next through whole periods, across a pole, or over a gap within tolerance.
fn assemble<'a>(
    brep: &'a ExactBRep,
    face: &Face<axiolid_brep::SurfaceId>,
    surface: &Surface,
    chart: &Chart,
    linear: Scalar,
) -> Result<Boundary<'a>, ExactMeasureError> {
    let topology = brep.topology();
    let mut pieces = Vec::new();
    let mut winding = [0_i64; 2];
    let mut wraps = [false; 2];
    let mut anchor = None;

    for bound in &face.bounds {
        let wire = topology
            .loops()
            .get(bound.loop_id.index())
            .ok_or(ExactMeasureError::DanglingReference)?;
        let reversed = bound.orientation == Orientation::Reversed;

        // Each use's pcurve interval already runs in the use's traversal
        // (ADR 0024); a reversed bound walks the loop backwards.
        let mut uses = Vec::with_capacity(wire.edges.len());
        for (index, use_) in wire.edges.iter().enumerate() {
            let pcurve = use_.pcurve.ok_or(ExactMeasureError::DanglingReference)?;
            let curve = brep
                .curves2()
                .get(pcurve.index())
                .ok_or(ExactMeasureError::DanglingReference)?;
            let interval = brep
                .pcurve_interval(bound.loop_id, index)
                .ok_or(ExactMeasureError::DanglingReference)?;
            let (start, end) = if reversed {
                (interval.end, interval.start)
            } else {
                (interval.start, interval.end)
            };
            uses.push((curve, start, end));
        }
        if reversed {
            uses.reverse();
        }
        let Some(&(first_curve, first_start, _)) = uses.first() else {
            continue;
        };

        let loop_start =
            evaluate2(first_curve, first_start).map_err(|_| ExactMeasureError::Evaluation)?;
        anchor.get_or_insert(loop_start);
        let mut offset = Vec2::ZERO;
        let mut cursor = loop_start;
        for (curve, start, end) in uses {
            let raw = evaluate2(curve, start).map_err(|_| ExactMeasureError::Evaluation)?;
            let (shift, bridge) = join(surface, chart, cursor, raw + offset, linear)?;
            offset += shift;
            if let Some(segment) = bridge {
                pieces.push(segment);
            }
            let piece = Piece::Curve {
                curve,
                start,
                end,
                offset,
            };
            cursor = piece.end_point()?;
            pieces.push(piece);
        }

        // Close the loop back onto its own start.
        let (shift, bridge) = join(surface, chart, cursor, loop_start, linear)?;
        if let Some(segment) = bridge {
            pieces.push(segment);
        }
        // `shift` moved the start onto the end: the loop ended that many
        // periods on from where it began.
        let turns = [
            periods(shift.x, chart.u_period),
            periods(shift.y, chart.v_period),
        ];
        for axis in 0..2 {
            winding[axis] += turns[axis];
            wraps[axis] |= turns[axis] != 0;
        }
    }

    let anchor = anchor.ok_or(ExactMeasureError::Degenerate)?;
    Ok(Boundary {
        pieces,
        winding,
        wraps,
        anchor,
    })
}

/// Whole periods in `shift`, which is already a multiple of `period`.
fn periods(shift: Scalar, period: Option<Scalar>) -> i64 {
    period.map_or(0, |period| (shift / period).round() as i64)
}

/// Join a loop that has reached `from` to the next use starting at `to`.
///
/// Returns the whole-period shift to apply to `to` and everything after it,
/// and the closing segment, if one is needed. The two must be the same point
/// on the surface; the parameter gap between them may be whole periods (a
/// seam), a stretch along a pole, or a residue within tolerance.
fn join<'a>(
    surface: &Surface,
    chart: &Chart,
    from: Point2,
    to: Point2,
    linear: Scalar,
) -> Result<(Vec2, Option<Piece<'a>>), ExactMeasureError> {
    let there = point(surface, from)?;
    let here = point(surface, to)?;
    if (there - here).length() > linear {
        return Err(ExactMeasureError::ParameterDomain(
            "consecutive pcurves of a face loop do not meet on the surface",
        ));
    }

    // Along a pole `u` is free: keep the stretch as it is, so its `H du`
    // is integrated, rather than folding it into periods.
    let on_pole = chart.poles.iter().any(|pole| {
        (from.y - pole).abs() <= parameter_slack(*pole)
            && (to.y - pole).abs() <= parameter_slack(*pole)
    });
    let shift = if on_pole {
        Vec2::ZERO
    } else {
        let wrap = |gap: Scalar, period: Option<Scalar>| {
            period.map_or(0.0, |period| (gap / period).round() * period)
        };
        Vec2::new(
            wrap(from.x - to.x, chart.u_period),
            wrap(from.y - to.y, chart.v_period),
        )
    };
    let landed = to + shift;
    let bridge = (landed != from).then_some(Piece::Segment { from, to: landed });
    Ok((shift, bridge))
}

fn parameter_slack(value: Scalar) -> Scalar {
    1e-9 * value.abs().max(1.0)
}

fn point(surface: &Surface, at: Point2) -> Result<Point3, ExactMeasureError> {
    evaluate(surface, at.x, at.y).map_err(|_| ExactMeasureError::Evaluation)
}

/// `H(u, v)`: the density integrated in `v` from the reference line.
fn inner_along_v(
    surface: &Surface,
    reference: Scalar,
    at: Point2,
    floor: &Sums,
) -> Result<Sums, ExactMeasureError> {
    adaptive(reference, at.y, floor, &mut |v| density(surface, at.x, v))
}

/// `G(u, v)`: the density integrated in `u` from the reference line.
fn inner_along_u(
    surface: &Surface,
    reference: Scalar,
    at: Point2,
    floor: &Sums,
) -> Result<Sums, ExactMeasureError> {
    adaptive(reference, at.x, floor, &mut |u| density(surface, u, at.y))
}

/// Every component's density at `(u, v)`: `|N|`, then the cone fields
/// dotted with `N = S_u x S_v`.
fn density(surface: &Surface, u: Scalar, v: Scalar) -> Result<Sums, ExactMeasureError> {
    let p = evaluate(surface, u, v).map_err(|_| ExactMeasureError::Evaluation)?;
    let (su, sv) = partials(surface, u, v).map_err(|_| ExactMeasureError::Evaluation)?;
    let n: Vec3 = su.cross(sv);
    let w = p.dot(n);
    Ok([
        n.length(),
        w / 3.0,
        p.x * w / 4.0,
        p.y * w / 4.0,
        p.z * w / 4.0,
        p.x * p.x * w / 5.0,
        p.y * p.y * w / 5.0,
        p.z * p.z * w / 5.0,
    ])
}

/// Absolute error below which a component is rounding noise: a component
/// that is zero by symmetry never meets a purely relative bound.
fn noise_floor(scale: Scalar) -> Sums {
    let mut floor = [0.0; COMPONENTS];
    for (slot, dimension) in floor.iter_mut().zip(DIMENSION) {
        *slot = 1e-15 * scale.powi(dimension);
    }
    floor
}

// Gauss-Kronrod 7/15 abscissae and weights (QUADPACK `qk15`).
const XGK: [Scalar; 8] = [
    0.991_455_371_120_812_6,
    0.949_107_912_342_758_5,
    0.864_864_423_359_769_1,
    0.741_531_185_599_394_4,
    0.586_087_235_467_691_1,
    0.405_845_151_377_397_2,
    0.207_784_955_007_898_5,
    0.0,
];
const WGK: [Scalar; 8] = [
    0.022_935_322_010_529_22,
    0.063_092_092_629_978_55,
    0.104_790_010_322_250_2,
    0.140_653_259_715_525_9,
    0.169_004_726_639_267_9,
    0.190_350_578_064_785_4,
    0.204_432_940_075_298_9,
    0.209_482_141_084_727_8,
];
const WG: [Scalar; 4] = [
    0.129_484_966_168_869_7,
    0.279_705_391_489_276_7,
    0.381_830_050_505_118_9,
    0.417_959_183_673_469_4,
];

/// One Gauss-Kronrod rule over `[a, b]`: the K15 value, its error estimate,
/// and the integral of the integrand's magnitude (the scale that estimate is
/// judged against).
fn kronrod(
    a: Scalar,
    b: Scalar,
    f: &mut dyn FnMut(Scalar) -> Result<Sums, ExactMeasureError>,
) -> Result<(Sums, Sums, Sums), ExactMeasureError> {
    let centre = 0.5 * (a + b);
    let half = 0.5 * (b - a);
    let mut k = [0.0; COMPONENTS];
    let mut g = [0.0; COMPONENTS];
    let mut magnitude = [0.0; COMPONENTS];
    let mut add = |value: Sums, kw: Scalar, gw: Scalar| {
        for c in 0..COMPONENTS {
            k[c] += kw * value[c];
            g[c] += gw * value[c];
            magnitude[c] += kw * value[c].abs();
        }
    };
    add(f(centre)?, WGK[7], WG[3]);
    for j in 0..7 {
        let dx = half * XGK[j];
        // Odd indices are the embedded Gauss nodes.
        let gw = if j % 2 == 1 { WG[j / 2] } else { 0.0 };
        add(f(centre - dx)?, WGK[j], gw);
        add(f(centre + dx)?, WGK[j], gw);
    }
    let mut value = [0.0; COMPONENTS];
    let mut error = [0.0; COMPONENTS];
    let mut scale = [0.0; COMPONENTS];
    for c in 0..COMPONENTS {
        value[c] = half * k[c];
        scale[c] = (half * magnitude[c]).abs();
        let raw = (half * (k[c] - g[c])).abs();
        // QUADPACK's estimate: the raw G7/K15 difference is G7's error, which
        // overstates K15's own by orders of magnitude on smooth integrands.
        error[c] = if scale[c] > 0.0 && raw > 0.0 {
            scale[c] * (200.0 * raw / scale[c]).powf(1.5).min(1.0)
        } else {
            raw
        };
    }
    Ok((value, error, scale))
}

/// Adaptive bisection until every piece meets the relative bound.
fn adaptive(
    a: Scalar,
    b: Scalar,
    floor: &Sums,
    f: &mut dyn FnMut(Scalar) -> Result<Sums, ExactMeasureError>,
) -> Result<Sums, ExactMeasureError> {
    let mut total = [0.0; COMPONENTS];
    if a == b {
        return Ok(total);
    }
    let mut pending = vec![(a, b)];
    let mut pieces = 0;
    while let Some((lo, hi)) = pending.pop() {
        pieces += 1;
        if pieces > MAX_PIECES {
            return Err(ExactMeasureError::NotConverged);
        }
        let (value, error, scale) = kronrod(lo, hi, f)?;
        let accepted = (0..COMPONENTS).all(|c| error[c] <= (RELATIVE * scale[c]).max(floor[c]));
        if accepted {
            for c in 0..COMPONENTS {
                total[c] += value[c];
            }
        } else {
            let mid = 0.5 * (lo + hi);
            if mid <= lo.min(hi) || mid >= lo.max(hi) {
                return Err(ExactMeasureError::NotConverged);
            }
            pending.push((lo, mid));
            pending.push((mid, hi));
        }
    }
    if total.iter().all(|value| value.is_finite()) {
        Ok(total)
    } else {
        Err(ExactMeasureError::NotConverged)
    }
}

#[cfg(test)]
mod tests {
    //! Faces whose domain reaches a pole, or fails to bound anything. No
    //! constructor builds these yet (a revolved circle is still refused), so
    //! they are assembled by hand: one circular edge shared by a curved face
    //! and a planar disc.

    use crate::exact::{exact_properties, ExactMeasureError};
    use axiolid_brep::{ExactBRep, ExactBRepBuilder};
    use axiolid_core::{Frame2, Frame3, Interval, Point3, Tolerance, Vec2, Vec3};
    use axiolid_curve::{Circle2, Circle3, Curve2, Curve3, KnotSpec, Line2};
    use axiolid_surface::{BSplineSurface, Cone, Cylinder, Plane, Sphere, Surface, Torus};
    use axiolid_topology::{
        Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
    };
    use core::f64::consts::{PI, TAU};

    const WORLD: Frame3 = Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    };

    /// A solid bounded by `surface` above (or below) the circle of radius
    /// `r` in `z = 0`, closed by the disc in that plane. `upward` says the
    /// curved face lies above the circle, so its loop runs `+u` and the
    /// disc faces down.
    fn capped(surface: Surface, r: f64, upward: bool) -> ExactBRep {
        let disc = Surface::Plane(Plane { frame: WORLD });
        let round = Circle2 {
            frame: Frame2 {
                origin: Vec2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: r,
        };
        capped_with(surface, r, upward, disc, round)
    }

    /// [`capped`], with the disc on `disc` and its rim at `round` in the
    /// disc's parameters.
    fn capped_with(
        surface: Surface,
        r: f64,
        upward: bool,
        disc: Surface,
        round: Circle2,
    ) -> ExactBRep {
        let mut b = ExactBRepBuilder::default();
        let vertex = b.topology_mut().add_vertex(Vertex {
            position: Point3::new(r, 0.0, 0.0),
        });
        let circle = b.add_curve3(Curve3::Circle(Circle3 {
            frame: WORLD,
            radius: r,
        }));
        let edge = b.topology_mut().add_edge(Edge {
            start: vertex,
            end: vertex,
            curve: Some(circle),
        });
        b.set_edge_interval(edge, Interval::new(0.0, TAU));

        let (wall_use, disc_use, forward, backward) = if upward {
            (
                Orientation::Forward,
                Orientation::Reversed,
                Interval::new(0.0, TAU),
                Interval::new(TAU, 0.0),
            )
        } else {
            (
                Orientation::Reversed,
                Orientation::Forward,
                Interval::new(TAU, 0.0),
                Interval::new(0.0, TAU),
            )
        };

        let rim = b.add_curve2(Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::X,
        }));
        let wall_loop = b.topology_mut().add_loop(Loop {
            edges: vec![EdgeUse {
                edge,
                orientation: wall_use,
                pcurve: Some(rim),
            }],
        });
        b.set_pcurve_interval(wall_loop, 0, forward);

        let round = b.add_curve2(Curve2::Circle(round));
        let disc_loop = b.topology_mut().add_loop(Loop {
            edges: vec![EdgeUse {
                edge,
                orientation: disc_use,
                pcurve: Some(round),
            }],
        });
        b.set_pcurve_interval(disc_loop, 0, backward);

        let wall_surface = b.add_surface(surface);
        let plane = b.add_surface(disc);
        let mut faces = Vec::new();
        for (surface, loop_id) in [(wall_surface, wall_loop), (plane, disc_loop)] {
            faces.push((
                b.topology_mut().add_face(Face {
                    surface: Some(surface),
                    bounds: vec![FaceBound {
                        loop_id,
                        orientation: Orientation::Forward,
                        outer: true,
                    }],
                    orientation: Orientation::Forward,
                }),
                Orientation::Forward,
            ));
        }
        let outer = b.topology_mut().add_shell(Shell {
            faces,
            closed: true,
        });
        b.topology_mut().add_solid(Solid {
            outer,
            voids: Vec::new(),
        });
        b.finish().expect("a valid capped solid")
    }

    fn close(what: &str, got: f64, expected: f64) {
        assert!(
            (got - expected).abs() <= 1e-11 * expected.abs().max(1.0),
            "{what}: expected {expected}, got {got}"
        );
    }

    #[test]
    fn a_hemisphere_reaches_its_north_pole() {
        let r = 1.5;
        let solid = capped(
            Surface::Sphere(Sphere {
                frame: WORLD,
                radius: r,
            }),
            r,
            true,
        );
        let props = exact_properties(&solid, Tolerance::METRE).expect("measurable");
        close("volume", props.signed_volume, 2.0 / 3.0 * PI * r.powi(3));
        close("area", props.area, 3.0 * PI * r * r);
        close("centroid z", props.centroid.z, 3.0 * r / 8.0);
        // int z^2 dV over the upper half ball is (2/15) pi r^5.
        close(
            "second moment z",
            props.second_moment_diagonal.z,
            2.0 / 15.0 * PI * r.powi(5),
        );
    }

    #[test]
    fn a_lower_hemisphere_reaches_its_south_pole() {
        // The loop runs -u, so the domain lies below it: choosing the north
        // pole would measure the complementary half and flip the sign.
        let r = 0.75;
        let solid = capped(
            Surface::Sphere(Sphere {
                frame: WORLD,
                radius: r,
            }),
            r,
            false,
        );
        let props = exact_properties(&solid, Tolerance::METRE).expect("measurable");
        close("volume", props.signed_volume, 2.0 / 3.0 * PI * r.powi(3));
        close("centroid z", props.centroid.z, -3.0 * r / 8.0);
    }

    #[test]
    fn a_cone_reaches_its_apex() {
        // Radius r at v = 0, shrinking to the apex at height h.
        let (r, h) = (2.0, 3.0);
        let solid = capped(
            Surface::Cone(Cone {
                frame: WORLD,
                radius: r,
                semi_angle: (-r / h).atan(),
            }),
            r,
            true,
        );
        let props = exact_properties(&solid, Tolerance::METRE).expect("measurable");
        close("volume", props.signed_volume, PI * r * r * h / 3.0);
        close(
            "area",
            props.area,
            PI * r * r + PI * r * (r * r + h * h).sqrt(),
        );
        close("centroid z", props.centroid.z, h / 4.0);
    }

    #[test]
    fn a_b_spline_face_is_integrated_like_any_other() {
        // The hemisphere again, closed by a bilinear B-spline patch lying in
        // z = 0 instead of a plane: parameters (u, v) in [0, 1]^2 map to
        // (-r + 2 r u, -r + 2 r v, 0), so the rim is the circle of radius
        // 1/2 about (1/2, 1/2) in the patch's parameters.
        let r = 1.25;
        let corner = |x: f64, y: f64| Point3::new(x, y, 0.0);
        let patch = BSplineSurface {
            u_degree: 1,
            v_degree: 1,
            control_points: vec![
                vec![corner(-r, -r), corner(-r, r)],
                vec![corner(r, -r), corner(r, r)],
            ],
            u_knots: vec![0.0, 1.0],
            u_multiplicities: vec![2, 2],
            v_knots: vec![0.0, 1.0],
            v_multiplicities: vec![2, 2],
            weights: None,
            u_closed: false,
            v_closed: false,
            knot_spec: KnotSpec::Unspecified,
            self_intersect: None,
        };
        let rim = Circle2 {
            frame: Frame2 {
                origin: Vec2::new(0.5, 0.5),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 0.5,
        };
        let solid = capped_with(
            Surface::Sphere(Sphere {
                frame: WORLD,
                radius: r,
            }),
            r,
            true,
            Surface::BSpline(patch),
            rim,
        );
        let props = exact_properties(&solid, Tolerance::METRE).expect("measurable");
        close("volume", props.signed_volume, 2.0 / 3.0 * PI * r.powi(3));
        close("area", props.area, 3.0 * PI * r * r);
    }

    /// Half a torus, `0 <= u <= pi`, closed by the two meridian discs. The
    /// torus face's loops are the meridians, which wind round the TUBE (in
    /// `v`), so this is the face integrated in `u` first.
    fn half_torus(major: f64, minor: f64) -> ExactBRep {
        let mut b = ExactBRepBuilder::default();
        let meridian = |b: &mut ExactBRepBuilder, x: f64| {
            let frame = Frame3 {
                origin: Point3::new(x * major, 0.0, 0.0),
                x: Vec3::X * x,
                y: Vec3::Z,
                z: (Vec3::X * x).cross(Vec3::Z),
            };
            let vertex = b.topology_mut().add_vertex(Vertex {
                position: Point3::new(x * (major + minor), 0.0, 0.0),
            });
            let circle = b.add_curve3(Curve3::Circle(Circle3 {
                frame,
                radius: minor,
            }));
            let edge = b.topology_mut().add_edge(Edge {
                start: vertex,
                end: vertex,
                curve: Some(circle),
            });
            b.set_edge_interval(edge, Interval::new(0.0, TAU));
            edge
        };
        let start = meridian(&mut b, 1.0);
        let end = meridian(&mut b, -1.0);

        let up = |b: &mut ExactBRepBuilder, u: f64| {
            b.add_curve2(Curve2::Line(Line2 {
                origin: Vec2::new(u, 0.0),
                direction: Vec2::Y,
            }))
        };
        let (at_start, at_end) = (up(&mut b, 0.0), up(&mut b, PI));
        // Anticlockwise in (u, v): up the u = pi side, down the u = 0 side.
        let tube_end = b.topology_mut().add_loop(Loop {
            edges: vec![EdgeUse {
                edge: end,
                orientation: Orientation::Forward,
                pcurve: Some(at_end),
            }],
        });
        b.set_pcurve_interval(tube_end, 0, Interval::new(0.0, TAU));
        // Stored running UP the u = 0 side and bound Reversed, so the face
        // walks it down: a reversed bound must reverse the pcurve too.
        let tube_start = b.topology_mut().add_loop(Loop {
            edges: vec![EdgeUse {
                edge: start,
                orientation: Orientation::Forward,
                pcurve: Some(at_start),
            }],
        });
        b.set_pcurve_interval(tube_start, 0, Interval::new(0.0, TAU));

        // Both end discs face -y. At u = 0 the plane's (x, y) are world
        // (x, z), which the meridian runs anticlockwise; at u = pi the
        // meridian's own x is world -x, so it runs clockwise and the disc
        // walks it backwards.
        let disc_frame = |x: f64| Frame3 {
            origin: Point3::new(x * major, 0.0, 0.0),
            x: Vec3::X,
            y: Vec3::Z,
            z: -Vec3::Y,
        };
        let rim = |b: &mut ExactBRepBuilder, x: f64| {
            b.add_curve2(Curve2::Circle(Circle2 {
                frame: Frame2 {
                    origin: Vec2::ZERO,
                    x: Vec2::new(x, 0.0),
                    y: Vec2::Y,
                },
                radius: minor,
            }))
        };
        let (rim_start, rim_end) = (rim(&mut b, 1.0), rim(&mut b, -1.0));
        let disc_start = b.topology_mut().add_loop(Loop {
            edges: vec![EdgeUse {
                edge: start,
                orientation: Orientation::Forward,
                pcurve: Some(rim_start),
            }],
        });
        b.set_pcurve_interval(disc_start, 0, Interval::new(0.0, TAU));
        let disc_end = b.topology_mut().add_loop(Loop {
            edges: vec![EdgeUse {
                edge: end,
                orientation: Orientation::Reversed,
                pcurve: Some(rim_end),
            }],
        });
        b.set_pcurve_interval(disc_end, 0, Interval::new(TAU, 0.0));

        let torus = b.add_surface(Surface::Torus(Torus {
            frame: WORLD,
            major_radius: major,
            minor_radius: minor,
        }));
        let plane_start = b.add_surface(Surface::Plane(Plane {
            frame: disc_frame(1.0),
        }));
        let plane_end = b.add_surface(Surface::Plane(Plane {
            frame: disc_frame(-1.0),
        }));
        let bound = |loop_id, outer| FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer,
        };
        let reversed = FaceBound {
            loop_id: tube_start,
            orientation: Orientation::Reversed,
            outer: false,
        };
        let mut faces = Vec::new();
        for (surface, bounds) in [
            (torus, vec![bound(tube_end, true), reversed]),
            (plane_start, vec![bound(disc_start, true)]),
            (plane_end, vec![bound(disc_end, true)]),
        ] {
            let face = b.topology_mut().add_face(Face {
                surface: Some(surface),
                bounds,
                orientation: Orientation::Forward,
            });
            faces.push((face, Orientation::Forward));
        }
        let outer = b.topology_mut().add_shell(Shell {
            faces,
            closed: true,
        });
        b.topology_mut().add_solid(Solid {
            outer,
            voids: Vec::new(),
        });
        b.finish().expect("a valid half torus")
    }

    #[test]
    fn a_face_winding_round_the_tube_is_integrated_along_it() {
        let (major, minor) = (3.0, 1.0);
        let props =
            exact_properties(&half_torus(major, minor), Tolerance::METRE).expect("measurable");
        close(
            "volume",
            props.signed_volume,
            PI * PI * major * minor * minor,
        );
        close(
            "area",
            props.area,
            2.0 * PI * PI * major * minor + 2.0 * PI * minor * minor,
        );
        // int y dV = int_0^pi sin u du * int rho^2 dA over the tube section.
        let moment = 2.0 * (PI * minor * minor * major * major + PI * minor.powi(4) / 4.0);
        close(
            "centroid y",
            props.centroid.y,
            moment / (PI * PI * major * minor * minor),
        );
    }

    #[test]
    fn the_adaptive_rule_meets_its_bound_where_one_panel_cannot() {
        // Every elementary face integrand is a low-order trigonometric
        // polynomial that a single K15 panel already integrates exactly, so
        // the solids never exercise the error bound. Runge's function does:
        // one panel is off in the fourth digit.
        use super::{adaptive, kronrod, COMPONENTS};
        let runge = |x: f64| 1.0 / (1.0 + 100.0 * x * x);
        let exact = 2.0 * 10.0_f64.atan() / 10.0;
        let floor = [0.0; COMPONENTS];
        let mut f = |x: f64| Ok([runge(x); COMPONENTS]);
        let (one_panel, _, _) = kronrod(-1.0, 1.0, &mut f).expect("finite");
        assert!((one_panel[0] - exact).abs() > 1e-6, "{}", one_panel[0]);
        let value = adaptive(-1.0, 1.0, &floor, &mut f).expect("converges");
        for component in value {
            assert!(
                (component - exact).abs() < 1e-14,
                "expected {exact}, got {component}"
            );
        }
    }

    #[test]
    fn a_lone_circle_on_a_cylinder_bounds_nothing() {
        // A cylinder has no pole: one loop round it leaves the domain open
        // towards infinity. Measuring it would invent a top.
        let solid = capped(
            Surface::Cylinder(Cylinder {
                frame: WORLD,
                radius: 1.0,
            }),
            1.0,
            true,
        );
        let error = exact_properties(&solid, Tolerance::METRE).expect_err("unbounded");
        assert!(
            matches!(error, ExactMeasureError::ParameterDomain(_)),
            "got {error:?}"
        );
    }
}
