//! Faces that are translates of each other, whatever chart their B-rep
//! trims them in (#227).
//!
//! The matched bound of [`super`] needs one parameter domain standing for a
//! face on both surfaces. Faces with the very same trims have that. A face
//! trimmed in world coordinates does not: a cap built on the world `xy`
//! frame and moved by `t` keeps its plane's axes, while its pcurves move by
//! `(t . x, t . y)`. It is still one domain, shifted.
//!
//! # The shifted chart
//!
//! Faces `F_A` on `S_A` and `F_B` on `S_B` are matched with a parameter
//! shift `s` when the loops of `F_B` are those of `F_A` moved by `s`, use by
//! use and parameter by parameter, to within a measured `delta`:
//! `|c_B(t) - (c_A(t) + s)| <= delta` for every pcurve pair over its
//! interval (the intervals equal). `S_B` is then re-charted onto `A`'s
//! parameters as `S'(q) = S_B(q + s)`, a surface of the same family: for a
//! plane the origin moves by `s_u x + s_v y`; for a cylinder, an elliptical
//! cylinder or a cone `s` runs along the axis only (`s_u = 0`, a turn about
//! the axis is a rotation, not a translation), the origin moves by `s_v z`,
//! and a cone's radius at `v = 0` becomes `r + s_v tan(alpha)`.
//!
//! For a point `S_A(q)` with `q` in `A`'s domain, either `q + s` lies in
//! `B`'s domain, and `S'(q) = S_B(q + s)` is a point of `B`'s boundary, or it
//! does not. Both domains are bounded by the same pieces in the same order
//! (seams and pole stretches are joined alike, as `delta` is far below the
//! joins' slack), so the straight homotopy between `B`'s loops and `A`'s
//! moved by `s` moves no boundary point farther than `delta`: a point
//! farther than `delta` from `B`'s boundary curves has the same winding in
//! both and is in one domain exactly when it is in the other. So `q + s`
//! is then within `delta` of a point `q'` of `B`'s boundary curves, whose
//! image is on `B`'s boundary, and `|S'(q) - S_B(q')| <= L delta` with `L` a
//! Lipschitz constant of `S'` over the patch grown by `delta`. Hence over a
//! patch of `A`'s domain
//!
//! `d(S_A(q), dB) <= |S_A - S'| + L delta + rounding`,
//!
//! with `|S_A - S'|` bounded as for any matched pair (same family, same
//! basis), and the rounding of building `S'` added.
//!
//! # A translation gives `|t|`
//!
//! When `F_B = F_A + t`, the re-charted surface is `S'(q) = S_A(q) + t`: for
//! a plane, `S_B(q + s) = o_B + s_u x + s_v y + u x + v y` with the same
//! axes, and `o_B + s_u x + s_v y = o_A + t` is how `B`'s builder placed the
//! translate of `A`'s origin; for a cylinder or cone the shift along the
//! axis does the same. Every basis coefficient agrees, the spread is zero,
//! and the bound is `|o_A - o_S'| = |t|` on every patch at once, up to
//! rounding: exact for the farthest point of a convex solid moved by `t`.
//!
//! # What is never matched
//!
//! The re-charted bound is sound for any two surfaces of one family, but
//! only a translate is matched here: the axes must agree, and so must the
//! shape (radius, semi-axes, angle, and a cone's apex after the shift), to
//! within [`GATE`] relative to the larger of one and the quantity; the
//! trims must agree within `GATE` relative to their size. A turned, resized
//! or re-trimmed face fails one of these and keeps the Lipschitz bound,
//! however small the difference. What the gate admits is still folded in:
//! the residual axis or radius difference into `|S_A - S'|` through the
//! coefficients, the trim residual through `L delta`. Spheres, tori and
//! B-splines have no translation that shifts their parameters; they match
//! only with identical trims (see [`super`]).

use axiolid_brep::ExactBRep;
use axiolid_core::{Frame3, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{Curve2, Curve3};
use axiolid_evaluate::evaluate3;
use axiolid_surface::{Cone, Cylinder, EllipticalCylinder, Plane, Surface};
use axiolid_topology::Face;

use crate::exact_distance::surface_of;

/// The largest axis, shape or trim difference, relative to the larger of
/// one and the quantity, that is still read as a translate. Rounding in a
/// builder that recomputes a translate from moved points stays many orders
/// below; a turn or a resize by a millionth is far above.
pub(super) const GATE: Scalar = 1e-9;

/// `B`'s face re-charted onto `A`'s parameters (see the module docs).
#[derive(Debug, Clone)]
pub(super) struct Shifted {
    /// `S'(q) = S_B(q + s)`.
    pub(super) surface: Surface,
    /// How far `B`'s trims are from `A`'s moved by `s`, in parameters.
    pub(super) delta: Scalar,
    /// Rounding in building `surface`, in length.
    pub(super) slack: Scalar,
}

impl Shifted {
    /// The bound `|S_A - S'| + L delta + slack` over the patch `[lo, hi]`,
    /// given `|S_A - S'|` from the matched bound.
    pub(super) fn bound(&self, matched: Scalar, lo: Point2, hi: Point2) -> Option<Scalar> {
        let lipschitz = lipschitz(&self.surface, lo, hi, self.delta)?;
        let bound = matched + lipschitz * self.delta + self.slack;
        bound.is_finite().then_some(bound)
    }
}

/// A Lipschitz constant of `surface` over `[lo, hi]` grown by `grow`: each
/// basis function of the planar, cylindrical and conical families has a
/// gradient of length at most `1 + |v|`.
fn lipschitz(surface: &Surface, lo: Point2, hi: Point2, grow: Scalar) -> Option<Scalar> {
    let (_, coefficients) = super::terms(surface)?;
    let reach_v = lo.y.abs().max(hi.y.abs()) + grow;
    let sum: Scalar = coefficients.iter().map(|a| a.length()).sum();
    Some(sum * (1.0 + reach_v))
}

/// Whether two quantities agree to within [`GATE`] relative to the larger
/// of one and either.
fn close(a: Scalar, b: Scalar) -> bool {
    (a - b).abs() <= GATE * a.abs().max(b.abs()).max(1.0)
}

fn same_axes(a: &Frame3, b: &Frame3) -> bool {
    [(a.x, b.x), (a.y, b.y), (a.z, b.z)]
        .iter()
        .all(|(p, q)| (*p - *q).length() <= GATE)
}

/// Which parameters a translation may shift, by family, or `None` for a
/// family or shape that is not a translate's.
fn shiftable(a: &Surface, b: &Surface) -> Option<[bool; 2]> {
    match (a, b) {
        (Surface::Plane(p), Surface::Plane(q)) => {
            same_axes(&p.frame, &q.frame).then_some([true, true])
        }
        (Surface::Cylinder(p), Surface::Cylinder(q)) => {
            (same_axes(&p.frame, &q.frame) && close(p.radius, q.radius)).then_some([false, true])
        }
        (Surface::EllipticalCylinder(p), Surface::EllipticalCylinder(q)) => {
            (same_axes(&p.frame, &q.frame)
                && close(p.semi_axis_x, q.semi_axis_x)
                && close(p.semi_axis_y, q.semi_axis_y))
            .then_some([false, true])
        }
        // The radius is compared after the shift, which moves it.
        (Surface::Cone(p), Surface::Cone(q)) => (same_axes(&p.frame, &q.frame)
            && close(p.semi_angle, q.semi_angle))
        .then_some([false, true]),
        _ => None,
    }
}

/// `S_B(q + s)` as a surface of `S_B`'s family, and the rounding of
/// building it.
fn recharted(b: &Surface, s: Vec2) -> Option<(Surface, Scalar)> {
    let moved = |f: &Frame3, by: Vec3| Frame3 {
        origin: f.origin + by,
        ..*f
    };
    let rounding = |f: &Frame3, extra: Scalar| {
        let reach = f.origin.length()
            + (s.x.abs() + s.y.abs()) * (f.x.length() + f.y.length() + f.z.length())
            + extra;
        1e-12 * reach
    };
    Some(match b {
        Surface::Plane(p) => (
            Surface::Plane(Plane {
                frame: moved(&p.frame, p.frame.x * s.x + p.frame.y * s.y),
            }),
            rounding(&p.frame, 0.0),
        ),
        Surface::Cylinder(c) => (
            Surface::Cylinder(Cylinder {
                frame: moved(&c.frame, c.frame.z * s.y),
                radius: c.radius,
            }),
            rounding(&c.frame, 0.0),
        ),
        Surface::EllipticalCylinder(c) => (
            Surface::EllipticalCylinder(EllipticalCylinder {
                frame: moved(&c.frame, c.frame.z * s.y),
                semi_axis_x: c.semi_axis_x,
                semi_axis_y: c.semi_axis_y,
            }),
            rounding(&c.frame, 0.0),
        ),
        Surface::Cone(c) => {
            let slope = c.semi_angle.tan();
            let radius = c.radius + s.y * slope;
            (
                Surface::Cone(Cone {
                    frame: moved(&c.frame, c.frame.z * s.y),
                    radius,
                    semi_angle: c.semi_angle,
                }),
                rounding(
                    &c.frame,
                    c.radius.abs() + radius.abs() + s.y.abs() * slope.abs(),
                ),
            )
        }
        _ => return None,
    })
}

/// The point a pcurve's shift is read from.
fn anchor(curve: &Curve2) -> Option<Point2> {
    match curve {
        Curve2::Line(l) => Some(l.origin),
        Curve2::Circle(c) => Some(c.frame.origin),
        Curve2::Ellipse(e) => Some(e.frame.origin),
        Curve2::Sinusoid(w) => Some(Point2::new(0.0, w.mean)),
        _ => None,
    }
}

/// `max over t in [t0, t1] of |c_B(t) - (c_A(t) + s)|` bounded above, and
/// the size of the pcurves for the gate; `None` for curves of different
/// kinds or a kind not compared.
fn deviation(a: &Curve2, b: &Curve2, s: Vec2, t0: Scalar, t1: Scalar) -> Option<(Scalar, Scalar)> {
    let reach_t = t0.abs().max(t1.abs());
    let (gap, size) = match (a, b) {
        (Curve2::Line(p), Curve2::Line(q)) => (
            (q.origin - p.origin - s).length() + reach_t * (q.direction - p.direction).length(),
            p.origin.length() + reach_t * p.direction.length(),
        ),
        (Curve2::Circle(p), Curve2::Circle(q)) => (
            (q.frame.origin - p.frame.origin - s).length()
                + (q.frame.x * q.radius - p.frame.x * p.radius).length()
                + (q.frame.y * q.radius - p.frame.y * p.radius).length(),
            p.frame.origin.length() + p.radius.abs(),
        ),
        (Curve2::Ellipse(p), Curve2::Ellipse(q)) => (
            (q.frame.origin - p.frame.origin - s).length()
                + (q.frame.x * q.semi_axis_x - p.frame.x * p.semi_axis_x).length()
                + (q.frame.y * q.semi_axis_y - p.frame.y * p.semi_axis_y).length(),
            p.frame.origin.length() + p.semi_axis_x.abs().max(p.semi_axis_y.abs()),
        ),
        // The graph `(t, mean + a cos t + b sin t)`.
        (Curve2::Sinusoid(p), Curve2::Sinusoid(q)) => (
            s.x.abs()
                + (q.mean - p.mean - s.y).abs()
                + (q.cosine - p.cosine).abs()
                + (q.sine - p.sine).abs(),
            reach_t + p.mean.abs() + p.cosine.abs() + p.sine.abs(),
        ),
        _ => return None,
    };
    // The difference itself is rounded, relative to what it is made of.
    let gap = gap + 1e-12 * (size + s.length());
    gap.is_finite().then_some((gap, size))
}

/// The point of `brep`'s edges farthest along `direction`, over straight,
/// circular and elliptic edges: a point certainly on the boundary, at a
/// line's end or where a conic's tangent is normal to `direction`.
///
/// Measured as a witness, it is what closes a translate's interval from
/// below. Every point of `A + t` lies at least `|t|` beyond the support
/// plane of `A` against `t`, so `A`'s support point `p` there is `|t|` from
/// all of `A + t`, whatever the shape: `h(dA, d(A + t)) = |t|`, and the
/// matched bound already holds every patch at `|t|`. On a face bounded by
/// lines and conics whose surface is ruled along its trims (planes,
/// cylinders, cones), a linear function is extreme on the trims, so the
/// support point of such a solid is on an edge.
///
/// Only the edges indexed by `edges` are searched: one item's, when the
/// B-rep holds several (#229).
pub(super) fn support_point(
    brep: &ExactBRep,
    direction: Vec3,
    edges: core::ops::Range<usize>,
) -> Option<(Point3, usize)> {
    let topology = brep.topology();
    let mut best: Option<(Scalar, Point3, usize)> = None;
    for (index, edge) in topology.edges().iter().enumerate() {
        if !edges.contains(&index) {
            continue;
        }
        let Some(curve) = edge.curve.and_then(|id| brep.curves3().get(id.index())) else {
            continue;
        };
        let Some(span) = topology
            .edge_id_at(index)
            .and_then(|id| brep.edge_interval(id))
        else {
            continue;
        };
        let (lo, hi) = (span.start.min(span.end), span.start.max(span.end));
        let mut at = vec![lo, hi];
        let conic = match curve {
            Curve3::Circle(c) => Some((&c.frame, c.radius, c.radius)),
            Curve3::Ellipse(e) => Some((&e.frame, e.semi_axis_x, e.semi_axis_y)),
            Curve3::Line(_) => None,
            _ => continue,
        };
        if let Some((frame, rx, ry)) = conic {
            // `direction . C(t) = k + A cos t + B sin t`, largest at
            // `atan2(B, A)` plus whole turns.
            let peak = (ry * direction.dot(frame.y)).atan2(rx * direction.dot(frame.x));
            let turn = std::f64::consts::TAU;
            let inside = peak + turn * ((lo - peak) / turn).ceil();
            if inside <= hi {
                at.push(inside);
            }
        }
        for t in at {
            let Ok(point) = evaluate3(curve, t) else {
                continue;
            };
            let height = direction.dot(point);
            if point.is_finite() && best.is_none_or(|(top, ..)| height > top) {
                best = Some((height, point, index));
            }
        }
    }
    best.map(|(_, point, edge)| (point, edge))
}

/// `B`'s face `fb` as a translate of `A`'s face `fa`, re-charted onto
/// `fa`'s parameters, or `None` when it is not one (see the module docs).
pub(super) fn translate(
    a: &ExactBRep,
    fa: &Face<axiolid_brep::SurfaceId>,
    b: &ExactBRep,
    fb: &Face<axiolid_brep::SurfaceId>,
) -> Option<Shifted> {
    let (Ok(sa), Ok(sb)) = (surface_of(a, fa.surface), surface_of(b, fb.surface)) else {
        return None;
    };
    let mask = shiftable(sa, sb)?;
    if fa.bounds.len() != fb.bounds.len() {
        return None;
    }
    let (ta, tb) = (a.topology(), b.topology());
    // Every pcurve pair over its interval, loop by loop and use by use.
    let mut pairs = Vec::new();
    for (ba, bb) in fa.bounds.iter().zip(&fb.bounds) {
        if ba.orientation != bb.orientation || ba.outer != bb.outer {
            return None;
        }
        let la = ta.loops().get(ba.loop_id.index())?;
        let lb = tb.loops().get(bb.loop_id.index())?;
        if la.edges.len() != lb.edges.len() {
            return None;
        }
        for (i, (ua, ub)) in la.edges.iter().zip(&lb.edges).enumerate() {
            let ca = super::pcurve(a, ua.pcurve)?;
            let cb = super::pcurve(b, ub.pcurve)?;
            let ia = a.pcurve_interval(ba.loop_id, i)?;
            if b.pcurve_interval(bb.loop_id, i)? != ia {
                return None;
            }
            pairs.push((ca, cb, ia));
        }
    }
    let &(first_a, first_b, _) = pairs.first()?;
    let raw = anchor(first_b)? - anchor(first_a)?;
    let s = Vec2::new(
        if mask[0] { raw.x } else { 0.0 },
        if mask[1] { raw.y } else { 0.0 },
    );
    let mut delta: Scalar = 0.0;
    let mut size: Scalar = 0.0;
    for (ca, cb, interval) in pairs {
        let (gap, scale) = deviation(ca, cb, s, interval.start, interval.end)?;
        delta = delta.max(gap);
        size = size.max(scale);
    }
    let (surface, slack) = recharted(sb, s)?;
    // A cone's apex is a pole of its chart, which a domain may reach: after
    // the shift it must sit where `A`'s does, and the residue is a trim
    // residue.
    if let (Surface::Cone(p), Surface::Cone(q)) = (sa, &surface) {
        if !close(p.radius, q.radius) {
            return None;
        }
        let (pa, pb) = (p.semi_angle.tan(), q.semi_angle.tan());
        if pa != 0.0 || pb != 0.0 {
            let pole = (-p.radius / pa - -q.radius / pb).abs();
            delta = delta.max(pole + 1e-12 * (p.radius / pa).abs());
        }
    }
    (delta.is_finite() && delta <= GATE * (1.0 + size + s.length())).then_some(Shifted {
        surface,
        delta,
        slack,
    })
}

#[cfg(test)]
mod tests {
    //! Which faces are translates, and the re-charted bound checked by
    //! dense sampling of `S_A(q)` against `S_B(q + s)`.

    use super::{translate, Shifted, GATE};
    use axiolid_brep::ExactBRep;
    use axiolid_core::{Frame2, Frame3, Interval, Point2, Point3, Vec2, Vec3};
    use axiolid_curve::{Circle2, Curve2, Line2};
    use axiolid_evaluate::surface::evaluate;
    use axiolid_surface::{Cone, Cylinder, Plane, Surface};

    /// A one-face sheet on `surface`, its loop of uses carrying `pcurves`
    /// over the unit interval.
    fn sheet(surface: Surface, pcurves: &[Curve2]) -> ExactBRep {
        use axiolid_topology::{Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Vertex};
        let mut builder = axiolid_brep::ExactBRepBuilder::default();
        let curve3 = builder.add_curve3(axiolid_curve::Curve3::Line(axiolid_curve::Line3 {
            origin: Point3::ZERO,
            direction: Vec3::X,
        }));
        let surface = builder.add_surface(surface);
        let ids: Vec<_> = pcurves
            .iter()
            .map(|c| builder.add_curve2(c.clone()))
            .collect();
        let count = pcurves.len();
        let topology = builder.topology_mut();
        let vertices: Vec<_> = (0..count)
            .map(|i| {
                topology.add_vertex(Vertex {
                    position: Point3::new(i as f64, 0.0, 0.0),
                })
            })
            .collect();
        let edges: Vec<_> = (0..count)
            .map(|i| {
                topology.add_edge(Edge {
                    start: vertices[i],
                    end: vertices[(i + 1) % count],
                    curve: Some(curve3),
                })
            })
            .collect();
        let loop_id = topology.add_loop(Loop {
            edges: edges
                .iter()
                .zip(&ids)
                .map(|(&edge, &pcurve)| EdgeUse {
                    edge,
                    orientation: Orientation::Forward,
                    pcurve: Some(pcurve),
                })
                .collect(),
        });
        topology.add_face(Face {
            surface: Some(surface),
            bounds: vec![FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation: Orientation::Forward,
        });
        for edge in edges {
            builder.set_edge_interval(edge, Interval::UNIT);
        }
        for use_index in 0..count {
            builder.set_pcurve_interval(loop_id, use_index, Interval::UNIT);
        }
        builder.finish().expect("a valid sheet")
    }

    fn matched(a: &ExactBRep, b: &ExactBRep) -> Option<Shifted> {
        translate(a, &a.topology().faces()[0], b, &b.topology().faces()[0])
    }

    /// A polygon's edges as line pcurves.
    fn polygon(points: &[Point2]) -> Vec<Curve2> {
        (0..points.len())
            .map(|i| {
                let (from, to) = (points[i], points[(i + 1) % points.len()]);
                Curve2::Line(Line2 {
                    origin: from,
                    direction: to - from,
                })
            })
            .collect()
    }

    fn moved(points: &[Point2], by: Vec2) -> Vec<Point2> {
        points.iter().map(|p| *p + by).collect()
    }

    fn plane(origin: Point3, turn: f64) -> Surface {
        let (s, c) = turn.sin_cos();
        Surface::Plane(Plane {
            frame: Frame3 {
                origin,
                x: Vec3::new(c, s, 0.0),
                y: Vec3::new(-s, c, 0.0),
                z: Vec3::Z,
            },
        })
    }

    fn upright(origin: Point3) -> Frame3 {
        Frame3 {
            origin,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        }
    }

    fn cylinder(origin: Point3, radius: f64) -> Surface {
        Surface::Cylinder(Cylinder {
            frame: upright(origin),
            radius,
        })
    }

    /// The largest `|S_A(q) - S_B(q + s)|` over a grid of `[lo, hi]`, the
    /// re-charted surface checked to be `S_B(q + s)` on the way.
    fn reach(a: &Surface, b: &Surface, shifted: &Shifted, s: Vec2, lo: Point2, hi: Point2) -> f64 {
        const N: usize = 16;
        let mut reach: f64 = 0.0;
        for i in 0..=N {
            for j in 0..=N {
                let q = Point2::new(
                    lo.x + (hi.x - lo.x) * i as f64 / N as f64,
                    lo.y + (hi.y - lo.y) * j as f64 / N as f64,
                );
                let pa = evaluate(a, q.x, q.y).expect("point");
                let pb = evaluate(b, q.x + s.x, q.y + s.y).expect("point");
                let recharted = evaluate(&shifted.surface, q.x, q.y).expect("point");
                assert!((recharted - pb).length() <= 1e-12 * (1.0 + pb.length()));
                reach = reach.max((pa - pb).length());
            }
        }
        reach
    }

    fn bound(a: &Surface, shifted: &Shifted, lo: Point2, hi: Point2) -> f64 {
        let matched = super::super::matched_bound(a, &shifted.surface, lo, hi).expect("bounded");
        shifted.bound(matched, lo, hi).expect("bounded")
    }

    const SECTION: [Point2; 3] = [
        Point2::new(1.0, 1.0),
        Point2::new(2.5, 1.25),
        Point2::new(1.5, 2.0),
    ];

    /// A cap on the world `xy` frame moved by `t`: its pcurves move by the
    /// plan part of `t`, its plane by the rest, and the bound is `|t|`.
    #[test]
    fn a_cap_trimmed_in_world_coordinates_is_a_translate_by_t() {
        let t = Vec3::new(0.3, -0.4, 1.2);
        let a = plane(Point3::new(0.0, 0.0, 2.0), 0.0);
        let b = plane(Point3::new(0.0, 0.0, 2.0 + t.z), 0.0);
        let s = Vec2::new(t.x, t.y);
        let sa = sheet(a.clone(), &polygon(&SECTION));
        let sb = sheet(b.clone(), &polygon(&moved(&SECTION, s)));
        let shifted = matched(&sa, &sb).expect("a translate");
        // Only the rounding of the comparison itself.
        assert!(shifted.delta <= 1e-10, "{}", shifted.delta);
        let (lo, hi) = (Point2::new(1.0, 1.0), Point2::new(2.5, 2.0));
        let bound = bound(&a, &shifted, lo, hi);
        let reach = reach(&a, &b, &shifted, s, lo, hi);
        assert!(reach <= bound, "{reach} > {bound}");
        assert!((bound - t.length()).abs() < 1e-9, "{bound} is not |t|");
    }

    /// The trims may differ by rounding: the residue is measured and folded
    /// in, and a residue above the gate is no translate.
    #[test]
    fn a_trim_residue_is_folded_in_and_a_larger_one_refused() {
        let s = Vec2::new(0.25, -0.5);
        let sa = sheet(plane(Point3::ZERO, 0.0), &polygon(&SECTION));
        for (nudge, accepted) in [(1e-13, true), (1e-6, false)] {
            let mut points = moved(&SECTION, s);
            points[1].x += nudge;
            let sb = sheet(plane(Point3::ZERO, 0.0), &polygon(&points));
            let shifted = matched(&sa, &sb);
            assert_eq!(shifted.is_some(), accepted, "nudge {nudge}");
            if let Some(shifted) = shifted {
                assert!(shifted.delta >= nudge, "{} < {nudge}", shifted.delta);
                assert!(shifted.delta <= GATE);
            }
        }
    }

    #[test]
    fn a_turned_cap_is_never_a_translate() {
        let sa = sheet(plane(Point3::ZERO, 0.0), &polygon(&SECTION));
        for turn in [1e-6_f64, 1e-3] {
            // Turned about the world origin: the plane's axes and the trim
            // turn together, a rigid motion but not a translation.
            let (sin, cos) = turn.sin_cos();
            let points: Vec<_> = SECTION
                .iter()
                .map(|p| Point2::new(cos * p.x - sin * p.y, sin * p.x + cos * p.y))
                .collect();
            let sb = sheet(plane(Point3::ZERO, turn), &polygon(&points));
            assert!(matched(&sa, &sb).is_none(), "turn {turn}");
            // The same turned trim on the unturned plane is another face.
            let sb = sheet(plane(Point3::ZERO, 0.0), &polygon(&points));
            assert!(matched(&sa, &sb).is_none(), "trim turned by {turn}");
            // The turned plane with the trim moved as a translate's would
            // be: the axes alone tell it apart.
            let sb = sheet(
                plane(Point3::ZERO, turn),
                &polygon(&moved(&SECTION, Vec2::new(0.5, 0.25))),
            );
            assert!(matched(&sa, &sb).is_none(), "plane turned by {turn}");
        }
    }

    /// A wall on a cylinder trimmed at world heights: moved by `t`, its
    /// axis moves by `t` across and its trims by `t.z` in `v`.
    fn wall(axis: Point3, radius: f64, lift: f64, turn: f64) -> ExactBRep {
        let (u0, u1) = (0.25 + turn, 1.75 + turn);
        let (v0, v1) = (3.0 + lift, 4.5 + lift);
        sheet(
            cylinder(axis, radius),
            &polygon(&[
                Point2::new(u0, v0),
                Point2::new(u1, v0),
                Point2::new(u1, v1),
                Point2::new(u0, v1),
            ]),
        )
    }

    #[test]
    fn a_cylinder_wall_moved_along_and_across_its_axis_is_a_translate() {
        let t = Vec3::new(0.1, -0.2, 0.3);
        let axis = Point3::new(1.0, 2.0, 0.0);
        let across = axis + Vec3::new(t.x, t.y, 0.0);
        let shifted =
            matched(&wall(axis, 0.5, 0.0, 0.0), &wall(across, 0.5, t.z, 0.0)).expect("a translate");
        let (sa, sb) = (cylinder(axis, 0.5), cylinder(across, 0.5));
        let (lo, hi) = (Point2::new(0.25, 3.0), Point2::new(1.75, 4.5));
        let bound = bound(&sa, &shifted, lo, hi);
        let reach = reach(&sa, &sb, &shifted, Vec2::new(0.0, t.z), lo, hi);
        assert!(reach <= bound, "{reach} > {bound}");
        assert!((bound - t.length()).abs() < 1e-9, "{bound} is not |t|");
    }

    #[test]
    fn a_resized_or_turned_cylinder_wall_is_never_a_translate() {
        let axis = Point3::new(1.0, 2.0, 0.0);
        let a = wall(axis, 0.5, 0.0, 0.0);
        // Another radius, by a millionth.
        assert!(matched(&a, &wall(axis, 0.5 + 1e-6, 0.0, 0.0)).is_none());
        // Turned about the axis: the trim moves in `u`.
        assert!(matched(&a, &wall(axis, 0.5, 0.0, 1e-6)).is_none());
        // Within rounding of the radius, still a translate.
        assert!(matched(&a, &wall(axis, 0.5 + 1e-15, 0.2, 0.0)).is_some());
    }

    #[test]
    fn a_cone_moved_along_its_axis_is_recharted_to_its_own_apex() {
        let (radius, angle) = (0.75_f64, 0.3_f64);
        let cone = |origin: Point3, radius: f64| {
            Surface::Cone(Cone {
                frame: upright(origin),
                radius,
                semi_angle: angle,
            })
        };
        let trim = |lift: f64| {
            polygon(&[
                Point2::new(0.5, 1.0 + lift),
                Point2::new(2.0, 1.0 + lift),
                Point2::new(2.0, 2.0 + lift),
                Point2::new(0.5, 2.0 + lift),
            ])
        };
        // `B` is `A` moved by `t`, its frame placed `d` further up the
        // axis: radius `r + d tan(angle)` there, and its trims `d` lower.
        let t = Vec3::new(0.2, 0.1, -0.4);
        let d = 0.6;
        let sa = cone(Point3::ZERO, radius);
        let sb = cone(Point3::ZERO + t + Vec3::Z * d, radius + d * angle.tan());
        let a = sheet(sa.clone(), &trim(0.0));
        let shifted = matched(&a, &sheet(sb.clone(), &trim(-d))).expect("a translate");
        let (lo, hi) = (Point2::new(0.5, 1.0), Point2::new(2.0, 2.0));
        let bound = bound(&sa, &shifted, lo, hi);
        let reach = reach(&sa, &sb, &shifted, Vec2::new(0.0, -d), lo, hi);
        assert!(reach <= bound, "{reach} > {bound}");
        assert!((bound - t.length()).abs() < 1e-9, "{bound} is not |t|");
        // A radius off by a millionth puts the apex elsewhere.
        let sc = cone(
            Point3::ZERO + t + Vec3::Z * d,
            radius + d * angle.tan() + 1e-6,
        );
        assert!(matched(&a, &sheet(sc, &trim(-d))).is_none());
        // Off by less than the gate, it is matched with the apex's move
        // folded into the trim residue.
        let nudge = 1e-10;
        let sc = cone(
            Point3::ZERO + t + Vec3::Z * d,
            radius + d * angle.tan() + nudge,
        );
        let shifted = matched(&a, &sheet(sc, &trim(-d))).expect("within the gate");
        assert!(
            shifted.delta >= 0.9 * nudge / angle.tan(),
            "{}",
            shifted.delta
        );
    }

    /// A cone of no angle has no apex: its radius is the only check.
    #[test]
    fn a_resized_cone_without_an_apex_is_never_a_translate() {
        let cone = |radius: f64| {
            Surface::Cone(Cone {
                frame: upright(Point3::ZERO),
                radius,
                semi_angle: 0.0,
            })
        };
        let trim = polygon(&[
            Point2::new(0.5, 1.0),
            Point2::new(2.0, 1.0),
            Point2::new(2.0, 2.0),
            Point2::new(0.5, 2.0),
        ]);
        let a = sheet(cone(0.5), &trim);
        assert!(matched(&a, &sheet(cone(0.5), &trim)).is_some());
        assert!(matched(&a, &sheet(cone(0.5 + 1e-6), &trim)).is_none());
    }

    #[test]
    fn a_circular_trim_moves_with_its_centre_only() {
        let circle = |centre: Point2, radius: f64| {
            Curve2::Circle(Circle2 {
                frame: Frame2 {
                    origin: centre,
                    x: Vec2::X,
                    y: Vec2::Y,
                },
                radius,
            })
        };
        let disc = |centre: Point2, radius: f64| {
            sheet(plane(Point3::ZERO, 0.0), &[circle(centre, radius)])
        };
        let a = disc(Point2::new(1.0, 1.0), 0.5);
        assert!(matched(&a, &disc(Point2::new(1.5, 0.75), 0.5)).is_some());
        assert!(matched(&a, &disc(Point2::new(1.5, 0.75), 0.5 + 1e-6)).is_none());
    }
}
