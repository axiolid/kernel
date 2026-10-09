//! Two boundaries shown to meet in space, for an early end to the distance
//! search (#273).
//!
//! Where two boundaries touch or cross, the distance is zero, but the
//! branch and bound can never close on it from its witnesses alone: a
//! patch centre or an edge midpoint lands on the contact only by chance, so
//! it splits the patch pairs round the contact until their spheres are
//! within the accuracy. In a building -- walls on slabs, windows in walls
//! -- that is the common case, and most of the time.
//!
//! So a pair the heap pops with a lower bound of zero is asked first for a
//! point of each element that the two share: the closest points of two
//! line edges, or a point of a line edge and the point of a planar face at
//! the same place -- where the line crosses the plane, or where an end or
//! the middle of the span stands on it. Each point is evaluated on its own
//! boundary exactly as any witness is (on the edge's curve; on the face's
//! surface at parameters inside a certified-inside patch, or that the
//! face's domain certifies), so the two are boundary points whatever they
//! are found to be. They count only when no more than the rounding margin
//! the lower bounds carry apart ([`within_rounding`]): no lower bound can
//! then rise above zero, so the interval `[0, d]` they give is what the
//! search would converge to, reached at once. A gap wider than that margin
//! is never taken for a touch; it is measured by the search as before.
//!
//! Curved faces and circle or ellipse edges show nothing here; such a
//! contact is still found by refinement.

use axiolid_core::{Point2, Point3, Scalar};
use axiolid_curve::Curve3;
use axiolid_evaluate::evaluate3;
use axiolid_surface::Surface;

use super::{point_on, surface_of, Element, Shape, Side};
use crate::exact::ExactMeasureError;

/// A point of `a` and a point of `b`, on their boundaries and within
/// rounding of each other, when the elements' families show one.
pub(super) fn touch(
    side_a: &Side<'_>,
    a: &Element,
    side_b: &Side<'_>,
    b: &Element,
) -> Result<Option<(Point3, Point3)>, ExactMeasureError> {
    Ok(match (a.shape, b.shape) {
        (Shape::Edge { .. }, Shape::Edge { .. }) => edge_edge(side_a, a, side_b, b)?,
        (Shape::Edge { .. }, Shape::Face { .. }) => edge_face(side_a, a, side_b, b)?,
        (Shape::Face { .. }, Shape::Edge { .. }) => {
            edge_face(side_b, b, side_a, a)?.map(|(on_edge, on_face)| (on_face, on_edge))
        }
        (Shape::Face { .. }, Shape::Face { .. }) => None,
    })
}

/// Whether two points are no farther apart than the rounding margin a
/// lower bound between elements about them carries (see `lower_bound`).
pub(super) fn within_rounding(p: Point3, q: Point3) -> bool {
    let d = (p - q).length();
    d <= 1e-12 * (p.length() + q.length())
}

/// The line carrying edge `edge`, when it is one, and its span in order.
fn line_span<'s>(
    side: &'s Side<'_>,
    element: &Element,
) -> Result<Option<(&'s Curve3, Scalar, Scalar)>, ExactMeasureError> {
    let Shape::Edge { edge, t0, t1 } = element.shape else {
        return Ok(None);
    };
    let curve = side.brep.topology().edges()[edge]
        .curve
        .and_then(|id| side.brep.curves3().get(id.index()))
        .ok_or(ExactMeasureError::DanglingReference)?;
    Ok(matches!(curve, Curve3::Line(_)).then_some((curve, t0.min(t1), t0.max(t1))))
}

fn on_curve(curve: &Curve3, t: Scalar) -> Result<Point3, ExactMeasureError> {
    evaluate3(curve, t).map_err(|_| crate::exact::EVALUATION)
}

/// The closest points of two line spans, evaluated on their curves.
fn edge_edge(
    side_a: &Side<'_>,
    a: &Element,
    side_b: &Side<'_>,
    b: &Element,
) -> Result<Option<(Point3, Point3)>, ExactMeasureError> {
    let (Some((ca, a0, a1)), Some((cb, b0, b1))) = (line_span(side_a, a)?, line_span(side_b, b)?)
    else {
        return Ok(None);
    };
    let (p0, p1) = (on_curve(ca, a0)?, on_curve(ca, a1)?);
    let (q0, q1) = (on_curve(cb, b0)?, on_curve(cb, b1)?);
    let (s, t) = segment_parameters(p0, p1, q0, q1);
    if !(s.is_finite() && t.is_finite()) {
        return Ok(None);
    }
    let p = on_curve(ca, a0 + (a1 - a0) * s)?;
    let q = on_curve(cb, b0 + (b1 - b0) * t)?;
    Ok(within_rounding(p, q).then_some((p, q)))
}

/// Parameters in `[0, 1]` of the closest points of segments `p0 p1` and
/// `q0 q1` (Ericson, *Real-Time Collision Detection*, 5.1.9); for parallel
/// segments, one closest pair of the several.
fn segment_parameters(p0: Point3, p1: Point3, q0: Point3, q1: Point3) -> (Scalar, Scalar) {
    let (d1, d2, r) = (p1 - p0, q1 - q0, p0 - q0);
    let (a, e, f) = (d1.dot(d1), d2.dot(d2), d2.dot(r));
    if a <= 0.0 && e <= 0.0 {
        return (0.0, 0.0);
    }
    if a <= 0.0 {
        return (0.0, (f / e).clamp(0.0, 1.0));
    }
    let c = d1.dot(r);
    if e <= 0.0 {
        return ((-c / a).clamp(0.0, 1.0), 0.0);
    }
    let b = d1.dot(d2);
    let denom = a * e - b * b;
    let mut s = if denom > 0.0 {
        ((b * f - c * e) / denom).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let mut t = (b * s + f) / e;
    if t < 0.0 {
        t = 0.0;
        s = (-c / a).clamp(0.0, 1.0);
    } else if t > 1.0 {
        t = 1.0;
        s = ((b - c) / a).clamp(0.0, 1.0);
    }
    (s, t)
}

/// A point of a line span and the point of a planar face patch at the same
/// place: where the line crosses the plane, else where an end or the middle
/// of the span stands on it.
fn edge_face(
    side_e: &Side<'_>,
    e: &Element,
    side_f: &Side<'_>,
    f: &Element,
) -> Result<Option<(Point3, Point3)>, ExactMeasureError> {
    let Some((curve, t0, t1)) = line_span(side_e, e)? else {
        return Ok(None);
    };
    let Curve3::Line(line) = curve else {
        return Ok(None);
    };
    let Shape::Face {
        face,
        lo,
        hi,
        inside,
    } = f.shape
    else {
        return Ok(None);
    };
    let surface = surface_of(side_f.brep, side_f.brep.topology().faces()[face].surface)?;
    let Surface::Plane(plane) = surface else {
        return Ok(None);
    };
    let frame = plane.frame;
    let (gxx, gxy, gyy) = (
        frame.x.dot(frame.x),
        frame.x.dot(frame.y),
        frame.y.dot(frame.y),
    );
    let gram = gxx * gyy - gxy * gxy;
    if gram.is_nan() || gram <= 0.0 {
        return Ok(None);
    }
    // Where the line crosses the plane, first, then the span's ends and
    // middle, for a line lying in the plane or standing on it.
    let normal = frame.x.cross(frame.y);
    let across = normal.dot(line.direction);
    let crossing = (across != 0.0)
        .then(|| normal.dot(frame.origin - line.origin) / across)
        .filter(|t| t.is_finite())
        .map(|t| t.clamp(t0, t1));
    let (u0, u1) = (lo.x.min(hi.x), lo.x.max(hi.x));
    let (v0, v1) = (lo.y.min(hi.y), lo.y.max(hi.y));
    let candidates = [crossing, Some(t0), Some(0.5 * (t0 + t1)), Some(t1)];
    for t in candidates.into_iter().flatten() {
        let p = on_curve(curve, t)?;
        let r = p - frame.origin;
        let (rx, ry) = (frame.x.dot(r), frame.y.dot(r));
        let u = ((gyy * rx - gxy * ry) / gram).clamp(u0, u1);
        let v = ((gxx * ry - gxy * rx) / gram).clamp(v0, v1);
        if !(u.is_finite() && v.is_finite()) {
            continue;
        }
        let at = Point2::new(u, v);
        let q = point_on(surface, at)?;
        if !within_rounding(p, q) {
            continue;
        }
        // The patch is in the face, or the face's domain says the point is.
        let on_face = inside
            || match &side_f.domains[face] {
                Some(domain) => domain.contains(at)? == Some(true),
                None => false,
            };
        if on_face {
            return Ok(Some((p, q)));
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    //! The early end against the search without it (#273): the same lower
    //! bound and a witness pair within rounding where the boxes touch or
    //! cross, in a handful of steps; the very same result, bit for bit,
    //! where they are apart -- by a whole gap or by `1e-9`.

    use axiolid_brep::{ExactBRep, ExactBRepBuilder};
    use axiolid_core::{Frame3, Interval, Point2, Point3, Tolerance, Vec2, Vec3};
    use axiolid_curve::{Curve2, Curve3, Line2, Line3};
    use axiolid_surface::{Plane, Surface};
    use axiolid_topology::{Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Vertex};

    use super::super::{search_with, Found, Metric, MAX_STEPS};
    use super::{segment_parameters, touch, within_rounding};
    use crate::exact_distance::Side;

    /// A box from corner `o` along the edge vectors `x`, `y`, `z`
    /// (orthogonal), with twelve shared line edges and six planar faces.
    fn cuboid(o: Point3, x: Vec3, y: Vec3, z: Vec3) -> ExactBRep {
        let mut builder = ExactBRepBuilder::default();
        let axes = [x, y, z];
        let corner = |i: usize| {
            o + x * (i & 1) as f64 + y * ((i >> 1) & 1) as f64 + z * ((i >> 2) & 1) as f64
        };
        let vertices: Vec<_> = (0..8)
            .map(|i| {
                builder.topology_mut().add_vertex(Vertex {
                    position: corner(i),
                })
            })
            .collect();
        // Edge from corner `i` along axis `k` (bit `k` of `i` clear).
        let mut edges = std::collections::HashMap::new();
        for i in 0..8usize {
            for (k, axis) in axes.iter().enumerate() {
                if i & (1 << k) != 0 {
                    continue;
                }
                let length = axis.length();
                let curve = builder.add_curve3(Curve3::Line(Line3 {
                    origin: corner(i),
                    direction: *axis / length,
                }));
                let edge = builder.topology_mut().add_edge(Edge {
                    start: vertices[i],
                    end: vertices[i | (1 << k)],
                    curve: Some(curve),
                });
                builder.set_edge_interval(edge, Interval::new(0.0, length));
                edges.insert((i, i | (1 << k)), edge);
            }
        }
        // Each face: base corner and its two axes; the third axis fixed.
        for (fixed, side) in [(2, 0), (2, 1), (1, 0), (1, 1), (0, 0), (0, 1)] {
            let (ka, kb) = match fixed {
                0 => (1, 2),
                1 => (0, 2),
                _ => (0, 1),
            };
            let base = side << fixed;
            let (la, lb) = (axes[ka].length(), axes[kb].length());
            let (ux, uy) = (axes[ka] / la, axes[kb] / lb);
            let surface = builder.add_surface(Surface::Plane(Plane {
                frame: Frame3 {
                    origin: corner(base),
                    x: ux,
                    y: uy,
                    z: ux.cross(uy),
                },
            }));
            // Counter-clockwise in (u, v): base, +a, +a+b, +b.
            let ring = [
                base,
                base | (1 << ka),
                base | (1 << ka) | (1 << kb),
                base | (1 << kb),
            ];
            let params = [
                Point2::new(0.0, 0.0),
                Point2::new(la, 0.0),
                Point2::new(la, lb),
                Point2::new(0.0, lb),
            ];
            let mut uses = Vec::new();
            let mut intervals = Vec::new();
            for m in 0..4 {
                let (from, to) = (ring[m], ring[(m + 1) % 4]);
                let (edge, orientation) = match edges.get(&(from, to)) {
                    Some(edge) => (*edge, Orientation::Forward),
                    None => (edges[&(to, from)], Orientation::Reversed),
                };
                let (p, q) = (params[m], params[(m + 1) % 4]);
                let length = (q - p).length();
                let pcurve = builder.add_curve2(Curve2::Line(Line2 {
                    origin: p,
                    direction: Vec2::new((q.x - p.x) / length, (q.y - p.y) / length),
                }));
                uses.push(EdgeUse {
                    edge,
                    orientation,
                    pcurve: Some(pcurve),
                });
                intervals.push(Interval::new(0.0, length));
            }
            let loop_id = builder.topology_mut().add_loop(Loop { edges: uses });
            for (index, interval) in intervals.into_iter().enumerate() {
                builder.set_pcurve_interval(loop_id, index, interval);
            }
            builder.topology_mut().add_face(Face {
                surface: Some(surface),
                bounds: vec![FaceBound {
                    loop_id,
                    orientation: Orientation::Forward,
                    outer: true,
                }],
                orientation: Orientation::Forward,
            });
        }
        builder.finish().expect("a valid box")
    }

    /// A one-face sheet on the plane `z = 0` bounded by the polygon
    /// `corners`, counter-clockwise.
    fn sheet(corners: &[Point2]) -> ExactBRep {
        let mut builder = ExactBRepBuilder::default();
        let surface = builder.add_surface(Surface::Plane(Plane {
            frame: Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
        }));
        let lift = |p: Point2| Point3::new(p.x, p.y, 0.0);
        let vertices: Vec<_> = corners
            .iter()
            .map(|p| {
                builder
                    .topology_mut()
                    .add_vertex(Vertex { position: lift(*p) })
            })
            .collect();
        let n = corners.len();
        let mut uses = Vec::new();
        let mut lengths = Vec::new();
        for i in 0..n {
            let (p, q) = (corners[i], corners[(i + 1) % n]);
            let length = (q - p).length();
            let direction = Vec2::new((q.x - p.x) / length, (q.y - p.y) / length);
            let curve = builder.add_curve3(Curve3::Line(Line3 {
                origin: lift(p),
                direction: Vec3::new(direction.x, direction.y, 0.0),
            }));
            let pcurve = builder.add_curve2(Curve2::Line(Line2 {
                origin: p,
                direction,
            }));
            let edge = builder.topology_mut().add_edge(Edge {
                start: vertices[i],
                end: vertices[(i + 1) % n],
                curve: Some(curve),
            });
            builder.set_edge_interval(edge, Interval::new(0.0, length));
            uses.push(EdgeUse {
                edge,
                orientation: Orientation::Forward,
                pcurve: Some(pcurve),
            });
            lengths.push(length);
        }
        let loop_id = builder.topology_mut().add_loop(Loop { edges: uses });
        for (index, length) in lengths.into_iter().enumerate() {
            builder.set_pcurve_interval(loop_id, index, Interval::new(0.0, length));
        }
        builder.topology_mut().add_face(Face {
            surface: Some(surface),
            bounds: vec![FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation: Orientation::Forward,
        });
        builder.finish().expect("a valid sheet")
    }

    fn unit_box(o: Point3) -> ExactBRep {
        cuboid(o, Vec3::X, Vec3::Y, Vec3::Z)
    }

    /// `v` turned about `axis` by `angle` (Rodrigues).
    fn turn(v: Vec3, axis: Vec3, angle: f64) -> Vec3 {
        let axis = axis.normalize();
        v * angle.cos() + axis.cross(v) * angle.sin() + axis * axis.dot(v) * (1.0 - angle.cos())
    }

    /// A box turned about `axis` by `angle`, its centre at `centre`.
    fn turned_box(centre: Point3, size: Vec3, axis: Vec3, angle: f64) -> ExactBRep {
        let (x, y, z) = (
            turn(Vec3::X * size.x, axis, angle),
            turn(Vec3::Y * size.y, axis, angle),
            turn(Vec3::Z * size.z, axis, angle),
        );
        cuboid(centre - (x + y + z) * 0.5, x, y, z)
    }

    fn tol() -> Tolerance {
        Tolerance::METRE
    }

    fn distance(a: &ExactBRep, b: &ExactBRep, accuracy: f64, touch: bool) -> Found {
        search_with(
            a,
            b,
            tol(),
            Metric::Space,
            MAX_STEPS,
            &mut |lower, upper| upper - lower <= accuracy,
            touch,
        )
        .expect("bounded")
    }

    /// The early end gives the same lower bound as the full search, a
    /// witness pair within rounding, and stops within `steps`.
    fn meets_early(a: &ExactBRep, b: &ExactBRep, steps: usize) {
        let accuracy = 1e-6;
        let early = distance(a, b, accuracy, true);
        let full = distance(a, b, accuracy, false);
        let bounds = &early.bounds;
        assert_eq!(bounds.lower, 0.0, "{bounds:?}");
        assert_eq!(full.bounds.lower, 0.0, "{:?}", full.bounds);
        // At least as tight as the full search, which may even run out of
        // budget before it reaches the accuracy.
        assert!(bounds.upper <= full.bounds.upper, "{:?}", full.bounds);
        assert!(bounds.upper <= accuracy, "{bounds:?}");
        eprintln!(
            "early: {} steps, upper {:e}; full: {} steps, upper {:e}",
            early.steps, bounds.upper, full.steps, full.bounds.upper
        );
        assert!(
            within_rounding(bounds.point_a, bounds.point_b),
            "{bounds:?}"
        );
        assert_eq!(bounds.upper, (bounds.point_a - bounds.point_b).length());
        assert!(
            early.steps <= steps,
            "{} steps, the full search {}",
            early.steps,
            full.steps
        );
        assert!(
            early.steps <= full.steps,
            "{} vs {}",
            early.steps,
            full.steps
        );
        // Deterministic: the same pair, the same points, every time.
        let again = distance(a, b, accuracy, true);
        assert_eq!(again.bounds, early.bounds);
        assert_eq!(again.steps, early.steps);
    }

    /// Apart, the early end never fires: the result is the full search's.
    fn unchanged(a: &ExactBRep, b: &ExactBRep, accuracy: f64) -> Found {
        let early = distance(a, b, accuracy, true);
        let full = distance(a, b, accuracy, false);
        assert_eq!(early.bounds, full.bounds);
        assert_eq!(early.steps, full.steps);
        early
    }

    #[test]
    fn two_boxes_sharing_a_face_meet_in_a_few_steps() {
        meets_early(
            &unit_box(Point3::ZERO),
            &unit_box(Point3::new(0.0, 0.0, 1.0)),
            8,
        );
    }

    #[test]
    fn two_turned_boxes_sharing_a_face_meet_in_a_few_steps() {
        // Off the axes, so the shared face is met only within rounding.
        let axis = Vec3::new(1.0, 2.0, 3.0);
        let size = Vec3::new(1.0, 0.7, 0.4);
        let a = turned_box(Point3::new(2.0, -1.0, 0.5), size, axis, 0.7);
        // Stacked, and slid along the shared face so the centres differ.
        let lift = turn(Vec3::new(0.3, 0.2, size.z), axis, 0.7);
        let b = turned_box(Point3::new(2.0, -1.0, 0.5) + lift, size, axis, 0.7);
        meets_early(&a, &b, 8);
    }

    #[test]
    fn a_box_standing_on_part_of_another_meets_in_a_few_steps() {
        meets_early(
            &unit_box(Point3::ZERO),
            &cuboid(
                Point3::new(0.3, 0.4, 1.0),
                Vec3::X * 0.5,
                Vec3::Y * 0.25,
                Vec3::Z * 2.0,
            ),
            8,
        );
    }

    #[test]
    fn two_crossing_boxes_meet_in_a_few_steps() {
        meets_early(
            &turned_box(
                Point3::new(0.5, 0.5, 0.5),
                Vec3::new(1.0, 1.0, 1.0),
                Vec3::Z,
                0.0,
            ),
            &turned_box(
                Point3::new(0.9, 0.7, 0.6),
                Vec3::new(1.5, 0.4, 0.3),
                Vec3::new(1.0, 2.0, 3.0),
                0.7,
            ),
            8,
        );
    }

    #[test]
    fn two_boxes_sharing_an_edge_meet_in_a_few_steps() {
        meets_early(
            &unit_box(Point3::ZERO),
            &unit_box(Point3::new(1.0, 0.0, 1.0)),
            8,
        );
    }

    #[test]
    fn two_edges_crossing_at_a_point_meet_in_a_few_steps() {
        // Two triangles whose long edges cross at one point, off both
        // midpoints, each lying on the far side of the other: the
        // boundaries meet only where the edges cross, which is on neither
        // face's certified interior, so only the edge pair shows it.
        let triangle = || {
            sheet(&[
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 0.0),
                Point2::new(0.0, 1.0),
            ])
        };
        let c = Point3::new(0.37, 0.63, 0.0);
        let out = Vec3::new(1.0, 1.0, 0.0).normalize();
        let along = (out + Vec3::Z).normalize();
        let inward = (out - Vec3::Z).normalize();
        // The triangle's long edge onto `along`, its inward normal onto
        // `inward`, its own normal onto their cross product.
        let local = axiolid_core::Mat3::from_cols(
            Vec3::new(-1.0, 1.0, 0.0).normalize(),
            Vec3::new(-1.0, -1.0, 0.0).normalize(),
            Vec3::Z,
        );
        let world = axiolid_core::Mat3::from_cols(along, inward, along.cross(inward));
        let turn = world * local.transpose();
        let place = axiolid_core::Transform3::from_mat3_translation(
            turn,
            c - turn * Vec3::new(0.3, 0.7, 0.0),
        );
        let tilted = triangle().transformed(&place).expect("a rigid placement");
        meets_early(&triangle(), &tilted, 8);
    }

    #[test]
    fn a_gap_of_a_nanometre_is_measured_as_a_gap() {
        let a = unit_box(Point3::ZERO);
        let b = unit_box(Point3::new(0.2, 0.1, 1.0 + 1e-9));
        let found = unchanged(&a, &b, 1e-12);
        assert!(found.bounds.lower > 0.0, "{:?}", found.bounds);
        assert!(found.bounds.lower <= 1e-9 && 1e-9 <= found.bounds.upper + 1e-15);
    }

    #[test]
    fn separated_boxes_measure_as_before() {
        let a = unit_box(Point3::ZERO);
        unchanged(&a, &unit_box(Point3::new(0.0, 0.0, 1.25)), 1e-9);
        unchanged(&a, &unit_box(Point3::new(1.5, 2.0, 0.5)), 1e-9);
        unchanged(
            &a,
            &turned_box(
                Point3::new(3.0, 0.5, 0.5),
                Vec3::new(1.0, 0.5, 2.0),
                Vec3::new(1.0, 2.0, 3.0),
                0.7,
            ),
            1e-6,
        );
    }

    #[test]
    fn a_touch_is_shown_only_where_the_elements_meet() {
        let a = unit_box(Point3::ZERO);
        let linear = tol().linear();
        let side_a = Side::new(&a, linear, Metric::Space).expect("a side");
        let shown = |b: &ExactBRep| {
            let side_b = Side::new(b, linear, Metric::Space).expect("a side");
            let mut count = 0;
            for ea in &side_a.elements {
                for eb in &side_b.elements {
                    if let Some((p, q)) = touch(&side_a, ea, &side_b, eb).expect("evaluates") {
                        assert!(within_rounding(p, q));
                        count += 1;
                    }
                }
            }
            count
        };
        // A slim bar through the top face, and the same bar lifted clear.
        let bar = |z: f64| {
            cuboid(
                Point3::new(0.4, 0.4, z),
                Vec3::X * 0.2,
                Vec3::Y * 0.2,
                Vec3::Z,
            )
        };
        assert!(shown(&bar(0.5)) > 0);
        assert_eq!(shown(&bar(1.0 + 1e-9)), 0);
        assert_eq!(shown(&bar(1.5)), 0);
        // Standing on it: the bar's bottom edges lie in the top face.
        assert!(shown(&bar(1.0)) > 0);
    }

    #[test]
    fn a_line_through_the_plane_outside_the_trim_is_not_a_touch() {
        // A triangle in z = 0, and a bar crossing that plane over the part
        // of the triangle's parameter box the triangle does not cover: its
        // edges pierce the plane, but not the face.
        let triangle = sheet(&[
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        ]);
        let bar = cuboid(
            Point3::new(0.7, 0.7, -0.5),
            Vec3::X * 0.2,
            Vec3::Y * 0.2,
            Vec3::Z,
        );
        let found = unchanged(&triangle, &bar, 1e-9);
        // The bar's corner (0.7, 0.7) is 0.4 / sqrt 2 from the hypotenuse.
        let gap = 0.4 / 2.0_f64.sqrt();
        assert!(found.bounds.lower > 0.0, "{:?}", found.bounds);
        assert!(found.bounds.lower <= gap && gap <= found.bounds.upper + 1e-12);
    }

    #[test]
    fn segment_parameters_find_crossing_and_overlapping_segments() {
        let (s, t) = segment_parameters(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.5, -1.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
        );
        assert_eq!((s, t), (0.25, 0.5));
        // Collinear and overlapping: a shared point.
        let (p0, p1) = (Point3::ZERO, Point3::new(2.0, 0.0, 0.0));
        let (q0, q1) = (Point3::new(3.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
        let (s, t) = segment_parameters(p0, p1, q0, q1);
        let p = p0 + (p1 - p0) * s;
        let q = q0 + (q1 - q0) * t;
        assert_eq!((p - q).length(), 0.0);
        // Skew and apart.
        let (s, t) = segment_parameters(
            Point3::ZERO,
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, -1.0, 1.0),
            Point3::new(0.5, 1.0, 1.0),
        );
        assert_eq!((s, t), (0.5, 0.5));
    }
}
