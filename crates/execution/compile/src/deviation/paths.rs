//! The deviation of each graph path the reference compiler meshes (#232).
//!
//! See the parent module for what each path's bound rests on.

use axiolid_construct::profile::{profile_deviation, ProfileDeviation};
use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_curve::Curve3;
use axiolid_mesh::TriMesh;
use axiolid_model::{GeometryGraph, GeometryNode, NodeId, SolidOperation};
use axiolid_primitive::Primitive;

use super::{Deviation, DeviationBound, DeviationPath};
use crate::certify::{certify, Cell, Patch};
use crate::directrix::DirectrixKind;

/// The deviation of a non-boolean solid the compiler built as `mesh`.
pub(crate) fn of_solid(
    graph: &GeometryGraph,
    operation: &SolidOperation,
    options: &ExecutionOptions,
    mesh: &TriMesh,
) -> GeomResult<Deviation> {
    let chord = crate::compiler::chord_error(options);
    let half = 0.5 * chord;
    let tolerance = options.tolerance();
    let profile = |id: NodeId, budget: Scalar| -> GeomResult<ProfileDeviation> {
        match graph.get(id) {
            Some(GeometryNode::Profile(shape)) => profile_deviation(shape, budget, tolerance),
            _ => Err(GeomError::InvalidInput(format!(
                "profile {id:?} is not a Profile node"
            ))),
        }
    };
    let proven = |deviation: ProfileDeviation, extra: Scalar| match deviation {
        ProfileDeviation::Bounded(d) => DeviationBound::Proven(d + extra),
        ProfileDeviation::Unbounded(reason) => DeviationBound::Unbounded(reason),
        _ => DeviationBound::Unbounded("profile deviation"),
    };
    let (path, detail, bound) = match operation {
        // Walls are the profile's chords moved along a straight line, caps
        // are its rings: both within the profile's bound.
        SolidOperation::Extrusion { profile: id, .. } => (
            DeviationPath::Extrusion,
            "",
            proven(profile(*id, chord)?, 0.0),
        ),
        // The profile gets half the budget, the turn the other half (#231).
        SolidOperation::Revolution { profile: id, .. } => (
            DeviationPath::Revolution,
            "",
            proven(profile(*id, half)?, half),
        ),
        SolidOperation::TaperedExtrusion { .. } => (
            DeviationPath::TaperedExtrusion,
            "",
            DeviationBound::Unbounded("tapered extrusion: its twisted walls are not bounded"),
        ),
        SolidOperation::TaperedRevolution {
            start_profile,
            end_profile,
            ..
        } => {
            let start = profile(*start_profile, half)?;
            let end = profile(*end_profile, half)?;
            let bound = match (start, end) {
                // Straight-edged rings are the exact profiles, vertex for
                // vertex, so the blended rings are the exact blend and
                // #231's turn bound is the whole story.
                (ProfileDeviation::Bounded(a), ProfileDeviation::Bounded(b))
                    if a == 0.0 && b == 0.0 =>
                {
                    DeviationBound::Proven(half)
                }
                _ => DeviationBound::Unbounded(
                    "tapered revolution between curved profiles: ring correspondence",
                ),
            };
            (DeviationPath::TaperedRevolution, "", bound)
        }
        SolidOperation::SweptDisk {
            directrix,
            radius,
            inner_radius,
            parameter_range,
            ..
        } => {
            // Lines, polylines (filleted or not) and chains of segments and
            // arcs are swept piece by piece by `axiolid_construct::pipe`,
            // which proves the budget (#232), mitred corners included (#245);
            // the compiler takes that path whenever the directrix resolves
            // to pieces, and refuses what it cannot prove.
            if crate::directrix::pieces(graph, *directrix, *parameter_range, options)?.is_some() {
                return Ok(Deviation::one(
                    DeviationPath::SweptDisk,
                    "segments and arcs",
                    DeviationBound::Proven(chord),
                ));
            }
            let kind = crate::directrix::kind(graph, *directrix, *parameter_range, options)?;
            match kind {
                DirectrixKind::Segment => (
                    DeviationPath::SweptDisk,
                    "segment",
                    DeviationBound::Proven(chord),
                ),
                DirectrixKind::Arc(_) => (
                    DeviationPath::SweptDisk,
                    "circular arc",
                    DeviationBound::Proven(chord),
                ),
                DirectrixKind::Smooth(curve, span) => {
                    let detail = match curve {
                        Curve3::Ellipse(_) => "ellipse",
                        _ => "B-spline",
                    };
                    let bound = tube_bound(
                        &curve,
                        (span.start, span.end),
                        *radius,
                        *inner_radius,
                        mesh,
                        chord,
                    )
                    .map_or(
                        DeviationBound::Unbounded("swept disk: the directrix could not be bounded"),
                        DeviationBound::Certified,
                    );
                    (DeviationPath::SweptDisk, detail, bound)
                }
                DirectrixKind::Other(name) => (
                    DeviationPath::SweptDisk,
                    name,
                    DeviationBound::Unbounded("swept disk along this directrix family"),
                ),
            }
        }
        SolidOperation::FixedReferenceSweep {
            profile: id,
            directrix,
            reference_direction,
            parameter_range,
        } => {
            let deviation = profile(*id, half)?;
            let kind = crate::directrix::kind(graph, *directrix, *parameter_range, options)?;
            match kind {
                // A fixed frame translated along a line: the walls are the
                // profile's chords, moved.
                DirectrixKind::Segment => (
                    DeviationPath::FixedReferenceSweep,
                    "segment",
                    proven(deviation, 0.0),
                ),
                // #231: with the reference normal to the arc's plane the
                // stations are rotations of one ring about the arc's axis.
                DirectrixKind::Arc(circle)
                    if circle
                        .frame
                        .x
                        .cross(circle.frame.y)
                        .normalize_or_zero()
                        .cross(reference_direction.normalize_or_zero())
                        .length()
                        <= 1e-9 =>
                {
                    (
                        DeviationPath::FixedReferenceSweep,
                        "circular arc",
                        proven(deviation, half),
                    )
                }
                _ => (
                    DeviationPath::FixedReferenceSweep,
                    "",
                    DeviationBound::Unbounded(
                        "fixed-reference sweep: frame law along this directrix is not bounded",
                    ),
                ),
            }
        }
        SolidOperation::SurfaceCurveSweep { .. } => (
            DeviationPath::SurfaceCurveSweep,
            "",
            DeviationBound::Unbounded("surface-curve sweep: frame law is not bounded"),
        ),
        SolidOperation::SectionedSpine { .. } => (
            DeviationPath::SectionedSpine,
            "",
            DeviationBound::Unbounded(
                "sectioned spine: the surface between sections is not defined exactly",
            ),
        ),
        SolidOperation::StationedSpine { .. } => (
            DeviationPath::StationedSpine,
            "",
            DeviationBound::Unbounded(STATIONED_SAMPLING),
        ),
        SolidOperation::BoundedHalfSpace { .. } => (
            DeviationPath::BoundedHalfSpace,
            "",
            DeviationBound::Unbounded(
                "bounded half-space: a finite stand-in for an unbounded solid",
            ),
        ),
        _ => (
            DeviationPath::Unreported,
            "",
            DeviationBound::Unbounded("solid family"),
        ),
    };
    Ok(Deviation::one(path, detail, bound))
}

/// Why station-placed sweeps are unbounded (#241): the sections between
/// stations are bisected until each placed point's midpoint test passes,
/// which samples the surface rather than bounding it.
const STATIONED_SAMPLING: &str =
    "stationed sections: the surface between stations is sampled, not bounded";

/// The deviation of a sectioned surface (#241); see [`STATIONED_SAMPLING`].
pub(crate) fn of_sectioned_surface() -> Deviation {
    Deviation::one(
        DeviationPath::SectionedSurface,
        "",
        DeviationBound::Unbounded(STATIONED_SAMPLING),
    )
}

/// The deviation of a CSG primitive meshed for chord budget `chord`.
///
/// Every curved primitive is proven within the budget by its tessellator
/// (`axiolid_reference::primitive`): spheres and tori by #231, cylinders and
/// cones by #232 (one way curved, the ring's sagitta is the whole distance;
/// a budget past the segment cap is refused, not clamped).
pub(crate) fn of_primitive(primitive: &Primitive, chord: Scalar) -> Deviation {
    let (detail, bound) = match *primitive {
        Primitive::Block { .. } => ("block", DeviationBound::Proven(0.0)),
        Primitive::Pyramid { .. } => ("pyramid", DeviationBound::Proven(0.0)),
        Primitive::Wedge { .. } => ("wedge", DeviationBound::Proven(0.0)),
        Primitive::Sphere { .. } => ("sphere", DeviationBound::Proven(chord)),
        Primitive::Torus { .. } => ("torus", DeviationBound::Proven(chord)),
        Primitive::Cylinder { .. } => ("cylinder", DeviationBound::Proven(chord)),
        Primitive::Cone { .. } => ("cone", DeviationBound::Proven(chord)),
        _ => ("primitive", DeviationBound::Unbounded("primitive family")),
    };
    Deviation::one(DeviationPath::Primitive, detail, bound)
}

/// The deviation of a curve-bounded plane: its boundaries are flattened by
/// the curve flattener, so the chord budget where every boundary family is
/// certified, plus the largest merge of near-duplicate points (`merged`),
/// scaled by the plane frame's stretch (the boundaries live in its
/// parameters).
pub(crate) fn of_curve_bounded(
    graph: &GeometryGraph,
    basis: NodeId,
    boundaries: &[NodeId],
    chord: Scalar,
    merged: Scalar,
) -> Deviation {
    let stretch = match graph.get(basis) {
        Some(GeometryNode::Surface(axiolid_surface::Surface::Plane(plane))) => {
            axiolid_reference::bound::frame_stretch3(&plane.frame)
        }
        _ => Scalar::INFINITY,
    };
    let mut worst = DeviationBound::Proven(0.0);
    for &id in boundaries {
        let bound = match graph.get(id) {
            Some(GeometryNode::Curve2(curve)) => {
                let straight = matches!(
                    curve,
                    axiolid_curve::Curve2::Line(_) | axiolid_curve::Curve2::Polyline(_)
                );
                if straight {
                    DeviationBound::Proven(0.0)
                } else if axiolid_reference::bound::certifies_flattening2(curve) {
                    DeviationBound::Proven(chord)
                } else {
                    DeviationBound::Unbounded("curve-bounded plane boundary family")
                }
            }
            Some(GeometryNode::Curve3(curve)) => {
                if matches!(curve, Curve3::Line(_) | Curve3::Polyline(_)) {
                    DeviationBound::Proven(0.0)
                } else if axiolid_reference::bound::certifies_flattening3(curve) {
                    DeviationBound::Proven(chord)
                } else {
                    DeviationBound::Unbounded("curve-bounded plane boundary family")
                }
            }
            _ => DeviationBound::Unbounded("curve-bounded plane boundary"),
        };
        worst = worst.worst(bound);
    }
    let bound = match worst {
        DeviationBound::Proven(d) if stretch.is_finite() => {
            DeviationBound::Proven((d + merged) * stretch)
        }
        DeviationBound::Proven(_) => DeviationBound::Unbounded("curve-bounded plane frame"),
        other => other,
    };
    Deviation::one(DeviationPath::CurveBoundedPlane, "", bound)
}

// --- certified tubes and round primitives ------------------------------------

/// A disk swept along a smooth curve: its wall (and inner wall), and the
/// two end caps square to the curve.
fn tube_bound(
    curve: &Curve3,
    span: (Scalar, Scalar),
    radius: Scalar,
    inner: Option<Scalar>,
    mesh: &TriMesh,
    target: Scalar,
) -> Option<Scalar> {
    let (lo, hi) = (span.0.min(span.1), span.0.max(span.1));
    if !(lo < hi && radius > 0.0) {
        return None;
    }
    let pieces = smooth_pieces(curve, lo, hi)?;
    let mut walls: Vec<TubeWall<'_>> = Vec::new();
    for &(a, b) in &pieces {
        let reference = perpendicular(unit_tangent(curve, 0.5 * (a + b))?);
        walls.push(TubeWall {
            curve,
            reference,
            radius,
        });
        if let Some(r) = inner {
            walls.push(TubeWall {
                curve,
                reference,
                radius: r,
            });
        }
    }
    let r0 = inner.unwrap_or(0.0);
    let caps: Vec<Disk> = [lo, hi]
        .iter()
        .map(|&s| {
            let point = axiolid_reference::curve::evaluate3(curve, s).ok()?;
            let t = unit_tangent(curve, s)?;
            let n = perpendicular(t);
            Some(Disk {
                centre: point,
                n,
                b: t.cross(n),
                radii: (r0, radius),
            })
        })
        .collect::<Option<_>>()?;
    let mut patches: Vec<&dyn Patch> = Vec::new();
    let mut cells = Vec::new();
    let per_piece = if inner.is_some() { 2 } else { 1 };
    for (k, wall) in walls.iter().enumerate() {
        let (a, b) = pieces[k / per_piece];
        patches.push(wall);
        for j in 0..8 {
            let t0 = core::f64::consts::TAU * j as Scalar / 8.0;
            let t1 = core::f64::consts::TAU * (j + 1) as Scalar / 8.0;
            cells.push(Cell {
                patch: patches.len() - 1,
                x: (a, b),
                y: (t0, t1),
            });
        }
    }
    for cap in &caps {
        patches.push(cap);
        for j in 0..8 {
            let t0 = core::f64::consts::TAU * j as Scalar / 8.0;
            let t1 = core::f64::consts::TAU * (j + 1) as Scalar / 8.0;
            cells.push(Cell {
                patch: patches.len() - 1,
                x: cap.radii,
                y: (t0, t1),
            });
        }
    }
    certify(&patches, cells, mesh, target)
}

/// `[lo, hi]` cut where the curve is not `C^2`, then halved until the
/// tangent turns by at most one radian over each piece, so a reference
/// normal to the tangent at a piece's middle stays well away from it.
fn smooth_pieces(curve: &Curve3, lo: Scalar, hi: Scalar) -> Option<Vec<(Scalar, Scalar)>> {
    let mut cuts: Vec<Scalar> = axiolid_reference::bound::continuity_breaks3(curve, 2)
        .into_iter()
        .filter(|&t| t > lo && t < hi)
        .collect();
    cuts.sort_by(Scalar::total_cmp);
    cuts.dedup();
    let mut stack = Vec::new();
    let mut from = lo;
    for cut in cuts.into_iter().chain(core::iter::once(hi)) {
        stack.push((from, cut, 0_u32));
        from = cut;
    }
    stack.reverse();
    let mut out = Vec::new();
    while let Some((a, b, depth)) = stack.pop() {
        let bounds = axiolid_reference::bound::curve_derivative_bounds3(curve, a, b)?;
        let speed = axiolid_reference::curve::derivative3(curve, 0.5 * (a + b))
            .ok()?
            .length();
        let slowest = speed - 0.5 * (b - a) * bounds.second;
        let turn = if slowest > 0.0 {
            (b - a) * bounds.second / slowest
        } else {
            Scalar::INFINITY
        };
        if turn <= 1.0 {
            out.push((a, b));
        } else if depth >= 24 {
            return None;
        } else {
            let m = 0.5 * (a + b);
            stack.push((m, b, depth + 1));
            stack.push((a, m, depth + 1));
        }
    }
    Some(out)
}

fn unit_tangent(curve: &Curve3, s: Scalar) -> Option<Vec3> {
    let d = axiolid_reference::curve::derivative3(curve, s).ok()?;
    let t = d.normalize_or_zero();
    (t != Vec3::ZERO && t.is_finite()).then_some(t)
}

/// A unit vector perpendicular to unit `t`, from the axis least along it.
fn perpendicular(t: Vec3) -> Vec3 {
    let axis = if t.x.abs() <= t.y.abs() && t.x.abs() <= t.z.abs() {
        Vec3::X
    } else if t.y.abs() <= t.z.abs() {
        Vec3::Y
    } else {
        Vec3::Z
    };
    (axis - t * t.dot(axis)).normalize()
}

/// The wall of a disk of `radius` swept along a smooth curve:
/// `T(s, theta) = c(s) + radius (cos theta n(s) + sin theta b(s))`, with
/// `n` the reference projected normal to the tangent `t` and normalised,
/// `b = t x n`. It covers every point at `radius` from `c(s)` in the
/// plane normal to the curve.
///
/// Second derivatives, with `D_k` bounding `|c^(k)|` and `sigma` the
/// smallest speed over the cell: `|t'| <= k1 = D_2 / sigma`, `|t''| <= k2 =
/// 2 D_3 / sigma + 6 D_2^2 / sigma^2`; for `P = a - (a.t) t`, `|P'| <= 2 k1`
/// and `|P''| <= 2 k2 + 2 k1^2`; for `n = P / |P|` with `|P| >= p`,
/// `|n'| <= m1 = |P'| / p` and `|n''| <= m2 = |P''| / p + 3 |P'|^2 / p^2`;
/// `|b'| <= k1 + m1`, `|b''| <= k2 + 2 k1 m1 + m2`. Then `|T_ss| <= D_2 +
/// r sqrt(m2^2 + |b''|^2)`, `|T_s theta| <= r sqrt(m1^2 + |b'|^2)` and
/// `|T_theta theta| = r`.
struct TubeWall<'c> {
    curve: &'c Curve3,
    reference: Vec3,
    radius: Scalar,
}

impl Patch for TubeWall<'_> {
    fn jet(&self, s: Scalar, theta: Scalar) -> Option<(Point3, Vec3, Vec3)> {
        let c = axiolid_reference::curve::evaluate3(self.curve, s).ok()?;
        let c1 = axiolid_reference::curve::derivative3(self.curve, s).ok()?;
        let c2 = axiolid_reference::curve::second_derivative3(self.curve, s).ok()?;
        let sigma = c1.length();
        if sigma <= 0.0 || sigma.is_nan() {
            return None;
        }
        let t = c1 / sigma;
        let t1 = (c2 - t * c2.dot(t)) / sigma;
        let a = self.reference;
        let at = a.dot(t);
        let p_vec = a - t * at;
        let p = p_vec.length();
        if p <= 0.0 || p.is_nan() {
            return None;
        }
        let n = p_vec / p;
        let p1 = -(t * a.dot(t1)) - t1 * at;
        let n1 = (p1 - n * n.dot(p1)) / p;
        let b = t.cross(n);
        let b1 = t1.cross(n) + t.cross(n1);
        let (sin, cos) = theta.sin_cos();
        let r = self.radius;
        let point = c + (n * cos + b * sin) * r;
        let ds = c1 + (n1 * cos + b1 * sin) * r;
        let dtheta = (b * cos - n * sin) * r;
        (point.is_finite() && ds.is_finite()).then_some((point, ds, dtheta))
    }

    fn second(&self, s: (Scalar, Scalar), _theta: (Scalar, Scalar)) -> Option<[Scalar; 3]> {
        let bounds = axiolid_reference::bound::curve_derivative_bounds3(self.curve, s.0, s.1)?;
        let (d2, d3) = (bounds.second, bounds.third);
        let mid = 0.5 * (s.0 + s.1);
        let half = 0.5 * (s.1 - s.0).abs();
        let c1 = axiolid_reference::curve::derivative3(self.curve, mid).ok()?;
        let sigma = c1.length() - half * d2;
        if sigma <= 0.0 || sigma.is_nan() {
            return None;
        }
        let k1 = d2 / sigma;
        let k2 = 2.0 * d3 / sigma + 6.0 * d2 * d2 / (sigma * sigma);
        let along = (self.reference.dot(c1 / c1.length()).abs() + k1 * half).min(1.0);
        let p = (1.0 - along * along).max(0.0).sqrt();
        if p <= 0.0 || p.is_nan() {
            return None;
        }
        let (dp1, dp2) = (2.0 * k1, 2.0 * k2 + 2.0 * k1 * k1);
        let m1 = dp1 / p;
        let m2 = dp2 / p + 3.0 * dp1 * dp1 / (p * p);
        let b1 = k1 + m1;
        let b2 = k2 + 2.0 * k1 * m1 + m2;
        let r = self.radius;
        let out = [
            d2 + r * (m2 * m2 + b2 * b2).sqrt(),
            r * (m1 * m1 + b1 * b1).sqrt(),
            r,
        ];
        out.iter().all(|v| v.is_finite()).then_some(out)
    }
}

/// A flat disk or annulus: `T(rho, theta) = centre + rho (cos theta n +
/// sin theta b)` for orthonormal `n, b`. `|T_rho rho| = 0`, `|T_rho theta|
/// = 1`, `|T_theta theta| = rho`.
struct Disk {
    centre: Point3,
    n: Vec3,
    b: Vec3,
    radii: (Scalar, Scalar),
}

impl Patch for Disk {
    fn jet(&self, rho: Scalar, theta: Scalar) -> Option<(Point3, Vec3, Vec3)> {
        let (sin, cos) = theta.sin_cos();
        let w = self.n * cos + self.b * sin;
        Some((
            self.centre + w * rho,
            w,
            (self.b * cos - self.n * sin) * rho,
        ))
    }

    fn second(&self, rho: (Scalar, Scalar), _theta: (Scalar, Scalar)) -> Option<[Scalar; 3]> {
        Some([0.0, 1.0, rho.0.abs().max(rho.1.abs())])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_curve::{BSplineCurve3, KnotSpec};

    /// Second partials by central differences of the exact first partials,
    /// over a lattice of the cell, never above the patch's bounds.
    fn assert_second_bounded(patch: &dyn Patch, x: (Scalar, Scalar), y: (Scalar, Scalar)) {
        let [a, b, c] = patch.second(x, y).unwrap();
        let h = 1e-6;
        for i in 0..=10 {
            for j in 0..=10 {
                let s = x.0 + (x.1 - x.0) * (0.05 + 0.9 * i as Scalar / 10.0);
                let t = y.0 + (y.1 - y.0) * (0.05 + 0.9 * j as Scalar / 10.0);
                let (_, xp, yp) = patch.jet(s + h, t).unwrap();
                let (_, xm, ym) = patch.jet(s - h, t).unwrap();
                let (_, _, yq) = patch.jet(s, t + h).unwrap();
                let (_, _, yr) = patch.jet(s, t - h).unwrap();
                let xx = ((xp - xm) / (2.0 * h)).length();
                let xy = ((yp - ym) / (2.0 * h)).length();
                let yy = ((yq - yr) / (2.0 * h)).length();
                let slack = 1e-4;
                assert!(
                    xx <= a * (1.0 + slack) + slack,
                    "xx {xx} > {a} at ({s}, {t})"
                );
                assert!(
                    xy <= b * (1.0 + slack) + slack,
                    "xy {xy} > {b} at ({s}, {t})"
                );
                assert!(
                    yy <= c * (1.0 + slack) + slack,
                    "yy {yy} > {c} at ({s}, {t})"
                );
            }
        }
    }

    #[test]
    fn a_tube_wall_bounds_its_second_partials() {
        let curve = Curve3::BSpline(BSplineCurve3 {
            degree: 3,
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.3, 0.0, 0.2),
                Point3::new(0.4, 0.4, 0.0),
                Point3::new(0.8, 0.3, 0.3),
                Point3::new(1.0, 0.0, 0.1),
            ],
            knots: vec![0.0, 0.5, 1.0],
            multiplicities: vec![4, 1, 4],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        });
        for (a, b) in smooth_pieces(&curve, 0.0, 1.0).unwrap() {
            let wall = TubeWall {
                curve: &curve,
                reference: perpendicular(unit_tangent(&curve, 0.5 * (a + b)).unwrap()),
                radius: 1.0,
            };
            assert_second_bounded(&wall, (a, b), (0.0, core::f64::consts::TAU));
        }
    }

    #[test]
    fn a_disk_bounds_its_second_partials() {
        let disk = Disk {
            centre: Point3::ZERO,
            n: Vec3::X,
            b: Vec3::Y,
            radii: (0.2, 0.7),
        };
        assert_second_bounded(&disk, (0.2, 0.7), (0.0, core::f64::consts::TAU));
    }
}
