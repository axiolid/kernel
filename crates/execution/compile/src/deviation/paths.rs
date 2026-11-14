//! The deviation of each graph path the reference compiler meshes (#232).
//!
//! See the parent module for what each path's bound rests on.

use axiolid_construct::pipe::{fixed_reference_plan, ring_reach, FixedReferencePlan};
use axiolid_construct::profile::{profile_deviation, profile_rings, ProfileDeviation};
use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_curve::{Curve2, Curve3, Elevated3, ElevationLaw};
use axiolid_mesh::TriMesh;
use axiolid_model::{GeometryGraph, GeometryNode, NodeId, SolidOperation};
use axiolid_primitive::Primitive;
use axiolid_reference::bound::{ElevatedPositionTable, SeamSide};

use super::{Deviation, DeviationBound, DeviationPath};
use crate::certify::{certify_thin_quads, Cell, Patch};
use crate::directrix::DirectrixKind;

/// The deviation of a non-boolean solid the compiler built as `mesh`,
/// against the budget of `options`; `built` holds the options it was built
/// with (a finer chord when the compiler rebuilt it, #252), which a pipe's
/// plan is read at.
pub(crate) fn of_solid(
    graph: &GeometryGraph,
    operation: &SolidOperation,
    options: &ExecutionOptions,
    built: &ExecutionOptions,
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
            fillet_radius,
        } => {
            // Lines, polylines (filleted or not, open or closed) and chains
            // of segments and arcs are swept piece by piece by
            // `axiolid_construct::pipe`, which proves the budget (#232),
            // mitred corners between segments and loops included (#245,
            // #248); a corner mitred beside an arc is certified against the
            // mesh (#248). The compiler takes that path whenever the
            // directrix resolves to pieces, and refuses what it cannot
            // bound.
            if let Some(path) =
                crate::directrix::pieces(graph, *directrix, *parameter_range, options)?
            {
                return super::mitre::of_pipe(
                    &path,
                    *radius,
                    *inner_radius,
                    *fillet_radius,
                    crate::compiler::chord_error(built),
                    chord,
                    mesh,
                );
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
                        Curve3::Elevated(_) => "elevated curve",
                        Curve3::Banked(_) => "banked curve",
                        _ => "B-spline",
                    };
                    let bound = if let Some(reason) = uncertified_points(&curve) {
                        DeviationBound::Unbounded(reason)
                    } else {
                        tube_bound(
                            &curve,
                            (span.start, span.end),
                            *radius,
                            *inner_radius,
                            mesh,
                            chord,
                        )
                        .map_or(
                            DeviationBound::Unbounded(
                                "swept disk: the directrix could not be bounded",
                            ),
                            DeviationBound::Certified,
                        )
                    };
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
            // A tangent-continuous chain of segments and arcs is swept
            // piece by piece (#122): proven when no arc twists the frame
            // past what the walls pay for.
            if let Some(path) = crate::directrix::fixed_reference_pieces(
                graph,
                *directrix,
                *parameter_range,
                options,
            )? {
                let reach = match graph.get(*id) {
                    Some(GeometryNode::Profile(shape)) => {
                        ring_reach(&profile_rings(shape, half, tolerance)?)
                    }
                    _ => 0.0,
                };
                let plan = fixed_reference_plan(
                    &path.pieces,
                    path.closed,
                    *reference_direction,
                    reach,
                    chord,
                    tolerance,
                )?;
                match plan {
                    FixedReferencePlan::Untwisted(chain) => (
                        DeviationPath::FixedReferenceSweep,
                        "segments and arcs",
                        proven(deviation, chain.walls()),
                    ),
                    _ => (
                        DeviationPath::FixedReferenceSweep,
                        "",
                        DeviationBound::Unbounded(
                            "fixed-reference sweep: the frame twists along an arc oblique to the reference",
                        ),
                    ),
                }
            } else {
                // Sampled directrices: a curve that is a segment in all but
                // spelling (a straight elevated curve) keeps its fixed frame;
                // the rest are not bounded.
                let kind = crate::directrix::kind(graph, *directrix, *parameter_range, options)?;
                match kind {
                    DirectrixKind::Segment => (
                        DeviationPath::FixedReferenceSweep,
                        "segment",
                        proven(deviation, 0.0),
                    ),
                    _ => (
                        DeviationPath::FixedReferenceSweep,
                        "",
                        DeviationBound::Unbounded(
                            "fixed-reference sweep: frame law along this directrix is not bounded",
                        ),
                    ),
                }
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
        SolidOperation::StationedSpine { .. } | SolidOperation::SectionsAtStations { .. } => {
            return match stationed_spine(graph, operation, options) {
                Some(run) => run.map(|(_, deviation)| deviation),
                None => Ok(Deviation::default()),
            };
        }
        SolidOperation::BoundedHalfSpace { .. } | SolidOperation::BoundedHalfSpaceAlongZ { .. } => {
            (
                DeviationPath::BoundedHalfSpace,
                "",
                DeviationBound::Unbounded(
                    "bounded half-space: a finite stand-in for an unbounded solid",
                ),
            )
        }
        _ => (
            DeviationPath::Unreported,
            "",
            DeviationBound::Unbounded("solid family"),
        ),
    };
    Ok(Deviation::one(path, detail, bound))
}

/// A run of sections at stations meshed with its deviation (#247): the
/// mesh `crate::station` samples and the bound it certifies for it, or
/// `None` for any other solid.
pub(crate) fn stationed_spine(
    graph: &GeometryGraph,
    operation: &SolidOperation,
    options: &ExecutionOptions,
) -> Option<GeomResult<(TriMesh, Deviation)>> {
    let run = match operation {
        SolidOperation::StationedSpine {
            directrix,
            sections,
            frame,
        } => {
            let sections: Vec<_> = sections.iter().copied().map(Into::into).collect();
            crate::station::stationed_spine_bounded(graph, *directrix, &sections, *frame, options)
        }
        SolidOperation::SectionsAtStations {
            directrix,
            sections,
            frame,
        } => crate::station::stationed_spine_bounded(graph, *directrix, sections, *frame, options),
        _ => return None,
    };
    Some(run.map(|(mesh, bound)| {
        (
            mesh,
            Deviation::one(DeviationPath::StationedSpine, "", bound),
        )
    }))
}

/// The deviation of a sectioned surface (#241, #247): the bound
/// `crate::station` certified for its mesh.
pub(crate) fn of_sectioned_surface(bound: DeviationBound) -> Deviation {
    Deviation::one(DeviationPath::SectionedSurface, "", bound)
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
/// parameters). A curve-relation boundary is bounded by its leaves (#255,
/// `crate::bounded::relation::bound`), its widest joint gap included.
pub(crate) fn of_curve_bounded(
    graph: &GeometryGraph,
    basis: NodeId,
    boundaries: &[NodeId],
    options: &ExecutionOptions,
    merged: Scalar,
) -> Deviation {
    // The boundaries were flattened to the budget shrunk by the stretch
    // (`crate::bounded::boundary_options`), so read them at that budget.
    let (stretch, options) = match graph.get(basis) {
        Some(GeometryNode::Surface(axiolid_surface::Surface::Plane(plane))) => (
            axiolid_reference::bound::frame_stretch3(&plane.frame),
            crate::bounded::boundary_options(&plane.frame, options),
        ),
        _ => (Scalar::INFINITY, options.clone()),
    };
    let options = &options;
    let chord = crate::compiler::chord_error(options);
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
            Some(GeometryNode::CurveRelation(_)) => {
                crate::bounded::relation::bound(graph, id, options)
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

/// Why a tube along an elevated or banked curve under an intrinsic profile
/// is unbounded: the points the certificate reads have no certified error.
pub(crate) const UNCERTIFIED_DIRECTRIX_POINTS: &str = "swept disk: the directrix's points have no certified error (an elevated or banked curve under an intrinsic profile)";

/// Share of the budget the tube's centre line may be off its exact points
/// ([`ElevatedPositionTable`]): the bound carries it, and the rest of the
/// budget is certified against the table's points.
const TUBE_POSITION_SHARE: Scalar = 1e-3;

/// The elevated curve under `curve`, or the one under a banked curve.
fn elevated_base(curve: &Curve3) -> Option<&Elevated3> {
    match curve {
        Curve3::Elevated(e) => Some(e),
        Curve3::Banked(b) => Some(&b.base),
        _ => None,
    }
}

/// Whether any piece of `law` is an intrinsic (vertical clothoid) profile.
fn has_intrinsic_profile(law: &ElevationLaw) -> bool {
    match law {
        ElevationLaw::Intrinsic { .. } => true,
        ElevationLaw::Piecewise { laws, .. } => laws.iter().any(has_intrinsic_profile),
        _ => false,
    }
}

/// [`UNCERTIFIED_DIRECTRIX_POINTS`] for a curve whose points have no
/// certified error: an elevated or banked curve under an intrinsic profile
/// (TODO(#481)).
fn uncertified_points(curve: &Curve3) -> Option<&'static str> {
    let base = elevated_base(curve)?;
    has_intrinsic_profile(&base.elevation).then_some(UNCERTIFIED_DIRECTRIX_POINTS)
}

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
    // The centre line of an elevated or banked curve over an intrinsic or a
    // chain plan is read from a certified table, whose error the bound
    // carries (and, over a chain plan, the error of its tangent and second
    // derivative: `wall_jet_charge`).
    let positions = match elevated_base(curve) {
        Some(base) if matches!(base.plan.as_ref(), Curve2::Intrinsic(_) | Curve2::Chain(_)) => {
            Some(ElevatedPositionTable::new(
                curve,
                TUBE_POSITION_SHARE * target,
                hi,
            )?)
        }
        _ => None,
    };
    let positions = positions.as_ref();
    let mut error = positions.map_or(0.0, ElevatedPositionTable::max_error);
    if error.is_nan() || error >= target {
        return None;
    }
    let mut jets = 0.0_f64;
    let mut walls: Vec<TubeWall<'_>> = Vec::new();
    for &(a, b) in &pieces {
        let reference = perpendicular(unit_tangent(curve, 0.5 * (a + b))?);
        if let Some(table) = positions {
            jets = jets.max(piece_jet_charge(curve, table, (a, b), radius)?);
        }
        walls.push(TubeWall {
            curve,
            positions,
            span: (a, b),
            reference,
            radius,
        });
        if let Some(r) = inner {
            walls.push(TubeWall {
                curve,
                positions,
                span: (a, b),
                reference,
                radius: r,
            });
        }
    }
    let r0 = inner.unwrap_or(0.0);
    // The end disks stand square to the table's tangent: a tangent off by
    // `2 e1 / |c'|` turns each disk, and moves its rim by `radius` times it.
    let mut caps_charge = 0.0_f64;
    let caps: Vec<Disk> = [(lo, SeamSide::Outgoing), (hi, SeamSide::Incoming)]
        .iter()
        .map(|&(s, side)| {
            let (point, t) = match positions {
                Some(table) => {
                    let (point, c1, _) = table.jet3_on(s, side)?;
                    let speed = c1.length();
                    if !(speed > 0.0 && speed.is_finite()) {
                        return None;
                    }
                    caps_charge = caps_charge.max(radius * 2.0 * table.jet_errors()[1] / speed);
                    (point, c1 / speed)
                }
                None => (
                    axiolid_reference::curve::evaluate3(curve, s).ok()?,
                    unit_tangent(curve, s)?,
                ),
            };
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
    error = tube_error(error, jets, caps_charge);
    if !error.is_finite() || error >= target {
        return None;
    }
    // A long tube's quads are long and thin and slightly twisted: joined
    // as flat regions up to a twentieth of the target thick, and cut
    // across their long edges, the search need not resolve every quad's
    // diagonal down to the bound's size (#252). It settles once within the
    // target: a long tube has too many cells to tighten further.
    certify_thin_quads(
        &patches,
        cells,
        mesh,
        target - error,
        target - error,
        TUBE_FLATNESS * target,
    )
    .map(|certified| (certified + error) * (1.0 + axiolid_reference::bound::ROUNDING))
}

/// What the tube's bound carries beyond the certificate: the centre's error
/// and the larger of the walls' jet charge and the end disks' tilt.
fn tube_error(centre: Scalar, jets: Scalar, caps: Scalar) -> Scalar {
    centre + jets.max(caps)
}

/// The charge for one wall piece `(a, b)` read from a table whose tangent and
/// second derivative carry errors (see [`wall_jet_charge`]): zero over an
/// intrinsic plan, whose derivatives are exact closed forms.
fn piece_jet_charge(
    curve: &Curve3,
    table: &ElevatedPositionTable,
    (a, b): (Scalar, Scalar),
    radius: Scalar,
) -> Option<Scalar> {
    let errors = table.jet_errors();
    if errors[1] == 0.0 && errors[2] == 0.0 {
        return Some(0.0);
    }
    let half = 0.5 * (b - a);
    let d2 = axiolid_reference::bound::curve_derivative_bounds3(curve, a, b)?.second;
    let (_, c1, _) = table.jet3_on(0.5 * (a + b), SeamSide::Outgoing)?;
    wall_jet_charge(errors, radius, d2, c1.length() - half * d2, half)
}

/// How far the wall `T(s, theta)` of radius `radius` built from the jets of
/// a table can lie from the exact wall, beyond the centre error `e0`, over
/// a piece of half-width `half`, when the table's first and second
/// derivatives are within `e1` and `e2` of the exact ones (`errors = [e0,
/// e1, e2]`); `d2` bounds `|c''|` and `sigma` the exact speed from below
/// before the table's error (`|c'| - half d2` at the middle, less `e1`).
///
/// With `sigma` the speed floor of the table's own tangent, `k1 = d2 /
/// sigma` bounds `|t'|` and the unit tangent is off by `et = 2 e1 / sigma`
/// (`|a/|a| - b/|b|| <= 2 |a - b| / |a|`), its derivative by `et1 = (2 e2 +
/// 2 d2 et + k1 e1) / sigma`. The reference's projection `P = a - (a.t) t`
/// then stays `p_min = cos tau` long (`tau` the turning of the tangent from
/// its middle value plus `et`), `n = P / |P|` is off by `4 et / p_min` and
/// `b = t x n` by `et + en`, and their derivatives by the same rules with
/// `|P'| <= 2 k1`, `m1 = 2 k1 / p_min`. The wall's point and its
/// `theta`-partial are off by `r (en + eb)`, its `s`-partial by `e1 + r
/// (en1 + eb1)`; a cell has `h_s <= half` and `h_theta <= pi`, and the factor 2
/// covers the products of two `e` (each below `1e-3`) the first-order
/// derivation drops. `None` where the speed floor, the turning or `et` leave
/// the region the derivation holds in.
fn wall_jet_charge(
    errors: [Scalar; 3],
    radius: Scalar,
    d2: Scalar,
    sigma: Scalar,
    half: Scalar,
) -> Option<Scalar> {
    let [_, e1, e2] = errors;
    if e1 == 0.0 && e2 == 0.0 {
        return Some(0.0);
    }
    let (point, s_partial) = wall_jet_errors(errors, radius, d2, sigma, half)?;
    let out = 2.0 * (point * (1.0 + core::f64::consts::PI) + s_partial * half);
    out.is_finite().then_some(out)
}

/// The first-order errors `wall_jet_charge` sums: how far the wall's point
/// and `theta`-partial (`r (en + eb)`) and its `s`-partial (`e1 + r (en1 +
/// eb1)`) built from jets within `errors` can lie from the exact ones.
fn wall_jet_errors(
    errors: [Scalar; 3],
    radius: Scalar,
    d2: Scalar,
    sigma: Scalar,
    half: Scalar,
) -> Option<(Scalar, Scalar)> {
    let [_, e1, e2] = errors;
    let sigma = sigma - e1;
    if !(sigma > 0.0 && sigma.is_finite()) {
        return None;
    }
    let k1 = d2 / sigma;
    let et = 2.0 * e1 / sigma;
    let tau = k1 * half + et;
    if !(tau < 1.5 && et <= 1e-3) {
        return None;
    }
    let p_min = tau.cos();
    let et1 = (2.0 * e2 + 2.0 * d2 * et + k1 * e1) / sigma;
    let en = 4.0 * et / p_min;
    let eb = et + en;
    let m1 = 2.0 * k1 / p_min;
    let ep1 = 2.0 * k1 * et + 2.0 * et1;
    let en1 = (2.0 * ep1 + 4.0 * k1 * en + 2.0 * m1 * et) / p_min;
    let eb1 = et1 + k1 * en + m1 * et + en1;
    let out = (radius * (en + eb), e1 + radius * (en1 + eb1));
    (out.0.is_finite() && out.1.is_finite()).then_some(out)
}

/// Fraction of the target a tube quad may be twisted and still be covered
/// as one flat region; its thickness is added to the bound.
const TUBE_FLATNESS: Scalar = 0.05;

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
    /// Certified centres within the table's `max_error` of the exact ones,
    /// for a curve over an intrinsic plan: the wall is the exact tube moved
    /// by at most that, since its frame is closed form. Over a chain plan
    /// the tangent and the second derivative are read from the table too,
    /// and `wall_jet_charge` charges their errors.
    positions: Option<&'c ElevatedPositionTable>,
    /// The plan distances the wall covers: a corner read at its end is read
    /// as the piece that ends there.
    span: (Scalar, Scalar),
    reference: Vec3,
    radius: Scalar,
}

impl TubeWall<'_> {
    /// The centre, its derivative and its second derivative at `s`.
    fn centre_jet(&self, s: Scalar) -> Option<(Point3, Vec3, Vec3)> {
        match self.positions {
            Some(table) => {
                let side = if s == self.span.1 && s != self.span.0 {
                    SeamSide::Incoming
                } else {
                    SeamSide::Outgoing
                };
                table.jet3_on(s, side)
            }
            None => Some((
                axiolid_reference::curve::evaluate3(self.curve, s).ok()?,
                axiolid_reference::curve::derivative3(self.curve, s).ok()?,
                axiolid_reference::curve::second_derivative3(self.curve, s).ok()?,
            )),
        }
    }
}

/// The wall's point and its partials in `s` and `theta` at `theta`, from the
/// centre's point `c`, derivative `c1` and second derivative `c2` and the
/// fixed `reference`.
fn wall_jet(
    reference: Vec3,
    radius: Scalar,
    (c, c1, c2): (Point3, Vec3, Vec3),
    theta: Scalar,
) -> Option<(Point3, Vec3, Vec3)> {
    let sigma = c1.length();
    if sigma <= 0.0 || sigma.is_nan() {
        return None;
    }
    let t = c1 / sigma;
    let t1 = (c2 - t * c2.dot(t)) / sigma;
    let a = reference;
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
    let r = radius;
    let point = c + (n * cos + b * sin) * r;
    let ds = c1 + (n1 * cos + b1 * sin) * r;
    let dtheta = (b * cos - n * sin) * r;
    (point.is_finite() && ds.is_finite()).then_some((point, ds, dtheta))
}

impl Patch for TubeWall<'_> {
    fn jet(&self, s: Scalar, theta: Scalar) -> Option<(Point3, Vec3, Vec3)> {
        wall_jet(self.reference, self.radius, self.centre_jet(s)?, theta)
    }

    fn second(&self, s: (Scalar, Scalar), _theta: (Scalar, Scalar)) -> Option<[Scalar; 3]> {
        let bounds = axiolid_reference::bound::curve_derivative_bounds3(self.curve, s.0, s.1)?;
        let (d2, d3) = (bounds.second, bounds.third);
        let mid = 0.5 * (s.0 + s.1);
        let half = 0.5 * (s.1 - s.0).abs();
        let (_, c1, _) = self.centre_jet(mid)?;
        // The table's tangent may be off by `e1`: the exact speed is above
        // `|c'| - e1`, and the exact tangent within `2 e1 / sigma` of it.
        let e1 = self.positions.map_or(0.0, |table| table.jet_errors()[1]);
        let sigma = c1.length() - half * d2 - e1;
        if sigma <= 0.0 || sigma.is_nan() {
            return None;
        }
        let k1 = d2 / sigma;
        let k2 = 2.0 * d3 / sigma + 6.0 * d2 * d2 / (sigma * sigma);
        let along =
            (self.reference.dot(c1 / c1.length()).abs() + k1 * half + 2.0 * e1 / sigma).min(1.0);
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
                positions: None,
                span: (a, b),
                reference: perpendicular(unit_tangent(&curve, 0.5 * (a + b)).unwrap()),
                radius: 1.0,
            };
            assert_second_bounded(&wall, (a, b), (0.0, core::f64::consts::TAU));
        }
    }

    #[test]
    fn a_tube_wall_over_a_clothoid_plan_reads_certified_points() {
        use axiolid_core::{Frame2, Point2, Vec2};
        use axiolid_curve::{CurvatureLaw, Elevated3, Intrinsic2};
        let frame = Frame2 {
            origin: Point2::new(20.0, 10.0),
            x: Vec2::new(0.8, 0.6),
            y: Vec2::new(-0.6, 0.8),
        };
        let plan = Curve2::Intrinsic(Intrinsic2::new(
            frame,
            CurvatureLaw::clothoid(0.0, 1.0 / 120.0, 60.0),
            60.0,
        ));
        let curve = Curve3::Elevated(Elevated3::new(
            plan,
            ElevationLaw::parabolic(50.0, 0.05, -0.03, 60.0),
        ));
        let table = ElevatedPositionTable::new(&curve, 1e-9, 60.0).unwrap();
        let wall = TubeWall {
            curve: &curve,
            positions: Some(&table),
            span: (40.0, 42.0),
            reference: perpendicular(unit_tangent(&curve, 41.0).unwrap()),
            radius: 0.0,
        };
        let read = wall.jet(41.0, 0.0).unwrap().0;
        assert_eq!(read, table.point3(41.0).unwrap().0);
        // Not vacuous: the evaluator's own point differs from the table's.
        let evaluated = axiolid_reference::curve::evaluate3(&curve, 41.0).unwrap();
        assert_ne!(read, evaluated);
        assert!(read.distance(evaluated) <= table.max_error() + 1e-9);
    }

    fn chain_frame() -> axiolid_core::Frame2 {
        use axiolid_core::{Point2, Vec2};
        axiolid_core::Frame2 {
            origin: Point2::new(20.0, 10.0),
            x: Vec2::new(0.8, 0.6),
            y: Vec2::new(-0.6, 0.8),
        }
    }

    /// A tube wall over an elevated harmonic chain plan (a straight 20 m, then
    /// `k = 0.01 sin 50 s`) reads the chain table's centre: the evaluator's
    /// point is 1.4e-2 off at 116.5 (#438).
    #[test]
    fn a_tube_wall_over_a_harmonic_chain_plan_reads_certified_points() {
        use axiolid_curve::{Chain2, ChainPiece2, CurvatureLaw, Elevated3};
        let plan = Curve2::Chain(Chain2::new(
            chain_frame(),
            vec![
                ChainPiece2::Intrinsic {
                    curvature: CurvatureLaw::Constant { curvature: 0.0 },
                    length: 20.0,
                },
                ChainPiece2::Intrinsic {
                    curvature: CurvatureLaw::Sinusoid {
                        mean: 0.0,
                        amplitude: 0.01,
                        angular_frequency: 50.0,
                        phase: 0.0,
                    },
                    length: 100.0,
                },
            ],
        ));
        let curve = Curve3::Elevated(Elevated3::new(
            plan,
            ElevationLaw::constant_grade(50.0, 0.04),
        ));
        let table = ElevatedPositionTable::new(&curve, 1e-6, 119.0).unwrap();
        let wall = TubeWall {
            curve: &curve,
            positions: Some(&table),
            span: (116.0, 117.0),
            reference: perpendicular(unit_tangent(&curve, 116.5).unwrap()),
            radius: 0.0,
        };
        let read = wall.jet(116.5, 0.0).unwrap().0;
        assert_eq!(read, table.point3(116.5).unwrap().0);
        let evaluated = axiolid_reference::curve::evaluate3(&curve, 116.5).unwrap();
        assert!(
            read.distance(evaluated) > 1e-2,
            "{}",
            read.distance(evaluated)
        );
    }

    /// `y = a x^3` by its binomial series: arc length from 0 to `x`.
    fn cubic_series(a: Scalar, x: Scalar) -> Scalar {
        let z = 9.0 * a * a * x.powi(4);
        let (mut binomial, mut power, mut sum) = (1.0, 1.0, 0.0);
        for n in 0..200 {
            let term = binomial * power / (4.0 * Scalar::from(n) + 1.0);
            sum += term;
            if term.abs() < 1e-30 {
                break;
            }
            binomial *= (0.5 - Scalar::from(n)) / (Scalar::from(n) + 1.0);
            power *= z;
        }
        x * sum
    }

    /// A wall built from a coarse table's jets over a cubic chain plan stays
    /// within the charge of the exact series wall, and the charge is not
    /// zero (#438).
    #[test]
    fn a_wall_from_a_coarse_chain_table_stays_within_its_charge() {
        use axiolid_core::{Point3, Vec2};
        use axiolid_curve::{BSplineCurve2, Chain2, ChainPiece2, Elevated3};
        let (reach, a) = (40.0, 1.0 / (6.0 * 300.0 * 40.0));
        let cubic = Curve2::BSpline(BSplineCurve2 {
            degree: 3,
            control_points: vec![
                axiolid_core::Point2::new(0.0, 0.0),
                axiolid_core::Point2::new(reach / 3.0, 0.0),
                axiolid_core::Point2::new(2.0 * reach / 3.0, 0.0),
                axiolid_core::Point2::new(reach, a * reach * reach * reach),
            ],
            knots: vec![0.0, 1.0],
            multiplicities: vec![4, 4],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        });
        let frame = chain_frame();
        let curve = Curve3::Elevated(Elevated3::new(
            Curve2::Chain(Chain2::new(
                frame,
                vec![ChainPiece2::Parametric {
                    curve: cubic,
                    start: 0.0,
                    length: 40.0,
                }],
            )),
            ElevationLaw::constant_grade(50.0, 0.04),
        ));
        // Coarse: the enclosures are wide, so the tangent and curvature
        // errors are large.
        let table = ElevatedPositionTable::new(&curve, 1e-3, 39.0).unwrap();
        let [e0, e1, e2] = table.jet_errors();
        assert!(e1 > 1e-9 && e2 > 1e-9, "{e1:e} {e2:e}");
        let (lo, hi, radius) = (10.0, 30.0, 0.3);
        let charge = piece_jet_charge(&curve, &table, (lo, hi), radius).unwrap();
        assert!(charge > 0.0);
        let wall = TubeWall {
            curve: &curve,
            positions: Some(&table),
            span: (lo, hi),
            reference: perpendicular(unit_tangent(&curve, 0.5 * (lo + hi)).unwrap()),
            radius,
        };
        // The exact centre jet at plan distance `s`, from the series.
        let exact = |s: Scalar| {
            let mut x = s;
            for _ in 0..60 {
                x -= (cubic_series(a, x) - s) / (1.0 + 9.0 * a * a * x.powi(4)).sqrt();
            }
            let slope = 3.0 * a * x * x;
            let heading = slope.atan();
            let kappa = 6.0 * a * x / (1.0 + slope * slope).powf(1.5);
            let world = |v: Vec2| frame.x * v.x + frame.y * v.y;
            let p = world(Vec2::new(x, a * x * x * x));
            let t = world(Vec2::new(heading.cos(), heading.sin()));
            let n = world(Vec2::new(-heading.sin(), heading.cos())) * kappa;
            (
                Point3::new(frame.origin.x + p.x, frame.origin.y + p.y, 50.0 + 0.04 * s),
                Vec3::new(t.x, t.y, 0.04),
                Vec3::new(n.x, n.y, 0.0),
            )
        };
        let mut worst = (0.0_f64, 0.0_f64, 0.0_f64);
        for i in 0..=20 {
            let s = lo + (hi - lo) * Scalar::from(i) / 20.0;
            for j in 0..8 {
                let theta = core::f64::consts::TAU * Scalar::from(j) / 8.0;
                let read = wall.jet(s, theta).unwrap();
                let truth = wall_jet(wall.reference, radius, exact(s), theta).unwrap();
                worst.0 = worst.0.max(read.0.distance(truth.0));
                worst.1 = worst.1.max((read.1 - truth.1).length());
                worst.2 = worst.2.max((read.2 - truth.2).length());
            }
        }
        assert!(
            worst.0 > 0.0 && worst.1 > 0.0,
            "{worst:?}: nothing to charge"
        );
        assert!(
            worst.0 <= e0 + charge,
            "point {:e} over {:e}",
            worst.0,
            e0 + charge
        );
        assert!(worst.1 <= charge, "s partial {:e} over {charge:e}", worst.1);
        assert!(
            worst.2 <= charge,
            "theta partial {:e} over {charge:e}",
            worst.2
        );
        // The bound carries the larger of the walls' and the disks' charge.
        let speed = table.jet3_on(hi, SeamSide::Incoming).unwrap().1.length();
        let caps = radius * 2.0 * e1 / speed;
        let carried = tube_error(e0, charge, caps);
        assert!(carried >= e0 + charge && carried >= e0 + caps && carried > e0);
        // Each term is the one carried where it is the larger.
        assert_eq!(tube_error(e0, 0.0, caps), e0 + caps);
        assert_eq!(tube_error(e0, charge, 0.0), e0 + charge);
        assert_eq!(tube_error(e0, 2.0, 1.0), e0 + 2.0);
        assert_eq!(tube_error(e0, 1.0, 2.0), e0 + 2.0);
    }

    /// Over an intrinsic plan the table's derivatives are exact closed forms:
    /// the wall is charged its centre's error only (#438).
    #[test]
    fn a_wall_over_an_intrinsic_plan_charges_only_its_centre() {
        assert_eq!(
            wall_jet_charge([1e-9, 0.0, 0.0], 0.3, 0.05, 0.9, 1.0),
            Some(0.0)
        );
    }

    /// The charge grows with each of the derivative errors, with the radius
    /// and with the piece's width, is first order in them, and refuses the
    /// regions the derivation does not hold in.
    #[test]
    fn a_wall_over_a_chain_plan_is_charged_its_jet_errors() {
        let charge = |e1, e2, r, half| wall_jet_charge([1e-9, e1, e2], r, 0.01, 0.9, half);
        let base = charge(1e-9, 1e-9, 0.3, 1.0).unwrap();
        assert!(base > 0.0 && base < 1e-6, "{base:e}");
        assert!(charge(2e-9, 1e-9, 0.3, 1.0).unwrap() > base);
        assert!(charge(1e-9, 2e-9, 0.3, 1.0).unwrap() > base);
        assert!(charge(1e-9, 1e-9, 0.6, 1.0).unwrap() > base);
        assert!(charge(1e-9, 1e-9, 0.3, 2.0).unwrap() > base);
        // Doubling both errors roughly doubles it.
        let double = charge(2e-9, 2e-9, 0.3, 1.0).unwrap();
        assert!((1.9..2.2).contains(&(double / base)), "{}", double / base);
        // A tangent error that leaves no speed, or a tangent known to
        // worse than a thousandth, is not charged: it is refused.
        assert!(charge(0.95, 1e-9, 0.3, 1.0).is_none());
        assert!(charge(1e-2, 1e-9, 0.3, 1.0).is_none());
        // A piece that turns the tangent past the derivation's range.
        assert!(wall_jet_charge([1e-9, 1e-9, 1e-9], 0.3, 2.0, 0.9, 1.0).is_none());
    }

    /// The first-order errors hold against jets perturbed by the worst the
    /// table may be off: adversarial directions of `c'` and `c''` (along and
    /// across the tangent, and pseudo-random), a quarter turn of a helix-like
    /// centre line, a reference not far from the tangent.
    #[test]
    fn a_wall_built_from_perturbed_jets_stays_within_its_charge() {
        let (e1, e2, radius) = (3e-6, 2e-5, 0.4);
        // c(s) = (cos 0.3 s / 0.3, sin 0.3 s / 0.3, 0.1 s): |c'| = 1.005,
        // |c''| = 0.3.
        let centre = |s: Scalar| {
            (
                Point3::new((0.3 * s).cos() / 0.3, (0.3 * s).sin() / 0.3, 0.1 * s),
                Vec3::new(-(0.3 * s).sin(), (0.3 * s).cos(), 0.1),
                Vec3::new(-0.3 * (0.3 * s).cos(), -0.3 * (0.3 * s).sin(), 0.0),
            )
        };
        let (a, b) = (1.0, 2.5);
        let half = 0.5 * (b - a);
        let d2 = 0.3;
        let speed = centre(0.5 * (a + b)).1.length();
        let errors = [0.0, e1, e2];
        let (point, s_partial) = wall_jet_errors(errors, radius, d2, speed - half * d2, half)
            .expect("inside the derivation");
        assert!(point > 0.0 && s_partial > e1);
        let reference = perpendicular(centre(0.5 * (a + b)).1.normalize());
        let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
        let mut random = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as Scalar / (1_u64 << 53) as Scalar * 2.0 - 1.0
        };
        let mut worst = (0.0_f64, 0.0_f64, 0.0_f64);
        for i in 0..=30 {
            let s = a + (b - a) * Scalar::from(i) / 30.0;
            let exact = centre(s);
            let t = exact.1.normalize();
            let across = t.cross(Vec3::Z).normalize();
            for (u, v) in [
                (t, t),
                (t, -t),
                (across, across),
                (across, -across),
                (-t, across),
                (
                    Vec3::new(random(), random(), random()).normalize(),
                    Vec3::new(random(), random(), random()).normalize(),
                ),
            ] {
                let read = (exact.0, exact.1 + u * e1, exact.2 + v * e2);
                for j in 0..=12 {
                    let theta = core::f64::consts::TAU * Scalar::from(j) / 12.0;
                    let x = wall_jet(reference, radius, exact, theta).unwrap();
                    let y = wall_jet(reference, radius, read, theta).unwrap();
                    worst.0 = worst.0.max(x.0.distance(y.0));
                    worst.1 = worst.1.max((x.2 - y.2).length());
                    worst.2 = worst.2.max((x.1 - y.1).length());
                }
            }
        }
        assert!(
            worst.0 > 0.0 && worst.2 > 0.0,
            "the perturbation did nothing"
        );
        assert!(worst.0 <= point, "point {:e} above {point:e}", worst.0);
        assert!(
            worst.1 <= point,
            "theta partial {:e} above {point:e}",
            worst.1
        );
        assert!(
            worst.2 <= s_partial,
            "s partial {:e} above {s_partial:e}",
            worst.2
        );
        // And the charge is not slack by orders of magnitude.
        assert!(point <= 1e3 * worst.0 && s_partial <= 1e3 * worst.2);
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
