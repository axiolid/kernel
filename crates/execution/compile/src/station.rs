//! Stations along a graph's curves, resolved (#241, ADR 0082).
//!
//! A [`CurveStation`] node, an [`OffsetByStations`] curve, a
//! [`StationedSpine`] solid and a [`SectionedSurface`] all name positions
//! by a distance along a basis curve. This module resolves them through
//! `axiolid_reference::station`, whose module documentation states the
//! distance and frame conventions and the accuracy every resolved point
//! carries (the arc-length inverse's `ARC_LENGTH_TOLERANCE` where the
//! measure is numerical).
//!
//! The basis must be an atomic 2D or 3D curve node; an instance or a curve
//! relation is refused by name, since its measure would be the relation's,
//! not its basis's.
//!
//! # Meshing between stations
//!
//! Everything between two stations is interpolated linearly in distance:
//! the offsets, and each section point with its counterpart. The meshed
//! sections stand at the stations and at distances bisected between them
//! until every placed point at an interval's midpoint is within the chord
//! budget of the midpoint of its chord. That test is a sample, not a
//! proof, so these paths report [`DeviationBound::Unbounded`] by name.
//!
//! [`OffsetByStations`]: axiolid_model::CurveRelation::OffsetByStations
//! [`StationedSpine`]: axiolid_model::SolidOperation::StationedSpine
//! [`SectionedSurface`]: axiolid_model::SurfaceRelation::SectionedSurface
//! [`DeviationBound::Unbounded`]: crate::DeviationBound::Unbounded

use axiolid_construct::loft::{loft_tapered, Station as LoftStation};
use axiolid_construct::profile::{profile_rings, Rings};
use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Frame3, Point2, Point3, Scalar, Tolerance};
use axiolid_curve::{Curve2, Curve3};
use axiolid_mesh::TriMesh;
use axiolid_model::{
    CurveRelation, CurveStation, GeometryGraph, GeometryNode, NodeId, Station, StationFrame,
    StationOffsets, StationedOpenSection, StationedSection,
};
use axiolid_reference::station::{station_section2, station_section3, SectionFrame};

/// Most sections one interval between two stations is bisected into.
const MAX_DEPTH: u32 = 12;

/// Most sections one stationed sweep or offset curve may mesh.
const MAX_SECTIONS: usize = 1 << 16;

/// A station resolved to a point and a frame.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResolvedStation {
    /// The station's point: the curve point moved by the offsets.
    pub point: Point3,
    /// The frame at the station, at `point`: `x` the tangent, `y` up,
    /// `z` to the right (`-lateral`), the curve-evaluation provider's
    /// layout.
    pub frame: Frame3,
    /// The basis curve's section frame at the distance, before the
    /// offsets, in the requested [`StationFrame`].
    pub section: SectionFrame,
}

/// Resolve a [`CurveStation`] node to its point and frame.
///
/// See the [module documentation](self) and
/// `axiolid_reference::station` for the conventions and the accuracy.
///
/// # Errors
///
/// A node that is not a [`GeometryNode::CurveStation`], a basis that is not
/// an atomic curve, a distance beyond the basis curve's length, and every
/// refusal of the curve evaluators, by name.
pub fn resolve(graph: &GeometryGraph, id: NodeId) -> GeomResult<ResolvedStation> {
    let Some(GeometryNode::CurveStation(CurveStation {
        basis,
        station,
        frame,
    })) = graph.get(id)
    else {
        return Err(GeomError::InvalidInput(format!(
            "node {id:?} is not a curve station"
        )));
    };
    let basis = Basis::of(graph, *basis)?;
    let section = basis.section(station.distance, *frame)?;
    let point = place(&section, &station.offsets, Point2::ZERO);
    let mut placed = section.frame();
    placed.origin = point;
    Ok(ResolvedStation {
        point,
        frame: placed,
        section,
    })
}

/// An atomic basis curve.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Basis<'g> {
    Two(&'g Curve2),
    Three(&'g Curve3),
}

impl<'g> Basis<'g> {
    pub(crate) fn of(graph: &'g GeometryGraph, id: NodeId) -> GeomResult<Self> {
        match graph.get(id) {
            Some(GeometryNode::Curve2(curve)) => Ok(Self::Two(curve)),
            Some(GeometryNode::Curve3(curve)) => Ok(Self::Three(curve)),
            Some(GeometryNode::Instance(_) | GeometryNode::CurveRelation(_)) => {
                Err(GeomError::UnsupportedInput {
                    backend: crate::BACKEND_ID,
                    operation: Operation::CurveEvaluation,
                    input: "a station along an instanced curve or a curve relation",
                })
            }
            Some(_) => Err(GeomError::InvalidInput(format!(
                "station basis {id:?} is not a curve"
            ))),
            None => Err(GeomError::InvalidInput(format!(
                "station basis {id:?} is outside the graph"
            ))),
        }
    }

    /// The section frame at `distance`, in `frame`.
    pub(crate) fn section(self, distance: Scalar, frame: StationFrame) -> GeomResult<SectionFrame> {
        let section = match self {
            Self::Two(curve) => station_section2(curve, distance)?,
            Self::Three(curve) => station_section3(curve, distance)?,
        };
        match frame {
            StationFrame::Section => Ok(section),
            StationFrame::Plan => section.plan(),
            _ => Err(GeomError::UnsupportedInput {
                backend: crate::BACKEND_ID,
                operation: Operation::CurveEvaluation,
                input: "a station frame this compiler does not know",
            }),
        }
    }
}

/// A section point `local` (profile `x` lateral, `y` up) placed with the
/// station offsets.
fn place(section: &SectionFrame, offsets: &StationOffsets, local: Point2) -> Point3 {
    section.place(
        local.x + offsets.lateral,
        local.y + offsets.vertical,
        offsets.longitudinal,
    )
}

fn lerp(a: Scalar, b: Scalar, u: Scalar) -> Scalar {
    a + (b - a) * u
}

fn lerp_offsets(a: &StationOffsets, b: &StationOffsets, u: Scalar) -> StationOffsets {
    StationOffsets::new(
        lerp(a.lateral, b.lateral, u),
        lerp(a.vertical, b.vertical, u),
        lerp(a.longitudinal, b.longitudinal, u),
    )
}

/// Whether the meshed span between two sections `lo` and `hi` passes the
/// sampling test: each placed point at the midpoint `mid` within `chord`
/// of its chord's midpoint, and each wall quad between neighbours `(i, j)`
/// within `chord` of flat ([`quad_spread`]), so the two triangles the loft
/// cuts it into stay near the bilinear patch.
fn passes(
    lo: &[Point3],
    hi: &[Point3],
    mid: &[Point3],
    edges: &[(usize, usize)],
    chord: Scalar,
) -> bool {
    let straight = lo
        .iter()
        .zip(hi)
        .zip(mid)
        .all(|((p, q), m)| (*m - 0.5 * (*p + *q)).length() <= chord);
    straight
        && edges
            .iter()
            .all(|&(i, j)| quad_spread(lo[i], lo[j], hi[i], hi[j]) <= chord)
}

/// How far the wall quad `a -> b` over `c -> d` strays from flat: the
/// smaller of its twist `|a - b - c + d| / 4` and the height of `d` above
/// the plane through `a, b, c`. A planar quad (a linear taper along a
/// straight directrix) needs no refinement.
fn quad_spread(a: Point3, b: Point3, c: Point3, d: Point3) -> Scalar {
    let twist = 0.25 * (a - b - c + d).length();
    let normal = (b - a).cross(c - a).normalize_or_zero();
    if normal == axiolid_core::Vec3::ZERO {
        twist
    } else {
        twist.min(normal.dot(d - a).abs())
    }
}

/// Distances at which an interval `[a, b]` is meshed, `a` included and `b`
/// excluded: bisected until [`passes`] holds, at most [`MAX_DEPTH`] deep.
fn refine<F>(
    a: Scalar,
    b: Scalar,
    chord: Scalar,
    edges: &[(usize, usize)],
    eval: &F,
    out: &mut Vec<Scalar>,
) -> GeomResult<()>
where
    F: Fn(Scalar) -> GeomResult<Vec<Point3>>,
{
    let start = eval(a)?;
    let end = eval(b)?;
    // Depth-first, left half first, so distances come out increasing.
    let mut stack = vec![(a, b, start, end, 0u32)];
    while let Some((lo, hi, at_lo, at_hi, depth)) = stack.pop() {
        let mid = 0.5 * (lo + hi);
        let at_mid = eval(mid)?;
        let straight = passes(&at_lo, &at_hi, &at_mid, edges, chord);
        if straight || depth >= MAX_DEPTH || !(mid > lo && mid < hi) {
            out.push(lo);
        } else {
            stack.push((mid, hi, at_mid.clone(), at_hi, depth + 1));
            stack.push((lo, mid, at_lo, at_mid, depth + 1));
        }
        if out.len() > MAX_SECTIONS {
            return Err(GeomError::BudgetExceeded {
                resource: "station sections",
            });
        }
    }
    Ok(())
}

/// Every distance a run of stations is meshed at, first and last included.
fn sample_run<F>(
    distances: &[Scalar],
    chord: Scalar,
    edges: &[(usize, usize)],
    eval: F,
) -> GeomResult<Vec<Scalar>>
where
    F: Fn(Scalar) -> GeomResult<Vec<Point3>>,
{
    let mut out = Vec::new();
    for pair in distances.windows(2) {
        refine(pair[0], pair[1], chord, edges, &eval, &mut out)?;
    }
    if let Some(last) = distances.last() {
        out.push(*last);
    }
    Ok(out)
}

/// The station interval `[i, i + 1]` containing `s` and the fraction
/// along it.
fn interval(distances: &[Scalar], s: Scalar) -> (usize, Scalar) {
    let last = distances.len() - 2;
    let i = distances
        .windows(2)
        .position(|pair| s <= pair[1])
        .unwrap_or(last)
        .min(last);
    let (a, b) = (distances[i], distances[i + 1]);
    (i, ((s - a) / (b - a)).clamp(0.0, 1.0))
}

/// The point of an offset curve by stations at distance `s`.
fn offset_point(
    basis: Basis<'_>,
    stations: &[Station],
    distances: &[Scalar],
    frame: StationFrame,
    s: Scalar,
) -> GeomResult<Point3> {
    let (i, u) = interval(distances, s);
    let offsets = lerp_offsets(&stations[i].offsets, &stations[i + 1].offsets, u);
    Ok(place(&basis.section(s, frame)?, &offsets, Point2::ZERO))
}

/// An [`CurveRelation::OffsetByStations`] sampled as a sweep directrix:
/// its parameter is the basis distance, so `range` narrows the run of
/// stations.
pub(crate) fn offset_curve_points(
    graph: &GeometryGraph,
    relation: &CurveRelation,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<Vec<Point3>> {
    let CurveRelation::OffsetByStations {
        basis,
        stations,
        frame,
    } = relation
    else {
        return Err(GeomError::InvalidInput(
            "not an offset curve by stations".into(),
        ));
    };
    if stations.len() < 2 {
        return Err(GeomError::InvalidInput(
            "an offset curve by stations needs at least two stations".into(),
        ));
    }
    let basis = Basis::of(graph, *basis)?;
    let distances: Vec<Scalar> = stations.iter().map(|station| station.distance).collect();
    let (first, last) = (distances[0], distances[distances.len() - 1]);
    let mut knots = distances.clone();
    if let Some((a, b)) = range {
        if !(a.is_finite() && b.is_finite()) || a.min(b) < first || a.max(b) > last {
            return Err(GeomError::InvalidInput(format!(
                "range ({a}, {b}) leaves the offset curve's span [{first}, {last}]"
            )));
        }
        let (lo, hi) = (a.min(b), a.max(b));
        knots.retain(|&d| d > lo && d < hi);
        knots.insert(0, lo);
        knots.push(hi);
    }
    let chord = crate::compiler::chord_error(options);
    let at = |s| offset_point(basis, stations, &distances, *frame, s);
    let samples = sample_run(&knots, chord, &[], |s| at(s).map(|p| vec![p]))?;
    let mut points = samples
        .into_iter()
        .map(at)
        .collect::<GeomResult<Vec<_>>>()?;
    if let Some((a, b)) = range {
        if a > b {
            points.reverse();
        }
    }
    Ok(points)
}

/// Linear interpolation of two ring sets of one structure.
fn lerp_rings(a: &Rings, b: &Rings, u: Scalar) -> Rings {
    let mix = |p: &Point2, q: &Point2| Point2::new(lerp(p.x, q.x, u), lerp(p.y, q.y, u));
    Rings {
        outer: a
            .outer
            .iter()
            .zip(&b.outer)
            .map(|(p, q)| mix(p, q))
            .collect(),
        holes: a
            .holes
            .iter()
            .zip(&b.holes)
            .map(|(ha, hb)| ha.iter().zip(hb).map(|(p, q)| mix(p, q)).collect())
            .collect(),
    }
}

fn same_structure(a: &Rings, b: &Rings) -> bool {
    a.outer.len() == b.outer.len()
        && a.holes.len() == b.holes.len()
        && a.holes
            .iter()
            .zip(&b.holes)
            .all(|(p, q)| p.len() == q.len())
}

/// Mesh a [`SolidOperation::StationedSpine`](axiolid_model::SolidOperation::StationedSpine).
///
/// The profiles are flattened to half the chord budget and the stations
/// sampled to the other half, as the other doubly curved sweeps split it.
pub(crate) fn stationed_spine(
    graph: &GeometryGraph,
    directrix: NodeId,
    sections: &[StationedSection],
    frame: StationFrame,
    options: &ExecutionOptions,
) -> GeomResult<TriMesh> {
    if sections.len() < 2 {
        return Err(GeomError::InvalidInput(
            "a station-placed spine needs at least two sections".into(),
        ));
    }
    let basis = Basis::of(graph, directrix)?;
    let half = 0.5 * crate::compiler::chord_error(options);
    let rings = sections
        .iter()
        .map(|section| match graph.get(section.profile) {
            Some(GeometryNode::Profile(shape)) => profile_rings(shape, half, options.tolerance()),
            _ => Err(GeomError::InvalidInput(format!(
                "spine section profile {:?} is not a Profile node",
                section.profile
            ))),
        })
        .collect::<GeomResult<Vec<_>>>()?;
    if rings.iter().any(|ring| !same_structure(ring, &rings[0])) {
        return Err(GeomError::InvalidInput(
            "station-placed spine sections must share their ring structure: vertices are \
             matched by ring and index"
                .into(),
        ));
    }
    let distances: Vec<Scalar> = sections.iter().map(|s| s.station.distance).collect();
    if distances.windows(2).any(|pair| pair[1] <= pair[0]) {
        return Err(GeomError::InvalidInput(
            "station-placed spine distances must increase strictly".into(),
        ));
    }
    let placed = |s: Scalar| -> GeomResult<(Rings, Vec<Point3>)> {
        let (i, u) = interval(&distances, s);
        let local = lerp_rings(&rings[i], &rings[i + 1], u);
        let offsets = lerp_offsets(
            &sections[i].station.offsets,
            &sections[i + 1].station.offsets,
            u,
        );
        let section = basis.section(s, frame)?;
        let points = local
            .outer
            .iter()
            .chain(local.holes.iter().flatten())
            .map(|p| place(&section, &offsets, *p))
            .collect();
        Ok((local, points))
    };
    // Wall edges: each ring closes on itself.
    let mut edges = Vec::new();
    let mut base = 0;
    for len in std::iter::once(rings[0].outer.len()).chain(rings[0].holes.iter().map(Vec::len)) {
        edges.extend((0..len).map(|k| (base + k, base + (k + 1) % len)));
        base += len;
    }
    let samples = sample_run(&distances, half, &edges, |s| {
        placed(s).map(|(_, points)| points)
    })?;
    let mut stations = Vec::with_capacity(samples.len());
    for s in samples {
        let (local, points) = placed(s)?;
        let mut it = points.into_iter();
        let mut loops = Vec::with_capacity(1 + local.holes.len());
        loops.push(it.by_ref().take(local.outer.len()).collect());
        for hole in &local.holes {
            loops.push(it.by_ref().take(hole.len()).collect());
        }
        stations.push(LoftStation { loops });
    }
    loft_tapered(&rings[0], &rings[rings.len() - 1], &stations)
}

/// The open polyline of a sectioned surface's section, checked against
/// its tags.
fn open_section(graph: &GeometryGraph, section: &StationedOpenSection) -> GeomResult<Vec<Point2>> {
    let path = match graph.get(section.profile) {
        Some(GeometryNode::OpenProfile(profile)) => profile.path,
        _ => {
            return Err(GeomError::InvalidInput(format!(
                "sectioned surface section {:?} is not an OpenProfile node",
                section.profile
            )))
        }
    };
    let points = match graph.get(path) {
        Some(GeometryNode::Curve2(Curve2::Polyline(polyline))) if !polyline.closed => {
            polyline.points.clone()
        }
        _ => {
            return Err(GeomError::UnsupportedInput {
                backend: crate::BACKEND_ID,
                operation: Operation::CurveEvaluation,
                input: "a sectioned surface section that is not an open polyline: its tags \
                        name polyline vertices",
            })
        }
    };
    if !section.tags.is_empty() && section.tags.len() != points.len() {
        return Err(GeomError::InvalidInput(format!(
            "a sectioned surface section has {} tags for {} vertices",
            section.tags.len(),
            points.len()
        )));
    }
    Ok(points)
}

/// Mesh a [`SurfaceRelation::SectionedSurface`](axiolid_model::SurfaceRelation::SectionedSurface)
/// as an open sheet.
pub(crate) fn sectioned_surface(
    graph: &GeometryGraph,
    directrix: NodeId,
    sections: &[StationedOpenSection],
    frame: StationFrame,
    options: &ExecutionOptions,
) -> GeomResult<TriMesh> {
    if sections.len() < 2 {
        return Err(GeomError::InvalidInput(
            "a sectioned surface needs at least two sections".into(),
        ));
    }
    let basis = Basis::of(graph, directrix)?;
    let polylines = sections
        .iter()
        .map(|section| open_section(graph, section))
        .collect::<GeomResult<Vec<_>>>()?;
    let count = polylines[0].len();
    if count < 2 || polylines.iter().any(|points| points.len() != count) {
        return Err(GeomError::InvalidInput(
            "sectioned surface sections must have one vertex per tag, at least two, and the \
             same count"
                .into(),
        ));
    }
    if sections
        .iter()
        .any(|section| section.tags != sections[0].tags)
    {
        return Err(GeomError::InvalidInput(
            "sectioned surface sections carry inconsistent tags".into(),
        ));
    }
    let distances: Vec<Scalar> = sections.iter().map(|s| s.station.distance).collect();
    if distances.windows(2).any(|pair| pair[1] <= pair[0]) {
        return Err(GeomError::InvalidInput(
            "sectioned surface distances must increase strictly".into(),
        ));
    }
    let placed = |s: Scalar| -> GeomResult<Vec<Point3>> {
        let (i, u) = interval(&distances, s);
        let offsets = lerp_offsets(
            &sections[i].station.offsets,
            &sections[i + 1].station.offsets,
            u,
        );
        let section = basis.section(s, frame)?;
        Ok(polylines[i]
            .iter()
            .zip(&polylines[i + 1])
            .map(|(p, q)| {
                place(
                    &section,
                    &offsets,
                    Point2::new(lerp(p.x, q.x, u), lerp(p.y, q.y, u)),
                )
            })
            .collect())
    };
    let chord = crate::compiler::chord_error(options);
    let edges: Vec<(usize, usize)> = (0..count - 1).map(|k| (k, k + 1)).collect();
    let samples = sample_run(&distances, chord, &edges, placed)?;
    let mut positions = Vec::with_capacity(samples.len() * count);
    for s in &samples {
        positions.extend(placed(*s)?);
    }
    let mut indices = Vec::with_capacity((samples.len() - 1) * (count - 1) * 6);
    for row in 0..samples.len() - 1 {
        let a = (row * count) as u32;
        let b = ((row + 1) * count) as u32;
        for k in 0..count as u32 - 1 {
            indices.extend([a + k, a + k + 1, b + k + 1, a + k, b + k + 1, b + k]);
        }
    }
    let mesh = TriMesh::new(positions, indices);
    reject_degenerate(&mesh, options.tolerance())?;
    Ok(mesh)
}

/// Refuse a sheet with a collapsed triangle: two sections that place a
/// tag at one point would leave a sliver no consumer can orient.
fn reject_degenerate(mesh: &TriMesh, tolerance: Tolerance) -> GeomResult<()> {
    let floor = tolerance.linear() * tolerance.linear();
    for triangle in mesh.indices.chunks_exact(3) {
        let [a, b, c] = [0, 1, 2].map(|k| mesh.positions[triangle[k] as usize]);
        if (b - a).cross(c - a).length() <= floor {
            return Err(GeomError::Degenerate(
                "a sectioned surface has a collapsed triangle between two sections".into(),
            ));
        }
    }
    Ok(())
}
