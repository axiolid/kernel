//! The sweep families that place a profile along a path.
//!
//! Each differs only in how the profile is carried: a tapered family blends
//! two profiles, a fixed-reference sweep keeps one direction, a
//! surface-curve sweep takes its up vector from a surface normal. The
//! stitching is shared with every other sweep in `loft`.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_mesh::TriMesh;

use crate::loft::{self, Frame, Station};
use crate::profile::Rings;

/// Blend two ring sets at parameter `t`.
///
/// Refuses a structural mismatch rather than truncating: two profiles with
/// different ring counts have no correspondence, and pairing them by index
/// would silently weld unrelated points.
fn blend_rings(a: &Rings, b: &Rings, t: Scalar) -> GeomResult<Rings> {
    if a.outer.len() != b.outer.len() || a.holes.len() != b.holes.len() {
        return Err(GeomError::InvalidInput(
            "tapered sweep profiles must share their ring structure".to_owned(),
        ));
    }
    let mut holes = Vec::with_capacity(a.holes.len());
    for (ha, hb) in a.holes.iter().zip(&b.holes) {
        if ha.len() != hb.len() {
            return Err(GeomError::InvalidInput(
                "tapered sweep holes must share their point count".to_owned(),
            ));
        }
        holes.push(
            ha.iter()
                .zip(hb)
                .map(|(p, q)| loft::blend(*p, *q, t))
                .collect(),
        );
    }
    Ok(Rings {
        outer: a
            .outer
            .iter()
            .zip(&b.outer)
            .map(|(p, q)| loft::blend(*p, *q, t))
            .collect(),
        holes,
    })
}

/// Extrude between two profiles, blending linearly along the direction.
///
/// Two stations suffice: the blend is linear, so intermediate ones would
/// add vertices without adding shape.
pub fn tapered_extrude(
    start: &Rings,
    end: &Rings,
    direction: Vec3,
    depth: Scalar,
) -> GeomResult<TriMesh> {
    if !depth.is_finite() || depth <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "extrusion depth must be positive and finite, got {depth}"
        )));
    }
    let d = direction.normalize_or_zero();
    if d == Vec3::ZERO {
        return Err(GeomError::InvalidInput(
            "extrusion direction must be a non-zero vector".to_owned(),
        ));
    }
    // Caps use the START rings, so the end cap is only correct when both
    // profiles share a structure. blend_rings enforces that.
    let far = blend_rings(start, end, 1.0)?;
    let s0 = loft::place(start, |p| Point3::new(p.x, p.y, 0.0));
    let s1 = loft::place(&far, |p| Point3::new(p.x, p.y, 0.0) + d * depth);
    loft::loft_tapered(start, &far, &[s0, s1])
}

/// Revolve between two profiles.
///
/// Unlike a plain revolution this can never close: a full turn would have
/// to meet the start profile with the end profile, which are different by
/// construction. It is always capped.
///
/// Every point of the surface the blended rings sweep lies within
/// `tolerance.linear()` of the mesh (#231). With `t` in `[0, 1]` along the
/// turn, a ring point moves as `f(t) = R(angle t) q(t)` with `q` affine in
/// `t` (the blend), so `f'' = angle^2 R'' q + 2 angle R' q'`, and
/// `|R'' v|` and `|R' v|` are the distance of `v` from the axis direction,
/// at most `|v|`. Hence `|f''| <= angle^2 rho_max + 2 |angle| d_max`, with
/// `rho_max` the largest distance of a start or end ring vertex from the
/// axis (a convex function of the blend, so largest at an end) and `d_max`
/// the largest distance a vertex travels between the two profiles. Taylor's
/// remainder bounds the distance of `f` from its chord over a step of
/// `1/n` by `|f''|_max / (8 n^2)`, and the chords between two stations span
/// the bilinear patch of each wall quad, which is within
/// `loft::quad_deviation` of its two triangles; the step count
/// grows until the sum of the two fits.
pub fn tapered_revolve(
    start: &Rings,
    end: &Rings,
    axis_origin: Point3,
    axis_direction: Vec3,
    angle: Scalar,
    tolerance: Tolerance,
) -> GeomResult<TriMesh> {
    if !angle.is_finite() || angle == 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "revolution angle must be finite and non-zero, got {angle}"
        )));
    }
    let dir = axis_direction.normalize_or_zero();
    if dir == Vec3::ZERO {
        return Err(GeomError::InvalidInput(
            "revolution axis must be a finite non-zero direction".to_owned(),
        ));
    }
    let far = blend_rings(start, end, 1.0)?;
    let mut max_r: Scalar = 0.0;
    let starts = start.outer.iter().chain(start.holes.iter().flatten());
    let ends = far.outer.iter().chain(far.holes.iter().flatten());
    let mut travel: Scalar = 0.0;
    for (p, q) in starts.zip(ends) {
        for p in [p, q] {
            let v = Point3::new(p.x, p.y, 0.0) - axis_origin;
            max_r = max_r.max((v - dir * dir.dot(v)).length());
        }
        travel = travel.max((*q - *p).length());
    }
    let budget = tolerance.linear();
    let curvature = angle * angle * max_r + 2.0 * angle.abs() * travel;
    let mut n = crate::revolve::steps(max_r, angle, budget)?;
    let stations = loop {
        let mut stations = Vec::with_capacity(n + 1);
        for s in (0..=n).rev() {
            let t = (s as Scalar) / (n as Scalar);
            let ring = blend_rings(start, end, t)?;
            let a = angle * t;
            stations.push(loft::place(&ring, |p| {
                crate::revolve::rotate(Point3::new(p.x, p.y, 0.0), axis_origin, dir, a)
            }));
        }
        let chord = curvature / (8.0 * (n * n) as Scalar);
        let twist = stations
            .windows(2)
            .map(|w| loft::span_deviation(&w[0], &w[1]))
            .fold(0.0, Scalar::max);
        if chord + twist <= budget {
            break stations;
        }
        n = crate::revolve::grow(n)?;
    };
    let mut mesh = loft::loft_tapered(&far, start, &stations)?;
    // Same cause and fix as `revolve::revolve` (#221): the station order
    // above only winds the walls outward for a positive angle, and `angle`
    // and `dir` matter only through their product, so a sign-of-`angle`
    // check alone cannot be made consistent with a sign flip of `dir` too.
    // Settle it from the built mesh's own orientation instead.
    if matches!(crate::extrude::outward_orientation(&mesh), Some(false)) {
        for triangle in mesh.indices.chunks_exact_mut(3) {
            triangle.swap(1, 2);
        }
    }
    Ok(mesh)
}

/// Sweep a profile along a sampled directrix with a fixed reference.
///
/// The reference direction is held constant, so the profile does not twist
/// with the path's torsion. That is what distinguishes this from a Frenet
/// sweep, whose frame rotates with the curve's binormal.
pub fn fixed_reference_sweep(
    rings: &Rings,
    path: &[Point3],
    reference: Vec3,
) -> GeomResult<TriMesh> {
    let frames = frames_along(path, None, |_| reference)?;
    let stations: Vec<Station> = frames
        .iter()
        .map(|f| loft::place(rings, |p| loft::at(f, p)))
        .collect();
    loft::loft(rings, &stations, false)
}

/// A directrix sampled into chords for a sweep that bounds its walls
/// (#231).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SampledPath {
    /// The samples, in the direction of travel.
    pub points: Vec<Point3>,
    /// The exact tangents at the first and the last sample, in the
    /// direction of travel, when the samples chord ONE smooth curve.
    ///
    /// `Some` places the end sections square to the curve rather than to
    /// its end chords, which lean by half the turn of a chord and tilt an
    /// end cap by up to `reach * sin(turn / 2)`, and opts the sweep into
    /// its chord bound. `None` marks a path whose chords are the path
    /// itself (a line or a polyline, whose corners have no curvature to
    /// refine) or one this crate cannot vouch is smooth, such as a
    /// composite; such a path is swept as given. A chain of segments and
    /// arcs is bounded piece by piece by [`crate::pipe`] instead (#232).
    pub end_tangents: Option<[Vec3; 2]>,
}

/// Rounds of halving the directrix budget before a sweep refuses.
const MAX_REFINEMENTS: usize = 16;

/// Stations for a sweep whose walls stay within `budget` of the surface
/// its rings trace along a smooth directrix (#231).
///
/// `frames_at(b)` samples the directrix to chord budget `b` and returns
/// its frames, and whether the sampled path is smooth
/// ([`SampledPath::end_tangents`]); a path that is not is swept once, as
/// given. For a smooth one the bound per span between two stations is
///
/// - the distance of the directrix from its chord: the smaller of `b` and
///   `(L/2) tan(phi/4)`, with `L` the chord and `phi` the turn of the
///   tangent across it, which is exactly the sagitta `R (1 - cos(phi/2))`
///   of a circular arc of radius `R`;
/// - plus `reach (1 - cos(psi/2))`, with `psi` the rotation between the two
///   frames and `reach` the furthest ring point from the directrix: the
///   extra sagitta of a profile point swept round with the frame;
/// - plus [`loft::quad_deviation`] over the span's wall quads.
///
/// Along a circular arc whose frames turn with it (a rotation-minimising
/// frame, or a fixed reference normal to the arc's plane) consecutive
/// stations are rotations of each other about the arc's axis, so the first
/// two terms are the sagitta `(R + reach)(1 - cos(phi/2))` of the circle
/// the furthest profile point runs on, and the bound is the one
/// [`crate::revolve`] proves. For other smooth curves the same terms are
/// the second-order estimate of the sampled curve's own deviation. The
/// budget is halved until every span fits, and refused with
/// `BudgetExceeded` after `MAX_REFINEMENTS` rounds.
fn bounded_stations(
    rings: &Rings,
    budget: Scalar,
    mut frames_at: impl FnMut(Scalar) -> GeomResult<(Vec<Frame>, bool)>,
) -> GeomResult<Vec<Station>> {
    if !(budget.is_finite() && budget > 0.0) {
        return Err(GeomError::InvalidInput(format!(
            "sweep chord budget must be positive and finite, got {budget}"
        )));
    }
    let reach = rings
        .outer
        .iter()
        .chain(rings.holes.iter().flatten())
        .map(|p| p.x.hypot(p.y))
        .fold(0.0, Scalar::max);
    let mut b = budget;
    for _ in 0..MAX_REFINEMENTS {
        let (frames, smooth) = frames_at(b)?;
        let stations: Vec<Station> = frames
            .iter()
            .map(|f| loft::place(rings, |p| loft::at(f, p)))
            .collect();
        if !smooth {
            return Ok(stations);
        }
        let worst = frames
            .windows(2)
            .zip(stations.windows(2))
            .map(|(f, s)| span_bound(&f[0], &f[1], b, reach) + loft::span_deviation(&s[0], &s[1]))
            .fold(0.0, Scalar::max);
        if worst <= budget {
            return Ok(stations);
        }
        b *= 0.5;
    }
    Err(GeomError::BudgetExceeded {
        resource: "sweep directrix refinement",
    })
}

/// The directrix and frame terms of [`bounded_stations`]'s span bound.
pub(crate) fn span_bound(from: &Frame, to: &Frame, budget: Scalar, reach: Scalar) -> Scalar {
    let (t0, t1) = (from.x.cross(from.y), to.x.cross(to.y));
    let turn = t0.dot(t1).clamp(-1.0, 1.0).acos();
    let chord = (to.origin - from.origin).length();
    let directrix = budget.min(0.5 * chord * (0.25 * turn).tan());
    let trace = from.x.dot(to.x) + from.y.dot(to.y) + t0.dot(t1);
    let rotation = (0.5 * (trace - 1.0)).clamp(-1.0, 1.0).acos();
    directrix + reach * (1.0 - (0.5 * rotation).cos())
}

/// [`fixed_reference_sweep`] whose walls stay within `budget` of the
/// surface its rings trace along a smooth directrix (#231).
///
/// `sample(b)` samples the directrix to chord budget `b`; it is called
/// again with a halved budget while the walls miss (see
/// [`SampledPath`]). The rings' own distance from a curved profile is the
/// caller's to add: the reference compiler flattens the profile to half
/// its chord budget and passes the other half here.
pub fn fixed_reference_sweep_within(
    rings: &Rings,
    mut sample: impl FnMut(Scalar) -> GeomResult<SampledPath>,
    reference: Vec3,
    budget: Scalar,
) -> GeomResult<TriMesh> {
    let stations = bounded_stations(rings, budget, |b| {
        let path = sample(b)?;
        let frames = frames_along(&path.points, path.end_tangents, |_| reference)?;
        Ok((frames, path.end_tangents.is_some()))
    })?;
    loft::loft(rings, &stations, false)
}

/// [`surface_curve_sweep`] whose walls stay within `budget` of the surface
/// its rings trace along a smooth directrix (#231).
///
/// `sample(b)` returns the directrix sampled to chord budget `b` and the
/// reference surface's normal at each sample; see
/// [`fixed_reference_sweep_within`].
pub fn surface_curve_sweep_within(
    rings: &Rings,
    mut sample: impl FnMut(Scalar) -> GeomResult<(SampledPath, Vec<Vec3>)>,
    budget: Scalar,
) -> GeomResult<TriMesh> {
    let stations = bounded_stations(rings, budget, |b| {
        let (path, normals) = sample(b)?;
        if normals.len() != path.points.len() {
            return Err(GeomError::InvalidInput(
                "a surface curve sweep needs one surface normal per directrix point".to_owned(),
            ));
        }
        let frames = frames_along(&path.points, path.end_tangents, |i| normals[i])?;
        Ok((frames, path.end_tangents.is_some()))
    })?;
    loft::loft(rings, &stations, false)
}

/// Build a frame at each path sample.
///
/// The tangent at an interior sample is the average of its two segment
/// directions, which keeps the profile from kinking at a corner. Endpoints
/// use the curve's exact tangents when `ends` gives them, else their single
/// adjacent segment.
fn frames_along(
    path: &[Point3],
    ends: Option<[Vec3; 2]>,
    up: impl Fn(usize) -> Vec3,
) -> GeomResult<Vec<Frame>> {
    let tangents = tangents_along(path, ends)?;
    let mut frames = Vec::with_capacity(path.len());
    for (i, (point, tangent)) in path.iter().zip(tangents).enumerate() {
        frames.push(Frame::from_reference(*point, tangent, up(i))?);
    }
    Ok(frames)
}

pub fn linear_extrusion_normals(path: &[Point3], direction: Vec3) -> GeomResult<Vec<Vec3>> {
    if path.len() < 2 {
        return Err(GeomError::InvalidInput(
            "a linear-extrusion surface needs at least two directrix points".into(),
        ));
    }
    let direction = direction.normalize_or_zero();
    if direction == Vec3::ZERO {
        return Err(GeomError::InvalidInput(
            "linear-extrusion direction must be finite and non-zero".into(),
        ));
    }
    let mut normals = Vec::with_capacity(path.len());
    for i in 0..path.len() {
        let tangent = if i == 0 {
            path[1] - path[0]
        } else if i + 1 == path.len() {
            path[i] - path[i - 1]
        } else {
            (path[i] - path[i - 1]).normalize_or_zero()
                + (path[i + 1] - path[i]).normalize_or_zero()
        };
        let normal = tangent.cross(direction).normalize_or_zero();
        if normal == Vec3::ZERO {
            return Err(GeomError::Degenerate(
                "directrix tangent is parallel to linear-extrusion direction".into(),
            ));
        }
        normals.push(normal);
    }
    Ok(normals)
}

/// Sweep a profile along a directrix lying on a reference surface.
///
/// The surface normal at each sample supplies the up direction, so the
/// profile stays oriented to the surface rather than to a global axis.
/// Callers pass the sampled normals because evaluating the surface belongs
/// to the surface provider, not to this sweep.
pub fn surface_curve_sweep(
    rings: &Rings,
    path: &[Point3],
    normals: &[Vec3],
) -> GeomResult<TriMesh> {
    if normals.len() != path.len() {
        return Err(GeomError::InvalidInput(
            "a surface curve sweep needs one surface normal per directrix point".to_owned(),
        ));
    }
    let frames = frames_along(path, None, |i| normals[i])?;
    let stations: Vec<Station> = frames
        .iter()
        .map(|f| loft::place(rings, |p| loft::at(f, p)))
        .collect();
    loft::loft(rings, &stations, false)
}

/// Sweep a disk along a directrix, optionally hollow.
///
/// `fillet_radius` is refused rather than ignored. The model's own docs say
/// a consumer that cannot round corners must refuse a `Some`, because
/// silently sharpening a pipe run produces geometry that builds, renders,
/// and is wrong.
pub fn swept_disk(
    path: &[Point3],
    radius: Scalar,
    inner_radius: Option<Scalar>,
    fillet_radius: Option<Scalar>,
    tolerance: Tolerance,
) -> GeomResult<TriMesh> {
    check_disk(radius, inner_radius, fillet_radius)?;
    // A disk is just a circular profile, so the sweep reuses the shared
    // loft. The section needs a frame that stays perpendicular to the path,
    // but which perpendicular does not matter for a circle. A single fixed
    // axis is NOT enough: a leg along that axis has no perpendicular
    // component and was refused, and a leg nearly along it projects to a
    // residue of arbitrary direction, rotating the ring between stations
    // so the loft connects vertex k to a rotated vertex k and the volume
    // collapses silently (axiolid/kernel#169). The reference is therefore
    // carried along the path by rotation-minimising frames.
    let rings = disk_rings(radius, inner_radius, tolerance.linear())?;
    let frames = rotation_minimising_frames(path, None)?;
    let stations: Vec<Station> = frames
        .iter()
        .map(|f| loft::place(&rings, |p| loft::at(f, p)))
        .collect();
    loft::loft(&rings, &stations, false)
}

/// [`swept_disk`] whose surface stays within `chord` of the exact tube
/// round a smooth directrix (#231).
///
/// Half the budget chords the disk: a ring of `n` points inscribed in a
/// circle of radius `r` is within `r (1 - cos(pi/n))` of it, and every
/// station places it rigidly, so every point of the exact tube is within
/// that of the surface the rings sweep. The other half bounds that surface
/// against its triangles, as [`fixed_reference_sweep_within`] does, with
/// the disk's radius as the reach: along a circular arc of radius `R` the
/// rotation-minimising frames turn with the arc, so the stations are
/// rotations of one ring about the arc's axis and the outer side of the
/// tube, `R + r` from it, sets the step. `sample(b)` samples the directrix
/// to chord budget `b`; see [`SampledPath`].
pub fn swept_disk_within(
    mut sample: impl FnMut(Scalar) -> GeomResult<SampledPath>,
    radius: Scalar,
    inner_radius: Option<Scalar>,
    fillet_radius: Option<Scalar>,
    chord: Scalar,
) -> GeomResult<TriMesh> {
    check_disk(radius, inner_radius, fillet_radius)?;
    let rings = disk_rings(radius, inner_radius, 0.5 * chord)?;
    let stations = bounded_stations(&rings, 0.5 * chord, |b| {
        let path = sample(b)?;
        let frames = rotation_minimising_frames(&path.points, path.end_tangents)?;
        Ok((frames, path.end_tangents.is_some()))
    })?;
    loft::loft(&rings, &stations, false)
}

/// Refuse a fillet and validate the disk's radii.
///
/// `fillet_radius` is refused rather than ignored. The model's own docs say
/// a consumer that cannot round corners must refuse a `Some`, because
/// silently sharpening a pipe run produces geometry that builds, renders,
/// and is wrong.
fn check_disk(
    radius: Scalar,
    inner_radius: Option<Scalar>,
    fillet_radius: Option<Scalar>,
) -> GeomResult<()> {
    if fillet_radius.is_some() {
        return Err(GeomError::Unsupported {
            backend: crate::BACKEND_ID,
            operation: axiolid_contracts::Operation::Sweep,
        });
    }
    check_radii(radius, inner_radius)
}

/// Validate a swept disk's outer and inner radii.
pub(crate) fn check_radii(radius: Scalar, inner_radius: Option<Scalar>) -> GeomResult<()> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "swept disk radius must be positive and finite, got {radius}"
        )));
    }
    if let Some(inner) = inner_radius {
        if !inner.is_finite() || inner <= 0.0 || inner >= radius {
            return Err(GeomError::InvalidInput(format!(
                "swept disk inner radius must be positive and below {radius}, got {inner}"
            )));
        }
    }
    Ok(())
}

/// Rotation-minimising frames along a sampled path, by double reflection.
///
/// Wang, Jüttler, Zheng and Liu, "Computation of Rotation Minimizing
/// Frames", ACM TOG 27(1), 2008: each step reflects the previous frame in
/// the bisector plane of the chord, then in the plane that maps the
/// reflected tangent onto the next tangent. Two reflections are a rotation,
/// so the reference stays unit length and perpendicular to the tangent at
/// every station by construction; it can never become parallel to it. The
/// frame has fourth-order accuracy in the step, and no twist beyond what
/// the path's own torsion forces.
///
/// Tangents match `frames_along` (averaged at interior samples, the exact
/// ones at the ends when `ends` gives them), so a corner is mitred the
/// same way as every other sweep family.
fn rotation_minimising_frames(path: &[Point3], ends: Option<[Vec3; 2]>) -> GeomResult<Vec<Frame>> {
    let seed = seed_reference(path)?;
    let tangents = tangents_along(path, ends)?;
    frames_carried(path, &tangents, seed)
}

/// Rotation-minimising frames along `path` with the given unit tangent at
/// each sample, the first frame's x being `seed` made perpendicular to the
/// first tangent; see [`rotation_minimising_frames`].
pub(crate) fn frames_carried(
    path: &[Point3],
    tangents: &[Vec3],
    seed: Vec3,
) -> GeomResult<Vec<Frame>> {
    let mut reference = seed;
    let mut frames = Vec::with_capacity(path.len());
    for i in 0..path.len() {
        frames.push(Frame::from_reference(path[i], tangents[i], reference)?);
        let Some(next) = path.get(i + 1) else { break };
        let chord = *next - path[i];
        let c1 = chord.dot(chord);
        if c1 == 0.0 {
            // A repeated sample: the frame does not move. `tangents_along`
            // has already refused a path whose tangent vanishes here.
            continue;
        }
        let reflected_ref = reference - chord * (2.0 / c1 * chord.dot(reference));
        let reflected_tan = tangents[i] - chord * (2.0 / c1 * chord.dot(tangents[i]));
        let fix = tangents[i + 1] - reflected_tan;
        let c2 = fix.dot(fix);
        reference = if c2 == 0.0 {
            reflected_ref
        } else {
            reflected_ref - fix * (2.0 / c2 * fix.dot(reflected_ref))
        };
    }
    Ok(frames)
}

/// Unit tangent at each sample, averaged at interior samples; at the two
/// ends the exact tangents when `ends` gives them, else the end chords.
fn tangents_along(path: &[Point3], ends: Option<[Vec3; 2]>) -> GeomResult<Vec<Vec3>> {
    if path.len() < 2 {
        return Err(GeomError::InvalidInput(format!(
            "a sweep directrix needs at least two points, got {}",
            path.len()
        )));
    }
    (0..path.len())
        .map(|i| {
            let raw = if i == 0 {
                ends.map_or(path[1] - path[0], |[start, _]| start)
            } else if i + 1 == path.len() {
                ends.map_or(path[i] - path[i - 1], |[_, end]| end)
            } else {
                (path[i] - path[i - 1]).normalize_or_zero()
                    + (path[i + 1] - path[i]).normalize_or_zero()
            };
            let tangent = raw.normalize_or_zero();
            if tangent == Vec3::ZERO {
                Err(GeomError::InvalidInput(
                    "sweep tangent must be a non-zero direction".to_owned(),
                ))
            } else {
                Ok(tangent)
            }
        })
        .collect()
}

/// A circular profile, hollow when `inner` is given, each ring within
/// `chord` of its circle.
pub(crate) fn disk_rings(
    radius: Scalar,
    inner: Option<Scalar>,
    chord: Scalar,
) -> GeomResult<Rings> {
    let circle = |r: Scalar, reverse: bool| -> GeomResult<Vec<Point2>> {
        let n = crate::revolve::steps(r, core::f64::consts::TAU, chord)?;
        let mut pts: Vec<Point2> = (0..n)
            .map(|k| {
                let a = core::f64::consts::TAU * (k as Scalar) / (n as Scalar);
                Point2::new(r * a.cos(), r * a.sin())
            })
            .collect();
        if reverse {
            pts.reverse();
        }
        Ok(pts)
    };
    // A hole ring runs opposite the outer ring so the triangulator reads it
    // as a void rather than a second island.
    Ok(Rings {
        outer: circle(radius, false)?,
        holes: match inner {
            Some(r) => vec![circle(r, true)?],
            None => Vec::new(),
        },
    })
}

/// A reference direction guaranteed not to be parallel to the first
/// segment.
///
/// A circular section has no preferred orientation, so any perpendicular
/// will do; what matters is that it is never degenerate.
fn seed_reference(path: &[Point3]) -> GeomResult<Vec3> {
    if path.len() < 2 {
        return Err(GeomError::InvalidInput(
            "a swept disk directrix needs at least two points".to_owned(),
        ));
    }
    let t = (path[1] - path[0]).normalize_or_zero();
    if t == Vec3::ZERO {
        return Err(GeomError::InvalidInput(
            "a swept disk directrix must not start with a zero-length segment".to_owned(),
        ));
    }
    let candidate = if t.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    Ok(candidate - t * t.dot(candidate))
}

/// Loft explicit sections placed along a spine.
///
/// The caller has already resolved each section's profile and placement, so
/// this only stitches them. Sections must share a ring structure: a spine
/// whose sections differ in topology has no vertex correspondence, and
/// pairing by index would weld unrelated points.
pub fn sectioned_spine(sections: &[(Rings, Vec<Point3>)]) -> GeomResult<TriMesh> {
    if sections.len() < 2 {
        return Err(GeomError::InvalidInput(format!(
            "a sectioned spine needs at least two sections, got {}",
            sections.len()
        )));
    }
    let stations: Vec<Station> = sections
        .iter()
        .map(|(rings, placed)| {
            let mut it = placed.iter().copied();
            let mut loops = Vec::with_capacity(1 + rings.holes.len());
            loops.push((0..rings.outer.len()).filter_map(|_| it.next()).collect());
            for hole in &rings.holes {
                loops.push((0..hole.len()).filter_map(|_| it.next()).collect());
            }
            Station { loops }
        })
        .collect();
    loft::loft(&sections[0].0, &stations, false)
}
