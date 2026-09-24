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
    for p in start.outer.iter().chain(far.outer.iter()) {
        let v = Point3::new(p.x, p.y, 0.0) - axis_origin;
        max_r = max_r.max((v - dir * dir.dot(v)).length());
    }
    let n = crate::revolve::steps(max_r, angle, tolerance.linear());
    let mut stations = Vec::with_capacity(n + 1);
    for s in 0..=n {
        let t = (s as Scalar) / (n as Scalar);
        let ring = blend_rings(start, end, t)?;
        let a = angle * t;
        stations.push(loft::place(&ring, |p| {
            crate::revolve::rotate(Point3::new(p.x, p.y, 0.0), axis_origin, dir, a)
        }));
    }
    let stations: Vec<_> = stations.into_iter().rev().collect();
    loft::loft_tapered(&far, start, &stations)
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
    let frames = frames_along(path, |_| reference)?;
    let stations: Vec<Station> = frames
        .iter()
        .map(|f| loft::place(rings, |p| loft::at(f, p)))
        .collect();
    loft::loft(rings, &stations, false)
}

/// Build a frame at each path sample.
///
/// The tangent at an interior sample is the average of its two segment
/// directions, which keeps the profile from kinking at a corner. Endpoints
/// use their single adjacent segment.
fn frames_along(path: &[Point3], up: impl Fn(usize) -> Vec3) -> GeomResult<Vec<Frame>> {
    if path.len() < 2 {
        return Err(GeomError::InvalidInput(format!(
            "a sweep directrix needs at least two points, got {}",
            path.len()
        )));
    }
    let mut frames = Vec::with_capacity(path.len());
    for i in 0..path.len() {
        let tangent = if i == 0 {
            path[1] - path[0]
        } else if i + 1 == path.len() {
            path[i] - path[i - 1]
        } else {
            (path[i] - path[i - 1]).normalize_or_zero()
                + (path[i + 1] - path[i]).normalize_or_zero()
        };
        frames.push(Frame::from_reference(path[i], tangent, up(i))?);
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
    let frames = frames_along(path, |i| normals[i])?;
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
    if fillet_radius.is_some() {
        return Err(GeomError::Unsupported {
            backend: crate::BACKEND_ID,
            operation: axiolid_contracts::Operation::Sweep,
        });
    }
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
    // A disk is just a circular profile, so the sweep reuses the shared
    // loft. The section needs a frame that stays perpendicular to the path,
    // but which perpendicular does not matter for a circle. A single fixed
    // axis is NOT enough: a leg along that axis has no perpendicular
    // component and was refused, and a leg nearly along it projects to a
    // residue of arbitrary direction, rotating the ring between stations
    // so the loft connects vertex k to a rotated vertex k and the volume
    // collapses silently (axiolid/kernel#169). The reference is therefore
    // carried along the path by rotation-minimising frames.
    let rings = disk_rings(radius, inner_radius, tolerance)?;
    let frames = rotation_minimising_frames(path)?;
    let stations: Vec<Station> = frames
        .iter()
        .map(|f| loft::place(&rings, |p| loft::at(f, p)))
        .collect();
    loft::loft(&rings, &stations, false)
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
/// Tangents match `frames_along` (averaged at interior samples), so a
/// corner is mitred the same way as every other sweep family.
fn rotation_minimising_frames(path: &[Point3]) -> GeomResult<Vec<Frame>> {
    let seed = seed_reference(path)?;
    let tangents = tangents_along(path)?;
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

/// Unit tangent at each sample, averaged at interior samples.
fn tangents_along(path: &[Point3]) -> GeomResult<Vec<Vec3>> {
    if path.len() < 2 {
        return Err(GeomError::InvalidInput(format!(
            "a sweep directrix needs at least two points, got {}",
            path.len()
        )));
    }
    (0..path.len())
        .map(|i| {
            let raw = if i == 0 {
                path[1] - path[0]
            } else if i + 1 == path.len() {
                path[i] - path[i - 1]
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

/// A circular profile, hollow when `inner` is given.
fn disk_rings(radius: Scalar, inner: Option<Scalar>, tolerance: Tolerance) -> GeomResult<Rings> {
    let circle = |r: Scalar, reverse: bool| -> Vec<Point2> {
        let n = crate::revolve::steps(r, core::f64::consts::TAU, tolerance.linear());
        let mut pts: Vec<Point2> = (0..n)
            .map(|k| {
                let a = core::f64::consts::TAU * (k as Scalar) / (n as Scalar);
                Point2::new(r * a.cos(), r * a.sin())
            })
            .collect();
        if reverse {
            pts.reverse();
        }
        pts
    };
    // A hole ring runs opposite the outer ring so the triangulator reads it
    // as a void rather than a second island.
    Ok(Rings {
        outer: circle(radius, false),
        holes: inner.map(|r| vec![circle(r, true)]).unwrap_or_default(),
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
