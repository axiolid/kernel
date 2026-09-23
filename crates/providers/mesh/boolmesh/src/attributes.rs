//! Carry attribute channels through a boolean (#116).
//!
//! Every result triangle lies inside ONE operand triangle: the boolean only
//! cuts faces, it never creates surface. The CSG core records which one
//! (`Tref`, lifted to the caller's numbering by `Manifold::face_src`), so a
//! channel's value at any result corner is the source triangle's value there,
//! read the same way whatever the channel's addressing.
//!
//! - A result corner that IS a source corner (bit-identical position) copies
//!   that corner's value. Nothing is derived; this is what makes a seam or a
//!   hard UV boundary survive exactly.
//! - Any other corner lies on a cut. Its value is derived in the source
//!   triangle under the channel's own [`Blend`]: barycentric for `Linear`,
//!   the dominant corner for `Nearest`, nothing for `None`.
//!
//! The output channel is corner-indexed: a cut vertex is shared by faces from
//! both operands, and by faces whose sources carry different values, so one
//! value per position would be wrong on at least one side.

use axiolid_core::{Point3, Scalar};
use axiolid_mesh::{AttributeChannel, AttributeFate, Blend, DropReason, TriMesh};
use std::collections::{HashMap, HashSet, VecDeque};

/// Which operand, and which of its triangles, one result face came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct FaceSource {
    /// `0` for the subject, `1` for the tool.
    pub operand: u8,
    /// Triangle index in that operand's `TriMesh`.
    pub triangle: usize,
}

/// How one result corner got its value.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Sample {
    /// Copied from a source corner.
    Copied,
    /// Derived inside the source triangle.
    Derived,
    /// No value: the operand has no such channel, or its triangle is unmapped.
    Unmapped,
}

/// Barycentric weights of `p` in triangle `abc`, projected onto the plane
/// where the triangle is largest.
///
/// A result corner lies on its source face up to the boolean's own rounding,
/// so the projection is the exact answer for an exact point and a stable one
/// for a rounded point. Degenerate triangles return `None`; the caller then
/// falls back to the nearest corner rather than dividing by zero.
fn barycentric(p: Point3, a: Point3, b: Point3, c: Point3) -> Option<[Scalar; 3]> {
    let n = (b - a).cross(c - a);
    let (ax, ay, az) = (n.x.abs(), n.y.abs(), n.z.abs());
    // Drop the axis the normal points along: that projection has the
    // largest area and so the best-conditioned solve.
    let proj = |q: Point3| {
        if ax >= ay && ax >= az {
            (q.y, q.z)
        } else if ay >= az {
            (q.z, q.x)
        } else {
            (q.x, q.y)
        }
    };
    let (p, a, b, c) = (proj(p), proj(a), proj(b), proj(c));
    let det = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
    if det == 0.0 || !det.is_finite() {
        return None;
    }
    let wb = ((p.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (p.1 - a.1)) / det;
    let wc = ((b.0 - a.0) * (p.1 - a.1) - (p.0 - a.0) * (b.1 - a.1)) / det;
    Some([1.0 - wb - wc, wb, wc])
}

/// Where a result corner's value is read from: an operand triangle, and the
/// barycentric weights of the corner in it.
struct Located {
    triangle: usize,
    weights: [Scalar; 3],
}

/// Weights within this tolerance of zero count as inside: a corner on a
/// shared edge is inside both triangles, and either gives the same value
/// for continuous data.
const INSIDE: Scalar = 1e-9;

/// Edge-adjacency of one operand, keyed by POSITION so it survives the
/// operand's own vertex duplication at seams: `(min, max)` position bits of
/// an edge to the triangles using it. Built lazily -- only a corner that
/// escapes its recorded triangle pays for it.
struct Adjacency {
    edges: HashMap<(PosKey, PosKey), Vec<usize>>,
}

type PosKey = [u64; 3];

fn pos_key(p: Point3) -> PosKey {
    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
}

impl Adjacency {
    fn new(mesh: &TriMesh) -> Self {
        let mut edges: HashMap<(PosKey, PosKey), Vec<usize>> = HashMap::new();
        for (t, tri) in mesh.indices.chunks_exact(3).enumerate() {
            for k in 0..3 {
                let a = pos_key(mesh.positions[tri[k] as usize]);
                let b = pos_key(mesh.positions[tri[(k + 1) % 3] as usize]);
                edges.entry((a.min(b), a.max(b))).or_default().push(t);
            }
        }
        Self { edges }
    }

    fn neighbours<'a>(&'a self, mesh: &'a TriMesh, t: usize) -> impl Iterator<Item = usize> + 'a {
        let tri = &mesh.indices[t * 3..t * 3 + 3];
        (0..3).flat_map(move |k| {
            let a = pos_key(mesh.positions[tri[k] as usize]);
            let b = pos_key(mesh.positions[tri[(k + 1) % 3] as usize]);
            self.edges
                .get(&(a.min(b), a.max(b)))
                .into_iter()
                .flatten()
                .copied()
        })
    }
}

/// Whether `a` and `b` of one operand lie in the same plane (unit normals
/// agree). Flood only crosses such edges: the result face is planar, so its
/// true source is coplanar with the recorded one.
fn coplanar(mesh: &TriMesh, a: usize, b: usize) -> bool {
    let n = |t: usize| {
        let q = |k: usize| mesh.positions[mesh.indices[t * 3 + k] as usize];
        (q(1) - q(0)).cross(q(2) - q(0)).normalize_or_zero()
    };
    n(a).dot(n(b)) > 1.0 - 1e-9
}

fn weights_in(mesh: &TriMesh, t: usize, p: Point3) -> Option<[Scalar; 3]> {
    let q = |k: usize| mesh.positions[mesh.indices[t * 3 + k] as usize];
    barycentric(p, q(0), q(1), q(2))
}

fn inside(w: &[Scalar; 3]) -> bool {
    w.iter().all(|&x| x >= -INSIDE)
}

/// Find the operand triangle a result corner lies in.
///
/// `p` is the corner, `probe` a point just inside the result face from it.
/// The recorded triangle wins when it contains `probe`. Otherwise the
/// coplanar region around it is searched breadth-first for the triangle
/// containing `probe` -- probe, not corner, so a corner on a seam picks the
/// side this face is on. Weights are then taken at `p` itself in that
/// triangle, where they are inside up to the probe's step.
///
/// Falls back to the recorded triangle when nothing contains the probe
/// (a sliver, or rounding at a region boundary): a near value there is
/// better than none, and the region search already rejected every better.
fn locate(
    mesh: &TriMesh,
    adjacency: &mut Option<Adjacency>,
    recorded: usize,
    p: Point3,
    probe: Point3,
) -> Option<Located> {
    let hit = |t: usize| {
        let w = weights_in(mesh, t, probe)?;
        inside(&w).then(|| weights_in(mesh, t, p)).flatten()
    };
    if let Some(weights) = hit(recorded) {
        return Some(Located {
            triangle: recorded,
            weights,
        });
    }
    let adjacency = adjacency.get_or_insert_with(|| Adjacency::new(mesh));
    let mut seen = HashSet::from([recorded]);
    let mut queue = VecDeque::from([recorded]);
    while let Some(t) = queue.pop_front() {
        for n in adjacency.neighbours(mesh, t) {
            if !seen.insert(n) || !coplanar(mesh, recorded, n) {
                continue;
            }
            if let Some(weights) = hit(n) {
                return Some(Located {
                    triangle: n,
                    weights,
                });
            }
            queue.push_back(n);
        }
    }
    weights_in(mesh, recorded, p).map(|weights| Located {
        triangle: recorded,
        weights,
    })
}

/// Sample one channel of one operand at one result corner.
///
/// Appends the tuple to `out` (when there is one) and says how it was got.
/// `at` is where `locate` put the corner.
fn sample(
    source: &TriMesh,
    channel: &AttributeChannel,
    at: &Located,
    p: Point3,
    out: &mut Vec<Scalar>,
) -> Sample {
    let base = at.triangle * 3;
    let Some(corners) = source.indices.get(base..base + 3) else {
        return Sample::Unmapped;
    };
    let value = |k: usize| channel.at_corner(&source.indices, base + k);
    // A triangle is mapped whole or not at all (validated in axiolid-mesh).
    let (Some(v0), Some(v1), Some(v2)) = (value(0), value(1), value(2)) else {
        return Sample::Unmapped;
    };
    let pos = |k: usize| source.positions[corners[k] as usize];
    let values = [v0, v1, v2];

    if let Some(k) = (0..3).find(|&k| pos(k) == p) {
        out.extend_from_slice(values[k]);
        return Sample::Copied;
    }
    let w = at.weights;
    match channel.blend {
        Blend::Linear => {
            for d in 0..channel.width {
                out.push(w[0] * values[0][d] + w[1] * values[1][d] + w[2] * values[2][d]);
            }
        }
        Blend::Nearest => {
            let k = (0..3).max_by(|&i, &j| w[i].total_cmp(&w[j])).unwrap_or(0);
            out.extend_from_slice(values[k]);
        }
        // Nothing may be derived; the caller drops the whole channel.
        Blend::None => {}
    }
    Sample::Derived
}

/// Build the result's channels and report each one's fate.
///
/// `result` is a boolean output, `sources[f]` the origin of its face `f`.
/// Channels are named by the SUBJECT; a tool channel of the same name,
/// width and blend supplies the tool's faces, and tool faces are unmapped
/// otherwise. Channels only the tool carries are not invented onto the
/// result -- the result is the subject, cut.
pub(crate) fn carry(
    subject: &TriMesh,
    tool: &TriMesh,
    result: &mut TriMesh,
    sources: &[FaceSource],
) -> Vec<(String, AttributeFate)> {
    let mut fates = Vec::with_capacity(subject.attributes.len());
    for channel in &subject.attributes {
        let tool_channel = tool.attributes.iter().find(|c| {
            c.name == channel.name && c.width == channel.width && c.blend == channel.blend
        });
        let operands = [(subject, Some(channel)), (tool, tool_channel)];
        match build(operands, channel, result, sources) {
            Ok((built, derived)) => {
                result.attributes.push(built);
                let fate = if derived {
                    AttributeFate::Interpolated
                } else {
                    AttributeFate::Preserved
                };
                fates.push((channel.name.clone(), fate));
            }
            Err(reason) => fates.push((channel.name.clone(), AttributeFate::Dropped(reason))),
        }
    }
    fates
}

/// One channel over the whole result. `Ok((channel, any_corner_derived))`.
fn build(
    operands: [(&TriMesh, Option<&AttributeChannel>); 2],
    channel: &AttributeChannel,
    result: &TriMesh,
    sources: &[FaceSource],
) -> Result<(AttributeChannel, bool), DropReason> {
    if sources.len() * 3 != result.indices.len() {
        // The provenance record does not describe this mesh. Refusing is
        // the only answer that cannot attach a value to the wrong face.
        return Err(DropReason::ProviderLimitation);
    }
    let mut adjacency: [Option<Adjacency>; 2] = [None, None];
    let mut values = Vec::with_capacity(result.indices.len() * channel.width);
    let mut corners = Vec::with_capacity(result.indices.len());
    let mut derived = false;
    for (face, source) in sources.iter().enumerate() {
        let side = usize::from(source.operand != 0);
        let (mesh, from) = operands[side];
        let tri = &result.indices[face * 3..face * 3 + 3];
        let q = |k: usize| result.positions[tri[k] as usize];
        let centroid = (q(0) + q(1) + q(2)) / 3.0;
        for k in 0..3 {
            let p = q(k);
            // Just inside this face from its corner: decides which side of a
            // source seam the corner belongs to for THIS face.
            let probe = p + (centroid - p) * PROBE_STEP;
            let before = values.len();
            let how = match from {
                Some(from) => match locate(mesh, &mut adjacency[side], source.triangle, p, probe) {
                    Some(at) => sample(mesh, from, &at, p, &mut values),
                    None => Sample::Unmapped,
                },
                None => Sample::Unmapped,
            };
            match how {
                Sample::Unmapped => corners.push(AttributeChannel::UNMAPPED),
                Sample::Derived if channel.blend == Blend::None => {
                    return Err(DropReason::NotBlendable);
                }
                Sample::Copied | Sample::Derived => {
                    derived |= how == Sample::Derived;
                    corners.push((before / channel.width) as u32);
                }
            }
        }
    }
    unmap_partial(&mut corners);
    let built = AttributeChannel::corner_indexed(
        channel.name.clone(),
        values,
        channel.width,
        channel.blend,
        corners,
    );
    Ok((built, derived))
}

/// Make every triangle mapped whole or not at all, as validation requires.
///
/// Corners are located independently, so one corner of a face can land in
/// an unmapped source triangle while its neighbours land in a mapped one
/// (the face straddles a partially mapped source region). Such a face has
/// no value the source defines over its whole extent: it is unmapped.
fn unmap_partial(corners: &mut [u32]) {
    for tri in corners.chunks_exact_mut(3) {
        if tri.contains(&AttributeChannel::UNMAPPED) {
            tri.fill(AttributeChannel::UNMAPPED);
        }
    }
}

/// Fraction of the way from a corner to its face's centroid that the probe
/// point moves. Small enough to stay in the triangle the corner touches,
/// large enough to clear f64 rounding at building coordinates.
const PROBE_STEP: Scalar = 1e-6;

/// Provenance for a result that was built rather than cut, by plane lookup.
///
/// For a result whose every face lies in some operand face and faces the
/// same way -- the analytic box path, which rebuilds the surface from a
/// grid -- each face's source is the operand triangle containing its
/// centroid among those whose plane and orientation match. `operands[i]`
/// carries a sign: `+1.0` if its faces appear as is (the subject of a
/// difference), `-1.0` if reversed (a cutter's faces line the hole).
///
/// `None` if any face has no such source: the result is then not what
/// this function assumes, and the caller must not attach values.
pub(crate) fn sources_by_plane(
    result: &TriMesh,
    operands: &[(&TriMesh, Scalar)],
) -> Option<Vec<FaceSource>> {
    let normal = |m: &TriMesh, t: usize| {
        let q = |k: usize| m.positions[m.indices[t * 3 + k] as usize];
        (q(1) - q(0)).cross(q(2) - q(0)).normalize_or_zero()
    };
    let mut out = Vec::with_capacity(result.triangle_count());
    // Plane-distance tolerance, relative to the result's size.
    let d = result.bounds().diagonal();
    let scale = d.x.max(d.y).max(d.z).max(Scalar::MIN_POSITIVE);
    for f in 0..result.triangle_count() {
        let q = |k: usize| result.positions[result.indices[f * 3 + k] as usize];
        let centroid = (q(0) + q(1) + q(2)) / 3.0;
        let n = normal(result, f);
        let found = operands.iter().enumerate().find_map(|(i, (m, sign))| {
            (0..m.triangle_count()).find_map(|t| {
                let facing = normal(m, t).dot(n) * sign > 1.0 - 1e-9;
                let a = m.positions[m.indices[t * 3] as usize];
                let on_plane = (centroid - a).dot(normal(m, t)).abs() <= scale * 1e-9;
                let w = (facing && on_plane)
                    .then(|| weights_in(m, t, centroid))
                    .flatten()?;
                inside(&w).then_some(FaceSource {
                    operand: i as u8,
                    triangle: t,
                })
            })
        })?;
        out.push(found);
    }
    Some(out)
}
