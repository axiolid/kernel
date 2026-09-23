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

/// Sample one channel of one operand at one result corner.
///
/// Appends the tuple to `out` (when there is one) and says how it was got.
fn sample(
    source: &TriMesh,
    channel: &AttributeChannel,
    triangle: usize,
    p: Point3,
    out: &mut Vec<Scalar>,
) -> Sample {
    let base = triangle * 3;
    let Some(corners) = source.indices.get(base..base + 3) else {
        return Sample::Unmapped;
    };
    let at = |k: usize| channel.at_corner(&source.indices, base + k);
    // Channel values for all three corners, or none: a triangle is mapped
    // whole or not at all (validated in axiolid-mesh).
    let (Some(v0), Some(v1), Some(v2)) = (at(0), at(1), at(2)) else {
        return Sample::Unmapped;
    };
    let pos = |k: usize| source.positions[corners[k] as usize];
    let values = [v0, v1, v2];

    if let Some(k) = (0..3).find(|&k| pos(k) == p) {
        out.extend_from_slice(values[k]);
        return Sample::Copied;
    }

    let weights = barycentric(p, pos(0), pos(1), pos(2));
    match (channel.blend, weights) {
        (Blend::Linear, Some(w)) => {
            for d in 0..channel.width {
                out.push(w[0] * values[0][d] + w[1] * values[1][d] + w[2] * values[2][d]);
            }
        }
        // Nearest: the corner with the largest weight. With no usable
        // weights (a sliver source) the nearest corner by distance is the
        // same question asked without the triangle's help.
        (Blend::Linear | Blend::Nearest, w) => {
            let k = match w {
                Some(w) => (0..3).max_by(|&i, &j| w[i].total_cmp(&w[j])).unwrap_or(0),
                None => (0..3)
                    .min_by(|&i, &j| {
                        (pos(i) - p)
                            .length_squared()
                            .total_cmp(&(pos(j) - p).length_squared())
                    })
                    .unwrap_or(0),
            };
            out.extend_from_slice(values[k]);
        }
        // No value may be derived. The caller drops the whole channel, so
        // pushing nothing here is never observed.
        (Blend::None, _) => {}
    }
    Sample::Derived
}

/// Build the result's channels and report each one's fate.
///
/// `result` is the boolean's output, `sources[f]` the origin of its face
/// `f`. Channels are named by the SUBJECT; a tool channel of the same name,
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
        match build(subject, tool, channel, tool_channel, result, sources) {
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
    subject: &TriMesh,
    tool: &TriMesh,
    channel: &AttributeChannel,
    tool_channel: Option<&AttributeChannel>,
    result: &TriMesh,
    sources: &[FaceSource],
) -> Result<(AttributeChannel, bool), DropReason> {
    if sources.len() * 3 != result.indices.len() {
        // The provenance record does not describe this mesh. Refusing is
        // the only answer that cannot attach a value to the wrong face.
        return Err(DropReason::ProviderLimitation);
    }
    let mut values = Vec::with_capacity(result.indices.len() * channel.width);
    let mut corners = Vec::with_capacity(result.indices.len());
    let mut derived = false;
    for (corner, &vertex) in result.indices.iter().enumerate() {
        let source = sources[corner / 3];
        let (mesh, from) = match source.operand {
            0 => (subject, Some(channel)),
            _ => (tool, tool_channel),
        };
        let p = result.positions[vertex as usize];
        let before = values.len();
        let how = match from {
            Some(from) => sample(mesh, from, source.triangle, p, &mut values),
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
    // A face may meet a source whose corner is mapped at one vertex and not
    // another only if its sources disagree, and a face has one source. So a
    // triangle is unmapped whole or mapped whole, as validation requires.
    Ok((
        AttributeChannel::corner_indexed(
            channel.name.clone(),
            values,
            channel.width,
            channel.blend,
            corners,
        ),
        derived,
    ))
}
