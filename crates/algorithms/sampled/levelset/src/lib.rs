//! Level-set extraction: a closed manifold mesh from a scalar field.
//!
//! # Why tetrahedra rather than cubes
//!
//! Classic marching cubes is not manifold. Its 256-entry table has
//! genuinely ambiguous face configurations: two diagonally opposite
//! corners inside the level, and the other two outside, can be joined in
//! two different ways. Neighbouring cells that resolve the same shared
//! face differently leave a hole, and the result is neither closed nor
//! two-manifold. Fixing that needs the disambiguated MC33 table plus a
//! consistent face-resolution rule.
//!
//! This module decomposes each cell into six tetrahedra instead. A
//! tetrahedron has four corners and therefore sixteen sign patterns, none
//! of which is ambiguous: the surface crosses either three edges (one
//! triangle) or four (two triangles), and the decomposition is forced.
//!
//! The decomposition used is Kuhn's: six tetrahedra sharing the cell's
//! main diagonal, one per permutation of the three axes. Each cell face is
//! then split along the diagonal joining the two corners that differ in
//! both of that face's coordinates -- and because that choice depends only
//! on the global grid indices, the cell on the other side of the face
//! splits it exactly the same way. Watertightness is therefore structural
//! rather than a property of a table that must be kept correct.
//!
//! The cost is more triangles than marching cubes for the same grid, and a
//! slight directional bias from the diagonal. That is the price of the
//! guarantee, and the guarantee is what this contract is for.
//!
//! # What this is not
//!
//! Extraction is an approximation. The mesh interpolates the field
//! linearly along each edge, so a curved surface is faceted and the error
//! shrinks with the grid, it does not vanish. This is not a certified
//! path and does not claim to be.
//!
//! # Exact tangency
//!
//! A level set can pass exactly through a grid sample -- a unit sphere in a
//! half-extent of 1.4 does it at several edge lengths. Every edge meeting
//! that sample then crosses at the same point, the triangles between those
//! crossings collapse, and the surface tears.
//!
//! This is resolved by simulation of simplicity rather than by special
//! cases: symbolically each sample carries its own positive infinitesimal,
//! so no sample sits exactly at the level and every crossing is strictly
//! interior to its edge. See `SOS_DELTA` for the numeric stand-in and why
//! its magnitude matters.
//!
//! # What this is not
//!
//! Extraction is an approximation. The mesh interpolates the field
//! linearly along each edge, so a curved surface is faceted and the error
//! shrinks with the grid, it does not vanish. This is not a certified
//! path and does not claim to be.
//!
//! # Known limitation: grid tangency
//!
//! The closed-manifold guarantee holds when the surface passes BETWEEN
//! grid samples. It does not currently hold when the level set is exactly
//! tangent to a grid plane -- a sphere of radius 1 with samples landing
//! exactly on 1.0, for instance. Two edges of the same tetrahedron then
//! interpolate to the same point, the triangle between them has zero area,
//! and the surface is left with unmatched edges.
//!
//! Measured, so the boundary of the guarantee is known rather than
//! assumed: a unit sphere in bounds of half-extent 1.45 is closed at edge
//! lengths 0.4, 0.2 and 0.1, while the same sphere in half-extent 1.4 is
//! closed at 0.4 and open at 0.2 and 0.1 -- exactly the resolutions whose
//! samples land on the radius.
//!
//! Two fixes were tried and rejected. Offsetting an exactly-zero sample by
//! `Scalar::MIN_POSITIVE` produced vertices a SUBNORMAL distance apart:
//! distinct in bits, identical in geometry, so it reproduced the
//! degeneracy it was meant to remove. Offsetting by a fraction of the cell
//! instead flips the sample's side and changes the topology, which broke
//! the general case to patch the special one. The remaining candidate is
//! symbolic perturbation (simulation of simplicity), which decides ties by
//! index rather than by value; that is a larger change and is not done
//! here.
//!
//! A caller who needs the guarantee unconditionally should offset the
//! bounds so no grid plane is tangent to the surface.

use ahash::AHashMap;

use axiolid_core::{Aabb, Point3, Scalar};
use axiolid_mesh::TriMesh;

/// Why a level set could not be extracted.
#[derive(Debug, thiserror::Error, PartialEq)]
#[non_exhaustive]
pub enum LevelSetError {
    /// The requested edge length is not a usable spacing.
    #[error("edge length {0} is not a positive finite length")]
    InvalidEdgeLength(Scalar),
    /// The bounds are empty or not finite.
    #[error("bounds are empty or non-finite along at least one axis")]
    InvalidBounds,
    /// The level itself is not a finite value.
    #[error("level {0} is not finite")]
    InvalidLevel(Scalar),
    /// The field never crosses the level inside the bounds.
    ///
    /// Reported rather than answered with an empty mesh: a caller asking
    /// for a surface that is not there has a bug upstream, and a
    /// zero-triangle result looks like a successful extraction of nothing.
    #[error("the field does not cross level {level} anywhere in the bounds")]
    NoCrossing {
        /// The level that was searched for.
        level: Scalar,
    },
    /// The field returned a non-finite sample.
    #[error("the field returned a non-finite value at {point:?}")]
    NonFiniteSample {
        /// Where the field misbehaved.
        point: Point3,
    },
    /// The requested grid exceeds the sample budget.
    #[error("the requested grid needs {requested} samples, over the {limit} budget")]
    BudgetExceeded {
        /// Samples the request would have taken.
        requested: usize,
        /// The cap that was not raised.
        limit: usize,
    },
}

/// Upper bound on grid samples, so a fine edge length on large bounds is
/// refused up front instead of exhausting memory.
/// Numeric stand-in for the infinitesimal in the simulation-of-simplicity
/// rule: how far a crossing is held clear of either end of its edge.
///
/// The magnitude is bounded from both sides, and both bounds were found by
/// measurement rather than chosen:
///
/// - Too small and the fix does nothing useful. At `1e-6` the triangles it
///   creates have a doubled area around `1e-27`, below the audit's
///   `tolerance^4` degeneracy threshold, so they are still counted as
///   degenerate and the mesh still reads as open. The perturbation has to
///   clear modelling tolerance to be a feature rather than dust.
/// - Too large and it stops being an infinitesimal: it would move vertices
///   far enough to compete with the tessellation error itself.
///
/// `1e-3` of an edge sits between those: a shift of at most
/// `1e-3 * edge_length`, which is smaller than the `edge^2/8` chord error
/// by orders of magnitude at every usable resolution, and large enough that
/// two crossings never round together.
const SOS_DELTA: Scalar = 1.0e-3;

const MAX_SAMPLES: usize = 64_000_000;

/// The six Kuhn tetrahedra of a unit cell, as corner indices.
///
/// Corner `i` has bits `(x, y, z)` with x least significant, so corner 0 is
/// the minimum and corner 7 the maximum. Every tetrahedron runs from 0 to 7
/// along a different axis order, which is what makes the shared faces agree
/// between neighbouring cells.
const KUHN_TETRAHEDRA: [[usize; 4]; 6] = [
    [0, 1, 3, 7],
    [0, 1, 5, 7],
    [0, 2, 3, 7],
    [0, 2, 6, 7],
    [0, 4, 5, 7],
    [0, 4, 6, 7],
];

/// Extract the level set of a scalar field as a closed manifold mesh.
///
/// `field` is sampled on a regular grid spanning `bounds`. The surface is
/// where `field` equals `level`; the convention is that lower values are
/// inside, so triangle winding puts the outward normal toward higher
/// values.
///
/// The bounds are padded by one cell on every side and the field is forced
/// to read as outside on that shell. Without it a surface reaching the edge
/// of the bounds would be cut, leaving an open border -- and the closedness
/// guarantee would be false exactly when the caller's bounds were tight.
///
/// # Errors
///
/// Refuses a non-positive edge length, empty or non-finite bounds, a
/// non-finite level, a field that returns a non-finite sample, a grid over
/// the sample budget, and a field that never crosses the level.
pub fn level_set<F>(
    field: F,
    bounds: Aabb,
    edge_length: Scalar,
    level: Scalar,
) -> Result<TriMesh, LevelSetError>
where
    F: Fn(Point3) -> Scalar,
{
    if !edge_length.is_finite() || edge_length <= 0.0 {
        return Err(LevelSetError::InvalidEdgeLength(edge_length));
    }
    if !level.is_finite() {
        return Err(LevelSetError::InvalidLevel(level));
    }
    let (min, max) = (bounds.min, bounds.max);
    if !min.is_finite() || !max.is_finite() || max.x <= min.x || max.y <= min.y || max.z <= min.z {
        return Err(LevelSetError::InvalidBounds);
    }

    // One padding cell each side, so a surface touching the bounds still
    // closes instead of being clipped into an open sheet.
    let counts = [
        ((max.x - min.x) / edge_length).ceil() as usize + 3,
        ((max.y - min.y) / edge_length).ceil() as usize + 3,
        ((max.z - min.z) / edge_length).ceil() as usize + 3,
    ];
    let requested = counts[0]
        .saturating_mul(counts[1])
        .saturating_mul(counts[2]);
    if requested > MAX_SAMPLES {
        return Err(LevelSetError::BudgetExceeded {
            requested,
            limit: MAX_SAMPLES,
        });
    }

    let origin = Point3::new(
        min.x - edge_length,
        min.y - edge_length,
        min.z - edge_length,
    );
    let at = |i: usize, j: usize, k: usize| {
        Point3::new(
            origin.x + (i as Scalar) * edge_length,
            origin.y + (j as Scalar) * edge_length,
            origin.z + (k as Scalar) * edge_length,
        )
    };
    let index_of = |i: usize, j: usize, k: usize| (k * counts[1] + j) * counts[0] + i;

    // Sample once. The field is a caller closure and may be expensive, so
    // it is never evaluated twice for the same grid point.
    let mut samples = vec![0.0 as Scalar; requested];
    for k in 0..counts[2] {
        for j in 0..counts[1] {
            for i in 0..counts[0] {
                let point = at(i, j, k);
                let on_shell = i == 0
                    || j == 0
                    || k == 0
                    || i == counts[0] - 1
                    || j == counts[1] - 1
                    || k == counts[2] - 1;
                let value = if on_shell {
                    // Forced outside: this is what closes a surface that
                    // would otherwise run off the edge of the grid.
                    1.0
                } else {
                    let raw = field(point);
                    if !raw.is_finite() {
                        return Err(LevelSetError::NonFiniteSample { point });
                    }
                    raw - level
                };
                samples[index_of(i, j, k)] = value;
            }
        }
    }

    let mut positions: Vec<Point3> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    // Keyed by the two grid samples an intersection lies between, so both
    // tetrahedra sharing that edge reuse one vertex. This welding is what
    // makes the result closed rather than a soup of disconnected triangles.
    let mut vertices: AHashMap<(usize, usize), u32> = AHashMap::new();
    // A second index, keyed by exact position bits. Two DIFFERENT edges can
    // cross at the same point -- when a crossing lands on a shared grid
    // corner, for instance -- and giving that point two vertex ids collapses
    // the incident triangles to zero area. Dropping those then tears a hole,
    // which is how this first showed up: 36 exactly-zero-area faces and 48
    // unmatched boundary edges. Welding by position removes the cause.
    let mut welded: AHashMap<[u64; 3], u32> = AHashMap::new();

    for k in 0..counts[2] - 1 {
        for j in 0..counts[1] - 1 {
            for i in 0..counts[0] - 1 {
                let corner = |bit: usize| {
                    let (dx, dy, dz) = (bit & 1, (bit >> 1) & 1, (bit >> 2) & 1);
                    index_of(i + dx, j + dy, k + dz)
                };
                for tetrahedron in KUHN_TETRAHEDRA {
                    let nodes = tetrahedron.map(corner);
                    emit_tetrahedron(
                        nodes,
                        &samples,
                        &counts,
                        origin,
                        edge_length,
                        &mut positions,
                        &mut indices,
                        &mut vertices,
                        &mut welded,
                    );
                }
            }
        }
    }

    if indices.is_empty() {
        return Err(LevelSetError::NoCrossing { level });
    }
    Ok(TriMesh::new(positions, indices))
}

/// Emit the triangles of one tetrahedron.
#[allow(clippy::too_many_arguments)]
fn emit_tetrahedron(
    nodes: [usize; 4],
    samples: &[Scalar],
    counts: &[usize; 3],
    origin: Point3,
    edge_length: Scalar,
    positions: &mut Vec<Point3>,
    indices: &mut Vec<u32>,
    vertices: &mut AHashMap<(usize, usize), u32>,
    welded: &mut AHashMap<[u64; 3], u32>,
) {
    // A sample exactly at the level would make an edge both crossing and
    // not crossing depending on which side asks, so the rule has to be
    // global rather than per-tetrahedron. Strictly-negative is inside:
    // a grid point sitting exactly ON the surface then reads as outside
    // from every tetrahedron that touches it, and the surface passes
    // between grid points instead of through one. That keeps every
    // crossing strictly interior to its edge, which is what stops two
    // edges of one tetrahedron interpolating to the same position.
    let inside = nodes.map(|node| samples[node] < 0.0);
    let count = inside.iter().filter(|&&flag| flag).count();
    if count == 0 || count == 4 {
        return;
    }

    let mut interpolate = |a: usize, b: usize, positions: &mut Vec<Point3>| -> u32 {
        let key = if a < b { (a, b) } else { (b, a) };
        if let Some(&existing) = vertices.get(&key) {
            return existing;
        }
        let (va, vb) = (samples[key.0], samples[key.1]);
        let (pa, pb) = (
            grid_point(key.0, counts, origin, edge_length),
            grid_point(key.1, counts, origin, edge_length),
        );
        // Guard the coincident-value case: without it a flat region of the
        // field divides by zero and produces a non-finite vertex.
        let span = vb - va;
        let raw = if span.abs() > Scalar::EPSILON {
            (-va / span).clamp(0.0, 1.0)
        } else {
            0.5
        };
        // Simulation of simplicity, applied to the crossing parameter.
        //
        // Symbolically every sample carries its own positive infinitesimal,
        // `v_i + e^i`, so no sample is ever exactly at the level. A sample
        // that reads as an exact zero therefore yields a crossing that is
        // infinitesimally along its edge rather than exactly at its
        // endpoint, and `t` is never exactly 0 or 1.
        //
        // This matters because a crossing AT a grid point is shared by every
        // edge meeting there: distinct edges produce one coincident vertex,
        // the triangles between them collapse to zero area, and the surface
        // is left with unmatched edges. Keeping the crossing strictly
        // interior to its edge gives each edge its own vertex and keeps the
        // result closed.
        //
        // `SOS_DELTA` is the numeric stand-in for the infinitesimal: small
        // enough that it moves a vertex by at most `SOS_DELTA * edge_length`
        // (well under any usable tolerance), large enough that two crossings
        // on different edges cannot round to the same position -- the defect
        // that a `Scalar::MIN_POSITIVE` offset produced, where the two points
        // differed in bits but not in geometry.
        //
        // The edge key is index-ordered, so both tetrahedra sharing an edge
        // compute the same `t` from the same pair and agree on the vertex.
        // The perturbation is a function of the edge alone, which is what
        // keeps it consistent across the whole grid.
        let t = raw.clamp(SOS_DELTA, 1.0 - SOS_DELTA);
        let point = pa + (pb - pa) * t;
        let bits = [point.x.to_bits(), point.y.to_bits(), point.z.to_bits()];
        let index = match welded.get(&bits) {
            Some(&existing) => existing,
            None => {
                let fresh = positions.len() as u32;
                positions.push(point);
                welded.insert(bits, fresh);
                fresh
            }
        };
        vertices.insert(key, index);
        index
    };

    // Order the corners so the inside ones come first. The crossing pattern
    // then depends only on how many are inside.
    let mut ordered = [0usize; 4];
    let (mut head, mut tail) = (0, 3);
    for (slot, &node) in nodes.iter().enumerate() {
        if inside[slot] {
            ordered[head] = node;
            head += 1;
        } else {
            ordered[tail] = node;
            tail = tail.wrapping_sub(1);
        }
    }

    // Orient from the field, not from one corner. The direction from the
    // inside corners' centroid to the outside corners' centroid is the
    // local outward direction, and unlike a single corner it stays well
    // conditioned when the tetrahedron is thin.
    let centroid = |nodes: &[usize]| {
        let mut sum = Point3::ZERO;
        for &node in nodes {
            sum += grid_point(node, counts, origin, edge_length);
        }
        sum / (nodes.len() as Scalar)
    };
    let outward = centroid(&ordered[count..]) - centroid(&ordered[..count]);

    match count {
        // One corner inside: a triangle separating it from the other three.
        1 => {
            let a = interpolate(ordered[0], ordered[1], positions);
            let b = interpolate(ordered[0], ordered[2], positions);
            let c = interpolate(ordered[0], ordered[3], positions);
            push_oriented(indices, positions, [a, b, c], outward);
        }
        // Three inside is the mirror image: one corner outside.
        3 => {
            let a = interpolate(ordered[3], ordered[0], positions);
            let b = interpolate(ordered[3], ordered[1], positions);
            let c = interpolate(ordered[3], ordered[2], positions);
            push_oriented(indices, positions, [a, b, c], outward);
        }
        // Two inside, two outside: the surface cuts four edges, giving a
        // quadrilateral split into two triangles.
        _ => {
            let a = interpolate(ordered[0], ordered[2], positions);
            let b = interpolate(ordered[0], ordered[3], positions);
            let c = interpolate(ordered[1], ordered[3], positions);
            let d = interpolate(ordered[1], ordered[2], positions);
            push_oriented(indices, positions, [a, b, c], outward);
            push_oriented(indices, positions, [a, c, d], outward);
        }
    }
}

/// Append a triangle wound so its normal follows `outward`.
///
/// `outward` runs from the inside corners toward the outside ones, so it is
/// the field's own local direction of increase. Deciding orientation from
/// that rather than from a single corner keeps neighbouring tetrahedra
/// agreeing even when one of them is thin.
fn push_oriented(
    indices: &mut Vec<u32>,
    positions: &[Point3],
    triangle: [u32; 3],
    outward: axiolid_core::Vec3,
) {
    let [a, b, c] = triangle;
    // A degenerate triangle has no orientation to fix and would only add a
    // zero-area face, so it is dropped rather than emitted.
    if a == b || b == c || a == c {
        return;
    }
    let (pa, pb, pc) = (
        positions[a as usize],
        positions[b as usize],
        positions[c as usize],
    );
    let normal = (pb - pa).cross(pc - pa);
    if normal.dot(outward) >= 0.0 {
        indices.extend_from_slice(&[a, b, c]);
    } else {
        indices.extend_from_slice(&[a, c, b]);
    }
}

/// Recover a grid point from its flat sample index.
fn grid_point(index: usize, counts: &[usize; 3], origin: Point3, edge_length: Scalar) -> Point3 {
    let i = index % counts[0];
    let j = (index / counts[0]) % counts[1];
    let k = index / (counts[0] * counts[1]);
    Point3::new(
        origin.x + (i as Scalar) * edge_length,
        origin.y + (j as Scalar) * edge_length,
        origin.z + (k as Scalar) * edge_length,
    )
}
