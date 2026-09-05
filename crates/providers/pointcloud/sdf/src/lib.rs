#![forbid(unsafe_code)]
//! Reference reconstruction: a signed-distance field from samples, extracted
//! as a level set.
//!
//! # Why this method
//!
//! The kernel already owns both halves of this: [`PointIndex`] answers
//! nearest-neighbour queries, and `axiolid-levelset` extracts a closed
//! manifold surface from any scalar field. Composing them gives a
//! reconstruction with **no external dependency** and no new geometric
//! machinery to verify.
//!
//! That matters more than raw quality. Poisson reconstruction produces a
//! smoother surface, but adopting it would mean vendoring a solver whose
//! numerics we cannot audit, to satisfy a contract whose whole purpose is
//! swappability. This provider exists so the contract is *verifiable* — a
//! better one can replace it without any consumer noticing.
//!
//! # How it works
//!
//! For a query point `p`, find the nearest samples and estimate the signed
//! distance to the surface they lie on:
//!
//! - **With normals**, project onto the neighbour's tangent plane. The sign
//!   is which side of that plane `p` falls on, so the surface passes exactly
//!   through the samples.
//! - **Without normals**, use unsigned distance offset by the sample
//!   spacing. This produces a surface *around* the points rather than
//!   through them, which is honest: with no orientation information there is
//!   no way to say which side is inside.
//!
//! The distinction is reported in the evidence, never hidden: a
//! positions-only reconstruction is a genuinely weaker result.
//!
//! # What it is not
//!
//! Not a hole filler. Where the capture has no data the field is
//! extrapolated from distant samples, and those triangles are counted in
//! `interpolated_triangles` so a caller can see how much of the surface is
//! inference rather than measurement.

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, Determinism, ExecutionOptions,
    ExecutionTarget, GeomResult, ScratchRequirement,
};
use axiolid_core::{Aabb, Point3, Scalar, Vec3};
use axiolid_levelset::level_set;
use axiolid_mesh::audit_mesh;
use axiolid_pointcloud::PointCloud;
use axiolid_pointcloud_reconstruction_contract::{
    PointcloudReconstruction, Reconstruction, ReconstructionEvidence, ReconstructionOutcome,
    ReconstructionRefusal, ReconstructionRequest, Resolution,
};
use axiolid_spatial::PointIndex;

/// How many neighbours the field estimate blends.
///
/// One neighbour makes the field a Voronoi surface with visible facets;
/// blending a few smooths it without washing out real detail. Six is enough
/// to span a local neighbourhood on a roughly uniform capture.
const BLEND_NEIGHBOURS: usize = 6;

/// Reference reconstruction provider.
#[derive(Debug, Default, Clone, Copy)]
pub struct SdfReconstruction;

impl SdfReconstruction {
    /// Backend identity.
    pub const ID: BackendId = BackendId::new("axiolid.pointcloud.sdf");

    /// Construct the provider.
    pub const fn new() -> Self {
        Self
    }
}

impl Backend for SdfReconstruction {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(Self::ID, ExecutionTarget::PortableCpu)
    }
}

impl PointcloudReconstruction for SdfReconstruction {
    fn scratch_requirement(&self) -> ScratchRequirement {
        // The grid dominates: bounded by the extent over the edge length.
        // Declared unbounded because that product is a function of the
        // request, not of this provider.
        ScratchRequirement::Unbounded
    }

    fn determinism(&self) -> Determinism {
        // Every step is order-independent: the grid is walked in index
        // order and neighbour queries break ties by point index, so the
        // same input yields the same floats.
        Determinism::Bitwise
    }

    fn cancellation_granularity(&self) -> CancellationGranularity {
        // The level-set extraction is a single opaque call.
        CancellationGranularity::None
    }

    fn minimum_points(&self) -> usize {
        // Fewer than four points cannot bound a volume, so no surface can
        // be estimated from them.
        4
    }

    fn reconstruct(
        &self,
        cloud: &PointCloud,
        request: &ReconstructionRequest,
        options: &ExecutionOptions,
    ) -> GeomResult<Reconstruction> {
        let points = cloud.points();
        if points.len() < self.minimum_points() {
            return Ok(Reconstruction::Refused(
                ReconstructionRefusal::TooFewPoints {
                    supplied: points.len(),
                    required: self.minimum_points(),
                },
            ));
        }

        let Some((min_corner, max_corner)) = cloud.bounds() else {
            return Ok(Reconstruction::Refused(
                ReconstructionRefusal::DegenerateExtent {
                    detail: "cloud has no finite bounds".to_owned(),
                },
            ));
        };

        // A set with no thickness on some axis has no surface to
        // reconstruct: fitting one would produce a zero-volume sheet and
        // present it as a solid.
        let span = max_corner - min_corner;
        let extent = span.x.max(span.y).max(span.z);
        if extent <= 0.0 {
            return Ok(Reconstruction::Refused(
                ReconstructionRefusal::DegenerateExtent {
                    detail: "all points are coincident".to_owned(),
                },
            ));
        }
        let thinnest = span.x.min(span.y).min(span.z);
        if thinnest <= extent * 1e-12 {
            return Ok(Reconstruction::Refused(
                ReconstructionRefusal::DegenerateExtent {
                    detail: format!(
                        "points are collinear or coplanar: extent {extent} but thinnest axis {thinnest}"
                    ),
                },
            ));
        }

        let index = PointIndex::build(points);
        let spacing = median_spacing(&index, points);
        // Explicit about NaN: a spacing that is not a positive number means
        // the samples carry no usable scale, whether it is zero or NaN.
        if !spacing.is_finite() || spacing <= 0.0 {
            return Ok(Reconstruction::Refused(
                ReconstructionRefusal::DegenerateExtent {
                    detail: "every sample is coincident with its neighbour".to_owned(),
                },
            ));
        }

        let edge_length = match request.resolution {
            Resolution::FromSampleSpacing => spacing,
            Resolution::TargetEdgeLength(requested) => {
                if !requested.is_finite() || requested <= 0.0 {
                    return Ok(Reconstruction::Refused(
                        ReconstructionRefusal::Unsupported {
                            detail: format!("edge length {requested} is not a usable length"),
                        },
                    ));
                }
                // Reconstructing finer than the samples resolve would
                // present interpolation as measurement. Refuse rather than
                // silently clamp, so the caller learns the data's limit.
                if requested < spacing * 0.5 {
                    return Ok(Reconstruction::Refused(
                        ReconstructionRefusal::ResolutionExceedsData {
                            requested,
                            sample_spacing: spacing,
                        },
                    ));
                }
                requested
            }
        };

        let normals = if request.use_normals {
            cloud.normals()
        } else {
            None
        };
        let used_normals = normals.is_some();

        // The influence radius must span several samples so the field is
        // continuous between them; too small and the surface breaks into
        // disconnected blobs around each point.
        let influence = spacing * 2.5;

        let field = |probe: Point3| -> Scalar {
            signed_distance(&index, points, normals, probe, influence, spacing)
        };

        // Pad the extraction volume so a surface reaching the capture's
        // edge still closes rather than being clipped open.
        let pad = Vec3::splat(influence + edge_length * 2.0);
        let mut padded = Aabb::default();
        padded.extend(min_corner - pad);
        padded.extend(max_corner + pad);

        options.check_cancelled()?;

        let mesh = match level_set(field, padded, edge_length, 0.0) {
            Ok(mesh) => mesh,
            Err(error) => {
                return Ok(Reconstruction::Refused(
                    ReconstructionRefusal::Unsupported {
                        detail: format!("level-set extraction refused: {error}"),
                    },
                ));
            }
        };

        let health = audit_mesh(&mesh, options.tolerance());
        let closed = health.is_closed_two_manifold();

        if request.require_closed && !closed {
            return Ok(Reconstruction::Refused(
                ReconstructionRefusal::CannotClose {
                    detail: format!(
                        "extraction left {} boundary edges; the capture does not cover the whole object",
                        health.boundary_edges
                    ),
                },
            ));
        }

        // Count how much of the result rests on measured data. A triangle
        // whose centroid is further from any sample than the influence
        // radius was inferred, not observed.
        let interpolated = count_interpolated(&mesh, &index, points, influence);

        let mut evidence = ReconstructionEvidence::measured();
        evidence.input_points = points.len();
        evidence.used_points = points.len() - index.rejected();
        evidence.output_triangles = mesh.indices.len() / 3;
        evidence.output_components = component_count(&mesh);
        evidence.closed = closed;
        evidence.achieved_edge_length = edge_length;
        evidence.sample_spacing = spacing;
        evidence.used_normals = used_normals;
        evidence.interpolated_triangles = interpolated;

        Ok(Reconstruction::Surface(Box::new(
            ReconstructionOutcome::new(mesh, evidence),
        )))
    }
}

/// Median distance from each sample to its nearest other sample.
///
/// The scale the capture actually resolves. Median rather than mean so a
/// handful of outliers cannot inflate it — an outlier sitting far from the
/// surface would otherwise raise the estimated spacing and coarsen the
/// whole reconstruction.
fn median_spacing(index: &PointIndex, points: &[Point3]) -> Scalar {
    // Sampling is enough for a median and keeps this linear on huge
    // captures. Deterministic stride, not random, so the estimate is
    // reproducible.
    let stride = (points.len() / 512).max(1);
    let mut distances: Vec<Scalar> = Vec::new();
    let mut hits = Vec::new();
    for point in points.iter().step_by(stride) {
        if !point.is_finite() {
            continue;
        }
        // Two nearest: the first is the point itself at distance zero.
        if index.nearest_into(*point, 2, &mut hits).is_err() {
            continue;
        }
        if let Some(hit) = hits.iter().find(|h| h.distance > 0.0) {
            distances.push(hit.distance);
        }
    }
    if distances.is_empty() {
        return 0.0;
    }
    distances.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    distances[distances.len() / 2]
}

/// Signed distance from `probe` to the surface the samples lie on.
///
/// Positive outside, negative inside, matching the convention
/// `axiolid-levelset` expects.
fn signed_distance(
    index: &PointIndex,
    points: &[Point3],
    normals: Option<&[Vec3]>,
    probe: Point3,
    influence: Scalar,
    spacing: Scalar,
) -> Scalar {
    let mut hits = Vec::new();
    if index
        .nearest_into(probe, BLEND_NEIGHBOURS, &mut hits)
        .is_err()
        || hits.is_empty()
    {
        // No usable samples: report "far outside" rather than zero, which
        // the extractor would read as a surface crossing and fabricate
        // geometry out of nothing.
        return influence.max(spacing);
    }

    match normals {
        Some(normals) => {
            // Signed distance to each neighbour's tangent plane, blended by
            // inverse-square weight. The surface passes through the samples
            // because a probe exactly on a sample's plane scores zero.
            let mut weighted = 0.0;
            let mut total = 0.0;
            for hit in &hits {
                let normal = normals[hit.index];
                let length = normal.length();
                if length <= 0.0 {
                    continue;
                }
                let plane_distance = (probe - points[hit.index]).dot(normal / length);
                // Offset by a small epsilon so a probe sitting exactly on a
                // sample still has a defined weight.
                let weight = 1.0 / (hit.distance * hit.distance + spacing * spacing * 1e-6);
                weighted += plane_distance * weight;
                total += weight;
            }
            if total > 0.0 {
                weighted / total
            } else {
                // Normals present but all degenerate: fall back to the
                // unsigned estimate rather than dividing by zero.
                hits[0].distance - spacing
            }
        }
        None => {
            // No orientation information exists, so no inside can be
            // determined. Offsetting the unsigned distance produces a
            // surface at `spacing` around the samples: a shrink-wrap, not a
            // fit. Weaker, and reported as such in the evidence.
            hits[0].distance - spacing
        }
    }
}

/// Triangles resting on inference rather than measurement.
///
/// A triangle whose centroid is further from every sample than the
/// influence radius was extrapolated across a gap in the capture.
fn count_interpolated(
    mesh: &axiolid_mesh::TriMesh,
    index: &PointIndex,
    points: &[Point3],
    influence: Scalar,
) -> usize {
    let _ = points;
    let mut count = 0;
    for triangle in mesh.indices.chunks_exact(3) {
        let centroid = (mesh.positions[triangle[0] as usize]
            + mesh.positions[triangle[1] as usize]
            + mesh.positions[triangle[2] as usize])
            / 3.0;
        match index.nearest(centroid) {
            Ok(Some(hit)) if hit.distance <= influence => {}
            _ => count += 1,
        }
    }
    count
}

/// Connected components of the result, by shared vertex index.
///
/// A capture with gaps often reconstructs into several shells; reporting the
/// count lets a caller detect that here rather than downstream.
fn component_count(mesh: &axiolid_mesh::TriMesh) -> usize {
    if mesh.positions.is_empty() {
        return 0;
    }
    let mut parent: Vec<usize> = (0..mesh.positions.len()).collect();
    fn find(parent: &mut [usize], mut node: usize) -> usize {
        while parent[node] != node {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        node
    }
    for triangle in mesh.indices.chunks_exact(3) {
        let a = find(&mut parent, triangle[0] as usize);
        let b = find(&mut parent, triangle[1] as usize);
        let c = find(&mut parent, triangle[2] as usize);
        parent[b] = a;
        parent[c] = a;
    }
    // Only vertices actually used by a triangle count as a component.
    let mut used = vec![false; mesh.positions.len()];
    for &corner in &mesh.indices {
        used[corner as usize] = true;
    }
    let mut roots = std::collections::BTreeSet::new();
    for (vertex, _) in used.iter().enumerate().filter(|(_, used)| **used) {
        let root = find(&mut parent, vertex);
        roots.insert(root);
    }
    roots.len()
}
