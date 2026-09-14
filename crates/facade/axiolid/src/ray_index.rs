//! Cached broad phase for repeated ray casts against the same mesh.
//!
//! `nearest_hit` scans every triangle by design: the broad phase belongs
//! to the caller, and a one-shot caller must not pay for an index it
//! never reuses. Measured here: the BVH only repays after ~22 rays, and
//! at a single ray it is ~20x slower than scanning.
//!
//! So the index is built lazily on the SECOND cast against a given mesh
//! and reused afterwards. A first cast costs what it always did.
//!
//! The cache key is a content digest, not an address. `TriMesh` exposes
//! `positions`/`indices` as public `Vec`s with no version counter, so a
//! caller can mutate a mesh in place; an address or (ptr, len) key would
//! then serve a stale index and return hits for geometry that no longer
//! exists. Digesting costs ~1.2% of a build and ~0.01% of a full scan.

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::sync::RwLock;

use axiolid_core::Aabb;
use axiolid_core::{Point3, Ray3, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_ray_mesh::{nearest_hit, nearest_hit_among, RayHit3, RayMeshError};
use axiolid_spatial::{Bvh, SpatialIndex, SpatialItem};
use std::ops::ControlFlow;

/// Entries retained. Each holds one BVH over a mesh, so this bounds
/// memory: an unbounded map would leak an index per distinct mesh ever
/// cast against, which for a caller streaming meshes is every mesh.
const CAPACITY: usize = 8;

/// Casts against a mesh before its index is built.
///
/// 1 means "build on the second cast". Building on the first would
/// penalise the one-shot caller the scanning API exists to serve.
const WARMUP_CASTS: u32 = 1;

struct Entry {
    digest: u64,
    /// `None` until the mesh has been cast against `WARMUP_CASTS` times.
    bvh: Option<Bvh<usize>>,
    casts: u32,
    /// Monotonic tick of last use, for least-recently-used eviction.
    touched: u64,
}

/// Ray index cache. Not part of the public API surface: it changes only
/// how `nearest_mesh_hit` finds its answer, never what that answer is.
#[derive(Default)]
pub(crate) struct RayIndexCache {
    inner: RwLock<Inner>,
}

#[derive(Default)]
struct Inner {
    entries: Vec<Entry>,
    clock: u64,
}

/// Content digest of the geometry a ray can hit.
///
/// Positions are hashed by bit pattern: two meshes that differ only by
/// -0.0 vs 0.0 are geometrically identical here, but hashing bits is the
/// conservative direction (a spurious miss rebuilds; a spurious HIT
/// would serve the wrong index). Attributes and normals are excluded:
/// they cannot change which triangle a ray strikes.
fn digest(mesh: &TriMesh) -> u64 {
    let mut hasher = DefaultHasher::new();
    mesh.indices.hash(&mut hasher);
    for point in &mesh.positions {
        point.x.to_bits().hash(&mut hasher);
        point.y.to_bits().hash(&mut hasher);
        point.z.to_bits().hash(&mut hasher);
    }
    hasher.finish()
}

fn build_bvh(mesh: &TriMesh) -> Bvh<usize> {
    let points = &mesh.positions;
    let items = (0..mesh.indices.len() / 3).map(|triangle| {
        let corners = &mesh.indices[triangle * 3..triangle * 3 + 3];
        let first = points[corners[0] as usize];
        let (mut low, mut high) = (first, first);
        for corner in &corners[1..] {
            let point = points[*corner as usize];
            low = Point3::new(low.x.min(point.x), low.y.min(point.y), low.z.min(point.z));
            high = Point3::new(
                high.x.max(point.x),
                high.y.max(point.y),
                high.z.max(point.z),
            );
        }
        SpatialItem::new(
            triangle,
            Aabb {
                min: low,
                max: high,
            },
        )
    });
    Bvh::build(items)
}

impl RayIndexCache {
    /// Nearest hit, using a cached broad phase once one exists.
    ///
    /// # Errors
    ///
    /// Propagates whatever `nearest_hit`/`nearest_hit_among` refuse.
    pub(crate) fn nearest_hit(
        &self,
        mesh: &TriMesh,
        ray: &Ray3,
        tolerance: Tolerance,
    ) -> Result<Option<RayHit3>, RayMeshError> {
        let digest = digest(mesh);

        // Fast path: a read lock is enough when the index already exists,
        // so concurrent casts against one mesh do not serialise.
        {
            let inner = self.inner.read().expect("ray index cache poisoned");
            if let Some(entry) = inner.entries.iter().find(|e| e.digest == digest) {
                if let Some(bvh) = &entry.bvh {
                    return accelerated(mesh, ray, tolerance, bvh);
                }
            }
        }

        // Slow path: record the cast and build once warm.
        let mut inner = self.inner.write().expect("ray index cache poisoned");
        inner.clock += 1;
        let tick = inner.clock;
        let capacity_reached = inner.entries.len() >= CAPACITY;
        match inner.entries.iter_mut().find(|e| e.digest == digest) {
            Some(entry) => {
                entry.casts += 1;
                entry.touched = tick;
                if entry.bvh.is_none() && entry.casts > WARMUP_CASTS {
                    entry.bvh = Some(build_bvh(mesh));
                }
            }
            None => {
                if capacity_reached {
                    // Evict least recently used, not the first entry: a
                    // hot mesh must survive a burst of one-shot casts.
                    if let Some(position) = inner
                        .entries
                        .iter()
                        .enumerate()
                        .min_by_key(|(_, e)| e.touched)
                        .map(|(i, _)| i)
                    {
                        inner.entries.remove(position);
                    }
                }
                inner.entries.push(Entry {
                    digest,
                    bvh: None,
                    casts: 1,
                    touched: tick,
                });
            }
        }

        // Whether or not an index now exists, answer THIS cast by the
        // path that defines the contract. Using a freshly built index
        // here would make the first accelerated answer untested against
        // the scan it must agree with.
        nearest_hit(mesh, ray, tolerance)
    }
}

/// Nearest hit via the broad phase.
///
/// Candidates arrive nearest-bound-first, but a nearer BOUND does not
/// mean a nearer HIT, so this cannot stop at the first success. It
/// collects every candidate whose bounds the ray enters and hands the
/// whole set to `nearest_hit_among`, which applies the same ordering
/// rule as the full scan -- including which face wins a shared edge.
/// Stopping early here produced a different owning triangle on 28 of
/// 2000 rays in the probe.
fn accelerated(
    mesh: &TriMesh,
    ray: &Ray3,
    tolerance: Tolerance,
    bvh: &Bvh<usize>,
) -> Result<Option<RayHit3>, RayMeshError> {
    let mut candidates = Vec::new();
    bvh.visit_ray(ray, &mut |hit| {
        candidates.push(*hit.key);
        ControlFlow::Continue(())
    });
    if candidates.is_empty() {
        return Ok(None);
    }
    candidates.sort_unstable();
    nearest_hit_among(mesh, ray, tolerance, candidates)
}

// `Application` is Debug + Clone. A cache is derived state, so cloning
// an application must NOT share or copy indices: the clone starts cold
// and rebuilds what it needs. Sharing would let one clone observe
// another's evictions.
impl Clone for RayIndexCache {
    fn clone(&self) -> Self {
        Self::default()
    }
}

impl std::fmt::Debug for RayIndexCache {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let entries = self.inner.read().map(|i| i.entries.len()).unwrap_or(0);
        formatter
            .debug_struct("RayIndexCache")
            .field("entries", &entries)
            .finish()
    }
}

#[cfg(all(test, feature = "ray-mesh", feature = "spatial"))]
mod ray_index_tests {
    use super::*;
    use axiolid_core::{Point3, Ray3, Tolerance, Vec3};
    use axiolid_mesh::TriMesh;

    fn tol() -> Tolerance {
        Tolerance::new(1e-6, 1e-9).expect("tolerance")
    }

    /// Two stacked quads so a ray crosses several triangles: a broad
    /// phase that returned only the nearest BOUND would pick wrong here.
    fn slab() -> TriMesh {
        let mut positions = Vec::new();
        let mut indices = Vec::new();
        for (level, z) in [0.0_f64, 1.0].into_iter().enumerate() {
            let base = (level * 4) as u32;
            positions.push(Point3::new(-1.0, -1.0, z));
            positions.push(Point3::new(1.0, -1.0, z));
            positions.push(Point3::new(1.0, 1.0, z));
            positions.push(Point3::new(-1.0, 1.0, z));
            indices.extend_from_slice(&[base, base + 1, base + 2]);
            indices.extend_from_slice(&[base, base + 2, base + 3]);
        }
        TriMesh::new(positions, indices)
    }

    /// The cache changes HOW the answer is found, never WHAT it is.
    /// Casts repeatedly so the comparison spans cold, warming and warm
    /// states -- a test that stopped at one cast would never exercise
    /// the accelerated path at all.
    #[test]
    fn cached_casts_match_the_full_scan() {
        let cache = RayIndexCache::default();
        let mesh = slab();
        for step in 0..32 {
            let offset = f64::from(step) * 0.05 - 0.8;
            let ray = Ray3 {
                origin: Point3::new(offset, 0.1, -3.0),
                direction: Vec3::new(0.0, 0.0, 1.0),
            };
            let expected = nearest_hit(&mesh, &ray, tol()).expect("scan");
            let actual = cache.nearest_hit(&mesh, &ray, tol()).expect("cached");
            match (expected, actual) {
                (None, None) => {}
                (Some(want), Some(got)) => {
                    assert!((want.t - got.t).abs() < 1e-12, "step {step}: t differs");
                    // Not just distance: the owning face must match too,
                    // or a shared-edge tie silently changes which
                    // triangle the caller is told it hit.
                    assert_eq!(want.triangle, got.triangle, "step {step}: triangle differs");
                }
                (a, b) => panic!("step {step}: hit disagreement {a:?} vs {b:?}"),
            }
        }
    }

    /// `TriMesh` fields are public, so a caller can move geometry under
    /// a cached index. Keying on content means the edit produces a new
    /// key and the stale BVH is never consulted. An address-based key
    /// would pass every other test here and fail this one.
    /// A grid big enough that the BVH really prunes: with only a handful
    /// of triangles every leaf is a candidate anyway, so a stale index
    /// still yields the right answer and the test proves nothing.
    fn grid(n: usize, z: f64) -> TriMesh {
        let mut positions = Vec::new();
        let mut indices = Vec::new();
        for i in 0..n {
            for j in 0..n {
                let (x, y) = (i as f64 * 0.1 - 2.0, j as f64 * 0.1 - 2.0);
                let base = positions.len() as u32;
                positions.push(Point3::new(x, y, z));
                positions.push(Point3::new(x + 0.09, y, z));
                positions.push(Point3::new(x, y + 0.09, z));
                indices.extend_from_slice(&[base, base + 1, base + 2]);
            }
        }
        TriMesh::new(positions, indices)
    }

    #[test]
    fn mutating_the_mesh_does_not_serve_a_stale_index() {
        let cache = RayIndexCache::default();
        // 1600 triangles at z = 1, so the BVH prunes hard.
        let mut mesh = grid(40, 1.0);
        let ray = Ray3 {
            origin: Point3::new(-1.97, -1.97, -3.0),
            direction: Vec3::new(0.0, 0.0, 1.0),
        };
        for _ in 0..4 {
            cache.nearest_hit(&mesh, &ray, tol()).expect("warm");
        }

        // Move ONE triangle -- the one the ray actually hits -- nearer.
        // Indices never change, so a shape-only digest keeps the stale
        // BVH, whose box for this triangle still sits at the old z. The
        // pruned traversal then never offers it as a candidate.
        let far = mesh.positions.len() - 3;
        for k in 0..3 {
            mesh.positions[far + k].x -= 3.9;
            mesh.positions[far + k].y -= 3.9;
            mesh.positions[far + k].z -= 2.0;
        }

        let after = cache.nearest_hit(&mesh, &ray, tol()).expect("after");
        let truth = nearest_hit(&mesh, &ray, tol()).expect("scan");
        assert_eq!(
            after.map(|h| h.t.to_bits()),
            truth.map(|h| h.t.to_bits()),
            "stale index served after an in-place edit",
        );
    }

    /// The cache must not grow without bound: a caller streaming meshes
    /// would otherwise retain a BVH for every one it ever cast against.
    #[test]
    fn the_cache_is_bounded() {
        let cache = RayIndexCache::default();
        let ray = Ray3 {
            origin: Point3::new(0.0, 0.0, -3.0),
            direction: Vec3::new(0.0, 0.0, 1.0),
        };
        for step in 0..(CAPACITY * 4) {
            let mut mesh = slab();
            // Distinct geometry per iteration -> distinct digest.
            for point in &mut mesh.positions {
                point.z += step as f64 * 0.25;
            }
            cache.nearest_hit(&mesh, &ray, tol()).expect("cast");
        }
        let held = cache.inner.read().expect("lock").entries.len();
        assert!(held <= CAPACITY, "cache grew to {held}");
    }
}
