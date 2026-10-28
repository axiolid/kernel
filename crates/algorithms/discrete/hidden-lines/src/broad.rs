//! The broad phase: which triangles can matter to an edge.
//!
//! A triangle is a candidate for an edge when the closed image rectangles of
//! their enclosures share a point and their depth ranges can order the right
//! way (a triangle wholly behind the edge cannot hide it; with a near plane,
//! a triangle wholly beyond it is removed). The boxes come from the sound
//! interval enclosures of the exact view-space vertices, so the filter never
//! drops a pair that the exact decision could keep.
//!
//! The boxes are indexed in an [`axiolid_spatial::Bvh`] over `(U, V, W)`; a
//! query is one box whose depth range encodes the depth test. Touching boxes
//! count (`Aabb::intersects` is closed), exactly as the linear scan does. The
//! answer is sorted ascending, so the candidates, their order and the pair
//! count are those of the linear scan. A box the hierarchy would refuse
//! (non-finite enclosure) makes the scan the answer; it never loses a pair.
//! `ScanHiddenLines` runs the scan throughout, for the benchmark baseline.

use axiolid_core::{Aabb, Point3};
use axiolid_spatial::{Bvh, SpatialItem};

use crate::edges::Tri;
use crate::view::Bounds;

/// A candidate query: the image rectangle and the depth range to keep.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Query {
    /// Only the image extent (`lo[0..2]`, `hi[0..2]`) is read.
    pub image: Bounds,
    /// A triangle needs `hi[2] >= w_lo`.
    pub w_lo: f64,
    /// A triangle needs `lo[2] <= w_hi`.
    pub w_hi: f64,
}

impl Query {
    fn admits(&self, tri: &Bounds) -> bool {
        tri.image_overlaps(&self.image) && tri.hi[2] >= self.w_lo && tri.lo[2] <= self.w_hi
    }
}

/// How candidates are found.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BroadPhase {
    /// The BVH (the provider's own).
    Tree,
    /// The linear scan, the definition the BVH must reproduce.
    Scan,
}

/// A hierarchy over the view-space boxes of the scene's triangles.
#[derive(Debug)]
pub(crate) struct Broad {
    bvh: Option<Bvh<usize>>,
}

impl Broad {
    pub(crate) fn build(tris: &[Tri], phase: BroadPhase) -> Self {
        if phase == BroadPhase::Scan {
            return Self { bvh: None };
        }
        let bvh = Bvh::build(tris.iter().enumerate().map(|(i, tri)| {
            let (lo, hi) = (tri.bounds.lo, tri.bounds.hi);
            SpatialItem::new(
                i,
                Aabb {
                    min: Point3::new(lo[0], lo[1], lo[2]),
                    max: Point3::new(hi[0], hi[1], hi[2]),
                },
            )
        }));
        Self {
            bvh: (bvh.rejected_items() == 0).then_some(bvh),
        }
    }

    /// The triangles `query` admits, ascending.
    pub(crate) fn candidates(&self, tris: &[Tri], query: &Query) -> Vec<usize> {
        let found = match &self.bvh {
            Some(bvh) => {
                let probe = Aabb {
                    min: Point3::new(query.image.lo[0], query.image.lo[1], query.w_lo),
                    max: Point3::new(query.image.hi[0], query.image.hi[1], query.w_hi),
                };
                let mut hits = Vec::new();
                bvh.query_aabb(&probe, &mut hits);
                // The hierarchy returns whole leaves: keep the boxes that
                // really overlap the probe.
                let mut keys: Vec<usize> = hits
                    .into_iter()
                    .filter_map(|i| bvh.item(i))
                    .filter(|item| query.admits(&tris[item.key].bounds))
                    .map(|item| item.key)
                    .collect();
                keys.sort_unstable();
                keys
            }
            _ => scan(tris, query),
        };
        #[cfg(test)]
        assert_eq!(
            found,
            scan(tris, query),
            "the hierarchy lost or added a pair"
        );
        found
    }
}

/// The linear scan: the definition the hierarchy must reproduce.
fn scan(tris: &[Tri], query: &Query) -> Vec<usize> {
    (0..tris.len())
        .filter(|&i| query.admits(&tris[i].bounds))
        .collect()
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use axiolid_contracts::ExecutionOptions;
    use axiolid_core::{Frame3, Point3, Tolerance, Vec3};
    use axiolid_mesh::TriMesh;
    use axiolid_mesh_hidden_lines_contract::{
        EdgeSelection, HiddenLineLimits, HiddenLineMesh, HiddenLineRequest, MeshHiddenLines,
        MeshId, MeshRole, OrthographicView,
    };

    use super::*;
    use crate::provider::{ScalarHiddenLines, ScanHiddenLines};

    fn tri(lo: [f64; 3], hi: [f64; 3]) -> Tri {
        Tri {
            v: [0, 0, 0],
            sigma: Ordering::Less,
            bounds: Bounds { lo, hi },
        }
    }

    fn query(lo: [f64; 3], hi: [f64; 3], w_lo: f64, w_hi: f64) -> Query {
        Query {
            image: Bounds { lo, hi },
            w_lo,
            w_hi,
        }
    }

    #[test]
    fn touching_boxes_are_candidates_wherever_the_leaves_fall() {
        let tris: Vec<Tri> = (0..40)
            .map(|i| {
                let x = f64::from(i);
                tri([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0])
            })
            .collect();
        let broad = Broad::build(&tris, BroadPhase::Tree);
        let mut touching = 0;
        // Integer x and leaf boundaries, y and w on both sides of every face.
        for x in (-1..=41).map(f64::from).chain([0.5, 7.5]) {
            for width in [0.0, 1.0] {
                for y in [-1.0, 0.0, 0.5, 1.0, 2.0] {
                    for w_lo in [-1.0, 0.0, 1.0, 2.0] {
                        for w_hi in [-1.0, 0.0, 1.0, 2.0, f64::INFINITY] {
                            let q = query([x, y, 0.0], [x + width, y, 0.0], w_lo, w_hi);
                            let found = broad.candidates(&tris, &q);
                            assert_eq!(found, scan(&tris, &q), "x {x} y {y} w {w_lo}..{w_hi}");
                            touching += usize::from(!found.is_empty());
                        }
                    }
                }
            }
        }
        assert!(touching > 1000);
        // The documented touching cases, spelled out.
        let q = query([5.0, 0.0, 0.0], [5.0, 0.0, 0.0], 1.0, 0.0);
        assert_eq!(broad.candidates(&tris, &q), vec![4, 5]);
        // A depth gap: nothing.
        let q = query([5.0, 0.0, 0.0], [5.0, 0.0, 0.0], 1.5, 5.0);
        assert!(broad.candidates(&tris, &q).is_empty());
    }

    #[test]
    fn a_non_finite_enclosure_falls_back_to_the_scan() {
        let tris = vec![
            tri([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            tri([f64::NEG_INFINITY, 0.0, 0.0], [f64::INFINITY, 1.0, 1.0]),
        ];
        let broad = Broad::build(&tris, BroadPhase::Tree);
        assert!(broad.bvh.is_none());
        let q = query([0.5, 0.5, 0.0], [0.5, 0.5, 0.0], 0.0, 1.0);
        assert_eq!(broad.candidates(&tris, &q), vec![0, 1]);
    }

    fn boxes(count: usize, seed: u64) -> Vec<TriMesh> {
        let mut state = seed | 1;
        let mut below = |n: u64| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 33) % n) as f64
        };
        (0..count)
            .map(|_| {
                let (x0, y0, z0) = (below(9) / 4.0, below(9) / 4.0, below(5) / 4.0);
                let (x1, y1, z1) = (
                    x0 + (1.0 + below(4)) / 4.0,
                    y0 + (1.0 + below(4)) / 4.0,
                    z0 + (1.0 + below(4)) / 4.0,
                );
                TriMesh::new(
                    vec![
                        Point3::new(x0, y0, z0),
                        Point3::new(x1, y0, z0),
                        Point3::new(x1, y1, z0),
                        Point3::new(x0, y1, z0),
                        Point3::new(x0, y0, z1),
                        Point3::new(x1, y0, z1),
                        Point3::new(x1, y1, z1),
                        Point3::new(x0, y1, z1),
                    ],
                    vec![
                        0, 3, 2, 0, 2, 1, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2,
                        3, 7, 2, 7, 6, 3, 0, 4, 3, 4, 7,
                    ],
                )
            })
            .collect()
    }

    fn draw(
        provider: &dyn MeshHiddenLines,
        meshes: &[TriMesh],
        frame: Frame3,
        near: Option<f64>,
    ) -> String {
        let entries: Vec<HiddenLineMesh<'_>> = meshes
            .iter()
            .enumerate()
            .map(|(i, m)| HiddenLineMesh::new(MeshId(i as u64 + 1), MeshRole::Solid, m))
            .collect();
        let mut view = OrthographicView::new(frame);
        if let Some(near) = near {
            view = view.with_near(near);
        }
        let limits = HiddenLineLimits::new(64, 4096, 4096, 10_000_000, 1_000_000);
        let request = HiddenLineRequest::new(&entries, view, EdgeSelection::new(Some(0.5)), limits);
        let options = ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"));
        let outcome = provider
            .hidden_lines(&request, &options)
            .expect("the scene draws");
        format!("{outcome:?}")
    }

    #[test]
    fn hierarchy_and_scan_draw_identical_outcomes() {
        let s6 = 1.0 / 6f64.sqrt();
        let s3 = 1.0 / 3f64.sqrt();
        let iso = Frame3 {
            origin: Point3::new(0.1, 0.2, 0.3),
            x: Vec3::new(
                std::f64::consts::FRAC_1_SQRT_2,
                -std::f64::consts::FRAC_1_SQRT_2,
                0.0,
            ),
            y: Vec3::new(s6, s6, -2.0 * s6),
            z: Vec3::new(s3, s3, s3),
        };
        let plan = Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        };
        // Looking along +Y (an elevation) and along a generic direction.
        let elevation = Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Z,
            z: -Vec3::Y,
        };
        let (a, b, c) = (0.3_f64, -0.5_f64, 0.8_f64);
        let norm = (a * a + b * b + c * c).sqrt();
        let z = Vec3::new(a / norm, b / norm, c / norm);
        let x_raw = Vec3::new(b, -a, 0.0);
        let xn = (a * a + b * b).sqrt();
        let x = Vec3::new(x_raw.x / xn, x_raw.y / xn, 0.0);
        let y = Vec3::new(
            z.y * x.z - z.z * x.y,
            z.z * x.x - z.x * x.z,
            z.x * x.y - z.y * x.x,
        );
        let oblique = Frame3 {
            origin: Point3::new(-0.4, 0.7, 0.1),
            x,
            y,
            z,
        };
        let mut checked = 0;
        for seed in 0..12 {
            let meshes = boxes(3 + seed as usize, seed);
            for frame in [iso, plan, elevation, oblique] {
                for near in [None, Some(0.55)] {
                    let scanned = draw(&ScanHiddenLines, &meshes, frame, near);
                    // The assertion in `candidates` also compares every query.
                    let indexed = draw(&ScalarHiddenLines::new(), &meshes, frame, near);
                    assert_eq!(scanned, indexed, "seed {seed}");
                    checked += 1;
                }
            }
        }
        assert_eq!(checked, 96);
    }
}
