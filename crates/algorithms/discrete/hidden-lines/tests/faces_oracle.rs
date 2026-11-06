//! An independent oracle for the visible-face regions: plain binary64 ray
//! casting on a sample grid.
//!
//! It shares no code with the provider. For each sample of a 64x64 grid over
//! the image it casts a ray through the scene in binary64 and records the
//! frontmost source: the frontmost kept triangle's mesh, or the cap of a solid
//! the near plane cuts (decided by the parity of the triangles of that solid
//! beyond the plane, as in `oracle.rs`). The sample must then lie, by an
//! ordinary point-in-polygon test, in exactly the regions of that source and
//! in none of the others.
//!
//! A sample is skipped as ambiguous when it lies within `1e-7` of a ring edge
//! of the output, or when the two frontmost sources are within `1e-7` of equal
//! depth (a tie owns both sources, which a binary64 depth cannot tell from a
//! near miss). At least 90% of the samples must be checked.
//!
//! Scenes are the seeded boxes of `oracle.rs`, some rotated by generic angles,
//! in the plan, an elevation, the isometric and the oblique views, with and
//! without a near plane. Every scene must be answered: a region that
//! collapses under rounding is snapped (ADR 0100), never refused. The
//! displacement the evidence reports must stay below a quarter of `1e-7`,
//! so the ambiguity margin above covers it, and every outcome must pass the
//! independent topology checker of `common::topology`.

mod common;

use axiolid_core::{Frame3, Point2};
use axiolid_hidden_lines::ScalarVisibleFaces;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    FaceGroup, FaceSource, OrthographicView, VisibleFaceOutcome,
};
use common::topology::assert_valid_topology;
use common::{scene, solids, try_faces, views};

const GRID: usize = 64;
const EPS: f64 = 1e-7;

type V3 = [f64; 3];
type Key = (usize, FaceSource);

fn to_view(mesh: &TriMesh, frame: &Frame3) -> Vec<V3> {
    mesh.positions
        .iter()
        .map(|p| {
            let d = *p - frame.origin;
            [d.dot(frame.x), d.dot(frame.y), d.dot(frame.z)]
        })
        .collect()
}

/// The image `(u, v)` and depth of the triangle `t` at `q`, when `q` is
/// inside it (closed, to rounding) and the triangle is not edge-on.
fn cover(v: &[V3], t: [u32; 3], q: [f64; 2]) -> Option<f64> {
    let [a, b, c] = t.map(|i| v[i as usize]);
    let area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    if area.abs() < 1e-12 {
        return None;
    }
    let w1 = ((q[0] - a[0]) * (c[1] - a[1]) - (q[1] - a[1]) * (c[0] - a[0])) / area;
    let w2 = ((b[0] - a[0]) * (q[1] - a[1]) - (b[1] - a[1]) * (q[0] - a[0])) / area;
    let w0 = 1.0 - w1 - w2;
    (w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0).then(|| w0 * a[2] + w1 * b[2] + w2 * c[2])
}

/// The frontmost source at `q`, or `None` when it is ambiguous (a near tie).
fn frontmost(
    meshes: &[(Vec<V3>, Vec<[u32; 3]>)],
    near: Option<f64>,
    cut: &[bool],
    q: [f64; 2],
) -> Result<Option<Key>, ()> {
    let mut found: Vec<(f64, Key)> = Vec::new();
    for (index, (v, tris)) in meshes.iter().enumerate() {
        let mut beyond = 0;
        for &t in tris {
            if let Some(w) = cover(v, t, q) {
                if near.is_some_and(|n| w > n) {
                    beyond += 1;
                } else {
                    found.push((w, (index, FaceSource::Group(FaceGroup(0)))));
                }
            }
        }
        if let (Some(n), true) = (near, cut[index]) {
            if beyond % 2 == 1 {
                found.push((n, (index, FaceSource::Cap)));
            }
        }
    }
    found.sort_by(|a, b| b.0.total_cmp(&a.0));
    let Some(&(depth, key)) = found.first() else {
        return Ok(None);
    };
    // Another source within the tolerance: a tie.
    if found.iter().any(|&(w, k)| k != key && depth - w < EPS) {
        return Err(());
    }
    Ok(Some(key))
}

fn distance_to_segment(p: Point2, a: Point2, b: Point2) -> f64 {
    let (ab, ap) = (b - a, p - a);
    let t = if ab.length_squared() == 0.0 {
        0.0
    } else {
        (ap.dot(ab) / ab.length_squared()).clamp(0.0, 1.0)
    };
    (p - (a + ab * t)).length()
}

fn inside(ring: &[Point2], p: Point2) -> bool {
    let mut odd = false;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        if (a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x) {
            odd = !odd;
        }
    }
    odd
}

fn near_a_ring_edge(outcome: &VisibleFaceOutcome, p: Point2) -> bool {
    outcome.regions.iter().any(|r| {
        std::iter::once(&r.outer).chain(&r.holes).any(|ring| {
            (0..ring.len())
                .any(|i| distance_to_segment(p, ring[i], ring[(i + 1) % ring.len()]) < EPS)
        })
    })
}

#[derive(Default)]
struct Tally {
    checked: usize,
    skipped: usize,
    occupied: usize,
    snapped: usize,
}

fn check(seed: u64, name: &str, frame: Frame3, with_near: bool, tally: &mut Tally) {
    let boxes = scene(seed);
    let meshes = solids(&boxes);
    let projected: Vec<(Vec<V3>, Vec<[u32; 3]>)> = boxes
        .iter()
        .map(|m| (to_view(m, &frame), m.indices.as_chunks::<3>().0.to_vec()))
        .collect();
    let depths = projected.iter().flat_map(|(v, _)| v.iter().map(|p| p[2]));
    let (lo, hi) = depths.fold((f64::MAX, f64::MIN), |(l, h), w| (l.min(w), h.max(w)));
    let near = with_near.then(|| (lo + hi) / 2.0);
    let cut: Vec<bool> = projected
        .iter()
        .map(|(v, _)| {
            near.is_some_and(|n| v.iter().any(|p| p[2] > n) && v.iter().any(|p| p[2] <= n))
        })
        .collect();
    let mut view = OrthographicView::new(frame);
    if let Some(n) = near {
        view = view.with_near(n);
    }
    let outcome = try_faces(&ScalarVisibleFaces::new(), &meshes, view, None)
        .unwrap_or_else(|error| panic!("seed {seed}, {name}, near {with_near}: {error}"));
    assert!(
        outcome.evidence.snap_displacement.get() * 4.0 < EPS,
        "seed {seed}, {name}: displacement {}",
        outcome.evidence.snap_displacement.get()
    );
    assert_valid_topology(&outcome);
    let evidence = &outcome.evidence;
    tally.snapped += usize::from(
        evidence.rerouted_edges + evidence.collapsed_regions + evidence.collapsed_holes > 0,
    );

    let (mut u_lo, mut u_hi, mut v_lo, mut v_hi) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for p in projected.iter().flat_map(|(v, _)| v) {
        (u_lo, u_hi) = (u_lo.min(p[0]), u_hi.max(p[0]));
        (v_lo, v_hi) = (v_lo.min(p[1]), v_hi.max(p[1]));
    }
    let (du, dv) = ((u_hi - u_lo) * 1.05, (v_hi - v_lo) * 1.05);
    for i in 0..GRID {
        for j in 0..GRID {
            // Cell centres of an irrational stride: no sample sits on a
            // grid-aligned box edge.
            let q = [
                u_lo - 0.025 * (u_hi - u_lo) + du * (i as f64 + 0.4142) / GRID as f64,
                v_lo - 0.025 * (v_hi - v_lo) + dv * (j as f64 + 0.7320) / GRID as f64,
            ];
            let p = Point2::new(q[0], q[1]);
            let Ok(expected) = frontmost(&projected, near, &cut, q) else {
                tally.skipped += 1;
                continue;
            };
            if near_a_ring_edge(&outcome, p) {
                tally.skipped += 1;
                continue;
            }
            tally.checked += 1;
            tally.occupied += usize::from(expected.is_some());
            // The sources whose regions hold the sample.
            let mut held: Vec<Key> = Vec::new();
            for region in &outcome.regions {
                let index = meshes
                    .iter()
                    .position(|m| m.id == region.mesh)
                    .expect("region of a request mesh");
                if inside(&region.outer, p) && !region.holes.iter().any(|h| inside(h, p)) {
                    held.push((index, region.source));
                }
            }
            held.sort();
            held.dedup();
            assert_eq!(
                held,
                expected.into_iter().collect::<Vec<_>>(),
                "seed {seed}, {name}, near {near:?}: sample {q:?}"
            );
        }
    }
}

#[test]
fn regions_agree_with_binary64_ray_casting() {
    let mut tally = Tally::default();
    for seed in 0..12 {
        for (name, frame) in views() {
            for with_near in [false, true] {
                check(seed, name, frame, with_near, &mut tally);
            }
        }
    }
    let total = tally.checked + tally.skipped;
    assert!(
        tally.checked * 10 >= total * 9,
        "checked {} of {total} samples",
        tally.checked
    );
    assert!(
        tally.occupied > 5_000,
        "only {} covered samples",
        tally.occupied
    );
    eprintln!(
        "faces oracle: {} checked ({} covered), {} skipped, {} scenes snapped of {}",
        tally.checked,
        tally.occupied,
        tally.skipped,
        tally.snapped,
        12 * views().len() * 2
    );
}

/// The two scenes that binary64 rounding used to refuse as `Degenerate`
/// (a region whose vertices round to one point) are answered, and the
/// answer passes the oracle like any other scene.
#[test]
fn the_known_collapses_are_snapped() {
    for (seed, view) in [(1, "oblique"), (9, "isometric")] {
        let (name, frame) = views()
            .into_iter()
            .find(|(name, _)| *name == view)
            .expect("the view");
        let mut tally = Tally::default();
        check(seed, name, frame, true, &mut tally);
        assert_eq!(tally.snapped, 1, "seed {seed}, {name}");
        assert!(tally.checked > 0, "seed {seed}, {name}");
    }
}
