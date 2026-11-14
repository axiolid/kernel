//! An independent oracle for see-through triangles (ADR 0111): plain
//! binary64 ray casting, sharing no code with the provider.
//!
//! Scenes are the seeded boxes of `oracle.rs`, with every box whose seeded
//! draw is 0 of 3 wholly see-through, in the plan, two elevations, the
//! isometric and the oblique views, with and without a near plane through
//! the middle of one box; the hidden lines also in a pinhole view.
//!
//! - Lines: each output segment is sampled at 1/4, 1/2 and 3/4 of its
//!   length. The ray through the sample meets the scene; the segment must be
//!   `Hidden` exactly when a kept opaque triangle or the cap of a cut solid
//!   lies strictly in front of the edge's point there, and otherwise count
//!   exactly the kept see-through triangles strictly in front.
//! - Faces: on a 64x64 grid over the image, the kept surfaces crossed by the
//!   ray are sorted front to back and cut after the first opaque one (a cap
//!   is opaque). The sample must lie in exactly one region of the frontmost
//!   source, and that region's `behind` must list the rest of the cut list,
//!   one layer per surface, when the frontmost is see-through, else nothing.
//!
//! A sample is skipped as ambiguous when it lies within `1e-7` of the image
//! edge of a triangle (or, for faces, of an output ring), or two depths that
//! decide it are within `1e-7` (a tie, which binary64 cannot tell from a near
//! miss). At least 90% of the face samples must be checked. Line samples lie
//! on an edge, so they are ambiguous wherever another triangle's image edge
//! runs along that edge's image, and how often that happens depends on the
//! view:
//!
//! - oblique and pinhole: the direction is generic, so only chance
//!   coincidences remain; 85% must be checked (about 89% and 91% are);
//! - isometric: two thirds of the boxes are axis-aligned on a 1/8 grid, and
//!   the isometric direction maps many grid edges of different boxes onto
//!   common image lines; 75% must be checked (about 77% to 80% are);
//! - the axis views: every box edge lies on the image line of a side of its
//!   own box seen edge-on, so about half the samples are ambiguous and only a
//!   minimum count is required.

mod common;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame3, Point2, Point3, Tolerance, Vec3};
use axiolid_hidden_lines::{ScalarHiddenLines, ScalarVisibleFaces};
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    BehindSurface, EdgeSelection, FaceGroup, FaceSource, HiddenLineMesh, HiddenLineOutcome,
    HiddenLineRequest, MeshHiddenLines, MeshVisibleFaces, Occlusion, View, ViewProjection,
    Visibility, VisibleFaceOutcome, VisibleFaceRequest,
};
use common::topology::assert_valid_topology;
use common::{scene, solids, views, Lcg, FACE_LIMITS, LIMITS};

const EPS: f64 = 1e-7;
const GRID: usize = 64;

type V3 = [f64; 3];

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"))
}

/// Which boxes of the scene `seed` are see-through: a seeded third.
fn glazed(seed: u64, count: usize) -> Vec<bool> {
    let mut rng = Lcg(seed.wrapping_mul(0xA24B_AED4_963E_E407) | 1);
    let mut flags: Vec<bool> = (0..count).map(|_| rng.below(3) == 0).collect();
    // At least one box see-through and, when there are two, one opaque.
    flags[0] = true;
    if count > 1 && flags.iter().all(|&f| f) {
        flags[1] = false;
    }
    flags
}

/// The scene in view space: per mesh the vertices `[U, V, W]`, the
/// triangles and whether they are see-through.
struct Projected {
    meshes: Vec<(Vec<V3>, Vec<[u32; 3]>, bool)>,
    perspective: bool,
    near: Option<f64>,
    /// Per mesh, whether the near plane cuts it (orthographic only).
    cut: Vec<bool>,
}

fn to_view(mesh: &TriMesh, frame: &Frame3) -> Vec<V3> {
    mesh.positions
        .iter()
        .map(|p| {
            let d = *p - frame.origin;
            [d.dot(frame.x), d.dot(frame.y), d.dot(frame.z)]
        })
        .collect()
}

impl Projected {
    fn new(boxes: &[TriMesh], glass: &[bool], view: &View) -> Self {
        let meshes: Vec<_> = boxes
            .iter()
            .zip(glass)
            .map(|(m, &g)| {
                (
                    to_view(m, &view.frame),
                    m.indices.as_chunks::<3>().0.to_vec(),
                    g,
                )
            })
            .collect();
        let cut = meshes
            .iter()
            .map(|(v, _, _)| {
                view.near
                    .is_some_and(|n| v.iter().any(|p| p[2] > n) && v.iter().any(|p| p[2] <= n))
            })
            .collect();
        Self {
            meshes,
            perspective: view.projection == ViewProjection::Perspective,
            near: view.near,
            cut,
        }
    }

    /// The image of a view-space point.
    fn image(&self, p: V3) -> [f64; 2] {
        if self.perspective {
            [p[0] / -p[2], p[1] / -p[2]]
        } else {
            [p[0], p[1]]
        }
    }
}

/// Where the ray through the image point `q` meets the triangle `t`: its
/// depth `W`, `None` when it misses or the triangle is seen edge-on, and
/// `Err` when `q` is within `EPS` of the triangle's image edges.
fn crossing(projected: &Projected, v: &[V3], t: [u32; 3], q: [f64; 2]) -> Result<Option<f64>, ()> {
    let corners = t.map(|i| v[i as usize]);
    let image = corners.map(|p| projected.image(p));
    for k in 0..3 {
        let (a, b) = (image[k], image[(k + 1) % 3]);
        if distance(q, a, b) < EPS {
            return Err(());
        }
    }
    let area = cross(image[0], image[1], image[2]);
    if area.abs() < 1e-12 {
        return Ok(None);
    }
    let w: Vec<f64> = (0..3)
        .map(|k| cross(image[(k + 1) % 3], image[(k + 2) % 3], q) / area)
        .collect();
    if w.iter().any(|&w| w < 0.0) {
        return Ok(None);
    }
    if projected.perspective {
        // Inverse depth is affine in the image.
        let inverse: f64 = (0..3).map(|k| w[k] / corners[k][2]).sum();
        Ok(Some(1.0 / inverse))
    } else {
        Ok(Some((0..3).map(|k| w[k] * corners[k][2]).sum()))
    }
}

fn cross(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

fn distance(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let (p, a, b) = (
        Point2::new(p[0], p[1]),
        Point2::new(a[0], a[1]),
        Point2::new(b[0], b[1]),
    );
    let (ab, ap) = (b - a, p - a);
    let t = if ab.length_squared() == 0.0 {
        0.0
    } else {
        (ap.dot(ab) / ab.length_squared()).clamp(0.0, 1.0)
    };
    (p - (a + ab * t)).length()
}

/// A kept surface the ray through a sample crosses.
#[derive(Debug, Clone, Copy)]
struct Hit {
    mesh: usize,
    cap: bool,
    see_through: bool,
    depth: f64,
}

/// Every kept surface the ray through `q` crosses, front to back, skipping
/// the triangles of `exclude` (an edge's own triangles) and those wholly
/// behind the depth `behind`, which cannot be in front of the point there;
/// `Err` when the sample is ambiguous.
fn hits(
    projected: &Projected,
    q: [f64; 2],
    exclude: Option<(usize, [u32; 2])>,
    behind: f64,
) -> Result<Vec<Hit>, ()> {
    let mut out = Vec::new();
    for (mesh, (v, tris, glass)) in projected.meshes.iter().enumerate() {
        let mut beyond = 0usize;
        for &t in tris {
            if exclude.is_some_and(|(m, [a, b])| m == mesh && t.contains(&a) && t.contains(&b))
                || t.iter().all(|&i| v[i as usize][2] < behind - EPS)
            {
                continue;
            }
            let Some(depth) = crossing(projected, v, t, q)? else {
                continue;
            };
            if let Some(near) = projected.near {
                if (depth - near).abs() < EPS {
                    return Err(());
                }
                if depth > near {
                    beyond += 1;
                    continue;
                }
            }
            out.push(Hit {
                mesh,
                cap: false,
                see_through: *glass,
                depth,
            });
        }
        // The cap: the point of the near plane over `q` is inside the solid
        // when an odd number of its triangles lie beyond the plane there.
        if projected.cut[mesh] && beyond % 2 == 1 {
            out.push(Hit {
                mesh,
                cap: true,
                see_through: false,
                depth: projected.near.expect("a cut has a near plane"),
            });
        }
    }
    out.sort_by(|a, b| b.depth.total_cmp(&a.depth));
    Ok(out)
}

#[derive(Default, Debug)]
struct Tally {
    checked: usize,
    skipped: usize,
    counted: usize,
}

// --- hidden lines -------------------------------------------------------------

fn check_lines(
    boxes: &[TriMesh],
    glass: &[bool],
    view: View,
    outcome: &HiddenLineOutcome,
    meshes: &[HiddenLineMesh<'_>],
    context: &str,
    tally: &mut Tally,
) {
    let projected = Projected::new(boxes, glass, &view);
    for segment in &outcome.segments {
        let mesh = meshes
            .iter()
            .position(|m| m.id == segment.source.mesh)
            .expect("a request mesh");
        let [a, b] = segment
            .source
            .vertices
            .map(|i| projected.meshes[mesh].0[i as usize]);
        let (start, end) = (segment.start, segment.end);
        for fraction in [0.25, 0.5, 0.75] {
            let q = [
                start.x + fraction * (end.x - start.x),
                start.y + fraction * (end.y - start.y),
            ];
            // The edge's point over `q`, from the image axis it moves along
            // most.
            let (ia, ib) = (projected.image(a), projected.image(b));
            let axis = usize::from((ib[1] - ia[1]).abs() > (ib[0] - ia[0]).abs());
            let s = if projected.perspective {
                let (d, dw) = (b[axis] - a[axis], b[2] - a[2]);
                (-q[axis] * a[2] - a[axis]) / (d + q[axis] * dw)
            } else {
                (q[axis] - ia[axis]) / (ib[axis] - ia[axis])
            };
            let depth = a[2] + s * (b[2] - a[2]);
            let Ok(hits) = hits(&projected, q, Some((mesh, segment.source.vertices)), depth) else {
                tally.skipped += 1;
                continue;
            };
            if hits.iter().any(|h| (h.depth - depth).abs() < EPS) {
                tally.skipped += 1;
                continue;
            }
            let front: Vec<&Hit> = hits.iter().filter(|h| h.depth > depth).collect();
            let hidden = front.iter().any(|h| !h.see_through);
            let count = if hidden { 0 } else { front.len() as u32 };
            tally.checked += 1;
            tally.counted += usize::from(count > 0);
            assert_eq!(
                (segment.visibility, segment.seen_through),
                (
                    if hidden {
                        Visibility::Hidden
                    } else {
                        Visibility::Visible
                    },
                    count
                ),
                "{context}: {segment:?} at {fraction}"
            );
        }
    }
}

fn flags_of(boxes: &[TriMesh], glass: &[bool]) -> Vec<Vec<Occlusion>> {
    boxes
        .iter()
        .zip(glass)
        .map(|(m, &g)| {
            let flag = if g {
                Occlusion::SeeThrough
            } else {
                Occlusion::Opaque
            };
            vec![flag; m.triangle_count()]
        })
        .collect()
}

/// A near plane through the middle of one seeded box.
fn mid_depth(seed: u64, boxes: &[TriMesh], frame: &Frame3) -> f64 {
    let mut rng = Lcg(seed.wrapping_mul(0xD1B5_4A32_D192_ED03) | 1);
    let mesh = &boxes[rng.below(boxes.len() as u64) as usize];
    let depths = mesh
        .positions
        .iter()
        .map(|p| (*p - frame.origin).dot(frame.z));
    let lo = depths.clone().fold(f64::INFINITY, f64::min);
    let hi = depths.fold(f64::NEG_INFINITY, f64::max);
    (lo + hi) / 2.0
}

fn pinhole() -> (&'static str, View) {
    // An eye above and beside the scene, looking at its middle.
    let eye = Point3::new(-1.3, -1.7, 3.1);
    let target = Point3::new(1.0, 1.0, 0.5);
    let z = (eye - target).normalize();
    let x = Vec3::Z.cross(z).normalize();
    let y = z.cross(x);
    let frame = Frame3 {
        origin: eye,
        x,
        y,
        z,
    };
    (
        "pinhole",
        View::perspective(frame, -0.01).expect("a pinhole"),
    )
}

#[test]
fn hidden_lines_agree_with_see_through_ray_casting() {
    let mut report = Vec::new();
    let mut cases: Vec<(&str, Frame3, bool)> = Vec::new();
    for (name, frame) in views() {
        cases.push((name, frame, false));
        cases.push((name, frame, true));
    }
    for (name, frame, with_near) in cases {
        let mut tally = Tally::default();
        for seed in 1..=40 {
            let boxes = scene(seed);
            let glass = glazed(seed, boxes.len());
            let flags = flags_of(&boxes, &glass);
            let meshes: Vec<HiddenLineMesh<'_>> = solids(&boxes)
                .into_iter()
                .zip(&flags)
                .map(|(m, f)| m.with_occlusion(f))
                .collect();
            let mut view = View::new(frame);
            if with_near {
                view = view.with_near(mid_depth(seed, &boxes, &frame));
            }
            let request =
                HiddenLineRequest::new(&meshes, view, EdgeSelection::new(Some(0.5)), LIMITS);
            let outcome = ScalarHiddenLines
                .hidden_lines(&request, &options())
                .unwrap_or_else(|error| panic!("seed {seed}, {name}: {error}"));
            let context = format!("seed {seed}, {name}, near {with_near}");
            check_lines(
                &boxes, &glass, view, &outcome, &meshes, &context, &mut tally,
            );
        }
        report.push((name, with_near, tally));
    }
    let (name, view) = pinhole();
    let mut tally = Tally::default();
    for seed in 1..=40 {
        let boxes = scene(seed);
        let glass = glazed(seed, boxes.len());
        let flags = flags_of(&boxes, &glass);
        let meshes: Vec<HiddenLineMesh<'_>> = solids(&boxes)
            .into_iter()
            .zip(&flags)
            .map(|(m, f)| m.with_occlusion(f))
            .collect();
        let request = HiddenLineRequest::new(&meshes, view, EdgeSelection::new(Some(0.5)), LIMITS);
        let outcome = ScalarHiddenLines
            .hidden_lines(&request, &options())
            .unwrap_or_else(|error| panic!("seed {seed}, {name}: {error}"));
        let context = format!("seed {seed}, {name}");
        check_lines(
            &boxes, &glass, view, &outcome, &meshes, &context, &mut tally,
        );
    }
    report.push((name, false, tally));
    for (name, near, tally) in &report {
        eprintln!("lines {name:>9} near {near}: {tally:?}");
        let total = tally.checked + tally.skipped;
        assert!(
            tally.checked >= 2000 && tally.counted >= 200,
            "{name}, near {near}: {tally:?}"
        );
        // The share each view must check; see the module docs for why.
        let percent = match *name {
            "oblique" | "pinhole" => 85,
            "isometric" => 75,
            _ => 0,
        };
        assert!(
            tally.checked * 100 >= total * percent,
            "{name}, near {near}: {tally:?}"
        );
    }
}

// --- visible faces ------------------------------------------------------------

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

fn near_a_ring_edge(outcome: &VisibleFaceOutcome, p: [f64; 2]) -> bool {
    outcome.regions.iter().any(|r| {
        std::iter::once(&r.outer).chain(&r.holes).any(|ring| {
            (0..ring.len()).any(|i| {
                let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                distance(p, [a.x, a.y], [b.x, b.y]) < EPS
            })
        })
    })
}

const G0: FaceSource = FaceSource::Group(FaceGroup(0));

#[test]
fn visible_faces_agree_with_see_through_ray_casting() {
    let mut report = Vec::new();
    for (name, frame) in views() {
        for with_near in [false, true] {
            let mut tally = Tally::default();
            for seed in 0..12 {
                let boxes = scene(seed);
                let glass = glazed(seed, boxes.len());
                let flags = flags_of(&boxes, &glass);
                let meshes: Vec<HiddenLineMesh<'_>> = solids(&boxes)
                    .into_iter()
                    .zip(&flags)
                    .map(|(m, f)| m.with_occlusion(f))
                    .collect();
                let mut view = View::new(frame);
                if with_near {
                    view = view.with_near(mid_depth(seed, &boxes, &frame));
                }
                let request = VisibleFaceRequest::new(
                    &meshes,
                    view,
                    FACE_LIMITS.with_max_behind_surfaces(1_000_000),
                );
                let outcome = ScalarVisibleFaces
                    .visible_faces(&request, &options())
                    .unwrap_or_else(|error| panic!("seed {seed}, {name}: {error}"));
                assert!(
                    outcome.evidence.snap_displacement.get() * 4.0 < EPS,
                    "seed {seed}, {name}"
                );
                assert_valid_topology(&outcome);
                let context = format!("seed {seed}, {name}, near {with_near}");
                check_faces(
                    &boxes, &glass, view, &outcome, &meshes, &context, &mut tally,
                );
            }
            report.push((name, with_near, tally));
        }
    }
    for (name, near, tally) in &report {
        eprintln!("faces {name:>9} near {near}: {tally:?}");
        let total = tally.checked + tally.skipped;
        assert!(tally.counted >= 1000, "{name}, near {near}: {tally:?}");
        assert!(
            tally.checked * 10 >= total * 9,
            "{name}, near {near}: {tally:?}"
        );
    }
}

fn check_faces(
    boxes: &[TriMesh],
    glass: &[bool],
    view: View,
    outcome: &VisibleFaceOutcome,
    meshes: &[HiddenLineMesh<'_>],
    context: &str,
    tally: &mut Tally,
) {
    let projected = Projected::new(boxes, glass, &view);
    let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
    for p in projected.meshes.iter().flat_map(|(v, _, _)| v) {
        for axis in 0..2 {
            lo[axis] = lo[axis].min(p[axis]);
            hi[axis] = hi[axis].max(p[axis]);
        }
    }
    for i in 0..GRID {
        for j in 0..GRID {
            // Irrational strides: no sample sits on a grid-aligned box edge.
            let q = [
                lo[0] - 0.025 * (hi[0] - lo[0])
                    + 1.05 * (hi[0] - lo[0]) * (i as f64 + 0.4142) / GRID as f64,
                lo[1] - 0.025 * (hi[1] - lo[1])
                    + 1.05 * (hi[1] - lo[1]) * (j as f64 + 0.7320) / GRID as f64,
            ];
            let Ok(hits) = hits(&projected, q, None, f64::NEG_INFINITY) else {
                tally.skipped += 1;
                continue;
            };
            if hits.windows(2).any(|w| w[0].depth - w[1].depth < EPS)
                || near_a_ring_edge(outcome, q)
            {
                tally.skipped += 1;
                continue;
            }
            tally.checked += 1;
            let p = Point2::new(q[0], q[1]);
            let holding: Vec<_> = outcome
                .regions
                .iter()
                .filter(|r| inside(&r.outer, p) && !r.holes.iter().any(|h| inside(h, p)))
                .collect();
            let Some(front) = hits.first() else {
                assert!(holding.is_empty(), "{context}: {q:?}");
                continue;
            };
            let source_of = |hit: &Hit| if hit.cap { FaceSource::Cap } else { G0 };
            let id = meshes[front.mesh].id;
            // No tie at the sample, so exactly one region holds it, and it is
            // the frontmost source's.
            assert_eq!(holding.len(), 1, "{context}: {q:?} in {holding:?}");
            let own = &holding;
            assert!(
                own[0].mesh == id && own[0].source == source_of(front),
                "{context}: {q:?} in {holding:?}"
            );
            let stop = hits
                .iter()
                .position(|h| !h.see_through)
                .map_or(hits.len(), |at| at + 1);
            let expected: Vec<BehindSurface> = if front.see_through {
                hits[1..stop]
                    .iter()
                    .enumerate()
                    .map(|(layer, h)| {
                        BehindSurface::new(meshes[h.mesh].id, source_of(h), layer as u32 + 1)
                    })
                    .collect()
            } else {
                Vec::new()
            };
            tally.counted += usize::from(!expected.is_empty());
            assert_eq!(own[0].behind, expected, "{context}: {q:?}");
        }
    }
}
