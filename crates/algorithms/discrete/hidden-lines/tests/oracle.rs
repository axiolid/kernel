//! An independent oracle for the provider: plain binary64 ray casting.
//!
//! It shares no code with the provider. For each output segment it samples
//! three points, lifts each to the source edge in 3D and asks, with ordinary
//! floating point, whether a kept triangle lies above it in the view. A
//! sample within `1e-9` of the scene size of a tie (an image boundary, an
//! equal depth) is skipped as ambiguous; every other sample must agree with
//! the provider. Coverage is checked per edge: the visible and hidden pieces
//! add up to the edge's image.
//!
//! A near plane removes what lies beyond it and caps the solids it cuts. The
//! oracle counts, with an upward ray from the sampled point's image at the
//! near depth, the crossings of each cut solid's triangles above the plane;
//! an odd count is inside the solid, so the cap covers the sample.
//!
//! Scenes are seeded boxes on a 1/8 grid, some rotated (by generic angles,
//! so no face is edge-on only up to rounding noise), in the plan view,
//! the four elevations, the isometric view and one oblique view.

mod common;

use std::cell::Cell;
use std::collections::BTreeMap;
use std::time::Instant;

use axiolid_contracts::GeomError;
use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{OrthographicView, Visibility};
use common::{
    box_mesh, draw, frame, isometric_frame, looking_along, plan_frame, scalar, solids, try_draw,
    BOX_EDGES,
};

type V3 = [f64; 3];

struct Lcg(u64);

impl Lcg {
    fn below(&mut self, n: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) % n
    }
}

/// 1 to 12 boxes on a 1/8 grid; a third of them rotated about their centre.
fn scene(seed: u64) -> Vec<TriMesh> {
    let mut rng = Lcg(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
    let count = 1 + rng.below(12) as usize;
    let about_z: [f64; 6] = [0.31, 0.57, 0.93, 1.21, 2.17, 2.6];
    let about_x: [f64; 3] = [0.2, 0.4, 0.7];
    let about_y: [f64; 3] = [0.15, 0.35, 0.55];
    (0..count)
        .map(|_| {
            let min = [
                rng.below(13) as f64 / 8.0,
                rng.below(13) as f64 / 8.0,
                rng.below(9) as f64 / 8.0,
            ];
            let size = [
                (1 + rng.below(8)) as f64 / 8.0,
                (1 + rng.below(8)) as f64 / 8.0,
                (1 + rng.below(8)) as f64 / 8.0,
            ];
            let max = [min[0] + size[0], min[1] + size[1], min[2] + size[2]];
            let mut mesh = box_mesh(min, max);
            if rng.below(3) == 0 {
                let centre = Point3::new(
                    (min[0] + max[0]) / 2.0,
                    (min[1] + max[1]) / 2.0,
                    (min[2] + max[2]) / 2.0,
                );
                let (sz, cz) = about_z[rng.below(6) as usize].sin_cos();
                let (sx, cx) = about_x[rng.below(3) as usize].sin_cos();
                let (sy, cy) = about_y[rng.below(3) as usize].sin_cos();
                for p in &mut mesh.positions {
                    let d = *p - centre;
                    // Generic angles about all three axes: no face is left
                    // parallel to an axis view, so none is edge-on up to noise.
                    let d = Vec3::new(cy * d.x + sy * d.z, d.y, -sy * d.x + cy * d.z);
                    let tilted = Vec3::new(d.x, cx * d.y - sx * d.z, sx * d.y + cx * d.z);
                    *p = centre
                        + Vec3::new(
                            cz * tilted.x - sz * tilted.y,
                            sz * tilted.x + cz * tilted.y,
                            tilted.z,
                        );
                }
            }
            mesh
        })
        .collect()
}

fn views() -> Vec<(&'static str, Frame3)> {
    let o = Point3::ZERO;
    vec![
        ("plan", plan_frame()),
        ("+X", frame(Vec3::Y, Vec3::Z, Vec3::X, o)),
        ("-X", frame(Vec3::Z, Vec3::Y, -Vec3::X, o)),
        ("+Y", frame(Vec3::Z, Vec3::X, Vec3::Y, o)),
        ("-Y", frame(Vec3::X, Vec3::Z, -Vec3::Y, o)),
        ("isometric", isometric_frame(Point3::new(0.1, 0.2, 0.3))),
        (
            "oblique",
            looking_along(Vec3::new(0.3, -0.5, 0.8), Point3::new(-0.4, 0.7, 0.1)),
        ),
    ]
}

// --- the oracle ------------------------------------------------------------

struct Oracle {
    /// Per mesh: view-space vertices.
    verts: Vec<Vec<V3>>,
    /// Per mesh: triangles.
    tris: Vec<Vec<[usize; 3]>>,
    margin: f64,
    flat: f64,
    near: Option<f64>,
    /// Per mesh: whether the near plane cuts it.
    cut: Vec<bool>,
    /// Whether the last verdict was a hidden one decided by a cap alone.
    by_cap: Cell<bool>,
}

/// The sampled point: image point `q` of the edge `a`-`b` with end vertices
/// `ends`.
struct Sample {
    a: usize,
    b: usize,
    ends: (V3, V3),
    q: [f64; 2],
}

#[derive(Debug, PartialEq)]
enum Verdict {
    Hidden,
    Visible,
    Ambiguous,
}

enum Cover {
    Out,
    In,
    Edge,
}

impl Oracle {
    fn new(meshes: &[TriMesh], frame: &Frame3, near: Option<f64>) -> Self {
        let verts: Vec<Vec<V3>> = meshes
            .iter()
            .map(|m| {
                m.positions
                    .iter()
                    .map(|p| {
                        let d = *p - frame.origin;
                        [d.dot(frame.x), d.dot(frame.y), d.dot(frame.z)]
                    })
                    .collect()
            })
            .collect();
        let tris = meshes
            .iter()
            .map(|m| {
                m.indices
                    .as_chunks::<3>()
                    .0
                    .iter()
                    .map(|t| t.map(|i| i as usize))
                    .collect()
            })
            .collect();
        let size = verts
            .iter()
            .flatten()
            .flat_map(|v| v.iter().map(|c| c.abs()))
            .fold(1.0f64, f64::max);
        let cut = verts
            .iter()
            .map(|vs| {
                near.is_some_and(|n| vs.iter().any(|v| v[2] > n) && vs.iter().any(|v| v[2] <= n))
            })
            .collect();
        Self {
            verts,
            tris,
            margin: 1e-9 * size,
            flat: 1e-12 * size * size,
            near,
            cut,
            by_cap: Cell::new(false),
        }
    }

    /// Hidden, visible or too close to a tie, at image point `q` of the edge
    /// `a`-`b` of mesh `mesh`.
    fn judge(&self, mesh: usize, a: usize, b: usize, q: [f64; 2]) -> Verdict {
        let m = self.margin;
        let (va, vb) = (self.verts[mesh][a], self.verts[mesh][b]);
        let d = [vb[0] - va[0], vb[1] - va[1]];
        let t = ((q[0] - va[0]) * d[0] + (q[1] - va[1]) * d[1]) / (d[0] * d[0] + d[1] * d[1]);
        let wp = va[2] + t * (vb[2] - va[2]);
        let sample = Sample {
            a,
            b,
            ends: (va, vb),
            q,
        };

        // The verdict for a candidate occluder point with depth `w` found
        // at cover class `cover`: Some(true) hides, Some(false) cannot.
        // None is a tie.
        let rule = |cover: Cover, w: f64| -> Option<bool> {
            match cover {
                Cover::Out => Some(false),
                Cover::In if w - wp > m => Some(true),
                Cover::In if w - wp < -m => Some(false),
                Cover::In => None,
                Cover::Edge if w - wp > -m => None,
                Cover::Edge => Some(false),
            }
        };

        let mut hidden = false;
        for (mi, triangles) in self.tris.iter().enumerate() {
            for tri in triangles {
                let own = mi == mesh && tri.contains(&a) && tri.contains(&b);
                let v = tri.map(|i| self.verts[mi][i]);
                let area2 = (v[1][0] - v[0][0]) * (v[2][1] - v[0][1])
                    - (v[1][1] - v[0][1]) * (v[2][0] - v[0][0]);
                let outcome = if area2.abs() > self.flat {
                    if own {
                        continue;
                    }
                    self.solid(&v, area2, sample.q, &rule)
                } else {
                    self.edge_on(tri, &v, own, &sample, &rule)
                };
                match outcome {
                    None => return Verdict::Ambiguous,
                    Some(true) => hidden = true,
                    Some(false) => {}
                }
            }
        }
        let by_triangles = hidden;
        if let Some(near) = self.near {
            for (mi, cut) in self.cut.iter().enumerate() {
                if !cut {
                    continue;
                }
                match self.cap_parity(mi, near, q) {
                    None => return Verdict::Ambiguous,
                    Some(false) => {}
                    // Equal depth does not occlude, and a sample at the
                    // plane is a tie.
                    Some(true) if wp < near - m => hidden = true,
                    Some(true) => return Verdict::Ambiguous,
                }
            }
        }
        self.by_cap.set(hidden && !by_triangles);
        if hidden {
            Verdict::Hidden
        } else {
            Verdict::Visible
        }
    }

    /// Whether `(q, near)` is inside mesh `mi`: the parity of the crossings
    /// of an upward ray with its triangles above the near plane. `None` when
    /// the ray grazes an image edge or a crossing is at the plane.
    fn cap_parity(&self, mi: usize, near: f64, q: [f64; 2]) -> Option<bool> {
        let m = self.margin;
        let mut crossings = 0;
        for tri in &self.tris[mi] {
            let v = tri.map(|i| self.verts[mi][i]);
            let area2 = (v[1][0] - v[0][0]) * (v[2][1] - v[0][1])
                - (v[1][1] - v[0][1]) * (v[2][0] - v[0][0]);
            if area2.abs() <= self.flat {
                continue;
            }
            let (cover, w) = self.cover_and_depth(&v, area2, q);
            match cover {
                Cover::Out => {}
                Cover::In if w > near + m => crossings += 1,
                Cover::In if w < near - m => {}
                Cover::Edge if w < near - m => {}
                Cover::In | Cover::Edge => return None,
            }
        }
        Some(crossings % 2 == 1)
    }

    /// A triangle with an image area: coverage by signed distances to its
    /// three image edges, depth by its plane.
    fn solid(
        &self,
        v: &[V3; 3],
        area2: f64,
        q: [f64; 2],
        rule: &dyn Fn(Cover, f64) -> Option<bool>,
    ) -> Option<bool> {
        let m = self.margin;
        let (cover, w) = self.cover_and_depth(v, area2, q);
        match (rule(cover, w), self.near) {
            // A point beyond the near plane is removed; one at it is a tie.
            (Some(true), Some(near)) if w > near + m => Some(false),
            (Some(true), Some(near)) if w >= near - m => None,
            (verdict, _) => verdict,
        }
    }

    /// Where `q` is against a triangle's image, and the depth of its plane
    /// there.
    fn cover_and_depth(&self, v: &[V3; 3], area2: f64, q: [f64; 2]) -> (Cover, f64) {
        let orient = |p: &V3, r: &V3, s: [f64; 2]| {
            (r[0] - p[0]) * (s[1] - p[1]) - (r[1] - p[1]) * (s[0] - p[0])
        };
        let mut nearest = f64::INFINITY;
        for k in 0..3 {
            let (p, r) = (&v[k], &v[(k + 1) % 3]);
            let length = ((r[0] - p[0]).powi(2) + (r[1] - p[1]).powi(2)).sqrt();
            nearest = nearest.min(orient(p, r, q) * area2.signum() / length);
        }
        let cover = if nearest < -self.margin {
            Cover::Out
        } else if nearest > self.margin {
            Cover::In
        } else {
            Cover::Edge
        };
        let l0 = orient(&v[1], &v[2], q) / orient(&v[1], &v[2], [v[0][0], v[0][1]]);
        let l1 = orient(&v[2], &v[0], q) / orient(&v[2], &v[0], [v[1][0], v[1][1]]);
        let l2 = 1.0 - l0 - l1;
        (cover, l0 * v[0][2] + l1 * v[1][2] + l2 * v[2][2])
    }

    /// A triangle seen edge-on: its edges along the sampled edge's line, each
    /// covering the image of its own extent at its own depth.
    fn edge_on(
        &self,
        tri: &[usize; 3],
        v: &[V3; 3],
        own: bool,
        sample: &Sample,
        rule: &dyn Fn(Cover, f64) -> Option<bool>,
    ) -> Option<bool> {
        let m = self.margin;
        let (a, b, q) = (sample.a, sample.b, sample.q);
        let (va, vb) = sample.ends;
        let mut any = false;
        // The shallowest depth the fibre reaches over `q`.
        let mut lowest = f64::INFINITY;
        for k in 0..3 {
            let (i, j) = (k, (k + 1) % 3);
            if own && [tri[i], tri[j]].iter().all(|x| *x == a || *x == b) {
                continue;
            }
            let (c, e) = (v[i], v[j]);
            let f = [e[0] - c[0], e[1] - c[1]];
            let length = (f[0] * f[0] + f[1] * f[1]).sqrt();
            if length <= m {
                continue;
            }
            let distance = |p: V3| ((p[0] - c[0]) * f[1] - (p[1] - c[1]) * f[0]).abs() / length;
            if distance(va) > m || distance(vb) > m {
                continue;
            }
            let s = ((q[0] - c[0]) * f[0] + (q[1] - c[1]) * f[1]) / (length * length);
            let slack = m / length;
            let cover = if s < -slack || s > 1.0 + slack {
                Cover::Out
            } else if s > slack && s < 1.0 - slack {
                Cover::In
            } else {
                Cover::Edge
            };
            let w = c[2] + s * (e[2] - c[2]);
            if !matches!(cover, Cover::Out) {
                lowest = lowest.min(w);
            }
            match rule(cover, w) {
                None => return None,
                Some(true) => any = true,
                Some(false) => {}
            }
        }
        // The fibre hides only through points the near plane keeps.
        match self.near {
            Some(near) if any && lowest >= near - m && lowest <= near + m => None,
            Some(near) if lowest > near + m => Some(false),
            _ => Some(any),
        }
    }
}

// --- the comparison --------------------------------------------------------

#[derive(Default)]
struct Tally {
    scenes: usize,
    segments: usize,
    robust: usize,
    ambiguous: usize,
    /// Robust samples hidden by a cap alone.
    capped: usize,
    /// Scenes the provider refused as a tie.
    refused: usize,
    seconds: f64,
}

/// The depth halfway through the extent of one of the scene's boxes, so
/// that the near plane cuts it.
fn mid_depth(seed: u64, meshes: &[TriMesh], view_frame: &Frame3) -> f64 {
    let mut rng = Lcg(seed.wrapping_mul(0xD1B5_4A32_D192_ED03) | 1);
    let mesh = &meshes[rng.below(meshes.len() as u64) as usize];
    let depths = mesh
        .positions
        .iter()
        .map(|p| (*p - view_frame.origin).dot(view_frame.z));
    let lo = depths.clone().fold(f64::INFINITY, f64::min);
    let hi = depths.fold(f64::NEG_INFINITY, f64::max);
    (lo + hi) / 2.0
}

fn check(seed: u64, name: &str, view_frame: Frame3, with_near: bool, tally: &mut Tally) {
    let meshes = scene(seed);
    let near = with_near.then(|| mid_depth(seed, &meshes, &view_frame));
    let oracle = Oracle::new(&meshes, &view_frame, near);
    let input = solids(&meshes);
    let mut view = OrthographicView::new(view_frame);
    if let Some(near) = near {
        view = view.with_near(near);
    }
    let started = Instant::now();
    let context = format!("scene {seed}, view {name}, near {near:?}");
    let outcome = if near.is_some() {
        match try_draw(&scalar(), &input, view) {
            Ok(outcome) => outcome,
            // A face in the near plane, or a cut line under a drawn edge, is
            // a tie the provider refuses.
            Err(error) if matches!(error.error, GeomError::Degenerate(_)) => {
                assert!(error.mesh.is_some(), "{context}: {error}");
                tally.refused += 1;
                tally.scenes += 1;
                return;
            }
            Err(error) => panic!("{context}: {error}"),
        }
    } else {
        draw(&scalar(), &input, view)
    };
    tally.seconds += started.elapsed().as_secs_f64();
    tally.scenes += 1;

    // Coverage: pieces of one edge add up to its image; only box edges, and
    // exactly those whose image is not a point, are drawn.
    let mut covered: BTreeMap<(u64, [u32; 2]), f64> = BTreeMap::new();
    for segment in &outcome.segments {
        let key = (segment.source.mesh.0, segment.source.vertices);
        *covered.entry(key).or_default() += (segment.end - segment.start).length();
    }
    for (index, mesh) in meshes.iter().enumerate() {
        for edge in BOX_EDGES {
            let image = |i: u32| {
                let d = mesh.positions[i as usize] - view_frame.origin;
                [d.dot(view_frame.x), d.dot(view_frame.y)]
            };
            let (p, r) = (image(edge[0]), image(edge[1]));
            let mut length = ((r[0] - p[0]).powi(2) + (r[1] - p[1]).powi(2)).sqrt();
            if let Some(near) = near {
                // Only the part of the edge not beyond the near plane.
                let depth =
                    |i: u32| (mesh.positions[i as usize] - view_frame.origin).dot(view_frame.z);
                let (wa, wb) = (depth(edge[0]), depth(edge[1]));
                let at = (near - wa) / (wb - wa);
                length *= match (wa <= near, wb <= near) {
                    (true, true) => 1.0,
                    (false, false) => 0.0,
                    (true, false) => at,
                    (false, true) => 1.0 - at,
                };
            }
            let got = covered.remove(&(index as u64 + 1, edge));
            if length > 1e-9 {
                let got = got.unwrap_or_else(|| panic!("{context}: edge {edge:?} not drawn"));
                assert!(
                    (got - length).abs() < 1e-9,
                    "{context}: edge {edge:?} covers {got} of {length}"
                );
            } else if near.is_some() {
                // Wholly beyond the plane, or a point in the image.
            } else {
                assert!(got.is_none(), "{context}: point edge {edge:?} drawn");
            }
        }
    }
    assert!(
        covered.is_empty(),
        "{context}: unexpected edges {covered:?}"
    );

    for segment in &outcome.segments {
        let (mesh, [a, b]) = (segment.source.mesh.0 as usize - 1, segment.source.vertices);
        for fraction in [0.25, 0.5, 0.75] {
            tally.segments += 1;
            let q = segment.start + (segment.end - segment.start) * fraction;
            let verdict = oracle.judge(mesh, a as usize, b as usize, [q.x, q.y]);
            let expected = match segment.visibility {
                Visibility::Visible => Verdict::Visible,
                Visibility::Hidden => Verdict::Hidden,
                _ => unreachable!("unknown visibility"),
            };
            match verdict {
                Verdict::Ambiguous => tally.ambiguous += 1,
                verdict => {
                    tally.robust += 1;
                    if oracle.by_cap.get() {
                        tally.capped += 1;
                    }
                    assert_eq!(verdict, expected, "{context}: {segment:?} at {fraction}");
                }
            }
        }
    }
}

#[test]
fn the_provider_agrees_with_an_independent_ray_cast() {
    const SEEDS: u64 = 40;
    let mut report = Vec::new();
    for (name, view_frame) in views() {
        let mut tally = Tally::default();
        for seed in 1..=SEEDS {
            check(seed, name, view_frame, false, &mut tally);
        }
        eprintln!(
            "{name:>9}: {} scenes, {} samples, {} robust, {} ambiguous, {:.3} s",
            tally.scenes, tally.segments, tally.robust, tally.ambiguous, tally.seconds
        );
        report.push((name, tally));
    }
    for (name, tally) in &report {
        assert!(
            tally.robust >= 100,
            "{name}: only {} robust samples",
            tally.robust
        );
    }
    for (name, tally) in report
        .iter()
        .filter(|(n, _)| ["isometric", "oblique"].contains(n))
    {
        assert!(
            tally.robust * 2 >= tally.segments,
            "{name}: too many ties: {} of {}",
            tally.ambiguous,
            tally.segments
        );
    }
}

#[test]
fn caps_of_cut_solids_agree_with_an_independent_ray_parity() {
    const SEEDS: u64 = 40;
    let (mut capped, mut scenes, mut refused) = (0, 0, 0);
    for (name, view_frame) in views() {
        let mut tally = Tally::default();
        for seed in 1..=SEEDS {
            check(seed, name, view_frame, true, &mut tally);
        }
        eprintln!(
            "{name:>9}: {} cut scenes ({} refused), {} samples, {} robust, {} by a cap, \
             {} ambiguous, {:.3} s",
            tally.scenes,
            tally.refused,
            tally.segments,
            tally.robust,
            tally.capped,
            tally.ambiguous,
            tally.seconds
        );
        assert!(
            tally.robust >= 50,
            "{name}: only {} robust samples",
            tally.robust
        );
        capped += tally.capped;
        scenes += tally.scenes;
        refused += tally.refused;
    }
    // The caps decide a good share of the samples, and ties are rare.
    assert!(capped >= 200, "only {capped} samples are hidden by a cap");
    assert!(
        refused * 10 <= scenes,
        "{refused} of {scenes} scenes refused"
    );
}

#[test]
fn a_cut_through_a_slab_at_a_rail_top_does_not_hide_the_rail() {
    // The slab [0,4]^2 x [0,2] is cut at depth 1, the top face of a rail
    // [3,6] x [1,1.5] x [0,1] that lies partly under the slab's image. The
    // rail is not cut: its top face is in the plane, at equal depth with the
    // slab's cap, which does not occlude it.
    let slab = box_mesh([0.0, 0.0, 0.0], [4.0, 4.0, 2.0]);
    let rail = box_mesh([3.0, 1.0, 0.0], [6.0, 1.5, 1.0]);
    let meshes = [slab, rail];
    let input = solids(&meshes);
    let outcome = draw(
        &scalar(),
        &input,
        OrthographicView::new(plan_frame()).with_near(1.0),
    );
    let top = |visibility| {
        outcome
            .segments
            .iter()
            .filter(|s| s.source.mesh.0 == 2 && s.source.vertices[0] >= 4)
            .filter(|s| s.visibility == visibility)
            .map(|s| (s.end - s.start).length())
            .sum::<f64>()
    };
    assert!((top(Visibility::Visible) - 7.0).abs() < 1e-12);
    assert_eq!(top(Visibility::Hidden), 0.0);
    // The slab's own bottom ring is under its cap.
    let slab_visible = outcome
        .segments
        .iter()
        .filter(|s| s.source.mesh.0 == 1 && s.visibility == Visibility::Visible)
        .count();
    assert_eq!(slab_visible, 0);
}
