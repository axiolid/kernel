//! An independent oracle for perspective hidden lines: binary64 rays cast
//! from the eye.
//!
//! It shares no code with the provider. For each output segment it samples
//! three image points, takes the point of the source edge on the eye ray
//! through each, and asks, with ordinary floating point, whether a kept
//! triangle is hit strictly nearer the eye. A sample within `1e-9` of the
//! scene size of a tie (a ray through a triangle's boundary, equal distance,
//! the near plane) is skipped as ambiguous; every other sample must agree
//! with the provider. Coverage is checked per edge: the pieces add up to the
//! image of the part of the edge in front of the near plane.
//!
//! A near plane removes what lies beyond it and caps the solids it cuts. The
//! oracle decides a cap by the parity of the crossings of the whole eye line
//! on the eye's side of the plane, the opposite direction to the provider's.
//!
//! Scenes are the seeded boxes of the orthographic oracle, seen from eyes
//! outside the scene and from eyes inside it, so that triangles lie behind
//! the eye and cross the near plane.

mod common;

use std::collections::BTreeMap;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Aabb2, Frame3, Point2, Point3, Tolerance, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    EdgeSelection, HiddenLineRequest, MeshHiddenLines, View, Visibility,
};
use common::{frame, looking_along, scalar, scene, solids, BOX_EDGES};

type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

struct Oracle {
    verts: Vec<Vec<V3>>,
    tris: Vec<Vec<[usize; 3]>>,
    margin: f64,
    near: f64,
    cut: Vec<bool>,
}

/// Where the eye line `mu * d` meets a triangle: `mu`, and the smallest
/// barycentric coordinate. `None` when the line is parallel to its plane.
fn line_hit(d: V3, v: [V3; 3]) -> Option<(f64, f64)> {
    let (e1, e2) = (sub(v[1], v[0]), sub(v[2], v[0]));
    let p = cross(d, e2);
    let det = dot(e1, p);
    let scale = dot(e1, e1).sqrt() * dot(e2, e2).sqrt() * dot(d, d).sqrt();
    if det.abs() <= 1e-12 * scale {
        return None;
    }
    // The origin is the eye: s = -v0.
    let s = [-v[0][0], -v[0][1], -v[0][2]];
    let u = dot(s, p) / det;
    let q = cross(s, e1);
    let w = dot(d, q) / det;
    let mu = dot(e2, q) / det;
    Some((mu, u.min(w).min(1.0 - u - w)))
}

impl Oracle {
    fn new(meshes: &[TriMesh], frame: &Frame3, near: f64) -> Self {
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
            .map(|vs| vs.iter().any(|v| v[2] > near) && vs.iter().any(|v| v[2] <= near))
            .collect();
        Self {
            verts,
            tris,
            margin: 1e-9 * size,
            near,
            cut,
        }
    }

    /// The point of the edge `a`-`b` on the eye ray through the image point
    /// `q` (at depth 1 in front of the eye), as its distance `lambda` along
    /// the ray `(q, -1)`.
    fn lambda(&self, mesh: usize, a: usize, b: usize, q: [f64; 2]) -> f64 {
        let (va, vb) = (self.verts[mesh][a], self.verts[mesh][b]);
        // U_A + t dU = q_x (-W_A - t dW), and likewise for V.
        let coefficient = |axis: usize| {
            (
                (vb[axis] - va[axis]) + q[axis] * (vb[2] - va[2]),
                -q[axis] * va[2] - va[axis],
            )
        };
        let (cx, rx) = coefficient(0);
        let (cy, ry) = coefficient(1);
        let t = if cx.abs() >= cy.abs() {
            rx / cx
        } else {
            ry / cy
        };
        -(va[2] + t * (vb[2] - va[2]))
    }

    /// Hidden, visible, or `None` for a tie.
    fn judge(&self, mesh: usize, a: usize, b: usize, q: [f64; 2]) -> Option<bool> {
        let m = self.margin;
        let d = [q[0], q[1], -1.0];
        let lambda = self.lambda(mesh, a, b, q);
        let near_distance = -self.near;
        let mut hidden = false;
        for (mi, triangles) in self.tris.iter().enumerate() {
            for tri in triangles {
                if mi == mesh && tri.contains(&a) && tri.contains(&b) {
                    continue;
                }
                let v = tri.map(|i| self.verts[mi][i]);
                // A line parallel to the plane meets it nowhere.
                let Some((mu, inside)) = line_hit(d, v) else {
                    continue;
                };
                // Only a hit between the near plane and the point covers.
                if mu < near_distance - m || mu > lambda + m || inside < -m {
                    continue;
                }
                if inside <= m {
                    return None;
                }
                if (mu - lambda).abs() <= m || (mu - near_distance).abs() <= m {
                    return None;
                }
                // Strictly nearer the eye than the point, and not removed
                // by the near plane.
                if mu < lambda && mu > near_distance {
                    hidden = true;
                }
            }
        }
        for (mi, cut) in self.cut.iter().enumerate() {
            if !*cut {
                continue;
            }
            let mut crossings = 0;
            for tri in &self.tris[mi] {
                let v = tri.map(|i| self.verts[mi][i]);
                let Some((mu, inside)) = line_hit(d, v) else {
                    continue;
                };
                // Only a crossing on the eye's side of the plane counts.
                if mu > near_distance + m || inside < -m {
                    continue;
                }
                if inside <= m || (mu - near_distance).abs() <= m {
                    return None;
                }
                // On the eye's side of the plane: nearer than the plane,
                // behind the eye included.
                if mu < near_distance {
                    crossings += 1;
                }
            }
            if crossings % 2 == 1 {
                if lambda - near_distance <= m {
                    return None;
                }
                hidden = true;
            }
        }
        Some(hidden)
    }
}

#[derive(Default)]
struct Tally {
    segments: usize,
    robust: usize,
    hidden: usize,
    ambiguous: usize,
}

/// Eyes: outside the scene looking in, and inside it.
fn eyes() -> Vec<(&'static str, Frame3)> {
    let centre = Point3::new(1.5, 1.5, 1.5);
    vec![
        (
            "above",
            looking_along(Vec3::new(0.05, 0.1, 1.0), Point3::new(1.4, 1.6, 8.0)),
        ),
        (
            "side",
            looking_along(Vec3::new(-1.0, 0.2, 0.15), Point3::new(-5.0, 1.1, 1.4)),
        ),
        (
            "corner",
            looking_along(Vec3::new(1.0, 1.0, 0.9), Point3::new(8.0, 8.0, 7.0)),
        ),
        (
            "inside-x",
            looking_along(Vec3::new(-1.0, 0.1, 0.05), Point3::new(1.3, 1.1, 1.2)),
        ),
        (
            "inside-y",
            looking_along(Vec3::new(0.2, -1.0, 0.3), Point3::new(1.3, 1.1, 1.2)),
        ),
        ("inside-z", looking_along(Vec3::new(0.1, 0.2, 1.0), centre)),
        (
            "plan",
            frame(Vec3::X, Vec3::Y, Vec3::Z, Point3::new(1.5, 1.5, 6.0)),
        ),
    ]
}

/// What a run limits the view with.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    /// A near plane at the eye's own scale, nothing cut.
    Plain,
    /// A near plane through the middle of a box, which caps it.
    Cut,
    /// A near plane, a far plane through the middle of the scene's depths,
    /// and a crop rectangle in image coordinates.
    Limited,
}

/// The crop of a limited run: an off-centre image rectangle. Its bounds are
/// no rational the scene's 1/8 grid reaches (an edge on a bound would lie
/// outside it by the rounding of the bound, not be a tie).
fn crop() -> Aabb2 {
    Aabb2 {
        min: Point2::new(-0.153_7, -0.251_1),
        max: Point2::new(0.201_3, 0.101_9),
    }
}

/// The part of the image segment `p`-`r` inside `region`, as its length
/// (Liang-Barsky).
fn clipped_length(p: [f64; 2], r: [f64; 2], region: Aabb2) -> f64 {
    let d = [r[0] - p[0], r[1] - p[1]];
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for (origin, delta, min, max) in [
        (p[0], d[0], region.min.x, region.max.x),
        (p[1], d[1], region.min.y, region.max.y),
    ] {
        if delta == 0.0 {
            if origin < min || origin > max {
                return 0.0;
            }
            continue;
        }
        let (t0, t1) = ((min - origin) / delta, (max - origin) / delta);
        lo = lo.max(t0.min(t1));
        hi = hi.min(t0.max(t1));
    }
    if hi <= lo {
        return 0.0;
    }
    (hi - lo) * (d[0] * d[0] + d[1] * d[1]).sqrt()
}

fn check(seed: u64, name: &str, view_frame: Frame3, mode: Mode, tally: &mut Tally) {
    let with_near = mode == Mode::Cut;
    let meshes = scene(seed);
    let near = if with_near {
        // The depth halfway through a box's extent, so the plane cuts it.
        let mesh = &meshes[(seed as usize) % meshes.len()];
        let depths: Vec<f64> = mesh
            .positions
            .iter()
            .map(|p| (*p - view_frame.origin).dot(view_frame.z))
            .collect();
        let lo = depths.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = depths.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        ((lo + hi) / 2.0).min(-0.05)
    } else {
        -0.05
    };
    let oracle = Oracle::new(&meshes, &view_frame, near);
    // A far plane halfway through the depths of the scene in front of the
    // near plane: it removes the far half of what the edges reach.
    let far = (mode == Mode::Limited)
        .then(|| {
            let depths = oracle
                .verts
                .iter()
                .flatten()
                .map(|v| v[2])
                .filter(|w| *w < near);
            let lo = depths.clone().fold(f64::INFINITY, f64::min);
            let hi = depths.fold(f64::NEG_INFINITY, f64::max);
            (lo.is_finite() && lo < hi).then_some((lo + hi) / 2.0)
        })
        .flatten();
    let region = (mode == Mode::Limited).then(crop);
    let input = solids(&meshes);
    let mut view = View::perspective(view_frame, near).expect("valid");
    if let Some(far) = far {
        view = view.with_far(far);
    }
    let context = format!("scene {seed}, eye {name}, near {near}, far {far:?}");
    let mut request =
        HiddenLineRequest::new(&input, view, EdgeSelection::new(Some(0.5)), common::LIMITS);
    if let Some(region) = region {
        request = request.with_region(region);
    }
    let options = ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"));
    let outcome = scalar()
        .hidden_lines(&request, &options)
        .unwrap_or_else(|error| panic!("{context}: {error}"));

    let mut covered: BTreeMap<(u64, [u32; 2]), f64> = BTreeMap::new();
    for segment in &outcome.segments {
        let key = (segment.source.mesh.0, segment.source.vertices);
        *covered.entry(key).or_default() += (segment.end - segment.start).length();
    }
    for index in 0..meshes.len() {
        for edge in BOX_EDGES {
            let [a, b] = edge.map(|i| oracle.verts[index][i as usize]);
            // The part of the edge at or in front of the near plane and not
            // beyond the far plane.
            let (mut ta, mut tb) = match (a[2] <= near, b[2] <= near) {
                (true, true) => (0.0, 1.0),
                (false, false) => (0.0, 0.0),
                (true, false) => (0.0, (near - a[2]) / (b[2] - a[2])),
                (false, true) => ((near - a[2]) / (b[2] - a[2]), 1.0),
            };
            if let Some(far) = far {
                match (a[2] >= far, b[2] >= far) {
                    (true, true) => {}
                    (false, false) => tb = ta,
                    (true, false) => tb = tb.min((far - a[2]) / (b[2] - a[2])),
                    (false, true) => ta = ta.max((far - a[2]) / (b[2] - a[2])),
                }
            }
            let at = |t: f64| {
                let p = [
                    a[0] + t * (b[0] - a[0]),
                    a[1] + t * (b[1] - a[1]),
                    a[2] + t * (b[2] - a[2]),
                ];
                [p[0] / -p[2], p[1] / -p[2]]
            };
            let length = if tb > ta {
                let (p, r) = (at(ta), at(tb));
                match region {
                    Some(region) => clipped_length(p, r, region),
                    None => ((r[0] - p[0]).powi(2) + (r[1] - p[1]).powi(2)).sqrt(),
                }
            } else {
                0.0
            };
            let got = covered.remove(&(index as u64 + 1, edge));
            if length > 1e-7 {
                let got = got.unwrap_or_else(|| panic!("{context}: edge {edge:?} not drawn"));
                assert!(
                    (got - length).abs() <= 1e-9 * length.max(1.0),
                    "{context}: edge {edge:?} covers {got} of {length}"
                );
            } else if length == 0.0 {
                // A vertex exactly at the plane may leave a rounding-sized
                // piece.
                assert!(
                    got.is_none_or(|length| length < 1e-12),
                    "{context}: removed edge {edge:?} drawn: {a:?} {b:?} {got:?}"
                );
            }
        }
    }
    assert!(
        covered.values().all(|length| *length < 1e-7),
        "{context}: unexpected edges {covered:?}"
    );

    for segment in &outcome.segments {
        let (mesh, [a, b]) = (segment.source.mesh.0 as usize - 1, segment.source.vertices);
        if let Some(region) = region {
            for p in [segment.start, segment.end] {
                assert!(
                    (region.min.x..=region.max.x).contains(&p.x)
                        && (region.min.y..=region.max.y).contains(&p.y),
                    "{context}: {segment:?} leaves the crop"
                );
            }
        }
        for fraction in [0.25, 0.5, 0.75] {
            tally.segments += 1;
            let q = segment.start + (segment.end - segment.start) * fraction;
            let expected = match segment.visibility {
                Visibility::Hidden => true,
                Visibility::Visible => false,
                _ => unreachable!("unknown visibility"),
            };
            match oracle.judge(mesh, a as usize, b as usize, [q.x, q.y]) {
                None => tally.ambiguous += 1,
                Some(hidden) => {
                    tally.robust += 1;
                    tally.hidden += usize::from(hidden);
                    assert_eq!(hidden, expected, "{context}: {segment:?} at {fraction}");
                }
            }
        }
    }
}

fn run(mode: Mode) {
    const SEEDS: u64 = 40;
    let mut total = Tally::default();
    for (name, view_frame) in eyes() {
        let mut tally = Tally::default();
        for seed in 1..=SEEDS {
            check(seed, name, view_frame, mode, &mut tally);
        }
        eprintln!(
            "{name:>9}: {}: {} samples, {} robust, {} ambiguous",
            match mode {
                Mode::Plain => "plain",
                Mode::Cut => "near cuts",
                Mode::Limited => "far and crop",
            },
            tally.segments,
            tally.robust,
            tally.ambiguous
        );
        total.segments += tally.segments;
        total.robust += tally.robust;
        total.hidden += tally.hidden;
        // The fixed crop sees little from an eye inside the scene, so a
        // limited run is held to its totals.
        if mode == Mode::Limited {
            continue;
        }
        assert!(
            tally.robust >= 50,
            "{name}: only {} robust samples",
            tally.robust
        );
        // Both verdicts are exercised.
        assert!(
            tally.hidden >= 20 && tally.robust - tally.hidden >= 20,
            "{name}: {} hidden of {} robust samples",
            tally.hidden,
            tally.robust
        );
        // Inside the scene most of what is drawn is a tie or a cap.
        let share = if name.starts_with("inside") { 5 } else { 2 };
        assert!(
            tally.robust * share >= tally.segments,
            "{name}: too many ties: {} of {}",
            tally.ambiguous,
            tally.segments
        );
    }
    if mode == Mode::Limited {
        assert!(
            total.robust >= 5000 && total.robust * 2 >= total.segments,
            "{} robust of {} samples",
            total.robust,
            total.segments
        );
        assert!(
            total.hidden >= 500 && total.robust - total.hidden >= 500,
            "{} hidden of {} robust samples",
            total.hidden,
            total.robust
        );
    }
}

#[test]
fn the_provider_agrees_with_an_independent_eye_ray_cast() {
    run(Mode::Plain);
}

#[test]
fn caps_of_cut_solids_agree_with_an_independent_eye_line_parity() {
    run(Mode::Cut);
}

#[test]
fn a_far_depth_and_a_crop_agree_with_the_eye_ray_cast() {
    run(Mode::Limited);
}
