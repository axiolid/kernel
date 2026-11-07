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

use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{View, Visibility};
use common::{draw, frame, looking_along, scalar, scene, solids, BOX_EDGES};

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

fn check(seed: u64, name: &str, view_frame: Frame3, with_near: bool, tally: &mut Tally) {
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
    let input = solids(&meshes);
    let view = View::perspective(view_frame, near).expect("valid");
    let context = format!("scene {seed}, eye {name}, near {near}");
    let outcome = draw(&scalar(), &input, view);

    let mut covered: BTreeMap<(u64, [u32; 2]), f64> = BTreeMap::new();
    for segment in &outcome.segments {
        let key = (segment.source.mesh.0, segment.source.vertices);
        *covered.entry(key).or_default() += (segment.end - segment.start).length();
    }
    for index in 0..meshes.len() {
        for edge in BOX_EDGES {
            let [a, b] = edge.map(|i| oracle.verts[index][i as usize]);
            // The part of the edge at or in front of the near plane.
            let (ta, tb) = match (a[2] <= near, b[2] <= near) {
                (true, true) => (0.0, 1.0),
                (false, false) => (0.0, 0.0),
                (true, false) => (0.0, (near - a[2]) / (b[2] - a[2])),
                (false, true) => ((near - a[2]) / (b[2] - a[2]), 1.0),
            };
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
                ((r[0] - p[0]).powi(2) + (r[1] - p[1]).powi(2)).sqrt()
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

fn run(with_near: bool) {
    const SEEDS: u64 = 40;
    for (name, view_frame) in eyes() {
        let mut tally = Tally::default();
        for seed in 1..=SEEDS {
            check(seed, name, view_frame, with_near, &mut tally);
        }
        eprintln!(
            "{name:>9}: near cuts {with_near}: {} samples, {} robust, {} ambiguous",
            tally.segments, tally.robust, tally.ambiguous
        );
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
}

#[test]
fn the_provider_agrees_with_an_independent_eye_ray_cast() {
    run(false);
}

#[test]
fn caps_of_cut_solids_agree_with_an_independent_eye_line_parity() {
    run(true);
}
