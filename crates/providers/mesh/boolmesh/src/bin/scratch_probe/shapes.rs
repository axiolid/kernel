//! Closed, outward meshes the scratch measurements run on.
//!
//! Shared by the `scratch_probe` binary and the tests that measure or pin
//! what the provider charges (`tests/scratch_bound.rs`,
//! `tests/batch_budget.rs`), so a fixture is defined once.

use axiolid_core::Point3;
use axiolid_mesh::{AttributeChannel, Blend, TriMesh};

/// An axis-aligned box from `min` to `max`.
pub fn box_at(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let positions: Vec<Point3> = vec![
        [x0, y0, z0].into(),
        [x1, y0, z0].into(),
        [x1, y1, z0].into(),
        [x0, y1, z0].into(),
        [x0, y0, z1].into(),
        [x1, y0, z1].into(),
        [x1, y1, z1].into(),
        [x0, y1, z1].into(),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(positions, indices)
}

/// Subdivide each triangle four ways, `levels` times, to grow triangle count.
pub fn subdivide(mesh: &TriMesh, levels: usize) -> TriMesh {
    let mut current = mesh.clone();
    for _ in 0..levels {
        let mut positions = current.positions.clone();
        let mut indices = Vec::new();
        for triangle in current.indices.as_chunks::<3>().0 {
            let [ia, ib, ic] = [triangle[0], triangle[1], triangle[2]];
            let (a, b, c) = (
                current.positions[ia as usize],
                current.positions[ib as usize],
                current.positions[ic as usize],
            );
            let base = positions.len() as u32;
            positions.push((a + b) * 0.5);
            positions.push((b + c) * 0.5);
            positions.push((c + a) * 0.5);
            let (ab, bc, ca) = (base, base + 1, base + 2);
            indices.extend_from_slice(&[ia, ab, ca]);
            indices.extend_from_slice(&[ab, ib, bc]);
            indices.extend_from_slice(&[ca, bc, ic]);
            indices.extend_from_slice(&[ab, bc, ca]);
        }
        current = TriMesh::new(positions, indices);
    }
    current
}

/// A closed prism over a counter-clockwise `profile`, from `z0` to `z1`,
/// with the cap triangulated by `cap` (counter-clockwise index triples into
/// the profile).
fn prism(profile: &[[f64; 2]], cap: &[[u32; 3]], z0: f64, z1: f64) -> TriMesh {
    let n = profile.len() as u32;
    let mut positions: Vec<Point3> = profile
        .iter()
        .map(|&[x, y]| Point3::new(x, y, z0))
        .collect();
    positions.extend(profile.iter().map(|&[x, y]| Point3::new(x, y, z1)));
    let mut indices = Vec::new();
    for &[a, b, c] in cap {
        indices.extend([a, c, b]);
        indices.extend([a + n, b + n, c + n]);
    }
    for i in 0..n {
        let j = (i + 1) % n;
        indices.extend([i, j, j + n, i, j + n, i + n]);
    }
    TriMesh::new(positions, indices)
}

/// A comb of `k` teeth along +y on a spine along x, extruded from `z0` to
/// `z1`: `16k + 12` triangles. The teeth cover x in `[0, 1]`, each half its
/// pitch wide and starting `phase` pitches in; the spine overhangs them.
fn comb(k: usize, phase: f64, length: f64, z0: f64, z1: f64) -> TriMesh {
    let pitch = 1.0 / k as f64;
    // Counter-clockwise: the spine's bottom edge, its right end, then the
    // top edge right to left with every tooth on it, then its left end.
    let mut profile = vec![[-0.2, -0.3], [1.2, -0.3], [1.2, 0.0]];
    for i in (0..k).rev() {
        let left = (i as f64 + phase) * pitch;
        let right = left + 0.5 * pitch;
        profile.extend([[right, 0.0], [right, length], [left, length], [left, 0.0]]);
    }
    profile.push([-0.2, 0.0]);
    // The spine as a fan from its first corner over every vertex on its
    // boundary, and each tooth as two triangles.
    let mut spine = vec![1u32, 2];
    let mut cap = Vec::new();
    for i in 0..k as u32 {
        let base = 3 + 4 * i;
        spine.extend([base, base + 3]);
        cap.push([base, base + 1, base + 2]);
        cap.push([base, base + 2, base + 3]);
    }
    spine.push(profile.len() as u32 - 1);
    cap.extend(spine.windows(2).map(|pair| [0, pair[0], pair[1]]));
    prism(&profile, &cap, z0, z1)
}

/// `mesh` mirrored across the plane `x = y`, kept outward.
fn mirrored(mut mesh: TriMesh) -> TriMesh {
    for p in &mut mesh.positions {
        std::mem::swap(&mut p.x, &mut p.y);
    }
    for corner in mesh.indices.as_chunks_mut::<3>().0 {
        corner.swap(1, 2);
    }
    mesh
}

/// Two perpendicular combs of `k` teeth whose teeth cross `k^2` times.
pub fn combs(k: usize) -> (TriMesh, TriMesh) {
    let subject = comb(k, 0.25, 1.2, 0.0, 1.0);
    let tool = mirrored(comb(k, 0.3, 1.23, 0.37, 1.41));
    (subject, tool)
}

/// A per-vertex channel of `width` scalars on `mesh`.
pub fn with_channel(mut mesh: TriMesh, width: usize) -> TriMesh {
    let values = mesh
        .positions
        .iter()
        .flat_map(|p| [p.x, p.y, p.z, p.x * p.y].into_iter().cycle().take(width))
        .collect();
    mesh.attributes
        .push(AttributeChannel::new("uv", values, width, Blend::Linear));
    mesh
}

/// A small deterministic generator, so the soups are the same every run.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// `n` parallel slabs, disjoint, rotated by `angle` about z and tilted by
/// `tilt` about x: one closed solid of `12n` triangles.
fn slabs(n: usize, gap: f64, angle: f64, tilt: f64) -> TriMesh {
    let mut positions = Vec::new();
    let mut indices = Vec::new();
    let (sa, ca) = angle.sin_cos();
    let (st, ct) = tilt.sin_cos();
    for i in 0..n {
        let x0 = -0.5 + i as f64 / n as f64;
        let x1 = x0 + gap / n as f64;
        let slab = box_at([x0, -0.7, -0.4], [x1, 0.7, 0.4]);
        let base = positions.len() as u32;
        positions.extend(slab.positions.iter().map(|p| {
            let (y, z) = (p.y * ct - p.z * st, p.y * st + p.z * ct);
            Point3::new(p.x * ca - y * sa, p.x * sa + y * ca, z)
        }));
        indices.extend(slab.indices.iter().map(|i| i + base));
    }
    TriMesh::new(positions, indices)
}

/// Two random soups of slabs crossing at random angles, from `seed`.
pub fn soup(seed: u64) -> (TriMesh, TriMesh) {
    let mut random = Lcg(seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1);
    let one = |random: &mut Lcg| {
        let n = 2 + (random.next() * 23.0) as usize;
        let gap = 0.2 + 0.6 * random.next();
        let angle = std::f64::consts::PI * random.next();
        let tilt = 0.6 * (random.next() - 0.5);
        slabs(n, gap, angle, tilt)
    };
    let subject = one(&mut random);
    let tool = one(&mut random);
    (subject, tool)
}

/// A slab pierced by an `m` by `m` grid of triangular pins: every pin
/// leaves a three-edged hole in both large faces, the most loops a face of
/// that many half-edges can have, so the most triangles per half-edge.
pub fn pins(m: usize) -> (TriMesh, TriMesh) {
    let subject = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 0.2]);
    let mut positions = Vec::new();
    let mut indices = Vec::new();
    let pitch = 1.0 / m as f64;
    let cap = [[0u32, 1, 2]];
    for i in 0..m {
        for j in 0..m {
            let (x, y) = ((i as f64 + 0.3) * pitch, (j as f64 + 0.35) * pitch);
            let r = 0.25 * pitch;
            let profile = [[x, y], [x + r, y + 0.1 * r], [x + 0.2 * r, y + r]];
            let pin = prism(&profile, &cap, -0.5, 0.7);
            let base = positions.len() as u32;
            positions.extend(pin.positions);
            indices.extend(pin.indices.iter().map(|i| i + base));
        }
    }
    (subject, TriMesh::new(positions, indices))
}
