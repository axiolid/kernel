//! Cost of carrying channels through a pairwise boolean (#116).
//! Best-of-5 wall clock per operand density, with and without a channel.

use std::time::Instant;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Tolerance};
use axiolid_mesh::{AttributeChannel, Blend, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;

/// A box whose six faces are each an `n` x `n` grid of quads.
fn dense_box(c: f64, s: f64, n: usize) -> TriMesh {
    let (mut ps, mut ix) = (Vec::new(), Vec::new());
    // (axis, sign): the face lies at coordinate `sign` along `axis`.
    for axis in 0..3 {
        for sign in [-1.0f64, 1.0] {
            let base = ps.len() as u32;
            let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
            for i in 0..=n {
                for j in 0..=n {
                    let mut q = [0.0; 3];
                    q[axis] = sign;
                    q[u] = -1.0 + 2.0 * i as f64 / n as f64;
                    q[v] = -1.0 + 2.0 * j as f64 / n as f64;
                    let h = s / 2.0;
                    ps.push(Point3::new(c + h * q[0], c + h * q[1], c + h * q[2]));
                }
            }
            let at = |i: usize, j: usize| base + (i * (n + 1) + j) as u32;
            for i in 0..n {
                for j in 0..n {
                    let (a, b) = (at(i, j), at(i + 1, j));
                    let (cc, d) = (at(i + 1, j + 1), at(i, j + 1));
                    // (u, v, axis) is right-handed, so a->b->cc is CCW seen
                    // from +axis: outward on the + face, flipped on the -.
                    if sign > 0.0 {
                        ix.extend_from_slice(&[a, b, cc, a, cc, d]);
                    } else {
                        ix.extend_from_slice(&[a, cc, b, a, d, cc]);
                    }
                }
            }
        }
    }
    TriMesh::new(ps, ix)
}

fn with_uv(mut m: TriMesh) -> TriMesh {
    let values = m.positions.iter().flat_map(|p| [p.x, p.y + p.z]).collect();
    m.attributes
        .push(AttributeChannel::new("uv", values, 2, Blend::Linear));
    m
}

fn best_ms(a: &TriMesh, b: &TriMesh) -> (f64, usize) {
    let opts = ExecutionOptions::new(Tolerance::METRE);
    let (mut best, mut tris) = (f64::MAX, 0);
    for _ in 0..5 {
        let t = Instant::now();
        let o = BoolmeshBoolean::new()
            .boolean(a, b, BooleanOperator::Difference, &opts)
            .expect("difference");
        best = best.min(t.elapsed().as_secs_f64() * 1e3);
        tris = o.mesh.triangle_count();
    }
    (best, tris)
}

fn main() {
    println!("n,input_tris,result_tris,plain_ms,uv_ms");
    for n in [8, 16, 32, 64] {
        let (a, b) = (dense_box(0.0, 2.0, n), dense_box(0.7, 2.0, n));
        let (plain, tris) = best_ms(&a, &b);
        let (uv, _) = best_ms(&with_uv(a.clone()), &with_uv(b.clone()));
        let input = a.triangle_count() + b.triangle_count();
        println!("{n},{input},{tris},{plain:.2},{uv:.2}");
    }
    println!("analytic_cutters,plain_ms,uv_ms");
    for k in [4, 16, 64] {
        println!(
            "{k},{:.3},{:.3}",
            analytic_ms(k, false),
            analytic_ms(k, true)
        );
    }
}

/// Analytic box path (wall minus a row of openings), best of 5, ms.
fn analytic_ms(n_cut: usize, uv: bool) -> f64 {
    let wall = dense_box(0.0, 2.0, 1);
    let opts = ExecutionOptions::new(Tolerance::METRE);
    let mut tools: Vec<TriMesh> = (0..n_cut)
        .map(|i| {
            let x = -0.9 + 1.8 * (i as f64 + 0.5) / n_cut as f64;
            dense_box_at([x, 0.0, 0.0], [0.9 / n_cut as f64, 3.0, 0.5])
        })
        .collect();
    let wall = if uv { with_uv(wall) } else { wall };
    if uv {
        tools = tools.into_iter().map(with_uv).collect();
    }
    let mut best = f64::MAX;
    for _ in 0..5 {
        let t = Instant::now();
        BoolmeshBoolean::new()
            .subtract_boxes_analytic(&wall, &tools, &opts, 10_000_000)
            .expect("ok")
            .expect("boxes");
        best = best.min(t.elapsed().as_secs_f64() * 1e3);
    }
    best
}

/// An axis-aligned box (12 triangles) with centre `c` and size `s`.
fn dense_box_at(c: [f64; 3], s: [f64; 3]) -> TriMesh {
    let mut m = dense_box(0.0, 2.0, 1);
    for p in &mut m.positions {
        *p = Point3::new(
            c[0] + p.x * s[0] / 2.0,
            c[1] + p.y * s[1] / 2.0,
            c[2] + p.z * s[2] / 2.0,
        );
    }
    m
}
