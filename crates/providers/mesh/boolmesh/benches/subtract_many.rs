//! `subtract_many` throughput against the ADR 0014 sequential baseline.
//!
//! Run with: `cargo bench -p axiolid-mesh-boolean-boolmesh`
//!
//! The question this answers is whether the batch override earns its
//! complexity. ADR 0014 recorded the sequential loop at n=16: 6.95 ms and
//! n=64: 48.68 ms. The `sequential` column reproduces that baseline in the
//! same process, so the speedup is a ratio of two measurements taken on the
//! same machine in the same run rather than against a remembered figure.

use std::hint::black_box;
use std::time::Instant;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;

/// Builds a subject plus its cutters for one benchmark case.
type Layout = fn(usize) -> (TriMesh, Vec<TriMesh>);

/// Axis-aligned box as a closed outward-wound solid.
fn box_solid(cx: f64, cy: f64, cz: f64, sx: f64, sy: f64, sz: f64) -> TriMesh {
    let (hx, hy, hz) = (sx / 2.0, sy / 2.0, sz / 2.0);
    let positions = vec![
        Point3::new(cx - hx, cy - hy, cz - hz),
        Point3::new(cx + hx, cy - hy, cz - hz),
        Point3::new(cx + hx, cy + hy, cz - hz),
        Point3::new(cx - hx, cy + hy, cz - hz),
        Point3::new(cx - hx, cy - hy, cz + hz),
        Point3::new(cx + hx, cy - hy, cz + hz),
        Point3::new(cx + hx, cy + hy, cz + hz),
        Point3::new(cx - hx, cy + hy, cz + hz),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(positions, indices)
}

/// Signed volume, so the two paths can be proven to agree.
fn volume(mesh: &TriMesh) -> f64 {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let (a, b, c) = (
                mesh.positions[t[0] as usize],
                mesh.positions[t[1] as usize],
                mesh.positions[t[2] as usize],
            );
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0
}

/// The IFC-dominant layout: a wall with N non-overlapping openings.
///
/// Windows and doors sit at distinct positions, so the overlap graph is empty
/// and every cutter fuses into one group. This is the case the override is
/// designed for.
fn wall_and_openings(n: usize) -> (TriMesh, Vec<TriMesh>) {
    let length = n as f64 + 1.0;
    let wall = box_solid(length / 2.0, 0.1, 1.5, length, 0.2, 3.0);
    let tools = (0..n)
        .map(|i| box_solid(0.75 + i as f64, 0.1, 1.0, 0.5, 0.5, 1.0))
        .collect();
    (wall, tools)
}

/// The worst case: every cutter overlaps every other.
///
/// Concentric cutters of increasing size, so the overlap graph is COMPLETE and
/// greedy colouring yields N groups of one. The override then performs exactly
/// the sequential work plus bounds construction. This measures the overhead
/// the batch path adds when its precondition fails outright, which is what
/// decides whether it is safe to enable unconditionally.
fn wall_and_nested_openings(n: usize) -> (TriMesh, Vec<TriMesh>) {
    let wall = box_solid(2.0, 0.1, 1.5, 4.0, 0.2, 3.0);
    let tools = (0..n)
        .map(|i| {
            let scale = 0.2 + 0.6 * (i as f64) / (n.max(2) - 1) as f64;
            box_solid(2.0, 0.1, 1.5, scale, 0.5, scale)
        })
        .collect();
    (wall, tools)
}

/// Time a closure once, after a warm-up run.
fn timed<F: FnMut() -> TriMesh>(mut f: F) -> (f64, f64) {
    let warm = f();
    let start = Instant::now();
    let result = black_box(f());
    let millis = start.elapsed().as_secs_f64() * 1e3;
    let (a, b) = (volume(&warm), volume(&result));
    assert!(
        (a - b).abs() <= 1e-9 * a.abs().max(1.0),
        "repeated runs must agree"
    );
    (millis, volume(&result))
}

/// The pre-override default: a serial loop of one-at-a-time booleans.
fn sequential(
    provider: &BoolmeshBoolean,
    subject: &TriMesh,
    tools: &[TriMesh],
    options: &ExecutionOptions,
) -> TriMesh {
    let mut current = subject.clone();
    for tool in tools {
        current = provider
            .boolean(&current, tool, BooleanOperator::Difference, options)
            .expect("sequential difference")
            .mesh;
    }
    current
}

fn main() {
    let provider = BoolmeshBoolean::new();
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let cases: [(&str, Layout); 2] = [
        ("disjoint openings (IFC-dominant layout)", wall_and_openings),
        (
            "nested openings (complete overlap graph)",
            wall_and_nested_openings,
        ),
    ];

    for (label, build) in cases {
        println!("\n{label}");
        println!(
            "{:>5}  {:>13}  {:>13}  {:>9}  {:>10}",
            "n", "sequential", "grouped", "speedup", "volumes"
        );
        for &n in &[1usize, 4, 16, 64] {
            let (subject, tools) = build(n);
            let (seq_ms, seq_vol) = timed(|| sequential(&provider, &subject, &tools, &options));
            let (grp_ms, grp_vol) = timed(|| {
                provider
                    .subtract_many(&subject, &tools, &options)
                    .expect("grouped difference")
                    .mesh
            });
            // The claim is an optimisation, not a behaviour change: if the
            // volumes disagree the speedup is meaningless.
            let gap = (seq_vol - grp_vol).abs();
            let agree = if gap <= 1e-9 * seq_vol.abs().max(1.0) {
                "equal".to_owned()
            } else {
                format!("{gap:.2e}")
            };
            println!(
                "{n:>5}  {seq_ms:>10.2} ms  {grp_ms:>10.2} ms  {:>8.2}x  {agree:>10}",
                seq_ms / grp_ms
            );
        }
    }
}
