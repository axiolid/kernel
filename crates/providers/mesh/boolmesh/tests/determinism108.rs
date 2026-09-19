//! TEMP probe for kernel#108: does one process produce a stable result?
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
fn boxm(c: [f64; 3], h: [f64; 3]) -> TriMesh {
    let v = (0..8)
        .map(|i| {
            let sx = if i & 1 == 0 { -h[0] } else { h[0] };
            let sy = if i & 2 == 0 { -h[1] } else { h[1] };
            let sz = if i & 4 == 0 { -h[2] } else { h[2] };
            Point3::new(c[0] + sx, c[1] + sy, c[2] + sz)
        })
        .collect();
    TriMesh::new(
        v,
        vec![
            0, 2, 1, 1, 2, 3, 4, 5, 6, 5, 7, 6, 0, 1, 4, 1, 5, 4, 2, 6, 3, 3, 6, 7, 0, 4, 2, 2, 4,
            6, 1, 3, 5, 3, 7, 5,
        ],
    )
}
#[test]
fn same_input_same_bytes_within_one_process() {
    let p = BoolmeshBoolean::new();
    let o = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let host = boxm([8.5, 0.1, 1.5], [8.5, 0.1, 1.5]);
    let tools: Vec<TriMesh> = (0..16)
        .map(|i| boxm([0.75 + i as f64, 0.1, 1.0], [0.25, 0.25, 0.5]))
        .collect();
    let fp = |m: &TriMesh| {
        (
            m.positions.len(),
            m.indices.len(),
            format!("{:?}", m.positions),
        )
    };
    let first = fp(&p.subtract_many(&host, &tools, &o).expect("sub").mesh);
    for _ in 0..30 {
        assert_eq!(
            fp(&p.subtract_many(&host, &tools, &o).expect("sub").mesh),
            first
        );
    }
}

/// The probe shape: ONE boolean against a fused multi-component tool.
///
/// Regression test for kernel#108. This permuted connectivity run to run
/// while `subtract_many` stayed stable: positions were byte-identical and
/// only `Half` ordering drifted. Two independent causes, both fixed —
/// `HashMap` iteration in `boolean45`, and `Rc::as_ptr` tie-breaks in the
/// ear-clip comparators. Removing either one alone makes this fail again.
#[test]
fn fused_multicomponent_tool_is_stable() {
    let p = BoolmeshBoolean::new();
    let o = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let host = boxm([8.5, 0.1, 1.5], [8.5, 0.1, 1.5]);
    let tools: Vec<TriMesh> = (0..16)
        .map(|i| boxm([0.75 + i as f64, 0.1, 1.0], [0.25, 0.25, 0.5]))
        .collect();
    // Concatenate every cutter into ONE disconnected mesh, then subtract
    // it in a single boolean -- exactly what the benchmark probe does.
    let mut pos = Vec::new();
    let mut idx: Vec<u32> = Vec::new();
    for t in &tools {
        let base = pos.len() as u32;
        pos.extend(t.positions.iter().copied());
        idx.extend(t.indices.iter().map(|i| i + base));
    }
    let fused = TriMesh::new(pos, idx);
    let op = BooleanOperator::Difference;
    let fp = |m: &TriMesh| format!("{:?}{:?}", m.positions, m.indices);
    let first = fp(&p.boolean(&host, &fused, op, &o).expect("fused").mesh);
    for _ in 0..40 {
        assert_eq!(
            fp(&p.boolean(&host, &fused, op, &o).expect("fused").mesh),
            first
        );
    }
}
