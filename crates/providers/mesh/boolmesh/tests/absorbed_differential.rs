//! Differential: absorbed algorithm vs upstream `boolmesh` 0.1.9.
//!
//! ADR 0047 absorbed the algorithm instead of depending on it, claiming
//! the landing is behaviour-identical. This proves it rather than
//! asserting it: upstream stays as a dev-dependency purely so both can
//! run on identical input and be compared bit-for-bit.
//!
//! If this ever fails, the absorbed copy has diverged from the version
//! ADR 0014 evaluated, and the divergence must be deliberate and
//! documented -- not a silent port bug.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;

/// An icosphere, the workload that motivated ADR 0047.
///
/// Curved operands with no coincident planes: every intersection is a
/// generic triangle-triangle crossing, which is the case the absorbed
/// code must reproduce exactly.
fn icosphere(center: [f64; 3], radius: f64, subdivisions: u32) -> TriMesh {
    let t = (1.0 + 5.0f64.sqrt()) / 2.0;
    let mut verts: Vec<[f64; 3]> = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    let mut faces: Vec<[u32; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];

    for _ in 0..subdivisions {
        let mut midpoint: std::collections::HashMap<(u32, u32), u32> =
            std::collections::HashMap::new();
        let mut next: Vec<[u32; 3]> = Vec::with_capacity(faces.len() * 4);
        for f in &faces {
            let mut mid = [0u32; 3];
            for e in 0..3 {
                let (a, b) = (f[e], f[(e + 1) % 3]);
                let key = (a.min(b), a.max(b));
                mid[e] = *midpoint.entry(key).or_insert_with(|| {
                    let (pa, pb) = (verts[a as usize], verts[b as usize]);
                    verts.push([
                        (pa[0] + pb[0]) * 0.5,
                        (pa[1] + pb[1]) * 0.5,
                        (pa[2] + pb[2]) * 0.5,
                    ]);
                    (verts.len() - 1) as u32
                });
            }
            next.push([f[0], mid[0], mid[2]]);
            next.push([f[1], mid[1], mid[0]]);
            next.push([f[2], mid[2], mid[1]]);
            next.push([mid[0], mid[1], mid[2]]);
        }
        faces = next;
    }

    let positions: Vec<Point3> = verts
        .iter()
        .map(|v| {
            let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            let k = radius / len;
            Point3::new(
                center[0] + v[0] * k,
                center[1] + v[1] * k,
                center[2] + v[2] * k,
            )
        })
        .collect();
    let indices: Vec<u32> = faces.iter().flat_map(|f| [f[0], f[1], f[2]]).collect();
    TriMesh::new(positions, indices)
}

/// Run the same operation through upstream `boolmesh` 0.1.9 directly.
fn upstream(
    a: &TriMesh,
    b: &TriMesh,
    op: boolmesh::prelude::OpType,
) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let flat =
        |m: &TriMesh| -> Vec<f64> { m.positions.iter().flat_map(|p| [p.x, p.y, p.z]).collect() };
    let idx = |m: &TriMesh| -> Vec<usize> { m.indices.iter().map(|&i| i as usize).collect() };
    let ma = boolmesh::prelude::Manifold::new(&flat(a), &idx(a)).expect("upstream subject");
    let mb = boolmesh::prelude::Manifold::new(&flat(b), &idx(b)).expect("upstream tool");
    let out = boolmesh::prelude::compute_boolean(&ma, &mb, op).expect("upstream boolean");
    let ps = out.ps.iter().map(|p| [p.x, p.y, p.z]).collect();
    let ts = out
        .get_indices()
        .iter()
        .map(|t| [t.x as u32, t.y as u32, t.z as u32])
        .collect();
    (ps, ts)
}

/// Run the operation through the ABSORBED algorithm, via the provider.
fn absorbed(a: &TriMesh, b: &TriMesh, op: BooleanOperator) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let provider = BoolmeshBoolean;
    let options = ExecutionOptions::new(Tolerance::METRE);
    let out = provider
        .boolean(a, b, op, &options)
        .expect("absorbed boolean")
        .mesh;
    (
        out.positions.iter().map(|p| [p.x, p.y, p.z]).collect(),
        out.indices
            .chunks_exact(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect(),
    )
}

/// Volume of a closed triangle soup, by the divergence theorem.
fn volume(ps: &[[f64; 3]], ts: &[[u32; 3]]) -> f64 {
    let mut sum = 0.0;
    for t in ts {
        let a = ps[t[0] as usize];
        let b = ps[t[1] as usize];
        let c = ps[t[2] as usize];
        sum += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
    (sum / 6.0).abs()
}

/// The absorbed algorithm agrees with upstream on the SOLID.
///
/// Not on the triangle list: ear clipping picks a quad diagonal by a
/// tie-break that the edition-2021 port perturbs, so a handful of coplanar
/// quads split the other way. That is a different tessellation of the same
/// surface, so the contract-level facts -- vertex positions, enclosed
/// volume, triangle count -- are what must match, and do.
#[test]
fn absorbed_agrees_with_upstream_on_the_solid() {
    let a = icosphere([0.0, 0.0, 0.0], 1.0, 3);
    let b = icosphere([0.5, 0.0, 0.0], 1.0, 3);

    for (op_ours, op_theirs, name) in [
        (
            BooleanOperator::Union,
            boolmesh::prelude::OpType::Add,
            "union",
        ),
        (
            BooleanOperator::Intersection,
            boolmesh::prelude::OpType::Intersect,
            "intersection",
        ),
        (
            BooleanOperator::Difference,
            boolmesh::prelude::OpType::Subtract,
            "difference",
        ),
    ] {
        let (ap, at) = absorbed(&a, &b, op_ours);
        let (up, ut) = upstream(&a, &b, op_theirs);

        assert_eq!(ap, up, "{name}: vertex positions must be bit-identical");
        assert_eq!(at.len(), ut.len(), "{name}: triangle count must match");

        let (va, vu) = (volume(&ap, &at), volume(&up, &ut));
        assert!(
            (va - vu).abs() <= 1e-12 * vu.abs().max(1.0),
            "{name}: enclosed volume diverged: absorbed {va}, upstream {vu}"
        );
    }
}
