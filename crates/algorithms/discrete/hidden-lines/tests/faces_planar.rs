//! `FaceGroups::Planar` (ADR 0106) gives the regions of the `PerTriangle`
//! request that carries the groups its definition names.
//!
//! The reference here is independent of the provider's resolver: normals in
//! dyadic arithmetic, every triangle pair that owns an edge by itself found
//! by comparing vertex pairs, and the components found by relaxing labels to
//! the smallest index until nothing changes.

mod common;

use axiolid_core::{Point2, Point3};
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;
use axiolid_hidden_lines::ScalarVisibleFaces;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    FaceGroup, FaceSource, HiddenLineMesh, MeshId, MeshRole, View, VisibleFaceOutcome,
};
use common::{box_mesh, plan_frame, scene, try_faces, views};

fn d(value: f64) -> Dyadic {
    Dyadic::from_f64(value)
}

/// `(b - a) × (c - a)`.
fn normal(mesh: &TriMesh, t: usize) -> [Dyadic; 3] {
    let p = |k: usize| mesh.positions[mesh.indices[3 * t + k] as usize];
    let (a, b, c) = (p(0), p(1), p(2));
    let diff = |p: Point3, q: Point3| {
        [
            d(p.x).sub(&d(q.x)),
            d(p.y).sub(&d(q.y)),
            d(p.z).sub(&d(q.z)),
        ]
    };
    let (u, v) = (diff(b, a), diff(c, a));
    [
        u[1].mul(&v[2]).sub(&u[2].mul(&v[1])),
        u[2].mul(&v[0]).sub(&u[0].mul(&v[2])),
        u[0].mul(&v[1]).sub(&u[1].mul(&v[0])),
    ]
}

fn is_zero(x: &Dyadic) -> bool {
    x.sign() == Some(Sign::Zero)
}

/// Whether the normals are parallel and point the same way.
fn flat(a: &[Dyadic; 3], b: &[Dyadic; 3]) -> bool {
    let cross = [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ];
    let dot = a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]));
    cross.iter().all(is_zero) && dot.sign() == Some(Sign::Positive)
}

/// The groups `FaceGroups::Planar` names for `mesh`.
fn reference(mesh: &TriMesh) -> Vec<FaceGroup> {
    let count = mesh.triangle_count();
    let edges = |t: usize| {
        let c = |k: usize| mesh.indices[3 * t + k];
        [(c(0), c(1)), (c(1), c(2)), (c(2), c(0))].map(|(i, j)| (i.min(j), i.max(j)))
    };
    let normals: Vec<_> = (0..count).map(|t| normal(mesh, t)).collect();
    let mut label: Vec<u32> = (0..count as u32).collect();
    let mut joins = Vec::new();
    for s in 0..count {
        for t in s + 1..count {
            for edge in edges(s) {
                let owners = (0..count).filter(|&u| edges(u).contains(&edge)).count();
                if owners == 2 && edges(t).contains(&edge) && flat(&normals[s], &normals[t]) {
                    joins.push((s, t));
                }
            }
        }
    }
    loop {
        let mut changed = false;
        for &(s, t) in &joins {
            let low = label[s].min(label[t]);
            for u in [s, t] {
                if label[u] != low {
                    label[u] = low;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    label.into_iter().map(FaceGroup).collect()
}

/// Equal up to `exact_escalations`, which depends on the path of the decisions.
fn assert_same(a: &VisibleFaceOutcome, b: &VisibleFaceOutcome, context: &str) {
    let mut evidence = b.evidence;
    evidence.exact_escalations = a.evidence.exact_escalations;
    assert_eq!(a.regions, b.regions, "{context}");
    assert_eq!(a.evidence, evidence, "{context}");
}

/// The regions and the evidence that does not count work. Reversing the
/// triangles of a large scene changes `candidate_pairs` by a few pairs with
/// any group mode (the broad phase is not winding-blind); the regions and
/// every count of the drawing stay equal.
fn assert_same_drawing(a: &VisibleFaceOutcome, b: &VisibleFaceOutcome, context: &str) {
    let mut evidence = b.evidence;
    evidence.exact_escalations = a.evidence.exact_escalations;
    evidence.candidate_pairs = a.evidence.candidate_pairs;
    assert_eq!(a.regions, b.regions, "{context}");
    assert_eq!(a.evidence, evidence, "{context}");
}

fn reversed(mesh: &TriMesh) -> TriMesh {
    let mut flipped = mesh.clone();
    for corner in flipped.indices.as_chunks_mut::<3>().0 {
        corner.swap(1, 2);
    }
    flipped
}

fn planar_meshes(boxes: &[TriMesh]) -> Vec<HiddenLineMesh<'_>> {
    boxes
        .iter()
        .enumerate()
        .map(|(i, m)| {
            HiddenLineMesh::new(MeshId(i as u64 + 1), MeshRole::Solid, m).with_planar_groups()
        })
        .collect()
}

#[test]
fn planar_equals_the_per_triangle_reference_on_seeded_scenes() {
    let provider = ScalarVisibleFaces::new();
    let (mut answered, mut merged, mut split) = (0, 0, 0);
    for seed in 0..20 {
        let boxes = scene(seed);
        let references: Vec<Vec<FaceGroup>> = boxes.iter().map(reference).collect();
        for groups in &references {
            let distinct = groups
                .iter()
                .collect::<std::collections::BTreeSet<_>>()
                .len();
            // A box is six faces; rotated ones fall apart where rounding
            // breaks exact coplanarity.
            if distinct == 6 {
                merged += 1;
            } else {
                split += 1;
            }
        }
        let planar: Vec<_> = boxes
            .iter()
            .enumerate()
            .map(|(i, m)| {
                HiddenLineMesh::new(MeshId(i as u64 + 1), MeshRole::Solid, m).with_planar_groups()
            })
            .collect();
        let explicit: Vec<_> = boxes
            .iter()
            .zip(&references)
            .enumerate()
            .map(|(i, (m, g))| {
                HiddenLineMesh::new(MeshId(i as u64 + 1), MeshRole::Solid, m).with_groups(g)
            })
            .collect();
        for (name, frame) in views()
            .into_iter()
            .filter(|(n, _)| ["plan", "-Y", "isometric", "oblique"].contains(n))
        {
            let view = View::new(frame);
            let context = format!("seed {seed}, {name}");
            let a = try_faces(&provider, &planar, view, None);
            let b = try_faces(&provider, &explicit, view, None);
            match (a, b) {
                (Ok(a), Ok(b)) => {
                    assert_same(&a, &b, &context);
                    answered += 1;
                }
                (a, b) => assert_eq!(a.err(), b.err(), "{context}"),
            }
        }
    }
    assert!(answered > 40, "{answered} scenes answered");
    assert!(
        merged > 0 && split > 0,
        "{merged} merged, {split} split boxes"
    );
}

#[test]
fn reversed_winding_gives_the_same_output() {
    let provider = ScalarVisibleFaces::new();
    for seed in 0..10 {
        let boxes = scene(seed);
        let flipped: Vec<_> = boxes.iter().map(reversed).collect();
        let (original, other) = (planar_meshes(&boxes), planar_meshes(&flipped));
        for (name, frame) in views()
            .into_iter()
            .filter(|(n, _)| ["plan", "oblique"].contains(n))
        {
            let view = View::new(frame);
            let context = format!("seed {seed}, {name}");
            match (
                try_faces(&provider, &original, view, None),
                try_faces(&provider, &other, view, None),
            ) {
                (Ok(a), Ok(b)) => assert_same_drawing(&a, &b, &context),
                (a, b) => assert_eq!(a.err(), b.err(), "{context}"),
            }
        }
    }
}

#[test]
fn the_plan_cube_shows_its_top_as_the_group_of_the_smaller_triangle() {
    let cube = box_mesh([0.0; 3], [2.0; 3]);
    let meshes = [HiddenLineMesh::new(MeshId(7), MeshRole::Solid, &cube).with_planar_groups()];
    let outcome = try_faces(
        &ScalarVisibleFaces::new(),
        &meshes,
        View::new(plan_frame()),
        None,
    )
    .expect("faces");
    assert_eq!(outcome.regions.len(), 1);
    let region = &outcome.regions[0];
    assert_eq!(region.source, FaceSource::Group(FaceGroup(2)));
    assert_eq!(region.triangles, [2, 3]);
}

#[test]
fn two_disjoint_coplanar_squares_are_two_groups() {
    let square = |x: f64, y: f64| {
        [
            Point3::new(x, y, 0.0),
            Point3::new(x + 1.0, y, 0.0),
            Point3::new(x + 1.0, y + 1.0, 0.0),
            Point3::new(x, y + 1.0, 0.0),
        ]
    };
    let mut positions = square(0.0, 0.0).to_vec();
    positions.extend(square(3.0, 0.0));
    let mesh = TriMesh::new(positions, vec![0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7]);
    let meshes = [HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh).with_planar_groups()];
    let outcome = try_faces(
        &ScalarVisibleFaces::new(),
        &meshes,
        View::new(plan_frame()),
        None,
    )
    .expect("faces");
    let found: Vec<_> = outcome
        .regions
        .iter()
        .map(|r| (r.source, r.triangles.clone(), r.outer.len()))
        .collect();
    assert_eq!(
        found,
        [
            (FaceSource::Group(FaceGroup(0)), vec![0, 1], 4),
            (FaceSource::Group(FaceGroup(2)), vec![2, 3], 4),
        ]
    );
    assert_eq!(outcome.regions[1].outer[0], Point2::new(3.0, 0.0));
}

#[test]
fn a_fold_keeps_parallel_opposite_facing_triangles_apart() {
    let mesh = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.5, 0.5, 0.0),
        ],
        vec![0, 1, 2, 1, 0, 3],
    );
    assert_eq!(reference(&mesh), [FaceGroup(0), FaceGroup(1)]);
    let meshes = [HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh).with_planar_groups()];
    let outcome = try_faces(
        &ScalarVisibleFaces::new(),
        &meshes,
        View::new(plan_frame()),
        None,
    )
    .expect("faces");
    let sources: Vec<_> = outcome.regions.iter().map(|r| r.source).collect();
    assert_eq!(
        sources,
        [
            FaceSource::Group(FaceGroup(0)),
            FaceSource::Group(FaceGroup(1))
        ]
    );
}
