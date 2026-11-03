//! How regions are assembled from the arrangement: the order of holes, and
//! the triangles of a region whose cluster of same-group triangles touches
//! none of its boundary.

mod common;

use axiolid_core::{Point2, Point3};
use axiolid_exact_arrangement::sequence_order;
use axiolid_hidden_lines::ScalarVisibleFaces;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{HiddenLineMesh, MeshId, MeshRole, OrthographicView};
use common::{faces, plan_frame};
use std::cmp::Ordering;

/// Four vertices and two triangles of an axis-parallel square at height `z`.
fn square(x0: f64, y0: f64, x1: f64, y1: f64, z: f64) -> (Vec<Point3>, Vec<u32>) {
    (
        vec![
            Point3::new(x0, y0, z),
            Point3::new(x1, y0, z),
            Point3::new(x1, y1, z),
            Point3::new(x0, y1, z),
        ],
        vec![0, 1, 2, 0, 2, 3],
    )
}

fn sheet(x0: f64, y0: f64, x1: f64, y1: f64, z: f64) -> TriMesh {
    let (vertices, indices) = square(x0, y0, x1, y1, z);
    TriMesh::new(vertices, indices)
}

fn permutations(items: &mut Vec<usize>, at: usize, out: &mut Vec<Vec<usize>>) {
    if at == items.len() {
        out.push(items.clone());
        return;
    }
    for k in at..items.len() {
        items.swap(at, k);
        permutations(items, at + 1, out);
        items.swap(at, k);
    }
}

/// The holes of a region come in ascending ring order whatever the order of
/// the occluders in the request.
#[test]
fn holes_are_sorted_whatever_the_request_order() {
    let spots = [(1.0, 1.0), (5.0, 1.0), (3.0, 4.0), (6.0, 6.0)];
    let mut orders = Vec::new();
    permutations(&mut vec![0, 1, 2, 3], 0, &mut orders);
    for order in orders {
        let mut inputs = vec![sheet(0.0, 0.0, 8.0, 8.0, 0.0)];
        for &k in &order {
            let (x, y) = spots[k];
            inputs.push(sheet(x, y, x + 1.0, y + 1.0, 1.0));
        }
        let meshes: Vec<HiddenLineMesh<'_>> = inputs
            .iter()
            .enumerate()
            .map(|(i, mesh)| HiddenLineMesh::new(MeshId(i as u64 + 1), MeshRole::Sheet, mesh))
            .collect();
        let outcome = faces(
            &ScalarVisibleFaces::new(),
            &meshes,
            OrthographicView::new(plan_frame()),
        );
        let big = outcome
            .regions
            .iter()
            .find(|r| r.mesh == MeshId(1))
            .expect("the back sheet is visible");
        assert_eq!(big.holes.len(), 4, "{order:?}");
        for pair in big.holes.windows(2) {
            assert_ne!(
                sequence_order(&pair[0], &pair[1]),
                Ordering::Greater,
                "{order:?}"
            );
        }
    }
}

/// A smaller square of the same group lying in front of a region, touching
/// none of its triangle edges, does not leave the region: its triangles are listed with the region around it,
/// also when that region is not the first one.
#[test]
fn a_cluster_inside_a_region_is_listed_with_that_region() {
    let mut vertices = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    // Triangles 0..2 A, 2..4 inside A, 4..6 B, 6..8 inside B.
    for (x0, y0, x1, y1, z) in [
        (0.0, 0.0, 2.0, 2.0, 0.0),
        (0.2, 1.2, 0.7, 1.7, 1.0),
        (10.0, 0.0, 12.0, 2.0, 0.0),
        (10.2, 1.2, 10.7, 1.7, 1.0),
    ] {
        let (v, i) = square(x0, y0, x1, y1, z);
        let base = vertices.len() as u32;
        vertices.extend(v);
        indices.extend(i.iter().map(|k| k + base));
    }
    let mesh = TriMesh::new(vertices, indices);
    let meshes = [HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh)];
    let outcome = faces(
        &ScalarVisibleFaces::new(),
        &meshes,
        OrthographicView::new(plan_frame()),
    );
    assert_eq!(outcome.regions.len(), 2);
    for (corner, expected) in [
        (Point2::new(0.0, 0.0), vec![0, 1, 2, 3]),
        (Point2::new(10.0, 0.0), vec![4, 5, 6, 7]),
    ] {
        let region = outcome
            .regions
            .iter()
            .find(|r| r.outer.contains(&corner))
            .expect("a region at the corner");
        assert_eq!(region.triangles, expected, "{corner:?}");
    }
}
