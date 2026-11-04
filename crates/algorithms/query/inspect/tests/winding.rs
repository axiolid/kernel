//! Certificate that a closed oriented cycle's winding number is 0 or 1 (#420).
//!
//! Every expected excess is a closed form: the volume where `w` leaves
//! {0, 1} is (integral of w^2 - integral of w) / 2 for these shells.

use axiolid_core::{Point3, Tolerance};
use axiolid_inspect::{certify_unit_winding, VolumeInterval, WindingError};
use axiolid_mesh::{TriMesh, TriangleMeshView};

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// An axis-aligned box, outward-wound.
fn cuboid(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let positions = vec![
        p(x0, y0, z0),
        p(x1, y0, z0),
        p(x1, y1, z0),
        p(x0, y1, z0),
        p(x0, y0, z1),
        p(x1, y0, z1),
        p(x1, y1, z1),
        p(x0, y1, z1),
    ];
    let triangles: [[u32; 3]; 12] = [
        [0, 2, 1],
        [0, 3, 2],
        [4, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [1, 2, 6],
        [1, 6, 5],
        [2, 3, 7],
        [2, 7, 6],
        [3, 0, 4],
        [3, 4, 7],
    ];
    TriMesh::new(positions, triangles.iter().flatten().copied().collect())
}

fn positions(m: &TriMesh) -> Vec<Point3> {
    (0..m.position_count()).map(|i| m.position(i)).collect()
}

/// Shells in one mesh, coincident positions welded to one vertex.
fn weld(shells: &[TriMesh]) -> TriMesh {
    let mut points: Vec<Point3> = Vec::new();
    let mut indices = Vec::new();
    for shell in shells {
        for t in 0..shell.triangle_count() {
            for i in shell.triangle(t) {
                let q = shell.position(i as usize);
                let at = points.iter().position(|&r| r == q).unwrap_or_else(|| {
                    points.push(q);
                    points.len() - 1
                });
                indices.push(at as u32);
            }
        }
    }
    TriMesh::new(points, indices)
}

/// The same surface wound the other way.
fn inverted(shell: &TriMesh) -> TriMesh {
    let mut indices = Vec::new();
    for t in 0..shell.triangle_count() {
        let [a, b, c] = shell.triangle(t);
        indices.extend([a, c, b].map(|i| i as u32));
    }
    TriMesh::new(positions(shell), indices)
}

fn certify(m: &TriMesh, budget: f64) -> Result<VolumeInterval, WindingError> {
    certify_unit_winding(m, Tolerance::ZERO, budget)
}

fn exceeds(m: &TriMesh, budget: f64, exact: f64) {
    match certify(m, budget) {
        Err(WindingError::Exceeds { excess }) => {
            assert!(excess.contains(exact), "{excess:?} misses {exact}");
            assert!(excess.width() < 1e-9, "{excess:?}");
            assert!(excess.lower > budget);
        }
        other => panic!("expected Exceeds, got {other:?}"),
    }
}

#[test]
fn unit_cube_is_certified() {
    let excess = certify(&cuboid([0.0; 3], [1.0; 3]), 1e-12).expect("certified");
    assert!(excess.upper <= 1e-12 && excess.lower == 0.0, "{excess:?}");
}

#[test]
fn cubes_sharing_an_edge_are_certified() {
    let m = weld(&[
        cuboid([0.0; 3], [1.0; 3]),
        cuboid([1.0, 1.0, 0.0], [2.0, 2.0, 1.0]),
    ]);
    certify(&m, 1e-12).expect("edge-touching shells");
}

#[test]
fn a_cavity_is_certified() {
    let m = weld(&[
        cuboid([0.0; 3], [2.0; 3]),
        inverted(&cuboid([0.5; 3], [1.5; 3])),
    ]);
    certify(&m, 1e-12).expect("w is 0 inside the cavity");
}

#[test]
fn a_duplicated_cube_exceeds() {
    let c = cuboid([0.0; 3], [1.0; 3]);
    exceeds(&weld(&[c.clone(), c]), 0.5, 1.0);
}

#[test]
fn a_triple_cube_has_excess_three() {
    // w = 3: E = 3 although the volume with w outside {0, 1} is 1.
    let c = cuboid([0.0; 3], [1.0; 3]);
    exceeds(&weld(&[c.clone(), c.clone(), c]), 0.5, 3.0);
}

#[test]
fn overlapping_cubes_exceed() {
    let m = weld(&[
        cuboid([0.0; 3], [1.0; 3]),
        cuboid([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
    ]);
    exceeds(&m, 0.1, 0.5);
}

#[test]
fn an_inverted_cube_exceeds() {
    exceeds(&inverted(&cuboid([0.0; 3], [1.0; 3])), 0.5, 1.0);
}

#[test]
fn a_budget_inside_the_bracket_is_undecided() {
    let c = cuboid([0.0; 3], [1.0; 3]);
    let m = weld(&[c.clone(), c]);
    let Err(WindingError::Exceeds { excess }) = certify(&m, 1e-12) else {
        panic!("expected Exceeds");
    };
    assert!(excess.width() > 0.0);
    let budget = excess.lower + excess.width() / 2.0;
    match certify(&m, budget) {
        Err(WindingError::Undecided { excess: e }) => assert_eq!(e, excess),
        other => panic!("expected Undecided, got {other:?}"),
    }
}

#[test]
fn an_open_box_is_refused() {
    let full = cuboid([0.0; 3], [1.0; 3]);
    let kept: Vec<u32> = (2..full.triangle_count())
        .flat_map(|t| full.triangle(t).map(|i| i as u32))
        .collect();
    let open = TriMesh::new(positions(&full), kept);
    assert!(matches!(
        certify(&open, 1.0),
        Err(WindingError::NotClosedOrientedCycle(h)) if h.unbalanced_edges > 0
    ));
}

#[test]
fn a_nan_position_is_refused() {
    let full = cuboid([0.0; 3], [1.0; 3]);
    let mut points = positions(&full);
    points[3] = p(f64::NAN, 0.0, 0.0);
    let indices = (0..full.triangle_count())
        .flat_map(|t| full.triangle(t).map(|i| i as u32))
        .collect();
    let m = TriMesh::new(points, indices);
    assert_eq!(certify(&m, 1.0), Err(WindingError::NonFinite));
}

#[test]
fn an_invalid_budget_is_refused() {
    let c = cuboid([0.0; 3], [1.0; 3]);
    for budget in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(certify(&c, budget), Err(WindingError::InvalidBudget));
    }
}
