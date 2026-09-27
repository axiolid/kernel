//! Certified line of sight (#185): hidden needs a complete argument,
//! visible a checked ray, and a graze is never inverted.

use axiolid_core::Point3;
use axiolid_inspect::{line_of_sight, line_of_sight_within, Sight, SightError};
use axiolid_mesh::TriMesh;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

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
    let indices = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(positions, indices)
}

/// A rectangle in the plane `x = at`, over `y` and `z` ranges, as two
/// triangles: a wall with no thickness.
fn wall(at: f64, y: (f64, f64), z: (f64, f64)) -> TriMesh {
    TriMesh::new(
        vec![
            p(at, y.0, z.0),
            p(at, y.1, z.0),
            p(at, y.1, z.1),
            p(at, y.0, z.1),
        ],
        vec![0, 1, 2, 0, 2, 3],
    )
}

const EYE: Point3 = Point3::new(0.0, 0.0, 0.0);

#[test]
fn a_target_behind_a_wall_is_hidden() {
    let target = cuboid([10.0, -1.0, -1.0], [11.0, 1.0, 1.0]);
    let w = wall(5.0, (-2.0, 2.0), (-2.0, 2.0));
    assert_eq!(
        line_of_sight(EYE, &target, &[&w]).unwrap(),
        Sight::Hidden { occluders: vec![0] }
    );
}

#[test]
fn a_target_in_front_of_the_wall_is_visible() {
    let target = cuboid([2.0, -0.5, -0.5], [3.0, 0.5, 0.5]);
    let w = wall(5.0, (-2.0, 2.0), (-2.0, 2.0));
    match line_of_sight(EYE, &target, &[&w]).unwrap() {
        Sight::Visible { through, .. } => assert!(through.x > 1.9 && through.x < 3.1),
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_target_is_seen_past_a_column_edge() {
    // The column hides the target's left part; its right part is in view.
    let target = cuboid([10.0, -1.0, -1.0], [11.0, 3.0, 1.0]);
    let column = cuboid([5.0, -1.0, -3.0], [6.0, 1.0, 3.0]);
    match line_of_sight(EYE, &target, &[&column]).unwrap() {
        Sight::Visible { through, triangle } => {
            // The witness ray passes the column on the right.
            let at_column = through * (5.0 / through.x);
            assert!(at_column.y > 1.0 || through.x < 5.0, "{through:?}");
            assert!(triangle < 12);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_column_seen_obliquely_hides_what_is_behind_it() {
    // From the eye the column shows its front and a side; rays to the
    // target enter through the side. Neither face alone covers the view:
    // the column as a solid does.
    let column = cuboid([4.0, 1.0, -1.0], [6.0, 3.0, 1.0]);
    let target = cuboid([9.0, 2.05, -0.1], [9.2, 2.3, 0.1]);
    assert_eq!(
        line_of_sight(EYE, &target, &[&column]).unwrap(),
        Sight::Hidden { occluders: vec![0] }
    );
}

#[test]
fn a_graze_is_never_inverted() {
    // The graze plane runs through the eye and the wall's top edge:
    // z = 0.4 x. A flat target at x = 10 reaching exactly to it (z = 4)
    // is hidden or undecided -- never visible.
    let w = wall(5.0, (-2.0, 2.0), (-2.0, 2.0));
    let flat = wall(10.0, (-1.0, 1.0), (-1.0, 4.0));
    let seen = line_of_sight(EYE, &flat, &[&w]).unwrap();
    assert!(!matches!(seen, Sight::Visible { .. }), "{seen:?}");
    // A hair above the plane is visible or undecided -- never hidden.
    let above = wall(10.0, (-1.0, 1.0), (-1.0, 4.0 + 1e-9));
    let seen = line_of_sight(EYE, &above, &[&w]).unwrap();
    assert!(!matches!(seen, Sight::Hidden { .. }), "{seen:?}");
}

#[test]
fn two_walls_meeting_at_a_seam_are_not_combined() {
    // Two separate wall meshes meet edge to edge in front of the target.
    // Together they hide it, but every sub-cone across the seam straddles
    // both: undecided, which is allowed. It must not be visible.
    let left = wall(5.0, (-2.0, 0.0), (-2.0, 2.0));
    let right = wall(5.0, (0.0, 2.0), (-2.0, 2.0));
    let target = cuboid([10.0, -1.0, -1.0], [11.0, 1.0, 1.0]);
    let seen = line_of_sight_within(EYE, &target, &[&left, &right], 2000).unwrap();
    assert!(!matches!(seen, Sight::Visible { .. }), "{seen:?}");
    // One wall mesh of both rectangles is merged into convex pieces and
    // covers it.
    let both = TriMesh::new(
        vec![
            p(5.0, -2.0, -2.0),
            p(5.0, 2.0, -2.0),
            p(5.0, 2.0, 2.0),
            p(5.0, -2.0, 2.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );
    assert_eq!(
        line_of_sight(EYE, &target, &[&both]).unwrap(),
        Sight::Hidden { occluders: vec![0] }
    );
}

#[test]
fn nothing_in_the_way() {
    let target = cuboid([10.0, -1.0, -1.0], [11.0, 1.0, 1.0]);
    let none: [&TriMesh; 0] = [];
    assert!(matches!(
        line_of_sight(EYE, &target, &none).unwrap(),
        Sight::Visible { .. }
    ));
}

#[test]
fn malformed_input_is_refused() {
    let target = cuboid([10.0, -1.0, -1.0], [11.0, 1.0, 1.0]);
    let none: [&TriMesh; 0] = [];
    assert_eq!(
        line_of_sight(p(f64::NAN, 0.0, 0.0), &target, &none),
        Err(SightError::NonFinite)
    );
    let empty = TriMesh::new(vec![p(1.0, 0.0, 0.0)], Vec::new());
    assert_eq!(
        line_of_sight(EYE, &empty, &none),
        Err(SightError::EmptyTarget)
    );
}

#[test]
fn a_sliver_in_view_is_not_called_hidden() {
    // A near plate hides all of the target plate but a strip 0.02 wide
    // along its edge, too thin for the sample rays. A far wall behind the
    // target covers it from the eye but cannot hide it: it is behind.
    let target = wall(10.0, (-1.0, 1.0), (-1.0, 1.0));
    let near = wall(5.0, (-2.0, 0.49), (-2.0, 2.0));
    let far = wall(20.0, (-10.0, 10.0), (-10.0, 10.0));
    let seen = line_of_sight_within(EYE, &target, &[&near, &far], 4000).unwrap();
    assert!(!matches!(seen, Sight::Hidden { .. }), "{seen:?}");
}

#[test]
fn a_column_seen_corner_on_hides_as_a_solid() {
    // Seen along its diagonal, the column's near faces meet at an edge in
    // the middle of the view, and so do its far faces: no single face
    // covers the cones across that plane, the column as a solid does.
    let column = cuboid([4.0, 4.0, -1.0], [6.0, 6.0, 1.0]);
    let target = cuboid([9.8, 9.8, -0.2], [10.2, 10.2, 0.2]);
    assert_eq!(
        line_of_sight(EYE, &target, &[&column]).unwrap(),
        Sight::Hidden { occluders: vec![0] }
    );
}
