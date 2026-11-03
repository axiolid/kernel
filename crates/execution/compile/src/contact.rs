//! Bodies of a union's result that lie within the tolerance (#356, #357).
//!
//! The mesh boolean has no tolerance (ADR 0047): it merges faces only where
//! they are exactly coplanar. Two operands that meet within the linear
//! tolerance, but not exactly, can therefore come back as two bodies, closed
//! and consistently wound, which the closure check labels a solid. The exact
//! compiler reads the same operands as one body (ADR 0080), so the two
//! compilers would disagree about a union the caller declared in contact.
//!
//! [`unmerged`] finds that outcome on the result: the first pair of index
//! components whose surfaces lie within the tolerance of each other. The
//! caller refuses the union by name instead of returning it. Reading the
//! result needs no point-in-mesh containment between the operands, and it
//! also sees a union that leaves a cavity open to the outside by less than
//! the tolerance.
//!
//! The cost is nothing for the common case, a result of one component.
//! Otherwise pairs of components whose boxes, grown by the tolerance, meet
//! are measured with [`axiolid_inspect::surface_gap`].
//!
//! What it does not see: one body whose faces cross with sub-rounding cracks
//! where contact faces met under different placements. That stays open
//! (`TODO(#391)`).

use axiolid_contracts::GeomError;
use axiolid_core::{Aabb, Point3, Scalar};
use axiolid_mesh::{decompose, TriMesh};

/// The start of the `GeomError::Degenerate` message a union is refused with
/// when its result holds two bodies within the linear tolerance of each other
/// (#356, #357): the operands meet within the tolerance, and the mesh boolean
/// left them apart.
pub const UNMERGED_UNION_CONTACT: &str =
    "union of operands that meet within the linear tolerance, left as separate bodies";

/// Two bodies of a union's result within `eps` of each other.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Unmerged {
    /// The distance between their surfaces.
    pub distance: Scalar,
    /// The first body's bounding box, as minimum and maximum corners.
    pub first: (Point3, Point3),
    /// The second body's bounding box.
    pub second: (Point3, Point3),
}

/// The first pair of index components of `mesh` (in
/// [`axiolid_mesh::decompose`] order, pairs `i < j` lexicographic) whose
/// surfaces lie within `eps`; `None` for fewer than two components or none
/// that close.
pub(crate) fn unmerged(mesh: &TriMesh, eps: Scalar) -> Option<Unmerged> {
    let parts = decompose(mesh);
    if parts.len() < 2 {
        return None;
    }
    let boxes: Vec<Aabb> = parts.iter().map(bounds).collect();
    for i in 0..parts.len() {
        for j in i + 1..parts.len() {
            if !close(&boxes[i], &boxes[j], eps) {
                continue;
            }
            if let Some(distance) = axiolid_inspect::surface_gap(&parts[i], &parts[j], eps) {
                return Some(Unmerged {
                    distance,
                    first: (boxes[i].min, boxes[i].max),
                    second: (boxes[j].min, boxes[j].max),
                });
            }
        }
    }
    None
}

fn bounds(mesh: &TriMesh) -> Aabb {
    let mut points = mesh.positions.iter();
    let mut bounds = Aabb::from_point(points.next().copied().unwrap_or(Point3::ZERO));
    for point in points {
        bounds.extend(*point);
    }
    bounds
}

/// Whether the boxes, grown by `eps`, meet.
fn close(a: &Aabb, b: &Aabb, eps: Scalar) -> bool {
    a.min.x - eps <= b.max.x
        && b.min.x - eps <= a.max.x
        && a.min.y - eps <= b.max.y
        && b.min.y - eps <= a.max.y
        && a.min.z - eps <= b.max.z
        && b.min.z - eps <= a.max.z
}

/// The refusal of a union that left `found` apart, within `eps`.
pub(crate) fn refusal(found: &Unmerged, eps: Scalar) -> GeomError {
    GeomError::Degenerate(format!(
        "{UNMERGED_UNION_CONTACT}: two bodies {:e} apart, within the linear tolerance {eps:e} \
         (boxes {:?} and {:?})",
        found.distance, found.first, found.second
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cube(lo: Point3, hi: Point3) -> TriMesh {
        let p = |x: bool, y: bool, z: bool| {
            Point3::new(
                if x { hi.x } else { lo.x },
                if y { hi.y } else { lo.y },
                if z { hi.z } else { lo.z },
            )
        };
        let positions = vec![
            p(false, false, false),
            p(true, false, false),
            p(true, true, false),
            p(false, true, false),
            p(false, false, true),
            p(true, false, true),
            p(true, true, true),
            p(false, true, true),
        ];
        let indices = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7,
            6, 3, 0, 4, 3, 4, 7,
        ];
        TriMesh::new(positions, indices)
    }

    fn unit(x: f64, y: f64, z: f64) -> TriMesh {
        cube(Point3::new(x, y, z), Point3::new(x + 1.0, y + 1.0, z + 1.0))
    }

    /// Two meshes as one, indices shifted, as a union result holds two bodies.
    fn together(a: &TriMesh, b: &TriMesh) -> TriMesh {
        let mut positions = a.positions.clone();
        let offset = positions.len() as u32;
        positions.extend(b.positions.iter().copied());
        let mut indices = a.indices.clone();
        indices.extend(b.indices.iter().map(|i| i + offset));
        TriMesh::new(positions, indices)
    }

    const EPS: f64 = 1e-6;

    #[test]
    fn one_body_is_never_unmerged() {
        assert_eq!(unmerged(&unit(0.0, 0.0, 0.0), EPS), None);
    }

    #[test]
    fn two_cubes_half_a_tolerance_apart_are_unmerged() {
        let mesh = together(&unit(0.0, 0.0, 0.0), &unit(1.0 + 0.5e-6, 0.0, 0.0));
        let found = unmerged(&mesh, EPS).expect("within the tolerance");
        assert!(
            (found.distance - 0.5e-6).abs() < 1e-15,
            "{}",
            found.distance
        );
        assert_eq!(found.first.1.x, 1.0);
        assert_eq!(found.second.0.x, 1.0 + 0.5e-6);
    }

    #[test]
    fn two_cubes_two_tolerances_apart_are_not() {
        let mesh = together(&unit(0.0, 0.0, 0.0), &unit(1.0 + 2e-6, 0.0, 0.0));
        assert_eq!(unmerged(&mesh, EPS), None);
    }

    #[test]
    fn cubes_touching_along_an_edge_are_unmerged() {
        let mesh = together(&unit(0.0, 0.0, 0.0), &unit(1.0, 1.0, 0.0));
        let found = unmerged(&mesh, EPS).expect("touching");
        assert_eq!(found.distance, 0.0);
    }

    #[test]
    fn a_body_in_a_cavity_further_than_the_tolerance_is_not() {
        let outer = cube(Point3::ZERO, Point3::new(10.0, 10.0, 10.0));
        let mut void = cube(Point3::new(2.0, 2.0, 2.0), Point3::new(8.0, 8.0, 8.0));
        void.indices
            .as_chunks_mut::<3>()
            .0
            .iter_mut()
            .for_each(|triangle| triangle.swap(1, 2));
        let island = cube(Point3::new(3.0, 3.0, 3.0), Point3::new(4.0, 4.0, 4.0));
        let mesh = together(&together(&outer, &void), &island);
        assert_eq!(unmerged(&mesh, EPS), None);
    }

    #[test]
    fn zero_tolerance_finds_only_touching() {
        let touching = together(&unit(0.0, 0.0, 0.0), &unit(1.0, 0.0, 0.0));
        assert!(unmerged(&touching, 0.0).is_some());
        let apart = together(&unit(0.0, 0.0, 0.0), &unit(1.0 + 1e-9, 0.0, 0.0));
        assert_eq!(unmerged(&apart, 0.0), None);
    }
}
