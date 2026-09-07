//! Coplanar overlap tests.
//!
//! Areas are computed by hand from each fixture. An overlap routine that
//! returns a plausible polygon with the wrong area is the failure mode
//! worth catching, so vertex counts alone are never the assertion.

use axiolid_core::Point3;
use axiolid_reference::coplanar::coplanar_overlap;

/// Area of a planar polygon, via the cross-product sum.
fn area(polygon: &[Point3]) -> f64 {
    if polygon.len() < 3 {
        return 0.0;
    }
    let mut total = Point3::new(0.0, 0.0, 0.0);
    for index in 1..polygon.len() - 1 {
        let u = polygon[index] - polygon[0];
        let v = polygon[index + 1] - polygon[0];
        total += u.cross(v);
    }
    total.length() / 2.0
}

/// Two triangles sharing a plane but no area contribute nothing.
///
/// This is the case that used to refuse whole operations: two walls in one
/// plane, metres apart. An empty overlap is the correct answer, not an error.
#[test]
fn coplanar_triangles_that_miss_each_other_have_no_overlap() {
    let left = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
    ];
    let right = [
        Point3::new(9.0, 0.0, 0.0),
        Point3::new(11.0, 0.0, 0.0),
        Point3::new(9.0, 2.0, 0.0),
    ];

    assert!(
        coplanar_overlap(left, right).is_empty(),
        "disjoint coplanar faces share no area"
    );
}

/// Identical triangles overlap in their whole area.
#[test]
fn identical_triangles_overlap_completely() {
    let triangle = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 0.0),
        Point3::new(0.0, 4.0, 0.0),
    ];

    let overlap = coplanar_overlap(triangle, triangle);

    // Legs of 4 and 4: area 8.
    assert!(
        (area(&overlap) - 8.0).abs() < 1e-9,
        "got {}",
        area(&overlap)
    );
}

/// A partial overlap has the area the geometry implies.
///
/// Both triangles are the half of a unit-ish square below its diagonal. The
/// second is shifted 2 along x, so the shared region is the triangle with
/// corners (2,0), (4,0), (2,2): legs of 2, area 2.
#[test]
fn a_partial_overlap_has_the_area_the_geometry_implies() {
    let first = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 0.0),
        Point3::new(0.0, 4.0, 0.0),
    ];
    let second = [
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(6.0, 0.0, 0.0),
        Point3::new(2.0, 4.0, 0.0),
    ];

    let overlap = coplanar_overlap(first, second);

    assert!(
        (area(&overlap) - 2.0).abs() < 1e-9,
        "got {}",
        area(&overlap)
    );
}

/// Triangles meeting along an edge share no AREA.
///
/// A shared edge is a line, not a region. Reporting it as an overlap would
/// push zero-area faces into the retriangulation downstream.
#[test]
fn edge_contact_is_not_an_area_overlap() {
    let left = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
    ];
    // Shares the edge x=2..0 along y=0 only by touching at x=2.
    let right = [
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 0.0),
        Point3::new(4.0, 2.0, 0.0),
    ];

    assert!(
        coplanar_overlap(left, right).is_empty(),
        "a shared edge or vertex is not an area"
    );
}

/// The overlap does not depend on which triangle is clipped.
///
/// Clipping is asymmetric in implementation but the shared region is not, so
/// a difference here would mean the winding normalisation is wrong.
#[test]
fn the_overlap_is_the_same_whichever_triangle_is_clipped() {
    let first = [
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(4.0, 0.0, 1.0),
        Point3::new(0.0, 4.0, 1.0),
    ];
    let second = [
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(5.0, 1.0, 1.0),
        Point3::new(1.0, 5.0, 1.0),
    ];

    let forward = area(&coplanar_overlap(first, second));
    let backward = area(&coplanar_overlap(second, first));

    assert!(forward > 0.0, "these triangles do overlap");
    assert!((forward - backward).abs() < 1e-9, "{forward} vs {backward}");
}

/// A triangle wholly inside another clips to itself.
#[test]
fn a_contained_triangle_clips_to_its_own_area() {
    let big = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(0.0, 10.0, 0.0),
    ];
    let small = [
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(3.0, 1.0, 0.0),
        Point3::new(1.0, 3.0, 0.0),
    ];

    let overlap = coplanar_overlap(small, big);

    // The small triangle has legs of 2: area 2.
    assert!(
        (area(&overlap) - 2.0).abs() < 1e-9,
        "got {}",
        area(&overlap)
    );
}
