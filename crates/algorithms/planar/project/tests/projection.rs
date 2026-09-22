//! Projection fixtures: union area, holes, degeneracy, prism clipping.
//!
//! Expected areas are computed from the geometry by hand, never recorded from
//! a previous run: a test that records its own output cannot detect a change.

use axiolid_core::{FrameError, Point2, Point3, Tolerance, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_overlay::{polygon_area, total_area, union_soup, Polygon, Ring};
use axiolid_project::{intersect_prism, project_mesh, Plane, ProjectionError};

fn tol() -> Tolerance {
    Tolerance::new(1e-9, 1e-9).expect("valid tolerance")
}

/// Two triangles forming a unit square in the z = 0 plane, lifted to `z`.
fn square_at(z: f64) -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, z),
            Point3::new(1.0, 0.0, z),
            Point3::new(1.0, 1.0, z),
            Point3::new(0.0, 1.0, z),
        ],
        vec![0, 1, 2, 0, 2, 3],
    )
}

#[test]
fn two_triangles_project_to_their_square() {
    // The two triangles tile the unit square exactly: area 1, one polygon.
    let result = project_mesh(&square_at(3.0), Plane::ground(), tol()).expect("valid projection");

    assert_eq!(result.evidence.input_triangles, 2);
    assert_eq!(result.evidence.degenerate_triangles, 0);
    assert_eq!(result.polygons.len(), 1);
    let area = total_area(&result.polygons);
    assert!((area - 1.0).abs() < 1e-9, "area was {area}");
}

#[test]
fn overlapping_triangles_are_unioned_not_summed() {
    // Two 2x2 squares at different heights, offset by (1, 1). Each has area
    // 4; they share a 1x1 overlap. Union = 4 + 4 - 1 = 7 by inclusion-exclusion,
    // computed independently of the implementation. Summing would give 8.
    let mesh = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 2.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
            Point3::new(1.0, 1.0, 5.0),
            Point3::new(3.0, 1.0, 5.0),
            Point3::new(3.0, 3.0, 5.0),
            Point3::new(1.0, 3.0, 5.0),
        ],
        vec![0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7],
    );

    let result = project_mesh(&mesh, Plane::ground(), tol()).expect("valid projection");
    let area = total_area(&result.polygons);
    assert!(
        (area - 7.0).abs() < 1e-9,
        "expected union area 7, got {area}"
    );
}

#[test]
fn a_through_hole_survives_projection() {
    // An annulus: outer 4x4 square, inner 2x2 hole, triangulated as four
    // trapezoidal bands. Area = 16 - 4 = 12, and the hole must be preserved
    // rather than filled by a hull or an outline.
    let mesh = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
            Point3::new(4.0, 4.0, 0.0),
            Point3::new(0.0, 4.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
            Point3::new(3.0, 3.0, 0.0),
            Point3::new(1.0, 3.0, 0.0),
        ],
        vec![
            0, 1, 5, 0, 5, 4, // bottom band
            1, 2, 6, 1, 6, 5, // right band
            2, 3, 7, 2, 7, 6, // top band
            3, 0, 4, 3, 4, 7, // left band
        ],
    );

    let result = project_mesh(&mesh, Plane::ground(), tol()).expect("valid projection");
    let area = total_area(&result.polygons);

    assert!(
        (area - 12.0).abs() < 1e-9,
        "expected annulus area 12, got {area}"
    );
    assert_eq!(result.evidence.output_holes, 1, "the hole must survive");
}

#[test]
fn edge_on_triangles_are_counted_not_silently_dropped() {
    // A vertical square seen from directly above projects to a line: zero
    // area. The distinction that matters is between "this mesh is edge-on"
    // and "this mesh was empty", so the count is reported.
    let wall = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 2.0),
            Point3::new(0.0, 0.0, 2.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );

    let result = project_mesh(&wall, Plane::ground(), tol()).expect("valid projection");

    assert_eq!(result.evidence.input_triangles, 2);
    assert_eq!(
        result.evidence.degenerate_triangles, 2,
        "both faces are edge-on"
    );
    assert!(
        result.polygons.is_empty(),
        "an edge-on wall has no footprint"
    );
}

#[test]
fn a_closed_solid_does_not_cancel_its_own_footprint() {
    // Top and bottom faces of a box: the bottom is wound the opposite way in
    // 3D. Orientation is not a planar fact, so both must contribute
    // positively. If back faces were subtracted the footprint would vanish.
    let slab = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 1.0),
        ],
        // bottom wound CW seen from above, top wound CCW
        vec![0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7],
    );

    let result = project_mesh(&slab, Plane::ground(), tol()).expect("valid projection");
    let area = total_area(&result.polygons);
    assert!(
        (area - 1.0).abs() < 1e-9,
        "expected footprint 1, got {area}"
    );
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon {
        outer: Ring {
            points: vec![
                Point2::new(x0, y0),
                Point2::new(x1, y0),
                Point2::new(x1, y1),
                Point2::new(x0, y1),
            ],
        },
        holes: Vec::new(),
    }
}

#[test]
fn a_prism_clips_the_footprint_to_its_region() {
    // The unit square footprint clipped by the prism over [0.5, 2] x [0.5, 2]
    // leaves [0.5, 1] x [0.5, 1]: area 0.25.
    let region = [rect(0.5, 0.5, 2.0, 2.0)];
    let result = intersect_prism(&square_at(0.0), Plane::ground(), &region, tol())
        .expect("valid prism intersection");

    let area = total_area(&result.polygons);
    assert!(
        (area - 0.25).abs() < 1e-9,
        "expected clipped area 0.25, got {area}"
    );
}

#[test]
fn a_disjoint_prism_yields_an_empty_footprint() {
    let region = [rect(10.0, 10.0, 12.0, 12.0)];
    let result = intersect_prism(&square_at(0.0), Plane::ground(), &region, tol())
        .expect("a disjoint prism is a valid query");
    assert!(result.polygons.is_empty());
}

/// A degenerate basis cannot reach `project_mesh` at all.
///
/// `Plane` is now a validated type, so the refusal happens at construction
/// rather than inside the projection. That is the improvement: an unvalidated
/// plane is no longer representable, so the failure cannot be forgotten by a
/// future caller.
#[test]
fn a_degenerate_basis_cannot_construct_a_plane() {
    let error = Plane::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol(),
    )
    .expect_err("two identical axes span a line, not a plane");
    assert_eq!(error, FrameError::Degenerate);
}

#[test]
fn an_out_of_range_index_is_refused() {
    let broken = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        vec![0, 1, 9],
    );
    let error = project_mesh(&broken, Plane::ground(), tol()).expect_err("index 9 does not exist");
    assert_eq!(error, ProjectionError::IndexOutOfRange);
}

/// Two conflicting "footprint" rules compose over one projection (ADR 0066).
///
/// The kernel deliberately ships no `Footprint` type, because the word does
/// not denote a fixed set of points: consumers disagree about whether an
/// overhang counts. This test stands in for those consumers and shows the
/// disagreement is expressible downstream, without either rule living in the
/// kernel and without re-running the geometry.
///
/// If a future change pushed one inclusion rule into `project_mesh`, one of
/// these two expectations would become unreachable -- which is exactly the
/// failure the ADR exists to prevent.
#[test]
fn conflicting_inclusion_rules_compose_over_one_projection() {
    // A 2x1 slab at z = 0 with a 1x1 overhang at z = 1, offset in +x.
    let base = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );
    let overhang = TriMesh::new(
        vec![
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(2.0, 0.0, 1.0),
            Point3::new(2.0, 1.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );

    let plane = Plane::ground();
    let base_projection = project_mesh(&base, plane, tol()).expect("base projects");
    let overhang_projection = project_mesh(&overhang, plane, tol()).expect("overhang projects");

    // Rule A -- "structure touching the ground": the overhang is excluded.
    let excluding = total_area(&base_projection.polygons);
    assert!(
        (excluding - 1.0).abs() < 1e-9,
        "grounded-only rule should see 1.0, saw {excluding}"
    );

    // Rule B -- "anything covering the site": the overhang is included. Both
    // rules read the same projections and differ only in what they combine.
    //
    // Unioned rather than area-summed. `total_area` adds polygon areas, so
    // summing would double-count wherever parts overlap and would report 2.0
    // even for an overhang sitting directly above the base. A footprint is the
    // area covered, not the area of material.
    let combined: Vec<Ring> = base_projection
        .polygons
        .iter()
        .chain(overhang_projection.polygons.iter())
        .map(|polygon| polygon.outer.clone())
        .collect();
    let merged = union_soup(&combined, tol()).expect("projections union");
    let including = total_area(&merged);
    assert!(
        (including - 2.0).abs() < 1e-9,
        "covering rule should see 2.0, saw {including}"
    );

    // The point of the ADR: the geometry is identical, the answers are not.
    assert!(
        including > excluding,
        "the two rules must be able to disagree"
    );

    // And an overhang stacked directly above the base adds no covered area,
    // which is what makes this a union and not a sum.
    let stacked = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 1.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );
    let stacked_projection = project_mesh(&stacked, plane, tol()).expect("stack projects");
    let overlapping: Vec<Ring> = base_projection
        .polygons
        .iter()
        .chain(stacked_projection.polygons.iter())
        .map(|polygon| polygon.outer.clone())
        .collect();
    let overlapped = union_soup(&overlapping, tol()).expect("overlap unions");
    let covered = total_area(&overlapped);
    assert!(
        (covered - 1.0).abs() < 1e-9,
        "stacked geometry covers 1.0, not 2.0; saw {covered}"
    );
}

/// A projection is reachable as `Polygon2`, and the conversion is lossy.
///
/// `project_mesh` returns overlay's `Polygon` rather than the foundation's
/// `Polygon2` because a projection can have holes and `Polygon2` models a
/// simple polygon. This pins both halves of that reasoning: the bridge exists,
/// and taking it costs the hole.
///
/// Written as a test rather than a doc line because "the conversion is lossy"
/// is the kind of claim that quietly stops being true.
#[test]
fn a_projection_converts_to_polygon2_but_an_outline_drops_the_hole() {
    // The annulus from `a_through_hole_survives_projection`: 4x4 outer with a
    // 2x2 hole, so covered area is 16 - 4 = 12 while the outline encloses 16.
    let mesh = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
            Point3::new(4.0, 4.0, 0.0),
            Point3::new(0.0, 4.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
            Point3::new(3.0, 3.0, 0.0),
            Point3::new(1.0, 3.0, 0.0),
        ],
        vec![
            0, 1, 5, 0, 5, 4, // bottom band
            1, 2, 6, 1, 6, 5, // right band
            2, 3, 7, 2, 7, 6, // top band
            3, 0, 4, 3, 4, 7, // left band
        ],
    );

    let result = project_mesh(&mesh, Plane::ground(), tol()).expect("valid projection");
    let polygon = result
        .polygons
        .first()
        .expect("the annulus projects to one polygon");
    assert!(polygon.has_holes(), "the fixture is meant to have a hole");

    // The polygon knows its true covered area, hole subtracted.
    assert!(
        (polygon_area(polygon) - 12.0).abs() < 1e-9,
        "covered area is 12, got {}",
        polygon_area(polygon)
    );

    // The outline is a Polygon2, and it has filled the hole in. 16, not 12.
    let outline = polygon.outline();
    assert!(
        (outline.area() - 16.0).abs() < 1e-9,
        "the outline encloses 16, got {}",
        outline.area()
    );
    assert!(
        outline.area() > polygon_area(polygon),
        "if these ever agree the conversion stopped being lossy and the \
         return type should be revisited"
    );

    // Round-tripping a ring through Polygon2 is, by contrast, exact.
    let restored = Ring::from(outline.clone());
    assert_eq!(restored, polygon.outer, "ring round trip must not drift");
}
