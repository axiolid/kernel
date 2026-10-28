//! The isometric cube with an offset origin: the view map is not exactly
//! orthonormal in binary64 and the shifted vertices are not exact in it, so
//! the interval filter cannot settle every sign; the exact tier does, and the
//! analytic drawing comes out.

mod common;

use axiolid_core::Point3;
use axiolid_mesh_hidden_lines_contract::{EdgeKind, OrthographicView, Visibility};
use common::{box_mesh, draw, isometric_frame, scalar, solids};

/// Escalations, the counts of (visible silhouettes, visible creases, hidden
/// creases) and the total number of segments.
fn analyse(origin: Point3) -> (usize, [usize; 3], usize) {
    let meshes = [box_mesh([0.0; 3], [1.0; 3])];
    let outcome = draw(
        &scalar(),
        &solids(&meshes),
        OrthographicView::new(isometric_frame(origin)),
    );
    let count = |kind, visibility| {
        outcome
            .segments
            .iter()
            .filter(|s| s.kind == kind && s.visibility == visibility)
            .count()
    };
    let side = (2.0f64 / 3.0).sqrt();
    for segment in &outcome.segments {
        let length = (segment.end - segment.start).length();
        assert!((length - side).abs() < 1e-12, "{segment:?}");
    }
    (
        outcome.evidence.exact_escalations,
        [
            count(EdgeKind::Silhouette, Visibility::Visible),
            count(EdgeKind::Crease, Visibility::Visible),
            count(EdgeKind::Crease, Visibility::Hidden),
        ],
        outcome.segments.len(),
    )
}

#[test]
fn the_offset_origin_cube_reproduces_the_analytic_drawing_with_exact_escalations() {
    let (escalations, counts, total) = analyse(Point3::new(0.1, 0.2, 0.3));
    assert_eq!(counts, [6, 3, 3]);
    assert_eq!(total, 12);
    assert!(escalations > 0, "the filter decided everything");
}

#[test]
fn the_origin_does_not_change_the_drawing() {
    let (_, shifted, _) = analyse(Point3::new(0.1, 0.2, 0.3));
    let (_, centred, _) = analyse(Point3::ZERO);
    assert_eq!(shifted, centred);
}
