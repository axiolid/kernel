// SPDX-License-Identifier: MPL-2.0

//! Conditioning of mass properties far from the world origin.
//!
//! The divergence-theorem sum evaluates `a . (b x c)` about the WORLD origin.
//! Every term scales with the cube of the coordinate magnitude, while the
//! answer stays the size of the object, so a small solid on a national grid
//! is computed as a difference of enormous nearly-equal numbers.
//!
//! Volume and centroid are translation-invariant in exact arithmetic, so
//! summing about a local origin is the same mathematics in better
//! conditioning. These tests pin that down.

use axiolid_core::{Point3, Tolerance};
use axiolid_measure::{surface_properties, volume_properties};
use axiolid_mesh::TriMesh;

/// An axis-aligned box with its minimum corner at `(ox, oy, oz)`.
fn boxx(ox: f64, oy: f64, oz: f64, s: f64) -> TriMesh {
    let positions = vec![
        Point3::new(ox, oy, oz),
        Point3::new(ox + s, oy, oz),
        Point3::new(ox + s, oy + s, oz),
        Point3::new(ox, oy + s, oz),
        Point3::new(ox, oy, oz + s),
        Point3::new(ox + s, oy, oz + s),
        Point3::new(ox + s, oy + s, oz + s),
        Point3::new(ox, oy + s, oz + s),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, // bottom
        4, 5, 6, 4, 6, 7, // top
        0, 1, 5, 0, 5, 4, // front
        1, 2, 6, 1, 6, 5, // right
        2, 3, 7, 2, 7, 6, // back
        3, 0, 4, 3, 4, 7, // left
    ];
    TriMesh::new(positions, indices)
}

/// Relative error of the measured volume against the analytic answer.
fn rel_error(ox: f64, s: f64) -> f64 {
    let mesh = boxx(ox, ox, 0.0, s);
    let measured = volume_properties(&mesh, Tolerance::METRE)
        .expect("a closed box is volume-usable")
        .signed_volume;
    let want = s * s * s;
    ((measured - want) / want).abs()
}

/// A metre-scale solid must measure correctly wherever the model sits.
///
/// Site coordinates in the millions are ordinary in infrastructure and
/// survey-derived models, so this is a realistic input, not a pathological
/// one.
#[test]
fn volume_is_accurate_at_survey_grid_coordinates() {
    for &origin in &[0.0f64, 1.0e3, 1.0e6, 1.0e7] {
        let err = rel_error(origin, 0.1);
        // The bound tracks the input's own representable grid, because that
        // is the real floor. At origin 1e7 adjacent f64 values are ~1.9e-9 m
        // apart, so a 0.1 m box's corners cannot be placed more precisely
        // than that and its volume inherits the error. Re-basing removed the
        // conditioning loss (25% -> 7e-9); what is left is the input, and
        // asserting a tighter bound would only be asserting that f64 has
        // more precision than it has.
        let ulp = if origin == 0.0 {
            f64::EPSILON
        } else {
            origin * f64::EPSILON
        };
        let floor = 10.0 * (ulp / 0.1);
        assert!(
            err < floor.max(1e-15),
            "a 0.1 m box at origin {origin:e} measured with relative error \
             {err:.3e}, above the {floor:.3e} representable-grid floor; \
             volume is translation-invariant, so this is conditioning loss"
        );
    }
}

/// The centroid must land on the solid, not somewhere else entirely.
#[test]
fn centroid_stays_on_the_solid_at_large_coordinates() {
    let origin = 1.0e7;
    let size = 0.1;
    let mesh = boxx(origin, origin, 0.0, size);
    let props = volume_properties(&mesh, Tolerance::METRE).expect("closed box");
    let want = Point3::new(origin + size / 2.0, origin + size / 2.0, size / 2.0);
    let drift = (props.centroid - want).length();
    assert!(
        drift < 1e-6,
        "centroid drifted {drift:.3e} m from the box centre at origin {origin:e}"
    );
}

/// The surface centroid must also survive large coordinates.
///
/// Honest status: this one PASSES even without re-basing, verified by
/// mutation (reverting `base` to the world origin leaves it green while the
/// two volume tests go red). Surface weights are areas (~1e-2 here), not
/// volumes, so the cancellation is far milder and f64 absorbs it at 1e7.
///
/// It is kept because the re-basing in `surface_properties` is real and
/// should not be silently reverted, and because the same sum degrades at
/// larger magnitudes than this fixture uses. It is NOT evidence that
/// `surface_properties` was broken -- it was not.
#[test]
fn surface_centroid_stays_on_the_solid_at_large_coordinates() {
    let origin = 1.0e7;
    let size = 0.1;
    let mesh = boxx(origin, origin, 0.0, size);
    let props = surface_properties(&mesh, Tolerance::METRE).expect("closed box");
    let want = Point3::new(origin + size / 2.0, origin + size / 2.0, size / 2.0);
    let drift = (props.centroid - want).length();
    assert!(
        drift < 1e-6,
        "surface centroid drifted {drift:.3e} m from the box centre at origin {origin:e}"
    );
    // Area is computed from edge differences, so it was always fine; assert
    // it anyway so a future "optimisation" cannot quietly break it.
    let want_area = 6.0 * size * size;
    let rel = ((props.area - want_area) / want_area).abs();
    assert!(rel < 1e-7, "surface area relative error {rel:.3e}");
}
