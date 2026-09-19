//! Morphology, clearance, and connected-component primitives.

mod fixtures;

use axiolid_core::{Interval, Vec3};
use axiolid_field_ops::{
    clearance_above, clearance_below, largest_free_span, sample_triangles_cpu, FieldChannel,
    LayeredFieldError, PlanarMask, Triangle3,
};

fn small_quad(w: f64, x0: f64, y0: f64, size: f64) -> [Triangle3; 2] {
    let p = |u: f64, v: f64| Vec3::new(x0 + u, y0 + v, w);
    [
        Triangle3::new(p(0.0, 0.0), p(size, 0.0), p(size, size)),
        Triangle3::new(p(0.0, 0.0), p(size, size), p(0.0, size)),
    ]
}

#[test]
fn mask_reflects_the_selected_channel() {
    let frame = fixtures::z_up_frame();
    let config = fixtures::config(frame, Vec3::new(3.0, 3.0, 6.0), 1.0);
    let field = sample_triangles_cpu(&config, &small_quad(2.0, 0.0, 0.0, 0.9)).unwrap();

    let surfaces = PlanarMask::from_field(&field, FieldChannel::SurfacePresence);
    assert_eq!(surfaces.dimensions(), (3, 3));
    assert_eq!(surfaces.count(), 1);
    assert_eq!(surfaces.get(0, 0), Some(true));
    assert_eq!(surfaces.get(2, 2), Some(false));
    assert_eq!(surfaces.get(3, 0), None);

    // No occupancy was constructed, so that channel is empty by construction.
    let occupancy = PlanarMask::from_field(&field, FieldChannel::OccupancyPresence);
    assert_eq!(occupancy.count(), 0);
}

#[test]
fn dilation_grows_by_a_metric_radius_and_erosion_reverses_it() {
    let frame = fixtures::z_up_frame();
    let config = fixtures::config(frame, Vec3::new(5.0, 5.0, 6.0), 1.0);
    let field = sample_triangles_cpu(&config, &small_quad(2.0, 2.0, 2.0, 0.9)).unwrap();
    let mask = PlanarMask::from_field(&field, FieldChannel::SurfacePresence);
    assert_eq!(mask.count(), 1);

    // Zero radius is the identity.
    assert_eq!(mask.dilate(&config, 0.0).unwrap(), mask);

    // One cell of reach: the Euclidean element gives the 4-neighbour cross.
    let grown = mask.dilate(&config, 1.0).unwrap();
    assert_eq!(grown.count(), 5);
    assert_eq!(grown.get(2, 1), Some(true));
    assert_eq!(grown.get(1, 1), Some(false), "diagonal exceeds radius 1");

    // Eroding the cross by the same radius returns the seed cell.
    let shrunk = grown.erode(&config, 1.0).unwrap();
    assert_eq!(shrunk.count(), 1);
    assert_eq!(shrunk.get(2, 2), Some(true));
}

#[test]
fn radius_is_metric_not_a_cell_count() {
    let frame = fixtures::z_up_frame();
    // Half-metre cells: a 1.0 radius must reach two cells, not one.
    let config = fixtures::config(frame, Vec3::new(5.0, 5.0, 6.0), 0.5);
    let field = sample_triangles_cpu(&config, &small_quad(2.0, 2.0, 2.0, 0.4)).unwrap();
    let mask = PlanarMask::from_field(&field, FieldChannel::SurfacePresence);
    let seed = mask.count();
    let grown = mask.dilate(&config, 1.0).unwrap();
    assert!(
        grown.count() > seed + 4,
        "a metric radius must scale with cell size, got {}",
        grown.count()
    );
}

#[test]
fn morphology_rejects_an_invalid_radius() {
    let frame = fixtures::z_up_frame();
    let config = fixtures::config(frame, Vec3::new(3.0, 3.0, 6.0), 1.0);
    let mask = PlanarMask::empty(3, 3).unwrap();
    assert_eq!(
        mask.dilate(&config, -1.0),
        Err(LayeredFieldError::InvalidEnvelope)
    );
    assert_eq!(
        mask.erode(&config, f64::NAN),
        Err(LayeredFieldError::InvalidEnvelope)
    );
}

#[test]
fn connected_components_are_labelled_deterministically() {
    let mut mask = PlanarMask::empty(5, 1).unwrap();
    mask.set(0, 0, true).unwrap();
    mask.set(1, 0, true).unwrap();
    // gap at x = 2
    mask.set(3, 0, true).unwrap();
    mask.set(4, 0, true).unwrap();

    let labels = mask.connected_components();
    assert_eq!(labels.count(), 2);
    assert!(labels.same_component((0, 0), (1, 0)));
    assert!(labels.same_component((3, 0), (4, 0)));
    assert!(!labels.same_component((1, 0), (3, 0)));
    assert_eq!(labels.label(2, 0), None);
    // Labels follow row-major discovery order, so they are stable.
    assert_eq!(labels.label(0, 0), Some(0));
    assert_eq!(labels.label(3, 0), Some(1));
}

#[test]
fn mask_set_operations_check_dimensions() {
    let a = PlanarMask::empty(2, 2).unwrap();
    let b = PlanarMask::empty(3, 2).unwrap();
    assert_eq!(a.intersect(&b), Err(LayeredFieldError::DimensionMismatch));
    assert_eq!(a.inverted().count(), 4);
    assert_eq!(a.inverted().intersect(&a).unwrap().count(), 0);
}

#[test]
fn clearance_reports_distance_to_the_next_layer() {
    let frame = fixtures::z_up_frame();
    let config = fixtures::config(frame, Vec3::new(2.0, 2.0, 8.0), 1.0);
    // Slabs at w = 1, 3, 5.
    let field = sample_triangles_cpu(&config, &fixtures::three_stacked_slabs(frame, 3.0)).unwrap();

    // Standing on the w = 1 slab: the next blocker above is w = 3.
    let above = clearance_above(&field, &config, 0, 0, 1.0).unwrap();
    assert!((above.distance - 2.0).abs() < 1e-9);
    assert_eq!(above.blocked_at, Some(3.0));
    assert!(!above.bounded_by_field);

    // Downward from the same layer: the next blocker is the field bound.
    let below = clearance_below(&field, &config, 0, 0, 1.0).unwrap();
    assert!(below.bounded_by_field);
    assert_eq!(below.blocked_at, None);
    assert!((below.distance - 1.0).abs() < 1e-9);

    // Above the topmost slab the span ends at the configured bound, not at
    // geometry: the field reports the limit rather than inventing a ceiling.
    let top = clearance_above(&field, &config, 0, 0, 5.0).unwrap();
    assert!(top.bounded_by_field);
    assert!((top.distance - 3.0).abs() < 1e-9);
}

#[test]
fn clearance_ignores_the_reference_surface_itself() {
    let frame = fixtures::z_up_frame();
    let config = fixtures::config(frame, Vec3::new(2.0, 2.0, 8.0), 1.0);
    let field = sample_triangles_cpu(&config, &fixtures::three_stacked_slabs(frame, 3.0)).unwrap();
    // Querying exactly at a slab must not report zero against that slab.
    let report = clearance_above(&field, &config, 0, 0, 3.0).unwrap();
    assert_eq!(report.blocked_at, Some(5.0));
    assert!(report.distance > 0.0);
}

#[test]
fn clearance_validates_its_inputs() {
    let frame = fixtures::z_up_frame();
    let config = fixtures::config(frame, Vec3::new(2.0, 2.0, 8.0), 1.0);
    let field = sample_triangles_cpu(&config, &fixtures::quad_at(frame, 2.0, 3.0)).unwrap();
    assert_eq!(
        clearance_above(&field, &config, 9, 9, 0.0),
        Err(LayeredFieldError::NodeOutsideField)
    );
    assert_eq!(
        clearance_above(&field, &config, 0, 0, f64::NAN),
        Err(LayeredFieldError::InvalidInterval)
    );
}

#[test]
fn largest_free_span_respects_occupancy_and_the_search_window() {
    let frame = fixtures::z_up_frame();
    let config = fixtures::config(frame, Vec3::new(2.0, 2.0, 10.0), 1.0);
    let closed = sample_triangles_cpu(&config, &fixtures::closed_slab(frame, 2.0, 4.0, 3.0))
        .unwrap()
        .derive_occupancy(axiolid_core::Tolerance::METRE)
        .unwrap();

    // Window 0..10 with 2..4 occupied: the larger free span is 4..10.
    let span = largest_free_span(&closed, 0, 0, Interval::new(0.0, 10.0))
        .unwrap()
        .unwrap();
    assert!((span.start - 4.0).abs() < 1e-9);
    assert!((span.end - 10.0).abs() < 1e-9);

    // A window entirely inside the occupied span has no free room at all.
    assert_eq!(
        largest_free_span(&closed, 0, 0, Interval::new(2.5, 3.5)).unwrap(),
        None
    );

    // An empty cell's whole window is free.
    let empty = sample_triangles_cpu(&config, &[]).unwrap();
    let free = largest_free_span(&empty, 0, 0, Interval::new(0.0, 3.0))
        .unwrap()
        .unwrap();
    assert!((free.length() - 3.0).abs() < 1e-9);
}

/// Direct transcription of the pre-optimization nested-window morphology.
///
/// Kept as the differential ORACLE for the row-span implementation in
/// `morphology.rs`. The fast path decomposes the Euclidean disc into one
/// contiguous span per row; this tests every offset explicitly. They must
/// agree bit for bit, because the span form is an exact rewrite of the
/// membership test, not an approximation of the disc.
fn windowed_reference(
    bits: &[bool],
    width: usize,
    height: usize,
    reach_squared: f64,
    steps: usize,
    dilate: bool,
) -> Vec<bool> {
    let mut out = vec![false; width * height];
    let steps = steps as isize;
    for y in 0..height {
        for x in 0..width {
            let mut value = !dilate;
            'window: for dy in -steps..=steps {
                for dx in -steps..=steps {
                    if (dx * dx + dy * dy) as f64 > reach_squared {
                        continue;
                    }
                    let nx = x as isize + dx;
                    let ny = y as isize + dy;
                    let neighbor =
                        if nx < 0 || ny < 0 || nx as usize >= width || ny as usize >= height {
                            false
                        } else {
                            bits[ny as usize * width + nx as usize]
                        };
                    if dilate && neighbor {
                        value = true;
                        break 'window;
                    }
                    if !dilate && !neighbor {
                        value = false;
                        break 'window;
                    }
                }
            }
            out[y * width + x] = value;
        }
    }
    out
}

/// The row-span morphology must equal the windowed oracle everywhere.
///
/// Sweeps grid size, radius, and pattern density, and checks BOTH
/// operations. Sizes and radii are chosen so discs land on, inside, and
/// across the border, since border handling is where the two forms could
/// most plausibly diverge.
#[test]
fn row_span_morphology_matches_the_windowed_reference() {
    let frame = fixtures::z_up_frame();
    let mut checked = 0;
    for n in [7usize, 16, 23] {
        for radius in [1.0f64, 1.5, 2.0, 3.0, 4.5] {
            for stride in [2usize, 3, 5] {
                for fill in [false, true] {
                    let config = fixtures::config(frame, Vec3::new(n as f64, n as f64, 4.0), 1.0);
                    let mut mask = PlanarMask::empty(n, n).unwrap();
                    let mut raw = vec![false; n * n];
                    // Two shapes per case. The sparse dot grid exercises dilation;
                    // the dense block exercises EROSION, which on a sparse mask
                    // short-circuits on the first unset neighbour and never
                    // reaches the border logic. `fill` inverts the roles so both
                    // branches are actually decided by the structuring element.
                    for y in 0..n {
                        for x in 0..n {
                            let set = if fill {
                                // Solid except a punched hole, so erosion has both
                                // an interior to keep and an edge to eat inward.
                                !(x >= n / 2 && y >= n / 2 && (x + y) % 5 == 0)
                            } else {
                                y % stride == 0 && x % (stride + 1) == 0
                            };
                            if set {
                                mask.set(x, y, true).unwrap();
                                raw[y * n + x] = true;
                            }
                        }
                    }
                    // Mirrors `radius_in_cells` for a unit cell size.
                    let steps = radius.ceil() as usize;
                    for dilate in [true, false] {
                        let got = if dilate {
                            mask.dilate(&config, radius)
                        } else {
                            mask.erode(&config, radius)
                        }
                        .unwrap();
                        let want = windowed_reference(&raw, n, n, radius * radius, steps, dilate);
                        for y in 0..n {
                            for x in 0..n {
                                assert_eq!(
                                got.get(x, y),
                                Some(want[y * n + x]),
                                "n={n} radius={radius} stride={stride} dilate={dilate} at ({x},{y})"
                            );
                            }
                        }
                        checked += 1;
                    }
                }
            }
        }
    }
    assert_eq!(checked, 180, "the sweep must actually run every case");
}
