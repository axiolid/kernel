//! Stepped union: bands against closed-form volume.
//!
//! A stepped union is refused by the single-prism path. The decomposition
//! is checked by TOTAL VOLUME against the inclusion-exclusion value
//! computed independently, so a dropped or duplicated band shows up as a
//! wrong number rather than a plausible-looking band list.

use axiolid_construct::boolean_exact::Prism;
use axiolid_construct::boolean_stepped::{union_prisms_stepped, Band};
use axiolid_core::{Point2, Tolerance};

fn square(half: f64) -> Vec<Point2> {
    vec![
        Point2::new(-half, -half),
        Point2::new(half, -half),
        Point2::new(half, half),
        Point2::new(-half, half),
    ]
}

fn prism(half: f64, bottom: f64, top: f64) -> Prism {
    Prism {
        rings: vec![square(half)],
        bottom,
        top,
    }
}

/// Shoelace area of a band's outer ring minus its holes.
fn band_area(band: &Band) -> f64 {
    let mut total = 0.0;
    for (index, ring) in band.rings.iter().enumerate() {
        let mut area = 0.0;
        for i in 0..ring.len() {
            let a = ring[i];
            let b = ring[(i + 1) % ring.len()];
            area += a.x * b.y - b.x * a.y;
        }
        area = 0.5 * area.abs();
        if index == 0 {
            total += area;
        } else {
            total -= area;
        }
    }
    total
}

fn volume(bands: &[Band]) -> f64 {
    bands
        .iter()
        .map(|band| band_area(band) * (band.top - band.bottom))
        .sum()
}

#[test]
fn a_stepped_union_has_the_inclusion_exclusion_volume() {
    // A wide short slab with a narrow tall tower standing on it: the
    // classic stepped shape the single-prism path refuses.
    let wide = prism(2.0, 0.0, 1.0);
    let narrow = prism(0.5, 0.0, 3.0);
    let bands = union_prisms_stepped(&wide, &narrow, Tolerance::METRE)
        .expect("a touching stepped union is representable as bands");

    // Two sections means two bands: the slab, then the tower alone.
    assert_eq!(bands.len(), 2, "expected one band per section change");

    // Volume by inclusion-exclusion, computed independently of the bands:
    //   wide 4x4x1 = 16, narrow 1x1x3 = 3, overlap 1x1x1 = 1.
    let expected = 16.0 + 3.0 - 1.0;
    let actual = volume(&bands);
    assert!(
        (actual - expected).abs() < 1e-9,
        "stepped union volume {actual} should equal {expected}"
    );
}

#[test]
fn bands_tile_the_full_height_without_gap_or_overlap() {
    let lower = prism(1.0, 0.0, 2.0);
    let upper = prism(1.5, 1.0, 4.0);
    let bands = union_prisms_stepped(&lower, &upper, Tolerance::METRE)
        .expect("overlapping spans are representable");

    // Three sections: lower alone, both, upper alone.
    assert_eq!(bands.len(), 3);
    assert!((bands[0].bottom - 0.0).abs() < 1e-12);
    assert!((bands[bands.len() - 1].top - 4.0).abs() < 1e-12);
    for pair in bands.windows(2) {
        assert!(
            (pair[0].top - pair[1].bottom).abs() < 1e-12,
            "a gap or overlap between bands loses or double-counts material"
        );
    }
    for band in &bands {
        assert!(band.top > band.bottom, "a band must have positive height");
    }
}

#[test]
fn equal_spans_collapse_to_a_single_band() {
    // This is the case the single-prism path already handles, so the
    // decomposition must not invent a step where there is none.
    let a = prism(1.0, 0.0, 2.0);
    let b = prism(1.5, 0.0, 2.0);
    let bands = union_prisms_stepped(&a, &b, Tolerance::METRE).expect("equal spans union");
    assert_eq!(bands.len(), 1, "equal spans are one prism, not a stack");
}

#[test]
fn prisms_that_do_not_meet_are_refused() {
    let low = prism(1.0, 0.0, 1.0);
    let high = prism(1.0, 5.0, 6.0);
    // Two separate solids: a band stack describes one, so this must refuse
    // rather than silently bridging the gap with material.
    assert!(union_prisms_stepped(&low, &high, Tolerance::METRE).is_err());
}

fn rect(hx: f64, hy: f64) -> Vec<Point2> {
    vec![
        Point2::new(-hx, -hy),
        Point2::new(hx, -hy),
        Point2::new(hx, hy),
        Point2::new(-hx, hy),
    ]
}

#[test]
fn an_overlap_band_is_the_union_of_both_sections_not_either_one() {
    // Cross-shaped overlap: neither bar contains the other, so a band that
    // copied just one operand's section would have the wrong area. The
    // earlier tests could not see this because the tower sat inside the slab.
    let across = Prism {
        rings: vec![rect(2.0, 0.5)],
        bottom: 0.0,
        top: 2.0,
    };
    let upright = Prism {
        rings: vec![rect(0.5, 2.0)],
        bottom: 1.0,
        top: 3.0,
    };
    let bands = union_prisms_stepped(&across, &upright, Tolerance::METRE)
        .expect("crossing bars meet and are representable");
    assert_eq!(bands.len(), 3);

    // Each bar is 4 x 1 = 4; their cross overlap is 1 x 1 = 1, so the
    // cross section area is 4 + 4 - 1 = 7 -- larger than either bar alone.
    let middle = &bands[1];
    let area = band_area(middle);
    assert!(
        (area - 7.0).abs() < 1e-9,
        "the overlap band must be the union of both sections, got area {area}"
    );

    // Total: bar1 4x1x2=8, bar2 4x1x2=8, shared cross 1x1x1=1.
    let expected = 8.0 + 8.0 - 1.0;
    let actual = volume(&bands);
    assert!(
        (actual - expected).abs() < 1e-9,
        "crossing stepped union volume {actual} should equal {expected}"
    );
}

#[test]
fn heights_within_tolerance_are_one_cut_not_a_sliver_band() {
    // The tops differ by far less than the tolerance, so they are the same
    // height. Treating them as two cuts would emit a band of ~1e-15
    // thickness: no solid can carry that, and it is not real geometry.
    let a = prism(1.0, 0.0, 2.0);
    let b = prism(1.5, 0.0, 2.0 + 1e-15);
    let bands = union_prisms_stepped(&a, &b, Tolerance::METRE).expect("near-equal spans");
    assert_eq!(
        bands.len(),
        1,
        "heights within tolerance must merge into one cut"
    );
    for band in &bands {
        assert!(
            band.top - band.bottom > 1e-9,
            "a sliver band of thickness {} is not representable",
            band.top - band.bottom
        );
    }
}
