//! Parameterised structural section extrusion (ADR 0057).
//!
//! Areas are checked against closed forms derived independently, including
//! the root-fillet contribution `r^2 (1 - pi/4)` per fillet. That term was
//! verified by Monte-Carlo integration before being used here.

use axiolid_brep_audit::geometric_audit;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::section_lower::section_contour;
use axiolid_core::{Tolerance, Vec3};
use axiolid_curve::Curve2;
use axiolid_profile::{Profile, SectionProfile};
use axiolid_surface::Surface;

const DEPTH: f64 = 2.0;

/// Area a single fillet adds over the sharp corner it replaces.
fn fillet_area(radius: f64) -> f64 {
    radius * radius * (1.0 - core::f64::consts::FRAC_PI_4)
}

/// Cross-sectional area via Green's theorem over the EXACT contour.
///
/// Integrating the contour directly keeps the arcs exact: a straight segment
/// contributes its trapezoid and a circular arc contributes its exact
/// sector-plus-triangle term, so the fillet material is measured rather than
/// chord-approximated. It shares no code with `exact_properties`, which
/// [`volume_of`] checks the built solid against.
fn contour_area(profile: &Profile) -> f64 {
    let section = match profile {
        Profile::Section(section) => section,
        _ => unreachable!("only sections here"),
    };
    let contour = section_contour(section).expect("a section lowers");
    let mut total = 0.0;
    for segment in &contour.outer.segments {
        let (from, to) = if segment.same_sense {
            (segment.domain.start, segment.domain.end)
        } else {
            (segment.domain.end, segment.domain.start)
        };
        match &segment.curve {
            Curve2::Line(line) => {
                let a = line.origin + line.direction * from;
                let b = line.origin + line.direction * to;
                total += a.perp_dot(b);
            }
            Curve2::Circle(circle) => {
                // Green's theorem over a circular arc, in closed form.
                let point = |t: f64| {
                    let (s, c) = t.sin_cos();
                    circle.frame.origin
                        + circle.frame.x * (circle.radius * c)
                        + circle.frame.y * (circle.radius * s)
                };
                let handedness = circle.frame.x.perp_dot(circle.frame.y);
                let sweep = (to - from) * handedness.signum();
                let a = point(from);
                let b = point(to);
                // Chord term plus the circular segment the arc bulges over.
                total += a.perp_dot(b);
                total += circle.radius * circle.radius * (sweep - sweep.sin());
            }
            other => panic!("unexpected segment kind {other:?}"),
        }
    }
    total / 2.0
}

fn volume_of(profile: &Profile) -> f64 {
    let solid = extrude_profile_exact(profile, Vec3::Z, DEPTH, Tolerance::METRE)
        .expect("a section extrudes");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "section solid must audit clean, found {:?}",
        health.defects()
    );
    // Area from the exact contour, depth from the extrusion: the solid is a
    // prism, so this IS its volume, and it works for curved walls too. The
    // built solid, fillet walls and all, must measure the same (#125).
    let expected = contour_area(profile) * DEPTH;
    let measured = axiolid_measure::exact_properties(&solid, Tolerance::METRE)
        .expect("a section solid is measurable")
        .signed_volume;
    assert!(
        (measured - expected).abs() <= 1e-11 * expected.abs(),
        "solid measures {measured}, contour gives {expected}"
    );
    expected
}

fn i_section(fillet: Option<f64>) -> Profile {
    Profile::Section(SectionProfile::I {
        depth: 0.4,
        width: 0.3,
        web_thickness: 0.011,
        flange_thickness: 0.019,
        fillet_radius: fillet,
        flange_edge_radius: None,
        flange_slope: None,
    })
}

/// The same I with an optional flange taper.
fn i_sloped(fillet: Option<f64>, slope: Option<f64>) -> Profile {
    Profile::Section(SectionProfile::I {
        depth: 0.4,
        width: 0.3,
        web_thickness: 0.011,
        flange_thickness: 0.019,
        fillet_radius: fillet,
        flange_edge_radius: None,
        flange_slope: slope,
    })
}

/// The `Profile::Section` payload, or a panic: these helpers only take
/// sections, so anything else is a test-authoring mistake, not a case to
/// handle.
fn section_of(profile: &Profile) -> &SectionProfile {
    match profile {
        Profile::Section(section) => section,
        _ => panic!("expected a section profile"),
    }
}

#[test]
fn a_sharp_i_section_matches_its_closed_form_area() {
    // 2*b*tf + (d - 2*tf)*tw, with no fillet material.
    let (d, b, tw, tf) = (0.4, 0.3, 0.011, 0.019);
    let expected = (2.0 * b * tf + (d - 2.0 * tf) * tw) * DEPTH;
    let got = volume_of(&i_section(None));
    assert!(
        (got - expected).abs() < 1e-12,
        "expected {expected}, got {got}"
    );
}

#[test]
fn the_root_fillets_carry_real_area() {
    // Measured: for this rolled section the four root fillets are 2.40% of
    // the cross-section. Dropping them leaves a shape that still looks like
    // an I and whose area, second moment and mass are all wrong.
    let (d, b, tw, tf, r) = (0.4, 0.3, 0.011, 0.019, 0.021);
    let sharp = (2.0 * b * tf + (d - 2.0 * tf) * tw) * DEPTH;
    let expected = sharp + 4.0 * fillet_area(r) * DEPTH;

    let got = volume_of(&i_section(Some(r)));
    assert!(
        (got - expected).abs() < 1e-12,
        "expected {expected}, got {got}"
    );

    // The claim that makes this worth testing at all.
    let share = (expected - sharp) / expected;
    assert!(
        (share - 0.024).abs() < 0.001,
        "root fillets should be about 2.4% of the section, got {share}"
    );
    assert!(got > sharp, "fillets must ADD material, not remove it");
}

#[test]
fn a_filleted_section_carries_exact_cylindrical_walls() {
    // Four root fillets become four genuine cylindrical walls, not a fan of
    // planar strips.
    let solid = extrude_profile_exact(&i_section(Some(0.021)), Vec3::Z, DEPTH, Tolerance::METRE)
        .expect("a filleted I extrudes");
    let radii: Vec<f64> = solid
        .surfaces()
        .iter()
        .filter_map(|s| match s {
            Surface::Cylinder(c) => Some(c.radius),
            _ => None,
        })
        .collect();
    assert_eq!(radii.len(), 4, "four root fillets, got {radii:?}");
    for radius in &radii {
        assert!(
            (radius - 0.021).abs() < 1e-12,
            "each fillet must carry the stated radius, got {radius}"
        );
    }
}

#[test]
fn a_t_section_matches_its_closed_form_area() {
    let (d, bf, tw, tf, r) = (0.3, 0.2, 0.01, 0.015, 0.012);
    let profile = Profile::Section(SectionProfile::T {
        depth: d,
        flange_width: bf,
        web_thickness: tw,
        flange_thickness: tf,
        fillet_radius: Some(r),
        flange_edge_radius: None,
        web_edge_radius: None,
        web_slope: None,
        flange_slope: None,
    });
    // Flange plus web, plus TWO root fillets.
    let expected = (bf * tf + (d - tf) * tw + 2.0 * fillet_area(r)) * DEPTH;
    let got = volume_of(&profile);
    assert!(
        (got - expected).abs() < 1e-12,
        "expected {expected}, got {got}"
    );
}

#[test]
fn a_u_section_matches_its_closed_form_area() {
    let (d, bf, tw, tf, r) = (0.3, 0.1, 0.0075, 0.0125, 0.01);
    let profile = Profile::Section(SectionProfile::U {
        depth: d,
        flange_width: bf,
        web_thickness: tw,
        flange_thickness: tf,
        fillet_radius: Some(r),
        edge_radius: None,
        flange_slope: None,
    });
    // Web over the full depth, plus the two flange outstands.
    let expected = (d * tw + 2.0 * (bf - tw) * tf + 2.0 * fillet_area(r)) * DEPTH;
    let got = volume_of(&profile);
    assert!(
        (got - expected).abs() < 1e-12,
        "expected {expected}, got {got}"
    );
}

#[test]
fn an_l_section_matches_its_closed_form_area() {
    let (d, w, t, r) = (0.15, 0.1, 0.012, 0.009);
    let profile = Profile::Section(SectionProfile::L {
        depth: d,
        width: Some(w),
        thickness: t,
        fillet_radius: Some(r),
        edge_radius: None,
        leg_slope: None,
    });
    // Two legs sharing the heel square, plus ONE root fillet.
    let expected = (w * t + (d - t) * t + fillet_area(r)) * DEPTH;
    let got = volume_of(&profile);
    assert!(
        (got - expected).abs() < 1e-12,
        "expected {expected}, got {got}"
    );
}

#[test]
fn an_absent_l_width_means_an_equal_angle() {
    // `None` means the source did not state a width, which for an angle
    // means equal legs -- not a zero-width leg.
    let (d, t) = (0.12, 0.01);
    let implied = Profile::Section(SectionProfile::L {
        depth: d,
        width: None,
        thickness: t,
        fillet_radius: None,
        edge_radius: None,
        leg_slope: None,
    });
    let stated = Profile::Section(SectionProfile::L {
        depth: d,
        width: Some(d),
        thickness: t,
        fillet_radius: None,
        edge_radius: None,
        leg_slope: None,
    });
    let a = volume_of(&implied);
    let b = volume_of(&stated);
    assert!((a - b).abs() < 1e-15, "{a} vs {b}");
}

#[test]
fn a_trapezium_matches_its_closed_form_area() {
    let (bottom, top, y, offset) = (0.4, 0.25, 0.2, 0.05);
    let profile = Profile::Section(SectionProfile::Trapezium {
        bottom_x: bottom,
        top_x: top,
        y,
        top_offset: offset,
    });
    // The offset shears the shape; area is unchanged by shear.
    let expected = (0.5 * (bottom + top) * y) * DEPTH;
    let got = volume_of(&profile);
    assert!(
        (got - expected).abs() < 1e-12,
        "expected {expected}, got {got}"
    );
}

#[test]
fn an_asymmetric_i_keeps_its_two_flange_widths() {
    // The whole reason the variant exists: collapsing the widths moves the
    // neutral axis. Different widths must give a different area.
    let make = |top: f64| {
        Profile::Section(SectionProfile::AsymmetricI {
            depth: 0.4,
            web_thickness: 0.011,
            bottom_flange_width: 0.3,
            bottom_flange_thickness: 0.019,
            bottom_fillet_radius: None,
            bottom_flange_edge_radius: None,
            bottom_flange_slope: None,
            top_flange_width: top,
            top_flange_thickness: Some(0.019),
            top_fillet_radius: None,
            top_flange_edge_radius: None,
            top_flange_slope: None,
        })
    };
    let narrow = volume_of(&make(0.2));
    let wide = volume_of(&make(0.3));
    // The difference is exactly the extra top-flange material.
    let expected = (0.3 - 0.2) * 0.019 * DEPTH;
    assert!(
        ((wide - narrow) - expected).abs() < 1e-12,
        "expected {expected}, got {}",
        wide - narrow
    );
}

#[test]
fn an_absent_top_thickness_means_the_bottom_thickness() {
    let make = |top: Option<f64>| {
        Profile::Section(SectionProfile::AsymmetricI {
            depth: 0.4,
            web_thickness: 0.011,
            bottom_flange_width: 0.3,
            bottom_flange_thickness: 0.019,
            bottom_fillet_radius: None,
            bottom_flange_edge_radius: None,
            bottom_flange_slope: None,
            top_flange_width: 0.25,
            top_flange_thickness: top,
            top_fillet_radius: None,
            top_flange_edge_radius: None,
            top_flange_slope: None,
        })
    };
    let a = volume_of(&make(None));
    let b = volume_of(&make(Some(0.019)));
    assert!(
        (a - b).abs() < 1e-15,
        "absent must mean the bottom value: {a} vs {b}"
    );
}

#[test]
fn a_z_section_matches_its_closed_form_area() {
    let (d, bf, tw, tf, r) = (0.2, 0.075, 0.008, 0.012, 0.008);
    let profile = Profile::Section(SectionProfile::Z {
        depth: d,
        flange_width: bf,
        web_thickness: tw,
        flange_thickness: tf,
        fillet_radius: Some(r),
        edge_radius: None,
    });
    // Full-depth web plus one outstand per flange, plus two root fillets.
    let expected = (d * tw + 2.0 * bf * tf + 2.0 * fillet_area(r)) * DEPTH;
    let got = volume_of(&profile);
    assert!(
        (got - expected).abs() < 1e-12,
        "expected {expected}, got {got}"
    );
}

#[test]
fn a_c_section_encloses_its_wall_not_its_envelope() {
    // A lipped channel is THIN-WALLED: the boundary follows the wall all the
    // way round, so the area is the material, not the full envelope. Getting
    // this wrong would over-report the section by a large factor.
    let (d, w, t, g) = (0.2, 0.075, 0.002, 0.02);
    let profile = Profile::Section(SectionProfile::C {
        depth: d,
        width: w,
        wall_thickness: t,
        girth: g,
        internal_fillet_radius: None,
    });
    let got = volume_of(&profile);
    let envelope = d * w * DEPTH;
    assert!(
        got < envelope * 0.5,
        "a thin wall must be far less than the envelope {envelope}, got {got}"
    );
    // Wall length times thickness, as a sanity bound: web + 2 flanges + 2
    // lips, each of thickness t.
    let wall = (d + 2.0 * (w - t) + 2.0 * (g - t)) * t * DEPTH;
    assert!(
        (got - wall).abs() < wall * 0.05,
        "expected about {wall}, got {got}"
    );
}

#[test]
fn a_declared_zero_slope_is_not_a_taper() {
    // `Some(0.0)` states a parallel flange explicitly; only a NON-ZERO slope
    // is a taper.
    let profile = Profile::Section(SectionProfile::I {
        depth: 0.4,
        width: 0.3,
        web_thickness: 0.011,
        flange_thickness: 0.019,
        fillet_radius: None,
        flange_edge_radius: None,
        flange_slope: Some(0.0),
    });
    assert!(extrude_profile_exact(&profile, Vec3::Z, DEPTH, Tolerance::METRE).is_ok());
}

#[test]
fn impossible_dimensions_are_refused() {
    // Flanges thicker than the depth leave no web; a web wider than the
    // flange is not an I at all.
    let no_web = Profile::Section(SectionProfile::I {
        depth: 0.03,
        width: 0.3,
        web_thickness: 0.011,
        flange_thickness: 0.019,
        fillet_radius: None,
        flange_edge_radius: None,
        flange_slope: None,
    });
    let fat_web = Profile::Section(SectionProfile::I {
        depth: 0.4,
        width: 0.01,
        web_thickness: 0.011,
        flange_thickness: 0.019,
        fillet_radius: None,
        flange_edge_radius: None,
        flange_slope: None,
    });
    for profile in [no_web, fat_web] {
        assert!(extrude_profile_exact(&profile, Vec3::Z, DEPTH, Tolerance::METRE).is_err());
    }
}

#[test]
fn a_fillet_too_large_for_its_edge_is_refused() {
    // The web is 0.011 wide, so a 0.5 root radius cannot fit between the two
    // fillets that share the web face.
    let profile = i_section(Some(0.5));
    let error = extrude_profile_exact(&profile, Vec3::Z, DEPTH, Tolerance::METRE)
        .expect_err("the radius does not fit");
    assert!(
        format!("{error:?}").contains("too large for the edge"),
        "got {error:?}"
    );
}

#[test]
fn a_tapered_flange_keeps_the_declared_mean_thickness() {
    // The taper pivots about the mid-point of the inner face, so the MEAN
    // flange thickness is unchanged and the area matches the parallel case.
    // Any other pivot silently changes the declared thickness.
    let parallel = i_sloped(None, None);
    let tapered = i_sloped(None, Some(8.0_f64.to_radians()));

    let flat = contour_area(&parallel);
    let sloped = contour_area(&tapered);
    assert!(
        (flat - sloped).abs() < 1e-12,
        "taper must preserve the mean thickness: {flat} vs {sloped}"
    );
}

#[test]
fn a_tapered_flange_actually_slopes() {
    // Guard against the taper being accepted and then ignored: the two ends
    // of the bottom flange's inner face must sit at different heights.
    let contour = section_contour(section_of(&i_sloped(None, Some(8.0_f64.to_radians()))))
        .expect("tapered lowers");
    let parallel = section_contour(section_of(&i_sloped(None, None))).expect("parallel lowers");

    // Compare corner heights directly: a line segment's origin is a corner.
    let heights = |profile: &axiolid_profile::ContourProfile| -> Vec<f64> {
        profile
            .outer
            .segments
            .iter()
            .filter_map(|segment| match &segment.curve {
                Curve2::Line(line) => Some(line.origin.y),
                _ => None,
            })
            .collect()
    };
    assert_ne!(
        heights(&contour),
        heights(&parallel),
        "a declared taper must change the outline"
    );
}

#[test]
fn a_tapered_fillet_stays_tangent_to_the_inclined_face() {
    // The whole point of the gap: the root fillet must touch the INCLINED
    // flange face, not a horizontal one. Tangency means the arc centre sits
    // exactly its own radius from both adjacent faces.
    let radius = 0.021;
    let tapered = i_sloped(Some(radius), Some(8.0_f64.to_radians()));
    let contour = section_contour(section_of(&tapered)).expect("tapered lowers");

    let arcs: Vec<_> = contour
        .outer
        .segments
        .iter()
        .filter_map(|segment| match &segment.curve {
            Curve2::Circle(circle) => Some(*circle),
            _ => None,
        })
        .collect();
    assert_eq!(arcs.len(), 4, "four root fillets, got {}", arcs.len());

    // Each fillet arc must carry the stated radius and be tangent to the two
    // straight segments it joins.
    for arc in &arcs {
        assert!(
            (arc.radius - radius).abs() < 1e-12,
            "fillet radius {} should be {radius}",
            arc.radius
        );
    }

    // Tangency, checked against the segments the arc actually JOINS -- not
    // the nearest line anywhere in the outline. An unrelated face can sit
    // closer than the radius, so a global minimum would fail on correct
    // geometry (it did: 0.0136 against a face the fillet never touches).
    let segments = &contour.outer.segments;
    for (index, segment) in segments.iter().enumerate() {
        let Curve2::Circle(arc) = &segment.curve else {
            continue;
        };
        let before = &segments[(index + segments.len() - 1) % segments.len()];
        let after = &segments[(index + 1) % segments.len()];
        for neighbour in [before, after] {
            let Curve2::Line(line) = &neighbour.curve else {
                continue;
            };
            let direction = line.direction.normalize();
            let to_centre = arc.frame.origin - line.origin;
            let distance = (to_centre - direction * to_centre.dot(direction)).length();
            assert!(
                (distance - radius).abs() < 1e-9,
                "fillet must sit one radius from its ADJACENT face, got {distance}"
            );
        }
    }
}

#[test]
fn an_absurdly_steep_slope_is_refused() {
    let section = i_sloped(None, Some(1.5));
    assert!(
        section_contour(section_of(&section)).is_err(),
        "a near-right-angle taper leaves no flange and must be refused"
    );
}
