//! Parametric profiles with decimal sizes lower and extrude at
//! `Tolerance::ZERO` (#250).
//!
//! Rolled-section tables state sizes like IPE 300's `tw = 0.0071`, none of
//! them dyadic, so a contour's segments meet only to the rounding of their
//! own evaluation: a line stored as `origin + t direction` ends at a rounded
//! `origin + direction`, an arc at the `cos`/`sin` of its sweep. Lowering
//! used to demand bit-equal joints at `Tolerance::ZERO` and refused every
//! such profile ("contour segments leave a gap of 2.6e-18"). Each ring
//! vertex is ONE value shared by the edge entering and the edge leaving it,
//! so the ring is closed by construction; the joint check now allows the
//! evaluation's own rounding and still refuses a contour open by more.
//!
//! Areas are closed forms derived independently of the kernel: a root
//! fillet adds `(1 - pi/4) r^2` over the sharp corner it replaces, a toe or
//! edge radius removes as much.

use axiolid_brep_audit::geometric_audit;
use axiolid_construct::contour_lower::{contour_to_arc_ring, orient_arc_ring};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::section_lower::{circle_contour, rectangle_contour, section_contour};
use axiolid_contracts::GeomError;
use axiolid_core::{Interval, Point2, Tolerance, Vec2, Vec3};
use axiolid_curve::{Curve2, Line2};
use axiolid_overlay::ArcRing;
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile,
    SectionProfile,
};

const PI: f64 = core::f64::consts::PI;
const LENGTH: f64 = 2.7;

/// Area one fillet adds over a sharp concave corner, or one edge radius
/// removes from a sharp convex one.
fn k(radius: f64) -> f64 {
    radius * radius * (1.0 - PI / 4.0)
}

/// A profile under test, its closed-form area, and a name for messages.
struct Case {
    name: String,
    profile: Profile,
    area: f64,
}

fn case(name: impl Into<String>, profile: Profile, area: f64) -> Case {
    Case {
        name: name.into(),
        profile,
        area,
    }
}

/// I sections from the IPE and HEA tables, with and without root fillets
/// and flange toe radii.
fn i_cases() -> Vec<Case> {
    // (name, h, b, tw, tf, r)
    let tables = [
        ("IPE 300", 0.3, 0.15, 0.0071, 0.0107, 0.015),
        ("IPE 200", 0.2, 0.1, 0.0056, 0.0085, 0.012),
        ("HEA 200", 0.19, 0.2, 0.0065, 0.01, 0.018),
        ("HEB 340", 0.34, 0.3, 0.012, 0.0215, 0.027),
    ];
    let mut cases = Vec::new();
    for (name, h, b, tw, tf, r) in tables {
        for (fillet, toe) in [(None, None), (Some(r), None), (Some(r), Some(0.0035))] {
            let sharp = 2.0 * b * tf + (h - 2.0 * tf) * tw;
            let area = sharp + 4.0 * fillet.map_or(0.0, k) - 4.0 * toe.map_or(0.0, k);
            cases.push(case(
                format!("{name} fillet {fillet:?} toe {toe:?}"),
                Profile::Section(SectionProfile::I {
                    depth: h,
                    width: b,
                    web_thickness: tw,
                    flange_thickness: tf,
                    fillet_radius: fillet,
                    flange_edge_radius: toe,
                    flange_slope: None,
                }),
                area,
            ));
        }
    }
    cases
}

/// Every other family, with decimal sizes and with and without radii.
fn other_cases() -> Vec<Case> {
    let mut cases = Vec::new();
    for radii in [false, true] {
        let pick = |r: f64| radii.then_some(r);
        let r = |value: f64| if radii { k(value) } else { 0.0 };

        // Asymmetric I: a welded plate girder with unequal flanges.
        let (h, tw, bb, tfb, bt, tft) = (0.43, 0.0093, 0.27, 0.0171, 0.19, 0.0127);
        cases.push(case(
            format!("asymmetric I radii {radii}"),
            Profile::Section(SectionProfile::AsymmetricI {
                depth: h,
                web_thickness: tw,
                bottom_flange_width: bb,
                bottom_flange_thickness: tfb,
                bottom_fillet_radius: pick(0.021),
                bottom_flange_edge_radius: pick(0.0043),
                bottom_flange_slope: None,
                top_flange_width: bt,
                top_flange_thickness: Some(tft),
                top_fillet_radius: pick(0.013),
                top_flange_edge_radius: pick(0.0031),
                top_flange_slope: None,
            }),
            bb * tfb + bt * tft + (h - tfb - tft) * tw + 2.0 * (r(0.021) + r(0.013))
                - 2.0 * (r(0.0043) + r(0.0031)),
        ));

        // T 100 (EN 10055).
        let (h, b, tw, tf) = (0.1, 0.1, 0.011, 0.011);
        cases.push(case(
            format!("T 100 radii {radii}"),
            Profile::Section(SectionProfile::T {
                depth: h,
                flange_width: b,
                web_thickness: tw,
                flange_thickness: tf,
                fillet_radius: pick(0.011),
                flange_edge_radius: pick(0.0055),
                web_edge_radius: pick(0.003),
                web_slope: None,
                flange_slope: None,
            }),
            b * tf + (h - tf) * tw + 2.0 * r(0.011) - 2.0 * r(0.0055) - 2.0 * r(0.003),
        ));

        // UPE 200: parallel-flange channel.
        let (h, b, tw, tf) = (0.2, 0.08, 0.006, 0.011);
        cases.push(case(
            format!("UPE 200 radii {radii}"),
            Profile::Section(SectionProfile::U {
                depth: h,
                flange_width: b,
                web_thickness: tw,
                flange_thickness: tf,
                fillet_radius: pick(0.013),
                edge_radius: pick(0.0045),
                flange_slope: None,
            }),
            2.0 * b * tf + (h - 2.0 * tf) * tw + 2.0 * r(0.013) - 2.0 * r(0.0045),
        ));

        // L 120 x 80 x 10, unequal legs.
        let (d, w, t) = (0.12, 0.08, 0.01);
        cases.push(case(
            format!("L 120x80x10 radii {radii}"),
            Profile::Section(SectionProfile::L {
                depth: d,
                width: Some(w),
                thickness: t,
                fillet_radius: pick(0.011),
                edge_radius: pick(0.0055),
                leg_slope: None,
            }),
            t * (d + w - t) + r(0.011) - 2.0 * r(0.0055),
        ));

        // Z section.
        let (h, b, tw, tf) = (0.2, 0.073, 0.0071, 0.0109);
        cases.push(case(
            format!("Z 200 radii {radii}"),
            Profile::Section(SectionProfile::Z {
                depth: h,
                flange_width: b,
                web_thickness: tw,
                flange_thickness: tf,
                fillet_radius: pick(0.011),
                edge_radius: pick(0.0055),
            }),
            h * tw + 2.0 * b * tf + 2.0 * r(0.011) - 2.0 * r(0.0055),
        ));

        // Cold-formed lipped channel.
        let (h, w, t, g) = (0.2, 0.075, 0.0025, 0.021);
        cases.push(case(
            format!("C 200 radii {radii}"),
            Profile::Section(SectionProfile::C {
                depth: h,
                width: w,
                wall_thickness: t,
                girth: g,
                internal_fillet_radius: pick(0.0037),
            }),
            t * (h + 2.0 * (w - t) + 2.0 * (g - t)) + 4.0 * r(0.0037),
        ));

        // Rectangles: rounded, hollow, and both (RHS 200 x 100 x 6.3).
        let (x, y, t, ro, ri) = (0.2, 0.1, 0.0063, 0.0158, 0.0095);
        let outer = x * y - (4.0 - PI) * if radii { ro * ro } else { 0.0 };
        let inner = (x - 2.0 * t) * (y - 2.0 * t) - (4.0 - PI) * if radii { ri * ri } else { 0.0 };
        cases.push(case(
            format!("rectangle radii {radii}"),
            Profile::Rectangle(RectangleProfile {
                x,
                y,
                thickness: None,
                outer_radius: pick(ro),
                inner_radius: None,
            }),
            outer,
        ));
        cases.push(case(
            format!("RHS 200x100x6.3 radii {radii}"),
            Profile::Rectangle(RectangleProfile {
                x,
                y,
                thickness: Some(t),
                outer_radius: pick(ro),
                inner_radius: pick(ri),
            }),
            outer - inner,
        ));
    }

    let (bottom, top, y) = (0.31, 0.17, 0.23);
    cases.push(case(
        "trapezium",
        Profile::Section(SectionProfile::Trapezium {
            bottom_x: bottom,
            top_x: top,
            y,
            top_offset: 0.07,
        }),
        (bottom + top) / 2.0 * y,
    ));

    // CHS 168.3 x 7.1.
    let (radius, t) = (0.08415, 0.0071);
    cases.push(case(
        "CHS 168.3x7.1",
        Profile::Circle(CircleProfile {
            radius,
            thickness: Some(t),
        }),
        PI * (radius * radius - (radius - t) * (radius - t)),
    ));
    cases
}

fn all_cases() -> Vec<Case> {
    let mut cases = i_cases();
    cases.extend(other_cases());
    cases
}

fn contour_of(profile: &Profile) -> ContourProfile {
    match profile {
        Profile::Section(section) => section_contour(section),
        Profile::Rectangle(rectangle) => rectangle_contour(rectangle),
        Profile::Circle(circle) => circle_contour(circle),
        other => panic!("not a parametric profile: {other:?}"),
    }
    .expect("a parametric profile lowers to a contour")
}

/// Every boundary of a contour profile, outer first.
fn boundaries(contour: &ContourProfile) -> Vec<&Contour> {
    core::iter::once(&contour.outer)
        .chain(&contour.holes)
        .collect()
}

/// The ring is the contour, joint for joint: one vertex per straight
/// segment, and each line's stored origin -- the point the router computed
/// once and shared with the segment before it -- is its vertex bit for bit.
fn assert_ring_matches(name: &str, contour: &Contour, ring: &ArcRing) {
    let mut index = 0;
    for segment in &contour.segments {
        match &segment.curve {
            Curve2::Line(line) => {
                assert!(segment.same_sense && segment.domain.start == 0.0);
                assert_eq!(
                    ring.vertices[index].point, line.origin,
                    "{name}: ring vertex {index} is not the shared corner"
                );
                assert_eq!(ring.vertices[index].bulge, 0.0, "{name}");
                index += 1;
            }
            Curve2::Circle(_) => {
                // A section arc turns a quarter at most: one bulged vertex.
                assert_ne!(ring.vertices[index].bulge, 0.0, "{name}: arc lost");
                index += 1;
            }
            other => panic!("{name}: unexpected segment {other:?}"),
        }
    }
    assert_eq!(index, ring.vertices.len(), "{name}: vertex count");
}

#[test]
fn every_parametric_family_lowers_at_zero_tolerance() {
    for case in all_cases() {
        let contour = contour_of(&case.profile);
        for boundary in boundaries(&contour) {
            let ring = contour_to_arc_ring(boundary, Tolerance::ZERO)
                .unwrap_or_else(|e| panic!("{}: {e}", case.name));
            assert_ring_matches(&case.name, boundary, &ring);
            // Lowering at a positive tolerance builds the very same ring:
            // the joint check decides whether a ring is built, never what.
            let at_metre = contour_to_arc_ring(boundary, Tolerance::METRE).expect("lowers");
            assert_eq!(ring, at_metre, "{}", case.name);
        }
    }
}

#[test]
fn every_parametric_family_extrudes_to_its_closed_form_at_zero_tolerance() {
    // An annular circle profile lowers (above) but its direct extrusion is
    // a separately refused kind at any tolerance; it is not a joint case.
    for case in all_cases()
        .into_iter()
        .filter(|case| !matches!(case.profile, Profile::Circle(_)))
    {
        let solid = extrude_profile_exact(&case.profile, Vec3::Z, LENGTH, Tolerance::ZERO)
            .unwrap_or_else(|e| panic!("{}: {e}", case.name));
        let health = geometric_audit(&solid, Tolerance::METRE);
        assert!(
            health.is_consistent(),
            "{}: {:?}",
            case.name,
            health.defects()
        );
        let topology = axiolid_topology::audit_brep(solid.topology());
        assert!(topology.is_closed_manifold(), "{}: {topology:?}", case.name);
        let volume = axiolid_measure::exact_properties(&solid, Tolerance::METRE)
            .expect("measurable")
            .signed_volume;
        let expected = case.area * LENGTH;
        assert!(
            (volume - expected).abs() <= 1e-12 * expected,
            "{}: volume {volume}, closed form {expected}",
            case.name
        );
    }
}

#[test]
fn tapered_flanges_lower_at_zero_tolerance() {
    // Tapers incline the faces the fillets meet, so the tangent points are
    // not on any axis and every joint rounds.
    let slope = 0.08_f64.atan();
    let profiles = [
        SectionProfile::I {
            depth: 0.3,
            width: 0.125,
            web_thickness: 0.0108,
            flange_thickness: 0.0162,
            fillet_radius: Some(0.0108),
            flange_edge_radius: Some(0.0065),
            flange_slope: Some(0.14_f64.atan()),
        },
        SectionProfile::U {
            depth: 0.2,
            flange_width: 0.075,
            web_thickness: 0.0085,
            flange_thickness: 0.0115,
            fillet_radius: Some(0.0115),
            edge_radius: Some(0.006),
            flange_slope: Some(slope),
        },
        SectionProfile::L {
            depth: 0.12,
            width: Some(0.08),
            thickness: 0.01,
            fillet_radius: Some(0.011),
            edge_radius: Some(0.0055),
            leg_slope: Some(0.03),
        },
        SectionProfile::T {
            depth: 0.1,
            flange_width: 0.1,
            web_thickness: 0.011,
            flange_thickness: 0.011,
            fillet_radius: Some(0.011),
            flange_edge_radius: Some(0.0055),
            web_edge_radius: Some(0.003),
            web_slope: Some(0.02),
            flange_slope: Some(0.02),
        },
    ];
    for section in profiles {
        let contour = section_contour(&section).expect("lowers");
        let ring = contour_to_arc_ring(&contour.outer, Tolerance::ZERO)
            .unwrap_or_else(|e| panic!("{section:?}: {e}"));
        assert_ring_matches("tapered", &contour.outer, &ring);
        let solid = extrude_profile_exact(
            &Profile::Section(section.clone()),
            Vec3::Z,
            1.5,
            Tolerance::ZERO,
        )
        .unwrap_or_else(|e| panic!("{section:?}: {e}"));
        // Same solid as at a positive tolerance, to the bit.
        let at_metre =
            extrude_profile_exact(&Profile::Section(section), Vec3::Z, 1.5, Tolerance::METRE)
                .expect("extrudes");
        let volume = |solid| {
            axiolid_measure::exact_properties(solid, Tolerance::METRE)
                .expect("measurable")
                .signed_volume
        };
        assert_eq!(volume(&solid), volume(&at_metre));
        assert!(volume(&solid) > 0.0);
    }
}

fn line(from: Point2, to: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}

/// A unit square whose second corner is reached `gap` short.
fn square_with_gap(gap: f64) -> Contour {
    let a = Point2::new(0.0, 0.0);
    let b = Point2::new(1.0, 0.0);
    let c = Point2::new(1.0, 1.0);
    let d = Point2::new(0.0, 1.0);
    Contour::new(vec![
        line(a, b - Vec2::new(gap, 0.0)),
        line(b, c),
        line(c, d),
        line(d, a),
    ])
}

/// A unit square whose last edge stops `gap` short of the start.
fn square_not_closing(gap: f64) -> Contour {
    let a = Point2::new(0.0, 0.0);
    let b = Point2::new(1.0, 0.0);
    let c = Point2::new(1.0, 1.0);
    let d = Point2::new(0.0, 1.0);
    Contour::new(vec![
        line(a, b),
        line(b, c),
        line(c, d),
        line(d, a + Vec2::new(0.0, gap)),
    ])
}

/// The rounding a joint may show is that of BOTH evaluations meeting there,
/// each scaled by what it is computed from -- not by where the joint is.
///
/// A line run in from far away ends to an ulp of its far origin, and a
/// big-radius arc starts to an ulp of its far centre, even where the joint
/// itself sits near the origin. Each case is checked to really round, so a
/// bound that ignored either side would refuse it at ZERO.
#[test]
fn a_joint_rounds_with_the_magnitudes_it_is_evaluated_from() {
    // A sliver triangle with one vertex 10 km out.
    let far = Point2::new(1e4, 0.3);
    let near = Point2::new(0.1, 0.3);
    let up = Point2::new(0.1, 0.7);
    assert_ne!(far + (near - far), near, "the long edge must round");
    let sliver = Contour::new(vec![line(far, near), line(near, up), line(up, far)]);
    contour_to_arc_ring(&sliver, Tolerance::ZERO).expect("closes to its rounding");

    // A shallow lens near the origin: an arc about a centre 10 km out, then
    // its chord back. The arc is the FIRST segment, so its start meets the
    // closing line's end.
    let start = Point2::new(0.13, 0.27);
    let centre = Point2::new(7000.3, -7100.7);
    let radius = (start - centre).length();
    let x = (start - centre).normalize();
    let y = Vec2::new(-x.y, x.x);
    let sweep: f64 = 1e-4;
    let circle = axiolid_curve::Circle2 {
        frame: axiolid_core::Frame2 {
            origin: centre,
            x,
            y,
        },
        radius,
    };
    assert_ne!(centre + x * radius, start, "the arc start must round");
    let end = centre + x * (radius * sweep.cos()) + y * (radius * sweep.sin());
    let lens = Contour::new(vec![
        ProfileSegment {
            curve: Curve2::Circle(circle),
            domain: Interval::new(0.0, sweep),
            same_sense: true,
        },
        line(end, start),
    ]);
    let ring = contour_to_arc_ring(&lens, Tolerance::ZERO).expect("closes to its rounding");
    assert_eq!(ring.vertices.len(), 2);

    // A 10 mm circle 10 km out as two arcs in different frames: the joint
    // rounds to an ulp of the CENTRE, a billion times the radius.
    let centre = Point2::new(7000.3, -7100.7);
    let radius = 0.01;
    let x = Vec2::new(1.0, 0.3).normalize();
    let y = Vec2::new(-x.y, x.x);
    let sweep: f64 = 1.1;
    let x2 = (x * sweep.cos() + y * sweep.sin()).normalize();
    let arc = |x: Vec2, to: f64| ProfileSegment {
        curve: Curve2::Circle(axiolid_curve::Circle2 {
            frame: axiolid_core::Frame2 {
                origin: centre,
                x,
                y: Vec2::new(-x.y, x.x),
            },
            radius,
        }),
        domain: Interval::new(0.0, to),
        same_sense: true,
    };
    let first_end = centre + x * (radius * sweep.cos()) + y * (radius * sweep.sin());
    assert_ne!(first_end, centre + x2 * radius, "the joint must round");
    let far_circle = Contour::new(vec![arc(x, sweep), arc(x2, core::f64::consts::TAU - sweep)]);
    contour_to_arc_ring(&far_circle, Tolerance::ZERO).expect("closes to its rounding");
}

#[test]
fn a_contour_open_by_more_than_its_rounding_is_still_refused_at_zero() {
    // 1e-13 is ~450 ulps of these coordinates: drawn, not rounded.
    for gap in [1e-13, 1e-9] {
        for contour in [square_with_gap(gap), square_not_closing(gap)] {
            let error = contour_to_arc_ring(&contour, Tolerance::ZERO)
                .expect_err("an open contour is refused at ZERO");
            assert!(
                matches!(&error, GeomError::InvalidInput(m)
                    if m.contains("gap") || m.contains("does not close")),
                "{error:?}"
            );
            // A tolerance that covers the gap welds it, as before.
            assert!(contour_to_arc_ring(&contour, Tolerance::METRE).is_ok());
        }
    }
    // An exactly shared square was always fine.
    let ring = contour_to_arc_ring(&square_with_gap(0.0), Tolerance::ZERO).expect("closed");
    assert_eq!(ring.vertices.len(), 4);
    // And orientation still reads the ring.
    assert!(orient_arc_ring(&ring, true).is_ok());
}
