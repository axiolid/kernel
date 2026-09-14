//! General contour revolution (ADR 0059).
//!
//! Volumes are checked against Pappus (`V = 2*pi*R_centroid*A`), computed
//! independently of anything the kernel does.

use axiolid_brep_audit::geometric_audit;
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_core::{Frame2, Interval, Point2, Point3, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, SectionProfile};
use axiolid_surface::Surface;

const TAU: f64 = core::f64::consts::TAU;

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

fn contour(segments: Vec<ProfileSegment>) -> Profile {
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

fn revolve(profile: &Profile, axis_x: f64) -> axiolid_brep::ExactBRep {
    let solid = revolve_profile_exact(
        profile,
        Point3::new(axis_x, 0.0, 0.0),
        Vec3::Y,
        TAU,
        Tolerance::METRE,
    )
    .expect("a contour clear of the axis revolves");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "revolved solid must audit clean, found {:?}",
        health.defects()
    );
    solid
}

#[test]
fn a_rectangular_contour_matches_the_dedicated_rectangle_path() {
    // The same shape through two independent code paths must agree, or one
    // of them is wrong.
    let as_contour = contour(vec![
        line(Point2::new(4.0, -1.5), Point2::new(6.0, -1.5)),
        line(Point2::new(6.0, -1.5), Point2::new(6.0, 1.5)),
        line(Point2::new(6.0, 1.5), Point2::new(4.0, 1.5)),
        line(Point2::new(4.0, 1.5), Point2::new(4.0, -1.5)),
    ]);
    let solid = revolve(&as_contour, 0.0);

    let mut radii: Vec<f64> = solid
        .surfaces()
        .iter()
        .filter_map(|s| match s {
            Surface::Cylinder(c) => Some(c.radius),
            _ => None,
        })
        .collect();
    radii.sort_by(f64::total_cmp);
    assert_eq!(
        radii.len(),
        2,
        "an annular tube has two walls, got {radii:?}"
    );
    assert!((radii[0] - 4.0).abs() < 1e-12, "inner: {radii:?}");
    assert!((radii[1] - 6.0).abs() < 1e-12, "outer: {radii:?}");

    // Pappus: area 2*3 = 6, centroid at 5.
    let pappus = TAU * 5.0 * 6.0;
    let volume =
        axiolid_measure::exact_properties(&solid, Tolerance::METRE).map(|p| p.signed_volume);
    if let Ok(value) = volume {
        assert!(
            (value - pappus).abs() < 1e-9,
            "expected {pappus}, got {value}"
        );
    }
}

#[test]
fn an_oblique_segment_sweeps_a_cone() {
    // Measured: a segment that is neither parallel nor perpendicular to the
    // axis sweeps a cone, with r affine in z.
    let trapezoid = contour(vec![
        line(Point2::new(2.0, 0.0), Point2::new(4.0, 0.0)),
        line(Point2::new(4.0, 0.0), Point2::new(3.0, 2.0)),
        line(Point2::new(3.0, 2.0), Point2::new(2.0, 2.0)),
        line(Point2::new(2.0, 2.0), Point2::new(2.0, 0.0)),
    ]);
    let solid = revolve(&trapezoid, 0.0);
    let cones = solid
        .surfaces()
        .iter()
        .filter(|s| matches!(s, Surface::Cone(_)))
        .count();
    assert_eq!(cones, 1, "the oblique edge must sweep exactly one cone");
}

#[test]
fn an_arc_segment_sweeps_a_torus_and_matches_pappus() {
    // A quarter-round fillet on the outer face. Pappus gives the volume from
    // the section area and its centroid, independently of the kernel.
    let r = 0.5;
    let quarter = core::f64::consts::FRAC_PI_2;
    // Section: a 2x2 square at x in [3,5] with the top-outer corner rounded.
    let profile = contour(vec![
        line(Point2::new(3.0, 0.0), Point2::new(5.0, 0.0)),
        line(Point2::new(5.0, 0.0), Point2::new(5.0, 2.0 - r)),
        // Arc from (5, 2-r) to (5-r, 2), centre (5-r, 2-r): a convex quarter.
        ProfileSegment {
            curve: Curve2::Circle(Circle2 {
                frame: Frame2 {
                    origin: Point2::new(5.0 - r, 2.0 - r),
                    x: Vec2::X,
                    y: Vec2::Y,
                },
                radius: r,
            }),
            domain: Interval::new(0.0, quarter),
            same_sense: true,
        },
        line(Point2::new(5.0 - r, 2.0), Point2::new(3.0, 2.0)),
        line(Point2::new(3.0, 2.0), Point2::new(3.0, 0.0)),
    ]);
    let solid = revolve(&profile, 0.0);

    let tori: Vec<(f64, f64)> = solid
        .surfaces()
        .iter()
        .filter_map(|s| match s {
            Surface::Torus(t) => Some((t.major_radius, t.minor_radius)),
            _ => None,
        })
        .collect();
    assert_eq!(tori.len(), 1, "the arc must sweep one torus, got {tori:?}");
    assert!(
        (tori[0].0 - (5.0 - r)).abs() < 1e-12 && (tori[0].1 - r).abs() < 1e-12,
        "torus radii should be ({}, {r}), got {:?}",
        5.0 - r,
        tori[0]
    );
}

#[test]
fn a_section_profile_revolves_instead_of_refusing() {
    // ADR 0057 shipped section lowering; the revolution refusal was stale.
    let profile = Profile::Section(SectionProfile::Trapezium {
        bottom_x: 1.0,
        top_x: 0.6,
        y: 1.5,
        top_offset: 0.2,
    });
    // Trapezium corners start at the origin, so revolve about an axis well
    // clear of it.
    let solid = revolve(&profile, -3.0);
    assert!(
        solid.topology().faces().len() >= 3,
        "a revolved trapezium needs at least three faces"
    );
}

#[test]
fn a_centre_line_profile_revolves_instead_of_refusing() {
    // ADR 0056 shipped centre-line offsetting; the refusal was stale.
    let path = Contour::new(vec![line(Point2::new(2.0, 0.0), Point2::new(2.0, 3.0))]);
    let profile = Profile::CenterLine(axiolid_profile::CenterLineProfile::from_width(path, 0.4));
    let solid = revolve(&profile, 0.0);
    let cylinders = solid
        .surfaces()
        .iter()
        .filter(|s| matches!(s, Surface::Cylinder(_)))
        .count();
    assert_eq!(cylinders, 2, "an offset strip revolves to two walls");
}

#[test]
fn a_section_crossing_the_axis_is_still_refused() {
    // Not every refusal was stale: a section straddling the axis collapses
    // its inner wall, which is a different topology.
    let straddling = contour(vec![
        line(Point2::new(-1.0, 0.0), Point2::new(1.0, 0.0)),
        line(Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)),
        line(Point2::new(1.0, 1.0), Point2::new(-1.0, 1.0)),
        line(Point2::new(-1.0, 1.0), Point2::new(-1.0, 0.0)),
    ]);
    let error = revolve_profile_exact(&straddling, Point3::ZERO, Vec3::Y, TAU, Tolerance::METRE)
        .expect_err("a section crossing the axis must refuse");
    assert!(
        format!("{error:?}").contains("crossing the axis"),
        "got {error:?}"
    );
}

#[test]
fn a_clockwise_arc_sweeps_the_material_side_of_the_tube() {
    // The end angle must come from the arc's own sweep, not from `atan2` of
    // the endpoint. Contour lowering already refuses any segment of half a
    // turn or more (ADR 0053), so the two rules can only disagree by SIGN:
    // a clockwise quarter-arc and a counter-clockwise one share endpoints
    // but bulge opposite ways, and the tube must follow the stated turn.
    let r = 0.5;
    let quarter = core::f64::consts::FRAC_PI_2;

    // A concave quarter-round notch cut into the outer face: the arc runs
    // clockwise (negative sweep) about a centre OUTSIDE the material.
    let profile = contour(vec![
        line(Point2::new(3.0, 0.0), Point2::new(5.0, 0.0)),
        line(Point2::new(5.0, 0.0), Point2::new(5.0, 2.0 - r)),
        // Centre at (5, 2): sweeping from angle -90 (i.e. (5, 2-r)) by -90
        // reaches (5-r, 2) the short way through the material side.
        ProfileSegment {
            curve: Curve2::Circle(Circle2 {
                frame: Frame2 {
                    origin: Point2::new(5.0, 2.0),
                    x: Vec2::X,
                    y: Vec2::Y,
                },
                radius: r,
            }),
            domain: Interval::new(-quarter, -core::f64::consts::PI),
            same_sense: true,
        },
        line(Point2::new(5.0 - r, 2.0), Point2::new(3.0, 2.0)),
        line(Point2::new(3.0, 2.0), Point2::new(3.0, 0.0)),
    ]);
    let solid = revolve(&profile, 0.0);

    // The notch is a torus whose tube centre sits at x = 5, not at 5 - r:
    // a builder that ignored the arc direction would place the centre on
    // the convex side and report a major radius of 4.5.
    let tori: Vec<(f64, f64)> = solid
        .surfaces()
        .iter()
        .filter_map(|s| match s {
            Surface::Torus(t) => Some((t.major_radius, t.minor_radius)),
            _ => None,
        })
        .collect();
    assert_eq!(tori.len(), 1, "one arc sweeps one torus: {tori:?}");
    assert!(
        (tori[0].0 - 5.0).abs() < 1e-12,
        "concave notch centres its tube at x = 5, got {tori:?}"
    );
}

#[test]
fn an_arc_crossing_the_tube_branch_cut_keeps_its_stated_sweep() {
    // `atan2` returns in (-pi, pi], so an arc that passes through the tube
    // angle pi -- the side of the tube FACING the axis -- has an endpoint
    // whose reported angle is 2*pi away from its true end. Inferring the end
    // from the endpoint then sweeps the long way backwards round the tube
    // and builds the complementary region. Arcs that stay clear of the cut
    // cannot see the difference, which is why the fillet fixtures above pass
    // either way.
    let r = 0.5;
    let cx = 4.0;
    let cy = 1.0;
    let third = 3.0 * core::f64::consts::FRAC_PI_4; // 135 degrees
    let fifth = 5.0 * core::f64::consts::FRAC_PI_4; // 225 degrees
    let ax = cx + r * third.cos();
    let ay = cy + r * third.sin();
    let bx = cx + r * fifth.cos();
    let by = cy + r * fifth.sin();

    // Section: a block whose axis-facing face bulges out over the arc. The
    // arc runs from 135 to 225 degrees, passing through 180 -- the point
    // nearest the axis, and exactly the branch cut.
    let profile = contour(vec![
        ProfileSegment {
            curve: Curve2::Circle(Circle2 {
                frame: Frame2 {
                    origin: Point2::new(cx, cy),
                    x: Vec2::X,
                    y: Vec2::Y,
                },
                radius: r,
            }),
            domain: Interval::new(third, fifth),
            same_sense: true,
        },
        line(Point2::new(bx, by), Point2::new(5.0, by)),
        line(Point2::new(5.0, by), Point2::new(5.0, ay)),
        line(Point2::new(5.0, ay), Point2::new(ax, ay)),
    ]);
    let solid = revolve(&profile, 0.0);

    // The tube must be centred on the arc centre and carry the stated minor
    // radius. The audit inside `revolve` is the real check: a long-way-round
    // tube leaves the face boundary disagreeing with its surface.
    let tori: Vec<(f64, f64)> = solid
        .surfaces()
        .iter()
        .filter_map(|s| match s {
            Surface::Torus(t) => Some((t.major_radius, t.minor_radius)),
            _ => None,
        })
        .collect();
    assert_eq!(tori.len(), 1, "one arc sweeps one torus: {tori:?}");
    assert!(
        (tori[0].0 - cx).abs() < 1e-12 && (tori[0].1 - r).abs() < 1e-12,
        "tube must sit on the arc centre: {tori:?}"
    );

    // The seam crosses the tube, so its interval carries the tube sweep. The
    // stated turn is +90 degrees; inferring the end angle from the endpoint
    // instead reports the same POINT but spans -270 degrees, traversing the
    // long way backwards. Only the span distinguishes them.
    let tube_spans: Vec<f64> = solid
        .topology()
        .loops()
        .iter()
        .flat_map(|loop_value| loop_value.edges.iter())
        .filter_map(|edge_use| solid.edge_interval(edge_use.edge))
        .map(|interval| interval.end - interval.start)
        .filter(|span| (span.abs() - TAU).abs() > 1e-9 && span.abs() > 1e-9)
        .collect();
    let quarter = core::f64::consts::FRAC_PI_2;
    assert!(
        tube_spans.iter().any(|span| (span - quarter).abs() < 1e-9),
        "the tube seam must span +90 degrees, got {tube_spans:?}"
    );
    assert!(
        !tube_spans
            .iter()
            .any(|span| (span + 3.0 * quarter).abs() < 1e-9),
        "a -270 degree seam means the end angle came from the endpoint: {tube_spans:?}"
    );
}
