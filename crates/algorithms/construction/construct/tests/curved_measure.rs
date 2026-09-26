//! Exact mass properties over curved faces (#125, ledger row C17).
//!
//! Every expectation is a closed form derived from the inputs -- `pi r^2 h`,
//! Pappus, a disc's parallel-axis moments -- never read back from the
//! kernel. Each solid mixes planar faces (fanned) with curved or
//! curve-bounded faces (integrated round their pcurves), so a mismatch
//! between the two paths' fields shows up as a wrong total.
//!
//! # Why this test lives in `axiolid-construct`
//!
//! Curved solids come from the constructors, and `axiolid-construct`
//! already depends on `axiolid-measure`; see `exact_measure.rs`.

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, clip_arc_prism_exact, ArcPrism};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_core::{
    BooleanOperator, Frame2, Interval, Plane3, Point2, Point3, Tolerance, Vec2, Vec3,
};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::{exact_properties, MassProperties};
use axiolid_overlay::ArcRing;
use axiolid_primitive::HalfSpace;
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, EllipseProfile, Profile, ProfileSegment,
    RectangleProfile,
};

const PI: f64 = std::f64::consts::PI;
const TAU: f64 = std::f64::consts::TAU;

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn measure(solid: &ExactBRep) -> MassProperties {
    let health = geometric_audit(solid, tol());
    assert!(
        health.is_consistent(),
        "fixture must audit clean: {:?}",
        health.defects()
    );
    exact_properties(solid, tol()).expect("a closed curved solid is measurable")
}

fn close(what: &str, got: f64, expected: f64) {
    let scale = expected.abs().max(1.0);
    assert!(
        (got - expected).abs() <= 1e-11 * scale,
        "{what}: expected {expected}, got {got} (off by {:e})",
        (got - expected).abs()
    );
}

fn extrude(profile: &Profile, depth: f64) -> ExactBRep {
    extrude_profile_exact(profile, Vec3::Z, depth, tol()).expect("extrudes exactly")
}

/// Revolve a profile a full turn about the line through `(axis_x, 0, 0)`
/// along `y`.
fn revolve(profile: &Profile, axis_x: f64) -> ExactBRep {
    revolve_profile_exact(profile, Point3::new(axis_x, 0.0, 0.0), Vec3::Y, TAU, tol())
        .expect("a profile clear of the axis revolves")
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

fn contour(segments: Vec<ProfileSegment>) -> Profile {
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

#[test]
fn a_cylinder_matches_its_closed_forms() {
    let (r, h) = (1.5, 2.0);
    let solid = extrude(
        &Profile::Circle(CircleProfile {
            radius: r,
            thickness: None,
        }),
        h,
    );
    let props = measure(&solid);

    let volume = PI * r * r * h;
    close("volume", props.signed_volume, volume);
    close("area", props.area, TAU * r * r + TAU * r * h);
    close("centroid x", props.centroid.x, 0.0);
    close("centroid y", props.centroid.y, 0.0);
    close("centroid z", props.centroid.z, h / 2.0);
    // int x^2 dV over a disc of radius r is pi r^4 / 4 per unit height.
    close(
        "second moment x",
        props.second_moment_diagonal.x,
        PI * r.powi(4) / 4.0 * h,
    );
    close(
        "second moment y",
        props.second_moment_diagonal.y,
        PI * r.powi(4) / 4.0 * h,
    );
    close(
        "second moment z",
        props.second_moment_diagonal.z,
        volume * h * h / 3.0,
    );
}

#[test]
fn a_plate_with_an_off_centre_round_hole_subtracts_the_bore() {
    // The bore is a cylinder wall used against its outward sense, and each
    // cap is planar with a circular hole: the curved wall's orientation and
    // the planar hole's winding must both subtract. Off-centre, so the hole
    // also moves the centroid.
    let (half, r, h) = (2.0, 0.75, 3.0);
    let (hx, hy) = (0.5, -0.25);
    let frame = Frame2 {
        origin: Point2::new(hx, hy),
        x: Vec2::X,
        y: Vec2::Y,
    };
    let quarter = std::f64::consts::FRAC_PI_2;
    let hole = Contour::new(
        (0..4)
            .map(|index| ProfileSegment {
                curve: Curve2::Circle(Circle2 { frame, radius: r }),
                domain: Interval::new(quarter * index as f64, quarter * (index + 1) as f64),
                same_sense: true,
            })
            .collect(),
    );
    let corners = [
        Point2::new(-half, -half),
        Point2::new(half, -half),
        Point2::new(half, half),
        Point2::new(-half, half),
    ];
    let outer = Contour::new(
        (0..4)
            .map(|i| line(corners[i], corners[(i + 1) % 4]))
            .collect(),
    );
    let solid = extrude(
        &Profile::Contour(ContourProfile {
            outer,
            holes: vec![hole],
        }),
        h,
    );
    let props = measure(&solid);

    let square = 4.0 * half * half;
    let disc = PI * r * r;
    close("volume", props.signed_volume, (square - disc) * h);
    close(
        "area",
        props.area,
        2.0 * (square - disc) + 8.0 * half * h + TAU * r * h,
    );
    close("centroid x", props.centroid.x, -disc * hx / (square - disc));
    close("centroid y", props.centroid.y, -disc * hy / (square - disc));
}

#[test]
fn an_elliptical_cylinder_matches_its_closed_forms() {
    let (a, b, h) = (3.0, 1.0, 2.0);
    let solid = extrude(
        &Profile::Ellipse(EllipseProfile {
            semi_axis_x: a,
            semi_axis_y: b,
        }),
        h,
    );
    let props = measure(&solid);
    close("volume", props.signed_volume, PI * a * b * h);
    // int x^2 over an ellipse is pi a^3 b / 4.
    close(
        "second moment x",
        props.second_moment_diagonal.x,
        PI * a.powi(3) * b / 4.0 * h,
    );
    close(
        "second moment y",
        props.second_moment_diagonal.y,
        PI * a * b.powi(3) / 4.0 * h,
    );
}

#[test]
fn an_off_axis_raised_disc_obeys_the_parallel_axis_theorem() {
    // Off the origin and off `z = 0`: every cone field sees non-zero
    // coordinates on every face, so no term vanishes by accident.
    let (cx, cy, r, bottom, top) = (3.0, -2.0, 1.25, 0.5, 2.5);
    let disc = |radius: f64| ArcPrism {
        section: ArcRing::circle(Point2::new(cx, cy), radius),
        bottom,
        top,
    };
    let solid = boolean_arc_prisms_exact(
        &disc(r),
        &disc(2.0 * r),
        BooleanOperator::Intersection,
        tol(),
    )
    .expect("a disc inside a disc intersects to itself");
    let props = measure(&solid);

    let h = top - bottom;
    let volume = PI * r * r * h;
    close("volume", props.signed_volume, volume);
    close("centroid x", props.centroid.x, cx);
    close("centroid y", props.centroid.y, cy);
    close("centroid z", props.centroid.z, (bottom + top) / 2.0);
    close(
        "second moment x",
        props.second_moment_diagonal.x,
        volume * (cx * cx + r * r / 4.0),
    );
    close(
        "second moment z",
        props.second_moment_diagonal.z,
        volume * (top.powi(3) - bottom.powi(3)) / (3.0 * h),
    );
}

#[test]
fn a_revolved_rectangle_matches_pappus() {
    // Two cylinder walls and two annular caps with circular edges.
    let section = contour(vec![
        line(Point2::new(4.0, -1.5), Point2::new(6.0, -1.5)),
        line(Point2::new(6.0, -1.5), Point2::new(6.0, 1.5)),
        line(Point2::new(6.0, 1.5), Point2::new(4.0, 1.5)),
        line(Point2::new(4.0, 1.5), Point2::new(4.0, -1.5)),
    ]);
    let props = measure(&revolve(&section, 0.0));
    close("volume", props.signed_volume, TAU * 5.0 * 6.0);
    close(
        "area",
        props.area,
        2.0 * PI * (36.0 - 16.0) + TAU * (4.0 + 6.0) * 3.0,
    );
    close("centroid y", props.centroid.y, 0.0);
}

#[test]
fn a_revolved_trapezoid_sweeps_a_cone_and_matches_pappus() {
    let corners = [
        Point2::new(2.0, 0.0),
        Point2::new(4.0, 0.0),
        Point2::new(3.0, 2.0),
        Point2::new(2.0, 2.0),
    ];
    let section = contour(
        (0..4)
            .map(|i| line(corners[i], corners[(i + 1) % 4]))
            .collect(),
    );
    let props = measure(&revolve(&section, 0.0));

    // Shoelace area and first moment about the axis (x = 0).
    let (mut area, mut moment) = (0.0, 0.0);
    for i in 0..4 {
        let (p, q) = (corners[i], corners[(i + 1) % 4]);
        let cross = p.x * q.y - q.x * p.y;
        area += cross / 2.0;
        moment += (p.x + q.x) * cross / 6.0;
    }
    assert!(area > 0.0);
    close("volume", props.signed_volume, TAU * moment);
    // Along the axis the solid's centroid weights the section by radius:
    // `int x y dA / int x dA`, not the section's own centroid.
    close(
        "centroid y",
        props.centroid.y,
        product_xy(&corners) / moment,
    );
}

/// `int x y dA` over a polygon, by the divergence theorem on its edges.
fn product_xy(corners: &[Point2]) -> f64 {
    let mut total = 0.0;
    for i in 0..corners.len() {
        let (p, q) = (corners[i], corners[(i + 1) % corners.len()]);
        let cross = p.x * q.y - q.x * p.y;
        total += cross * (p.x * (2.0 * p.y + q.y) + q.x * (p.y + 2.0 * q.y)) / 24.0;
    }
    total
}

#[test]
fn a_revolved_fillet_sweeps_a_torus_and_matches_pappus() {
    // A 2 x 2 square at x in [3, 5] with its top-outer corner rounded by r.
    let r = 0.5;
    let quarter = std::f64::consts::FRAC_PI_2;
    let section = contour(vec![
        line(Point2::new(3.0, 0.0), Point2::new(5.0, 0.0)),
        line(Point2::new(5.0, 0.0), Point2::new(5.0, 2.0 - r)),
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
    let props = measure(&revolve(&section, 0.0));

    // int x dA of the section: the square, less the r x r corner square,
    // plus back the quarter disc (centroid 4r / 3pi beyond its centre).
    let corner_square = r * r * (5.0 - r / 2.0);
    let quarter_disc = PI * r * r / 4.0 * (5.0 - r + 4.0 * r / (3.0 * PI));
    let moment = 4.0 * 4.0 - corner_square + quarter_disc;
    close("volume", props.signed_volume, TAU * moment);
}

#[test]
fn a_revolved_rounded_rectangle_sweeps_four_tori_and_matches_pappus() {
    let (x, y, r, centre) = (2.0, 3.0, 0.5, 5.0);
    let profile = Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: Some(r),
        inner_radius: None,
    });
    let solid = revolve_profile_exact(
        &profile,
        Point3::new(-centre, 0.0, 0.0),
        Vec3::Y,
        TAU,
        tol(),
    )
    .expect("a rounded rectangle clear of the axis revolves");
    let props = measure(&solid);
    // Symmetric section: its centroid is its centre, `centre` from the axis.
    let area = x * y - (4.0 - PI) * r * r;
    close("volume", props.signed_volume, TAU * centre * area);
}

#[test]
fn a_sloped_cut_through_a_column_matches_the_closed_form() {
    // The cut rim is an ellipse edge with a sinusoid pcurve on the wall;
    // the cap is a plane bounded by that ellipse. Volume is the section
    // area times the plane's height over the section centroid.
    let (cx, cy, r) = (2.0, -1.0, 1.5);
    let (h, gx, gy) = (8.0, 0.3, -0.4);
    let column = ArcPrism {
        section: ArcRing::circle(Point2::new(cx, cy), r),
        bottom: 1.0,
        top: 20.0,
    };
    let below = HalfSpace {
        boundary: Plane3 {
            origin: Point3::new(0.0, 0.0, h),
            normal: Vec3::new(-gx, -gy, 1.0),
        },
        agreement: false,
    };
    let solid = clip_arc_prism_exact(&column, &below, tol()).expect("representable");
    let props = measure(&solid);
    let mean = h + gx * cx + gy * cy;
    close("volume", props.signed_volume, PI * r * r * (mean - 1.0));
}

#[test]
fn exact_and_tessellated_measures_agree_on_a_curved_solid() {
    // Two implementations sharing no code: the exact integral must sit
    // inside the mesh's chord error, and much closer than the cruder mesh.
    use axiolid_construct::extrude::extrude_profile;
    use axiolid_construct::profile::profile_rings;
    use axiolid_measure::{Measure, MeshMeasure};

    let profile = Profile::Circle(CircleProfile {
        radius: 1.0,
        thickness: None,
    });
    let exact = measure(&extrude(&profile, 2.0)).signed_volume;
    let mut previous = f64::INFINITY;
    for chord in [1e-2, 1e-4, 1e-6] {
        let rings = profile_rings(&profile, chord, tol()).expect("rings");
        let mesh = extrude_profile(&rings, Vec3::Z, 2.0, tol()).expect("mesh");
        let meshed = MeshMeasure.measure(&mesh, tol()).expect("mesh measure");
        let gap = (exact - meshed.signed_volume).abs();
        assert!(
            gap < previous,
            "finer mesh must close in on the exact value"
        );
        previous = gap;
    }
    assert!(previous < 1e-5, "finest mesh still {previous} away");
}
