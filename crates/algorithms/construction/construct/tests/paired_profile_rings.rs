//! Profiles flattened together (#361): equal structure, vertices on the
//! exact curves, the interpolated section within the reported deviation,
//! orientation, a derived profile's stretch, the half step of circle-only
//! runs, and the families that do not pair.

use axiolid_construct::profile::paired_profile_rings;
use axiolid_core::{Frame2, Interval, Point2, Scalar, Tolerance, Transform2, Vec2};
use axiolid_curve::{BSplineCurve2, Circle2, Curve2, KnotSpec, Line2, Polyline2};
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, EllipseProfile, Profile, ProfileSegment,
    RectangleProfile,
};
use axiolid_reference::signed_area2;

const CHORD: Scalar = 0.005;
const TAU: Scalar = std::f64::consts::TAU;

fn tolerance() -> Tolerance {
    Tolerance::MILLIMETRE
}

fn disk(radius: Scalar) -> Profile {
    Profile::Circle(CircleProfile {
        radius,
        thickness: None,
    })
}

fn annulus(radius: Scalar, thickness: Scalar) -> Profile {
    Profile::Circle(CircleProfile {
        radius,
        thickness: Some(thickness),
    })
}

fn ellipse(a: Scalar, b: Scalar) -> Profile {
    Profile::Ellipse(EllipseProfile {
        semi_axis_x: a,
        semi_axis_y: b,
    })
}

fn rounded(radius: Scalar) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x: 2.0,
        y: 1.4,
        thickness: None,
        outer_radius: Some(radius),
        inner_radius: None,
    })
}

fn sharp(x: Scalar, y: Scalar) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

fn pair(profiles: &[&Profile]) -> axiolid_construct::profile::PairedRings {
    paired_profile_rings(profiles, CHORD, tolerance())
        .unwrap()
        .expect("the profiles pair")
}

fn distance_to_segment(p: Point2, a: Point2, b: Point2) -> Scalar {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0);
    (p - (a + ab * t)).length()
}

fn distance_to_ring(p: Point2, ring: &[Point2]) -> Scalar {
    (0..ring.len())
        .map(|k| distance_to_segment(p, ring[k], ring[(k + 1) % ring.len()]))
        .fold(Scalar::INFINITY, Scalar::min)
}

/// The polygon the loft makes between two matched rings at fraction `u`.
fn blended(a: &[Point2], b: &[Point2], u: Scalar) -> Vec<Point2> {
    a.iter().zip(b).map(|(p, q)| p.lerp(*q, u)).collect()
}

/// How far the exact section `exact(u, t)` (t in [0, 2 pi]) strays from the
/// blended rings, over a dense grid of `u` and `t`.
fn stray(
    paired: &axiolid_construct::profile::PairedRings,
    exact: impl Fn(Scalar, Scalar) -> Point2,
) -> Scalar {
    let mut worst: Scalar = 0.0;
    for i in 0..=20 {
        let u = i as Scalar / 20.0;
        let ring = blended(&paired.rings[0].outer, &paired.rings[1].outer, u);
        for k in 0..4000 {
            let t = TAU * k as Scalar / 4000.0;
            worst = worst.max(distance_to_ring(exact(u, t), &ring));
        }
    }
    worst
}

#[test]
fn equal_structure_and_vertices_on_the_exact_curves() {
    let (small, large) = (disk(0.5), disk(0.8));
    let paired = pair(&[&small, &large]);
    let [a, b] = [&paired.rings[0], &paired.rings[1]];
    assert_eq!(a.outer.len(), b.outer.len());
    assert!(a.holes.is_empty() && b.holes.is_empty());
    assert!(a.outer.iter().all(|p| (p.length() - 0.5).abs() < 1e-12));
    assert!(b.outer.iter().all(|p| (p.length() - 0.8).abs() < 1e-12));
    assert!(signed_area2(&a.outer) > 0.0);
    assert!(paired.deviation <= CHORD * 1.0001, "{}", paired.deviation);

    let ring = annulus(1.0, 0.2);
    let thick = annulus(1.5, 0.3);
    let paired = pair(&[&ring, &thick]);
    let [a, b] = [&paired.rings[0], &paired.rings[1]];
    assert_eq!(a.outer.len(), b.outer.len());
    assert_eq!(a.holes.len(), 1);
    assert_eq!(a.holes[0].len(), b.holes[0].len());
    assert!(signed_area2(&a.holes[0]) < 0.0);
    assert!(a.holes[0].iter().all(|p| (p.length() - 0.8).abs() < 1e-12));
    assert!(b.holes[0].iter().all(|p| (p.length() - 1.2).abs() < 1e-12));

    let flat = ellipse(1.0, 0.5);
    let round = disk(0.8);
    let paired = pair(&[&flat, &round]);
    assert!(paired.rings[0]
        .outer
        .iter()
        .all(|p| ((p.x / 1.0).powi(2) + (p.y / 0.5).powi(2) - 1.0).abs() < 1e-12));
}

#[test]
fn the_interpolated_section_stays_within_the_deviation() {
    let (small, large) = (disk(0.5), disk(0.8));
    let paired = pair(&[&small, &large]);
    let far = stray(&paired, |u, t| {
        let r = 0.5 + 0.3 * u;
        Point2::new(r * t.cos(), r * t.sin())
    });
    assert!(far <= paired.deviation, "{far} > {}", paired.deviation);
    // The bound is not slack by orders of magnitude.
    assert!(
        far > 0.2 * paired.deviation,
        "{far} vs {}",
        paired.deviation
    );

    let flat = ellipse(1.0, 0.5);
    let round = disk(0.8);
    let paired = pair(&[&flat, &round]);
    let far = stray(&paired, |u, t| {
        Point2::new((1.0 - 0.2 * u) * t.cos(), (0.5 + 0.3 * u) * t.sin())
    });
    assert!(far <= paired.deviation, "{far} > {}", paired.deviation);
}

fn mirrored(basis: Profile) -> Profile {
    Profile::Derived {
        basis: Box::new(basis),
        transform: Transform2::from_scale(Vec2::new(-1.0, 1.0)),
    }
}

#[test]
fn a_mirrored_profile_does_not_pair_with_an_unmirrored_one() {
    // The ring of the mirrored profile is reversed and its start point stays
    // mirrored, so every vertex would match a point of the opposite side.
    for (a, b) in [
        (mirrored(rounded(0.1)), rounded(0.3)),
        (mirrored(disk(0.5)), disk(0.8)),
        (mirrored(sharp(1.0, 1.0)), sharp(2.0, 2.0)),
    ] {
        assert!(paired_profile_rings(&[&a, &b], CHORD, tolerance())
            .unwrap()
            .is_none());
    }
    // Two mirrored profiles agree and pair, each ring counter-clockwise.
    let (a, b) = (mirrored(rounded(0.1)), mirrored(rounded(0.3)));
    let paired = pair(&[&a, &b]);
    for rings in &paired.rings {
        assert!(signed_area2(&rings.outer) > 0.0);
    }
    assert_eq!(paired.rings[0].outer.len(), paired.rings[1].outer.len());
    assert!(paired.rings[0]
        .outer
        .iter()
        .all(|p| p.x.abs() <= 1.0 + 1e-12 && p.y.abs() <= 0.7 + 1e-12));
}

#[test]
fn a_non_uniformly_scaled_derived_profile_pairs_and_charges_its_stretch() {
    let squashed = Profile::Derived {
        basis: Box::new(disk(0.5)),
        transform: Transform2::from_scale(Vec2::new(2.0, 1.0)),
    };
    let round = disk(0.8);
    let paired = pair(&[&squashed, &round]);
    assert!(paired.rings[0]
        .outer
        .iter()
        .all(|p| ((p.x / 1.0).powi(2) + (p.y / 0.5).powi(2) - 1.0).abs() < 1e-12));
    let far = stray(&paired, |u, t| {
        Point2::new(
            ((1.0 - u) * 1.0 + u * 0.8) * t.cos(),
            ((1.0 - u) * 0.5 + u * 0.8) * t.sin(),
        )
    });
    assert!(far <= paired.deviation, "{far} > {}", paired.deviation);
    assert!(paired.deviation <= CHORD * 1.0001);
}

/// A half annulus of outer radius 1 and the given width: the outer arc over
/// `[0, pi]`, a line along the axis, the inner arc back in four pieces, a
/// line out to the outer radius. Thin enough that a coarse sampling of the
/// ring gets the sign of its area wrong.
fn band(width: Scalar, reversed: bool) -> Profile {
    let inner = 1.0 - width;
    let pi = std::f64::consts::PI;
    let line = |from: Point2, to: Point2| ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::new(0.0, 1.0),
        same_sense: true,
    };
    let arc = |radius: Scalar, start: Scalar, end: Scalar, same_sense: bool| ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(0.0, 0.0),
                x: Vec2::new(1.0, 0.0),
                y: Vec2::new(0.0, 1.0),
            },
            radius,
        }),
        domain: Interval::new(start, end),
        same_sense,
    };
    let mut segments = vec![
        arc(1.0, 0.0, pi, true),
        line(Point2::new(-1.0, 0.0), Point2::new(-inner, 0.0)),
    ];
    for k in 0..4 {
        let (hi, lo) = (
            pi - k as Scalar * pi / 4.0,
            pi - (k + 1) as Scalar * pi / 4.0,
        );
        segments.push(arc(inner, lo, hi, false));
    }
    segments.push(line(Point2::new(inner, 0.0), Point2::new(1.0, 0.0)));
    if reversed {
        // The same band written clockwise.
        segments.reverse();
        for segment in &mut segments {
            segment.same_sense = !segment.same_sense;
        }
    }
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

#[test]
fn a_thin_band_is_oriented_by_its_exact_area() {
    for width in [0.01, 0.005] {
        for reversed in [false, true] {
            let one = band(width, reversed);
            let other = band(width, reversed);
            let paired = pair(&[&one, &other]);
            for rings in &paired.rings {
                // A fine polygon of the true area: positive and near pi w.
                assert!(signed_area2(&rings.outer) > 0.0, "{width} {reversed}");
            }
        }
        // Band to band: every vertex matches a nearby one, not the far side.
        let (a, b) = (band(width, false), band(width / 2.0, false));
        let paired = pair(&[&a, &b]);
        let [ra, rb] = [&paired.rings[0].outer, &paired.rings[1].outer];
        assert_eq!(ra.len(), rb.len());
        assert!(signed_area2(ra) > 0.0 && signed_area2(rb) > 0.0);
        assert!(ra.iter().zip(rb).all(|(p, q)| (*p - *q).length() <= width));
    }
}

#[test]
fn a_scaled_derived_profile_charges_its_stretch() {
    let scaled = Profile::Derived {
        basis: Box::new(disk(0.5)),
        transform: Transform2::from_scale(Vec2::splat(2.0)),
    };
    // The stretched circle sets the count: the smaller one needs fewer steps.
    let other = disk(0.3);
    let paired = pair(&[&scaled, &other]);
    let ring = &paired.rings[0].outer;
    assert!(ring.iter().all(|p| (p.length() - 1.0).abs() < 1e-12));
    // The chords of the radius-1 circle really sag this far, and the
    // reported deviation covers it.
    let sagitta = 1.0 - (std::f64::consts::PI / ring.len() as Scalar).cos();
    assert!(
        paired.deviation >= sagitta,
        "{} < {sagitta}",
        paired.deviation
    );
    assert!(paired.deviation <= CHORD * 1.0001);
}

#[test]
fn the_half_step_is_only_for_runs_of_circles() {
    let (small, large) = (disk(0.5), disk(0.8));
    let paired = pair(&[&small, &large]);
    let ring = &paired.rings[0].outer;
    let half = std::f64::consts::PI / ring.len() as Scalar;
    assert!((ring[0].x - 0.5 * half.cos()).abs() < 1e-12);
    assert!((ring[0].y - 0.5 * half.sin()).abs() < 1e-12);

    let flat = ellipse(1.0, 0.5);
    let mixed = pair(&[&small, &flat]);
    let first = mixed.rings[0].outer[0];
    assert!(
        (first.x - 0.5).abs() < 1e-12 && first.y.abs() < 1e-12,
        "{first:?}"
    );
}

#[test]
fn one_segment_against_four_does_not_pair() {
    let round = disk(0.8);
    let square = sharp(1.0, 1.0);
    assert!(paired_profile_rings(&[&round, &square], CHORD, tolerance())
        .unwrap()
        .is_none());
    // A ring against a ring and a hole.
    let ring = annulus(1.0, 0.2);
    assert!(paired_profile_rings(&[&round, &ring], CHORD, tolerance())
        .unwrap()
        .is_none());
}

fn square_with(edge: impl Fn([Point2; 4]) -> ProfileSegment) -> Profile {
    let corners = [
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(0.0, 1.0),
    ];
    let mut segments: Vec<ProfileSegment> = (0..3)
        .map(|k| ProfileSegment {
            curve: Curve2::Line(Line2 {
                origin: corners[k],
                direction: corners[k + 1] - corners[k],
            }),
            domain: Interval::new(0.0, 1.0),
            same_sense: true,
        })
        .collect();
    segments.push(edge(corners));
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

#[test]
fn polyline_and_spline_segments_do_not_pair() {
    let polyline = square_with(|c| ProfileSegment {
        curve: Curve2::Polyline(Polyline2 {
            points: vec![c[3], Point2::new(-0.1, 0.5), c[0]],
            closed: false,
        }),
        domain: Interval::new(0.0, 2.0),
        same_sense: true,
    });
    let spline = square_with(|c| ProfileSegment {
        curve: Curve2::BSpline(BSplineCurve2 {
            degree: 1,
            control_points: vec![c[3], c[0]],
            knots: vec![0.0, 1.0],
            multiplicities: vec![2, 2],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        }),
        domain: Interval::new(0.0, 1.0),
        same_sense: true,
    });
    for profile in [&polyline, &spline] {
        assert!(
            paired_profile_rings(&[profile, profile], CHORD, tolerance())
                .unwrap()
                .is_none()
        );
    }
    // Other families do not pair either.
    let composite = Profile::Composite(vec![disk(0.5)]);
    assert!(
        paired_profile_rings(&[&composite, &composite], CHORD, tolerance())
            .unwrap()
            .is_none()
    );
}

#[test]
fn a_chord_error_that_is_not_positive_is_refused() {
    let round = disk(0.5);
    for bad in [0.0, -1.0, Scalar::NAN] {
        assert!(paired_profile_rings(&[&round, &round], bad, tolerance()).is_err());
    }
}
