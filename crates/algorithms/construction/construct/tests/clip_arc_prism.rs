//! A curved prism cut by a sloped plane (#120 / ADR 0071).
//!
//! The claim: the cut is exact. Each cylindrical wall stays a `Cylinder`
//! whose sloped rim is an `Ellipse3` edge with a `Sinusoid2` pcurve, and the
//! solid passes the geometric audit, which checks every pcurve against its
//! 3D edge.
//!
//! The volume oracle is closed form, not read from the kernel: cutting a
//! vertical prism of section `S` with the plane `z = h + g . p` leaves
//! `area(S) * (top - h - g . centroid(S))` of material below the top cap,
//! because the height of the removed wedge is affine in plan position.
//! The measured side is the divergence theorem on the B-rep's own cap
//! loops -- the walls are vertical and contribute nothing to `V = int z n_z`.

use axiolid_brep::ExactBRep;
use axiolid_construct::boolean_exact::{clip_arc_prism_exact, ArcPrism};
use axiolid_contracts::GeomError;
use axiolid_core::{Plane3, Point2, Point3, Tolerance, Vec3};
use axiolid_curve::{Curve2, Curve3};
use axiolid_overlay::{ArcRing, ArcVertex};
use axiolid_primitive::HalfSpace;
use axiolid_surface::Surface;
use axiolid_topology::Orientation;

const PI: f64 = std::f64::consts::PI;

fn disc(cx: f64, cy: f64, r: f64) -> ArcRing {
    ArcRing::circle(Point2::new(cx, cy), r)
}

/// A `w` x `h` rectangle whose right side is a half-disc bulging outward:
/// straight AND curved walls in one section.
fn stadium_end(w: f64, h: f64) -> ArcRing {
    ArcRing {
        vertices: vec![
            ArcVertex::straight(Point2::new(0.0, 0.0)),
            ArcVertex::bulged(Point2::new(w, 0.0), 1.0),
            ArcVertex::straight(Point2::new(w, h)),
            ArcVertex::straight(Point2::new(0.0, h)),
        ],
    }
}

fn prism(section: ArcRing, bottom: f64, top: f64) -> ArcPrism {
    ArcPrism {
        section,
        bottom,
        top,
    }
}

/// The plane `z = height + gx x + gy y`, keeping the side below it.
fn below(height: f64, gx: f64, gy: f64) -> HalfSpace {
    HalfSpace {
        boundary: Plane3 {
            origin: Point3::new(0.0, 0.0, height),
            normal: Vec3::new(-gx, -gy, 1.0),
        },
        agreement: false,
    }
}

/// The same plane, keeping the side above it.
fn above(height: f64, gx: f64, gy: f64) -> HalfSpace {
    HalfSpace {
        agreement: true,
        ..below(height, gx, gy)
    }
}

/// `V = sum over cap faces of int z n_z dA`, from the cap loops.
///
/// For a planar face, `int z n_z dA = (n_z / |n|) * int z dA`, and with the
/// face as a polygon (curved edges sampled densely) both the area vector
/// and the first moment come from a fan. Walls are vertical: n_z = 0.
fn volume_from_caps(solid: &ExactBRep) -> f64 {
    let topo = solid.topology();
    let mut total = 0.0;
    for face in topo.faces() {
        let surface = &solid.surfaces()[face.surface.unwrap().index()];
        let Surface::Plane(plane) = surface else {
            continue;
        };
        if plane.frame.z.z.abs() < 1e-12 {
            continue; // a vertical planar wall
        }
        for bound in &face.bounds {
            let ring = loop_points(solid, bound.loop_id.index());
            // Fan about the first point: area vector and z-weighted area.
            let o = ring[0];
            let mut area_vec = Vec3::ZERO;
            let mut z_moment = 0.0;
            for w in ring.windows(2).skip(1) {
                let cross = (w[0] - o).cross(w[1] - o) * 0.5;
                area_vec += cross;
                z_moment += cross.z * (o.z + w[0].z + w[1].z) / 3.0;
            }
            let sense = match face.orientation {
                Orientation::Forward => 1.0,
                Orientation::Reversed => -1.0,
            };
            let _ = area_vec;
            total += sense * z_moment;
        }
    }
    total.abs()
}

/// Dense points around a loop, in traversal order, closing on itself.
fn loop_points(solid: &ExactBRep, loop_index: usize) -> Vec<Point3> {
    let topo = solid.topology();
    let lp = &topo.loops()[loop_index];
    let mut out = Vec::new();
    for use_ in &lp.edges {
        let edge = &topo.edges()[use_.edge.index()];
        let curve = &solid.curves3()[edge.curve.unwrap().index()];
        let span = solid.edge_interval(use_.edge).unwrap();
        let n = match curve {
            Curve3::Line(_) => 1,
            _ => 4096,
        };
        for i in 0..n {
            let f = i as f64 / n as f64;
            let f = match use_.orientation {
                Orientation::Forward => f,
                Orientation::Reversed => 1.0 - f,
            };
            let t = span.start + (span.end - span.start) * f;
            out.push(axiolid_evaluate::curve::evaluate3(curve, t).unwrap());
        }
    }
    let first = out[0];
    out.push(first);
    out
}

fn audit(solid: &ExactBRep) {
    let health = axiolid_brep_audit::geometric_audit(solid, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
}

fn count_surfaces(solid: &ExactBRep) -> (usize, usize) {
    let mut planes = 0;
    let mut cylinders = 0;
    for s in solid.surfaces() {
        match s {
            Surface::Plane(_) => planes += 1,
            Surface::Cylinder(_) => cylinders += 1,
            other => panic!("unexpected surface {other:?}"),
        }
    }
    (planes, cylinders)
}

fn sinusoid_pcurves(solid: &ExactBRep) -> usize {
    solid
        .curves2()
        .iter()
        .filter(|c| matches!(c, Curve2::Sinusoid(_)))
        .count()
}

fn ellipse_edges(solid: &ExactBRep) -> usize {
    solid
        .curves3()
        .iter()
        .filter(|c| matches!(c, Curve3::Ellipse(_)))
        .count()
}

#[test]
fn a_column_under_a_sloped_roof_keeps_its_cylinder_and_cuts_an_ellipse() {
    // Radius-1 column from 0 to 10, roof z = 6 + 0.5 x through its axis.
    let column = prism(disc(0.0, 0.0, 1.0), 0.0, 10.0);
    let solid = clip_arc_prism_exact(&column, &below(6.0, 0.5, 0.0), Tolerance::METRE)
        .expect("a roof cutting cleanly through the column is representable");
    audit(&solid);

    let (planes, cylinders) = count_surfaces(&solid);
    assert_eq!(planes, 2, "a flat bottom and a sloped top cap");
    assert!(
        cylinders >= 1,
        "the wall stays a cylinder, not a ruled face"
    );
    assert!(
        ellipse_edges(&solid) >= 1,
        "the roof cuts the wall in an ellipse"
    );
    assert!(
        sinusoid_pcurves(&solid) >= 1,
        "trimmed on the wall by the wave"
    );

    // Symmetric section: the centroid is the axis, so the mean height is 6.
    let expected = PI * 6.0;
    let got = volume_from_caps(&solid);
    assert!(
        (got - expected).abs() < 1e-5 * expected,
        "volume {got}, expected {expected}"
    );
    // Highest and lowest points of the cut are at x = +1 and x = -1.
    let zs: Vec<f64> = solid
        .topology()
        .vertices()
        .iter()
        .map(|v| v.position.z)
        .collect();
    assert!(zs.iter().all(|&z| (0.0..=6.5 + 1e-12).contains(&z)));
}

#[test]
fn an_off_axis_sloped_roof_matches_the_closed_form_volume() {
    // Centroid off the origin, and a gradient in both x and y, so the mean
    // height is NOT the height at the origin.
    let column = prism(disc(2.0, -1.0, 1.5), 1.0, 20.0);
    let (h, gx, gy) = (8.0, 0.3, -0.4);
    let solid =
        clip_arc_prism_exact(&column, &below(h, gx, gy), Tolerance::METRE).expect("representable");
    audit(&solid);
    let area = PI * 1.5 * 1.5;
    let mean = h + gx * 2.0 - gy;
    let expected = area * (mean - 1.0);
    let got = volume_from_caps(&solid);
    assert!(
        (got - expected).abs() < 1e-5 * expected,
        "volume {got}, expected {expected}"
    );
}

#[test]
fn a_mixed_section_gets_sloped_lines_on_planar_walls_and_ellipses_on_curved_ones() {
    // 4 x 2 rectangle with a half-disc (radius 1, centre (4,1)) on the right.
    let section = stadium_end(4.0, 2.0);
    let solid = clip_arc_prism_exact(
        &prism(section, 0.0, 10.0),
        &below(5.0, 0.2, 0.1),
        Tolerance::METRE,
    )
    .expect("representable");
    audit(&solid);
    assert_eq!(ellipse_edges(&solid), 1, "one arc edge, one ellipse cut");
    assert_eq!(sinusoid_pcurves(&solid), 1, "one sloped cylinder rim");
    // area = 8 + pi/2; centroid from the rectangle (2, 1) and the half-disc
    // (4 + 4/(3 pi), 1), area-weighted.
    let (a_rect, a_half) = (8.0, PI / 2.0);
    let cx = (a_rect * 2.0 + a_half * (4.0 + 4.0 / (3.0 * PI))) / (a_rect + a_half);
    let cy = 1.0;
    let expected = (a_rect + a_half) * (5.0 + 0.2 * cx + 0.1 * cy);
    let got = volume_from_caps(&solid);
    assert!(
        (got - expected).abs() < 1e-5 * expected,
        "volume {got}, expected {expected}"
    );
}

#[test]
fn keeping_the_upper_side_replaces_the_bottom_cap() {
    let column = prism(disc(0.0, 0.0, 1.0), 0.0, 10.0);
    let solid = clip_arc_prism_exact(&column, &above(3.0, 0.5, 0.0), Tolerance::METRE)
        .expect("a sloped floor through the column is representable");
    audit(&solid);
    assert!(sinusoid_pcurves(&solid) >= 1);
    let expected = PI * (10.0 - 3.0);
    let got = volume_from_caps(&solid);
    assert!(
        (got - expected).abs() < 1e-5 * expected,
        "volume {got}, expected {expected}"
    );
}

#[test]
fn a_downward_normal_selects_the_same_side_as_its_flipped_twin() {
    // The same plane written with the opposite normal, kept on the side that
    // normal points into, is the same half-space. The results must agree.
    let column = prism(disc(0.0, 0.0, 1.0), 0.0, 10.0);
    let flipped = HalfSpace {
        boundary: Plane3 {
            origin: Point3::new(0.0, 0.0, 6.0),
            normal: Vec3::new(0.5, 0.0, -1.0),
        },
        agreement: true,
    };
    let a = clip_arc_prism_exact(&column, &below(6.0, 0.5, 0.0), Tolerance::METRE).unwrap();
    let b = clip_arc_prism_exact(&column, &flipped, Tolerance::METRE).unwrap();
    let (va, vb) = (volume_from_caps(&a), volume_from_caps(&b));
    assert!((va - vb).abs() < 1e-9 * va, "{va} vs {vb}");
}

#[test]
fn the_cut_edge_lies_on_both_the_cylinder_and_the_plane() {
    // Independent of the audit: sample every ellipse edge along its span
    // and check each point is on the radius-1 cylinder and on the roof.
    let column = prism(disc(0.0, 0.0, 1.0), 0.0, 10.0);
    let solid = clip_arc_prism_exact(&column, &below(6.0, 0.5, 0.25), Tolerance::METRE).unwrap();
    let topo = solid.topology();
    let mut checked = 0;
    for (i, edge) in topo.edges().iter().enumerate() {
        let curve = &solid.curves3()[edge.curve.unwrap().index()];
        if !matches!(curve, Curve3::Ellipse(_)) {
            continue;
        }
        let span = solid.edge_interval(topo.edge_id_at(i).unwrap()).unwrap();
        for k in 0..=64 {
            let t = span.start + (span.end - span.start) * k as f64 / 64.0;
            let p = axiolid_evaluate::curve::evaluate3(curve, t).unwrap();
            let r = (p.x * p.x + p.y * p.y).sqrt();
            assert!((r - 1.0).abs() < 1e-12, "off the cylinder: r = {r}");
            let roof = 6.0 + 0.5 * p.x + 0.25 * p.y;
            assert!((p.z - roof).abs() < 1e-12, "off the roof by {}", p.z - roof);
            checked += 1;
        }
        // The edge must run from one of its vertices to the other.
        let start = axiolid_evaluate::curve::evaluate3(curve, span.start).unwrap();
        let end = axiolid_evaluate::curve::evaluate3(curve, span.end).unwrap();
        let vs = topo.vertices();
        assert!((start - vs[edge.start.index()].position).length() < 1e-12);
        assert!((end - vs[edge.end.index()].position).length() < 1e-12);
    }
    assert!(checked > 0, "no ellipse edge was sampled");
}

#[test]
fn a_plane_crossing_a_cap_inside_the_section_is_refused_by_name() {
    // z = 0.5 + x crosses z = 0 at x = -0.5, inside the unit disc.
    let column = prism(disc(0.0, 0.0, 1.0), 0.0, 10.0);
    let err = clip_arc_prism_exact(&column, &below(0.5, 1.0, 0.0), Tolerance::METRE)
        .expect_err("a plane that crosses the bottom cap leaves two kinds of cap");
    assert!(
        matches!(err, GeomError::UnsupportedInput { input, .. } if input.contains("crosses a cap")),
        "{err:?}"
    );
}

#[test]
fn a_vertical_plane_is_refused_by_name() {
    let column = prism(disc(0.0, 0.0, 1.0), 0.0, 10.0);
    let wall = HalfSpace {
        boundary: Plane3 {
            origin: Point3::ZERO,
            normal: Vec3::X,
        },
        agreement: true,
    };
    let err = clip_arc_prism_exact(&column, &wall, Tolerance::METRE).unwrap_err();
    assert!(
        matches!(err, GeomError::UnsupportedInput { input, .. } if input.contains("parallel")),
        "{err:?}"
    );
}

#[test]
fn a_plane_that_keeps_everything_or_nothing_is_decided_by_the_true_range() {
    // Radius 1 around the origin, gradient 0.5: the plane varies by +-0.5
    // over the disc. At height 10.4 it dips to 9.9 at x = -1 -- inside the
    // prism -- so it must cut, not be treated as clear of the top.
    let column = prism(disc(0.0, 0.0, 1.0), 0.0, 10.0);
    let err = clip_arc_prism_exact(&column, &below(10.4, 0.5, 0.0), Tolerance::METRE)
        .expect_err("dips through the top cap");
    assert!(matches!(err, GeomError::UnsupportedInput { .. }), "{err:?}");
    // At 10.6 the lowest point is 10.1: clear, the whole prism is kept.
    let whole = clip_arc_prism_exact(&column, &below(10.6, 0.5, 0.0), Tolerance::METRE)
        .expect("clear of the prism");
    assert_eq!(sinusoid_pcurves(&whole), 0);
    let got = volume_from_caps(&whole);
    assert!((got - PI * 10.0).abs() < 1e-5 * PI * 10.0, "{got}");
    // Entirely below the prism on the kept side: nothing is left.
    let err =
        clip_arc_prism_exact(&column, &below(-0.6, 0.5, 0.0), Tolerance::METRE).expect_err("empty");
    assert!(matches!(err, GeomError::Degenerate(_)), "{err:?}");
}

#[test]
fn only_the_uncut_cap_keeps_its_name() {
    use axiolid_brep::{FaceName, SweptFace};
    let column = prism(disc(0.0, 0.0, 1.0), 0.0, 10.0);
    let solid = clip_arc_prism_exact(&column, &below(6.0, 0.5, 0.0), Tolerance::METRE).unwrap();
    let faces = solid.topology().faces();
    let names: Vec<_> = (0..faces.len())
        .filter_map(|i| solid.face_name(solid.topology().face_id_at(i).unwrap()))
        .collect();
    assert!(names.contains(&&FaceName::swept(SweptFace::StartCap)));
    assert!(
        !names.contains(&&FaceName::swept(SweptFace::EndCap)),
        "the cut cap is the roof's, not the column's end cap"
    );
}

#[test]
fn the_range_includes_extremes_inside_an_arc_not_only_its_vertices() {
    // `ArcRing::circle` puts its two vertices at (+-1, 0). A plane sloped
    // along y is level at both vertices but dips 0.5 below them at (0, -1)
    // and rises 0.5 above at (0, 1): only the in-arc extremes show that a
    // plane at 10.4 dips through the top cap at 10.
    let column = prism(disc(0.0, 0.0, 1.0), 0.0, 10.0);
    let err = clip_arc_prism_exact(&column, &below(10.4, 0.0, 0.5), Tolerance::METRE)
        .expect_err("dips to 9.9 inside the arc");
    assert!(
        matches!(err, GeomError::UnsupportedInput { input, .. } if input.contains("crosses a cap")),
        "{err:?}"
    );
    // The same on the kept-above side: level 0.4 at the vertices, but it
    // dips to -0.1 inside the arc, through the bottom cap.
    let err = clip_arc_prism_exact(&column, &above(0.4, 0.0, 0.5), Tolerance::METRE)
        .expect_err("dips to -0.1 inside the arc");
    assert!(matches!(err, GeomError::UnsupportedInput { .. }), "{err:?}");
}

#[test]
fn keeping_the_upper_side_refuses_a_plane_through_the_bottom_cap() {
    // z = 0.3 + 0.5 x spans -0.2 .. 0.8 over the disc: it crosses z = 0.
    let column = prism(disc(0.0, 0.0, 1.0), 0.0, 10.0);
    let err = clip_arc_prism_exact(&column, &above(0.3, 0.5, 0.0), Tolerance::METRE)
        .expect_err("crosses the bottom cap");
    assert!(
        matches!(err, GeomError::UnsupportedInput { input, .. } if input.contains("crosses a cap")),
        "{err:?}"
    );
}

#[test]
fn keeping_the_upper_side_names_only_the_top_cap() {
    use axiolid_brep::{FaceName, SweptFace};
    let column = prism(disc(0.0, 0.0, 1.0), 0.0, 10.0);
    let solid = clip_arc_prism_exact(&column, &above(3.0, 0.5, 0.0), Tolerance::METRE).unwrap();
    let names: Vec<_> = (0..solid.topology().faces().len())
        .filter_map(|i| solid.face_name(solid.topology().face_id_at(i).unwrap()))
        .collect();
    assert!(names.contains(&&FaceName::swept(SweptFace::EndCap)));
    assert!(
        !names.contains(&&FaceName::swept(SweptFace::StartCap)),
        "the cut floor is the half-space's, not the column's start cap"
    );
}

#[test]
fn a_concave_arc_edge_is_cut_in_the_right_direction() {
    // A 4 x 2 block with a half-disc notch (radius 1, centre (2, 2)) bitten
    // out of its top edge. The notch is a CLOCKWISE arc inside a
    // counter-clockwise ring, so its sweep is negative and the cut ellipse
    // is walked backwards; the audit and the edge-endpoint check catch a
    // reversed span.
    let section = ArcRing {
        vertices: vec![
            ArcVertex::straight(Point2::new(0.0, 0.0)),
            ArcVertex::straight(Point2::new(4.0, 0.0)),
            ArcVertex::straight(Point2::new(4.0, 2.0)),
            ArcVertex::bulged(Point2::new(3.0, 2.0), -1.0),
            ArcVertex::straight(Point2::new(1.0, 2.0)),
            ArcVertex::straight(Point2::new(0.0, 2.0)),
        ],
    };
    let solid = clip_arc_prism_exact(
        &prism(section, 0.0, 10.0),
        &below(5.0, 0.2, 0.1),
        Tolerance::METRE,
    )
    .expect("a notched block under a sloped roof is representable");
    audit(&solid);
    assert_eq!(ellipse_edges(&solid), 1);
    // area = 8 - pi/2; the notch's centroid sits 4/(3 pi) below y = 2.
    let (a_rect, a_notch) = (8.0, PI / 2.0);
    let area = a_rect - a_notch;
    let cx = 2.0;
    let cy = (a_rect - a_notch * (2.0 - 4.0 / (3.0 * PI))) / area;
    let expected = area * (5.0 + 0.2 * cx + 0.1 * cy);
    let got = volume_from_caps(&solid);
    assert!(
        (got - expected).abs() < 1e-5 * expected,
        "volume {got}, expected {expected}"
    );
}
