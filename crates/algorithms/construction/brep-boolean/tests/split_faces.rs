//! Faces split along their section edges (#167, ADR 0075 stage 1).
//!
//! Oracles are areas from the inputs: a plane face's area is its area in
//! its own orthonormal parameters, a cylinder face's is the radius times
//! it. Regions must cover the face exactly (their areas sum to the face's)
//! and match the closed forms of the pieces the section cuts off.

use axiolid_brep::ExactBRep;
use axiolid_brep_boolean::{section_edges, split_face, Region, SectionEdge};
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, clip_arc_prism_exact, ArcPrism};
use axiolid_core::{BooleanOperator, Plane3, Point2, Point3, Tolerance, Vec3};
use axiolid_evaluate::curve::evaluate2;
use axiolid_overlay::ArcRing;
use axiolid_primitive::HalfSpace;
use axiolid_surface::Surface;
use axiolid_topology::FaceId;

const PI: f64 = std::f64::consts::PI;

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn prism(section: ArcRing, bottom: f64, top: f64) -> ArcPrism {
    ArcPrism {
        section,
        bottom,
        top,
    }
}

fn solid(section: ArcRing, bottom: f64, top: f64) -> ExactBRep {
    boolean_arc_prisms_exact(
        &prism(section.clone(), bottom, top),
        &prism(section, bottom, top),
        BooleanOperator::Intersection,
        tol(),
    )
    .expect("a solid")
}

fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> ArcRing {
    ArcRing::from_points(&[
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ])
}

fn surface_of(brep: &ExactBRep, face: FaceId) -> Surface {
    brep.surfaces()[brep.topology().faces()[face.index()]
        .surface
        .unwrap()
        .index()]
    .clone()
}

/// Signed parameter area of a loop, densely sampled along its exact pcurves.
fn loop_area(pieces: &[axiolid_brep_boolean::Piece]) -> f64 {
    let mut pts = Vec::new();
    for piece in pieces {
        let n = 4000;
        for i in 0..n {
            let p = piece.pspan.start + (piece.pspan.end - piece.pspan.start) * i as f64 / n as f64;
            pts.push(evaluate2(&piece.pcurve, p).unwrap());
        }
    }
    let mut area = 0.0;
    for i in 0..pts.len() {
        let (p, q) = (pts[i], pts[(i + 1) % pts.len()]);
        area += p.x * q.y - q.x * p.y;
    }
    0.5 * area
}

/// A region's area on its surface.
fn area(region: &Region, surface: &Surface) -> f64 {
    let scale = match surface {
        Surface::Plane(_) => 1.0,
        Surface::Cylinder(c) => c.radius,
        other => panic!("no area scale for {other:?}"),
    };
    let mut total = loop_area(&region.outer);
    for hole in &region.holes {
        total += loop_area(hole);
    }
    total * scale
}

/// Split every face of `brep` that carries sections; `other` is the second
/// operand, `first` whether `brep` is the first operand of `edges`.
fn split_all(
    brep: &ExactBRep,
    other: &ExactBRep,
    edges: &[SectionEdge],
    first: bool,
) -> Vec<(FaceId, Vec<Region>)> {
    let mut out = Vec::new();
    for index in 0..brep.topology().faces().len() {
        let face = brep.topology().face_id_at(index).unwrap();
        let mine: Vec<SectionEdge> = edges
            .iter()
            .filter(|e| {
                if first {
                    e.face_a == face
                } else {
                    e.face_b == face
                }
            })
            .cloned()
            .collect();
        let others: Vec<Surface> = mine
            .iter()
            .map(|e| surface_of(other, if first { e.face_b } else { e.face_a }))
            .collect();
        let regions = split_face(brep, face, &mine, &others, tol()).expect("split");
        out.push((face, regions));
    }
    out
}

fn close(what: &str, got: f64, want: f64, eps: f64) {
    assert!(
        (got - want).abs() <= eps,
        "{what}: expected {want}, got {got}"
    );
}

#[test]
fn a_box_roof_pierced_by_a_pipe_splits_into_a_ring_and_a_disc() {
    let block = solid(square(-1.0, -1.0, 1.0, 1.0), 0.0, 2.0);
    let r = 0.5;
    let pipe = solid(ArcRing::circle(Point2::new(0.2, -0.1), r), -1.0, 3.0);
    let edges = section_edges(&block, &pipe, tol()).expect("sections");
    for (face, regions) in split_all(&block, &pipe, &edges, true) {
        let surface = surface_of(&block, face);
        let total: f64 = regions.iter().map(|r| area(r, &surface)).sum();
        let pierced = edges.iter().any(|e| e.face_a == face);
        if pierced {
            // Floor or roof: 4 = (4 - pi r^2) + pi r^2.
            assert_eq!(regions.len(), 2, "{regions:#?}");
            let mut areas: Vec<f64> = regions.iter().map(|r| area(r, &surface)).collect();
            areas.sort_by(f64::total_cmp);
            close("disc", areas[0], PI * r * r, 1e-5);
            close("ring", areas[1], 4.0 - PI * r * r, 1e-5);
            assert!(
                regions.iter().any(|r| r.holes.len() == 1),
                "the ring has the hole"
            );
        } else {
            assert_eq!(regions.len(), 1);
        }
        // Split or not, the regions cover the 2 x 2 face exactly.
        close("face total", total, 4.0, 1e-5);
    }
    // The pipe's two half-walls each split into three bands, of heights
    // 1, 2 and 1.
    for (face, regions) in split_all(&pipe, &block, &edges, false) {
        let surface = surface_of(&pipe, face);
        if !matches!(surface, Surface::Cylinder(_)) {
            continue;
        }
        let mut areas: Vec<f64> = regions.iter().map(|g| area(g, &surface)).collect();
        areas.sort_by(f64::total_cmp);
        assert_eq!(areas.len(), 3, "{areas:?}");
        close("low band", areas[0], PI * r * 1.0, 1e-5);
        close("high band", areas[1], PI * r * 1.0, 1e-5);
        close("middle band", areas[2], PI * r * 2.0, 1e-5);
    }
}

#[test]
fn a_sloped_roof_and_the_pipe_through_it_split_exactly() {
    let column = prism(square(-2.0, -2.0, 2.0, 2.0), 0.0, 10.0);
    let normal = Vec3::new(-0.4, 0.2, 1.0);
    let roof = HalfSpace {
        boundary: Plane3 {
            origin: Point3::new(0.0, 0.0, 3.0),
            normal,
        },
        agreement: false,
    };
    let block = clip_arc_prism_exact(&column, &roof, tol()).expect("a sloped block");
    let (cx, cy, r) = (0.3, 0.2, 0.6);
    let pipe = solid(ArcRing::circle(Point2::new(cx, cy), r), -1.0, 12.0);
    let edges = section_edges(&block, &pipe, tol()).expect("sections");
    // On the roof: an ellipse of area pi r^2 / cos(theta).
    let cos = normal.normalize().z;
    for (face, regions) in split_all(&block, &pipe, &edges, true) {
        let surface = surface_of(&block, face);
        let Surface::Plane(plane) = &surface else {
            continue;
        };
        if plane.frame.z.normalize().z.abs() > 0.99 {
            continue;
        }
        if !edges.iter().any(|e| e.face_a == face) {
            continue;
        }
        let mut areas: Vec<f64> = regions.iter().map(|g| area(g, &surface)).collect();
        areas.sort_by(f64::total_cmp);
        assert_eq!(areas.len(), 2);
        close("ellipse", areas[0], PI * r * r / cos, 1e-5);
    }
    // The pipe's walls below the roof: the roof height over the pipe's
    // circle averages its height at the centre, 3 + 0.4 cx - 0.2 cy, and
    // the wall starts at z = -1.
    let mean = 3.0 + 0.4 * cx - 0.2 * cy;
    let mut below = 0.0;
    for (face, regions) in split_all(&pipe, &block, &edges, false) {
        let surface = surface_of(&pipe, face);
        if !matches!(surface, Surface::Cylinder(_)) {
            continue;
        }
        for region in &regions {
            // A region below the roof has its boundary, on average, under
            // the roof's mean height.
            let mid_v = region
                .outer
                .iter()
                .map(|p| {
                    evaluate2(&p.pcurve, 0.5 * (p.pspan.start + p.pspan.end))
                        .unwrap()
                        .y
                })
                .sum::<f64>()
                / region.outer.len() as f64;
            if mid_v < mean {
                below += area(region, &surface);
            }
        }
    }
    close(
        "wall below the roof",
        below,
        2.0 * PI * r * (mean + 1.0),
        1e-4,
    );
}

#[test]
fn a_box_roof_cut_at_its_corner_keeps_its_area() {
    let block = solid(square(-1.0, -1.0, 1.0, 1.0), 0.0, 2.0);
    let (cx, cy, r) = (1.0, 0.9, 0.5);
    let pipe = solid(ArcRing::circle(Point2::new(cx, cy), r), -1.0, 3.0);
    let edges = section_edges(&block, &pipe, tol()).expect("sections");
    for (face, regions) in split_all(&block, &pipe, &edges, true) {
        let surface = surface_of(&block, face);
        let total: f64 = regions.iter().map(|g| area(g, &surface)).sum();
        let Surface::Plane(plane) = &surface else {
            continue;
        };
        if plane.frame.z.z.abs() > 0.99 && edges.iter().any(|e| e.face_a == face) {
            assert_eq!(regions.len(), 2, "{regions:#?}");
            // The corner piece: the part of the disc over the box, by a
            // fine grid over its bounding square.
            let n = 2000;
            let mut cut = 0.0;
            for i in 0..n {
                for j in 0..n {
                    let x = cx - r + 2.0 * r * (i as f64 + 0.5) / n as f64;
                    let y = cy - r + 2.0 * r * (j as f64 + 0.5) / n as f64;
                    if (x - cx).hypot(y - cy) <= r && x <= 1.0 && y <= 1.0 {
                        cut += (2.0 * r / n as f64).powi(2);
                    }
                }
            }
            let mut areas: Vec<f64> = regions.iter().map(|g| area(g, &surface)).collect();
            areas.sort_by(f64::total_cmp);
            close("corner piece", areas[0], cut, 2e-3);
        }
        close("face total", total, 4.0, 1e-5);
    }
}
