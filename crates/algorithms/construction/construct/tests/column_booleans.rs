//! Stepped coaxial booleans: results that are not one prism (#120).
//!
//! A union of prisms with different heights, or a difference whose tool
//! stops partway up, leaves a solid whose cross-section changes with
//! height. These are built exactly as vertical columns over one planar
//! arrangement (ADR 0072): a curved wall stays a `Cylinder`, and a ledge
//! is a planar face at the height where the section changes.
//!
//! Oracles are closed forms from the inputs, never read back from the
//! kernel. Planar results are measured by `exact_properties`; curved ones
//! by the divergence theorem on the caps (walls are vertical, so only the
//! caps carry `int z n_z dA`).

use axiolid_brep::{ExactBRep, FaceName, Operand, SweptFace};
use axiolid_construct::boolean_exact::{
    boolean_arc_prisms_exact, boolean_arc_prisms_exact_solids, boolean_prisms_exact,
    boolean_prisms_exact_solids, ArcPrism, Prism,
};
use axiolid_contracts::GeomError;
use axiolid_core::{BooleanOperator, Point2, Point3, Tolerance};
use axiolid_curve::Curve3;
use axiolid_overlay::ArcRing;
use axiolid_surface::Surface;
use axiolid_topology::Orientation;

const PI: f64 = std::f64::consts::PI;

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Point2> {
    vec![
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ]
}

fn prism(ring: Vec<Point2>, bottom: f64, top: f64) -> Prism {
    Prism {
        rings: vec![ring],
        bottom,
        top,
    }
}

fn disc(cx: f64, r: f64, bottom: f64, top: f64) -> ArcPrism {
    ArcPrism {
        section: ArcRing::circle(Point2::new(cx, 0.0), r),
        bottom,
        top,
    }
}

fn planar_volume(solid: &ExactBRep) -> f64 {
    axiolid_measure::exact_properties(solid, tol())
        .expect("all-planar solid is measurable")
        .signed_volume
}

/// `V = sum over caps of int z n_z dA`, curved cap edges sampled densely.
fn cap_volume(solid: &ExactBRep) -> f64 {
    let topo = solid.topology();
    let mut total = 0.0;
    for face in topo.faces() {
        let Surface::Plane(plane) = &solid.surfaces()[face.surface.unwrap().index()] else {
            continue;
        };
        if plane.frame.z.z.abs() < 1e-12 {
            continue;
        }
        for bound in &face.bounds {
            let lp = &topo.loops()[bound.loop_id.index()];
            let mut ring: Vec<Point3> = Vec::new();
            for use_ in &lp.edges {
                let edge = &topo.edges()[use_.edge.index()];
                let curve = &solid.curves3()[edge.curve.unwrap().index()];
                let span = solid.edge_interval(use_.edge).unwrap();
                let n = if matches!(curve, Curve3::Line(_)) {
                    1
                } else {
                    4096
                };
                for i in 0..n {
                    let f = i as f64 / n as f64;
                    let f = match use_.orientation {
                        Orientation::Forward => f,
                        Orientation::Reversed => 1.0 - f,
                    };
                    let t = span.start + (span.end - span.start) * f;
                    ring.push(axiolid_evaluate::curve::evaluate3(curve, t).unwrap());
                }
            }
            // Loops are wound in the plane's own frame (counter-clockwise
            // seen from +z); a face used `Reversed` points down, so its
            // contribution flips. These builders use every face forward in
            // its shell, so the face orientation is the whole story.
            let sense = match face.orientation {
                Orientation::Forward => 1.0,
                Orientation::Reversed => -1.0,
            };
            let o = ring[0];
            let mut moment = 0.0;
            for w in ring[1..].windows(2) {
                let cross = (w[0] - o).cross(w[1] - o) * 0.5;
                moment += cross.z * (o.z + w[0].z + w[1].z) / 3.0;
            }
            total += sense * moment;
        }
    }
    total
}

fn audit(solid: &ExactBRep) {
    let health = axiolid_brep_audit::geometric_audit(solid, tol());
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topo = axiolid_topology::audit_brep(solid.topology());
    assert!(topo.is_closed_manifold(), "{topo:?}");
}

fn close(got: f64, want: f64, rel: f64) {
    assert!(
        (got - want).abs() <= rel * want.abs(),
        "got {got}, want {want}"
    );
}

fn cylinders(solid: &ExactBRep) -> usize {
    solid
        .surfaces()
        .iter()
        .filter(|s| matches!(s, Surface::Cylinder(_)))
        .count()
}

fn names(solid: &ExactBRep) -> Vec<FaceName> {
    let topo = solid.topology();
    (0..topo.faces().len())
        .filter_map(|i| solid.face_name(topo.face_id_at(i).unwrap()).cloned())
        .collect()
}

/// Lens area of two unit discs whose centres are `d` apart.
fn lens(d: f64) -> f64 {
    2.0 * (d / 2.0).acos() - (d / 2.0) * (4.0 - d * d).sqrt()
}

#[test]
fn a_slab_and_a_tower_unite_into_one_stepped_solid() {
    // 4 x 4 x 1 slab, 2 x 4 x 5 tower over its right half:
    // V = 16 * 1 + 8 * (5 - 1) = 48.
    let slab = prism(rect(0.0, 0.0, 4.0, 4.0), 0.0, 1.0);
    let tower = prism(rect(2.0, 0.0, 4.0, 4.0), 0.0, 5.0);
    let solid = boolean_prisms_exact(&slab, &tower, BooleanOperator::Union, tol())
        .expect("a stepped union is one solid");
    audit(&solid);
    close(planar_volume(&solid), 48.0, 1e-12);
    // The ledge at z = 1 is the slab's own top cap, cut back by the tower.
    assert!(names(&solid).contains(&FaceName::swept(SweptFace::EndCap).fragment(Operand::Subject)));
    assert!(names(&solid).contains(&FaceName::swept(SweptFace::EndCap).fragment(Operand::Tool)));
}

#[test]
fn a_short_tool_cuts_a_step_out_of_one_end() {
    // 10 x 4 x 3 minus a 2 x 2 x 1.5 corner notch: 120 - 6 = 114.
    let block = prism(rect(0.0, 0.0, 10.0, 4.0), 0.0, 3.0);
    let notch = prism(rect(0.0, 0.0, 2.0, 2.0), 0.0, 1.5);
    let solid = boolean_prisms_exact(&block, &notch, BooleanOperator::Difference, tol())
        .expect("a notch leaves one solid");
    audit(&solid);
    close(planar_volume(&solid), 114.0, 1e-12);
    // The notch's ceiling is the tool's top cap, facing down into the notch.
    assert!(names(&solid).contains(&FaceName::swept(SweptFace::EndCap).fragment(Operand::Tool)));
}

#[test]
fn a_blind_pocket_from_the_top_leaves_a_floor() {
    // Pocket 8 x 2, from z = 2 up through the top at 3: 120 - 16 = 104.
    let block = prism(rect(0.0, 0.0, 10.0, 4.0), 0.0, 3.0);
    let pocket = prism(rect(1.0, 1.0, 9.0, 3.0), 2.0, 5.0);
    let solid = boolean_prisms_exact(&block, &pocket, BooleanOperator::Difference, tol())
        .expect("a pocket leaves one solid");
    audit(&solid);
    close(planar_volume(&solid), 104.0, 1e-12);
    // The top cap has a hole where the pocket opens.
    let top_bounds = solid
        .topology()
        .faces()
        .iter()
        .filter(|f| f.orientation == Orientation::Forward)
        .map(|f| f.bounds.len())
        .max()
        .unwrap();
    assert!(top_bounds >= 2, "the pocket must open through the top cap");
}

#[test]
fn the_multi_solid_path_agrees_with_the_single_solid_path() {
    let slab = prism(rect(0.0, 0.0, 4.0, 4.0), 0.0, 1.0);
    let tower = prism(rect(2.0, 0.0, 4.0, 4.0), 0.0, 5.0);
    let one = boolean_prisms_exact(&slab, &tower, BooleanOperator::Union, tol()).unwrap();
    let many = boolean_prisms_exact_solids(&slab, &tower, BooleanOperator::Union, tol()).unwrap();
    assert_eq!(many.len(), 1);
    close(planar_volume(&many[0]), planar_volume(&one), 1e-15);
}

#[test]
fn a_short_round_column_joined_to_a_tall_one_keeps_both_cylinders() {
    // Unit discs 0.5 apart, heights 1 and 5:
    // V = 5 pi (tall) + (pi - lens) * 1 (the short one's own part).
    let short = disc(0.0, 1.0, 0.0, 1.0);
    let tall = disc(0.5, 1.0, 0.0, 5.0);
    let solid = boolean_arc_prisms_exact(&short, &tall, BooleanOperator::Union, tol())
        .expect("a stepped curved union is one solid");
    audit(&solid);
    close(cap_volume(&solid), 5.0 * PI + (PI - lens(0.5)), 1e-6);
    // Two walls on the tall disc (one below the ledge, one above) would be
    // fine; zero would mean the wall was faceted.
    assert!(cylinders(&solid) >= 2, "both discs keep cylindrical walls");
}

#[test]
fn a_counterbore_is_a_stepped_round_hole() {
    // Radius-2 column, height 3, minus a radius-1 bore from z = 2 upward:
    // V = 4 pi * 3 - pi * 1 = 11 pi.
    let column = disc(0.0, 2.0, 0.0, 3.0);
    let bore = disc(0.0, 1.0, 2.0, 4.0);
    let solid = boolean_arc_prisms_exact(&column, &bore, BooleanOperator::Difference, tol())
        .expect("a counterbore is one solid");
    audit(&solid);
    close(cap_volume(&solid), 11.0 * PI, 1e-6);
    // Outer wall of the column, and the bore's own wall.
    let radii: Vec<f64> = solid
        .surfaces()
        .iter()
        .filter_map(|s| match s {
            Surface::Cylinder(c) => Some(c.radius),
            _ => None,
        })
        .collect();
    assert!(radii.iter().any(|r| (r - 2.0).abs() < 1e-12), "{radii:?}");
    assert!(radii.iter().any(|r| (r - 1.0).abs() < 1e-12), "{radii:?}");
}

#[test]
fn a_stepped_result_that_falls_apart_is_every_piece() {
    // A tall disc minus a strip that only reaches z = 1: below z = 1 two
    // D-shaped halves, above it the whole disc. One solid, not two.
    let column = disc(0.0, 1.0, 0.0, 3.0);
    let strip = ArcPrism {
        section: ArcRing::new(
            rect(-0.25, -2.0, 0.25, 2.0)
                .into_iter()
                .map(axiolid_overlay::ArcVertex::straight)
                .collect(),
        ),
        bottom: 0.0,
        top: 1.0,
    };
    let pieces =
        boolean_arc_prisms_exact_solids(&column, &strip, BooleanOperator::Difference, tol())
            .expect("stepped");
    assert_eq!(pieces.len(), 1, "the upper part joins both legs");
    audit(&pieces[0]);
    // Strip over the disc: |x| <= 0.25 inside the unit disc,
    //   area = 2 int_{-0.25}^{0.25} sqrt(1 - x^2) dx.
    let a = 0.25f64;
    let strip_area = 2.0 * (a * (1.0 - a * a).sqrt() + a.asin());
    close(cap_volume(&pieces[0]), 3.0 * PI - strip_area, 1e-6);
}

#[test]
fn solids_touching_only_along_an_edge_are_refused_by_name() {
    // Squares meeting at one corner in plan, different heights: the union
    // shares a single vertical edge, which is not a manifold solid.
    let a = prism(rect(0.0, 0.0, 1.0, 1.0), 0.0, 1.0);
    let b = prism(rect(1.0, 1.0, 2.0, 2.0), 0.0, 2.0);
    let error = boolean_prisms_exact_solids(&a, &b, BooleanOperator::Union, tol())
        .expect_err("non-manifold");
    assert!(
        matches!(error, GeomError::UnsupportedInput { input, .. } if input.contains("touch along an edge")),
        "{error:?}"
    );
}

#[test]
fn equal_spans_still_take_the_single_prism_path() {
    // Unchanged behaviour: same spans, one prism, named after the subject.
    let a = prism(rect(0.0, 0.0, 2.0, 1.0), 0.0, 1.0);
    let b = prism(rect(1.0, 0.0, 3.0, 1.0), 0.0, 1.0);
    let solid = boolean_prisms_exact(&a, &b, BooleanOperator::Union, tol()).unwrap();
    close(planar_volume(&solid), 3.0, 1e-12);
    // Both caps are the subject's own (same planes, subject first).
    let names = names(&solid);
    assert!(names.contains(&FaceName::swept(SweptFace::StartCap).fragment(Operand::Subject)));
    assert!(names.contains(&FaceName::swept(SweptFace::EndCap).fragment(Operand::Subject)));
}

#[test]
fn blocks_meeting_along_a_horizontal_edge_are_refused_by_name() {
    // Side by side in plan, one below z = 1 and one above it: they share
    // only the edge x = 1, z = 1. Walls on the shared piece change sides
    // there, and the union is not a manifold solid.
    let low = prism(rect(0.0, 0.0, 1.0, 1.0), 0.0, 1.0);
    let high = prism(rect(1.0, 0.0, 2.0, 1.0), 1.0, 2.0);
    let error = boolean_prisms_exact_solids(&low, &high, BooleanOperator::Union, tol())
        .expect_err("non-manifold");
    assert!(
        matches!(error, GeomError::UnsupportedInput { input, .. } if input.contains("touch along an edge")),
        "{error:?}"
    );
}

#[test]
fn hole_walls_are_numbered_after_the_outer_ring() {
    // A 4 x 4 frame with a 2 x 2 opening, joined to a taller bar on its
    // right side. The opening's four walls are the subject's walls 4..8:
    // its own extrusion numbers the outer ring first, then the hole.
    let hole = {
        let mut ring = rect(1.0, 1.0, 3.0, 3.0);
        ring.reverse();
        ring
    };
    let frame = Prism {
        rings: vec![rect(0.0, 0.0, 4.0, 4.0), hole],
        bottom: 0.0,
        top: 1.0,
    };
    let bar = prism(rect(4.0, 0.0, 5.0, 4.0), 0.0, 2.0);
    let solid = boolean_prisms_exact(&frame, &bar, BooleanOperator::Union, tol())
        .expect("a frame and a bar touching along a face are one solid");
    audit(&solid);
    // 16 - 4 = 12 for the frame, 4 * 2 = 8 for the bar.
    close(planar_volume(&solid), 20.0, 1e-12);
    let names = names(&solid);
    for ordinal in 4..8 {
        let want = FaceName::swept(SweptFace::Side(ordinal)).fragment(Operand::Subject);
        assert!(names.contains(&want), "missing {want:?} in {names:?}");
    }
}
