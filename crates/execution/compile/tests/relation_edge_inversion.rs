//! Edges without pcurves on B-rep faces over surface relations (#454, ADR
//! 0110).
//!
//! Every face here has edges with their 3D curves only, inverted onto the
//! relation by the certified search along its swept parameter. Each case
//! reads the mesh back as `relation_faces.rs` does: a dense sample of the
//! exact surface within the reported bound, and closure, watertightness
//! and volume for a closed wall. Where the relation's derivative bounds are
//! an atomic surface's (a cylinder's), so is its bound. The adversarial
//! cases at the end each end in a named refusal or in a report whose bound
//! holds: never in a certified bound below the deviation measured.

mod common;

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use axiolid_contracts::GeomResult;
use axiolid_core::{Frame3, Interval, Point3, Scalar, Vec3};
use axiolid_curve::{BSplineCurve3, Curve3, KnotSpec};
use axiolid_mesh_compile::DeviationBound;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId, RuledSide, SurfaceRelation};
use axiolid_surface::{Cylinder, Surface, Torus};
use axiolid_topology::{
    BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};

use common::*;

/// The lens's report where the surface cannot be shown invertible.
const LENS: &str = "the surface is not certifiably invertible along an edge without a pcurve";

/// The refusals of the inverter, as their text begins or reads.
const AXIS: &str = "runs into the axis of its surface of revolution";
const APEX: &str = "runs into the apex of its ruled surface";
const COLLAPSED: &str = "runs into a collapsed ruling of its ruled surface";
const NOT_UNIQUE: &str = "has no unique parameters on its";
const OVER_BUDGET: &str = "certified inversion of an edge sample onto a surface relation";

/// An inverted edge's span term, `d + e` plus the lens, is sampled to the
/// budget, and the face adds its triangles' term on top: on an atomic
/// surface too the bound lands a little above the budget (a cylinder sheet
/// without pcurves reports 1.09e-3 at 1e-3, on its atomic cylinder as on
/// its revolution). The cases here check that the bound holds, and that it
/// stays within this many budgets.
const INVERTED_CEILING: Scalar = 2.0;

fn without_pcurves(sheet: Sheet<'_>) -> Sheet<'_> {
    Sheet {
        pcurves: [false; 4],
        ..sheet
    }
}

/// `check_sheet` without the budget, within [`INVERTED_CEILING`] of it;
/// the bound.
fn assert_inverted(name: &str, support: Support<'_>, sheet: &Sheet<'_>, family: &str) -> Scalar {
    let (_, _, report) = check_sheet(name, support, sheet, family, false);
    let bound = report.bound.unwrap();
    assert!(
        bound <= INVERTED_CEILING * BUDGET,
        "{name}: bound {bound:e} beyond {INVERTED_CEILING} budgets"
    );
    bound
}

/// `check_closed_wall` without the budget, within [`INVERTED_CEILING`] of
/// it; the bound.
fn assert_inverted_wall(name: &str, wall: &Wall<'_>, expect: &Solid3<'_>) -> Scalar {
    let bound = check_closed_wall(name, wall, expect, false);
    assert!(
        bound <= INVERTED_CEILING * BUDGET,
        "{name}: bound {bound:e} beyond {INVERTED_CEILING} budgets"
    );
    bound
}

/// The relation's bound is the atomic surface's, to rounding.
fn assert_same_bound(name: &str, relation: Scalar, atomic: Scalar) {
    assert!(
        (relation - atomic).abs() <= 1e-9 * atomic,
        "{name}: {relation:e} on the relation against {atomic:e} on the atomic surface"
    );
}

fn cylinder_control(b: &mut GeometryGraphBuilder) -> NodeId {
    atomic(
        b,
        Surface::Cylinder(Cylinder {
            frame: frame_at(0.0),
            radius: RADIUS,
        }),
    )
}

fn extrusion(b: &mut GeometryGraphBuilder, curve: Curve3, direction: Vec3) -> NodeId {
    let swept_curve = push_curve3(b, curve);
    b.push(GeometryNode::SurfaceRelation(
        SurfaceRelation::LinearExtrusion {
            swept_curve,
            direction,
        },
    ))
    .unwrap()
}

/// The cylinder patch over `u` with its arcs and straight rulings.
fn cylinder_patch<'a>(u: (Scalar, Scalar)) -> Sheet<'a> {
    Sheet {
        u,
        ..cylinder_sheet(false, None)
    }
}

// --- linear extrusion ---------------------------------------------------------------

#[test]
fn an_extruded_cylinder_patch_without_pcurves_meshes_as_the_atomic_one() {
    let support = |b: &mut GeometryGraphBuilder| extruded_circle(b, Vec3::Z);
    // The second runs across the circle's seam at u = 0, which the cover
    // joins.
    for (name, u) in [("patch", (0.3, 2.0)), ("patch across u = 0", (-0.8, 0.9))] {
        let sheet = cylinder_patch(u);
        let relation = assert_inverted(&format!("extruded {name}"), &support, &sheet, EXTRUSION);
        let atomic = assert_inverted(
            &format!("cylinder {name}"),
            &cylinder_control,
            &sheet,
            "cylinder",
        );
        assert_same_bound(name, relation, atomic);
    }
}

#[test]
fn an_extruded_cylinder_wall_without_pcurves_is_a_watertight_solid() {
    let support = |b: &mut GeometryGraphBuilder| extruded_circle(b, Vec3::Z);
    let wall = |support| Wall {
        pcurves: false,
        ..cylinder_wall(support)
    };
    let relation = assert_inverted_wall(
        "extruded wall",
        &wall(&support),
        &Solid3 {
            family: EXTRUSION,
            ..cylinder_solid()
        },
    );
    let atomic = assert_inverted_wall(
        "cylinder wall",
        &wall(&cylinder_control),
        &Solid3 {
            family: "cylinder",
            ..cylinder_solid()
        },
    );
    assert_same_bound("extruded wall", relation, atomic);
}

/// A wall of one face bounded by its two rims, with no seam edge.
fn ring_wall(support: Support<'_>) -> GeomResult<Compiled> {
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let support = support(&mut builder);
    let low = brep.add_vertex(Vertex {
        position: Point3::new(RADIUS, 0.0, 0.0),
    });
    let high = brep.add_vertex(Vertex {
        position: Point3::new(RADIUS, 0.0, HEIGHT),
    });
    let mut rim = |vertex, z| {
        let curve = push_curve3(&mut builder, horizontal_circle(0.0, 0.0, z, RADIUS));
        brep.add_edge(Edge {
            start: vertex,
            end: vertex,
            curve: Some(curve),
        })
    };
    let (bottom, top) = (rim(low, 0.0), rim(high, HEIGHT));
    let ring = |edge, orientation| Loop {
        edges: vec![EdgeUse {
            edge,
            orientation,
            pcurve: None,
        }],
    };
    let bottom = brep.add_loop(ring(bottom, Orientation::Forward));
    let top = brep.add_loop(ring(top, Orientation::Reversed));
    let face = brep.add_face(Face {
        surface: Some(support),
        bounds: vec![
            FaceBound {
                loop_id: bottom,
                orientation: Orientation::Forward,
                outer: true,
            },
            FaceBound {
                loop_id: top,
                orientation: Orientation::Forward,
                outer: false,
            },
        ],
        orientation: Orientation::Forward,
    });
    brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    compile(builder, brep)
}

/// A full turn bounded by two rims and no seam edge is not triangulated
/// on an atomic cylinder either: its rings do not close in the parameter
/// plane. The relation comes to the same refusal, not to another one.
#[test]
fn a_full_turn_wall_bounded_by_its_two_rims_is_refused_as_the_atomic_one_is() {
    let support = |b: &mut GeometryGraphBuilder| extruded_circle(b, Vec3::Z);
    let refusal = |support| match ring_wall(support) {
        Ok((_, _, report)) => panic!("expected a refusal, got {report:?}"),
        Err(error) => error.to_string(),
    };
    let (relation, atomic) = (refusal(&support), refusal(&cylinder_control));
    assert!(atomic.contains("cannot be triangulated"), "{atomic}");
    assert_eq!(relation, atomic);
}

#[test]
fn an_obliquely_extruded_circle_without_pcurves_meshes_within_its_bound() {
    let direction = Vec3::new(0.4, 0.0, 2.0);
    let support = |b: &mut GeometryGraphBuilder| extruded_circle(b, direction);
    let point = |u: Scalar, v: Scalar| cylinder_sheet_point(u, 0.0) + direction * v;
    let sheet = Sheet {
        u: (0.2, 2.9),
        v: (0.0, 0.25),
        point: &point,
        curves: [
            Some(horizontal_circle(0.0, 0.0, 0.0, RADIUS)),
            None,
            Some(horizontal_circle(0.1, 0.0, HEIGHT, RADIUS)),
            None,
        ],
        pcurves: [false; 4],
        hole: None,
        window: None,
    };
    assert_inverted("oblique extrusion", &support, &sheet, EXTRUSION);
}

/// The cubic Bezier profile at `u`, by de Casteljau.
fn casteljau(points: &[Point3], u: Scalar) -> Point3 {
    let mut p = points.to_vec();
    while p.len() > 1 {
        p = p.windows(2).map(|w| w[0].lerp(w[1], u)).collect();
    }
    p[0]
}

/// The control points of a cubic Bezier's piece over `a..b`.
fn bezier_piece(points: [Point3; 4], a: Scalar, b: Scalar) -> [Point3; 4] {
    let split = |p: [Point3; 4], t: Scalar| {
        let l1 = [p[0].lerp(p[1], t), p[1].lerp(p[2], t), p[2].lerp(p[3], t)];
        let l2 = [l1[0].lerp(l1[1], t), l1[1].lerp(l1[2], t)];
        let m = l2[0].lerp(l2[1], t);
        ([p[0], l1[0], l2[0], m], [m, l2[1], l1[2], p[3]])
    };
    split(split(points, b).0, a / b).1
}

#[test]
fn an_extruded_cubic_profile_without_pcurves_meshes_within_its_bound() {
    let direction = Vec3::new(0.1, 0.8, 0.0);
    let support = |b: &mut GeometryGraphBuilder| extrusion(b, bezier(bezier_points()), direction);
    let point = |u: Scalar, v: Scalar| casteljau(&bezier_points(), u) + direction * v;
    let lifted = bezier_points().map(|p| p + direction * 0.5);
    let sheet = Sheet {
        u: (0.0, 1.0),
        v: (0.0, 0.5),
        point: &point,
        curves: [
            Some(bezier(bezier_points())),
            None,
            Some(bezier(lifted)),
            None,
        ],
        pcurves: [false; 4],
        hole: None,
        window: None,
    };
    assert_inverted("extruded cubic", &support, &sheet, EXTRUSION);
}

// --- revolution ---------------------------------------------------------------------

fn frustum_revolution(b: &mut GeometryGraphBuilder) -> NodeId {
    revolution(
        b,
        line3(Point3::new(0.2, 0.0, 0.0), Vec3::new(0.2, 0.0, 1.0)),
        Point3::ZERO,
        Vec3::Z,
    )
}

fn frustum_sheet_point(u: Scalar, v: Scalar) -> Point3 {
    let r = 0.2 + 0.2 * v;
    Point3::new(r * u.cos(), r * u.sin(), v)
}

fn frustum_patch<'a>(u: (Scalar, Scalar)) -> Sheet<'a> {
    Sheet {
        u,
        v: (0.0, HEIGHT),
        point: &frustum_sheet_point,
        curves: [
            Some(horizontal_circle(0.0, 0.0, 0.0, 0.2)),
            None,
            Some(horizontal_circle(0.0, 0.0, HEIGHT, 0.3)),
            None,
        ],
        pcurves: [false; 4],
        hole: None,
        window: None,
    }
}

#[test]
fn a_revolved_slanted_line_without_pcurves_meshes_as_its_frustum() {
    let family = "surface of revolution";
    assert_inverted(
        "frustum patch",
        &frustum_revolution,
        &frustum_patch((0.4, 2.2)),
        family,
    );
    // Across u = pi, where the angle read back wraps.
    assert_inverted(
        "frustum patch across u = pi",
        &frustum_revolution,
        &frustum_patch((2.5, 3.9)),
        family,
    );
    let wall = Wall {
        pcurves: false,
        ..frustum_wall(&frustum_revolution)
    };
    assert_inverted_wall("frustum wall", &wall, &frustum_solid());
}

/// The revolved wall shares each rim with its planar cap: neither use has
/// a pcurve, both invert the rim's one set of samples.
#[test]
fn a_revolved_line_wall_and_its_caps_share_their_rims_without_pcurves() {
    let wall = |support| Wall {
        pcurves: false,
        ..cylinder_wall(support)
    };
    let relation = assert_inverted_wall(
        "revolved line wall",
        &wall(&cylinder_revolution),
        &cylinder_solid(),
    );
    let atomic = assert_inverted_wall(
        "cylinder wall",
        &wall(&cylinder_control),
        &Solid3 {
            family: "cylinder",
            ..cylinder_solid()
        },
    );
    assert_same_bound("revolved line wall", relation, atomic);
}

fn meridian(u: Scalar) -> Curve3 {
    circle_in(
        Point3::new(1.0 + 0.4 * u.cos(), 2.0 + 0.4 * u.sin(), 0.1),
        Vec3::new(u.cos(), u.sin(), 0.0),
        Vec3::Z,
        0.1,
    )
}

/// The torus patch over `u` by `v`, its rims at constant `v` and its
/// meridians at constant `u`.
fn torus_sheet<'a>(u: (Scalar, Scalar), v: (Scalar, Scalar)) -> Sheet<'a> {
    let rim = |v: Scalar| horizontal_circle(1.0, 2.0, 0.1 + 0.1 * v.sin(), 0.4 + 0.1 * v.cos());
    Sheet {
        u,
        v,
        point: &torus_point,
        curves: [
            Some(rim(v.0)),
            Some(meridian(u.1)),
            Some(rim(v.1)),
            Some(meridian(u.0)),
        ],
        pcurves: [false; 4],
        hole: None,
        window: None,
    }
}

fn torus_revolution(b: &mut GeometryGraphBuilder) -> NodeId {
    revolution(b, torus_profile(), Point3::new(1.0, 2.0, 0.0), Vec3::Z)
}

fn torus_control(b: &mut GeometryGraphBuilder) -> NodeId {
    atomic(
        b,
        Surface::Torus(Torus {
            frame: Frame3 {
                origin: Point3::new(1.0, 2.0, 0.1),
                ..frame_at(0.0)
            },
            major_radius: 0.4,
            minor_radius: 0.1,
        }),
    )
}

#[test]
fn a_revolved_circle_patch_without_pcurves_meshes_as_its_torus() {
    let family = "surface of revolution";
    // The second runs across both seams: the turn's at u = pi and the
    // profile's at v = 0.
    for (name, sheet) in [
        ("torus patch", without_pcurves(torus_patch(&torus_point))),
        (
            "torus patch across both seams",
            torus_sheet((2.6, 3.7), (-0.6, 0.5)),
        ),
    ] {
        // The relation's derivative bounds are not the torus's own, so
        // neither is its bound; the atomic control checks the fixture.
        assert_inverted(name, &torus_revolution, &sheet, family);
        assert_inverted(name, &torus_control, &sheet, "torus");
    }
}

#[test]
fn a_revolved_cubic_profile_without_pcurves_meshes_within_its_bound() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(b, bezier(bezier_points()), Point3::ZERO, Vec3::Z)
    };
    let turned = bezier_points().map(|p| Point3::new(-p.x, -p.y, p.z));
    let sheet = Sheet {
        u: (0.0, PI),
        v: (0.0, 1.0),
        point: &bezier_point,
        curves: [
            Some(horizontal_circle(0.0, 0.0, 0.0, 0.3)),
            Some(bezier(turned)),
            Some(horizontal_circle(0.0, 0.0, 0.3, 0.3)),
            Some(bezier(bezier_points())),
        ],
        pcurves: [false; 4],
        hole: None,
        window: None,
    };
    assert_inverted("revolved cubic", &support, &sheet, "surface of revolution");
}

#[test]
fn a_revolved_sheet_with_a_window_without_pcurves_meshes_round_it() {
    let sheet = Sheet {
        window: Some(Window {
            u: (1.2, 1.9),
            v: (0.15, 0.35),
            curves: [
                Some(horizontal_circle(0.0, 0.0, 0.15, RADIUS)),
                None,
                Some(horizontal_circle(0.0, 0.0, 0.35, RADIUS)),
                None,
            ],
        }),
        ..cylinder_sheet(false, None)
    };
    let (mesh, _, _) = check_sheet(
        "revolved sheet with a window",
        &cylinder_revolution,
        &sheet,
        "surface of revolution",
        false,
    );
    for t in mesh.indices.as_chunks::<3>().0 {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        let m = (a + b + c) / 3.0;
        let (u, z) = (m.y.atan2(m.x), m.z);
        let inside = u > 1.2 + 0.01 && u < 1.9 - 0.01 && z > 0.15 + 0.01 && z < 0.35 - 0.01;
        assert!(!inside, "a triangle centred in the window at ({u}, {z})");
    }
}

#[test]
fn a_ring_mixing_uses_with_and_without_pcurves_meshes_within_its_bound() {
    for pcurves in [[true, false, true, false], [false, true, false, true]] {
        let sheet = Sheet {
            pcurves,
            ..cylinder_sheet(true, None)
        };
        assert_inverted(
            &format!("mixed ring {pcurves:?}"),
            &cylinder_revolution,
            &sheet,
            "surface of revolution",
        );
    }
}

#[test]
fn the_same_face_with_and_without_pcurves_is_certified_both_ways() {
    let family = "surface of revolution";
    assert_sheet(
        "with pcurves",
        &cylinder_revolution,
        &cylinder_sheet(true, None),
        family,
    );
    let sheet = cylinder_sheet(false, None);
    let relation = assert_inverted("without pcurves", &cylinder_revolution, &sheet, family);
    let atomic = assert_inverted("atomic without", &cylinder_control, &sheet, "cylinder");
    assert_same_bound("without pcurves", relation, atomic);
}

#[test]
fn the_same_face_meshes_the_same_on_every_call() {
    let sheet = without_pcurves(torus_patch(&torus_point));
    let (first, _, first_report) = open_sheet(&torus_revolution, &sheet).unwrap();
    for _ in 0..2 {
        let (again, _, report) = open_sheet(&torus_revolution, &sheet).unwrap();
        assert_eq!(first.positions, again.positions);
        assert_eq!(first.indices, again.indices);
        assert_eq!(first_report.bound, report.bound);
    }
}

// --- ruled surfaces -------------------------------------------------------------------

fn saddle_point(u: Scalar, v: Scalar) -> Point3 {
    Point3::new(u * (1.0 - v), v, u * v)
}

#[test]
fn a_ruled_saddle_with_a_window_without_pcurves_meshes_within_its_bound() {
    let support = |b: &mut GeometryGraphBuilder| {
        ruled(
            b,
            [
                (line3(Point3::ZERO, Vec3::X), 0.0, 1.0),
                (line3(Point3::new(0.0, 1.0, 0.0), Vec3::Z), 0.0, 1.0),
            ],
        )
    };
    let sheet = Sheet {
        u: (0.0, 1.0),
        v: (0.0, 1.0),
        point: &saddle_point,
        curves: [None, None, None, None],
        pcurves: [false; 4],
        hole: None,
        window: Some(Window {
            u: (0.3, 0.6),
            v: (0.3, 0.6),
            curves: [None, None, None, None],
        }),
    };
    assert_inverted("ruled saddle", &support, &sheet, RULED);
}

#[test]
fn an_arc_ruled_against_a_line_without_pcurves_meshes_within_its_bound() {
    let support = |b: &mut GeometryGraphBuilder| arc_line_support(b, false);
    let sheet = without_pcurves(unit_sheet(
        &arc_line_point,
        [
            Some(horizontal_circle(0.0, 0.0, 0.0, 1.0)),
            None,
            Some(line3(Point3::new(1.0, 0.0, 1.0), Vec3::new(-1.0, 1.0, 0.0))),
            None,
        ],
    ));
    assert_inverted("arc against a line", &support, &sheet, RULED);
}

// --- adversarial ----------------------------------------------------------------------

/// What an adversarial case may come to.
enum Outcome {
    /// A certified bound, which must hold.
    Certified,
    /// A report unbounded for this reason.
    Unbounded(&'static str),
    /// A refusal whose text holds this.
    Refused(&'static str),
}

/// The case comes to one of `allowed`; whatever it reports, a certified
/// bound is never below the deviation measured at `samples`.
fn honest(name: &str, result: GeomResult<Compiled>, samples: &[Point3], allowed: &[Outcome]) {
    match result {
        Err(error) => {
            let text = error.to_string();
            eprintln!("{name}: refused, {text}");
            assert!(
                allowed
                    .iter()
                    .any(|o| matches!(o, Outcome::Refused(s) if text.contains(s))),
                "{name}: refused by an outcome not listed: {error:?}"
            );
        }
        Ok((mesh, _, report)) => match report.bound {
            Some(bound) => {
                let worst = farthest(&mesh, samples);
                eprintln!("{name}: certified {bound:e}, measured {worst:e}");
                // The measurement's own rounding: on a plane patch
                // certified exact, `closest_point_on_triangle` lands 1.8e-13
                // off a sample inside a long thin triangle.
                let scale = samples
                    .iter()
                    .map(|p| p.abs().max_element())
                    .fold(1.0, Scalar::max);
                assert!(
                    worst <= bound * (1.0 + 1e-9) + 1e-12 * scale,
                    "{name}: certified {bound:e} below the measured {worst:e}"
                );
                assert!(
                    allowed.iter().any(|o| matches!(o, Outcome::Certified)),
                    "{name}: certified {bound:e}, an outcome not listed"
                );
            }
            None => {
                let reasons: Vec<&str> = report
                    .contributions
                    .iter()
                    .filter_map(|c| match c.bound {
                        DeviationBound::Unbounded(reason) => Some(reason),
                        _ => None,
                    })
                    .collect();
                assert!(
                    allowed.iter().any(|o| matches!(
                        o,
                        Outcome::Unbounded(s) if reasons.iter().any(|r| r.contains(s))
                    )),
                    "{name}: unbounded for a reason not listed: {reasons:?}"
                );
            }
        },
    }
}

fn honest_sheet(name: &str, support: Support<'_>, sheet: &Sheet<'_>, allowed: &[Outcome]) {
    honest(
        name,
        open_sheet(support, sheet),
        &sheet_samples(sheet),
        allowed,
    );
}

#[test]
fn a_frustum_ten_tolerances_wide_at_its_end_is_certified_or_unbounded_by_name() {
    let tau = BUDGET;
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(
            b,
            line3(Point3::new(10.0 * tau, 0.0, 0.0), Vec3::new(0.2, 0.0, 1.0)),
            Point3::ZERO,
            Vec3::Z,
        )
    };
    let point = |u: Scalar, v: Scalar| {
        let r = 10.0 * tau + 0.2 * v;
        Point3::new(r * u.cos(), r * u.sin(), v)
    };
    let sheet = Sheet {
        point: &point,
        curves: [
            Some(horizontal_circle(0.0, 0.0, 0.0, 10.0 * tau)),
            None,
            Some(horizontal_circle(0.0, 0.0, HEIGHT, 10.0 * tau + 0.1)),
            None,
        ],
        ..frustum_patch((0.4, 2.2))
    };
    honest_sheet(
        "narrow frustum",
        &support,
        &sheet,
        &[Outcome::Certified, Outcome::Unbounded(LENS)],
    );
}

#[test]
fn an_edge_through_the_apex_of_a_revolved_line_is_refused_at_the_axis() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(
            b,
            line3(Point3::ZERO, Vec3::new(0.2, 0.0, 1.0)),
            Point3::ZERO,
            Vec3::Z,
        )
    };
    let point = |u: Scalar, v: Scalar| Point3::new(0.2 * v * u.cos(), 0.2 * v * u.sin(), v);
    let sheet = Sheet {
        u: (0.4, 2.2),
        v: (-0.5, 0.5),
        point: &point,
        curves: [
            Some(circle_in(
                Point3::new(0.0, 0.0, -0.5),
                -Vec3::X,
                -Vec3::Y,
                0.1,
            )),
            None,
            Some(horizontal_circle(0.0, 0.0, 0.5, 0.1)),
            None,
        ],
        pcurves: [false; 4],
        hole: None,
        window: None,
    };
    honest_sheet(
        "through the apex",
        &support,
        &sheet,
        &[Outcome::Refused(AXIS)],
    );
}

/// `c(s) = (u^3, 0, u)` with `u = 2s - 1`, scaled by 0.3: its tangent runs
/// along `z` at `s = 1/2`.
fn cubic_fold() -> [Point3; 4] {
    [
        (-1.0, -1.0),
        (1.0, -1.0 / 3.0),
        (-1.0, 1.0 / 3.0),
        (1.0, 1.0),
    ]
    .map(|(x, z)| Point3::new(0.3 * x, 0.0, 0.3 * z))
}

#[test]
fn a_profile_tangent_to_the_direction_is_unbounded_or_refused_by_name() {
    let support = |b: &mut GeometryGraphBuilder| extrusion(b, bezier(cubic_fold()), Vec3::Z);
    let point = |u: Scalar, v: Scalar| casteljau(&cubic_fold(), u) + Vec3::Z * v;
    let (a, b) = (0.3, 0.7);
    let low = bezier_piece(cubic_fold(), a, b);
    let sheet = Sheet {
        u: (a, b),
        v: (0.0, HEIGHT),
        point: &point,
        curves: [
            Some(bezier(low)),
            None,
            Some(bezier(low.map(|p| p + Vec3::Z * HEIGHT))),
            None,
        ],
        pcurves: [false; 4],
        hole: None,
        window: None,
    };
    honest_sheet(
        "tangent to the direction",
        &support,
        &sheet,
        &[
            Outcome::Unbounded(LENS),
            Outcome::Refused(NOT_UNIQUE),
            Outcome::Refused(OVER_BUDGET),
        ],
    );
}

#[test]
fn a_profile_extruded_in_its_own_plane_is_refused_as_not_unique() {
    let support = |b: &mut GeometryGraphBuilder| {
        extrusion(b, circle_in(Point3::ZERO, Vec3::X, Vec3::Z, 0.3), Vec3::X)
    };
    let point = |u: Scalar, v: Scalar| Point3::new(0.3 * u.cos() + v, 0.0, 0.3 * u.sin());
    let sheet = Sheet {
        u: (0.3, 1.2),
        v: (0.0, HEIGHT),
        point: &point,
        curves: [
            Some(circle_in(Point3::ZERO, Vec3::X, Vec3::Z, 0.3)),
            None,
            Some(circle_in(
                Point3::new(HEIGHT, 0.0, 0.0),
                Vec3::X,
                Vec3::Z,
                0.3,
            )),
            None,
        ],
        pcurves: [false; 4],
        hole: None,
        window: None,
    };
    honest_sheet(
        "in its own plane",
        &support,
        &sheet,
        &[Outcome::Refused(NOT_UNIQUE)],
    );
}

#[test]
fn a_full_circle_revolved_about_its_diameter_is_refused_as_not_unique() {
    let support =
        |b: &mut GeometryGraphBuilder| revolution(b, zone_circle(), Point3::ZERO, Vec3::Z);
    let wall = Wall {
        pcurves: false,
        ..zone_wall(&support)
    };
    honest(
        "sphere covered twice",
        closed_wall(&wall),
        &[],
        &[Outcome::Refused(NOT_UNIQUE), Outcome::Refused(AXIS)],
    );
}

#[test]
fn a_line_perpendicular_to_the_axis_is_refused_as_not_unique() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(
            b,
            line3(Point3::new(0.0, 0.0, 0.2), Vec3::X),
            Point3::ZERO,
            Vec3::Z,
        )
    };
    let point = |u: Scalar, v: Scalar| Point3::new(v * u.cos(), v * u.sin(), 0.2);
    let sheet = Sheet {
        u: (0.3, 1.2),
        v: (0.1, 0.3),
        point: &point,
        curves: [
            Some(horizontal_circle(0.0, 0.0, 0.2, 0.1)),
            None,
            Some(horizontal_circle(0.0, 0.0, 0.2, 0.3)),
            None,
        ],
        pcurves: [false; 4],
        hole: None,
        window: None,
    };
    honest_sheet(
        "plane covered twice",
        &support,
        &sheet,
        &[Outcome::Refused(NOT_UNIQUE)],
    );
}

/// A cubic Bezier in the `xy` plane crossing itself on `x = 0`, symmetric
/// under `x -> -x` with `s -> 1 - s`.
fn figure_eight() -> [Point3; 4] {
    [(-1.0, 0.0), (2.0, 1.5), (-2.0, 1.5), (1.0, 0.0)].map(|(x, y)| Point3::new(x, y, 0.0))
}

#[test]
fn a_self_crossing_profile_is_refused_on_its_crossing_line() {
    // The crossing is at `s` and `1 - s` with `x(s) = 0`, `s` in (0, 1/4).
    let (mut lo, mut hi) = (0.0, 0.25);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if casteljau(&figure_eight(), mid).x < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let s = 0.5 * (lo + hi);
    let support = |b: &mut GeometryGraphBuilder| extrusion(b, bezier(figure_eight()), Vec3::Z);
    let point = |u: Scalar, v: Scalar| casteljau(&figure_eight(), u) + Vec3::Z * v;
    let low = bezier_piece(figure_eight(), s, s + 0.15);
    let sheet = Sheet {
        u: (s, s + 0.15),
        v: (0.0, HEIGHT),
        point: &point,
        curves: [
            Some(bezier(low)),
            None,
            Some(bezier(low.map(|p| p + Vec3::Z * HEIGHT))),
            None,
        ],
        pcurves: [false; 4],
        hole: None,
        window: None,
    };
    honest_sheet(
        "figure eight",
        &support,
        &sheet,
        &[Outcome::Refused(NOT_UNIQUE)],
    );
}

#[test]
fn an_edge_ten_tolerances_off_its_surface_is_refused_with_a_floor() {
    let support = |b: &mut GeometryGraphBuilder| extruded_circle(b, Vec3::Z);
    let off = RADIUS + 10.0 * BUDGET;
    let point = |u: Scalar, v: Scalar| {
        let r = if v == HEIGHT { off } else { RADIUS };
        Point3::new(r * u.cos(), r * u.sin(), v)
    };
    let sheet = Sheet {
        point: &point,
        curves: [
            Some(horizontal_circle(0.0, 0.0, 0.0, RADIUS)),
            None,
            Some(horizontal_circle(0.0, 0.0, HEIGHT, off)),
            None,
        ],
        ..cylinder_patch((0.3, 2.0))
    };
    let text = refused(open_sheet(&support, &sheet)).to_string();
    let floor: Scalar = text
        .split("leaves its surface of linear extrusion by at least ")
        .nth(1)
        .and_then(|rest| rest.split(':').next())
        .unwrap_or_else(|| panic!("{text}"))
        .parse()
        .unwrap();
    assert!(floor > BUDGET && floor <= 10.0 * BUDGET, "{floor:e}");
}

/// A closed cubic B-spline circle of `RADIUS` in the `xy` plane, its seam at
/// `(RADIUS, 0, 0)` and its knots at the quarters.
fn closed_spline_circle() -> Curve3 {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let p = |x: Scalar, y: Scalar| Point3::new(x * RADIUS, y * RADIUS, 0.0);
    Curve3::BSpline(BSplineCurve3 {
        degree: 2,
        control_points: vec![
            p(1.0, 0.0),
            p(1.0, 1.0),
            p(0.0, 1.0),
            p(-1.0, 1.0),
            p(-1.0, 0.0),
            p(-1.0, -1.0),
            p(0.0, -1.0),
            p(1.0, -1.0),
            p(1.0, 0.0),
        ],
        knots: vec![0.0, 0.25, 0.5, 0.75, 1.0],
        multiplicities: vec![3, 2, 2, 2, 3],
        weights: Some(vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0]),
        closed: true,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

#[test]
fn a_seam_edge_of_a_closed_spline_profile_is_refused_as_not_unique() {
    let support = |b: &mut GeometryGraphBuilder| extrusion(b, closed_spline_circle(), Vec3::Z);
    // Exact at the knots 0 and 1/4, the patch's corners.
    let point = |u: Scalar, v: Scalar| cylinder_sheet_point(TAU * u, v);
    let sheet = Sheet {
        u: (0.0, 0.25),
        point: &point,
        ..cylinder_patch((0.0, FRAC_PI_2))
    };
    honest_sheet(
        "closed spline seam",
        &support,
        &sheet,
        &[Outcome::Refused(NOT_UNIQUE)],
    );
}

/// One face over `support` bounded by `corners`, edge `i` running from
/// corner `i` to the next with 3D curve `curves[i]`, no pcurves.
fn polygon_face(
    support: Support<'_>,
    corners: &[Point3],
    curves: Vec<Option<Curve3>>,
) -> GeomResult<Compiled> {
    let mut builder = GeometryGraphBuilder::new();
    let mut brep: BRep<NodeId> = BRep::default();
    let support = support(&mut builder);
    let vertices: Vec<_> = corners
        .iter()
        .map(|&position| brep.add_vertex(Vertex { position }))
        .collect();
    let edges = curves
        .into_iter()
        .enumerate()
        .map(|(i, curve)| {
            let curve = curve.map(|c| push_curve3(&mut builder, c));
            let edge = brep.add_edge(Edge {
                start: vertices[i],
                end: vertices[(i + 1) % vertices.len()],
                curve,
            });
            EdgeUse {
                edge,
                orientation: Orientation::Forward,
                pcurve: None,
            }
        })
        .collect();
    let loop_id = brep.add_loop(Loop { edges });
    let face = brep.add_face(Face {
        surface: Some(support),
        bounds: vec![FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    compile(builder, brep)
}

#[test]
fn an_edge_into_a_ruled_apex_is_refused_at_the_apex() {
    let apex = Point3::new(0.0, 0.0, 1.0);
    let support = |b: &mut GeometryGraphBuilder| {
        let curve = push_curve3(b, horizontal_circle(0.0, 0.0, 0.0, 1.0));
        b.push(GeometryNode::SurfaceRelation(SurfaceRelation::Ruled {
            first: RuledSide::Curve {
                curve,
                interval: Interval {
                    start: 0.0,
                    end: FRAC_PI_2,
                },
            },
            second: RuledSide::Apex(apex),
        }))
        .unwrap()
    };
    let rim = |angle: Scalar| Point3::new(angle.cos(), angle.sin(), 0.0);
    let corners = [rim(0.2), rim(1.2), apex];
    let result = polygon_face(
        &support,
        &corners,
        vec![Some(horizontal_circle(0.0, 0.0, 0.0, 1.0)), None, None],
    );
    honest("into the apex", result, &[], &[Outcome::Refused(APEX)]);
}

#[test]
fn an_edge_into_a_collapsed_ruling_is_refused_there() {
    let support = |b: &mut GeometryGraphBuilder| {
        ruled(
            b,
            [
                (line3(Point3::ZERO, Vec3::X), 0.0, 1.0),
                (
                    line3(Point3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, -1.0)),
                    0.0,
                    1.0,
                ),
            ],
        )
    };
    let corners = [Point3::ZERO, Point3::X, Point3::Z];
    let result = polygon_face(&support, &corners, vec![None, None, None]);
    honest(
        "collapsed ruling",
        result,
        &[],
        &[Outcome::Refused(COLLAPSED)],
    );
}

/// An edge crossing the seam of a closed swept spline off a sample: the
/// spline has no period, so its samples' parameters jump by the domain
/// there and the edge is refused, not read by continuity.
#[test]
fn an_edge_across_a_closed_spline_s_seam_off_a_sample_is_refused_as_not_unique() {
    let support = |b: &mut GeometryGraphBuilder| extrusion(b, closed_spline_circle(), Vec3::Z);
    let sheet = cylinder_patch((-0.4, 0.7));
    honest_sheet(
        "across the closed spline's seam",
        &support,
        &sheet,
        &[Outcome::Refused(NOT_UNIQUE)],
    );
}

fn flat_ellipse(z: Scalar, a: Scalar, b: Scalar) -> Curve3 {
    Curve3::Ellipse(axiolid_curve::Ellipse3 {
        frame: frame_at(z),
        semi_axis_x: a,
        semi_axis_y: b,
    })
}

/// Over a flat ellipse's vertex the generators move at the short
/// semi-axis's speed, not the long one's: edges a fraction of the tolerance
/// off the surface are no fold, and the face is certified against the exact
/// surface.
#[test]
fn a_flat_ellipse_extruded_over_its_slow_vertex_is_no_fold() {
    for (a, b, off) in [(3.0, 0.1, 0.1), (1.2, 0.1, 0.7)] {
        let off = off * BUDGET;
        let support =
            move |bld: &mut GeometryGraphBuilder| extrusion(bld, flat_ellipse(0.0, a, b), Vec3::Z);
        let given =
            move |u: Scalar, v: Scalar| Point3::new((a + off) * u.cos(), (b + off) * u.sin(), v);
        let exact = move |u: Scalar, v: Scalar| Point3::new(a * u.cos(), b * u.sin(), v);
        let sheet = Sheet {
            point: &given,
            curves: [
                Some(flat_ellipse(0.0, a + off, b + off)),
                None,
                Some(flat_ellipse(HEIGHT, a + off, b + off)),
                None,
            ],
            ..cylinder_patch((-0.3, 0.3))
        };
        let inside = Sheet {
            point: &exact,
            ..cylinder_patch((-0.29, 0.29))
        };
        honest(
            &format!("ellipse {a} x {b}, edges {off:e} off"),
            open_sheet(&support, &sheet),
            &sheet_samples(&inside),
            &[Outcome::Certified],
        );
    }
}

/// Quarter arcs of radii 1 and 0.03 ruled a unit apart: near the small rim
/// the ruling pivots, and the edges' uniform samples must be as close as
/// the lens needs there all along. The face is certified, at a cost the
/// same face with pcurves does not pay (ADR 0110); this pins that cost
/// until edges are sampled by their local speed.
#[test]
fn a_ruling_pivoting_about_a_slow_rim_meshes_at_the_cost_of_uniform_samples() {
    let (big, small) = (1.0, 0.03);
    let support = move |bld: &mut GeometryGraphBuilder| {
        ruled(
            bld,
            [
                (horizontal_circle(0.0, 0.0, 0.0, big), 0.0, FRAC_PI_2),
                (horizontal_circle(0.0, 0.0, 1.0, small), 0.0, FRAC_PI_2),
            ],
        )
    };
    let exact = move |u: Scalar, v: Scalar| {
        let angle = FRAC_PI_2 * u;
        Point3::new(angle.cos(), angle.sin(), 0.0) * (big * (1.0 - v) + small * v)
            + Point3::new(0.0, 0.0, v)
    };
    let curves = || {
        [
            Some(horizontal_circle(0.0, 0.0, 0.0, big)),
            None,
            Some(horizontal_circle(0.0, 0.0, 1.0, small)),
            None,
        ]
    };
    let with = check_sheet(
        "pivoting ruling, pcurves",
        &support,
        &unit_sheet(&exact, curves()),
        RULED,
        true,
    );
    let without = without_pcurves(unit_sheet(&exact, curves()));
    let (mesh, _, report) = check_sheet("pivoting ruling", &support, &without, RULED, false);
    let bound = report.bound.unwrap();
    assert!(
        bound <= INVERTED_CEILING * BUDGET,
        "{bound:e} beyond {INVERTED_CEILING} budgets"
    );
    let (cheap, dear) = (with.0.indices.len() / 3, mesh.indices.len() / 3);
    assert!(
        dear > 10 * cheap,
        "{dear} triangles without pcurves against {cheap} with them: update ADR 0110"
    );
}

/// Two straight pieces meeting at `knot` with a corner: `(0, 0)` to `(3, 0)`
/// at speed `3 / knot`, then `(3, 0)` to `(3, 0.07)` a hundred times
/// slower.
fn kinked(degree: u16, knot: Scalar) -> Curve3 {
    let p = |x: Scalar, y: Scalar| Point3::new(x, y, 0.0);
    let (control_points, multiplicities) = match degree {
        1 => (vec![p(0.0, 0.0), p(3.0, 0.0), p(3.0, 0.07)], vec![2, 1, 2]),
        _ => (
            vec![
                p(0.0, 0.0),
                p(1.5, 0.0),
                p(3.0, 0.0),
                p(3.0, 0.035),
                p(3.0, 0.07),
            ],
            vec![3, 2, 3],
        ),
    };
    Curve3::BSpline(BSplineCurve3 {
        degree,
        control_points,
        knots: vec![0.0, knot, 1.0],
        multiplicities,
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// A plane patch on the fast side of a profile's corner, its edge a
/// hundredth of the corner's piece from it: no cell of the cover spans the
/// corner, where the speed jumps, so the slow side's speed never stands for
/// the fast side's and the minimiser is never dropped.
#[test]
fn a_patch_beside_a_profile_s_corner_meshes_without_pcurves() {
    for degree in [1, 2] {
        for knot in [0.3, 0.5] {
            let speed = 3.0 / knot;
            let support =
                move |b: &mut GeometryGraphBuilder| extrusion(b, kinked(degree, knot), Vec3::Z);
            let point = move |u: Scalar, v: Scalar| Point3::new(speed * u, 0.0, v);
            for (x0, x1) in [(1.0, 2.0), (2.0, 2.985)] {
                let sheet = Sheet {
                    u: (x0 / speed, x1 / speed),
                    v: (0.0, HEIGHT),
                    point: &point,
                    curves: [None, None, None, None],
                    pcurves: [false; 4],
                    hole: None,
                    window: None,
                };
                honest(
                    &format!("degree {degree}, corner at {knot}, x in [{x0}, {x1}]"),
                    open_sheet(&support, &sheet),
                    &sheet_samples(&sheet),
                    &[Outcome::Certified],
                );
            }
        }
    }
}
