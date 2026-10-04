//! Warped B-rep faces that declare no surface (#257).
//!
//! A faceted B-rep face without a surface is its boundary polygon. When its
//! corners are not coplanar it has no single flat surface, so its deviation
//! report must not read exact: it is measured as a warped authored polygon
//! face is, by the slab its corners span about their fit plane (#261), and
//! certified by name. Every expected value here is worked out by hand from
//! the corners.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Tolerance, Transform3, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, DeviationReport, ReferenceMeshCompiler};
use axiolid_model::{GeometryGraphBuilder, GeometryNode, Instance, NodeId};
use axiolid_topology::{
    BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};
use std::collections::HashMap;

const WARPED: &str = "non-planar face without a surface";

/// A solid of faces without surfaces over shared, welded vertices. Each
/// face is its outer ring of vertex indices plus holes; every edge is
/// added once and used in both directions.
fn faceted(corners: &[[f64; 3]], faces: &[(&[usize], &[&[usize]])]) -> BRep<NodeId> {
    let mut brep = BRep::default();
    let verts: Vec<_> = corners
        .iter()
        .map(|c| {
            brep.add_vertex(Vertex {
                position: Vec3::new(c[0], c[1], c[2]),
            })
        })
        .collect();
    let mut edges = HashMap::new();
    let mut ids = Vec::new();
    for &(outer, holes) in faces {
        let mut bounds = Vec::new();
        for (index, ring) in std::iter::once(outer)
            .chain(holes.iter().copied())
            .enumerate()
        {
            let mut uses = Vec::new();
            for i in 0..ring.len() {
                let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                let key = (a.min(b), a.max(b));
                let edge = *edges.entry(key).or_insert_with(|| {
                    brep.add_edge(Edge {
                        start: verts[key.0],
                        end: verts[key.1],
                        curve: None,
                    })
                });
                uses.push(EdgeUse {
                    edge,
                    orientation: if a == key.0 {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    },
                    pcurve: None,
                });
            }
            bounds.push(FaceBound {
                loop_id: brep.add_loop(Loop { edges: uses }),
                orientation: Orientation::Forward,
                outer: index == 0,
            });
        }
        ids.push(brep.add_face(Face {
            surface: None,
            bounds,
            orientation: Orientation::Forward,
        }));
    }
    let shell = brep.add_shell(Shell {
        faces: ids.iter().map(|&f| (f, Orientation::Forward)).collect(),
        closed: faces.len() > 1,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    brep
}

/// The unit box with its bottom corners (0..4, over (0,0), (1,0), (1,1),
/// (0,1)) at `bottom` and its top corners (4..8) at `top`. Every side face
/// stays planar (each lies in `x` or `y` constant).
fn box_with(bottom: [f64; 4], top: [f64; 4]) -> BRep<NodeId> {
    let plan = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let corners: Vec<[f64; 3]> = (0..8)
        .map(|i| {
            let z = if i < 4 { bottom[i] } else { top[i - 4] };
            [plan[i % 4][0], plan[i % 4][1], z]
        })
        .collect();
    let quads: [&[usize]; 6] = [
        &[0, 3, 2, 1],
        &[4, 5, 6, 7],
        &[0, 1, 5, 4],
        &[1, 2, 6, 5],
        &[2, 3, 7, 6],
        &[3, 0, 4, 7],
    ];
    let faces: Vec<(&[usize], &[&[usize]])> = quads.iter().map(|&q| (q, &[][..])).collect();
    faceted(&corners, &faces)
}

fn box_with_top(top: [f64; 4]) -> BRep<NodeId> {
    box_with([0.0; 4], top)
}

/// Compile under `options`, optionally through a uniform scale, with the
/// deviation report.
fn compile_with(
    brep: BRep<NodeId>,
    scale: Option<f64>,
    options: ExecutionOptions,
) -> (TriMesh, DeviationReport) {
    let mut builder = GeometryGraphBuilder::new();
    let mut root = builder.push(GeometryNode::BRep(brep)).expect("push");
    if let Some(factor) = scale {
        root = builder
            .push(GeometryNode::Instance(Instance {
                source: root,
                transform: Transform3::from_scale(Vec3::splat(factor)),
            }))
            .expect("push");
    }
    let graph = builder.finish(vec![root]).expect("finish");
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options)
        .expect("a faceted B-rep compiles");
    (outcome.mesh, report)
}

/// [`compile_with`] at a millimetre, the budget being the tolerance.
fn compile(brep: BRep<NodeId>, scale: Option<f64>) -> (TriMesh, DeviationReport) {
    compile_with(brep, scale, ExecutionOptions::new(Tolerance::MILLIMETRE))
}

/// The warped-face contribution of a report, if any.
fn warp_of(report: &DeviationReport) -> Option<DeviationBound> {
    report
        .contributions
        .iter()
        .find(|c| c.path == DeviationPath::BRepFace && c.detail == WARPED)
        .map(|c| c.bound)
}

fn certified(report: &DeviationReport) -> f64 {
    match warp_of(report) {
        Some(DeviationBound::Certified(warp)) => warp,
        _ => panic!("a warped face is certified by name, never proven or absent: {report:?}"),
    }
}

/// #257's done-when: a box whose top face is a +-5 cm saddle has no single
/// flat top, so it is not reported exact. Its slab about the fit plane
/// z = 1 is 10 cm, certified by name; the planar faces still say `0`.
#[test]
fn a_box_with_one_warped_face_reports_its_slab() {
    let (mesh, report) = compile(box_with_top([1.05, 0.95, 1.05, 0.95]), None);
    assert_eq!(mesh.indices.len(), 36, "every quad gives two triangles");
    let warp = certified(&report);
    assert!((0.1..0.1 + 1e-12).contains(&warp), "warp {warp}");
    assert_eq!(report.bound, Some(warp), "{report:?}");
    assert!(!report.meets_requested(), "{report:?}");
    assert!(
        report
            .contributions
            .iter()
            .any(|c| c.path == DeviationPath::BRepFace
                && c.detail == "plane"
                && c.bound == DeviationBound::Proven(0.0)),
        "{report:?}"
    );
    assert!(
        !report
            .contributions
            .iter()
            .any(|c| matches!(c.bound, DeviationBound::Proven(v) if v > 0.0)),
        "no proven claim covers the warp: {report:?}"
    );
}

/// One corner lifted 5 cm: the top is `z = 1 + h x y`, its fit plane takes
/// the linear part and leaves the corners `h/4` alternately above and
/// below it, a slab of `h/2`, the gap between its two diagonal
/// triangulations at the centre. Along that plane's Newell normal
/// `(-h, -h, 2)` it is a fraction of a percent narrower.
#[test]
fn a_lifted_corner_reports_half_its_lift() {
    let h = 0.05;
    let (_, report) = compile(box_with_top([1.0, 1.0, 1.0 + h, 1.0]), None);
    let want = 0.5 * h * 2.0 / (4.0 + 2.0 * h * h).sqrt();
    let warp = certified(&report);
    assert!(
        (want..want + 1e-12).contains(&warp),
        "warp {warp}, want {want}"
    );
}

/// The regression axioval found (#257): a faceted box with one corner
/// lifted 5 cm reported every face `Proven(0)`. It now reports the
/// measured slab, about 2.5 cm (see above), and whether the mesh meets the
/// requested budget follows from that bound: not at a millimetre, yes at
/// a 3 cm chord budget, no at 2 cm.
#[test]
fn a_box_with_one_corner_lifted_five_centimetres_is_not_exact() {
    let lifted = || box_with_top([1.0, 1.0, 1.05, 1.0]);
    let (_, report) = compile(lifted(), None);
    assert!(
        report.bound.is_some_and(|b| b > 0.02),
        "the warp is measured, not 0: {report:?}"
    );
    let warp = certified(&report);
    assert_eq!(report.bound, Some(warp), "{report:?}");
    assert!(!report.meets_requested(), "not within 1 mm: {report:?}");
    for (budget, meets) in [(0.03, true), (0.02, false)] {
        let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
            .with_chord_error(budget)
            .expect("a valid budget");
        let (_, report) = compile_with(lifted(), None, options);
        assert_eq!(report.bound, Some(warp), "{report:?}");
        assert_eq!(
            report.meets_requested(),
            meets,
            "budget {budget}: {report:?}"
        );
    }
}

/// A planar box, and one whose top is warped within the linear tolerance,
/// stay exact: no warped-face contribution and every bound `Proven(0)`.
#[test]
fn a_planar_box_stays_exact() {
    for top in [[1.0; 4], [1.0, 1.0, 1.0 + 9e-4, 1.0]] {
        let (_, report) = compile(box_with_top(top), None);
        assert_eq!(warp_of(&report), None, "{report:?}");
        assert_eq!(report.bound, Some(0.0), "{report:?}");
        assert!(
            report
                .contributions
                .iter()
                .all(|c| c.bound == DeviationBound::Proven(0.0)),
            "{report:?}"
        );
    }
}

/// The warp is a world length: an instance scaled by two doubles it.
#[test]
fn a_scaled_instance_doubles_the_warp() {
    let (_, report) = compile(box_with_top([1.05, 0.95, 1.05, 0.95]), Some(2.0));
    let warp = certified(&report);
    assert!((0.2..0.2 + 1e-12).contains(&warp), "warp {warp}");
}

/// A hole off its face's plane warps the face: every corner is measured.
/// The flat outer ring fixes the fit plane at z = 0, so the slab runs to
/// the hole corner's 5 cm.
#[test]
fn a_hole_off_the_plane_of_its_face_is_measured() {
    let corners = [
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [4.0, 4.0, 0.0],
        [0.0, 4.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0, 2.0, 0.0],
        [2.0, 2.0, 0.05],
        [2.0, 1.0, 0.0],
    ];
    let (_, report) = compile(
        faceted(&corners, &[(&[0, 1, 2, 3], &[&[4, 5, 6, 7]])]),
        None,
    );
    let warp = certified(&report);
    assert!((0.05..0.05 + 1e-12).contains(&warp), "warp {warp}");
}

/// Of two warped faces the worst sets the bound, whichever is met first:
/// a +-2 cm saddle and a +-5 cm one, bottom and top in both orders.
#[test]
fn the_worst_warped_face_sets_the_bound() {
    let saddle = |h: f64, z: f64| [z + h, z - h, z + h, z - h];
    for (bottom, top) in [(0.02, 0.05), (0.05, 0.02)] {
        let (_, report) = compile(box_with(saddle(bottom, 0.0), saddle(top, 1.0)), None);
        let warp = certified(&report);
        assert!((0.1..0.1 + 1e-12).contains(&warp), "warp {warp}");
    }
}
