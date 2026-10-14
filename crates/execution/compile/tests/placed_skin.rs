//! A wall cut in its own frame and placed 5600 km out by one `Instance`
//! (#312): its doors a residue `r` short of one face, and as far past the
//! other, as an IFC export leaves them through composed placements.
//!
//! The snap reached only the rounding of the operands' own coordinates
//! (`16 f64::EPSILON * 3.67`, about `1.3e-14`), so a skin of `2e-14` to
//! `1e-11` was kept as two faces that thin; placed where a coordinate rounds
//! to about `1e-9`, the faces crossed and the volume was refused as
//! self-intersecting. The snap now reaches the rounding of the coordinates
//! the result is placed at (about `2e-8` there): every such residue closes
//! placed, at the origin the operands' own rounding still decides, and a
//! micrometre stays a real skin.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Scalar, Tolerance, Transform2, Transform3, Vec2, Vec3};
use axiolid_inspect::{enclosed_volume, VolumeInterval};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::deviation::SNAPPED_OPERANDS;
use axiolid_mesh_compile::{DeviationBound, DeviationReport, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::{MeshClosure, MeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{Profile, RectangleProfile};

const DOOR_AREA: Scalar = 1.01 * 2.26;

/// More than the wall's surface area, which bounds how far moving every
/// point of it by `d` changes its volume (`area * d`).
const AREA: Scalar = 40.0;

/// Residues the snap closes once placed.
const RESIDUES: [Scalar; 7] = [0.0, 5e-15, 2e-14, 1e-13, 1e-12, 1e-11, 1e-10];

/// The wall less both doors cut through it.
fn net() -> Scalar {
    3.65 * 0.25 * 3.67 - 2.0 * DOOR_AREA * 0.25
}

fn rectangle(x: Scalar, y: Scalar, transform: Transform2) -> GeometryNode {
    GeometryNode::Profile(Profile::Derived {
        basis: Box::new(Profile::Rectangle(RectangleProfile {
            x,
            y,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        })),
        transform,
    })
}

/// Turned a quarter and then to grid north, 5600 km out.
fn site() -> Transform3 {
    let (s, c) = (0.039_608_610_093_419_1_f64, 0.999_215_271_103_513_f64);
    let norm = s.hypot(c);
    let (s, c) = (s / norm, c / norm);
    Transform3::from_cols(
        Vec3::new(-s, c, 0.0),
        Vec3::new(-c, -s, 0.0),
        Vec3::Z,
        Vec3::new(600_014.096_183_109_1, 5_599_984.596_242_446, 204.3),
    )
}

/// The rounding of a coordinate placed `size` from the origin.
fn ulp(size: Scalar) -> Scalar {
    size * Scalar::EPSILON
}

/// How far placing the wall at the site may move its volume: each placed
/// coordinate rounds by about an ulp of the site's size.
fn placement_slack() -> Scalar {
    AREA * 4.0 * ulp(5.7e6)
}

/// The wall `3.65 x 0.25 x 3.67`, its faces at `y = +-0.125`, less two door
/// openings extruded `0.25` along `-y` from `y = 0.125 - residue`; and the
/// root placing it by `world`. Returns the graph, the root and the net.
fn wall(residue: Scalar, world: Transform3) -> (GeometryGraph, NodeId, NodeId) {
    nested_wall(residue, &[world])
}

/// The wall placed by one instance per entry of `worlds`, each inside the
/// next.
fn nested_wall(residue: Scalar, worlds: &[Transform3]) -> (GeometryGraph, NodeId, NodeId) {
    let mut builder = GeometryGraphBuilder::new();
    let mut push = |node| builder.push(node).expect("a valid node");
    let profile = push(rectangle(
        3.65,
        0.25,
        Transform2::from_cols(Vec2::NEG_X, Vec2::NEG_Y, Vec2::new(1.825, 0.0)),
    ));
    let mut net = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: 3.67,
    }));
    for x in [3.325, 1.35] {
        let profile = push(rectangle(
            1.01,
            2.26,
            Transform2::from_cols(Vec2::Y, Vec2::NEG_X, Vec2::new(1.13, 0.505)),
        ));
        let solid = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 0.25,
        }));
        let placed = push(GeometryNode::Instance(Instance {
            source: solid,
            transform: Transform3::from_cols(
                Vec3::Z,
                Vec3::NEG_X,
                Vec3::NEG_Y,
                Vec3::new(x, 0.125 - residue, 0.2),
            ),
        }));
        net = push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: net,
            right: placed,
            operator: BooleanOperator::Difference,
        }));
    }
    let mut root = net;
    for &transform in worlds {
        root = push(GeometryNode::Instance(Instance {
            source: root,
            transform,
        }));
    }
    (
        builder.finish(vec![root]).expect("a valid graph"),
        root,
        net,
    )
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn compile(residue: Scalar, world: Transform3) -> (TriMesh, DeviationReport) {
    let (graph, root, _) = wall(residue, world);
    let (outcome, report) = compiler()
        .compile_mesh_with_deviation(&graph, root, &options())
        .unwrap_or_else(|e| panic!("residue {residue:e}: {e}"));
    assert_eq!(outcome.closure, MeshClosure::Solid, "residue {residue:e}");
    (outcome.mesh, report)
}

/// The largest snap the report names, `0` when nothing was snapped.
fn snapped(report: &DeviationReport) -> Scalar {
    report
        .contributions
        .iter()
        .filter(|c| c.detail == SNAPPED_OPERANDS)
        .map(|c| match c.bound {
            DeviationBound::Certified(moved) => moved,
            other => panic!("a snap is certified: {other:?}"),
        })
        .fold(0.0, Scalar::max)
}

/// `expected` within `slack` of the certified interval.
fn assert_volume(volume: VolumeInterval, expected: Scalar, slack: Scalar, what: &str) {
    assert!(
        volume.lower - slack <= expected && expected <= volume.upper + slack,
        "{what}: [{}, {}] does not hold {expected} within {slack:e}",
        volume.lower,
        volume.upper
    );
}

#[test]
fn a_residue_below_the_placed_rounding_measures_placed_and_at_the_origin() {
    for residue in RESIDUES {
        let (origin, origin_report) = compile(residue, Transform3::IDENTITY);
        let (placed, placed_report) = compile(residue, site());
        let at_origin = enclosed_volume(&origin)
            .unwrap_or_else(|e| panic!("residue {residue:e} at the origin: {e:?}"));
        let at_site = enclosed_volume(&placed)
            .unwrap_or_else(|e| panic!("residue {residue:e} placed: {e:?}"));
        let (moved_origin, moved_site) = (snapped(&origin_report), snapped(&placed_report));
        // Placed, the skin closes: the doors moved onto the wall's faces by
        // the residue, never beyond the placed rounding, and reported.
        if residue > 0.0 {
            assert!(
                moved_site >= 0.5 * residue && moved_site <= 16.0 * ulp(8.4e6),
                "residue {residue:e}: moved {moved_site:e}"
            );
            assert!(
                placed_report.bound.is_some_and(|b| b >= moved_site),
                "{placed_report:?}"
            );
        } else {
            assert_eq!(moved_site, 0.0, "flush needs no snap");
        }
        assert_volume(
            at_site,
            net(),
            AREA * moved_site + placement_slack(),
            &format!("residue {residue:e} placed"),
        );
        // At the origin the operands' own rounding decides, as before: the
        // exported residue closes, a thicker one stays as its skin.
        if residue > 16.0 * ulp(3.67) {
            assert_eq!(moved_origin, 0.0, "residue {residue:e}: {origin_report:?}");
        }
        // Unsnapped, the operands are cut as given: the mesh lies within the
        // residue of the snapped one, and the skin is kept.
        let rounding = AREA * ulp(4.0);
        assert_volume(
            at_origin,
            net(),
            AREA * (moved_origin + residue) + rounding,
            &format!("residue {residue:e} at the origin"),
        );
        if moved_origin == 0.0 && residue > 0.0 {
            assert!(
                at_origin.lower > net(),
                "residue {residue:e}: the skin is kept"
            );
        }
        // Both agree within the residue, the reported snaps and the
        // placement's rounding, and within the reported deviations.
        let slack = AREA * (moved_origin + moved_site + residue) + rounding + placement_slack();
        let agree = |slack: Scalar| {
            at_origin.lower - slack <= at_site.upper && at_site.lower - slack <= at_origin.upper
        };
        assert!(
            agree(slack),
            "residue {residue:e}: [{}, {}] at the origin, [{}, {}] placed",
            at_origin.lower,
            at_origin.upper,
            at_site.lower,
            at_site.upper
        );
        let reported =
            origin_report.bound.expect("bounded") + placed_report.bound.expect("bounded");
        assert!(agree(AREA * reported), "residue {residue:e}");
    }
}

#[test]
fn a_skin_above_the_placed_rounding_stays_a_skin() {
    // Fifty times the placed reach: authored, kept, and measured.
    let residue = 1e-6;
    let skin = 2.0 * DOOR_AREA * residue;
    for (name, world) in [("origin", Transform3::IDENTITY), ("site", site())] {
        let (mesh, report) = compile(residue, world);
        assert_eq!(snapped(&report), 0.0, "{name}: {report:?}");
        let volume = enclosed_volume(&mesh).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_volume(volume, net() + skin, placement_slack(), name);
        assert!(
            volume.upper + placement_slack() < net() + skin + 0.5 * skin
                && volume.lower - placement_slack() > net() + 0.5 * skin,
            "{name}: the skin is kept: [{}, {}]",
            volume.lower,
            volume.upper
        );
    }
}

#[test]
fn one_source_placed_near_and_far_compiles_once_per_scale() {
    // The unplaced net and the same net placed at the site, in one batch
    // sharing one cache: the placed one must not reuse the net compiled
    // for the origin, whose skin is kept and crosses itself once placed.
    let residue = 1e-12;
    let (graph, root, unplaced) = wall(residue, site());
    let meshes = compiler()
        .compile_mesh_batch(&graph, &[unplaced, root], &options())
        .expect("both compile");
    let near = enclosed_volume(&meshes[0]).expect("the net measures");
    assert_volume(near, net_volume_with_skin(residue), AREA * ulp(4.0), "near");
    let far = enclosed_volume(&meshes[1]).expect("the placed net measures");
    assert_volume(
        far,
        net(),
        AREA * 16.0 * ulp(8.4e6) + placement_slack(),
        "far",
    );
    // And the other way round.
    let meshes = compiler()
        .compile_mesh_batch(&graph, &[root, unplaced], &options())
        .expect("both compile");
    enclosed_volume(&meshes[0]).expect("the placed net measures");
    let near = enclosed_volume(&meshes[1]).expect("the net measures");
    assert_volume(
        near,
        net_volume_with_skin(residue),
        AREA * ulp(4.0),
        "near, second",
    );
}

#[test]
fn a_placement_split_over_nested_instances_reaches_as_far() {
    // The site's turn inside, its translation outside: the boolean's placed
    // scale is the outer instance's, carried through the inner one.
    for residue in [1e-12, 1e-11] {
        let site = site();
        let (graph, root, _) = nested_wall(
            residue,
            &[
                Transform3::from_mat3(site.matrix3),
                Transform3::from_translation(site.translation),
            ],
        );
        let (outcome, report) = compiler()
            .compile_mesh_with_deviation(&graph, root, &options())
            .unwrap_or_else(|e| panic!("residue {residue:e}: {e}"));
        let moved = snapped(&report);
        assert!(moved >= 0.5 * residue, "residue {residue:e}: {report:?}");
        let volume = enclosed_volume(&outcome.mesh)
            .unwrap_or_else(|e| panic!("residue {residue:e} nested: {e:?}"));
        assert_volume(volume, net(), AREA * moved + placement_slack(), "nested");
    }
}

fn net_volume_with_skin(residue: Scalar) -> Scalar {
    net() + 2.0 * DOOR_AREA * residue
}
