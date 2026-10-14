//! Windowed stations along an offset by distances beside an ellipse
//! (#307), the deck of openbimrs/ifc#443: a sectioned solid along the
//! ellipse of semi-axes 10 and 6 about `(30, 0)`, offset by `0.5` to
//! station 2, widening to `0.75` at station 20 and held to the end, its
//! rectangles of area 2 at stations 0 and 5. Its windowed sections need
//! the varying piece's length certified to a window's certification; they
//! resolve, and the windowed and window-free decks agree.

use axiolid_contracts::{ExecutionOptions, GeomResult};
use axiolid_core::{Frame2, Point2, Scalar, Tolerance, Vec2};
use axiolid_curve::{Curve2, Ellipse2};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::station::{curve_path, resolve};
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    CurveRelation, CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId,
    RunStation, SeamSide, SectionAtStation, SolidOperation, Station, StationFrame, StationOffsets,
};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_reference::station::CompositeBasis;

/// A window whose certification, `1e-3` of it, is `1e-8`.
const WINDOW: Scalar = 1e-5;

fn lateral(distance: Scalar) -> StationOffsets {
    StationOffsets::new(distance, 0.0, 0.0)
}

/// The deck's directrix: the ellipse offset by distances at stations.
fn offset(b: &mut GeometryGraphBuilder) -> NodeId {
    let ellipse = b
        .push_value(Curve2::Ellipse(Ellipse2 {
            frame: Frame2 {
                origin: Point2::new(30.0, 0.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            semi_axis_x: 10.0,
            semi_axis_y: 6.0,
        }))
        .unwrap();
    b.push(GeometryNode::CurveRelation(
        CurveRelation::OffsetsAtStations {
            basis: ellipse,
            stations: vec![
                RunStation::at_basis_start(lateral(0.5)),
                RunStation::new(Station::new(2.0, lateral(0.5))),
                RunStation::new(Station::new(20.0, lateral(0.75))),
                RunStation::at_basis_end(lateral(0.75)),
            ],
            frame: StationFrame::Section,
        },
    ))
    .unwrap()
}

/// The deck with its sections windowed by `window`, and the stations at
/// `along` on its directrix with the same window: the graph, the solid and
/// the stations.
fn deck(window: Scalar, along: &[Scalar]) -> (GeometryGraph, NodeId, Vec<NodeId>) {
    let mut b = GeometryGraphBuilder::new();
    let directrix = offset(&mut b);
    let rectangle = b
        .push_value(Profile::Rectangle(RectangleProfile {
            x: 2.0,
            y: 1.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }))
        .unwrap();
    let sections = [0.0, 5.0]
        .iter()
        .map(|&d| SectionAtStation::new(rectangle, Station::at(d)).with_seam_window(window))
        .collect();
    let solid = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix,
                sections,
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let stations: Vec<NodeId> = along
        .iter()
        .map(|&d| {
            b.push_value(
                CurveStation::new(directrix, Station::at(d))
                    .with_seam_side(SeamSide::Incoming)
                    .with_seam_window(window),
            )
            .unwrap()
        })
        .collect();
    let mut roots = vec![solid];
    roots.extend(&stations);
    (b.finish(roots).unwrap(), solid, stations)
}

fn compile(graph: &GeometryGraph, root: NodeId) -> GeomResult<TriMesh> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(
        graph,
        root,
        &ExecutionOptions::new(Tolerance::MILLIMETRE),
    )
}

/// The kernel's own distance of the directrix's joint at station 20: the
/// offset's length to it.
fn joint_at_20() -> Scalar {
    let mut b = GeometryGraphBuilder::new();
    let directrix = offset(&mut b);
    let graph = b.finish(vec![directrix]).unwrap();
    let path = curve_path(&graph, directrix).unwrap();
    let basis = CompositeBasis::from_path(&path).unwrap();
    let pieces = basis.pieces();
    pieces[0].length() + pieces[1].length()
}

#[test]
fn the_windowed_deck_resolves_and_agrees_with_the_window_free_one() {
    let joint = joint_at_20();
    let along = [30.0, joint - 0.4 * WINDOW, joint];
    let (windowed, solid, stations) = deck(WINDOW, &along);
    let (plain, plain_solid, plain_stations) = deck(0.0, &along);
    let mesh = compile(&windowed, solid).unwrap();
    // No seam lies within a window of the sections at 0 and 5: they read
    // where they stand, bitwise as without a window.
    let free = compile(&plain, plain_solid).unwrap();
    assert_eq!(mesh.indices, free.indices);
    assert_eq!(mesh.positions, free.positions);
    // Pappus: the rectangles' area, 2, times the 5 between them.
    let volume = volume_properties(&mesh, Tolerance::MILLIMETRE)
        .unwrap()
        .signed_volume;
    assert!((volume - 10.0).abs() <= 1e-3, "{volume}");
    // A station 30 along, past the varying piece, reads where it stands.
    assert_eq!(
        resolve(&windowed, stations[0]).unwrap(),
        resolve(&plain, plain_stations[0]).unwrap()
    );
    // One within its window of the joint at station 20 stands on it.
    let on = resolve(&plain, plain_stations[2]).unwrap();
    assert_eq!(resolve(&windowed, stations[1]).unwrap(), on);
    assert_ne!(resolve(&plain, plain_stations[1]).unwrap(), on);
}
