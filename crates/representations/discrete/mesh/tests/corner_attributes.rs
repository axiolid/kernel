//! Corner-indexed attribute channels (#112).
//!
//! The fixture is the IFC4 specification's own worked example for
//! `IfcIndexedTriangleTextureMap` (IFC4 ADD2 TC1, Figure 415): a 1x1x2 box,
//! 8 positions, 12 triangles, 8 texture vertices. Indices are the spec's
//! 1-based `CoordIndex` / `TexCoordIndex`, shifted to 0-based here. Every
//! one of the 8 positions carries 3 different UVs, so no per-vertex channel
//! can hold this data without splitting positions.

use axiolid_core::{Point3, Tolerance};
use axiolid_mesh::{audit_mesh, AttributeChannel, Blend, MeshValidationError, TriMesh};

const COORD_INDEX: [[u32; 3]; 12] = [
    [1, 6, 5],
    [1, 2, 6],
    [6, 2, 7],
    [7, 2, 3],
    [7, 8, 6],
    [6, 8, 5],
    [5, 8, 1],
    [1, 8, 4],
    [4, 2, 1],
    [2, 4, 3],
    [4, 8, 7],
    [7, 3, 4],
];
const TEX_COORD_INDEX: [[u32; 3]; 12] = [
    [1, 4, 3],
    [1, 2, 4],
    [3, 1, 4],
    [4, 1, 2],
    [8, 7, 6],
    [6, 7, 5],
    [4, 3, 2],
    [2, 3, 1],
    [5, 8, 7],
    [8, 5, 6],
    [2, 4, 3],
    [3, 1, 2],
];
const TEX_COORDS: [[f64; 2]; 8] = [
    [0.0, -0.5],
    [1.0, -0.5],
    [0.0, 1.5],
    [1.0, 1.5],
    [0.0, 0.0],
    [0.0, 1.0],
    [1.0, 0.0],
    [1.0, 1.0],
];

fn one_based(rows: &[[u32; 3]]) -> Vec<u32> {
    rows.iter().flatten().map(|&i| i - 1).collect()
}

fn figure_415() -> TriMesh {
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    let mut mesh = TriMesh::new(
        vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(1.0, 1.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, 0.0, 2.0),
            p(1.0, 0.0, 2.0),
            p(1.0, 1.0, 2.0),
            p(0.0, 1.0, 2.0),
        ],
        one_based(&COORD_INDEX),
    );
    mesh.attributes.push(AttributeChannel::corner_indexed(
        "uv",
        TEX_COORDS.iter().flatten().copied().collect(),
        2,
        Blend::Linear,
        one_based(&TEX_COORD_INDEX),
    ));
    mesh
}

#[test]
fn the_spec_box_validates_and_stays_closed_with_the_channel_attached() {
    let mesh = figure_415();
    mesh.validate_structure()
        .expect("corner-indexed channel validates");
    let health = audit_mesh(&mesh, Tolerance::MILLIMETRE);
    assert!(
        health.is_closed_two_manifold(),
        "UV seams must not change adjacency: {health:?}"
    );
    assert_eq!(mesh.positions.len(), 8, "no position was split for a seam");
}

#[test]
fn every_corner_round_trips_its_own_value() {
    let mesh = figure_415();
    let channel = &mesh.attributes[0];
    for (triangle, (coords, texs)) in COORD_INDEX.iter().zip(TEX_COORD_INDEX.iter()).enumerate() {
        for k in 0..3 {
            let corner = triangle * 3 + k;
            assert_eq!(mesh.indices[corner], coords[k] - 1);
            let want = TEX_COORDS[(texs[k] - 1) as usize];
            assert_eq!(
                channel.at_corner(&mesh.indices, corner),
                Some(&want[..]),
                "triangle {triangle} corner {k}"
            );
        }
    }
}

#[test]
fn one_position_really_carries_three_values() {
    // The property that makes a per-vertex channel insufficient, pinned so
    // a fixture edit that loses it cannot pass silently.
    let mesh = figure_415();
    let channel = &mesh.attributes[0];
    for position in 0..8u32 {
        let mut seen: Vec<[u64; 2]> = (0..mesh.indices.len())
            .filter(|&c| mesh.indices[c] == position)
            .map(|c| {
                let uv = channel.at_corner(&mesh.indices, c).expect("mapped");
                [uv[0].to_bits(), uv[1].to_bits()]
            })
            .collect();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(seen.len(), 3, "position {position}");
    }
}

#[test]
fn a_per_vertex_channel_reads_through_at_corner_too() {
    let mesh = figure_415();
    let ids: Vec<f64> = (0..8).map(f64::from).collect();
    let channel = AttributeChannel::new("id", ids, 1, Blend::Nearest);
    for (corner, &position) in mesh.indices.iter().enumerate() {
        assert_eq!(
            channel.at_corner(&mesh.indices, corner),
            Some(&[f64::from(position)][..])
        );
    }
}

fn with_corners(corners: Vec<u32>) -> TriMesh {
    let mut mesh = figure_415();
    mesh.attributes[0].corner_indices = Some(corners);
    mesh
}

#[test]
fn an_unmapped_triangle_is_valid_and_reads_as_no_value() {
    // IfcTrafficSignLibrary maps 12 of 48 triangles: the rest carry nothing.
    let mut corners = one_based(&TEX_COORD_INDEX);
    corners[..3].fill(AttributeChannel::UNMAPPED);
    let mesh = with_corners(corners);
    mesh.validate_structure()
        .expect("an unmapped triangle is representable");
    for corner in 0..3 {
        assert_eq!(mesh.attributes[0].at_corner(&mesh.indices, corner), None);
    }
    assert!(mesh.attributes[0].at_corner(&mesh.indices, 3).is_some());
}

#[test]
fn malformed_corner_channels_are_refused_by_name() {
    let mut partial = one_based(&TEX_COORD_INDEX);
    partial[4] = AttributeChannel::UNMAPPED;
    let mut out_of_range = one_based(&TEX_COORD_INDEX);
    out_of_range[7] = 8;
    let mut short = one_based(&TEX_COORD_INDEX);
    short.pop();
    let mut ragged = figure_415();
    ragged.attributes[0].values.pop();

    let cases = [
        (with_corners(partial), "partial"),
        (with_corners(out_of_range), "out of range"),
        (with_corners(short), "short"),
        (ragged, "ragged"),
    ];
    for (mesh, what) in cases {
        let error = mesh.validate_structure().expect_err(what);
        let named = match &error {
            MeshValidationError::AttributePartiallyMapped { name, triangle } => {
                assert_eq!((what, *triangle), ("partial", 1));
                name
            }
            MeshValidationError::AttributeCornerIndexOutOfRange {
                name,
                index,
                value_count,
            } => {
                assert_eq!((what, *index, *value_count), ("out of range", 8, 8));
                name
            }
            MeshValidationError::AttributeCornerCount {
                name,
                expected,
                actual,
            } => {
                assert_eq!((what, *expected, *actual), ("short", 36, 35));
                name
            }
            MeshValidationError::AttributeRaggedValues {
                name,
                values,
                width,
            } => {
                assert_eq!((what, *values, *width), ("ragged", 15, 2));
                name
            }
            other => panic!("{what}: unexpected {other:?}"),
        };
        assert_eq!(named, "uv", "{what}");
    }
}
