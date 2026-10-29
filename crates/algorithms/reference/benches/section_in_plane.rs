//! In-plane faces of a section request (ADR 0093, #366): what the exact
//! arrangement costs against the grid of faces it merges.
//!
//! Four cases, each timed on its own wall clock, best of a few runs:
//!
//! - `grid`: a sheet of 64 x 64 cells, 8,192 triangles, all in the plane and
//!   conforming. Nothing crosses, nothing overlaps.
//! - `box on a slab`: a slab and a box standing on it, 24 triangles, run
//!   many times: the footprint is a hole of the slab top and a `Both` region.
//! - `footprint on a grid`: the 8,192-triangle grid with one square of the
//!   opposite facing over its middle.
//! - `skew grids`: two sheets of 24 x 24 cells in one plane, the second turned
//!   by 7 degrees and facing the other way: about three thousand crossings.
//!
//! Run: `cargo bench -p axiolid-reference --bench section_in_plane`

use std::hint::black_box;
use std::time::Instant;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame3, Point3, Tolerance, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::{
    MeshPlaneSection, RegionLimits, SectionLimits, SectionRequest, SectionRole,
};
use axiolid_reference::section::ScalarSection;

/// A conforming grid of `cells` x `cells` squares of side `size`, each split
/// on its diagonal, in the plane z = 0, turned by `angle` about the origin
/// and moved by `shift`. A `flipped` grid winds clockwise.
fn grid(cells: usize, size: f64, angle: f64, shift: (f64, f64), flipped: bool) -> TriMesh {
    let (s, c) = angle.sin_cos();
    let row = cells + 1;
    let mut positions = Vec::with_capacity(row * row);
    for j in 0..row {
        for i in 0..row {
            let (x, y) = (i as f64 * size, j as f64 * size);
            positions.push(Point3::new(
                c * x - s * y + shift.0,
                s * x + c * y + shift.1,
                0.0,
            ));
        }
    }
    let mut indices = Vec::new();
    for j in 0..cells {
        for i in 0..cells {
            let a = (j * row + i) as u32;
            let (b, d) = (a + 1, a + row as u32);
            let e = d + 1;
            if flipped {
                indices.extend_from_slice(&[a, e, b, a, d, e]);
            } else {
                indices.extend_from_slice(&[a, b, e, a, e, d]);
            }
        }
    }
    TriMesh::new(positions, indices)
}

/// Several meshes as separate shells.
fn combine(parts: &[TriMesh]) -> TriMesh {
    let mut positions = Vec::new();
    let mut indices = Vec::new();
    for part in parts {
        let offset = positions.len() as u32;
        positions.extend(part.positions.iter().copied());
        indices.extend(part.indices.iter().map(|index| index + offset));
    }
    TriMesh::new(positions, indices)
}

/// A box `x` by `y` by `z` as 8 vertices and 12 outward triangles.
fn cuboid(x: (f64, f64), y: (f64, f64), z: (f64, f64)) -> TriMesh {
    let positions = vec![
        Point3::new(x.0, y.0, z.0),
        Point3::new(x.1, y.0, z.0),
        Point3::new(x.1, y.1, z.0),
        Point3::new(x.0, y.1, z.0),
        Point3::new(x.0, y.0, z.1),
        Point3::new(x.1, y.0, z.1),
        Point3::new(x.1, y.1, z.1),
        Point3::new(x.0, y.1, z.1),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(positions, indices)
}

fn request(z: f64, role: SectionRole) -> SectionRequest {
    let frame = Frame3 {
        origin: Point3::new(0.0, 0.0, z),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    };
    SectionRequest::new(
        frame,
        SectionLimits::new(1_000_000, 1_000_000, 1_000_000, 1_000_000),
        RegionLimits::new(1_000_000, 10_000_000, 1_000_000_000),
    )
    .with_role(role)
}

/// Best of `runs` wall-clock times of `repeat` sections, per section, with
/// what the section answered.
fn time(mesh: &TriMesh, request: &SectionRequest, runs: usize, repeat: usize) -> (f64, String) {
    let options = ExecutionOptions::new(Tolerance::METRE);
    let mut best = f64::INFINITY;
    let mut answer = String::new();
    for _ in 0..runs {
        let start = Instant::now();
        let mut last = None;
        for _ in 0..repeat {
            last = Some(black_box(ScalarSection::new().section_request(
                black_box(mesh),
                black_box(request),
                &options,
            )));
        }
        let elapsed = start.elapsed().as_secs_f64() * 1e3 / repeat as f64;
        best = best.min(elapsed);
        answer = match last {
            Some(Ok(section)) => format!(
                "{} regions, {} crossings, {} pairs",
                section.regions.len(),
                section.evidence.in_plane_crossings,
                section.evidence.in_plane_pairs
            ),
            Some(Err(error)) => format!("refused: {error:?}"),
            None => String::new(),
        };
    }
    (best, answer)
}

fn main() {
    println!("in-plane faces of a section request\n");
    println!(
        "{:<24} {:>10} {:>12}  answer",
        "case", "triangles", "ms/section"
    );
    let big = grid(64, 1.0, 0.0, (0.0, 0.0), false);
    let footprint = grid(1, 8.0, 0.0, (28.0, 28.0), true);
    let cases: Vec<(&str, TriMesh, SectionRequest, usize, usize)> = vec![
        ("grid", big.clone(), request(0.0, SectionRole::Sheet), 5, 1),
        (
            "box on a slab",
            combine(&[
                cuboid((0.0, 4.0), (0.0, 4.0), (0.0, 1.0)),
                cuboid((1.0, 2.0), (1.0, 2.0), (1.0, 2.0)),
            ]),
            request(1.0, SectionRole::Solid),
            5,
            200,
        ),
        (
            "footprint on a grid",
            combine(&[big, footprint]),
            request(0.0, SectionRole::Sheet),
            5,
            1,
        ),
        (
            "skew grids",
            combine(&[
                grid(24, 1.0, 0.0, (0.0, 0.0), false),
                grid(24, 1.0, 7.0_f64.to_radians(), (3.0, 1.0), true),
            ]),
            request(0.0, SectionRole::Sheet),
            3,
            1,
        ),
    ];
    for (name, mesh, request, runs, repeat) in cases {
        let (ms, answer) = time(&mesh, &request, runs, repeat);
        println!(
            "{:<24} {:>10} {:>12.3}  {answer}",
            name,
            mesh.indices.len() / 3,
            ms
        );
    }
}
