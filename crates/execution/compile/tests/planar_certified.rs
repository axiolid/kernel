//! Planar faces through the certified ear clipper (#260).
//!
//! Authored polygon faces and planar B-rep faces used earcut, which left
//! the T-junction #253 found in profile caps: for two holes in one band a
//! triangle edge ran along the band's bottom line past both holes' inner
//! corners. Here the #253 layouts become the top and bottom faces of
//! prisms, authored and as B-reps, which must compile to closed
//! two-manifolds of volume area x depth. Faces are surface patches, so
//! rings touching at a single point triangulate there (#262) with their
//! exact area, while the extrusion of the same rings, a solid, is refused
//! by name. Planar faces keep their deviation contribution: exact.
//!
//! Coordinates are dyadic, so areas and volumes are exact sums.

use std::collections::HashMap;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Interval, Point2, Point3, Tolerance, Vec3};
use axiolid_curve::{Curve2, Line2};
use axiolid_measure::volume_properties;
use axiolid_mesh::{PolygonFace, PolygonMesh, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId, SolidOperation};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};
use axiolid_topology::{
    BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};

const DEPTH: f64 = 2.0;

/// A face as rings of corner indices: the outer ring and the holes.
type FaceRings = (Vec<usize>, Vec<Vec<usize>>);

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Point2> {
    vec![p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)]
}

fn twice_ring_area(ring: &[Point2]) -> f64 {
    (0..ring.len())
        .map(|i| {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            a.x * b.y - b.x * a.y
        })
        .sum()
}

fn twice_area(mesh: &TriMesh) -> f64 {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            (b - a).cross(c - a).length()
        })
        .sum()
}

fn close(got: f64, want: f64, what: &str) {
    assert!(
        (got - want).abs() <= 1e-12 * want.abs().max(1.0),
        "{what}: got {got}, want {want}"
    );
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

fn compile_node(node: GeometryNode) -> Result<TriMesh, GeomError> {
    let mut b = GeometryGraphBuilder::new();
    let root = b.push(node).unwrap();
    let graph = b.finish(vec![root]).unwrap();
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(&graph, root, &options())
}

/// A prism over `outer` (counter-clockwise) and `holes` (clockwise):
/// corners bottom then top, faces as rings of corner indices (bottom,
/// top, then one wall per ring edge), each facing out.
fn prism(outer: &[Point2], holes: &[Vec<Point2>]) -> (Vec<Point3>, Vec<FaceRings>) {
    let rings: Vec<&[Point2]> = core::iter::once(outer)
        .chain(holes.iter().map(Vec::as_slice))
        .collect();
    let count: usize = rings.iter().map(|r| r.len()).sum();
    let mut corners = Vec::with_capacity(2 * count);
    for z in [0.0, DEPTH] {
        for ring in &rings {
            corners.extend(ring.iter().map(|q| Point3::new(q.x, q.y, z)));
        }
    }
    let mut ranges = Vec::new();
    let mut start = 0;
    for ring in &rings {
        ranges.push(start..start + ring.len());
        start += ring.len();
    }
    let reversed = |r: &core::ops::Range<usize>| -> Vec<usize> { r.clone().rev().collect() };
    let lifted =
        |r: &core::ops::Range<usize>| -> Vec<usize> { r.clone().map(|i| i + count).collect() };
    let mut faces = vec![
        (
            reversed(&ranges[0]),
            ranges[1..].iter().map(reversed).collect(),
        ),
        (lifted(&ranges[0]), ranges[1..].iter().map(lifted).collect()),
    ];
    for range in &ranges {
        let n = range.len();
        for k in 0..n {
            let (a, b) = (range.start + k, range.start + (k + 1) % n);
            faces.push((vec![a, b, b + count, a + count], Vec::new()));
        }
    }
    (corners, faces)
}

fn authored(corners: &[Point3], faces: &[FaceRings]) -> PolygonMesh {
    let index = |v: &Vec<usize>| v.iter().map(|&i| i as u32).collect::<Vec<u32>>();
    PolygonMesh {
        positions: corners.to_vec(),
        faces: faces
            .iter()
            .map(|(outer, holes)| PolygonFace {
                outer: index(outer),
                holes: holes.iter().map(index).collect(),
            })
            .collect(),
    }
}

/// Faces without surfaces over shared, welded vertices; every edge added
/// once and used in both directions (as in `brep_warped_faces.rs`). A
/// solid when `solid`, else a surface model (a shell and no solid).
fn faceted(corners: &[Point3], faces: &[FaceRings], solid: bool) -> BRep<NodeId> {
    let mut brep = BRep::default();
    let verts: Vec<_> = corners
        .iter()
        .map(|&position| brep.add_vertex(Vertex { position }))
        .collect();
    let mut edges = HashMap::new();
    let mut ids = Vec::new();
    for (outer, holes) in faces {
        let mut bounds = Vec::new();
        for (index, ring) in core::iter::once(outer).chain(holes).enumerate() {
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
        closed: solid,
    });
    if solid {
        brep.add_solid(Solid {
            outer: shell,
            voids: Vec::new(),
        });
    }
    brep
}

/// Both prisms of the layout close with volume area x depth, and the B-rep
/// reports its planar faces exact.
fn check_prism(outer: Vec<Point2>, holes: Vec<Vec<Point2>>, what: &str) {
    let holes: Vec<Vec<Point2>> = holes
        .into_iter()
        .map(|mut h| {
            if twice_ring_area(&h) > 0.0 {
                h.reverse();
            }
            h
        })
        .collect();
    let twice = twice_ring_area(&outer) + holes.iter().map(|h| twice_ring_area(h)).sum::<f64>();
    let (corners, faces) = prism(&outer, &holes);
    let mesh = compile_node(GeometryNode::PolygonMesh(authored(&corners, &faces)))
        .unwrap_or_else(|e| panic!("{what}: authored prism refused: {e:?}"));
    let volume = volume_properties(&mesh, Tolerance::METRE)
        .unwrap_or_else(|e| panic!("{what}: authored prism is not closed: {e:?}"))
        .signed_volume;
    close(
        volume,
        twice / 2.0 * DEPTH,
        &format!("{what}: authored volume"),
    );

    let mut b = GeometryGraphBuilder::new();
    let root = b
        .push(GeometryNode::BRep(faceted(&corners, &faces, true)))
        .unwrap();
    let graph = b.finish(vec![root]).unwrap();
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options())
        .unwrap_or_else(|e| panic!("{what}: B-rep prism refused: {e:?}"));
    let volume = volume_properties(&outcome.mesh, Tolerance::METRE)
        .unwrap_or_else(|e| panic!("{what}: B-rep prism is not closed: {e:?}"))
        .signed_volume;
    close(
        volume,
        twice / 2.0 * DEPTH,
        &format!("{what}: B-rep volume"),
    );
    assert!(
        report
            .contributions
            .iter()
            .all(|c| c.path == DeviationPath::BRepFace && c.bound == DeviationBound::Proven(0.0)),
        "{what}: {report:?}"
    );
}

#[test]
fn two_holes_in_one_band_close_both_prisms() {
    // The #253 profile: earcut left a T-junction along the band.
    check_prism(
        rect(-2.0, -2.0, 2.0, 2.0),
        vec![rect(-1.5, -0.5, -0.5, 0.5), rect(0.5, -0.5, 1.5, 0.5)],
        "side by side",
    );
    check_prism(
        rect(-3.0, -2.0, 3.0, 2.0),
        vec![
            rect(-2.5, -0.5, -1.5, 0.5),
            rect(-0.5, -0.5, 0.5, 0.5),
            rect(1.5, -0.5, 2.5, 0.5),
        ],
        "three on one ray",
    );
}

/// A small deterministic generator, as in construct's `profile_holes.rs`.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// The #253 hole shapes: rectangle, triangle, L, diamond, rectangle with
/// collinear midpoints, in a `w x h` box at the origin.
fn shape(kind: u64, w: f64, h: f64) -> Vec<Point2> {
    match kind {
        0 => rect(0.0, 0.0, w, h),
        1 => vec![p(0.0, 0.0), p(w, 0.0), p(w / 2.0, h)],
        2 => vec![
            p(0.0, 0.0),
            p(w, 0.0),
            p(w, h / 2.0),
            p(w / 2.0, h / 2.0),
            p(w / 2.0, h),
            p(0.0, h),
        ],
        3 => vec![
            p(w / 2.0, 0.0),
            p(w, h / 2.0),
            p(w / 2.0, h),
            p(0.0, h / 2.0),
        ],
        _ => vec![
            p(0.0, 0.0),
            p(w / 2.0, 0.0),
            p(w, 0.0),
            p(w, h / 2.0),
            p(w, h),
            p(w / 2.0, h),
            p(0.0, h),
            p(0.0, h / 2.0),
        ],
    }
}

#[test]
fn random_hole_layouts_close_both_prisms() {
    // #253's on-grid property, through the mesh compiler's planar paths:
    // holes sharing rows, columns and edge lines all the time.
    let mut rng = Lcg(0x260);
    for round in 0..60 {
        let q = 0.25;
        let mut boxes: Vec<(f64, f64, f64, f64)> = Vec::new();
        let target = 1 + rng.below(9) as usize;
        for _ in 0..target * 8 {
            if boxes.len() == target {
                break;
            }
            let w = (2 + rng.below(16)) as f64 * q;
            let h = (2 + rng.below(16)) as f64 * q;
            let x = (1 + rng.below(63)) as f64 * q - 8.0;
            let y = (1 + rng.below(63)) as f64 * q - 8.0;
            if x + w > 8.0 - q || y + h > 8.0 - q {
                continue;
            }
            let clear = boxes.iter().all(|&(bx, by, bw, bh)| {
                x >= bx + bw + q || bx >= x + w + q || y >= by + bh + q || by >= y + h + q
            });
            if clear {
                boxes.push((x, y, w, h));
            }
        }
        let holes = boxes
            .iter()
            .map(|&(x, y, w, h)| {
                shape(rng.below(5), w, h)
                    .into_iter()
                    .map(|c| p(c.x + x, c.y + y))
                    .collect()
            })
            .collect();
        let outer: Vec<Point2> = if round % 3 == 0 {
            shape(4, 16.0, 16.0)
                .into_iter()
                .map(|c| p(c.x - 8.0, c.y - 8.0))
                .collect()
        } else {
            rect(-8.0, -8.0, 8.0, 8.0)
        };
        check_prism(outer, holes, &format!("round {round}"));
    }
}

/// One face, authored and as a single-face B-rep surface, with exactly
/// the given twice-area.
fn check_patch(outer: Vec<Point2>, holes: Vec<Vec<Point2>>, twice: f64, what: &str) {
    let rings: Vec<&Vec<Point2>> = core::iter::once(&outer).chain(&holes).collect();
    let mut corners = Vec::new();
    let mut face: FaceRings = (Vec::new(), Vec::new());
    for (r, ring) in rings.iter().enumerate() {
        let start = corners.len();
        corners.extend(ring.iter().map(|q| Point3::new(q.x, q.y, 0.0)));
        let indices: Vec<usize> = (start..corners.len()).collect();
        if r == 0 {
            face.0 = indices;
        } else {
            face.1.push(indices);
        }
    }
    let faces = vec![face];
    let mesh = compile_node(GeometryNode::PolygonMesh(authored(&corners, &faces)))
        .unwrap_or_else(|e| panic!("{what}: authored face refused: {e:?}"));
    assert_eq!(twice_area(&mesh), twice, "{what}: authored area");
    for t in mesh.indices.chunks_exact(3) {
        let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
        assert!((b - a).cross(c - a).z > 0.0, "{what}: {t:?} faces down");
    }
    let mesh = compile_node(GeometryNode::BRep(faceted(&corners, &faces, false)))
        .unwrap_or_else(|e| panic!("{what}: B-rep face refused: {e:?}"));
    assert_eq!(twice_area(&mesh), twice, "{what}: B-rep area");
}

#[test]
fn pinched_faces_triangulate_as_surface_patches() {
    // Two footprints touching at a corner, one ring through it twice.
    let union = vec![
        p(0.0, 0.0),
        p(2.0, 0.0),
        p(2.0, 1.0),
        p(4.0, 1.0),
        p(4.0, 2.0),
        p(2.0, 2.0),
        p(2.0, 1.0),
        p(0.0, 1.0),
    ];
    check_patch(union, vec![], 8.0, "rooms at a corner");
    // A hole touching the outer ring at its corner, one inside its edge,
    // and two holes touching each other at a corner.
    let square = || rect(0.0, 0.0, 4.0, 4.0);
    check_patch(
        square(),
        vec![vec![p(0.0, 0.0), p(1.0, 2.0), p(2.0, 1.0)]],
        29.0,
        "hole at the outer corner",
    );
    check_patch(
        square(),
        vec![vec![p(2.0, 0.0), p(1.0, 1.0), p(3.0, 1.0)]],
        30.0,
        "hole on the outer edge",
    );
    check_patch(
        square(),
        vec![rect(1.0, 1.0, 2.0, 2.0), rect(2.0, 2.0, 3.0, 3.0)],
        28.0,
        "holes at a corner",
    );
}

#[test]
fn a_doubled_corner_inside_a_ring_is_dropped() {
    // Corner 2 repeats corner 1 exactly; it is left out and every
    // triangle still names the authored corners it was built from.
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(1.0, 2.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
    ];
    let mesh = compile_node(GeometryNode::PolygonMesh(PolygonMesh {
        positions,
        faces: vec![PolygonFace {
            outer: (0..7).collect(),
            holes: Vec::new(),
        }],
    }))
    .expect("a doubled corner is dropped, not refused");
    assert_eq!(twice_area(&mesh), 6.0);
    assert_eq!(
        mesh.indices.len(),
        12,
        "an L of six corners: four triangles"
    );
    assert!(!mesh.indices.contains(&2), "{:?}", mesh.indices);
    for t in mesh.indices.chunks_exact(3) {
        let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
        assert!((b - a).cross(c - a).z > 0.0, "{t:?} faces down");
    }
}

fn line(from: Point2, to: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}

fn polygon(points: &[Point2]) -> Contour {
    Contour::new(
        (0..points.len())
            .map(|i| line(points[i], points[(i + 1) % points.len()]))
            .collect(),
    )
}

#[test]
fn extruding_the_same_rings_is_refused_by_name() {
    // A solid over a pinched profile would share one wall edge between
    // four faces: the mesh extrusion refuses it.
    let mut b = GeometryGraphBuilder::new();
    let profile = b
        .push(GeometryNode::Profile(Profile::Contour(ContourProfile {
            outer: polygon(&rect(0.0, 0.0, 4.0, 4.0)),
            holes: vec![polygon(&[p(0.0, 0.0), p(1.0, 2.0), p(2.0, 1.0)])],
        })))
        .unwrap();
    let root = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: DEPTH,
        }))
        .unwrap();
    let graph = b.finish(vec![root]).unwrap();
    let error = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(&graph, root, &options())
        .expect_err("a pinched solid is refused");
    assert!(
        matches!(&error, GeomError::InvalidInput(m) if m.contains("touches")),
        "{error:?}"
    );
}
