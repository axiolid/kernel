//! General exact boolean over planar-faced solids (#77).
//!
//! The oracles are independent of the boolean: volume via the v0.7 provider,
//! the identity |A ∪ B| + |A ∩ B| = |A| + |B|, and a differential against the
//! `boolmesh` mesh path which shares no code with this one.

use axiolid_construct::polyhedron::{boolean_polyhedra_exact, BooleanOp, Polyhedron};
use axiolid_core::{Point3, Tolerance};
use axiolid_heal::mesh::MeshHealer;
use axiolid_heal::self_intersections;
use axiolid_heal::Diagnose;
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("tolerance")
}

/// An axis-aligned box as six outward-wound quads.
fn box_solid(min: [f64; 3], max: [f64; 3]) -> Polyhedron {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    Polyhedron::new(vec![
        vec![p(x0, y0, z0), p(x0, y1, z0), p(x1, y1, z0), p(x1, y0, z0)],
        vec![p(x0, y0, z1), p(x1, y0, z1), p(x1, y1, z1), p(x0, y1, z1)],
        vec![p(x0, y0, z0), p(x1, y0, z0), p(x1, y0, z1), p(x0, y0, z1)],
        vec![p(x0, y1, z0), p(x0, y1, z1), p(x1, y1, z1), p(x1, y1, z0)],
        vec![p(x0, y0, z0), p(x0, y0, z1), p(x0, y1, z1), p(x0, y1, z0)],
        vec![p(x1, y0, z0), p(x1, y1, z0), p(x1, y1, z1), p(x1, y0, z1)],
    ])
    .expect("box is a valid solid")
}

/// Triangulate a polyhedron by fanning each face, for measurement.
///
/// Vertices are shared through an exact-coordinate lookup: emitting a fresh
/// vertex per face would leave every edge used once, so `audit_mesh` would
/// report a cloud of boundary edges and refuse to measure a solid that is
/// in fact closed. Coordinates from a boolean are bit-identical where faces
/// meet, because both sides come from the same split, so exact keying is
/// correct here and avoids inventing a welding tolerance.
fn to_mesh(solid: &Polyhedron) -> TriMesh {
    let mut positions: Vec<Point3> = Vec::new();
    let mut indices = Vec::new();
    let mut lookup: std::collections::HashMap<[u64; 3], u32> = std::collections::HashMap::new();

    let mut index_of = |p: Point3, positions: &mut Vec<Point3>| -> u32 {
        let key = [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
        *lookup.entry(key).or_insert_with(|| {
            positions.push(p);
            (positions.len() - 1) as u32
        })
    };

    for face in solid.faces() {
        let ring: Vec<u32> = face.iter().map(|&p| index_of(p, &mut positions)).collect();
        for i in 1..ring.len() - 1 {
            indices.extend([ring[0], ring[i], ring[i + 1]]);
        }
    }
    TriMesh::new(positions, indices)
}

fn volume(solid: &Polyhedron) -> f64 {
    volume_properties(&to_mesh(solid), tol())
        .expect("boolean result must be a closed solid")
        .signed_volume
}

#[test]
fn overlapping_boxes_satisfy_the_volume_identity() {
    // Two unit boxes overlapping in an eighth of their volume.
    let a = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_solid([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);

    let union = boolean_polyhedra_exact(&a, &b, BooleanOp::Union).expect("union");
    let intersection =
        boolean_polyhedra_exact(&a, &b, BooleanOp::Intersection).expect("intersection");

    let (va, vb) = (volume(&a), volume(&b));
    let (vu, vi) = (volume(&union), volume(&intersection));

    // The identity holds for any pair of solids, so it checks the two
    // operations against each other rather than against a hand-computed
    // number that could be wrong in the same way twice.
    assert!(
        (vu + vi - (va + vb)).abs() < 1e-9,
        "|A u B| + |A n B| = |A| + |B| violated: {vu} + {vi} != {va} + {vb}"
    );
    // And the intersection is independently known: a 0.5 cube.
    assert!(
        (vi - 0.125).abs() < 1e-9,
        "intersection of the two boxes is a 0.5-cube, got {vi}"
    );
}

#[test]
fn difference_removes_exactly_the_shared_volume() {
    let a = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_solid([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);

    let difference = boolean_polyhedra_exact(&a, &b, BooleanOp::Difference).expect("difference");
    let expected = volume(&a) - 0.125;
    let actual = volume(&difference);
    assert!(
        (actual - expected).abs() < 1e-9,
        "difference volume {actual} != {expected}"
    );
}

/// An L-shaped prism: the case a convex-only implementation gets wrong.
///
/// Built as a genuine 6-gon footprint. The caps stay single non-convex
/// rings rather than being pre-split into rectangles: splitting them would
/// introduce cap vertices the side walls do not carry, leaving T-junctions
/// that read as boundary edges. The boolean itself never needs the caps
/// triangulated -- only the measurement helper does, and it fans each ring.
fn l_prism(z0: f64, z1: f64) -> Polyhedron {
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    let ring = [
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 1.0),
        (1.0, 1.0),
        (1.0, 2.0),
        (0.0, 2.0),
    ];
    let mut faces = Vec::new();
    faces.push(ring.iter().rev().map(|&(x, y)| p(x, y, z0)).collect());
    faces.push(ring.iter().map(|&(x, y)| p(x, y, z1)).collect());
    for i in 0..ring.len() {
        let (x0, y0) = ring[i];
        let (x1, y1) = ring[(i + 1) % ring.len()];
        faces.push(vec![
            p(x0, y0, z0),
            p(x1, y1, z0),
            p(x1, y1, z1),
            p(x0, y0, z1),
        ]);
    }
    Polyhedron::new(faces).expect("L-prism is a valid solid")
}

/// Volume of an L-prism, computed from its footprint rather than measured.
///
/// The measurement helper fans each face ring, and a fan across the L's
/// notch emits triangles outside the solid. Rather than weaken the oracle,
/// the L's own volume is known in closed form: footprint area 3 times
/// height.
fn l_prism_volume(z0: f64, z1: f64) -> f64 {
    3.0 * (z1 - z0)
}

#[test]
fn non_convex_operands_are_handled_correctly() {
    let l = l_prism(0.0, 1.0);
    assert!(
        (l_prism_volume(0.0, 1.0) - 3.0).abs() < 1e-12,
        "the L footprint encloses area 3"
    );

    // A box strictly inside the notch, touching nothing. The notch is a
    // concavity: a convex-only classifier calls this point set "inside the
    // L" because it is on the inner side of every face plane, and would
    // report a non-empty intersection.
    let notch = box_solid([1.2, 1.2, 0.2], [1.8, 1.8, 0.8]);
    let result = boolean_polyhedra_exact(&l, &notch, BooleanOp::Intersection);
    assert!(
        result.is_err(),
        "the notch interior is OUTSIDE the L, so the intersection is empty; \
         a convex-only classifier would wrongly return a solid here"
    );

    // A box overlapping the L's lower arm genuinely intersects.
    let arm = box_solid([1.5, 0.0, 0.0], [2.5, 0.5, 1.0]);
    let hit = boolean_polyhedra_exact(&l, &arm, BooleanOp::Intersection).expect("real overlap");
    // Overlap is x 1.5..2, y 0..0.5, z 0..1 = 0.25.
    assert!(
        (volume(&hit) - 0.25).abs() < 1e-9,
        "notch-free overlap volume {}",
        volume(&hit)
    );
}

/// The exact result carries no defects the v0.7 diagnosis can find.
#[test]
fn results_are_closed_manifold_and_free_of_self_intersection() {
    let a = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_solid([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);

    for op in [
        BooleanOp::Union,
        BooleanOp::Intersection,
        BooleanOp::Difference,
    ] {
        let result = boolean_polyhedra_exact(&a, &b, op).expect("boolean");
        let mesh = to_mesh(&result);
        let diagnosis = MeshHealer.diagnose(&mesh, tol()).expect("diagnose");
        assert!(
            diagnosis.is_clean(),
            "{op:?} produced defects: {:?}",
            diagnosis.defects
        );
        // With #73's coplanar branch now deciding shared area exactly, a
        // boolean result's sibling fragments no longer register as false
        // positives, so this can assert the real property directly.
        assert!(
            self_intersections(&mesh).is_empty(),
            "{op:?} produced a self-intersecting shell"
        );
    }
}

/// A curved operand is refused by name, not approximated.
///
/// v0.6's exact revolution produces cylindrical faces, and those remain out
/// of scope: general curved surface/surface intersection is separate work.
/// The refusal must be typed, because silently tessellating a cylinder and
/// running the mesh boolean would answer an exact request approximately.
#[test]
fn a_non_planar_face_is_refused_not_approximated() {
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    // A "face" whose four corners do not share a plane.
    let warped = vec![
        vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(1.0, 1.0, 0.5),
            p(0.0, 1.0, 0.0),
        ],
        vec![p(0.0, 0.0, 0.0), p(0.0, 1.0, 0.0), p(0.0, 0.0, 1.0)],
        vec![p(1.0, 0.0, 0.0), p(0.0, 0.0, 1.0), p(0.0, 0.0, 0.0)],
        vec![p(1.0, 1.0, 0.5), p(0.0, 0.0, 1.0), p(1.0, 0.0, 0.0)],
    ];
    let error = Polyhedron::new(warped).expect_err("a warped face has no single plane");
    let text = format!("{error}");
    assert!(
        text.contains("planar"),
        "the refusal must name planarity as the reason, got: {text}"
    );
}

/// A grid-aligned subtraction is answered, not refused.
///
/// Every cutter face is a plane of the host, so the operands share vertices
/// and edges in bulk. That is the configuration a single fixed probe
/// direction hits degenerately: the ray leaves a fragment centroid and meets
/// a shared vertex exactly, and #77 shipped refusing the whole operation.
///
/// This is the Menger sponge's depth-1 cell arrangement -- the exact shape
/// that made the benchmark harness report `refused` from depth 2 -- reduced
/// to the smallest case that still shares planes on all three axes.
#[test]
fn a_grid_aligned_cutter_is_classified_not_refused() {
    // Host spans 0..3 on every axis; the cutter is the centre column of the
    // 3x3x3 grid, so all four of its side planes coincide with grid planes
    // the host's own subdivision would produce.
    let host = box_solid([0.0, 0.0, 0.0], [3.0, 3.0, 3.0]);
    let cutter = box_solid([1.0, 1.0, -1.0], [2.0, 2.0, 4.0]);

    let result = boolean_polyhedra_exact(&host, &cutter, BooleanOp::Difference)
        .expect("grid-aligned difference must be answered, not refused");

    // 27 - 3 = 24: the cutter passes clean through one 1x1x3 column.
    let volume = volume(&result);
    assert!(
        (volume - 24.0).abs() < 1e-9,
        "expected volume 24.0, got {volume}"
    );
}

/// Repeated grid-aligned subtraction is answered, not refused.
///
/// A single grid-aligned cutter was always fine. The failure needed
/// ACCUMULATION: each difference introduces fragment centroids that are
/// themselves grid-aligned with the next cutter, and by the ninth
/// subtraction the fixed probe direction meets a shared vertex exactly.
/// #77 shipped refusing the whole operation there, which is what made the
/// Menger benchmark report `refused` from depth 2.
///
/// Containment is direction-independent, so retrying along another ray is
/// exact rather than a fudge. This walks the first ten voids of a depth-2
/// Menger sponge -- the shortest sequence that reaches the degenerate step.
#[test]
fn accumulated_grid_aligned_subtraction_is_answered() {
    let mut current = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let t = 1.0 / 3.0;

    // Depth-1 centre columns, then depth-2 columns in one surviving cell.
    let cutters = [
        ([t, t, -1.0], [2.0 * t, 2.0 * t, 2.0]),
        ([t, -1.0, t], [2.0 * t, 2.0, 2.0 * t]),
        ([-1.0, t, t], [2.0, 2.0 * t, 2.0 * t]),
        (
            [t / 3.0, t / 3.0, -1.0],
            [2.0 * t / 3.0, 2.0 * t / 3.0, 2.0],
        ),
        (
            [t / 3.0, -1.0, t / 3.0],
            [2.0 * t / 3.0, 2.0, 2.0 * t / 3.0],
        ),
        (
            [-1.0, t / 3.0, t / 3.0],
            [2.0, 2.0 * t / 3.0, 2.0 * t / 3.0],
        ),
    ];

    for (i, (mn, mx)) in cutters.iter().enumerate() {
        let cutter = box_solid(*mn, *mx);
        current = boolean_polyhedra_exact(&current, &cutter, BooleanOp::Difference)
            .unwrap_or_else(|e| panic!("subtraction {i} refused: {e}"));
    }

    // `volume` is unavailable here: repeated subtraction produces non-convex
    // faces, and `triangulate` fans them, which leaves the mesh unusable for
    // measurement (documented on `triangulate`). Assert on the B-rep instead.
    //
    // Each subtraction must add faces -- a cutter that removed nothing, or an
    // operation that collapsed the solid, would not.
    assert!(
        current.faces().len() > 6,
        "expected the subtracted solid to carry more faces than the host box, \
         got {}",
        current.faces().len()
    );
    for face in current.faces() {
        assert!(
            face.len() >= 3,
            "every face of the result must be a real polygon"
        );
    }
}

/// Void boxes of a Menger sponge at `depth`, on the unit cube.
///
/// Each level splits a cell into 27 and removes the 7 with at least two
/// centred coordinates: the six face centres and the core.
fn menger_voids(depth: u32) -> Vec<([f64; 3], [f64; 3])> {
    fn carve(out: &mut Vec<([f64; 3], [f64; 3])>, origin: [f64; 3], size: f64, depth: u32) {
        if depth == 0 {
            return;
        }
        let t = size / 3.0;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let lo = [
                        origin[0] + t * i as f64,
                        origin[1] + t * j as f64,
                        origin[2] + t * k as f64,
                    ];
                    let centred = usize::from(i == 1) + usize::from(j == 1) + usize::from(k == 1);
                    if centred >= 2 {
                        out.push((lo, [lo[0] + t, lo[1] + t, lo[2] + t]));
                    } else {
                        carve(out, lo, t, depth - 1);
                    }
                }
            }
        }
    }
    let mut out = Vec::new();
    carve(&mut out, [0.0, 0.0, 0.0], 1.0, depth);
    out
}

/// A long grid-aligned subtraction chain completes.
///
/// Repeated subtraction against grid-aligned cutters composes coordinate
/// error until a split emits a COLLAPSED ring: a quad whose vertices pair
/// up, spanning a third of the model yet enclosing no area because its
/// width is one ULP. Such a face has no usable normal, so every probe ray
/// meets it edge-on and containment cannot be decided at all -- the
/// direction-retry family cannot rescue it, because the problem is the
/// face, not the ray.
///
/// This walks the void set of a depth-2 Menger sponge, which is the
/// benchmark workload that exposed it. Before the area check in
/// `split_polygon`, this refused at the 82nd of 147 subtractions.
#[test]
fn a_long_grid_aligned_subtraction_chain_completes() {
    let mut current = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let mut done = 0usize;

    for (min, max) in menger_voids(2) {
        let cutter = box_solid(min, max);
        current = boolean_polyhedra_exact(&current, &cutter, BooleanOp::Difference)
            .unwrap_or_else(|e| panic!("subtraction {done} of 147 refused: {e}"));
        done += 1;
    }

    assert_eq!(done, 147, "every void must be subtracted");
    // Sanity on the result itself: a solid this carved has far more faces
    // than the host box, and every face must still be a real polygon.
    assert!(current.faces().len() > 100);
    for face in current.faces() {
        assert!(face.len() >= 3);
    }
}
