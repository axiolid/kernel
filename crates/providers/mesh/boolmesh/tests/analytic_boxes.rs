//! The opt-in analytic box-subtraction path.
//!
//! Two obligations are tested here, and the second matters more than the first:
//!
//! 1. when the analytic path answers, it agrees with the general solver;
//! 2. when it cannot be exact, it declines instead of guessing.
//!
//! A fast path that returns a plausible wrong answer is worse than no fast path,
//! so the refusal cases are not an afterthought.

mod support;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
use support::{boxx, volume};

/// Generous ceiling: these fixtures induce grids far below it, so a `None` in
/// these tests means a genuine refusal, never a budget trip.
const CELLS: usize = 1 << 20;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

/// A wall with `n` evenly spaced openings, the shape this path exists for.
fn wall_with_openings(n: usize) -> (TriMesh, Vec<TriMesh>) {
    let wall = boxx(0.0, 0.0, 0.0, 10.0, 0.3, 3.0, 0.0);
    let tools = (0..n)
        .map(|i| {
            let x = -4.0 + (i as f64) * (8.0 / (n.max(2) - 1) as f64);
            boxx(x, 0.0, 0.9, 0.5, 1.0, 1.2, 0.0)
        })
        .collect();
    (wall, tools)
}

/// A vertex position as exact bits, so coincident corners compare equal without
/// float tolerance games. Positions here are produced by the same arithmetic on
/// both sides of a shared face, so bitwise equality is the right test.
type PosKey = (u64, u64, u64);

/// Every directed edge must appear exactly once, and its reverse exactly once.
///
/// This is the oracle signed volume cannot provide. A crack, a T-junction, or a
/// duplicated face leaves an unpaired edge, but can still integrate to the right
/// volume through cancellation -- so a volume-only check would pass a mesh that
/// is not a closed solid.
fn edge_pairing_defects(mesh: &TriMesh) -> usize {
    use std::collections::HashMap;
    // Key on POSITION, not vertex index. Index-keyed pairing silently passes a
    // mesh whose coincident corners are duplicate vertices -- exactly the defect
    // (cracks, unmerged interior faces) this check exists to find.
    let key = |i: u32| -> PosKey {
        let p = mesh.positions[i as usize];
        (p.x.to_bits(), p.y.to_bits(), p.z.to_bits())
    };
    let mut counts: HashMap<(PosKey, PosKey), i32> = HashMap::new();
    for tri in mesh.indices.as_chunks::<3>().0 {
        for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            let (ka, kb) = (key(a), key(b));
            let (k, delta) = if ka < kb {
                ((ka, kb), 1)
            } else {
                ((kb, ka), -1)
            };
            *counts.entry(k).or_insert(0) += delta;
        }
    }
    // Zero means each undirected edge was traversed once in each direction.
    counts.values().filter(|v| **v != 0).count()
}

/// Count faces that appear more than once at the same location.
///
/// Two coincident triangles with opposite winding are the signature of an
/// interior face pair that should have cancelled. They are invisible to both
/// signed volume (they sum to zero) and edge pairing (each edge still balances),
/// so this is the only one of the three oracles that detects them. A "closed"
/// mesh full of internal walls is not the solid the caller asked for.
fn duplicate_faces(mesh: &TriMesh) -> usize {
    use std::collections::HashMap;
    let mut seen: HashMap<[PosKey; 3], usize> = HashMap::new();
    for tri in mesh.indices.as_chunks::<3>().0 {
        let mut key: Vec<PosKey> = tri
            .iter()
            .map(|&i| {
                let p = mesh.positions[i as usize];
                (p.x.to_bits(), p.y.to_bits(), p.z.to_bits())
            })
            .collect();
        // Sort so winding does not distinguish a face from its opposite twin.
        key.sort_unstable();
        *seen.entry([key[0], key[1], key[2]]).or_insert(0) += 1;
    }
    seen.values().filter(|c| **c > 1).count()
}

#[test]
fn agrees_with_the_general_solver_on_a_wall_with_openings() {
    let provider = BoolmeshBoolean::new();
    let opts = options();

    for n in [1usize, 2, 4, 8] {
        let (wall, tools) = wall_with_openings(n);

        let analytic = provider
            .subtract_boxes_analytic(&wall, &tools, &opts, CELLS)
            .expect("analytic path must not error on boxes")
            .unwrap_or_else(|| panic!("analytic path declined {n} axis-aligned openings"));

        let exact = provider
            .subtract_many(&wall, &tools, &opts)
            .expect("general solver");

        // Volumes agree: the two paths computed the same solid.
        let (va, ve) = (volume(&analytic.mesh), volume(&exact.mesh));
        assert!(
            (va - ve).abs() <= 1e-9 * ve.abs().max(1.0),
            "n={n}: analytic volume {va} != general volume {ve}"
        );

        // ...and the analytic result is genuinely closed, which volume alone
        // would not establish.
        assert_eq!(
            edge_pairing_defects(&analytic.mesh),
            0,
            "n={n}: analytic result is not a closed surface"
        );

        // The caller can tell which machinery ran.
        assert!(analytic.evidence.analytic_path);
        assert!(!exact.evidence.analytic_path);
    }
}

#[test]
fn output_is_byte_identical_across_repeated_runs() {
    // The general solver is nondeterministic here (upstream `boolmesh` dedups
    // vertices through a randomly-seeded `HashMap`). The analytic path uses an
    // ordered map specifically so it is not, and that claim is worth a test.
    let provider = BoolmeshBoolean::new();
    let opts = options();
    let (wall, tools) = wall_with_openings(6);

    let first = provider
        .subtract_boxes_analytic(&wall, &tools, &opts, CELLS)
        .expect("analytic")
        .expect("recognised");

    for _ in 0..8 {
        let again = provider
            .subtract_boxes_analytic(&wall, &tools, &opts, CELLS)
            .expect("analytic")
            .expect("recognised");
        assert_eq!(
            first.mesh.indices, again.mesh.indices,
            "index order drifted"
        );
        assert_eq!(
            first.mesh.positions.len(),
            again.mesh.positions.len(),
            "vertex count drifted"
        );
        for (a, b) in first.mesh.positions.iter().zip(&again.mesh.positions) {
            assert_eq!((a.x, a.y, a.z), (b.x, b.y, b.z), "vertex position drifted");
        }
    }
}

#[test]
fn declines_a_rotated_subject() {
    // 30 degrees about Z. The bounding box is still a box, which is exactly the
    // trap: a detector that accepted on extent alone would cut the wrong solid.
    let provider = BoolmeshBoolean::new();
    let wall = boxx(0.0, 0.0, 0.0, 10.0, 0.3, 3.0, std::f64::consts::FRAC_PI_6);
    let tools = vec![boxx(0.0, 0.0, 0.9, 0.5, 1.0, 1.2, 0.0)];
    assert!(provider
        .subtract_boxes_analytic(&wall, &tools, &options(), CELLS)
        .expect("refusal is not an error")
        .is_none());
}

#[test]
fn declines_a_rotated_tool() {
    let provider = BoolmeshBoolean::new();
    let wall = boxx(0.0, 0.0, 0.0, 10.0, 0.3, 3.0, 0.0);
    let tools = vec![boxx(0.0, 0.0, 0.9, 0.5, 1.0, 1.2, 0.4)];
    assert!(provider
        .subtract_boxes_analytic(&wall, &tools, &options(), CELLS)
        .expect("refusal is not an error")
        .is_none());
}

#[test]
fn declines_a_non_box_subject() {
    // A tetrahedron: 4 triangles, so it fails on structure rather than on any
    // coordinate coincidence.
    let provider = BoolmeshBoolean::new();
    let tet = TriMesh::new(
        vec![
            axiolid_core::Point3::new(0.0, 0.0, 0.0),
            axiolid_core::Point3::new(2.0, 0.0, 0.0),
            axiolid_core::Point3::new(0.0, 2.0, 0.0),
            axiolid_core::Point3::new(0.0, 0.0, 2.0),
        ],
        vec![0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3],
    );
    let tools = vec![boxx(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0)];
    assert!(provider
        .subtract_boxes_analytic(&tet, &tools, &options(), CELLS)
        .expect("refusal is not an error")
        .is_none());
}

#[test]
fn declines_a_box_with_one_corner_pulled_off_the_lattice() {
    // 12 triangles and every face still assigned to a plane, but one vertex is
    // displaced. This isolates the lattice check: a detector that trusted the
    // bounding box would accept this and cut a solid that does not exist.
    let provider = BoolmeshBoolean::new();
    let mut dented = boxx(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 0.0);
    let p = dented.positions[0];
    dented.positions[0] = axiolid_core::Point3::new(p.x + 0.9, p.y, p.z);
    let tools = vec![boxx(0.0, 0.0, 0.0, 1.0, 1.0, 8.0, 0.0)];
    assert!(provider
        .subtract_boxes_analytic(&dented, &tools, &options(), CELLS)
        .expect("refusal is not an error")
        .is_none());
}

#[test]
fn result_is_a_closed_surface_when_the_grid_is_subdivided() {
    // A single centred opening splits the wall into a 3x3x3 grid, so the solid
    // cells share many interior faces. Those faces must cancel: if the emitter
    // wrote them, the surface would be self-intersecting internally while still
    // integrating to the correct volume. Edge pairing is what detects that;
    // volume alone cannot.
    let provider = BoolmeshBoolean::new();
    let wall = boxx(0.0, 0.0, 0.0, 6.0, 1.0, 6.0, 0.0);
    let tools = vec![boxx(0.0, 0.0, 2.0, 2.0, 4.0, 2.0, 0.0)];

    let out = provider
        .subtract_boxes_analytic(&wall, &tools, &options(), CELLS)
        .expect("analytic")
        .expect("recognised");

    assert_eq!(
        edge_pairing_defects(&out.mesh),
        0,
        "interior faces were emitted: surface is not closed"
    );

    // Edge pairing and volume are both blind to a cancelled interior face pair:
    // the edges still balance and the volumes still sum to zero. Only a
    // coincident-face check sees them, so the emitter's cancellation rule is
    // pinned here.
    assert_eq!(
        duplicate_faces(&out.mesh),
        0,
        "interior faces were emitted: duplicate coincident triangles"
    );

    // The exact triangle count for this fixture: a 3x3x3 grid with the centre
    // column removed has a known boundary, and any extra face changes it.
    assert_eq!(
        out.mesh.indices.len() / 3,
        64,
        "unexpected triangle count for a single centred opening"
    );

    // Independent confirmation that the cut really happened, so the test cannot
    // be satisfied by returning the uncut wall. The opening is deeper than the
    // wall is thick, so the removed volume is the CLIPPED overlap (the wall's
    // 1.0 thickness), not the tool's own volume.
    let removed = 2.0 * 1.0 * 2.0;
    let expected = volume(&wall) - removed;
    let got = volume(&out.mesh);
    assert!(
        (got - expected).abs() <= 1e-9 * expected,
        "volume {got} != expected {expected}"
    );
}

#[test]
fn declines_a_box_with_a_face_swapped_for_a_duplicate() {
    // 12 triangles and every vertex on the lattice, but one face plane carries
    // 4 triangles while another carries none: an open box with a doubled face.
    // This isolates the per-plane count check -- the count and lattice checks
    // both pass. Cutting an open surface would produce nonsense.
    let provider = BoolmeshBoolean::new();
    let mut open = boxx(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 0.0);
    // Overwrite the first face's two triangles with copies of the second face's.
    let (a, b) = (open.indices[6], open.indices[7]);
    let c = open.indices[8];
    let (d, e, f) = (open.indices[9], open.indices[10], open.indices[11]);
    open.indices[0] = a;
    open.indices[1] = b;
    open.indices[2] = c;
    open.indices[3] = d;
    open.indices[4] = e;
    open.indices[5] = f;
    let tools = vec![boxx(0.0, 0.0, 0.0, 1.0, 1.0, 8.0, 0.0)];
    assert!(provider
        .subtract_boxes_analytic(&open, &tools, &options(), CELLS)
        .expect("refusal is not an error")
        .is_none());
}

/// 12 valid box triangles plus two stray indices: 38, not a multiple of 3.
///
/// `chunks_exact(3)` silently drops a trailing partial chunk, so the face-plane
/// loop sees a perfect box and would accept. Only the explicit index-count check
/// notices the malformed buffer. Verified by mutation: disabling that check
/// makes this input accepted.
fn box_with_trailing_partial_triangle() -> TriMesh {
    let base = boxx(0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 0.0);
    let mut indices = base.indices.clone();
    indices.push(0);
    indices.push(1);
    TriMesh::new(base.positions, indices)
}

#[test]
fn declines_a_malformed_index_buffer() {
    let provider = BoolmeshBoolean::new();
    let tools = vec![boxx(0.0, 0.0, 0.0, 1.0, 1.0, 8.0, 0.0)];
    assert!(provider
        .subtract_boxes_analytic(
            &box_with_trailing_partial_triangle(),
            &tools,
            &options(),
            CELLS
        )
        .expect("refusal is not an error")
        .is_none());
}

/// A box carrying one extra vertex that no triangle references.
///
/// This is the only input where the lattice rule and the face-plane rule
/// disagree: the plane check inspects referenced vertices only, so an unused
/// off-lattice position is invisible to it. Accepting such a mesh would mean
/// the analytic path silently ignores geometry the caller supplied -- and if
/// that vertex is later referenced by a caller-side edit, the "recognised box"
/// was never a box at all.
fn box_with_stray_vertex() -> TriMesh {
    let base = boxx(0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 0.0);
    let mut positions = base.positions.clone();
    positions.push(axiolid_core::Point3::new(0.37, 0.41, 0.53));
    TriMesh::new(positions, base.indices)
}

#[test]
fn declines_a_box_with_an_off_lattice_stray_vertex() {
    let provider = BoolmeshBoolean::new();
    let tools = vec![boxx(0.0, 0.0, 0.0, 1.0, 1.0, 8.0, 0.0)];
    assert!(provider
        .subtract_boxes_analytic(&box_with_stray_vertex(), &tools, &options(), CELLS)
        .expect("refusal is not an error")
        .is_none());
}

#[test]
fn declines_when_the_grid_would_exceed_the_budget() {
    // The grid is O(n^3) in cutter count; a caller-supplied ceiling is what
    // stops a pathological input from allocating without bound.
    let provider = BoolmeshBoolean::new();
    let (wall, tools) = wall_with_openings(8);
    assert!(provider
        .subtract_boxes_analytic(&wall, &tools, &options(), 4)
        .expect("refusal is not an error")
        .is_none());
}

#[test]
fn declines_when_no_tool_meets_the_subject() {
    // Nothing to cut is not a cut of nothing: hand it back so the caller keeps
    // one code path for "unchanged", rather than inventing a result here.
    let provider = BoolmeshBoolean::new();
    let wall = boxx(0.0, 0.0, 0.0, 10.0, 0.3, 3.0, 0.0);
    let far = vec![boxx(100.0, 100.0, 100.0, 1.0, 1.0, 1.0, 0.0)];
    assert!(provider
        .subtract_boxes_analytic(&wall, &far, &options(), CELLS)
        .expect("refusal is not an error")
        .is_none());
    assert!(provider
        .subtract_boxes_analytic(&wall, &[], &options(), CELLS)
        .expect("refusal is not an error")
        .is_none());
}

#[test]
fn result_survives_the_orientation_contract() {
    // The analytic construction is held to the same outward-orientation rule as
    // the general path: an inside-out solid poisons every later subtraction.
    let provider = BoolmeshBoolean::new();
    let (wall, tools) = wall_with_openings(3);
    let out = provider
        .subtract_boxes_analytic(&wall, &tools, &options(), CELLS)
        .expect("analytic")
        .expect("recognised");
    assert!(volume(&out.mesh) > 0.0, "result is inside-out");

    // A wall with three openings is still one connected piece.
    assert_eq!(out.evidence.output_components, 1);

    // And it really did remove material.
    let solid = provider
        .boolean(&wall, &wall, BooleanOperator::Union, &options())
        .expect("union");
    assert!(volume(&out.mesh) < volume(&solid.mesh));
}
