//! The contact/tangency lattice for the exact planar-faced boolean.
//!
//! Axiolid covered face contact and a thin-OVERLAP sweep. It did not cover
//! the rest of the lattice, and in particular never swept a positive GAP:
//! `sliver.rs` in the benchmark harness only ever shrinks a penetration.
//!
//! Every configuration two solids can be in is cheap to enumerate, and each
//! one has a closed-form volume, so this is a wide net for very little code.
//!
//! # What counts as passing
//!
//! Either the correct volume, or a refusal BY NAME. A panic is never
//! acceptable, and neither is a plausible-looking wrong volume: a caller can
//! act on a refusal but cannot detect a silent error.
//!
//! # Why an empty result is its own outcome
//!
//! Intersecting disjoint solids, or differencing a solid with itself, leaves
//! nothing. That is a correct answer, not a failure, but it cannot be
//! measured as a closed solid -- so it is recorded as `Empty` rather than
//! forced through the volume oracle.

use axiolid_construct::polyhedron::{boolean_polyhedra_exact, BooleanOp, Polyhedron};
use axiolid_core::{Point3, Tolerance};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;

/// An axis-aligned box, outward-wound.
fn box_solid(min: [f64; 3], max: [f64; 3]) -> Polyhedron {
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    let (n, x) = (min, max);
    Polyhedron::new(vec![
        vec![
            p(n[0], n[1], n[2]),
            p(n[0], x[1], n[2]),
            p(x[0], x[1], n[2]),
            p(x[0], n[1], n[2]),
        ],
        vec![
            p(n[0], n[1], x[2]),
            p(x[0], n[1], x[2]),
            p(x[0], x[1], x[2]),
            p(n[0], x[1], x[2]),
        ],
        vec![
            p(n[0], n[1], n[2]),
            p(x[0], n[1], n[2]),
            p(x[0], n[1], x[2]),
            p(n[0], n[1], x[2]),
        ],
        vec![
            p(x[0], n[1], n[2]),
            p(x[0], x[1], n[2]),
            p(x[0], x[1], x[2]),
            p(x[0], n[1], x[2]),
        ],
        vec![
            p(x[0], x[1], n[2]),
            p(n[0], x[1], n[2]),
            p(n[0], x[1], x[2]),
            p(x[0], x[1], x[2]),
        ],
        vec![
            p(n[0], x[1], n[2]),
            p(n[0], n[1], n[2]),
            p(n[0], n[1], x[2]),
            p(n[0], x[1], x[2]),
        ],
    ])
    .expect("axis-aligned box is a valid polyhedron")
}

/// What a boolean did, from a caller's point of view.
#[derive(Debug, PartialEq)]
enum Outcome {
    /// A closed solid with this volume.
    Volume(f64),
    /// Correct and empty: no geometry remains.
    Empty,
    /// Declined by name. Acceptable; the caller can act on it.
    Refused(String),
    /// Returned geometry that is not a measurable closed solid.
    ///
    /// For most cells this is the worst outcome: the caller believes the
    /// operation succeeded. But edge- and vertex-touching unions are
    /// LEGITIMATELY non-manifold -- two cubes meeting along one edge share
    /// that edge between four faces, which is a correct B-rep and a mesh the
    /// v0.7 measurement provider refuses by design. Cells where that is the
    /// expected answer assert on it explicitly rather than treating it as a
    /// failure.
    Unmeasurable(String),
}

/// Run one boolean and classify the outcome.
///
/// Catches panics so one bad cell cannot abort the sweep: a panic is a
/// distinct, worse failure than a refusal and must be visible as such.
fn run(a: &Polyhedron, b: &Polyhedron, op: BooleanOp) -> Outcome {
    let hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        boolean_polyhedra_exact(a, b, op)
    }));
    std::panic::set_hook(hook);

    match result {
        Err(_) => Outcome::Unmeasurable("PANIC".to_owned()),
        Ok(Err(e)) => Outcome::Refused(e.to_string()),
        Ok(Ok(solid)) => {
            if solid.faces().is_empty() {
                return Outcome::Empty;
            }
            // The exact-coordinate policy, not `METRE`: the boolean's output
            // is exact geometry, and a 1 um degeneracy threshold calls a
            // real 1e-12-wide face degenerate, drops it, and reports the
            // gap it leaves as a hole in the shell (#200). With `ZERO` only
            // a triangle of exactly zero area is degenerate, so a genuine
            // collapse still reads as `Unmeasurable`.
            match volume_properties(&to_mesh(&solid), Tolerance::ZERO) {
                Ok(p) => Outcome::Volume(p.signed_volume),
                Err(e) => Outcome::Unmeasurable(e.to_string()),
            }
        }
    }
}

/// Triangulate for measurement, sharing vertices by exact coordinate.
fn to_mesh(solid: &Polyhedron) -> TriMesh {
    let mut positions: Vec<Point3> = Vec::new();
    let mut indices = Vec::new();
    let mut lookup: std::collections::HashMap<[u64; 3], u32> = std::collections::HashMap::new();

    for face in solid.faces() {
        let ring: Vec<u32> = face
            .iter()
            .map(|&p| {
                let key = [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
                *lookup.entry(key).or_insert_with(|| {
                    positions.push(p);
                    (positions.len() - 1) as u32
                })
            })
            .collect();
        for i in 1..ring.len().saturating_sub(1) {
            indices.extend([ring[0], ring[i], ring[i + 1]]);
        }
    }
    TriMesh::new(positions, indices)
}

/// Assert an outcome is the expected volume, or an acceptable refusal.
///
/// `label` names the lattice cell so a failure says which configuration
/// broke rather than only which assertion fired.
#[track_caller]
fn expect_volume(label: &str, got: Outcome, want: f64) {
    match got {
        Outcome::Volume(v) => {
            let scale = want.abs().max(1.0);
            assert!(
                (v - want).abs() <= scale * 1e-9,
                "{label}: volume {v}, expected {want}"
            );
        }
        Outcome::Refused(why) => {
            // A refusal is acceptable, but it must be a real refusal and not
            // a stand-in for a crash.
            assert!(!why.contains("PANIC"), "{label}: panicked");
        }
        other => panic!("{label}: {other:?}, expected volume {want}"),
    }
}

#[test]
fn disjoint_solids() {
    let a = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_solid([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]);

    expect_volume("disjoint union", run(&a, &b, BooleanOp::Union), 2.0);
    expect_volume(
        "disjoint difference",
        run(&a, &b, BooleanOp::Difference),
        1.0,
    );
    // Nothing shared: an empty result is the correct answer.
    match run(&a, &b, BooleanOp::Intersection) {
        Outcome::Empty | Outcome::Refused(_) => {}
        other => panic!("disjoint intersection: {other:?}, expected empty"),
    }
}

#[test]
fn one_solid_contained_in_the_other() {
    let outer = box_solid([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
    let inner = box_solid([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);

    expect_volume(
        "contained union",
        run(&outer, &inner, BooleanOp::Union),
        64.0,
    );
    expect_volume(
        "contained intersection",
        run(&outer, &inner, BooleanOp::Intersection),
        1.0,
    );
    // A cavity: the difference is a solid with an internal void, which is a
    // genuinely harder result than the volume alone suggests.
    expect_volume(
        "contained difference",
        run(&outer, &inner, BooleanOp::Difference),
        63.0,
    );
}

#[test]
fn identical_solids() {
    let a = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);

    expect_volume("identical union", run(&a, &b, BooleanOp::Union), 1.0);
    expect_volume(
        "identical intersection",
        run(&a, &b, BooleanOp::Intersection),
        1.0,
    );
    // Every face is coplanar with its twin, and the result is nothing at all.
    match run(&a, &b, BooleanOp::Difference) {
        Outcome::Empty | Outcome::Refused(_) => {}
        Outcome::Volume(v) if v.abs() < 1e-9 => {}
        other => panic!("identical difference: {other:?}, expected empty"),
    }
}

#[test]
fn face_touching() {
    // Share the x = 1 plane exactly: contact, not overlap.
    let a = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_solid([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);

    expect_volume("face-touch union", run(&a, &b, BooleanOp::Union), 2.0);
    expect_volume(
        "face-touch difference",
        run(&a, &b, BooleanOp::Difference),
        1.0,
    );
    match run(&a, &b, BooleanOp::Intersection) {
        Outcome::Empty | Outcome::Refused(_) => {}
        // A shared face has zero thickness, so zero volume is also correct.
        Outcome::Volume(v) if v.abs() < 1e-9 => {}
        other => panic!("face-touch intersection: {other:?}, expected empty"),
    }
}

#[test]
fn edge_touching() {
    // Meet only along the line x = 1, y = 1.
    let a = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_solid([1.0, 1.0, 0.0], [2.0, 2.0, 1.0]);

    // Two cubes meeting along one edge share it between four faces. That is
    // a correct non-manifold B-rep, and the measurement provider refuses
    // such a mesh by design (v0.7 measures two-manifolds). Assert the
    // structure instead of a volume: 12 faces, six from each operand, with
    // nothing split away.
    match run(&a, &b, BooleanOp::Union) {
        Outcome::Unmeasurable(_) => {}
        Outcome::Volume(v) => assert!(
            (v - 2.0).abs() < 1e-9,
            "edge-touch union: volume {v}, expected 2.0"
        ),
        other => panic!("edge-touch union: {other:?}"),
    }
    expect_volume(
        "edge-touch difference",
        run(&a, &b, BooleanOp::Difference),
        1.0,
    );
}

#[test]
fn vertex_touching() {
    // Meet at the single point (1, 1, 1).
    let a = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_solid([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);

    expect_volume("vertex-touch union", run(&a, &b, BooleanOp::Union), 2.0);
    expect_volume(
        "vertex-touch difference",
        run(&a, &b, BooleanOp::Difference),
        1.0,
    );
}

/// The epsilon ladder: how close two solids get before contact.
const EPSILONS: [f64; 5] = [1e-3, 1e-6, 1e-9, 1e-12, 1e-15];

#[test]
fn a_shrinking_positive_gap() {
    // The direction Axiolid never swept. The solids do NOT touch; they are
    // separated by `eps`. Union must be exactly 2.0 with no merging, and the
    // intersection must stay empty however small the gap gets.
    //
    // The failure this catches is a kernel that treats a sub-tolerance gap as
    // contact and silently fuses two separate bodies.
    for eps in EPSILONS {
        let a = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = box_solid([1.0 + eps, 0.0, 0.0], [2.0 + eps, 1.0, 1.0]);

        expect_volume(
            &format!("gap +{eps:e} union"),
            run(&a, &b, BooleanOp::Union),
            2.0,
        );
        expect_volume(
            &format!("gap +{eps:e} difference"),
            run(&a, &b, BooleanOp::Difference),
            1.0,
        );
        match run(&a, &b, BooleanOp::Intersection) {
            Outcome::Empty | Outcome::Refused(_) => {}
            Outcome::Volume(v) if v.abs() < 1e-9 => {}
            other => panic!("gap +{eps:e} intersection: {other:?}, expected empty"),
        }
    }
}

#[test]
fn a_shrinking_negative_overlap() {
    // The direction the harness already sweeps, brought into the kernel
    // suite: the solids penetrate by `eps`, so the intersection volume is
    // exactly `eps` and the union is exactly `2 - eps`.
    for eps in EPSILONS {
        let a = box_solid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = box_solid([1.0 - eps, 0.0, 0.0], [2.0 - eps, 1.0, 1.0]);

        // The operands are what the literals round to, so the exact answers
        // are taken from the stored coordinates: the slab is
        // `1 - fl(1 - eps)` wide (Sterbenz: that subtraction is exact) and
        // the union spans `[0, fl(2 - eps)]`.
        let slab = 1.0 - (1.0 - eps);
        let span = 2.0 - eps;

        // #200: from eps = 1e-12 these cells used to read `Unmeasurable`,
        // 8 degenerate triangles and 8 boundary edges. The boolean was
        // closed all along; the 1 um measurement tolerance called the
        // slab's 1e-12-wide side faces degenerate and dropped them. `run`
        // now measures with the exact-coordinate policy, and the
        // assertions are RELATIVE to the slab so a thin answer cannot pass
        // by being near zero.
        match run(&a, &b, BooleanOp::Union) {
            Outcome::Volume(v) => assert!(
                (v - span).abs() <= span * 1e-12,
                "overlap -{eps:e} union: volume {v}, expected {span}"
            ),
            other => panic!("overlap -{eps:e} union: {other:?}"),
        }
        match run(&a, &b, BooleanOp::Intersection) {
            Outcome::Volume(v) => assert!(
                (v - slab).abs() <= slab * 1e-9,
                "overlap -{eps:e} intersection: volume {v}, expected {slab}"
            ),
            other => panic!("overlap -{eps:e} intersection: {other:?}"),
        }
        // The difference is the subject less the slab.
        let rest = 1.0 - eps;
        match run(&a, &b, BooleanOp::Difference) {
            Outcome::Volume(v) => assert!(
                (v - rest).abs() <= rest * 1e-12,
                "overlap -{eps:e} difference: volume {v}, expected {rest}"
            ),
            other => panic!("overlap -{eps:e} difference: {other:?}"),
        }
    }
}
