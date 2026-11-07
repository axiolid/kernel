"""Mutation probe for the mesh-to-exact direction of a measured boolean (#329).

Each mutant weakens one decision of `reverse::mesh_to_exact` or its wiring in
`of_boolean`; the `deviation::boolean` unit tests or `tests/boolean_deviation.rs`
must fail on it. A build error also counts as killed.

Not probed, because no fixture separates them (each drops a term of the
error the bound is widened by, and on every fixture that term is below 1e-13,
under the search's own rounding margin of about 1e-11, so no bound a test can
read changes): the recentring of the mesh in the bound, the witness error in
the final widening, and the weld's move in the witness error. They are the
argument in `reverse.rs`'s docs (steps 4 to 7), not something a value pins;
the witness's own error is pinned by
`the_witness_error_holds_the_slack_and_the_rounding_of_every_face`. Every
other soundness decision is killed.
"""
import probe_lib

B = "crates/execution/compile/src/deviation/boolean.rs"
R = "crates/execution/compile/src/deviation/boolean/reverse.rs"
UNIT = ["-p", "axiolid-mesh-compile", "--lib", "deviation::boolean"]
INTEGRATION = ["-p", "axiolid-mesh-compile", "--test", "boolean_deviation"]
TESTS = [UNIT, INTEGRATION]

MUTANTS = [
    ("the other direction never measured", B,
     "            if faces.iter().all(|f| f.domain.affine) {",
     "            if false && faces.iter().all(|f| f.domain.affine) {"),
    ("a curved result measured in reverse", B,
     "            if faces.iter().all(|f| f.domain.affine) {",
     "            if faces.iter().any(|f| f.domain.affine) {"),
    ("the other direction not widened by the readings and the floor", B,
     "            covered(bound, report.rounding_floor(), reading_distance(&report)),",
     "            bound,"),
    ("the reverse search settles at the whole budget", R,
     "    let settle = if widening < budget {\n        budget - widening\n    } else {\n        budget\n    };",
     "    let settle = budget;"),
    ("the readings shown to the search as well as widened (inflated)", R,
     "settle, flattened, max_cells)",
     "settle, widening, max_cells)"),
    ("the mesh and the witness not recentred", R,
     "    let origin = (lo + hi) * 0.5;",
     "    let origin = Vec3::ZERO * (lo + hi);"),
    ("the witness's rounding term dropped", R,
     "        let rounding = 9.0 * Scalar::EPSILON * (rel.abs().max_element() + reach);",
     "        let rounding = 0.0 * (rel.abs().max_element() + reach);"),
    ("the slack term dropped from the witness error", R,
     "        let off = 2.0 * domain.slack * (du.length() + dv.length());",
     "        let off = 0.0 * domain.slack * (du.length() + dv.length());"),
    ("a mesh cell inside across the hypotenuse", R,
     "        if x.1 + y.1 <= 1.0 {",
     "        if x.0 + y.0 <= 1.0 {"),
    ("the hypotenuse clip keeps the far side", R,
     "            Vec2::new(-1.0, -1.0) / core::f64::consts::SQRT_2,",
     "            Vec2::new(1.0, 1.0) / core::f64::consts::SQRT_2,"),
    ("witness faces not turned outward", R,
     "            indices.extend(if face.outward { [a, b, c] } else { [a, c, b] });",
     "            indices.extend([a, b, c]);"),
    ("the witness not welded", R,
     "    let (mesh, moved) = weld(&positions, &indices, scale);",
     "    let (mesh, moved) = (TriMesh::new(positions.clone(), indices.clone()), 0.0);"),
    ("holes dropped from the witness", R,
     "            .filter(|&(k, r)| k != outer && r.len() >= 3)",
     "            .filter(|&(k, r)| k != outer && r.len() >= 3 && false)"),
    ("a hole taken as the outer ring", B,
     "        domain.outer = face.bounds.iter().position(|bound| bound.outer);",
     "        domain.outer = face.bounds.iter().position(|bound| !bound.outer);"),
    ("a face without an outer ring refused instead of skipped", R,
     "    let outer = outer.filter(|&k| rings.get(k).is_some_and(|r| r.len() >= 3))?;",
     "    let outer = outer.filter(|&k| rings.get(k).is_some_and(|r| r.len() >= 3)).unwrap();"),
]

probe_lib.main(MUTANTS, tests=TESTS, timeout=1800, reject_uncompiled=True)
