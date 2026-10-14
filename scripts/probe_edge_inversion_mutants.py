"""Mutation probe for curved faces whose edges carry only 3D curves (#310).

Each mutant weakens one decision of `crates/execution/compile/src/brep/
inversion.rs` or its caller in `brep.rs`: the refusal at a cone apex or
sphere pole, the unwrapping and its quarter-turn guard, the strict pass of
the certified B-spline inversion, the polish inside its enclosure, the
domain clamp, the check that an edge curve runs between its vertices, and
the terms the deviation report adds per span (the parameter chord's image,
the lens and the containment that sizes it).
`tests/edge_inversion.rs` or the module's unit tests must fail on it. A
build error also counts as killed.

Not listed, as equivalent at test resolution:
- the parameter chord's image left out of a span's term: the report adds
  the face's triangle bound to its trim terms, and that bound already
  covers the surface over every triangle edge, the chord included, so no
  sample can tell the image from it.
- the lens term left out: it is of second order in a width the edge's
  sampling keeps far below the budget, so dropping it moves no bound
  past a measured sample.
  The containment check that sizes it is listed.
- `align` and the B-spline parameter check: both act only on a periodic
  B-spline face, which no fixture has.

First run: 8 / 13. The survivors were the two mutants above, the
quarter-turn guard, the polish and the span terms; the guard and the
polish got unit tests, and the flat patch trimmed by a 3D circle (whose
bound is its trim terms alone) was added for the span terms. 11 / 11
killed after.
"""
import probe_lib

I = "crates/execution/compile/src/brep/inversion.rs"
B = "crates/execution/compile/src/brep.rs"
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--test", "edge_inversion"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::inversion"],
]

MUTANTS = [
    ("edges without a pcurve refused again", B,
     "            let Some(pcurve) = use_.pcurve else {\n",
     "            let Some(pcurve) = use_.pcurve else {\n"
     "                return Err(GeomError::Unsupported {\n"
     "                    backend: crate::BACKEND_ID,\n"
     "                    operation: axiolid_contracts::Operation::Tessellation,\n"
     "                });\n"
     "                #[allow(unreachable_code)]\n"),
    ("no refusal at an apex or pole", I,
     "        if radial.is_nan() || radial <= self.tolerance.linear() {",
     "        if radial.is_nan() {"),
    ("edge parameters not unwrapped", I,
     "            localized_parameter(params[i - 1].x, params[i].x, pu),",
     "            params[i].x + 0.0 * localized_parameter(params[i - 1].x, params[i].x, pu),"),
    ("quarter-turn guard reads half a period", I,
     "period.is_none_or(|p| d.abs() < 0.25 * p)",
     "period.is_none_or(|p| d.abs() < 0.5 * p)"),
    ("B-spline inverted at the model tolerance only", I,
     "        let outcome = match run(search.strict)? {",
     "        let outcome = match run(search.model)? {"),
    ("certified B-spline parameters not polished", I,
     "        for _ in 0..POLISH_STEPS {",
     "        for _ in 0..0 {"),
    ("polish not clamped to the domain", I,
     "        let (u_lo, u_hi) = (u_lo.max(u_domain.0), u_hi.min(u_domain.1));",
     "        let (u_lo, u_hi) = (u_lo, u_hi);"),
    ("edge curve need not start at its vertex", I,
     "        if gap.is_nan() || gap > linear {",
     "        if gap.is_nan() {"),
    ("span terms left out of the report", I,
     "            self.trim = self.trim.max(near + lens);",
     "            self.trim = self.trim.max(lens);"),
    ("lens containment not checked", I,
     "            if on_box >= sigma {",
     "            if true {"),
    ("sampling ignores the certified terms", I,
     "            if certified.is_none_or(|worst| worst <= budget) {",
     "            if certified.is_none_or(|worst| worst <= budget) || true {"),
]

probe_lib.main(MUTANTS, tests=TARGETS, command=probe_lib.CARGO_TEST + ("--all-features",), timeout=1800)
