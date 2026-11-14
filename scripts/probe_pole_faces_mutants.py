"""Mutation probe for B-rep faces that reach a pole of their surface (#463, ADR 0112).

Each mutant weakens one decision of the tessellator's pole handling: the
pole predicate of `crates/execution/compile/src/brep/support.rs`, the pole
end of an edge and its bound along the meridians in `brep/inversion.rs`,
the chains, rows, refusals, fans and their demotion in `brep/pole.rs`, and
the uneven lattice of a ring with a pole row in `brep.rs`. The fixtures of
`tests/pole_faces.rs`, the edge inversion tests at a cone apex and a ruled
apex (`tests/edge_inversion.rs`, `tests/relation_edge_inversion.rs`),
`brep/exact_tests.rs` and the unit tests of `brep::pole` must fail on it. A
build error also counts as killed.

Not listed, as equivalent at test resolution:
- the one-span edge between two poles left finished: an edge from pole to
  pole on a sphere or a revolution is a meridian of half a turn, whose
  single span's midpoint misses any budget the fixtures set, so it is
  halved before the guard is read.
- a row narrower than the angular tolerance inserted anyway: every fixture
  with a pole junction has a row a whole cell wide.
- the box along the meridians not clamped at the pole a span ends at: it
  only widens the box, so the bound grows and stays a bound; on every
  fixture the pole is also an end of the support's range of `v` (a
  sphere's latitudes, a cone's apex, a profile's domain), which clamps it
  the same.
- the meridian term `dvv dv^2 / 8` left out of a span into a pole
  (`pole_span`), and the meridian spread left out of a span near one
  (`meridian_lens`): the argument needs them, no measurement shows them.
  They are zero on a cone and along a ruled surface's rulings. On tilted
  circles into a sphere's pole and on spans near and into the pole of the
  canopy's curved profile (`brep::inversion` unit tests) the farthest
  points measured are the edge's own: 2.21e-3 against the chord bound
  2.23e-3 into the canopy's pole, 6.7e-3 against `d + e` = 1.8e-2 near
  it, so either term left out still bounds every sample. Those tests kill
  the chord bound left out instead.

First run: 12 / 15. The survivors were the chord bound and the meridian
term of a span into a pole, and the fans' corner offsets: no fixture's
pole span dominated its face's bound, and every fixture's vertices lie on
the surface. Unit tests against dense sampling of a span's lens and of a
fan with lifted corners were added, and the meridian term is listed above.
"""
import probe_lib

S = "crates/execution/compile/src/brep/support.rs"
I = "crates/execution/compile/src/brep/inversion.rs"
P = "crates/execution/compile/src/brep/pole.rs"
B = "crates/execution/compile/src/brep.rs"
TESTS = [
    ["-p", "axiolid-mesh-compile", "--test", "pole_faces"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::"],
    ["-p", "axiolid-mesh-compile", "--test", "relation_edge_inversion"],
    ["-p", "axiolid-mesh-compile", "--test", "edge_inversion"],
]

MUTANTS = [
    # The pole predicate and the pole end.
    ("is_pole always false", S,
     "            None => false,\n            Some(PoleKind::Rotational) => self\n",
     "            _ if true => false,\n            None => false,\n            Some(PoleKind::Rotational) => self\n"),
    ("a pole end takes u = 0", I,
     [("        params[0].x = params[1].x;\n", "        params[0].x = 0.0;\n"),
      ("        params[n].x = params[n - 1].x;\n", "        params[n].x = 0.0;\n")],
     TESTS),
    ("pole() applied to interior samples", I,
     "            raw.push(inverter.invert(mid)?);\n",
     "            raw.push(match inverter.pole(mid)? {\n"
     "                Some(v) => Point2::new(0.0, v),\n"
     "                None => inverter.invert(mid)?,\n"
     "            });\n"),
    ("a span into a pole halved onto the pole", I,
     "                crowded |= !matches!(inverter.pole(mid), Ok(None));\n",
     "                crowded |= false;\n"),
    # The bound along the meridians (I6).
    ("pole_span without span.chord", I,
     "        let near = span.chord + offset;\n",
     "        let near = offset;\n"),
    ("a span across a pole bounded along the meridians", I,
     "        if clear.is_nan() || clear <= arc {\n",
     "        if false {\n"),
    # Chains, rows and refusals (I2 to I4).
    ("closure row reversed", P,
     "        if us.windows(2).all(|w| w[1] > w[0]) != (start.x > end.x) {\n",
     "        if us.windows(2).all(|w| w[1] > w[0]) == (start.x > end.x) {\n"),
    ("a full turn read from two distinct edges", P,
     "        if k != 2 || !one_seam(&chains[0], &chains[1]) {\n",
     "        if k != 2 {\n"),
    ("P2 replaced by the positive lune", P,
     "    if k >= 2 && !chains.iter().all(|chain| chain.anchored) {\n",
     "    if false {\n"),
    ("full-turn row set to 0", P,
     "        chains[1].shift_by((turn / period).round() * period);\n",
     "        chains[1].shift_by(0.0 * (turn / period).round() * period);\n"),
    ("row not welded to the pole vertex", P,
     "        shared.extend(core::iter::repeat_n(pole, row.len()));\n",
     "        shared.extend(core::iter::repeat_n(pole + 1, row.len()));\n"),
    # Fans (I7).
    ("fan Q without duu", P,
     "            0.125 * (bounds.duu * du * du + 2.0 * bounds.duv * du * dv + bounds.dvv * dv * dv)\n",
     "            0.125 * (0.0 * du * du + 2.0 * bounds.duv * du * dv + bounds.dvv * dv * dv)\n"),
    ("fan without corner offsets", P,
     "        let rest = across + offset;\n",
     "        let rest = across;\n"),
    ("demotion removed", P,
     "                Some(terms) if terms.rest <= budget => row.push(triangle),\n",
     "                Some(_) => row.push(triangle),\n"),
    # The uneven lattice of a ring with a pole row.
    ("a pole ring with uneven samples left to the clipper", B,
     "    if uniform.is_some() || !has_pole_row(boundary) {\n",
     "    if true {\n"),
]

probe_lib.main(MUTANTS, tests=TESTS, command=probe_lib.CARGO_TEST + ("--all-features",), timeout=1800)
