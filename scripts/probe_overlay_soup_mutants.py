"""Mutation probe for the soup speed-ups in the exact overlay (#218): each mutant must be caught.

Every shortcut keeps a fallback, and each mutant below removes or bends
one: a stop taken as a simple crossing when it is not, a segment taken as
carried by its own ring alone when another boundary runs along it, the
rings carried across a stop updated the wrong way, a lazy point's box no
longer around the point, a filter deciding what it cannot. The tests run
with debug assertions off, so the cross-checks are off and a kill comes from a
wrong result, not from an assertion comparing the shortcut with the count.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "crates/algorithms/planar/overlay/src"
TESTS = [["--config", "profile.dev.debug-assertions=false", "-p", "axiolid-overlay"]]
ARR = "exact_arc/arrangement.rs"
POINT = "exact_arc/point.rs"

MUTANTS = [
    # Simple stops: one crossing, interior to both segments, met once.
    ("every crossing taken as simple", ARR,
     [("add(&mut cuts[index], x, simple.then_some(near));",
       "add(&mut cuts[index], x, Some(near));"),
      ("add(&mut cuts[near], x, simple.then_some(index));",
       "add(&mut cuts[near], x, Some(index));")], TESTS[0]),
    # Not listed: a check that the point is at no end of either segment.
    # It is implied: such a point is met again (the end is a stop, or the
    # next edge of that ring meets it too), which clears the tag; checked,
    # removing it survives.
    ("a stop met twice keeps its crossing", ARR,
     "            Some(held) => held.1 = None,",
     "            Some(_) => {}"),
    # Lone segments: carried by their own ring alone.
    ("overlapping segments taken as lone", ARR,
     "            if straight && xs.len() > 1 {",
     "            if false {"),
    ("arcs taken as lone", ARR,
     "        let lone = !edge.is_arc() && !overlapped[index];",
     "        let lone = !overlapped[index];"),
    # Rings carried across a stop.
    ("ring left but kept", ARR,
     "                            inside.remove(at);\n",
     ""),
    ("ring entered but not added", ARR,
     "                        Err(at) => inside.insert(at, other),",
     "                        Err(_) => {}"),
    ("rings carried across every stop", ARR,
     "            let carried = match (&before, w[0].1) {\n"
     "                (Some(previous), Some(crossing)) => {",
     "            let carried = match (&before, w[0].1.or(Some(0))) {\n"
     "                (Some(previous), Some(crossing)) => {"),
    # Identity of copies.
    ("box-overlapping points taken as copies", POINT,
     "    if Arc::ptr_eq(&a.exact, &b.exact) {",
     "    if true {"),
    ("every point interned as the first", ARR,
     "        if let Some(&index) = by_id.get(&x.id()) {",
     "        if let Some(&index) = by_id.values().next() {"),
    # Lazy points: their boxes must hold them.
    ("crossing box along y from x", POINT,
     "(i(p0.x).add(&t.mul(&dx)), i(p0.y).add(&t.mul(&dy_)))",
     "(i(p0.x).add(&t.mul(&dx)), i(p0.y).add(&t.mul(&dx)))"),
    ("sample box misses the offset", POINT,
     "            a_.add(&Interval::point(u).mul(&Interval::point(b).sub(&a_)))",
     "            a_.add(&Interval::point(u).mul(&Interval::point(b)))"),
    # Filters: each decides only what its boxes prove.
    ("rounding filter puts the lower midpoint above", POINT,
     "                -0.5 * (r - lo)",
     "                0.5 * (r - lo)"),
    ("fixed tangents' cross filtered as their dot", POINT,
     "                u.bx.mul(&v.by).sub(&u.by.mul(&v.bx))",
     "                u.bx.mul(&v.bx).add(&u.by.mul(&v.by))"),
    ("reversed direction keeps its boxes", POINT,
     "            bx: self.bx.neg(),",
     "            bx: self.bx,"),
    ("equal x coordinates taken as one point", "exact_arc.rs",
     "    if a == b || c == a || c == b {",
     "    if a == b || c == a || c.x == b.x {"),
    # Buffers and the snap's grid.
    ("query keeps the previous hits", "exact_arc/boxes.rs",
     "        out.clear();\n        out.extend_from_slice(&self.always);",
     "        out.extend_from_slice(&self.always);"),
    ("snap grid misses its top row", "exact_overlay.rs",
     "            for cy in y0..=y1 {",
     "            for cy in y0..y1 {"),
    ("snap band in y without the tolerance", "exact_overlay.rs",
     "            let slack = tolerance + 1e-9",
     "            let slack = 1e-9"),
]

probe_lib.main(MUTANTS, tests=TESTS, root=SRC, timeout=1200)
