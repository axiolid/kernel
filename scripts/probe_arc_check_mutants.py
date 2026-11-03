"""Mutation probe for the exact arc-polygon check: each mutant must be caught."""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "crates/algorithms/planar/overlay/src"

MUTANTS = [
    ("adjacency exemption dropped", "arc_check.rs",
     "            let next = j == i + 1;",
     "            let next = false;"),
    ("every point of an adjacent pair exempted", "arc_check.rs",
     "!(next && same_point(x, &edges[i].p1))",
     "!next"),
    ("closing pair not exempted", "arc_check.rs",
     "            let closing = i == 0 && j == n - 1;",
     "            let closing = false;"),
    ("closing pair not checked", "arc_check.rs",
     "        for &j in hits.iter().filter(|&&j| j > i) {",
     "        for &j in hits.iter().filter(|&&j| j > i && !(i == 0 && j == n - 1)) {"),
    ("box broad phase reads false", "arc_check.rs",
     "|b| b.overlaps(&edges[i].bounds)",
     "|_b| false"),
    ("holes not nested", "arc_check.rs",
     "    if simple {\n        nesting(&rings, tolerance, &mut check);",
     "    if false {\n        nesting(&rings, tolerance, &mut check);"),
    ("holes overlapping not reported", "arc_check.rs",
     "crowded.insert([only.min(other) - 1, only.max(other) - 1]);",
     "let _ = (only, other);"),
    ("outer containment not read", "arc_check.rs",
     "                if !edge.inside_left(0) {",
     "                if false {"),
    ("overlap not reported", "arc_check.rs",
     "overlaps.insert([a, b]);",
     "if false {\n                        overlaps.insert([a, b]);\n                    }"),
    ("contacts reported as defects", "arc_check.rs",
     "    pub fn is_valid(&self) -> bool {\n        self.defects.is_empty()",
     "    pub fn is_valid(&self) -> bool {\n        self.defects.is_empty() && self.contacts.is_empty()"),
    ("nesting run on a non-simple ring", "arc_check.rs",
     "                at,\n            });\n            simple = false;",
     "                at,\n            });"),
    ("zero-area ring skips the crossing test", "arc_check.rs",
     "            if !zero_area {",
     "            if true {"),
    ("outer-only pieces treated as holes", "arc_check.rs",
     "            [] | [0] => {}",
     "            [] => {}"),
    ("a ring meets itself at a vertex", "arc_check.rs",
     "                if !at_vertex[v].contains(&ring) {",
     "                if true {"),
    ("closing contact cycle not detected", "arc_check.rs",
     "            if a == b {\n                closing.push(ring);",
     "            if false {\n                closing.push(ring);"),
    ("all contact points merged into one node", "arc_check.rs",
     "        let point = ring_count + vertex;",
     "        let point = ring_count;"),
    ("disconnected interior judged beside other defects", "arc_check.rs",
     "    let disconnected = if overlaps.is_empty() && outside.is_empty() && crowded.is_empty() {",
     "    let disconnected = if true {"),
]

probe_lib.main(MUTANTS, tests=[["--release", "-p", "axiolid-overlay",
         "--test", "arc_check"]], root=SRC, timeout=900)
