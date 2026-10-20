"""Mutation probe for grown soups and the overlay's memory budget (#338,
`crates/algorithms/planar/overlay/src/exact_overlay.rs`,
`src/exact_arc/arrangement.rs` and `src/budget.rs`).

Each mutant weakens one step -- the chain reduction group by group (the
fallback, the groups, a group's simplicity check), the pairs of later twins
left unasked, the meetings at an edge's own ends left unrecorded, the
point equality read off boxes, and each charge of the budget and its
refusal -- and must turn a test red.

Equivalent mutants, deliberately not listed:
- joining groups by any shared edge rather than an edge given each way:
  larger groups cancel the same edges, and a group whose cycles cross
  falls back to its rings as given, so only the cost changes;
- leaving out, on the later of two edges, a point at an end of the earlier
  one: no scene was found where it changes the subdivision. Such a point
  is either where the earlier ring leaves the later edge's line, which the
  leaving edge's own crossing reports, or inside a stretch the earlier
  ring carries, whose pieces the earlier ring keeps (the case the other
  way round, a point at the later edge's end, is listed);
- the charges' constants: `tests/grown_soup.rs` measures that they bound
  the heap with room to spare, so a smaller constant shows only where its
  room is used up (the charges themselves, removed, are listed).
"""
import probe_lib

R = "crates/algorithms/planar/overlay/src/exact_overlay.rs"
A = "crates/algorithms/planar/overlay/src/exact_arc/arrangement.rs"
P = "crates/algorithms/planar/overlay/src/exact_arc/point.rs"
B = "crates/algorithms/planar/overlay/src/budget.rs"
L = "crates/algorithms/planar/overlay/src/lib.rs"
UNIT = ["-p", "axiolid-overlay", "--lib", "exact_overlay::"]
GROWN = ["-p", "axiolid-overlay", "--test", "grown_soup"]
SNAP = ["-p", "axiolid-overlay", "--test", "soup_snapshots"]
ALL = [UNIT, GROWN, SNAP]
ARR = ["-p", "axiolid-overlay", "--test", "arrangement"]

MUTANTS = [
    ('a crossing cycle keeps every ring as given', R,
     '        return cycles;\n    }\n    grouped(rings)',
     '        return cycles;\n    }\n    rings', ALL),
    ('groups not joined by edges given each way', R,
     '            if forward && back {',
     '            if false {', ALL),
    ("a group's crossing cycles kept", R,
     '            Some(cycles) if cycles.iter().all(|cycle| simple(cycle)) => out.extend(cycles),',
     '            Some(cycles) => out.extend(cycles),', ALL),
    ('a first twin left unasked too', A,
     '            if twin[index] == twin[near] && twin[index] != index {',
     '            if twin[index] == twin[near] {', [SNAP, GROWN]),
    ('a later twin asked about no edge', A,
     '            if twin[index] == twin[near] && twin[index] != index {',
     '            if twin[index] != index && repeated[near] {', [SNAP, GROWN]),
    ("a meeting at the other edge's end left out", A,
     '                if !edge.plainly_ends_at(x) {',
     '                if !theirs.plainly_ends_at(x) {', [ARR, SNAP, GROWN]),
    ('points with exact boxes plainly the same', P,
     '        || (exact(a) && exact(b) && a.bx.lo() == b.bx.lo() && a.by.lo() == b.by.lo())',
     '        || (exact(a) && exact(b) && a.bx.lo() == b.bx.lo())', [ARR, SNAP, GROWN]),
    ('a budget never refuses', B,
     '            Some(budget) if self.used > budget => Err(OverlayRefusal::BudgetExceeded {',
     '            Some(budget) if false => Err(OverlayRefusal::BudgetExceeded {', GROWN),
    ('a budget refuses at its own size', B,
     '            Some(budget) if self.used > budget => Err(OverlayRefusal::BudgetExceeded {',
     '            Some(budget) if self.used >= budget => Err(OverlayRefusal::BudgetExceeded {',
     GROWN),
    ('no budget refuses', L,
     '        &mut Charge::within(execution.memory_budget_bytes()),\n    )\n}\n\n/// [`overlay_with`], charged',
     '        &mut Charge::within(Some(execution.memory_budget_bytes().unwrap_or(0))),\n    )\n}\n\n/// [`overlay_with`], charged',
     GROWN),
    ('the operands not charged', L,
     '    charge.add(vertices(subject) + vertices(clip), budget::PER_VERTEX)?;',
     '    charge.add(0, budget::PER_VERTEX)?;', GROWN),
    ("a soup's rings not charged", L,
     '    charge.add(\n        rings.iter().map(|ring| ring.points.len()).sum(),\n        budget::PER_VERTEX,\n    )?;',
     '    charge.add(0, budget::PER_VERTEX)?;', GROWN),
    ('the edges not charged', A,
     '        rings.iter().map(|ring| ring.vertices.len()).sum(),\n        budget::PER_EDGE,',
     '        0,\n        budget::PER_EDGE,', GROWN),
    ('the meetings not charged', A,
     '            charge.add(recorded, budget::PER_MEETING)?;',
     '            charge.add(0 * recorded, budget::PER_MEETING)?;', GROWN),
    ('the pieces not charged', A,
     '    charge.add(count, budget::PER_PIECE)?;',
     '    charge.add(0 * count, budget::PER_PIECE)?;', GROWN),
]

probe_lib.main(MUTANTS, timeout=900)
