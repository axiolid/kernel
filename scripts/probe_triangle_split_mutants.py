"""Mutation probe for the shared triangle split (#422).

`axiolid_mesh::split_triangle_at_side_corners`
(`crates/representations/discrete/mesh/src/split.rs`) is the one side-bitmask
ear clipping that `axiolid-mesh-compile` and `axiolid-heal` call. Each mutant
breaks one decision of it and must turn the mesh crate's `split` unit tests
red (every mix of 0..=4 corners per side must be a proper, equal-area,
same-winding cover; the piece order is fixed; a failed split leaves the
output as it was), or the consumers' T-junction and planar tests.

Not listed: the refusal when no ear exists. A ring of a triangle with corners
on its sides always has one (the exhaustive test shows it), so the branch is
a defensive refusal that no input reaches.

Runs on `scripts/probe_lib.py`: each mutant in its own process group under a
timeout that kills the group, the source always restored, and a missing or
ambiguous anchor aborts the probe. A mutant that does not compile counts as
BROKEN, not killed.
"""
import probe_lib

S = "crates/representations/discrete/mesh/src/split.rs"
MESH = [["-p", "axiolid-mesh", "--lib", "split"]]
HEAL = [["-p", "axiolid-heal", "--test", "split_t_junctions"]]
COMPILE = [["-p", "axiolid-mesh-compile", "--test", "authored_polygons"]]
ALL = MESH + HEAL + COMPILE
TIMEOUT = 1800

MUTANTS = [
    ('the triangle corners on the wrong side masks', S,
     "ring.push((t[i], (1 << i) | (1 << ((i + 2) % 3))));",
     "ring.push((t[i], (1 << i) | (1 << ((i + 1) % 3))));", MESH),
    ('a side corner shares every side', S,
     "ring.extend(sides[i].iter().map(|&c| (c, 1u8 << i)));",
     "ring.extend(sides[i].iter().map(|&c| (c, 0b111u8)));", MESH),
    ('collinearity needs the three masks to agree on two sides', S,
     "masks.fold(0b111u8, |acc, m| acc & m) != 0;",
     "masks.fold(0b111u8, |acc, m| acc & m).count_ones() > 1;", MESH),
    ('collinearity read as any shared side of two points', S,
     "masks.fold(0b111u8, |acc, m| acc & m) != 0;",
     "masks.fold(0u8, |acc, m| acc | m) != 0;", MESH),
    ('the ear leaves the ring unremoved', S,
     "        ring.remove(k);\n", "        ring.swap_remove(k);\n", MESH),
    ('the last ear taken instead of the first', S,
     "let ear = (0..n).find(|&k| {", "let ear = (0..n).rev().find(|&k| {", MESH),
    ('the ear emitted with reversed winding', S,
     "out.extend_from_slice(&[ring[prev].0, ring[k].0, ring[next].0]);",
     "out.extend_from_slice(&[ring[next].0, ring[k].0, ring[prev].0]);", ALL),
    ('the previous neighbour taken as the next', S,
     "let (prev, next) = ((k + n - 1) % n, (k + 1) % n);",
     "let (prev, next) = ((k + 1) % n, (k + 1) % n);", MESH),
    ('the final triangle dropped', S,
     "    out.extend_from_slice(&[ring[0].0, ring[1].0, ring[2].0]);\n    true",
     "    true", MESH),
    ('a side corner takes the place of a triangle corner', S,
     "ring.push((t[i], ", "ring.push((t[(i + 1) % 3], ", ALL),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
