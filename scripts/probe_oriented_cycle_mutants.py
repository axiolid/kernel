"""Mutation probe for volume over oriented edge cycles (#416, ADR 0099).

Each mutant breaks one rule of the oriented-cycle audit, of the volume
measure built on it, or of the T-junction split that feeds it, and must turn
a test of its crate red. All 21 are killed.
"""
import probe_lib

A = "crates/representations/discrete/mesh/src/audit.rs"
M = "crates/algorithms/query/measure/src/mesh.rs"
T = "crates/algorithms/repair/heal/src/tjunction.rs"
H = "crates/algorithms/repair/heal/src/mesh.rs"
C = "crates/algorithms/repair/heal/src/carry.rs"

MESH = ["-p", "axiolid-mesh"]
MEASURE = ["-p", "axiolid-measure", "--test", "cycle_volume"]
HEAL = ["-p", "axiolid-heal", "--test", "split_t_junctions"]

MUTANTS = [
    ('vec sink counts no imbalance', A, '            if winding != 0 {\n                summary.unbalanced += 1;', '            if false && winding != 0 {\n                summary.unbalanced += 1;', MESH),
    ('map sink counts only boundary', A, '.filter(|(_, &(_, winding))| winding != 0);', '.filter(|(_, &(count, winding))| count == 1 && winding != 0);', MESH),
    ('first unbalanced edge replaced by the last', A, 'summary.first_unbalanced.get_or_insert(key);', 'summary.first_unbalanced = Some(key);', MESH),
    ('cycle predicate drops the degenerate check', A, '            && self.health.degenerate_triangles == 0\n            && self.unbalanced_edges == 0', '            && self.unbalanced_edges == 0', MESH),
    ('cycle predicate drops the usable check', A, '        self.health.is_surface_usable()\n            && self.health.degenerate_triangles == 0\n            && self.unbalanced_edges == 0', '        self.health.degenerate_triangles == 0\n            && self.unbalanced_edges == 0', MESH),
    ('measure gates on two-manifold', M, '    if !audit.is_closed_oriented_cycle() {', '    if !audit.health.is_closed_two_manifold() {', MEASURE),
    ('measure does not gate at all', M, '    if !audit.is_closed_oriented_cycle() {', '    if false && !audit.is_closed_oriented_cycle() {', MEASURE),
    ('helper sums about the world origin', M, '    let base = mesh.position(0);\n    let (mut volume, mut weighted)', '    let base = Point3::ZERO;\n    let (mut volume, mut weighted)', MEASURE),
    ('centroid divides by three', M, 'weighted += (a + b + c) * (tetra / 4.0);', 'weighted += (a + b + c) * (tetra / 3.0);', MEASURE),
    ('split never fires', T, '                    && distance <= linear', '                    && distance < 0.0', HEAL),
    ('endpoint exclusion removed', T, '                    && (pm - pa).length() > linear\n                    && (pm - pb).length() > linear\n', '', HEAL),
    ('unbalanced filter removed', T, '        .filter(|&(_, &winding)| winding != 0)', '        .filter(|_| true)', HEAL),
    ('points sorted reversed', T, 'found.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));', 'found.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)));', HEAL),
    ('max deviation last instead of max', T, 'report.max_deviation = report.max_deviation.max(distance);', 'report.max_deviation = distance;', HEAL),
    ('orientation guard removed', T, 'normal.length_squared() > limit && normal.dot(parent) > 0.0', 'normal.length_squared() > limit', HEAL),
    ('degenerate guard removed', T, 'normal.length_squared() > limit && normal.dot(parent) > 0.0', 'normal.dot(parent) > 0.0', HEAL),
    ('corner-value decline removed', T, 'sides[s].is_empty() || carry::corners_agree(mesh, t * 3 + s, t * 3 + (s + 1) % 3)', 'sides[s].is_empty() || true', HEAL),
    ('all pieces appended', T, '            if slot < 3 {', '            if false {', HEAL),
    ('per-vertex data check removed', T, '.is_some_and(|&at| carry::vertex_data_agrees(mesh, low, high, m, at))', '.is_some_and(|_| true)', HEAL),
    ('linear channel compared bitwise to the start', C, 'a[k] + (b[k] - a[k]) * t == p[k]', 'a[k] == p[k]', HEAL),
    ('repair always reports applied', H, '                    summary.split_triangles > 0', '                    true', HEAL),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
