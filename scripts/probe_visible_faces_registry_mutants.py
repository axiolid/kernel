"""Mutation probe for the visible-face registry (#383, ADR 0094).

Each mutant removes one check the registry owes the contract, loosens a
limit by one, or lets a refusal lose the mesh it is about, and must turn a
registry test red.
"""
import probe_lib

R = "crates/execution/dispatch/src/visible_faces.rs"
TESTS = [
    "-p", "axiolid-dispatch", "--no-default-features", "--features", "mesh-hidden-lines",
    "--test", "mesh_visible_faces_contract",
]

MUTANTS = [
    ('group count unchecked', R, '            if groups.len() != entry.mesh.triangle_count() {', '            if false && groups.len() != entry.mesh.triangle_count() {', TESTS),
    ('group count refusal loses its mesh id', R, '                return Err(HiddenLineError::mesh(\n                    entry.id,\n                    GeomError::InvalidInput(format!(\n                        "{LABEL} mesh has', '                return Err(HiddenLineError::request(\n                    GeomError::InvalidInput(format!(\n                        "{LABEL} mesh has', TESTS),
    ('view change accepted', R, '    if outcome.view != request.view {', '    if false && outcome.view != request.view {', TESTS),
    ('region of a foreign mesh accepted', R, '            .position(|entry| entry.id == region.mesh)', '            .position(|entry| entry.id == region.mesh || true)', TESTS),
    ('short ring accepted', R, '            if ring.len() < 3 {', '            if false && ring.len() < 3 {', TESTS),
    ('repeated ring point accepted', R, '                .any(|(a, b)| a == b)', '                .any(|(a, b)| a == b && false)', TESTS),
    ('ring wrap-around pair ignored', R, '                .zip(ring.iter().cycle().skip(1))', '                .zip(ring.iter().skip(1))', TESTS),
    ('non-finite ring point accepted', R, '            if ring.iter().any(|p| !p.is_finite()) {', '            if false && ring.iter().any(|p| !p.is_finite()) {', TESTS),
    ('holes not counted as region vertices', R, '        for ring in std::iter::once(&region.outer).chain(&region.holes) {', '        for ring in std::iter::once(&region.outer) {', TESTS),
    ('descending triangles accepted', R, '        if region.triangles.windows(2).any(|w| w[0] >= w[1]) {', '        if false && region.triangles.windows(2).any(|w| w[0] >= w[1]) {', TESTS),
    ('repeated triangle accepted', R, '        if region.triangles.windows(2).any(|w| w[0] >= w[1]) {', '        if region.triangles.windows(2).any(|w| w[0] > w[1]) {', TESTS),
    ('triangle out of range accepted', R, '        if region.triangles.iter().any(|&t| t as usize >= count) {', '        if false && region.triangles.iter().any(|&t| t as usize >= count) {', TESTS),
    ('triangle of another group accepted', R, '                    .any(|&t| entry.groups.of(t as usize) != group)', '                    .any(|&t| entry.groups.of(t as usize) != group && false)', TESTS),
    ('cap of a sheet accepted', R, '                if entry.role != MeshRole::Solid || request.view.near.is_none() {', '                if request.view.near.is_none() {', TESTS),
    ('cap without a near plane accepted', R, '                if entry.role != MeshRole::Solid || request.view.near.is_none() {', '                if entry.role != MeshRole::Solid {', TESTS),
    ('cap with triangles accepted', R, '                if !region.triangles.is_empty() {', '                if false && !region.triangles.is_empty() {', TESTS),
    ('caps not counted', R, '                caps += 1;', '                caps += 0;', TESTS),
    ('mesh order unchecked', R, '            let order = last_position\n                .cmp(&position)\n', '            let order = Ordering::Equal\n', TESTS),
    ('source order unchecked', R, '                .then(last_source.cmp(&region.source))', '                .then(Ordering::Less)', TESTS),
    ('outer ring order unchecked', R, '                .then_with(|| sequence_order(last_outer, &region.outer));', '                .then_with(|| Ordering::Less);', TESTS),
    ('order not enforced', R, '            if order.is_ge() {', '            if false && order.is_ge() {', TESTS),
    ('equal regions accepted', R, '            if order.is_ge() {', '            if order.is_gt() {', TESTS),
    ('region limit unchecked', R, '    if outcome.regions.len() > limits.max_regions {', '    if false && outcome.regions.len() > limits.max_regions {', TESTS),
    ('region limit exclusive', R, '    if outcome.regions.len() > limits.max_regions {', '    if outcome.regions.len() >= limits.max_regions {', TESTS),
    ('region vertex limit unchecked', R, '    if vertices > limits.max_region_vertices {', '    if false && vertices > limits.max_region_vertices {', TESTS),
    ('region vertex limit exclusive', R, '    if vertices > limits.max_region_vertices {', '    if vertices >= limits.max_region_vertices {', TESTS),
    ('candidate limit unchecked', R, '    if evidence.candidate_pairs > limits.max_candidate_pairs {', '    if false && evidence.candidate_pairs > limits.max_candidate_pairs {', TESTS),
    ('candidate limit exclusive', R, '    if evidence.candidate_pairs > limits.max_candidate_pairs {', '    if evidence.candidate_pairs >= limits.max_candidate_pairs {', TESTS),
    ('source meshes unchecked', R, '        || evidence.source_meshes != request.meshes.len()\n', '', TESTS),
    ('source triangles unchecked', R, '        || Some(evidence.source_triangles) != triangles\n', '', TESTS),
    ('region count unchecked', R, '        || evidence.regions != outcome.regions.len()\n', '', TESTS),
    ('cap count unchecked', R, '        || evidence.cap_regions != caps\n', '', TESTS),
    ('refusal naming a foreign mesh accepted', R, '                if !request.meshes.iter().any(|entry| entry.id == *id) {', '                if false && !request.meshes.iter().any(|entry| entry.id == *id) {', TESTS),
    ('fallback on any error', R, '                Err(error) => return Err(error),', '                Err(error) => {\n                    last_retryable = Some(error);\n                }', TESTS),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
