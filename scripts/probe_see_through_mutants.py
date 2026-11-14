"""Mutation probe for see-through triangles (#496, ADR 0111).

Each mutant breaks one part of how see-through triangles are kept from a
provider that does not honour them, checked, counted or stacked: the default
claim, the suite's claim inference, the registries' routing and outcome
checks, the hidden-line count, the visible-face stacks, the mixed-group and
budget refusals, and must turn a named test red.
"""
import probe_lib

C = "crates/contracts/operations/mesh-hidden-lines/src/contract.rs"
F = "crates/contracts/operations/mesh-hidden-lines/src/faces.rs"
S = "crates/contracts/operations/mesh-hidden-lines/src/conformance/see_through.rs"
HR = "crates/execution/dispatch/src/hidden_lines.rs"
VF = "crates/execution/dispatch/src/visible_faces.rs"
VC = "crates/execution/dispatch/src/view_checks.rs"
HL = "crates/algorithms/discrete/hidden-lines/src/"
OCC = HL + "occlusion.rs"
EDG = HL + "edges.rs"
PRO = HL + "provider.rs"
GRP = HL + "groups.rs"
MOD = HL + "faces.rs"
LAB = HL + "faces/label.rs"
REG = HL + "faces/regions.rs"
FPR = HL + "faces/provider.rs"

SUITE = ["-p", "axiolid-mesh-hidden-lines-contract", "--lib", "see_through"]
DISPATCH = ["-p", "axiolid-dispatch", "--no-default-features", "--features", "mesh-hidden-lines"]
LINES = DISPATCH + ["--test", "mesh_hidden_lines_contract"]
FACES = DISPATCH + ["--test", "mesh_visible_faces_contract"]
PROVIDERS = ["-p", "axiolid-hidden-lines", "--test", "see_through"]
LINE_SUITE = PROVIDERS + ["the_hidden_line_providers_pass_the_see_through_suite"]
FACE_SUITE = PROVIDERS + ["the_visible_face_providers_pass_the_see_through_suite"]
OCC_LIB = ["-p", "axiolid-hidden-lines", "--lib", "occlusion::tests"]
LABEL_LIB = ["-p", "axiolid-hidden-lines", "--lib", "label::tests"]
GROUP_LIB = ["-p", "axiolid-hidden-lines", "--lib", "groups::tests"]
ORACLE = ["-p", "axiolid-hidden-lines", "--test", "see_through_oracle"]

DEFAULT_CLAIM = '        occlusion == Occlusion::Opaque\n    }\n'
CLAIM_LINES = (
    'pub fn run(provider: &dyn MeshHiddenLines) -> SeeThroughConformanceReport {\n'
    '        let claimed = provider.supports_occlusion(Occlusion::SeeThrough);'
)
CLAIM_FACES = (
    'pub fn run_faces(provider: &dyn MeshVisibleFaces) -> SeeThroughConformanceReport {\n'
    '        let claimed = provider.supports_occlusion(Occlusion::SeeThrough);'
)
FIRST_SITE = '                && (!see_through || entry.provider.supports_occlusion(Occlusion::SeeThrough))\n'
SECOND_SITE = '                || (see_through && !entry.provider.supports_occlusion(Occlusion::SeeThrough))\n'
ADJACENT = '            && previous_visibility == Some((segment.visibility, segment.seen_through))'
PREVIOUS = '        previous_visibility = Some((segment.visibility, segment.seen_through));'
REGISTER = '        if provider.supports_occlusion(Occlusion::SeeThrough) {'

MUTANTS = [
    # 1. The default claim.
    ('hidden lines: the default claims see-through triangles', C, DEFAULT_CLAIM, '        true\n    }\n', SUITE),
    ('visible faces: the default claims see-through triangles', F, DEFAULT_CLAIM, '        true\n    }\n', SUITE),
    # 2, 3. The suite infers the claim.
    ('the suite never infers a hidden-line claim', S, CLAIM_LINES, CLAIM_LINES.replace('provider.supports_occlusion(Occlusion::SeeThrough)', 'false'), SUITE),
    ('the suite never infers a visible-face claim', S, CLAIM_FACES, CLAIM_FACES.replace('provider.supports_occlusion(Occlusion::SeeThrough)', 'false'), SUITE),
    # 4 to 7. Registry routing.
    ('hidden lines: an unclaiming provider counts as matching', HR, FIRST_SITE, '', LINES),
    ('hidden lines: an unclaiming provider is called', HR, SECOND_SITE, '', LINES),
    ('visible faces: an unclaiming provider counts as matching', VF, FIRST_SITE, '', FACES),
    ('visible faces: an unclaiming provider is called', VF, SECOND_SITE, '', FACES),
    # 8. The occlusion length check, in the registry and in the provider.
    ('the registry ignores the occlusion length', VC, '        TriangleOcclusion::PerTriangle(flags) if flags.len() != entry.mesh.triangle_count() => {', '        TriangleOcclusion::PerTriangle(flags) if false && flags.len() != entry.mesh.triangle_count() => {', LINES),
    ('the provider ignores the occlusion length', EDG, '            if flags.len() != count {', '            if false && flags.len() != count {', PROVIDERS),
    # 9, 10. Hidden-line outcome checks.
    ('the adjacency check ignores the count', HR, [(ADJACENT, ADJACENT.replace('segment.seen_through', '0')), (PREVIOUS, PREVIOUS.replace('segment.seen_through', '0'))], LINES),
    ('a hidden piece with a count is accepted', HR, '            if segment.visibility == Visibility::Hidden {\n                return Err(violation("hidden segment counts see-through triangles"));', '            if false {\n                return Err(violation("hidden segment counts see-through triangles"));', LINES),
    ('a count without a see-through triangle is accepted', HR, '            if !see_through {', '            if false {', LINES),
    # 11, 12. Visible-face outcome checks.
    ('the behind layer order is not checked', VF, '        if !in_order {', '        if false && !in_order {', FACES),
    ('the registry ignores the behind budget', VF, '        if behind > limits.max_behind_surfaces {', '        if false && behind > limits.max_behind_surfaces {', FACES),
    ('regions are not ordered by what is behind them', VF, '            let order = order.then_with(|| behind_order(request, last_behind, &region.behind));', '            let order = order;', FACES),
    ('the registry accepts a mixed single group', VF, '                if flags.contains(&Occlusion::Opaque) {', '                if false && flags.contains(&Occlusion::Opaque) {', FACES),
    ('the registry accepts a region listing its own source', VF, '        if surface.layer == 0 && (surface.mesh, surface.source) == (region.mesh, region.source) {', '        if false && surface.layer == 0 && (surface.mesh, surface.source) == (region.mesh, region.source) {', FACES),
    ('the registry accepts an occluder before the last layer', VF, '        if occludes && surface.layer != last.layer {', '        if false && occludes && surface.layer != last.layer {', FACES),
    ('the registry accepts a cap behind a mesh that has none', VF, '        if surface.source == FaceSource::Cap\n', '        if false && surface.source == FaceSource::Cap\n', FACES),
    ('the registry accepts a stack behind an opaque region', VF, '    if region.source == FaceSource::Cap || !entry.occlusion.has_see_through() || opaque_group {', '    if false {', FACES),
    ('visible faces: behind surfaces are not checked', VF, '        || evidence.behind_surfaces != behind\n', '', FACES),
    ('visible faces: see-through triangles are not checked', VF, '        || evidence.see_through_triangles != see_through_triangles(request.meshes)\n', '', FACES),
    ('hidden lines: seen-through segments are not checked', HR, '        || evidence.seen_through_segments != seen_through\n', '', LINES),
    ('hidden lines: see-through triangles are not checked', HR, '        || evidence.see_through_triangles != see_through_triangles(request.meshes)\n', '', LINES),
    ('a short occlusion slice is routed first', HR, '        validate_occlusion_lengths(request.meshes)?;\n', '', LINES),
    ('the registry accepts a mixed face group', VF, '        if let Some(pair) = pairs.windows(2).find(|pair| pair[0].0 == pair[1].0) {', '        if let Some(pair) = pairs.windows(2).find(|pair| false && pair[0].0 == pair[1].0) {', FACES),
    # Registration verifies a claim.
    ('hidden lines: register_conformant ignores a see-through claim', HR, REGISTER, '        if false && provider.supports_occlusion(Occlusion::SeeThrough) {', PROVIDERS + ["register_conformant"]),
    ('visible faces: register_conformant ignores a see-through claim', VF, REGISTER, '        if false && provider.supports_occlusion(Occlusion::SeeThrough) {', PROVIDERS + ["register_conformant"]),
    # 13 to 17. The hidden-line count.
    ('a see-through triangle hides', OCC, '        if tri.see_through {', '        if false && tri.see_through {', LINE_SUITE),
    ('a see-through triangle is not counted', OCC, '                glass.extend(see_through(solver, scene, tri, kept)?);', '                glass.extend(see_through(solver, scene, tri, kept)?.filter(|_| false));', LINE_SUITE),
    ('the count uses the closed image clip', OCC, '        match solver.clip_perturbed(span, &side)? {', '        match solver.clip(span, &side)? {', LINE_SUITE),
    ('equal counts are not merged', OCC, '            Some(last) if ours && last.seen_through == seen_through => last.hi = to,', '            Some(last) if false && ours && last.seen_through == seen_through => last.hi = to,', LINE_SUITE),
    ('edge-on glass occludes', OCC, '        if tri.see_through {', '        if tri.see_through && tri.sigma != Ordering::Equal {', OCC_LIB),
    ('seen-through segments are not counted', PRO, '            evidence.seen_through_segments += 1;', '            evidence.seen_through_segments += 0;', LINE_SUITE),
    ('the count includes glass at equal depth', OCC, '    let in_front = scene.constraint(Kind::PlaneAbove { t, neg, lift }, true);', '    let in_front = scene.constraint(Kind::PlaneAbove { t, neg, lift }, false);', LINE_SUITE),
    ('the count includes glass beyond the near plane', OCC, '            solver.clip(span, &not_cut)\n', '            let _ = not_cut;\n            Ok(Some(span))\n', [LINE_SUITE, PROVIDERS, ORACLE]),
    ('the sweep counts only the ends before a gap', OCC, '            if solver.cmp(end, at)? == Ordering::Greater {', '            if solver.cmp(end, at)? != Ordering::Less {', [LINE_SUITE, OCC_LIB]),
    ('collapsed neighbours of different counts are joined', PRO, '                if last.visibility == visibility && last.seen_through == piece.seen_through {', '                if last.visibility == visibility {', PROVIDERS + ["a_collapsed_piece"]),
    # 18. The prune.
    ('the faces scene lets glass prune', OCC, '        if tri.see_through {', '        if tri.see_through && scene.count_see_through {', FACE_SUITE),
    # 19 to 21. The visible-face stacks.
    ('the stack does not stop at the first occluder', LAB, '                    Ordering::Less => {}', '                    Ordering::Less => stop.push(surface),', [FACE_SUITE, LABEL_LIB]),
    ('glass behind the occluder stays in the stack', LAB, '                    Ordering::Less => continue,', '                    Ordering::Less => {}', LABEL_LIB),
    ('the stack omits the occluder', LAB, '        if !stop.is_empty() {', '        if false && !stop.is_empty() {', FACE_SUITE),
    ('a tie behind the front is not counted', LAB, '        label.tied |= layer > 0 && keys.len() > 1;', '        label.tied |= false && layer > 0 && keys.len() > 1;', LABEL_LIB),
    ('the region boundary ignores the stack', LAB, '        self.keys == other.keys && self.behind == other.behind', '        self.keys == other.keys', FACE_SUITE),
    ('behind surfaces are not counted', FPR, '    evidence.behind_surfaces = regions.iter().map(|r| r.behind.len()).sum();', '    evidence.behind_surfaces = 0;', FACE_SUITE),
    # 22 to 24. Refusals and planar groups.
    ('the provider accepts a mixed face group', MOD, '    match pairs.windows(2).find(|pair| pair[0].0 == pair[1].0) {', '    match pairs.windows(2).find(|pair| false && pair[0].0 == pair[1].0) {', FACE_SUITE),
    ('planar groups ignore the flags', GRP, '            || flags.is_some_and(|flags| flags[t1 as usize] != flags[t2 as usize])\n', '', [FACE_SUITE, GROUP_LIB]),
    ("the provider ignores the behind budget", REG, '            if behind_total > limits.max_behind_surfaces {', '            if false && behind_total > limits.max_behind_surfaces {', FACE_SUITE),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
