"""Mutation probe for the projection routing and refusal (#387, ADR 0107).

Each mutant removes one part of how a perspective view is kept from a provider
that cannot draw it: the registry's routing by `supports_projection`, the view
validation, the default claim, the providers' refusal by name, or the claim
check of the perspective suite, and must turn a test red.
"""
import probe_lib

C = "crates/contracts/operations/mesh-hidden-lines/src/contract.rs"
F = "crates/contracts/operations/mesh-hidden-lines/src/faces.rs"
S = "crates/contracts/operations/mesh-hidden-lines/src/conformance/perspective.rs"
HR = "crates/execution/dispatch/src/hidden_lines.rs"
VF = "crates/execution/dispatch/src/visible_faces.rs"
VC = "crates/execution/dispatch/src/view_checks.rs"
VW = "crates/algorithms/discrete/hidden-lines/src/view.rs"
HP = "crates/algorithms/discrete/hidden-lines/src/provider.rs"
FP = "crates/algorithms/discrete/hidden-lines/src/faces/provider.rs"

CONTRACT = ["-p", "axiolid-mesh-hidden-lines-contract", "--lib"]
LINES = [
    "-p", "axiolid-dispatch", "--no-default-features", "--features", "mesh-hidden-lines",
    "--test", "mesh_hidden_lines_contract",
]
FACES = [
    "-p", "axiolid-dispatch", "--no-default-features", "--features", "mesh-hidden-lines",
    "--test", "mesh_visible_faces_contract",
]
PROVIDERS = ["-p", "axiolid-hidden-lines", "--test", "perspective"]

FIRST_SITE = (
    '            if matches_device(options.device(), descriptor.id, descriptor.target)\n'
    '                && entry.provider.supports_projection(request.view.projection)\n'
)
FIRST_SITE_OFF = '            if matches_device(options.device(), descriptor.id, descriptor.target)\n'
SECOND_SITE = '                || !entry.provider.supports_projection(request.view.projection)\n'
DEFAULT_CLAIM = '        projection == ViewProjection::Orthographic\n    }\n'
CLAIM_LINES = (
    'pub fn run(provider: &dyn MeshHiddenLines) -> PerspectiveConformanceReport {\n'
    '        let claimed = provider.supports_projection(ViewProjection::Perspective);'
)
CLAIM_FACES = (
    'pub fn run_faces(provider: &dyn MeshVisibleFaces) -> PerspectiveConformanceReport {\n'
    '        let claimed = provider.supports_projection(ViewProjection::Perspective);'
)

MUTANTS = [
    # Registry routing.
    ('hidden lines: unsupported views still count as a matching provider', HR, FIRST_SITE, FIRST_SITE_OFF, LINES),
    ('hidden lines: an unsupporting provider is called', HR, SECOND_SITE, '', LINES),
    ('visible faces: unsupported views still count as a matching provider', VF, FIRST_SITE, FIRST_SITE_OFF, FACES),
    ('visible faces: an unsupporting provider is called', VF, SECOND_SITE, '', FACES),
    # View validation.
    ('a perspective view without a near depth accepted', VC, '            if !view.near.is_some_and(|near| near < 0.0) {', '            if false && !view.near.is_some_and(|near| near < 0.0) {', LINES),
    ('a perspective view with near zero accepted', VC, 'view.near.is_some_and(|near| near < 0.0)', 'view.near.is_some_and(|near| near <= 0.0)', LINES),
    ('a perspective view with a positive near accepted (faces)', VC, 'view.near.is_some_and(|near| near < 0.0)', 'view.near.is_some_and(|near| near < 1.0)', FACES),
    # Constructor.
    ('the constructor accepts near zero', C, '        if !near.is_finite() || near >= 0.0 {', '        if !near.is_finite() || near > 0.0 {', CONTRACT),
    ('the constructor accepts a positive near', C, '        if !near.is_finite() || near >= 0.0 {', '        if !near.is_finite() || near > 1.0 {', CONTRACT),
    ('the constructor accepts a non-finite near', C, '        if !near.is_finite() || near >= 0.0 {', '        if near >= 0.0 {', CONTRACT),
    ('the constructor builds an orthographic view', C, '            projection: ViewProjection::Perspective,\n        })', '            projection: ViewProjection::Orthographic,\n        })', CONTRACT),
    # Default claim.
    ('hidden lines: the default claims every projection', C, DEFAULT_CLAIM, '        true\n    }\n', PROVIDERS),
    ('visible faces: the default claims every projection', F, DEFAULT_CLAIM, '        true\n    }\n', PROVIDERS),
    ('hidden lines: the default claims no projection', C, DEFAULT_CLAIM, '        false\n    }\n', PROVIDERS),
    # Providers refuse by name.
    ('the providers draw a perspective view', VW, '    if view.projection == ViewProjection::Orthographic {\n        return Ok(());', '    if true {\n        return Ok(());', PROVIDERS),
    ('the providers refuse an orthographic view', VW, '    if view.projection == ViewProjection::Orthographic {\n        return Ok(());', '    if view.projection != ViewProjection::Orthographic {\n        return Ok(());', PROVIDERS),
    # The suite.
    ('the suite never infers a hidden-line claim', S, CLAIM_LINES, CLAIM_LINES.replace('provider.supports_projection(ViewProjection::Perspective)', 'false'), CONTRACT),
    ('the suite never infers a visible-face claim', S, CLAIM_FACES, CLAIM_FACES.replace('provider.supports_projection(ViewProjection::Perspective)', 'false'), CONTRACT),
    ('the suite accepts a provider that draws', S, '        Ok(()) => report\n            .failures\n            .push(PerspectiveConformanceFailure::UnclaimedButDrawn { operation }),', '        Ok(()) => {}', CONTRACT),
    ('the suite accepts any refusal', S, '            error: GeomError::Unsupported { .. } | GeomError::UnsupportedInput { .. },', '            error: _,', CONTRACT),
    # Registration.
    ('hidden lines: registration ignores the perspective report', HR, '                .map(PerspectiveConformanceFailure::HiddenLinesBase),\n        );\n        if !report.is_success() {', '                .map(PerspectiveConformanceFailure::HiddenLinesBase),\n        );\n        if false && !report.is_success() {', LINES),
    ('visible faces: registration ignores the perspective report', VF, '                .map(PerspectiveConformanceFailure::VisibleFacesBase),\n        );\n        if !report.is_success() {', '                .map(PerspectiveConformanceFailure::VisibleFacesBase),\n        );\n        if false && !report.is_success() {', FACES),

    # The providers' own refusal, one mutant per provider.
    ('ScalarHiddenLines draws a perspective view', HP, '        require_orthographic(&request.view, Self::ID, Operation::MeshHiddenLines)?;\n        draw(BroadPhase::Tree, request, options)', '        draw(BroadPhase::Tree, request, options)', PROVIDERS),
    ('ScanHiddenLines draws a perspective view', HP, '        require_orthographic(\n            &request.view,\n            self.descriptor().id,\n            Operation::MeshHiddenLines,\n        )?;\n        draw(BroadPhase::Scan, request, options)', '        draw(BroadPhase::Scan, request, options)', PROVIDERS),
    ('ScalarVisibleFaces draws a perspective view', FP, '        require_orthographic(&request.view, Self::ID, Operation::MeshVisibleFaces)?;\n        run(BroadPhase::Tree, true, request, options)', '        run(BroadPhase::Tree, true, request, options)', PROVIDERS),
    ('ScanVisibleFaces draws a perspective view', FP, '        require_orthographic(\n            &request.view,\n            self.descriptor().id,\n            Operation::MeshVisibleFaces,\n        )?;\n        run(BroadPhase::Scan, true, request, options)', '        run(BroadPhase::Scan, true, request, options)', PROVIDERS),
    ('UnprunedVisibleFaces draws a perspective view', FP, '        require_orthographic(\n            &request.view,\n            self.descriptor().id,\n            Operation::MeshVisibleFaces,\n        )?;\n        run(BroadPhase::Tree, false, request, options)', '        run(BroadPhase::Tree, false, request, options)', PROVIDERS),
    # The default claim of visible faces.
    ('visible faces: the default claims no projection', F, DEFAULT_CLAIM, '        false\n    }\n', PROVIDERS),
    # Registration.
    ('hidden lines: perspective registration skips the base suite', HR, '                .map(PerspectiveConformanceFailure::HiddenLinesBase),', '                .filter(|_| false)\n                .map(PerspectiveConformanceFailure::HiddenLinesBase),', LINES),
    ('visible faces: perspective registration skips the base suite', VF, '                .map(PerspectiveConformanceFailure::VisibleFacesBase),', '                .filter(|_| false)\n                .map(PerspectiveConformanceFailure::VisibleFacesBase),', FACES),
    ('hidden lines: register_conformant ignores a perspective claim', HR, '        if provider.supports_projection(ViewProjection::Perspective) {', '        if false && provider.supports_projection(ViewProjection::Perspective) {', PROVIDERS),
    ('visible faces: register_conformant ignores a perspective claim', VF, '        if provider.supports_projection(ViewProjection::Perspective) {', '        if false && provider.supports_projection(ViewProjection::Perspective) {', PROVIDERS),
    ('hidden lines: register_conformant accepts a failed perspective suite', HR, '            if !perspective.is_success() {', '            if false && !perspective.is_success() {', PROVIDERS),
    ('visible faces: register_conformant accepts a failed perspective suite', VF, '            if !perspective.is_success() {', '            if false && !perspective.is_success() {', PROVIDERS),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
