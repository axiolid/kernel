"""Mutation probe for the drawings on `Application` (#371, ADR 0105).

Each mutant drops or crosses one thing the facade owes the drawing bundle (a
feature, a builder call, a descriptor row, the provider named in a refusal,
the memoised verdict) and must turn the closure gate or a test of the facade
red.

Not observable, and so not a mutant: removing the verdict check in
`Application::hidden_lines` / `visible_faces` (the `verified(..).map_err(..)?`
step). The portable providers pass their suites, so the check never fires in
a test; that is why `verified` is generic and tested directly (its mutant
below), and the call sites are two lines each.
"""
import probe_lib

DRAWING = "crates/facade/axiolid/src/application/drawing.rs"
APPLICATION = "crates/facade/axiolid/src/application.rs"
INTEGRATION = "crates/facade/axiolid/src/integration.rs"
MANIFEST = "crates/facade/axiolid/Cargo.toml"

FACADE = ["cargo", "test", "-q", "-p", "axiolid", "--features", "application"]
INTEGRATION_TEST = FACADE + ["--test", "application"]
LIB_TEST = FACADE + ["--lib"]
CLOSURE = ["cargo", "xtask", "architecture", "closure", "check"]

MUTANTS = [
    ('provider feature dropped from application', MANIFEST,
     '    "dispatch-mesh-hidden-lines",\n    "hidden-lines-provider",\n    "portable-provider",',
     '    "dispatch-mesh-hidden-lines",\n    "portable-provider",', [CLOSURE, INTEGRATION_TEST]),
    ('contract dropped from discrete', MANIFEST,
     '    "mesh-section",\n    "mesh-hidden-lines",\n    "generate",',
     '    "mesh-section",\n    "generate",', [CLOSURE]),
    ('descriptor provider ids swapped', DRAWING,
     [('            Operation::MeshHiddenLines,\n            hidden_lines,\n',
       '            Operation::MeshHiddenLines,\n            visible_faces,\n'),
      ('            Operation::MeshVisibleFaces,\n            visible_faces,\n',
       '            Operation::MeshVisibleFaces,\n            hidden_lines,\n')],
     [INTEGRATION_TEST]),
    ('descriptor exactness exact', DRAWING,
     '                exactness: Exactness::ToleranceBounded,',
     '                exactness: Exactness::Exact,', [INTEGRATION_TEST]),
    ('portable omits hidden lines', APPLICATION,
     '            .with_portable_hidden_lines()\n', '', [INTEGRATION_TEST]),
    ('portable omits visible faces', APPLICATION,
     '            .with_portable_visible_faces()\n', '', [INTEGRATION_TEST]),
    ('refusal loses the selected provider', DRAWING,
     '        let provider = self.hidden_lines_provider.unwrap_or(APPLICATION_ID);',
     '        let provider = APPLICATION_ID;', [INTEGRATION_TEST]),
    ('refusal source is a plain geometry error', DRAWING,
     '        source: Box::new(ApplicationErrorSource::HiddenLines(source)),',
     '        source: Box::new(ApplicationErrorSource::Geometry(source.error)),', [INTEGRATION_TEST]),
    ('audit not forwarded to the registry', DRAWING,
     '        self.hidden_lines.audit(request, options).map_err(|source| {',
     '        self.hidden_lines.audit(request, options).map(|_| Vec::new()).map_err(|source| {', [INTEGRATION_TEST]),
    ('failing verdict accepted', DRAWING,
     '        if report.is_success() {\n            Ok(())',
     '        if true {\n            Ok(())', [LIB_TEST]),
    ('face refusal reports the line operation', DRAWING,
     "|source| drawing_error(Operation::MeshVisibleFaces, provider, options, source))\n    }",
     "|source| drawing_error(Operation::MeshHiddenLines, provider, options, source))\n    }", [INTEGRATION_TEST]),
    ('face refusal loses the selected provider', DRAWING,
     '        let provider = self.visible_faces_provider.unwrap_or(APPLICATION_ID);',
     '        let provider = APPLICATION_ID;', [INTEGRATION_TEST]),
    ('audit refusal loses the selected provider', DRAWING,
     '                self.hidden_lines_provider.unwrap_or(APPLICATION_ID),\n                options,',
     '                APPLICATION_ID,\n                options,', [INTEGRATION_TEST]),
    ('face report always a success', APPLICATION,
     '            Self::MeshVisibleFaces(report) => report.is_success(),',
     '            Self::MeshVisibleFaces(_) => true,', [LIB_TEST]),
    ('line report always a success', APPLICATION,
     '            Self::MeshHiddenLines(report) => report.is_success(),',
     '            Self::MeshHiddenLines(_) => true,', [LIB_TEST]),
    ('enabled features omit the provider', INTEGRATION,
     '    compiled_feature!(features, "hidden-lines-provider");\n', '', [INTEGRATION_TEST]),
]

probe_lib.main(MUTANTS, timeout=1800, command=())
