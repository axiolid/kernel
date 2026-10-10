"""Mutation probe for the geometry graph wire format (#267, ADR 0085).

Each mutant breaks one rule the format's promise rests on: the version
check (another major, a newer minor, a missing or unparsed version, another
format name, all checked before the graph), the refusal of an unknown kind
or field by name, the rebuild through the builder's validation (a refused
node skipped instead of refusing the payload, roots dropped, references
branded with a foreign graph, a node reference read outside a graph), the
refusal of non-finite numbers on read and write and of trailing bytes, and
the frozen wire names (a variant, a node kind and a field renamed on the
wire). Each must turn `axiolid-model`'s `tests/wire.rs` red.

Not probed, because no test can tell them apart: the path text inside a
`NonFinite` refusal of a map key that is not text (the writer never writes
one), and the `quoted_after` fallback for a message without a closing
clause (serde always writes one).

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import os, pathlib, signal, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
W = "crates/representations/modeling/graph/src/wire.rs"
I = "crates/representations/modeling/graph/src/id.rs"
N = "crates/representations/modeling/graph/src/node.rs"
S = "crates/representations/analytic/primitive/src/solid.rs"
TARGET = ["-p", "axiolid-model", "--features", "serde", "--test", "wire"]
TIMEOUT = 1800

DENY = (
    '#[cfg_attr(\n'
    '    feature = "serde",\n'
    '    derive(serde::Serialize, serde::Deserialize),\n'
    '    serde(deny_unknown_fields)\n'
    ')]\n'
    'pub enum Primitive {'
)

MUTANTS = [
    ('a newer minor is read', W,
     'payload.major == self.major && payload.minor <= self.minor',
     'payload.major == self.major'),
    ('another major is read', W,
     'payload.major == self.major && payload.minor <= self.minor',
     'payload.minor <= self.minor'),
    ('the version is not checked', W,
     'Some(version) if FORMAT_VERSION.reads(version) => Ok(version),',
     'Some(version) => Ok(version),'),
    ('a missing version reads as the current one', W,
     'let found = header.version.ok_or(WireError::MissingVersion)?;',
     'let found = header.version.unwrap_or_else(|| FORMAT_VERSION.to_string());'),
    ('any format name is read', W,
     'Some(format) if format == FORMAT_NAME => {}',
     'Some(_) => {}'),
    ('the JSON version is checked after the graph', W,
     '    let header: Header = serde_json::from_str(text).map_err(|error| json_error(&error, None))?;\n'
     '    let version = check_header(header)?;',
     '    let header: Header = serde_json::from_str(text).map_err(|error| json_error(&error, None))?;\n'
     '    let version = check_header(header).unwrap_or(FORMAT_VERSION);'),
    ('the CBOR version is checked after the graph', W,
     '    let header: Header = value.deserialized().map_err(malformed)?;\n'
     '    let version = check_header(header)?;',
     '    let header: Header = value.deserialized().map_err(malformed)?;\n'
     '    let version = check_header(header).unwrap_or(FORMAT_VERSION);'),
    ('an unknown kind is reported as malformed', W,
     'return WireError::UnknownKind { kind, version };',
     'return WireError::Malformed { detail: kind };'),
    ('an unknown field is reported as malformed', W,
     'return WireError::UnknownField { field, version };',
     'return WireError::Malformed { detail: field };'),
    ('the envelope tolerates unknown fields', W,
     '#[derive(Deserialize)]\n#[serde(deny_unknown_fields)]\nstruct EnvelopeIn {',
     '#[derive(Deserialize)]\nstruct EnvelopeIn {'),
    ('a value tolerates unknown fields', S,
     DENY,
     DENY.replace('    serde(deny_unknown_fields)\n', '')
         .replace('derive(serde::Serialize, serde::Deserialize),', 'derive(serde::Serialize, serde::Deserialize)')),
    ('a refused node is skipped, not refused', W,
     '            .push(node)\n'
     '            .map_err(|error| WireError::InvalidGraph {\n'
     '                node: Some(index),\n'
     '                error,\n'
     '            })?;',
     '            .push(node)\n'
     '            .ok();\n'
     '        let _ = index;'),
    ('the roots are dropped', W,
     '.finish(graph.roots)',
     '.finish(Vec::new())'),
    ('JSON references take a foreign brand', W,
     '        let _brand = WireBrand::enter(owner);\n        serde_json::from_str(text)',
     '        let _brand = WireBrand::enter(GraphId::fresh());\n        serde_json::from_str(text)'),
    ('a node reference reads outside a graph', I,
     '            None => Err(serde::de::Error::custom(\n'
     '                "a node reference is only meaningful inside a geometry graph payload",\n'
     '            )),',
     '            None => Ok(Self {\n'
     '                graph: GraphId::fresh(),\n'
     '                index,\n'
     '            }),'),
    ('a CBOR NaN is read', W,
     '    check_finite(&value, &mut String::new())?;\n    let owner',
     '    let owner'),
    ('a NaN is written', W,
     '    checked_value(graph)?;\n    serde_json::to_string',
     '    serde_json::to_string'),
    ('a JSON overflow is reported as malformed', W,
     'if message.starts_with("number out of range") {',
     'if false {'),
    ('trailing CBOR bytes are ignored', W,
     'if !rest.is_empty() {',
     'if false {'),
    ('a variant renamed on the wire', S,
     '    Block {',
     '    #[cfg_attr(feature = "serde", serde(rename = "Cuboid"))]\n    Block {'),
    ('a node kind renamed on the wire', N,
     '    BoundingBox(Aabb),',
     '    #[cfg_attr(feature = "serde", serde(rename = "Aabb"))]\n    BoundingBox(Aabb),'),
    ('a field renamed on the wire', N,
     '    /// Reused source node.\n    pub source: NodeId,',
     '    /// Reused source node.\n    #[cfg_attr(feature = "serde", serde(rename = "src"))]\n    pub source: NodeId,'),
]


def run(target):
    """Run the tests; a timeout kills the whole process group. Returns the
    exit code and whether the mutant failed to compile."""
    process = subprocess.Popen(
        ["cargo", "test", "-q", *target],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        start_new_session=True,
    )
    try:
        output, _ = process.communicate(timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.communicate()
        return -1, False
    return process.returncode, "could not compile" in output


survivors = []
broken = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        code, uncompiled = run(TARGET)
    finally:
        path.write_text(original)
    # A mutant that does not compile proves nothing about the tests.
    status = "BROKEN" if uncompiled else "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if uncompiled:
        broken.append(name)
    elif code == 0:
        survivors.append(name)
killed = len(MUTANTS) - len(survivors) - len(broken)
print(f"{killed}/{len(MUTANTS)} killed")
sys.exit(1 if survivors or broken else 0)
