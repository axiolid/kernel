# 0085 — The geometry graph has a versioned wire format

- **Status:** Accepted
- **Date:** 2026-10-10
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #267. A source adapter lowers its geometry into `axiolid-model`'s
`GeometryGraph`, and its language bindings (JavaScript/WASM, Python, C,
.NET) want to hand that neutral representation to a host as a value: a
viewer or a kernel in another language could then evaluate the exact
geometry itself instead of receiving only triangles. That makes the
graph a wire format, and until now it was only a Rust data structure.
No encoding was documented, and nothing said which changes to
`GeometryNode`, `Primitive`, `Profile` or a curve or surface variant
break a serialised graph, as opposed to breaking the Rust API.

The maintainer decisions on #267 bind this record: commit now; a JSON
and a binary encoding behind an optional `serde` feature, off by
default; the binary encoding is CBOR; every payload carries an explicit
format version; additions are minor, renames, removals and changes of
meaning or units are major; and a reader that meets content it does not
know refuses the whole payload by name, never skips it.

ADR 0011 declined a flat serialised view until a real consumer existed.
One exists now. ADR 0039 asks that any serialised protocol be
explicitly versioned and preserve unknown fields or fail closed; this
format fails closed.

## Decision

We will publish one versioned wire format for `GeometryGraph`, named
`axiolid-geometry-graph`, starting at version `1.0`, with two encodings
of the same serde data model: JSON (RFC 8259, via `serde_json`) and
CBOR (RFC 8949, via `ciborium`). `axiolid-model` owns the format, its
version and its reader; the representation crates whose values the
graph carries derive `Serialize`/`Deserialize` behind their own
optional `serde` features and own no policy.

### Envelope

Every payload is one map with exactly three entries, in this order:

```json
{
  "format": "axiolid-geometry-graph",
  "version": "1.0",
  "graph": { "nodes": [ ... ], "roots": [ 0, 3 ] }
}
```

- `format` is the literal name. Any other value, or none, is refused.
- `version` is the string `"MAJOR.MINOR"`, two decimal integers.
- `graph.nodes` lists the nodes in insertion order, which is
  topological; a node reference is the zero-based index of an earlier
  node. `graph.roots` lists root indices in the caller's order.

### Data model

- **Enums are externally tagged.** A unit variant is its name as a
  string (`"Union"`); any other variant is a one-entry map from its name
  to its content (`{"Primitive": {"Sphere": {"radius": 1.0}}}`). A
  newtype variant's content is the wrapped value, a struct variant's
  content is a map of its fields, and a tuple field is an array.
  Rust enum discriminants never reach the wire.
- **Names are the Rust names as of format 1.0**: node kinds, variants
  and fields keep the spelling they had when this format was published.
  A later Rust rename must keep the wire name with
  `#[serde(rename = "...")]`; the committed golden payloads
  (`crates/representations/modeling/graph/tests/wire/`) fail the build
  if it does not.
- **Vectors and transforms are arrays** in the layout `glam`'s own
  serde support writes: a 2D or 3D point or vector is `[x, y]` or
  `[x, y, z]`; a 3D affine transform is twelve numbers, the three
  columns of its linear part then the translation; a 2D affine transform
  is six numbers in the same column order.
- **An absent optional value is written as `null`**; a reader treats a
  missing optional field as absent as well.
- **Integers** (node and topology indices, degrees, multiplicities) are
  unsigned integers. **Real numbers** are IEEE-754 doubles and
  round-trip bit-exactly: JSON writes the shortest decimal that reads
  back to the same double, `-0.0` included, and CBOR writes the
  shortest IEEE float (half, single or double) that holds the value
  exactly.
- **An integer reads as a real when it is exactly a double.** JavaScript
  and Python CBOR encoders routinely write whole-number doubles as
  integer items, and JSON has one number type, so where a real is
  expected both readers accept an integer `i` with
  `i as f64 as i128 == i` (every `|i| <= 2^53`, and larger integers that
  are doubles); a lossless conversion cannot produce wrong geometry. A
  CBOR integer that is not exactly a double is refused as `Malformed`,
  naming its path, wherever it stands: no index or count of the format
  reaches `2^53`, and rounding one would change the geometry silently.
  Integer zero has no sign and reads as `+0.0`; a payload that needs
  `-0.0` writes it as a float, as this writer does. The CBOR reader
  maps the decoded item tree onto the JSON data model before
  deserialising, so the two encodings are read by one deserialiser; it
  unwraps the self-describe tag (55799) and refuses every other tag, byte
  strings, simple values, non-text map keys and repeated keys as
  `Malformed`.
- **Numbers are finite.** JSON cannot carry a NaN or an infinity, so
  neither encoding does: a writer refuses a graph holding one, naming
  where it is, and a reader refuses a payload holding one. Lengths and
  angles keep the units the Rust values document (the caller's length
  unit, radians).

### Version rules

| Change | Version |
| --- | --- |
| Add a node kind, an enum variant, or an optional field (absent means what the payload meant before) | minor |
| Rename or remove a kind, variant or field; change a field's type, meaning or units; make an optional field required | major |
| Change the envelope | major |

A payload of an older minor of the same major stays readable by every
later reader of that major. The committed golden payloads are the
regression guard: each published version keeps its own, and today's
reader must read each one to the identical graph.

### Refusal rules

A reader returns a typed `WireError` and no graph whenever it cannot
guarantee the whole payload:

- **Version first.** Before reading the graph, the reader checks the
  envelope. A wrong or missing `format` is `UnknownFormat`, a missing
  `version` is `MissingVersion`, and a version it cannot parse, of
  another major, or of a newer minor than it supports is
  `UnsupportedVersion`, naming the version found and the version
  supported. A newer minor is refused even though it might contain
  nothing new, because the reader cannot know that without the
  vocabulary it lacks.
- **Unknown kinds and variants.** A node kind or variant the reader
  does not know, in a payload whose version it supports, is
  `UnknownKind`, naming the kind and the payload's version. It is never
  skipped: skipping would return incomplete geometry silently.
- **Unknown fields.** A field the reader does not know is
  `UnknownField`. Adding a field is a minor change, so a genuine newer
  payload is already refused by its version; an unknown field in a
  supported version is a malformed or hand-edited payload.
- **Validation on read.** A decoded graph is rebuilt node by node
  through `GeometryGraphBuilder::push`, the same validation as any
  other construction: references to earlier nodes of the right family,
  dimension rules, station rules. A failure is `InvalidGraph`, naming
  the node index and the `GraphError`. There is no unchecked
  construction path; a node reference deserialises only inside a graph
  payload, where it takes the new graph's brand.
- **Everything else** that does not decode (bad syntax, a wrong type,
  a missing required field, trailing bytes, exceeding the decoders'
  nesting limit) is `Malformed`, with the decoder's message. A
  non-finite number (a CBOR NaN or infinity, a JSON literal past the
  double range) is `NonFinite`, with its path or position.

### Dependencies

The encoders are optional and outside the default build: `serde`
(derive), `serde_json` and `ciborium` are enabled only by the `serde`
feature. All are pure Rust and permissively licensed (`serde`,
`serde_json`: MIT or Apache-2.0; `ciborium`: Apache-2.0, its `half`
dependency MIT or Apache-2.0); `ciborium` was already in the lockfile
through the benchmark harness. `glam`'s `serde` feature supplies the
vector and transform impls. No internal dependency edge and no closure
profile changes: serde is external, and every allowlisted internal
edge stays as it was.

`axiolid-core`, `axiolid-linear`, `axiolid-curve`, `axiolid-surface`,
`axiolid-profile`, `axiolid-primitive`, `axiolid-mesh`,
`axiolid-topology` and `axiolid-model` each gain a `serde` feature; the
facade's `serde` feature forwards to whichever of them are in the
build.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Skip unknown node kinds | Returns incomplete geometry silently, the failure this project exists to refuse. |
| Accept a newer minor and refuse only if an unknown name appears | Sound only if every minor change is a new name. A minor change may also add an optional field whose absence an older reader would misread; the version is the only reliable signal. |
| Tolerate unknown fields | An older reader would drop data it was handed; failing closed matches ADR 0039. |
| A Rust-specific binary layout (`bincode`, `postcard`) | No mature reader in JavaScript, Python, C or .NET; CBOR has one in each. |
| A schema language (Protocol Buffers, FlatBuffers) | Pulls a code generator and a second source of truth into a format-neutral crate; the serde data model already is the schema, and the golden payloads pin it. |
| Explicit `#[serde(rename)]` on every variant and field now | Hundreds of attributes repeating the Rust names. The golden payloads catch a rename that changes the wire, and the rename attribute is added then. |
| Hand-written wire mirror types | A second copy of every representation type to keep in step; the derive plus golden files gives the same guarantee. |
| Serialise non-finite numbers in CBOR only | The two encodings would no longer carry the same set of graphs. |
| Require float items for reals in CBOR | Refuses what common JavaScript and Python encoders write for whole numbers, though the conversion is lossless; JSON already accepts integers. |
| Round any CBOR integer to the nearest double | Changes the geometry silently past `2^53`. |

## Consequences

**Positive**

- A host in any language can read and write the exact graph with an
  off-the-shelf JSON or CBOR library.
- A breaking wire change cannot land unnoticed: the golden payloads
  and the variant-coverage test fail.
- A reader never returns a partial graph.

**Negative / costs**

- A wire name is now frozen beside each Rust name; renaming a variant
  costs a `serde(rename)` attribute.
- Reading parses the envelope header before the graph, and the CBOR
  reader decodes to a value tree and maps it onto the JSON data model
  first (checking finiteness and integer exactness), so reading costs
  more than one pass.
- An integer `-0` does not exist, so a host that writes `-0.0` as an
  integer zero reads back `+0.0`.
- A graph holding a non-finite number, which the builder accepts, cannot
  be written.

**Follow-ups / risks to watch**

- A `glam` upgrade that changed its serde layout would change the wire;
  the golden payloads catch it, and the old layout is then kept by hand.
- Very deep nesting (a long `Profile::Derived` chain, a deeply nested
  piecewise law) is limited by the decoders' recursion limits and
  refused as `Malformed`.

## Relation to existing code

- `crates/representations/modeling/graph/src/wire.rs`: the envelope,
  `FORMAT_NAME`, `FORMAT_VERSION`, `WireError`, `to_json`,
  `from_json`, `to_cbor`, `from_cbor`.
- `crates/representations/modeling/graph/src/id.rs`: node references
  serialise as indices and deserialise only under a graph payload.
- `crates/representations/modeling/graph/tests/wire.rs` and
  `tests/wire/`: round trips over every node kind and variant, golden
  payloads, refusals.
- The `serde` features of the representation crates listed above.
