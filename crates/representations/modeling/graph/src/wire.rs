//! The versioned wire format of a [`GeometryGraph`] (ADR 0085).
//!
//! One serde data model, two encodings: JSON (RFC 8259) and CBOR
//! (RFC 8949). Every payload is an envelope naming the format and its
//! version around the graph:
//!
//! ```json
//! {
//!   "format": "axiolid-geometry-graph",
//!   "version": "1.0",
//!   "graph": { "nodes": [{ "Primitive": { "Sphere": { "radius": 1.0 } } }], "roots": [0] }
//! }
//! ```
//!
//! `graph.nodes` lists the nodes in insertion order, which is topological,
//! and a node reference is the zero-based index of an earlier node.
//! Enums are externally tagged by their variant names, fields by their
//! names, a point or vector is an array of its coordinates and a transform
//! the columns of its linear part then its translation. Real numbers
//! round-trip bit-exactly and are finite.
//!
//! # Versions
//!
//! Adding a node kind, an enum variant or an optional field is a minor
//! version; renaming or removing anything, or changing its type, meaning or
//! units, is a major version. A reader reads its own major up to its own
//! minor, so every older minor of [`FORMAT_VERSION`] stays readable.
//!
//! # Refusals
//!
//! A reader never returns part of a payload. It refuses, by a named
//! [`WireError`]: a missing or other format name, a missing version, a
//! version it does not support (another major, or a newer minor, checked
//! before the graph is read), a node kind or variant it does not know
//! ([`WireError::UnknownKind`], never skipped), a field it does not know, a
//! non-finite number, anything that does not decode, and a graph that the
//! builder's validation refuses: a decoded graph is rebuilt node by node
//! through [`GeometryGraphBuilder::push`], so a payload cannot build a graph
//! the API could not.
//!
//! A writer refuses a graph holding a non-finite number, which JSON cannot
//! carry, rather than write a payload that would read back differently.

use core::fmt;

use ciborium::Value;
use serde::{Deserialize, Serialize};

use crate::id::{GraphId, WireBrand};
use crate::{GeometryGraph, GeometryGraphBuilder, GeometryNode, GraphError, NodeId};

/// The format name every payload carries in its `format` entry.
pub const FORMAT_NAME: &str = "axiolid-geometry-graph";

/// The format version this crate writes, and the newest it reads.
pub const FORMAT_VERSION: FormatVersion = FormatVersion::new(1, 0);

/// A `MAJOR.MINOR` wire format version.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FormatVersion {
    /// Incremented by any change an older reader would misread.
    pub major: u32,
    /// Incremented by additions: a node kind, a variant, an optional field.
    pub minor: u32,
}

impl FormatVersion {
    /// The version `major.minor`.
    #[must_use]
    pub const fn new(major: u32, minor: u32) -> Self {
        Self { major, minor }
    }

    /// Parse `"MAJOR.MINOR"`, two decimal integers; `None` for anything else.
    #[must_use]
    pub fn parse(text: &str) -> Option<Self> {
        let (major, minor) = text.split_once('.')?;
        let number = |part: &str| {
            if part.is_empty() || !part.bytes().all(|byte| byte.is_ascii_digit()) {
                None
            } else {
                part.parse().ok()
            }
        };
        Some(Self::new(number(major)?, number(minor)?))
    }

    /// Whether a reader of this version reads a payload of `payload`: the
    /// same major, and a minor no newer than its own.
    #[must_use]
    pub const fn reads(self, payload: Self) -> bool {
        payload.major == self.major && payload.minor <= self.minor
    }
}

impl fmt::Display for FormatVersion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}.{}", self.major, self.minor)
    }
}

/// Why a graph could not be written or a payload could not be read.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WireError {
    /// The payload's `format` entry is missing or names another format.
    UnknownFormat {
        /// The format named, if any.
        found: Option<String>,
    },
    /// The payload has no `version` entry.
    MissingVersion,
    /// The payload's version is not `MAJOR.MINOR`, is of another major, or
    /// is a newer minor than this reader supports.
    UnsupportedVersion {
        /// The version as the payload states it.
        found: String,
        /// The newest version this reader reads.
        supported: FormatVersion,
    },
    /// A node kind or enum variant this reader does not know. The payload
    /// is refused whole: skipping the kind would lose geometry silently.
    UnknownKind {
        /// The kind or variant name.
        kind: String,
        /// The payload's format version.
        version: FormatVersion,
    },
    /// A field this reader does not know.
    UnknownField {
        /// The field name.
        field: String,
        /// The payload's format version.
        version: FormatVersion,
    },
    /// A number that is NaN or infinite, which neither encoding carries.
    NonFinite {
        /// Where it is: a path into the payload for a graph or a CBOR
        /// payload, the decoder's position for a JSON payload.
        path: String,
    },
    /// The decoded graph fails the validation every graph construction
    /// runs.
    InvalidGraph {
        /// The index of the node refused, or `None` when the roots are.
        node: Option<usize>,
        /// Why.
        error: GraphError,
    },
    /// The payload does not decode: bad syntax, a wrong type, a missing
    /// field, trailing bytes, or nesting past the decoder's limit.
    Malformed {
        /// The decoder's message.
        detail: String,
    },
}

impl fmt::Display for WireError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownFormat { found: Some(found) } => {
                write!(f, "payload format `{found}` is not `{FORMAT_NAME}`")
            }
            Self::UnknownFormat { found: None } => {
                write!(f, "payload names no format; expected `{FORMAT_NAME}`")
            }
            Self::MissingVersion => write!(f, "payload names no format version"),
            Self::UnsupportedVersion { found, supported } => write!(
                f,
                "format version `{found}` is not supported; this reader reads \
                 {major}.0 to {supported}",
                major = supported.major
            ),
            Self::UnknownKind { kind, version } => write!(
                f,
                "unknown kind `{kind}` in a format {version} payload; refusing the whole payload"
            ),
            Self::UnknownField { field, version } => {
                write!(f, "unknown field `{field}` in a format {version} payload")
            }
            Self::NonFinite { path } => write!(f, "non-finite number at {path}"),
            Self::InvalidGraph {
                node: Some(node),
                error,
            } => write!(f, "invalid graph at node {node}: {error}"),
            Self::InvalidGraph { node: None, error } => write!(f, "invalid graph roots: {error}"),
            Self::Malformed { detail } => write!(f, "malformed payload: {detail}"),
        }
    }
}

impl std::error::Error for WireError {}

#[derive(Serialize)]
struct EnvelopeOut<'a> {
    format: &'static str,
    version: String,
    graph: GraphOut<'a>,
}

#[derive(Serialize)]
struct GraphOut<'a> {
    nodes: &'a [GeometryNode],
    roots: &'a [NodeId],
}

/// The envelope's identifying entries, read before the graph so that the
/// version is checked before any vocabulary is.
#[derive(Deserialize)]
struct Header {
    format: Option<String>,
    version: Option<String>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct EnvelopeIn {
    #[allow(dead_code)]
    format: String,
    #[allow(dead_code)]
    version: String,
    graph: GraphIn,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GraphIn {
    nodes: Vec<GeometryNode>,
    roots: Vec<NodeId>,
}

fn envelope(graph: &GeometryGraph) -> EnvelopeOut<'_> {
    EnvelopeOut {
        format: FORMAT_NAME,
        version: FORMAT_VERSION.to_string(),
        graph: GraphOut {
            nodes: graph.nodes(),
            roots: graph.roots(),
        },
    }
}

fn malformed(error: impl fmt::Display) -> WireError {
    WireError::Malformed {
        detail: error.to_string(),
    }
}

/// The value tree of `graph`'s payload, refused if it holds a non-finite
/// number.
fn checked_value(graph: &GeometryGraph) -> Result<Value, WireError> {
    let value = Value::serialized(&envelope(graph)).map_err(malformed)?;
    check_finite(&value, &mut String::new())?;
    Ok(value)
}

fn check_finite(value: &Value, path: &mut String) -> Result<(), WireError> {
    let at = path.len();
    let result = match value {
        Value::Float(number) if !number.is_finite() => Err(WireError::NonFinite {
            path: if path.is_empty() {
                "the payload root".to_owned()
            } else {
                path.clone()
            },
        }),
        Value::Array(items) => items.iter().enumerate().try_for_each(|(index, item)| {
            path.truncate(at);
            path.push_str(&format!("[{index}]"));
            check_finite(item, path)
        }),
        Value::Map(entries) => entries.iter().try_for_each(|(key, item)| {
            path.truncate(at);
            match key {
                Value::Text(name) if at == 0 => path.push_str(name),
                Value::Text(name) => {
                    path.push('.');
                    path.push_str(name);
                }
                _ => {
                    check_finite(key, path)?;
                    path.truncate(at);
                    path.push_str("[key]");
                }
            }
            check_finite(item, path)
        }),
        Value::Tag(_, inner) => check_finite(inner, path),
        _ => Ok(()),
    };
    path.truncate(at);
    result
}

/// The version a payload's header states, if this reader reads it.
fn check_header(header: Header) -> Result<FormatVersion, WireError> {
    match header.format {
        Some(format) if format == FORMAT_NAME => {}
        found => return Err(WireError::UnknownFormat { found }),
    }
    let found = header.version.ok_or(WireError::MissingVersion)?;
    match FormatVersion::parse(&found) {
        Some(version) if FORMAT_VERSION.reads(version) => Ok(version),
        _ => Err(WireError::UnsupportedVersion {
            found,
            supported: FORMAT_VERSION,
        }),
    }
}

/// The name serde quotes after `prefix` in a decoder message, up to the
/// clause that follows it.
fn quoted_after(message: &str, prefix: &str) -> Option<String> {
    let start = message.find(prefix)? + prefix.len();
    let rest = &message[start..];
    let end = ["`, expected", "`, there are no"]
        .iter()
        .filter_map(|clause| rest.rfind(clause))
        .max()
        .or_else(|| rest.find('`'))?;
    Some(rest[..end].to_owned())
}

/// Name a decoder failure: an unknown kind or field by its name, anything
/// else as malformed.
fn decode_error(message: String, version: FormatVersion) -> WireError {
    if let Some(kind) = quoted_after(&message, "unknown variant `") {
        return WireError::UnknownKind { kind, version };
    }
    if let Some(field) = quoted_after(&message, "unknown field `") {
        return WireError::UnknownField { field, version };
    }
    WireError::Malformed { detail: message }
}

/// A JSON decoder failure: a number literal past the double range is a
/// non-finite number.
fn json_error(error: &serde_json::Error, version: Option<FormatVersion>) -> WireError {
    let message = error.to_string();
    if message.starts_with("number out of range") {
        return WireError::NonFinite {
            path: format!("line {} column {}", error.line(), error.column()),
        };
    }
    match version {
        Some(version) => decode_error(message, version),
        None => WireError::Malformed { detail: message },
    }
}

/// Rebuild the decoded graph through the builder's validation.
fn build(owner: GraphId, graph: GraphIn) -> Result<GeometryGraph, WireError> {
    let mut builder = GeometryGraphBuilder::with_owner(owner);
    for (index, node) in graph.nodes.into_iter().enumerate() {
        builder
            .push(node)
            .map_err(|error| WireError::InvalidGraph {
                node: Some(index),
                error,
            })?;
    }
    builder
        .finish(graph.roots)
        .map_err(|error| WireError::InvalidGraph { node: None, error })
}

/// Write `graph` as a JSON payload.
///
/// # Errors
///
/// [`WireError::NonFinite`] if the graph holds a NaN or an infinity.
pub fn to_json(graph: &GeometryGraph) -> Result<String, WireError> {
    checked_value(graph)?;
    serde_json::to_string(&envelope(graph)).map_err(malformed)
}

/// Read a graph from a JSON payload.
///
/// # Errors
///
/// A named [`WireError`] for every payload this reader cannot read whole;
/// see the [module documentation](self).
pub fn from_json(text: &str) -> Result<GeometryGraph, WireError> {
    let header: Header = serde_json::from_str(text).map_err(|error| json_error(&error, None))?;
    let version = check_header(header)?;
    let owner = GraphId::fresh();
    let envelope: EnvelopeIn = {
        let _brand = WireBrand::enter(owner);
        serde_json::from_str(text).map_err(|error| json_error(&error, Some(version)))?
    };
    build(owner, envelope.graph)
}

/// Write `graph` as a CBOR payload.
///
/// # Errors
///
/// [`WireError::NonFinite`] if the graph holds a NaN or an infinity.
pub fn to_cbor(graph: &GeometryGraph) -> Result<Vec<u8>, WireError> {
    let value = checked_value(graph)?;
    let mut bytes = Vec::new();
    ciborium::into_writer(&value, &mut bytes).map_err(malformed)?;
    Ok(bytes)
}

/// Read a graph from a CBOR payload: one data item, nothing after it.
///
/// # Errors
///
/// A named [`WireError`] for every payload this reader cannot read whole;
/// see the [module documentation](self).
pub fn from_cbor(bytes: &[u8]) -> Result<GeometryGraph, WireError> {
    let mut rest = bytes;
    let value: Value = ciborium::from_reader(&mut rest).map_err(malformed)?;
    if !rest.is_empty() {
        return Err(WireError::Malformed {
            detail: format!("{} trailing bytes after the payload", rest.len()),
        });
    }
    let header: Header = value.deserialized().map_err(malformed)?;
    let version = check_header(header)?;
    check_finite(&value, &mut String::new())?;
    let owner = GraphId::fresh();
    let envelope: EnvelopeIn = {
        let _brand = WireBrand::enter(owner);
        value
            .deserialized()
            .map_err(|error| decode_error(error.to_string(), version))?
    };
    build(owner, envelope.graph)
}

impl GeometryGraph {
    /// This graph as a JSON payload of the wire format; see [`to_json`].
    ///
    /// # Errors
    ///
    /// [`WireError::NonFinite`] if the graph holds a NaN or an infinity.
    pub fn to_json(&self) -> Result<String, WireError> {
        to_json(self)
    }

    /// A graph read from a JSON payload of the wire format; see
    /// [`from_json`].
    ///
    /// # Errors
    ///
    /// A named [`WireError`] for every payload that cannot be read whole.
    pub fn from_json(text: &str) -> Result<Self, WireError> {
        from_json(text)
    }

    /// This graph as a CBOR payload of the wire format; see [`to_cbor`].
    ///
    /// # Errors
    ///
    /// [`WireError::NonFinite`] if the graph holds a NaN or an infinity.
    pub fn to_cbor(&self) -> Result<Vec<u8>, WireError> {
        to_cbor(self)
    }

    /// A graph read from a CBOR payload of the wire format; see
    /// [`from_cbor`].
    ///
    /// # Errors
    ///
    /// A named [`WireError`] for every payload that cannot be read whole.
    pub fn from_cbor(bytes: &[u8]) -> Result<Self, WireError> {
        from_cbor(bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn versions_parse_strictly() {
        assert_eq!(FormatVersion::parse("1.0"), Some(FormatVersion::new(1, 0)));
        assert_eq!(
            FormatVersion::parse("12.34"),
            Some(FormatVersion::new(12, 34))
        );
        for bad in ["1", "1.", ".0", "1.0.0", "+1.0", "1.-0", "v1.0", " 1.0", ""] {
            assert_eq!(FormatVersion::parse(bad), None, "{bad:?}");
        }
    }

    #[test]
    fn a_reader_reads_its_own_major_up_to_its_own_minor() {
        let reader = FormatVersion::new(1, 2);
        assert!(reader.reads(FormatVersion::new(1, 0)));
        assert!(reader.reads(FormatVersion::new(1, 2)));
        assert!(!reader.reads(FormatVersion::new(1, 3)));
        assert!(!reader.reads(FormatVersion::new(2, 0)));
        assert!(!reader.reads(FormatVersion::new(0, 9)));
    }

    #[test]
    fn decoder_messages_name_the_unknown_kind_or_field() {
        let version = FORMAT_VERSION;
        assert_eq!(
            decode_error(
                "unknown variant `Hyperboloid`, expected one of `Plane`, `Cylinder` at line 1 column 9"
                    .to_owned(),
                version
            ),
            WireError::UnknownKind {
                kind: "Hyperboloid".to_owned(),
                version
            }
        );
        assert_eq!(
            decode_error(
                "unknown field `colour`, expected `radius`".to_owned(),
                version
            ),
            WireError::UnknownField {
                field: "colour".to_owned(),
                version
            }
        );
        assert_eq!(
            decode_error(
                "unknown variant `A`b`, there are no variants".to_owned(),
                version
            ),
            WireError::UnknownKind {
                kind: "A`b".to_owned(),
                version
            }
        );
        assert!(matches!(
            decode_error("invalid type: string, expected f64".to_owned(), version),
            WireError::Malformed { .. }
        ));
    }
}
