//! Backend identity metadata. Operation traits remain the sole capability truth.

use core::fmt;

/// Maximum identifier length in bytes.
///
/// Sized so driver-enumerated accelerator names (`cuda:0`, `hip:1`,
/// `vulkan:discrete:0`) fit without allocating, while keeping [`BackendId`]
/// small enough to stay `Copy` inside every error variant.
const IDENTIFIER_CAPACITY: usize = 47;

/// Stable provider identifier for logs and explicit selection.
///
/// Identity is stored inline as fixed-capacity UTF-8 rather than
/// `&'static str`, because accelerator backends (CUDA, HIP, Vulkan) enumerate
/// their devices at runtime and cannot produce `'static` text without leaking.
/// Keeping the value inline preserves `Copy`, so `BackendId` can continue to
/// live inside `GeomError` variants and `DevicePreference` without forcing an
/// allocation or a lifetime onto the error type.
#[derive(Clone, Copy)]
pub struct BackendId {
    bytes: [u8; IDENTIFIER_CAPACITY],
    len: u8,
}

/// An identifier was longer than [`BackendId::CAPACITY`] bytes.
///
/// Rejecting is deliberate: a silently truncated identity would make two
/// distinct devices compare equal, which would corrupt provider selection and
/// error attribution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BackendIdTooLong {
    /// Length of the rejected identifier, in bytes.
    pub len: usize,
}

impl fmt::Display for BackendIdTooLong {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "backend identifier is {} bytes, limit is {}",
            self.len,
            BackendId::CAPACITY
        )
    }
}

impl core::error::Error for BackendIdTooLong {}

impl BackendId {
    /// Maximum identifier length in bytes.
    pub const CAPACITY: usize = IDENTIFIER_CAPACITY;

    /// Construct an identifier known at compile time.
    ///
    /// # Panics
    ///
    /// Panics if `value` exceeds [`BackendId::CAPACITY`] bytes. This is a
    /// `const fn`, so an over-long literal fails the build rather than a test
    /// run. Use [`BackendId::try_new`] for runtime-derived text.
    pub const fn new(value: &str) -> Self {
        match Self::build(value.as_bytes()) {
            Some(id) => id,
            None => panic!("backend identifier exceeds BackendId::CAPACITY"),
        }
    }

    /// Construct an identifier from runtime-derived text, such as a
    /// driver-enumerated device name.
    ///
    /// # Errors
    ///
    /// Returns [`BackendIdTooLong`] when the text exceeds
    /// [`BackendId::CAPACITY`] bytes. The identity is never truncated.
    pub fn try_new(value: &str) -> Result<Self, BackendIdTooLong> {
        Self::build(value.as_bytes()).ok_or(BackendIdTooLong { len: value.len() })
    }

    /// Shared inline copy used by both constructors.
    ///
    /// Written as an index loop because `copy_from_slice` is not `const`.
    const fn build(source: &[u8]) -> Option<Self> {
        if source.len() > IDENTIFIER_CAPACITY {
            return None;
        }
        let mut bytes = [0_u8; IDENTIFIER_CAPACITY];
        let mut index = 0;
        while index < source.len() {
            bytes[index] = source[index];
            index += 1;
        }
        Some(Self {
            bytes,
            len: source.len() as u8,
        })
    }

    /// Identifier text.
    pub fn as_str(&self) -> &str {
        // The only constructors copy from a `&str`, so the populated prefix is
        // always valid UTF-8. Zero padding beyond `len` is never included.
        core::str::from_utf8(&self.bytes[..self.len as usize])
            .expect("identifier bytes originate from &str")
    }
}

// Comparison and hashing use the semantic prefix, never the zero padding, so a
// runtime-built identity equals a compile-time one with the same text.
impl PartialEq for BackendId {
    fn eq(&self, other: &Self) -> bool {
        self.as_str() == other.as_str()
    }
}

impl Eq for BackendId {}

impl PartialOrd for BackendId {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for BackendId {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.as_str().cmp(other.as_str())
    }
}

impl core::hash::Hash for BackendId {
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.as_str().hash(state);
    }
}

impl fmt::Debug for BackendId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "BackendId({:?})", self.as_str())
    }
}

impl fmt::Display for BackendId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Broad execution target. Specific ISA/device features stay in provider crates.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ExecutionTarget {
    /// Portable scalar CPU implementation.
    PortableCpu,
    /// Runtime-selected CPU implementation.
    OptimizedCpu,
    /// General-purpose GPU compute.
    Gpu,
    /// Other accelerator supplied downstream.
    Accelerator,
}

/// Operation name used for diagnostics only. Implementing an operation trait is
/// the capability proof; this enum never drives capability discovery.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Operation {
    CurveEvaluation,
    SurfaceEvaluation,
    ProfileTriangulation,
    Sweep,
    Tessellation,
    MeshBoolean,
    MeshPlaneSection,
    PointcloudReconstruction,
    SpatialQuery,
    Measurement,
    Healing,
    GraphCompilation,
}

/// Provider identity only. It deliberately contains no operation booleans.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BackendDescriptor {
    /// Stable implementation identity.
    pub id: BackendId,
    /// Hardware class used by execution policy.
    pub target: ExecutionTarget,
}

impl BackendDescriptor {
    /// Construct provider identity metadata.
    pub const fn new(id: BackendId, target: ExecutionTarget) -> Self {
        Self { id, target }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A driver-enumerated accelerator identity (`cuda:0`, `hip:1`) is only
    /// known at runtime, so identifiers must not require `&'static str`.
    #[test]
    fn identifiers_accept_runtime_owned_text() {
        let ordinal = 1_u32;
        let runtime = format!("cuda:{ordinal}");
        let id = BackendId::try_new(&runtime).expect("short runtime identity");
        assert_eq!(id.as_str(), "cuda:1");
        assert_eq!(id.to_string(), "cuda:1");
    }

    #[test]
    fn runtime_and_const_identifiers_compare_equal() {
        const STATIC: BackendId = BackendId::new("cuda:0");
        let runtime = BackendId::try_new(&String::from("cuda:0")).expect("identity");
        assert_eq!(STATIC, runtime);
        assert_eq!(STATIC.cmp(&runtime), core::cmp::Ordering::Equal);
    }

    #[test]
    fn identifiers_stay_copy_and_orderable() {
        fn assert_copy<T: Copy + Ord + core::hash::Hash>() {}
        assert_copy::<BackendId>();
        let a = BackendId::try_new("aa").expect("identity");
        let b = BackendId::try_new("aab").expect("identity");
        assert!(a < b, "zero padding must not invert lexicographic order");
    }

    #[test]
    fn over_long_identifiers_are_rejected_not_truncated() {
        let too_long = "x".repeat(BackendId::CAPACITY + 1);
        assert!(BackendId::try_new(&too_long).is_err());
        let at_limit = "x".repeat(BackendId::CAPACITY);
        assert_eq!(
            BackendId::try_new(&at_limit).expect("identity").as_str(),
            at_limit
        );
    }

    /// Two devices whose names share a prefix must stay distinct. A truncating
    /// implementation would alias them and misroute both selection and blame.
    #[test]
    fn long_shared_prefix_devices_do_not_alias() {
        let base = "x".repeat(BackendId::CAPACITY - 1);
        let first = BackendId::try_new(&format!("{base}0")).expect("identity");
        let second = BackendId::try_new(&format!("{base}1")).expect("identity");
        assert_ne!(first, second);
    }

    #[test]
    fn hashing_matches_equality_across_construction_paths() {
        use std::collections::HashSet;
        let mut seen = HashSet::new();
        seen.insert(BackendId::new("hip:0"));
        assert!(seen.contains(&BackendId::try_new("hip:0").expect("identity")));
    }

    #[test]
    fn multibyte_identifiers_round_trip() {
        let id = BackendId::try_new("gpu-µ-0").expect("identity");
        assert_eq!(id.as_str(), "gpu-µ-0");
    }
}
