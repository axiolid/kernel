//! Cache and core-topology detection for tuning, not for capability gating.
//!
//! # Why this is separate from `CpuFeatures`
//!
//! `CpuFeatures` answers *what a machine can execute* -- a wrong answer there
//! is a crash. This module answers *what shape the machine is* -- a wrong
//! answer here is only a bad tuning choice. The two failure modes differ by
//! orders of magnitude, so they stay separate types.
//!
//! # Why every field is optional
//!
//! Detection reads sysfs, which is absent on some targets and restricted in
//! some containers. An unknown cache size is reported as `None`, never as a
//! plausible default: a fabricated 32 KiB that is really 64 KiB produces a
//! tuning decision made on a number nobody measured. Callers must supply
//! their own fallback explicitly, so the guess is visible at the call site.
//!
//! # Evidence this matters
//!
//! A radix sort in the mesh audit executed 18.9% fewer instructions than the
//! comparison sort it replaced and ran 2.2x slower, because scattering into
//! 256 buckets missed L1 four times as often. Instruction count could not see
//! it; cache size explains it. That is the class of decision this module
//! exists to inform.

/// One cache level's measured geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CacheLevel {
    /// Level number: 1, 2, or 3.
    pub level: u8,
    /// Total size in bytes.
    pub bytes: usize,
    /// Cache line size in bytes.
    pub line_bytes: usize,
    /// How many logical CPUs share this level.
    ///
    /// A shared last level means per-thread working sets contend; a private
    /// L2 means they do not. Work-splitting decisions need this, not just
    /// the size.
    pub shared_by: usize,
}

/// The measured shape of the host machine.
///
/// Every field is optional because an unmeasured value must stay visibly
/// absent rather than silently defaulted.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct CpuTopology {
    /// Logical CPUs usable by this process, honouring CPU affinity.
    pub logical_cpus: Option<usize>,
    /// Data and unified cache levels, ascending by level.
    pub caches: Vec<CacheLevel>,
    /// Whether cores have differing capability (P/E cores, Snapdragon X).
    ///
    /// `None` means undetermined, which is not the same as `Some(false)`.
    /// Even work splitting is wrong on a heterogeneous machine, so the
    /// distinction matters.
    pub heterogeneous_cores: Option<bool>,
}

impl CpuTopology {
    /// Detect the host's shape. Unavailable values stay `None`.
    pub fn detect() -> Self {
        Self {
            logical_cpus: detect_logical_cpus(),
            caches: detect_caches(),
            heterogeneous_cores: detect_heterogeneous(),
        }
    }

    /// Bytes in the given cache level, if measured.
    pub fn cache_bytes(&self, level: u8) -> Option<usize> {
        self.caches
            .iter()
            .find(|cache| cache.level == level)
            .map(|cache| cache.bytes)
    }

    /// Bytes in the largest measured cache level.
    pub fn last_level_bytes(&self) -> Option<usize> {
        self.caches.iter().map(|cache| cache.bytes).max()
    }

    /// Cache line size, if measured.
    pub fn line_bytes(&self) -> Option<usize> {
        self.caches.first().map(|cache| cache.line_bytes)
    }

    /// Whether a working set of `bytes` fits in the given level.
    ///
    /// `None` when the level was not measured -- the caller must decide what
    /// to do without the number rather than receive a fabricated `false`.
    pub fn fits_in_cache(&self, bytes: usize, level: u8) -> Option<bool> {
        self.cache_bytes(level).map(|capacity| bytes <= capacity)
    }
}

/// Logical CPUs available to this process.
///
/// `available_parallelism` honours cgroup limits and CPU affinity, so a
/// container pinned to 2 cores of a 20-core host reports 2. Reading the raw
/// core count would oversubscribe such a machine badly.
fn detect_logical_cpus() -> Option<usize> {
    std::thread::available_parallelism().ok().map(Into::into)
}

/// Parse a sysfs size string such as `32K`, `4096K`, or `16M`.
fn parse_size(text: &str) -> Option<usize> {
    let text = text.trim();
    let (digits, multiplier) = match text.as_bytes().last()? {
        b'K' => (&text[..text.len() - 1], 1024),
        b'M' => (&text[..text.len() - 1], 1024 * 1024),
        b'G' => (&text[..text.len() - 1], 1024 * 1024 * 1024),
        _ => (text, 1),
    };
    digits.parse::<usize>().ok()?.checked_mul(multiplier)
}

/// Count entries in a sysfs CPU list such as `0-19` or `0,4,8`.
fn parse_cpu_list(text: &str) -> Option<usize> {
    let mut total = 0usize;
    for part in text.trim().split(',') {
        if part.is_empty() {
            continue;
        }
        match part.split_once('-') {
            Some((low, high)) => {
                let low: usize = low.trim().parse().ok()?;
                let high: usize = high.trim().parse().ok()?;
                total += high.checked_sub(low)?.checked_add(1)?;
            }
            None => total += 1,
        }
    }
    (total > 0).then_some(total)
}

/// Read data and unified cache levels from sysfs.
///
/// Instruction caches are skipped: they say nothing about a data working
/// set, which is what tuning decisions are about.
#[cfg(target_os = "linux")]
fn detect_caches() -> Vec<CacheLevel> {
    use std::fs::read_to_string;

    let mut caches = Vec::new();
    for index in 0..16 {
        let base = format!("/sys/devices/system/cpu/cpu0/cache/index{index}");
        let Ok(kind) = read_to_string(format!("{base}/type")) else {
            break;
        };
        let kind = kind.trim();
        if kind != "Data" && kind != "Unified" {
            continue;
        }
        let level = read_to_string(format!("{base}/level"))
            .ok()
            .and_then(|text| text.trim().parse::<u8>().ok());
        let bytes = read_to_string(format!("{base}/size"))
            .ok()
            .and_then(|text| parse_size(&text));
        let line_bytes = read_to_string(format!("{base}/coherency_line_size"))
            .ok()
            .and_then(|text| text.trim().parse::<usize>().ok());
        let shared_by = read_to_string(format!("{base}/shared_cpu_list"))
            .ok()
            .and_then(|text| parse_cpu_list(&text));
        // A partially-read level is dropped rather than completed with
        // invented values.
        if let (Some(level), Some(bytes), Some(line_bytes), Some(shared_by)) =
            (level, bytes, line_bytes, shared_by)
        {
            caches.push(CacheLevel {
                level,
                bytes,
                line_bytes,
                shared_by,
            });
        }
    }
    caches.sort_unstable_by_key(|cache| cache.level);
    caches.dedup_by_key(|cache| cache.level);
    caches
}

/// Cache geometry is unavailable off Linux; callers see an empty list.
#[cfg(not(target_os = "linux"))]
fn detect_caches() -> Vec<CacheLevel> {
    Vec::new()
}

/// Whether the machine has cores of differing capability.
///
/// Two independent signals, because neither alone covers both vendors:
///
/// - `cpu/types/` exists on Intel hybrid parts (P/E cores).
/// - Differing `cpu_capacity` values cover Arm big.LITTLE, which is how a
///   Snapdragon X presents its heterogeneous cores.
///
/// Returns `None` when neither signal is readable, because "undetermined"
/// and "homogeneous" lead to different splitting choices.
#[cfg(target_os = "linux")]
fn detect_heterogeneous() -> Option<bool> {
    use std::fs::read_to_string;

    if std::path::Path::new("/sys/devices/system/cpu/types").is_dir() {
        return Some(true);
    }

    let mut capacities = Vec::new();
    for cpu in 0..256 {
        let path = format!("/sys/devices/system/cpu/cpu{cpu}/cpu_capacity");
        if !std::path::Path::new(&path).exists() {
            // CPU numbering is dense from 0; the first gap ends the scan.
            if cpu == 0 {
                return None;
            }
            break;
        }
        if let Some(value) = read_to_string(&path)
            .ok()
            .and_then(|text| text.trim().parse::<u32>().ok())
        {
            capacities.push(value);
        }
    }
    if capacities.is_empty() {
        return None;
    }
    let first = capacities[0];
    Some(capacities.iter().any(|value| *value != first))
}

/// Core heterogeneity is unavailable off Linux.
#[cfg(not(target_os = "linux"))]
fn detect_heterogeneous() -> Option<bool> {
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizes_parse_with_and_without_suffixes() {
        assert_eq!(parse_size("32K"), Some(32 * 1024));
        assert_eq!(parse_size("4096K"), Some(4096 * 1024));
        assert_eq!(parse_size("16M"), Some(16 * 1024 * 1024));
        assert_eq!(parse_size("512"), Some(512));
        assert_eq!(parse_size(""), None);
        assert_eq!(parse_size("garbage"), None);
    }

    #[test]
    fn cpu_lists_count_ranges_and_singletons() {
        assert_eq!(parse_cpu_list("0-19"), Some(20));
        assert_eq!(parse_cpu_list("0"), Some(1));
        assert_eq!(parse_cpu_list("0,4,8"), Some(3));
        assert_eq!(parse_cpu_list("0-3,8-11"), Some(8));
        assert_eq!(parse_cpu_list(""), None);
    }

    /// A descending range must not silently wrap into a huge count.
    #[test]
    fn a_malformed_range_refuses_rather_than_wrapping() {
        assert_eq!(parse_cpu_list("19-0"), None);
    }

    /// An unmeasured level must report absence, never a plausible default.
    #[test]
    fn an_unmeasured_cache_level_answers_none() {
        let topology = CpuTopology::default();
        assert_eq!(topology.cache_bytes(1), None);
        assert_eq!(topology.last_level_bytes(), None);
        assert_eq!(topology.fits_in_cache(1024, 1), None);
    }

    #[test]
    fn fit_queries_compare_against_the_measured_capacity() {
        let topology = CpuTopology {
            logical_cpus: Some(8),
            caches: vec![CacheLevel {
                level: 1,
                bytes: 32 * 1024,
                line_bytes: 64,
                shared_by: 1,
            }],
            heterogeneous_cores: Some(false),
        };
        assert_eq!(topology.fits_in_cache(32 * 1024, 1), Some(true));
        assert_eq!(topology.fits_in_cache(32 * 1024 + 1, 1), Some(false));
        assert_eq!(topology.fits_in_cache(1024, 2), None);
    }

    /// Detection must never invent values on the host it runs on.
    ///
    /// This asserts internal consistency rather than specific sizes, so it
    /// stays true on any machine including CI containers.
    #[test]
    fn detection_is_self_consistent_on_this_host() {
        let topology = CpuTopology::detect();
        if let Some(cpus) = topology.logical_cpus {
            assert!(cpus >= 1, "a usable machine has at least one CPU");
        }
        let mut previous = 0usize;
        for cache in &topology.caches {
            assert!(cache.bytes > 0, "a reported cache has a size");
            assert!(cache.line_bytes > 0, "a reported cache has a line size");
            assert!(cache.shared_by >= 1, "a cache is shared by >= 1 CPU");
            assert!(
                cache.bytes >= previous,
                "cache levels grow with level: L{} is {} bytes after {previous}",
                cache.level,
                cache.bytes,
            );
            previous = cache.bytes;
        }
    }
}
