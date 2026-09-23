use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

pub type Result<T> = std::result::Result<T, String>;

pub const LEDGER: &str = "architecture/capability-ledger.toml";
pub const PACKAGES: &str = "architecture/reference-packages.toml";

/// Every package path that exists in the pinned OCCT and CGAL trees.
///
/// Generated once per pin from the trees themselves, so the gate can refuse
/// a mistyped reference path without either tree being present.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferencePackages {
    pub occt: std::collections::BTreeSet<String>,
    pub cgal: std::collections::BTreeSet<String>,
}

pub fn load_packages(root: &Path) -> Result<ReferencePackages> {
    let path = root.join(PACKAGES);
    let text = std::fs::read_to_string(&path)
        .map_err(|error| format!("read {}: {error}", path.display()))?;
    toml::from_str(&text).map_err(|error| format!("{PACKAGES}: {error}"))
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Ledger {
    pub reference: BTreeMap<String, Reference>,
    #[serde(default)]
    pub issue: Vec<Issue>,
    pub capability: Vec<Capability>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Reference {
    pub repository: String,
    #[serde(default)]
    pub tag: Option<String>,
    pub commit: String,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Issue {
    pub key: String,
    pub title: String,
    pub priority: String,
    pub effort: String,
    /// GitHub issue number; 0 until filed.
    pub number: u32,
    /// Issue keys that must close first (mirrors GitHub "blocked by").
    #[serde(default)]
    pub blocked_by: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Level {
    Absent,
    Narrow,
    Implemented,
}

impl Level {
    pub fn as_str(self) -> &'static str {
        match self {
            Level::Absent => "absent",
            Level::Narrow => "narrow",
            Level::Implemented => "implemented",
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Capability {
    pub id: String,
    pub area: String,
    pub name: String,
    pub level: Level,
    pub summary: String,
    #[serde(default)]
    pub evidence: Vec<String>,
    #[serde(default)]
    pub occt: Vec<String>,
    #[serde(default)]
    pub cgal: Vec<String>,
    #[serde(default)]
    pub issue: Option<String>,
}

impl Ledger {
    pub fn issue(&self, key: &str) -> Option<&Issue> {
        self.issue.iter().find(|issue| issue.key == key)
    }

    /// An issue by key, or by GitHub number with or without a leading `#`.
    pub fn find_issue(&self, target: &str) -> Option<&Issue> {
        let number = target.trim_start_matches('#').parse::<u32>().ok();
        self.issue
            .iter()
            .find(|issue| issue.key == target || (number.is_some() && Some(issue.number) == number))
    }
}

pub fn workspace_root() -> Result<PathBuf> {
    // xtask lives at tools/xtask; the workspace root is two levels up.
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    manifest
        .ancestors()
        .nth(2)
        .map(Path::to_path_buf)
        .ok_or_else(|| "cannot locate workspace root".to_owned())
}

pub fn load(root: &Path) -> Result<Ledger> {
    let path = root.join(LEDGER);
    let text = std::fs::read_to_string(&path)
        .map_err(|error| format!("read {}: {error}", path.display()))?;
    parse(&text)
}

pub fn parse(text: &str) -> Result<Ledger> {
    toml::from_str(text).map_err(|error| format!("{LEDGER}: {error}"))
}
