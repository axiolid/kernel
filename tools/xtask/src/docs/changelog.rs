//! Release sections of the per-crate changelogs (ADR 0067), and
//! `docs/reference/changelog.md`, which collects them.
//!
//! A release is a dated `## [x.y.z] - YYYY-MM-DD` section. The release commit
//! dates it, so what a page calls "released" is a property of the commit, not
//! of when the docs are built: reading git tags instead would make every
//! release turn `main`'s docs stale, since the tag is pushed after the
//! release commit merges. `[Unreleased]` is never a release, so adding to it
//! never makes a generated page stale.
//!
//! Crates published only as part of a workspace-wide release before ADR 0067
//! have no dated section of their own. For them the release is the dated
//! section of `docs/CHANGELOG.md` whose version the manifest still carries.

use std::cmp::Ordering;

use super::workspace::{Crate, Workspace};
use super::{Output, Result, BANNER, REPO};

const PAGE: &str = "docs/reference/changelog.md";
const WORKSPACE_CHANGELOG: &str = "docs/CHANGELOG.md";

/// One dated release section.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Section {
    pub(super) version: String,
    pub(super) date: String,
    pub(super) body: String,
}

/// What a crate has released, and where the notes for it live.
pub(super) struct Release {
    pub(super) version: String,
    pub(super) date: String,
    /// The notes of this release, links made absolute for the site. `None`
    /// when the release was workspace-wide and the crate has no own section.
    pub(super) notes: Option<String>,
}

fn read(workspace: &Workspace, rel: &str) -> Result<String> {
    std::fs::read_to_string(workspace.root.join(rel)).map_err(|error| format!("{rel}: {error}"))
}

fn crate_changelog(workspace: &Workspace, krate: &Crate) -> Result<String> {
    read(workspace, &format!("{}/CHANGELOG.md", krate.dir))
}

/// The newest release of a crate, if it has one.
pub(super) fn latest(workspace: &Workspace, krate: &Crate) -> Result<Option<Release>> {
    let own = sections(&crate_changelog(workspace, krate)?);
    if let Some(newest) = own
        .into_iter()
        .max_by(|a, b| compare_versions(&a.version, &b.version))
    {
        let notes = (!newest.body.is_empty()).then(|| absolutise(&newest.body, &krate.dir));
        return Ok(Some(Release {
            version: newest.version,
            date: newest.date,
            notes,
        }));
    }
    let shared = sections(&read(workspace, WORKSPACE_CHANGELOG)?);
    Ok(shared
        .into_iter()
        .find(|section| section.version == krate.version)
        .map(|section| Release {
            version: section.version,
            date: section.date,
            notes: None,
        }))
}

/// `docs/reference/changelog.md`: every publishable crate's releases,
/// grouped by crate, newest first within each.
pub(super) fn page(workspace: &Workspace) -> Result<Output> {
    let mut crates = Vec::new();
    for krate in &workspace.crates {
        crates.push((krate.name.as_str(), crate_changelog(workspace, krate)?));
    }
    Ok(Output::new(workspace, PAGE, render(&crates)))
}

fn render(crates: &[(&str, String)]) -> String {
    let mut lines = vec![
        BANNER.to_owned(),
        String::new(),
        "# Per-crate changelog".to_owned(),
        String::new(),
        "Every publishable crate versions and publishes independently \
         ([ADR 0067](/adr/0067-crates-version-independently)); this page \
         collects each crate's own `CHANGELOG.md`, newest release first per \
         crate. Workspace-wide narrative — breaking bumps and coordinated \
         releases — stays in the [top-level changelog](/CHANGELOG)."
            .to_owned(),
        String::new(),
    ];
    let mut any_release = false;
    for (name, text) in crates {
        let sections = sections(text);
        if sections.is_empty() {
            continue;
        }
        any_release = true;
        lines.push(format!("## {name}"));
        lines.push(String::new());
        for section in sections {
            lines.push(format!("### {} - {}", section.version, section.date));
            lines.push(String::new());
            if !section.body.is_empty() {
                lines.push(section.body);
                lines.push(String::new());
            }
        }
        lines.push(String::new());
    }
    if !any_release {
        lines.push(
            "No crate has a dated release yet under independent versioning; \
             every crate's history to date lives in the top-level changelog."
                .to_owned(),
        );
        lines.push(String::new());
    }
    format!("{}\n", lines.join("\n").trim_end_matches('\n'))
}

/// Every dated section in order of appearance. A section's body runs to the
/// next dated heading; `[Unreleased]` is not one, so its notes, which come
/// first, belong to no section.
pub(super) fn sections(text: &str) -> Vec<Section> {
    // (heading start, body start, version, date)
    let mut headings = Vec::new();
    let mut offset = 0;
    for line in text.split_inclusive('\n') {
        if let Some((version, date)) = dated_heading(line) {
            headings.push((offset, offset + line.len(), version, date));
        }
        offset += line.len();
    }
    let mut out = Vec::new();
    for (index, (_, body_start, version, date)) in headings.iter().enumerate() {
        let end = headings.get(index + 1).map_or(text.len(), |next| next.0);
        out.push(Section {
            version: version.clone(),
            date: date.clone(),
            body: text[*body_start..end].trim_matches('\n').to_owned(),
        });
    }
    out
}

/// `## [0.3.1] - 2026-09-27` → (version, date). Only a dated heading counts.
fn dated_heading(line: &str) -> Option<(String, String)> {
    let rest = line.strip_prefix("## [")?;
    let close = rest.find(']')?;
    let version = &rest[..close];
    if version.is_empty() {
        return None;
    }
    let date = rest[close + 1..].strip_prefix(" - ")?.trim_end();
    let shape = date.len() == 10
        && date.bytes().enumerate().all(|(i, b)| match i {
            4 | 7 => b == b'-',
            _ => b.is_ascii_digit(),
        });
    shape.then(|| (version.to_owned(), date.to_owned()))
}

/// Versions compared numerically, piece by piece.
pub(super) fn compare_versions(a: &str, b: &str) -> Ordering {
    let key = |version: &str| -> Vec<u64> {
        version
            .split(['.', '-', '+'])
            .map(|piece| piece.parse().unwrap_or(0))
            .collect()
    };
    key(a).cmp(&key(b))
}

/// Rewrite repository-relative `[text](target)` links to GitHub URLs.
///
/// A changelog is read in the repository, where `../x.md` resolves from the
/// crate directory, and quoted on the site, where it does not and VitePress
/// fails the build on the dead link. Site links (`/adr/…`), anchors and URLs
/// stay as they are.
pub(super) fn absolutise(body: &str, dir: &str) -> String {
    let mut out = String::with_capacity(body.len());
    let mut rest = body;
    while let Some(open) = rest.find('[') {
        out.push_str(&rest[..open]);
        match link_at(&rest[open..]) {
            Some((text, target, length)) => {
                let keep = target.contains("://")
                    || target.starts_with('#')
                    || target.starts_with('/')
                    || target.starts_with("mailto:");
                if keep {
                    out.push_str(&rest[open..open + length]);
                } else {
                    out.push_str(&format!("[{text}]({REPO}/blob/main/{})", join(dir, target)));
                }
                rest = &rest[open + length..];
            }
            None => {
                out.push('[');
                rest = &rest[open + 1..];
            }
        }
    }
    out.push_str(rest);
    out
}

/// `dir/target` with `.` and `..` resolved.
fn join(dir: &str, target: &str) -> String {
    let mut parts: Vec<&str> = dir.split('/').filter(|p| !p.is_empty()).collect();
    for piece in target.split('/') {
        match piece {
            "" | "." => {}
            ".." => {
                parts.pop();
            }
            _ => parts.push(piece),
        }
    }
    parts.join("/")
}

/// `[text](target)` at the start of `s`: (text, target, matched length).
fn link_at(s: &str) -> Option<(&str, &str, usize)> {
    let inner = s.strip_prefix('[')?;
    let close = inner.find(']')?;
    if close == 0 {
        return None;
    }
    let after = inner[close + 1..].strip_prefix('(')?;
    let end = after.find(')')?;
    if end == 0 {
        return None;
    }
    Some((&inner[..close], &after[..end], 1 + close + 2 + end + 1))
}

#[cfg(test)]
mod tests {
    use super::*;

    const TWO_RELEASES: &str = "# Changelog\n\n\
        ## [Unreleased]\n\n- pending, must not appear\n\n\
        ## [0.3.1] - 2026-09-20\n\n### Added\n\n- b\n\n\
        ## [0.3.0] - 2026-09-01\n\n### Added\n\n- a\n";

    #[test]
    fn sections_skip_unreleased_and_keep_order() {
        let found = sections(TWO_RELEASES);
        let versions: Vec<&str> = found.iter().map(|s| s.version.as_str()).collect();
        let dates: Vec<&str> = found.iter().map(|s| s.date.as_str()).collect();
        assert_eq!(versions, ["0.3.1", "0.3.0"]);
        assert_eq!(dates, ["2026-09-20", "2026-09-01"]);
        assert_eq!(found[0].body, "### Added\n\n- b");
        assert!(found.iter().all(|s| !s.body.contains("pending")));
    }

    #[test]
    fn only_unreleased_means_no_sections() {
        assert!(sections("# Changelog\n\n## [Unreleased]\n\n- pending\n").is_empty());
    }

    #[test]
    fn headings_need_a_date() {
        assert_eq!(
            dated_heading("## [0.2.1] - 2026-09-23\n"),
            Some(("0.2.1".into(), "2026-09-23".into()))
        );
        assert_eq!(dated_heading("## [Unreleased]"), None);
        assert_eq!(dated_heading("## [0.1.0] - soon"), None);
        assert_eq!(dated_heading("### [0.1.0] - 2026-09-23"), None);
    }

    #[test]
    fn page_skips_crates_without_a_release() {
        let unreleased = "# Changelog\n\n## [Unreleased]\n\n- pending\n".to_owned();
        let released =
            "# Changelog\n\n## [Unreleased]\n\n## [0.3.1] - 2026-09-20\n\n- shipped\n".to_owned();
        let page = render(&[
            ("axiolid-unreleased", unreleased.clone()),
            ("axiolid-released", released),
        ]);
        assert!(!page.contains("axiolid-unreleased"));
        assert!(page.contains("## axiolid-released\n\n### 0.3.1 - 2026-09-20\n\n- shipped\n"));
        assert!(page.starts_with(BANNER));

        let empty = render(&[("axiolid-probe", unreleased)]);
        assert!(empty.contains("No crate has a dated release yet"));
    }

    #[test]
    fn versions_compare_numerically() {
        assert_eq!(compare_versions("0.10.0", "0.9.0"), Ordering::Greater);
        assert_eq!(compare_versions("0.2.0", "0.2.0"), Ordering::Equal);
    }

    #[test]
    fn relative_links_resolve_against_the_crate() {
        assert_eq!(
            absolutise(
                "see [adr](../../docs/x.md), [site](/adr/1), [u](https://a.b) and [y](#z)",
                "crates/a/b"
            ),
            format!(
                "see [adr]({REPO}/blob/main/crates/docs/x.md), [site](/adr/1), [u](https://a.b) and [y](#z)"
            )
        );
        assert_eq!(absolutise("[a] (b) [] (c)", "x"), "[a] (b) [] (c)");
    }
}
