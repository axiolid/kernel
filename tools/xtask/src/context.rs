//! Guard where the repository keeps its context (ADR 0078).
//!
//! The root `AGENTS.md` is the one file every contributor reads first: short,
//! stable, and the only one of its name. Each crate carries a `README.md` (its
//! crates.io page) with its purpose and the design notes no test, ADR or
//! module doc already holds; the reasoning behind a module lives in its `//!`
//! docs. Open work is not context: it lives in GitHub issues and `TODO(#N)`
//! markers, never in a checked-in plan.
//!
//! `check` keeps that shape: a nested `AGENTS.md`, a plan file or a stray
//! root-level note cannot regrow, a crate cannot ship without a README, a
//! README cannot grow into a manual or a checklist, a code marker cannot
//! float free of an issue, and a pointer to a README cannot dangle after a
//! file moves.

use std::path::{Path, PathBuf};
use std::process::Command;

use cargo_metadata::MetadataCommand;

type Result<T> = std::result::Result<T, String>;

/// The root `AGENTS.md` is read in full before any change, so it stays short.
const ROOT_AGENTS_MAX_LINES: usize = 120;

/// A crate README is a crates.io page and a place for a few design notes, not
/// a manual: the API belongs in rustdoc, the site in `docs/`.
const README_MAX_LINES: usize = 150;

/// The workspace crate count the README scan must at least reach, so a
/// layout change that makes the filter match nothing fails instead of passing.
const MIN_CRATES: usize = 50;

/// Files exempt from the marker and pointer scans. Dated records may name
/// files that no longer exist, and rewriting them would falsify history.
fn is_exempt(path: &str) -> bool {
    path.starts_with("docs/adr/")
        || path.ends_with("CHANGELOG.md")
        || path == "docs/reference/changelog.md"
        || path.starts_with("docs/research/")
        // This checker names the files it forbids.
        || path == "tools/xtask/src/context.rs"
}

pub fn check() -> Result<()> {
    let metadata = MetadataCommand::new()
        .no_deps()
        .exec()
        .map_err(|error| format!("cargo metadata failed: {error}"))?;
    let root = PathBuf::from(metadata.workspace_root.as_str());
    let files = repository_files(&root)?;
    let mut problems = Vec::new();

    root_agents(&root, &mut problems);
    retired_context_files(&files, &mut problems);
    root_files(&files, &mut problems);

    let mut crates = 0;
    for package in &metadata.packages {
        let manifest = PathBuf::from(package.manifest_path.as_str());
        let dir = manifest.parent().expect("a manifest has a directory");
        let rel = relative(&root, dir);
        if !(rel.starts_with("crates/") || rel.starts_with("tools/")) {
            continue;
        }
        crates += 1;
        crate_readme(&root, dir, &mut problems);
        let publishable = package.publish.as_ref().is_none_or(|r| !r.is_empty());
        // cargo reports the path as written, relative to the manifest.
        let declared = package.readme.as_ref().map(|p| p.as_str());
        if publishable && declared != Some("README.md") {
            problems.push(format!(
                "{rel}/Cargo.toml: declare `readme = \"README.md\"` so crates.io shows the crate's own page"
            ));
        }
    }
    if crates < MIN_CRATES {
        problems.push(format!(
            "found {crates} crates under crates/ and tools/, expected at least {MIN_CRATES}: the scan no longer matches the layout"
        ));
    }

    for file in &files {
        if file.ends_with(".rs") && !is_exempt(file) {
            todo_markers(&root, file, &mut problems);
        }
        if (file.ends_with(".rs") || file.ends_with(".md")) && !is_exempt(file) {
            pointers(&root, file, &mut problems);
        }
    }

    if problems.is_empty() {
        println!(
            "context: ok ({crates} crate READMEs, {} files scanned)",
            files.len()
        );
        Ok(())
    } else {
        Err(format!(
            "{} problem(s):\n  {}",
            problems.len(),
            problems.join("\n  ")
        ))
    }
}

/// Tracked files plus untracked ones not ignored, so a new plan file is caught
/// before it is committed.
fn repository_files(root: &Path) -> Result<Vec<String>> {
    let output = Command::new("git")
        .args(["ls-files", "--cached", "--others", "--exclude-standard"])
        .current_dir(root)
        .output()
        .map_err(|error| format!("git ls-files failed: {error}"))?;
    if !output.status.success() {
        return Err("git ls-files failed".into());
    }
    let mut files: Vec<String> = String::from_utf8_lossy(&output.stdout)
        .lines()
        .filter(|path| !path.starts_with("docs/node_modules/") && root.join(path).is_file())
        .map(str::to_owned)
        .collect();
    files.sort();
    files.dedup();
    Ok(files)
}

fn relative(root: &Path, path: &Path) -> String {
    path.strip_prefix(root)
        .unwrap_or(path)
        .to_string_lossy()
        .replace('\\', "/")
}

fn has_checkbox(text: &str) -> Option<usize> {
    text.lines().position(|line| {
        let line = line.trim_start();
        ["- [ ]", "- [x]", "- [X]", "* [ ]", "* [x]", "* [X]"]
            .iter()
            .any(|box_| line.starts_with(box_))
    })
}

fn root_agents(root: &Path, problems: &mut Vec<String>) {
    let Ok(text) = std::fs::read_to_string(root.join("AGENTS.md")) else {
        problems.push("AGENTS.md: the root AGENTS.md is missing".into());
        return;
    };
    let lines = text.lines().count();
    if lines > ROOT_AGENTS_MAX_LINES {
        problems.push(format!(
            "AGENTS.md: {lines} lines, limit {ROOT_AGENTS_MAX_LINES}; move detail into a README, module doc or ADR"
        ));
    }
    if let Some(line) = has_checkbox(&text) {
        problems.push(format!(
            "AGENTS.md:{}: task checkbox; open work lives in issues",
            line + 1
        ));
    }
}

/// Nested agent context and checked-in plans are retired: their content lives
/// in module docs, READMEs, ADRs, tests and issues.
fn retired_context_files(files: &[String], problems: &mut Vec<String>) {
    for file in files {
        let name = file.rsplit('/').next().unwrap_or(file).to_ascii_lowercase();
        let nested_agents = (name == "agents.md" || name == "claude.md") && file != "AGENTS.md";
        let plan = name.ends_with(".md")
            && (name == "plan.md"
                || name.starts_with("plan-")
                || name.starts_with("plan_")
                || name.ends_with("-plan.md"));
        if nested_agents {
            problems.push(format!(
                "{file}: only the root AGENTS.md exists; put crate context in its README.md or module docs"
            ));
        } else if plan || file.starts_with("docs/plans/") {
            problems.push(format!(
                "{file}: plans are not checked in; open work lives in GitHub issues (ADR 0078)"
            ));
        }
    }
}

/// The repository root holds the workspace manifest and its front pages.
/// Session plans, findings and notes do not belong there: their content goes
/// to issues, ADRs, tests or module docs.
const ROOT_FILES: &[&str] = &[
    "AGENTS.md",
    "Cargo.lock",
    "Cargo.toml",
    "LICENSE",
    "README.md",
    "rust-toolchain.toml",
];

fn root_files(files: &[String], problems: &mut Vec<String>) {
    for file in files {
        if !file.contains('/') && !file.starts_with('.') && !ROOT_FILES.contains(&file.as_str()) {
            problems.push(format!(
                "{file}: the root holds only {}; move it into an issue, ADR, test or module doc",
                ROOT_FILES.join(", ")
            ));
        }
    }
}

fn crate_readme(root: &Path, dir: &Path, problems: &mut Vec<String>) {
    let readme = dir.join("README.md");
    let rel = relative(root, &readme);
    let Ok(text) = std::fs::read_to_string(&readme) else {
        problems.push(format!("{rel}: every crate has a README.md"));
        return;
    };
    let lines = text.lines().count();
    if lines > README_MAX_LINES {
        problems.push(format!(
            "{rel}: {lines} lines, limit {README_MAX_LINES}; the API belongs in rustdoc, the guide in docs/"
        ));
    }
    if let Some(line) = has_checkbox(&text) {
        problems.push(format!(
            "{rel}:{}: task checkbox; open work lives in issues",
            line + 1
        ));
    }
}

/// A work marker in code names its issue: `TODO(#123)`. Anything else is a
/// plan nobody tracks.
fn todo_markers(root: &Path, file: &str, problems: &mut Vec<String>) {
    let Ok(text) = std::fs::read_to_string(root.join(file)) else {
        return;
    };
    for (index, line) in text.lines().enumerate() {
        for marker in ["TODO", "FIXME", "XXX"] {
            let mut rest = line;
            while let Some(at) = rest.find(marker) {
                let before = rest[..at].chars().next_back();
                let after = &rest[at + marker.len()..];
                let is_word = before.is_none_or(|c| !c.is_alphanumeric() && c != '_')
                    && after
                        .chars()
                        .next()
                        .is_none_or(|c| !c.is_alphanumeric() && c != '_');
                if is_word && !(marker == "TODO" && names_issue(after)) {
                    problems.push(format!(
                        "{file}:{}: `{marker}` without an issue; write `TODO(#N)`",
                        index + 1
                    ));
                }
                rest = after;
            }
        }
    }
}

fn names_issue(after: &str) -> bool {
    let Some(inner) = after.strip_prefix("(#") else {
        return false;
    };
    let digits = inner.chars().take_while(char::is_ascii_digit).count();
    digits > 0 && inner[digits..].starts_with(')')
}

/// Path-like tokens ending in a context file name, e.g. `../README.md` or
/// `crates/foo/AGENTS.md`.
fn path_tokens<'a>(text: &'a str, name: &str) -> Vec<(usize, &'a str)> {
    let mut tokens = Vec::new();
    for (index, line) in text.lines().enumerate() {
        let mut from = 0;
        while let Some(at) = line[from..].find(name) {
            let end = from + at + name.len();
            let start = line[..from + at]
                .rfind(|c: char| !(c.is_alphanumeric() || "._-/".contains(c)))
                .map_or(0, |i| i + 1);
            tokens.push((index + 1, &line[start..end]));
            from = end;
        }
    }
    tokens
}

fn pointers(root: &Path, file: &str, problems: &mut Vec<String>) {
    let Ok(text) = std::fs::read_to_string(root.join(file)) else {
        return;
    };
    for (line, token) in path_tokens(&text, "PLAN.md") {
        problems.push(format!(
            "{file}:{line}: `{token}` points at a retired plan; link the issue instead"
        ));
    }
    for (line, token) in path_tokens(&text, "AGENTS.md") {
        if token.contains('/') {
            problems.push(format!(
                "{file}:{line}: `{token}` points at a nested AGENTS.md; only the root one exists"
            ));
        }
    }
    let dir = Path::new(file).parent().unwrap_or(Path::new(""));
    for (line, token) in path_tokens(&text, "README.md") {
        if !token.contains('/') || token.starts_with("http") || token.contains("://") {
            continue;
        }
        let resolves = root.join(dir).join(token).is_file() || root.join(token).is_file();
        if !resolves {
            problems.push(format!("{file}:{line}: `{token}` does not resolve"));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn issue_markers_need_a_number() {
        assert!(names_issue("(#12) fix"));
        assert!(!names_issue("(#) fix"));
        assert!(!names_issue(": fix"));
        assert!(!names_issue("(12)"));
    }

    #[test]
    fn tokens_stop_at_punctuation() {
        let tokens = path_tokens("see `../foo/README.md`, and [x](a/README.md)", "README.md");
        assert_eq!(tokens, vec![(1, "../foo/README.md"), (1, "a/README.md")]);
    }

    #[test]
    fn checkboxes_are_found() {
        assert_eq!(has_checkbox("a\n  - [ ] b"), Some(1));
        assert_eq!(has_checkbox("- [link](x)"), None);
    }
}
