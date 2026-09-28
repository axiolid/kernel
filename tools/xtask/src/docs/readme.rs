//! A crate's README: the prose its reference page shows, and the lint that
//! keeps it a crates.io page.
//!
//! The README is the one hand-written summary of a crate (ADR 0078), so the
//! page takes its Overview and its `##` sections (the design notes) from
//! there rather than from the `//!` docs, which would say the same twice.
//!
//! Each publishable crate's README is its crates.io page, so it links to the
//! three places the crate is documented and says nothing that turns false at
//! the next release.
//!
//! - The links: `docs.rs/<crate>`, the generated reference page
//!   `<site>/reference/crates/<crate>`, and the repository. Each names the
//!   README's own crate, so a README copied from a sibling fails.
//! - No publication claim ("not yet published", "unpublished", …): crates.io
//!   copies the README at publish time, where such a line is false on
//!   arrival.
//! - No version claim: a three-part version (`0.3.1`), a `v`-prefixed one
//!   (`v0.4`) or a `version = "…"` pin. crates.io shows the version; a
//!   README that repeats it drifts at the next bump.

use super::changelog::absolutise;
use super::workspace::{Crate, Workspace};
use super::{Result, REPO, SITE};

/// The README's path relative to the workspace root.
fn readme_path(krate: &Crate) -> String {
    format!("{}/{}", krate.dir, "README.md")
}

/// The parts of a README a reference page shows.
pub(super) struct Prose {
    /// The introduction: everything between the title and the install
    /// snippet or link list.
    pub(super) overview: String,
    /// Every `##` section, verbatim.
    pub(super) notes: String,
}

pub(super) fn prose(workspace: &Workspace, krate: &Crate) -> Result<Prose> {
    let rel = readme_path(krate);
    let text = std::fs::read_to_string(workspace.root.join(&rel))
        .map_err(|error| format!("{rel}: {error}"))?;
    let (overview, notes) = split(&text);
    let publish = |part: &str| escape(&absolutise(part.trim(), &krate.dir));
    Ok(Prose {
        overview: publish(&overview),
        notes: publish(&notes),
    })
}

/// (introduction, `##` sections) of a README.
fn split(text: &str) -> (String, String) {
    let mut overview = Vec::new();
    let mut notes = Vec::new();
    let mut in_intro = true;
    let mut in_notes = false;
    for line in text.lines().skip_while(|l| !l.starts_with("# ")).skip(1) {
        if line.starts_with("## ") {
            in_intro = false;
            in_notes = true;
        } else if in_intro && (line.starts_with("```") || line.starts_with("- ")) {
            in_intro = false;
        }
        if in_intro {
            overview.push(line);
        } else if in_notes {
            notes.push(line);
        }
    }
    (overview.join("\n"), notes.join("\n"))
}

/// Escape `{{` outside code: VitePress compiles each page as a Vue template,
/// where it would open an interpolation.
fn escape(text: &str) -> String {
    let mut out = Vec::new();
    let mut fenced = false;
    for line in text.lines() {
        if line.trim_start().starts_with("```") {
            fenced = !fenced;
        }
        if fenced || !line.contains("{{") {
            out.push(line.to_owned());
            continue;
        }
        let mut escaped = String::with_capacity(line.len());
        let mut in_code = false;
        let mut chars = line.chars().peekable();
        while let Some(c) = chars.next() {
            match c {
                '`' => {
                    in_code = !in_code;
                    escaped.push(c);
                }
                '{' if !in_code && chars.peek() == Some(&'{') => escaped.push_str("&#123;"),
                _ => escaped.push(c),
            }
        }
        out.push(escaped);
    }
    out.join("\n")
}

/// `README:line: why` for every rule a publishable README breaks.
pub(super) fn problems(workspace: &Workspace) -> Vec<String> {
    let mut out = Vec::new();
    for krate in &workspace.crates {
        let rel = readme_path(krate);
        match std::fs::read_to_string(workspace.root.join(&rel)) {
            Ok(text) => out.extend(check(&krate.name, &rel, &text)),
            Err(error) => out.push(format!("{rel}: {error}")),
        }
    }
    out
}

/// The README lines every publishable crate carries, in this order.
pub(super) fn links(name: &str) -> [String; 3] {
    [
        format!("- API documentation: [docs.rs/{name}](https://docs.rs/{name})"),
        format!(
            "- Reference page: [{}]({SITE}/reference/crates/{name})",
            SITE.trim_start_matches("https://")
        ),
        format!(
            "- Source and issues: [{}]({REPO})",
            REPO.trim_start_matches("https://")
        ),
    ]
}

fn check(name: &str, rel: &str, text: &str) -> Vec<String> {
    let mut out = Vec::new();
    for line in links(name) {
        if !text.lines().any(|l| l.trim_end() == line) {
            out.push(format!("{rel}: missing the line `{line}`"));
        }
    }
    for (index, line) in text.lines().enumerate() {
        let at = index + 1;
        if claims_unpublished(line) {
            out.push(format!(
                "{rel}:{at}: publication claim; crates.io copies the README when it publishes, so drop it"
            ));
        }
        if claims_version(line) {
            out.push(format!(
                "{rel}:{at}: version claim; crates.io shows the version, so the README does not repeat it"
            ));
        }
    }
    out
}

/// Phrases that assert a package is not (yet) on a registry.
const UNPUBLISHED: &[&str] = &[
    "not published",
    "not yet published",
    "unpublished",
    "not been published",
    "not released yet",
    "not yet released",
    "not on crates.io",
    "not yet on crates.io",
];

fn claims_unpublished(line: &str) -> bool {
    let lower = line.to_ascii_lowercase();
    UNPUBLISHED.iter().any(|phrase| lower.contains(phrase))
}

/// A three-part version, a `v`-prefixed version, or a `version = "` pin.
fn claims_version(line: &str) -> bool {
    if line.contains("version = \"") {
        return true;
    }
    let bytes = line.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let starts_token =
            i == 0 || !(bytes[i - 1].is_ascii_alphanumeric() || bytes[i - 1] == b'.');
        if !(starts_token && bytes[i].is_ascii_alphanumeric()) {
            i += 1;
            continue;
        }
        let prefixed = bytes[i] == b'v';
        let start = if prefixed { i + 1 } else { i };
        let mut parts = 0;
        let mut j = start;
        loop {
            let digits = bytes[j..].iter().take_while(|b| b.is_ascii_digit()).count();
            if digits == 0 {
                break;
            }
            parts += 1;
            j += digits;
            if j + 1 < bytes.len() && bytes[j] == b'.' && bytes[j + 1].is_ascii_digit() {
                j += 1;
            } else {
                break;
            }
        }
        if parts >= 3 || (prefixed && parts >= 2) {
            return true;
        }
        i = j.max(i + 1);
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    fn clean(name: &str) -> String {
        format!("# {name}\n\nWhat it is.\n\n{}\n", links(name).join("\n"))
    }

    #[test]
    fn a_readme_with_its_own_links_passes() {
        assert!(check("axiolid-x", "R", &clean("axiolid-x")).is_empty());
    }

    #[test]
    fn links_must_name_the_readmes_own_crate() {
        let problems = check("axiolid-x", "R", &clean("axiolid-y"));
        assert_eq!(problems.len(), 2, "{problems:?}");
        assert!(problems.iter().any(|p| p.contains("docs.rs/axiolid-x")));
        assert!(problems
            .iter()
            .any(|p| p.contains("reference/crates/axiolid-x")));
    }

    #[test]
    fn publication_and_version_claims_are_rejected() {
        let text = format!(
            "{}\nIt is not yet published.\nSee v0.4 and 0.3.1.\naxiolid = {{ version = \"0.3\" }}\n",
            clean("axiolid-x")
        );
        let problems = check("axiolid-x", "R", &text);
        let lines: Vec<&str> = problems
            .iter()
            .filter_map(|p| p.split(':').nth(1))
            .collect();
        assert_eq!(lines, ["9", "10", "11"], "{problems:?}");
    }

    #[test]
    fn prose_is_the_intro_and_the_sections() {
        let text = "# axiolid-x\n\nWhat it is,\nand is not.\n\n```bash\ncargo add axiolid-x\n```\n\n- API documentation: x\n\n## Design notes\n\nWhy.\n\n```rust\nlet a = {{ b }};\n```\n";
        let (overview, notes) = split(text);
        assert_eq!(overview.trim(), "What it is,\nand is not.");
        assert_eq!(
            notes,
            "## Design notes\n\nWhy.\n\n```rust\nlet a = {{ b }};\n```"
        );
        assert_eq!(escape(&notes), notes);
        assert_eq!(escape("a {{ b }} `{{ c }}`"), "a &#123;{ b }} `{{ c }}`");
    }

    #[test]
    fn plain_numbers_are_not_versions() {
        assert!(!claims_version("ADR 0078, a 0.5 mm gap, 1e-9, kernel#9"));
        assert!(!claims_version("orient2d and Point3 values"));
        assert!(claims_version("since 1.2.3"));
        assert!(claims_version("(v0.4 ABI)"));
    }
}
