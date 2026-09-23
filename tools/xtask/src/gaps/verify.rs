use super::model::{self, Ledger, Level, Result};
use std::collections::BTreeSet;
use std::path::Path;

const AREAS: &[&str] = &[
    "foundations",
    "curves-surfaces",
    "brep",
    "mesh",
    "triangulation",
    "planar",
    "point-sets",
    "spatial",
];

pub fn check() -> Result<()> {
    let root = model::workspace_root()?;
    let ledger = model::load(&root)?;
    let packages = model::load_packages(&root)?;
    let errors = validate(&ledger, &packages, &root);
    if !errors.is_empty() {
        return Err(format!("{}:\n  {}", model::LEDGER, errors.join("\n  ")));
    }
    let open = ledger
        .capability
        .iter()
        .filter(|row| row.level.is_open())
        .count();
    let scoped = ledger
        .capability
        .iter()
        .filter(|row| row.level == Level::Scoped)
        .count();
    println!(
        "capability ledger ok: {} capabilities ({open} open, {scoped} scoped), {} issues",
        ledger.capability.len(),
        ledger.issue.len()
    );
    Ok(())
}

/// Every mechanically checkable property; returns all violations, not the first.
pub fn validate(ledger: &Ledger, packages: &model::ReferencePackages, root: &Path) -> Vec<String> {
    let mut errors = Vec::new();
    for row in &ledger.capability {
        for (library, paths, known) in [
            ("occt", &row.occt, &packages.occt),
            ("cgal", &row.cgal, &packages.cgal),
        ] {
            for path in paths.iter().filter(|path| !known.contains(*path)) {
                errors.push(format!(
                    "{}: {library} path `{path}` is not a package in the pinned tree ({})",
                    row.id,
                    model::PACKAGES
                ));
            }
        }
    }
    for (name, reference) in &ledger.reference {
        if reference.commit.len() != 40 || !reference.commit.bytes().all(|b| b.is_ascii_hexdigit())
        {
            errors.push(format!(
                "reference.{name}: commit must be a full 40-hex SHA"
            ));
        }
    }
    if !ledger.reference.contains_key("occt") || !ledger.reference.contains_key("cgal") {
        errors.push("reference.occt and reference.cgal must both be pinned".to_owned());
    }
    let mut issue_keys = BTreeSet::new();
    for issue in &ledger.issue {
        if !issue_keys.insert(issue.key.as_str()) {
            errors.push(format!("issue `{}` declared twice", issue.key));
        }
        for (field, value) in [("priority", &issue.priority), ("effort", &issue.effort)] {
            if !["Urgent", "High", "Medium", "Low"].contains(&value.as_str()) {
                errors.push(format!(
                    "issue `{}`: {field} `{value}` is not a board option",
                    issue.key
                ));
            }
        }
    }
    for issue in &ledger.issue {
        for blocker in &issue.blocked_by {
            if !issue_keys.contains(blocker.as_str()) || blocker == &issue.key {
                errors.push(format!(
                    "issue `{}`: blocked_by `{blocker}` is not another issue key",
                    issue.key
                ));
            }
        }
    }
    let mut ids = BTreeSet::new();
    let mut used_issues = BTreeSet::new();
    for row in &ledger.capability {
        let id = &row.id;
        if !ids.insert(id.as_str()) {
            errors.push(format!("{id}: declared twice"));
        }
        if !AREAS.contains(&row.area.as_str()) {
            errors.push(format!("{id}: unknown area `{}`", row.area));
        }
        if row.summary.trim().is_empty() {
            errors.push(format!("{id}: empty summary"));
        }
        // A claim of capability must point at code. `absent` needs none;
        // `scoped` may cite what exists but need not (out-of-scope rows).
        if matches!(row.level, Level::Narrow | Level::Implemented) && row.evidence.is_empty() {
            errors.push(format!(
                "{id}: level `{}` cites no evidence",
                row.level.as_str()
            ));
        }
        for cite in &row.evidence {
            if let Some(problem) = evidence_problem(root, cite) {
                errors.push(format!("{id}: evidence `{cite}`: {problem}"));
            }
        }
        if row.occt.is_empty() && row.cgal.is_empty() {
            errors.push(format!(
                "{id}: no OCCT or CGAL reference; not a comparison row"
            ));
        }
        match &row.issue {
            Some(key) if !issue_keys.contains(key.as_str()) => {
                errors.push(format!("{id}: issue `{key}` is not declared in [[issue]]"));
            }
            Some(key) => {
                used_issues.insert(key.as_str());
                if !row.level.is_open() {
                    errors.push(format!(
                        "{id}: {} but still assigned to issue `{key}`; remove `issue`",
                        row.level.as_str()
                    ));
                }
            }
            None => {}
        }
        let rationale = row
            .scope_rationale
            .as_deref()
            .is_some_and(|text| !text.trim().is_empty());
        match (row.level == Level::Scoped, rationale) {
            (true, false) => errors.push(format!(
                "{id}: scoped needs a non-empty scope_rationale saying why it is not built"
            )),
            (false, true) => errors.push(format!(
                "{id}: scope_rationale on a `{}` row; only scoped rows carry one",
                row.level.as_str()
            )),
            _ => {}
        }
    }
    for key in issue_keys.difference(&used_issues) {
        errors.push(format!("issue `{key}` tracks no capability row"));
    }
    errors
}

/// `path` or `path::symbol`. The file must exist; the symbol, when given,
/// must be defined in it (a rename then fails the gate instead of rotting).
fn evidence_problem(root: &Path, cite: &str) -> Option<String> {
    let (path, symbol) = match cite.split_once("::") {
        Some((path, symbol)) => (path, Some(symbol)),
        None => (cite, None),
    };
    let text = match std::fs::read_to_string(root.join(path)) {
        Ok(text) => text,
        Err(_) => return Some("file does not exist".to_owned()),
    };
    let symbol = symbol?;
    let name = symbol.rsplit("::").next().unwrap_or(symbol);
    defines(&text, name).then_some(()).map_or(
        Some(format!("`{name}` is not defined in this file")),
        |()| None,
    )
}

/// Whether an item named `name` is declared: `fn`, `struct`, `enum`,
/// `trait`, `type`, `const`, `static` or `mod`, at any visibility.
fn defines(text: &str, name: &str) -> bool {
    const KEYWORDS: &[&str] = &[
        "fn", "struct", "enum", "trait", "type", "const", "static", "mod",
    ];
    text.lines().any(|line| {
        let mut words = line.split(|c: char| !(c.is_alphanumeric() || c == '_'));
        let mut previous = "";
        words.any(|word| {
            let hit = word == name && KEYWORDS.contains(&previous);
            if !word.is_empty() {
                previous = word;
            }
            hit
        })
    })
}
