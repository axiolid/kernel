//! Minimal-dependency-closure verification.
//!
//! The ownership checker in `checks.rs` validates DECLARED metadata without
//! resolving dependencies. That cannot prove what a downstream application
//! actually compiles. This module resolves an isolated consumer fixture with
//! Cargo and compares the real internal package set against a declared profile.
//!
//! `cargo tree` is documented as close to, but not identical with, the build
//! graph, so a graph assertion alone is not proof: `check` also runs a real
//! `cargo check` for each fixture.

use super::model::Result;
use cargo_metadata::MetadataCommand;
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const PROFILES: &str = "docs/architecture/closure-profiles.toml";

#[derive(Debug, Clone)]
pub struct ClosureProfile {
    pub name: String,
    pub manifest: String,
    pub description: String,
    pub features: Vec<String>,
    pub expected_internal: BTreeSet<String>,
    pub forbidden_internal: BTreeSet<String>,
}

/// Resolved facts about one profile.
#[derive(Debug)]
pub struct ResolvedClosure {
    pub actual_internal: BTreeSet<String>,
    pub external: BTreeSet<String>,
}

pub fn load_profiles(root: &Path) -> Result<Vec<ClosureProfile>> {
    let path = root.join(PROFILES);
    let text =
        fs::read_to_string(&path).map_err(|error| format!("read {}: {error}", path.display()))?;
    parse_profiles(&text)
}

/// Minimal TOML array-of-tables reader for this file's fixed shape.
///
/// A dedicated parser dependency is not justified for one declaration file the
/// repository fully controls; the format is validated strictly instead.
fn parse_profiles(text: &str) -> Result<Vec<ClosureProfile>> {
    let mut profiles: Vec<ClosureProfile> = Vec::new();
    let mut current: Option<ClosureProfile> = None;
    let mut list_key: Option<String> = None;

    for raw in text.lines() {
        let line = strip_comment(raw).trim().to_string();
        if line.is_empty() {
            continue;
        }
        if line == "[[profile]]" {
            if let Some(profile) = current.take() {
                profiles.push(profile);
            }
            current = Some(ClosureProfile {
                name: String::new(),
                manifest: String::new(),
                description: String::new(),
                features: Vec::new(),
                expected_internal: BTreeSet::new(),
                forbidden_internal: BTreeSet::new(),
            });
            continue;
        }
        let profile = current
            .as_mut()
            .ok_or_else(|| format!("{PROFILES}: value outside a [[profile]] table: {line}"))?;

        if let Some(key) = list_key.clone() {
            if line == "]" {
                list_key = None;
                continue;
            }
            let value = unquote(line.trim_end_matches(','))?;
            match key.as_str() {
                "features" => profile.features.push(value),
                "expected_internal" => {
                    profile.expected_internal.insert(value);
                }
                "forbidden_internal" => {
                    profile.forbidden_internal.insert(value);
                }
                other => return Err(format!("{PROFILES}: unknown list `{other}`")),
            }
            continue;
        }

        let (key, value) = line
            .split_once('=')
            .ok_or_else(|| format!("{PROFILES}: expected `key = value`, got: {line}"))?;
        let key = key.trim();
        let value = value.trim();
        if value == "[" {
            list_key = Some(key.to_owned());
            continue;
        }
        if value.starts_with('[') && value.ends_with(']') {
            let inner = &value[1..value.len() - 1];
            for item in inner.split(',').map(str::trim).filter(|i| !i.is_empty()) {
                let item = unquote(item)?;
                match key {
                    "features" => profile.features.push(item),
                    "expected_internal" => {
                        profile.expected_internal.insert(item);
                    }
                    "forbidden_internal" => {
                        profile.forbidden_internal.insert(item);
                    }
                    other => return Err(format!("{PROFILES}: unknown list `{other}`")),
                }
            }
            continue;
        }
        match key {
            "name" => profile.name = unquote(value)?,
            "manifest" => profile.manifest = unquote(value)?,
            "description" => profile.description = unquote(value)?,
            other => return Err(format!("{PROFILES}: unknown key `{other}`")),
        }
    }
    if let Some(profile) = current.take() {
        profiles.push(profile);
    }
    if profiles.is_empty() {
        return Err(format!("{PROFILES}: no profiles declared"));
    }
    for profile in &profiles {
        if profile.name.is_empty() || profile.manifest.is_empty() {
            return Err(format!(
                "{PROFILES}: a profile is missing `name` or `manifest`"
            ));
        }
        if profile.expected_internal.is_empty() {
            return Err(format!(
                "{}: profile `{}` declares no expected packages",
                PROFILES, profile.name
            ));
        }
        // A package cannot be simultaneously required and banned; that would
        // make the profile unfalsifiable in one direction.
        let contradiction: Vec<_> = profile
            .expected_internal
            .intersection(&profile.forbidden_internal)
            .cloned()
            .collect();
        if !contradiction.is_empty() {
            return Err(format!(
                "{}: profile `{}` both expects and forbids: {}",
                PROFILES,
                profile.name,
                contradiction.join(", ")
            ));
        }
    }
    Ok(profiles)
}

fn strip_comment(line: &str) -> &str {
    let mut in_string = false;
    for (index, character) in line.char_indices() {
        match character {
            '"' => in_string = !in_string,
            '#' if !in_string => return &line[..index],
            _ => {}
        }
    }
    line
}

fn unquote(value: &str) -> Result<String> {
    let value = value.trim();
    value
        .strip_prefix('"')
        .and_then(|rest| rest.strip_suffix('"'))
        .map(ToOwned::to_owned)
        .ok_or_else(|| format!("{PROFILES}: expected a quoted string, got: {value}"))
}

/// Resolve a profile's real closure with Cargo.
pub fn resolve(root: &Path, profile: &ClosureProfile) -> Result<ResolvedClosure> {
    let manifest = root.join(&profile.manifest);
    if !manifest.exists() {
        return Err(format!(
            "profile `{}`: fixture manifest {} is missing",
            profile.name,
            manifest.display()
        ));
    }
    let mut command = MetadataCommand::new();
    command.manifest_path(&manifest);
    // `--offline`, not `--frozen`: the fixture is its own workspace root and
    // legitimately needs to write its own lock file on first resolution.
    command.other_options(vec!["--offline".to_owned()]);
    if profile.features.is_empty() {
        command.features(cargo_metadata::CargoOpt::NoDefaultFeatures);
    } else {
        command.features(cargo_metadata::CargoOpt::SomeFeatures(
            profile.features.clone(),
        ));
    }
    let metadata = command
        .exec()
        .map_err(|error| format!("profile `{}`: cargo metadata failed: {error}", profile.name))?;

    // The fixture's own root package is the consumer, not part of the Axiolid
    // closure. Resolve it from the metadata rather than matching a hardcoded
    // name, so a new profile cannot silently fail this exclusion.
    let root: Option<String> = metadata
        .resolve
        .as_ref()
        .and_then(|resolve| resolve.root.as_ref())
        .and_then(|id| {
            metadata
                .packages
                .iter()
                .find(|package| &package.id == id)
                .map(|package| package.name.to_string())
        });

    let mut internal = BTreeSet::new();
    let mut external = BTreeSet::new();
    for package in &metadata.packages {
        let name = package.name.to_string();
        if root.as_deref() == Some(name.as_str()) {
            continue;
        }
        if name == "axiolid" || name.starts_with("axiolid-") {
            internal.insert(name);
        } else {
            external.insert(name);
        }
    }
    Ok(ResolvedClosure {
        actual_internal: internal,
        external,
    })
}

const DOC: &str = "docs/architecture/closure-profiles.md";

/// Render the declared profiles and their measured closures as documentation.
///
/// Generated, not hand-maintained: a closure table that drifts from the
/// verified profiles is worse than none, because it invites trust.
fn render(root: &Path) -> Result<String> {
    let profiles = load_profiles(root)?;
    let mut out = String::new();
    out.push_str("# Minimal compilation closures\n\n");
    out.push_str("Generated by `cargo xtask architecture closure docs`. Do not edit by hand.\n\n");
    out.push_str(
        "Each profile is a compatibility promise verified by `cargo xtask architecture closure check`, \
which resolves an isolated downstream fixture and then compiles it. \
See [ADR 0036](../adr/0036-use-case-specific-compilation-closures.md).\n\n",
    );

    out.push_str("| Profile | Internal packages | Purpose |\n| --- | --- | --- |\n");
    for profile in &profiles {
        out.push_str(&format!(
            "| `{}` | {} | {} |\n",
            profile.name,
            profile.expected_internal.len(),
            profile.description
        ));
    }

    for profile in &profiles {
        out.push_str(&format!(
            "\n## {}\n\n{}\n\n",
            profile.name, profile.description
        ));
        out.push_str(&format!("Fixture: `{}`\n\n", profile.manifest));
        out.push_str("Compiles exactly these internal packages:\n\n");
        for package in &profile.expected_internal {
            out.push_str(&format!("- `{package}`\n"));
        }
        out.push_str(&format!(
            "\nAnd must never acquire any of {} forbidden packages, including `{}`.\n",
            profile.forbidden_internal.len(),
            profile
                .forbidden_internal
                .iter()
                .take(3)
                .cloned()
                .collect::<Vec<_>>()
                .join("`, `")
        ));
    }
    Ok(out)
}

pub fn docs(root: &Path) -> Result<()> {
    let rendered = render(root)?;
    let path = root.join(DOC);
    fs::write(&path, rendered).map_err(|error| format!("write {DOC}: {error}"))?;
    println!("generated {DOC}");
    Ok(())
}

/// Fail when the committed documentation no longer matches the profiles.
pub fn check_docs(root: &Path) -> Result<()> {
    let rendered = render(root)?;
    let path = root.join(DOC);
    let current = fs::read_to_string(&path).unwrap_or_default();
    if current == rendered {
        Ok(())
    } else {
        Err(format!(
            "{DOC} is stale; run `cargo xtask architecture closure docs`"
        ))
    }
}

pub fn check(root: &Path) -> Result<()> {
    let profiles = load_profiles(root)?;
    let mut errors = Vec::new();

    for profile in &profiles {
        let resolved = resolve(root, profile)?;

        let unexpected: Vec<_> = resolved
            .actual_internal
            .difference(&profile.expected_internal)
            .cloned()
            .collect();
        let missing: Vec<_> = profile
            .expected_internal
            .difference(&resolved.actual_internal)
            .cloned()
            .collect();
        let forbidden: Vec<_> = resolved
            .actual_internal
            .intersection(&profile.forbidden_internal)
            .cloned()
            .collect();

        if !forbidden.is_empty() {
            errors.push(format!(
                "closure `{}` contains forbidden packages: {}",
                profile.name,
                forbidden.join(", ")
            ));
        }
        if !unexpected.is_empty() {
            errors.push(format!(
                "closure `{}` grew: {} (update the profile deliberately if intended)",
                profile.name,
                unexpected.join(", ")
            ));
        }
        if !missing.is_empty() {
            errors.push(format!(
                "closure `{}` no longer contains: {}",
                profile.name,
                missing.join(", ")
            ));
        }

        // A resolved graph is not a compiled program. Build the fixture too.
        let manifest = root.join(&profile.manifest);
        let mut cargo = Command::new(std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()));
        cargo
            .arg("check")
            .arg("--manifest-path")
            .arg(&manifest)
            .arg("--quiet");
        if profile.features.is_empty() {
            cargo.arg("--no-default-features");
        } else {
            cargo.arg("--features").arg(profile.features.join(","));
        }
        let status = cargo
            .status()
            .map_err(|error| format!("profile `{}`: cargo check failed: {error}", profile.name))?;
        if !status.success() {
            errors.push(format!("closure `{}` does not compile", profile.name));
        }
    }

    if let Err(error) = check_docs(root) {
        errors.push(error);
    }

    if errors.is_empty() {
        println!("closure check passed: {} profile(s)", profiles.len());
        Ok(())
    } else {
        errors.sort();
        Err(format!(
            "{} closure violation(s):\n- {}",
            errors.len(),
            errors.join("\n- ")
        ))
    }
}

pub fn explain(root: &Path, name: &str) -> Result<()> {
    let profiles = load_profiles(root)?;
    let profile = profiles
        .iter()
        .find(|profile| profile.name == name)
        .ok_or_else(|| {
            let available: Vec<_> = profiles.iter().map(|p| p.name.clone()).collect();
            format!(
                "unknown profile `{name}`; available: {}",
                available.join(", ")
            )
        })?;
    let resolved = resolve(root, profile)?;

    println!("profile: {}", profile.name);
    println!("{}", profile.description);
    println!();
    println!("fixture manifest: {}", profile.manifest);
    println!();
    println!("features: {:?}", profile.features);
    println!();
    println!(
        "resolved internal packages ({}):",
        resolved.actual_internal.len()
    );
    for package in &resolved.actual_internal {
        println!("    {package}");
    }
    println!();
    println!("external packages ({}):", resolved.external.len());
    for package in &resolved.external {
        println!("    {package}");
    }
    println!();
    println!(
        "forbidden packages checked: {}",
        profile.forbidden_internal.len()
    );
    println!();
    println!("verify with:");
    println!("    cargo xtask architecture closure check");
    println!(
        "    cargo tree --manifest-path {} --no-default-features --edges normal,build",
        profile.manifest
    );
    Ok(())
}

pub fn workspace_root() -> Result<PathBuf> {
    let metadata = MetadataCommand::new()
        .no_deps()
        .exec()
        .map_err(|error| format!("cargo metadata failed: {error}"))?;
    Ok(PathBuf::from(metadata.workspace_root.as_str()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_the_declared_profile_shape() {
        let profiles = parse_profiles(
            r#"
[[profile]]
name = "demo"
manifest = "tests/consumers/demo/Cargo.toml"
description = "d"
features = []
expected_internal = [
    "axiolid-core",
]
forbidden_internal = ["axiolid-mesh"]
"#,
        )
        .expect("valid profile");
        assert_eq!(profiles.len(), 1);
        assert!(profiles[0].expected_internal.contains("axiolid-core"));
        assert!(profiles[0].forbidden_internal.contains("axiolid-mesh"));
    }

    /// A profile that both requires and bans a package could never fail one of
    /// those checks honestly, so it is rejected at load time.
    #[test]
    fn rejects_a_contradictory_profile() {
        let error = parse_profiles(
            r#"
[[profile]]
name = "bad"
manifest = "m"
expected_internal = ["axiolid-core"]
forbidden_internal = ["axiolid-core"]
"#,
        )
        .expect_err("contradiction must be refused");
        assert!(error.contains("both expects and forbids"), "{error}");
    }

    #[test]
    fn rejects_a_profile_with_no_expectations() {
        let error = parse_profiles(
            r#"
[[profile]]
name = "empty"
manifest = "m"
expected_internal = []
"#,
        )
        .expect_err("an empty expectation proves nothing");
        assert!(error.contains("no expected packages"), "{error}");
    }

    #[test]
    fn comment_stripping_preserves_hashes_inside_strings() {
        assert_eq!(
            strip_comment(r#"name = "a#b" # trailing"#).trim(),
            r#"name = "a#b""#
        );
    }
}
