use super::model::{Architecture, Result};
use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

const LAYERS: &[(&str, &str, &str)] = &[
    ("foundation", "foundation.", "crates/foundation/"),
    (
        "representations",
        "representation.",
        "crates/representations/",
    ),
    ("contracts", "contract.", "crates/contracts/"),
    ("algorithms", "algorithm.", "crates/algorithms/"),
    ("providers", "provider.", "crates/providers/"),
    ("execution", "execution.", "crates/execution/"),
    ("facade", "facade", "crates/facade/"),
    ("tools", "tool.", "tools/"),
];

pub fn validate(architecture: &Architecture) -> Result<()> {
    let mut errors = Vec::new();
    let workspace_names: BTreeSet<_> = architecture.packages.keys().cloned().collect();
    check_registered_manifests(architecture, &mut errors)?;

    for package in architecture.packages.values() {
        let Some((_, role_prefix, path_prefix)) =
            LAYERS.iter().find(|(layer, _, _)| *layer == package.layer)
        else {
            errors.push(format!(
                "{}: unknown layer `{}`",
                package.name, package.layer
            ));
            continue;
        };
        if !package.role.starts_with(role_prefix) {
            errors.push(format!(
                "{}: role `{}` does not agree with layer `{}`",
                package.name, package.role, package.layer
            ));
        }
        if !package.path.starts_with(path_prefix) {
            errors.push(format!(
                "{}: path `{}` must be under `{}` for layer `{}`",
                package.name, package.path, path_prefix, package.layer
            ));
        }
        if package.domain.ends_with(".temporary") || package.domain == "common.temporary" {
            // Temporary migration ownership is explicit and visible in generated docs.
        }
        for declared in &package.allowed_internal_dependencies {
            if !workspace_names.contains(declared) {
                errors.push(format!(
                    "{}: allowed dependency `{declared}` is not a workspace package",
                    package.name
                ));
            }
        }
        if package.actual_internal_dependencies != package.allowed_internal_dependencies {
            let undeclared: Vec<_> = package
                .actual_internal_dependencies
                .difference(&package.allowed_internal_dependencies)
                .cloned()
                .collect();
            let stale: Vec<_> = package
                .allowed_internal_dependencies
                .difference(&package.actual_internal_dependencies)
                .cloned()
                .collect();
            if !undeclared.is_empty() {
                errors.push(format!(
                    "{}: undeclared internal dependencies: {}",
                    package.name,
                    undeclared.join(", ")
                ));
            }
            if !stale.is_empty() {
                errors.push(format!(
                    "{}: stale allowed internal dependencies: {}",
                    package.name,
                    stale.join(", ")
                ));
            }
        }
        for dependency in &package.production_internal_dependencies {
            let dependency_package = &architecture.packages[dependency];
            if !layer_edge_allowed(&package.layer, &dependency_package.layer) {
                errors.push(format!(
                    "{} ({}) must not depend on {} ({})",
                    package.name, package.layer, dependency, dependency_package.layer
                ));
            }
        }
        if package.public && package.layer == "tools" {
            errors.push(format!(
                "{}: tools must not be public packages",
                package.name
            ));
        }
        if !package.format_neutral && package.layer != "tools" {
            errors.push(format!(
                "{}: Axiolid production packages must remain format-neutral",
                package.name
            ));
        }
    }

    super::naming::validate(architecture, &mut errors)?;
    super::source_checks::validate(architecture, &mut errors)?;
    super::source_checks::validate_no_product_names(&architecture.root, &mut errors)?;
    if errors.is_empty() {
        Ok(())
    } else {
        errors.sort();
        Err(format!(
            "{} architecture violation(s):\n- {}",
            errors.len(),
            errors.join("\n- ")
        ))
    }
}

fn check_registered_manifests(architecture: &Architecture, errors: &mut Vec<String>) -> Result<()> {
    let registered: BTreeSet<_> = architecture
        .packages
        .values()
        .map(|package| format!("{}/Cargo.toml", package.path))
        .collect();
    let mut discovered = Vec::new();
    for root in ["crates", "tools"] {
        collect_manifests(&architecture.root.join(root), &mut discovered)?;
    }
    for manifest in discovered {
        let relative = manifest
            .strip_prefix(&architecture.root)
            .map_err(|_| format!("{} is outside workspace root", manifest.display()))?
            .to_string_lossy()
            .replace('\\', "/");
        if !registered.contains(&relative) {
            errors.push(format!(
                "unregistered package manifest `{relative}`; add an explicit workspace member and architecture metadata or remove it"
            ));
        }
    }
    Ok(())
}

fn collect_manifests(directory: &Path, manifests: &mut Vec<std::path::PathBuf>) -> Result<()> {
    if !directory.exists() {
        return Ok(());
    }
    for entry in
        fs::read_dir(directory).map_err(|error| format!("read {}: {error}", directory.display()))?
    {
        let path = entry
            .map_err(|error| format!("read directory entry: {error}"))?
            .path();
        if path.is_dir() {
            collect_manifests(&path, manifests)?;
        } else if path.file_name().and_then(|name| name.to_str()) == Some("Cargo.toml") {
            manifests.push(path);
        }
    }
    Ok(())
}

fn layer_edge_allowed(from: &str, to: &str) -> bool {
    match from {
        "foundation" => false,
        "representations" => matches!(to, "foundation" | "representations"),
        "contracts" => matches!(to, "foundation" | "representations" | "contracts"),
        "algorithms" => matches!(
            to,
            "foundation" | "representations" | "contracts" | "algorithms"
        ),
        "providers" => matches!(
            to,
            "foundation" | "representations" | "contracts" | "algorithms" | "providers"
        ),
        "execution" => !matches!(to, "facade" | "tools"),
        "facade" => to != "tools",
        "tools" => true,
        _ => false,
    }
}

#[cfg(test)]
mod invariant_tests {
    use super::layer_edge_allowed;

    #[test]
    fn role_dag_rejects_upward_dependencies() {
        assert!(!layer_edge_allowed("foundation", "representations"));
        assert!(!layer_edge_allowed("representations", "algorithms"));
        assert!(!layer_edge_allowed("contracts", "providers"));
        assert!(!layer_edge_allowed("algorithms", "execution"));
        assert!(!layer_edge_allowed("providers", "execution"));
        assert!(!layer_edge_allowed("execution", "facade"));
    }

    #[test]
    fn role_dag_allows_declared_cross_role_shapes() {
        assert!(layer_edge_allowed("representations", "foundation"));
        assert!(layer_edge_allowed("contracts", "representations"));
        assert!(layer_edge_allowed("algorithms", "contracts"));
        assert!(layer_edge_allowed("providers", "algorithms"));
        assert!(layer_edge_allowed("execution", "providers"));
        assert!(layer_edge_allowed("facade", "execution"));
    }
}

#[cfg(test)]
mod boundary_tests {
    use super::layer_edge_allowed;

    #[test]
    fn foundation_is_a_true_root() {
        assert!(!layer_edge_allowed("foundation", "representations"));
        assert!(!layer_edge_allowed("foundation", "foundation"));
    }

    #[test]
    fn contracts_cannot_reach_implementations() {
        assert!(!layer_edge_allowed("contracts", "algorithms"));
        assert!(!layer_edge_allowed("contracts", "providers"));
        assert!(!layer_edge_allowed("contracts", "execution"));
    }

    #[test]
    fn algorithms_cannot_select_providers_or_execution() {
        assert!(!layer_edge_allowed("algorithms", "providers"));
        assert!(!layer_edge_allowed("algorithms", "execution"));
        assert!(!layer_edge_allowed("algorithms", "facade"));
    }
}
