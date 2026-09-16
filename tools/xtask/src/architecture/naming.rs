//! Crate names derived from architecture metadata.
//!
//! ADR 0004 chose machine-checkable metadata over a naming convention, and
//! naming then drifted precisely because nothing read the names. This module
//! closes that gap: the expected name is DERIVED from `role` + `domain`, so a
//! name and its metadata cannot disagree silently.
//!
//! See ADR 0064.

use super::model::{Architecture, PackageArchitecture, Result};

/// Names that predate the rule and cannot change without a breaking release.
///
/// Each entry is `(actual, expected, reason)`. Deliberately NOT tied to a
/// version: a rename is not worth manufacturing a breaking release for, so
/// these land in the next release that is breaking for an independent
/// reason.
///
/// The list is CLOSED: a
/// violation not listed here fails, and an entry that no longer violates
/// ALSO fails, so a fixed exception cannot linger as permanent amnesty.
const GRANDFATHERED: &[(&str, &str, &str)] = &[
    (
        "axiolid-tessellation-contract",
        "axiolid-tessellate-contract",
        "published; noun for a verb domain; rename at the next breaking release",
    ),
    (
        "axiolid-exact-compile-contract",
        "axiolid-brep-compile-contract",
        "published; domain is brep.compile; rename at the next breaking release",
    ),
    (
        "axiolid-mesh-compile",
        "axiolid-graph-compile",
        "published; orchestrator impersonating a contract's provider; rename at the next breaking release",
    ),
];

/// `mesh.boolean` -> `mesh-boolean`.
fn domain_slug(domain: &str) -> String {
    domain.replace('.', "-")
}

/// The name a package's own metadata implies, or `None` where the rule
/// constrains the name without fully determining it.
///
/// A provider's engine suffix is a free choice -- `boolmesh` names an
/// upstream library, not a domain -- so providers are checked by PREFIX
/// rather than equality.
fn expected_name(package: &PackageArchitecture) -> Option<String> {
    if package.role == "contract.operation" {
        return Some(format!("axiolid-{}-contract", domain_slug(&package.domain)));
    }
    None
}

/// Check every package name against the rule derived from its metadata.
pub fn validate(architecture: &Architecture, errors: &mut Vec<String>) -> Result<()> {
    // Reserved: a contract's name with `-contract` stripped. Nobody may take
    // it, because it reads as THE implementation of that capability.
    let reserved: Vec<(String, String)> = architecture
        .packages
        .values()
        .filter(|p| p.role == "contract.operation")
        .map(|p| {
            (
                format!("axiolid-{}", domain_slug(&p.domain)),
                p.name.clone(),
            )
        })
        .collect();

    let mut unused_exceptions: Vec<&str> = GRANDFATHERED
        .iter()
        .filter(|(actual, _, _)| architecture.packages.contains_key(*actual))
        .map(|(actual, _, _)| *actual)
        .collect();

    for package in architecture.packages.values() {
        if package.layer == "tools" {
            continue;
        }
        let exception = GRANDFATHERED.iter().find(|(a, _, _)| *a == package.name);
        let mut violation: Option<String> = None;

        // Rule 1: a contract's name is fully determined by its domain.
        if let Some(expected) = expected_name(package) {
            if package.name != expected {
                violation = Some(format!(
                    "{}: role `{}` with domain `{}` implies name `{}`",
                    package.name, package.role, package.domain, expected
                ));
            }
        }

        // Rule 2: a provider carries its capability domain AND an engine suffix.
        if violation.is_none() && package.role.starts_with("provider.") {
            let prefix = format!("axiolid-{}-", domain_slug(&package.domain));
            if !package.name.starts_with(&prefix) || package.name.len() <= prefix.len() {
                violation = Some(format!(
                    "{}: provider for domain `{}` must be named `{}<engine>`",
                    package.name, package.domain, prefix
                ));
            }
        }

        // Rule 3: only a contract may claim `-contract`, and nobody may take
        // the bare capability name a contract reserves.
        if violation.is_none() && package.role != "contract.operation" {
            if package.name.ends_with("-contract") {
                violation = Some(format!(
                    "{}: only `contract.operation` packages may be named `-contract` (role is `{}`)",
                    package.name, package.role
                ));
            } else if let Some((_, owner)) = reserved.iter().find(|(bare, _)| *bare == package.name)
            {
                violation = Some(format!(
                    "{}: name is reserved by contract `{}`; an implementation needs an engine suffix",
                    package.name, owner
                ));
            }
        }

        match (violation, exception) {
            // A known exception, still violating: allowed, and recorded as used.
            (Some(_), Some((actual, _, _))) => {
                unused_exceptions.retain(|name| name != actual);
            }
            // A new violation.
            (Some(message), None) => errors.push(message),
            // Listed as an exception but no longer violating: the list is stale.
            (None, Some((actual, expected, _))) => {
                unused_exceptions.retain(|name| name != actual);
                errors.push(format!(
                    "{actual}: listed as a grandfathered naming exception (expected `{expected}`) but now complies; remove it from GRANDFATHERED"
                ));
            }
            (None, None) => {}
        }
    }

    for name in unused_exceptions {
        errors.push(format!(
            "{name}: grandfathered naming exception is unused; remove it from GRANDFATHERED"
        ));
    }
    Ok(())
}
