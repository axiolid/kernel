//! The publishable crates, as the docs generators see them.
//!
//! Layer, role and internal dependencies come from the architecture model,
//! so the docs and `cargo xtask architecture check` read the same
//! `[package.metadata.axiolid]` through the same code.

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

use cargo_metadata::MetadataCommand;

use super::Result;
use crate::architecture::model::Architecture;

/// The layers a publishable crate may declare, in reading order, with the
/// heading each gets on the index. A crate in any other layer fails the run,
/// so a new layer cannot silently drop out of the reference.
pub(super) const LAYERS: &[(&str, &str)] = &[
    ("foundation", "Foundation"),
    ("representations", "Representations"),
    ("contracts", "Contracts"),
    ("algorithms", "Algorithms"),
    ("providers", "Providers"),
    ("execution", "Execution"),
    ("facade", "Facade"),
];

/// The facade package, whose features expose the other crates.
pub(super) const FACADE: &str = "axiolid";

pub(super) struct Crate {
    pub(super) name: String,
    pub(super) description: String,
    /// The version on `main`, which may be ahead of the latest release.
    pub(super) version: String,
    /// The crate directory relative to the workspace root.
    pub(super) dir: String,
    pub(super) layer: String,
    pub(super) role: String,
    /// The library's Rust name, which names its rustdoc directory.
    pub(super) lib: Option<String>,
    pub(super) features: BTreeMap<String, Vec<String>>,
    /// Workspace crates this one depends on outside dev-dependencies.
    pub(super) internal_deps: BTreeSet<String>,
}

pub(super) struct Workspace {
    pub(super) root: PathBuf,
    /// Publishable crates, sorted by name.
    pub(super) crates: Vec<Crate>,
    pub(super) architecture: Architecture,
}

impl Workspace {
    pub(super) fn load() -> Result<Self> {
        let metadata = MetadataCommand::new()
            .no_deps()
            .exec()
            .map_err(|error| format!("cargo metadata failed: {error}"))?;
        let architecture = Architecture::from_metadata(&metadata)?;
        let mut crates = Vec::new();
        for package in metadata.workspace_packages() {
            let name = package.name.to_string();
            let Some(model) = architecture.packages.get(&name) else {
                continue;
            };
            if !model.publish {
                continue;
            }
            if !LAYERS.iter().any(|(layer, _)| *layer == model.layer) {
                return Err(format!(
                    "{name}: layer `{}` has no section in the crate reference; add it to docs::workspace::LAYERS",
                    model.layer
                ));
            }
            let description = package
                .description
                .as_deref()
                .map(str::trim)
                .filter(|text| !text.is_empty())
                .ok_or_else(|| format!("{name}: a published crate needs a `description`"))?
                .to_owned();
            let lib = package
                .targets
                .iter()
                .find(|target| {
                    target.is_lib()
                        || target.is_rlib()
                        || target.is_cdylib()
                        || target.is_staticlib()
                })
                .map(|target| target.name.replace('-', "_"));
            crates.push(Crate {
                name,
                description,
                version: package.version.to_string(),
                dir: model.path.clone(),
                layer: model.layer.clone(),
                role: model.role.clone(),
                lib,
                features: package.features.clone(),
                internal_deps: model.production_internal_dependencies.clone(),
            });
        }
        crates.sort_by(|a, b| a.name.cmp(&b.name));
        if !crates.iter().any(|c| c.name == FACADE) {
            return Err(format!("the facade `{FACADE}` is not a published crate"));
        }
        Ok(Self {
            root: architecture.root.clone(),
            crates,
            architecture,
        })
    }

    pub(super) fn get(&self, name: &str) -> Option<&Crate> {
        self.crates.iter().find(|c| c.name == name)
    }

    /// Publishable crates in one layer, sorted by name.
    pub(super) fn layer(&self, layer: &str) -> Vec<&Crate> {
        self.crates.iter().filter(|c| c.layer == layer).collect()
    }
}
