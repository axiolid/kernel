//! `docs/.vitepress/data/facts.json`: the crate facts the site config reads.
//!
//! The Reference sidebar is built from `crates.groups`, so a new crate or
//! layer reaches the sidebar with no edit to `config.ts`.

use super::changelog;
use super::workspace::{Workspace, LAYERS};
use super::{Output, Result};

const TARGET: &str = "docs/.vitepress/data/facts.json";

pub(super) fn generate(workspace: &Workspace) -> Result<Output> {
    let mut per_crate = serde_json::Map::new();
    for krate in &workspace.crates {
        let release = changelog::latest(workspace, krate)?;
        per_crate.insert(
            krate.name.clone(),
            serde_json::json!({
                "description": krate.description,
                "layer": krate.layer,
                "role": krate.role,
                "version": krate.version,
                "released": release.as_ref().map(|r| r.version.clone()),
                "released_date": release.as_ref().map(|r| r.date.clone()),
                "path": krate.dir,
            }),
        );
    }
    let groups: Vec<serde_json::Value> = LAYERS
        .iter()
        .map(|(key, title)| {
            let members: Vec<&str> = workspace
                .layer(key)
                .into_iter()
                .map(|c| c.name.as_str())
                .collect();
            serde_json::json!({ "key": key, "title": title, "crates": members })
        })
        .collect();
    let facts = serde_json::json!({
        "crates": {
            "total": workspace.crates.len(),
            "groups": groups,
        },
        "crate": per_crate,
    });
    let mut json = serde_json::to_string_pretty(&facts).map_err(|error| error.to_string())?;
    json.push('\n');
    Ok(Output::new(workspace, TARGET, json))
}
