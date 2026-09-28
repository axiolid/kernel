//! One generated reference page per publishable crate, and the index that
//! groups them by layer.

use std::collections::BTreeMap;

use super::changelog::{self, Release};
use super::readme;
use super::workspace::{Crate, Workspace, FACADE, LAYERS};
use super::{Output, Result, BANNER, REPO};

pub(super) fn generate(workspace: &Workspace) -> Result<Vec<Output>> {
    let facade = workspace
        .get(FACADE)
        .ok_or_else(|| format!("the facade `{FACADE}` is not a published crate"))?;
    let mut releases = BTreeMap::new();
    for krate in &workspace.crates {
        releases.insert(krate.name.as_str(), changelog::latest(workspace, krate)?);
    }
    let mut outputs = vec![Output::new(
        workspace,
        "docs/reference/index.md",
        index(workspace, &releases),
    )];
    for krate in &workspace.crates {
        let page = page(
            workspace,
            krate,
            facade,
            releases[krate.name.as_str()].as_ref(),
        )?;
        outputs.push(Output::new(
            workspace,
            format!("docs/reference/crates/{}.md", krate.name),
            page,
        ));
    }
    Ok(outputs)
}

fn released(release: Option<&Release>) -> String {
    release.map_or_else(|| "not released".to_owned(), |r| r.version.clone())
}

fn index(workspace: &Workspace, releases: &BTreeMap<&str, Option<Release>>) -> String {
    let mut out = vec![
        BANNER.to_owned(),
        String::new(),
        "# Crate reference".to_owned(),
        String::new(),
        format!(
            "Axiolid publishes {} crates, each versioned and released on its own. The \
             [`{FACADE}`](./crates/{FACADE}) facade re-exports them behind features; every \
             crate can also be used directly. [Selecting a package](./selecting-packages) \
             says which to start from.",
            workspace.crates.len()
        ),
        String::new(),
        "Every page here is generated from the crate itself: its manifest, its crate \
         documentation and its changelog. Sections follow the layers of the \
         [crate map](/architecture/crate-map)."
            .to_owned(),
    ];
    for (layer, title) in LAYERS {
        let members = workspace.layer(layer);
        if members.is_empty() {
            continue;
        }
        out.extend([
            String::new(),
            format!("## {title}"),
            String::new(),
            "| Crate | Role | Latest release | Description |".to_owned(),
            "| --- | --- | --- | --- |".to_owned(),
        ]);
        for krate in members {
            out.push(format!(
                "| [`{name}`](./crates/{name}) | `{}` | {} | {} |",
                krate.role,
                released(releases[krate.name.as_str()].as_ref()),
                krate.description,
                name = krate.name
            ));
        }
    }
    out.push(String::new());
    out.join("\n")
}

fn page(
    workspace: &Workspace,
    krate: &Crate,
    facade: &Crate,
    release: Option<&Release>,
) -> Result<String> {
    let mut out = vec![
        BANNER.to_owned(),
        String::new(),
        format!("# {}", krate.name),
        String::new(),
        format!("{}.", krate.description.trim_end_matches('.')),
        String::new(),
        "| | |".to_owned(),
        "| --- | --- |".to_owned(),
    ];
    match release {
        Some(release) => {
            out.push(format!(
                "| Latest release | {} ({}) |",
                release.version, release.date
            ));
            if release.version != krate.version {
                out.push(format!("| On `main` | {} (unreleased) |", krate.version));
            }
            out.push(format!(
                "| crates.io | [`{0}`](https://crates.io/crates/{0}) |",
                krate.name
            ));
        }
        None => out.push(format!(
            "| Latest release | not released (`main` is {}) |",
            krate.version
        )),
    }
    let exposed = facade_features(facade, &krate.name);
    if !exposed.is_empty() {
        out.push(format!(
            "| Facade | [`{FACADE}`](./{FACADE}) feature {} |",
            exposed.join(", ")
        ));
    }
    out.push(format!("| Layer | {} (`{}`) |", krate.layer, krate.role));
    let mut api = Vec::new();
    if let Some(lib) = &krate.lib {
        // Built into the site by the docs workflow; VitePress adds the base.
        api.push(format!("[rustdoc](/api/rustdoc/{lib}/index.html)"));
    }
    if release.is_some() {
        api.push(format!("[docs.rs](https://docs.rs/{})", krate.name));
    }
    if !api.is_empty() {
        out.push(format!("| API documentation | {} |", api.join(" · ")));
    }
    out.push(format!(
        "| Source | [`{dir}/`]({REPO}/tree/main/{dir}) |",
        dir = krate.dir
    ));

    let prose = readme::prose(workspace, krate)?;
    if !prose.overview.is_empty() {
        out.extend([
            String::new(),
            "## Overview".to_owned(),
            String::new(),
            prose.overview,
        ]);
    }
    if !prose.notes.is_empty() {
        out.extend([String::new(), prose.notes]);
    }
    features(workspace, krate, &mut out);
    if !krate.internal_deps.is_empty() {
        out.extend([String::new(), "## Depends on".to_owned(), String::new()]);
        for dep in &krate.internal_deps {
            out.push(format!("- {}", crate_link(workspace, dep)));
        }
    }
    out.extend([String::new(), "## Changes".to_owned(), String::new()]);
    match release {
        Some(Release {
            version,
            date,
            notes: Some(notes),
        }) => {
            out.push(format!("Latest release, {version} ({date}):"));
            out.push(String::new());
            out.push(notes.clone());
        }
        Some(Release { version, date, .. }) => out.push(format!(
            "Released in the workspace-wide {version} release ({date}), before crates \
             versioned independently; its notes are in the [workspace changelog](/CHANGELOG)."
        )),
        None => out.push("No release yet.".to_owned()),
    }
    out.push(String::new());
    out.push(format!(
        "Full history: [`{dir}/CHANGELOG.md`]({REPO}/blob/main/{dir}/CHANGELOG.md)",
        dir = krate.dir
    ));
    out.push(String::new());
    Ok(out.join("\n"))
}

/// The facade features that enable `name` directly, as inline code.
fn facade_features(facade: &Crate, name: &str) -> Vec<String> {
    let dep = format!("dep:{name}");
    facade
        .features
        .iter()
        .filter(|(_, values)| values.contains(&dep))
        .map(|(feature, _)| format!("`{feature}`"))
        .collect()
}

/// The crate's own features, with what each turns on.
fn features(workspace: &Workspace, krate: &Crate, out: &mut Vec<String>) {
    let defaults = krate.features.get("default").cloned().unwrap_or_default();
    let listed: Vec<(&String, &Vec<String>)> = krate
        .features
        .iter()
        .filter(|(feature, _)| feature.as_str() != "default")
        .collect();
    if listed.is_empty() {
        return;
    }
    out.extend([
        String::new(),
        "## Features".to_owned(),
        String::new(),
        format!(
            "Default: {}.",
            if defaults.is_empty() {
                "none".to_owned()
            } else {
                defaults
                    .iter()
                    .map(|f| format!("`{f}`"))
                    .collect::<Vec<_>>()
                    .join(", ")
            }
        ),
        String::new(),
        "| Feature | Enables |".to_owned(),
        "| --- | --- |".to_owned(),
    ]);
    for (feature, values) in listed {
        let enables: Vec<String> = values
            .iter()
            .map(|value| match value.strip_prefix("dep:") {
                Some(dep) => crate_link(workspace, dep),
                None => format!("`{value}`"),
            })
            .collect();
        let enables = if enables.is_empty() {
            "—".to_owned()
        } else {
            enables.join(", ")
        };
        out.push(format!("| `{feature}` | {enables} |"));
    }
}

/// A workspace crate links to its page; anything else is plain code.
fn crate_link(workspace: &Workspace, name: &str) -> String {
    if workspace.get(name).is_some() {
        format!("[`{name}`](./{name})")
    } else {
        format!("`{name}`")
    }
}
