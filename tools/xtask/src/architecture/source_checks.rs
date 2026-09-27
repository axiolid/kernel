use super::model::{Architecture, Result};
use rustc_lexer::{tokenize, TokenKind};
use std::fs;
use std::path::{Path, PathBuf};

// Production source is intentionally format-neutral. Tests may name external
// formats when proving that adapters and architecture gates reject leakage.
const FORBIDDEN_FORMAT_TERMS: &[&str] = &[
    "ifc", "revit", "solibri", "cset", "pkl", "protobuf", "prost",
];
const FORBIDDEN_NATIVE_DEPENDENCIES: &[&str] = &["bindgen", "cxx", "autocxx", "cmake"];

pub fn validate(architecture: &Architecture, errors: &mut Vec<String>) -> Result<()> {
    for package in architecture.packages.values() {
        if package.layer == "tools" {
            continue;
        }
        let src = architecture.root.join(&package.path).join("src");
        if !src.exists() {
            errors.push(format!("{}: missing src directory", package.name));
            continue;
        }
        let files = rust_files(&src)?;
        for file in &files {
            let source = fs::read_to_string(file)
                .map_err(|error| format!("read {}: {error}", file.display()))?;
            // Comments are explanatory prose, not executable coupling. Lex first so
            // identifiers and string literals remain enforceable without rejecting docs.
            let code = source_without_comments(&source);
            let lower = code.to_ascii_lowercase();
            if package.format_neutral {
                for term in FORBIDDEN_FORMAT_TERMS {
                    if lower.contains(term) {
                        errors.push(format!(
                            "{}: {} contains forbidden source-format term `{term}`",
                            package.name,
                            relative(&architecture.root, file)
                        ));
                    }
                }
            }
            for placeholder in ["todo!(", "unimplemented!("] {
                if code.contains(placeholder) {
                    errors.push(format!(
                        "{}: {} contains panicking placeholder `{placeholder}`",
                        package.name,
                        relative(&architecture.root, file)
                    ));
                }
            }
            validate_declared_module(&src, file, &source, errors)?;
        }
        // Use Cargo's resolved declaration metadata, not manifest key spelling: a
        // dependency renamed with `package = "prost"` must remain detectable.
        for dependency in &package.declared_dependency_packages {
            let lower = dependency.to_ascii_lowercase();
            if FORBIDDEN_NATIVE_DEPENDENCIES.contains(&lower.as_str()) {
                errors.push(format!(
                    "{}: forbidden native bridge dependency `{dependency}`",
                    package.name
                ));
            }
            if package.format_neutral {
                for term in FORBIDDEN_FORMAT_TERMS {
                    if lower.contains(term) {
                        errors.push(format!(
                            "{}: forbidden source-format dependency `{dependency}` contains `{term}`",
                            package.name
                        ));
                    }
                }
            }
        }
        let unsafe_boundary = package.name == "axiolid-capi";
        let forbids_unsafe = (matches!(
            package.layer.as_str(),
            "foundation" | "representations" | "contracts" | "facade"
        ) || package.name == "axiolid-backend-gpu")
            && !unsafe_boundary;
        let lib = src.join("lib.rs");
        if forbids_unsafe {
            let source = fs::read_to_string(&lib)
                .map_err(|error| format!("read {}: {error}", lib.display()))?;
            if !source.contains("#![forbid(unsafe_code)]") {
                errors.push(format!(
                    "{}: {} must forbid unsafe code",
                    package.name,
                    relative(&architecture.root, &lib)
                ));
            }
        } else if unsafe_boundary {
            let source = fs::read_to_string(&lib)
                .map_err(|error| format!("read {}: {error}", lib.display()))?;
            if !source.contains("#![deny(unsafe_op_in_unsafe_fn)]") {
                errors.push(format!(
                    "{}: {} must deny unsafe operations outside explicit unsafe blocks",
                    package.name,
                    relative(&architecture.root, &lib)
                ));
            }
        }
    }
    Ok(())
}

fn source_without_comments(source: &str) -> String {
    let mut code = String::with_capacity(source.len());
    let mut offset = 0;
    for token in tokenize(source) {
        let end = offset + token.len;
        match token.kind {
            TokenKind::LineComment | TokenKind::BlockComment { .. } => code.push(' '),
            _ => code.push_str(&source[offset..end]),
        }
        offset = end;
    }
    debug_assert_eq!(offset, source.len());
    code
}

fn rust_files(root: &Path) -> Result<Vec<PathBuf>> {
    fn visit(dir: &Path, files: &mut Vec<PathBuf>) -> Result<()> {
        for entry in
            fs::read_dir(dir).map_err(|error| format!("read {}: {error}", dir.display()))?
        {
            let path = entry
                .map_err(|error| format!("read directory entry: {error}"))?
                .path();
            if path.is_dir() {
                visit(&path, files)?;
            } else if path.extension().and_then(|value| value.to_str()) == Some("rs") {
                files.push(path);
            }
        }
        Ok(())
    }
    let mut files = Vec::new();
    visit(root, &mut files)?;
    files.sort();
    Ok(files)
}

fn validate_declared_module(
    src: &Path,
    file: &Path,
    _source: &str,
    errors: &mut Vec<String>,
) -> Result<()> {
    if file.file_name().and_then(|value| value.to_str()) == Some("lib.rs")
        || file
            .parent()
            .and_then(Path::file_name)
            .and_then(|value| value.to_str())
            == Some("bin")
        || file.components().any(|part| part.as_os_str() == "tests")
    {
        return Ok(());
    }
    let stem = file
        .file_stem()
        .and_then(|value| value.to_str())
        .ok_or_else(|| format!("{} has no UTF-8 module stem", file.display()))?;
    let parent = file
        .parent()
        .ok_or_else(|| format!("{} has no parent", file.display()))?;
    let in_bin_directory = parent
        .parent()
        .and_then(Path::file_name)
        .and_then(|value| value.to_str())
        == Some("bin");
    if in_bin_directory && stem == "main" {
        return Ok(());
    }
    let candidate = if parent == src {
        src.join("lib.rs")
    } else if in_bin_directory {
        // A multi-file binary, `src/bin/<name>/main.rs`, declares its
        // sibling modules from `main.rs`.
        parent.join("main.rs")
    } else {
        parent.with_extension("rs")
    };
    let declaration = if candidate.exists() {
        candidate
    } else {
        parent.join("mod.rs")
    };
    let declaring_source = fs::read_to_string(&declaration)
        .map_err(|error| format!("read {}: {error}", declaration.display()))?;
    if !declaring_source.contains(&format!("mod {stem};"))
        && !declaring_source.contains(&format!("pub mod {stem};"))
    {
        errors.push(format!(
            "{} is not declared by {}",
            file.display(),
            declaration.display()
        ));
    }
    Ok(())
}

fn relative(root: &Path, path: &Path) -> String {
    path.strip_prefix(root)
        .unwrap_or(path)
        .display()
        .to_string()
}

// Downstream product names must not appear anywhere in the repository, not just
// in Rust sources. Docs, PLAN files and scripts are published too, and naming a
// specific end-user product in a public infrastructure repo leaks a commercial
// relationship the kernel is not supposed to model.
//
// The two files that carry the banned list itself are exempt: they are the
// machinery that enforces this rule.
const FORBIDDEN_PRODUCT_TERMS: &[&str] = &["solibri"];

const PRODUCT_TERM_EXEMPT_SUFFIXES: &[&str] = &[
    "tools/xtask/src/architecture/source_checks.rs",
    "crates/algorithms/sampled/field/tests/neutrality.rs",
];

pub fn validate_no_product_names(root: &Path, errors: &mut Vec<String>) -> Result<()> {
    let tracked = std::process::Command::new("git")
        .args(["ls-files", "-z"])
        .current_dir(root)
        .output()
        .map_err(|error| format!("git ls-files: {error}"))?;
    if !tracked.status.success() {
        return Err("git ls-files failed".to_owned());
    }
    let listing = String::from_utf8_lossy(&tracked.stdout);
    for entry in listing.split('\0').filter(|entry| !entry.is_empty()) {
        let normalized = entry.replace('\\', "/");
        if PRODUCT_TERM_EXEMPT_SUFFIXES
            .iter()
            .any(|exempt| normalized.ends_with(exempt) || normalized == *exempt)
        {
            continue;
        }
        let path = root.join(entry);
        let Ok(source) = fs::read_to_string(&path) else {
            continue; // binary or unreadable: nothing textual to leak
        };
        let haystack = format!("{}\n{normalized}", source.to_ascii_lowercase());
        for term in FORBIDDEN_PRODUCT_TERMS {
            if haystack.contains(term) {
                errors.push(format!(
                    "{normalized}: contains downstream product name `{term}`"
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::source_without_comments;

    #[test]
    fn removes_line_doc_and_nested_block_comments() {
        let source = r#"// protobuf
//! prost docs
const OK: &str = "neutral"; /* outer prost /* nested protobuf */ */"#;
        let code = source_without_comments(source);
        assert!(!code.contains("protobuf"));
        assert!(!code.contains("prost"));
        assert!(code.contains("const OK"));
    }

    #[test]
    fn preserves_identifiers_and_string_literals() {
        let source = r##"const PROTOBUF: &str = r#"prost // not a comment"#;"##;
        let code = source_without_comments(source);
        assert!(code.contains("PROTOBUF"));
        assert!(code.contains("prost // not a comment"));
    }
}
