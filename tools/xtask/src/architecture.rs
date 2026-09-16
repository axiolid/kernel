mod checks;
mod closure;
mod model;
mod naming;
mod render;
mod source_checks;

use model::{Architecture, Result};

pub fn check() -> Result<()> {
    let architecture = Architecture::load()?;
    checks::validate(&architecture)?;
    render::check_docs(&architecture)?;
    println!(
        "architecture check passed: {} packages",
        architecture.packages.len()
    );
    Ok(())
}

pub fn list() -> Result<()> {
    let architecture = Architecture::load()?;
    checks::validate(&architecture)?;
    print!("{}", render::list(&architecture));
    Ok(())
}

pub fn graph() -> Result<()> {
    let architecture = Architecture::load()?;
    checks::validate(&architecture)?;
    print!("{}", render::graph(&architecture));
    Ok(())
}

pub fn closure_check() -> Result<()> {
    let root = closure::workspace_root()?;
    closure::check(&root)
}

pub fn closure_docs() -> Result<()> {
    let root = closure::workspace_root()?;
    closure::docs(&root)
}

pub fn closure_explain(name: &str) -> Result<()> {
    let root = closure::workspace_root()?;
    closure::explain(&root, name)
}

pub fn docs() -> Result<()> {
    let architecture = Architecture::load()?;
    checks::validate(&architecture)?;
    render::write_docs(&architecture)?;
    println!("generated architecture documentation");
    Ok(())
}
