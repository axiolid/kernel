//! Stamp a SONAME (ELF) / install name (Mach-O) onto the cdylib.
//!
//! Rust's `crate-type = ["cdylib"]` emits a shared object with no SONAME.
//! Without one, a consumer's `DT_NEEDED` entry records the path the
//! linker happened to see -- often relative, e.g. `./libaxiolid_capi.so`.
//! Per ELF semantics a `NEEDED` value containing `/` is treated as a
//! literal filesystem path and is NEVER searched for via RPATH/RUNPATH,
//! so the consumer fails to load from any directory but the build tree
//! (axiolid/kernel#70).
//!
//! Stamping the bare filename makes the dependency relocatable: the
//! loader resolves it through the normal search path, which is what
//! every distro-packaged shared library relies on.

fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    // Only the cdylib carries a soname; rlib/staticlib are unaffected.
    let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    match target_os.as_str() {
        // ELF: -soname is read by the consumer's linker and copied
        // verbatim into its DT_NEEDED entry.
        "linux" | "android" | "freebsd" | "netbsd" | "openbsd" | "dragonfly" => {
            println!("cargo:rustc-cdylib-link-arg=-Wl,-soname,libaxiolid_capi.so");
        }
        // Mach-O: the install name plays the same role. `@rpath/` keeps
        // it relocatable instead of baking an absolute build path.
        "macos" | "ios" => {
            println!("cargo:rustc-cdylib-link-arg=-Wl,-install_name,@rpath/libaxiolid_capi.dylib");
        }
        // Windows has no soname concept: the import library carries the
        // DLL name already, so there is nothing to stamp.
        _ => {}
    }
}
