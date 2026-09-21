//! `default = []` is a deliberate 1.0 decision, not an accident (kernel#9).
//!
//! The facade compiles no geometry unless the consumer names a capability.
//! A `compile_error!` in `lib.rs` explains that when nothing is selected --
//! but a compile failure cannot be asserted from inside the crate that
//! fails, so this test guards the MANIFEST invariants the diagnostic
//! depends on. If someone restores a non-empty default, or drops the
//! `standard` escape hatch, the decision breaks silently without this.

use std::path::Path;

fn manifest() -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    std::fs::read_to_string(&path).expect("the facade manifest is readable")
}

/// The default feature set is empty.
#[test]
fn the_default_feature_set_is_empty() {
    let text = manifest();
    assert!(
        text.contains("default = []"),
        "the facade must default to no capabilities: a non-empty default makes \\
         every consumer pay for geometry they did not ask for (kernel#9)"
    );
}

/// `standard` reproduces the pre-0.4 default exactly.
///
/// The `compile_error!` tells upgrading consumers to add this one feature.
/// If its contents drift from the old default, that instruction silently
/// becomes wrong and the migration stops being one line.
#[test]
fn standard_reproduces_the_previous_default() {
    let text = manifest();
    assert!(
        text.contains(r#"standard = ["mesh", "cpu", "integration"]"#),
        "`standard` must stay byte-equal to the pre-0.4 default, because the \\
         migration note in the compile_error names it as the drop-in"
    );
}

/// Every feature the diagnostic advertises is a real feature.
///
/// The message is aimed at coding agents, which will paste these names
/// straight into a manifest. A name that no longer exists would send them
/// into a second, worse error, so the advice is checked against the
/// manifest rather than trusted to stay current by hand.
#[test]
fn the_diagnostic_only_recommends_features_that_exist() {
    let text = manifest();
    let advertised = [
        "standard",
        "mesh",
        "linear-intersection",
        "measure",
        "ray-mesh",
        "application",
        "nurbs",
        "brep",
        "tessellation",
        "pointcloud-provider",
        "field-ops",
        "heal",
        "discrete",
        "parametric",
        "advanced",
        "full",
        "mesh-boolean",
    ];
    for name in advertised {
        assert!(
            text.contains(&format!("{name} = [")),
            "the compile_error recommends `{name}`, but the manifest defines no \\
             such feature -- an agent following that advice would hit a second error"
        );
    }
}
