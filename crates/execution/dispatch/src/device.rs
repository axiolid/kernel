//! Device preference matching, shared by every operation registry.
//!
//! One implementation so a new operation cannot drift from the established
//! routing semantics: adding a registry that quietly treats `Cpu` as
//! portable-only would silently exclude optimised backends.

use axiolid_contracts::{BackendId, DevicePreference, ExecutionTarget};

/// Whether a backend satisfies the caller's device preference.
pub(crate) fn matches_device(
    preference: DevicePreference,
    id: BackendId,
    target: ExecutionTarget,
) -> bool {
    match preference {
        DevicePreference::Auto => true,
        DevicePreference::Cpu => {
            matches!(
                target,
                ExecutionTarget::PortableCpu | ExecutionTarget::OptimizedCpu
            )
        }
        DevicePreference::Gpu => matches!(target, ExecutionTarget::Gpu),
        DevicePreference::Backend(required) => required == id,
    }
}
