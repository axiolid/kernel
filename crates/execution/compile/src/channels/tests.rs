//! Unit tests for fate bookkeeping (#115).
//!
//! The graph-level behaviour is pinned end to end in
//! `tests/graph_channels.rs`; these cover the two combination rules that
//! the integration tests reach only indirectly.

use axiolid_mesh::{AttributeFate, DropReason, TriMesh};

use super::{after_boolean, Fates};

fn dropped(reason: DropReason) -> AttributeFate {
    AttributeFate::Dropped(reason)
}

/// Parallel inputs combine worst-of, in either order: a channel dropped in
/// one Collection member is untrustworthy as a whole.
#[test]
fn parallel_fates_keep_the_worst_in_either_order() {
    let worst = dropped(DropReason::IncompatibleChannels);
    for order in [
        [AttributeFate::Preserved, worst.clone()],
        [worst.clone(), AttributeFate::Preserved],
    ] {
        let mut fates = Fates::default();
        for fate in order {
            fates.record("uv", fate);
        }
        assert_eq!(fates.into_vec(), vec![("uv".to_owned(), worst.clone())]);
    }
    let mut fates = Fates::default();
    fates.record("uv", AttributeFate::Interpolated);
    fates.record("uv", AttributeFate::Preserved);
    assert_eq!(
        fates.into_vec(),
        vec![("uv".to_owned(), AttributeFate::Interpolated)],
        "interpolated outranks preserved"
    );
}

/// Sequential composition keeps the FIRST drop reason: the step that lost
/// the data. A boolean reporting a later reason must not overwrite it.
#[test]
fn a_boolean_keeps_an_upstream_drop_reason() {
    let mut subject = Fates::default();
    subject.record("uv", dropped(DropReason::IncompatibleChannels));
    let built = after_boolean(
        TriMesh::default(),
        &subject,
        None,
        vec![("uv".to_owned(), dropped(DropReason::NotBlendable))],
    );
    assert_eq!(
        built.fates.into_vec(),
        vec![("uv".to_owned(), dropped(DropReason::IncompatibleChannels))]
    );
}

/// A channel tracked upstream that the provider no longer mentions was lost
/// in the boolean. Silence is not preservation.
#[test]
fn a_channel_the_boolean_omits_is_reported_dropped() {
    let mut subject = Fates::default();
    subject.record("uv", AttributeFate::Preserved);
    let built = after_boolean(TriMesh::default(), &subject, None, Vec::new());
    assert_eq!(
        built.fates.into_vec(),
        vec![("uv".to_owned(), dropped(DropReason::ProviderLimitation))]
    );
}

/// The tool's history counts too: a channel the tool had dropped upstream
/// stays dropped even though the subject preserved it.
#[test]
fn the_tools_upstream_fates_are_folded_in() {
    let mut subject = Fates::default();
    subject.record("uv", AttributeFate::Preserved);
    let mut tool = Fates::default();
    tool.record("uv", dropped(DropReason::IncompatibleChannels));
    let built = after_boolean(
        TriMesh::default(),
        &subject,
        Some(&tool),
        vec![("uv".to_owned(), AttributeFate::Preserved)],
    );
    assert_eq!(
        built.fates.into_vec(),
        vec![("uv".to_owned(), dropped(DropReason::IncompatibleChannels))]
    );
}
