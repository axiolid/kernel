//! A boolean must report what happened to attribute channels.
//!
//! Before this, a channel simply vanished: the caller compared `is_some()`
//! before and after to discover the loss, and got no reason for it. The cut
//! genuinely cannot preserve per-vertex data -- new vertices along the seam
//! have no preimage in either operand -- so the contract is not preservation
//! but disclosure.

mod support;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_mesh::{AttributeChannel, AttributeFate, Blend, DropReason};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
use support::boxx;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

/// A dropped channel is named, with the reason it could not survive.
#[test]
fn a_dropped_channel_is_reported_rather_than_silently_lost() {
    let mut subject = boxx(0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 0.0);
    let tool = boxx(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 0.0);

    let vertices = subject.positions.len();
    subject.attributes.push(AttributeChannel::new(
        "source_entity",
        (0..vertices).map(|i| i as f64).collect(),
        1,
        // A source-entity handle is a label: averaging two handles would
        // name an entity that does not exist.
        Blend::Nearest,
    ));

    let outcome = BoolmeshBoolean::new()
        .boolean(&subject, &tool, BooleanOperator::Difference, &options())
        .expect("difference");

    assert_eq!(
        outcome.evidence.attribute_fates,
        vec![(
            "source_entity".to_owned(),
            AttributeFate::Dropped(DropReason::ProviderLimitation)
        )],
        "the channel must be named as dropped, with a reason"
    );
}

/// A channel that forbids derivation reports that, not a provider limit.
///
/// The distinction matters to a caller deciding whether to retry: a
/// provider limitation might be lifted by a better backend, whereas
/// `Blend::None` is the data itself saying no value is derivable.
#[test]
fn a_non_blendable_channel_names_the_data_not_the_provider() {
    let mut subject = boxx(0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 0.0);
    let tool = boxx(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 0.0);

    let vertices = subject.positions.len();
    subject.attributes.push(AttributeChannel::new(
        "opaque",
        (0..vertices).map(|i| i as f64).collect(),
        1,
        Blend::None,
    ));

    let outcome = BoolmeshBoolean::new()
        .boolean(&subject, &tool, BooleanOperator::Difference, &options())
        .expect("difference");

    assert_eq!(
        outcome.evidence.attribute_fates,
        vec![(
            "opaque".to_owned(),
            AttributeFate::Dropped(DropReason::NotBlendable)
        )]
    );
}

/// A subject with no channels reports an empty list, not a fabricated one.
#[test]
fn a_subject_without_channels_reports_nothing() {
    let subject = boxx(0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 0.0);
    let tool = boxx(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 0.0);

    let outcome = BoolmeshBoolean::new()
        .boolean(&subject, &tool, BooleanOperator::Difference, &options())
        .expect("difference");

    assert!(outcome.evidence.attribute_fates.is_empty());
}

/// A corner-indexed channel is named like any other (#112).
///
/// Addressing is not a reason to go quiet: the cut drops it for the same
/// reason it drops a per-vertex channel, and says so.
#[test]
fn a_corner_indexed_channel_is_reported_too() {
    let mut subject = boxx(0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 0.0);
    let tool = boxx(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 0.0);

    let corners = subject.indices.len() as u32;
    subject.attributes.push(AttributeChannel::corner_indexed(
        "uv",
        (0..corners).flat_map(|c| [f64::from(c), 0.0]).collect(),
        2,
        Blend::Linear,
        (0..corners).collect(),
    ));
    subject.validate_structure().expect("valid operand");

    let outcome = BoolmeshBoolean::new()
        .boolean(&subject, &tool, BooleanOperator::Difference, &options())
        .expect("difference");

    assert_eq!(
        outcome.evidence.attribute_fates,
        vec![(
            "uv".to_owned(),
            AttributeFate::Dropped(DropReason::ProviderLimitation)
        )]
    );
}
