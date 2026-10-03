//! `TrimSelector::ArcLength` (#239): a parameter-kind selector, finite or
//! refused, and stored exactly.

use axiolid_core::Vec2;
use axiolid_curve::{Curve2, Line2};
use axiolid_model::{
    CurveRelation, GeometryGraphBuilder, GraphError, OpenProfile, TrimSelector, TrimmingPreference,
};

fn line() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Vec2::ZERO,
        direction: Vec2::X,
    })
}

fn open_trim(
    start: Vec<TrimSelector>,
    end: Vec<TrimSelector>,
    preference: TrimmingPreference,
) -> Result<(), GraphError> {
    let mut builder = GeometryGraphBuilder::new();
    let basis = builder.push_value(line()).unwrap();
    let trimmed = builder
        .push_value(CurveRelation::Trimmed {
            basis,
            start,
            end,
            sense_agreement: true,
            preference,
        })
        .unwrap();
    builder.push_value(OpenProfile::new(trimmed)).map(|_| ())
}

#[test]
fn an_arc_length_trim_is_an_open_profile_under_every_preference_that_admits_it() {
    for preference in [
        TrimmingPreference::Parameter,
        TrimmingPreference::Unspecified,
    ] {
        open_trim(
            vec![TrimSelector::ArcLength(0.0)],
            vec![TrimSelector::ArcLength(12.5)],
            preference,
        )
        .unwrap_or_else(|error| panic!("{preference:?}: {error:?}"));
    }
    // Mixed with a parameter selector, as a source may state both.
    open_trim(
        vec![TrimSelector::Parameter(0.0)],
        vec![TrimSelector::ArcLength(3.0)],
        TrimmingPreference::Parameter,
    )
    .unwrap();
}

#[test]
fn an_arc_length_does_not_satisfy_a_cartesian_preference() {
    assert!(open_trim(
        vec![TrimSelector::ArcLength(0.0)],
        vec![TrimSelector::ArcLength(1.0)],
        TrimmingPreference::Cartesian,
    )
    .is_err());
}

#[test]
fn a_non_finite_arc_length_is_refused() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(open_trim(
            vec![TrimSelector::ArcLength(0.0)],
            vec![TrimSelector::ArcLength(value)],
            TrimmingPreference::Parameter,
        )
        .is_err());
    }
}

#[test]
fn equal_arc_lengths_are_an_empty_trim() {
    assert!(open_trim(
        vec![TrimSelector::ArcLength(2.0)],
        vec![TrimSelector::ArcLength(2.0)],
        TrimmingPreference::Parameter,
    )
    .is_err());
}

#[test]
fn the_selector_keeps_its_value_exactly() {
    let value = 0.1 + 0.2;
    let selector = TrimSelector::ArcLength(value);
    assert_eq!(selector, TrimSelector::ArcLength(0.300_000_000_000_000_04));
    assert_ne!(selector, TrimSelector::Parameter(value));
}
