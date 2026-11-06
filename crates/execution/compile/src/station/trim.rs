//! The parameter a trim of a curve relation is read in (#424, ADR 0082).
//!
//! A trim reads its selectors in the parameter of the curve it trims, and
//! every relation has one:
//!
//! - a trim has its basis's parameter (ISO 10303-42), an offset its basis's,
//!   a placed curve its source's;
//! - an atomic curve has its own.
//!
//! This module reads the span `[lo, hi]` of that parameter a trim keeps (the
//! trim's own span), and narrows it to a window a trim above it keeps in
//! turn. A closed conic is a periodic carrier: a trim of it runs from its
//! start selector in its sense, across the seam if it must, and a window
//! on a trim of it is read modulo the period into the trim's arc, as a sweep
//! range is on a directrix (`crate::directrix`).
//!
//! The polyline compiler (`crate::polyline`) and the station basis
//! (`super::basis`) read a trim of a trim through the same span, so a
//! sampled directrix, a polyline and a station on it agree.

use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::Scalar;
use axiolid_curve::{Curve2, Curve3};
use axiolid_model::{
    CurveRelation, GeometryGraph, GeometryNode, NodeId, TrimSelector, TrimmingPreference,
};
use axiolid_reference::curve::{domain2, domain3};
use axiolid_reference::station::StationCurve;

/// Most relations and placements a trim's carrier is read through.
const MAX_DEPTH: usize = 256;

/// Why a point or arc-length selector, or a Cartesian preference, on a trim
/// of a relation is refused. TODO(#432): read them on a relation.
const TRIM_SELECTOR_ON_RELATION: &str = "a point or arc-length trim selector, or a Cartesian \
     preference, on a trim of a curve relation or of a placed curve: only a parameter selector \
     is read in a relation's parameter";

/// Why a trim whose parameters leave the trimmed curve's is refused.
const TRIM_OUTSIDE: &str = "a trim's parameters lie outside the trimmed curve it trims";

/// Why a trim of a composite that runs against its sense is refused.
/// TODO(#433): decide whether it crosses the closing joint.
const TRIM_ACROSS_CLOSING_JOINT: &str = "a trim of a closed composite curve whose selectors run \
     against its sense: whether it crosses the closing joint is not read";

/// The curve a trim cuts, as far as its parameter goes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Carrier {
    /// A closed conic: its parameter runs round with this period.
    Periodic(Scalar),
    /// A curve over `[lo, hi]` of its parameter; one that repeats with a
    /// `period` is a trim of a closed conic, and a `closed` one meets itself
    /// where it ends.
    Interval {
        lo: Scalar,
        hi: Scalar,
        period: Option<Scalar>,
        closed: bool,
    },
    /// An unbounded line.
    Unbounded,
}

impl Carrier {
    /// The period the parameter repeats with, if it does.
    fn period(self) -> Option<Scalar> {
        match self {
            Self::Periodic(period) => Some(period),
            Self::Interval { period, .. } => period,
            Self::Unbounded => None,
        }
    }
}

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: crate::BACKEND_ID,
        operation: Operation::CurveEvaluation,
        input,
    }
}

/// The carrier of an atomic 2D curve.
fn carrier2(curve: &Curve2) -> Carrier {
    match curve {
        Curve2::Circle(_) | Curve2::Ellipse(_) => Carrier::Periodic(core::f64::consts::TAU),
        Curve2::Line(_) => Carrier::Unbounded,
        _ => {
            let domain = domain2(curve);
            interval(domain.start, domain.end)
        }
    }
}

/// The carrier of an atomic 3D curve.
fn carrier3(curve: &Curve3) -> Carrier {
    match curve {
        Curve3::Circle(_) | Curve3::Ellipse(_) => Carrier::Periodic(core::f64::consts::TAU),
        Curve3::Line(_) => Carrier::Unbounded,
        _ => {
            let domain = domain3(curve);
            interval(domain.start, domain.end)
        }
    }
}

fn interval(a: Scalar, b: Scalar) -> Carrier {
    Carrier::Interval {
        lo: a.min(b),
        hi: a.max(b),
        period: None,
        closed: false,
    }
}

/// The carrier of the curve `id`, through placements, profiles, offsets and
/// trims.
///
/// # Errors
///
/// A node with no parameter a trim could be read in, by name; the nesting
/// is bounded.
pub(crate) fn carrier(
    graph: &GeometryGraph,
    id: NodeId,
    slack: Scalar,
    depth: usize,
) -> GeomResult<Carrier> {
    if depth > MAX_DEPTH {
        return Err(GeomError::BudgetExceeded {
            resource: "relations and placements nested in a trim",
        });
    }
    match graph.get(id) {
        Some(GeometryNode::Curve2(curve)) => Ok(carrier2(curve)),
        Some(GeometryNode::Curve3(curve)) => Ok(carrier3(curve)),
        Some(GeometryNode::Instance(instance)) => carrier(graph, instance.source, slack, depth + 1),
        Some(GeometryNode::InstanceAtStation(at)) => carrier(graph, at.source, slack, depth + 1),
        Some(GeometryNode::OpenProfile(profile)) => carrier(graph, profile.path, slack, depth + 1),
        Some(GeometryNode::CurveRelation(CurveRelation::Offset { basis, .. })) => {
            carrier(graph, *basis, slack, depth + 1)
        }
        Some(GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis,
            start,
            end,
            sense_agreement,
            preference,
        })) => {
            let span = trimmed_span(
                graph,
                *basis,
                (start, end),
                (*sense_agreement, *preference),
                None,
                slack,
                depth + 1,
            )?;
            let period = carrier(graph, *basis, slack, depth + 1)?.period();
            Ok(Carrier::Interval {
                lo: span.0,
                hi: span.1,
                period,
                closed: false,
            })
        }
        Some(GeometryNode::CurveRelation(CurveRelation::Composite { .. })) => Err(unsupported(
            "a trim of a curve relation: a composite is not trimmed to a polyline yet",
        )),
        _ => Err(unsupported(
            "a trim of a node with no curve parameter to read it in: a collection, a surface \
             curve, an offset by stations or a node that is not a curve",
        )),
    }
}

/// The span `[lo, hi]` of its basis's parameter a trim of `basis` keeps,
/// narrowed to the `incoming` window of a trim above it, when one is given.
///
/// The basis's own curve reads its selectors as a station basis does
/// ([`super::trim_parameters`]: a parameter, an arc length or a point); any
/// other basis reads a parameter selector in its carrier. `slack` is the
/// tolerance a parameter may lie beyond an end by.
///
/// # Errors
///
/// A selector other than a parameter on a relation, a span outside the
/// carrier's or empty, and a trim of a closed composite running against its
/// sense, by name.
pub(crate) fn trimmed_span(
    graph: &GeometryGraph,
    basis: NodeId,
    (start, end): (&[TrimSelector], &[TrimSelector]),
    (sense, preference): (bool, TrimmingPreference),
    incoming: Option<(Scalar, Scalar)>,
    slack: Scalar,
    depth: usize,
) -> GeomResult<(Scalar, Scalar)> {
    let carrier = carrier(graph, basis, slack, depth)?;
    let own = match graph.get(basis) {
        Some(GeometryNode::Curve2(curve)) => {
            super::trim_parameters(StationCurve::Two(curve), start, end, sense, preference)?
        }
        Some(GeometryNode::Curve3(curve)) => {
            super::trim_parameters(StationCurve::Three(curve), start, end, sense, preference)?
        }
        _ => {
            let a = relation_parameter(start, preference)?;
            let b = relation_parameter(end, preference)?;
            cut(carrier, (a, b), sense, slack)?
        }
    };
    match incoming {
        None => Ok(own),
        Some(window) => narrowed(own, carrier.period(), window, slack),
    }
}

/// A trim selector list read as a parameter of a relation: its first
/// parameter-kind selector, finite.
///
/// # Errors
///
/// An arc-length selector, a list with no parameter selector (a point only)
/// and a Cartesian preference, by name; a parameter that is not finite.
fn relation_parameter(
    selectors: &[TrimSelector],
    preference: TrimmingPreference,
) -> GeomResult<Scalar> {
    if preference == TrimmingPreference::Cartesian {
        return Err(unsupported(TRIM_SELECTOR_ON_RELATION));
    }
    match selectors
        .iter()
        .find(|selector| crate::directrix::is_parameter_kind(selector))
    {
        Some(TrimSelector::Parameter(value)) if value.is_finite() => Ok(*value),
        Some(TrimSelector::Parameter(_)) => Err(GeomError::InvalidInput(
            "a trim's parameter is not finite".into(),
        )),
        _ => Err(unsupported(TRIM_SELECTOR_ON_RELATION)),
    }
}

/// The span `[lo, hi]` a trim between parameters `a` and `b` of a relation
/// keeps, in the sense it runs by.
fn cut(
    carrier: Carrier,
    (a, b): (Scalar, Scalar),
    sense: bool,
    slack: Scalar,
) -> GeomResult<(Scalar, Scalar)> {
    if a == b {
        return Err(GeomError::Degenerate(
            "a trimmed curve has an empty parameter interval".into(),
        ));
    }
    match carrier {
        Carrier::Periodic(period) => Ok(crate::directrix::periodic_arc(period, a, b, sense)),
        Carrier::Unbounded => Ok((a.min(b), a.max(b))),
        Carrier::Interval {
            lo,
            hi,
            period,
            closed,
        } => {
            if closed && (b > a) != sense {
                return Err(unsupported(TRIM_ACROSS_CLOSING_JOINT));
            }
            narrowed((lo, hi), period, (a, b), slack)
        }
    }
}

/// The part of the span `(lo, hi)` the `window` names: its ends read modulo
/// `period` into the span where the carrier repeats, within `slack` of the
/// span's ends otherwise.
///
/// # Errors
///
/// A window outside the span, and an empty one, by name.
pub(crate) fn narrowed(
    (lo, hi): (Scalar, Scalar),
    period: Option<Scalar>,
    (a, b): (Scalar, Scalar),
    slack: Scalar,
) -> GeomResult<(Scalar, Scalar)> {
    if let Some(period) = period {
        return crate::directrix::range_on_periodic_arc((lo, hi), period, (a, b), slack).map_err(
            |error| match error {
                GeomError::InvalidInput(_) => GeomError::InvalidInput(TRIM_OUTSIDE.into()),
                other => other,
            },
        );
    }
    let (x, y) = (a.min(b), a.max(b));
    if x < lo - slack || y > hi + slack {
        return Err(GeomError::InvalidInput(TRIM_OUTSIDE.into()));
    }
    let (x, y) = (x.max(lo), y.min(hi));
    if x >= y {
        return Err(GeomError::Degenerate(
            "a trimmed curve has an empty parameter interval".into(),
        ));
    }
    Ok((x, y))
}
