//! The curve a station is measured along, read from the graph (#241,
//! #264, #285).
//!
//! An atomic 2D or 3D curve is read as it is, unbounded on a line. A curve
//! placed at a station is its source carried by the placement. A curve
//! relation -- a composite, a trim, a surface curve whose 3D curve
//! governs, nested and placed in any combination -- is flattened into
//! `axiolid_reference::station::StationPiece`s, the spans of atomic curves
//! it runs along, and measured as one
//! `axiolid_reference::station::CompositeBasis`, whose module
//! documentation states the distance convention, the joint rule and the
//! exactness a composite carries. [`curve_path`] hands the same pieces out
//! as the owned, neutral [`CurvePath`] (#290), which a caller reads through
//! the curve-evaluation contract's `path_*` queries exactly as a station on
//! the relation is read here; the pieces resolved for a station borrow the
//! graph's curves, so resolving one copies none.
//!
//! The flattening reads a relation as a composite directrix does
//! (`crate::directrix`): a composite's segments in order, a segment whose
//! sense disagrees reversed; an untrimmed line inside a composite is its
//! parameter domain `[0, 1]`; a trim of an atomic curve between the
//! basis parameters its selectors name (a point selector inverted, an arc
//! length measured, a closed conic's trim possibly across its parameter
//! seam), reversed where its sense disagrees; a trim of a relation between
//! parameter selectors read as distances along it, which needs the
//! relation measured in arc length (its parameter); a curve placed at a
//! station carried by its placement, a trim of it in its source's
//! parameter.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Scalar, Tolerance, Transform3};
use axiolid_curve::{Curve2, Curve3, CurvePath, PathPiece};
use axiolid_model::{
    CurveRelation, GeometryGraph, GeometryNode, MasterRepresentation, NodeId, SeamSide,
    StationFrame, TrimSelector, TrimmingPreference,
};
use axiolid_reference::station::{
    exact_station_seams2, exact_station_seams3, station_frame_is_exact2, station_frame_is_exact3,
    station_seams2, station_seams3, station_section2_on, station_section3_on, CompositeBasis,
    DistanceConvention, SectionFrame, StationCurve, StationPiece, StationSeam, JOINT_TOLERANCE,
};

/// Most relations and placements a station basis may be nested through.
const MAX_BASIS_DEPTH: usize = 64;

/// The refusal of a basis no station can be measured along.
pub(crate) const UNSUPPORTED_BASIS: &str =
    "a station along an instanced curve, or along a curve relation other than a composite, a \
     trim, a surface curve whose 3D curve governs, or a curve placed at a station";

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: crate::BACKEND_ID,
        operation: Operation::CurveEvaluation,
        input,
    }
}

/// A graph curve flattened: one atomic curve whole (possibly placed), or
/// spans of atomic curves end to end.
enum Flat<'g> {
    Atomic {
        curve: StationCurve<'g>,
        placed: Option<Transform3>,
        /// Whether the placement, if any, is exact.
        exact: bool,
    },
    Pieces(Vec<StationPiece<'g>>),
}

impl<'g> Flat<'g> {
    /// As pieces: an atomic curve whole (a line its domain `[0, 1]`).
    fn into_pieces(self) -> GeomResult<Vec<StationPiece<'g>>> {
        match self {
            Self::Atomic {
                curve,
                placed,
                exact,
            } => {
                let piece = StationPiece::whole(curve)?;
                Ok(vec![match placed {
                    Some(rigid) => piece.placed(rigid, exact),
                    None => piece,
                }])
            }
            Self::Pieces(pieces) => Ok(pieces),
        }
    }
}

/// The curve `id` a distance runs along, flattened into the neutral
/// [`CurvePath`] (#290): the pieces a station on it is measured along, its
/// curves copied, to be read through the curve-evaluation contract's
/// `path_*` queries.
///
/// A relation (a composite, a trim, a surface curve whose 3D curve
/// governs, a curve placed at a station, nested in any combination) is
/// flattened as a station basis is (see the module documentation); an
/// atomic curve is one piece, whole, an untrimmed line its parameter domain
/// `[0, 1]` as inside a composite (a station along the atomic curve itself
/// reads a line unbounded, so query an atomic basis through the per-curve
/// methods). Placements at stations are resolved here, each exact only on
/// an exactly read frame. The path is not measured: the contract's
/// provider refuses one whose pieces do not meet or measure alike, by name.
///
/// # Errors
///
/// A node that is not a curve, an instanced curve or a relation no
/// distance runs along, and a trim or a placement that does not resolve,
/// each by name.
pub fn curve_path(graph: &GeometryGraph, id: NodeId) -> GeomResult<CurvePath> {
    Ok(flatten(graph, id, 0)?
        .into_pieces()?
        .into_iter()
        .map(PathPiece::from)
        .collect())
}

/// Pieces traversed the other way: in reverse order, each reversed.
fn reversed(pieces: Vec<StationPiece<'_>>) -> Vec<StationPiece<'_>> {
    pieces
        .into_iter()
        .rev()
        .map(StationPiece::reversed)
        .collect()
}

/// A station's basis: an atomic curve, possibly placed at stations (#264),
/// or a curve relation measured as a composite (#285).
#[derive(Debug, Clone)]
pub(crate) enum Basis<'g> {
    Atomic {
        curve: StationCurve<'g>,
        /// The composed placements, when the curve is placed.
        placed: Option<Transform3>,
        /// Whether a frame read on it is exact (a line, placed in exact
        /// frames).
        exact: bool,
    },
    Composite(CompositeBasis<'g>),
}

impl<'g> Basis<'g> {
    pub(crate) fn of(graph: &'g GeometryGraph, id: NodeId) -> GeomResult<Self> {
        Self::of_depth(graph, id, 0)
    }

    pub(super) fn of_depth(graph: &'g GeometryGraph, id: NodeId, depth: usize) -> GeomResult<Self> {
        Ok(match flatten(graph, id, depth)? {
            Flat::Atomic {
                curve,
                placed,
                exact,
            } => Basis::Atomic {
                curve,
                placed,
                exact: exact
                    && match curve {
                        StationCurve::Two(curve) => station_frame_is_exact2(curve),
                        StationCurve::Three(curve) => station_frame_is_exact3(curve),
                        _ => false,
                    },
            },
            Flat::Pieces(pieces) => Basis::Composite(CompositeBasis::new(pieces)?),
        })
    }

    /// The section frame at `distance`, in `frame`, read from `side` on a
    /// seam or a joint.
    pub(crate) fn section_on(
        &self,
        distance: Scalar,
        frame: StationFrame,
        side: SeamSide,
    ) -> GeomResult<SectionFrame> {
        let section = match self {
            Self::Atomic { curve, placed, .. } => {
                let section = match curve {
                    StationCurve::Two(curve) => station_section2_on(curve, distance, side)?,
                    StationCurve::Three(curve) => station_section3_on(curve, distance, side)?,
                    _ => {
                        return Err(unsupported(
                            "a station basis curve this compiler does not know",
                        ))
                    }
                };
                match placed {
                    Some(rigid) => section.carried(*rigid, curve.convention())?,
                    None => section,
                }
            }
            Self::Composite(composite) => composite.section_on(distance, side)?,
        };
        match frame {
            StationFrame::Section => Ok(section),
            StationFrame::Plan => section.plan(),
            _ => Err(unsupported("a station frame this compiler does not know")),
        }
    }

    /// Whether the frame at `distance`, read from `side`, is exact,
    /// rounding aside.
    pub(super) fn frame_is_exact(&self, distance: Scalar, side: SeamSide) -> bool {
        match self {
            Self::Atomic { exact, .. } => *exact,
            Self::Composite(composite) => composite.frame_is_exact_at(distance, side),
        }
    }

    /// Every seam (a B-spline's corner knots by quadrature); a placement
    /// moves none of them.
    pub(super) fn seams(&self) -> GeomResult<Vec<StationSeam>> {
        match self {
            Self::Atomic { curve, .. } => match curve {
                StationCurve::Two(curve) => station_seams2(curve),
                StationCurve::Three(curve) => station_seams3(curve),
                _ => Err(unsupported(
                    "a station basis curve this compiler does not know",
                )),
            },
            Self::Composite(composite) => composite.seams(),
        }
    }

    /// Every seam, each exact or refused.
    pub(super) fn exact_seams(&self) -> GeomResult<Vec<StationSeam>> {
        match self {
            Self::Atomic { curve, .. } => match curve {
                StationCurve::Two(curve) => exact_station_seams2(curve),
                StationCurve::Three(curve) => exact_station_seams3(curve),
                _ => Err(unsupported(
                    "a station basis curve this compiler does not know",
                )),
            },
            Self::Composite(composite) => composite.exact_seams(),
        }
    }
}

fn flatten<'g>(graph: &'g GeometryGraph, id: NodeId, depth: usize) -> GeomResult<Flat<'g>> {
    if depth > MAX_BASIS_DEPTH {
        return Err(GeomError::BudgetExceeded {
            resource: "relations and placements nested in a station basis",
        });
    }
    match graph.get(id) {
        Some(GeometryNode::Curve2(curve)) => Ok(Flat::Atomic {
            curve: StationCurve::Two(curve),
            placed: None,
            exact: true,
        }),
        Some(GeometryNode::Curve3(curve)) => Ok(Flat::Atomic {
            curve: StationCurve::Three(curve),
            placed: None,
            exact: true,
        }),
        // A curve placed at a station (#264): its source carried by the
        // placement, whatever the source is.
        Some(GeometryNode::InstanceAtStation(placed)) => {
            let source = flatten(graph, placed.source, depth + 1)?;
            let placement = super::placement_at(graph, &placed.station, depth + 1)?;
            Ok(match source {
                Flat::Atomic {
                    curve,
                    placed: inner,
                    exact,
                } => Flat::Atomic {
                    curve,
                    placed: Some(match inner {
                        Some(inner) => placement.transform * inner,
                        None => placement.transform,
                    }),
                    exact: exact && placement.exact,
                },
                Flat::Pieces(pieces) => Flat::Pieces(
                    pieces
                        .into_iter()
                        .map(|piece| piece.placed(placement.transform, placement.exact))
                        .collect(),
                ),
            })
        }
        // A composite (#285): its segments end to end, a segment whose
        // sense disagrees traversed backwards.
        Some(GeometryNode::CurveRelation(CurveRelation::Composite { segments })) => {
            let mut out = Vec::new();
            for segment in segments {
                let pieces = flatten(graph, segment.curve, depth + 1)?.into_pieces()?;
                out.extend(match segment.same_sense {
                    true => pieces,
                    false => reversed(pieces),
                });
            }
            Ok(Flat::Pieces(out))
        }
        Some(GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis,
            start,
            end,
            sense_agreement,
            preference,
        })) => {
            let pieces = match flatten(graph, *basis, depth + 1)? {
                Flat::Atomic {
                    curve,
                    placed,
                    exact,
                } => {
                    let (lo, hi) =
                        trim_parameters(curve, start, end, *sense_agreement, *preference)?;
                    let piece = StationPiece::between_parameters(curve, lo, hi)?;
                    vec![match placed {
                        Some(rigid) => piece.placed(rigid, exact),
                        None => piece,
                    }]
                }
                Flat::Pieces(pieces) => {
                    let composite = CompositeBasis::new(pieces)?;
                    if composite.convention() != DistanceConvention::ArcLength3d {
                        return Err(unsupported(
                            "a trim of a plan-measured curve relation as a station basis: its \
                             parameter is its arc length, not the plan distance it is measured \
                             in",
                        ));
                    }
                    let a = crate::directrix::parameter_only(start, *preference, "start")?;
                    let b = crate::directrix::parameter_only(end, *preference, "end")?;
                    composite.pieces_between(a.min(b), a.max(b))?
                }
            };
            Ok(Flat::Pieces(if *sense_agreement {
                pieces
            } else {
                reversed(pieces)
            }))
        }
        Some(GeometryNode::CurveRelation(CurveRelation::SurfaceCurve {
            curve_3d,
            master: MasterRepresentation::Curve3d,
            ..
        })) => flatten(graph, *curve_3d, depth + 1),
        Some(GeometryNode::Instance(_) | GeometryNode::CurveRelation(_)) => {
            Err(unsupported(UNSUPPORTED_BASIS))
        }
        Some(_) => Err(GeomError::InvalidInput(format!(
            "station basis {id:?} is not a curve"
        ))),
        None => Err(GeomError::InvalidInput(format!(
            "station basis {id:?} is outside the graph"
        ))),
    }
}

/// The tolerance a point selector is inverted to: the joint tolerance at
/// the scale of the points named.
fn point_tolerance(selectors: &[&[TrimSelector]]) -> Tolerance {
    let scale = selectors
        .iter()
        .flat_map(|list| list.iter())
        .map(|selector| match selector {
            TrimSelector::Point2(p) => p.abs().max_element(),
            TrimSelector::Point3(p) => p.abs().max_element(),
            _ => 0.0,
        })
        .fold(1.0, Scalar::max);
    Tolerance::new(JOINT_TOLERANCE * scale, Tolerance::METRE.angular()).unwrap_or(Tolerance::METRE)
}

/// The basis parameters `lo < hi` a trim of an atomic curve spans, as a
/// composite directrix reads them: a closed conic from the start
/// selector in the trim's sense (possibly across its parameter seam), any
/// other curve between its two selectors.
fn trim_parameters(
    curve: StationCurve<'_>,
    start: &[TrimSelector],
    end: &[TrimSelector],
    sense: bool,
    preference: TrimmingPreference,
) -> GeomResult<(Scalar, Scalar)> {
    let tolerance = point_tolerance(&[start, end]);
    let (a, b) = match curve {
        StationCurve::Three(curve) => (
            crate::directrix::parameter(start, preference, "start", curve, sense, tolerance)?,
            crate::directrix::parameter(end, preference, "end", curve, sense, tolerance)?,
        ),
        StationCurve::Two(curve) => (
            parameter2(start, preference, "start", curve, sense, tolerance)?,
            parameter2(end, preference, "end", curve, sense, tolerance)?,
        ),
        _ => {
            return Err(unsupported(
                "a station basis curve this compiler does not know",
            ))
        }
    };
    if a == b {
        return Err(GeomError::Degenerate(
            "a trimmed station basis has an empty interval".into(),
        ));
    }
    let closed = matches!(
        curve,
        StationCurve::Two(Curve2::Circle(_) | Curve2::Ellipse(_))
            | StationCurve::Three(Curve3::Circle(_) | Curve3::Ellipse(_))
    );
    if closed {
        return crate::directrix::periodic_trim_interval(
            core::f64::consts::TAU,
            a,
            b,
            sense,
            None,
            &ExecutionOptions::new(tolerance),
        );
    }
    Ok((a.min(b), a.max(b)))
}

/// A trim selector of a 2D curve read as its parameter, by the rules a 3D
/// directrix's trims are read by (`crate::directrix`): a parameter as is,
/// an arc length measured from parameter `0` in the trim's sense, a point
/// inverted exactly or refused.
fn parameter2(
    selectors: &[TrimSelector],
    preference: TrimmingPreference,
    label: &str,
    curve: &Curve2,
    sense: bool,
    tolerance: Tolerance,
) -> GeomResult<Scalar> {
    let measured = |selector: &TrimSelector| -> GeomResult<Option<Scalar>> {
        match selector {
            TrimSelector::Parameter(value) => Ok(Some(*value)),
            TrimSelector::ArcLength(length) => {
                if !length.is_finite() {
                    return Err(GeomError::InvalidInput(format!(
                        "trimmed station basis {label} arc length must be finite"
                    )));
                }
                let signed = if sense { *length } else { -*length };
                axiolid_reference::arc_parameter::parameter_at_arc_length2(curve, 0.0, signed)
                    .map(Some)
                    .map_err(|error| {
                        GeomError::InvalidInput(format!(
                            "trimmed station basis {label} arc length {length} does not resolve \
                             on its basis: {error}"
                        ))
                    })
            }
            _ => Ok(None),
        }
    };
    let is_measure = |selector: &&TrimSelector| {
        matches!(
            selector,
            TrimSelector::Parameter(_) | TrimSelector::ArcLength(_)
        )
    };
    let first_measure = || -> GeomResult<Option<Scalar>> {
        selectors.iter().find(is_measure).map_or(Ok(None), measured)
    };
    let invert = || {
        selectors.iter().find_map(|selector| match selector {
            TrimSelector::Point2(point) => {
                axiolid_reference::curve::invert2(curve, *point, tolerance).ok()
            }
            _ => None,
        })
    };
    let selected = match preference {
        TrimmingPreference::Parameter => first_measure()?,
        TrimmingPreference::Unspecified => match selectors.first() {
            Some(selector) if is_measure(&selector) => measured(selector)?,
            _ => invert(),
        },
        TrimmingPreference::Cartesian => match invert() {
            Some(value) => Some(value),
            None => first_measure()?,
        },
    };
    selected.filter(|value| value.is_finite()).ok_or_else(|| {
        GeomError::InvalidInput(format!(
            "trimmed station basis {label} needs a parameter selector, or a point selector \
             that lies on the basis curve"
        ))
    })
}
