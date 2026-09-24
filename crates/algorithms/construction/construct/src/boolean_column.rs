//! Coaxial booleans and plane cuts whose result is not one prism (#120).
//!
//! Both reduce to a column description ([`crate::column`]): the plan is cut
//! by every operand ring (and, for a plane cut, by the lines where the plane
//! meets the two cap heights), and above each cell the solid is a stack of
//! height intervals that depends only on which rings contain the cell.
//!
//! - A stepped boolean: two prisms spanning different heights. Over a cell
//!   the solid is the boolean of the two operands' height intervals there.
//! - A plane crossing a cap: over a cell the kept part of the prism runs
//!   from a cap or the plane to a cap or the plane, whichever is nearer.

use axiolid_brep::{ExactBRep, FaceName, Operand, SweptFace};
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{BooleanOperator, Point2, Scalar, Tolerance, Vec2};
use axiolid_overlay::{ArcArrangement, ArcRing, ArcVertex, EdgeSource};

use crate::column::{build_columns, Block, Columns};
use crate::extrude_arc::Level;

/// One operand of a coaxial boolean, as rings in the arrangement.
pub(crate) struct ColumnOperand {
    /// Outer ring first, then holes; any winding.
    pub(crate) rings: Vec<ArcRing>,
    pub(crate) bottom: Scalar,
    pub(crate) top: Scalar,
}

/// Where an operand's rings and heights landed.
struct Placed {
    first: usize,
    count: usize,
    bottom: usize,
    top: usize,
}

impl Placed {
    /// Inside the outer ring and outside every hole.
    fn contains(&self, mask: &[bool]) -> bool {
        mask[self.first] && !(self.first + 1..self.first + self.count).any(|ring| mask[ring])
    }

    fn spans(&self, gap: usize) -> bool {
        self.bottom <= gap && gap < self.top
    }
}

/// Distinct heights, merged within the tolerance, lowest first.
fn distinct_heights(values: &[Scalar], tolerance: Tolerance) -> Vec<Scalar> {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let mut out: Vec<Scalar> = Vec::with_capacity(sorted.len());
    for value in sorted {
        match out.last() {
            Some(&last) if value - last <= tolerance.linear() => {}
            _ => out.push(value),
        }
    }
    out
}

fn height_index(heights: &[Scalar], value: Scalar, tolerance: Tolerance) -> usize {
    heights
        .iter()
        .position(|&h| (h - value).abs() <= tolerance.linear())
        .unwrap_or(0)
}

/// Merge the solid gaps `0..gaps` into blocks of consecutive heights.
fn blocks(gaps: usize, solid: impl Fn(usize) -> bool) -> Vec<Block> {
    let mut out = Vec::new();
    let mut gap = 0;
    while gap < gaps {
        if !solid(gap) {
            gap += 1;
            continue;
        }
        let start = gap;
        while gap < gaps && solid(gap) {
            gap += 1;
        }
        out.push((start, gap));
    }
    out
}

/// The wall name of a piece: the first operand ring that carries it.
///
/// Wall ordinals count edges across an operand's rings in order, the same
/// numbering the operand's own extrusion uses, so a fragment names the
/// exact input wall it came from.
fn wall_name(sources: &[EdgeSource], placed: &[(Operand, &Placed, &[usize])]) -> Option<FaceName> {
    let source = sources.iter().min_by_key(|s| s.ring)?;
    for (operand, place, counts) in placed {
        if source.ring >= place.first && source.ring < place.first + place.count {
            let before: usize = counts[..source.ring - place.first].iter().sum();
            let ordinal = u32::try_from(before + source.edge).ok()?;
            return Some(FaceName::swept(SweptFace::Side(ordinal)).fragment(*operand));
        }
    }
    None
}

/// Exact coaxial boolean of two operands, any spans, one B-rep per piece.
///
/// Solids are ordered by their lowest vertex, `x` then `y`.
pub(crate) fn coaxial_columns(
    subject: &ColumnOperand,
    tool: &ColumnOperand,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> GeomResult<Vec<ExactBRep>> {
    let heights = distinct_heights(
        &[subject.bottom, subject.top, tool.bottom, tool.top],
        tolerance,
    );
    let planes: Vec<Level> = heights.iter().map(|&h| Level::flat(h)).collect();
    let mut rings = Vec::with_capacity(subject.rings.len() + tool.rings.len());
    rings.extend(subject.rings.iter().cloned());
    rings.extend(tool.rings.iter().cloned());
    let arrangement = ArcArrangement::new(&rings, tolerance)
        .map_err(|error| GeomError::InvalidInput(format!("coaxial boolean section: {error:?}")))?;
    let at = |value| height_index(&heights, value, tolerance);
    let s = Placed {
        first: 0,
        count: subject.rings.len(),
        bottom: at(subject.bottom),
        top: at(subject.top),
    };
    let t = Placed {
        first: subject.rings.len(),
        count: tool.rings.len(),
        bottom: at(tool.bottom),
        top: at(tool.top),
    };
    let solid_in = |mask: &[bool], gap: usize| {
        let a = s.contains(mask) && s.spans(gap);
        let b = t.contains(mask) && t.spans(gap);
        match operator {
            BooleanOperator::Union => a || b,
            BooleanOperator::Intersection => a && b,
            _ => a && !b,
        }
    };
    if !matches!(
        operator,
        BooleanOperator::Union | BooleanOperator::Intersection | BooleanOperator::Difference
    ) {
        return Err(crate::boolean_exact::unsupported(
            "unknown exact prism boolean operator",
        ));
    }
    let gaps = heights.len() - 1;
    let stack = |mask: &[bool]| blocks(gaps, |gap| solid_in(mask, gap));

    // A cap is named after the operand whose own cap lies on that plane
    // facing the same way, subject first; failing that, after the tool's
    // opposite cap (a ledge or pocket floor left by a difference).
    //
    // The subject's opposite cap can never bound the result: facing up on
    // the subject's bottom plane (or down on its top) needs solid on the
    // side where the subject has none, which union, intersection and
    // difference only get from the tool -- whose own cap then ends there
    // and is matched first. So that case is not a branch here.
    let cap_name = |plane: usize, up: bool| {
        let own = |p: &Placed| {
            if up {
                p.top == plane
            } else {
                p.bottom == plane
            }
        };
        let other = |p: &Placed| {
            if up {
                p.bottom == plane
            } else {
                p.top == plane
            }
        };
        let face = |at_top: bool| {
            FaceName::swept(if at_top {
                SweptFace::EndCap
            } else {
                SweptFace::StartCap
            })
        };
        if own(&s) {
            Some(face(up).fragment(Operand::Subject))
        } else if own(&t) {
            Some(face(up).fragment(Operand::Tool))
        } else if other(&t) {
            Some(face(!up).fragment(Operand::Tool))
        } else {
            None
        }
    };
    let counts_s: Vec<usize> = subject.rings.iter().map(|r| r.vertices.len()).collect();
    let counts_t: Vec<usize> = tool.rings.iter().map(|r| r.vertices.len()).collect();
    let placed = [
        (Operand::Subject, &s, counts_s.as_slice()),
        (Operand::Tool, &t, counts_t.as_slice()),
    ];
    let name_wall = |sources: &[EdgeSource], _side: &[bool]| wall_name(sources, &placed);
    let solids = build_columns(&Columns {
        arrangement: &arrangement,
        planes: &planes,
        stack: &stack,
        cap_name: &cap_name,
        wall_name: &name_wall,
        tolerance,
    })?;
    ordered(solids, tolerance)
}

/// A prism between two flat heights cut by `level`, keeping one side,
/// when the plane may cross either cap inside the section.
///
/// The plan is additionally cut where the plane reaches the top height and
/// where it reaches the bottom height; between those lines the plane is the
/// kept part's cap, beyond them the prism's own cap is.
pub(crate) fn clip_columns(
    section: &ArcRing,
    (bottom, top): (Scalar, Scalar),
    level: Level,
    keeps_above: bool,
    tolerance: Tolerance,
) -> GeomResult<Vec<ExactBRep>> {
    let mut rings = vec![section.clone()];
    // Half-planes where the plane is at or above the top, and at or below
    // the bottom, as boxes around the section clipped by a line. Each is
    // present only if the line actually crosses the box.
    let bounds = plan_box(section);
    let above_top = half_plane_ring(&bounds, level, top, true);
    let below_bottom = half_plane_ring(&bounds, level, bottom, false);
    let above_index = above_top.map(|ring| {
        rings.push(ring);
        rings.len() - 1
    });
    let below_index = below_bottom.map(|ring| {
        rings.push(ring);
        rings.len() - 1
    });
    let arrangement = ArcArrangement::new(&rings, tolerance)
        .map_err(|error| GeomError::InvalidInput(format!("arc prism clip section: {error:?}")))?;
    let planes = [Level::flat(bottom), Level::flat(top), level];
    let (flat_bottom, flat_top, cut) = (0, 1, 2);
    let stack = |mask: &[bool]| -> Vec<Block> {
        if !mask[0] {
            return Vec::new();
        }
        let over_top = above_index.is_some_and(|i| mask[i]);
        let under_bottom = below_index.is_some_and(|i| mask[i]);
        if keeps_above {
            if over_top {
                Vec::new()
            } else {
                vec![(if under_bottom { flat_bottom } else { cut }, flat_top)]
            }
        } else if under_bottom {
            Vec::new()
        } else {
            vec![(flat_bottom, if over_top { flat_top } else { cut })]
        }
    };
    let cap_name = |plane: usize, _up: bool| match plane {
        0 => Some(FaceName::swept(SweptFace::StartCap)),
        1 => Some(FaceName::swept(SweptFace::EndCap)),
        _ => None,
    };
    let name_wall = |sources: &[EdgeSource], _side: &[bool]| {
        let source = sources.iter().find(|s| s.ring == 0)?;
        Some(FaceName::swept(SweptFace::Side(
            u32::try_from(source.edge).ok()?,
        )))
    };
    let solids = build_columns(&Columns {
        arrangement: &arrangement,
        planes: &planes,
        stack: &stack,
        cap_name: &cap_name,
        wall_name: &name_wall,
        tolerance,
    })?;
    ordered(solids, tolerance)
}

/// Axis-aligned box around a ring, padded so a clip line through the box
/// never runs along the section's own boundary.
fn plan_box(ring: &ArcRing) -> (Point2, Point2) {
    let mut lo = Point2::new(Scalar::INFINITY, Scalar::INFINITY);
    let mut hi = Point2::new(Scalar::NEG_INFINITY, Scalar::NEG_INFINITY);
    let count = ring.vertices.len();
    for index in 0..count {
        let vertex = ring.vertices[index];
        let next = ring.vertices[(index + 1) % count].point;
        // An arc stays within its sagitta of the chord's box.
        let chord = (next - vertex.point).length();
        let pad = vertex.bulge.abs() * chord / 2.0;
        for p in [vertex.point, next] {
            lo = lo.min(p - Vec2::splat(pad));
            hi = hi.max(p + Vec2::splat(pad));
        }
    }
    let margin = (hi - lo).max_element().max(1.0);
    (lo - Vec2::splat(margin), hi + Vec2::splat(margin))
}

/// The part of the box where `level` is at or above `height` (`above`) or
/// at or below it, as a counter-clockwise ring; `None` if that part is
/// empty or the whole box.
fn half_plane_ring(
    (lo, hi): &(Point2, Point2),
    level: Level,
    height: Scalar,
    above: bool,
) -> Option<ArcRing> {
    let side = |p: Point2| {
        let d = level.at(p) - height;
        if above {
            d
        } else {
            -d
        }
    };
    let corners = [
        Point2::new(lo.x, lo.y),
        Point2::new(hi.x, lo.y),
        Point2::new(hi.x, hi.y),
        Point2::new(lo.x, hi.y),
    ];
    let values = corners.map(side);
    if values.iter().all(|&v| v >= 0.0) || values.iter().all(|&v| v <= 0.0) {
        return None;
    }
    // Sutherland-Hodgman against one line: keep corners on the wanted side,
    // add the crossing point on every edge that changes side.
    let mut points = Vec::with_capacity(5);
    for index in 0..4 {
        let (p, q) = (corners[index], corners[(index + 1) % 4]);
        let (a, b) = (values[index], values[(index + 1) % 4]);
        if a >= 0.0 {
            points.push(p);
        }
        if (a > 0.0 && b < 0.0) || (a < 0.0 && b > 0.0) {
            let t = a / (a - b);
            points.push(p + (q - p) * t);
        }
    }
    points.dedup();
    (points.len() >= 3).then(|| ArcRing::new(points.into_iter().map(ArcVertex::straight).collect()))
}

/// Order solids by their lowest vertex, `x` then `y`, as the other
/// multi-solid entry points do; gate each through the geometric audit.
fn ordered(solids: Vec<ExactBRep>, tolerance: Tolerance) -> GeomResult<Vec<ExactBRep>> {
    let mut keyed = solids
        .into_iter()
        .map(|solid| {
            let key = solid
                .topology()
                .vertices()
                .iter()
                .map(|v| (v.position.x, v.position.y))
                .min_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)))
                .unwrap_or((0.0, 0.0));
            (key, solid)
        })
        .collect::<Vec<_>>();
    keyed.sort_by(|a, b| a.0 .0.total_cmp(&b.0 .0).then(a.0 .1.total_cmp(&b.0 .1)));
    keyed
        .into_iter()
        .map(|(_, solid)| crate::boolean_exact::gate_geometry(solid, tolerance))
        .collect()
}

/// Whether a coaxial boolean of these spans is stepped: not one prism.
pub(crate) fn is_stepped(
    subject: (Scalar, Scalar),
    tool: (Scalar, Scalar),
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> bool {
    match operator {
        BooleanOperator::Union => {
            !tolerance.eq(subject.0, tool.0) || !tolerance.eq(subject.1, tool.1)
        }
        BooleanOperator::Difference => {
            // Only a tool covering the subject's whole span leaves one
            // prism. A tool reaching in from one end, stopping inside at
            // both ends, or missing the span altogether does not.
            !(tool.0 <= subject.0 + tolerance.linear() && tool.1 >= subject.1 - tolerance.linear())
        }
        _ => false,
    }
}
