//! Edge samples inverted onto a surface relation (#454, ADR 0110).
//!
//! An edge without a pcurve reaches a face over a relation as 3D samples,
//! and the tessellator needs each as `(u, v)`. Every supported relation is
//! a one-parameter family of generators along a swept parameter `t`: the
//! circle about the axis through `c(v)` (revolution, `t = v`), the line
//! `c(u) + s d` (linear extrusion, `t = u`) and the ruling from `A(u)` to
//! `B(u)` (ruled, `t = u`). A sample `P` is inverted in three steps.
//!
//! 1. **Locate `t`.** `g(t)`, the distance from `P` to the generator at
//!    `t`, changes no faster than the swept curve moves, so over a cell
//!    `[a, b]` with midpoint `m` no point is nearer than `g(m) - L (b - a) /
//!    2`, `L` the curve's certified speed bound over the cell
//!    ([`curve_derivative_bounds3`]; the larger of the sides' for a ruled
//!    surface; `|e|` for a line of direction `e`). A branch-and-bound over
//!    the bracket of `t` drops a cell whose lower bound exceeds `min(tau,
//!    best + slack)`, `best` the least `g` evaluated, `slack =`
//!    [`MINIMISER_SLACK`] `tau` and `tau` the linear tolerance, and keeps a
//!    cell as a leaf once `L (b - a) / 2 <= slack`. The leaves cover every
//!    `t` with `g(t) <= min g + slack` and lie within `min g + 3 slack`.
//! 2. **Require one minimiser.** The leaves must merge into one interval,
//!    or on a periodic bracket into two touching both of its ends, which
//!    meet across the seam. Otherwise two generators pass within the slack
//!    of the nearest one: the sample has no unique parameters (a closed
//!    swept curve's seam, a surface covered twice or one meeting itself),
//!    and it is refused, never resolved by a guess.
//! 3. **Read and check.** Gauss-Newton on `t` sharpens it inside its
//!    interval, the other parameter is read in closed form, and the
//!    residual `|S(u, v) - P|`, measured by [`RelationSurface::evaluate`],
//!    must be within `tau`.
//!
//! | Relation | `g(t)` | other parameter |
//! |---|---|---|
//! | Revolution | `hypot(rho_q - rho_p, h_q - h_p)`, `q = c(t) - O`, `p = P - O`, `h` along the axis and `rho` from it | `u = atan2((a x q_perp).p_perp, q_perp.p_perp)` |
//! | Extrusion | `|w - (w.d) d / |d|^2|`, `w = c(t) - P` | `v = (P - c(t)).d / |d|^2` |
//! | Ruled | distance from `P` to the segment `A(t)`-`B(t)` | `v = (P - A).(B - A) / |B - A|^2` clamped to `[0, 1]` |
//!
//! The bracket of `t` is a swept B-spline's domain, `[0, 2 pi]` for a
//! circle or an ellipse (periodic) and `[0, 1]` for a ruled surface, each
//! cut into [`COARSE_CELLS`] cells up front with their speed bounds and
//! midpoint generators. A swept line is unbounded, so it is handled
//! exactly instead: extruded, it sweeps a plane, inverted in closed form;
//! revolved, its bracket is cut per sample to where both `|h_q - h_p|` and
//! `rho_q - rho_p` can be at most `tau`.
//!
//! Refused by name, before any search: a sample within `tau` of the axis
//! of revolution or of a ruled surface's apex, where every `u` names the
//! point, and after it, a sample on a ruling shorter than `tau`, where
//! every `v` does (#463 owns those poles). A sample that leaves the surface
//! is refused with a certified lower bound on its distance. A cover that
//! needs more than [`RELATION_COVER_WORK`] evaluations, as on a fold of the
//! surface, is refused as over budget.

use super::{RelationSurface, SweptCurve};
use crate::ruled::Side;
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::Curve3;
use axiolid_reference::bound::{curve_derivative_bounds3, ROUNDING};
use std::f64::consts::TAU;

/// Cells a fixed bracket of the swept parameter is cut into up front.
const COARSE_CELLS: usize = 64;

/// The cover keeps every swept parameter within this fraction of the
/// linear tolerance of the nearest generator's distance.
const MINIMISER_SLACK: Scalar = 1e-3;

/// Distance evaluations one sample's cover may take.
const RELATION_COVER_WORK: usize = 100_000;

/// Gauss-Newton steps sharpening the swept parameter.
const POLISH_STEPS: usize = 8;

/// An off-surface sample's distance floor is refined until it is this
/// fraction of the least distance evaluated, within the work cap.
const FLOOR_SHARPNESS: Scalar = 0.9;

/// A line extruded along a direction within this fraction of its own
/// sweeps no surface.
const PARALLEL: Scalar = 1e-12;

const BUDGET: &str = "certified inversion of an edge sample onto a surface relation";

/// The generator at one swept parameter.
#[derive(Debug, Clone, Copy)]
enum Generator {
    /// The swept curve's point, which the generating circle or line passes
    /// through.
    Point(Point3),
    /// The ruling from `A(u)` to `B(u)`.
    Ruling(Point3, Point3),
}

/// A cell of the swept parameter with its generator at the midpoint.
#[derive(Debug, Clone, Copy)]
struct Coarse {
    a: Scalar,
    b: Scalar,
    speed: Scalar,
    generator: Generator,
}

/// A cell of the cover and the sample's distance from its midpoint
/// generator.
#[derive(Debug, Clone, Copy)]
struct Cell {
    a: Scalar,
    b: Scalar,
    speed: Scalar,
    distance: Scalar,
}

impl Cell {
    fn middle(&self) -> Scalar {
        0.5 * (self.a + self.b)
    }

    /// How far the distance may fall inside the cell.
    fn reach(&self) -> Scalar {
        0.5 * self.speed * (self.b - self.a)
    }

    fn lower(&self) -> Scalar {
        self.distance - self.reach()
    }
}

/// The range the swept parameter is searched over.
#[derive(Debug, Clone, Copy)]
struct Bracket {
    start: Scalar,
    end: Scalar,
    periodic: bool,
}

/// How the swept parameter is found.
#[derive(Debug, Clone)]
enum Sweep {
    /// Searched over a fixed bracket, cut into cells up front.
    Fixed(Bracket, Vec<Coarse>),
    /// A line revolved: searched over a bracket cut per sample. `origin` is
    /// the line's origin less the axis origin.
    RevolvedLine { origin: Vec3, direction: Vec3 },
    /// A line extruded: the plane it sweeps, inverted in closed form.
    ExtrudedLine { origin: Point3, direction: Vec3 },
}

/// A sample and what each distance from it reuses.
struct Probe {
    point: Point3,
    /// For a revolution, the part of `P - O` normal to the axis, its
    /// length and its height along the axis.
    perp: Vec3,
    rho: Scalar,
    height: Scalar,
}

/// Inverts edge samples onto one surface relation.
#[derive(Debug, Clone)]
pub(crate) struct RelationInverter<'g> {
    surface: RelationSurface<'g>,
    tau: Scalar,
    sweep: Sweep,
}

impl<'g> RelationInverter<'g> {
    /// The inverter of `surface` at the tolerance's linear part.
    pub(crate) fn new(surface: RelationSurface<'g>, tolerance: Tolerance) -> GeomResult<Self> {
        let tau = tolerance.linear();
        let line = match &surface {
            RelationSurface::Revolution { curve, .. }
            | RelationSurface::Extrusion { curve, .. } => match curve.curve {
                Curve3::Line(line) => Some(line),
                _ => None,
            },
            RelationSurface::Ruled { .. } => None,
        };
        let sweep = match (&surface, line) {
            (RelationSurface::Revolution { origin, .. }, Some(line)) => {
                swept_line_direction(line.direction)?;
                Sweep::RevolvedLine {
                    origin: line.origin - *origin,
                    direction: line.direction,
                }
            }
            (RelationSurface::Extrusion { direction, .. }, Some(line)) => {
                let e = swept_line_direction(line.direction)?;
                let normal = across(*direction, e);
                if normal.length() <= PARALLEL * e.length() {
                    return Err(GeomError::Degenerate(
                        "a surface of linear extrusion sweeps a line along its own direction: \
                         it is a line, not a surface"
                            .to_owned(),
                    ));
                }
                Sweep::ExtrudedLine {
                    origin: line.origin,
                    direction: e,
                }
            }
            (RelationSurface::Revolution { curve, .. }, None)
            | (RelationSurface::Extrusion { curve, .. }, None) => {
                fixed(&surface, curve_bracket(curve))?
            }
            (RelationSurface::Ruled { .. }, _) => fixed(
                &surface,
                Bracket {
                    start: 0.0,
                    end: 1.0,
                    periodic: false,
                },
            )?,
        };
        Ok(Self {
            surface,
            tau,
            sweep,
        })
    }

    /// The parameters `(u, v)` naming `point`, or a refusal by name.
    pub(crate) fn invert(&self, point: Point3) -> GeomResult<Point2> {
        if !point.is_finite() {
            return Err(GeomError::InvalidInput(
                "an edge sample to invert onto a surface relation is not finite".to_owned(),
            ));
        }
        let probe = self.probe(point)?;
        let (bracket, seeds, mut work) = match &self.sweep {
            Sweep::ExtrudedLine { origin, direction } => {
                return self.invert_plane(&probe, *origin, *direction);
            }
            Sweep::Fixed(bracket, cells) => {
                let seeds = cells
                    .iter()
                    .map(|cell| Cell {
                        a: cell.a,
                        b: cell.b,
                        speed: cell.speed,
                        distance: self.distance(&probe, cell.generator),
                    })
                    .collect();
                (*bracket, seeds, cells.len())
            }
            Sweep::RevolvedLine { origin, direction } => {
                let bracket = self
                    .line_bracket(&probe, *origin, *direction)
                    .map_err(|floor| self.off_surface(floor))?;
                let speed = direction.length() * (1.0 + ROUNDING);
                let seeds = cut(bracket)
                    .map(|(a, b)| {
                        let middle = generator(&self.surface, 0.5 * (a + b))?;
                        Ok(Cell {
                            a,
                            b,
                            speed,
                            distance: self.distance(&probe, middle),
                        })
                    })
                    .collect::<GeomResult<Vec<_>>>()?;
                (bracket, seeds, COARSE_CELLS)
            }
        };
        let (kept, floor) = self.cover(&probe, seeds, &mut work)?;
        if kept.is_empty() {
            return Err(self.off_surface(floor));
        }
        let (start, piece) = unique_piece(&kept, bracket).ok_or_else(|| self.ambiguous())?;
        let t = self.polish(&probe, start, piece)?;
        let t = if bracket.periodic {
            t.rem_euclid(bracket.end - bracket.start) + bracket.start
        } else {
            t
        };
        let uv = self.complete(&probe, t)?;
        let off = (self.surface.evaluate(uv.x, uv.y)? - point).length();
        if off.is_nan() || off > self.tau {
            return Err(self.off_surface(floor));
        }
        Ok(uv)
    }

    /// The sample, refused at a pole of its surface.
    fn probe(&self, point: Point3) -> GeomResult<Probe> {
        let mut probe = Probe {
            point,
            perp: Vec3::ZERO,
            rho: 0.0,
            height: 0.0,
        };
        match &self.surface {
            RelationSurface::Revolution { origin, axis, .. } => {
                let p = point - *origin;
                probe.height = p.dot(*axis);
                probe.perp = p - *axis * probe.height;
                probe.rho = probe.perp.length();
                if probe.rho.is_nan() || probe.rho <= self.tau {
                    return Err(GeomError::Degenerate(
                        "an edge without a p-curve runs into the axis of its surface of \
                         revolution, where every u names the point: its parameters there are \
                         not determined"
                            .to_owned(),
                    ));
                }
            }
            RelationSurface::Ruled { sides } => {
                for side in sides {
                    if let Side::Apex(apex) = side {
                        let near = (point - *apex).length();
                        if near.is_nan() || near <= self.tau {
                            return Err(GeomError::Degenerate(
                                "an edge without a p-curve runs into the apex of its ruled \
                                 surface, where every u names the point: its parameters there \
                                 are not determined"
                                    .to_owned(),
                            ));
                        }
                    }
                }
            }
            RelationSurface::Extrusion { .. } => {}
        }
        Ok(probe)
    }

    /// The sample's distance from a generator.
    fn distance(&self, probe: &Probe, generator: Generator) -> Scalar {
        match (&self.surface, generator) {
            (RelationSurface::Revolution { origin, axis, .. }, Generator::Point(c)) => {
                let q = c - *origin;
                let height = q.dot(*axis);
                let rho = (q - *axis * height).length();
                (rho - probe.rho).hypot(height - probe.height)
            }
            (RelationSurface::Extrusion { direction, .. }, Generator::Point(c)) => {
                across(*direction, c - probe.point).length()
            }
            (_, Generator::Ruling(a, b)) => {
                let d = b - a;
                let length = d.dot(d);
                let s = if length > 0.0 {
                    ((probe.point - a).dot(d) / length).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                (a + d * s - probe.point).length()
            }
            _ => Scalar::NAN,
        }
    }

    /// One half of `cell`, its distance evaluated.
    fn half(&self, probe: &Probe, cell: &Cell, a: Scalar, b: Scalar) -> GeomResult<Cell> {
        Ok(Cell {
            a,
            b,
            speed: cell.speed,
            distance: self.distance(probe, generator(&self.surface, 0.5 * (a + b))?),
        })
    }

    /// The branch-and-bound of the module notes: the leaves kept, in
    /// ascending order, and a lower bound on the sample's distance from
    /// every generator over the seeds.
    fn cover(
        &self,
        probe: &Probe,
        seeds: Vec<Cell>,
        work: &mut usize,
    ) -> GeomResult<(Vec<Cell>, Scalar)> {
        let slack = MINIMISER_SLACK * self.tau;
        let mut best = seeds
            .iter()
            .map(|cell| cell.distance)
            .fold(Scalar::INFINITY, Scalar::min);
        let mut pending = seeds;
        pending.reverse();
        let (mut kept, mut dropped) = (Vec::new(), Vec::new());
        while let Some(cell) = pending.pop() {
            if cell.lower().is_nan() || cell.lower() > self.tau.min(best + slack) {
                dropped.push(cell);
                continue;
            }
            let m = cell.middle();
            if cell.reach() <= slack || !(m > cell.a && m < cell.b) {
                kept.push(cell);
                continue;
            }
            *work += 2;
            if *work > RELATION_COVER_WORK {
                return Err(GeomError::BudgetExceeded { resource: BUDGET });
            }
            let low = self.half(probe, &cell, cell.a, m)?;
            let high = self.half(probe, &cell, m, cell.b)?;
            best = best.min(low.distance).min(high.distance);
            pending.push(high);
            pending.push(low);
        }
        let threshold = self.tau.min(best + slack);
        let (kept, late): (Vec<Cell>, Vec<Cell>) =
            kept.into_iter().partition(|cell| cell.lower() <= threshold);
        let floor = if kept.is_empty() {
            dropped.extend(late);
            self.sharpen(probe, dropped, best, work)
        } else {
            kept.iter()
                .map(Cell::lower)
                .fold(Scalar::INFINITY, Scalar::min)
        };
        Ok((kept, floor.max(0.0)))
    }

    /// A lower bound on the sample's distance from the generators over
    /// `cells`, refined until it is [`FLOOR_SHARPNESS`] of the least
    /// distance evaluated; past the work cap, the bound found so far.
    fn sharpen(
        &self,
        probe: &Probe,
        cells: Vec<Cell>,
        mut best: Scalar,
        work: &mut usize,
    ) -> Scalar {
        let mut floor = Scalar::INFINITY;
        let mut pending = cells;
        while let Some(cell) = pending.pop() {
            let m = cell.middle();
            let settled = cell.lower() >= FLOOR_SHARPNESS * best
                || !(m > cell.a && m < cell.b)
                || *work + 2 > RELATION_COVER_WORK;
            let halves = if settled {
                None
            } else {
                *work += 2;
                self.half(probe, &cell, cell.a, m)
                    .and_then(|low| Ok((low, self.half(probe, &cell, m, cell.b)?)))
                    .ok()
            };
            match halves {
                Some((low, high)) => {
                    best = best.min(low.distance).min(high.distance);
                    pending.push(high);
                    pending.push(low);
                }
                None => floor = floor.min(cell.lower()),
            }
        }
        floor
    }

    /// The parameters at the swept parameter `t`, the other one read in
    /// closed form.
    fn complete(&self, probe: &Probe, t: Scalar) -> GeomResult<Point2> {
        match &self.surface {
            RelationSurface::Revolution {
                curve,
                origin,
                axis,
            } => {
                let q = curve.point(t)? - *origin;
                let q_perp = q - *axis * q.dot(*axis);
                let u = axis
                    .cross(q_perp)
                    .dot(probe.perp)
                    .atan2(q_perp.dot(probe.perp));
                Ok(Point2::new(u, t))
            }
            RelationSurface::Extrusion { curve, direction } => {
                let c = curve.point(t)?;
                Ok(Point2::new(
                    t,
                    (probe.point - c).dot(*direction) / direction.dot(*direction),
                ))
            }
            RelationSurface::Ruled { sides } => {
                let (a, b) = (sides[0].at(t)?, sides[1].at(t)?);
                let d = b - a;
                let length = d.length();
                if length.is_nan() || length <= self.tau {
                    return Err(GeomError::Degenerate(
                        "an edge without a p-curve runs into a collapsed ruling of its ruled \
                         surface, where every v names the point: its parameters there are not \
                         determined"
                            .to_owned(),
                    ));
                }
                let v = ((probe.point - a).dot(d) / (length * length)).clamp(0.0, 1.0);
                Ok(Point2::new(t, v))
            }
        }
    }

    /// Gauss-Newton on `|S - P|^2` along the swept parameter from `start`,
    /// kept inside `piece`; the better parameter is kept at each step.
    fn polish(
        &self,
        probe: &Probe,
        start: Scalar,
        (lo, hi): (Scalar, Scalar),
    ) -> GeomResult<Scalar> {
        let residual = |uv: Point2| -> GeomResult<Vec3> {
            Ok(self.surface.evaluate(uv.x, uv.y)? - probe.point)
        };
        let mut t = start;
        let mut uv = self.complete(probe, t)?;
        let mut best = residual(uv)?.length();
        for _ in 0..POLISH_STEPS {
            if best == 0.0 {
                break;
            }
            let r = residual(uv)?;
            let (su, sv) = self.surface.tangents(uv.x, uv.y)?;
            let (a, b, c) = (su.dot(su), su.dot(sv), sv.dot(sv));
            let det = a * c - b * b;
            if !(det.is_finite() && det > 0.0) {
                break;
            }
            let (gu, gv) = (su.dot(r), sv.dot(r));
            let step = match self.surface {
                RelationSurface::Revolution { .. } => (a * gv - b * gu) / det,
                _ => (c * gu - b * gv) / det,
            };
            let next = (t - step).clamp(lo, hi);
            if !next.is_finite() {
                break;
            }
            let next_uv = self.complete(probe, next)?;
            let next_best = residual(next_uv)?.length();
            if next_best < best {
                (t, uv, best) = (next, next_uv, next_best);
            } else {
                break;
            }
        }
        Ok(t)
    }

    /// The range of a revolved line's parameter where the sample can be
    /// within `tau` of its circle, or a lower bound on the sample's
    /// distance from the surface when there is none.
    fn line_bracket(
        &self,
        probe: &Probe,
        origin: Vec3,
        direction: Vec3,
    ) -> Result<Bracket, Scalar> {
        let RelationSurface::Revolution { axis, .. } = &self.surface else {
            return Err(0.0);
        };
        let tau = self.tau;
        let everything = (Scalar::NEG_INFINITY, Scalar::INFINITY);
        let (h0, e_a) = (origin.dot(*axis), direction.dot(*axis));
        let (q_perp, e_perp) = (origin - *axis * h0, direction - *axis * e_a);
        // Every point of the line is at height `h0 + t e_a`.
        let height = if e_a != 0.0 {
            let (centre, half) = ((probe.height - h0) / e_a, tau / e_a.abs());
            (centre - half, centre + half)
        } else if (h0 - probe.height).abs() <= tau {
            everything
        } else {
            return Err((h0 - probe.height).abs());
        };
        // ... and at `|q_perp + t e_perp|` from the axis, at least `m`.
        let ee = e_perp.dot(e_perp);
        let reach = if ee > 0.0 {
            let t_m = -q_perp.dot(e_perp) / ee;
            let m = (q_perp + e_perp * t_m).length();
            let outer = probe.rho + tau;
            if m > outer {
                return Err(m - probe.rho);
            }
            let half = (outer * outer - m * m).sqrt() / ee.sqrt();
            (t_m - half, t_m + half)
        } else {
            let m = q_perp.length();
            if (m - probe.rho).abs() <= tau {
                everything
            } else {
                return Err((m - probe.rho).abs());
            }
        };
        let (start, end) = (height.0.max(reach.0), height.1.min(reach.1));
        if start.is_nan() || end.is_nan() || start > end {
            return Err(tau);
        }
        if !(start.is_finite() && end.is_finite()) {
            return Err(0.0);
        }
        Ok(Bracket {
            start,
            end,
            periodic: false,
        })
    }

    /// A sample on the plane an extruded line sweeps, in closed form.
    fn invert_plane(&self, probe: &Probe, origin: Point3, e: Vec3) -> GeomResult<Point2> {
        let RelationSurface::Extrusion { direction, .. } = &self.surface else {
            return Err(self.ambiguous());
        };
        let (e_perp, w_perp) = (
            across(*direction, e),
            across(*direction, origin - probe.point),
        );
        let t = -w_perp.dot(e_perp) / e_perp.dot(e_perp);
        let c = origin + e * t;
        let v = (probe.point - c).dot(*direction) / direction.dot(*direction);
        let off = (self.surface.evaluate(t, v)? - probe.point).length();
        if off.is_nan() || off > self.tau {
            let floor = (w_perp + e_perp * t).length();
            return Err(self.off_surface(if floor.is_nan() { 0.0 } else { floor }));
        }
        Ok(Point2::new(t, v))
    }

    fn off_surface(&self, floor: Scalar) -> GeomError {
        GeomError::Degenerate(format!(
            "an edge without a p-curve leaves its {} by at least {floor:e}: its samples must lie \
             on the surface",
            self.surface.family()
        ))
    }

    fn ambiguous(&self) -> GeomError {
        GeomError::Degenerate(format!(
            "an edge sample without a p-curve has no unique parameters on its {} (a seam of a \
             closed swept curve, a surface covered twice or a surface meeting itself)",
            self.surface.family()
        ))
    }
}

/// `x` less its part along `d`.
fn across(d: Vec3, x: Vec3) -> Vec3 {
    x - d * (x.dot(d) / d.dot(d))
}

/// A swept line's direction, refused when it has none.
fn swept_line_direction(direction: Vec3) -> GeomResult<Vec3> {
    if !(direction.is_finite() && direction.length() > 0.0) {
        return Err(GeomError::InvalidInput(
            "a swept line's direction is zero or not finite".to_owned(),
        ));
    }
    Ok(direction)
}

/// A swept B-spline's domain, or a conic's turn.
fn curve_bracket(curve: &SweptCurve<'_>) -> Bracket {
    match curve.domain {
        Some(domain) => Bracket {
            start: domain.start.min(domain.end),
            end: domain.start.max(domain.end),
            periodic: false,
        },
        None => Bracket {
            start: 0.0,
            end: TAU,
            periodic: true,
        },
    }
}

/// The bracket's [`COARSE_CELLS`] cells; neighbours share their ends.
fn cut(bracket: Bracket) -> impl Iterator<Item = (Scalar, Scalar)> {
    let at = move |i: usize| {
        if i == COARSE_CELLS {
            bracket.end
        } else {
            bracket.start + (bracket.end - bracket.start) * (i as Scalar / COARSE_CELLS as Scalar)
        }
    };
    (0..COARSE_CELLS).map(move |i| (at(i), at(i + 1)))
}

/// The generator at the swept parameter `t`.
fn generator(surface: &RelationSurface<'_>, t: Scalar) -> GeomResult<Generator> {
    match surface {
        RelationSurface::Revolution { curve, .. } | RelationSurface::Extrusion { curve, .. } => {
            Ok(Generator::Point(curve.point(t)?))
        }
        RelationSurface::Ruled { sides } => Ok(Generator::Ruling(sides[0].at(t)?, sides[1].at(t)?)),
    }
}

/// The fixed bracket's coarse cells, each with its certified speed bound
/// and its midpoint generator.
fn fixed(surface: &RelationSurface<'_>, bracket: Bracket) -> GeomResult<Sweep> {
    let cells = cut(bracket)
        .map(|(a, b)| {
            let speed = match surface {
                RelationSurface::Revolution { curve, .. }
                | RelationSurface::Extrusion { curve, .. } => {
                    curve_derivative_bounds3(curve.curve, a, b).map(|d| d.first)
                }
                RelationSurface::Ruled { sides } => sides[0]
                    .speeds(a, b)
                    .zip(sides[1].speeds(a, b))
                    .map(|((first, _), (second, _))| first.max(second)),
            };
            let speed = speed.filter(|s| s.is_finite()).ok_or_else(|| {
                GeomError::Degenerate(format!(
                    "the swept curve of a {} has no certified speed bound",
                    surface.family()
                ))
            })?;
            Ok(Coarse {
                a,
                b,
                speed: speed * (1.0 + ROUNDING),
                generator: generator(surface, 0.5 * (a + b))?,
            })
        })
        .collect::<GeomResult<Vec<_>>>()?;
    Ok(Sweep::Fixed(bracket, cells))
}

/// The start of the polish and the connected interval of the cover holding
/// it, when the leaves form one: on a periodic bracket two pieces touching
/// both ends meet across the seam, and the piece is read below the start.
/// `None` when the cover has separated parts.
fn unique_piece(kept: &[Cell], bracket: Bracket) -> Option<(Scalar, (Scalar, Scalar))> {
    let mut sorted = kept.to_vec();
    sorted.sort_by(|x, y| x.a.total_cmp(&y.a));
    let mut parts: Vec<(Scalar, Scalar)> = Vec::new();
    for cell in &sorted {
        match parts.last_mut() {
            Some(last) if cell.a <= last.1 => last.1 = last.1.max(cell.b),
            _ => parts.push((cell.a, cell.b)),
        }
    }
    let nearest = sorted
        .iter()
        .fold(None::<&Cell>, |best, cell| match best {
            Some(best) if cell.distance.is_nan() || cell.distance >= best.distance => Some(best),
            _ => Some(cell),
        })?
        .middle();
    match parts.as_slice() {
        [part] => Some((nearest, *part)),
        [first, last] if bracket.periodic && first.0 <= bracket.start && last.1 >= bracket.end => {
            let period = bracket.end - bracket.start;
            let start = if nearest > first.1 {
                nearest - period
            } else {
                nearest
            };
            Some((start, (last.0 - period, first.1)))
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::{Frame3, Interval};
    use axiolid_curve::{BSplineCurve, Circle3, Ellipse3, KnotSpec, Line3};
    use axiolid_model::{GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, RuledSide};
    use std::f64::consts::{FRAC_PI_2, PI};

    const TAU_TEST: Scalar = 1e-4;

    fn tolerance() -> Tolerance {
        Tolerance::new(TAU_TEST, 1e-9).unwrap()
    }

    fn graph_of(curves: Vec<Curve3>) -> (GeometryGraph, Vec<NodeId>) {
        let mut builder = GeometryGraphBuilder::new();
        let ids: Vec<NodeId> = curves
            .into_iter()
            .map(|curve| builder.push(GeometryNode::Curve3(curve)).unwrap())
            .collect();
        (builder.finish(vec![ids[0]]).unwrap(), ids)
    }

    fn line(origin: Point3, direction: Vec3) -> Curve3 {
        Curve3::Line(Line3 { origin, direction })
    }

    fn circle(origin: Point3, x: Vec3, y: Vec3, radius: Scalar) -> Curve3 {
        Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin,
                x,
                y,
                z: x.cross(y),
            },
            radius,
        })
    }

    fn ellipse() -> Curve3 {
        Curve3::Ellipse(Ellipse3 {
            frame: Frame3 {
                origin: Point3::new(0.6, 0.1, 0.2),
                x: Vec3::X,
                y: Vec3::new(0.0, 0.6, 0.8),
                z: Vec3::new(0.0, -0.8, 0.6),
            },
            semi_axis_x: 0.3,
            semi_axis_y: 0.15,
        })
    }

    fn cubic(points: &[Point3]) -> Curve3 {
        Curve3::BSpline(BSplineCurve {
            degree: 3,
            control_points: points.to_vec(),
            knots: vec![0.0, 1.0],
            multiplicities: vec![4, 4],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        })
    }

    fn open_cubic() -> Curve3 {
        cubic(&[
            Point3::new(0.4, 0.0, 0.0),
            Point3::new(0.6, 0.1, 0.3),
            Point3::new(0.3, -0.1, 0.6),
            Point3::new(0.5, 0.0, 0.9),
        ])
    }

    /// A far axis origin, a non-unit axis.
    fn far_axis() -> (Point3, Vec3) {
        (Point3::new(0.0, 0.0, -40.0), Vec3::new(0.0, 0.0, 3.0))
    }

    fn revolved(graph: &GeometryGraph, id: NodeId) -> RelationSurface<'_> {
        let (origin, axis) = far_axis();
        RelationSurface::revolution(graph, id, origin, axis).unwrap()
    }

    fn extruded(graph: &GeometryGraph, id: NodeId) -> RelationSurface<'_> {
        RelationSurface::extrusion(graph, id, Vec3::new(0.3, -0.2, 2.5)).unwrap()
    }

    fn sides(ids: &[NodeId], intervals: [(Scalar, Scalar); 2]) -> [RuledSide; 2] {
        [0, 1].map(|i| RuledSide::Curve {
            curve: ids[i],
            interval: Interval {
                start: intervals[i].0,
                end: intervals[i].1,
            },
        })
    }

    /// Whether `got` names `want`'s point, its parameters equal modulo the
    /// surface's periods.
    fn assert_round_trip(name: &str, surface: &RelationSurface<'_>, u: Scalar, v: Scalar) {
        let inverter = RelationInverter::new(*surface, tolerance()).unwrap();
        let point = surface.evaluate(u, v).unwrap();
        let got = inverter
            .invert(point)
            .unwrap_or_else(|e| panic!("{name} at ({u}, {v}): {e}"));
        let off = (surface.evaluate(got.x, got.y).unwrap() - point).length();
        let scale = point.length().max(1.0);
        assert!(
            off <= 1e-9 * scale,
            "{name} at ({u}, {v}): residual {off:e}"
        );
        let (pu, pv) = surface.periods();
        let same = |a: Scalar, b: Scalar, period: Option<Scalar>| {
            let d = a - b;
            let d = period.map_or(d, |p| d - p * (d / p).round());
            d.abs() <= 1e-6
        };
        assert!(
            same(got.x, u, pu) && same(got.y, v, pv),
            "{name}: ({u}, {v}) came back as {got:?}"
        );
    }

    fn grid(us: (Scalar, Scalar), vs: (Scalar, Scalar)) -> impl Iterator<Item = (Scalar, Scalar)> {
        (0..6).flat_map(move |i| {
            (0..5).map(move |j| {
                (
                    us.0 + (us.1 - us.0) * (i as Scalar + 0.3) / 6.0,
                    vs.0 + (vs.1 - vs.0) * (j as Scalar + 0.2) / 5.0,
                )
            })
        })
    }

    #[test]
    fn revolved_samples_round_trip_for_every_swept_family() {
        let torus_tube = circle(Point3::new(0.4, 0.0, 0.1), Vec3::X, Vec3::Z, 0.1);
        let slanted = line(Point3::new(0.2, 0.1, 0.0), Vec3::new(0.4, 0.0, 2.0));
        let (graph, ids) = graph_of(vec![slanted, torus_tube, ellipse(), open_cubic()]);
        let cases = [
            ("line", (-1.0, 1.0)),
            ("torus", (0.0, TAU)),
            ("ellipse", (0.0, TAU)),
            ("cubic", (0.0, 1.0)),
        ];
        for (id, (name, vs)) in ids.iter().zip(cases) {
            let surface = revolved(&graph, *id);
            for (u, v) in grid((-PI + 0.01, PI), vs) {
                assert_round_trip(name, &surface, u, v);
            }
        }
    }

    #[test]
    fn extruded_samples_round_trip_for_every_swept_family() {
        let skew = line(Point3::new(0.3, 0.1, 0.2), Vec3::new(1.5, -0.4, 0.2));
        let round = circle(Point3::new(0.1, 0.2, 0.0), Vec3::X, Vec3::Y, 0.3);
        let (graph, ids) = graph_of(vec![skew, round, ellipse(), open_cubic()]);
        let cases = [
            ("line", (-2.0, 2.0)),
            ("circle", (0.0, TAU)),
            ("ellipse", (0.0, TAU)),
            ("cubic", (0.0, 1.0)),
        ];
        for (id, (name, us)) in ids.iter().zip(cases) {
            let surface = extruded(&graph, *id);
            for (u, v) in grid(us, (-3.0, 4.0)) {
                assert_round_trip(name, &surface, u, v);
            }
        }
    }

    #[test]
    fn ruled_samples_round_trip_for_every_pair_of_sides() {
        let (graph, ids) = graph_of(vec![
            line(Point3::ZERO, Vec3::new(2.0, 0.0, 0.0)),
            line(Point3::new(0.0, 1.0, 1.0), Vec3::new(0.0, 0.0, 1.5)),
            circle(Point3::ZERO, Vec3::X, Vec3::Y, 1.0),
            line(Point3::new(1.0, 0.0, 1.0), Vec3::new(-1.0, 1.0, 0.0)),
            open_cubic(),
        ]);
        let skew = sides(&ids[0..2], [(0.0, 0.5), (0.0, 1.0)]);
        let arc_line = sides(&ids[2..4], [(0.0, FRAC_PI_2), (0.0, 1.0)]);
        let cases: [(&str, RelationSurface<'_>); 3] = [
            (
                "line/line",
                RelationSurface::ruled(&graph, &skew[0], &skew[1]).unwrap(),
            ),
            (
                "arc/line",
                RelationSurface::ruled(&graph, &arc_line[0], &arc_line[1]).unwrap(),
            ),
            (
                "curve/apex",
                RelationSurface::ruled(
                    &graph,
                    &RuledSide::Curve {
                        curve: ids[4],
                        interval: Interval {
                            start: 0.0,
                            end: 1.0,
                        },
                    },
                    &RuledSide::Apex(Point3::new(0.0, 0.0, 2.0)),
                )
                .unwrap(),
            ),
        ];
        for (name, surface) in cases {
            for (u, v) in grid((0.0, 1.0), (0.0, 0.9)) {
                assert_round_trip(name, &surface, u, v);
            }
        }
    }

    fn torus() -> (GeometryGraph, Vec<NodeId>) {
        graph_of(vec![circle(
            Point3::new(0.4, 0.0, 0.1),
            Vec3::X,
            Vec3::Z,
            0.1,
        )])
    }

    fn normal(surface: &RelationSurface<'_>, u: Scalar, v: Scalar) -> Vec3 {
        let (su, sv) = surface.tangents(u, v).unwrap();
        su.cross(sv).normalize()
    }

    fn message(error: GeomError) -> String {
        match error {
            GeomError::Degenerate(text) => text,
            other => panic!("not a named degeneracy: {other:?}"),
        }
    }

    #[test]
    fn a_sample_half_the_tolerance_off_still_inverts() {
        let (graph, ids) = torus();
        let surface = revolved(&graph, ids[0]);
        let inverter = RelationInverter::new(surface, tolerance()).unwrap();
        let (u, v) = (0.7, 2.0);
        let point = surface.evaluate(u, v).unwrap() + normal(&surface, u, v) * (0.5 * TAU_TEST);
        let got = inverter.invert(point).unwrap();
        assert!(
            (got.x - u).abs() <= 1e-6 && (got.y - v).abs() <= 1e-6,
            "{got:?}"
        );
    }

    #[test]
    fn a_sample_five_tolerances_off_is_refused_with_a_floor() {
        let (graph, ids) = torus();
        let surface = revolved(&graph, ids[0]);
        let inverter = RelationInverter::new(surface, tolerance()).unwrap();
        let (u, v) = (0.7, 2.0);
        let point = surface.evaluate(u, v).unwrap() + normal(&surface, u, v) * (5.0 * TAU_TEST);
        let text = message(inverter.invert(point).unwrap_err());
        let floor: Scalar = text
            .strip_prefix("an edge without a p-curve leaves its surface of revolution by at least ")
            .and_then(|rest| rest.split(':').next())
            .unwrap_or_else(|| panic!("{text}"))
            .parse()
            .unwrap();
        assert!(
            (4.0 * TAU_TEST..=5.0 * TAU_TEST).contains(&floor),
            "{floor:e}"
        );
    }

    #[test]
    fn a_sample_just_before_a_ruling_reads_v_zero() {
        let (graph, ids) = graph_of(vec![
            line(Point3::ZERO, Vec3::X),
            line(Point3::new(0.0, 1.0, 1.0), Vec3::Z),
        ]);
        let pair = sides(&ids, [(0.0, 1.0), (0.0, 1.0)]);
        let surface = RelationSurface::ruled(&graph, &pair[0], &pair[1]).unwrap();
        let inverter = RelationInverter::new(surface, tolerance()).unwrap();
        let u = 0.4;
        let (a, b) = (
            surface.evaluate(u, 0.0).unwrap(),
            surface.evaluate(u, 1.0).unwrap(),
        );
        let point = a - (b - a).normalize() * (0.5 * TAU_TEST);
        let got = inverter.invert(point).unwrap();
        assert_eq!(got.y, 0.0);
        assert!((got.x - u).abs() <= 1e-3, "{got:?}");
    }

    fn refusal(surface: RelationSurface<'_>, point: Point3) -> String {
        let inverter = RelationInverter::new(surface, tolerance()).unwrap();
        message(inverter.invert(point).unwrap_err())
    }

    #[test]
    fn r1_a_sample_on_the_axis_is_refused() {
        let (graph, ids) = torus();
        let text = refusal(revolved(&graph, ids[0]), Point3::new(0.0, 0.0, 0.05));
        assert!(
            text.contains("runs into the axis of its surface of revolution"),
            "{text}"
        );
    }

    #[test]
    fn r2_a_sample_at_the_apex_is_refused() {
        let (graph, ids) = graph_of(vec![open_cubic()]);
        let apex = Point3::new(0.0, 0.0, 2.0);
        let surface = RelationSurface::ruled(
            &graph,
            &RuledSide::Curve {
                curve: ids[0],
                interval: Interval {
                    start: 0.0,
                    end: 1.0,
                },
            },
            &RuledSide::Apex(apex),
        )
        .unwrap();
        let text = refusal(surface, apex + Vec3::X * (0.5 * TAU_TEST));
        assert!(
            text.contains("runs into the apex of its ruled surface"),
            "{text}"
        );
    }

    #[test]
    fn r3_a_sample_on_a_collapsed_ruling_is_refused() {
        let (graph, ids) = graph_of(vec![
            line(Point3::ZERO, Vec3::X),
            line(Point3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, -1.0)),
        ]);
        let pair = sides(&ids, [(0.0, 1.0), (0.0, 1.0)]);
        let surface = RelationSurface::ruled(&graph, &pair[0], &pair[1]).unwrap();
        let text = refusal(surface, Point3::new(1.0, 0.0, 0.0));
        assert!(
            text.contains("runs into a collapsed ruling of its ruled surface"),
            "{text}"
        );
    }

    #[test]
    fn r4_a_revolved_line_s_sample_off_its_bracket_is_refused_with_a_floor() {
        let (graph, ids) = graph_of(vec![line(Point3::new(0.3, 0.0, 0.0), Vec3::Z)]);
        let surface = RelationSurface::revolution(&graph, ids[0], Point3::ZERO, Vec3::Z).unwrap();
        let text = refusal(surface, Point3::new(0.0, 0.5, 0.2));
        assert!(
            text.starts_with(
                "an edge without a p-curve leaves its surface of revolution by at least 2e-1"
            ),
            "{text}"
        );
    }

    #[test]
    fn r5_a_surface_covered_twice_is_refused() {
        let (graph, ids) = graph_of(vec![circle(Point3::ZERO, Vec3::X, Vec3::Z, 0.3)]);
        let surface = RelationSurface::revolution(&graph, ids[0], Point3::ZERO, Vec3::Z).unwrap();
        let point = surface.evaluate(0.4, 0.3).unwrap();
        let text = refusal(surface, point);
        assert!(
            text.contains("has no unique parameters on its surface of revolution"),
            "{text}"
        );
    }

    #[test]
    fn r6_a_profile_running_along_the_direction_exhausts_the_cover() {
        let profile = Curve3::BSpline(BSplineCurve {
            degree: 1,
            control_points: vec![
                Point3::ZERO,
                Point3::new(0.0, 0.0, 1.0),
                Point3::new(1.0, 0.0, 1.0),
            ],
            knots: vec![0.0, 0.5, 1.0],
            multiplicities: vec![2, 1, 2],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        });
        let (graph, ids) = graph_of(vec![profile]);
        let surface = RelationSurface::extrusion(&graph, ids[0], Vec3::Z).unwrap();
        let inverter = RelationInverter::new(surface, tolerance()).unwrap();
        let error = inverter.invert(Point3::new(0.0, 0.0, 0.5)).unwrap_err();
        assert!(
            matches!(error, GeomError::BudgetExceeded { resource } if resource == BUDGET),
            "{error:?}"
        );
    }

    #[test]
    fn r7_a_line_extruded_along_itself_is_refused() {
        let (graph, ids) = graph_of(vec![line(Point3::ZERO, Vec3::new(0.0, 0.0, 2.0))]);
        let surface = RelationSurface::extrusion(&graph, ids[0], Vec3::Z).unwrap();
        let text = message(RelationInverter::new(surface, tolerance()).unwrap_err());
        assert!(
            text.contains("sweeps a line along its own direction"),
            "{text}"
        );
    }

    #[test]
    fn a_swept_line_without_a_direction_is_invalid() {
        let (graph, ids) = graph_of(vec![line(Point3::X, Vec3::ZERO)]);
        let turned = RelationSurface::revolution(&graph, ids[0], Point3::ZERO, Vec3::Z).unwrap();
        let swept = RelationSurface::extrusion(&graph, ids[0], Vec3::Z).unwrap();
        for surface in [turned, swept] {
            let error = RelationInverter::new(surface, tolerance()).unwrap_err();
            assert!(matches!(error, GeomError::InvalidInput(_)), "{error:?}");
        }
    }

    #[test]
    fn a_seam_sample_of_a_periodic_bracket_joins_across_it() {
        let (graph, ids) = torus();
        let surface = revolved(&graph, ids[0]);
        for v in [0.0, 1e-9, TAU - 1e-9] {
            assert_round_trip("torus seam", &surface, 0.3, v);
        }
    }

    #[test]
    fn the_same_sample_inverts_to_the_same_parameters() {
        let (graph, ids) = graph_of(vec![open_cubic()]);
        let surface = extruded(&graph, ids[0]);
        let inverter = RelationInverter::new(surface, tolerance()).unwrap();
        let point = surface.evaluate(0.37, 1.2).unwrap();
        let first = inverter.invert(point).unwrap();
        assert_eq!(first, inverter.invert(point).unwrap());
    }
}
