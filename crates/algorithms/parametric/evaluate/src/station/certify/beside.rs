//! Certified lengths of offsets beside a base in closed form: a line, a
//! circle, an ellipse or a B-spline (#294, #307).
//!
//! The offset `P(t) = B(t) + D(t)` is measured by chords at the base's
//! native parameter ([`super::panels`]). `P''` and `P'''` are bounded from
//! the base's derivative suprema `|B'|` to `|B''''|` ([`crate::bound`])
//! and the frame's: each axis of the displacement is a unit vector `w /
//! |w|` with `w` the base's tangent, its plan projection or `V x B'`
//! ([`super::panels::unit_jet`]), and each component `c_k` of the law is
//! a factor, so `D^(j)` is bounded by Leibniz from the suprema of `c_k`
//! and of the axes.
//!
//! # A law that varies along a base measured by quadrature (#307)
//!
//! An offset by distances interpolates its law linearly in the base's
//! measure. On a line or a circle that measure is the native parameter
//! times a constant. On an ellipse or a B-spline it is the base's arc
//! length, a quadrature, so the law at a native parameter `t` is read
//! through the base's own certified length: the base is measured by
//! chords first (a table of panels), and `lambda(t)`, the base's arc
//! length from the piece's first end to `t`, is known to an interval
//! (the panels before `t` and a sub-panel to it). The measured curve is
//! the offset whose law runs linearly in that true arc length between the
//! points the kernel reads as the piece's ends: `v(t) = va + (vb - va)
//! lambda(t) / Lambda`, `Lambda` the base's true length between them, so
//! it passes exactly through both ends. A point read at the interval's
//! middle is off that curve by at most `|law rate| (vb - va)` times the
//! ratio's half-width, an error the chain of chords carries like rounding.
//! Its derivatives: `v' = (vb - va) |B'| / Lambda`, `v'' = (vb - va) (u .
//! B'') / Lambda`, `|v'''| <= (vb - va) (|B''|^2 / |B'| + |B'''|) /
//! Lambda`, with `Lambda` bounded below by the table.

use std::cell::Cell;

use axiolid_contracts::GeomResult;
use axiolid_core::{Scalar, Vec3};
use axiolid_curve::{OffsetFrame, OffsetLaw};

use super::panels::{
    product, unit_jet, unit_rate, Compensated, Jet, Measured, Panel, PANEL_SHARE, VELOCITY_ROUNDING,
};
use super::{uncertified, DistanceInterval, Native};
use crate::station::offset::StationOffset;
use crate::station::reference_up_frame;

/// The total turning of a chain of chords the base table's tolerance is
/// chosen for, with the `2` of its ends ([`super::panels`]): the table is
/// made fine enough that a point's error, carried along a chain that turns
/// this much, takes a small share of the tolerance. A chain that turns
/// more is still bounded exactly; it only spends more of the tolerance.
const TURNING: Scalar = 2.0 + 2.0 * core::f64::consts::PI;

/// The base's certified arc length at every native parameter of the
/// piece, from its chord panels.
pub(super) struct BaseTable {
    panels: Vec<Panel>,
    /// Lower and upper bounds of the length at each panel's start, and at
    /// the last panel's end.
    lower: Vec<Scalar>,
    upper: Vec<Scalar>,
    /// The rounding of the chain and of its prefix sums.
    slack: Scalar,
}

impl BaseTable {
    pub(super) fn build(
        base: &Measured<'_>,
        t0: Scalar,
        t1: Scalar,
        breaks: &[Scalar],
        target: Scalar,
    ) -> GeomResult<Self> {
        let panels = base.chords(t0, t1, breaks, target)?.panels;
        if panels.is_empty() {
            return Err(uncertified());
        }
        let (mut lower, mut upper) = (vec![0.0], vec![0.0]);
        let (mut low, mut high) = (Compensated::default(), Compensated::default());
        for panel in &panels {
            low.add(panel.chord + panel.low);
            high.add(panel.chord + panel.high);
            lower.push(low.value());
            upper.push(high.value());
        }
        let (chain_low, _) = base.chain(&panels);
        let rounding = (low.value() - chain_low).max(0.0);
        let slack = rounding + 8.0 * Scalar::EPSILON * high.value();
        Ok(Self {
            panels,
            lower,
            upper,
            slack,
        })
    }

    /// Bounds of the base's whole length between the piece's ends.
    pub(super) fn total(&self) -> (Scalar, Scalar) {
        let n = self.panels.len();
        (self.lower[n] - self.slack, self.upper[n] + self.slack)
    }

    /// Bounds of `lambda(t)`, the base's arc length from the first panel's
    /// start to `t`.
    pub(super) fn at(&self, base: &Measured<'_>, t: Scalar) -> GeomResult<(Scalar, Scalar)> {
        let n = self.panels.len();
        if t <= self.panels[0].a.t {
            return Ok((0.0, 0.0));
        }
        if t >= self.panels[n - 1].b.t {
            return Ok(self.total());
        }
        let i = self.panels.partition_point(|p| p.b.t <= t).min(n - 1);
        let panel = self.panels[i];
        if t == panel.a.t {
            return Ok((self.lower[i] - self.slack, self.upper[i] + self.slack));
        }
        // Monotone: the length to `t` lies between the panel's ends'.
        let whole = (self.lower[i] - self.slack, self.upper[i + 1] + self.slack);
        let part = base.panel(panel.a, base.node(t)?)?;
        let rounding = 2.0 * base.point_error(&panel.a).max(base.point_error(&part.b))
            + 4.0 * Scalar::EPSILON * part.chord;
        let lo = self.lower[i] + part.chord + part.low - rounding - self.slack;
        let hi = self.upper[i] + part.chord + part.high + rounding + self.slack;
        Ok(if lo.is_finite() && hi.is_finite() {
            (lo.max(whole.0), hi.min(whole.1))
        } else {
            whole
        })
    }
}

/// Where along its base piece a law is read at a native parameter.
enum Progress {
    /// A constant law: anywhere.
    Fixed,
    /// The base distance is the native parameter times this (a line's
    /// `|d|`, a circle's radius).
    Rate(Scalar),
    /// Through the base's certified arc length, and the bounds of the
    /// base's length between the piece's ends.
    Table(BaseTable, (Scalar, Scalar)),
}

/// The certified length of an offset of a closed-form base between its
/// base measures `va < vb`, by chords; see the
/// [module documentation](self).
pub(super) fn chord_offset_length(
    offset: &StationOffset<'_>,
    va: Scalar,
    vb: Scalar,
    tolerance: Scalar,
) -> GeomResult<DistanceInterval> {
    let Some(beside) = Beside::new(offset, va, vb, tolerance)? else {
        return Ok(DistanceInterval::point(0.0));
    };
    beside.measured(|curve| curve.chord_turning(beside.t0, beside.t1, &beside.breaks, tolerance))
}

/// An offset of a closed-form base between its base measures `va < vb`
/// (into the base piece), as the curve its length and its frames are
/// certified on (#307, #316): the native parameters of its ends, and its
/// law with how it is read along the base.
pub(super) struct Beside<'c> {
    offset: StationOffset<'c>,
    pub(super) native: Native<'c>,
    progress: Progress,
    /// The native parameters of the ends, ascending.
    pub(super) t0: Scalar,
    pub(super) t1: Scalar,
    /// Where the measured curve's third derivative may jump.
    pub(super) breaks: Vec<Scalar>,
    va: Scalar,
    vb: Scalar,
    /// Each law component's rate per unit base distance, and its size.
    rates: [Scalar; 3],
    sizes: [Scalar; 3],
}

impl<'c> Beside<'c> {
    /// The offset between `va < vb`, its base measured by chords for
    /// `tolerance` where its law varies along a base measured by
    /// quadrature; `None` for an empty span.
    pub(super) fn new(
        offset: &StationOffset<'c>,
        va: Scalar,
        vb: Scalar,
        tolerance: Scalar,
    ) -> GeomResult<Option<Self>> {
        let base = offset.base();
        let measure = |v: Scalar| {
            if base.reversed {
                base.end - v
            } else {
                base.start + v
            }
        };
        let (ma, mb) = (measure(va), measure(vb));
        let native = Native::of(base.curve, ma.min(mb), ma.max(mb))?;
        let (ta, tb) = (native.parameter(ma)?, native.parameter(mb)?);
        let (t0, t1) = (ta.min(tb), ta.max(tb));
        if t1.is_nan() || t1 <= t0 {
            return Ok(None);
        }
        let mut breaks = native.breaks(3);
        breaks.extend(native.breaks(1));
        let law = offset.law();
        let length = offset.base_length();
        let pairs = match law {
            OffsetLaw::Linear { start, end, .. } => [
                (start.lateral, end.lateral),
                (start.vertical, end.vertical),
                (start.longitudinal, end.longitudinal),
            ],
            _ => [(0.0, 0.0); 3],
        };
        // Each component's rate per unit base distance, and its size.
        let rates = pairs.map(|(a, b)| (b - a) / length);
        let sizes = pairs.map(|(a, b)| a.abs().max(b.abs()));
        let rsum: Scalar = rates.iter().map(|r| r.abs()).sum();
        let progress = if law.is_constant() {
            Progress::Fixed
        } else if let Some(rate) = native.rate() {
            Progress::Rate(rate)
        } else {
            native.measured(|base_curve| {
                let target = (PANEL_SHARE * tolerance).min(tolerance / (16.0 * rsum * TURNING));
                let table = BaseTable::build(base_curve, t0, t1, &breaks, target)?;
                let total = table.total();
                if !(total.0 > 0.0 && total.1.is_finite()) {
                    return Err(uncertified());
                }
                Ok(Progress::Table(table, total))
            })?
        };
        Ok(Some(Self {
            offset: *offset,
            native,
            progress,
            t0,
            t1,
            breaks,
            va,
            vb,
            rates,
            sizes,
        }))
    }

    /// The direction of travel along the base's native parameter: `-1` on
    /// a reversed base.
    pub(super) fn travel(&self) -> Scalar {
        if self.offset.base().reversed {
            -1.0
        } else {
            1.0
        }
    }

    /// `run` given the offset measured by chords at its base's native
    /// parameter, its law read as the [module documentation](self) says.
    pub(super) fn measured<R>(
        &self,
        run: impl FnOnce(&Measured<'_>) -> GeomResult<R>,
    ) -> GeomResult<R> {
        self.native.measured(|base_curve| {
            offset_chords(
                &self.offset,
                &self.native,
                base_curve,
                &self.progress,
                (self.t0, self.t1),
                (self.va, self.vb, self.rates, self.sizes),
                run,
            )
        })
    }
}

/// The offset itself measured by chords, its law read by `progress`, handed
/// to `run`.
fn offset_chords<R>(
    offset: &StationOffset<'_>,
    native: &Native<'_>,
    base_curve: &Measured<'_>,
    progress: &Progress,
    (t0, t1): (Scalar, Scalar),
    (va, vb, rates, sizes): (Scalar, Scalar, [Scalar; 3], [Scalar; 3]),
    run: impl FnOnce(&Measured<'_>) -> GeomResult<R>,
) -> GeomResult<R> {
    let base = offset.base();
    let law = offset.law();
    let travel = if base.reversed { -1.0 } else { 1.0 };
    let placement = base.placement;
    let length = offset.base_length();
    let span = vb - va;
    let rsum: Scalar = rates.iter().map(|r| r.abs()).sum();
    let move_point = |p| placement.map_or(p, |rigid| rigid.transform_point3(p));
    let move_vector = |v: Vec3| placement.map_or(v, |rigid| rigid.transform_vector3(v));
    // The base distance into the piece at `t`, and the bound of its error.
    let memo = Cell::new((Scalar::NAN, (0.0, 0.0)));
    let read = |t: Scalar| -> GeomResult<(Scalar, Scalar)> {
        let (at, value) = memo.get();
        if at == t {
            return Ok(value);
        }
        let value = match progress {
            Progress::Fixed => (0.0, 0.0),
            Progress::Rate(rate) => {
                let m = t * rate;
                (
                    if base.reversed {
                        base.end - m
                    } else {
                        m - base.start
                    },
                    0.0,
                )
            }
            Progress::Table(table, (whole_lo, whole_hi)) => {
                let ratio = if t <= t0 {
                    (0.0, 0.0)
                } else if t >= t1 {
                    (1.0, 1.0)
                } else {
                    let (lo, hi) = table.at(base_curve, t)?;
                    (
                        (lo / whole_hi).clamp(0.0, 1.0),
                        (hi / whole_lo).clamp(0.0, 1.0),
                    )
                };
                let middle = 0.5 * (ratio.0 + ratio.1);
                let v = if base.reversed {
                    vb - span * middle
                } else {
                    va + span * middle
                };
                (v, 0.5 * span * (ratio.1 - ratio.0))
            }
        };
        memo.set((t, value));
        Ok(value)
    };
    // `dv/dt` at `t`, `d` the base's velocity there, and its error.
    let pace = |d: Vec3| -> (Scalar, Scalar) {
        match progress {
            Progress::Fixed => (0.0, 0.0),
            Progress::Rate(rate) => (travel * rate, 0.0),
            Progress::Table(_, (lo, hi)) => {
                let speed = d.length();
                (
                    travel * span * speed * 2.0 / (lo + hi),
                    0.5 * span * speed * (1.0 / lo - 1.0 / hi),
                )
            }
        }
    };
    let point = |t: Scalar| -> GeomResult<(axiolid_core::Point3, Scalar)> {
        let p = move_point(native.point(t)?);
        let d = move_vector(native.derivative(t)?);
        let section = reference_up_frame(p, travel * d)?;
        let (v, error) = read(t)?;
        Ok((
            section.point + offset.displacement(&section, v)?,
            rsum * error,
        ))
    };
    // `P' = B' + D'`, the displacement's frame differentiated through the
    // unit vectors it is built of.
    let velocity = |t: Scalar| -> GeomResult<(Vec3, Scalar)> {
        let d = move_vector(native.derivative(t)?);
        let (w, dw) = (
            travel * d,
            travel * move_vector(native.second_derivative(t)?),
        );
        let (tangent, dtangent) = unit_rate(w, dw)?;
        let (lateral, dlateral) = unit_rate(Vec3::Z.cross(w), Vec3::Z.cross(dw))?;
        let up = (
            tangent.cross(lateral),
            dtangent.cross(lateral) + tangent.cross(dlateral),
        );
        let mut slip = 0.0;
        let displacement = match law {
            OffsetLaw::Planar { distance } => distance * dlateral,
            OffsetLaw::Directed {
                distance,
                reference_direction,
            } => {
                let v = reference_direction.normalize();
                distance * unit_rate(v.cross(w), v.cross(dw))?.1
            }
            OffsetLaw::Linear { start, end, frame } => {
                let axes = match frame {
                    OffsetFrame::Section => [(lateral, dlateral), up, (tangent, dtangent)],
                    OffsetFrame::Plan => {
                        let flat = |v: Vec3| Vec3::new(v.x, v.y, 0.0);
                        let (along, dalong) = unit_rate(flat(w), flat(dw))?;
                        [
                            (Vec3::Z.cross(along), Vec3::Z.cross(dalong)),
                            (Vec3::Z, Vec3::ZERO),
                            (along, dalong),
                        ]
                    }
                    _ => return Err(uncertified()),
                };
                let (v, error) = read(t)?;
                let (rate, rate_error) = pace(d);
                let u = (v / length).clamp(0.0, 1.0);
                let pairs = [
                    (start.lateral, end.lateral),
                    (start.vertical, end.vertical),
                    (start.longitudinal, end.longitudinal),
                ];
                let mut sum = Vec3::ZERO;
                for (k, (&(a, b), (axis, daxis))) in pairs.iter().zip(axes).enumerate() {
                    sum += (b - a) / length * rate * axis + (a + (b - a) * u) * daxis;
                    slip += rates[k].abs() * (rate_error + error * daxis.length());
                }
                sum
            }
            _ => return Err(uncertified()),
        };
        let out = d + displacement;
        Ok((out, slip + VELOCITY_ROUNDING * out.length()))
    };
    let bounds = |a: Scalar, b: Scalar| -> GeomResult<(Scalar, Option<Scalar>)> {
        let derivatives = native.bounds(a, b)?;
        let fourth = native.fourth(a, b).unwrap_or(Scalar::INFINITY);
        let (w1, w2, w3) = (derivatives.first, derivatives.second, derivatives.third);
        let h = b - a;
        let d = move_vector(native.derivative(0.5 * (a + b))?);
        let reach = 0.5 * h * w2;
        let flat = Vec3::new(d.x, d.y, 0.0).length();
        let floor = d.length() - reach;
        let tangent = unit_jet(floor, w2, w3, fourth);
        let lateral = unit_jet(flat - reach, w2, w3, fourth);
        let up = product(tangent, lateral);
        // A law component's suprema: its size and the derivatives of its
        // reading along the base.
        let component = |k: usize| -> Jet<4> {
            let r = rates[k].abs();
            match progress {
                Progress::Fixed => [sizes[k], 0.0, 0.0, 0.0],
                Progress::Rate(rate) => [sizes[k], r * rate, 0.0, 0.0],
                Progress::Table(_, (lo, _)) => {
                    let rho = r * span / lo;
                    let bend = if floor > 0.0 {
                        w2 * w2 / floor + w3
                    } else {
                        Scalar::INFINITY
                    };
                    [sizes[k], rho * w1, rho * w2, rho * bend]
                }
            }
        };
        let displacement = match law {
            OffsetLaw::Planar { distance } => product([distance.abs(), 0.0, 0.0, 0.0], lateral),
            OffsetLaw::Directed {
                distance,
                reference_direction,
            } => {
                let v = reference_direction.normalize();
                let r = v.cross(d).length() - reach;
                product([distance.abs(), 0.0, 0.0, 0.0], unit_jet(r, w2, w3, fourth))
            }
            OffsetLaw::Linear { frame, .. } => {
                let axes = match frame {
                    OffsetFrame::Section => [lateral, up, tangent],
                    OffsetFrame::Plan => [lateral, [1.0, 0.0, 0.0, 0.0], lateral],
                    _ => return Err(uncertified()),
                };
                (0..3).fold([0.0; 4], |acc, k| {
                    let term = product(component(k), axes[k]);
                    core::array::from_fn(|j| acc[j] + term[j])
                })
            }
            _ => return Err(uncertified()),
        };
        let third = w3 + displacement[3];
        Ok((w2 + displacement[2], (third.is_finite()).then_some(third)))
    };
    let reach = sizes.iter().sum::<Scalar>()
        + match law {
            OffsetLaw::Planar { distance } | OffsetLaw::Directed { distance, .. } => distance.abs(),
            _ => 0.0,
        };
    let shift = placement.map_or(0.0, |rigid| rigid.translation.abs().max_element());
    run(&Measured {
        point: &point,
        velocity: &velocity,
        bounds: &bounds,
        scale: native.extent() + reach + shift,
    })
}
