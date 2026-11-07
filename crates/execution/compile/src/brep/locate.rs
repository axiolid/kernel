//! Locating a shared edge's samples on one use's pcurve (#413).
//!
//! The edge pass ([`super::edges`]) samples an edge once, through one
//! pcurve use. Another use of the edge, on the neighbouring face, runs its
//! own pcurve over the same 3D curve, and nothing makes the two walk it at
//! the same pace: a traced section, an IFC or Revit pcurve and the pcurve of
//! a cap whose circle starts elsewhere than the wall's rim all parameterise
//! one edge differently. Reading the other pcurve at the same uniform
//! parameters then lands off the shared samples.
//!
//! So each cached sample `X_i` is located on the other pcurve: the
//! parameter `t_i` minimising `|S(p(t)) - X_i|` is searched with Brent's
//! method (the squared distance, which is smooth where the distance has a
//! corner) over `[t_{i-1}, guess]` only. `t_{i-1}` is the sample before,
//! already located, so the parameters stay strictly monotone; `guess` is
//! where the uniform spacing of the samples left would put the sample after
//! this one. A minimum at either end of that bracket is not a minimum of the
//! bracket, so it refuses. The two ends of the use are not searched: they
//! are the use's own span ends.
//!
//! The vertices are the shared samples whatever the search found, so the
//! result is a bound, not an assumption: a sample's measured distance from
//! the surface at its located parameters is refused above the chord budget
//! (an end above the budget the caller admits), and below it is the vertex
//! offset the triangles' bound of the face already carries.
//!
//! A located parameter within a thousandth of a step of the sample's uniform
//! parameter is put on it: a pcurve at the pace of the sampling use then
//! triangulates exactly as it would at the uniform parameters, not at the
//! search's rounding of them.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Scalar, Vec3};

use super::support::Support;
use super::view::Pcurve;

/// Why a shared sample has no place on a neighbouring pcurve.
pub(super) const NOT_LOCATED: &str = "a shared edge sample is not located on a neighbouring pcurve";

/// How far from either end of its bracket a located parameter must stay,
/// as a fraction of the bracket: a minimum nearer is the bracket's end.
const END_MARGIN: Scalar = 1e-6;
/// How near its uniform parameter, as a fraction of the step from the sample
/// before, a located parameter is put on it.
const SNAP: Scalar = 1e-3;
/// The relative and absolute parameter tolerance of the search.
const REL_TOL: Scalar = 1e-12;
const ABS_TOL: Scalar = 1e-15;
const MAX_ITERATIONS: usize = 100;
/// How often a bracket whose minimum is at its far end is doubled.
const MAX_WIDENINGS: usize = 4;
/// The golden section's smaller part, `(3 - sqrt 5) / 2`.
const GOLDEN: Scalar = 0.381_966_011_250_105_1;

/// The samples of one edge use, located on its pcurve.
#[derive(Debug, Clone)]
pub(super) struct Located {
    /// The pcurve parameter of each sample, both ends included, strictly
    /// monotone from the span's start to its end.
    pub ts: Vec<Scalar>,
    /// The pcurve's point at each parameter.
    pub params: Vec<Point2>,
}

/// How far a located sample may sit from its pcurve's image.
#[derive(Debug, Clone, Copy)]
pub(super) struct Limits {
    /// For a sample inside the edge: the chord budget.
    pub interior: Scalar,
    /// For the use's two ends, which an exact boolean closes only within its
    /// reading (#334).
    pub end: Scalar,
}

/// The pcurve parameter at fraction `s` of its span.
fn at(pcurve: &Pcurve<'_>, s: Scalar) -> Scalar {
    pcurve.span.start + (pcurve.span.end - pcurve.span.start) * s
}

/// The `i`th of `n` uniform parameters, as [`super::trim_samples`] takes them.
fn uniform(pcurve: &Pcurve<'_>, i: usize, n: usize) -> Scalar {
    pcurve.span.start + (pcurve.span.end - pcurve.span.start) * (i as Scalar) / (n as Scalar)
}

/// The pcurve's point at `t` and the surface point over it.
fn image(pcurve: &Pcurve<'_>, surface: Support<'_>, t: Scalar) -> GeomResult<(Point2, Vec3)> {
    let p = axiolid_reference::curve::evaluate2(pcurve.curve, t)?;
    Ok((p, surface.evaluate(p.x, p.y)?))
}

/// Locate `targets` (both ends included, in the use's direction) on
/// `pcurve` over `surface`, or refuse by name.
pub(super) fn locate_samples(
    pcurve: &Pcurve<'_>,
    surface: Support<'_>,
    targets: &[Vec3],
    limits: Limits,
) -> GeomResult<Located> {
    let n = match targets.len().checked_sub(1) {
        Some(n) if n >= 1 => n,
        _ => {
            return Err(GeomError::Degenerate(
                "curved edge needs both endpoint samples".to_owned(),
            ))
        }
    };
    let mut ts = Vec::with_capacity(n + 1);
    let mut params = Vec::with_capacity(n + 1);
    // The fraction of the span the last located sample sits at.
    let mut previous = 0.0;
    for (i, &target) in targets.iter().enumerate() {
        let (t, s) = if i == 0 || i == n {
            (uniform(pcurve, i, n), if i == 0 { 0.0 } else { 1.0 })
        } else {
            let t = uniform(pcurve, i, n);
            let s = i as Scalar / n as Scalar;
            let (_, point) = image(pcurve, surface, t)?;
            // The sample the pcurve itself produced (the source use) sits
            // exactly on its uniform parameter.
            if s > previous && point == target {
                (t, s)
            } else {
                let found = search(pcurve, surface, (targets[i - 1], target), previous, (i, n))?;
                let located = at(pcurve, found);
                let before = ts.last().copied().unwrap_or(pcurve.span.start);
                let ahead = (pcurve.span.end - pcurve.span.start).signum();
                if (t - before) * ahead > 0.0
                    && (located - t).abs() <= SNAP * (located - before).abs()
                {
                    // A pcurve at the pace of the sampling use: its uniform
                    // parameters, not the search's rounding of them.
                    (t, s)
                } else {
                    (located, found)
                }
            }
        };
        let (p, point) = image(pcurve, surface, t)?;
        let miss = (point - target).length();
        let limit = if i == 0 || i == n {
            limits.end
        } else {
            limits.interior
        };
        if !miss.is_finite() || miss > limit {
            return Err(GeomError::Degenerate(NOT_LOCATED.to_owned()));
        }
        if let Some(&before) = ts.last() {
            let ahead = (pcurve.span.end - pcurve.span.start).signum();
            if (t - before) * ahead <= 0.0 {
                return Err(GeomError::Degenerate(NOT_LOCATED.to_owned()));
            }
        }
        ts.push(t);
        params.push(p);
        previous = s;
    }
    Ok(Located { ts, params })
}

/// The fraction of the span nearest to the later of `targets`, searched
/// between the sample before (`previous`) and a guess of where the sample
/// after sample `i` of `n` lies.
///
/// The guess is the uniform spacing of the samples left, scaled by how much
/// faster the pcurve walks the edge than that spacing does where the sample
/// before sits (its 3D travel over one uniform step, against the distance
/// to travel): a pcurve a few times slower or faster than the sampling use
/// is still bracketed, one that changes pace more than that within a step
/// is not.
fn search(
    pcurve: &Pcurve<'_>,
    surface: Support<'_>,
    (before, target): (Vec3, Vec3),
    previous: Scalar,
    (i, n): (usize, usize),
) -> GeomResult<Scalar> {
    let uniform = (1.0 - previous) / (n - i + 1) as Scalar;
    let (_, from) = image(pcurve, surface, at(pcurve, previous))?;
    let (_, ahead) = image(pcurve, surface, at(pcurve, previous + uniform))?;
    let walked = (ahead - from).length();
    let to_walk = (target - before).length();
    let step = if walked > 0.0 && to_walk.is_finite() {
        uniform * to_walk / walked
    } else {
        uniform
    };
    let cost = |s: Scalar| -> GeomResult<Scalar> {
        let (_, point) = image(pcurve, surface, at(pcurve, s))?;
        Ok((point - target).length_squared())
    };
    // A minimum at the far end of the bracket is not the bracket's: widen
    // it, a few times, before refusing.
    let lo = previous;
    let mut width = 2.0 * step;
    for _ in 0..=MAX_WIDENINGS {
        let hi = (lo + width).min(1.0);
        let found = minimise(cost, lo, hi, lo + 0.5 * (hi - lo))?;
        let margin = END_MARGIN * (hi - lo);
        if found - lo <= margin {
            break;
        }
        if hi - found > margin {
            return Ok(found);
        }
        if hi >= 1.0 {
            break;
        }
        width *= 2.0;
    }
    Err(GeomError::Degenerate(NOT_LOCATED.to_owned()))
}

/// Brent's minimiser of `f` over `(lo, hi)`, started at `start`.
fn minimise(
    mut f: impl FnMut(Scalar) -> GeomResult<Scalar>,
    lo: Scalar,
    hi: Scalar,
    start: Scalar,
) -> GeomResult<Scalar> {
    let (mut a, mut b) = (lo, hi);
    let (mut x, mut w, mut v) = (start, start, start);
    let mut fx = f(x)?;
    let (mut fw, mut fv) = (fx, fx);
    let (mut d, mut e): (Scalar, Scalar) = (0.0, 0.0);
    for _ in 0..MAX_ITERATIONS {
        let mid = 0.5 * (a + b);
        let tol = REL_TOL * x.abs() + ABS_TOL;
        if (x - mid).abs() <= 2.0 * tol - 0.5 * (b - a) {
            break;
        }
        let mut golden = true;
        if e.abs() > tol {
            let r = (x - w) * (fx - fv);
            let mut q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            }
            q = q.abs();
            let before = e;
            e = d;
            if !(p.abs() >= (0.5 * q * before).abs() || p <= q * (a - x) || p >= q * (b - x)) {
                d = p / q;
                let u = x + d;
                if u - a < 2.0 * tol || b - u < 2.0 * tol {
                    d = tol.copysign(mid - x);
                }
                golden = false;
            }
        }
        if golden {
            e = if x >= mid { a - x } else { b - x };
            d = GOLDEN * e;
        }
        let u = if d.abs() >= tol {
            x + d
        } else {
            x + tol.copysign(d)
        };
        let fu = f(u)?;
        if fu <= fx {
            if u >= x {
                a = x;
            } else {
                b = x;
            }
            (v, w, x) = (w, x, u);
            (fv, fw, fx) = (fw, fx, fu);
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || w == x {
                (v, w) = (w, u);
                (fv, fw) = (fw, fu);
            } else if fu <= fv || v == x || v == w {
                (v, fv) = (u, fu);
            }
        }
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::Support;
    use super::*;
    use axiolid_core::{Frame3, Interval, Point3};
    use axiolid_curve::{Curve2, Polyline2};
    use axiolid_surface::{Plane, Surface};

    fn plane() -> Surface {
        Surface::Plane(Plane {
            frame: Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
        })
    }

    /// The line from the origin to `(1, 0)`, on the plane that reads `(u, v)`
    /// as `(x, y)`.
    fn line() -> Curve2 {
        Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)],
            closed: false,
        })
    }

    fn pcurve(curve: &Curve2) -> Pcurve<'_> {
        Pcurve {
            curve,
            span: axiolid_reference::curve::domain2(curve),
            key: 0,
        }
    }

    fn targets(xs: &[Scalar]) -> Vec<Vec3> {
        xs.iter().map(|&x| Vec3::new(x, 0.0, 0.0)).collect()
    }

    fn locate(xs: &[Scalar], limit: Scalar) -> GeomResult<Located> {
        let curve = line();
        locate_samples(
            &pcurve(&curve),
            Support::Surface(&plane()),
            &targets(xs),
            Limits {
                interior: limit,
                end: limit,
            },
        )
    }

    #[test]
    fn samples_at_a_different_pace_are_found_where_they_lie() {
        // A quadratic pace: the uniform guess is far from every sample.
        let xs: Vec<Scalar> = (0..=8).map(|i| (i as Scalar / 8.0).powi(2)).collect();
        let located = locate(&xs, 1e-9).unwrap();
        for (p, x) in located.params.iter().zip(&xs) {
            assert!((p.x - x).abs() < 1e-9, "{p:?} for {x}");
        }
        assert!(located.ts.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn the_pcurve_own_samples_stay_on_their_uniform_parameters() {
        let xs: Vec<Scalar> = (0..=4).map(|i| i as Scalar / 4.0).collect();
        let located = locate(&xs, 1e-9).unwrap();
        assert_eq!(located.ts, [0.0, 0.25, 0.5, 0.75, 1.0]);
    }

    #[test]
    fn a_sample_off_the_surface_is_refused_above_the_limit() {
        let curve = line();
        let mut targets = targets(&[0.0, 0.25, 0.5, 0.75, 1.0]);
        targets[2].y = 3e-3;
        targets[1].y = 1e-3;
        let limits = Limits {
            interior: 1e-2,
            end: 1e-2,
        };
        locate_samples(
            &pcurve(&curve),
            Support::Surface(&plane()),
            &targets,
            limits,
        )
        .unwrap();
        let error = locate_samples(
            &pcurve(&curve),
            Support::Surface(&plane()),
            &targets,
            Limits {
                interior: 2e-3,
                end: 1e-2,
            },
        )
        .unwrap_err();
        assert!(
            matches!(&error, GeomError::Degenerate(text) if text == NOT_LOCATED),
            "{error:?}"
        );
    }

    #[test]
    fn an_end_may_sit_further_than_the_interior_by_its_own_limit() {
        let curve = line();
        let mut targets = targets(&[0.0, 0.5, 1.0]);
        targets[2].y = 5e-3;
        let at = |interior, end| {
            locate_samples(
                &pcurve(&curve),
                Support::Surface(&plane()),
                &targets,
                Limits { interior, end },
            )
        };
        assert!(at(1e-3, 1e-2).is_ok());
        assert!(at(1e-3, 1e-3).is_err());
        assert!(at(1e-2, 1e-3).is_err());
    }

    #[test]
    fn a_sample_before_its_predecessor_is_refused_not_clamped() {
        // The third sample lies behind the second; the generous limit would
        // admit the second's own parameter as its location.
        let error = locate(&[0.0, 0.5, 0.45, 1.0], 1.0).unwrap_err();
        assert!(
            matches!(&error, GeomError::Degenerate(text) if text == NOT_LOCATED),
            "{error:?}"
        );
    }

    #[test]
    fn a_sample_beyond_the_bracket_is_refused_not_found() {
        // The pcurve dwells at x = 0.5 for a hundred parameter units out of
        // 102: its pace there gives the samples left no guess of how far the
        // next one lies, and it is past every widening of the bracket.
        let mut points = vec![Point2::new(0.0, 0.0)];
        points.extend((0..=100).map(|_| Point2::new(0.5, 0.0)));
        points.push(Point2::new(1.0, 0.0));
        let curve = Curve2::Polyline(Polyline2 {
            points,
            closed: false,
        });
        let xs: Vec<Scalar> = (0..=200).map(|i| i as Scalar / 200.0).collect();
        let error = locate_samples(
            &pcurve(&curve),
            Support::Surface(&plane()),
            &targets(&xs),
            Limits {
                interior: 1.0,
                end: 1.0,
            },
        )
        .unwrap_err();
        assert!(
            matches!(&error, GeomError::Degenerate(text) if text == NOT_LOCATED),
            "{error:?}"
        );
    }

    #[test]
    fn a_pace_a_hundred_times_off_the_uniform_one_is_followed() {
        // The pcurve spends half its parameters on the first thousandth of
        // the edge: the first sample lies half the domain away from its
        // uniform parameter, far past what widening the bracket reaches, and
        // only the pace of the pcurve at the start tells where to look.
        let curve = Curve2::Polyline(Polyline2 {
            points: vec![
                Point2::new(0.0, 0.0),
                Point2::new(0.001, 0.0),
                Point2::new(1.0, 0.0),
            ],
            closed: false,
        });
        let xs: Vec<Scalar> = (0..=400).map(|i| i as Scalar / 400.0).collect();
        let located = locate_samples(
            &pcurve(&curve),
            Support::Surface(&plane()),
            &targets(&xs),
            Limits {
                interior: 1e-9,
                end: 1e-9,
            },
        )
        .unwrap();
        for (p, x) in located.params.iter().zip(&xs) {
            assert!((p.x - x).abs() < 1e-9, "{p:?} for {x}");
        }
    }

    #[test]
    fn a_reversed_span_is_located_from_its_start() {
        let curve = line();
        let reversed = Pcurve {
            curve: &curve,
            span: Interval {
                start: 1.0,
                end: 0.0,
            },
            key: 0,
        };
        let xs = [1.0, 0.6, 0.3, 0.0];
        let located = locate_samples(
            &reversed,
            Support::Surface(&plane()),
            &targets(&xs),
            Limits {
                interior: 1e-9,
                end: 1e-9,
            },
        )
        .unwrap();
        assert!(
            located.ts.windows(2).all(|w| w[0] > w[1]),
            "{:?}",
            located.ts
        );
        for (p, x) in located.params.iter().zip(xs) {
            assert!((p.x - x).abs() < 1e-9, "{p:?} for {x}");
        }
    }

    #[test]
    fn a_sample_at_a_corner_is_put_on_the_corner_exactly() {
        let curve = Curve2::Polyline(Polyline2 {
            points: vec![
                Point2::new(0.0, 0.0),
                Point2::new(0.37, 0.0),
                Point2::new(1.0, 0.0),
            ],
            closed: false,
        });
        let xs = [0.0, 0.1, 0.2, 0.3, 0.37, 0.5, 0.7, 0.9, 1.0];
        let located = locate_samples(
            &pcurve(&curve),
            Support::Surface(&plane()),
            &targets(&xs),
            Limits {
                interior: 1e-9,
                end: 1e-9,
            },
        )
        .unwrap();
        assert!(located.ts.contains(&1.0), "{:?}", located.ts);
    }

    #[test]
    fn a_pcurve_at_the_pace_of_the_samples_keeps_their_uniform_parameters() {
        // The samples come from another pcurve, at the same pace: their
        // images are a rounding off this one's, and the search finds
        // parameters a rounding off the uniform ones. The uniform ones are
        // kept, as the face across has them.
        let curve = Curve2::Polyline(Polyline2 {
            points: vec![
                Point2::new(0.0, 0.0),
                Point2::new(0.3, 0.0),
                Point2::new(1.0, 0.0),
            ],
            closed: false,
        });
        let n = 49;
        let xs: Vec<Scalar> = (0..=n).map(|i| i as Scalar / n as Scalar).collect();
        // A pace of its own: x(t) is piecewise linear in t, not uniform, so
        // the uniform parameters of the samples are not where they lie.
        let located = locate_samples(
            &pcurve(&curve),
            Support::Surface(&plane()),
            &targets(&xs),
            Limits {
                interior: 1e-9,
                end: 1e-9,
            },
        )
        .unwrap();
        assert_ne!(located.ts[10], 2.0 * 10.0 / 49.0);
        // The samples taken from the pcurve itself are its uniform ones.
        let own: Vec<Vec3> = (0..=n)
            .map(|i| {
                let t = 2.0 * i as Scalar / n as Scalar;
                let p = axiolid_reference::curve::evaluate2(&curve, t).unwrap();
                // Nudged by a rounding, as another pcurve's image would be.
                Vec3::new(p.x + 1e-16 * (i % 3) as Scalar, 0.0, 0.0)
            })
            .collect();
        let located = locate_samples(
            &pcurve(&curve),
            Support::Surface(&plane()),
            &own,
            Limits {
                interior: 1e-9,
                end: 1e-9,
            },
        )
        .unwrap();
        for (i, t) in located.ts.iter().enumerate() {
            assert_eq!(*t, 2.0 * i as Scalar / n as Scalar, "sample {i}");
        }
    }

    #[test]
    fn fewer_than_two_samples_have_no_edge() {
        assert!(locate(&[0.0], 1.0).is_err());
        assert!(locate(&[], 1.0).is_err());
    }
}
