//! Intrinsic (natural) plane curve data: curvature as a function of arc length.
//!
//! A plane curve is fixed up to rigid motion by its curvature law k(s). Anchoring
//! that law to a start frame fixes it absolutely. This is the *natural equation*
//! (Cesaro/Whewell) rather than a parametric map, and it is the honest way to
//! carry a clothoid: the spiral has no elementary parametric form, so a kernel
//! that only speaks parametrically has to approximate it before it has even
//! stored it.
//!
//! # Exact and symbolic, deliberately not evaluated
//!
//! Everything here is closed-form symbolic data and closed-form symbolic
//! operations on it: differentiate the law, negate it, measure the total turning
//! it accumulates. No module in this crate turns a law into points.
//!
//! That boundary is not laziness, it is the mathematics. Recovering position from
//! k(s) requires integrating the tangent angle and then integrating the tangent,
//! which for a linear law is the Fresnel integral -- not an elementary function.
//! Any point you get is a quadrature result carrying a tolerance, so it belongs
//! to an evaluator that can state that tolerance, never to the representation.
//! Storing the law exactly and refusing to fake points is what keeps a clothoid
//! a clothoid instead of a polyline that used to be one.

use axiolid_core::{Frame2, Scalar};

/// One sinusoidal term of a curvature law:
/// `amplitude * sin(angular_frequency * s + phase)`.
///
/// A term is deliberately not a law: it carries no mean. The constant part of
/// a composite law lives in its polynomial, so a given function has one
/// representation instead of many that differ only in where the mean was put.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Harmonic {
    /// Peak deviation contributed by this term.
    pub amplitude: Scalar,
    /// Radians of phase per unit arc length.
    pub angular_frequency: Scalar,
    /// Phase offset at `s = 0`, in radians.
    pub phase: Scalar,
}

/// Curvature as a closed-form function of arc length.
///
/// Sign follows the usual plane convention: positive curvature turns the
/// tangent counter-clockwise. `s` is arc length measured from the curve start,
/// so a law is meaningful on `[0, length]` of the curve that carries it.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum CurvatureLaw {
    /// `k(s) = curvature`.
    ///
    /// Zero is a straight line; any other value is a circular arc of radius
    /// `1 / curvature`.
    Constant {
        /// The constant curvature.
        curvature: Scalar,
    },
    /// `k(s) = coefficients[0] + coefficients[1] * s + coefficients[2] * s^2 + ...`
    ///
    /// Degree 1 is the clothoid (Euler spiral), whose curvature is linear in
    /// arc length; that is the transition spiral used between straight and
    /// circular track so lateral acceleration ramps linearly instead of
    /// stepping. Higher degrees cover the Bloss and cubic-parabola families.
    ///
    /// An empty coefficient list is the zero polynomial: a straight line.
    Polynomial {
        /// Coefficients in ascending powers of arc length.
        coefficients: Vec<Scalar>,
    },
    /// `k(s) = mean + amplitude * sin(angular_frequency * s + phase)`.
    ///
    /// The sinusoidal transition family (Klein, cosine ramps). Kept distinct
    /// from `Polynomial` because it is exactly representable this way and a
    /// truncated series would not be.
    Sinusoid {
        /// Curvature the oscillation is centred on.
        mean: Scalar,
        /// Peak deviation from `mean`.
        amplitude: Scalar,
        /// Radians of phase per unit arc length.
        angular_frequency: Scalar,
        /// Phase offset at `s = 0`, in radians.
        phase: Scalar,
    },
    /// `k(s) = sum c_i s^i + sum A_j sin(w_j s + p_j)`.
    ///
    /// A polynomial and any number of harmonic terms at once. Neither
    /// `Polynomial` nor `Sinusoid` can hold a law with both a secular trend and
    /// an oscillation, so a curve of that shape previously had to be refused or
    /// approximated; this variant stores it exactly.
    ///
    /// The shape is flat and additive rather than a recursive `Sum(Vec<Self>)`.
    /// A recursive sum would let the same function be written in unboundedly
    /// many ways, would make `is_constant` a search over arbitrary trees, and
    /// would admit nested sums that mean nothing extra. Flattening keeps one
    /// canonical slot per kind of term, keeps the family closed under
    /// differentiation and integration, and keeps the structural predicates a
    /// finite check over two lists.
    ///
    /// Empty `harmonics` is exactly the polynomial law; an empty polynomial with
    /// empty harmonics is the zero law. Both are legal, so a caller assembling
    /// terms never has to special-case the empty stage.
    Composite {
        /// Coefficients in ascending powers of arc length.
        polynomial: Vec<Scalar>,
        /// Additive sinusoidal terms.
        harmonics: Vec<Harmonic>,
    },
    /// Pieces laid end to end along arc length, each with its own law.
    ///
    /// `breaks` holds the INTERIOR seam positions in arc length from the
    /// curve start, so `laws.len() == breaks.len() + 1` and piece `i` spans
    /// `breaks[i - 1] .. breaks[i]`, the first starting at `0` and the last
    /// ending at the carrying curve's length.
    ///
    /// Each piece's law is written in its OWN arc length, restarting at zero
    /// at its seam, so a piece does not depend on where it sits and moving
    /// one never rewrites its coefficients.
    ///
    /// This variant exists because a piecewise profile genuinely cannot be
    /// decomposed into several `Intrinsic2` values. Every piece after the
    /// first would need an absolute start frame whose origin is the position
    /// at the seam, and that position is the non-elementary integral this
    /// crate refuses to compute. Holding the pieces in ONE curve keeps a
    /// single absolute frame at the start and anchors the interior purely by
    /// arc length, so no interior position is ever required.
    ///
    /// Unlike a summed law, pieces are disjoint and ordered: the seams are
    /// observable data, not a redundant re-encoding of one function. That is
    /// why nesting is meaningful here and was not for `Composite` -- a piece
    /// may itself be piecewise, expressing refinement.
    ///
    /// Mismatched lengths stay representable, as everywhere else in this
    /// crate; the operations report `None`/`false` rather than guessing.
    Piecewise {
        /// Interior seam positions in arc length, ascending.
        breaks: Vec<Scalar>,
        /// One law per piece; `laws.len() == breaks.len() + 1`.
        laws: Vec<CurvatureLaw>,
    },
}

impl CurvatureLaw {
    /// The zero law: a straight line.
    #[must_use]
    pub const fn straight() -> Self {
        Self::Constant { curvature: 0.0 }
    }

    /// A circular arc of the given signed curvature.
    #[must_use]
    pub const fn circular(curvature: Scalar) -> Self {
        Self::Constant { curvature }
    }

    /// A clothoid whose curvature runs from `start` to `end` over `length`.
    ///
    /// This is the transition-spiral constructor: the rate is derived rather
    /// than asked for, because the two endpoint curvatures and the length are
    /// what an alignment actually specifies.
    ///
    /// A non-finite or zero `length` cannot define a rate, so the result is a
    /// constant `start` law -- the honest degenerate answer, not a division by
    /// zero smuggled into a coefficient.
    #[must_use]
    pub fn clothoid(start: Scalar, end: Scalar, length: Scalar) -> Self {
        if !length.is_finite() || length == 0.0 {
            return Self::Constant { curvature: start };
        }
        Self::Polynomial {
            coefficients: vec![start, (end - start) / length],
        }
    }

    /// A transition whose linear ramp carries one full sine correction over
    /// its length: `k(s) = start + (d/L) s - (d / 2pi) sin(2 pi s / L)`, where
    /// `d = end - start`.
    ///
    /// The sine term removes the curvature-rate step a plain clothoid has at
    /// each end, so the rate starts and ends at zero instead of jumping. The
    /// mean rate is still `d / L`, so the total turning is unchanged from the
    /// clothoid's `(start + end) / 2 * L`.
    ///
    /// As with `clothoid`, a non-finite or zero `length` cannot define a rate,
    /// so the result degrades to a constant `start` law.
    #[must_use]
    pub fn sine_corrected_transition(start: Scalar, end: Scalar, length: Scalar) -> Self {
        if !length.is_finite() || length == 0.0 {
            return Self::Constant { curvature: start };
        }
        let delta = end - start;
        let turn = core::f64::consts::TAU;
        Self::Composite {
            polynomial: vec![start, delta / length],
            harmonics: vec![Harmonic {
                amplitude: -delta / turn,
                angular_frequency: turn / length,
                phase: 0.0,
            }],
        }
    }
    /// Pieces laid end to end, each carrying its own law.
    ///
    /// The seams are interior positions in ascending arc length; the caller
    /// supplies one more law than seam. A mismatch is storable and reported
    /// by `is_well_formed`, not rejected here.
    #[must_use]
    pub fn piecewise(breaks: Vec<Scalar>, laws: Vec<CurvatureLaw>) -> Self {
        Self::Piecewise { breaks, laws }
    }

    /// Whether the stored shape is internally consistent.
    ///
    /// Only `Piecewise` can be malformed: it carries two lists whose lengths
    /// must agree and seams that must ascend. Every other variant is
    /// well-formed by construction, so this is `true` for them.
    ///
    /// Structural, like the other predicates: it inspects stored data and
    /// never evaluates the law. Non-finite or descending seams are reported
    /// rather than silently sorted, because reordering would change which
    /// piece owns which arc length.
    #[must_use]
    pub fn is_well_formed(&self) -> bool {
        match self {
            Self::Piecewise { breaks, laws } => {
                // n pieces need n-1 interior seams. The empty law is the
                // one exception: zero pieces carry zero seams, and it is a
                // legitimate value (an alignment with nothing in it yet)
                // rather than a broken one, so it is well formed and turns
                // nothing.
                if laws.is_empty() {
                    return breaks.is_empty();
                }
                laws.len() == breaks.len() + 1
                    && breaks.iter().all(|b| b.is_finite())
                    && breaks.windows(2).all(|w| w[0] < w[1])
                    && laws.iter().all(Self::is_well_formed)
            }
            _ => true,
        }
    }
    /// Whether the law is identically zero, i.e. a straight line.
    ///
    /// Exact: this is a structural test on the stored coefficients, not a
    /// sampled one.
    #[must_use]
    pub fn is_straight(&self) -> bool {
        match self {
            Self::Constant { curvature } => *curvature == 0.0,
            Self::Polynomial { coefficients } => coefficients.iter().all(|c| *c == 0.0),
            Self::Sinusoid {
                mean,
                amplitude,
                angular_frequency,
                ..
            } => {
                // A zero frequency freezes the sine at its phase value, so the
                // law is constant but not necessarily zero; that case is only
                // straight when the frozen value cancels the mean, which needs
                // evaluation. Report false rather than guess.
                *mean == 0.0 && *amplitude == 0.0 && angular_frequency.is_finite()
            }
            Self::Composite {
                polynomial,
                harmonics,
            } => {
                // Every polynomial coefficient must vanish, and every harmonic
                // must contribute nothing. A harmonic contributes nothing only
                // when its amplitude is zero: a zero-frequency term freezes at
                // A*sin(p), which cancels only for particular phases, and
                // deciding that needs evaluation. Report false rather than
                // guess, matching the Sinusoid precedent.
                polynomial.iter().all(|c| *c == 0.0)
                    && harmonics
                        .iter()
                        .all(|h| h.amplitude == 0.0 && h.angular_frequency.is_finite())
            }
            // Straight overall exactly when every piece is straight. The
            // seams are irrelevant to this question: a union of zero-curvature
            // pieces is zero-curvature whatever the break positions.
            Self::Piecewise { laws, .. } => laws.iter().all(Self::is_straight),
        }
    }

    /// Whether the law is constant in arc length.
    #[must_use]
    pub fn is_constant(&self) -> bool {
        match self {
            Self::Constant { .. } => true,
            Self::Polynomial { coefficients } => coefficients.iter().skip(1).all(|c| *c == 0.0),
            Self::Sinusoid { amplitude, .. } => *amplitude == 0.0,
            Self::Composite {
                polynomial,
                harmonics,
            } => {
                // Constant in s: no polynomial term above degree 0 survives, and
                // no harmonic actually oscillates. A zero-frequency harmonic is
                // frozen at A*sin(p) and so IS constant, unlike the straightness
                // case where its value would also have to cancel.
                polynomial.iter().skip(1).all(|c| *c == 0.0)
                    && harmonics
                        .iter()
                        .all(|h| h.amplitude == 0.0 || h.angular_frequency == 0.0)
            }
            // Constant across the WHOLE curve needs every piece constant and
            // every piece equal to its neighbours: a staircase of differing
            // constants is piecewise-constant but not constant.
            //
            // Structural equality is the honest test available. Two pieces can
            // be equal in value while differing in form (`Constant { 0 }` versus
            // an empty `Polynomial`), and settling that needs evaluation, so
            // report false rather than guess -- the Sinusoid precedent.
            Self::Piecewise { laws, .. } => {
                self.is_well_formed()
                    && laws.iter().all(Self::is_constant)
                    && laws.windows(2).all(|w| w[0] == w[1])
            }
        }
    }

    /// The law's value, when it does not vary with arc length.
    ///
    /// Returns `None` for a law that varies, so a caller that needs a
    /// single number cannot silently read one off a varying law. The
    /// variants that `is_constant` accepts are exactly the ones answered
    /// here: a frozen harmonic contributes `A * sin(p)`, which is why a
    /// zero-FREQUENCY term is constant without being zero.
    #[must_use]
    pub fn constant_value(&self) -> Option<Scalar> {
        if !self.is_constant() {
            return None;
        }
        match self {
            Self::Constant { curvature } => Some(*curvature),
            Self::Polynomial { coefficients } => Some(coefficients.first().copied().unwrap_or(0.0)),
            Self::Sinusoid { mean, .. } => Some(*mean),
            Self::Composite {
                polynomial,
                harmonics,
            } => {
                let constant = polynomial.first().copied().unwrap_or(0.0);
                let frozen: Scalar = harmonics.iter().map(|h| h.amplitude * h.phase.sin()).sum();
                Some(constant + frozen)
            }
            // `is_constant` already required every piece equal, so the first
            // piece speaks for the whole law.
            Self::Piecewise { laws, .. } => laws.first().and_then(Self::constant_value),
        }
    }
    /// Interior seam positions strictly inside `(0, span)`, ascending.
    ///
    /// A piecewise law is only piecewise-smooth: its value can jump at a
    /// seam. Gauss-Legendre quadrature assumes the integrand is smooth
    /// across a panel, so a panel straddling a seam loses most of its
    /// accuracy -- measured at 3.0e-3 on a joined curve, against 1e-9
    /// once panels break at the seam. An integrator must therefore split
    /// its panels here rather than spreading them uniformly.
    ///
    /// Nested piecewise laws report their inner seams too, in absolute
    /// arc length from this law's origin, because a piece may itself be
    /// piecewise and its seams are just as discontinuous.
    #[must_use]
    pub fn seams_within(&self, span: Scalar) -> Vec<Scalar> {
        let mut out = Vec::new();
        self.collect_seams(0.0, span, &mut out);
        out.sort_by(Scalar::total_cmp);
        out.dedup();
        out
    }

    fn collect_seams(&self, origin: Scalar, span: Scalar, out: &mut Vec<Scalar>) {
        let Self::Piecewise { breaks, laws } = self else {
            return;
        };
        if laws.len() != breaks.len() + 1 {
            return;
        }
        let mut start = 0.0;
        for (index, piece) in laws.iter().enumerate() {
            let end = breaks.get(index).copied().unwrap_or(Scalar::INFINITY);
            let absolute = origin + start;
            if absolute > 0.0 && absolute < span {
                out.push(absolute);
            }
            piece.collect_seams(absolute, span, out);
            if !end.is_finite() || origin + end >= span {
                break;
            }
            start = end;
        }
    }
    /// The derivative law `dk/ds`, in closed form.
    ///
    /// Exact symbolic differentiation. The sharpness of a clothoid is the
    /// constant this returns.
    #[must_use]
    pub fn derivative(&self) -> Self {
        match self {
            Self::Constant { .. } => Self::Constant { curvature: 0.0 },
            Self::Polynomial { coefficients } => Self::Polynomial {
                coefficients: coefficients
                    .iter()
                    .enumerate()
                    .skip(1)
                    .map(|(power, c)| *c * power as Scalar)
                    .collect(),
            },
            // d/ds [m + A sin(w s + p)] = A w cos(w s + p)
            //                           = A w sin(w s + p + pi/2)
            // Folding the cosine back into a sine keeps the family closed under
            // differentiation, so no new variant is needed.
            Self::Sinusoid {
                amplitude,
                angular_frequency,
                phase,
                ..
            } => Self::Sinusoid {
                mean: 0.0,
                amplitude: amplitude * angular_frequency,
                angular_frequency: *angular_frequency,
                phase: phase + core::f64::consts::FRAC_PI_2,
            },
            // Differentiation is linear, so each part differentiates in place:
            // the polynomial drops a degree, and each harmonic keeps its
            // frequency while gaining a factor w and a quarter-turn of phase.
            // The result is again a Composite, so the family stays closed.
            Self::Composite {
                polynomial,
                harmonics,
            } => Self::Composite {
                polynomial: polynomial
                    .iter()
                    .enumerate()
                    .skip(1)
                    .map(|(power, c)| *c * power as Scalar)
                    .collect(),
                harmonics: harmonics
                    .iter()
                    .map(|h| Harmonic {
                        amplitude: h.amplitude * h.angular_frequency,
                        angular_frequency: h.angular_frequency,
                        phase: h.phase + core::f64::consts::FRAC_PI_2,
                    })
                    .collect(),
            },
            // Differentiate each piece in its own arc length. Seams are
            // untouched: a piece restarts at zero at its seam, so its
            // derivative is again a law in the same local parameter.
            //
            // dk/ds is generally DISCONTINUOUS at a seam. That is the correct
            // answer, not a defect: a profile assembled from pieces jumps
            // wherever the pieces disagree, and nothing here smooths it.
            Self::Piecewise { breaks, laws } => Self::Piecewise {
                breaks: breaks.clone(),
                laws: laws.iter().map(Self::derivative).collect(),
            },
        }
    }

    /// The same function re-written in a coordinate that starts at `a`.
    ///
    /// Returns the law `g` with `g(u) = self(a + u)`. This is what a TRIM
    /// needs: restricting a curve to `[a, b]` does not approximate the
    /// shape, it re-anchors the same law, so the trimmed curve is exactly
    /// the original one on that span (ADR 0062).
    ///
    /// The family is closed under this operation, which is why trimming is
    /// exact rather than a refit:
    ///
    /// - a polynomial shifts by the binomial expansion of `(a + u)^i`;
    /// - a sinusoid shifts purely in PHASE, `p -> p + w * a`, because the
    ///   amplitude and frequency do not depend on where the window starts;
    /// - a piecewise law drops the pieces that end before `a`, rebases the
    ///   one containing `a`, and keeps the rest with their seams moved back.
    ///
    /// Returns `None` when `a` is not finite, or when a piecewise law is
    /// malformed, since the shifted law would then be a guess.
    #[must_use]
    pub fn shifted(&self, a: Scalar) -> Option<Self> {
        if !a.is_finite() {
            return None;
        }
        match self {
            Self::Constant { curvature } => Some(Self::Constant {
                curvature: *curvature,
            }),
            Self::Polynomial { coefficients } => Some(Self::Polynomial {
                coefficients: shift_polynomial(coefficients, a),
            }),
            Self::Sinusoid {
                mean,
                amplitude,
                angular_frequency,
                phase,
            } => Some(Self::Sinusoid {
                mean: *mean,
                amplitude: *amplitude,
                angular_frequency: *angular_frequency,
                phase: phase + angular_frequency * a,
            }),
            Self::Composite {
                polynomial,
                harmonics,
            } => Some(Self::Composite {
                polynomial: shift_polynomial(polynomial, a),
                harmonics: harmonics
                    .iter()
                    .map(|h| Harmonic {
                        amplitude: h.amplitude,
                        angular_frequency: h.angular_frequency,
                        phase: h.phase + h.angular_frequency * a,
                    })
                    .collect(),
            }),
            Self::Piecewise { breaks, laws } => shift_piecewise(breaks, laws, a),
        }
    }

    /// The mirrored law `-k(s)`: the same curve reflected.
    #[must_use]
    pub fn reversed_orientation(&self) -> Self {
        match self {
            Self::Constant { curvature } => Self::Constant {
                curvature: -curvature,
            },
            Self::Polynomial { coefficients } => Self::Polynomial {
                coefficients: coefficients.iter().map(|c| -c).collect(),
            },
            Self::Sinusoid {
                mean,
                amplitude,
                angular_frequency,
                phase,
            } => Self::Sinusoid {
                mean: -mean,
                amplitude: -amplitude,
                angular_frequency: *angular_frequency,
                phase: *phase,
            },
            // Negation is linear too: negate every coefficient and every
            // amplitude. Frequencies and phases are untouched, so mirroring
            // twice returns the original law exactly.
            Self::Composite {
                polynomial,
                harmonics,
            } => Self::Composite {
                polynomial: polynomial.iter().map(|c| -c).collect(),
                harmonics: harmonics
                    .iter()
                    .map(|h| Harmonic {
                        amplitude: -h.amplitude,
                        angular_frequency: h.angular_frequency,
                        phase: h.phase,
                    })
                    .collect(),
            },
            // Mirroring negates curvature everywhere, so it negates each piece
            // in place. The seams are positions in arc length, not curvature
            // values, so they are unchanged -- this mirrors the curve about its
            // start tangent rather than reversing its direction of travel.
            Self::Piecewise { breaks, laws } => Self::Piecewise {
                breaks: breaks.clone(),
                laws: laws.iter().map(Self::reversed_orientation).collect(),
            },
        }
    }
}

/// Integral over `[0, s]` of `sum c_i x^i`, i.e. `sum c_i s^(i+1) / (i+1)`.
/// Turning accumulated by `law` over `span` arc length from its own start.
///
/// Closed form for every variant, including nested pieces: each piece is
/// integrated over its OWN subinterval and the results are added. Recursion
/// terminates because a piece's span is strictly shorter than its parent's.
///
/// Returns `None` when the shape cannot define an integral: a malformed
/// piecewise law, or seams that fall outside `[0, span]`. Reporting the
/// refusal beats inventing a clamp the caller did not ask for.
fn turning_over(law: &CurvatureLaw, span: Scalar) -> Option<Scalar> {
    match law {
        CurvatureLaw::Constant { curvature } => Some(curvature * span),
        CurvatureLaw::Polynomial { coefficients } => Some(polynomial_turning(coefficients, span)),
        CurvatureLaw::Sinusoid {
            mean,
            amplitude,
            angular_frequency,
            phase,
        } => Some(
            mean * span
                + harmonic_turning(
                    &Harmonic {
                        amplitude: *amplitude,
                        angular_frequency: *angular_frequency,
                        phase: *phase,
                    },
                    span,
                ),
        ),
        CurvatureLaw::Composite {
            polynomial,
            harmonics,
        } => Some(
            polynomial_turning(polynomial, span)
                + harmonics
                    .iter()
                    .map(|h| harmonic_turning(h, span))
                    .sum::<Scalar>(),
        ),
        CurvatureLaw::Piecewise { breaks, laws } => {
            if !law.is_well_formed() {
                return None;
            }
            // A seam at or before zero would mean the pieces do not tile the
            // span from its start, which is a malformed law rather than a
            // short window.
            if breaks.iter().any(|b| *b <= 0.0) {
                return None;
            }
            // Seams BEYOND the span are fine: integrating over `[0, span]`
            // simply stops inside whichever piece contains `span`, and the
            // later pieces are never reached. Rejecting them would make a
            // partial evaluation of a piecewise curve impossible, which is
            // exactly what trimming and mid-curve evaluation need.
            let mut total = 0.0;
            let mut start = 0.0;
            for (index, piece) in laws.iter().enumerate() {
                let end = breaks.get(index).copied().unwrap_or(span).min(span);
                if end <= start {
                    break;
                }
                total += turning_over(piece, end - start)?;
                start = end;
            }
            Some(total)
        }
    }
}
/// Coefficients of `p(a + u)` in ascending powers of `u`.
///
/// Binomial expansion: the `j`-th shifted coefficient collects
/// `c_i * C(i, j) * a^(i-j)` over every `i >= j`. Exact in the family:
/// a degree-`n` polynomial shifts to a degree-`n` polynomial.
fn shift_polynomial(coefficients: &[Scalar], a: Scalar) -> Vec<Scalar> {
    let n = coefficients.len();
    let mut out = vec![0.0; n];
    for (i, c) in coefficients.iter().enumerate() {
        // Pascal's triangle row `i`, built incrementally so no factorial
        // overflows and no combinatorial function is needed.
        let mut binomial = 1.0;
        for (j, slot) in out.iter_mut().enumerate().take(i + 1) {
            *slot += c * binomial * a.powi((i - j) as i32);
            // C(i, j+1) = C(i, j) * (i - j) / (j + 1)
            binomial = binomial * ((i - j) as Scalar) / ((j + 1) as Scalar);
        }
    }
    out
}

/// Restrict a piecewise law to the window starting at `a`.
///
/// Pieces wholly before `a` are dropped; the piece containing `a` is
/// rebased into its own coordinate; later pieces keep their laws and move
/// their seams back by `a`. Each piece is already written in its own arc
/// length, so only the piece straddling `a` is rewritten.
fn shift_piecewise(breaks: &[Scalar], laws: &[CurvatureLaw], a: Scalar) -> Option<CurvatureLaw> {
    if laws.len() != breaks.len() + 1 {
        return None;
    }
    if breaks.windows(2).any(|w| w[1] <= w[0]) || breaks.iter().any(|b| !b.is_finite()) {
        return None;
    }
    // Which piece contains `a`? Pieces are [0, b0), [b0, b1), ...
    let index = breaks.iter().take_while(|b| **b <= a).count();
    let piece_start = if index == 0 { 0.0 } else { breaks[index - 1] };
    let head = laws[index].shifted(a - piece_start)?;
    if index == breaks.len() {
        // `a` lies in the final piece: the window is that piece alone.
        return Some(head);
    }
    let mut kept = Vec::with_capacity(laws.len() - index);
    kept.push(head);
    kept.extend_from_slice(&laws[index + 1..]);
    let moved: Vec<Scalar> = breaks[index..].iter().map(|b| b - a).collect();
    Some(CurvatureLaw::Piecewise {
        breaks: moved,
        laws: kept,
    })
}
fn polynomial_turning(coefficients: &[Scalar], s: Scalar) -> Scalar {
    coefficients
        .iter()
        .enumerate()
        .map(|(power, c)| c * s.powi(power as i32 + 1) / (power as Scalar + 1.0))
        .sum()
}

/// Integral over `[0, s]` of `A sin(w x + p)`.
///
/// That antiderivative is `-(A/w)[cos(w s + p) - cos(p)]`, which is undefined
/// at `w == 0`. The integrand is then the constant `A sin(p)`, so the integral
/// is `A sin(p) s` -- the honest limit, not a special case invented to avoid a
/// division.
fn harmonic_turning(harmonic: &Harmonic, s: Scalar) -> Scalar {
    let Harmonic {
        amplitude,
        angular_frequency,
        phase,
    } = *harmonic;
    if angular_frequency == 0.0 {
        return amplitude * phase.sin() * s;
    }
    -(amplitude / angular_frequency) * ((angular_frequency * s + phase).cos() - phase.cos())
}
/// A plane curve given by its natural equation: a curvature law anchored to a
/// start frame and run for a finite arc length.
///
/// The frame supplies the rigid motion that k(s) alone cannot: its origin is the
/// curve start, and its `x` axis is the start tangent direction.
///
/// Like every other value in this crate, dirty imported data stays
/// representable -- a non-positive length is storable, and naming it is the
/// job of a validator, not of the type.
#[derive(Debug, Clone, PartialEq)]
pub struct Intrinsic2 {
    /// Start frame: origin at the curve start, `x` along the start tangent.
    pub start: Frame2,
    /// Curvature as a function of arc length from `start`.
    pub curvature: CurvatureLaw,
    /// Arc length the law is defined over.
    pub length: Scalar,
}

impl Intrinsic2 {
    /// Anchor a curvature law to a start frame over an arc length.
    #[must_use]
    pub const fn new(start: Frame2, curvature: CurvatureLaw, length: Scalar) -> Self {
        Self {
            start,
            curvature,
            length,
        }
    }

    /// Whether the curve is a straight segment.
    #[must_use]
    pub fn is_straight(&self) -> bool {
        self.curvature.is_straight()
    }

    /// Tangent heading at arc length `s` from the start, in radians.
    ///
    /// The integral of `k` over `[0, s]`, which every law in the family has in
    /// closed form -- so the heading is exact even though the POSITION is not
    /// and has to be quadratured. Measured from the start tangent, so the
    /// absolute heading is this plus the start frame's rotation.
    ///
    /// Returns `None` when `s` is not finite or a piecewise law does not tile
    /// `[0, s]`, matching [`Self::total_turning`].
    #[must_use]
    pub fn heading_at(&self, s: Scalar) -> Option<Scalar> {
        if !s.is_finite() {
            return None;
        }
        turning_over(&self.curvature, s)
    }

    /// Total tangent turning over the curve, in radians, in closed form.
    ///
    /// This is the integral of k(s) over `[0, length]` -- exact for every law in
    /// the family, because each one has an elementary antiderivative. It is the
    /// *angle* that integrates in closed form; position does not, which is why
    /// this method exists and an `evaluate` does not.
    ///
    /// Returns `None` when the length is not finite, since the integral is then
    /// undefined rather than merely large.
    #[must_use]
    pub fn total_turning(&self) -> Option<Scalar> {
        if !self.length.is_finite() {
            return None;
        }
        // Over the DECLARED length, the pieces must tile the whole domain: a
        // seam at or beyond the end means the stored law is not the one being
        // asked about, so refuse rather than clamp. This is stricter than
        // `heading_at`, which asks for a partial span and legitimately stops
        // inside whichever piece contains the sample.
        if let CurvatureLaw::Piecewise { breaks, .. } = &self.curvature {
            if breaks.iter().any(|b| *b >= self.length) {
                return None;
            }
        }
        turning_over(&self.curvature, self.length)
    }

    /// Upper bound on the total variation of heading over `[0, s]`, in radians.
    ///
    /// `total_turning` is the SIGNED integral of `k`; over a zero-mean
    /// oscillation it is zero however violently the curve wiggles. A quadrature
    /// budget derived from it would then spend one panel on a curve that needs
    /// many, so this returns the integral of `|k|` instead -- or a bound on it.
    ///
    /// Exact for `Constant`. For the other families an exact total variation
    /// needs the sign changes of `k`, so this returns a cheap upper bound from
    /// the triangle inequality: a bound is what a panel budget wants, and
    /// overestimating costs panels while underestimating costs correctness.
    ///
    /// `None` when `s` is not finite, matching `heading_at`.
    #[must_use]
    pub fn turning_variation_bound(&self, s: Scalar) -> Option<Scalar> {
        if !s.is_finite() {
            return None;
        }
        variation_bound(&self.curvature, s)
    }
}

/// Upper bound on the integral of `|k|` over `[0, span]`.
fn variation_bound(law: &CurvatureLaw, span: Scalar) -> Option<Scalar> {
    let span_abs = span.abs();
    match law {
        // Exact: |k| is constant, so the integral is just |k| * span.
        CurvatureLaw::Constant { curvature } => Some(curvature.abs() * span_abs),
        // sum |c_i| s^(i+1) / (i+1) >= |int sum c_i s^i|, term by term.
        CurvatureLaw::Polynomial { coefficients } => {
            Some(polynomial_variation(coefficients, span_abs))
        }
        // |mean| * span + |amplitude| * span bounds both parts; the harmonic
        // term integrates to at most its amplitude times the span.
        CurvatureLaw::Sinusoid {
            mean, amplitude, ..
        } => Some((mean.abs() + amplitude.abs()) * span_abs),
        CurvatureLaw::Composite {
            polynomial,
            harmonics,
        } => Some(
            polynomial_variation(polynomial, span_abs)
                + harmonics
                    .iter()
                    .map(|h| h.amplitude.abs() * span_abs)
                    .sum::<Scalar>(),
        ),
        // Pieces tile the span, so their variations add. Each piece is written
        // in its own arc length, exactly as `turning_over` treats them.
        CurvatureLaw::Piecewise { breaks, laws } => {
            if !law.is_well_formed() {
                return None;
            }
            if breaks.iter().any(|b| *b <= 0.0) {
                return None;
            }
            // As in `turning_over`: a seam past the span just means the later
            // pieces are not reached by this window.
            let mut total = 0.0;
            let mut start = 0.0;
            for (index, piece) in laws.iter().enumerate() {
                let end = breaks.get(index).copied().unwrap_or(span_abs).min(span_abs);
                if end <= start {
                    break;
                }
                total += variation_bound(piece, end - start)?;
                start = end;
            }
            Some(total)
        }
    }
}

/// Term-by-term bound on the integral of |polynomial| over `[0, span]`.
fn polynomial_variation(coefficients: &[Scalar], span: Scalar) -> Scalar {
    coefficients
        .iter()
        .enumerate()
        .map(|(power, c)| c.abs() * span.powi(power as i32 + 1) / (power as Scalar + 1.0))
        .sum()
}
