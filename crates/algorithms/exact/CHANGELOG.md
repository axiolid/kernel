# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- First release (ADR 0068, #154). Filtered exact arithmetic for
  constructions: `Interval` (outward-rounded `f64`, the fast tier),
  `Dyadic` (big-integer mantissa times a power of two, the exact tier),
  and `certify`, which runs one `SignExpr` in the first and falls back to
  the second only when the interval cannot decide.
- `Root2`, `sign_root`, `sign_two_roots`: signs and comparisons of
  `(a + b*sqrt(c)) / d`, including across different radicands, by squaring
  with case analysis. Exact zeros are reported as zeros.
- Constructions: `crossing_orientation` (which side of a line two lines
  cross), `line_circle_hits` (missed, tangent or two hits, with tangency
  decided exactly), `LineHit::{cmp_param, orientation}` and
  `compare_along`.
