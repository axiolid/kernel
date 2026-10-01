#![forbid(unsafe_code)]

//! Certified geometric predicates: the shared exact-arithmetic substrate.
//!
//! A predicate answers a *sign* question, and a sign drives topology, so a
//! plausible answer is not good enough. Every public predicate here escalates
//! to exact arithmetic rather than comparing against a global epsilon.
//!
//! This package is deliberately narrow. It carries no curve, surface, mesh,
//! B-rep, provider, or execution dependency, so a consumer that needs only
//! certified signs — linear intersection, polygon orientation, NURBS root
//! isolation, topology classification — pays for arithmetic and nothing else.
//! That includes big integers: the exact tier here is expansion arithmetic
//! (plus a fixed-size dyadic fallback for `orient3d`), and exact
//! *constructions* over big integers live in `axiolid-exact` instead. The
//! broad `axiolid-reference` oracle re-exports these items unchanged
//! (ADR 0036).
//!
//! A `*_filter` function may return `Certified::Uncertain`; a public
//! predicate never does.

pub mod arithmetic;
pub mod expansion;
pub mod orient3;
mod orient3_dyadic;
pub mod orientation;
pub mod scene;
pub mod sphere;
pub mod static_filter;

pub use expansion::{two_diff, two_product, two_sum};
pub use orient3::{orient3d, orient3d_filter};
pub use orientation::{orient2d, orient2d_filter};
pub use sphere::{
    in_diametral_sphere, in_diametral_sphere_filter, incircle, incircle_filter, insphere,
    insphere_filter,
};
pub use static_filter::StaticFilter;
