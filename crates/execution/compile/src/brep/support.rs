//! The surface a B-rep face lies on, as the tessellator reads it.
//!
//! The tessellator never matches a surface itself. It asks a [`Support`]
//! for a point, its parameter periods, its family name and whether two
//! faces share it, and a [`SupportOracle`] for the derivative bounds over
//! a parameter box. A support is an atomic [`axiolid_surface::Surface`] or
//! a surface relation read through a [`RelationSurface`] (#411, ADR 0103);
//! the view keeps the tessellator closed over one `Copy` type, so a further
//! support lands as one more arm here.

use crate::surface_relation::{RelationOracle, RelationSurface};
use axiolid_contracts::GeomResult;
use axiolid_core::{Point3, Scalar};
use axiolid_model::SurfaceRelation;
use axiolid_reference::bound::{SurfaceBoundOracle, SurfaceBounds};
use axiolid_surface::Surface;

/// What a face lies on.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Support<'a> {
    /// An atomic surface value.
    Surface(&'a Surface),
    /// A surface relation, with the node it was read from.
    Relation {
        node: &'a SurfaceRelation,
        surface: RelationSurface<'a>,
    },
}

impl<'a> Support<'a> {
    /// The atomic surface, when the support is one.
    pub(crate) fn atomic(self) -> Option<&'a Surface> {
        match self {
            Self::Surface(surface) => Some(surface),
            Self::Relation { .. } => None,
        }
    }

    /// The support's point at `(u, v)`.
    pub(crate) fn evaluate(self, u: Scalar, v: Scalar) -> GeomResult<Point3> {
        match self {
            Self::Surface(surface) => axiolid_reference::surface::evaluate(surface, u, v),
            Self::Relation { surface, .. } => surface.evaluate(u, v),
        }
    }

    /// The support's parameter periods.
    pub(crate) fn periods(self) -> (Option<Scalar>, Option<Scalar>) {
        match self {
            Self::Surface(surface) => super::surface_periods(surface),
            Self::Relation { surface, .. } => surface.periods(),
        }
    }

    /// The support's family name, for a report.
    pub(crate) fn family(self) -> &'static str {
        match self {
            Self::Surface(surface) => super::surface_family(surface),
            Self::Relation { surface, .. } => surface.family(),
        }
    }

    /// Whether both are the same node, by identity.
    pub(crate) fn same(self, other: Self) -> bool {
        match (self, other) {
            (Self::Surface(a), Self::Surface(b)) => std::ptr::eq(a, b),
            (Self::Relation { node: a, .. }, Self::Relation { node: b, .. }) => std::ptr::eq(a, b),
            _ => false,
        }
    }
}

/// The derivative bounds of a [`Support`] over parameter boxes.
#[derive(Debug, Clone)]
pub(crate) enum SupportOracle<'a> {
    /// An atomic surface's oracle.
    Surface(SurfaceBoundOracle<'a>),
    /// A surface relation's oracle.
    Relation(RelationOracle<'a>),
}

impl<'a> SupportOracle<'a> {
    /// The oracle of `support`, or `None` when it has no derivative bound.
    pub(crate) fn new(support: Support<'a>) -> Option<Self> {
        match support {
            Support::Surface(surface) => SurfaceBoundOracle::new(surface).map(Self::Surface),
            Support::Relation { surface, .. } => Some(Self::Relation(RelationOracle::new(surface))),
        }
    }

    /// The bounds over the parameter box, or `None` when it has none.
    pub(crate) fn bounds(&self, u: (Scalar, Scalar), v: (Scalar, Scalar)) -> Option<SurfaceBounds> {
        match self {
            Self::Surface(oracle) => oracle.bounds(u, v),
            Self::Relation(oracle) => oracle.bounds(u, v),
        }
    }
}
