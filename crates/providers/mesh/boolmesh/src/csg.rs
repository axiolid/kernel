//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//--- This Source Code Form is subject to the terms of the Mozilla Public
//--- License v.2.0.

//! The mesh boolean, absorbed from `boolmesh` 0.1.9 (ADR 0047).
//!
//! Absorbed rather than depended upon so that the defects and the hot
//! paths are reachable: over 99% of a boolean's runtime is inside
//! [`compute_boolean`], which upstream exposes as a single opaque call.
//! See `docs/adr/0047-absorb-mesh-boolean.md` for the full reasoning and
//! the measurements behind it.
//!
//! # Provenance
//!
//! Upstream is <https://github.com/komietty/boolmesh> by Saki Komikado,
//! MPL-2.0, a from-scratch Rust implementation of the algorithm described
//! by Elalish's Manifold. Every file in this module keeps its original
//! copyright header. Axiolid is MPL-2.0 itself, so absorbing adds no
//! licence obligation the project has not already accepted -- but those
//! headers must survive refactoring.
//!
//! # What was not absorbed
//!
//! Upstream's `compose` module (cube, sphere, torus, cone, cylinder,
//! extrude, fractal) is deliberately absent. Axiolid has its own
//! primitives and profile extrusion; a second set would be a second
//! answer to one question.

pub(crate) mod boolean03;
pub(crate) mod boolean45;
pub(crate) mod common;
pub(crate) mod manifold;
pub(crate) mod simplification;
pub(crate) mod triangulation;

// Absorbed modules import these through the module root, exactly as they did
// through upstream's crate root. Keeping the re-export means the absorbed
// files need no import rewriting beyond `crate::` -> `crate::csg::`.
pub(crate) use common::*;
pub(crate) use manifold::*;

use boolean03::boolean03;
use boolean45::boolean45;
use manifold::cleanup_unused_verts;
use simplification::simplify_topology;
use triangulation::triangulate;

pub(crate) use common::{OpType, Vec3u};
pub(crate) use manifold::Manifold;

/// Boolean of two closed, oriented manifolds.
///
/// Preserved verbatim from upstream `boolmesh::compute_boolean` so that
/// absorption is behaviour-identical; the conformance and differential
/// suites are the gate on that claim. Optimisation is separate work.
pub(crate) fn compute_boolean(
    mp: &Manifold,
    mq: &Manifold,
    op: OpType,
) -> Result<Manifold, String> {
    let eps = mp.eps.max(mq.eps);
    let tol = mp.tol.max(mq.tol);

    let b03 = boolean03(mp, mq, &op);
    let mut b45 = boolean45(mp, mq, &b03, &op);
    let mut trg = triangulate(mp, mq, &b45, eps)?;

    simplify_topology(
        &mut trg.hs,
        &mut b45.ps,
        &mut trg.ns,
        &mut trg.rs,
        b45.nv_from_p,
        b45.nv_from_q,
        eps,
    );

    cleanup_unused_verts(&mut b45.ps, &mut trg.hs);

    Manifold::new_impl(
        b45.ps,
        trg.hs
            .chunks(3)
            .map(|hs| Vec3u::new(hs[0].tail, hs[1].tail, hs[2].tail))
            .collect(),
        Some(eps),
        Some(tol),
    )
}
