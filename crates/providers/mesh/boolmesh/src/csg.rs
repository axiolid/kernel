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
use manifold::{cleanup_unused_verts, halfedges_are_two_manifold};
use simplification::simplify_topology;
use triangulation::triangulate;

pub(crate) use common::{OpType, Vec3u};
pub(crate) use manifold::Manifold;

/// The boolean's result as plain geometry.
///
/// `compute_boolean` used to end by feeding its result back through
/// `Manifold::new_impl`, which sorts faces into Morton order, builds a
/// half-edge mesh, a BVH, and a coplanar-face index, then validates
/// two-manifoldness. The provider reads none of that: `from_manifold`
/// takes positions and triangles and nothing else, and the provider
/// re-checks orientation itself.
///
/// Returning the geometry directly skips that whole rebuild. The
/// validity it used to assert is not lost -- see `compute_boolean`.
pub(crate) struct BooleanMesh {
    pub ps: Vec<Vec3>,
    pub tris: Vec<Vec3u>,
}

/// Boolean of two closed, oriented manifolds.
///
/// Returns plain geometry rather than a rebuilt [`Manifold`]. Upstream
/// ended by calling `Manifold::new_impl` on the result, which sorts the
/// faces into Morton order, rebuilds a half-edge mesh, builds a BVH and
/// a coplanar-face index, and validates two-manifoldness. The provider
/// consumes none of it: `from_manifold` read only positions and
/// triangles, so every one of those structures was discarded on the
/// next line.
///
/// Two behaviours of that call did matter and are kept explicitly:
///
/// - An empty result was signalled by `new_impl` failing on an empty
///   position matrix, which the provider matches on to return the empty
///   solid. The same error is raised directly here.
/// - `new_impl` rejected a non-two-manifold result. `simplify_topology`
///   can in principle leave one, so the check is retained -- but on the
///   half-edge data the boolean already has, rather than on a fresh
///   half-edge mesh built solely to ask the question.
pub(crate) fn compute_boolean(
    mp: &Manifold,
    mq: &Manifold,
    op: OpType,
) -> Result<BooleanMesh, String> {
    let eps = mp.eps.max(mq.eps);

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

    // Validate BEFORE cleanup: `cleanup_unused_verts` renumbers `tail` and
    // `head` but leaves `pair` addressing the pre-cleanup half-edge order,
    // so afterwards `pair` can point past the end of the array. Upstream
    // never hit that because it validated a freshly rebuilt half-edge mesh;
    // here the check has to happen while the indices are still coherent.
    if !halfedges_are_two_manifold(&trg.hs) {
        return Err("The input mesh is not manifold".into());
    }

    cleanup_unused_verts(&mut b45.ps, &mut trg.hs);

    // Preserves the signal the provider matches on: upstream reached this
    // through `edge_topology`, which refuses an empty position matrix.
    if b45.ps.is_empty() {
        return Err("empty pos matrix".into());
    }
    if trg.hs.is_empty() {
        return Err("empty idx matrix".into());
    }

    Ok(BooleanMesh {
        tris: trg
            .hs
            .chunks(3)
            .map(|hs| Vec3u::new(hs[0].tail, hs[1].tail, hs[2].tail))
            .collect(),
        ps: b45.ps,
    })
}
