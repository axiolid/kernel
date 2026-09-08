//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

pub mod kernel01;
pub mod kernel02;
pub mod kernel03;
pub mod kernel11;
pub mod kernel12;
use crate::csg::boolean03::kernel03::winding03;
use crate::csg::boolean03::kernel12::intersect12;
use crate::csg::common::{OpType, Vec3};
use crate::csg::manifold::Manifold;

pub struct Boolean03 {
    pub p1q2: Vec<[usize; 2]>,
    pub p2q1: Vec<[usize; 2]>,
    pub x12: Vec<i32>,
    pub x21: Vec<i32>,
    pub w03: Vec<i32>,
    pub w30: Vec<i32>,
    pub v12: Vec<Vec3>,
    pub v21: Vec<Vec3>,
}

pub fn boolean03(mp: &Manifold, mq: &Manifold, op: &OpType) -> Boolean03 {
    let e = if op == &OpType::Add { 1. } else { -1. };
    let mut p1q2 = vec![];
    let mut p2q1 = vec![];
    let x12;
    let v12;
    let x21;
    let v21;
    let w03;
    let w30;

    #[cfg(feature = "parallel")]
    {
        (((x12, v12), w03), ((x21, v21), w30)) = rayon::join(
            || {
                rayon::join(
                    || intersect12(mp, mq, &mut p1q2, e, true),
                    || winding03(mp, mq, e, true),
                )
            },
            || {
                rayon::join(
                    || intersect12(mp, mq, &mut p2q1, e, false),
                    || winding03(mp, mq, e, false),
                )
            },
        );
    }

    #[cfg(not(feature = "parallel"))]
    {
        ((x12, v12), w03) = (
            intersect12(mp, mq, &mut p1q2, e, true),
            winding03(mp, mq, e, true),
        );
        ((x21, v21), w30) = (
            intersect12(mp, mq, &mut p2q1, e, false),
            winding03(mp, mq, e, false),
        );
    }

    Boolean03 {
        p1q2,
        p2q1,
        x12,
        x21,
        w03,
        w30,
        v12,
        v21,
    }
}

/// Alternative to [`boolean03`] using [`kernel03::winding03_fast`] for the
/// winding-number classification.
///
/// Not run under the same `rayon::join` as `intersect12`: the fast path
/// needs `intersect12`'s `p1q2` before it can start, so this sequences
/// where `boolean03` parallelised. The two `fwd` directions (`p1q2` vs
/// `p2q1`) are still independent of each other and stay parallel.
pub fn boolean03_fast(mp: &Manifold, mq: &Manifold, op: &OpType) -> Boolean03 {
    let e = if op == &OpType::Add { 1. } else { -1. };
    let mut p1q2 = vec![];
    let mut p2q1 = vec![];

    #[cfg(feature = "parallel")]
    let (((x12, v12), w03), ((x21, v21), w30)) = rayon::join(
        || {
            let r = intersect12(mp, mq, &mut p1q2, e, true);
            let w = kernel03::winding03_fast(mp, mq, e, true, &p1q2);
            (r, w)
        },
        || {
            let r = intersect12(mp, mq, &mut p2q1, e, false);
            let w = kernel03::winding03_fast(mp, mq, e, false, &p2q1);
            (r, w)
        },
    );

    #[cfg(not(feature = "parallel"))]
    let (((x12, v12), w03), ((x21, v21), w30)) = {
        let a = {
            let r = intersect12(mp, mq, &mut p1q2, e, true);
            let w = kernel03::winding03_fast(mp, mq, e, true, &p1q2);
            (r, w)
        };
        let b = {
            let r = intersect12(mp, mq, &mut p2q1, e, false);
            let w = kernel03::winding03_fast(mp, mq, e, false, &p2q1);
            (r, w)
        };
        (a, b)
    };

    Boolean03 {
        p1q2,
        p2q1,
        x12,
        x21,
        w03,
        w30,
        v12,
        v21,
    }
}
