//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

use super::kernel02::Kernel02;
use crate::csg::bounds::BPos;
use crate::csg::{Manifold, Real, Vec2};

pub fn winding03(mp: &Manifold, mq: &Manifold, expand: Real, fwd: bool) -> Vec<i32> {
    let ma = if fwd { mp } else { mq };
    let mb = if fwd { mq } else { mp };

    let mut w03 = vec![0; ma.nv];
    let k02 = Kernel02 {
        ps_p: &ma.ps,
        ps_q: &mb.ps,
        hs_q: &mb.hs,
        ns: &mp.vert_normals,
        expand,
        fwd,
    };

    // The xy grid: this query never reads z (see `BPos::overlaps_node`),
    // and is run once per operand, so an O(n) structure beats a tree here
    // -- there is no second call to amortize a tree build against.
    mb.planar_grid.collision(
        &ma.ps
            .iter()
            .enumerate()
            .map(|(i, p)| BPos {
                id: Some(i),
                pos: Vec2::new(p.x, p.y),
            })
            .collect::<Vec<_>>(),
        &mut |a, b| {
            if let Some((s, _)) = k02.op(a, b) {
                w03[a] += s * if fwd { 1 } else { -1 };
            }
        },
    );

    w03
}

fn find(parent: &mut [usize], mut node: usize) -> usize {
    while parent[node] != node {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    node
}

fn union(parent: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (find(parent, a), find(parent, b));
    if ra != rb {
        parent[rb] = ra;
    }
}

/// Alternative to [`winding03`]: one query per connected component of
/// unbroken edges, not one per vertex.
///
/// An edge of `a` whose endpoints are both untouched by any intersection
/// (absent from `p1q2`) cannot cross `b`'s boundary, so its two endpoints
/// share the same winding number. Union-finding on that fact, then
/// querying once per component and flood-filling, gives the same answer
/// as [`winding03`] with far fewer queries when few components exist --
/// which is the common case away from a thin intersection band.
///
/// `p1q2` must be the exact set this call's own [`super::intersect12`]
/// produced for the same `mp`/`mq`/`fwd`: both classifiers share the same
/// `Kernel02` predicate, so the two are consistent by construction, but a
/// `p1q2` from a different tolerance or a different mesh pair breaks that
/// invariant silently.
///
/// # Trade-off against [`winding03`]
///
/// Not an approximation -- the flood-filled answer is provably identical
/// to the per-vertex one, given a correct `p1q2`. What differs is fault
/// containment: a bug in edge-break detection mislabels an entire
/// component here, where the same bug under [`winding03`] would corrupt
/// only the vertices it directly touches. Opt-in, not a default, until
/// this has run against a wider correctness corpus than sphere unions.
pub fn winding03_fast(
    mp: &Manifold,
    mq: &Manifold,
    expand: Real,
    fwd: bool,
    p1q2: &[[usize; 2]],
) -> Vec<i32> {
    let ma = if fwd { mp } else { mq };
    let mb = if fwd { mq } else { mp };

    // p1q2's column layout depends on fwd: intersect12(fwd=true) pushes
    // [mp_edge, mq_face], but intersect12(fwd=false) pushes
    // [mp_face, mq_edge] -- see kernel12.rs's `rec` closure, which
    // swaps to [b, a] in the fwd=false branch. `ma`'s edge id (the one
    // this function needs to mark "broken") is therefore column 0 when
    // fwd, column 1 when not. Column selection here MUST track the same
    // fwd this function was called with, or every edge is misread.
    let col = if fwd { 0 } else { 1 };
    let mut broken: Vec<usize> = p1q2.iter().map(|pair| pair[col]).collect();
    broken.sort_unstable();
    broken.dedup();

    let mut parent: Vec<usize> = (0..ma.nv).collect();
    for (i, h) in ma.hs.iter().enumerate() {
        if !h.is_forward() {
            continue;
        }
        if broken.binary_search(&i).is_ok() {
            continue;
        }
        union(&mut parent, h.tail, h.head);
    }

    let k02 = Kernel02 {
        ps_p: &ma.ps,
        ps_q: &mb.ps,
        hs_q: &mb.hs,
        ns: &mp.vert_normals,
        expand,
        fwd,
    };

    // One representative vertex per component: its root, first-seen.
    let mut reps: Vec<usize> = Vec::new();
    let mut root_of_rep: Vec<usize> = vec![usize::MAX; ma.nv];
    for v in 0..ma.nv {
        let r = find(&mut parent, v);
        if root_of_rep[r] == usize::MAX {
            root_of_rep[r] = reps.len();
            reps.push(v);
        }
    }

    let mut w03_rep = vec![0i32; reps.len()];
    let query = reps
        .iter()
        .enumerate()
        .map(|(i, &v)| BPos {
            id: Some(i),
            pos: Vec2::new(ma.ps[v].x, ma.ps[v].y),
        })
        .collect::<Vec<_>>();
    mb.planar_grid.collision(&query, &mut |i, b| {
        if let Some((s, _)) = k02.op(reps[i], b) {
            w03_rep[i] += s * if fwd { 1 } else { -1 };
        }
    });

    let mut w03 = vec![0; ma.nv];
    for (v, slot) in w03.iter_mut().enumerate() {
        let r = find(&mut parent, v);
        *slot = w03_rep[root_of_rep[r]];
    }
    w03
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::csg::boolean03::kernel12::intersect12;

    /// A crude icosphere: enough facets for a real intersection band, no
    /// external dependency.
    fn icosphere(cx: f64, cy: f64, cz: f64, r: f64, subdiv: u32) -> (Vec<f64>, Vec<usize>) {
        let t = (1.0 + 5.0f64.sqrt()) / 2.0;
        let mut verts: Vec<[f64; 3]> = vec![
            [-1.0, t, 0.0],
            [1.0, t, 0.0],
            [-1.0, -t, 0.0],
            [1.0, -t, 0.0],
            [0.0, -1.0, t],
            [0.0, 1.0, t],
            [0.0, -1.0, -t],
            [0.0, 1.0, -t],
            [t, 0.0, -1.0],
            [t, 0.0, 1.0],
            [-t, 0.0, -1.0],
            [-t, 0.0, 1.0],
        ];
        let mut faces: Vec<[u32; 3]> = vec![
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ];
        for _ in 0..subdiv {
            let mut mid: std::collections::HashMap<(u32, u32), u32> =
                std::collections::HashMap::new();
            let mut next: Vec<[u32; 3]> = Vec::with_capacity(faces.len() * 4);
            for f in &faces {
                let mut m = [0u32; 3];
                for e in 0..3 {
                    let (a, b) = (f[e], f[(e + 1) % 3]);
                    let key = (a.min(b), a.max(b));
                    m[e] = *mid.entry(key).or_insert_with(|| {
                        let (pa, pb) = (verts[a as usize], verts[b as usize]);
                        verts.push([
                            (pa[0] + pb[0]) * 0.5,
                            (pa[1] + pb[1]) * 0.5,
                            (pa[2] + pb[2]) * 0.5,
                        ]);
                        (verts.len() - 1) as u32
                    });
                }
                next.push([f[0], m[0], m[2]]);
                next.push([f[1], m[1], m[0]]);
                next.push([f[2], m[2], m[1]]);
                next.push([m[0], m[1], m[2]]);
            }
            faces = next;
        }
        let mut pos = Vec::with_capacity(verts.len() * 3);
        for v in &verts {
            let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            let k = r / len;
            pos.push(cx + v[0] * k);
            pos.push(cy + v[1] * k);
            pos.push(cz + v[2] * k);
        }
        let idx = faces
            .iter()
            .flat_map(|f| [f[0] as usize, f[1] as usize, f[2] as usize])
            .collect();
        (pos, idx)
    }

    fn assert_winding_fast_matches_slow(
        mp: &Manifold,
        mq: &Manifold,
        expand: Real,
        fwd: bool,
        case: &str,
    ) {
        let mut p1q2 = vec![];
        intersect12(mp, mq, &mut p1q2, expand, fwd);
        let slow = winding03(mp, mq, expand, fwd);
        let fast = winding03_fast(mp, mq, expand, fwd, &p1q2);
        assert_eq!(slow, fast, "{case}: fwd={fwd}");
    }

    #[test]
    fn matches_slow_path_on_overlapping_spheres() {
        let (pa, ia) = icosphere(0.0, 0.0, 0.0, 1.0, 3);
        let (pb, ib) = icosphere(1.2, 0.0, 0.0, 1.0, 3);
        let mp = Manifold::new(&pa, &ia).expect("sphere a");
        let mq = Manifold::new(&pb, &ib).expect("sphere b");
        for fwd in [true, false] {
            assert_winding_fast_matches_slow(&mp, &mq, 1.0, fwd, "overlapping spheres");
        }
    }

    #[test]
    fn matches_slow_path_on_deeply_nested_spheres() {
        // b almost entirely inside a: exercises the "whole mesh minus a
        // tiny band is one component" case the optimisation targets.
        let (pa, ia) = icosphere(0.0, 0.0, 0.0, 2.0, 3);
        let (pb, ib) = icosphere(0.0, 0.0, 0.0, 1.0, 3);
        let mp = Manifold::new(&pa, &ia).expect("outer sphere");
        let mq = Manifold::new(&pb, &ib).expect("inner sphere");
        for fwd in [true, false] {
            assert_winding_fast_matches_slow(&mp, &mq, 1.0, fwd, "nested spheres");
        }
    }

    #[test]
    fn matches_slow_path_on_barely_touching_spheres() {
        // Radius 1 + radius 1, centres 1.999 apart: a hair from tangency,
        // the thinnest possible intersection band -- most vertices are far
        // from it, few are near, stressing the epsilon boundary between
        // "broken" and "not broken".
        let (pa, ia) = icosphere(0.0, 0.0, 0.0, 1.0, 3);
        let (pb, ib) = icosphere(1.999, 0.0, 0.0, 1.0, 3);
        let mp = Manifold::new(&pa, &ia).expect("sphere a");
        let mq = Manifold::new(&pb, &ib).expect("sphere b");
        for fwd in [true, false] {
            assert_winding_fast_matches_slow(&mp, &mq, 1.0, fwd, "barely touching spheres");
        }
    }

    #[test]
    fn matches_slow_path_at_different_subdivisions_on_each_side() {
        // Asymmetric density: b's edges are much larger than a's, so a's
        // component near the cut can straddle several of b's faces --
        // exercising a union-find over a's edges against a broad phase
        // built from a differently-scaled b.
        let (pa, ia) = icosphere(0.0, 0.0, 0.0, 1.0, 4);
        let (pb, ib) = icosphere(1.2, 0.0, 0.0, 1.0, 1);
        let mp = Manifold::new(&pa, &ia).expect("dense sphere a");
        let mq = Manifold::new(&pb, &ib).expect("coarse sphere b");
        for fwd in [true, false] {
            assert_winding_fast_matches_slow(&mp, &mq, 1.0, fwd, "asymmetric subdivision");
        }
    }

    #[test]
    fn matches_slow_path_when_disjoint() {
        // No intersection at all: p1q2 is empty, so every vertex of a is
        // one giant component. This is the degenerate case the whole
        // optimisation is built around -- must not panic on an empty
        // broken-edge set.
        let (pa, ia) = icosphere(0.0, 0.0, 0.0, 1.0, 2);
        let (pb, ib) = icosphere(10.0, 0.0, 0.0, 1.0, 2);
        let mp = Manifold::new(&pa, &ia).expect("sphere a");
        let mq = Manifold::new(&pb, &ib).expect("sphere b");
        for fwd in [true, false] {
            assert_winding_fast_matches_slow(&mp, &mq, 1.0, fwd, "disjoint spheres");
        }
    }
}
