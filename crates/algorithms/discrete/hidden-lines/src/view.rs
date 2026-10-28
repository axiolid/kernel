//! The view's affine map and the per-vertex view-space data.
//!
//! A model point `p` maps to `U = (p - o)·x`, `V = (p - o)·y`, depth
//! `W = (p - o)·z`: exact polynomials of the binary64 inputs. Every vertex
//! keeps its exact value (the exact tier) and an enclosure (bounding boxes
//! only; a decision never reads the enclosure directly).

use axiolid_core::{Frame3, Point3};
use axiolid_exact::{Arith, Dyadic, Interval};

/// One vertex in view space.
#[derive(Debug, Clone)]
pub(crate) struct Vert {
    /// Exact `[U, V, W]`.
    pub exact: [Dyadic; 3],
    /// Sound enclosures of `[U, V, W]`.
    pub enc: [Interval; 3],
}

/// `[U, V, W]` of `p` evaluated in arithmetic `T`.
pub(crate) fn project<T: Arith>(p: Point3, frame: &Frame3) -> [T; 3] {
    let d = [
        T::from_f64(p.x).sub(&T::from_f64(frame.origin.x)),
        T::from_f64(p.y).sub(&T::from_f64(frame.origin.y)),
        T::from_f64(p.z).sub(&T::from_f64(frame.origin.z)),
    ];
    let along = |axis: axiolid_core::Vec3| {
        d[0].mul(&T::from_f64(axis.x))
            .add(&d[1].mul(&T::from_f64(axis.y)))
            .add(&d[2].mul(&T::from_f64(axis.z)))
    };
    [along(frame.x), along(frame.y), along(frame.z)]
}

impl Vert {
    /// The view-space vertex for a finite model point.
    pub(crate) fn new(p: Point3, frame: &Frame3) -> Self {
        let exact = project::<Dyadic>(p, frame);
        let enc = [
            exact[0].enclosure(),
            exact[1].enclosure(),
            exact[2].enclosure(),
        ];
        Self { exact, enc }
    }
}

/// A closed box in view space, enclosing a set of vertices.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Bounds {
    pub lo: [f64; 3],
    pub hi: [f64; 3],
}

impl Bounds {
    /// The box enclosing `ids`; `ids` is never empty.
    pub(crate) fn of(verts: &[Vert], ids: &[u32]) -> Self {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &id in ids {
            for axis in 0..3 {
                let e = verts[id as usize].enc[axis];
                lo[axis] = lo[axis].min(e.lo());
                hi[axis] = hi[axis].max(e.hi());
            }
        }
        Self { lo, hi }
    }

    /// Whether the image rectangles share a point (closed).
    pub(crate) fn image_overlaps(&self, other: &Self) -> bool {
        self.lo[0] <= other.hi[0]
            && other.lo[0] <= self.hi[0]
            && self.lo[1] <= other.hi[1]
            && other.lo[1] <= self.hi[1]
    }
}
