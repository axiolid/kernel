//! The view's affine map and the per-vertex view-space data.
//!
//! A model point `p` maps to `U = (p - o)·x`, `V = (p - o)·y`, depth
//! `W = (p - o)·z`: exact polynomials of the binary64 inputs. Every vertex
//! keeps its exact value (the exact tier) and an enclosure (bounding boxes
//! only; a decision never reads the enclosure directly).

use axiolid_contracts::{BackendId, GeomError, Operation};
use axiolid_core::{Frame3, Point3};
use axiolid_exact::{Arith, Dyadic, Interval};
use axiolid_mesh_hidden_lines_contract::{HiddenLineError, View, ViewProjection};

/// The refusal of a view this crate does not draw, by name. The providers
/// leave `supports_projection` at its default, so a registry does not reach
/// them with such a view (ADR 0107); a direct call is refused here.
pub(crate) fn require_orthographic(
    view: &View,
    backend: BackendId,
    operation: Operation,
) -> Result<(), HiddenLineError> {
    if view.projection == ViewProjection::Orthographic {
        return Ok(());
    }
    Err(GeomError::UnsupportedInput {
        backend,
        operation,
        input: "perspective view",
    }
    .into())
}

/// One vertex in view space.
///
/// The exact coordinates are homogeneous (ADR 0107): `[U, V, W, h]` stands
/// for the point `(U / h, V / h, W / h)`. An orthographic view has `h = 1`,
/// so the first three are the view coordinates themselves; a perspective view
/// will carry a depth-dependent `h`.
#[derive(Debug, Clone)]
pub(crate) struct Vert {
    /// Exact `[U, V, W, h]`.
    pub exact: [Dyadic; 4],
    /// Sound enclosures of `[U, V, W]`.
    pub enc: [Interval; 3],
}

/// The homogeneous weight of an orthographic vertex.
pub(crate) fn one() -> Dyadic {
    Dyadic::from_parts(1.into(), 0)
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
    /// Whether the weight `h` is exactly one, as in every orthographic view.
    pub(crate) fn is_affine(&self) -> bool {
        self.exact[3].exponent() == 0 && *self.exact[3].mantissa() == 1.into()
    }

    /// The coordinate `axis` (0 to 2) as an affine value `U`, `V` or `W`.
    ///
    /// TODO(#387): a perspective vertex has `h != 1`; every caller of this
    /// reads the coordinate as if it were, and has to use the homogeneous
    /// form instead. The assertion makes a missed site fail in a debug build.
    pub(crate) fn affine(&self, axis: usize) -> &Dyadic {
        debug_assert!(self.is_affine(), "a site that assumes h == 1 met h != 1");
        &self.exact[axis]
    }

    /// The view-space vertex for a finite model point.
    pub(crate) fn new(p: Point3, frame: &Frame3) -> Self {
        let [u, v, w] = project::<Dyadic>(p, frame);
        let enc = [u.enclosure(), v.enclosure(), w.enclosure()];
        Self {
            exact: [u, v, w, one()],
            enc,
        }
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
        // TODO(#387): `enc` encloses `U, V, W`; a perspective box needs the
        // enclosure of the image `U / h, V / h`.
        debug_assert!(ids.iter().all(|&id| verts[id as usize].is_affine()));
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
