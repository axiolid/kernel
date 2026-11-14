//! The view's affine map and the per-vertex view-space data.
//!
//! A model point `p` maps to `U = (p - o)·x`, `V = (p - o)·y`, depth
//! `W = (p - o)·z`: exact polynomials of the binary64 inputs. Every vertex
//! keeps its exact value (the exact tier) and an enclosure (bounding boxes
//! only; a decision never reads the enclosure directly). The image of an
//! orthographic vertex is `(U, V)`; that of a perspective vertex, the eye
//! being the origin, is `(U, V) / -W` (ADR 0107).

use axiolid_contracts::{BackendId, GeomError, Operation};
use axiolid_core::{Frame3, Point3};
use axiolid_exact::{Arith, Dyadic, Interval};
use axiolid_mesh_hidden_lines_contract::{HiddenLineError, HiddenLineMesh, View, ViewProjection};

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

/// The refusal of see-through triangles by a provider that does not honour
/// them, by name. Such a provider leaves `supports_occlusion` at its default,
/// so a registry does not reach it with them (ADR 0111); a direct call is
/// refused here.
pub(crate) fn require_opaque(
    meshes: &[HiddenLineMesh<'_>],
    backend: BackendId,
    operation: Operation,
) -> Result<(), HiddenLineError> {
    if !meshes.iter().any(|entry| entry.occlusion.has_see_through()) {
        return Ok(());
    }
    Err(GeomError::UnsupportedInput {
        backend,
        operation,
        input: "see-through triangles",
    }
    .into())
}

/// One vertex in view space.
///
/// The exact coordinates are homogeneous (ADR 0107): `[U, V, W, h]` stands
/// for the point `(U / h, V / h, W / h)`. The first three are always the view
/// coordinates themselves. An orthographic view has `h = 1`; a perspective
/// view has `h = -W`, so the point is its central projection onto the plane
/// `W = -1` and `U / h`, `V / h` are the image coordinates.
#[derive(Debug, Clone)]
pub(crate) struct Vert {
    /// Exact `[U, V, W, h]`.
    pub exact: [Dyadic; 4],
    /// Sound enclosures of the image coordinates and the depth: `[U, V, W]`
    /// of an orthographic vertex, `[U / h, V / h, W]` of a perspective one
    /// (the whole line for the image when `h` is not certainly positive).
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
    /// TODO(#387): the visible-face path reads the image coordinates as if
    /// `h` were one (P4). The assertion makes a perspective vertex fail in a
    /// debug build there, and the hidden-line path never calls this for one.
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

    /// The view-space vertex for a finite model point of a perspective view:
    /// `h = -W`, which is positive exactly where the point is in front of the
    /// eye.
    pub(crate) fn perspective(p: Point3, frame: &Frame3) -> Self {
        let [u, v, w] = project::<Dyadic>(p, frame);
        let h = w.neg();
        let (w_enc, h_enc) = (w.enclosure(), h.enclosure());
        let image = |x: &Dyadic| {
            if h_enc.lo() > 0.0 {
                x.enclosure().quotient(h_enc)
            } else {
                Interval::WHOLE
            }
        };
        let enc = [image(&u), image(&v), w_enc];
        Self {
            exact: [u, v, w, h],
            enc,
        }
    }

    /// The vertex of `p` in `view`.
    pub(crate) fn in_view(p: Point3, view: &View) -> Self {
        match view.projection {
            ViewProjection::Perspective => Self::perspective(p, &view.frame),
            _ => Self::new(p, &view.frame),
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
