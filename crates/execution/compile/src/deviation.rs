//! How far the exact geometry may lie from a compiled mesh (#232).
//!
//! [`ReferenceMeshCompiler::compile_mesh_with_deviation`] returns, next to
//! the mesh, a [`DeviationReport`]: a certified upper bound on the distance
//! from any point of the exact surface the graph describes to the nearest
//! mesh triangle, which paths contributed to it, and whether it is within
//! the chord budget the caller asked for. This is the quantity #231 proves
//! for revolutions, spheres, tori and sweeps along one arc; for the paths
//! that have no proof by construction it is computed for the mesh at hand.
//!
//! # A bound, never an estimate
//!
//! Every number in a report is an upper bound: either a construction's
//! proof ([`DeviationBound::Proven`]) or a computation that only uses
//! certified derivative bounds and exact distances
//! ([`DeviationBound::Certified`]). Where neither exists the path is named
//! [`DeviationBound::Unbounded`] with its reason, and the report's overall
//! bound is `None`; nothing is filled in with a sampled or heuristic value.
//! The bounds are one-sided, exact surface to mesh: a mesh may carry extra
//! area (a cap the exact solid does not have) without changing them.
//!
//! # How each path is bounded
//!
//! - Authored meshes are their own exact surface: `0`.
//! - Profiles ([`axiolid_construct::profile::profile_deviation`]): the
//!   chord budget for every family the flattener certifies (lines, arcs,
//!   circles, ellipses, B-splines), plus any merged near-duplicate point,
//!   scaled by a derived profile's stretch. An extrusion's walls and caps
//!   are within that of the exact solid.
//! - Revolutions: the profile's bound plus the half budget the turn gets
//!   (#231). Tapered revolutions only between straight-edged profiles,
//!   whose rings correspond exactly; tapered extrusions are unbounded
//!   (their twisted walls are not measured).
//! - Swept disks: along lines, polylines (filleted) and chains of
//!   segments and arcs by `axiolid_construct::pipe`'s proof (#232), along
//!   one circular arc by #231's; along an ellipse or a B-spline, certified against the exact tube
//!   (walls, bore and square end caps) by branch and bound, below.
//!   Fixed-reference sweeps along a segment, or along a circle whose plane
//!   the reference is normal to, by #231; other frame laws, polyline and
//!   composite directrices are unbounded here.
//! - Primitives: spheres and tori by #231, cylinders and cones by #232
//!   (`axiolid_reference::primitive`), flat-faced ones exactly.
//! - B-rep faces: planar faces with straight edges exactly; curved faces
//!   per triangle from the surface's second derivatives over the
//!   triangle, plus the trim lens of each pcurve chord (the B-rep
//!   tessellator's module notes give the argument).
//!
//! # Branch and bound against an exact parameterisation
//!
//! Where the mesh's vertices sit at no known surface parameters (a tube's
//! stations stand on frames read off its sampled path, a boolean's cut
//! curve wherever the mesh boolean put it), the exact surface
//! `T` is covered by parameter cells. A cell with centre `x_c` and
//! half-widths `h` maps into `T(x_c) + DT(x_c) [-h, h]` grown by
//! `e2 = 1/2 (A h_x^2 + 2 B h_x h_y + C h_y^2)`, with `A, B, C` certified
//! bounds on `T`'s second partials over the cell (for a tube, from the
//! directrix's derivative bounds through its normalised frame). Distance
//! to one triangle is convex, so over a piece of that parallelogram it is
//! at most the largest of the piece's corners'; the best triangle per piece
//! and the worst piece, plus `e2`, bound the cell. Cells are split worst
//! first until the worst bound is within 10% of the largest exact distance
//! sampled at a cell centre, or a work cap is reached; the reported value
//! is always the worst bound over all cells.
//! - Instances scale the bound by their transform's largest stretch;
//!   collections take the worst member.
//! - Booleans (#235): cutting two meshes moves the intersection curve, and
//!   a one-sided bound on each operand does not bound it, so nothing is
//!   derived from the operands. Where `ReferenceExactCompiler` builds the
//!   exact result of the same node (differences of placed extrusions,
//!   #228, and half-space clips, #234), the mesh is certified against that
//!   exact B-rep by the branch and bound above, over each face's trimmed
//!   parameter domain; a boolean it refuses is unbounded with the
//!   refusal's name. The bound is relative to the exact compiler's result,
//!   the exact boolean of operands moved by at most the tolerance (#228),
//!   not to the boolean of the unperturbed operands. It runs only for a
//!   report, and only for the booleans whose result is emitted, not for
//!   the inner booleans of a chain; the search stops once within the
//!   requested budget.
//!
//! [`ReferenceMeshCompiler::compile_mesh_with_deviation`]: crate::ReferenceMeshCompiler::compile_mesh_with_deviation

use axiolid_core::{Scalar, Transform3};

mod boolean;
mod paths;
pub(crate) use boolean::{emitted_booleans, of_boolean};
pub(crate) use paths::{of_curve_bounded, of_primitive, of_sectioned_surface, of_solid};

/// The construction path a part of the mesh came from.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DeviationPath {
    /// An authored triangle or polygon mesh.
    AuthoredMesh,
    /// A linear extrusion of a profile.
    Extrusion,
    /// A revolution of a profile.
    Revolution,
    /// An extrusion between two profiles.
    TaperedExtrusion,
    /// A revolution between two profiles.
    TaperedRevolution,
    /// A disk swept along a directrix.
    SweptDisk,
    /// A profile swept with a fixed reference direction.
    FixedReferenceSweep,
    /// A profile swept along a curve on a surface.
    SurfaceCurveSweep,
    /// Sections placed along a spine.
    SectionedSpine,
    /// Sections standing at stations along a directrix (#241).
    StationedSpine,
    /// Open sections at stations joined by tag (#241).
    SectionedSurface,
    /// A half-space bounded by a polygon.
    BoundedHalfSpace,
    /// A CSG primitive.
    Primitive,
    /// A face of a B-rep.
    BRepFace,
    /// A plane trimmed by boundary curves.
    CurveBoundedPlane,
    /// A boolean of two meshes.
    Boolean,
    /// A path that does not report its deviation.
    Unreported,
}

/// A bound on how far the exact surface may lie from the mesh.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DeviationBound {
    /// Guaranteed by the construction's own proof, without inspecting the
    /// mesh (#231): typically the chord budget itself.
    Proven(Scalar),
    /// Computed for this mesh from certified derivative bounds and exact
    /// distances (#232).
    Certified(Scalar),
    /// No certified bound exists for this path; the reason names it.
    Unbounded(&'static str),
}

impl DeviationBound {
    /// The bound's value, `None` when unbounded.
    #[must_use]
    pub const fn value(self) -> Option<Scalar> {
        match self {
            Self::Proven(value) | Self::Certified(value) => Some(value),
            Self::Unbounded(_) => None,
        }
    }

    /// The worse of two bounds: unbounded wins, else the larger value
    /// (keeping its kind).
    #[must_use]
    pub fn worst(self, other: Self) -> Self {
        match (self.value(), other.value()) {
            (None, _) => self,
            (_, None) => other,
            (Some(a), Some(b)) => {
                if b > a {
                    other
                } else {
                    self
                }
            }
        }
    }

    fn scaled(self, factor: Scalar) -> Self {
        match self {
            Self::Proven(value) => Self::Proven(value * factor),
            Self::Certified(value) => Self::Certified(value * factor),
            unbounded => unbounded,
        }
    }
}

/// One path's share of a report.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeviationContribution {
    /// The construction path.
    pub path: DeviationPath,
    /// What within the path: a surface family, a directrix family, a
    /// primitive kind. Empty when the path says it all.
    pub detail: &'static str,
    /// The worst bound over every part of the mesh this path made.
    pub bound: DeviationBound,
}

/// A certified bound on the distance from the exact surface to a compiled
/// mesh (#232).
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct DeviationReport {
    /// The chord budget the compilation was asked for, in world units:
    /// `ExecutionOptions::chord_error`, else the linear tolerance.
    pub requested: Scalar,
    /// An upper bound on the distance from every point of the exact
    /// surface to the mesh's triangles, in world units; `None` when any
    /// contribution is unbounded.
    pub bound: Option<Scalar>,
    /// Each path that contributed, with the worst bound it reached, in
    /// first-seen order.
    pub contributions: Vec<DeviationContribution>,
}

impl DeviationReport {
    /// Whether the mesh is certified within the requested chord budget.
    ///
    /// `false` means "not certified", not "missed": an unbounded path, or a
    /// bound above the budget whose true deviation may still be below it.
    #[must_use]
    pub fn meets_requested(&self) -> bool {
        self.bound.is_some_and(|bound| bound <= self.requested)
    }

    /// The contributions with no certified bound.
    pub fn unbounded(&self) -> impl Iterator<Item = &DeviationContribution> {
        self.contributions
            .iter()
            .filter(|c| matches!(c.bound, DeviationBound::Unbounded(_)))
    }
}

/// A node's deviation while the graph is compiled: its contributions,
/// merged by path and detail.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Deviation(Vec<DeviationContribution>);

impl Default for Deviation {
    /// Unreported: a path that never states its deviation must not read as
    /// exact.
    fn default() -> Self {
        Self::one(
            DeviationPath::Unreported,
            "",
            DeviationBound::Unbounded("this path does not report its deviation"),
        )
    }
}

impl Deviation {
    pub(crate) fn one(path: DeviationPath, detail: &'static str, bound: DeviationBound) -> Self {
        Self(vec![DeviationContribution {
            path,
            detail,
            bound,
        }])
    }

    /// No contribution yet; [`Self::add`] fills it.
    pub(crate) const fn empty() -> Self {
        Self(Vec::new())
    }

    /// Fold one contribution in, keeping the worst per path and detail.
    pub(crate) fn add(&mut self, path: DeviationPath, detail: &'static str, bound: DeviationBound) {
        match self
            .0
            .iter_mut()
            .find(|c| c.path == path && c.detail == detail)
        {
            Some(existing) => existing.bound = existing.bound.worst(bound),
            None => self.0.push(DeviationContribution {
                path,
                detail,
                bound,
            }),
        }
    }

    /// Fold another node's contributions in.
    pub(crate) fn absorb(&mut self, other: &Self) {
        for c in &other.0 {
            self.add(c.path, c.detail, c.bound);
        }
    }

    /// The same contributions after a transform: lengths grow by at most
    /// its largest stretch.
    pub(crate) fn transformed(&self, transform: Transform3) -> Self {
        let factor = stretch(transform);
        Self(
            self.0
                .iter()
                .map(|c| DeviationContribution {
                    bound: c.bound.scaled(factor),
                    ..*c
                })
                .collect(),
        )
    }

    pub(crate) fn report(&self, requested: Scalar) -> DeviationReport {
        let bound = if self.0.is_empty() {
            None
        } else {
            self.0.iter().try_fold(0.0, |worst: Scalar, c| {
                c.bound.value().map(|value| worst.max(value))
            })
        };
        DeviationReport {
            requested,
            bound,
            contributions: self.0.clone(),
        }
    }
}

/// An upper bound on the largest factor a transform's linear part stretches
/// a length by: the largest column length when the columns are orthogonal,
/// else the Frobenius norm. The same rule the compiler uses to shrink a
/// chord budget into an instance's local space.
pub(crate) fn stretch(transform: Transform3) -> Scalar {
    let m = transform.matrix3;
    let (sx, sy, sz) = (m.x_axis.length(), m.y_axis.length(), m.z_axis.length());
    let eps = 32.0 * f64::EPSILON;
    let orthogonal = m.x_axis.dot(m.y_axis).abs() <= eps * sx * sy
        && m.x_axis.dot(m.z_axis).abs() <= eps * sx * sz
        && m.y_axis.dot(m.z_axis).abs() <= eps * sy * sz;
    if orthogonal {
        sx.max(sy).max(sz) * (1.0 + 3.0 * eps)
    } else {
        (sx * sx + sy * sy + sz * sz).sqrt()
    }
}
