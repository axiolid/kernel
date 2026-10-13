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
//! - Authored meshes are their own exact surface: `0`. An authored polygon
//!   face within the linear tolerance of its plane counts as planar and
//!   adds nothing. A face warped beyond it has no unique surface (#254):
//!   it is triangulated in its fit plane and reported
//!   [`DeviationBound::Certified`] under the detail
//!   `"non-planar authored face"`, never `Proven`. The value is the width
//!   of the slab the face's corners (holes included) span about that fit
//!   plane, largest minus smallest signed distance (#261): every reading
//!   of the face that lies in that slab (either diagonal triangulation,
//!   the face flattened onto the plane, a bilinear patch) is within it of
//!   the mesh wherever it lies over the face, and the mesh within it of
//!   each reading that covers the face (`crate::planar` gives the
//!   argument); the worst face of the mesh is reported.
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
//!   one circular arc by #231's; along an ellipse, a B-spline, or an
//!   elevated or banked curve without a grade break (#252, its derivative
//!   bounds in plan distance from `axiolid_reference::elevated`), certified
//!   against the exact tube (walls, bore and square end caps) by branch and
//!   bound, below; an elevated curve whose grade breaks is unbounded by
//!   name. A straight plan under a constant grade is a segment.
//!   Fixed-reference sweeps along a segment, or along a circle whose plane
//!   the reference is normal to, by #231; other frame laws, polyline and
//!   composite directrices are unbounded here.
//! - Primitives: spheres and tori by #231, cylinders and cones by #232
//!   (`axiolid_reference::primitive`), flat-faced ones exactly.
//! - B-rep faces: planar faces with straight edges exactly; a face that
//!   declares no surface and whose corners lie off their fit plane beyond
//!   the linear tolerance (#257) as a warped authored face is, by its
//!   slab width, [`DeviationBound::Certified`] under the detail
//!   `"non-planar face without a surface"`; curved faces
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
//! is always the worst bound over all cells. A tube (#252) settles once
//! within the budget, joins its slightly twisted quads as flat regions up
//! to a twentieth of the budget thick (that thickness is part of the
//! bound), and cuts each cell the way that lowers its halves more, so a
//! long tube's long, thin quads do not force every cell down to the width
//! of their diagonals.
//! - Instances scale the bound by an upper bound on their transform's
//!   largest stretch (`stretch`): the largest column length grown by
//!   `1 + 96 eps` when the columns are orthogonal, else the Frobenius
//!   norm. The margin keeps it above one even for an identity or a
//!   rotation, since an `f64` rotation is orthonormal only to rounding;
//!   the source is meshed to the budget shrunk by that factor, down to the
//!   largest value whose rounded product with it is within the budget
//!   (`shrunk_budget`, as #255 shrinks a curve-bounded plane's), so a
//!   rigid or scaled instance of a body meshed to exactly its budget
//!   reports at most the requested budget, nested instances included
//!   (#258). Uniform, non-uniform and sheared instances grow the bound by
//!   their full stretch and stay sound. The rounding of moving each vertex
//!   (a few ulps of its coordinates, so of the placement's distance from
//!   the origin) is not part of the bound, as the rounding of computing a
//!   construction's own vertices is not.
//!   Collections take the worst member.
//! - Runs of sections at stations (#247): stationed spines and sectioned
//!   surfaces are bounded cell by cell as they are sampled
//!   (`crate::station`'s `bound` module): the vertex paths' chord bound
//!   from certified derivative bounds of the frame law, breakpoints at
//!   seams inside a cell, and the wall quad's
//!   `axiolid_construct::loft::quad_deviation`, plus a spine's profile
//!   flattening, [`DeviationBound::Certified`]; a numerically framed basis,
//!   a cell beside a mitre on a curved piece or with a changing section,
//!   and different curved profiles are [`DeviationBound::Unbounded`] by
//!   name.
//! - Instances at stations (#264) move their source rigidly: on a line
//!   basis the frame is exact and the source's bound stands; on any other
//!   the frame is the exact curve's at a distance read by an estimate (the
//!   arc-length inverse, quadrature), so the placement is
//!   [`DeviationBound::Unbounded`] under the path
//!   [`DeviationPath::StationPlacement`].
//! - Booleans (#235): cutting two meshes moves the intersection curve, and
//!   a one-sided bound on each operand does not bound it, so nothing is
//!   derived from the operands. Where `ReferenceExactCompiler` builds the
//!   exact result of the same node (differences of placed extrusions,
//!   #228, and half-space clips, #234; for the report also booleans with
//!   a mesh operand, any operator, below), the mesh is certified against that
//!   exact B-rep by the branch and bound above, over each face's trimmed
//!   parameter domain; a boolean it refuses is unbounded with the
//!   refusal's name. The bound is relative to the exact compiler's result,
//!   the exact boolean of operands moved by at most the tolerance (#228),
//!   not to the boolean of the unperturbed operands. It runs only for a
//!   report, and only for the booleans whose result is emitted, not for
//!   the inner booleans of a chain; the search stops once within the
//!   requested budget.
//! - Operands snapped together (#276, #291): before the mesh boolean,
//!   operand coordinates a rounding residue (16 `f64::EPSILON` of the
//!   operands' largest coordinate) from the other operand's axis-aligned
//!   faces are set onto them exactly (`crate::snap`), so a skin that thin
//!   is not cut as one; a thicker skin is kept. The largest move is a
//!   [`DeviationBound::Certified`] contribution of
//!   [`DeviationPath::Boolean`] under the detail [`SNAPPED_OPERANDS`],
//!   carried through instances and enclosing booleans, the measured ones
//!   included; it never exceeds the tolerance. A boolean whose snapped
//!   operands are refused is cut as given and reports no snap.
//! - Warped operands (#266): an operand that is an authored mesh, or a
//!   B-rep whose faces are all planes with straight edges or warped faces
//!   without a surface, is measured as the exact polyhedron of its
//!   triangle mesh: a warped face flattened to its triangles, within its
//!   slab width `w` of every reading of the face (`crate::planar`), both
//!   ways. That is how an operand's certified bound composes into a
//!   boolean's: not as a one-sided bound (above), but as a two-sided
//!   perturbation of the operand. The boolean's mesh is measured against
//!   the exact boolean of the flattened operands (ADR 0080: the exact
//!   boolean of operands perturbed within the reported deviation), by `d`,
//!   and its contribution is [`DeviationBound::Certified`] `d + w`, `w`
//!   the largest slab width of any operand flattened beneath it, under
//!   the detail "measured against the exact boolean of operands flattened
//!   within their slab width, plus it". The search settles within the
//!   budget less `w`, so a boolean whose slab widths fit the budget
//!   reports within it, and `meets_requested` follows. Each warped
//!   operand's own contribution ("non-planar authored face", "non-planar
//!   face without a surface") is carried as the snap is, so the report
//!   names it. Transversal to the cut, the boolean of the unflattened operands
//!   is within `d + w` of the mesh too; near tangency the guarantee is
//!   ADR 0080's, relative to the flattened operands. Where the exact
//!   boolean cannot take the flattened operand (a mesh that is not a
//!   closed two-manifold, a B-rep with a curved face, a tool that is
//!   itself a boolean, a configuration the general boolean refuses) the
//!   boolean is [`DeviationBound::Unbounded`] with that refusal's name,
//!   and the warp is still named. An unwarped mesh operand (`w = 0`) is
//!   measured like any exact operand, against the exact boolean of the
//!   operands as given when no decision within tolerance fired.
//!
//! [`ReferenceMeshCompiler::compile_mesh_with_deviation`]: crate::ReferenceMeshCompiler::compile_mesh_with_deviation

use axiolid_core::{Scalar, Transform3};

mod boolean;
mod paths;
pub(crate) use boolean::{emitted_booleans, of_boolean};
pub(crate) use paths::{
    of_curve_bounded, of_primitive, of_sectioned_surface, of_solid, stationed_spine,
};

/// The [`DeviationContribution::detail`] of authored polygon faces warped
/// beyond the linear tolerance (#254).
pub(crate) const WARPED_AUTHORED_FACE: &str = "non-planar authored face";

/// The [`DeviationContribution::detail`] of a node placed at a station
/// whose frame is read at a numerically resolved distance (#264).
pub(crate) const NUMERICAL_STATION_FRAME: &str = "frame at a numerically resolved station";

/// The [`DeviationContribution::detail`] of B-rep faces that declare no
/// surface and are warped beyond the linear tolerance (#257).
pub(crate) const WARPED_BREP_FACE: &str = "non-planar face without a surface";

/// The [`DeviationContribution::detail`] of a boolean whose operands were
/// moved onto each other's faces by a rounding residue, within the linear
/// tolerance, before the mesh boolean cut them (#276, #291): its bound is
/// the largest distance any operand vertex moved, so the mesh is the
/// boolean of operands perturbed by at most that much. It travels with the boolean's result through
/// instances (scaled) and enclosing booleans, also where an enclosing
/// boolean is measured against its exact result.
pub const SNAPPED_OPERANDS: &str = "operands snapped onto each other within tolerance";

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
    /// A half-space bounded by a polygon or a profile.
    BoundedHalfSpace,
    /// A CSG primitive.
    Primitive,
    /// A face of a B-rep.
    BRepFace,
    /// A plane trimmed by boundary curves.
    CurveBoundedPlane,
    /// A boolean of two meshes, measured against its exact result (#235).
    ///
    /// An operand's own bound does not compose into it as a one-sided
    /// bound: a boolean is measured, not derived. A warped operand (#266)
    /// enters as a two-sided perturbation: the exact result is the boolean
    /// of the operands flattened to their triangles, each within its slab
    /// width `w` of every reading, and the contribution is the measured
    /// distance plus the largest `w` (see the module's "Warped operands").
    Boolean,
    /// A path that does not report its deviation.
    Unreported,
    /// Sections standing at stations along a directrix (#241).
    StationedSpine,
    /// Open sections at stations joined by tag (#241).
    SectionedSurface,
    /// A node placed in the frame of a curve station (#264), on a basis
    /// whose frame is not exact.
    StationPlacement,
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

    /// Fold in `other`'s operand perturbations: its [`SNAPPED_OPERANDS`]
    /// contributions (#276) and its warped faces' slab widths
    /// ([`WARPED_AUTHORED_FACE`], [`WARPED_BREP_FACE`], #266). They stay
    /// reported when a boolean replaces the rest of an operand's
    /// deviation, so the report names every operand that was moved or
    /// flattened, through any chain of booleans.
    pub(crate) fn carry_perturbations(&mut self, other: &Self) {
        for c in &other.0 {
            let carried = matches!(
                (c.path, c.detail),
                (DeviationPath::Boolean, SNAPPED_OPERANDS)
                    | (DeviationPath::AuthoredMesh, WARPED_AUTHORED_FACE)
                    | (DeviationPath::BRepFace, WARPED_BREP_FACE)
            );
            if carried {
                self.add(c.path, c.detail, c.bound);
            }
        }
    }

    /// The largest slab width of a B-rep's warped faces when its mesh is
    /// the polyhedron it describes up to them (#266): every contribution
    /// an exact plane with straight edges, or a warped face without a
    /// surface. `None` when any face is curved, chorded or off its plane.
    pub(crate) fn polyhedral_warp(&self) -> Option<Scalar> {
        self.0
            .iter()
            .try_fold(0.0, |widest: Scalar, c| match (c.path, c.detail, c.bound) {
                (DeviationPath::BRepFace, "plane", DeviationBound::Proven(0.0)) => Some(widest),
                (DeviationPath::BRepFace, WARPED_BREP_FACE, DeviationBound::Certified(width)) => {
                    Some(widest.max(width))
                }
                _ => None,
            })
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

    /// The worst certified bound above `budget`, if any: what a mesh
    /// rebuilt finer could still bring within it. A proven bound is the
    /// construction's own and does not move with a finer mesh.
    pub(crate) fn certified_over(&self, budget: Scalar) -> Option<Scalar> {
        self.0
            .iter()
            .filter_map(|c| match c.bound {
                DeviationBound::Certified(value) if value > budget => Some(value),
                _ => None,
            })
            .reduce(Scalar::max)
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
/// else the Frobenius norm. The same factor the compiler shrinks a chord
/// budget by into an instance's local space ([`shrunk_budget`]).
///
/// The `1 + 96 eps` margin on the orthogonal branch covers what the column
/// lengths do not: columns orthogonal only to `32 eps` (an `f64` rotation
/// matrix is never exactly orthonormal) raise the largest singular value
/// above the largest column length by at most about `32 eps` of it, and
/// each computed length is itself rounded. It is why a rigid transform's
/// factor is a hair above one, never exactly one: no `f64` rotation can be
/// shown to stretch nothing, so none is read as stretch 1 (#258). The
/// instance's source is instead meshed to a budget shrunk by this factor,
/// so the grown bound still fits the request.
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

/// The largest budget `b` from `budget / stretch` down whose product with
/// `stretch`, rounded as [`Deviation::transformed`] rounds it, is within
/// `budget` (#255, #258).
///
/// `budget / stretch` alone is not enough: at some factors (a scale of
/// 0.33 or 1.19 at a budget of `1e-4`) the quotient multiplied back
/// rounds one ulp above `budget`, and a source meshed to exactly its local
/// budget would report just over the request. Rounded multiplication is
/// monotone, so any source bound `d <= b` then reports `d * stretch <=
/// budget`, and nested instances keep that through every level.
/// `budget` itself when `stretch` is not a positive finite factor (the
/// callers refuse those first), and the plain quotient when no positive
/// value fits (a quotient in the subnormal range, which no mesh reaches).
pub(crate) fn shrunk_budget(budget: Scalar, stretch: Scalar) -> Scalar {
    if !(stretch.is_finite() && stretch > 0.0) {
        return budget;
    }
    let mut shrunk = budget / stretch;
    while shrunk > 0.0 && shrunk * stretch > budget {
        shrunk = shrunk.next_down();
    }
    if shrunk > 0.0 {
        shrunk
    } else {
        budget / stretch
    }
}

#[cfg(test)]
mod tests {
    use axiolid_core::{Mat3, Vec3};

    use super::*;

    /// Columns orthogonal only to `30 eps` pass as orthogonal, and the
    /// largest singular value then exceeds the largest column length: the
    /// margin must cover it (#258).
    #[test]
    fn the_stretch_covers_columns_orthogonal_only_to_rounding() {
        let skew = 30.0 * f64::EPSILON;
        let m = Mat3::from_cols(Vec3::X, Vec3::new(skew, 1.0, 0.0), Vec3::Z);
        // Both columns have length one in f64; the largest singular value
        // of [[1, skew], [0, 1]] is above 1 + skew / 2.
        assert_eq!(m.y_axis.length(), 1.0);
        let factor = stretch(Transform3::from_mat3(m));
        assert!(factor >= 1.0 + 0.5 * skew, "{factor}");
        assert!(factor < 1.0 + 200.0 * f64::EPSILON, "{factor}");
        assert!(stretch(Transform3::IDENTITY) > 1.0);
    }

    /// The shrunk budget times the stretch is within the budget, and the
    /// next value up would not be (unless the quotient already fit).
    #[test]
    fn the_shrunk_budget_fits_after_rounding() {
        for budget in [1e-4, 1e-3, 3.7e-5, 0.125] {
            for factor in [0.33, 1.0, 1.19, 2.38, 2.510_038_280_268_978, 1e6] {
                let s = factor * (1.0 + 96.0 * f64::EPSILON);
                let b = shrunk_budget(budget, s);
                assert!(b * s <= budget, "{budget} {factor}");
                assert!(b <= budget / s);
                assert!(b == budget / s || b.next_up() * s > budget);
            }
        }
        // The quotient alone rounds over at this factor.
        let s = 0.33 * (1.0 + 96.0 * f64::EPSILON);
        assert!((1e-4 / s) * s > 1e-4);
        assert!(shrunk_budget(1e-4, s) < 1e-4 / s);
    }
}
