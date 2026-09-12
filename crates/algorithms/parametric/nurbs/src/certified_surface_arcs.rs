//! Certified curved surface/surface analysis with explicit accounting.
//!
//! A whole-patch transversality bound cannot certify a curved pair: the
//! surface normal sweeps, so its interval hull straddles the other normal
//! and `|n1 x n2|` has lower bound zero. Affine patches escape this only
//! because a constant normal makes the bound sharp.
//!
//! So transversality is certified PER CELL. Cells that cannot be certified
//! are not dropped and not subdivided forever -- they are returned as
//! located [`RegionKind::Tangential`] or [`RegionKind::BudgetExhausted`]
//! regions, so a caller always learns which part of the domain is still
//! unproven.
//!
//! The leaves tile the root parameter box exactly. See
//! [`audit_coverage`] for the machine-checkable statement of that.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Scalar;
use axiolid_surface::BSplineSurface;

use crate::certified_refinement::RefinementBudget;
use crate::certified_surface_bezier::{piecewise_bezier_patches, Patch};
use crate::certified_surface_surface_intersection::{
    normal_cross_squared_lower_bound, patches_are_disjoint, SurfaceSurfaceParameterBox,
};
use crate::ParameterInterval;

/// Largest subdivision depth the coverage audit can represent exactly.
///
/// Each subdivision splits all four parameter axes, so a leaf at depth `d`
/// covers `16^-d` of its root box. The audit sums `16^(MAX_AUDIT_DEPTH - d)`
/// as `u128`, needing `4 * MAX_AUDIT_DEPTH <= 127` to stay exact. A larger
/// bound would silently saturate, and a saturated sum can mask a real gap.
pub const MAX_AUDIT_DEPTH: u32 = 31;

/// What was certified about one leaf region of the parameter domain.
///
/// Every leaf carries exactly one kind. A caller that treats
/// [`Self::Transversal`] as "the answer" and ignores the rest is
/// reading an incomplete result, which is why the unproven kinds carry
/// their own located boxes rather than being folded into a count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegionKind {
    /// Proven to contain no intersection: the patch bounding boxes are
    /// disjoint in at least one coordinate.
    Empty,
    /// Proven transversal: `|n1 x n2|` is bounded strictly away from zero,
    /// so any intersection here is a regular curve, never a surface patch
    /// or an isolated tangential touch.
    Transversal,
    /// Not proven transversal, and shrinking the cell will not help.
    ///
    /// The surfaces are tangent, near-tangent, or coincident somewhere in
    /// this box. This is a located refusal, not a failure: the region is
    /// returned so the caller can decide, refine by other means, or
    /// report it. It is never silently discarded.
    Tangential,
    /// Subdivision stopped on policy (depth or work budget) before the
    /// region could be classified either way.
    ///
    /// Distinct from [`Self::Tangential`]: that is a statement about the
    /// geometry, this is a statement about the budget. Raising the budget
    /// may resolve it; raising it cannot resolve a tangency.
    BudgetExhausted,
}

/// One classified leaf of the certified subdivision.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CertifiedRegion {
    /// The parameter box this statement applies to.
    pub box_: SurfaceSurfaceParameterBox,
    /// What was proven about it.
    pub kind: RegionKind,
    /// Bisection depth, where the root pair is depth 0.
    ///
    /// Each step splits all four parameter axes, so a leaf at depth `d`
    /// covers `2^-d` of the root box. The coverage audit relies on this.
    pub depth: u32,
    /// Certified lower bound on `|n1 x n2|^2` over the box.
    ///
    /// Strictly positive exactly when `kind` is [`RegionKind::Transversal`].
    pub normal_separation_lower_bound: Scalar,
}

/// Result of certified curved surface/surface analysis.
///
/// Deliberately not an `Option`-like "answer or nothing": a curved pair
/// is routinely *partly* provable, and collapsing that to a single
/// verdict is what loses geometry. The caller gets the transversal
/// regions and the unproven ones together, and can always ask whether
/// the analysis was complete.
#[derive(Debug, Clone, PartialEq)]
pub struct CertifiedSurfaceArcs3 {
    /// Every leaf, in deterministic order, tiling the root box exactly.
    pub regions: Vec<CertifiedRegion>,
    /// Patch pairs examined.
    pub visited_patch_pairs: u32,
    /// Deepest leaf produced.
    pub max_depth_reached: u32,
}

impl CertifiedSurfaceArcs3 {
    /// Whether every leaf was proven either empty or transversal.
    ///
    /// `false` means part of the domain is still unproven; the regions
    /// say which part. Callers that need a total answer must branch on
    /// this rather than assuming the arcs are exhaustive.
    pub fn is_fully_certified(&self) -> bool {
        self.regions
            .iter()
            .all(|region| matches!(region.kind, RegionKind::Empty | RegionKind::Transversal))
    }

    /// Regions that remain unproven, in deterministic order.
    pub fn unproven(&self) -> impl Iterator<Item = &CertifiedRegion> {
        self.regions.iter().filter(|region| {
            matches!(
                region.kind,
                RegionKind::Tangential | RegionKind::BudgetExhausted
            )
        })
    }
}

/// Why a set of regions failed to account for the whole domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoverageFault {
    /// The leaves cover less than the root box: geometry was dropped.
    ///
    /// This is the failure that silently loses an intersection, so it is
    /// reported as a fault rather than tolerated.
    Gap,
    /// The leaves cover more than the root box: a region was emitted twice
    /// or a parent was kept alongside its children.
    Overlap,
    /// A leaf is deeper than the audit can represent exactly.
    ///
    /// Reported instead of silently saturating, because a saturated sum
    /// could mask a real gap.
    DepthOverflow,
}

/// Verify that the regions tile the root parameter box exactly.
///
/// # The invariant
///
/// Subdivision bisects all four parameter axes, so a leaf at depth `d`
/// covers exactly `16^-d` of its root box -- a dyadic rational, never a
/// rounded quantity. Summing `16^(MAX_AUDIT_DEPTH - d)` over the leaves
/// gives exactly `16^MAX_AUDIT_DEPTH` per root box if and only if the
/// leaves tile every root with no gap and no overlap.
///
/// The arithmetic is integer `u128` throughout: no tolerance, no
/// accumulated float error, no judgement call. A dropped leaf makes the
/// sum too small, a duplicated one makes it too large, and either way
/// this returns a fault.
///
/// # Why this check outlives the code it checks
///
/// It constrains only *completeness*, never *content*. Future work may
/// add refusal kinds, sharpen bounds, or change how tangency is handled;
/// none of that is allowed to stop accounting for the domain. So this
/// keeps catching the invisible bug class without needing a rewrite.
pub fn audit_coverage(
    regions: &[CertifiedRegion],
    root_pair_count: u32,
) -> Result<(), CoverageFault> {
    let mut total: u128 = 0;
    for region in regions {
        if region.depth > MAX_AUDIT_DEPTH {
            return Err(CoverageFault::DepthOverflow);
        }
        // 16^(MAX - depth), exact: each level splits all four axes.
        total += 1u128 << (4 * (MAX_AUDIT_DEPTH - region.depth));
    }
    let per_root: u128 = 1u128 << (4 * MAX_AUDIT_DEPTH);
    let expected = per_root * u128::from(root_pair_count);
    match total.cmp(&expected) {
        std::cmp::Ordering::Equal => Ok(()),
        std::cmp::Ordering::Less => Err(CoverageFault::Gap),
        std::cmp::Ordering::Greater => Err(CoverageFault::Overlap),
    }
}

/// Explicit policy for certified curved analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CertifiedSurfaceArcsOptions {
    max_depth: u32,
}

impl CertifiedSurfaceArcsOptions {
    /// Validate and construct a policy.
    ///
    /// `max_depth` is capped at [`MAX_AUDIT_DEPTH`] so the coverage audit
    /// stays exact by construction rather than by caller discipline.
    pub fn new(max_depth: u32) -> GeomResult<Self> {
        if max_depth == 0 || max_depth > MAX_AUDIT_DEPTH {
            return Err(GeomError::InvalidInput(format!(
                "certified arc max_depth must be in 1..={MAX_AUDIT_DEPTH}"
            )));
        }
        Ok(Self { max_depth })
    }
}

impl Default for CertifiedSurfaceArcsOptions {
    fn default() -> Self {
        Self { max_depth: 8 }
    }
}

/// Certify the transversality structure of a curved surface pair.
///
/// Subdivides the parameter domain until each leaf is provably empty,
/// provably transversal, or provably not worth subdividing further, and
/// returns every leaf. Unlike the affine path this never collapses to a
/// single verdict, because a curved pair is routinely part-provable and
/// collapsing is what loses geometry.
///
/// The returned regions always tile the domain exactly; see
/// [`audit_coverage`].
pub fn certify_surface_arcs(
    first: &BSplineSurface,
    second: &BSplineSurface,
    options: CertifiedSurfaceArcsOptions,
) -> GeomResult<CertifiedSurfaceArcs3> {
    let mut budget = RefinementBudget::new(1_000_000, "certified surface arc refinement");
    let first_patches = piecewise_bezier_patches(first, &mut budget)?;
    let second_patches = piecewise_bezier_patches(second, &mut budget)?;
    let pair_count = first_patches
        .len()
        .checked_mul(second_patches.len())
        .ok_or(GeomError::BudgetExceeded {
            resource: "certified surface arc pair count",
        })?;
    let visited_patch_pairs = u32::try_from(pair_count).map_err(|_| GeomError::BudgetExceeded {
        resource: "certified surface arc pair count",
    })?;

    let mut regions = Vec::new();
    let mut max_depth_reached = 0;
    for first_patch in &first_patches {
        for second_patch in &second_patches {
            classify_recursive(
                first_patch,
                second_patch,
                0,
                options.max_depth,
                &mut regions,
                &mut max_depth_reached,
            )?;
        }
    }
    Ok(CertifiedSurfaceArcs3 {
        regions,
        visited_patch_pairs,
        max_depth_reached,
    })
}

/// Classify one cell pair, subdividing only when that can settle it.
///
/// Every control-flow path emits exactly one leaf or recurses into four
/// children. That is what makes exact coverage a property of the
/// structure rather than something a test has to hope for.
fn classify_recursive(
    first: &Patch,
    second: &Patch,
    depth: u32,
    max_depth: u32,
    regions: &mut Vec<CertifiedRegion>,
    max_depth_reached: &mut u32,
) -> GeomResult<()> {
    *max_depth_reached = (*max_depth_reached).max(depth);
    let box_ = pair_box(first, second);

    // Disjoint bounding boxes prove emptiness outright, at any depth.
    if patches_are_disjoint(first, second)? {
        push_region(regions, box_, RegionKind::Empty, depth, 0.0)?;
        return Ok(());
    }

    // A strictly positive normal separation proves the intersection here
    // is a regular curve. This is the same bound the affine path uses;
    // only the scope changed, from whole-patch to per-cell.
    let separation = normal_cross_squared_lower_bound(first, second)?;
    if separation > 0.0 {
        push_region(regions, box_, RegionKind::Transversal, depth, separation)?;
        return Ok(());
    }

    if depth >= max_depth {
        // Out of budget, not out of geometry. Distinguishing the two lets
        // a caller tell "raise the budget" from "this will never resolve".
        let kind = if is_persistently_tangential(first, second)? {
            RegionKind::Tangential
        } else {
            RegionKind::BudgetExhausted
        };
        push_region(regions, box_, kind, depth, 0.0)?;
        return Ok(());
    }

    let (first_children, second_children) = (split_patch(first)?, split_patch(second)?);
    for first_child in &first_children {
        for second_child in &second_children {
            classify_recursive(
                first_child,
                second_child,
                depth + 1,
                max_depth,
                regions,
                max_depth_reached,
            )?;
        }
    }
    Ok(())
}

/// Whether a cell pair looks tangential rather than merely under-refined.
///
/// Probes the four quadrant children: if none of them recovers a positive
/// normal separation, the obstruction is geometric (a tangency curve runs
/// through the cell) rather than a matter of cell size. This is a
/// heuristic for *labelling* an already-refused region, never a
/// certificate -- both labels mean "unproven", and the distinction only
/// tells the caller whether raising the budget could help.
fn is_persistently_tangential(first: &Patch, second: &Patch) -> GeomResult<bool> {
    let first_children = split_patch(first)?;
    let second_children = split_patch(second)?;
    for first_child in &first_children {
        for second_child in &second_children {
            if patches_are_disjoint(first_child, second_child)? {
                continue;
            }
            if normal_cross_squared_lower_bound(first_child, second_child)? > 0.0 {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Split a patch at its parameter midpoint in both directions.
fn split_patch(patch: &Patch) -> GeomResult<[Patch; 4]> {
    let u_mid = patch.u_start * 0.5 + patch.u_end * 0.5;
    let v_mid = patch.v_start * 0.5 + patch.v_end * 0.5;
    Ok([
        patch.restrict(patch.u_start, u_mid, patch.v_start, v_mid)?,
        patch.restrict(u_mid, patch.u_end, patch.v_start, v_mid)?,
        patch.restrict(patch.u_start, u_mid, v_mid, patch.v_end)?,
        patch.restrict(u_mid, patch.u_end, v_mid, patch.v_end)?,
    ])
}

fn pair_box(first: &Patch, second: &Patch) -> SurfaceSurfaceParameterBox {
    SurfaceSurfaceParameterBox {
        first_u: ParameterInterval {
            start: first.u_start,
            end: first.u_end,
        },
        first_v: ParameterInterval {
            start: first.v_start,
            end: first.v_end,
        },
        second_u: ParameterInterval {
            start: second.u_start,
            end: second.u_end,
        },
        second_v: ParameterInterval {
            start: second.v_start,
            end: second.v_end,
        },
    }
}

fn push_region(
    regions: &mut Vec<CertifiedRegion>,
    box_: SurfaceSurfaceParameterBox,
    kind: RegionKind,
    depth: u32,
    normal_separation_lower_bound: Scalar,
) -> GeomResult<()> {
    if regions.len() == regions.capacity() {
        regions
            .try_reserve(1)
            .map_err(|_| GeomError::BudgetExceeded {
                resource: "certified surface arc region allocation",
            })?;
    }
    regions.push(CertifiedRegion {
        box_,
        kind,
        depth,
        normal_separation_lower_bound,
    });
    Ok(())
}
