//! Certified curved surface/surface analysis.
//!
//! These tests pin the property that matters for trusting the result:
//! the classified regions account for the WHOLE parameter domain. A
//! result that quietly loses a region would still look plausible, so
//! completeness is asserted structurally rather than by eyeballing counts.

use axiolid_core::Point3;
use axiolid_curve::KnotSpec;
use axiolid_nurbs::{
    audit_coverage, certify_surface_arcs, CertifiedRegion, CertifiedSurfaceArcsOptions,
    CoverageFault, RegionKind, MAX_AUDIT_DEPTH,
};
use axiolid_surface::BSplineSurface;

/// A quarter cylinder of radius 1 about the origin, extruded along z.
///
/// Rational quadratic: the surface normal sweeps 90 degrees across the
/// patch, which is exactly what defeats a whole-patch transversality bound.
fn quarter_cylinder() -> BSplineSurface {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    BSplineSurface {
        u_degree: 2,
        v_degree: 1,
        control_points: vec![
            vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 2.0)],
            vec![Point3::new(1.0, 1.0, 0.0), Point3::new(1.0, 1.0, 2.0)],
            vec![Point3::new(0.0, 1.0, 0.0), Point3::new(0.0, 1.0, 2.0)],
        ],
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![3, 3],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![2, 2],
        weights: Some(vec![vec![1.0, 1.0], vec![w, w], vec![1.0, 1.0]]),
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    }
}

/// A flat plane at constant y, spanning the cylinder.
fn plane_at_y(y: f64) -> BSplineSurface {
    BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        control_points: vec![
            vec![Point3::new(-2.0, y, -1.0), Point3::new(-2.0, y, 3.0)],
            vec![Point3::new(2.0, y, -1.0), Point3::new(2.0, y, 3.0)],
        ],
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![2, 2],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![2, 2],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    }
}

/// A region with a real parameter box, used for the pure measure tests.
///
/// The box content is irrelevant to coverage -- the audit reasons about
/// depths -- so this borrows a genuine one from a trivial analysis rather
/// than fabricating a value the type deliberately does not allow.
fn region(kind: RegionKind, depth: u32) -> CertifiedRegion {
    let arcs = certify_surface_arcs(
        &plane_at_y(0.0),
        &plane_at_y(9.0),
        CertifiedSurfaceArcsOptions::new(1).expect("valid policy"),
    )
    .expect("analysis terminates");
    let box_ = arcs.regions.first().expect("at least one region").box_;
    CertifiedRegion {
        box_,
        kind,
        depth,
        normal_separation_lower_bound: 0.0,
    }
}

/// The audit accepts exactly one whole root box.
#[test]
fn a_single_undivided_region_covers_its_root() {
    assert_eq!(audit_coverage(&[region(RegionKind::Empty, 0)], 1), Ok(()));
}

/// Sixteen depth-1 children tile one root exactly.
///
/// Pins the measure arithmetic itself: each subdivision splits all four
/// parameter axes, so a child is 1/16 of its parent, not 1/2 or 1/4.
#[test]
fn sixteen_children_tile_one_root() {
    let children: Vec<_> = (0..16).map(|_| region(RegionKind::Empty, 1)).collect();
    assert_eq!(audit_coverage(&children, 1), Ok(()));
}

/// A dropped region is caught as a gap.
///
/// This is the whole point of the audit. Losing one child of sixteen
/// leaves a result that still looks structurally fine -- correct kinds,
/// plausible depths, no panic -- while an intersection has gone missing.
/// Only the measure sum reveals it.
#[test]
fn a_dropped_region_is_reported_as_a_gap() {
    let short: Vec<_> = (0..15).map(|_| region(RegionKind::Empty, 1)).collect();
    assert_eq!(audit_coverage(&short, 1), Err(CoverageFault::Gap));
}

/// A duplicated region is caught as an overlap.
///
/// The mirror failure: emitting a parent alongside its children, or
/// visiting a cell pair twice, double-counts geometry.
#[test]
fn a_duplicated_region_is_reported_as_an_overlap() {
    let long: Vec<_> = (0..17).map(|_| region(RegionKind::Empty, 1)).collect();
    assert_eq!(audit_coverage(&long, 1), Err(CoverageFault::Overlap));
}

/// Keeping a parent next to its own children is an overlap, not a gap.
///
/// A plausible refactoring mistake -- push the region, then also recurse --
/// and one that silently doubles part of the domain.
#[test]
fn a_parent_kept_beside_its_children_is_an_overlap() {
    let mut regions = vec![region(RegionKind::Tangential, 0)];
    regions.extend((0..16).map(|_| region(RegionKind::Empty, 1)));
    assert_eq!(audit_coverage(&regions, 1), Err(CoverageFault::Overlap));
}

/// Depth beyond the exact range is refused rather than saturated.
///
/// A shift past the integer width would wrap or saturate, and a saturated
/// total can mask a real gap. Refusing keeps the audit trustworthy.
#[test]
fn depth_beyond_the_exact_range_is_refused() {
    let deep = vec![region(RegionKind::Empty, MAX_AUDIT_DEPTH + 1)];
    assert_eq!(audit_coverage(&deep, 1), Err(CoverageFault::DepthOverflow));
}

/// A curved pair is analysed, and the analysis accounts for everything.
///
/// The plane cuts the quarter cylinder, so the result must contain real
/// transversal regions -- and whatever it could not prove must still be
/// present as located regions rather than missing.
#[test]
fn a_curved_pair_is_classified_with_complete_coverage() {
    let arcs = certify_surface_arcs(
        &quarter_cylinder(),
        &plane_at_y(0.5),
        CertifiedSurfaceArcsOptions::new(3).expect("valid policy"),
    )
    .expect("analysis terminates");

    assert_eq!(
        audit_coverage(&arcs.regions, arcs.visited_patch_pairs),
        Ok(()),
        "classified regions must account for the whole domain"
    );

    let transversal = arcs
        .regions
        .iter()
        .filter(|region| region.kind == RegionKind::Transversal)
        .count();
    assert!(
        transversal > 0,
        "a plane cutting a cylinder has transversal regions"
    );
}

/// Every transversal region carries a strictly positive separation bound.
///
/// The label and the evidence must agree: a region may only be called
/// transversal on the strength of a positive certified bound, never on
/// the strength of having reached some depth.
#[test]
fn transversal_regions_carry_positive_evidence() {
    let arcs = certify_surface_arcs(
        &quarter_cylinder(),
        &plane_at_y(0.5),
        CertifiedSurfaceArcsOptions::new(3).expect("valid policy"),
    )
    .expect("analysis terminates");

    for region in &arcs.regions {
        if region.kind == RegionKind::Transversal {
            assert!(
                region.normal_separation_lower_bound > 0.0,
                "transversal region without positive evidence: {region:?}"
            );
        }
    }
}

/// A tangency is reported as a located region, not hidden or dropped.
///
/// The plane y=1 touches the cylinder crest along a line. The honest
/// outcome is an explicit unproven region -- not a confident empty
/// result, and not subdividing forever chasing a bound that cannot exist.
#[test]
fn a_tangency_is_surfaced_as_an_unproven_region() {
    let arcs = certify_surface_arcs(
        &quarter_cylinder(),
        &plane_at_y(1.0),
        CertifiedSurfaceArcsOptions::new(3).expect("valid policy"),
    )
    .expect("analysis terminates");

    assert_eq!(
        audit_coverage(&arcs.regions, arcs.visited_patch_pairs),
        Ok(()),
        "a tangential pair must still account for the whole domain"
    );
    assert!(
        arcs.unproven().count() > 0,
        "tangency must be reported, not silently absorbed"
    );
    assert!(
        !arcs.is_fully_certified(),
        "a tangential pair must not claim full certification"
    );
}

/// Disjoint surfaces are fully certified with no unproven remainder.
///
/// Guards the opposite failure: a system that always reports something
/// unproven is as useless as one that never does.
#[test]
fn disjoint_surfaces_are_fully_certified() {
    let arcs = certify_surface_arcs(
        &quarter_cylinder(),
        &plane_at_y(50.0),
        CertifiedSurfaceArcsOptions::new(3).expect("valid policy"),
    )
    .expect("analysis terminates");

    assert_eq!(
        audit_coverage(&arcs.regions, arcs.visited_patch_pairs),
        Ok(())
    );
    assert!(
        arcs.is_fully_certified(),
        "a clearly disjoint pair must be fully certified"
    );
    assert_eq!(arcs.unproven().count(), 0);
}

/// The depth policy is validated rather than trusted.
#[test]
fn arc_policy_rejects_depths_outside_the_exact_range() {
    assert!(CertifiedSurfaceArcsOptions::new(0).is_err());
    assert!(CertifiedSurfaceArcsOptions::new(MAX_AUDIT_DEPTH + 1).is_err());
    assert!(CertifiedSurfaceArcsOptions::new(1).is_ok());
}
