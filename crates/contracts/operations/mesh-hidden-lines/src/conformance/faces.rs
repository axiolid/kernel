//! Shared conformance checks for portable visible-face providers.
//!
//! Like the hidden-line suite, it calls the provider directly, so a provider
//! must refuse bad input itself. Each exact case states the regions it must
//! return, bit for bit; the area cases state the area per source where the
//! coordinates are not binary64. Every outcome is also checked for shape
//! (ring orientation and start, sort order, triangle lists, caps, evidence).
//! The scenes and the arithmetic behind them are in the (private) fixtures,
//! checked by their own unit tests.

use core::cmp::Ordering;
use core::fmt;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Point2, Tolerance};

use crate::face_fixtures::{
    area_cases, cases, face_refusals, EvidenceExpect, FaceFixture, AREA_TOLERANCE,
};
use crate::fixtures::reversed;
use crate::{
    FaceSource, HiddenLineError, HiddenLineMesh, MeshId, MeshRole, MeshVisibleFaces,
    OrthographicView, VisibleFaceLimits, VisibleFaceOutcome, VisibleFaceRequest, VisibleRegion,
};

#[derive(Debug, Clone, PartialEq)]
pub enum FaceConformanceFailure {
    /// The provider refused a request it must answer, or refused a request
    /// that it must refuse with a different error.
    ReturnedError {
        case: &'static str,
        error: HiddenLineError,
    },
    /// The provider answered a request it must refuse.
    MissingRefusal { case: &'static str },
    /// The refusal names the wrong mesh, or none.
    WrongRefusalMesh {
        case: &'static str,
        expected: Option<MeshId>,
        actual: Option<MeshId>,
    },
    /// The regions differ from the expected ones.
    WrongRegions {
        case: &'static str,
        expected: Vec<VisibleRegion>,
        actual: Vec<VisibleRegion>,
    },
    /// The total area of the regions of one source is wrong.
    WrongArea {
        case: &'static str,
        mesh: MeshId,
        source: FaceSource,
        expected: f64,
        actual: f64,
    },
    /// The outcome has regions of sources other than the expected ones.
    WrongRegionCount {
        case: &'static str,
        expected: usize,
        actual: usize,
    },
    /// Evidence disagrees with the request, the outcome or the scene.
    IncorrectEvidence { case: &'static str },
    /// Two identical requests returned different outcomes.
    NonDeterministic { case: &'static str },
    /// Reversing every triangle changed the outcome.
    OrientationDependent { case: &'static str },
    /// A region breaks the shape the contract fixes.
    InvalidShape { case: &'static str, detail: String },
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct FaceConformanceReport {
    pub failures: Vec<FaceConformanceFailure>,
}

impl FaceConformanceReport {
    pub fn is_success(&self) -> bool {
        self.failures.is_empty()
    }
}

impl fmt::Display for FaceConformanceReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_success() {
            write!(f, "conformant")
        } else {
            write!(f, "{:?}", self.failures)
        }
    }
}

pub struct FaceConformanceSuite;

/// What the generic checks need to know about the request.
struct Scene<'a> {
    name: &'static str,
    meshes: &'a [FaceFixture],
    view: OrthographicView,
}

impl FaceConformanceSuite {
    pub fn run(provider: &dyn MeshVisibleFaces) -> FaceConformanceReport {
        let mut report = FaceConformanceReport::default();
        let options = ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("valid tolerance"));

        for case in cases() {
            let scene = Scene {
                name: case.name,
                meshes: &case.meshes,
                view: case.view,
            };
            let first = answer(
                provider,
                &case.meshes,
                case.view,
                case.region,
                case.limits,
                &options,
            );
            let second = answer(
                provider,
                &case.meshes,
                case.view,
                case.region,
                case.limits,
                &options,
            );
            match (first, second) {
                (Ok(a), Ok(b)) => {
                    check_shape(&mut report, &scene, &a);
                    if a.regions != case.regions {
                        report.failures.push(FaceConformanceFailure::WrongRegions {
                            case: case.name,
                            expected: case.regions.clone(),
                            actual: a.regions.clone(),
                        });
                    }
                    check_expected_evidence(&mut report, case.name, &case.evidence, &a);
                    if a != b {
                        report
                            .failures
                            .push(FaceConformanceFailure::NonDeterministic { case: case.name });
                    }
                    check_orientation(
                        provider,
                        &mut report,
                        &scene,
                        case.region,
                        case.limits,
                        &a,
                        &options,
                    );
                }
                (Err(error), _) | (_, Err(error)) => {
                    report.failures.push(FaceConformanceFailure::ReturnedError {
                        case: case.name,
                        error,
                    })
                }
            }
        }

        for case in area_cases() {
            let scene = Scene {
                name: case.name,
                meshes: &case.meshes,
                view: case.view,
            };
            let first = answer(
                provider,
                &case.meshes,
                case.view,
                None,
                case.limits,
                &options,
            );
            let second = answer(
                provider,
                &case.meshes,
                case.view,
                None,
                case.limits,
                &options,
            );
            match (first, second) {
                (Ok(a), Ok(b)) => {
                    check_shape(&mut report, &scene, &a);
                    check_areas(&mut report, &case.areas, case.name, &a);
                    if a != b {
                        report
                            .failures
                            .push(FaceConformanceFailure::NonDeterministic { case: case.name });
                    }
                    check_orientation(
                        provider,
                        &mut report,
                        &scene,
                        None,
                        case.limits,
                        &a,
                        &options,
                    );
                }
                (Err(error), _) | (_, Err(error)) => {
                    report.failures.push(FaceConformanceFailure::ReturnedError {
                        case: case.name,
                        error,
                    })
                }
            }
        }

        for refusal in face_refusals() {
            match answer(
                provider,
                &refusal.meshes,
                refusal.view,
                None,
                refusal.limits,
                &options,
            ) {
                Ok(_) => report
                    .failures
                    .push(FaceConformanceFailure::MissingRefusal { case: refusal.name }),
                Err(error) if !(refusal.accepts)(&error.error) => {
                    report.failures.push(FaceConformanceFailure::ReturnedError {
                        case: refusal.name,
                        error,
                    })
                }
                Err(error) if error.mesh != refusal.mesh => {
                    report
                        .failures
                        .push(FaceConformanceFailure::WrongRefusalMesh {
                            case: refusal.name,
                            expected: refusal.mesh,
                            actual: error.mesh,
                        })
                }
                Err(_) => {}
            }
        }

        report
    }
}

fn answer(
    provider: &dyn MeshVisibleFaces,
    fixtures: &[FaceFixture],
    view: OrthographicView,
    region: Option<axiolid_core::Aabb2>,
    limits: VisibleFaceLimits,
    options: &ExecutionOptions,
) -> Result<VisibleFaceOutcome, HiddenLineError> {
    let meshes: Vec<HiddenLineMesh<'_>> = fixtures
        .iter()
        .map(|f| {
            let mesh = HiddenLineMesh::new(f.id, f.role, &f.mesh);
            if f.groups.is_empty() {
                mesh
            } else {
                mesh.with_groups(&f.groups)
            }
        })
        .collect();
    let mut request = VisibleFaceRequest::new(&meshes, view, limits);
    if let Some(region) = region {
        request = request.with_region(region);
    }
    provider.visible_faces(&request, options)
}

/// Equal up to `exact_escalations`, which may differ between two otherwise
/// identical requests.
fn same_result(a: &VisibleFaceOutcome, b: &VisibleFaceOutcome) -> bool {
    let mut evidence = b.evidence;
    evidence.exact_escalations = a.evidence.exact_escalations;
    a.view == b.view && a.regions == b.regions && a.evidence == evidence
}

fn check_orientation(
    provider: &dyn MeshVisibleFaces,
    report: &mut FaceConformanceReport,
    scene: &Scene<'_>,
    region: Option<axiolid_core::Aabb2>,
    limits: VisibleFaceLimits,
    expected: &VisibleFaceOutcome,
    options: &ExecutionOptions,
) {
    let flipped: Vec<FaceFixture> = scene
        .meshes
        .iter()
        .map(|f| FaceFixture {
            id: f.id,
            role: f.role,
            mesh: reversed(&f.mesh),
            groups: f.groups.clone(),
        })
        .collect();
    match answer(provider, &flipped, scene.view, region, limits, options) {
        Ok(outcome) if same_result(&outcome, expected) => {}
        Ok(_) => report
            .failures
            .push(FaceConformanceFailure::OrientationDependent { case: scene.name }),
        Err(error) => report.failures.push(FaceConformanceFailure::ReturnedError {
            case: scene.name,
            error,
        }),
    }
}

fn signed_area(ring: &[Point2]) -> f64 {
    let n = ring.len();
    (0..n)
        .map(|i| {
            let (a, b) = (ring[i], ring[(i + 1) % n]);
            a.x * b.y - b.x * a.y
        })
        .sum::<f64>()
        / 2.0
}

fn point_order(a: &Point2, b: &Point2) -> Ordering {
    a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y))
}

fn sequence_order(a: &[Point2], b: &[Point2]) -> Ordering {
    a.iter()
        .zip(b)
        .map(|(p, q)| point_order(p, q))
        .find(|o| o.is_ne())
        .unwrap_or_else(|| a.len().cmp(&b.len()))
}

/// Whether the near plane cuts a solid: some vertex is beyond it (`W >
/// near`) and some is not.
fn cut_by_near(fixture: &FaceFixture, view: &OrthographicView) -> bool {
    let Some(near) = view.near else {
        return false;
    };
    let depth = |p: &axiolid_core::Point3| (*p - view.frame.origin).dot(view.frame.z);
    let beyond = fixture
        .mesh
        .positions
        .iter()
        .filter(|p| depth(p) > near)
        .count();
    fixture.role == MeshRole::Solid && beyond > 0 && beyond < fixture.mesh.positions.len()
}

fn check_ring(detail: &mut Option<String>, what: &str, ring: &[Point2], outer: bool) {
    let mut fail = |text: String| {
        detail.get_or_insert(text);
    };
    if ring.len() < 3 {
        return fail(format!("{what} has fewer than three points"));
    }
    if ring.iter().any(|p| !p.x.is_finite() || !p.y.is_finite()) {
        return fail(format!("{what} has a non-finite point"));
    }
    if ring.iter().any(|p| point_order(p, &ring[0]).is_lt()) {
        fail(format!("{what} does not start at its smallest point"));
    }
    if ring.first() == ring.last() {
        fail(format!("{what} repeats its first point"));
    }
    let area = signed_area(ring);
    if outer && area <= 0.0 {
        fail(format!("{what} is not counter-clockwise"));
    }
    if !outer && area >= 0.0 {
        fail(format!("{what} is not clockwise"));
    }
}

/// The checks every outcome must pass, whatever the scene.
fn check_shape(
    report: &mut FaceConformanceReport,
    scene: &Scene<'_>,
    outcome: &VisibleFaceOutcome,
) {
    let mut invalid = |detail: String| {
        report.failures.push(FaceConformanceFailure::InvalidShape {
            case: scene.name,
            detail,
        });
    };
    let position = |id: MeshId| scene.meshes.iter().position(|m| m.id == id);
    let mut previous: Option<(usize, FaceSource, &[Point2])> = None;
    let mut caps = 0;
    for region in &outcome.regions {
        let Some(index) = position(region.mesh) else {
            invalid(format!(
                "region of mesh {} that is not in the request",
                region.mesh.0
            ));
            continue;
        };
        let fixture = &scene.meshes[index];
        let mut detail = None;
        check_ring(&mut detail, "outer ring", &region.outer, true);
        for hole in &region.holes {
            check_ring(&mut detail, "hole", hole, false);
        }
        if region
            .holes
            .windows(2)
            .any(|w| sequence_order(&w[0], &w[1]).is_ge())
        {
            detail.get_or_insert("holes are not strictly sorted".into());
        }
        if region.triangles.windows(2).any(|w| w[0] >= w[1]) {
            detail.get_or_insert("triangles are not strictly ascending".into());
        }
        let count = fixture.mesh.triangle_count();
        if region.triangles.iter().any(|&t| t as usize >= count) {
            detail.get_or_insert("a triangle index is out of range".into());
        }
        match region.source {
            FaceSource::Group(group) => {
                if region
                    .triangles
                    .iter()
                    .any(|&t| (t as usize) < count && fixture.group(t as usize) != group)
                {
                    detail.get_or_insert("a triangle is not in the region's group".into());
                }
            }
            FaceSource::Cap => {
                caps += 1;
                if !cut_by_near(fixture, &scene.view) {
                    detail.get_or_insert("cap of a mesh the near plane does not cut".into());
                }
                if !region.triangles.is_empty() {
                    detail.get_or_insert("cap lists triangles".into());
                }
            }
        }
        if let Some((prev_mesh, prev_source, prev_outer)) = previous {
            let order = prev_mesh
                .cmp(&index)
                .then(prev_source.cmp(&region.source))
                .then_with(|| sequence_order(prev_outer, &region.outer));
            if order.is_ge() {
                detail.get_or_insert("regions are not strictly sorted".into());
            }
        }
        previous = Some((index, region.source, &region.outer));
        if let Some(detail) = detail {
            invalid(format!("mesh {}: {detail}", region.mesh.0));
        }
    }

    let evidence = &outcome.evidence;
    let triangles: usize = scene.meshes.iter().map(|m| m.mesh.triangle_count()).sum();
    if !evidence.is_derived_from_input_mesh()
        || evidence.source_meshes != scene.meshes.len()
        || evidence.source_triangles != triangles
        || evidence.regions != outcome.regions.len()
        || evidence.cap_regions != caps
    {
        report
            .failures
            .push(FaceConformanceFailure::IncorrectEvidence { case: scene.name });
    }
}

fn check_expected_evidence(
    report: &mut FaceConformanceReport,
    case: &'static str,
    expect: &EvidenceExpect,
    outcome: &VisibleFaceOutcome,
) {
    let evidence = &outcome.evidence;
    let wrong = (expect.intersections && evidence.intersection_segments == 0)
        || (expect.ties && evidence.tied_sides == 0)
        || expect
            .edge_on
            .is_some_and(|n| evidence.edge_on_triangles != n);
    if wrong {
        report
            .failures
            .push(FaceConformanceFailure::IncorrectEvidence { case });
    }
}

fn check_areas(
    report: &mut FaceConformanceReport,
    expected: &[(MeshId, FaceSource, f64)],
    case: &'static str,
    outcome: &VisibleFaceOutcome,
) {
    let known = |r: &VisibleRegion| {
        expected
            .iter()
            .any(|&(m, s, _)| m == r.mesh && s == r.source)
    };
    let unexpected = outcome.regions.iter().filter(|r| !known(r)).count();
    if unexpected > 0 {
        report
            .failures
            .push(FaceConformanceFailure::WrongRegionCount {
                case,
                expected: expected.len(),
                actual: unexpected + expected.len(),
            });
    }
    for &(mesh, source, area) in expected {
        let actual: f64 = outcome
            .regions
            .iter()
            .filter(|r| r.mesh == mesh && r.source == source)
            .map(|r| signed_area(&r.outer) + r.holes.iter().map(|h| signed_area(h)).sum::<f64>())
            .sum();
        // Written so that a NaN area fails.
        let close = (actual - area).abs() <= AREA_TOLERANCE * area;
        if !close {
            report.failures.push(FaceConformanceFailure::WrongArea {
                case,
                mesh,
                source,
                expected: area,
                actual,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use axiolid_contracts::{Backend, BackendDescriptor, BackendId, ExecutionTarget, GeomError};

    use super::*;
    use crate::VisibleFaceEvidence;

    const REFUSAL_CASES: [&str; 7] = [
        "three triangles on one edge",
        "inconsistent winding",
        "sheet declared solid",
        "group count mismatch",
        "region limit",
        "region vertex limit",
        "candidate limit",
    ];

    /// Refuses every request: the offending mesh (the second of two) for a
    /// per-mesh refusal, or the whole request when `blame_offender` is off,
    /// and the budget the limits name.
    #[derive(Debug)]
    struct Refuser {
        blame_offender: bool,
    }

    impl Backend for Refuser {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(
                BackendId::new("face-conformance-refuser"),
                ExecutionTarget::PortableCpu,
            )
        }
    }

    impl MeshVisibleFaces for Refuser {
        fn visible_faces(
            &self,
            request: &VisibleFaceRequest<'_>,
            _options: &ExecutionOptions,
        ) -> Result<VisibleFaceOutcome, HiddenLineError> {
            if request.meshes.len() == 2 && request.meshes[1].id == MeshId(2) {
                let mut error = GeomError::NotManifold("refused".into());
                let groups = request.meshes[1].groups;
                if matches!(groups, crate::FaceGroups::PerTriangle(g) if g.len() == 4) {
                    error = GeomError::InvalidInput("refused".into());
                }
                if request.meshes[1].mesh.triangle_count() == 12
                    && matches!(groups, crate::FaceGroups::Single)
                {
                    // The two-box scenes carry budgets; refuse by limits below.
                    return self.budget(request);
                }
                return Err(if self.blame_offender {
                    HiddenLineError::mesh(request.meshes[1].id, error)
                } else {
                    HiddenLineError::request(error)
                });
            }
            Err(GeomError::Cancelled.into())
        }
    }

    impl Refuser {
        fn budget(
            &self,
            request: &VisibleFaceRequest<'_>,
        ) -> Result<VisibleFaceOutcome, HiddenLineError> {
            let resource = if request.limits.max_regions == 1 {
                "visible-face regions"
            } else if request.limits.max_region_vertices == 5 {
                "visible-face region vertices"
            } else if request.limits.max_candidate_pairs == 1 {
                "visible-face candidate pairs"
            } else {
                return Err(GeomError::Cancelled.into());
            };
            Err(GeomError::BudgetExceeded { resource }.into())
        }
    }

    fn refusal_failures(report: &FaceConformanceReport) -> Vec<&FaceConformanceFailure> {
        report
            .failures
            .iter()
            .filter(|failure| match failure {
                FaceConformanceFailure::ReturnedError { case, .. }
                | FaceConformanceFailure::MissingRefusal { case }
                | FaceConformanceFailure::WrongRefusalMesh { case, .. } => {
                    REFUSAL_CASES.contains(case)
                }
                _ => false,
            })
            .collect()
    }

    #[test]
    fn a_provider_that_names_the_offending_mesh_passes_the_refusal_cases() {
        let report = FaceConformanceSuite::run(&Refuser {
            blame_offender: true,
        });
        assert!(refusal_failures(&report).is_empty(), "{report}");
        // It answers nothing, so every region case fails.
        assert!(!report.is_success());
    }

    #[test]
    fn a_refusal_that_loses_the_mesh_id_fails_exactly_the_mesh_cases() {
        let report = FaceConformanceSuite::run(&Refuser {
            blame_offender: false,
        });
        let wrong = refusal_failures(&report);
        assert_eq!(wrong.len(), 4, "{wrong:?}");
        for failure in wrong {
            assert!(matches!(
                failure,
                FaceConformanceFailure::WrongRefusalMesh { actual: None, .. }
            ));
        }
    }

    /// Answers every request with no regions and honest evidence.
    #[derive(Debug)]
    struct Empty;

    impl Backend for Empty {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(
                BackendId::new("face-conformance-empty"),
                ExecutionTarget::PortableCpu,
            )
        }
    }

    impl MeshVisibleFaces for Empty {
        fn visible_faces(
            &self,
            request: &VisibleFaceRequest<'_>,
            _options: &ExecutionOptions,
        ) -> Result<VisibleFaceOutcome, HiddenLineError> {
            let triangles = request.meshes.iter().map(|m| m.mesh.triangle_count()).sum();
            Ok(VisibleFaceOutcome::new(
                request.view,
                Vec::new(),
                VisibleFaceEvidence::input_mesh(request.meshes.len(), triangles),
            ))
        }
    }

    #[test]
    fn an_empty_answer_and_missing_refusals_are_reported() {
        let report = FaceConformanceSuite::run(&Empty);
        let cases = cases().len();
        let missing = report
            .failures
            .iter()
            .filter(|f| matches!(f, FaceConformanceFailure::MissingRefusal { .. }))
            .count();
        let regions = report
            .failures
            .iter()
            .filter(|f| matches!(f, FaceConformanceFailure::WrongRegions { .. }))
            .count();
        let areas = report
            .failures
            .iter()
            .filter(|f| matches!(f, FaceConformanceFailure::WrongArea { .. }))
            .count();
        assert_eq!(missing, REFUSAL_CASES.len(), "{report}");
        assert_eq!(regions, cases, "{report}");
        assert_eq!(areas, 3, "{report}");
    }

    /// Returns the expected regions of one exact case, whatever it is asked,
    /// with evidence that is deliberately wrong.
    #[derive(Debug)]
    struct Fixed(Vec<VisibleRegion>);

    impl Backend for Fixed {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(
                BackendId::new("face-conformance-fixed"),
                ExecutionTarget::PortableCpu,
            )
        }
    }

    impl MeshVisibleFaces for Fixed {
        fn visible_faces(
            &self,
            request: &VisibleFaceRequest<'_>,
            _options: &ExecutionOptions,
        ) -> Result<VisibleFaceOutcome, HiddenLineError> {
            Ok(VisibleFaceOutcome::new(
                request.view,
                self.0.clone(),
                VisibleFaceEvidence::input_mesh(0, 0),
            ))
        }
    }

    #[test]
    fn evidence_that_disagrees_with_the_request_is_reported() {
        let first = cases().remove(0);
        let report = FaceConformanceSuite::run(&Fixed(first.regions.clone()));
        let hit = report.failures.iter().any(|f| {
            matches!(f, FaceConformanceFailure::IncorrectEvidence { case } if *case == first.name)
        });
        assert!(hit, "{report}");
        // The regions themselves are right for that case.
        assert!(!report.failures.iter().any(|f| {
            matches!(f, FaceConformanceFailure::WrongRegions { case, .. } if *case == first.name)
        }));
    }

    #[test]
    fn the_expected_regions_pass_the_shape_checks() {
        for case in cases() {
            let triangles: usize = case.meshes.iter().map(|m| m.mesh.triangle_count()).sum();
            let caps = case
                .regions
                .iter()
                .filter(|r| r.source == FaceSource::Cap)
                .count();
            let mut evidence = crate::VisibleFaceEvidence::input_mesh(case.meshes.len(), triangles);
            evidence.regions = case.regions.len();
            evidence.cap_regions = caps;
            let outcome = VisibleFaceOutcome::new(case.view, case.regions.clone(), evidence);
            let scene = Scene {
                name: case.name,
                meshes: &case.meshes,
                view: case.view,
            };
            let mut report = FaceConformanceReport::default();
            check_shape(&mut report, &scene, &outcome);
            assert!(report.is_success(), "{}: {report}", case.name);
        }
    }

    #[test]
    fn shape_checks_catch_each_kind_of_bad_region() {
        let case = cases().remove(0);
        let scene = Scene {
            name: case.name,
            meshes: &case.meshes,
            view: case.view,
        };
        let triangles = 24;
        let mut evidence = crate::VisibleFaceEvidence::input_mesh(2, triangles);
        evidence.regions = 2;
        let run = |regions: Vec<VisibleRegion>| {
            let mut e = evidence;
            e.regions = regions.len();
            let mut report = FaceConformanceReport::default();
            check_shape(
                &mut report,
                &scene,
                &VisibleFaceOutcome::new(case.view, regions, e),
            );
            report
        };
        let good = case.regions.clone();
        assert!(run(good.clone()).is_success());

        let mutate = |f: &dyn Fn(&mut Vec<VisibleRegion>)| {
            let mut regions = good.clone();
            f(&mut regions);
            run(regions)
        };
        let shaped = |r: &FaceConformanceReport| {
            r.failures
                .iter()
                .any(|f| matches!(f, FaceConformanceFailure::InvalidShape { .. }))
        };
        // Clockwise outer ring.
        assert!(shaped(&mutate(&|r| r[0].outer.reverse())));
        // Ring that does not start at its smallest point.
        assert!(shaped(&mutate(&|r| r[0].outer.rotate_left(1))));
        // Regions out of order.
        assert!(shaped(&mutate(&|r| r.reverse())));
        // Unsorted or out-of-range triangles.
        assert!(shaped(&mutate(&|r| r[0].triangles = vec![3, 2])));
        assert!(shaped(&mutate(&|r| r[0].triangles = vec![2, 99])));
        // A triangle of another group; a cap of an uncut solid.
        assert!(shaped(&mutate(&|r| r[0].source = FaceSource::Cap)));
        // A hole that is counter-clockwise.
        assert!(shaped(&mutate(&|r| {
            r[0].holes = vec![vec![
                Point2::new(0.5, 0.5),
                Point2::new(1.0, 0.5),
                Point2::new(1.0, 1.0),
            ]]
        })));
        // A mesh that is not in the request.
        assert!(shaped(&mutate(&|r| r[0].mesh = MeshId(99))));
    }
}
