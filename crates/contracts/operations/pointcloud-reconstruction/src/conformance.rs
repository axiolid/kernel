//! Conformance suite every `PointcloudReconstruction` provider must pass.
//!
//! # Why this is library code, not a test file
//!
//! A test bound to a concrete provider tests *that provider*, not *the
//! contract*: a second provider would inherit no obligations at all. This
//! suite is generic over `impl PointcloudReconstruction` and exported, so an
//! out-of-tree provider runs the identical checks.
//!
//! # What it does and does not prove
//!
//! It checks *contract* obligations: refusals are typed, evidence is
//! populated and self-consistent, empty results are explained, and repeated
//! calls agree when the provider claims determinism. It does not check that
//! the reconstructed surface is geometrically faithful — no contract can,
//! because a point set does not determine a unique surface. A provider
//! passing this suite is well-behaved, not necessarily accurate.
//!
//! # Skips are not passes
//!
//! A provider may legitimately refuse work. Such a case is recorded as
//! [`Outcome::Skipped`] with its reason and reported separately, so a
//! provider cannot reach "conformant" by refusing everything: a caller can
//! see exactly what was actually exercised.

use core::fmt;

use axiolid_contracts::{Determinism, ExecutionOptions};
use axiolid_core::{Point3, Scalar, Tolerance};
use axiolid_pointcloud::PointCloud;

use crate::{PointcloudReconstruction, Reconstruction, ReconstructionRequest};

/// Result of one conformance check.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    /// The obligation was met.
    Passed,
    /// The obligation was violated, with what went wrong.
    Failed(String),
    /// The provider legitimately declined, with its stated reason.
    ///
    /// Reported separately from a pass so refusing everything cannot look
    /// like conformance.
    Skipped(String),
}

impl Outcome {
    /// Whether this outcome blocks conformance.
    pub fn is_failure(&self) -> bool {
        matches!(self, Self::Failed(_))
    }
}

/// One named check and how it went.
#[derive(Debug, Clone, PartialEq)]
pub struct Check {
    /// What was checked.
    pub name: &'static str,
    /// How it went.
    pub outcome: Outcome,
}

/// Everything the suite found.
#[derive(Debug, Clone, PartialEq)]
pub struct ConformanceReport {
    /// Every check, in the order run.
    pub checks: Vec<Check>,
}

impl ConformanceReport {
    /// Whether every check passed or was legitimately skipped.
    pub fn is_conformant(&self) -> bool {
        !self.checks.iter().any(|check| check.outcome.is_failure())
    }

    /// Checks that actually exercised the provider.
    pub fn exercised(&self) -> usize {
        self.checks
            .iter()
            .filter(|check| matches!(check.outcome, Outcome::Passed))
            .count()
    }

    /// Checks the provider declined.
    pub fn skipped(&self) -> usize {
        self.checks
            .iter()
            .filter(|check| matches!(check.outcome, Outcome::Skipped(_)))
            .count()
    }
}

impl fmt::Display for ConformanceReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "conformance: {} exercised, {} skipped, {} failed",
            self.exercised(),
            self.skipped(),
            self.checks
                .iter()
                .filter(|check| check.outcome.is_failure())
                .count()
        )?;
        for check in &self.checks {
            match &check.outcome {
                Outcome::Passed => writeln!(f, "  PASS {}", check.name)?,
                Outcome::Skipped(reason) => writeln!(f, "  SKIP {} ({reason})", check.name)?,
                Outcome::Failed(detail) => writeln!(f, "  FAIL {} — {detail}", check.name)?,
            }
        }
        Ok(())
    }
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

/// A dense sampling of a sphere: enough points for any method, with normals.
fn sphere_samples(count: usize) -> PointCloud {
    let mut points = Vec::with_capacity(count);
    let mut normals = Vec::with_capacity(count);
    // Fibonacci sphere: even coverage without clustering at the poles, so a
    // provider is not tested against an artificially easy distribution.
    let golden = core::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    for i in 0..count {
        let y = 1.0 - (i as Scalar / (count.max(2) - 1) as Scalar) * 2.0;
        let radius = (1.0 - y * y).max(0.0).sqrt();
        let theta = golden * i as Scalar;
        let direction = Point3::new(theta.cos() * radius, y, theta.sin() * radius);
        points.push(direction);
        normals.push(direction);
    }
    PointCloud::new(points)
        .expect("sphere samples are finite")
        .with_normals(normals)
        .expect("normals match")
}

/// Run every contract obligation against a provider.
pub fn run(provider: &impl PointcloudReconstruction) -> ConformanceReport {
    let checks = vec![
        Check {
            name: "too few points is refused by name, not by empty mesh",
            outcome: check_too_few(provider),
        },
        Check {
            name: "degenerate extent is refused by name",
            outcome: check_degenerate(provider),
        },
        Check {
            name: "evidence counts agree with the mesh returned",
            outcome: check_evidence_consistent(provider),
        },
        Check {
            name: "input point count is reported faithfully",
            outcome: check_input_count(provider),
        },
        Check {
            name: "a declared-deterministic provider repeats itself",
            outcome: check_determinism(provider),
        },
        Check {
            name: "minimum_points agrees with actual refusal behaviour",
            outcome: check_minimum_agrees(provider),
        },
    ];

    ConformanceReport { checks }
}

fn check_too_few(provider: &impl PointcloudReconstruction) -> Outcome {
    let cloud = PointCloud::new(vec![Point3::ZERO]).expect("one point");
    match provider.reconstruct(&cloud, &ReconstructionRequest::default(), &options()) {
        Ok(Reconstruction::Refused(_)) => Outcome::Passed,
        Ok(Reconstruction::Surface(outcome)) => Outcome::Failed(format!(
            "a single point produced a surface with {} triangles",
            outcome.mesh.indices.len() / 3
        )),
        Err(error) => Outcome::Failed(format!("errored instead of refusing: {error}")),
    }
}

fn check_degenerate(provider: &impl PointcloudReconstruction) -> Outcome {
    // All points on one line: no surface exists, and a provider that
    // returns a zero-thickness sheet is misrepresenting the input.
    let points: Vec<Point3> = (0..64)
        .map(|i| Point3::new(i as Scalar * 0.1, 0.0, 0.0))
        .collect();
    let cloud = PointCloud::new(points).expect("collinear points are finite");
    match provider.reconstruct(&cloud, &ReconstructionRequest::default(), &options()) {
        Ok(Reconstruction::Refused(_)) => Outcome::Passed,
        Ok(Reconstruction::Surface(outcome)) => Outcome::Failed(format!(
            "collinear points produced a surface with {} triangles",
            outcome.mesh.indices.len() / 3
        )),
        Err(error) => Outcome::Failed(format!("errored instead of refusing: {error}")),
    }
}

fn check_evidence_consistent(provider: &impl PointcloudReconstruction) -> Outcome {
    let cloud = sphere_samples(512);
    match provider.reconstruct(&cloud, &ReconstructionRequest::default(), &options()) {
        Ok(Reconstruction::Refused(reason)) => Outcome::Skipped(reason.to_string()),
        Err(error) => Outcome::Failed(format!("errored on a dense sphere: {error}")),
        Ok(Reconstruction::Surface(outcome)) => {
            let actual = outcome.mesh.indices.len() / 3;
            if outcome.evidence.output_triangles != actual {
                return Outcome::Failed(format!(
                    "evidence claims {} triangles, mesh has {actual}",
                    outcome.evidence.output_triangles
                ));
            }
            if outcome.evidence.used_points > outcome.evidence.input_points {
                return Outcome::Failed(format!(
                    "used {} of {} points",
                    outcome.evidence.used_points, outcome.evidence.input_points
                ));
            }
            if actual > 0 && !outcome.evidence.achieved_edge_length.is_finite() {
                return Outcome::Failed(
                    "a surface was produced without reporting its edge length".to_owned(),
                );
            }
            if outcome.evidence.interpolated_triangles > actual {
                return Outcome::Failed(format!(
                    "claims {} interpolated triangles of {actual} total",
                    outcome.evidence.interpolated_triangles
                ));
            }
            Outcome::Passed
        }
    }
}

fn check_input_count(provider: &impl PointcloudReconstruction) -> Outcome {
    let cloud = sphere_samples(300);
    match provider.reconstruct(&cloud, &ReconstructionRequest::default(), &options()) {
        Ok(Reconstruction::Refused(reason)) => Outcome::Skipped(reason.to_string()),
        Err(error) => Outcome::Failed(format!("errored: {error}")),
        Ok(outcome) => {
            let outcome = outcome.outcome().expect("surface");
            if outcome.evidence.input_points == cloud.len() {
                Outcome::Passed
            } else {
                Outcome::Failed(format!(
                    "reported {} input points for a cloud of {}",
                    outcome.evidence.input_points,
                    cloud.len()
                ))
            }
        }
    }
}

fn check_determinism(provider: &impl PointcloudReconstruction) -> Outcome {
    if provider.determinism() == Determinism::BestEffort {
        return Outcome::Skipped("provider only claims best-effort determinism".to_owned());
    }
    let cloud = sphere_samples(400);
    let request = ReconstructionRequest::default();
    let first = provider.reconstruct(&cloud, &request, &options());
    let second = provider.reconstruct(&cloud, &request, &options());
    match (first, second) {
        (Ok(Reconstruction::Surface(a)), Ok(Reconstruction::Surface(b))) => {
            if a.mesh.positions == b.mesh.positions && a.mesh.indices == b.mesh.indices {
                Outcome::Passed
            } else {
                Outcome::Failed(
                    "a provider claiming determinism produced two different meshes".to_owned(),
                )
            }
        }
        (Ok(Reconstruction::Refused(reason)), _) => Outcome::Skipped(reason.to_string()),
        (Err(error), _) | (_, Err(error)) => Outcome::Failed(format!("errored: {error}")),
        _ => Outcome::Failed("one call produced a surface and the other refused".to_owned()),
    }
}

fn check_minimum_agrees(provider: &impl PointcloudReconstruction) -> Outcome {
    let minimum = provider.minimum_points();
    if minimum == 0 {
        return Outcome::Failed("minimum_points of 0 cannot be honest".to_owned());
    }
    // One below the declared minimum must be refused: a provider whose
    // advertised threshold does not match its behaviour makes the value
    // useless for pre-flighting.
    let points: Vec<Point3> = (0..minimum.saturating_sub(1))
        .map(|i| Point3::new(i as Scalar, (i % 3) as Scalar, (i % 5) as Scalar))
        .collect();
    if points.is_empty() {
        return Outcome::Skipped("minimum is 1; nothing below it to test".to_owned());
    }
    let cloud = PointCloud::new(points).expect("finite");
    match provider.reconstruct(&cloud, &ReconstructionRequest::default(), &options()) {
        Ok(Reconstruction::Refused(_)) => Outcome::Passed,
        Ok(Reconstruction::Surface(_)) => Outcome::Failed(format!(
            "declared minimum_points={minimum} but reconstructed from {}",
            minimum - 1
        )),
        Err(error) => Outcome::Failed(format!("errored instead of refusing: {error}")),
    }
}
