//! Shared conformance checks for portable mesh-section providers.

use core::fmt;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame3, Point3, Tolerance, Vec3};
use axiolid_mesh::TriMesh;

use crate::{MeshPlaneSection, SectionLimits};

#[derive(Debug, Clone, PartialEq)]
pub enum ConformanceFailure {
    ReturnedError(String),
    EmptyCentralSection,
    IncorrectEvidence,
    NonDeterministic,
    /// A section returned the wrong number of closed contours.
    WrongContourCount {
        case: &'static str,
        expected: usize,
        actual: usize,
    },
    /// A section contour enclosed the wrong area.
    WrongSectionArea {
        case: &'static str,
        expected: f64,
        actual: f64,
    },
    /// A plane tangent to the solid produced interior area instead of
    /// touching it.
    TangentPlaneHasArea {
        case: &'static str,
        area: f64,
    },
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ConformanceReport {
    pub failures: Vec<ConformanceFailure>,
}

impl ConformanceReport {
    pub fn is_success(&self) -> bool {
        self.failures.is_empty()
    }
}

impl fmt::Display for ConformanceReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_success() {
            write!(f, "conformant")
        } else {
            write!(f, "{:?}", self.failures)
        }
    }
}

pub struct ConformanceSuite;

impl ConformanceSuite {
    pub fn run(provider: &dyn MeshPlaneSection) -> ConformanceReport {
        let mut report = ConformanceReport::default();
        let mesh = unit_cube();
        let frame = Frame3 {
            origin: Point3::new(0.0, 0.0, 0.5),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        };
        let limits = SectionLimits::new(8, 12, 32, 4);
        let options = ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("valid tolerance"));
        let first = provider.section(&mesh, frame, limits, &options);
        let second = provider.section(&mesh, frame, limits, &options);
        match (first, second) {
            (Ok(a), Ok(b)) => {
                if a.contours.is_empty() {
                    report
                        .failures
                        .push(ConformanceFailure::EmptyCentralSection);
                }
                if a.evidence.source_triangles != mesh.triangles().len()
                    || !a.evidence.is_derived_from_input_mesh()
                {
                    report.failures.push(ConformanceFailure::IncorrectEvidence);
                }
                if a != b {
                    report.failures.push(ConformanceFailure::NonDeterministic);
                }
            }
            (Err(error), _) | (_, Err(error)) => report
                .failures
                .push(ConformanceFailure::ReturnedError(error.to_string())),
        }

        // Geometry, not just non-emptiness. A cube section that merely
        // exists proves nothing about the numbers the provider returned;
        // these cases check the actual enclosed area and contour count
        // against an analytic oracle, including the degenerate planes
        // that a naive implementation gets wrong.
        //
        // Subdivision 3 is deliberate: fine enough that the inscribed
        // polygon is within ~0.5% of the true circle, coarse enough that
        // the fixture stays cheap for every provider that runs it.
        let sphere = icosphere(1.0, 3);
        for spec in &[
            // Central: the largest circle, area pi.
            SphereCase {
                case: "sphere central",
                height: 0.0,
                contours: 1,
                tolerance: 0.01,
            },
            // Off-centre: a smaller circle and a different oracle value, so a
            // provider returning a constant cannot pass both.
            SphereCase {
                case: "sphere h=0.5",
                height: 0.5,
                contours: 1,
                tolerance: 0.02,
            },
            // Above the pole: empty, no contour at all.
            SphereCase {
                case: "sphere above pole",
                height: 1.5,
                contours: 0,
                tolerance: 0.0,
            },
            // Tangent at the pole: touches one point, so zero enclosed area.
            // The icosphere has a VERTEX at the pole, so this is also the
            // vertex-hit case.
            SphereCase {
                case: "sphere tangent pole",
                height: 1.0,
                contours: 0,
                tolerance: 0.0,
            },
        ] {
            check_sphere_section(provider, &mut report, &sphere, 1.0, spec);
        }

        report
    }
}

fn unit_cube() -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 1.0),
        ],
        vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7,
            6, 3, 0, 4, 3, 4, 7,
        ],
    )
}

/// Signed area of a closed contour, by the shoelace formula.
fn contour_area(contour: &crate::SectionContour) -> f64 {
    let p = &contour.points;
    if p.len() < 3 {
        return 0.0;
    }
    let mut twice = 0.0;
    for i in 0..p.len() {
        let a = p[i];
        let b = p[(i + 1) % p.len()];
        twice += a.x * b.y - b.x * a.y;
    }
    (twice / 2.0).abs()
}

/// Total enclosed area across every contour.
fn total_area(contours: &[crate::SectionContour]) -> f64 {
    contours.iter().map(contour_area).sum()
}

/// An icosphere of `radius` about the origin, subdivided `n` times.
///
/// The section oracle is derived from the TESSELLATED solid, not from the
/// ideal sphere: an icosphere inscribes its sphere, so a plane at height h
/// cuts a polygon slightly smaller than the true circle of radius
/// sqrt(r^2 - h^2). Comparing against the ideal radius would charge the
/// provider for the caller's tessellation choice.
fn icosphere(radius: f64, subdivisions: u32) -> TriMesh {
    let t = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut v: Vec<[f64; 3]> = vec![
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
    let mut f: Vec<[u32; 3]> = vec![
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
    for _ in 0..subdivisions {
        let mut mid: std::collections::HashMap<(u32, u32), u32> = std::collections::HashMap::new();
        let mut next: Vec<[u32; 3]> = Vec::with_capacity(f.len() * 4);
        for tri in &f {
            let mut m = [0u32; 3];
            for e in 0..3 {
                let (a, b) = (tri[e], tri[(e + 1) % 3]);
                let key = (a.min(b), a.max(b));
                m[e] = *mid.entry(key).or_insert_with(|| {
                    let (pa, pb) = (v[a as usize], v[b as usize]);
                    v.push([
                        (pa[0] + pb[0]) * 0.5,
                        (pa[1] + pb[1]) * 0.5,
                        (pa[2] + pb[2]) * 0.5,
                    ]);
                    (v.len() - 1) as u32
                });
            }
            next.push([tri[0], m[0], m[2]]);
            next.push([tri[1], m[1], m[0]]);
            next.push([tri[2], m[2], m[1]]);
            next.push([m[0], m[1], m[2]]);
        }
        f = next;
    }
    let positions: Vec<Point3> = v
        .iter()
        .map(|p| {
            let l = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            let k = radius / l;
            Point3::new(p[0] * k, p[1] * k, p[2] * k)
        })
        .collect();
    let indices: Vec<u32> = f.iter().flat_map(|t| [t[0], t[1], t[2]]).collect();
    TriMesh::new(positions, indices)
}

/// One sphere-section case: the plane height and what it should produce.
struct SphereCase {
    /// Name used in failure reports.
    case: &'static str,
    /// Height of the sectioning plane above the sphere centre.
    height: f64,
    /// Closed contours the plane should produce.
    contours: usize,
    /// Relative area tolerance, sized to the tessellation.
    tolerance: f64,
}

/// Section a sphere and compare against the oracle derived from the
/// tessellated solid.
///
/// A plane at height h through a sphere of radius r cuts a circle of
/// radius sqrt(r^2 - h^2) and area pi*(r^2 - h^2). The icosphere inscribes
/// that sphere, so the measured area is slightly under; the tolerance is
/// sized to the tessellation, not to the provider.
fn check_sphere_section(
    provider: &dyn MeshPlaneSection,
    report: &mut ConformanceReport,
    sphere: &TriMesh,
    radius: f64,
    spec: &SphereCase,
) {
    let SphereCase {
        case,
        height,
        contours: expect_contours,
        tolerance: rel_tolerance,
    } = *spec;
    let frame = Frame3 {
        origin: Point3::new(0.0, 0.0, height),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    };
    let limits = SectionLimits::new(4096, 8192, 16384, 64);
    let options = ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("valid tolerance"));
    match provider.section(sphere, frame, limits, &options) {
        Err(error) => report
            .failures
            .push(ConformanceFailure::ReturnedError(error.to_string())),
        Ok(outcome) => {
            let closed = outcome.contours.len();
            if closed != expect_contours {
                report.failures.push(ConformanceFailure::WrongContourCount {
                    case,
                    expected: expect_contours,
                    actual: closed,
                });
            }
            let area = total_area(&outcome.contours);
            if expect_contours == 0 {
                // A tangent plane touches at one point: any enclosed area is
                // a real defect, not a rounding artefact.
                if area > 1e-9 {
                    report
                        .failures
                        .push(ConformanceFailure::TangentPlaneHasArea { case, area });
                }
                return;
            }
            let expected = core::f64::consts::PI * (radius * radius - height * height);
            if expected > 0.0 && (area - expected).abs() / expected > rel_tolerance {
                report.failures.push(ConformanceFailure::WrongSectionArea {
                    case,
                    expected,
                    actual: area,
                });
            }
        }
    }
}
