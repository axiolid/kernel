#![forbid(unsafe_code)]
//! Shared adversarial and degenerate geometry fixtures.
//!
//! # Why a crate rather than a directory of files
//!
//! A degenerate case is usually a NUMBER, not a file: a 2e-9 plane tilt, a
//! sliver a fraction of a millimetre wide, two vertices that coincide. Storing
//! those as mesh files invites silent corruption -- an exporter rounds a
//! coordinate and the fixture stops being degenerate while still passing.
//!
//! Constructing them in code keeps the exact bit pattern under version
//! control and makes the reproduction steps the fixture itself.
//!
//! # Provenance
//!
//! Every fixture carries a [`Provenance`] naming where it came from and under
//! what licence. Fixtures here are ORIGINAL: constructed from published bug
//! descriptions and geometric first principles, not copied from any corpus.
//! That keeps the licence question trivial and the repository redistributable.
//!
//! # Adding a fixture
//!
//! Add a constructor returning [`Fixture`], fill in every [`Provenance`]
//! field, and state in `expectation` what an implementation must do -- not what
//! it currently does. A fixture that records present behaviour cannot detect a
//! regression, because the regression becomes the new expectation.

use axiolid_core::Point3;
use axiolid_mesh::TriMesh;

/// Where a fixture came from and under what terms it may be redistributed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Provenance {
    /// Where the case came from: an issue, a specification, or first
    /// principles. Specific enough that a reader can go and check it.
    pub source: &'static str,
    /// Licence covering redistribution of this fixture's data.
    ///
    /// Every fixture so far is original work under the repository licence. A
    /// fixture taken from an external corpus must name that corpus's licence
    /// here.
    pub licence: &'static str,
    /// What an implementation must do with it, stated as a requirement.
    pub expectation: &'static str,
}

/// A named mesh fixture with its provenance.
#[derive(Debug, Clone, PartialEq)]
pub struct Fixture {
    /// Stable identifier, usable in a test name or a failure message.
    pub name: &'static str,
    /// The geometry.
    pub mesh: TriMesh,
    /// Where it came from and what it demands.
    pub provenance: Provenance,
}

/// A well-formed unit cube. The control: every operation must handle it.
#[must_use]
pub fn unit_cube() -> Fixture {
    Fixture {
        name: "unit_cube",
        mesh: box_mesh(1.0, 1.0, 1.0),
        provenance: Provenance {
            source: "First principles: the simplest closed two-manifold solid.",
            licence: "Original work, same licence as the repository.",
            expectation: "Volume 1, closed, every operation succeeds.",
        },
    }
}

/// A sliver triangle: three nearly collinear vertices.
///
/// Area is ~5e-11, far below any sane linear tolerance squared, so a normal
/// computed by cross product is dominated by rounding.
#[must_use]
pub fn sliver_triangle() -> Fixture {
    Fixture {
        name: "sliver_triangle",
        mesh: TriMesh::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.5, 1.0e-10, 0.0),
            ],
            vec![0, 1, 2],
        ),
        provenance: Provenance {
            source: "First principles: the classic degenerate-normal case.",
            licence: "Original work, same licence as the repository.",
            expectation: "Refuse as degenerate, or handle without producing NaN.",
        },
    }
}

/// A triangle with two coincident vertices: zero area, not merely small.
#[must_use]
pub fn duplicate_vertex_triangle() -> Fixture {
    Fixture {
        name: "duplicate_vertex_triangle",
        mesh: TriMesh::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            vec![0, 1, 2],
        ),
        provenance: Provenance {
            source: "First principles: exact degeneracy, no tolerance can rescue it.",
            licence: "Original work, same licence as the repository.",
            expectation: "Refuse as degenerate. It has no normal and no area.",
        },
    }
}

/// An open box: the top face is missing, so the shell is not closed.
///
/// Volume is undefined for an open shell. A provider that reports a plausible
/// number here is guessing, which is exactly the failure this fixture catches.
#[must_use]
pub fn open_shell() -> Fixture {
    let mut mesh = box_mesh(1.0, 1.0, 1.0);
    // Drop the last two triangles: one face of the cube.
    mesh.indices.truncate(mesh.indices.len() - 6);
    Fixture {
        name: "open_shell",
        mesh,
        provenance: Provenance {
            source: "First principles: volume needs a closed boundary.",
            licence: "Original work, same licence as the repository.",
            expectation: "Report not-closed; refuse volume rather than guess.",
        },
    }
}

/// Two cubes whose sizes differ by nine orders of magnitude.
///
/// Catches absolute epsilons: a tolerance tuned for the metre-scale box
/// swallows the nanometre box whole, and a bounds routine that adds the two
/// loses the small one to rounding entirely.
#[must_use]
pub fn scale_disparity() -> Fixture {
    let mut mesh = box_mesh(1.0e3, 1.0e3, 1.0e3);
    let small = box_mesh(1.0e-6, 1.0e-6, 1.0e-6);
    let base = u32::try_from(mesh.positions.len()).expect("fixture is small");
    mesh.positions.extend(small.positions.iter().copied());
    mesh.indices
        .extend(small.indices.iter().map(|index| index + base));
    Fixture {
        name: "scale_disparity",
        mesh,
        provenance: Provenance {
            source: "First principles: absolute epsilons fail across scales.",
            licence: "Original work, same licence as the repository.",
            expectation: "Bounds must contain both boxes; no relative feature is lost.",
        },
    }
}

/// The ADR 0014 near-degenerate half-space column, in millimetres.
///
/// A 250 x 250 x 11940 mm column clipped by a plane tilted 2e-9 off axis. The
/// historical failure was a flyaway: output escaping the input bounds by
/// kilometres. See `docs/adr/0014-adopt-boolmesh-mesh-boolean.md`.
#[must_use]
pub fn millimetre_column() -> Fixture {
    Fixture {
        name: "millimetre_column",
        mesh: column_mesh(),
        provenance: Provenance {
            source: "ADR 0014, upstream issue 1155: a half-space clip flyaway.",
            licence: "Original reconstruction from the published bug description.",
            expectation: "Clipped output stays within the column bounds, or refuses.",
        },
    }
}

/// Two cubes sharing exactly one face plane: coplanar boolean contact.
///
/// Coplanar faces are the hardest boolean case: the classifier must decide
/// whether the shared plane is inside, outside, or on the boundary, and a
/// tolerance-based answer flips with the ordering of the operands.
#[must_use]
pub fn coplanar_contact() -> (Fixture, Fixture) {
    let mut second = box_mesh(1.0, 1.0, 1.0);
    for position in &mut second.positions {
        position.x += 1.0;
    }
    (
        Fixture {
            name: "coplanar_contact_left",
            mesh: box_mesh(1.0, 1.0, 1.0),
            provenance: COPLANAR,
        },
        Fixture {
            name: "coplanar_contact_right",
            mesh: second,
            provenance: COPLANAR,
        },
    )
}

const COPLANAR: Provenance = Provenance {
    source: "First principles: the canonical coplanar-boolean ambiguity.",
    licence: "Original work, same licence as the repository.",
    expectation: "Union volume is 2 exactly; intersection has zero volume.",
};

/// Every single-mesh fixture in the corpus.
///
/// Iterating this is how a differential test picks up new fixtures without
/// being edited: add a constructor here and every consumer covers it.
#[must_use]
pub fn corpus() -> Vec<Fixture> {
    vec![
        unit_cube(),
        sliver_triangle(),
        duplicate_vertex_triangle(),
        open_shell(),
        scale_disparity(),
        millimetre_column(),
    ]
}

/// An axis-aligned box with a corner at the origin, wound outward.
#[must_use]
pub fn box_mesh(sx: f64, sy: f64, sz: f64) -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(sx, 0.0, 0.0),
            Point3::new(sx, sy, 0.0),
            Point3::new(0.0, sy, 0.0),
            Point3::new(0.0, 0.0, sz),
            Point3::new(sx, 0.0, sz),
            Point3::new(sx, sy, sz),
            Point3::new(0.0, sy, sz),
        ],
        vec![
            0, 2, 1, 0, 3, 2, // bottom, wound outward (downward)
            4, 5, 6, 4, 6, 7, // top
            0, 1, 5, 0, 5, 4, // front
            1, 2, 6, 1, 6, 5, // right
            2, 3, 7, 2, 7, 6, // back
            3, 0, 4, 3, 4, 7, // left
        ],
    )
}

/// The ADR 0014 column: 250 x 250 mm in plan, from z = 11940 to z = 23880.
fn column_mesh() -> TriMesh {
    let mut mesh = box_mesh(250.0, 250.0, 11_940.0);
    for position in &mut mesh.positions {
        position.x -= 125.0;
        position.y -= 125.0;
        position.z += 11_940.0;
    }
    mesh
}
