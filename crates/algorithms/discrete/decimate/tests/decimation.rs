//! Decimation preserves what it claims to preserve (#75).
//!
//! The oracles are independent of the decimator: volume via the v0.7
//! measurement provider, and defect classes via the v0.7 diagnosis. A
//! decimator that damages a mesh must be caught by the same machinery that
//! finds damage from any other source.

use axiolid_core::{Point3, Tolerance};
use axiolid_decimate::{decimate, DecimateTarget};
use axiolid_heal::mesh::MeshHealer;
use axiolid_heal::{DefectKind, Diagnose};
use axiolid_measure::volume_properties;
use axiolid_mesh::{audit_mesh, TriMesh};

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("tolerance")
}

/// An icosphere-ish closed mesh: subdivided octahedron, normalised.
///
/// Dense enough that decimation has real work to do, and closed so volume
/// is defined.
fn sphere(subdivisions: usize) -> TriMesh {
    let mut positions: Vec<Point3> = vec![
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, -1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(0.0, 0.0, -1.0),
    ];
    let mut faces: Vec<[u32; 3]> = vec![
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];

    for _ in 0..subdivisions {
        let mut next = Vec::with_capacity(faces.len() * 4);
        let mut cache: std::collections::BTreeMap<(u32, u32), u32> =
            std::collections::BTreeMap::new();
        for f in &faces {
            let mut mid = |a: u32, b: u32| -> u32 {
                let key = (a.min(b), a.max(b));
                if let Some(&hit) = cache.get(&key) {
                    return hit;
                }
                let p = (positions[a as usize] + positions[b as usize]) / 2.0;
                let n = p / p.length();
                let index = u32::try_from(positions.len()).expect("fits");
                positions.push(n);
                cache.insert(key, index);
                index
            };
            let (a, b, c) = (f[0], f[1], f[2]);
            let (ab, bc, ca) = (mid(a, b), mid(b, c), mid(c, a));
            next.push([a, ab, ca]);
            next.push([ab, b, bc]);
            next.push([ca, bc, c]);
            next.push([ab, bc, ca]);
        }
        faces = next;
    }

    let indices = faces.iter().flat_map(|f| f.iter().copied()).collect();
    TriMesh::new(positions, indices)
}

#[test]
fn decimation_reduces_and_reports_its_deviation() {
    let mesh = sphere(3);
    let before = mesh.indices.len() / 3;

    let (out, report) = decimate(&mesh, DecimateTarget::MaxDeviation(0.15), tol())
        .expect("a closed sphere is decimatable");

    assert!(
        report.collapses > 0,
        "nothing collapsed on a dense sphere: {report:?}"
    );
    assert!(
        out.indices.len() / 3 < before,
        "triangle count did not fall: {before} -> {}",
        out.indices.len() / 3
    );
    assert!(
        report.max_deviation <= 0.15,
        "reported deviation {} exceeds the bound it promised",
        report.max_deviation
    );
}

#[test]
fn the_result_carries_no_new_defects() {
    let mesh = sphere(3);
    let (out, report) =
        decimate(&mesh, DecimateTarget::MaxDeviation(0.15), tol()).expect("decimatable");
    assert!(
        report.collapses > 0,
        "the test proves nothing if nothing ran"
    );

    // The v0.7 diagnosis is the oracle: a decimator that inverts triangles
    // or opens the shell is caught by the same machinery that finds damage
    // from any other source.
    let diagnosis = MeshHealer.diagnose(&out, tol()).expect("diagnosable");
    let classes: Vec<DefectKind> = diagnosis.defects.iter().map(|d| d.kind).collect();
    assert!(
        !classes.contains(&DefectKind::SelfIntersection),
        "decimation introduced self-intersection: {classes:?}"
    );
    assert!(
        !classes.contains(&DefectKind::NonManifoldEdge),
        "decimation introduced a non-manifold edge: {classes:?}"
    );
    assert!(
        !classes.contains(&DefectKind::OpenShell),
        "decimation opened the shell: {classes:?}"
    );
    assert!(
        !classes.contains(&DefectKind::InconsistentOrientation),
        "decimation inverted a triangle: {classes:?}"
    );
}

#[test]
fn volume_is_preserved_within_the_stated_bound() {
    let mesh = sphere(3);
    let original = volume_properties(&mesh, tol())
        .expect("closed")
        .signed_volume;

    let bound = 0.15;
    let (out, report) =
        decimate(&mesh, DecimateTarget::MaxDeviation(bound), tol()).expect("decimatable");
    assert!(report.collapses > 0);

    let reduced = volume_properties(&out, tol())
        .expect("still closed")
        .signed_volume;
    // A surface perturbed by at most `bound` changes volume by at most
    // roughly area * bound. The sphere has area ~4*pi, so this is a
    // generous but non-vacuous ceiling derived from the bound itself
    // rather than from the observed answer.
    let ceiling = 4.0 * std::f64::consts::PI * bound;
    assert!(
        (original - reduced).abs() < ceiling,
        "volume moved {} which exceeds the {ceiling} implied by a {bound} bound",
        (original - reduced).abs()
    );
}

#[test]
fn repeated_runs_are_identical() {
    let mesh = sphere(2);
    let (a, ra) = decimate(&mesh, DecimateTarget::MaxDeviation(0.1), tol()).expect("ok");
    let (b, rb) = decimate(&mesh, DecimateTarget::MaxDeviation(0.1), tol()).expect("ok");
    assert_eq!(a.indices, b.indices, "index buffers differ between runs");
    assert_eq!(a.positions, b.positions, "positions differ between runs");
    assert_eq!(ra, rb, "reports differ between runs");
}

#[test]
fn an_impossible_bound_is_a_reported_noop() {
    let mesh = sphere(2);
    // Far below any edge length: every collapse must be refused, and the
    // mesh must come back untouched rather than damaged.
    let (out, report) =
        decimate(&mesh, DecimateTarget::MaxDeviation(1e-12), tol()).expect("valid request");
    assert!(
        report.is_noop(),
        "collapses happened under an impossible bound: {report:?}"
    );
    assert_eq!(out.indices.len(), mesh.indices.len());
    assert!(
        report.rejected_deviation > 0,
        "refusals must be reported, not silent: {report:?}"
    );
}

#[test]
fn a_ragged_index_buffer_is_refused() {
    let bad = TriMesh::new(vec![Point3::new(0.0, 0.0, 0.0)], vec![0, 0]);
    assert!(decimate(&bad, DecimateTarget::MaxDeviation(0.1), tol()).is_err());
}

/// An unsafe collapse is refused rather than performed.
///
/// The sphere fixture never rejects anything: on a smooth convex surface
/// every candidate collapse is safe, so `rejected_unsafe` stays 0 and the
/// guards are never exercised. This fixture triggers one deliberately.
///
/// A tall thin pyramid whose apex sits almost above one base corner. The
/// short apex edge is cheap to collapse, and its endpoints share three
/// neighbours rather than two, so the link condition refuses it: performing
/// it would weld the fan into a non-manifold configuration.
///
/// Verified by mutation: removing the link condition makes this test fail.
/// The normal-inversion guard alongside it is NOT verified by this fixture
/// -- the link condition rejects first, so the inversion branch never runs.
/// TODO(#201): add a fixture that passes the link condition and still flips a
/// normal.
#[test]
fn an_unsafe_collapse_is_refused_not_performed() {
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 0.0),
        Point3::new(4.0, 4.0, 0.0),
        Point3::new(0.0, 4.0, 0.0),
        // Apex almost exactly above corner 0: the edge 0-4 is short, so
        // collapsing it is cheap, and it drags the apex across the base.
        Point3::new(0.05, 0.05, 3.0),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, // base, wound downward
        0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4,
    ];
    let spike = TriMesh::new(positions, indices);

    let original = volume_properties(&spike, tol())
        .expect("closed")
        .signed_volume;

    // Generous enough to admit the short apex edge, and every other edge is
    // 4 units long so nothing else is a cheap candidate.
    let (out, report) =
        decimate(&spike, DecimateTarget::MaxDeviation(2.0), tol()).expect("valid request");

    assert!(
        report.rejected_unsafe > 0,
        "the inversion guard never fired, so this fixture proves nothing: {report:?}"
    );

    let health = audit_mesh(&out, tol());
    assert_eq!(
        health.inconsistent_winding_edges, 0,
        "decimation produced inconsistent winding: {health:?}"
    );
    if let Ok(props) = volume_properties(&out, tol()) {
        assert!(
            props.signed_volume > 0.0,
            "the result encloses negative volume, so a face was inverted \
             (input was {original}): {}",
            props.signed_volume
        );
    }
}
