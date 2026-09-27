//! Torus and wedge primitives (#142, ledger row C3) against closed forms.
//!
//! The torus mesh is a grid of planar trapezoids, n steps round the axis by
//! m round the tube, so its volume has a closed form of its own: each of the
//! n sectors is the tube polygon swept between two half-planes with a
//! constant Jacobian `sin(2 pi / n) rho`, giving
//! `V = n sin(2 pi / n) * R * (m / 2) r^2 sin(2 pi / m)`. That is exact for
//! the mesh, and tends to `2 pi^2 R r^2` as n and m grow. The wedge's faces
//! are all planar, so its volume is the prismatoid formula exactly.

use std::collections::{HashMap, HashSet};

use axiolid_core::{Scalar, Tolerance};
use axiolid_measure::{surface_properties, volume_properties};
use axiolid_mesh::TriMesh;
use axiolid_primitive::Primitive;
use axiolid_reference::primitive::tessellate_primitive;
use core::f64::consts::{PI, TAU};

fn tol_for(chord: Scalar) -> Tolerance {
    Tolerance::new((chord * 1e-3).max(1e-12), 1e-9).expect("tolerance")
}

fn torus(major_radius: Scalar, minor_radius: Scalar) -> Primitive {
    Primitive::Torus {
        major_radius,
        minor_radius,
    }
}

fn wedge(
    x: Scalar,
    y: Scalar,
    height: Scalar,
    top_x: [Scalar; 2],
    top_y: [Scalar; 2],
) -> Primitive {
    Primitive::Wedge {
        x,
        y,
        height,
        top_x_min: top_x[0],
        top_x_max: top_x[1],
        top_y_min: top_y[0],
        top_y_max: top_y[1],
    }
}

/// Closed, consistently oriented, one component; returns V - E + F.
fn euler_characteristic(mesh: &TriMesh) -> i64 {
    let mut directed: HashSet<(u32, u32)> = HashSet::new();
    for t in mesh.indices.chunks_exact(3) {
        assert!(
            t[0] != t[1] && t[1] != t[2] && t[2] != t[0],
            "degenerate {t:?}"
        );
        for k in 0..3 {
            let edge = (t[k], t[(k + 1) % 3]);
            assert!(directed.insert(edge), "edge {edge:?} used twice one way");
        }
    }
    for &(a, b) in &directed {
        assert!(directed.contains(&(b, a)), "edge ({a}, {b}) is a boundary");
    }
    // One component: flood the vertex graph.
    let mut adjacent: HashMap<u32, Vec<u32>> = HashMap::new();
    for &(a, b) in &directed {
        adjacent.entry(a).or_default().push(b);
    }
    let mut seen = HashSet::from([mesh.indices[0]]);
    let mut stack = vec![mesh.indices[0]];
    while let Some(v) = stack.pop() {
        for &w in &adjacent[&v] {
            if seen.insert(w) {
                stack.push(w);
            }
        }
    }
    assert_eq!(
        seen.len(),
        mesh.positions.len(),
        "one component, no stray vertices"
    );
    let (v, e, f) = (
        mesh.positions.len() as i64,
        directed.len() as i64 / 2,
        mesh.indices.len() as i64 / 3,
    );
    v - e + f
}

/// Grid counts (round the axis, round the tube) read back from the mesh:
/// the vertices at angle 0 round the axis are the one tube ring with y = 0
/// and x > 0.
fn grid(mesh: &TriMesh) -> (usize, usize) {
    let m = mesh
        .positions
        .iter()
        .filter(|p| p.y == 0.0 && p.x > 0.0)
        .count();
    assert_eq!(mesh.positions.len() % m, 0);
    (mesh.positions.len() / m, m)
}

#[test]
fn a_torus_mesh_has_its_trapezoid_grids_exact_volume_and_area() {
    let (big, r) = (3.0, 1.0);
    let chord = 1e-3;
    let mesh = tessellate_primitive(&torus(big, r), Tolerance::new(chord, 1e-9).unwrap()).unwrap();
    let (n, m) = grid(&mesh);
    assert!(n > m && m >= 3, "n {n}, m {m}");
    assert_eq!(mesh.indices.len(), 6 * n * m);
    // Every vertex is on the torus.
    for p in &mesh.positions {
        let rho = (p.x * p.x + p.y * p.y).sqrt();
        let off = ((rho - big).powi(2) + p.z * p.z).sqrt() - r;
        assert!(off.abs() < 1e-12, "{p:?} is {off} off the torus");
    }
    // A torus has genus 1.
    assert_eq!(euler_characteristic(&mesh), 0);

    let (nf, mf) = (n as Scalar, m as Scalar);
    let exact = nf * (TAU / nf).sin() * big * 0.5 * mf * r * r * (TAU / mf).sin();
    let v = volume_properties(&mesh, tol_for(chord))
        .unwrap()
        .signed_volume;
    assert!((v - exact).abs() <= 1e-12 * exact, "{v} vs {exact}");

    // Each cell's trapezoid: parallel chords 2 rho sin(pi / n) apart by
    // the distance between their lines.
    let (s, c) = ((PI / nf).sin(), (PI / nf).cos());
    let ring = |j: usize| {
        let phi = TAU * (j % m) as Scalar / mf;
        (big + r * phi.cos(), r * phi.sin())
    };
    let area: Scalar = (0..m)
        .map(|j| {
            let ((ra, za), (rb, zb)) = (ring(j), ring(j + 1));
            nf * s * (ra + rb) * ((zb - za).powi(2) + ((rb - ra) * c).powi(2)).sqrt()
        })
        .sum();
    let got = surface_properties(&mesh, tol_for(chord)).unwrap().area;
    assert!((got - area).abs() <= 1e-12 * area, "{got} vs {area}");
}

#[test]
fn a_torus_converges_on_two_pi_squared_r_r_squared() {
    let (big, r) = (5.0, 2.0);
    let want = 2.0 * PI * PI * big * r * r;
    let want_area = 4.0 * PI * PI * big * r;
    let mut previous = Scalar::INFINITY;
    for chord in [1e-2, 1e-3, 1e-4] {
        let mesh =
            tessellate_primitive(&torus(big, r), Tolerance::new(chord, 1e-9).unwrap()).unwrap();
        let v = volume_properties(&mesh, tol_for(chord))
            .unwrap()
            .signed_volume;
        let a = surface_properties(&mesh, tol_for(chord)).unwrap().area;
        // Inscribed: both fall short, by about the chord budget per unit
        // of surface for the volume.
        assert!(
            v < want && want - v <= 2.0 * chord * want_area,
            "{chord}: {v} vs {want}"
        );
        assert!(a < want_area, "{chord}: {a} vs {want_area}");
        let error = want - v;
        assert!(error < previous, "error must shrink: {error} vs {previous}");
        previous = error;
    }
}

#[test]
fn a_torus_faces_outward_everywhere() {
    // Positive volume alone would pass with a few flipped cells; every
    // triangle's normal must point away from the tube's centre circle.
    let (big, r) = (3.0, 1.0);
    let mesh = tessellate_primitive(&torus(big, r), Tolerance::new(1e-2, 1e-9).unwrap()).unwrap();
    for t in mesh.indices.chunks_exact(3) {
        let [a, b, c] = [0, 1, 2].map(|k| mesh.positions[t[k] as usize]);
        let centre = (a + b + c) / 3.0;
        let rho = (centre.x * centre.x + centre.y * centre.y).sqrt();
        let spine = axiolid_core::Point3::new(centre.x / rho * big, centre.y / rho * big, 0.0);
        assert!((b - a).cross(c - a).dot(centre - spine) > 0.0);
    }
}

#[test]
fn horn_and_spindle_tori_and_bad_radii_are_refused() {
    for (big, r, needle) in [
        (1.0, 1.0, "horn torus"),
        (1.0, 1.5, "spindle torus"),
        (0.0, 1.0, "major radius"),
        (2.0, -1.0, "minor radius"),
        (f64::NAN, 1.0, "major radius"),
        (2.0, f64::INFINITY, "minor radius"),
    ] {
        match tessellate_primitive(&torus(big, r), Tolerance::MILLIMETRE) {
            Err(axiolid_contracts::GeomError::InvalidInput(message)) => {
                assert!(message.contains(needle), "{big}, {r}: {message}");
            }
            other => panic!("{big}, {r}: {other:?}"),
        }
    }
}

/// Prismatoid: `h / 6 (A0 + 4 Am + A1)`, exact for planar faces between
/// parallel rectangles.
fn prismatoid(x: Scalar, y: Scalar, h: Scalar, top_x: [Scalar; 2], top_y: [Scalar; 2]) -> Scalar {
    let (tx, ty) = (top_x[1] - top_x[0], top_y[1] - top_y[0]);
    h / 6.0 * (x * y + (x + tx) * (y + ty) + tx * ty)
}

#[test]
fn an_occt_ltx_wedge_has_its_closed_form_volume_and_faces() {
    // OCCT MakeWedge(dx = 4, dy = 2 (our height), dz = 3 (our y), ltx = 1).
    let (x, y, h, ltx) = (4.0, 3.0, 2.0, 1.0);
    let mesh =
        tessellate_primitive(&wedge(x, y, h, [0.0, ltx], [0.0, y]), Tolerance::MILLIMETRE).unwrap();
    assert_eq!(euler_characteristic(&mesh), 2);
    assert_eq!(mesh.positions.len(), 8);
    assert_eq!(mesh.indices.len(), 36, "six quads");
    let v = volume_properties(&mesh, Tolerance::MILLIMETRE)
        .unwrap()
        .signed_volume;
    // OCCT's ltx form: dy dz (dx + ltx) / 2.
    assert_eq!(v, h * y * (x + ltx) / 2.0);
    // Faces: base 12, top 3, two trapezoids of 5, the upright side 6 and
    // the slanted one 3 x sqrt(3^2 + 2^2).
    let area = surface_properties(&mesh, Tolerance::MILLIMETRE)
        .unwrap()
        .area;
    let want = 12.0 + 3.0 + 5.0 + 5.0 + 6.0 + 3.0 * 13.0_f64.sqrt();
    assert!((area - want).abs() < 1e-12, "{area} vs {want}");
}

#[test]
fn a_general_wedge_matches_the_prismatoid_formula() {
    // A top that is narrower in both directions and shifted, part of it
    // overhanging the base.
    let (x, y, h, tx, ty) = (4.0, 3.0, 2.5, [1.0, 4.5], [0.5, 2.0]);
    let mesh = tessellate_primitive(&wedge(x, y, h, tx, ty), Tolerance::MILLIMETRE).unwrap();
    assert_eq!(euler_characteristic(&mesh), 2);
    let v = volume_properties(&mesh, Tolerance::MILLIMETRE)
        .unwrap()
        .signed_volume;
    let want = prismatoid(x, y, h, tx, ty);
    assert!((v - want).abs() < 1e-12, "{v} vs {want}");
}

#[test]
fn a_top_collapsed_to_an_edge_or_a_point_stays_a_closed_solid() {
    // A wedge proper: the top is a ridge along y at x = 2 (a gable).
    let ridge = tessellate_primitive(
        &wedge(4.0, 3.0, 2.0, [2.0, 2.0], [0.0, 3.0]),
        Tolerance::MILLIMETRE,
    )
    .unwrap();
    assert_eq!(euler_characteristic(&ridge), 2);
    assert_eq!(ridge.positions.len(), 6);
    // Base, two roof quads, two gable triangles.
    assert_eq!(ridge.indices.len() / 3, 2 + 2 + 2 + 1 + 1);
    let v = volume_properties(&ridge, Tolerance::MILLIMETRE)
        .unwrap()
        .signed_volume;
    assert_eq!(v, 0.5 * 4.0 * 2.0 * 3.0);

    // A point: a pyramid with its apex over a corner.
    let apex = tessellate_primitive(
        &wedge(4.0, 3.0, 2.0, [0.0, 0.0], [0.0, 0.0]),
        Tolerance::MILLIMETRE,
    )
    .unwrap();
    assert_eq!(euler_characteristic(&apex), 2);
    assert_eq!(apex.positions.len(), 5);
    // The two sides through the apex's corner are triangles, as are the
    // other two; with the base, six triangles.
    assert_eq!(apex.indices.len() / 3, 6);
    let v = volume_properties(&apex, Tolerance::MILLIMETRE)
        .unwrap()
        .signed_volume;
    assert!((v - 4.0 * 3.0 * 2.0 / 3.0).abs() < 1e-12, "{v}");
}

#[test]
fn degenerate_wedges_are_refused() {
    for (primitive, needle) in [
        (wedge(0.0, 3.0, 2.0, [0.0, 1.0], [0.0, 3.0]), "wedge x"),
        (wedge(4.0, -3.0, 2.0, [0.0, 1.0], [0.0, 3.0]), "wedge y"),
        (wedge(4.0, 3.0, 0.0, [0.0, 1.0], [0.0, 3.0]), "wedge height"),
        (wedge(4.0, 3.0, 2.0, [1.0, 0.5], [0.0, 3.0]), "min <= max"),
        (wedge(4.0, 3.0, 2.0, [0.0, 1.0], [2.0, 1.0]), "min <= max"),
        (
            wedge(4.0, 3.0, 2.0, [f64::NAN, 1.0], [0.0, 3.0]),
            "top x min",
        ),
        (
            wedge(4.0, 3.0, 2.0, [0.0, 1.0], [0.0, f64::INFINITY]),
            "top y max",
        ),
    ] {
        match tessellate_primitive(&primitive, Tolerance::MILLIMETRE) {
            Err(axiolid_contracts::GeomError::InvalidInput(message)) => {
                assert!(message.contains(needle), "{primitive:?}: {message}");
            }
            other => panic!("{primitive:?}: {other:?}"),
        }
    }
}
