//! Isotropic remeshing: edge lengths converge on the target, valences
//! improve, and topology, orientation, features and the surface are kept.

use std::f64::consts::{PI, TAU};

use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_mesh::{HalfedgeBuildError, HalfedgeMesh, TriMesh};
use axiolid_refine::remesh::{remesh, RemeshError, RemeshOptions, RemeshReport};

// ---- fixtures -------------------------------------------------------------

/// Latitude-longitude sphere, outward wound. The poles have valence
/// `slices`, far from the optimal six.
fn uv_sphere(radius: Scalar, stacks: u32, slices: u32) -> TriMesh {
    let mut positions = vec![Point3::new(0.0, 0.0, radius)];
    for i in 1..stacks {
        let theta = PI * Scalar::from(i) / Scalar::from(stacks);
        for j in 0..slices {
            let phi = TAU * Scalar::from(j) / Scalar::from(slices);
            positions.push(Point3::new(
                radius * theta.sin() * phi.cos(),
                radius * theta.sin() * phi.sin(),
                radius * theta.cos(),
            ));
        }
    }
    positions.push(Point3::new(0.0, 0.0, -radius));
    let south = positions.len() as u32 - 1;
    let ring = |i: u32, j: u32| 1 + (i - 1) * slices + j % slices;
    let mut indices = Vec::new();
    for j in 0..slices {
        indices.extend([0, ring(1, j), ring(1, j + 1)]);
        indices.extend([south, ring(stacks - 1, j + 1), ring(stacks - 1, j)]);
    }
    for i in 1..stacks - 1 {
        for j in 0..slices {
            let (a, b, c, d) = (
                ring(i, j),
                ring(i + 1, j),
                ring(i + 1, j + 1),
                ring(i, j + 1),
            );
            indices.extend([a, b, c, a, c, d]);
        }
    }
    TriMesh::new(positions, indices)
}

/// Torus around the z axis, outward wound, its quads split by diagonals
/// alternating in a checkerboard so valences are four and eight.
fn torus(major: Scalar, minor: Scalar, around: u32, tube: u32) -> TriMesh {
    let mut positions = Vec::new();
    for i in 0..around {
        let u = TAU * Scalar::from(i) / Scalar::from(around);
        for j in 0..tube {
            let v = TAU * Scalar::from(j) / Scalar::from(tube);
            let r = major + minor * v.cos();
            positions.push(Point3::new(r * u.cos(), r * u.sin(), minor * v.sin()));
        }
    }
    let at = |i: u32, j: u32| (i % around) * tube + j % tube;
    let mut indices = Vec::new();
    for i in 0..around {
        for j in 0..tube {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            if (i + j) % 2 == 0 {
                indices.extend([a, b, c, a, c, d]);
            } else {
                indices.extend([a, b, d, b, c, d]);
            }
        }
    }
    TriMesh::new(positions, indices)
}

fn torus_normal(major: Scalar, p: Point3) -> Vec3 {
    let radial = Vec3::new(p.x, p.y, 0.0).normalize() * major;
    (p - radial).normalize()
}

/// The unit square in the xy plane, fanned from its centre to `n` points
/// per side: one vertex of valence `4 n`, every border vertex of valence 3.
fn fan_square(n: u32) -> TriMesh {
    let mut positions = vec![Point3::new(0.5, 0.5, 0.0)];
    let step = 1.0 / Scalar::from(n);
    for k in 0..n {
        positions.push(Point3::new(Scalar::from(k) * step, 0.0, 0.0));
    }
    for k in 0..n {
        positions.push(Point3::new(1.0, Scalar::from(k) * step, 0.0));
    }
    for k in 0..n {
        positions.push(Point3::new(1.0 - Scalar::from(k) * step, 1.0, 0.0));
    }
    for k in 0..n {
        positions.push(Point3::new(0.0, 1.0 - Scalar::from(k) * step, 0.0));
    }
    let border = 4 * n;
    let mut indices = Vec::new();
    for k in 0..border {
        indices.extend([0, 1 + k, 1 + (k + 1) % border]);
    }
    TriMesh::new(positions, indices)
}

/// Box `[0, sx] x [0, sy] x [0, sz]`, outward wound, two triangles a face.
fn box_mesh(sx: Scalar, sy: Scalar, sz: Scalar) -> TriMesh {
    let p = |x: Scalar, y: Scalar, z: Scalar| Point3::new(x * sx, y * sy, z * sz);
    let positions = vec![
        p(0.0, 0.0, 0.0),
        p(1.0, 0.0, 0.0),
        p(1.0, 1.0, 0.0),
        p(0.0, 1.0, 0.0),
        p(0.0, 0.0, 1.0),
        p(1.0, 0.0, 1.0),
        p(1.0, 1.0, 1.0),
        p(0.0, 1.0, 1.0),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, // bottom
        4, 5, 6, 4, 6, 7, // top
        0, 1, 5, 0, 5, 4, // front
        1, 2, 6, 1, 6, 5, // right
        2, 3, 7, 2, 7, 6, // back
        3, 0, 4, 3, 4, 7, // left
    ];
    TriMesh::new(positions, indices)
}

/// Closed cylinder of the given radius around the z axis from `z = 0` to
/// `height`, outward wound, each cap fanned from its centre. The two rims
/// are vertices `0..segments` (bottom) and `segments..2 segments` (top).
fn capped_cylinder(radius: Scalar, height: Scalar, segments: u32) -> TriMesh {
    let mut positions = Vec::new();
    for z in [0.0, height] {
        for j in 0..segments {
            let phi = TAU * Scalar::from(j) / Scalar::from(segments);
            positions.push(Point3::new(radius * phi.cos(), radius * phi.sin(), z));
        }
    }
    positions.push(Point3::new(0.0, 0.0, 0.0));
    positions.push(Point3::new(0.0, 0.0, height));
    let (bottom, top) = (2 * segments, 2 * segments + 1);
    let mut indices = Vec::new();
    for j in 0..segments {
        let k = (j + 1) % segments;
        let (a, b, c, d) = (j, k, segments + k, segments + j);
        indices.extend([a, b, c, a, c, d]);
        indices.extend([bottom, k, j]);
        indices.extend([top, segments + j, segments + k]);
    }
    TriMesh::new(positions, indices)
}

/// L-shaped prism: the L `(0,0) (2,0) (2,1) (1,1) (1,2) (0,2)` in the xy
/// plane extruded to `z = 1`, outward wound. Its edge at `(1, 1)` is a
/// concave (reflex) crease.
fn l_prism() -> TriMesh {
    let outline = [
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 1.0),
        (1.0, 1.0),
        (1.0, 2.0),
        (0.0, 2.0),
    ];
    let n = outline.len() as u32;
    let mut positions = Vec::new();
    for z in [0.0, 1.0] {
        for &(x, y) in &outline {
            positions.push(Point3::new(x, y, z));
        }
    }
    // The outline runs counter-clockwise seen from +z; each cap is a fan
    // from corner 0, which sees every other corner.
    let cap = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5]];
    let mut indices = Vec::new();
    for [a, b, c] in cap {
        indices.extend([a, c, b]); // bottom faces -z
        indices.extend([n + a, n + b, n + c]); // top faces +z
    }
    for j in 0..n {
        let k = (j + 1) % n;
        indices.extend([j, k, n + k, j, n + k, n + j]);
    }
    TriMesh::new(positions, indices)
}

// ---- measurements ---------------------------------------------------------

fn normal(mesh: &TriMesh, t: [u32; 3]) -> Vec3 {
    let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
    (b - a).cross(c - a)
}

fn centroid(mesh: &TriMesh, t: [u32; 3]) -> Point3 {
    let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
    (a + b + c) / 3.0
}

/// Structure of a mesh that remeshing must keep.
fn topology(mesh: &TriMesh) -> (i64, usize) {
    let he = HalfedgeMesh::from_tri_mesh(mesh).expect("manifold, consistently wound");
    he.validate().expect("halfedge invariants");
    (he.euler_characteristic(), he.boundary_loops().len())
}

fn edge_lengths(mesh: &TriMesh) -> Vec<Scalar> {
    let he = HalfedgeMesh::from_tri_mesh(mesh).expect("manifold");
    he.edges()
        .map(|e| {
            let [a, b] = he.edge_vertices(e);
            (he.position(a) - he.position(b)).length()
        })
        .collect()
}

fn band_fraction(mesh: &TriMesh, target: Scalar) -> Scalar {
    let lengths = edge_lengths(mesh);
    let inside = lengths
        .iter()
        .filter(|&&l| (0.8 * target..=4.0 / 3.0 * target).contains(&l))
        .count();
    inside as Scalar / lengths.len() as Scalar
}

/// Largest distance of a mesh from a smooth surface, sampled on a
/// barycentric grid of every triangle; `distance` is the surface's distance
/// function.
fn chord_error(mesh: &TriMesh, distance: impl Fn(Point3) -> Scalar) -> Scalar {
    const STEPS: u32 = 8;
    let mut worst: Scalar = 0.0;
    for t in mesh.triangles() {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        for i in 0..=STEPS {
            for j in 0..=STEPS - i {
                let u = Scalar::from(i) / Scalar::from(STEPS);
                let v = Scalar::from(j) / Scalar::from(STEPS);
                worst = worst.max(distance(a + (b - a) * u + (c - a) * v));
            }
        }
    }
    worst
}

fn hausdorff(a: &TriMesh, b: &TriMesh) -> Scalar {
    axiolid_measure::hausdorff_distance(a, b, 1e-5)
        .expect("measurable")
        .distance
        .upper
}

/// Run twice and require bit-identical output, then check what every
/// remeshing keeps: manifoldness, Euler characteristic, boundary loops, and
/// that edits only ever removed or added well-shaped triangles.
fn remesh_checked(mesh: &TriMesh, options: RemeshOptions) -> (TriMesh, RemeshReport) {
    let (out, report) = remesh(mesh, options).expect("remeshable");
    let (again, again_report) = remesh(mesh, options).expect("remeshable");
    assert_eq!(out, again, "remeshing is deterministic");
    assert_eq!(report, again_report);
    assert_eq!(topology(&out), topology(mesh), "topology preserved");
    assert_eq!(report.output_triangles, out.triangle_count());
    eprintln!("{report:#?}");
    (out, report)
}

// ---- convergence on closed and open surfaces ------------------------------

#[test]
fn sphere_converges_to_the_target_and_stays_on_the_sphere() {
    let radius = 1.0;
    let input = uv_sphere(radius, 12, 24);
    let target = 0.15;
    let (out, report) = remesh_checked(&input, RemeshOptions::new(target));

    let band = band_fraction(&out, target);
    assert!(band >= 0.9, "{band} of the edges in [4/5, 4/3] L");
    assert!((report.edges_in_band - band).abs() < 1e-12);
    assert!(
        report.valence_deviation_after < report.valence_deviation_before,
        "valence {} -> {}",
        report.valence_deviation_before,
        report.valence_deviation_after
    );
    assert!(report.valence_deviation_after < 0.5);
    assert!(
        report.min_angle > 20f64.to_radians(),
        "min angle {}",
        report.min_angle.to_degrees()
    );
    for t in out.triangles() {
        assert!(normal(&out, t).dot(centroid(&out, t)) > 0.0, "outward");
    }
    // Every vertex is projected onto the input.
    assert!(
        report.max_vertex_deviation < 1e-12,
        "{}",
        report.max_vertex_deviation
    );
    // Both meshes are inscribed in the sphere, so the remeshed surface is
    // within the two chord errors of the input: the sphere's own sagitta
    // over each mesh's edges.
    let sphere = |p: Point3| (radius - p.length()).abs();
    let bound = chord_error(&input, sphere) + chord_error(&out, sphere);
    let distance = hausdorff(&input, &out);
    eprintln!("sphere hausdorff {distance}, bound {bound}");
    assert!(distance <= bound, "Hausdorff {distance} over {bound}");
}

#[test]
fn torus_converges_and_keeps_its_genus() {
    let (major, minor) = (2.0, 0.6);
    let input = torus(major, minor, 40, 14);
    let target = 0.2;
    let (out, report) = remesh_checked(&input, RemeshOptions::new(target));
    assert_eq!(topology(&out).0, 0, "genus one");

    let band = band_fraction(&out, target);
    assert!(band >= 0.9, "{band} of the edges in [4/5, 4/3] L");
    assert!(report.valence_deviation_after < report.valence_deviation_before);
    assert!(
        report.min_angle > 20f64.to_radians(),
        "min angle {}",
        report.min_angle.to_degrees()
    );
    for t in out.triangles() {
        let c = centroid(&out, t);
        assert!(normal(&out, t).dot(torus_normal(major, c)) > 0.0, "outward");
    }
    assert!(
        report.max_vertex_deviation < 1e-12,
        "{}",
        report.max_vertex_deviation
    );
    // As on the sphere: within the two meshes' chord errors of the torus.
    let surface = |p: Point3| {
        let radial = Vec3::new(p.x, p.y, 0.0).normalize() * major;
        ((p - radial).length() - minor).abs()
    };
    let bound = chord_error(&input, surface) + chord_error(&out, surface);
    let distance = hausdorff(&input, &out);
    eprintln!("torus hausdorff {distance}, bound {bound}");
    assert!(distance <= bound, "Hausdorff {distance} over {bound}");
}

#[test]
fn plane_patch_keeps_its_border_and_its_plane() {
    let input = fan_square(4);
    let target = 0.1;
    // The border is protected whether or not sharp features are.
    for feature_angle in [Some(RemeshOptions::DEFAULT_FEATURE_ANGLE), None] {
        let options = RemeshOptions::new(target).with_feature_angle(feature_angle);
        check_square_patch(&input, options);
    }
}

fn check_square_patch(input: &TriMesh, options: RemeshOptions) {
    let target = options.target_edge_length;
    let (out, report) = remesh_checked(input, options);

    let band = band_fraction(&out, target);
    assert!(band >= 0.9, "{band} of the edges in [4/5, 4/3] L");
    // Border vertices aim at valence 4, interior ones at 6.
    assert!(
        report.valence_deviation_after < 0.3,
        "{}",
        report.valence_deviation_after
    );
    assert!(
        report.min_angle > 20f64.to_radians(),
        "min angle {}",
        report.min_angle.to_degrees()
    );
    for t in out.triangles() {
        assert!(normal(&out, t).z > 0.0, "upward, as the input");
    }
    for p in &out.positions {
        assert_eq!(p.z, 0.0, "in the plane");
        assert!((0.0..=1.0).contains(&p.x) && (0.0..=1.0).contains(&p.y));
    }
    // The four corners are kept bit-identical, and the border is still the
    // square's: its vertices lie on the square's sides.
    for corner in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
        assert!(
            out.positions.iter().any(|p| (p.x, p.y) == corner),
            "{corner:?} kept"
        );
    }
    let he = HalfedgeMesh::from_tri_mesh(&out).expect("manifold");
    for h in he.boundary_halfedges() {
        let p = he.position(he.target(h));
        assert!(
            p.x == 0.0 || p.x == 1.0 || p.y == 0.0 || p.y == 1.0,
            "border vertex {p:?} on a side"
        );
    }
    assert!(
        hausdorff(input, &out) <= 1e-4,
        "certified upper bound of a zero distance"
    );
}

#[test]
fn sharp_box_keeps_its_corners_and_edges() {
    let (sx, sy, sz) = (2.0, 1.0, 1.0);
    let input = box_mesh(sx, sy, sz);
    let target = 0.2;
    let (out, report) = remesh_checked(&input, RemeshOptions::new(target));

    let band = band_fraction(&out, target);
    assert!(band >= 0.9, "{band} of the edges in [4/5, 4/3] L");
    assert!(
        report.min_angle > 20f64.to_radians(),
        "min angle {}",
        report.min_angle.to_degrees()
    );
    // All eight corners kept exactly.
    for corner in &input.positions {
        assert!(out.positions.contains(corner), "corner {corner:?} kept");
    }
    // Every vertex on a face of the box; vertices on two faces lie on an
    // edge of it, and each triangle faces outwards from the box face it
    // lies on.
    let on = |v: Scalar, hi: Scalar| v == 0.0 || v == hi;
    for p in &out.positions {
        let faces = [on(p.x, sx), on(p.y, sy), on(p.z, sz)];
        assert!(faces.iter().any(|&f| f), "{p:?} on the box");
    }
    let middle = Point3::new(sx / 2.0, sy / 2.0, sz / 2.0);
    for t in out.triangles() {
        let n = normal(&out, t);
        let axis = n.abs().max_element();
        assert!(axis > 0.999 * n.length(), "triangle within one box face");
        assert!(n.dot(centroid(&out, t) - middle) > 0.0, "outward");
    }
    assert!(report.max_vertex_deviation < 1e-12);
    assert!(
        hausdorff(&input, &out) <= 1e-4,
        "certified upper bound of a zero distance"
    );
    // Feature edges: 12 box edges, each now cut into pieces near L.
    let perimeter = 4.0 * (sx + sy + sz);
    assert!(report.protected_edges as Scalar >= perimeter / (4.0 / 3.0 * target));
}

#[test]
fn without_feature_protection_the_box_still_never_folds() {
    // Only the (empty) boundary is protected: the creases are free to be
    // flipped across and smoothed over, so the shape may round off, but the
    // result stays a valid, consistently oriented closed surface.
    let input = box_mesh(1.0, 1.0, 1.0);
    let options = RemeshOptions::new(0.25).with_feature_angle(None);
    let (out, report) = remesh_checked(&input, options);
    assert_eq!(report.protected_edges, 0);
    let middle = Point3::new(0.5, 0.5, 0.5);
    for t in out.triangles() {
        assert!(
            normal(&out, t).dot(centroid(&out, t) - middle) > 0.0,
            "outward"
        );
    }
}

#[test]
fn curved_feature_lines_keep_every_input_vertex() {
    // A capped cylinder: its rims are 90-degree creases that are polygons,
    // not straight lines, so no rim vertex may go -- removing one would cut
    // the corner off the feature. The rim edges are shorter than 4/5 L and
    // stay so; everything else converges.
    let (radius, height, segments) = (1.0, 1.5, 64);
    let input = capped_cylinder(radius, height, segments);
    let target = 0.25;
    let (out, report) = remesh_checked(&input, RemeshOptions::new(target));
    for rim in &input.positions[..2 * segments as usize] {
        assert!(out.positions.contains(rim), "rim vertex {rim:?} kept");
    }
    assert_eq!(report.protected_edges, 2 * segments as usize);
    assert!(report.edges_in_band >= 0.75, "{}", report.edges_in_band);
    assert!(
        report.min_angle > 10f64.to_radians(),
        "{}",
        report.min_angle.to_degrees()
    );
    for t in out.triangles() {
        let (n, c) = (normal(&out, t), centroid(&out, t));
        let outward = if c.z == 0.0 {
            -n.z
        } else if c.z == height {
            n.z
        } else {
            n.dot(Vec3::new(c.x, c.y, 0.0))
        };
        assert!(outward > 0.0, "outward at {c:?}");
    }
    assert!(report.max_vertex_deviation < 1e-12);
}

/// Every output triangle lies in one face plane of an axis-aligned
/// polyhedron: its normal is axis-aligned and every corner shares that
/// coordinate with an input triangle facing the same way.
fn assert_on_axis_faces(input: &TriMesh, out: &TriMesh) {
    let faces: Vec<(usize, Scalar, Scalar)> = input
        .triangles()
        .map(|t| {
            let n = normal(input, t);
            let axis = (0..3)
                .max_by(|&i, &j| n[i].abs().total_cmp(&n[j].abs()))
                .unwrap_or(0);
            (axis, n[axis].signum(), input.positions[t[0] as usize][axis])
        })
        .collect();
    for t in out.triangles() {
        let n = normal(out, t);
        let axis = (0..3)
            .max_by(|&i, &j| n[i].abs().total_cmp(&n[j].abs()))
            .unwrap_or(0);
        assert!(n[axis].abs() > 0.999 * n.length(), "axis-aligned triangle");
        let corners = t.map(|i| out.positions[i as usize][axis]);
        assert!(
            faces.iter().any(|&(a, sign, level)| a == axis
                && sign == n[axis].signum()
                && corners.iter().all(|&c| c == level)),
            "triangle {t:?} lies on an input face, facing its way"
        );
    }
}

#[test]
fn concave_crease_keeps_vertices_on_their_own_side() {
    // At a reflex edge the closest point of the whole surface to a vertex
    // sliding over one face can lie on the other face; projection is
    // restricted to the vertex's own patch so it never jumps there.
    let input = l_prism();
    let target = 0.15;
    let (out, report) = remesh_checked(&input, RemeshOptions::new(target));
    for corner in &input.positions {
        assert!(out.positions.contains(corner), "corner {corner:?} kept");
    }
    assert_on_axis_faces(&input, &out);
    assert!(report.max_vertex_deviation < 1e-12);
    assert!(report.edges_in_band >= 0.9, "{}", report.edges_in_band);
    assert!(hausdorff(&input, &out) <= 1e-4);
}

#[test]
fn a_thin_slab_keeps_both_rims_of_its_narrow_faces() {
    // The narrow faces are thinner than 4/5 L: their cross edges join two
    // feature vertices and must not collapse, or one rim would be dragged
    // onto the other.
    let input = box_mesh(1.0, 1.0, 0.1);
    let (out, report) = remesh_checked(&input, RemeshOptions::new(0.2));
    for corner in &input.positions {
        assert!(out.positions.contains(corner), "corner {corner:?} kept");
    }
    assert_on_axis_faces(&input, &out);
    assert!(report.max_vertex_deviation < 1e-12);
    assert!(hausdorff(&input, &out) <= 1e-4);
}

/// Unit square of `n x n` cells with interior vertices jittered by up to a
/// quarter cell (a fixed linear congruential sequence), each cell split
/// along the diagonal that gives the fatter pair of triangles.
fn jittered_square(n: u32, seed: u64) -> TriMesh {
    let mut state = seed;
    let mut jitter = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 11) as Scalar / (1u64 << 53) as Scalar - 0.5) * 0.5
    };
    let h = 1.0 / Scalar::from(n);
    let mut positions = Vec::new();
    for i in 0..=n {
        for j in 0..=n {
            let inside = |k: u32| k > 0 && k < n;
            let (dx, dy) = (jitter(), jitter());
            let x = Scalar::from(j) * h + if inside(i) && inside(j) { dx * h } else { 0.0 };
            let y = Scalar::from(i) * h + if inside(i) && inside(j) { dy * h } else { 0.0 };
            positions.push(Point3::new(x, y, 0.0));
        }
    }
    let at = |i: u32, j: u32| i * (n + 1) + j;
    let area = |p: &[Point3], [a, b, c]: [u32; 3]| {
        let (a, b, c) = (p[a as usize], p[b as usize], p[c as usize]);
        (b - a).cross(c - a).z
    };
    let mut indices = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let (a, b, c, d) = (at(i, j), at(i, j + 1), at(i + 1, j + 1), at(i + 1, j));
            let first = [[a, b, c], [a, c, d]];
            let second = [[a, b, d], [b, c, d]];
            let worst =
                |pair: [[u32; 3]; 2]| area(&positions, pair[0]).min(area(&positions, pair[1]));
            let pair = if worst(first) >= worst(second) {
                first
            } else {
                second
            };
            for t in pair {
                assert!(area(&positions, t) > 0.0, "valid jittered input");
                indices.extend(t);
            }
        }
    }
    TriMesh::new(positions, indices)
}

#[test]
fn irregular_planar_input_never_folds() {
    // Jittered grids, coarsened and refined: every guarded edit has
    // candidates that would fold or flatten a triangle.
    for seed in 1..=4 {
        let input = jittered_square(14, seed);
        for target in [0.03, 0.12, 0.25] {
            let (out, report) = remesh_checked(&input, RemeshOptions::new(target));
            for t in out.triangles() {
                assert!(
                    normal(&out, t).z > 0.0,
                    "seed {seed}, target {target}: upward"
                );
            }
            assert!(
                report.min_angle > 10f64.to_radians(),
                "seed {seed}, target {target}: min angle {}",
                report.min_angle.to_degrees()
            );
            assert!(report.edges_in_band >= 0.8, "{}", report.edges_in_band);
        }
    }
}

fn ridge_height(x: Scalar) -> Scalar {
    0.2 * (0.5 - (x - 0.5).abs())
}

/// The unit square folded into a shallow roof along `x = 0.5`, as a
/// `4 x 4` grid: its normals turn by about 22.6 degrees across the ridge.
fn roof() -> TriMesh {
    let mut positions = Vec::new();
    for i in 0..=4 {
        for j in 0..=4 {
            let x = Scalar::from(j) / 4.0;
            positions.push(Point3::new(x, Scalar::from(i) / 4.0, ridge_height(x)));
        }
    }
    let at = |i: u32, j: u32| i * 5 + j;
    let mut indices = Vec::new();
    for i in 0..4 {
        for j in 0..4 {
            let (a, b, c, d) = (at(i, j), at(i, j + 1), at(i + 1, j + 1), at(i + 1, j));
            indices.extend([a, b, c, a, c, d]);
        }
    }
    TriMesh::new(positions, indices)
}

#[test]
fn a_ridge_meeting_the_border_keeps_both() {
    // The ridge is a feature (10 degree threshold) that ends on the
    // protected border. An edge between a ridge vertex and a border vertex
    // joins two feature lines and must never collapse: either the ridge or
    // the border would be cut short.
    let input = roof();
    let options = RemeshOptions::new(0.08).with_feature_angle(Some(10f64.to_radians()));
    let (out, report) = remesh_checked(&input, options);
    for p in &out.positions {
        assert!((p.z - ridge_height(p.x)).abs() < 1e-12, "{p:?} on the roof");
    }
    let ridge = out.positions.iter().filter(|p| p.x == 0.5).count();
    assert!(
        ridge as Scalar >= 1.0 / (4.0 / 3.0 * 0.08),
        "{ridge} ridge vertices"
    );
    for t in out.triangles() {
        assert!(normal(&out, t).z > 0.0, "upward");
    }
    assert!(report.edges_in_band >= 0.8, "{}", report.edges_in_band);
    assert!(hausdorff(&input, &out) <= 1e-4);
}

#[test]
fn coarsening_collapses_towards_the_target() {
    // A target well above the input's spacing: collapses dominate.
    let input = uv_sphere(1.0, 40, 80);
    let target = 0.3;
    let (out, report) = remesh_checked(&input, RemeshOptions::new(target));
    assert!(report.output_triangles * 4 < report.input_triangles);
    assert!(report.collapses > report.splits);
    assert!(
        band_fraction(&out, target) >= 0.9,
        "{}",
        report.edges_in_band
    );
    // The input is regular almost everywhere; the output need only be near
    // regular, not better than it.
    assert!(report.valence_deviation_after < 0.5);
    assert!(
        report.min_angle > 20f64.to_radians(),
        "{}",
        report.min_angle.to_degrees()
    );
    for t in out.triangles() {
        assert!(normal(&out, t).dot(centroid(&out, t)) > 0.0, "outward");
    }
}

#[test]
fn zero_iterations_return_the_input_connectivity() {
    let input = uv_sphere(1.0, 6, 8);
    let (out, report) = remesh(&input, RemeshOptions::new(0.3).with_iterations(0)).expect("ok");
    assert_eq!(out, input);
    assert_eq!(
        report.splits + report.collapses + report.flips + report.moves,
        0
    );
}

// ---- refusals -------------------------------------------------------------

fn triangle_pair() -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        vec![0, 1, 2, 2, 1, 3],
    )
}

#[test]
fn invalid_targets_are_refused() {
    let mesh = triangle_pair();
    for target in [0.0, -1.0, Scalar::NAN, Scalar::INFINITY] {
        let error = remesh(&mesh, RemeshOptions::new(target)).expect_err("refused");
        assert!(matches!(error, RemeshError::InvalidTarget(_)), "{error}");
    }
    for angle in [-0.1, 4.0, Scalar::NAN] {
        let options = RemeshOptions::new(0.5).with_feature_angle(Some(angle));
        let error = remesh(&mesh, options).expect_err("refused");
        assert!(
            matches!(error, RemeshError::InvalidFeatureAngle(_)),
            "{error}"
        );
    }
}

#[test]
fn non_manifold_input_is_refused_by_name() {
    // A third triangle on edge 1-2.
    let mut fin = triangle_pair();
    fin.positions.push(Point3::new(0.5, 0.5, 1.0));
    fin.indices.extend([1, 2, 4]);
    let error = remesh(&fin, RemeshOptions::new(0.5)).expect_err("refused");
    assert!(
        matches!(
            error,
            RemeshError::InvalidMesh(HalfedgeBuildError::NonManifoldEdge { .. })
        ),
        "{error}"
    );

    let mut flipped = triangle_pair();
    flipped.indices = vec![0, 1, 2, 1, 2, 3];
    let error = remesh(&flipped, RemeshOptions::new(0.5)).expect_err("refused");
    assert!(
        matches!(
            error,
            RemeshError::InvalidMesh(HalfedgeBuildError::InconsistentOrientation { .. })
        ),
        "{error}"
    );

    let mut ragged = triangle_pair();
    ragged.indices.pop();
    let error = remesh(&ragged, RemeshOptions::new(0.5)).expect_err("refused");
    assert!(matches!(error, RemeshError::InvalidMesh(_)), "{error}");
}

#[test]
fn non_finite_and_degenerate_input_is_refused() {
    let mut nan = triangle_pair();
    nan.positions[3].x = Scalar::NAN;
    let error = remesh(&nan, RemeshOptions::new(0.5)).expect_err("refused");
    assert_eq!(error, RemeshError::NonFinitePosition { vertex: 3 });

    let mut flat = triangle_pair();
    flat.positions[3] = Point3::new(-1.0, 2.0, 0.0); // on the line 1-2
    let error = remesh(&flat, RemeshOptions::new(0.5)).expect_err("refused");
    assert_eq!(error, RemeshError::DegenerateTriangle { triangle: 1 });
}

#[test]
fn an_impossible_resolution_is_refused_before_any_work() {
    let error = remesh(&uv_sphere(1000.0, 6, 8), RemeshOptions::new(1e-3)).expect_err("refused");
    assert!(
        matches!(error, RemeshError::BudgetExceeded { .. }),
        "{error}"
    );
}

#[test]
fn attribute_channels_are_reported_dropped() {
    let mut mesh = triangle_pair();
    mesh.attributes.push(axiolid_mesh::AttributeChannel::new(
        "u",
        vec![0.0, 1.0, 2.0, 3.0],
        1,
        axiolid_mesh::Blend::Linear,
    ));
    let (out, report) = remesh(&mesh, RemeshOptions::new(0.3)).expect("ok");
    assert!(out.attributes.is_empty());
    assert_eq!(report.attribute_fates.len(), 1);
    assert!(matches!(
        report.attribute_fates[0].1,
        axiolid_mesh::AttributeFate::Dropped(_)
    ));
}
