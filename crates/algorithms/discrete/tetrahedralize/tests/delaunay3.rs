//! `Delaunay3` against brute-force checkers.
//!
//! Every triangulation here is checked for: positive orientation of every
//! tetrahedron; the empty-circumsphere property against every input point;
//! symmetric adjacency with each shared face separating its two apexes;
//! the Euler characteristic of a ball; hull faces that are exactly the
//! brute-force convex hull; and coverage of the hull without overlap, by
//! counting the tetrahedra that contain sampled query points.
//!
//! Integer inputs (grids, lattice spheres, coplanar layers) are judged by an
//! `i128` oracle independent of the expansion arithmetic in
//! `axiolid-predicates`; other inputs by the certified predicates.

use std::collections::{BTreeSet, HashMap};

use axiolid_core::Point3;
use axiolid_guarantees::Sign;
use axiolid_predicates::{insphere, orient3d};
use axiolid_tetrahedralize::{Delaunay3, Delaunay3Error, Insertion, Tetrahedra};

// ---------------------------------------------------------------- oracles

fn sign_of(v: i128) -> i32 {
    v.signum() as i32
}

fn sign_value(s: Sign) -> i32 {
    match s {
        Sign::Positive => 1,
        Sign::Negative => -1,
        _ => 0,
    }
}

fn as_int(p: Point3) -> Option<[i128; 3]> {
    let ok = |v: f64| v.fract() == 0.0 && v.abs() < 1.0e6;
    (ok(p.x) && ok(p.y) && ok(p.z)).then_some([p.x as i128, p.y as i128, p.z as i128])
}

/// Sign of `(b - a) . ((c - a) x (d - a))`.
fn volume_sign(a: Point3, b: Point3, c: Point3, d: Point3) -> i32 {
    if let (Some(a), Some(b), Some(c), Some(d)) = (as_int(a), as_int(b), as_int(c), as_int(d)) {
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let w = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
        return sign_of(
            u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0]),
        );
    }
    sign_value(orient3d(a, b, d, c).sign().expect("certified"))
}

/// Positive when `e` is strictly inside the circumsphere of the positively
/// oriented `a, b, c, d`.
fn in_sphere_sign(a: Point3, b: Point3, c: Point3, d: Point3, e: Point3) -> i32 {
    if let (Some(a), Some(b), Some(c), Some(d), Some(e)) =
        (as_int(a), as_int(b), as_int(c), as_int(d), as_int(e))
    {
        let row = |p: [i128; 3]| {
            let q = [p[0] - e[0], p[1] - e[1], p[2] - e[2]];
            [q[0], q[1], q[2], q[0] * q[0] + q[1] * q[1] + q[2] * q[2]]
        };
        let m = [row(a), row(b), row(c), row(d)];
        // 4x4 determinant by cofactors along the first row.
        let det3 = |r: [usize; 3], col: [usize; 3]| {
            let x = |i: usize, j: usize| m[r[i]][col[j]];
            x(0, 0) * (x(1, 1) * x(2, 2) - x(1, 2) * x(2, 1))
                - x(0, 1) * (x(1, 0) * x(2, 2) - x(1, 2) * x(2, 0))
                + x(0, 2) * (x(1, 0) * x(2, 1) - x(1, 1) * x(2, 0))
        };
        let rows = [1, 2, 3];
        let det = m[0][0] * det3(rows, [1, 2, 3]) - m[0][1] * det3(rows, [0, 2, 3])
            + m[0][2] * det3(rows, [0, 1, 3])
            - m[0][3] * det3(rows, [0, 1, 2]);
        // With rows (p - e, |p - e|^2) the determinant is negative inside
        // for a right-handed tetrahedron.
        return -sign_of(det);
    }
    sign_value(insphere(a, b, d, c, e).sign().expect("certified"))
}

// ---------------------------------------------------------------- checker

/// Everything a Delaunay tetrahedralization of `points` must satisfy.
fn check(points: &[Point3], delaunay: &Delaunay3) -> Tetrahedra {
    let mesh = delaunay.tetrahedra();
    let tets = &mesh.tetrahedra;
    assert_eq!(delaunay.tetrahedron_count(), tets.len());
    assert_eq!(delaunay.points().len(), points.len());

    // Distinct vertices: one per distinct coordinate triple.
    let mut first: HashMap<[u64; 3], usize> = HashMap::new();
    for (i, p) in points.iter().enumerate() {
        let key = [p.x + 0.0, p.y + 0.0, p.z + 0.0].map(f64::to_bits);
        first.entry(key).or_insert(i);
        let vertex = delaunay.vertex_of(i);
        assert_eq!(
            key,
            [
                points[vertex].x + 0.0,
                points[vertex].y + 0.0,
                points[vertex].z + 0.0
            ]
            .map(f64::to_bits),
            "point {i} maps to a vertex with other coordinates"
        );
        assert_eq!(delaunay.vertex_of(vertex), vertex);
    }
    let vertices: BTreeSet<usize> = (0..points.len()).map(|i| delaunay.vertex_of(i)).collect();
    let used: BTreeSet<usize> = tets.iter().flatten().copied().collect();
    assert_eq!(
        used, vertices,
        "every distinct point is a vertex, nothing else is"
    );

    let p = |v: usize| points[v];

    // Orientation and empty circumspheres.
    for (t, tet) in tets.iter().enumerate() {
        let [a, b, c, d] = tet.map(p);
        assert_eq!(
            volume_sign(a, b, c, d),
            1,
            "tetrahedron {t} {tet:?} not positive"
        );
        for &v in &vertices {
            if tet.contains(&v) {
                continue;
            }
            assert!(
                in_sphere_sign(a, b, c, d, p(v)) <= 0,
                "point {v} strictly inside the circumsphere of {tet:?}"
            );
        }
    }

    check_adjacency(points, &mesh);
    check_euler(&mesh);
    check_hull(points, delaunay, &mesh, &vertices);
    mesh
}

fn face(tet: [usize; 4], i: usize) -> [usize; 3] {
    let mut f: Vec<usize> = (0..4).filter(|&k| k != i).map(|k| tet[k]).collect();
    f.sort_unstable();
    [f[0], f[1], f[2]]
}

fn check_adjacency(points: &[Point3], mesh: &Tetrahedra) {
    let p = |v: usize| points[v];
    let mut faces: HashMap<[usize; 3], Vec<(usize, usize)>> = HashMap::new();
    for (t, tet) in mesh.tetrahedra.iter().enumerate() {
        for i in 0..4 {
            faces.entry(face(*tet, i)).or_default().push((t, i));
        }
    }
    for (t, tet) in mesh.tetrahedra.iter().enumerate() {
        for i in 0..4 {
            let shared = &faces[&face(*tet, i)];
            match mesh.neighbors[t][i] {
                None => assert_eq!(shared.len(), 1, "hull face of {t} is shared"),
                Some(n) => {
                    assert_eq!(shared.len(), 2, "face of {t} not shared by exactly two");
                    let j = mesh.neighbors[n]
                        .iter()
                        .position(|&m| m == Some(t))
                        .expect("adjacency is symmetric");
                    assert_eq!(face(mesh.tetrahedra[n], j), face(*tet, i));
                    // The two apexes lie on opposite sides of the face.
                    let f = face(*tet, i).map(p);
                    let s1 = volume_sign(f[0], f[1], f[2], p(tet[i]));
                    let s2 = volume_sign(f[0], f[1], f[2], p(mesh.tetrahedra[n][j]));
                    assert_eq!(s1 * s2, -1, "tetrahedra {t} and {n} overlap");
                }
            }
        }
    }
}

/// A triangulated ball: V - E + F - T = 1.
fn check_euler(mesh: &Tetrahedra) {
    let mut v = BTreeSet::new();
    let mut e = BTreeSet::new();
    let mut f = BTreeSet::new();
    for tet in &mesh.tetrahedra {
        for (k, &a) in tet.iter().enumerate() {
            v.insert(a);
            for &b in &tet[k + 1..] {
                e.insert((a.min(b), a.max(b)));
            }
        }
        for i in 0..4 {
            f.insert(face(*tet, i));
        }
    }
    let chi = v.len() as i64 - e.len() as i64 + f.len() as i64 - mesh.tetrahedra.len() as i64;
    assert_eq!(chi, 1, "Euler characteristic of the tetrahedralization");
}

fn check_hull(
    points: &[Point3],
    delaunay: &Delaunay3,
    mesh: &Tetrahedra,
    vertices: &BTreeSet<usize>,
) {
    let p = |v: usize| points[v];
    let hull = delaunay.hull_triangles();

    // Same faces as the tetrahedra's unshared faces.
    let from_tets: BTreeSet<[usize; 3]> = mesh
        .tetrahedra
        .iter()
        .zip(&mesh.neighbors)
        .flat_map(|(tet, n)| (0..4).filter(|&i| n[i].is_none()).map(|i| face(*tet, i)))
        .collect();
    let sorted = |t: [usize; 3]| {
        let mut s = t;
        s.sort_unstable();
        s
    };
    let from_hull: BTreeSet<[usize; 3]> = hull.iter().map(|&t| sorted(t)).collect();
    assert_eq!(
        from_hull, from_tets,
        "hull triangles are the unshared faces"
    );
    assert_eq!(from_hull.len(), hull.len());

    // Closed, consistently oriented surface: every directed edge once.
    let mut directed = BTreeSet::new();
    for t in &hull {
        for k in 0..3 {
            assert!(directed.insert((t[k], t[(k + 1) % 3])), "edge repeated");
        }
    }
    for &(a, b) in &directed {
        assert!(directed.contains(&(b, a)), "hull is not closed at {a}-{b}");
    }

    // Every hull triangle supports the point set from outside.
    for t in &hull {
        let [a, b, c] = t.map(p);
        for &v in vertices {
            assert!(
                volume_sign(a, b, c, p(v)) <= 0,
                "point {v} beyond hull triangle {t:?}"
            );
        }
    }

    // Brute force: every supporting plane through three vertices carries at
    // least one hull triangle. Run only where n^4 is small.
    let list: Vec<usize> = vertices.iter().copied().collect();
    if list.len() <= 60 {
        for (x, &i) in list.iter().enumerate() {
            for (y, &j) in list.iter().enumerate().skip(x + 1) {
                for &k in list.iter().skip(y + 1) {
                    let (a, b, c) = (p(i), p(j), p(k));
                    if is_collinear(a, b, c) {
                        continue;
                    }
                    let sides: BTreeSet<i32> =
                        list.iter().map(|&v| volume_sign(a, b, c, p(v))).collect();
                    if sides.contains(&1) && sides.contains(&-1) {
                        continue;
                    }
                    assert!(
                        hull.iter()
                            .any(|t| t.iter().all(|&v| volume_sign(a, b, c, p(v)) == 0)),
                        "supporting plane through {i} {j} {k} has no hull triangle"
                    );
                }
            }
        }
    }
}

fn is_collinear(a: Point3, b: Point3, c: Point3) -> bool {
    let u = [b.x - a.x, b.y - a.y, b.z - a.z];
    let v = [c.x - a.x, c.y - a.y, c.z - a.z];
    // Exact for the small integer inputs the brute-force hull runs on.
    u[1] * v[2] - u[2] * v[1] == 0.0
        && u[2] * v[0] - u[0] * v[2] == 0.0
        && u[0] * v[1] - u[1] * v[0] == 0.0
}

/// Sampled coverage: a query strictly inside one tetrahedron is inside no
/// other, and a query is inside some tetrahedron iff it is inside the hull.
fn check_coverage(
    points: &[Point3],
    delaunay: &Delaunay3,
    mesh: &Tetrahedra,
    seed: u64,
    count: usize,
) {
    let p = |v: usize| points[v];
    let hull = delaunay.hull_triangles();
    let (lo, hi) = bounds(points);
    let mut rng = Rng(seed);
    for _ in 0..count {
        let q = Point3::new(
            rng.range(lo[0] - 1.0, hi[0] + 1.0),
            rng.range(lo[1] - 1.0, hi[1] + 1.0),
            rng.range(lo[2] - 1.0, hi[2] + 1.0),
        );
        let containing = mesh
            .tetrahedra
            .iter()
            .filter(|tet| {
                let [a, b, c, d] = tet.map(p);
                volume_sign(q, b, c, d) >= 0
                    && volume_sign(a, q, c, d) >= 0
                    && volume_sign(a, b, q, d) >= 0
                    && volume_sign(a, b, c, q) >= 0
            })
            .count();
        let strictly_outside = hull.iter().any(|t| {
            let [a, b, c] = t.map(p);
            volume_sign(a, b, c, q) > 0
        });
        if strictly_outside {
            assert_eq!(
                containing, 0,
                "query outside the hull lies in a tetrahedron"
            );
        } else {
            assert!(
                containing >= 1,
                "query inside the hull lies in no tetrahedron"
            );
            // Random float queries are never on a face, so exactly one.
            assert_eq!(containing, 1, "tetrahedra overlap at {q:?}");
        }
    }
}

fn bounds(points: &[Point3]) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in points {
        for (k, v) in [p.x, p.y, p.z].into_iter().enumerate() {
            lo[k] = lo[k].min(v);
            hi[k] = hi[k].max(v);
        }
    }
    (lo, hi)
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
    fn int(&mut self, n: i64) -> i64 {
        (self.next() % n as u64) as i64
    }
    fn shuffle<T>(&mut self, items: &mut [T]) {
        for i in (1..items.len()).rev() {
            let j = (self.next() % (i as u64 + 1)) as usize;
            items.swap(i, j);
        }
    }
}

fn grid(n: i64) -> Vec<Point3> {
    let mut points = Vec::new();
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                points.push(Point3::new(x as f64, y as f64, z as f64));
            }
        }
    }
    points
}

/// Tetrahedra as sets of coordinates, comparable across index orders.
fn canonical(points: &[Point3], mesh: &Tetrahedra) -> BTreeSet<Vec<[u64; 3]>> {
    mesh.tetrahedra
        .iter()
        .map(|tet| {
            let mut key: Vec<[u64; 3]> = tet
                .iter()
                .map(|&v| [points[v].x, points[v].y, points[v].z].map(f64::to_bits))
                .collect();
            key.sort_unstable();
            key
        })
        .collect()
}

/// Insert `points` one by one in a seeded random order.
fn incremental(points: &[Point3], seed: u64) -> (Vec<Point3>, Delaunay3) {
    let mut order: Vec<Point3> = points.to_vec();
    Rng(seed).shuffle(&mut order);
    let mut delaunay = Delaunay3::new();
    for &q in &order {
        delaunay.insert(q).expect("valid point");
    }
    (order, delaunay)
}

// ---------------------------------------------------------------- tests

#[test]
fn single_tetrahedron() {
    let points = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
    ];
    let delaunay = Delaunay3::from_points(&points).expect("full-dimensional");
    let mesh = check(&points, &delaunay);
    assert_eq!(mesh.tetrahedra.len(), 1);
    assert_eq!(mesh.neighbors[0], [None; 4]);
    assert_eq!(delaunay.hull_triangles().len(), 4);
}

#[test]
fn random_points_match_the_brute_force_checker() {
    for seed in 1..=12u64 {
        let mut rng = Rng(seed * 0x9E37_79B9);
        let n = 20 + 5 * seed as usize;
        let points: Vec<Point3> = (0..n)
            .map(|_| {
                Point3::new(
                    rng.range(-5.0, 5.0),
                    rng.range(-5.0, 5.0),
                    rng.range(-5.0, 5.0),
                )
            })
            .collect();
        let delaunay = Delaunay3::from_points(&points).expect("random points span 3D");
        let mesh = check(&points, &delaunay);
        check_coverage(&points, &delaunay, &mesh, seed, 300);
    }
}

#[test]
fn integer_points_match_the_brute_force_checker() {
    // Small integer coordinates collide often: duplicates, coplanar and
    // cospherical subsets all occur, judged by the i128 oracle.
    for seed in 1..=12u64 {
        let mut rng = Rng(seed * 0x5851_F42D);
        let points: Vec<Point3> = (0..60)
            .map(|_| Point3::new(rng.int(5) as f64, rng.int(5) as f64, rng.int(4) as f64))
            .collect();
        let delaunay = Delaunay3::from_points(&points).expect("spans 3D");
        let mesh = check(&points, &delaunay);
        check_coverage(&points, &delaunay, &mesh, seed, 200);
    }
}

#[test]
fn cubic_grid_is_fully_cospherical_and_still_valid() {
    // Every unit cube of a grid has eight cospherical corners, and every
    // hull face holds a planar grid of cocircular points.
    let points = grid(5);
    let delaunay = Delaunay3::from_points(&points).expect("grid spans 3D");
    let mesh = check(&points, &delaunay);
    check_coverage(&points, &delaunay, &mesh, 7, 300);
    // 6 faces, each a 4x4 grid of squares split in two.
    assert_eq!(delaunay.hull_triangles().len(), 6 * 16 * 2);
    // Total volume 64 (six times it, exact in integers).
    let six_volume: f64 = mesh
        .tetrahedra
        .iter()
        .map(|t| {
            let [a, b, c, d] = t.map(|v| points[v]);
            let u = [b.x - a.x, b.y - a.y, b.z - a.z];
            let v = [c.x - a.x, c.y - a.y, c.z - a.z];
            let w = [d.x - a.x, d.y - a.y, d.z - a.z];
            u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0])
        })
        .sum();
    assert_eq!(six_volume, 6.0 * 64.0);
}

#[test]
fn the_result_does_not_depend_on_insertion_order() {
    // The perturbation ranks points by coordinates only, so the perturbed
    // Delaunay triangulation is unique: every order gives the same one,
    // even on the maximally degenerate grid.
    let sets = [
        grid(4),
        {
            let mut rng = Rng(99);
            (0..80)
                .map(|_| Point3::new(rng.int(4) as f64, rng.int(4) as f64, rng.int(3) as f64))
                .collect()
        },
        lattice_sphere(5),
        quadruple_sphere(4097.0),
    ];
    for points in &sets {
        let reference = Delaunay3::from_points(points).expect("spans 3D");
        let expected = canonical(points, &reference.tetrahedra());
        for seed in 1..=6u64 {
            let (order, delaunay) = incremental(points, seed);
            let mesh = check(&order, &delaunay);
            assert_eq!(canonical(&order, &mesh), expected, "order {seed} differs");
        }
    }
}

/// Integer points with x^2 + y^2 + z^2 = r^2, plus the centre.
fn lattice_sphere(r: i64) -> Vec<Point3> {
    let mut points = vec![Point3::new(0.0, 0.0, 0.0)];
    for x in -r..=r {
        for y in -r..=r {
            for z in -r..=r {
                if x * x + y * y + z * z == r * r {
                    points.push(Point3::new(x as f64, y as f64, z as f64));
                }
            }
        }
    }
    points
}

#[test]
fn points_on_one_sphere() {
    // r = 5: 30 lattice points, every one cospherical with every other.
    let points = lattice_sphere(5);
    assert!(points.len() > 30);
    let delaunay = Delaunay3::from_points(&points).expect("spans 3D");
    let mesh = check(&points, &delaunay);
    check_coverage(&points, &delaunay, &mesh, 3, 300);
    // Without the centre: all points on the hull and on one sphere.
    let shell = &points[1..];
    let delaunay = Delaunay3::from_points(shell).expect("spans 3D");
    check(shell, &delaunay);
}

/// The Pythagorean quadruple `(m^2, 2, 2m; m^2 + 2)`: every sign and
/// permutation of it and the six axis points at distance `m^2 + 2` lie on
/// one sphere; the centre is added.
fn quadruple_sphere(m: f64) -> Vec<Point3> {
    let r = m * m + 2.0;
    let mut points = vec![Point3::new(0.0, 0.0, 0.0)];
    let base = [m * m, 2.0, 2.0 * m];
    for perm in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        for signs in 0..8u32 {
            let s = |k: u32| if (signs >> k) & 1 == 1 { -1.0 } else { 1.0 };
            points.push(Point3::new(
                s(0) * base[perm[0]],
                s(1) * base[perm[1]],
                s(2) * base[perm[2]],
            ));
        }
    }
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            let mut c = [0.0; 3];
            c[axis] = sign * r;
            points.push(Point3::new(c[0], c[1], c[2]));
        }
    }
    points
}

#[test]
fn cospherical_points_with_wide_coordinate_differences() {
    // With m = 4097 the coordinate differences span 25 bits, beyond the
    // in-sphere test's narrow-grid integer tier, so every tie is decided by
    // the expansion arithmetic. Four of these points on one circle make a
    // thin tetrahedron on which the old `insphere` filter certified a wrong
    // sign, and the cavity stopped being star-shaped.
    let points = quadruple_sphere(4097.0);
    let delaunay = Delaunay3::from_points(&points).expect("spans 3D");
    let mesh = check(&points, &delaunay);
    check_coverage(&points, &delaunay, &mesh, 17, 200);
}

#[test]
fn coplanar_layers_and_points_in_hull_face_planes() {
    // Three parallel planes of random lattice points: whole layers are
    // coplanar, and the outer two are hull facets with points inside them.
    let mut rng = Rng(4242);
    let mut points = Vec::new();
    for z in [0.0, 3.0, 7.0] {
        for _ in 0..25 {
            points.push(Point3::new(rng.int(9) as f64, rng.int(9) as f64, z));
        }
    }
    let delaunay = Delaunay3::from_points(&points).expect("spans 3D");
    let mesh = check(&points, &delaunay);
    check_coverage(&points, &delaunay, &mesh, 11, 300);

    // Points arriving later in the plane of an existing hull face, on and
    // off its circumcircle.
    let mut delaunay = Delaunay3::new();
    let base = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 0.0),
        Point3::new(0.0, 4.0, 0.0),
        Point3::new(0.0, 0.0, 4.0),
        Point3::new(4.0, 4.0, 0.0), // on the circumcircle of the z = 0 face
        Point3::new(6.0, -1.0, 0.0), // in the plane, outside the circle
        Point3::new(1.0, 1.0, 0.0), // in the plane, inside the hull face
        Point3::new(-1.0, 2.0, 0.0), // in the plane, outside the hull
    ];
    for q in base {
        delaunay.insert(q).expect("valid");
    }
    let mesh = check(&base, &delaunay);
    check_coverage(&base, &delaunay, &mesh, 5, 300);
}

#[test]
fn duplicates_are_merged() {
    let mut points = grid(3);
    let copies: Vec<Point3> = points.iter().step_by(4).copied().collect();
    points.extend(copies);
    points.push(Point3::new(-0.0, 0.0, -0.0)); // equal to (0, 0, 0)
    let delaunay = Delaunay3::from_points(&points).expect("spans 3D");
    check(&points, &delaunay);
    assert_eq!(delaunay.vertex_of(27), 0);
    assert_eq!(delaunay.vertex_of(points.len() - 1), 0);

    let mut delaunay = Delaunay3::new();
    let a = Point3::new(1.0, 2.0, 3.0);
    assert_eq!(delaunay.insert(a), Ok(Insertion::Vertex(0)));
    assert_eq!(delaunay.insert(a), Ok(Insertion::Duplicate(0)));
    for q in grid(2) {
        delaunay.insert(q).expect("valid");
    }
    assert_eq!(
        delaunay.insert(Point3::new(1.0, 1.0, 1.0)),
        Ok(Insertion::Duplicate(9))
    );
    assert_eq!(delaunay.insert(a), Ok(Insertion::Duplicate(0)));
}

#[test]
fn incremental_insertion_waits_for_the_first_tetrahedron() {
    let mut delaunay = Delaunay3::new();
    assert_eq!(delaunay.dimension(), None);
    let mut points = Vec::new();
    // A long coplanar prefix, including collinear points and duplicates.
    for i in 0..10 {
        points.push(Point3::new(f64::from(i), 0.0, 0.0));
    }
    points.push(Point3::new(3.0, 0.0, 0.0));
    for i in 0..10 {
        points.push(Point3::new(f64::from(i % 4), f64::from(i / 4 + 1), 0.0));
    }
    for (k, &q) in points.iter().enumerate() {
        delaunay.insert(q).expect("valid");
        let expected = if k == 0 {
            0
        } else if k < 11 {
            1
        } else {
            2
        };
        assert_eq!(delaunay.dimension(), Some(expected));
        assert_eq!(delaunay.tetrahedron_count(), 0);
    }
    let apex = Point3::new(1.5, 1.5, 2.0);
    points.push(apex);
    delaunay.insert(apex).expect("valid");
    assert_eq!(delaunay.dimension(), Some(3));
    check(&points, &delaunay);
    points.push(Point3::new(1.5, 1.5, -2.0));
    delaunay.insert(points[points.len() - 1]).expect("valid");
    let mesh = check(&points, &delaunay);
    check_coverage(&points, &delaunay, &mesh, 13, 200);
}

#[test]
fn low_dimensional_input_is_refused_by_name() {
    let refuse = |points: &[Point3]| Delaunay3::from_points(points).err();
    assert_eq!(
        refuse(&[]),
        Some(Delaunay3Error::NotFullDimensional { dimension: None })
    );
    let p = Point3::new(1.0, 1.0, 1.0);
    assert_eq!(
        refuse(&[p, p, p]),
        Some(Delaunay3Error::NotFullDimensional { dimension: Some(0) })
    );
    let line: Vec<Point3> = (0..5)
        .map(|i| Point3::new(f64::from(i), 2.0 * f64::from(i), 1.0))
        .collect();
    assert_eq!(
        refuse(&line),
        Some(Delaunay3Error::NotFullDimensional { dimension: Some(1) })
    );
    let plane: Vec<Point3> = grid(4)
        .into_iter()
        .map(|q| Point3::new(q.x, q.y, q.x + q.y))
        .collect();
    assert_eq!(
        refuse(&plane),
        Some(Delaunay3Error::NotFullDimensional { dimension: Some(2) })
    );
}

#[test]
fn unrepresentable_coordinates_are_refused_by_name() {
    let mut points = grid(2);
    points.push(Point3::new(0.5, f64::NAN, 0.5));
    assert_eq!(
        Delaunay3::from_points(&points).err(),
        Some(Delaunay3Error::NonFinite { index: 8 })
    );
    points[8] = Point3::new(0.5, 1.0e31, 0.5);
    assert_eq!(
        Delaunay3::from_points(&points).err(),
        Some(Delaunay3Error::OutOfRange { index: 8 })
    );
    points[8] = Point3::new(0.5, 1.0e-31, 0.5);
    assert_eq!(
        Delaunay3::from_points(&points).err(),
        Some(Delaunay3Error::OutOfRange { index: 8 })
    );
    let mut delaunay = Delaunay3::from_points(&points[..8]).expect("cube");
    assert_eq!(
        delaunay.insert(Point3::new(f64::INFINITY, 0.0, 0.0)),
        Err(Delaunay3Error::NonFinite { index: 8 })
    );
    assert_eq!(delaunay.points().len(), 8, "a refused point is not stored");
    assert_eq!(axiolid_tetrahedralize::MIN_COORDINATE, 2f64.powi(-100));
    assert_eq!(axiolid_tetrahedralize::MAX_COORDINATE, 2f64.powi(100));
}

#[test]
fn extreme_but_accepted_magnitudes() {
    // Coordinates at the ends of the accepted range: the exact predicates
    // must still decide, and the result must be a valid Delaunay
    // triangulation (judged by the certified predicates).
    for scale in [2f64.powi(-90), 2f64.powi(90)] {
        let points: Vec<Point3> = grid(3)
            .into_iter()
            .map(|q| Point3::new(q.x * scale, q.y * scale, q.z * scale))
            .collect();
        let delaunay = Delaunay3::from_points(&points).expect("in range");
        check(&points, &delaunay);
    }
}

#[test]
fn larger_random_sets_are_locally_delaunay() {
    // Too large for the all-pairs checker: local Delaunay on every interior
    // face (which implies global Delaunay), orientation, adjacency, Euler.
    for (seed, n) in [(1u64, 5_000usize), (2, 5_000)] {
        let mut rng = Rng(seed);
        let points: Vec<Point3> = (0..n)
            .map(|_| Point3::new(rng.unit(), rng.unit(), rng.unit()))
            .collect();
        let delaunay = Delaunay3::from_points(&points).expect("spans 3D");
        let mesh = delaunay.tetrahedra();
        check_adjacency(&points, &mesh);
        check_euler(&mesh);
        for (t, tet) in mesh.tetrahedra.iter().enumerate() {
            let [a, b, c, d] = tet.map(|v| points[v]);
            assert_eq!(volume_sign(a, b, c, d), 1);
            for i in 0..4 {
                if let Some(n) = mesh.neighbors[t][i] {
                    let j = mesh.neighbors[n]
                        .iter()
                        .position(|&m| m == Some(t))
                        .expect("symmetric");
                    let apex = points[mesh.tetrahedra[n][j]];
                    assert!(
                        in_sphere_sign(a, b, c, d, apex) <= 0,
                        "face {t}/{i} not locally Delaunay"
                    );
                }
            }
        }
        let used: BTreeSet<usize> = mesh.tetrahedra.iter().flatten().copied().collect();
        assert_eq!(used.len(), n);
    }
}
