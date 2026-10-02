//! Tessellation of CSG primitives into closed solids.
//!
//! Every variant is analytic, so each is emitted directly. A block has no
//! curvature to sample; a sphere's is exactly known. The chord budget still
//! decides the radial segment count so a primitive and a swept face of the
//! same radius agree on what a tolerance means.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_primitive::Primitive;

/// Radial segments needed for `radius` under a chord budget.
///
/// Same sagitta rule the curve flattener uses: r(1 - cos(pi/n)) <= tol.
/// Clamped low so a degenerate tolerance cannot ask for an unbounded mesh.
fn segments(radius: Scalar, tolerance: Scalar) -> usize {
    if !(radius.is_finite() && radius > 0.0 && tolerance.is_finite() && tolerance > 0.0) {
        return 3;
    }
    let ratio = 1.0 - (tolerance / radius).min(1.0);
    let n = (core::f64::consts::PI / ratio.acos().max(1e-9)).ceil();
    (n as usize).clamp(3, 4096)
}

/// [`segments`] for one direction of a doubly curved surface, refused
/// rather than clamped when the budget needs more than the clamp allows:
/// a clamped count would leave the surface outside the bound the caller
/// asked for (#231).
fn segments_within(radius: Scalar, tolerance: Scalar) -> GeomResult<usize> {
    if !(tolerance.is_finite() && tolerance > 0.0) {
        return Err(GeomError::InvalidInput(format!(
            "chord budget must be positive and finite, got {tolerance}"
        )));
    }
    let n = segments(radius, tolerance);
    if n == 4096 && radius * (1.0 - (core::f64::consts::PI / 4096.0).cos()) > tolerance {
        return Err(GeomError::BudgetExceeded {
            resource: "primitive surface segments",
        });
    }
    Ok(n)
}

/// Tessellate one CSG primitive into a closed, outward-wound solid.
///
/// A sphere and a torus are doubly curved: every point of the exact
/// surface lies within `tolerance.linear()` of the mesh (#231). Each is the
/// revolution of an inscribed polygon (a meridian polygon, a tube circle
/// polygon) about its axis, which lies in the polygon's plane. A point of
/// the exact surface is the rotation of a point within the polygon's
/// sagitta `s_p` of the polygon, and a point of the revolved polygon is on
/// an arc of radius at most `rho_max` (the widest vertex) whose chord lies
/// in a planar trapezoid of the mesh, so within `rho_max (1 - cos(h/2))`
/// of it for a step `h` round the axis. Each direction gets half the
/// budget, so the two add to at most the whole.
///
/// Outward winding is not decoration: `axiolid-mesh-boolean-boolmesh` and the clash
/// containment test both read signed volume, and an inverted primitive
/// silently produces negative volume and wrong verdicts.
pub fn tessellate_primitive(primitive: &Primitive, tolerance: Tolerance) -> GeomResult<TriMesh> {
    let tol = tolerance.linear();
    match primitive {
        Primitive::Block { x, y, z } => block(*x, *y, *z),
        Primitive::Sphere { radius } => sphere(*radius, tol),
        Primitive::Cylinder { radius, height } => cylinder(*radius, *height, tol),
        Primitive::Cone { radius, height } => cone(*radius, *height, tol),
        Primitive::Pyramid { x, y, height } => pyramid(*x, *y, *height),
        Primitive::Torus {
            major_radius,
            minor_radius,
        } => torus(*major_radius, *minor_radius, tol),
        Primitive::Wedge {
            x,
            y,
            height,
            top_x_min,
            top_x_max,
            top_y_min,
            top_y_max,
        } => wedge(
            [*x, *y, *height],
            [*top_x_min, *top_x_max],
            [*top_y_min, *top_y_max],
        ),
        // The enum is non_exhaustive: a new family is unsupported, never
        // silently approximated by the nearest one.
        _ => Err(GeomError::Unsupported {
            backend: crate::ScalarBoolean::ID,
            operation: axiolid_contracts::Operation::Tessellation,
        }),
    }
}

/// Validate a positive finite extent, naming the offender.
fn positive(value: Scalar, what: &str) -> GeomResult<Scalar> {
    if !value.is_finite() || value <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "{what} must be positive and finite, got {value}"
        )));
    }
    Ok(value)
}

/// Axis-aligned block centred on the local origin.
fn block(x: Scalar, y: Scalar, z: Scalar) -> GeomResult<TriMesh> {
    let (hx, hy, hz) = (
        positive(x, "block x")? / 2.0,
        positive(y, "block y")? / 2.0,
        positive(z, "block z")? / 2.0,
    );
    let p = vec![
        Point3::new(-hx, -hy, -hz),
        Point3::new(hx, -hy, -hz),
        Point3::new(hx, hy, -hz),
        Point3::new(-hx, hy, -hz),
        Point3::new(-hx, -hy, hz),
        Point3::new(hx, -hy, hz),
        Point3::new(hx, hy, hz),
        Point3::new(-hx, hy, hz),
    ];
    // Outward winding, verified by the positive-volume test rather than by
    // reading the index list.
    let i = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    Ok(TriMesh::new(p, i))
}

/// Rectangular pyramid: base on z = 0, apex on +z.
fn pyramid(x: Scalar, y: Scalar, height: Scalar) -> GeomResult<TriMesh> {
    let (hx, hy) = (
        positive(x, "pyramid x")? / 2.0,
        positive(y, "pyramid y")? / 2.0,
    );
    let h = positive(height, "pyramid height")?;
    let p = vec![
        Point3::new(-hx, -hy, 0.0),
        Point3::new(hx, -hy, 0.0),
        Point3::new(hx, hy, 0.0),
        Point3::new(-hx, hy, 0.0),
        Point3::new(0.0, 0.0, h),
    ];
    let i = vec![0, 2, 1, 0, 3, 2, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4];
    Ok(TriMesh::new(p, i))
}

/// Wedge: base `[0, x] x [0, y]` at z = 0, top `[x0, x1] x [y0, y1]` at
/// z = `height`.
///
/// Every face is planar: the x sides join edges parallel to y, the y sides
/// edges parallel to x. A top collapsed to a segment or a point shares
/// vertices, so the faces that lose area are dropped and the rest meet at
/// the shared vertices rather than along zero-length edges.
fn wedge(
    [x, y, height]: [Scalar; 3],
    [x0, x1]: [Scalar; 2],
    [y0, y1]: [Scalar; 2],
) -> GeomResult<TriMesh> {
    let (x, y, h) = (
        positive(x, "wedge x")?,
        positive(y, "wedge y")?,
        positive(height, "wedge height")?,
    );
    for (value, what) in [
        (x0, "wedge top x min"),
        (x1, "wedge top x max"),
        (y0, "wedge top y min"),
        (y1, "wedge top y max"),
    ] {
        if !value.is_finite() {
            return Err(GeomError::InvalidInput(format!(
                "{what} must be finite, got {value}"
            )));
        }
    }
    if x0 > x1 || y0 > y1 {
        return Err(GeomError::InvalidInput(format!(
            "wedge top must have min <= max, got x {x0}..{x1}, y {y0}..{y1}"
        )));
    }
    let corners = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(x, 0.0, 0.0),
        Point3::new(x, y, 0.0),
        Point3::new(0.0, y, 0.0),
        Point3::new(x0, y0, h),
        Point3::new(x1, y0, h),
        Point3::new(x1, y1, h),
        Point3::new(x0, y1, h),
    ];
    // Outward faces over the corners above: bottom, top, then the sides at
    // y = 0, x = x, y = y and x = 0.
    const FACES: [[usize; 4]; 6] = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ];
    let mut p: Vec<Point3> = Vec::with_capacity(8);
    let index: Vec<u32> = corners
        .iter()
        .map(|&c| match p.iter().position(|&q| q == c) {
            Some(k) => k as u32,
            None => {
                p.push(c);
                (p.len() - 1) as u32
            }
        })
        .collect();
    let mut i = Vec::with_capacity(36);
    for face in FACES {
        let mut ring: Vec<u32> = Vec::with_capacity(4);
        for corner in face {
            let v = index[corner];
            if ring.last() != Some(&v) {
                ring.push(v);
            }
        }
        if ring.len() > 1 && ring.first() == ring.last() {
            ring.pop();
        }
        if ring.len() < 3 {
            continue;
        }
        // Each face is convex, so a fan triangulates it.
        for k in 1..ring.len() - 1 {
            i.extend([ring[0], ring[k], ring[k + 1]]);
        }
    }
    Ok(TriMesh::new(p, i))
}

/// Ring torus about +z, tube centre circle of radius `major` in z = 0.
///
/// A grid of `n` steps round the axis by `m` round the tube, each sized by
/// the same chord rule as the other curved primitives at HALF the budget
/// (#231): `n` for the outer equator, the largest circle round the axis.
/// Each grid cell is a planar trapezoid (its two edges round the axis are
/// parallel chords), split in two, so the tube polygon's sagitta and the
/// equator's are the only deviations and their sum fits the budget.
fn torus(major: Scalar, minor: Scalar, tol: Scalar) -> GeomResult<TriMesh> {
    let big = positive(major, "torus major radius")?;
    let r = positive(minor, "torus minor radius")?;
    if r >= big {
        let kind = if r == big {
            "a horn torus (minor radius equal to major)"
        } else {
            "a spindle torus (minor radius above major)"
        };
        return Err(GeomError::InvalidInput(format!(
            "torus minor radius {r} must be below major radius {big}: {kind} \
             does not bound a two-manifold solid"
        )));
    }
    let n = segments_within(big + r, 0.5 * tol)?;
    let m = segments_within(r, 0.5 * tol)?;
    let mut p = Vec::with_capacity(n * m);
    for i in 0..n {
        let theta = core::f64::consts::TAU * (i as Scalar) / (n as Scalar);
        for j in 0..m {
            let phi = core::f64::consts::TAU * (j as Scalar) / (m as Scalar);
            let rho = big + r * phi.cos();
            p.push(Point3::new(
                rho * theta.cos(),
                rho * theta.sin(),
                r * phi.sin(),
            ));
        }
    }
    let at = |i: usize, j: usize| ((i % n) * m + (j % m)) as u32;
    let mut idx = Vec::with_capacity(n * m * 6);
    for i in 0..n {
        for j in 0..m {
            // Round the axis, then round the tube: outward, since the
            // axis tangent crossed with the tube tangent points away from
            // the tube's centre.
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            idx.extend([a, b, c]);
            idx.extend([a, c, d]);
        }
    }
    Ok(TriMesh::new(p, idx))
}

/// A ring of `n` points at `radius`, height `z`.
fn ring(radius: Scalar, z: Scalar, n: usize) -> Vec<Point3> {
    (0..n)
        .map(|k| {
            let a = core::f64::consts::TAU * (k as Scalar) / (n as Scalar);
            Point3::new(radius * a.cos(), radius * a.sin(), z)
        })
        .collect()
}

/// Cylinder along +z, base on z = 0.
fn cylinder(radius: Scalar, height: Scalar, tol: Scalar) -> GeomResult<TriMesh> {
    let r = positive(radius, "cylinder radius")?;
    let h = positive(height, "cylinder height")?;
    let n = segments(r, tol);
    let mut p = ring(r, 0.0, n);
    p.extend(ring(r, h, n));
    p.push(Point3::new(0.0, 0.0, 0.0));
    p.push(Point3::new(0.0, 0.0, h));
    let (bc, tc) = (2 * n, 2 * n + 1);
    let mut i = Vec::with_capacity(n * 12);
    for k in 0..n {
        let (a, b) = (k, (k + 1) % n);
        // Side quad, then the two caps. The base fan is wound opposite to
        // the top so both face away from the enclosed volume.
        i.extend([a as u32, b as u32, (b + n) as u32]);
        i.extend([a as u32, (b + n) as u32, (a + n) as u32]);
        i.extend([bc as u32, b as u32, a as u32]);
        i.extend([tc as u32, (a + n) as u32, (b + n) as u32]);
    }
    Ok(TriMesh::new(p, i))
}

/// Cone along +z: base ring on z = 0, apex at height.
fn cone(radius: Scalar, height: Scalar, tol: Scalar) -> GeomResult<TriMesh> {
    let r = positive(radius, "cone radius")?;
    let h = positive(height, "cone height")?;
    let n = segments(r, tol);
    let mut p = ring(r, 0.0, n);
    p.push(Point3::new(0.0, 0.0, 0.0));
    p.push(Point3::new(0.0, 0.0, h));
    let (base, apex) = (n, n + 1);
    let mut i = Vec::with_capacity(n * 6);
    for k in 0..n {
        let (a, b) = (k as u32, ((k + 1) % n) as u32);
        i.extend([base as u32, b, a]);
        i.extend([a, b, apex as u32]);
    }
    Ok(TriMesh::new(p, i))
}

/// Sphere centred on the local origin, as a UV mesh.
///
/// Half the budget goes round the axis and half down the meridian (#231).
fn sphere(radius: Scalar, tol: Scalar) -> GeomResult<TriMesh> {
    let r = positive(radius, "sphere radius")?;
    let n = segments_within(r, 0.5 * tol)?;
    // Half as many stacks as segments: the polar direction spans PI, not TAU,
    // so equal counts would oversample it by 2x for the same chord error.
    // Rounded UP: an odd `n` halved down would make the polar step longer
    // than `TAU / n` and its chords deeper than the budget.
    let stacks = n.div_ceil(2).max(2);
    let mut p = Vec::with_capacity((stacks + 1) * n);
    for i in 0..=stacks {
        let v = core::f64::consts::PI * (i as Scalar) / (stacks as Scalar);
        for j in 0..n {
            let u = core::f64::consts::TAU * (j as Scalar) / (n as Scalar);
            p.push(Point3::new(
                r * v.sin() * u.cos(),
                r * v.sin() * u.sin(),
                r * v.cos(),
            ));
        }
    }
    let mut idx = Vec::with_capacity(stacks * n * 6);
    for i in 0..stacks {
        for j in 0..n {
            let jn = (j + 1) % n;
            // A pole row's n entries are the same point, so they must map
            // to ONE index or no triangle around the pole shares an edge.
            let row = |r: usize, k: usize| -> u32 {
                if r == 0 || r == stacks {
                    (r * n) as u32
                } else {
                    (r * n + k) as u32
                }
            };
            let a = row(i, j);
            let b = row(i, jn);
            let c = row(i + 1, j);
            let d = row(i + 1, jn);
            // Pole rows collapse to a point; skip the degenerate half.
            if i > 0 {
                idx.extend([a, c, b]);
            }
            if i + 1 < stacks {
                idx.extend([b, c, d]);
            }
        }
    }
    Ok(TriMesh::new(p, idx))
}
