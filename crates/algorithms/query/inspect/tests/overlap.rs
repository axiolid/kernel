//! Certified shared and difference volumes of closed meshes (#183).
//!
//! Every expected volume is a rational closed form; containment is checked
//! exactly (in dyadics), so a bracket that misses the true value by an ulp
//! fails.

use axiolid_core::Point3;
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;
use axiolid_inspect::{
    difference_volume, enclosed_volume, intersection_volume, Operand, OverlapError, VolumeInterval,
};
use axiolid_mesh::TriMesh;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn mesh(positions: Vec<Point3>, triangles: &[[u32; 3]]) -> TriMesh {
    TriMesh::new(positions, triangles.iter().flatten().copied().collect())
}

/// An axis-aligned box, outward-wound.
fn cuboid(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let positions = vec![
        p(x0, y0, z0),
        p(x1, y0, z0),
        p(x1, y1, z0),
        p(x0, y1, z0),
        p(x0, y0, z1),
        p(x1, y0, z1),
        p(x1, y1, z1),
        p(x0, y1, z1),
    ];
    mesh(
        positions,
        &[
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ],
    )
}

/// A prism over a counter-clockwise triangle, from `z0` to `z1`.
fn prism(base: [(f64, f64); 3], z0: f64, z1: f64) -> TriMesh {
    let mut positions: Vec<Point3> = base.iter().map(|&(x, y)| p(x, y, z0)).collect();
    positions.extend(base.iter().map(|&(x, y)| p(x, y, z1)));
    let mut faces = vec![[0, 2, 1], [3, 4, 5]];
    for i in 0..3u32 {
        let j = (i + 1) % 3;
        faces.push([i, j, j + 3]);
        faces.push([i, j + 3, i + 3]);
    }
    mesh(positions, &faces)
}

/// The corner tetrahedron x, y, z >= 0, x + y + z <= s.
fn corner(s: f64) -> TriMesh {
    mesh(
        vec![
            p(0.0, 0.0, 0.0),
            p(s, 0.0, 0.0),
            p(0.0, s, 0.0),
            p(0.0, 0.0, s),
        ],
        &[[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
    )
}

fn turned(m: &TriMesh) -> TriMesh {
    let mut out = m.clone();
    for t in out.indices.chunks_exact_mut(3) {
        t.swap(1, 2);
    }
    out
}

fn moved(m: &TriMesh, by: [f64; 3]) -> TriMesh {
    let mut out = m.clone();
    for q in &mut out.positions {
        *q = p(q.x + by[0], q.y + by[1], q.z + by[2]);
    }
    out
}

/// `lower <= num / den <= upper`, decided exactly.
fn brackets(v: VolumeInterval, num: f64, den: f64) {
    let (n, d) = (Dyadic::from_f64(num), Dyadic::from_f64(den));
    let below = Dyadic::from_f64(v.lower).mul(&d).sub(&n).sign();
    let above = Dyadic::from_f64(v.upper).mul(&d).sub(&n).sign();
    assert!(
        below != Some(Sign::Positive),
        "{v:?} lower above {num}/{den}"
    );
    assert!(
        above != Some(Sign::Negative),
        "{v:?} upper below {num}/{den}"
    );
    assert!(v.lower >= 0.0, "{v:?} negative");
    assert!(v.width() <= 1e-9 * (1.0 + num / den), "{v:?} too wide");
}

#[test]
fn a_box_encloses_its_volume_either_way_round() {
    let a = cuboid([0.0, 0.0, 0.0], [2.0, 3.0, 0.7]);
    brackets(enclosed_volume(&a).unwrap(), 4.2, 1.0);
    brackets(enclosed_volume(&turned(&a)).unwrap(), 4.2, 1.0);
    brackets(enclosed_volume(&corner(1.0)).unwrap(), 1.0, 6.0);
}

#[test]
fn overlapping_boxes_share_their_overlap() {
    let a = cuboid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let b = cuboid([1.0, 0.5, -1.0], [3.0, 1.5, 0.5]);
    // [1, 2] x [0.5, 1.5] x [0, 0.5].
    brackets(intersection_volume(&a, &b).unwrap(), 0.5, 1.0);
    brackets(intersection_volume(&b, &a).unwrap(), 0.5, 1.0);
    brackets(difference_volume(&a, &b).unwrap(), 7.5, 1.0);
    brackets(difference_volume(&b, &a).unwrap(), 2.5, 1.0);
}

#[test]
fn coplanar_faces_are_counted_once() {
    // Sharing the top, bottom, front and back planes.
    let a = cuboid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let b = cuboid([1.0, 0.0, 0.0], [3.0, 2.0, 2.0]);
    brackets(intersection_volume(&a, &b).unwrap(), 4.0, 1.0);
    brackets(difference_volume(&a, &b).unwrap(), 4.0, 1.0);
}

#[test]
fn a_duplicate_shares_everything_and_leaves_nothing() {
    let a = cuboid([0.0, 0.0, 0.0], [2.0, 1.0, 3.0]);
    brackets(intersection_volume(&a, &a).unwrap(), 6.0, 1.0);
    brackets(difference_volume(&a, &a).unwrap(), 0.0, 1.0);
    // Wound the other way, the same solid.
    brackets(intersection_volume(&a, &turned(&a)).unwrap(), 6.0, 1.0);
}

#[test]
fn one_inside_the_other() {
    let outer = cuboid([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
    let inner = corner(1.0);
    let inner = moved(&inner, [1.0, 1.0, 1.0]);
    brackets(intersection_volume(&outer, &inner).unwrap(), 1.0, 6.0);
    brackets(difference_volume(&inner, &outer).unwrap(), 0.0, 1.0);
    brackets(difference_volume(&outer, &inner).unwrap(), 383.0, 6.0);
}

#[test]
fn touching_bodies_share_nothing_never_less() {
    let a = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    // Face to face, edge to edge, corner to corner, and apart.
    for b in [
        cuboid([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]),
        cuboid([1.0, 1.0, 0.0], [2.0, 2.0, 1.0]),
        cuboid([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]),
        cuboid([0.0, 0.0, 1.0], [1.0, 1.0, 2.0]),
        cuboid([3.0, 0.0, 0.0], [4.0, 1.0, 1.0]),
    ] {
        let v = intersection_volume(&a, &b).unwrap();
        brackets(v, 0.0, 1.0);
        brackets(difference_volume(&a, &b).unwrap(), 1.0, 1.0);
    }
}

#[test]
fn a_slanted_face_splits_the_shadow() {
    // The corner tetrahedron against the half-size cube: the cube less
    // the small corner beyond x + y + z = 1, 1/8 - 1/48 = 5/48.
    let tetra = corner(1.0);
    let cube = cuboid([0.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
    brackets(intersection_volume(&tetra, &cube).unwrap(), 5.0, 48.0);
    brackets(difference_volume(&tetra, &cube).unwrap(), 3.0, 48.0);
    brackets(difference_volume(&cube, &tetra).unwrap(), 1.0, 48.0);
    // Shifted so no vertex coordinate is shared: the box [a, a + 1/2] with
    // a = (1/4, 1/8, 1/16). Below the plane x + y + z = 1 is the corner
    // simplex from the box's low corner, of side 1 - 7/16 = 9/16, less the
    // three pieces past the box's far faces, each a simplex of side
    // 9/16 - 1/2 = 1/16: (729 - 3) / (6 * 4096).
    let cube = cuboid([0.25, 0.125, 0.0625], [0.75, 0.625, 0.5625]);
    brackets(intersection_volume(&tetra, &cube).unwrap(), 726.0, 24576.0);
}

#[test]
fn prisms_at_an_angle() {
    // A triangular prism across a box: the box's square [1, 3]^2 less the
    // corner past x + y = 4 (area 2), one unit tall.
    let wedge = prism([(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)], 0.0, 1.0);
    let cube = cuboid([1.0, 1.0, -1.0], [3.0, 3.0, 2.0]);
    brackets(intersection_volume(&wedge, &cube).unwrap(), 2.0, 1.0);
    brackets(difference_volume(&cube, &wedge).unwrap(), 10.0, 1.0);
    // Far from the origin, the same.
    let far = [1.0e6, -2.0e6, 3.0e5];
    let v = intersection_volume(&moved(&wedge, far), &moved(&cube, far)).unwrap();
    assert!(v.contains(2.0), "{v:?}");
    assert!(v.width() < 1e-6, "{v:?}");
}

#[test]
fn refusals_name_the_operand() {
    let a = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let mut open = a.clone();
    open.indices.truncate(33);
    assert_eq!(
        intersection_volume(&a, &open),
        Err(OverlapError::NotClosed {
            operand: Operand::Second
        })
    );
    // Two boxes through each other, as one mesh.
    let b = cuboid([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);
    let mut both = a.clone();
    let shift = both.positions.len() as u32;
    both.positions.extend(b.positions.iter().copied());
    both.indices.extend(b.indices.iter().map(|i| i + shift));
    assert!(matches!(
        intersection_volume(&both, &a),
        Err(OverlapError::SelfIntersecting {
            operand: Operand::First,
            ..
        })
    ));
    let mut far = a.clone();
    far.positions[0] = p(f64::NAN, 0.0, 0.0);
    assert_eq!(
        enclosed_volume(&far),
        Err(OverlapError::NonFinite {
            operand: Operand::First
        })
    );
}

#[test]
fn random_boxes_match_their_closed_form() {
    let mut s = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        // Multiples of 1/16, so the closed form is exact.
        ((s >> 40) % 64) as f64 / 16.0
    };
    for _ in 0..40 {
        let mut corner_pair = || {
            let (a, b) = (next(), next());
            (a.min(b), a.max(b) + 0.0625)
        };
        let (ax, ay, az) = (corner_pair(), corner_pair(), corner_pair());
        let (bx, by, bz) = (corner_pair(), corner_pair(), corner_pair());
        let a = cuboid([ax.0, ay.0, az.0], [ax.1, ay.1, az.1]);
        let b = cuboid([bx.0, by.0, bz.0], [bx.1, by.1, bz.1]);
        let overlap =
            |(p0, p1): (f64, f64), (q0, q1): (f64, f64)| (p1.min(q1) - p0.max(q0)).max(0.0);
        let exact = overlap(ax, bx) * overlap(ay, by) * overlap(az, bz);
        brackets(intersection_volume(&a, &b).unwrap(), exact, 1.0);
        let whole = (ax.1 - ax.0) * (ay.1 - ay.0) * (az.1 - az.0);
        brackets(difference_volume(&a, &b).unwrap(), whole - exact, 1.0);
    }
}

/// An extrusion of a counter-clockwise outline, with its cap triangles.
fn extrusion(outline: &[(f64, f64)], caps: &[[u32; 3]], z0: f64, z1: f64) -> TriMesh {
    let n = outline.len() as u32;
    let mut positions: Vec<Point3> = outline.iter().map(|&(x, y)| p(x, y, z0)).collect();
    positions.extend(outline.iter().map(|&(x, y)| p(x, y, z1)));
    let mut faces = Vec::new();
    for &[a, b, c] in caps {
        faces.push([a, c, b]);
        faces.push([a + n, b + n, c + n]);
    }
    for i in 0..n {
        let j = (i + 1) % n;
        faces.push([i, j, j + n]);
        faces.push([i, j + n, i + n]);
    }
    mesh(positions, &faces)
}

#[test]
fn a_non_convex_solid() {
    // An L: the unit square's notch beyond (1, 1) removed from [0, 2]^2.
    let outline = [
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 1.0),
        (1.0, 1.0),
        (1.0, 2.0),
        (0.0, 2.0),
    ];
    let caps = [[3, 4, 5], [3, 5, 0], [3, 0, 1], [3, 1, 2]];
    let l = extrusion(&outline, &caps, 0.0, 1.0);
    brackets(enclosed_volume(&l).unwrap(), 3.0, 1.0);
    // A square straddling the notch's corner: three quarters of it in the L.
    let cube = cuboid([0.5, 0.5, -1.0], [1.5, 1.5, 2.0]);
    brackets(intersection_volume(&l, &cube).unwrap(), 0.75, 1.0);
    brackets(difference_volume(&cube, &l).unwrap(), 2.25, 1.0);
    // Standing on its side, so the notch is seen edge-on from above.
    let side: Vec<Point3> = l.positions.iter().map(|q| p(q.x, q.z, q.y)).collect();
    let side = mesh(
        side,
        &l.indices
            .chunks_exact(3)
            .map(|t| [t[0], t[1], t[2]])
            .collect::<Vec<_>>(),
    );
    let cube = cuboid([0.5, -1.0, 0.5], [1.5, 2.0, 1.5]);
    brackets(intersection_volume(&side, &cube).unwrap(), 0.75, 1.0);
}

#[test]
fn two_slanted_solids() {
    // The corner simplex against itself moved by 1/4 along each axis:
    // x, y, z >= 1/4 and x + y + z <= 1, a simplex of side 1/4.
    let a = corner(1.0);
    let b = moved(&a, [0.25, 0.25, 0.25]);
    brackets(intersection_volume(&a, &b).unwrap(), 1.0, 384.0);
    brackets(difference_volume(&a, &b).unwrap(), 63.0, 384.0);
}
