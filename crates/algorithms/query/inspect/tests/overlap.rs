//! Certified shared and difference volumes of closed meshes (#183).
//!
//! Every expected volume is a rational closed form; containment is checked
//! exactly (in dyadics), so a bracket that misses the true value by an ulp
//! fails. The one exception is the exported cabinet in its room (#313),
//! whose volume is pinned by its own certified bracket.

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

// Flush contact within rounding (#313).

/// `v` holds the dyadic `exact` and is no wider than `width`.
fn holds(v: VolumeInterval, exact: &Dyadic, width: f64) {
    holds_sixth(v, &exact.mul(&d(6.0)), width);
}

fn d(v: f64) -> Dyadic {
    Dyadic::from_f64(v)
}

/// Six times the volume a closed, outward-wound mesh encloses, exactly:
/// the sum of its triangles' determinants.
fn six_volume(m: &TriMesh) -> Dyadic {
    m.indices.chunks_exact(3).fold(d(0.0), |sum, t| {
        let [a, b, c] = [0, 1, 2].map(|k| m.positions[t[k] as usize]);
        let (a, b, c) = (
            [d(a.x), d(a.y), d(a.z)],
            [d(b.x), d(b.y), d(b.z)],
            [d(c.x), d(c.y), d(c.z)],
        );
        let cross = [
            b[1].mul(&c[2]).sub(&b[2].mul(&c[1])),
            b[2].mul(&c[0]).sub(&b[0].mul(&c[2])),
            b[0].mul(&c[1]).sub(&b[1].mul(&c[0])),
        ];
        sum.add(
            &a[0]
                .mul(&cross[0])
                .add(&a[1].mul(&cross[1]))
                .add(&a[2].mul(&cross[2])),
        )
    })
}

/// `v` holds `six / 6` and is no wider than `width`.
fn holds_sixth(v: VolumeInterval, six: &Dyadic, width: f64) {
    let below = d(v.lower).mul(&d(6.0)).sub(six).sign();
    let above = d(v.upper).mul(&d(6.0)).sub(six).sign();
    assert!(
        below != Some(Sign::Positive),
        "{v:?} lower above the volume"
    );
    assert!(
        above != Some(Sign::Negative),
        "{v:?} upper below the volume"
    );
    assert!(v.lower >= 0.0, "{v:?} negative");
    assert!(v.width() <= width, "{v:?} wider than {width:e}");
}

/// Moved by `steps` floats along each axis.
fn nudged(q: Point3, steps: [i32; 3]) -> Point3 {
    let step = |mut v: f64, n: i32| {
        for _ in 0..n.unsigned_abs() {
            v = if n > 0 { v.next_up() } else { v.next_down() };
        }
        v
    };
    p(
        step(q.x, steps[0]),
        step(q.y, steps[1]),
        step(q.z, steps[2]),
    )
}

fn box_volume(min: [f64; 3], max: [f64; 3]) -> Dyadic {
    (0..3).fold(d(1.0), |v, k| v.mul(&d(max[k]).sub(&d(min[k]))))
}

#[test]
fn exactly_coplanar_contact_is_tight() {
    // A cabinet on the floor and against the wall of a room.
    let room = cuboid([3.2, 5.0, 0.05], [8.15, 8.8, 2.25]);
    let (min, max) = ([7.55, 5.5, 0.05], [8.15, 7.1, 0.95]);
    let cabinet = cuboid(min, max);
    let v = box_volume(min, max);
    holds(intersection_volume(&cabinet, &room).unwrap(), &v, 1e-12);
    holds(difference_volume(&cabinet, &room).unwrap(), &d(0.0), 1e-12);
}

#[test]
fn contact_within_a_small_gap_is_tight() {
    let room = cuboid([3.2, 5.0, 0.05], [8.15, 8.8, 2.25]);
    for eps in [1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9] {
        // Lifted off the floor and drawn back from the wall: inside.
        let (min, max) = ([7.55 - eps, 5.5, 0.05 + eps], [8.15 - eps, 7.1, 0.95]);
        let cabinet = cuboid(min, max);
        let v = box_volume(min, max);
        holds(intersection_volume(&cabinet, &room).unwrap(), &v, 1e-12);
        holds(difference_volume(&cabinet, &room).unwrap(), &d(0.0), 1e-12);
        // Sunk into the floor: the part below it is outside.
        let (min, max) = ([7.55, 5.5, 0.05 - eps], [8.15, 7.1, 0.95]);
        let cabinet = cuboid(min, max);
        let below = box_volume(min, [max[0], max[1], 0.05]);
        let whole = box_volume(min, max);
        holds(
            intersection_volume(&cabinet, &room).unwrap(),
            &whole.sub(&below),
            1e-12,
        );
        holds(difference_volume(&cabinet, &room).unwrap(), &below, 1e-12);
    }
}

#[test]
fn near_vertical_slivers_stay_tight() {
    // The cabinet's corners at x1 are off by a few floats, as exported
    // vertical faces are, so its face there casts a needle of shadow --
    // a plane of slope near 1e30 -- which the floor's diagonal, from
    // (3.2, 5) to (9, 8.8), crosses near y = 8.24.
    let room = cuboid([3.2, 5.0, 0.05], [9.0, 8.8, 2.25]);
    for nudge in [
        [[-1, 0, 0], [0, 0, 0], [0, 0, 0], [0, 1, 0]],
        [[-2, 0, 0], [0, 0, 0], [0, 0, 0], [0, -1, 0]],
        [[1, 0, 0], [0, 0, 0], [0, 0, 0], [0, 1, 0]],
        [[0, 0, 0], [1, 0, 0], [0, 0, 0], [1, 1, 0]],
        [[0, 0, 0], [0, 1, 0], [-1, 0, 0], [0, 0, 0]],
    ] {
        let mut cabinet = cuboid([7.55, 7.6, 0.05], [8.15, 8.6, 0.95]);
        // Corners 1 and 2 at the bottom, 5 and 6 at the top.
        for (corner, steps) in [1, 2, 5, 6].into_iter().zip(nudge) {
            cabinet.positions[corner] = nudged(cabinet.positions[corner], steps);
        }
        let six = six_volume(&cabinet);
        holds_sixth(enclosed_volume(&cabinet).unwrap(), &six, 1e-12);
        holds_sixth(intersection_volume(&cabinet, &room).unwrap(), &six, 1e-12);
        holds_sixth(intersection_volume(&room, &cabinet).unwrap(), &six, 1e-12);
        holds(difference_volume(&cabinet, &room).unwrap(), &d(0.0), 1e-12);
    }
}

// A sink cabinet standing in a room, from the public buildingSMART IFC4
// building-architecture sample (IfcFurniture #176 in IfcSpace #89), as
// meshed by ifc-geometry. Its back lies within about 5e-14 of the room's
// wall, its feet within about 3e-15 of the floor, and its vertical faces
// lean by rounding residues of about 1e-15.
const CABINET_POSITIONS: &[[f64; 3]] = &[
    [8.000000000000096, 6.374999999999978, 0.9000000000000488],
    [8.000000000000098, 6.824999999999976, 0.7000000000000488],
    [8.000000000000096, 6.374999999999978, 0.7000000000000488],
    [8.000000000000098, 6.824999999999976, 0.9000000000000486],
    [8.150000000000098, 7.099999999999977, 0.9000000000000482],
    [8.150000000000096, 5.499999999999945, 4.92527347391429e-14],
    [8.150000000000098, 7.0999999999999766, 4.818758725377847e-14],
    [8.150000000000096, 5.499999999999945, 0.9000000000000488],
    [7.600000000000097, 6.374999999999979, 0.7000000000000487],
    [7.600000000000098, 6.824999999999977, 0.7000000000000486],
    [7.550000000000098, 7.099999999999978, 0.8500000000000483],
    [7.550000000000098, 7.099999999999978, 0.9000000000000483],
    [7.5700000000001, 7.099999999999978, 0.8500000000000483],
    [7.5700000000001, 7.099999999999977, 4.8363654576841063e-14],
    [7.550000000000097, 5.499999999999946, 0.900000000000048],
    [7.550000000000097, 5.499999999999946, 0.8500000000000485],
    [7.570000000000098, 5.499999999999946, 4.849199783150134e-14],
    [7.600000000000097, 6.374999999999979, 0.9000000000000486],
    [7.600000000000098, 6.824999999999977, 0.9000000000000484],
    [7.570000000000094, 5.499999999999946, 0.8500000000000485],
];
const CABINET_INDICES: &[u32] = &[
    0, 1, 2, 1, 0, 3, 4, 5, 6, 5, 4, 7, 8, 1, 9, 1, 8, 2, 4, 10, 11, 10, 4, 12, 12, 4, 13, 13, 4,
    6, 14, 10, 15, 10, 14, 11, 5, 13, 6, 13, 5, 16, 0, 8, 17, 8, 0, 2, 1, 18, 9, 18, 1, 3, 11, 18,
    4, 18, 11, 14, 18, 14, 17, 17, 14, 0, 4, 3, 7, 3, 4, 18, 7, 3, 0, 7, 0, 14, 19, 13, 16, 13, 19,
    12, 9, 17, 8, 17, 9, 18, 19, 14, 15, 14, 19, 7, 7, 19, 16, 7, 16, 5, 19, 10, 12, 10, 19, 15,
];
const ROOM_POSITIONS: &[[f64; 3]] = &[
    [7.70000000000012, 7.599999999999978, 5.0624038294699857e-14],
    [7.70000000000012, 8.29999999999998, 5.0624038294699857e-14],
    [8.150000000000073, 8.29999999999998, 5.0624038294699857e-14],
    [8.150000000000148, 8.800000000000098, 5.0624038294699857e-14],
    [3.200000000000032, 8.8000000000001, 5.0624038294699857e-14],
    [
        3.2000000000000317,
        4.999999999999939,
        5.0624038294699857e-14,
    ],
    [8.150000000000148, 4.999999999999938, 5.0624038294699857e-14],
    [8.150000000000073, 7.599999999999978, 5.0624038294699857e-14],
    [7.700000000000079, 7.599999999999906, 2.2000000000000934],
    [7.700000000000079, 8.299999999999907, 2.2000000000000934],
    [8.150000000000032, 8.299999999999907, 2.2000000000000934],
    [8.150000000000107, 8.800000000000026, 2.2000000000000934],
    [3.1999999999999913, 8.800000000000027, 2.2000000000000934],
    [3.1999999999999913, 4.999999999999867, 2.2000000000000934],
    [8.150000000000107, 4.999999999999867, 2.2000000000000934],
    [8.150000000000032, 7.599999999999906, 2.2000000000000934],
];
const ROOM_INDICES: &[u32] = &[
    9, 10, 11, 1, 3, 2, 11, 12, 13, 3, 5, 4, 13, 14, 15, 5, 7, 6, 9, 11, 13, 1, 5, 3, 13, 15, 8, 5,
    0, 7, 8, 9, 13, 0, 5, 1, 0, 1, 9, 0, 9, 8, 1, 2, 10, 1, 10, 9, 2, 3, 11, 2, 11, 10, 3, 4, 12,
    3, 12, 11, 4, 5, 13, 4, 13, 12, 5, 6, 14, 5, 14, 13, 6, 7, 15, 6, 15, 14, 7, 0, 8, 7, 8, 15,
];

fn sample(positions: &[[f64; 3]], indices: &[u32], shift: [f64; 3]) -> TriMesh {
    TriMesh::new(
        positions
            .iter()
            .map(|q| p(q[0] + shift[0], q[1] + shift[1], q[2] + shift[2]))
            .collect(),
        indices.to_vec(),
    )
}

#[test]
fn a_cabinet_flush_in_its_room_is_tight() {
    let room = sample(ROOM_POSITIONS, ROOM_INDICES, [0.0; 3]);
    // As exported, and a micrometre clear of the wall and the floor.
    for shift in [[0.0; 3], [-1e-6, 0.0, 1e-6]] {
        let cabinet = sample(CABINET_POSITIONS, CABINET_INDICES, shift);
        let whole = enclosed_volume(&cabinet).unwrap();
        let both = intersection_volume(&cabinet, &room).unwrap();
        let rest = difference_volume(&cabinet, &room).unwrap();
        assert!(both.width() < 1e-12, "{both:?}");
        assert!(both.lower > whole.lower - 1e-12, "{both:?} {whole:?}");
        assert!(both.upper <= whole.upper, "{both:?} {whole:?}");
        assert!((both.lower - 0.8008).abs() < 1e-12, "{both:?}");
        assert!(rest.upper < 1e-12, "{rest:?}");
    }
}
