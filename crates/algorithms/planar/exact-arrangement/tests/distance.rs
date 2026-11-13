//! Exact distance predicates through the public API.

use std::cmp::Ordering;

use axiolid_contracts::{GeomError, Sign};
use axiolid_core::Point2;
use axiolid_exact::Dyadic;
use axiolid_exact_arrangement::{ExactPoints, Nearest, PointId};

fn flat(points: &mut ExactPoints, x: f64, y: f64) -> PointId {
    points.push_flat(Point2::new(x, y)).unwrap()
}

/// The point `(x / w, y / w)` of integer homogeneous coordinates.
fn rational(points: &mut ExactPoints, x: i64, y: i64, w: i64) -> PointId {
    let dyadic = |value: i64| Dyadic::try_from_f64(value as f64).unwrap();
    points.push([dyadic(x), dyadic(y), dyadic(w)]).unwrap()
}

#[test]
fn the_nearest_point_is_an_end_or_the_foot() {
    let mut points = ExactPoints::new();
    let a = flat(&mut points, 0.0, 0.0);
    let b = flat(&mut points, 4.0, 0.0);
    let behind = flat(&mut points, -1.0, 3.0);
    let at_start = flat(&mut points, 0.0, 3.0);
    let inside = flat(&mut points, 1.0, 3.0);
    let at_end = flat(&mut points, 4.0, 3.0);
    let beyond = flat(&mut points, 5.0, 3.0);
    let on_line = flat(&mut points, 2.0, 0.0);
    let nearest = |point| points.nearest_on_segment(point, a, b).unwrap();
    assert_eq!(nearest(behind), Nearest::Start);
    assert_eq!(nearest(at_start), Nearest::Start);
    assert_eq!(nearest(inside), Nearest::Interior);
    assert_eq!(nearest(at_end), Nearest::End);
    assert_eq!(nearest(beyond), Nearest::End);
    assert_eq!(nearest(on_line), Nearest::Interior);
    // The ends themselves are the nearest points to themselves.
    assert_eq!(nearest(a), Nearest::Start);
    assert_eq!(nearest(b), Nearest::End);
}

#[test]
fn the_foot_of_a_perpendicular_is_exact_and_not_always_binary64() {
    let mut points = ExactPoints::new();
    let a = flat(&mut points, 0.0, 0.0);
    let b = flat(&mut points, 4.0, 0.0);
    let above = flat(&mut points, 1.0, 3.0);
    let foot = points.foot(above, a, b).unwrap();
    assert_eq!(points.round(foot).unwrap(), Point2::new(1.0, 0.0));
    // The foot of (0, 5) on the line through the origin and (3, 4) is
    // (12 / 5, 16 / 5): on the line, at distance 3, and no binary64 point.
    let (c, d) = (flat(&mut points, 0.0, 0.0), flat(&mut points, 3.0, 4.0));
    let v = flat(&mut points, 0.0, 5.0);
    let foot = points.foot(v, c, d).unwrap();
    assert!(points.flat(foot).is_none());
    assert_eq!(points.orient(c, d, foot).unwrap(), Sign::Zero);
    assert_eq!(
        points.compare_distance(v, foot, 3.0).unwrap(),
        Ordering::Equal
    );
    assert_eq!(
        points.compare_distance(v, foot, 3.0_f64.next_up()).unwrap(),
        Ordering::Less
    );
    assert_eq!(
        points
            .compare_distance(v, foot, 3.0_f64.next_down())
            .unwrap(),
        Ordering::Greater
    );
    assert_eq!(
        points.nearest_on_segment(v, c, d).unwrap(),
        Nearest::Interior
    );
}

#[test]
fn a_distance_equal_to_the_length_is_equal() {
    let mut points = ExactPoints::new();
    let a = flat(&mut points, 0.0, 0.0);
    let b = flat(&mut points, 3.0, 4.0);
    assert_eq!(points.compare_distance(a, b, 5.0).unwrap(), Ordering::Equal);
    assert_eq!(
        points.compare_distance(a, b, 5.0_f64.next_up()).unwrap(),
        Ordering::Less
    );
    assert_eq!(
        points.compare_distance(a, b, 5.0_f64.next_down()).unwrap(),
        Ordering::Greater
    );
    assert_eq!(points.compare_distance(a, a, 0.0).unwrap(), Ordering::Equal);
    assert_eq!(
        points.compare_distance(a, a, 1e-300).unwrap(),
        Ordering::Less
    );
}

#[test]
fn the_distance_to_a_line_is_exact() {
    let mut points = ExactPoints::new();
    let a = flat(&mut points, 0.0, 0.0);
    let b = flat(&mut points, 3.0, 4.0);
    // (0, 5) is at distance 3 from the line, wherever the foot is.
    let v = flat(&mut points, 0.0, 5.0);
    assert_eq!(
        points.compare_line_distance(v, a, b, 3.0).unwrap(),
        Ordering::Equal
    );
    assert_eq!(
        points
            .compare_line_distance(v, a, b, 3.0_f64.next_up())
            .unwrap(),
        Ordering::Less
    );
    assert_eq!(
        points
            .compare_line_distance(v, a, b, 3.0_f64.next_down())
            .unwrap(),
        Ordering::Greater
    );
    // A point beyond the end is still measured to the line.
    let far = flat(&mut points, 30.0, 40.0);
    assert_eq!(
        points.compare_line_distance(far, a, b, 0.0).unwrap(),
        Ordering::Equal
    );
}

#[test]
fn tiny_and_huge_coordinates_are_decided() {
    let mut points = ExactPoints::new();
    let a = flat(&mut points, 1e-300, 0.0);
    let b = flat(&mut points, 3e-300, 0.0);
    let v = flat(&mut points, 2e-300, 1e-300);
    assert_eq!(
        points.nearest_on_segment(v, a, b).unwrap(),
        Nearest::Interior
    );
    assert_eq!(
        points.compare_line_distance(v, a, b, 1e-300).unwrap(),
        Ordering::Equal
    );
    let (c, d) = (flat(&mut points, 1e150, 0.0), flat(&mut points, 3e150, 0.0));
    let w = flat(&mut points, 2e150, 1e150);
    assert_eq!(
        points.compare_line_distance(w, c, d, 1e150).unwrap(),
        Ordering::Equal
    );
}

#[test]
fn bad_input_is_refused() {
    let mut points = ExactPoints::new();
    let a = flat(&mut points, 0.0, 0.0);
    let same = flat(&mut points, 0.0, 0.0);
    let b = flat(&mut points, 1.0, 0.0);
    let invalid = |result: Result<_, GeomError>| matches!(result, Err(GeomError::InvalidInput(_)));
    assert!(invalid(points.nearest_on_segment(b, a, a).map(drop)));
    assert!(invalid(points.nearest_on_segment(b, a, same).map(drop)));
    assert!(invalid(points.foot(b, a, same).map(drop)));
    assert!(invalid(
        points.compare_line_distance(b, a, same, 1.0).map(drop)
    ));
    for length in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(invalid(points.compare_distance(a, b, length).map(drop)));
        assert!(invalid(
            points.compare_line_distance(a, a, b, length).map(drop)
        ));
    }
}

/// A seeded xorshift generator.
struct Xorshift(u64);

impl Xorshift {
    fn below(&mut self, bound: u64) -> i64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % bound) as i64
    }
}

type Point = [i128; 3];

fn delta(p: Point, q: Point) -> [i128; 2] {
    [q[0] * p[2] - p[0] * q[2], q[1] * p[2] - p[1] * q[2]]
}

fn dot(u: [i128; 2], v: [i128; 2]) -> i128 {
    u[0] * v[0] + u[1] * v[1]
}

fn order(value: i128) -> Ordering {
    value.cmp(&0)
}

#[test]
fn a_thousand_seeded_cases_agree_with_integer_arithmetic() {
    let mut random = Xorshift(0xD1B5_4A32_D192_ED03);
    let (mut starts, mut interiors, mut ends, mut ties) = (0, 0, 0, 0);
    for _ in 0..1000 {
        let mut points = ExactPoints::new();
        let mut draw = |points: &mut ExactPoints| {
            let w = if random.below(2) == 0 {
                1
            } else {
                1 + random.below(5)
            };
            let (x, y) = (random.below(17) - 8, random.below(17) - 8);
            (rational(points, x, y, w), [x as i128, y as i128, w as i128])
        };
        let (a, ha) = draw(&mut points);
        let (b, hb) = draw(&mut points);
        let (v, hv) = draw(&mut points);
        let length = random.below(12) as f64;
        let g = length as i128;
        if points.same(a, b).unwrap() {
            continue;
        }
        let (u, w) = (delta(ha, hb), delta(ha, hv));
        let (along, uu) = (dot(u, w), dot(u, u));
        let expected = if along <= 0 {
            Nearest::Start
        } else if along * hb[2] >= uu * hv[2] {
            Nearest::End
        } else {
            Nearest::Interior
        };
        assert_eq!(points.nearest_on_segment(v, a, b).unwrap(), expected);
        match expected {
            Nearest::Start => starts += 1,
            Nearest::Interior => interiors += 1,
            Nearest::End => ends += 1,
        }
        let cross = u[0] * w[1] - u[1] * w[0];
        let line = order(cross * cross - g * g * ha[2] * ha[2] * hv[2] * hv[2] * uu);
        assert_eq!(points.compare_line_distance(v, a, b, length).unwrap(), line);
        let apart = delta(ha, hv);
        let pair = order(dot(apart, apart) - g * g * (ha[2] * hv[2]).pow(2));
        assert_eq!(points.compare_distance(a, v, length).unwrap(), pair);
        ties += usize::from(line == Ordering::Equal || pair == Ordering::Equal);
        if expected == Nearest::Interior {
            // The foot lies on the line and its connection to the point is
            // perpendicular, which `compare_distance` and the line distance
            // agree on.
            let foot = points.foot(v, a, b).unwrap();
            assert_eq!(points.orient(a, b, foot).unwrap(), Sign::Zero);
            assert_eq!(
                points.compare_distance(v, foot, length).unwrap(),
                points.compare_line_distance(v, a, b, length).unwrap()
            );
        }
    }
    assert!(
        starts > 50 && interiors > 50 && ends > 50,
        "{starts} {interiors} {ends}"
    );
    assert!(ties > 0);
}
