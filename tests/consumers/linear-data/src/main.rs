//! Linear data modelling with no query layer.
//!
//! Alignments, routes, and centrelines are often STORED and exchanged
//! long before anything intersects them. Such a consumer needs the
//! linear representation and core values only -- no intersection
//! algorithm, and certainly no mesh or provider.

use axiolid_core::Point3;
use axiolid_linear::{Polyline, Segment};

fn main() {
    let centreline = Polyline {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(3.0, 4.0, 0.0),
        ],
        closed: false,
    };
    assert_eq!(centreline.points.len(), 3);

    // Length is a property of the stored geometry, not a query against
    // other geometry, so it belongs in this closure.
    let total: f64 = centreline
        .points
        .windows(2)
        .map(|w| (w[1] - w[0]).length())
        .sum();
    assert_eq!(total, 7.0);

    let span = Segment {
        start: Point3::new(0.0, 0.0, 0.0),
        end: Point3::new(0.0, 0.0, 2.5),
    };
    assert_eq!((span.end - span.start).length(), 2.5);

    println!("linear-data ok");
}
