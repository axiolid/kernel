//! A proximity rule checker over POINTS, with no discrete geometry.
//!
//! Clash and spacing rules are commonly stated over sensor or node
//! positions rather than solids. Such an application needs a spatial
//! index and core values -- and must NOT drag in meshes, B-rep, or a
//! provider to ask "is anything too close to anything else?".

use axiolid_core::Point3;
use axiolid_spatial::PointIndex;

const MINIMUM_SPACING: f64 = 1.0;

fn main() {
    let nodes = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(5.0, 0.0, 0.0),
        Point3::new(5.4, 0.0, 0.0),
    ];
    let index = PointIndex::build(&nodes);
    assert_eq!(index.len(), 3);

    // The rule: no two nodes closer than the minimum spacing. The pair at
    // 5.0 and 5.4 violates it; the pair at 0.0 and 5.0 does not.
    let mut violations = 0;
    for (i, node) in nodes.iter().enumerate() {
        index
            .for_each_within(*node, MINIMUM_SPACING, |hit| {
                if hit.index != i {
                    violations += 1;
                }
            })
            .expect("the radius is finite and positive");
    }

    // Counted once from each side.
    assert_eq!(violations, 2, "the 5.0/5.4 pair breaches minimum spacing");

    println!("spatial-rule-checker ok");
}
