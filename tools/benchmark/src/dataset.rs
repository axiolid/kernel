//! Small regression datasets: inputs that once broke something.
//!
//! # Why these are code and not files
//!
//! A binary fixture in git is opaque -- a reviewer cannot see what makes it
//! interesting, and it rots silently when the format changes. These are
//! constructed in code with the reason written next to them, so a dataset
//! carries its own provenance.
//!
//! Anything genuinely large belongs in the sibling `benchmarks` repo, which
//! owns bigger corpora. The bar for living here: small enough to read, and
//! tied to a specific defect or a specific measurement.

use axiolid_core::Point3;
use axiolid_mesh::TriMesh;

use crate::workload::{self, Workload};

/// A dataset with the reason it exists.
pub struct Dataset {
    /// Stable identifier used in benchmark row names.
    pub name: &'static str,
    /// Why this input is worth keeping.
    pub rationale: &'static str,
    /// The geometry itself.
    pub mesh: TriMesh,
}

/// A sliver triangle: three near-collinear points.
///
/// The case where a floating-point orientation filter is least trustworthy,
/// and where `orient3d` escalates to exact arithmetic. Kept because the cost
/// of that escalation is the number worth watching.
pub fn sliver() -> Dataset {
    let mesh = TriMesh {
        positions: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0e-13, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ],
        indices: vec![0, 1, 2],
        ..TriMesh::default()
    };
    Dataset {
        name: "sliver",
        rationale: "near-collinear: forces predicate escalation",
        mesh,
    }
}

/// Two boxes sharing exactly one face plane.
///
/// Coplanar contact is the classic boolean degeneracy: the answer depends on
/// how ties are broken, so it is worth measuring separately from the general
/// case where operands overlap cleanly.
pub fn coplanar_contact() -> (Workload, Workload) {
    let left = workload::box_mesh(1.0, 1.0, 1.0);
    let mut right = workload::box_mesh(1.0, 1.0, 1.0);
    for position in &mut right.mesh.positions {
        position.x += 1.0;
    }
    (left, right)
}

/// Every dataset, for iterating in a benchmark or a test.
pub fn corpus() -> Vec<Dataset> {
    vec![sliver()]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_dataset_states_why_it_exists() {
        // A dataset without a rationale is a file nobody dares delete.
        for dataset in corpus() {
            assert!(
                !dataset.rationale.is_empty(),
                "{} carries no rationale",
                dataset.name
            );
        }
    }

    #[test]
    fn coplanar_boxes_touch_without_overlapping() {
        let (left, right) = coplanar_contact();
        let left_max = left
            .mesh
            .positions
            .iter()
            .map(|p| p.x)
            .fold(f64::NEG_INFINITY, f64::max);
        let right_min = right
            .mesh
            .positions
            .iter()
            .map(|p| p.x)
            .fold(f64::INFINITY, f64::min);
        // Exactly touching: the degeneracy only exists if these are equal.
        assert_eq!(left_max, right_min);
    }
}
