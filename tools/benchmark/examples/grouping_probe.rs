//! How much independent work the wall workload actually exposes.
//!
//! `subtract_grouped` fuses mutually non-intersecting tools into groups and
//! processes groups SEQUENTIALLY, because each group cuts the result of the
//! previous one. Parallelism is therefore bounded by group SIZE, not group
//! count. This prints both so the ceiling is measured, not assumed.

use axiolid_benchmark::workload;
use axiolid_core::Aabb;

/// Reimplements the provider's grouping rule to observe it from outside.
fn disjoint_groups(bounds: &[Aabb]) -> Vec<Vec<usize>> {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut group_bounds: Vec<Vec<Aabb>> = Vec::new();
    'tool: for (index, bound) in bounds.iter().enumerate() {
        for (group, members) in group_bounds.iter_mut().enumerate() {
            if members.iter().all(|existing| !existing.intersects(bound)) {
                groups[group].push(index);
                members.push(*bound);
                continue 'tool;
            }
        }
        groups.push(vec![index]);
        group_bounds.push(vec![*bound]);
    }
    groups
}

fn main() {
    for openings in [4_usize, 16, 64, 256] {
        let wall = workload::wall_with_openings(openings);
        let bounds: Vec<Aabb> = wall
            .tools
            .iter()
            .map(axiolid_mesh::TriMesh::bounds)
            .collect();
        let groups = disjoint_groups(&bounds);
        let largest = groups.iter().map(Vec::len).max().unwrap_or(0);
        println!(
            "openings {openings:>4}: {} group(s), largest {largest}, sequential steps {}",
            groups.len(),
            groups.len(),
        );
    }
}
