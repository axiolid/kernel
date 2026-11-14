//! Face groups resolved from a mesh (ADR 0106).
//!
//! `Single` and `PerTriangle` groups are read as given. `Planar` is derived:
//! triangles are joined across every edge that has exactly two triangles
//! whose model-space normals `n1` and `n2` are parallel and point the same
//! way, decided exactly as `n1 × n2 = 0` in all three components and
//! `n1 · n2 > 0` (the interval filter first, dyadic arithmetic when it
//! cannot decide; no tolerance). A group is an edge-connected component and
//! its id is the smallest triangle index in it, found by union-find that
//! links the larger root under the smaller, so the id depends on neither the
//! view nor the winding of the whole mesh.
//!
//! The pass polls cancellation every 1024 triangles and edge runs, and a
//! failed allocation is the `memory` budget refusal.

use std::cmp::Ordering;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::Point3;
use axiolid_exact::{Arith, SignExpr};
use axiolid_guarantees::Sign;
use axiolid_mesh_hidden_lines_contract::{FaceGroup, FaceGroups, HiddenLineError, HiddenLineMesh};

use crate::decide::decide;
use crate::edges::{normal, occlusion_flags, triangle_points};

/// The group of every triangle of `entry`.
///
/// `escalations` counts the decisions the interval filter could not settle.
pub(crate) fn resolve_groups(
    entry: &HiddenLineMesh<'_>,
    options: &ExecutionOptions,
    escalations: &mut usize,
) -> Result<Vec<FaceGroup>, HiddenLineError> {
    let count = entry.mesh.triangle_count();
    let memory = || HiddenLineError::from(GeomError::BudgetExceeded { resource: "memory" });
    let refuse = |error: GeomError| HiddenLineError::mesh(entry.id, error);
    match entry.groups {
        FaceGroups::Single => {
            let mut groups = Vec::new();
            groups.try_reserve_exact(count).map_err(|_| memory())?;
            groups.resize(count, FaceGroup(0));
            Ok(groups)
        }
        FaceGroups::PerTriangle(given) => {
            if given.len() != count {
                return Err(refuse(GeomError::InvalidInput(format!(
                    "{} face groups for {count} triangles",
                    given.len()
                ))));
            }
            let mut groups = Vec::new();
            groups.try_reserve_exact(count).map_err(|_| memory())?;
            groups.extend_from_slice(given);
            Ok(groups)
        }
        FaceGroups::Planar => planar(entry, options, escalations),
        _ => Err(refuse(GeomError::InvalidInput(
            "unknown face group mode".into(),
        ))),
    }
}

/// Which relation of two normals a [`Relation`] signs.
#[derive(Debug, Clone, Copy)]
enum Part {
    /// One component of `n1 × n2`.
    Cross(usize),
    /// `n1 · n2`.
    Dot,
}

/// A component of `n1 × n2`, or `n1 · n2`, of two triangles' normals.
struct Relation<'a> {
    first: &'a [Point3; 3],
    second: &'a [Point3; 3],
    part: Part,
}

impl SignExpr for Relation<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let (a, b) = (normal::<T>(self.first), normal::<T>(self.second));
        match self.part {
            Part::Cross(axis) => {
                let (i, j) = ((axis + 1) % 3, (axis + 2) % 3);
                a[i].mul(&b[j]).sub(&a[j].mul(&b[i])).sign()
            }
            Part::Dot => a[0]
                .mul(&b[0])
                .add(&a[1].mul(&b[1]))
                .add(&a[2].mul(&b[2]))
                .sign(),
        }
    }
}

/// Whether the normals are parallel and point the same way, exactly.
fn same_plane_and_facing(
    first: &[Point3; 3],
    second: &[Point3; 3],
    escalations: &mut usize,
) -> Result<bool, GeomError> {
    for part in [Part::Cross(0), Part::Cross(1), Part::Cross(2)] {
        let relation = Relation {
            first,
            second,
            part,
        };
        if decide(&relation, escalations)? != Ordering::Equal {
            return Ok(false);
        }
    }
    let dot = Relation {
        first,
        second,
        part: Part::Dot,
    };
    Ok(decide(&dot, escalations)? == Ordering::Greater)
}

/// The root of `t`, halving the path on the way.
fn find(parent: &mut [u32], mut t: u32) -> u32 {
    while parent[t as usize] != t {
        let grand = parent[parent[t as usize] as usize];
        parent[t as usize] = grand;
        t = grand;
    }
    t
}

fn planar(
    entry: &HiddenLineMesh<'_>,
    options: &ExecutionOptions,
    escalations: &mut usize,
) -> Result<Vec<FaceGroup>, HiddenLineError> {
    let memory = || HiddenLineError::from(GeomError::BudgetExceeded { resource: "memory" });
    let refuse = |error: GeomError| HiddenLineError::mesh(entry.id, error);
    let mesh = entry.mesh;
    let count = mesh.triangle_count();
    // Triangles of different occlusion are never joined (ADR 0111).
    let flags = occlusion_flags(entry)?;
    let mut records: Vec<(u32, u32, u32)> = Vec::new();
    records.try_reserve_exact(3 * count).map_err(|_| memory())?;
    let mut parent: Vec<u32> = Vec::new();
    parent.try_reserve_exact(count).map_err(|_| memory())?;

    for (t, corner) in mesh.indices.as_chunks::<3>().0.iter().enumerate() {
        if t % 1024 == 0 {
            options.check_cancelled()?;
        }
        parent.push(t as u32);
        for k in 0..3 {
            let (i, j) = (corner[k], corner[(k + 1) % 3]);
            records.push((i.min(j), i.max(j), t as u32));
        }
    }
    records.sort_unstable();

    let points = |t: u32| triangle_points(mesh, &mesh.indices[3 * t as usize..][..3]);
    for (n, run) in records
        .chunk_by(|a, b| (a.0, a.1) == (b.0, b.1))
        .enumerate()
    {
        if n % 1024 == 0 {
            options.check_cancelled()?;
        }
        let &[(_, _, t1), (_, _, t2)] = run else {
            continue;
        };
        if t1 == t2
            || flags.is_some_and(|flags| flags[t1 as usize] != flags[t2 as usize])
            || !same_plane_and_facing(&points(t1), &points(t2), escalations).map_err(refuse)?
        {
            continue;
        }
        let (r1, r2) = (find(&mut parent, t1), find(&mut parent, t2));
        parent[r1.max(r2) as usize] = r1.min(r2);
    }

    let mut groups = Vec::new();
    groups.try_reserve_exact(count).map_err(|_| memory())?;
    for t in 0..count as u32 {
        groups.push(FaceGroup(find(&mut parent, t)));
    }
    Ok(groups)
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_contracts::CancellationToken;
    use axiolid_core::Tolerance;
    use axiolid_mesh::TriMesh;
    use axiolid_mesh_hidden_lines_contract::{MeshId, MeshRole};

    fn square() -> TriMesh {
        TriMesh::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![0, 1, 2, 0, 2, 3],
        )
    }

    fn options() -> ExecutionOptions {
        ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"))
    }

    #[test]
    fn a_planar_square_is_one_group_named_by_its_first_triangle() {
        let mesh = square();
        let entry = HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh).with_planar_groups();
        let groups = resolve_groups(&entry, &options(), &mut 0).expect("groups");
        assert_eq!(groups, [FaceGroup(0), FaceGroup(0)]);
    }

    #[test]
    fn given_groups_are_returned_and_a_wrong_length_names_the_mesh() {
        let mesh = square();
        let given = [FaceGroup(4), FaceGroup(9)];
        let entry = HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh).with_groups(&given);
        assert_eq!(
            resolve_groups(&entry, &options(), &mut 0).expect("groups"),
            given
        );
        let short = HiddenLineMesh::new(MeshId(5), MeshRole::Sheet, &mesh).with_groups(&given[..1]);
        let error = resolve_groups(&short, &options(), &mut 0).expect_err("length");
        assert_eq!(error.mesh, Some(MeshId(5)));
        assert!(matches!(error.error, GeomError::InvalidInput(_)));
        let single = HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh);
        assert_eq!(
            resolve_groups(&single, &options(), &mut 0).expect("groups"),
            [FaceGroup(0); 2]
        );
    }

    #[test]
    fn a_cancelled_request_stops_the_planar_pass() {
        let mesh = square();
        let entry = HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh).with_planar_groups();
        let token = CancellationToken::new();
        token.cancel();
        let error = resolve_groups(&entry, &options().with_cancellation(token), &mut 0)
            .expect_err("cancelled");
        assert_eq!(error.error, GeomError::Cancelled);
    }

    #[test]
    fn planar_groups_join_only_triangles_of_one_occlusion() {
        use axiolid_mesh_hidden_lines_contract::Occlusion::{Opaque, SeeThrough};
        let mesh = square();
        for (flags, expected) in [
            ([Opaque, SeeThrough], [FaceGroup(0), FaceGroup(1)]),
            ([SeeThrough, SeeThrough], [FaceGroup(0), FaceGroup(0)]),
            ([Opaque, Opaque], [FaceGroup(0), FaceGroup(0)]),
        ] {
            let entry = HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh)
                .with_planar_groups()
                .with_occlusion(&flags);
            let groups = resolve_groups(&entry, &options(), &mut 0).expect("groups");
            assert_eq!(groups, expected, "{flags:?}");
        }
    }
}
