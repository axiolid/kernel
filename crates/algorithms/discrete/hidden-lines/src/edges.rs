//! Triangles, edge topology and edge classification.
//!
//! Per mesh: every triangle gets its image orientation (exactly), every edge
//! its incident triangles, and the selected edges their kind. An edge with
//! one incident triangle is a boundary; with two it is a silhouette when the
//! triangles' image orientations differ (edge-on is its own orientation),
//! otherwise a crease when the normals' cosine is strictly below the
//! threshold, decided exactly. The kind follows from the geometry; the
//! selection only decides whether an edge of that kind is drawn.

use std::cmp::Ordering;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::Point3;
use axiolid_exact::{sign_root, Arith, SignExpr};
use axiolid_guarantees::Sign;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    EdgeKind, EdgeSelection, HiddenLineError, MeshId, MeshRole,
};

use crate::decide::decide;
use crate::recipe::Orient2;
use crate::view::{Bounds, Vert};

/// One input triangle in view space.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Tri {
    /// Vertex ids into the scene's vertex list.
    pub v: [u32; 3],
    /// Image orientation: counter-clockwise, clockwise or edge-on.
    pub sigma: Ordering,
    pub bounds: Bounds,
}

/// A selected edge, oriented from its lower vertex to its higher one.
#[derive(Debug, Clone, Copy)]
pub(crate) struct DrawEdge {
    /// Position of the mesh in the request.
    pub mesh: usize,
    /// Position indices within the mesh, ascending.
    pub local: [u32; 2],
    /// Vertex ids into the scene's vertex list.
    pub ends: [u32; 2],
    pub kind: EdgeKind,
    pub bounds: Bounds,
}

/// `(b - a) × (c - a)` in model space, in arithmetic `T`.
fn normal<T: Arith>(p: &[Point3; 3]) -> [T; 3] {
    let d = |i: usize, j: usize| {
        let (to, from) = (p[i].to_array(), p[j].to_array());
        [0, 1, 2].map(|k| T::from_f64(to[k]).sub(&T::from_f64(from[k])))
    };
    let (u, v) = (d(1, 0), d(2, 0));
    let comp = |i: usize, j: usize| u[i].mul(&v[j]).sub(&u[j].mul(&v[i]));
    [comp(1, 2), comp(2, 0), comp(0, 1)]
}

/// One component of a triangle's model-space normal.
struct NormalComponent<'a> {
    points: &'a [Point3; 3],
    axis: usize,
}

impl SignExpr for NormalComponent<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        normal::<T>(self.points)[self.axis].sign()
    }
}

/// `n1·n2 - cosine·|n1|·|n2|` as the sign of `d + (-cosine)·sqrt(s)`.
struct Crease<'a> {
    first: &'a [Point3; 3],
    second: &'a [Point3; 3],
    cosine: f64,
}

impl SignExpr for Crease<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let (n1, n2) = (normal::<T>(self.first), normal::<T>(self.second));
        let dot =
            |a: &[T; 3], b: &[T; 3]| a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]));
        let d = dot(&n1, &n2);
        let s = dot(&n1, &n1).mul(&dot(&n2, &n2));
        sign_root(&d, &T::from_f64(-self.cosine), &s)
    }
}

/// What one mesh contributes.
pub(crate) struct MeshEdges<'a> {
    pub index: usize,
    pub id: MeshId,
    pub role: MeshRole,
    pub mesh: &'a TriMesh,
    /// Scene vertex id of the mesh's first position.
    pub base: u32,
}

/// Append the mesh's triangles to `tris` and its selected edges to `edges`.
pub(crate) fn collect(
    input: &MeshEdges<'_>,
    verts: &[Vert],
    selection: &EdgeSelection,
    options: &ExecutionOptions,
    escalations: &mut usize,
    tris: &mut Vec<Tri>,
    edges: &mut Vec<DrawEdge>,
) -> Result<(), HiddenLineError> {
    let refuse = |error: GeomError| HiddenLineError::mesh(input.id, error);
    let memory = || refuse(GeomError::BudgetExceeded { resource: "memory" });
    let first_tri = tris.len();
    let count = input.mesh.triangle_count();
    tris.try_reserve(count).map_err(|_| memory())?;
    let mut records: Vec<(u32, u32, u32, bool)> = Vec::new();
    records.try_reserve_exact(3 * count).map_err(|_| memory())?;

    for (t, corner) in input.mesh.indices.as_chunks::<3>().0.iter().enumerate() {
        if t % 1024 == 0 {
            options.check_cancelled()?;
        }
        let points = triangle_points(input.mesh, corner);
        for axis in 0..3 {
            let sign = decide(
                &NormalComponent {
                    points: &points,
                    axis,
                },
                escalations,
            )
            .map_err(refuse)?;
            if sign != Ordering::Equal {
                break;
            }
            if axis == 2 {
                return Err(refuse(GeomError::Degenerate(format!(
                    "triangle {t} is exactly collinear"
                ))));
            }
        }
        let v = [0, 1, 2].map(|k| input.base + corner[k]);
        let sigma = decide(
            &Orient2 {
                verts,
                a: v[0],
                b: v[1],
                c: v[2],
            },
            escalations,
        )
        .map_err(refuse)?;
        tris.push(Tri {
            v,
            sigma,
            bounds: Bounds::of(verts, &v),
        });
        for k in 0..3 {
            let (i, j) = (corner[k], corner[(k + 1) % 3]);
            records.push(if i < j {
                (i, j, t as u32, true)
            } else {
                (j, i, t as u32, false)
            });
        }
    }
    records.sort_unstable();

    let mut start = 0;
    while start < records.len() {
        let (a, b) = (records[start].0, records[start].1);
        let end = start
            + records[start..]
                .iter()
                .take_while(|r| (r.0, r.1) == (a, b))
                .count();
        let group = &records[start..end];
        start = end;
        options.check_cancelled()?;
        let kind = match group {
            [_] if input.role == MeshRole::Solid => {
                return Err(refuse(GeomError::NotManifold(format!(
                    "a solid must be closed; edge ({a}, {b}) has one triangle"
                ))));
            }
            [_] => selection.boundary.then_some(EdgeKind::Boundary),
            [first, second] if first.3 == second.3 => {
                return Err(refuse(GeomError::NotManifold(format!(
                    "edge ({a}, {b}) is traversed twice in the same direction"
                ))));
            }
            [first, second] => {
                let (t1, t2) = (first_tri + first.2 as usize, first_tri + second.2 as usize);
                let triangle = |t: usize| {
                    triangle_points(input.mesh, &input.mesh.indices[3 * (t - first_tri)..][..3])
                };
                if tris[t1].sigma != tris[t2].sigma {
                    // A silhouette is never a crease, selected or not.
                    selection.silhouette.then_some(EdgeKind::Silhouette)
                } else if let Some(cosine) = selection.crease_cosine {
                    let crease = Crease {
                        first: &triangle(t1),
                        second: &triangle(t2),
                        cosine,
                    };
                    (decide(&crease, escalations).map_err(refuse)? == Ordering::Less)
                        .then_some(EdgeKind::Crease)
                } else {
                    None
                }
            }
            _ => {
                return Err(refuse(GeomError::NotManifold(format!(
                    "edge ({a}, {b}) has {} triangles",
                    group.len()
                ))));
            }
        };
        if let Some(kind) = kind {
            let ends = [input.base + a, input.base + b];
            edges.try_reserve(1).map_err(|_| memory())?;
            edges.push(DrawEdge {
                mesh: input.index,
                local: [a, b],
                ends,
                kind,
                bounds: Bounds::of(verts, &ends),
            });
        }
    }
    Ok(())
}

fn triangle_points(mesh: &TriMesh, corner: &[u32]) -> [Point3; 3] {
    [0, 1, 2].map(|k| mesh.positions[corner[k] as usize])
}
