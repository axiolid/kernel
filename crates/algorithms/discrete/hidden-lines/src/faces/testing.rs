//! Scenes for the unit tests of the face operation.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Aabb2, Frame3, Point2, Point3, Tolerance, Vec3};
use axiolid_exact_arrangement::{PairBudget, Segment};
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{HiddenLineMesh, MeshId, MeshRole, OrthographicView};

use super::points::Points;
use super::segments::PAIRS;
use super::{Faces, Work};
use crate::broad::{Broad, BroadPhase};
use crate::prepare::{prepare, Prepared, ResourceNames, SharedLimits};

pub(crate) fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"))
}

/// Looking down the world z axis: image `(x, y)`, nearer is higher.
pub(crate) fn plan() -> OrthographicView {
    OrthographicView::new(Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    })
}

pub(crate) fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// A horizontal quad `[x0,x1] x [y0,y1]` at height `z`, facing up.
pub(crate) fn sheet(x0: f64, y0: f64, x1: f64, y1: f64, z: f64) -> TriMesh {
    TriMesh::new(
        vec![p(x0, y0, z), p(x1, y0, z), p(x1, y1, z), p(x0, y1, z)],
        vec![0, 1, 2, 0, 2, 3],
    )
}

/// An outward-wound box.
pub(crate) fn solid(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let (x0, y0, z0, x1, y1, z1) = (min[0], min[1], min[2], max[0], max[1], max[2]);
    TriMesh::new(
        vec![
            p(x0, y0, z0),
            p(x1, y0, z0),
            p(x1, y1, z0),
            p(x0, y1, z0),
            p(x0, y0, z1),
            p(x1, y0, z1),
            p(x1, y1, z1),
            p(x0, y1, z1),
        ],
        vec![
            0, 3, 2, 0, 2, 1, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7,
            6, 3, 0, 4, 3, 4, 7,
        ],
    )
}

/// A prepared scene; build [`Faces`] from it with [`Setup::faces`].
pub(crate) struct Setup {
    pub prepared: Prepared,
    pub broad: Broad,
    pub meshes: Vec<(MeshId, MeshRole)>,
    pub view: OrthographicView,
    pub region: Option<Aabb2>,
    pub inputs: Vec<TriMesh>,
}

impl Setup {
    pub(crate) fn new(
        parts: &[(MeshRole, TriMesh)],
        view: OrthographicView,
        region: Option<Aabb2>,
    ) -> Self {
        let inputs: Vec<TriMesh> = parts.iter().map(|(_, m)| m.clone()).collect();
        let meshes: Vec<(MeshId, MeshRole)> = parts
            .iter()
            .enumerate()
            .map(|(i, (role, _))| (MeshId(i as u64 + 1), *role))
            .collect();
        let entries = Self::entries(&meshes, &inputs);
        let prepared = prepare(
            &entries,
            &view,
            region,
            SharedLimits {
                max_meshes: 100,
                max_source_vertices: 10_000,
                max_source_triangles: 10_000,
            },
            &ResourceNames::VISIBLE_FACES,
            &options(),
        )
        .expect("prepared");
        let broad = Broad::build(&prepared.tris, BroadPhase::Tree);
        Self {
            prepared,
            broad,
            meshes,
            view,
            region,
            inputs,
        }
    }

    fn entries<'a>(
        meshes: &[(MeshId, MeshRole)],
        inputs: &'a [TriMesh],
    ) -> Vec<HiddenLineMesh<'a>> {
        meshes
            .iter()
            .zip(inputs)
            .map(|(&(id, role), mesh)| HiddenLineMesh::new(id, role, mesh))
            .collect()
    }

    pub(crate) fn work(&self, max_pairs: usize) -> Work {
        Work {
            budget: PairBudget::new(max_pairs, PAIRS),
            escalations: 0,
        }
    }

    pub(crate) fn faces(&self, work: &mut Work) -> Faces<'_> {
        let entries = Self::entries(&self.meshes, &self.inputs);
        Faces::new(
            &self.prepared,
            &self.broad,
            &entries,
            (self.view.near, self.view.far, self.region),
            work,
            &options(),
        )
        .expect("faces")
    }

    pub(crate) fn points(&self) -> Points {
        Points::new(self.prepared.verts.len())
    }
}

/// The ends of a segment as flat coordinates.
pub(crate) fn flat(points: &Points, segment: &Segment) -> [(f64, f64); 2] {
    segment.ends.map(|end| {
        let p: Point2 = points.store.flat(end).expect("a binary64 point");
        (p.x, p.y)
    })
}
