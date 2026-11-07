//! Exact image points for the arrangement: mesh vertices, edge splits and
//! plane crossings, all as homogeneous dyadic triples.

use axiolid_contracts::GeomResult;
use axiolid_exact::{Arith, Dyadic};
use axiolid_exact_arrangement::{ExactPoints, PointId};

use crate::param::Param;
use crate::view::Vert;

/// The point store, with each mesh vertex stored once.
pub(crate) struct Points {
    pub store: ExactPoints,
    vertices: Vec<Option<PointId>>,
}

impl Points {
    pub(crate) fn new(vertex_count: usize) -> Self {
        Self {
            store: ExactPoints::new(),
            vertices: vec![None; vertex_count],
        }
    }

    /// The homogeneous image `(U, V, h)` of a scene vertex.
    pub(crate) fn vertex(&mut self, verts: &[Vert], vertex: u32) -> GeomResult<PointId> {
        if let Some(point) = self.vertices[vertex as usize] {
            return Ok(point);
        }
        let exact = &verts[vertex as usize].exact;
        let point = self
            .store
            .push([exact[0].clone(), exact[1].clone(), exact[3].clone()])?;
        self.vertices[vertex as usize] = Some(point);
        Ok(point)
    }

    /// The image point of the edge `start -> end` at `param`.
    pub(crate) fn at_param(
        &mut self,
        verts: &[Vert],
        start: u32,
        end: u32,
        param: Param,
    ) -> GeomResult<PointId> {
        match param {
            Param::Zero => self.vertex(verts, start),
            Param::One => self.vertex(verts, end),
            Param::Root(root) => {
                let a = root.constraint.value::<Dyadic>(verts, start);
                let b = root.constraint.value::<Dyadic>(verts, end);
                self.crossing(verts, [start, end], [&a, &b])
            }
        }
    }

    /// The point of the edge `ends` where a function with the values
    /// `values` at its ends (of opposite signs) vanishes.
    pub(crate) fn crossing(
        &mut self,
        verts: &[Vert],
        ends: [u32; 2],
        values: [&Dyadic; 2],
    ) -> GeomResult<PointId> {
        let [p, q] = ends.map(|end| &verts[end as usize].exact);
        let [a, b] = values;
        let coordinate = |axis: usize| a.mul(&q[axis]).sub(&b.mul(&p[axis]));
        // TODO(#387): a hint is a coordinate that does not move along the
        // edge; with `h != 1` it is `U / h`.
        debug_assert!(ends.iter().all(|&end| verts[end as usize].is_affine()));
        let hints = [0, 1].map(|axis| fixed(&p[axis], &q[axis]));
        // The third component is the weight, `h`, of the homogeneous point.
        self.store
            .push_hinted([coordinate(0), coordinate(1), coordinate(3)], hints)
    }
}

/// The coordinate of an edge that does not move in it, as a binary64 hint.
fn fixed(from: &Dyadic, to: &Dyadic) -> Option<f64> {
    (from == to).then(|| from.to_f64())
}
