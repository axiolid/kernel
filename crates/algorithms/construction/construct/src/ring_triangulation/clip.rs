//! Ear clipping of the bridged, weakly simple ring.
//!
//! A node `v` between `a` and `b` is an ear when `a v b` turns strictly
//! left and no other node lies in the closed triangle, except nodes at the
//! triangle's own corners (the far ends of bridges). Testing the closed
//! triangle is the point: a vertex on the diagonal `a b` blocks the ear,
//! so no triangle edge ever runs past a vertex, which is exactly the
//! T-junction earcut produced. A polygon with holes always has a
//! triangulation on its own vertices, and the leaves of its dual tree are
//! ears in this sense, so a ring with no ear is not a valid polygon; it is
//! refused rather than returned half-cut.
//!
//! Only nodes that do not turn strictly left can block an ear, so only
//! those are filed in the search grid. If the closed triangle holds any
//! node off its corners, take one farthest from `a b`: the part of the
//! triangle beyond it is free of boundary, hence inside the polygon, so
//! the node at that point whose sector faces it has both neighbours on the
//! near side and cannot turn strictly left. Cutting an ear only narrows its
//! neighbours' angles, so a node that turns left never stops doing so.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Point2;

use super::bridge::Polygon;
use super::orient;

/// Cut `polygon` into counter-clockwise triangles of vertex indices.
pub(super) fn clip_ears(points: &[Point2], polygon: &mut Polygon) -> GeomResult<Vec<[u32; 3]>> {
    let mut triangles = Vec::with_capacity(polygon.len.saturating_sub(2));
    let mut grid = Grid::new(points, polygon);
    let mut v = polygon.start;
    // Nodes visited since the last cut; a full lap without one is a stall.
    let mut since_cut = 0;
    while polygon.len > 3 {
        let node = polygon.nodes[v];
        if is_ear(points, polygon, &grid, v) {
            let (a, b) = (node.prev, node.next);
            triangles.push([
                polygon.nodes[a].vertex,
                node.vertex,
                polygon.nodes[b].vertex,
            ]);
            polygon.nodes[a].next = b;
            polygon.nodes[b].prev = a;
            polygon.len -= 1;
            grid.remove(v);
            since_cut = 0;
            // Skip past the corner just changed: cutting there next would
            // fan every triangle from one vertex, each wider than the last.
            v = polygon.nodes[b].next;
        } else {
            since_cut += 1;
            if since_cut > polygon.len {
                return Err(GeomError::Degenerate(format!(
                    "profile triangulation found no ear among {} remaining vertices",
                    polygon.len
                )));
            }
            v = node.next;
        }
    }
    // The last three nodes are the last triangle. A flat one would mean a
    // wrong cut earlier; the certificate refuses it with every other defect.
    let node = polygon.nodes[v];
    triangles.push([
        polygon.nodes[node.prev].vertex,
        node.vertex,
        polygon.nodes[node.next].vertex,
    ]);
    Ok(triangles)
}

/// Whether node `v` is an ear of `polygon`.
fn is_ear(points: &[Point2], polygon: &Polygon, grid: &Grid, v: usize) -> bool {
    let node = polygon.nodes[v];
    let at = |vertex: u32| points[vertex as usize];
    let (a, p, b) = (
        at(polygon.nodes[node.prev].vertex),
        at(node.vertex),
        at(polygon.nodes[node.next].vertex),
    );
    if orient(a, p, b) <= 0 {
        return false;
    }
    let low = Point2::new(a.x.min(p.x).min(b.x), a.y.min(p.y).min(b.y));
    let high = Point2::new(a.x.max(p.x).max(b.x), a.y.max(p.y).max(b.y));
    let blocked = grid.any_in(low, high, |other| {
        let q = at(polygon.nodes[other].vertex);
        if q.x < low.x || q.x > high.x || q.y < low.y || q.y > high.y {
            return false;
        }
        // The corners themselves, and the far ends of bridges at them.
        if q == a || q == p || q == b {
            return false;
        }
        orient(a, p, q) >= 0 && orient(p, b, q) >= 0 && orient(b, a, q) >= 0
    });
    !blocked
}

/// A uniform grid over the nodes, so an ear test reads only the nodes near
/// its triangle. Cells hold the indices of the nodes that can block an ear;
/// a node leaves its cell when it is cut.
struct Grid {
    origin: Point2,
    cell: f64,
    columns: usize,
    rows: usize,
    cells: Vec<Vec<usize>>,
    /// Per node, its cell and its position there, while it is filed.
    slots: Vec<Option<(usize, usize)>>,
}

impl Grid {
    fn new(points: &[Point2], polygon: &Polygon) -> Self {
        let mut low = Point2::new(f64::INFINITY, f64::INFINITY);
        let mut high = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
        for node in &polygon.nodes {
            let p = points[node.vertex as usize];
            low = Point2::new(low.x.min(p.x), low.y.min(p.y));
            high = Point2::new(high.x.max(p.x), high.y.max(p.y));
        }
        let span = (high.x - low.x).max(high.y - low.y);
        // About one node per cell along the longer side's share.
        let per_side = (polygon.nodes.len() as f64).sqrt().ceil().max(1.0);
        let cell = if span > 0.0 && span.is_finite() {
            span / per_side
        } else {
            1.0
        };
        let count = |extent: f64| ((extent / cell).floor() as usize).saturating_add(1);
        let (columns, rows) = (count(high.x - low.x), count(high.y - low.y));
        let mut grid = Self {
            origin: low,
            cell,
            columns,
            rows,
            cells: vec![Vec::new(); columns * rows],
            slots: vec![None; polygon.nodes.len()],
        };
        for (index, node) in polygon.nodes.iter().enumerate() {
            let p = points[node.vertex as usize];
            let (prev, next) = (
                points[polygon.nodes[node.prev].vertex as usize],
                points[polygon.nodes[node.next].vertex as usize],
            );
            if orient(prev, p, next) > 0 {
                continue;
            }
            let cell = grid.row(p.y) * columns + grid.column(p.x);
            grid.slots[index] = Some((cell, grid.cells[cell].len()));
            grid.cells[cell].push(index);
        }
        grid
    }

    /// Take a cut node out of its cell, if it was filed.
    fn remove(&mut self, node: usize) {
        let Some((cell, position)) = self.slots[node].take() else {
            return;
        };
        self.cells[cell].swap_remove(position);
        if let Some(&moved) = self.cells[cell].get(position) {
            self.slots[moved] = Some((cell, position));
        }
    }

    fn column(&self, x: f64) -> usize {
        (((x - self.origin.x) / self.cell).floor().max(0.0) as usize).min(self.columns - 1)
    }

    fn row(&self, y: f64) -> usize {
        (((y - self.origin.y) / self.cell).floor().max(0.0) as usize).min(self.rows - 1)
    }

    /// Whether `hit` holds for any node in a cell meeting the box; `hit`
    /// decides containment exactly. Every step of `column` and `row` is
    /// monotone under rounding, so a node inside the box is never filed
    /// outside its cell range.
    fn any_in(&self, low: Point2, high: Point2, mut hit: impl FnMut(usize) -> bool) -> bool {
        let (c0, c1) = (self.column(low.x), self.column(high.x));
        let (r0, r1) = (self.row(low.y), self.row(high.y));
        for r in r0..=r1 {
            for c in c0..=c1 {
                if self.cells[r * self.columns + c].iter().any(|&n| hit(n)) {
                    return true;
                }
            }
        }
        false
    }
}
