// SPDX-License-Identifier: MPL-2.0

//! Where point location starts its walk.
//!
//! A walk from the last triangle crosses O(sqrt n) triangles for points in
//! random order, which made the walk the dominant cost once adjacency was
//! kept locally (#268). A coarse grid over the input's bounding box
//! remembers, per cell, the last vertex inserted there; the walk starts
//! from a triangle at that vertex, a few triangles from the target.
//!
//! The hint only changes where the walk starts, never what it finds when
//! the point lies strictly inside a triangle: that triangle is the only one
//! containing it. A point on an edge or a vertex lies in more than one
//! closed triangle, and which one the walk reaches decides the output's
//! triangle order. For those the build repeats the walk from the last
//! triangle, as before, so the output stays identical to the unhinted
//! build.

use axiolid_core::Point2;

/// No vertex recorded for a cell.
const EMPTY: u32 = u32::MAX;

/// The last vertex inserted in each cell of a grid over the input.
#[derive(Debug, Clone)]
pub(crate) struct CellHints {
    min: Point2,
    /// Cells per unit length.
    scale: f64,
    /// Cells along each axis.
    side: usize,
    cells: Vec<u32>,
}

impl CellHints {
    /// A grid of about one cell per four points over `[min, max]`.
    pub(crate) fn new(min: Point2, max: Point2, points: usize) -> Self {
        let side = ((points / 4) as f64).sqrt().ceil().max(1.0) as usize;
        let span = (max.x - min.x).max(max.y - min.y);
        let scale = if span > 0.0 && span.is_finite() {
            side as f64 / span
        } else {
            0.0
        };
        Self {
            min,
            scale,
            side,
            cells: vec![EMPTY; side * side],
        }
    }

    /// The cell `p` falls in, clamped to the grid.
    fn cell(&self, p: Point2) -> usize {
        // `as usize` saturates, and maps NaN to 0.
        let ix = (((p.x - self.min.x) * self.scale) as usize).min(self.side - 1);
        let iy = (((p.y - self.min.y) * self.scale) as usize).min(self.side - 1);
        iy * self.side + ix
    }

    /// The last vertex recorded in `p`'s cell.
    pub(crate) fn get(&self, p: Point2) -> Option<u32> {
        #[cfg(test)]
        if tests::UNHINTED.with(std::cell::Cell::get) {
            return None;
        }
        let v = self.cells[self.cell(p)];
        (v != EMPTY).then_some(v)
    }

    /// Record `v`, at `p`, as its cell's last vertex.
    pub(crate) fn set(&mut self, p: Point2, v: u32) {
        let cell = self.cell(p);
        self.cells[cell] = v;
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use std::cell::Cell;

    use axiolid_core::Point2;

    use crate::{triangulate, Constraint};

    thread_local! {
        /// Triangles visited by point-location walks on this thread.
        pub(crate) static STEPS: Cell<usize> = const { Cell::new(0) };
        /// Every walk starts from the last triangle, as before #268.
        pub(crate) static UNHINTED: Cell<bool> = const { Cell::new(false) };
    }

    /// Seeded uniform points in `[0, side)^2`.
    fn random(count: usize, side: f64, seed: u64) -> Vec<Point2> {
        let mut s = seed;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            side * (s >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..count).map(|_| Point2::new(next(), next())).collect()
    }

    /// The triangle list, and the walk steps it took, with or without hints.
    fn run(points: &[Point2], constraints: &[Constraint], hinted: bool) -> (Vec<u32>, usize) {
        UNHINTED.with(|c| c.set(!hinted));
        STEPS.with(|c| c.set(0));
        let tri = triangulate(points, constraints).expect("triangulates");
        UNHINTED.with(|c| c.set(false));
        (tri.triangles, STEPS.with(Cell::get))
    }

    #[test]
    fn hints_shorten_the_walk() {
        let points = random(1500, 1.0, 0x2684);
        let (_, hinted) = run(&points, &[], true);
        let (_, unhinted) = run(&points, &[], false);
        assert!(
            hinted * 3 < unhinted,
            "hinted walks took {hinted} steps, unhinted {unhinted}"
        );
        assert!(hinted < 12 * points.len(), "{hinted} steps for 1500 points");
    }

    #[test]
    fn hints_leave_the_output_unchanged() {
        // Grids put points on edges, where two triangles hold the point and
        // the walk's choice decides the triangle order; duplicates put them
        // on vertices.
        let mut cases: Vec<(Vec<Point2>, Vec<Constraint>)> = Vec::new();
        for side in [5usize, 13, 30] {
            let grid: Vec<Point2> = (0..side * side)
                .map(|k| Point2::new((k % side) as f64, (k / side) as f64))
                .collect();
            let mut reversed = grid.clone();
            reversed.reverse();
            let mut doubled = grid.clone();
            doubled.extend_from_slice(&grid[..side]);
            // In row order most points land on the hull, whose edges have
            // one finite triangle; shuffled, they land on interior edges,
            // where the two triangles beside them differ (#299).
            let mut shuffled = grid.clone();
            let order = random(grid.len(), 1.0, 0x2686 + side as u64);
            for k in (1..shuffled.len()).rev() {
                let j = (order[k].x * (k + 1) as f64) as usize;
                shuffled.swap(k, j);
            }
            cases.extend([grid, reversed, doubled, shuffled].map(|g| (g, Vec::new())));
        }
        let mut points = random(800, 10.0, 0x2685);
        cases.push((points.clone(), Vec::new()));
        let a = points.len() as u32;
        points.push(Point2::new(-0.5, 2.0));
        points.push(Point2::new(10.5, 7.0));
        cases.push((points, vec![Constraint::new(a, a + 1)]));
        for (k, (points, constraints)) in cases.iter().enumerate() {
            assert_eq!(
                run(points, constraints, true).0,
                run(points, constraints, false).0,
                "case {k}"
            );
        }
    }
}
