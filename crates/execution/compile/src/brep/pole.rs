//! Faces that reach a pole of their surface (#463, ADR 0112).
//!
//! A pole is a line of parameters, `v = v_pole`, that one point carries for
//! every `u`: a sphere's poles, a cone's apex, a revolution's profile on its
//! axis, a ruled surface's apex side ([`super::support::PoleKind`]). The
//! B-rep has no edge there (ADR 0059): the pole is a vertex, and a face's
//! loop leaves the pole line with a gap in parameters that is no gap on the
//! surface. A sphere or a revolution reaching its axis also needs a seam; a
//! trim that circles a pole with no seam edge stays as it was.
//!
//! **Rows.** Where two uses of a ring meet at a pole (both an end at it,
//! read from the use's pcurve by [`Support::is_pole`] or from the inverter
//! for a use without one), the ring gets a row of samples `(u_k, v_pole)`
//! between them, each carrying the pole vertex's mesh index, so the row's
//! segments are ring edges of no length on the surface (I1). The ring is cut
//! into chains at its pole junctions, each localised by continuity and,
//! where a use carries a pcurve, shifted by whole periods so its first
//! pcurve use keeps its parameters (I2). With one pole the row runs from the
//! chain's end back to its start; with several, every chain anchored, the
//! rows run as the pcurves place them, so between two poles the pcurves
//! name the side of the meridians the face covers (the far meridian of a
//! second half at `u + pi`, not back at `0`); two poles without that, on a
//! periodic `u`, with the chains a whole number of periods apart, are one
//! full turn; anything else is refused (P2). A row narrower than the
//! angular tolerance is no row, one wider than a period is invalid (P3), and
//! a hole reaching a pole is refused (P4) (I3). A row copies the `u` of the
//! ring's one run of samples at one `v` spanning it (a rim), bit for bit,
//! so the face stays a lattice; otherwise it is uniform, no step above an
//! eighth of a period or the surface's bend along `u` allows (I4).
//!
//! **Fans.** A ring that is a rectangle with a pole row is meshed as a grid
//! ([`super::recognise_grid`], whose lattice may then be uneven), and each
//! cell of a pole row is one triangle into the pole, a [`Fan`], bounded over
//! the cell's whole parameter quad by the bilinear interpolation bound
//! ([`FaceBound::fan`]): second order, where the triangle bound is first
//! order at a pole. A fan is split only at the middle of its edge opposite
//! the pole, into two fans, while the term along the row is what keeps it
//! above the budget. A row whose cells' term across it and corner offsets
//! already miss the budget is demoted, all of it, to the nominal triangles
//! (I7). A ring that is no grid, or a demoted row, is meshed by the
//! clipper and bounded by the triangle bound, whose first-order term near a
//! pole needs fine triangles there and can leave the face above the budget
//! (I8).
//!
//! **Supports without a pole kind.** A face whose trim closes across a
//! collapsed edge of any other surface (a B-spline's degenerate edge) is
//! refused by name (P1).

use super::support::Support;
use super::{
    localized_parameter, refinement_midpoint, surface_triangle, CurvedMeshState, FaceBound,
    FaceContext, GridPatch, LocalEdgeKey, SurfaceTriangle, SurfaceVertex, MAX_CURVED_EDGE_SEGMENTS,
};
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Scalar};
use axiolid_mesh::TriMesh;
use std::collections::{HashMap, HashSet};

/// Whether a trim jumps from `end` to `start` across a collapsed edge of
/// `surface`: both are points where `u` does not move the point but not the
/// same parameter point. Read on supports without a [`support::PoleKind`]
/// only, whose poles the tessellator does not close.
fn across_pole(support: Support<'_>, end: Point2, start: Point2) -> bool {
    // A surface relation without a pole kind has no poles of its own.
    let Some(surface) = support.atomic() else {
        return false;
    };
    let pole = |p: Point2| match axiolid_reference::surface::partials(surface, p.x, p.y) {
        Ok((du, dv)) => du.length() <= 1e-9 * (1.0 + dv.length()),
        Err(_) => false,
    };
    end != start && pole(end) && pole(start)
}

/// P1: the refusal for a face [`across_pole`] names. `TODO(#499)`.
fn pole_refusal() -> GeomError {
    GeomError::UnsupportedInput {
        backend: crate::BACKEND_ID,
        operation: axiolid_contracts::Operation::Tessellation,
        input: "a curved B-rep face whose trim closes across a collapsed edge of a surface that \
                is not a sphere, a cone, a surface of revolution or a ruled surface with an apex",
    }
}

/// One edge use's samples on a ring of a curved face, both ends included.
pub(super) struct UseSamples {
    pub(super) params: Vec<Point2>,
    pub(super) shared: Vec<u32>,
    /// Whether the use carries a pcurve, whose parameters are as authored.
    pub(super) pcurve: bool,
    /// Whether its first and its last sample are a pole of the surface.
    pub(super) at_pole: [bool; 2],
}

impl UseSamples {
    fn first(&self) -> Point2 {
        self.params[0]
    }

    fn last(&self) -> Point2 {
        self.params[self.params.len() - 1]
    }

    /// Move the use by the whole periods that bring `from`, one of its own
    /// parameters, nearest `to`.
    fn shift(&mut self, from: Point2, to: Point2, periods: (Option<Scalar>, Option<Scalar>)) {
        let shift = Point2::new(
            localized_parameter(to.x, from.x, periods.0),
            localized_parameter(to.y, from.y, periods.1),
        ) - from;
        for p in &mut self.params {
            *p += shift;
        }
    }
}

/// The uses between two pole junctions of a ring, in ring order.
struct Chain {
    uses: Vec<UseSamples>,
    /// Whether a use of it carries a pcurve: its position is then authored.
    anchored: bool,
}

impl Chain {
    fn first(&self) -> Point2 {
        self.uses[0].first()
    }

    fn last(&self) -> Point2 {
        self.uses[self.uses.len() - 1].last()
    }

    /// Localise the uses by continuity, from the first pcurve use, which
    /// keeps its parameters, or from the first use.
    fn localise(&mut self, periods: (Option<Scalar>, Option<Scalar>)) {
        let anchor = self.uses.iter().position(|u| u.pcurve).unwrap_or(0);
        for i in anchor + 1..self.uses.len() {
            let (from, to) = (self.uses[i].first(), self.uses[i - 1].last());
            self.uses[i].shift(from, to, periods);
        }
        for i in (0..anchor).rev() {
            let (from, to) = (self.uses[i].last(), self.uses[i + 1].first());
            self.uses[i].shift(from, to, periods);
        }
    }

    fn shift_by(&mut self, du: Scalar) {
        for use_ in &mut self.uses {
            for p in &mut use_.params {
                p.x += du;
            }
        }
    }
}

/// P2: two or more poles without the pcurves that place the meridians
/// between them. `TODO(#500)`.
fn sides_refusal() -> GeomError {
    GeomError::UnsupportedInput {
        backend: crate::BACKEND_ID,
        operation: axiolid_contracts::Operation::Tessellation,
        input: "a curved B-rep face between two poles whose edge uses carry no pcurves: the side \
                of its meridians it covers is not determined",
    }
}

/// One ring of a curved face from its uses' samples (#463, ADR 0112): each
/// use's samples but its last, and where two uses meet at a pole of the
/// surface, a row of samples along the pole line, all welded to the pole's
/// vertex. Gaps between uses that do not meet at a pole are bounded as
/// before.
///
/// The ring is cut into chains at its pole junctions, each localised by
/// continuity and anchored by its first pcurve use. With one pole the row
/// runs from the chain's end back to its start; with several, every chain
/// anchored, the rows run as the pcurves place them; two poles without
/// that, on a periodic `u`, close only as one full turn; anything else is
/// refused (P2). A row wider than one period is invalid (P3), and a hole
/// reaching a pole is refused (P4).
pub(super) fn pole_ring(
    state: &mut CurvedMeshState<'_>,
    ctx: &FaceContext<'_>,
    surface: Support<'_>,
    bound: &mut FaceBound<'_>,
    uses: Vec<UseSamples>,
    outer: bool,
) -> GeomResult<(Vec<Point2>, Vec<u32>)> {
    let count = uses.len();
    let junctions: Vec<bool> = (0..count)
        .map(|j| uses[j].at_pole[1] && uses[(j + 1) % count].at_pole[0])
        .collect();
    if surface.pole_kind().is_none() {
        for j in 0..count {
            if across_pole(surface, uses[j].last(), uses[(j + 1) % count].first()) {
                return Err(pole_refusal());
            }
        }
    }
    let Some(first) = junctions.iter().position(|&pole| pole) else {
        let mut uv = Vec::new();
        let mut shared = Vec::new();
        for (j, use_) in uses.iter().enumerate() {
            bound.gap(use_.last(), uses[(j + 1) % count].first(), surface);
            let n = use_.params.len() - 1;
            uv.extend_from_slice(&use_.params[..n]);
            shared.extend_from_slice(&use_.shared[..n]);
        }
        return Ok((uv, shared));
    };
    if !outer {
        // TODO(#501): a hole whose trim reaches a pole.
        return Err(GeomError::UnsupportedInput {
            backend: crate::BACKEND_ID,
            operation: axiolid_contracts::Operation::Tessellation,
            input: "a hole of a curved B-rep face whose trim reaches a pole",
        });
    }
    // Chains in ring order, starting after the first pole junction.
    let mut chains: Vec<Chain> = Vec::new();
    let mut uses: Vec<Option<UseSamples>> = uses.into_iter().map(Some).collect();
    let mut current: Vec<UseSamples> = Vec::new();
    for k in 1..=count {
        let j = (first + k) % count;
        let use_ = uses[j].take().ok_or_else(|| {
            GeomError::Degenerate("a curved ring's use was taken twice".to_owned())
        })?;
        current.push(use_);
        if junctions[j] {
            let anchored = current.iter().any(|u| u.pcurve);
            chains.push(Chain {
                uses: core::mem::take(&mut current),
                anchored,
            });
        }
    }
    let periods = surface.periods();
    for chain in &mut chains {
        chain.localise(periods);
        for pair in chain.uses.windows(2) {
            bound.gap(pair[0].last(), pair[1].first(), surface);
        }
    }
    let k = chains.len();
    let angular = ctx.tolerance.angular();
    let row = |chains: &[Chain], c: usize| chains[(c + 1) % k].first().x - chains[c].last().x;
    if k >= 2 && !chains.iter().all(|chain| chain.anchored) {
        // Two poles, the chains on one meridian: one full turn.
        let period = periods.0.ok_or_else(sides_refusal)?;
        let whole = |du: Scalar| (du - (du / period).round() * period).abs() <= angular;
        if k != 2 || !(whole(row(&chains, 0)) && whole(row(&chains, 1))) {
            return Err(sides_refusal());
        }
        let turn = chains[0].last().x + period - chains[1].first().x;
        chains[1].shift_by((turn / period).round() * period);
    }
    let (mut uv, mut shared) = (Vec::new(), Vec::new());
    let runs = level_runs(&chains);
    let samples: Vec<Point2> = chains
        .iter()
        .flat_map(|chain| chain.uses.iter())
        .flat_map(|use_| use_.params.iter().copied())
        .collect();
    let ring_box = parameter_box(&samples);
    for c in 0..k {
        for use_ in &chains[c].uses {
            let n = use_.params.len() - 1;
            uv.extend_from_slice(&use_.params[..n]);
            shared.extend_from_slice(&use_.shared[..n]);
        }
        let (end, start) = (chains[c].last(), chains[(c + 1) % k].first());
        let last = &chains[c].uses[chains[c].uses.len() - 1];
        let pole = last.shared[last.shared.len() - 1];
        if chains[(c + 1) % k].uses[0].shared[0] != pole {
            return Err(GeomError::Degenerate(
                "a curved ring's uses meet at a pole on two different vertices".to_owned(),
            ));
        }
        let du = start.x - end.x;
        if du.abs() <= angular {
            bound.gap(end, start, surface);
            continue;
        }
        if periods.0.is_some_and(|period| du.abs() > period + angular) {
            return Err(GeomError::InvalidInput(
                "a pole row of a curved B-rep face spans more than one period".to_owned(),
            ));
        }
        uv.push(end);
        shared.push(pole);
        let row = row_samples(ctx, surface, bound, &runs, (end, start), ring_box)?;
        state.reserve_face_vertices(row.len())?;
        shared.extend(core::iter::repeat_n(pole, row.len()));
        uv.extend(row);
    }
    Ok((uv, shared))
}

/// The maximal runs of a ring's samples at one `v`, strictly monotone in
/// `u`, as `(u values, v)`: a rim the pole rows may copy (#463).
fn level_runs(chains: &[Chain]) -> Vec<(Vec<Scalar>, Scalar)> {
    // A pole's samples are no rim: each breaks a run (`None`).
    let samples: Vec<Option<Point2>> = chains
        .iter()
        .flat_map(|chain| chain.uses.iter())
        .flat_map(|use_| {
            let last = use_.params.len() - 1;
            use_.params.iter().enumerate().map(move |(k, &p)| {
                let pole = (k == 0 && use_.at_pole[0]) || (k == last && use_.at_pole[1]);
                (!pole).then_some(p)
            })
        })
        .collect();
    // Relative to the ring's spans, as `recognise_grid` reads its border.
    let all = chains
        .iter()
        .flat_map(|chain| chain.uses.iter())
        .flat_map(|use_| use_.params.iter());
    let (lo, hi) = all
        .clone()
        .fold((Scalar::INFINITY, Scalar::NEG_INFINITY), |(lo, hi), p| {
            (lo.min(p.y), hi.max(p.y))
        });
    let (left, right) = all.fold((Scalar::INFINITY, Scalar::NEG_INFINITY), |(lo, hi), p| {
        (lo.min(p.x), hi.max(p.x))
    });
    let same = |a: Scalar, b: Scalar| (a - b).abs() <= 1e-9 * (hi - lo);
    // Two uses' shared end, read once.
    let repeated =
        |a: Point2, b: Point2| same(a.y, b.y) && (a.x - b.x).abs() <= 1e-9 * (right - left);
    let mut runs: Vec<(Vec<Scalar>, Scalar)> = Vec::new();
    let mut current: Vec<Point2> = Vec::new();
    let mut close = |current: &mut Vec<Point2>| {
        let monotone = current.windows(2).all(|w| w[1].x > w[0].x)
            || current.windows(2).all(|w| w[1].x < w[0].x);
        if current.len() >= 2 && monotone {
            runs.push((current.iter().map(|p| p.x).collect(), current[0].y));
        }
        current.clear();
    };
    for p in samples {
        let Some(p) = p else {
            close(&mut current);
            continue;
        };
        if let Some(&previous) = current.last() {
            if repeated(previous, p) {
                continue;
            }
            if !same(previous.y, p.y) {
                close(&mut current);
            }
        }
        current.push(p);
    }
    close(&mut current);
    runs
}

/// The samples strictly inside a pole row from `end` to `start` (I4): a
/// run's `u` values bit for bit where exactly one run spans the row's
/// ends, else uniform, as fine as a period over eight and as the surface's
/// bend along `u` over the ring asks.
fn row_samples(
    ctx: &FaceContext<'_>,
    surface: Support<'_>,
    bound: &FaceBound<'_>,
    runs: &[(Vec<Scalar>, Scalar)],
    (end, start): (Point2, Point2),
    ring_box: Option<(Point2, Point2)>,
) -> GeomResult<Vec<Point2>> {
    let (lo, hi) = (end.x.min(start.x), end.x.max(start.x));
    let slack = 1e-9 * (hi - lo);
    let spanning: Vec<&Vec<Scalar>> = runs
        .iter()
        .filter(|(us, _)| {
            let (a, b) = (us[0].min(us[us.len() - 1]), us[0].max(us[us.len() - 1]));
            a <= lo + slack && b >= hi - slack
        })
        .map(|(us, _)| us)
        .collect();
    let along = |u: Scalar| {
        let t = (u - end.x) / (start.x - end.x);
        Point2::new(u, end.y + (start.y - end.y) * t)
    };
    if let [run] = spanning.as_slice() {
        let mut us: Vec<Scalar> = run
            .iter()
            .copied()
            .filter(|&u| u > lo + slack && u < hi - slack)
            .collect();
        if us.windows(2).all(|w| w[1] > w[0]) != (start.x > end.x) {
            us.reverse();
        }
        return Ok(us.into_iter().map(along).collect());
    }
    let du = start.x - end.x;
    let mut step = surface.periods().0.map_or(Scalar::INFINITY, |p| p / 8.0);
    if let Some(box_) = ring_box {
        if let Some(b) = bound
            .oracle
            .as_ref()
            .and_then(|oracle| oracle.bounds((box_.0.x, box_.1.x), (box_.0.y, box_.1.y)))
        {
            if b.duu > 0.0 {
                step = step.min((8.0 * ctx.chord_error / b.duu).sqrt());
            }
        }
    }
    let cells = (du.abs() / step).ceil();
    let cells = if cells.is_finite() && cells >= 1.0 {
        cells
    } else {
        1.0
    };
    if cells > MAX_CURVED_EDGE_SEGMENTS as Scalar {
        return Err(GeomError::BudgetExceeded {
            resource: "curved edge samples",
        });
    }
    let n = (cells as usize).next_power_of_two();
    if n > MAX_CURVED_EDGE_SEGMENTS {
        return Err(GeomError::BudgetExceeded {
            resource: "curved edge samples",
        });
    }
    Ok((1..n)
        .map(|i| along(end.x + du * (i as Scalar / n as Scalar)))
        .collect())
}

/// The parameter box of a ring's samples.
fn parameter_box(ring: &[Point2]) -> Option<(Point2, Point2)> {
    let first = *ring.first()?;
    Some(
        ring.iter()
            .fold((first, first), |(lo, hi), p| (lo.min(*p), hi.max(*p))),
    )
}

/// A triangle into a pole that stands for a whole quad of a grid's pole row
/// (#463, ADR 0112). The quad's corners are the pole corner
/// `vertices[corner]`, the far corners `vertices[corner + 1]` (on the pole
/// corner's meridian) and `vertices[corner + 2]`, and `other`, the pole
/// corner on the second far corner's meridian; both pole corners carry the
/// pole's mesh vertex, so the quad's bilinear interpolant of the mesh lies
/// on the triangle.
#[derive(Debug, Clone, Copy)]
pub(super) struct Fan {
    corner: u8,
    other: Point2,
}

impl Fan {
    /// The quad's corners in order: pole, first far, second far, other
    /// pole.
    fn quad(self, vertices: [SurfaceVertex; 3]) -> [Point2; 4] {
        let p = usize::from(self.corner);
        [
            vertices[p].uv,
            vertices[(p + 1) % 3].uv,
            vertices[(p + 2) % 3].uv,
            self.other,
        ]
    }

    /// The edge opposite the pole, the only one a fan splits.
    pub(super) fn opposite(self) -> usize {
        (usize::from(self.corner) + 1) % 3
    }
}

/// A fan's bound ([`FaceBound::fan`]), split by what splitting reduces.
pub(super) struct FanTerms {
    /// The term along the pole row, which splitting the far edge quarters.
    pub(super) along: Scalar,
    /// The term across the row plus the corners' offsets.
    pub(super) rest: Scalar,
}

/// The fans of a grid's bottom (`v` lowest) and top pole rows, or `None`
/// for a row that is no pole row or is demoted (#463, ADR 0112). A border
/// lattice row whose vertices are all one mesh vertex is a pole row; each
/// of its cells is one fan, from the nominal triangle that does not
/// collapse. The row is demoted, all of it, when a cell's bound without
/// its term along the row, `dvv dv^2 / 8` and the corners' offsets,
/// misses the budget or has no smooth bound: splitting the edge opposite
/// the pole cannot reduce those.
pub(super) fn fan_rows(
    (mesh, budget): (&TriMesh, Scalar),
    surface: Support<'_>,
    bound: &FaceBound<'_>,
    patch: &GridPatch,
    cell: impl Fn(usize, usize) -> [SurfaceVertex; 4],
    boundary_edges: &HashSet<LocalEdgeKey>,
) -> GeomResult<[Option<Vec<SurfaceTriangle>>; 2]> {
    let pole_row = |iv: usize, top: bool| {
        (0..patch.nu).all(|iu| {
            let [a, b, c, d] = cell(iv, iu);
            if top {
                c.mesh == d.mesh && d.mesh == cell(iv, 0)[3].mesh
            } else {
                a.mesh == b.mesh && a.mesh == cell(iv, 0)[0].mesh
            }
        })
    };
    let rows = [pole_row(0, false), pole_row(patch.nv - 1, true)];
    let mut fans = [None, None];
    for (side, slot) in fans.iter_mut().enumerate() {
        // One cell high between two poles: no cell has a far edge.
        if !rows[side] || (patch.nv == 1 && rows == [true, true]) {
            continue;
        }
        let mut row = Vec::with_capacity(patch.nu);
        for iu in 0..patch.nu {
            let triangle = if side == 0 {
                let [a, b, c, d] = cell(0, iu);
                fan_triangle([b, c, d], 0, a.uv, boundary_edges)
            } else {
                let [a, b, c, d] = cell(patch.nv - 1, iu);
                fan_triangle([d, a, b], 0, c.uv, boundary_edges)
            };
            match bound.fan(mesh, surface, &triangle)? {
                Some(terms) if terms.rest <= budget => row.push(triangle),
                _ => break,
            }
        }
        if row.len() == patch.nu {
            *slot = Some(row);
        }
    }
    Ok(fans)
}

fn fan_triangle(
    vertices: [SurfaceVertex; 3],
    corner: u8,
    other: Point2,
    boundary_edges: &HashSet<LocalEdgeKey>,
) -> SurfaceTriangle {
    SurfaceTriangle {
        fan: Some(Fan { corner, other }),
        ..surface_triangle(vertices, boundary_edges)
    }
}

/// Split a fan at the middle of its far edge into two fans: the bilinear
/// quad's line `s = 1/2` cuts it into two quads, whose pole corners meet at
/// the middle of the pole row's span, the same mesh vertex.
pub(super) fn split_fan(
    state: &mut CurvedMeshState<'_>,
    surface: Support<'_>,
    triangle: SurfaceTriangle,
    cache: &mut HashMap<LocalEdgeKey, SurfaceVertex>,
    next_local: &mut u32,
    out: &mut Vec<SurfaceTriangle>,
) -> GeomResult<()> {
    let Some(fan) = triangle.fan else {
        return Err(GeomError::Degenerate(
            "a nominal curved triangle split as a fan".to_owned(),
        ));
    };
    let p = usize::from(fan.corner);
    let (pole, near, far) = (
        triangle.vertices[p],
        triangle.vertices[(p + 1) % 3],
        triangle.vertices[(p + 2) % 3],
    );
    let middle = refinement_midpoint(state, surface, near, far, cache, next_local)?;
    let split_pole = SurfaceVertex {
        uv: (pole.uv + fan.other) * 0.5,
        ..pole
    };
    let edge = |i: usize| triangle.boundary[(p + i) % 3];
    let depth = triangle.depth + 1;
    out.extend([
        SurfaceTriangle {
            vertices: [pole, near, middle],
            boundary: [edge(0), false, false],
            depth,
            fan: Some(Fan {
                corner: 0,
                other: split_pole.uv,
            }),
        },
        SurfaceTriangle {
            vertices: [split_pole, middle, far],
            boundary: [false, false, edge(2)],
            depth,
            fan: Some(Fan {
                corner: 0,
                other: fan.other,
            }),
        },
    ]);
    Ok(())
}

impl FaceBound<'_> {
    /// The bound of a fan over its parameter quad (#463, ADR 0112), or
    /// `None` where the quad's box has no smooth second-derivative bound.
    ///
    /// With `A, B` the quad's pole corners and `D, C` its far corners on
    /// their meridians, map the unit square onto the quad bilinearly,
    /// `phi(s, t) = (1 - t) ((1 - s) A + s B) + t ((1 - s) D + s C)`, and
    /// let `F = S o phi`. Interpolating `F` linearly in `s` and then in `t`
    /// errs by at most `(sup|F_ss| + sup|F_tt|) / 8`: each linear step's
    /// error is `s (1 - s) / 2` times a second derivative, and the second
    /// step only blends the first's errors on its two edges. `phi` is
    /// linear along each of `s` and `t`, so `F_ss = S''[phi_s, phi_s]`,
    /// `phi_s` a blend of `B - A` and `C - D` (`phi_t` of `D - A` and `C -
    /// B`), and `|F_ss| <= A d_u^2 + 2 B d_u d_v + C d_v^2`, with `A, B, C`
    /// bounding `|S_uu|, |S_uv|, |S_vv|` over the quad's box and `d` the
    /// larger component of the two. Replacing the corners' surface points by
    /// the mesh vertices moves the interpolant by at most their largest
    /// offset. Both pole corners carry the pole's vertex `P`, so the
    /// interpolant of the mesh is `(1 - t) P + t ((1 - s) X_D + s X_C)`:
    /// every point of the quad's surface is within the bound of the fan's
    /// triangle, and every point of the triangle is reached. The term along
    /// the row is `along`; `rest` is the term across it plus the offsets,
    /// which splitting the edge opposite the pole does not reduce.
    pub(super) fn fan(
        &self,
        mesh: &TriMesh,
        surface: Support<'_>,
        triangle: &SurfaceTriangle,
    ) -> GeomResult<Option<FanTerms>> {
        let (Some(fan), Some(oracle)) = (triangle.fan, self.oracle.as_ref()) else {
            return Ok(None);
        };
        let [a, d, c, b] = fan.quad(triangle.vertices);
        let (lo, hi) = (a.min(b).min(c).min(d), a.max(b).max(c).max(d));
        let Some(bounds) = oracle.bounds((lo.x, hi.x), (lo.y, hi.y)) else {
            return Ok(None);
        };
        if !bounds.smooth {
            return Ok(None);
        }
        let reach = |x: Point2, y: Point2| (x.x.abs().max(y.x.abs()), x.y.abs().max(y.y.abs()));
        let q = |(du, dv): (Scalar, Scalar)| {
            0.125 * (bounds.duu * du * du + 2.0 * bounds.duv * du * dv + bounds.dvv * dv * dv)
        };
        let along = q(reach(b - a, c - d));
        let across = q(reach(d - a, c - b));
        let p = usize::from(fan.corner);
        let pole = mesh.positions[triangle.vertices[p].mesh as usize];
        let mut offset: Scalar = 0.0;
        for (uv, at) in [
            (a, pole),
            (b, pole),
            (
                d,
                mesh.positions[triangle.vertices[(p + 1) % 3].mesh as usize],
            ),
            (
                c,
                mesh.positions[triangle.vertices[(p + 2) % 3].mesh as usize],
            ),
        ] {
            offset = offset.max((surface.evaluate(uv.x, uv.y)? - at).length());
        }
        let rest = across + offset;
        if !(along.is_finite() && rest.is_finite()) {
            return Ok(None);
        }
        Ok(Some(FanTerms { along, rest }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep::{refine_curved_face, EdgeSamples, WeldedVertices};
    use crate::deviation::DeviationBound;
    use axiolid_core::{Frame3, Point3, Vec3};
    use axiolid_measure::proximity::closest_point_on_triangle;
    use axiolid_surface::{Cone, Sphere, Surface};
    use std::f64::consts::{FRAC_PI_2, TAU};

    fn frame() -> Frame3 {
        Frame3 {
            origin: Point3::new(0.1, -0.2, 0.3),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        }
    }

    /// A sphere with its north pole at `v = pi / 2`, and a cone with its apex
    /// at `v = 0.8`.
    fn surfaces() -> [(Surface, Scalar); 2] {
        let (radius, height) = (0.5, 0.8);
        [
            (
                Surface::Sphere(Sphere {
                    frame: frame(),
                    radius,
                }),
                FRAC_PI_2,
            ),
            (
                Surface::Cone(Cone {
                    frame: frame(),
                    radius,
                    semi_angle: -(radius / height).atan(),
                }),
                height,
            ),
        ]
    }

    fn at(surface: &Surface, u: Scalar, v: Scalar) -> Point3 {
        axiolid_reference::surface::evaluate(surface, u, v).unwrap()
    }

    /// The distance from `p` to the triangle.
    fn distance(p: Point3, triangle: [Point3; 3]) -> Scalar {
        closest_point_on_triangle(p, triangle)
            .map(|q| p.distance(q))
            .unwrap()
    }

    /// I7 against dense closed-form sampling: every surface point over a
    /// fan's quad lies within its bound of the fan's triangle, for quads
    /// next to a sphere's pole and a cone's apex, small and large, square
    /// and slanted; and the bound is of second order, a quarter of it for
    /// half the quad.
    #[test]
    fn a_fan_bounds_its_whole_quad() {
        for (surface, pole_v) in surfaces() {
            let support = Support::Surface(&surface);
            let bound = FaceBound::new(support);
            let mut previous: Option<Scalar> = None;
            for (du, dv) in [
                (1.6, 0.6),
                (0.8, 0.3),
                (0.4, 0.15),
                (0.2, 0.075),
                (0.05, 0.4),
            ] {
                let u0 = 0.3;
                let (a, b) = (Point2::new(u0, pole_v), Point2::new(u0 + du, pole_v));
                let (d, c) = (
                    Point2::new(u0, pole_v - dv),
                    Point2::new(u0 + du, pole_v - dv),
                );
                let mesh = TriMesh {
                    positions: vec![
                        at(&surface, d.x, d.y),
                        at(&surface, c.x, c.y),
                        at(&surface, a.x, a.y),
                    ],
                    indices: Vec::new(),
                    normals: None,
                    attributes: Vec::new(),
                };
                let vertex = |uv: Point2, mesh: u32| SurfaceVertex {
                    uv,
                    mesh,
                    local: mesh,
                };
                let fan = SurfaceTriangle {
                    vertices: [vertex(a, 2), vertex(d, 0), vertex(c, 1)],
                    boundary: [false; 3],
                    depth: 0,
                    fan: Some(Fan {
                        corner: 0,
                        other: b,
                    }),
                };
                let terms = bound.fan(&mesh, support, &fan).unwrap().unwrap();
                let q = terms.along + terms.rest;
                let triangle = [2, 0, 1].map(|i| mesh.positions[i]);
                let mut measured: Scalar = 0.0;
                for i in 0..=60 {
                    for j in 0..=60 {
                        let (s, t) = (i as Scalar / 60.0, j as Scalar / 60.0);
                        let p = at(&surface, u0 + du * s, pole_v - dv * t);
                        measured = measured.max(distance(p, triangle));
                    }
                }
                assert!(
                    measured <= q,
                    "{surface:?} {du} {dv}: measured {measured:e} above {q:e}"
                );
                if let (Some(previous), true) = (previous, du < 1.0 && dv < 0.4) {
                    assert!(q <= 0.3 * previous, "{du} {dv}: {q:e} after {previous:e}");
                }
                previous = Some(q);
            }
        }
    }

    /// A fan's corners off the surface (a far corner, the pole vertex) widen
    /// its bound by their offsets, which the measured deviation needs.
    #[test]
    fn a_fan_bound_holds_its_corners_off_the_surface() {
        let (surface, pole_v) = surfaces()[0].clone();
        let support = Support::Surface(&surface);
        let (u0, du, dv) = (0.3, 0.2, 0.075);
        let (a, b) = (Point2::new(u0, pole_v), Point2::new(u0 + du, pole_v));
        let (d, c) = (
            Point2::new(u0, pole_v - dv),
            Point2::new(u0 + du, pole_v - dv),
        );
        for (corner, lift) in [(0, 0.01), (2, 0.004)] {
            let mut positions = vec![
                at(&surface, d.x, d.y),
                at(&surface, c.x, c.y),
                at(&surface, a.x, a.y),
            ];
            positions[corner] += Vec3::new(0.3, -0.2, 0.9).normalize() * lift;
            let mesh = TriMesh {
                positions,
                indices: Vec::new(),
                normals: None,
                attributes: Vec::new(),
            };
            let vertex = |uv: Point2, mesh: u32| SurfaceVertex {
                uv,
                mesh,
                local: mesh,
            };
            let fan = SurfaceTriangle {
                vertices: [vertex(a, 2), vertex(d, 0), vertex(c, 1)],
                boundary: [false; 3],
                depth: 0,
                fan: Some(Fan {
                    corner: 0,
                    other: b,
                }),
            };
            let terms = FaceBound::new(support)
                .fan(&mesh, support, &fan)
                .unwrap()
                .unwrap();
            let triangle = [2, 0, 1].map(|i| mesh.positions[i]);
            let mut measured: Scalar = 0.0;
            for i in 0..=60 {
                for j in 0..=60 {
                    let (s, t) = (i as Scalar / 60.0, j as Scalar / 60.0);
                    let p = at(&surface, u0 + du * s, pole_v - dv * t);
                    measured = measured.max(distance(p, triangle));
                }
            }
            let q = terms.along + terms.rest;
            assert!(measured <= q, "corner {corner}: {measured:e} above {q:e}");
            assert!(
                measured > q - lift,
                "corner {corner}: the offset is not needed, {measured:e} against {q:e}"
            );
        }
    }

    /// A pole row whose cells' term across the row misses the budget is
    /// demoted, all of it, to the nominal triangles (I7), which are still
    /// bounded, and certified, by the triangle bound (I8); with a budget it
    /// meets, the row is fans.
    #[test]
    fn a_demoted_pole_row_is_still_certified() {
        let (surface, pole_v) = surfaces()[0].clone();
        let support = Support::Surface(&surface);
        let (nu, dv) = (4, 0.3);
        let du = TAU / 8.0 / nu as Scalar;
        let rim = pole_v - dv;
        let mut positions: Vec<Vec3> = (0..=nu)
            .map(|i| at(&surface, du * i as Scalar, rim))
            .collect();
        positions.push(at(&surface, 0.0, pole_v));
        let pole = nu as u32 + 1;
        let lattice: Vec<SurfaceVertex> = (0..=nu)
            .map(|i| SurfaceVertex {
                uv: Point2::new(du * i as Scalar, rim),
                mesh: i as u32,
                local: i as u32,
            })
            .chain((0..=nu).map(|i| SurfaceVertex {
                uv: Point2::new(du * i as Scalar, pole_v),
                mesh: pole,
                local: (nu + 1 + i) as u32,
            }))
            .collect();
        let patch = GridPatch {
            nu,
            nv: 1,
            u_start: 0.0,
            v_start: rim,
            du,
            dv,
            lattice: None,
        };
        let cell = |iv: usize, iu: usize| {
            let row = |j: usize, i: usize| lattice[j * (nu + 1) + i];
            [
                row(iv, iu),
                row(iv, iu + 1),
                row(iv + 1, iu + 1),
                row(iv + 1, iu),
            ]
        };
        let mut mesh = TriMesh {
            positions,
            indices: Vec::new(),
            normals: None,
            attributes: Vec::new(),
        };
        let bound = FaceBound::new(support);
        let edges = HashSet::new();
        let tight = fan_rows((&mesh, 1e-3), support, &bound, &patch, cell, &edges).unwrap();
        assert!(tight[0].is_none() && tight[1].is_none(), "demoted at 1e-3");
        let loose = fan_rows((&mesh, 1e-2), support, &bound, &patch, cell, &edges).unwrap();
        assert_eq!(loose[1].as_ref().map(Vec::len), Some(nu), "fans at 1e-2");

        let triangles: Vec<SurfaceTriangle> = (0..nu)
            .flat_map(|iu| {
                let [a, b, c, d] = cell(0, iu);
                [
                    surface_triangle([a, b, d], &edges),
                    surface_triangle([b, c, d], &edges),
                ]
            })
            .collect();
        let mut welded = WeldedVertices::new();
        let mut edge_samples = EdgeSamples::new();
        let mut records = lattice.len();
        let mut state = CurvedMeshState {
            mesh: &mut mesh,
            welded: &mut welded,
            edge_samples: &mut edge_samples,
            chord_error: 1e-3,
            endpoint_limit: 1e-3,
            face_local_vertices: lattice.len(),
            total_curved_records: &mut records,
            new_welded: Vec::new(),
            new_edge_samples: Vec::new(),
        };
        let mut bound = FaceBound::new(support);
        let next_local = lattice.len() as u32;
        refine_curved_face(
            &mut state, support, triangles, next_local, false, &mut bound,
        )
        .unwrap();
        let DeviationBound::Certified(certified) = bound.finish() else {
            panic!("the nominal triangles of a demoted row are bounded");
        };
        let triangles: Vec<[Point3; 3]> = mesh
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .map(|t| t.map(|i| mesh.positions[i as usize]))
            .collect();
        let mut measured: Scalar = 0.0;
        for i in 0..=80 {
            for j in 0..=40 {
                let p = at(
                    &surface,
                    du * nu as Scalar * i as Scalar / 80.0,
                    rim + dv * j as Scalar / 40.0,
                );
                let near = triangles
                    .iter()
                    .map(|&t| distance(p, t))
                    .fold(Scalar::INFINITY, Scalar::min);
                measured = measured.max(near);
            }
        }
        assert!(
            measured <= certified,
            "measured {measured:e} above {certified:e}"
        );
    }

    /// A ring's rim run, its pole samples apart, is what a row copies (I4).
    #[test]
    fn a_row_copies_its_rim_and_unwrapping_leaves_it() {
        let uses = vec![
            UseSamples {
                params: (0..=4).map(|i| Point2::new(i as Scalar, 0.0)).collect(),
                shared: (0..=4).collect(),
                pcurve: true,
                at_pole: [false, false],
            },
            UseSamples {
                params: vec![
                    Point2::new(4.0, 0.0),
                    Point2::new(4.0, 0.5),
                    Point2::new(4.0, 1.0),
                ],
                shared: vec![4, 5, 9],
                pcurve: true,
                at_pole: [false, true],
            },
            UseSamples {
                params: vec![
                    Point2::new(0.0, 1.0),
                    Point2::new(0.0, 0.5),
                    Point2::new(0.0, 0.0),
                ],
                shared: vec![9, 6, 0],
                pcurve: true,
                at_pole: [true, false],
            },
        ];
        let chains = [Chain {
            uses,
            anchored: true,
        }];
        let runs = level_runs(&chains);
        assert_eq!(runs.len(), 1, "{runs:?}");
        assert_eq!(runs[0].0, vec![0.0, 1.0, 2.0, 3.0, 4.0]);
        // The row from the chain's end back to its start, unwrapped by
        // continuity on a periodic support, keeps its `u` values.
        let mut ring: Vec<Point2> = [4.0, 3.0, 2.0, 1.0, 0.0]
            .iter()
            .map(|&u| Point2::new(u, 1.0))
            .collect();
        let before = ring.clone();
        let sphere = &surfaces()[0].0;
        super::super::unwrap_parameter_ring(Support::Surface(sphere), &mut ring, None).unwrap();
        assert_eq!(ring, before);
    }
}
