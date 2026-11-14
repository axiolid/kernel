//! The sources frontmost beside a piece.
//!
//! For one side of a piece the candidates are the triangles that cover the
//! side point, are not removed and have a kept depth there
//! (`far <= W <= near`, a tie with the bound counts), and the caps of the
//! solids the near plane cuts whose ray parity at the point is odd. The
//! frontmost have the greatest depth; every candidate of that depth is kept,
//! because equal depth never occludes.
//!
//! When a frontmost triangle is see-through (ADR 0111), the side also gets
//! the stack seen through it: the candidates grouped into depth layers by
//! exactly equal depth, front to back, up to and including the first layer
//! that holds an occluding surface (an opaque triangle or a cap). A frontmost
//! layer that already holds one is the whole stack. Otherwise the frontmost
//! occluding layer among the other candidates is found by the same linear
//! pass as the frontmost, the see-through candidates strictly behind it are
//! dropped, those at its depth join it, and the rest are inserted into layers
//! by binary search. Each of these extra depth comparisons counts against the
//! pair budget. A see-through front source gets the stack without its own
//! layer-0 entry, interned in [`Stacks`]; every other source, and every side
//! without a see-through front, gets the empty stack 0, so a request without
//! see-through triangles makes no extra comparison.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_exact_arrangement::{ExactPoints, PointId};
use axiolid_mesh_hidden_lines_contract::FaceSource;

use super::side::{sign_at, SidePoint};
use super::{Faces, Work};
use crate::broad::Query;
use crate::view::Bounds;

/// A source of a drawing: a mesh (by position in the request) and what part
/// of it.
pub(crate) type Key = (usize, FaceSource);

/// A source and the id of the stack seen through it in [`Stacks`]: what a
/// region is a region of.
pub(crate) type Class = (Key, u32);

/// One surface of a stack: the mesh's position, its source and its layer.
pub(crate) type Entry = (usize, FaceSource, u32);

/// The frontmost sources beside one side of a piece.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Label {
    /// The sources, ascending.
    pub keys: Vec<Key>,
    /// Per source, the frontmost triangles (scene indices, ascending); empty
    /// for a cap.
    pub fronts: Vec<Vec<usize>>,
    /// Per source, the id of the stack seen through it; empty when every
    /// stack is the empty one, so a label has one form.
    pub behind: Vec<u32>,
    /// Whether a layer behind the frontmost holds two or more sources.
    pub tied: bool,
}

impl Label {
    /// The front triangles of `key`.
    pub(crate) fn fronts_of(&self, key: &Key) -> &[usize] {
        self.keys
            .binary_search(key)
            .map_or(&[][..], |at| &self.fronts[at])
    }

    /// The stack of the source at `at`.
    pub(crate) fn stack(&self, at: usize) -> u32 {
        self.behind.get(at).copied().unwrap_or(0)
    }

    /// The classes, ascending.
    pub(crate) fn classes(&self) -> impl Iterator<Item = Class> + '_ {
        self.keys
            .iter()
            .enumerate()
            .map(|(at, &key)| (key, self.stack(at)))
    }

    /// Whether `class` is one of the classes.
    pub(crate) fn has(&self, class: &Class) -> bool {
        self.keys
            .binary_search(&class.0)
            .is_ok_and(|at| self.stack(at) == class.1)
    }

    /// Whether both labels have the same classes.
    pub(crate) fn same_classes(&self, other: &Self) -> bool {
        self.keys == other.keys && self.behind == other.behind
    }
}

/// The distinct stacks of one request, each named by an id; id 0 is the
/// empty stack. Allocation failures are the `memory` budget refusal.
#[derive(Debug)]
pub(crate) struct Stacks {
    lists: Vec<Vec<Entry>>,
    /// The ids ordered by their lists.
    sorted: Vec<u32>,
}

impl Default for Stacks {
    fn default() -> Self {
        Self {
            lists: vec![Vec::new()],
            sorted: vec![0],
        }
    }
}

impl Stacks {
    /// The entries of the stack `id`.
    pub(crate) fn get(&self, id: u32) -> &[Entry] {
        &self.lists[id as usize]
    }

    /// The id of `stack`, added if new.
    fn intern(&mut self, stack: Vec<Entry>) -> GeomResult<u32> {
        let memory = || GeomError::BudgetExceeded { resource: "memory" };
        let lists = &self.lists;
        match self
            .sorted
            .binary_search_by(|&id| lists[id as usize].cmp(&stack))
        {
            Ok(at) => Ok(self.sorted[at]),
            Err(at) => {
                let id = u32::try_from(self.lists.len()).map_err(|_| memory())?;
                self.lists.try_reserve(1).map_err(|_| memory())?;
                self.sorted.try_reserve(1).map_err(|_| memory())?;
                self.lists.push(stack);
                self.sorted.insert(at, id);
                Ok(id)
            }
        }
    }
}

/// A surface that can be frontmost.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Surface {
    Triangle(usize),
    /// The cap of the solid `faces.caps[index]`, at depth `near`.
    Cap(usize),
}

/// The candidates around the start of a piece, shared by its two sides.
struct Around {
    triangles: Vec<usize>,
    /// Triangles that may be beyond `near`, for cap parity.
    beyond: Vec<usize>,
}

fn around(faces: &Faces<'_>, points: &ExactPoints, at: PointId) -> Around {
    let [x_lo, y_lo, x_hi, y_hi] = points.bounds(at);
    let image = Bounds {
        lo: [x_lo, y_lo, 0.0],
        hi: [x_hi, y_hi, 0.0],
    };
    let triangles = faces.broad.candidates(
        faces.tris,
        &Query {
            image,
            w_lo: faces.far.unwrap_or(f64::NEG_INFINITY),
            w_hi: faces.near.unwrap_or(f64::INFINITY),
        },
    );
    let beyond = match faces.near {
        Some(near) if !faces.caps.is_empty() => faces.broad.candidates(
            faces.tris,
            &Query {
                image,
                w_lo: near,
                w_hi: f64::INFINITY,
            },
        ),
        _ => Vec::new(),
    };
    Around { triangles, beyond }
}

/// The labels beside the piece `from -> to`: its left and its right.
pub(crate) fn label_piece(
    faces: &Faces<'_>,
    points: &ExactPoints,
    ends: [PointId; 2],
    work: &mut Work,
) -> GeomResult<[Label; 2]> {
    let near = around(faces, points, ends[0]);
    let left = label_side(
        faces,
        &SidePoint::new(points, ends[0], ends[1], true),
        &near,
        work,
    )?;
    let right = label_side(
        faces,
        &SidePoint::new(points, ends[0], ends[1], false),
        &near,
        work,
    )?;
    Ok([left, right])
}

/// Whether the side point lies outside the crop box.
fn cropped(faces: &Faces<'_>, side: &SidePoint, work: &mut Work) -> GeomResult<bool> {
    for form in &faces.crop {
        if sign_at(side, form, &mut work.escalations)? == Ordering::Less {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether the triangle `t` covers the side point.
fn covers(faces: &Faces<'_>, side: &SidePoint, t: usize, work: &mut Work) -> GeomResult<bool> {
    let clockwise = faces.tris[t].sigma == Ordering::Less;
    let Some(plane) = faces.planes[t].as_ref() else {
        return Ok(false);
    };
    for form in &plane.sides {
        let sign = sign_at(side, form, &mut work.escalations)?;
        let sign = if clockwise { sign.reverse() } else { sign };
        if sign != Ordering::Greater {
            return Ok(false);
        }
    }
    Ok(true)
}

fn label_side(
    faces: &Faces<'_>,
    side: &SidePoint,
    around: &Around,
    work: &mut Work,
) -> GeomResult<Label> {
    if cropped(faces, side, work)? {
        return Ok(Label::default());
    }
    let mut candidates: Vec<Surface> = Vec::new();
    for &t in &around.triangles {
        if !faces.covers(t) {
            continue;
        }
        work.budget.count()?;
        if !covers(faces, side, t, work)? {
            continue;
        }
        let plane = faces.planes[t]
            .as_ref()
            .expect("a covering triangle has a plane");
        let mut kept = true;
        if let Some(near) = faces.near {
            kept &= sign_at(side, &plane.below(near), &mut work.escalations)? != Ordering::Less;
        }
        if let Some(far) = faces.far {
            kept &= sign_at(side, &plane.above(far), &mut work.escalations)? != Ordering::Less;
        }
        if kept {
            candidates.push(Surface::Triangle(t));
        }
    }
    for index in 0..faces.caps.len() {
        if cap_parity(faces, side, index, &around.beyond, work)? {
            candidates.push(Surface::Cap(index));
        }
    }

    let mut best: Vec<Surface> = Vec::new();
    for &surface in &candidates {
        let Some(&first) = best.first() else {
            best.push(surface);
            continue;
        };
        match depth_order(faces, side, surface, first, work)? {
            Ordering::Greater => {
                best.clear();
                best.push(surface);
            }
            Ordering::Equal => best.push(surface),
            Ordering::Less => {}
        }
    }

    let mut grouped: BTreeMap<Key, Vec<usize>> = BTreeMap::new();
    for &surface in &best {
        match surface {
            Surface::Triangle(t) => grouped
                .entry((faces.mesh_of[t], FaceSource::Group(faces.group_of[t])))
                .or_default()
                .push(t),
            Surface::Cap(index) => {
                let mesh = faces.mesh_of[faces.caps[index].tris.start];
                grouped.entry((mesh, FaceSource::Cap)).or_default();
            }
        }
    }
    let mut label = Label::default();
    for (key, mut fronts) in grouped {
        fronts.sort_unstable();
        fronts.dedup();
        label.keys.push(key);
        label.fronts.push(fronts);
    }
    if best.iter().any(|&surface| see_through(faces, surface)) {
        stack(faces, side, &candidates, &best, &mut label, work)?;
    }
    Ok(label)
}

fn see_through(faces: &Faces<'_>, surface: Surface) -> bool {
    matches!(surface, Surface::Triangle(t) if faces.tris[t].see_through)
}

fn key_of(faces: &Faces<'_>, surface: Surface) -> Key {
    match surface {
        Surface::Triangle(t) => (faces.mesh_of[t], FaceSource::Group(faces.group_of[t])),
        Surface::Cap(index) => (faces.mesh_of[faces.caps[index].tris.start], FaceSource::Cap),
    }
}

/// The depth order of two surfaces, counted against the pair budget: only
/// stacks ask for it.
fn counted_order(
    faces: &Faces<'_>,
    side: &SidePoint,
    x: Surface,
    y: Surface,
    work: &mut Work,
) -> GeomResult<Ordering> {
    work.budget.count()?;
    depth_order(faces, side, x, y, work)
}

/// The stacks of `label`, whose frontmost surfaces `best` hold a
/// see-through triangle; `candidates` are all surfaces at the side point.
fn stack(
    faces: &Faces<'_>,
    side: &SidePoint,
    candidates: &[Surface],
    best: &[Surface],
    label: &mut Label,
    work: &mut Work,
) -> GeomResult<()> {
    // Layers front to back; the frontmost is `best`.
    let mut layers: Vec<Vec<Surface>> = vec![best.to_vec()];
    if best.iter().all(|&surface| see_through(faces, surface)) {
        let others = || candidates.iter().copied().filter(|s| !best.contains(s));
        // The frontmost occluding layer behind `best`.
        let mut stop: Vec<Surface> = Vec::new();
        for surface in others().filter(|&s| !see_through(faces, s)) {
            match stop.first() {
                None => stop.push(surface),
                Some(&first) => match counted_order(faces, side, surface, first, work)? {
                    Ordering::Greater => stop = vec![surface],
                    Ordering::Equal => stop.push(surface),
                    Ordering::Less => {}
                },
            }
        }
        let mut between: Vec<Vec<Surface>> = Vec::new();
        for surface in others().filter(|&s| see_through(faces, s)) {
            if let Some(&first) = stop.first() {
                match counted_order(faces, side, surface, first, work)? {
                    Ordering::Less => continue,
                    Ordering::Equal => {
                        stop.push(surface);
                        continue;
                    }
                    Ordering::Greater => {}
                }
            }
            let (mut low, mut high) = (0, between.len());
            let mut joined = false;
            while low < high {
                let mid = low + (high - low) / 2;
                match counted_order(faces, side, surface, between[mid][0], work)? {
                    Ordering::Greater => high = mid,
                    Ordering::Less => low = mid + 1,
                    Ordering::Equal => {
                        between[mid].push(surface);
                        joined = true;
                        break;
                    }
                }
            }
            if !joined {
                between.insert(low, vec![surface]);
            }
        }
        layers.extend(between);
        if !stop.is_empty() {
            layers.push(stop);
        }
    }

    // One entry per source and layer, ordered by layer, mesh and source.
    let mut entries: Vec<Entry> = Vec::new();
    for (layer, surfaces) in layers.iter().enumerate() {
        let mut keys: Vec<Key> = surfaces.iter().map(|&s| key_of(faces, s)).collect();
        keys.sort_unstable();
        keys.dedup();
        label.tied |= layer > 0 && keys.len() > 1;
        let layer =
            u32::try_from(layer).map_err(|_| GeomError::BudgetExceeded { resource: "memory" })?;
        entries.extend(keys.into_iter().map(|(mesh, source)| (mesh, source, layer)));
    }
    let mut behind = vec![0; label.keys.len()];
    for (at, key) in label.keys.iter().enumerate() {
        // A group is all see-through or all opaque, and a cap occludes.
        let glazed = label.fronts[at].iter().any(|&t| faces.tris[t].see_through);
        if glazed {
            let own: Vec<Entry> = entries
                .iter()
                .copied()
                .filter(|&(mesh, source, layer)| (layer, mesh, source) != (0, key.0, key.1))
                .collect();
            behind[at] = work.stacks.intern(own)?;
        }
    }
    if behind.iter().any(|&id| id != 0) {
        label.behind = behind;
    }
    Ok(())
}

/// Whether the point of the near plane above the side point lies inside the
/// solid `faces.caps[index]`: an odd number of its triangles lie strictly
/// beyond the plane there. A triangle in the plane is not a crossing, and a
/// triangle seen edge-on covers nothing.
fn cap_parity(
    faces: &Faces<'_>,
    side: &SidePoint,
    index: usize,
    beyond: &[usize],
    work: &mut Work,
) -> GeomResult<bool> {
    let Some(near) = faces.near else {
        return Ok(false);
    };
    let range = &faces.caps[index].tris;
    let first = beyond.partition_point(|&t| t < range.start);
    let last = beyond.partition_point(|&t| t < range.end);
    let mut odd = false;
    for &t in &beyond[first..last] {
        let Some(plane) = faces.planes[t].as_ref() else {
            continue;
        };
        work.budget.count()?;
        if covers(faces, side, t, work)?
            && sign_at(side, &plane.below(near), &mut work.escalations)? == Ordering::Less
        {
            odd = !odd;
        }
    }
    Ok(odd)
}

/// `Greater` when `x` is in front of `y` at the side point.
fn depth_order(
    faces: &Faces<'_>,
    side: &SidePoint,
    x: Surface,
    y: Surface,
    work: &mut Work,
) -> GeomResult<Ordering> {
    let plane = |t: usize| {
        faces.planes[t]
            .as_ref()
            .ok_or_else(|| GeomError::Degenerate("an edge-on triangle has no depth".into()))
    };
    let near = || {
        faces
            .near
            .ok_or_else(|| GeomError::Degenerate("a cap without a near plane".into()))
    };
    let escalations = &mut work.escalations;
    match (x, y) {
        (Surface::Triangle(a), Surface::Triangle(b)) => {
            sign_at(side, &plane(a)?.in_front_of(plane(b)?), escalations)
        }
        // `near - W_t`: positive when the cap is in front.
        (Surface::Cap(_), Surface::Triangle(t)) => {
            sign_at(side, &plane(t)?.below(near()?), escalations)
        }
        (Surface::Triangle(t), Surface::Cap(_)) => {
            sign_at(side, &plane(t)?.below(near()?), escalations).map(Ordering::reverse)
        }
        (Surface::Cap(_), Surface::Cap(_)) => Ok(Ordering::Equal),
    }
}

#[cfg(test)]
mod tests {
    use axiolid_core::{Aabb2, Point2};
    use axiolid_mesh_hidden_lines_contract::{FaceGroup, MeshRole};

    use super::*;
    use crate::caps::Cap;
    use crate::faces::testing::{plan, sheet, solid, Setup};

    const GROUP: FaceSource = FaceSource::Group(FaceGroup(0));

    /// The labels beside the piece `from -> to`.
    fn beside(setup: &Setup, from: (f64, f64), to: (f64, f64)) -> [Label; 2] {
        let mut work = setup.work(usize::MAX);
        let faces = setup.faces(&mut work);
        let mut points = setup.points();
        let mut push = |(x, y)| points.store.push_flat(Point2::new(x, y)).expect("point");
        let ends = [push(from), push(to)];
        label_piece(&faces, &points.store, ends, &mut work).expect("labels")
    }

    fn keys(label: &Label) -> Vec<Key> {
        label.keys.clone()
    }

    #[test]
    fn a_sheet_labels_what_it_covers_and_leaves_the_rest_empty() {
        let setup = Setup::new(
            &[(MeshRole::Sheet, sheet(0.0, 0.0, 4.0, 4.0, 0.0))],
            plan(),
            None,
        );
        let [left, right] = beside(&setup, (1.0, 1.0), (2.0, 1.0));
        assert_eq!(keys(&left), vec![(0, GROUP)]);
        assert_eq!(left, right);
        assert!(!left.fronts_of(&(0, GROUP)).is_empty());
        let [left, right] = beside(&setup, (5.0, 5.0), (6.0, 5.0));
        assert!(left.keys.is_empty() && right.keys.is_empty());
    }

    #[test]
    fn the_side_point_beside_a_triangle_edge_is_inside_on_one_side_only() {
        let setup = Setup::new(
            &[(MeshRole::Sheet, sheet(0.0, 0.0, 4.0, 4.0, 0.0))],
            plan(),
            None,
        );
        // Up along the left edge x = 0: left is x < 0, right is the sheet.
        let [left, right] = beside(&setup, (0.0, 1.0), (0.0, 2.0));
        assert!(left.keys.is_empty());
        assert_eq!(keys(&right), vec![(0, GROUP)]);
        // Along the diagonal both sides are inside, in different triangles.
        let [left, right] = beside(&setup, (1.0, 1.0), (2.0, 2.0));
        assert_eq!(keys(&left), keys(&right));
        assert_ne!(left.fronts, right.fronts);
    }

    #[test]
    fn labelling_counts_its_candidate_pairs_against_the_budget() {
        let setup = Setup::new(
            &[(MeshRole::Sheet, sheet(0.0, 0.0, 4.0, 4.0, 0.0))],
            plan(),
            None,
        );
        let mut work = setup.work(0);
        let faces = setup.faces(&mut work);
        let mut points = setup.points();
        let mut push = |(x, y)| points.store.push_flat(Point2::new(x, y)).expect("point");
        let ends = [push((1.0, 1.0)), push((2.0, 1.0))];
        let error = label_piece(&faces, &points.store, ends, &mut work).expect_err("over budget");
        assert!(matches!(error, GeomError::BudgetExceeded { .. }), "{error}");
    }

    #[test]
    fn the_nearer_sheet_wins_and_equal_depth_keeps_both() {
        let low = sheet(0.0, 0.0, 4.0, 4.0, 0.0);
        let setup = Setup::new(
            &[
                (MeshRole::Sheet, low.clone()),
                (MeshRole::Sheet, sheet(1.0, 1.0, 3.0, 3.0, 1.0)),
            ],
            plan(),
            None,
        );
        assert_eq!(
            keys(&beside(&setup, (2.0, 2.0), (2.5, 2.0))[0]),
            vec![(1, GROUP)]
        );
        assert_eq!(
            keys(&beside(&setup, (0.5, 0.5), (1.0, 0.5))[0]),
            vec![(0, GROUP)]
        );
        let tie = Setup::new(
            &[
                (MeshRole::Sheet, low),
                (MeshRole::Sheet, sheet(1.0, 1.0, 3.0, 3.0, 0.0)),
            ],
            plan(),
            None,
        );
        assert_eq!(
            keys(&beside(&tie, (2.0, 2.0), (2.5, 2.0))[0]),
            vec![(0, GROUP), (1, GROUP)]
        );
    }

    #[test]
    fn depth_at_the_near_plane_counts_as_kept() {
        let at = |z| {
            Setup::new(
                &[(MeshRole::Sheet, sheet(0.0, 0.0, 4.0, 4.0, z))],
                plan().with_near(1.0),
                None,
            )
        };
        assert_eq!(
            keys(&beside(&at(1.0), (1.0, 1.0), (2.0, 1.0))[0]),
            vec![(0, GROUP)]
        );
        assert!(beside(&at(2.0), (1.0, 1.0), (2.0, 1.0))[0].keys.is_empty());
        let far = Setup::new(
            &[(MeshRole::Sheet, sheet(0.0, 0.0, 4.0, 4.0, -1.0))],
            plan().with_far(-1.0),
            None,
        );
        assert_eq!(
            keys(&beside(&far, (1.0, 1.0), (2.0, 1.0))[0]),
            vec![(0, GROUP)]
        );
    }

    #[test]
    fn a_cut_solid_gives_a_cap_not_its_top_face() {
        // The top at z = 2 is beyond near = 1: it counts for the parity of the
        // ray and is never a candidate.
        let setup = Setup::new(
            &[(MeshRole::Solid, solid([0.0; 3], [2.0, 2.0, 2.0]))],
            plan().with_near(1.0),
            None,
        );
        let [left, _] = beside(&setup, (1.0, 1.0), (1.5, 1.0));
        assert_eq!(keys(&left), vec![(0, FaceSource::Cap)]);
        assert!(left.fronts_of(&(0, FaceSource::Cap)).is_empty());
        assert!(beside(&setup, (3.0, 1.0), (3.5, 1.0))[0].keys.is_empty());
    }

    #[test]
    fn a_triangle_in_the_near_plane_is_not_a_crossing() {
        // Force the cap of a box whose top lies exactly in the plane.
        let setup = Setup::new(
            &[(MeshRole::Solid, solid([0.0; 3], [2.0, 2.0, 1.0]))],
            plan().with_near(1.0),
            None,
        );
        let mut work = setup.work(usize::MAX);
        let mut faces = setup.faces(&mut work);
        let forced = [Cap { tris: 0..12 }];
        faces.caps = &forced;
        let mut points = setup.points();
        let mut push = |x, y| points.store.push_flat(Point2::new(x, y)).expect("point");
        let ends = [push(1.0, 1.0), push(1.5, 1.0)];
        let [left, _] = label_piece(&faces, &points.store, ends, &mut work).expect("labels");
        // The top is the only candidate; there is no cap.
        assert_eq!(keys(&left), vec![(0, GROUP)]);
    }

    #[test]
    fn the_crop_empties_what_lies_outside() {
        let region = Aabb2 {
            min: Point2::new(0.0, 0.0),
            max: Point2::new(1.0, 1.0),
        };
        let setup = Setup::new(
            &[(MeshRole::Sheet, sheet(0.0, 0.0, 4.0, 4.0, 0.0))],
            plan(),
            Some(region),
        );
        assert!(beside(&setup, (2.0, 2.0), (3.0, 2.0))[0].keys.is_empty());
        assert_eq!(
            keys(&beside(&setup, (0.25, 0.5), (0.75, 0.5))[0]),
            vec![(0, GROUP)]
        );
        // Along the crop's own side: inside on one side only.
        let [left, right] = beside(&setup, (1.0, 0.25), (1.0, 0.75));
        assert_eq!(keys(&left), vec![(0, GROUP)]);
        assert!(right.keys.is_empty());
    }

    /// Sheets `(x0, y0, x1, y1, z, see-through)` in plan; the left label
    /// beside the piece `from -> to`, the stack of each of its sources and the
    /// pairs counted.
    fn stacked(
        sheets: &[(f64, f64, f64, f64, f64, bool)],
        from: (f64, f64),
        to: (f64, f64),
    ) -> (Label, Vec<Vec<Entry>>, usize) {
        let parts: Vec<_> = sheets
            .iter()
            .map(|&(x0, y0, x1, y1, z, _)| (MeshRole::Sheet, sheet(x0, y0, x1, y1, z)))
            .collect();
        let mut setup = Setup::new(&parts, plan(), None);
        for (index, &(.., glass)) in sheets.iter().enumerate() {
            let range = setup.prepared.ranges[index].triangles.clone();
            for tri in &mut setup.prepared.tris[range] {
                tri.see_through = glass;
            }
        }
        let mut work = setup.work(usize::MAX);
        let faces = setup.faces(&mut work);
        let mut points = setup.points();
        let mut push = |(x, y)| points.store.push_flat(Point2::new(x, y)).expect("point");
        let ends = [push(from), push(to)];
        let [left, _] = label_piece(&faces, &points.store, ends, &mut work).expect("labels");
        let stacks = (0..left.keys.len())
            .map(|at| work.stacks.get(left.stack(at)).to_vec())
            .collect();
        (left, stacks, work.budget.used)
    }

    const FLOOR: (f64, f64, f64, f64, f64, bool) = (0.0, 0.0, 4.0, 4.0, 0.0, false);

    #[test]
    fn an_opaque_front_has_no_stack() {
        let (label, stacks, _) = stacked(
            &[FLOOR, (1.0, 1.0, 3.0, 3.0, 1.0, false)],
            (2.0, 2.0),
            (2.5, 2.0),
        );
        assert_eq!(keys(&label), vec![(1, GROUP)]);
        assert!(label.behind.is_empty() && !label.tied);
        assert_eq!(stacks, [Vec::<Entry>::new()]);
    }

    #[test]
    fn a_see_through_front_lists_down_to_the_first_occluder() {
        let (label, stacks, _) = stacked(
            &[FLOOR, (1.0, 1.0, 3.0, 3.0, 1.0, true)],
            (2.0, 2.0),
            (2.5, 2.0),
        );
        assert_eq!(keys(&label), vec![(1, GROUP)]);
        assert_eq!(stacks, [vec![(0, GROUP, 1)]]);
        // A see-through layer between, an opaque sheet ending the stack and a
        // see-through sheet behind it, which is dropped.
        let scene = [
            FLOOR,
            (0.0, 0.0, 4.0, 4.0, 0.5, true),
            (1.0, 1.0, 3.0, 3.0, 1.0, false),
            (0.0, 0.0, 4.0, 4.0, 1.5, true),
            (0.0, 0.0, 4.0, 4.0, 2.0, true),
        ];
        let (_, stacks, _) = stacked(&scene, (2.0, 2.0), (2.5, 2.0));
        assert_eq!(stacks, [vec![(3, GROUP, 1), (2, GROUP, 2)]]);
        let (_, stacks, _) = stacked(&scene, (0.25, 0.25), (0.5, 0.25));
        assert_eq!(stacks, [vec![(3, GROUP, 1), (1, GROUP, 2), (0, GROUP, 3)]]);
    }

    #[test]
    fn a_see_through_front_tied_with_an_opaque_one_lists_the_tie_only() {
        let (label, stacks, _) = stacked(
            &[
                (0.0, 0.0, 4.0, 4.0, -1.0, false),
                (0.0, 0.0, 4.0, 4.0, 0.0, false),
                (1.0, 1.0, 3.0, 3.0, 0.0, true),
            ],
            (2.0, 2.0),
            (2.5, 2.0),
        );
        assert_eq!(keys(&label), vec![(1, GROUP), (2, GROUP)]);
        assert_eq!(stacks, [vec![], vec![(1, GROUP, 0)]]);
        assert!(!label.tied);
    }

    #[test]
    fn a_tie_behind_the_front_is_counted_and_layers_sort_by_depth() {
        // Two see-through sheets tied at depth 1 and one at 2 behind the
        // front at 3, inserted out of depth order.
        let (label, stacks, _) = stacked(
            &[
                FLOOR,
                (0.0, 0.0, 4.0, 4.0, 1.0, true),
                (0.0, 0.0, 4.0, 4.0, 2.0, true),
                (0.0, 0.0, 4.0, 4.0, 1.0, true),
                (0.0, 0.0, 4.0, 4.0, 3.0, true),
            ],
            (2.0, 2.0),
            (2.5, 2.0),
        );
        assert!(label.tied);
        assert_eq!(
            stacks,
            [vec![
                (2, GROUP, 1),
                (1, GROUP, 2),
                (3, GROUP, 2),
                (0, GROUP, 3)
            ]]
        );
    }

    #[test]
    fn the_stack_comparisons_count_against_the_budget() {
        let glass = [
            FLOOR,
            (1.0, 1.0, 3.0, 3.0, 1.0, true),
            (1.0, 1.0, 3.0, 3.0, 0.5, true),
        ];
        let opaque = glass.map(|(x0, y0, x1, y1, z, _)| (x0, y0, x1, y1, z, false));
        let (_, _, with_stack) = stacked(&glass, (2.0, 2.0), (2.5, 2.0));
        let (_, _, without) = stacked(&opaque, (2.0, 2.0), (2.5, 2.0));
        // The stop search compares nothing (one occluder); the glass at 0.5
        // is compared with the stop, on each side of the piece.
        assert_eq!(with_stack, without + 2);
    }
}
