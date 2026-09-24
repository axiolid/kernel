//! Exact solids made of vertical columns over a planar arrangement (#120).
//!
//! A coaxial boolean with differing spans, and a prism cut by a plane that
//! crosses a cap, are the same shape: the plan is cut into cells, and above
//! each cell the solid occupies a stack of height intervals, each bounded
//! below and above by a plane (flat, or sloped as in ADR 0071). The cells
//! come from one exact [`ArcArrangement`], so every face shares its vertices
//! with every other face by index.
//!
//! # Faces
//!
//! - A cap lies on one plane, facing up or down, over the cells whose stack
//!   has an interval ending on that plane from that side. Its boundary is
//!   whatever the arrangement links for that cell set.
//! - A wall stands on one arrangement piece, over the heights where the two
//!   cells either side disagree about being solid. It is a `Cylinder` over an
//!   arc piece and a plane over a straight one, and it faces the empty side.
//! - A vertical edge stands at an arrangement vertex, between consecutive
//!   heights that any face meets there, so walls meeting at one vertex share
//!   the same vertical edges rather than overlapping ones.
//!
//! # Scope
//!
//! Interval ends that lie on one plane over a cell must stay apart from the
//! other ends over that cell (the caller refuses crossings before building).
//! Each connected piece of solid becomes one `ExactBRep`. Two pieces that
//! touch along an edge are refused, since that edge would bound four faces
//! and the result would not be a manifold solid. A result enclosing a
//! cavity is refused too: void shells are not tessellated in this kernel.

use std::collections::BTreeMap;

use axiolid_brep::{ExactBRep, ExactBRepBuilder, FaceName};
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Frame3, Interval, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Curve3, Ellipse2, Line2, Line3, Sinusoid2};
use axiolid_overlay::ArcRing;
use axiolid_overlay::{ArcArrangement, ArrangementEdgeUse, EdgeSource};
use axiolid_surface::{Cylinder, Plane, Surface};
use axiolid_topology::{
    audit_brep, Edge, EdgeId, EdgeUse, Face, FaceBound, FaceId, Orientation, Shell, ShellId, Solid,
    Vertex, VertexId,
};

use crate::extrude_arc::{
    arc_geometry, circle2_frame, circle_of, sloped_arc_edge, ArcGeometry, Level,
};
use crate::extrude_exact::{add_loop, identity_frame3};
use crate::BACKEND_ID;

/// Refusal of a shape this builder cannot represent, named for the caller.
fn unsupported(input: &'static str) -> GeomError {
    crate::boolean_exact::unsupported(input)
}

/// A broken internal invariant: the builder produced something it should not.
fn contract(detail: impl Into<String>) -> GeomError {
    GeomError::BackendContractViolation {
        backend: BACKEND_ID,
        detail: detail.into(),
    }
}

/// Disjoint sets with path halving; roots are the smallest member.
struct Sets(Vec<usize>);

impl Sets {
    fn new(count: usize) -> Self {
        Self((0..count).collect())
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }

    fn join(&mut self, a: usize, b: usize) {
        let (a, b) = (self.find(a), self.find(b));
        if a != b {
            self.0[a.max(b)] = a.min(b);
        }
    }
}

/// Plan geometry of one arrangement piece, in its own direction.
struct Piece {
    from: Point2,
    to: Point2,
    bulge: Scalar,
    /// Start, middle and end: where two levels are compared along it.
    samples: [Point2; 3],
}

fn piece(arrangement: &ArcArrangement, index: usize) -> GeomResult<Piece> {
    let edge = &arrangement.edges()[index];
    let from = arrangement.vertices()[edge.from];
    let to = arrangement.vertices()[edge.to];
    let mid = if edge.bulge == 0.0 {
        (from + to) * 0.5
    } else {
        let arc = arc_geometry(from, to, edge.bulge)?;
        let radial = from - arc.centre;
        let angle = radial.y.atan2(radial.x) + 0.5 * arc.sweep;
        arc.centre + Vec2::new(angle.cos(), angle.sin()) * arc.radius
    };
    Ok(Piece {
        from,
        to,
        bulge: edge.bulge,
        samples: [from, mid, to],
    })
}

/// The distinct heights over one piece, bottom to top.
///
/// Planes that agree within the tolerance at every sample of the piece are
/// one class: over that piece they are the same surface, so they share one
/// rim edge. Three samples settle it for an arc (an affine height that
/// vanishes at three points of a circle vanishes everywhere) and two for a
/// straight piece.
struct Classes {
    /// Representative plane of each class, lowest class first.
    reps: Vec<usize>,
    /// Plane index -> class index.
    of: BTreeMap<usize, usize>,
}

fn classes(
    planes: &[Level],
    used: impl IntoIterator<Item = usize>,
    piece: &Piece,
    tolerance: Tolerance,
) -> GeomResult<Classes> {
    let mut members: Vec<usize> = used.into_iter().collect();
    members.sort_unstable();
    members.dedup();
    let [_, mid, _] = piece.samples;
    members.sort_by(|&a, &b| {
        planes[a]
            .at(mid)
            .total_cmp(&planes[b].at(mid))
            .then(a.cmp(&b))
    });
    let same = |a: usize, b: usize| {
        piece
            .samples
            .iter()
            .all(|&p| tolerance.eq(planes[a].at(p), planes[b].at(p)))
    };
    let mut reps: Vec<usize> = Vec::new();
    let mut of = BTreeMap::new();
    for plane in members {
        // Members sorted by mid height: a coincident plane is next to the
        // class it joins, or within the run of planes equal at the middle.
        let joined = reps
            .iter()
            .rposition(|&rep| same(rep, plane))
            .filter(|&class| tolerance.eq(planes[reps[class]].at(mid), planes[plane].at(mid)));
        let class = match joined {
            Some(class) => class,
            None => {
                reps.push(plane);
                reps.len() - 1
            }
        };
        of.insert(plane, class);
    }
    // Distinct classes must not swap order along the piece: the caller
    // splits cells wherever two bounding planes cross.
    for pair in reps.windows(2) {
        let (low, high) = (&planes[pair[0]], &planes[pair[1]]);
        for &p in &piece.samples {
            if low.at(p) > high.at(p) + tolerance.linear() {
                return Err(contract(
                    "two bounding planes cross inside one arrangement piece",
                ));
            }
        }
    }
    Ok(Classes { reps, of })
}

/// One height interval of the solid over a cell: plane indices, low first.
pub(crate) type Block = (usize, usize);

/// Names a wall from the input edges under its piece and the solid side's
/// ring membership.
pub(crate) type WallNamer<'a> = &'a dyn Fn(&[EdgeSource], &[bool]) -> Option<FaceName>;

/// What a column solid is made of.
pub(crate) struct Columns<'a> {
    /// The plan subdivision; every face boundary is made of its pieces.
    pub(crate) arrangement: &'a ArcArrangement,
    /// Bounding planes, referred to by index.
    pub(crate) planes: &'a [Level],
    /// The solid over a cell, from which rings contain the cell: disjoint
    /// intervals, lowest first. Empty where there is no solid.
    pub(crate) stack: &'a dyn Fn(&[bool]) -> Vec<Block>,
    /// Name of a cap on `plane`; `up` is true for a cap facing up.
    pub(crate) cap_name: &'a dyn Fn(usize, bool) -> Option<FaceName>,
    /// Name of a wall on a piece with these sources, whose solid lies on
    /// the side where the rings in the given mask are inside.
    pub(crate) wall_name: WallNamer<'a>,
    /// Distances below this are one height.
    pub(crate) tolerance: Tolerance,
}

/// A topological edge, named before any builder exists.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum Key {
    /// Where the plane class represented by `rep` runs along `piece`.
    Rim { piece: usize, rep: usize },
    /// The vertical edge at arrangement vertex `vertex`, from height
    /// cluster `k` up to cluster `k + 1`.
    Vert { vertex: usize, k: usize },
}

/// What a face stands on.
enum Kind {
    /// A cap on `plane`; `plan` holds its region rings, outer first.
    Cap {
        plane: usize,
        up: bool,
        plan: Vec<Vec<ArrangementEdgeUse>>,
    },
    /// A wall on `piece`, with the solid on its left when `left`.
    Wall { piece: usize, left: bool },
}

/// A face before any builder exists: its loops as keyed edge uses.
struct Sketch {
    kind: Kind,
    /// Loops, outer first; `true` walks the key's edge forward.
    bounds: Vec<Vec<(Key, bool)>>,
    name: Option<FaceName>,
}

/// Everything about one arrangement piece the faces need.
struct PieceData {
    geo: Piece,
    classes: Classes,
    left: Vec<bool>,
    right: Vec<bool>,
    left_stack: Vec<Block>,
    right_stack: Vec<Block>,
}

impl PieceData {
    /// The class representative of `plane` on this piece.
    fn rep(&self, plane: usize) -> GeomResult<usize> {
        self.classes
            .of
            .get(&plane)
            .map(|&class| self.classes.reps[class])
            .ok_or_else(|| contract("a face uses a plane its boundary piece does not carry"))
    }

    /// Whether `stack` is solid in the gap above class `gap`.
    fn solid(&self, stack: &[Block], gap: usize) -> bool {
        stack
            .iter()
            .any(|&(lo, hi)| self.classes.of[&lo] <= gap && self.classes.of[&hi] > gap)
    }
}

/// Heights met at each arrangement vertex, clustered: `(lowest, highest)`
/// of each cluster, bottom to top.
struct Heights(Vec<Vec<(Scalar, Scalar)>>);

impl Heights {
    /// Cluster raw heights per vertex. A cluster spans at most the
    /// tolerance; two heights closer than that but split across clusters
    /// would make the vertex ambiguous, so that is refused.
    fn new(raw: Vec<Vec<Scalar>>, tolerance: Tolerance) -> GeomResult<Self> {
        let mut out = Vec::with_capacity(raw.len());
        for mut values in raw {
            values.sort_by(|a, b| a.total_cmp(b));
            let mut clusters: Vec<(Scalar, Scalar)> = Vec::new();
            for value in values {
                match clusters.last_mut() {
                    Some(last) if value - last.0 <= tolerance.linear() => last.1 = value,
                    Some(last) if value - last.1 <= tolerance.linear() => {
                        return Err(unsupported(
                            "column heights at one corner closer than the tolerance \
                             but not equal within it",
                        ))
                    }
                    _ => clusters.push((value, value)),
                }
            }
            out.push(clusters);
        }
        Ok(Self(out))
    }

    /// Which cluster at `vertex` holds `z`, which must be a height that
    /// was clustered (same computation, so an exact match).
    fn index(&self, vertex: usize, z: Scalar) -> GeomResult<usize> {
        self.0[vertex]
            .iter()
            .position(|&(lo, hi)| lo <= z && z <= hi)
            .ok_or_else(|| contract("a face corner height was never clustered"))
    }

    /// The height a cluster's vertex sits at.
    fn z(&self, vertex: usize, k: usize) -> Scalar {
        let (lo, hi) = self.0[vertex][k];
        0.5 * (lo + hi)
    }
}

/// A wall before its vertical edges are known.
struct WallPlan {
    piece: usize,
    left: bool,
    low: usize,
    high: usize,
}

/// Faces of a column description, keyed but not yet built.
struct Plan {
    data: Vec<PieceData>,
    heights: Heights,
    sketches: Vec<Sketch>,
}

fn plan(c: &Columns<'_>) -> GeomResult<Plan> {
    let arrangement = c.arrangement;
    let rings = arrangement.ring_count();
    let mut data = Vec::with_capacity(arrangement.edges().len());
    for (index, edge) in arrangement.edges().iter().enumerate() {
        let left: Vec<bool> = (0..rings).map(|r| edge.inside_left(r)).collect();
        let right: Vec<bool> = (0..rings).map(|r| edge.inside_right(r)).collect();
        let left_stack = (c.stack)(&left);
        let right_stack = (c.stack)(&right);
        let geo = piece(arrangement, index)?;
        let used = left_stack
            .iter()
            .chain(&right_stack)
            .flat_map(|&(lo, hi)| [lo, hi]);
        let classes = classes(c.planes, used, &geo, c.tolerance)?;
        data.push(PieceData {
            geo,
            classes,
            left,
            right,
            left_stack,
            right_stack,
        });
    }

    // Walls: per piece, maximal runs of gaps where exactly one side is
    // solid, and on the same side throughout.
    let mut walls = Vec::new();
    for (index, d) in data.iter().enumerate() {
        let gaps = d.classes.reps.len().saturating_sub(1);
        let mut gap = 0;
        while gap < gaps {
            let l = d.solid(&d.left_stack, gap);
            let r = d.solid(&d.right_stack, gap);
            if l == r {
                gap += 1;
                continue;
            }
            let start = gap;
            while gap < gaps
                && d.solid(&d.left_stack, gap) == l
                && d.solid(&d.right_stack, gap) == r
            {
                gap += 1;
            }
            walls.push(WallPlan {
                piece: index,
                left: l,
                low: d.classes.reps[start],
                high: d.classes.reps[gap],
            });
        }
    }

    // Caps: per plane and facing, the regions where a block ends there.
    let mut caps = Vec::new();
    for plane in 0..c.planes.len() {
        for up in [false, true] {
            let ends = |mask: &[bool]| {
                (c.stack)(mask)
                    .iter()
                    .any(|&(lo, hi)| if up { hi == plane } else { lo == plane })
            };
            let regions = arrangement
                .regions(ends)
                .map_err(|error| contract(format!("a cap boundary did not link: {error:?}")))?;
            for region in regions {
                let mut plan = vec![region.outer];
                plan.extend(region.holes);
                caps.push((plane, up, plan));
            }
        }
    }

    // Heights every rim reaches at every vertex.
    let vertices = arrangement.vertices();
    let mut raw = vec![Vec::new(); vertices.len()];
    let mut reach = |piece: usize, rep: usize| {
        let edge = &arrangement.edges()[piece];
        for vertex in [edge.from, edge.to] {
            raw[vertex].push(c.planes[rep].at(vertices[vertex]));
        }
    };
    for wall in &walls {
        reach(wall.piece, wall.low);
        reach(wall.piece, wall.high);
    }
    for (plane, _, plan) in &caps {
        for uses in plan {
            for u in uses {
                reach(u.edge, data[u.edge].rep(*plane)?);
            }
        }
    }
    let heights = Heights::new(raw, c.tolerance)?;

    let mut sketches = Vec::with_capacity(walls.len() + caps.len());
    for wall in walls {
        let edge = &arrangement.edges()[wall.piece];
        let (from, to) = if wall.left {
            (edge.from, edge.to)
        } else {
            (edge.to, edge.from)
        };
        let k =
            |vertex: usize, rep: usize| heights.index(vertex, c.planes[rep].at(vertices[vertex]));
        let (to_low, to_high) = (k(to, wall.low)?, k(to, wall.high)?);
        let (from_low, from_high) = (k(from, wall.low)?, k(from, wall.high)?);
        let rim = |rep: usize| Key::Rim {
            piece: wall.piece,
            rep,
        };
        let mut uses = vec![(rim(wall.low), wall.left)];
        uses.extend((to_low..to_high).map(|k| (Key::Vert { vertex: to, k }, true)));
        uses.push((rim(wall.high), !wall.left));
        uses.extend(
            (from_low..from_high)
                .rev()
                .map(|k| (Key::Vert { vertex: from, k }, false)),
        );
        let d = &data[wall.piece];
        let side = if wall.left { &d.left } else { &d.right };
        sketches.push(Sketch {
            name: (c.wall_name)(&edge.sources, side),
            kind: Kind::Wall {
                piece: wall.piece,
                left: wall.left,
            },
            bounds: vec![uses],
        });
    }
    for (plane, up, plan) in caps {
        let mut bounds = Vec::with_capacity(plan.len());
        for uses in &plan {
            let mut keyed = Vec::with_capacity(uses.len());
            for u in uses {
                let rep = data[u.edge].rep(plane)?;
                keyed.push((Key::Rim { piece: u.edge, rep }, !u.reversed));
            }
            bounds.push(keyed);
        }
        sketches.push(Sketch {
            name: (c.cap_name)(plane, up),
            kind: Kind::Cap { plane, up, plan },
            bounds,
        });
    }
    Ok(Plan {
        data,
        heights,
        sketches,
    })
}

/// `integral of (height + gradient . p) dA` over a plan ring (signed by
/// winding).
///
/// Used only for its SIGN, to tell an outer shell from a cavity, so arcs
/// are sampled rather than integrated in closed form: a shell's volume is
/// never near zero, and the sampling error is a small fraction of it.
fn plane_integral(level: &Level, ring: &ArcRing) -> Scalar {
    const STEPS: usize = 32;
    let count = ring.vertices.len();
    let mut points = Vec::with_capacity(count * STEPS);
    for index in 0..count {
        let vertex = ring.vertices[index];
        let next = ring.vertices[(index + 1) % count].point;
        points.push(vertex.point);
        if vertex.bulge == 0.0 {
            continue;
        }
        if let Ok(arc) = arc_geometry(vertex.point, next, vertex.bulge) {
            let radial = vertex.point - arc.centre;
            let start = radial.y.atan2(radial.x);
            for step in 1..STEPS {
                let t = start + arc.sweep * step as Scalar / STEPS as Scalar;
                points.push(arc.centre + Vec2::new(t.cos(), t.sin()) * arc.radius);
            }
        }
    }
    let (mut area, mut mx, mut my) = (0.0, 0.0, 0.0);
    for index in 0..points.len() {
        let p = points[index];
        let q = points[(index + 1) % points.len()];
        let cross = p.x * q.y - q.x * p.y;
        area += cross / 2.0;
        mx += (p.x + q.x) * cross / 6.0;
        my += (p.y + q.y) * cross / 6.0;
    }
    level.height * area + level.gradient.x * mx + level.gradient.y * my
}

/// Build the solids a column description encloses.
///
/// One `ExactBRep` per connected piece of solid. A result with an enclosed
/// cavity is refused (see below). Faces are grouped into shells by the edges they share;
/// an edge shared by more than two faces means two pieces touch along it,
/// which is not a manifold solid and is refused.
pub(crate) fn build_columns(c: &Columns<'_>) -> GeomResult<Vec<ExactBRep>> {
    let plan = plan(c)?;
    let mut users: BTreeMap<Key, Vec<usize>> = BTreeMap::new();
    for (face, sketch) in plan.sketches.iter().enumerate() {
        for uses in &sketch.bounds {
            for (key, _) in uses {
                users.entry(*key).or_default().push(face);
            }
        }
    }
    let mut sets = Sets::new(plan.sketches.len());
    for faces in users.values() {
        match faces.len() {
            2 => sets.join(faces[0], faces[1]),
            1 => return Err(contract("a column face edge has no neighbour")),
            _ => {
                return Err(unsupported(
                    "column solids that touch along an edge (not a manifold)",
                ))
            }
        }
    }
    let mut shells: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for face in 0..plan.sketches.len() {
        shells.entry(sets.find(face)).or_default().push(face);
    }

    // Outer shells enclose positive volume; a cavity's shell faces into
    // the cavity, so its volume comes out negative. Walls are vertical,
    // so only caps contribute (divergence theorem with the field (0,0,z)).
    let mut outers = Vec::new();
    let mut voids = Vec::new();
    for faces in shells.into_values() {
        let mut volume = 0.0;
        for &face in &faces {
            if let Kind::Cap {
                plane,
                up,
                plan: rings,
            } = &plan.sketches[face].kind
            {
                let sign = if *up { 1.0 } else { -1.0 };
                for uses in rings {
                    let ring = c.arrangement.ring(uses);
                    volume += sign * plane_integral(&c.planes[*plane], &ring);
                }
            }
        }
        if volume > 0.0 {
            outers.push(faces);
        } else if volume < 0.0 {
            voids.push(faces);
        } else {
            return Err(contract("a column shell encloses no volume"));
        }
    }
    // An enclosed cavity would be a void shell, and this kernel reads void
    // shells as boolean intent, not geometry: the mesh compiler tessellates
    // only the outer shell (`void_shells_do_not_add_surface`). Returning one
    // would lose the cavity silently downstream, so it is refused by name.
    if !voids.is_empty() {
        return Err(unsupported(
            "exact prism boolean leaving an enclosed cavity (void shells are not geometry here)",
        ));
    }

    let mut solids = Vec::with_capacity(outers.len());
    for outer in outers {
        let mut emit = Emit::new(c, &plan);
        let outer = emit.shell(&outer)?;
        solids.push(emit.finish(outer)?);
    }
    Ok(solids)
}

/// Builds one solid's faces from sketches, creating each vertex and edge
/// the first time a face uses it.
struct Emit<'a> {
    c: &'a Columns<'a>,
    plan: &'a Plan,
    builder: ExactBRepBuilder,
    vertices: BTreeMap<(usize, usize), VertexId>,
    edges: BTreeMap<Key, EdgeId>,
}

impl<'a> Emit<'a> {
    fn new(c: &'a Columns<'a>, plan: &'a Plan) -> Self {
        Self {
            c,
            plan,
            builder: ExactBRepBuilder::default(),
            vertices: BTreeMap::new(),
            edges: BTreeMap::new(),
        }
    }

    fn point(&self, vertex: usize, k: usize) -> Point3 {
        let p = self.c.arrangement.vertices()[vertex];
        Point3::new(p.x, p.y, self.plan.heights.z(vertex, k))
    }

    fn vertex(&mut self, vertex: usize, k: usize) -> VertexId {
        if let Some(&id) = self.vertices.get(&(vertex, k)) {
            return id;
        }
        let position = self.point(vertex, k);
        let id = self.builder.topology_mut().add_vertex(Vertex { position });
        self.vertices.insert((vertex, k), id);
        id
    }

    /// Start and end corners of an edge, as `(vertex, cluster)`.
    fn ends(&self, key: Key) -> GeomResult<((usize, usize), (usize, usize))> {
        match key {
            Key::Rim { piece, rep } => {
                let edge = &self.c.arrangement.edges()[piece];
                let at = |vertex: usize| -> GeomResult<(usize, usize)> {
                    let p = self.c.arrangement.vertices()[vertex];
                    let k = self.plan.heights.index(vertex, self.c.planes[rep].at(p))?;
                    Ok((vertex, k))
                };
                Ok((at(edge.from)?, at(edge.to)?))
            }
            Key::Vert { vertex, k } => Ok(((vertex, k), (vertex, k + 1))),
        }
    }

    /// The 3D support of an edge and its native span, start to end.
    fn curve(&self, key: Key) -> GeomResult<(Curve3, Interval)> {
        let ((v0, k0), (v1, k1)) = self.ends(key)?;
        let line = |a: Point3, b: Point3| {
            (
                Curve3::Line(Line3 {
                    origin: a,
                    direction: b - a,
                }),
                Interval::UNIT,
            )
        };
        match key {
            Key::Vert { .. } => Ok(line(self.point(v0, k0), self.point(v1, k1))),
            Key::Rim { piece, rep } => {
                let geo = &self.plan.data[piece].geo;
                if geo.bulge == 0.0 {
                    return Ok(line(self.point(v0, k0), self.point(v1, k1)));
                }
                let arc = arc_geometry(geo.from, geo.to, geo.bulge)?;
                let level = self.c.planes[rep];
                if level.gradient == Vec2::ZERO {
                    let start = Point3::new(geo.from.x, geo.from.y, level.height);
                    Ok((
                        Curve3::Circle(circle_of(&arc, level.height, start)?),
                        Interval::new(0.0, arc.sweep),
                    ))
                } else {
                    sloped_arc_edge(&arc, geo.from, level, 0.0)
                }
            }
        }
    }

    fn edge(&mut self, key: Key) -> GeomResult<EdgeId> {
        if let Some(&id) = self.edges.get(&key) {
            return Ok(id);
        }
        let ((v0, k0), (v1, k1)) = self.ends(key)?;
        if (v0, k0) == (v1, k1) {
            return Err(contract("a column edge starts and ends at one corner"));
        }
        let (curve, interval) = self.curve(key)?;
        let start = self.vertex(v0, k0);
        let end = self.vertex(v1, k1);
        let curve = self.builder.add_curve3(curve);
        let id = self.builder.topology_mut().add_edge(Edge {
            start,
            end,
            curve: Some(curve),
        });
        self.builder.set_edge_interval(id, interval);
        self.edges.insert(key, id);
        Ok(id)
    }

    fn shell(&mut self, faces: &[usize]) -> GeomResult<ShellId> {
        let mut ids = Vec::with_capacity(faces.len());
        for &face in faces {
            ids.push((self.face(&self.plan.sketches[face])?, Orientation::Forward));
        }
        Ok(self.builder.topology_mut().add_shell(Shell {
            faces: ids,
            closed: true,
        }))
    }

    fn finish(mut self, outer: ShellId) -> GeomResult<ExactBRep> {
        self.builder.topology_mut().add_solid(Solid {
            outer,
            voids: Vec::new(),
        });
        let exact = self
            .builder
            .finish()
            .map_err(|error| contract(format!("column assembly failed: {error}")))?;
        let health = audit_brep(exact.topology());
        if !health.is_closed_manifold() {
            return Err(contract(format!(
                "column solid is not a closed manifold: {health:?}"
            )));
        }
        Ok(exact)
    }
}

/// The surface a face lies on, and how to express an edge in it.
enum Support {
    /// A horizontal cap: parameters are plan coordinates.
    Flat,
    /// A sloped cap: parameters are coordinates in its frame.
    Sloped(Frame3),
    /// A straight wall: origin, along-piece axis, and up.
    Planar { origin: Point3, x: Vec3 },
    /// An arc wall on a cylinder anchored at the directed piece's start,
    /// so `u` is the angle from that start and `v` the absolute height.
    Round {
        arc: ArcGeometry,
        frame: Frame3,
        /// Arrangement vertex the directed piece starts at.
        start: usize,
    },
}

impl Emit<'_> {
    /// Pcurve of `key` on `support`, in the edge's own direction with the
    /// span that matches its 3D curve's native span proportionally.
    fn pcurve(&self, support: &Support, key: Key) -> GeomResult<(Curve2, Interval)> {
        let ((v0, k0), (v1, k1)) = self.ends(key)?;
        let line = |a: Vec2, b: Vec2| {
            (
                Curve2::Line(Line2 {
                    origin: a,
                    direction: b - a,
                }),
                Interval::UNIT,
            )
        };
        match support {
            Support::Planar { origin, x } => {
                let uv = |p: Point3| {
                    let d = p - *origin;
                    Vec2::new(d.dot(*x), d.z)
                };
                Ok(line(uv(self.point(v0, k0)), uv(self.point(v1, k1))))
            }
            Support::Round { arc, frame, start } => {
                let u = |vertex: usize| if vertex == *start { 0.0 } else { arc.sweep };
                match key {
                    Key::Vert { vertex, .. } => {
                        let a = self.point(v0, k0).z;
                        let b = self.point(v1, k1).z;
                        Ok(line(Vec2::new(u(vertex), a), Vec2::new(u(vertex), b)))
                    }
                    Key::Rim { rep, .. } => {
                        let (u0, u1) = (u(v0), u(v1));
                        let level = self.c.planes[rep];
                        if level.gradient == Vec2::ZERO {
                            Ok(line(
                                Vec2::new(u0, level.height),
                                Vec2::new(u1, level.height),
                            ))
                        } else {
                            let x = Vec2::new(frame.x.x, frame.x.y);
                            let y = Vec2::new(frame.y.x, frame.y.y);
                            Ok((
                                Curve2::Sinusoid(Sinusoid2 {
                                    mean: level.at(arc.centre),
                                    cosine: arc.radius * level.gradient.dot(x),
                                    sine: arc.radius * level.gradient.dot(y),
                                }),
                                Interval::new(u0, u1),
                            ))
                        }
                    }
                }
            }
            Support::Flat => {
                let Key::Rim { piece, .. } = key else {
                    return Err(contract("a cap uses a vertical edge"));
                };
                let geo = &self.plan.data[piece].geo;
                if geo.bulge == 0.0 {
                    return Ok(line(geo.from, geo.to));
                }
                let arc = arc_geometry(geo.from, geo.to, geo.bulge)?;
                Ok((
                    Curve2::Circle(Circle2 {
                        frame: circle2_frame(&arc, geo.from),
                        radius: arc.radius,
                    }),
                    Interval::new(0.0, arc.sweep),
                ))
            }
            Support::Sloped(frame) => {
                let point = |p: Point3| {
                    let d = p - frame.origin;
                    Vec2::new(d.dot(frame.x), d.dot(frame.y))
                };
                let vector = |v: Vec3| Vec2::new(v.dot(frame.x), v.dot(frame.y));
                let (curve, interval) = self.curve(key)?;
                let pcurve = match curve {
                    Curve3::Line(l) => Curve2::Line(Line2 {
                        origin: point(l.origin),
                        direction: vector(l.direction),
                    }),
                    Curve3::Ellipse(e) => Curve2::Ellipse(Ellipse2 {
                        frame: Frame2 {
                            origin: point(e.frame.origin),
                            x: vector(e.frame.x),
                            y: vector(e.frame.y),
                        },
                        semi_axis_x: e.semi_axis_x,
                        semi_axis_y: e.semi_axis_y,
                    }),
                    _ => {
                        return Err(contract(
                            "a sloped cap edge is neither a line nor an ellipse",
                        ))
                    }
                };
                Ok((pcurve, interval))
            }
        }
    }

    fn face(&mut self, sketch: &Sketch) -> GeomResult<FaceId> {
        let (support, surface, orientation) = match &sketch.kind {
            Kind::Wall { piece, left } => {
                let edge = &self.c.arrangement.edges()[*piece];
                let geo = &self.plan.data[*piece].geo;
                // The directed piece has the solid on its left, like a
                // counter-clockwise ring edge in the prism path.
                let (start, from, to, bulge) = if *left {
                    (edge.from, geo.from, geo.to, geo.bulge)
                } else {
                    (edge.to, geo.to, geo.from, -geo.bulge)
                };
                if bulge == 0.0 {
                    let plan = Vec3::new(to.x - from.x, to.y - from.y, 0.0);
                    let length = plan.length();
                    if length == 0.0 || !length.is_finite() {
                        return Err(contract("a column wall stands on a zero-length piece"));
                    }
                    let x = plan / length;
                    let origin = Point3::new(from.x, from.y, 0.0);
                    let surface = Surface::Plane(Plane {
                        frame: Frame3 {
                            origin,
                            x,
                            y: Vec3::Z,
                            z: x.cross(Vec3::Z),
                        },
                    });
                    (Support::Planar { origin, x }, surface, Orientation::Forward)
                } else {
                    let arc = arc_geometry(from, to, bulge)?;
                    let circle = circle_of(&arc, 0.0, Point3::new(from.x, from.y, 0.0))?;
                    let surface = Surface::Cylinder(Cylinder {
                        frame: circle.frame,
                        radius: arc.radius,
                    });
                    (
                        Support::Round {
                            arc,
                            frame: circle.frame,
                            start,
                        },
                        surface,
                        Orientation::Forward,
                    )
                }
            }
            Kind::Cap { plane, up, .. } => {
                let level = self.c.planes[*plane];
                let orientation = if *up {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                if level.gradient == Vec2::ZERO {
                    let frame = identity_frame3(Point3::new(0.0, 0.0, level.height));
                    (Support::Flat, Surface::Plane(Plane { frame }), orientation)
                } else {
                    let frame = level.sloped_frame()?;
                    (
                        Support::Sloped(frame),
                        Surface::Plane(Plane { frame }),
                        orientation,
                    )
                }
            }
        };
        let surface = self.builder.add_surface(surface);
        let mut bounds = Vec::with_capacity(sketch.bounds.len());
        for (index, keyed) in sketch.bounds.iter().enumerate() {
            let mut uses = Vec::with_capacity(keyed.len());
            let mut intervals = Vec::with_capacity(keyed.len());
            for &(key, forward) in keyed {
                let edge = self.edge(key)?;
                let (pcurve, native) = self.pcurve(&support, key)?;
                let pcurve = self.builder.add_curve2(pcurve);
                uses.push(EdgeUse {
                    edge,
                    orientation: if forward {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    },
                    pcurve: Some(pcurve),
                });
                intervals.push(if forward {
                    native
                } else {
                    Interval::new(native.end, native.start)
                });
            }
            bounds.push(FaceBound {
                loop_id: add_loop(&mut self.builder, uses, intervals),
                orientation: Orientation::Forward,
                outer: index == 0,
            });
        }
        let face = self.builder.topology_mut().add_face(Face {
            surface: Some(surface),
            bounds,
            orientation,
        });
        if let Some(name) = &sketch.name {
            self.builder.set_face_name(face, name.clone());
        }
        Ok(face)
    }
}
