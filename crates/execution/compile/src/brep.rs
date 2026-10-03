//! Faceted B-rep tessellation.
//!
//! A brep face is a planar polygon in an arbitrary plane, possibly
//! concave or holed. A triangle fan is wrong for both, so each face is
//! projected to its own plane, triangulated with the same earcut path
//! profiles use, and lifted back. Shared vertices stay shared: the loop
//! indices already reference interned topology vertices.
//!
//! # Deviation (#232)
//!
//! Every face also reports a certified bound on how far its exact surface
//! lies from its triangles, as a [`Deviation`] per surface family:
//!
//! - A planar face with straight edges is its own polygon: the bound is
//!   how far its vertices sit off its declared plane (zero without one). A
//!   curved edge on a planar face is chorded vertex to vertex here, so it
//!   is unbounded by name, as is a face skipped for having no area.
//! - A curved face is bounded per triangle in its parameter domain. With
//!   `L` the linear interpolant of the surface at the triangle's parameter
//!   corners, Taylor's theorem gives `|S(x) - L(x)| <= 1/2 sum_i l_i
//!   q(y_i - x)`, `l_i` the barycentric weights and `q(d) = A d_u^2 + 2 B
//!   |d_u d_v| + C d_v^2 <= (A + B) d_u^2 + (C + B) d_v^2`, where `A, B, C`
//!   bound `|S_uu|, |S_uv|, |S_vv|` over the triangle's parameter box
//!   ([`axiolid_reference::bound::SurfaceBoundOracle`]). In coordinates
//!   scaled by `sqrt(A + B)` and `sqrt(C + B)` the sum is `R^2 - |z - c|^2`
//!   for any centre `c` and the corners' largest distance `R` from it, so
//!   at most the squared radius of the triangle's smallest enclosing
//!   circle there. Where a knot line of full multiplicity crosses the box
//!   the first-order bound `max_edges (|S_u| |d_u| + |S_v| |d_v|)` is used
//!   instead. Each mesh vertex's distance from the surface at its own
//!   parameters (a seam vertex shared from the neighbouring face) is added.
//! - The trim: the exact domain differs from the triangulated polygon by
//!   lenses between each pcurve span and its chord. A lens point is within
//!   the span's certified chord bound `d` of the chord
//!   ([`axiolid_reference::bound::chord_bound2`]), which is an edge of a
//!   triangle, so its surface point is within `G d` of that triangle's
//!   bound, `G = sqrt(|S_u|^2 + |S_v|^2)` over the lens. A gap between
//!   consecutive pcurves adds the same way.
//!
//! A triangle the mesh drops because two of its corners share a vertex (a
//! pole) still covers its parameter triangle; its bound counts only if the
//! segment it collapses to is an edge of the face's emitted triangles.
//! For a rational B-spline the interpolation coefficients come from its
//! homogeneous numerator and weight, not from the quotient rule
//! ([`axiolid_reference::bound::SurfaceBounds::interpolation`]); both sets
//! are tried and the smaller bound kept.
//!
//! The bound also steers the mesh where that is cheap: a triangle whose
//! certified bound misses the budget has its widest free edge split, and
//! a trim edge is sampled until its certified chord bound fits, both up to
//! caps beyond which the bound is reported as it stands.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Scalar, Vec3};
use axiolid_reference::bound::SurfaceBoundOracle;

use crate::deviation::{Deviation, DeviationBound, DeviationPath};
use axiolid_mesh::TriMesh;
use axiolid_mesh_compile_contract::MeshClosure;
use axiolid_model::NodeId;
use axiolid_topology::{BRep, Orientation};
use std::collections::{HashMap, HashSet};

use crate::planar::{earcut_planar_face, earcut_projected, newell_normal, plane_axes};

const MAX_BREP_FACES: usize = 1 << 16;
const MAX_BREP_TOPOLOGY_ITEMS: usize = 1 << 20;
const MAX_BREP_EDGE_USES: usize = 1 << 20;
const MAX_TOTAL_CURVED_RECORDS: usize = 1 << 22;
const MAX_MESH_POSITIONS: usize = 1 << 22;
const MAX_MESH_INDICES: usize = 1 << 24;

fn check_tessellation_input_budget(brep: &BRep<NodeId>) -> GeomResult<()> {
    for (count, resource) in [
        (brep.vertices().len(), "B-rep vertices"),
        (brep.edges().len(), "B-rep edges"),
        (brep.loops().len(), "B-rep loops"),
        (brep.shells().len(), "B-rep shells"),
        (brep.solids().len(), "B-rep solids"),
    ] {
        if count > MAX_BREP_TOPOLOGY_ITEMS {
            return Err(GeomError::BudgetExceeded { resource });
        }
    }
    if brep.faces().len() > MAX_BREP_FACES {
        return Err(GeomError::BudgetExceeded {
            resource: "B-rep faces",
        });
    }
    let edge_uses = brep.loops().iter().try_fold(0_usize, |total, wire| {
        total
            .checked_add(wire.edges.len())
            .ok_or(GeomError::BudgetExceeded {
                resource: "B-rep edge uses",
            })
    })?;
    if edge_uses > MAX_BREP_EDGE_USES {
        return Err(GeomError::BudgetExceeded {
            resource: "B-rep edge uses",
        });
    }
    let face_bounds = brep.faces().iter().try_fold(0_usize, |total, face| {
        total
            .checked_add(face.bounds.len())
            .ok_or(GeomError::BudgetExceeded {
                resource: "B-rep face bounds",
            })
    })?;
    let shell_faces = brep.shells().iter().try_fold(0_usize, |total, shell| {
        total
            .checked_add(shell.faces.len())
            .ok_or(GeomError::BudgetExceeded {
                resource: "B-rep shell face uses",
            })
    })?;
    for (count, resource) in [
        (face_bounds, "B-rep face bounds"),
        (shell_faces, "B-rep shell face uses"),
    ] {
        if count > MAX_BREP_EDGE_USES {
            return Err(GeomError::BudgetExceeded { resource });
        }
    }
    Ok(())
}

fn consume_tessellation_work(total: &mut usize, count: usize) -> GeomResult<()> {
    let next = total.checked_add(count).ok_or(GeomError::BudgetExceeded {
        resource: "expanded B-rep tessellation work",
    })?;
    if next > MAX_BREP_EDGE_USES {
        return Err(GeomError::BudgetExceeded {
            resource: "expanded B-rep tessellation work",
        });
    }
    *total = next;
    Ok(())
}

fn checked_output_len(
    current: usize,
    add: usize,
    limit: usize,
    resource: &'static str,
) -> GeomResult<usize> {
    current
        .checked_add(add)
        .filter(|&next| next <= limit)
        .ok_or(GeomError::BudgetExceeded { resource })
}

fn extend_output<T: Copy>(
    target: &mut Vec<T>,
    values: &[T],
    limit: usize,
    resource: &'static str,
) -> GeomResult<()> {
    checked_output_len(target.len(), values.len(), limit, resource)?;
    target.extend_from_slice(values);
    Ok(())
}

fn check_mesh_budget(mesh: &TriMesh) -> GeomResult<()> {
    if mesh.positions.len() > MAX_MESH_POSITIONS {
        return Err(GeomError::BudgetExceeded {
            resource: "tessellated mesh positions",
        });
    }
    if mesh.indices.len() > MAX_MESH_INDICES {
        return Err(GeomError::BudgetExceeded {
            resource: "tessellated mesh indices",
        });
    }
    Ok(())
}

/// Tessellate one faceted B-rep into a triangle mesh.
///
/// With solids, every solid's outer shell and void shells contribute
/// surface and the result is [`MeshClosure::Solid`]; several disjoint solids
/// mesh as several closed components. A void is a cavity: its
/// triangles face into it, away from the material, so the mesh's signed
/// volume is the outer volume less every cavity. A void can only remove
/// material, so one authored facing the other way (the STEP convention
/// stores a void shell outward and reverses it on use) is turned round
/// rather than added. A void shell that is not closed is refused: it would
/// leave the mesh open.
///
/// Without a solid -- a surface model, e.g. IFC `IfcShellBasedSurfaceModel`
/// -- every shell is tessellated as authored and the result is
/// [`MeshClosure::Surface`], open or closed alike (#161). Promoting a shell
/// to a solid here would claim a volume the source never declared, so the
/// flag travels with the mesh and volume readers refuse it.
pub fn tessellate(
    brep: &BRep<NodeId>,
    graph: &axiolid_model::GeometryGraph,
    tolerance: axiolid_core::Tolerance,
    chord_error: Scalar,
) -> GeomResult<(TriMesh, MeshClosure, Deviation)> {
    check_tessellation_input_budget(brep)?;
    // Structure before geometry. A dangling handle or an open loop
    // produces a mesh that looks plausible and is wrong, so the
    // topology is audited before a single triangle is emitted.
    let health = axiolid_topology::audit_brep(brep);
    if !health.is_tessellable() {
        return Err(GeomError::InvalidInput(format!(
            "brep topology is not tessellable: {health:?}"
        )));
    }

    // Each shell with whether it is a void. Every solid contributes: a
    // B-rep of several disjoint solids (a composite profile whose members
    // do not touch) meshes all of them.
    let (shells, closure): (Vec<(&axiolid_topology::Shell, bool)>, MeshClosure) =
        if brep.solids().is_empty() {
            if brep.shells().is_empty() {
                return Err(GeomError::InvalidInput(
                    "brep has neither a solid nor a shell".to_string(),
                ));
            }
            (
                brep.shells().iter().map(|shell| (shell, false)).collect(),
                MeshClosure::Surface,
            )
        } else {
            let mut shells = Vec::new();
            for solid in brep.solids() {
                let shell = brep
                    .shells()
                    .get(solid.outer.index())
                    .ok_or_else(|| GeomError::InvalidInput("outer shell missing".to_string()))?;
                shells.push((shell, false));
                for void in &solid.voids {
                    let shell = brep
                        .shells()
                        .get(void.index())
                        .ok_or_else(|| GeomError::InvalidInput("void shell missing".to_string()))?;
                    check_void_closed(brep, shell)?;
                    shells.push((shell, true));
                }
            }
            (shells, MeshClosure::Solid)
        };
    if shells.iter().any(|(shell, _)| shell.faces.is_empty()) {
        return Err(GeomError::InvalidInput(
            "a tessellated shell has no faces".to_string(),
        ));
    }
    let shell_faces = || shells.iter().flat_map(|(shell, _)| shell.faces.iter());

    let mut expanded_work = 0_usize;
    for &(face_id, _) in shell_faces() {
        consume_tessellation_work(&mut expanded_work, 1)?;
        let face = brep
            .faces()
            .get(face_id.index())
            .ok_or_else(|| GeomError::InvalidInput("face missing".to_owned()))?;
        consume_tessellation_work(&mut expanded_work, face.bounds.len())?;
        for bound in &face.bounds {
            let wire = brep
                .loops()
                .get(bound.loop_id.index())
                .ok_or_else(|| GeomError::InvalidInput("loop missing".to_owned()))?;
            consume_tessellation_work(&mut expanded_work, wire.edges.len())?;
        }
    }

    let mut mesh = TriMesh::default();
    let ctx = FaceContext {
        brep,
        graph,
        tolerance,
        chord_error,
        closure,
    };
    let mut edge_cache: EdgeSamples = EdgeSamples::new();
    let mut welded: std::collections::HashMap<axiolid_topology::VertexId, u32> =
        std::collections::HashMap::new();
    let mut total_curved_records = 0_usize;
    let mut deviation = Deviation::empty();
    for &(shell, void) in &shells {
        let first_index = mesh.indices.len();
        for &(face_id, shell_sense) in &shell.faces {
            let face = brep
                .faces()
                .get(face_id.index())
                .ok_or_else(|| GeomError::InvalidInput("face missing".to_string()))?;
            let flip = (shell_sense == Orientation::Reversed)
                ^ (face.orientation == Orientation::Reversed);
            append_face(
                &mut mesh,
                &ctx,
                face,
                flip,
                &mut welded,
                &mut edge_cache,
                &mut total_curved_records,
                &mut deviation,
            )?;
            check_mesh_budget(&mesh)?;
        }
        if void {
            face_into_cavity(&mut mesh, first_index);
        }
    }
    Ok((mesh, closure, deviation))
}

/// A void shell must be closed on its own: every edge its faces use is
/// used exactly twice within it.
fn check_void_closed(brep: &BRep<NodeId>, shell: &axiolid_topology::Shell) -> GeomResult<()> {
    let mut uses: std::collections::HashMap<axiolid_topology::EdgeId, usize> =
        std::collections::HashMap::new();
    for &(face_id, _) in &shell.faces {
        let face = brep
            .faces()
            .get(face_id.index())
            .ok_or_else(|| GeomError::InvalidInput("face missing".to_owned()))?;
        for bound in &face.bounds {
            let wire = brep
                .loops()
                .get(bound.loop_id.index())
                .ok_or_else(|| GeomError::InvalidInput("loop missing".to_owned()))?;
            for use_ in &wire.edges {
                *uses.entry(use_.edge).or_default() += 1;
            }
        }
    }
    if shell.faces.is_empty() || uses.values().any(|&count| count != 2) {
        return Err(GeomError::InvalidInput(
            "a void shell is not closed: a cavity must be bounded on every side".to_owned(),
        ));
    }
    Ok(())
}

/// Turn the triangles from `first_index` on so they enclose negative
/// volume: a cavity removes material whichever way it was authored.
fn face_into_cavity(mesh: &mut TriMesh, first_index: usize) {
    let positions = &mesh.positions;
    let volume: Scalar = mesh.indices[first_index..]
        .chunks_exact(3)
        .map(|t| {
            let (a, b, c) = (
                positions[t[0] as usize],
                positions[t[1] as usize],
                positions[t[2] as usize],
            );
            a.dot(b.cross(c))
        })
        .sum();
    if volume > 0.0 {
        for triangle in mesh.indices[first_index..].chunks_exact_mut(3) {
            triangle.swap(1, 2);
        }
    }
}

/// Triangulate one face and append it to the mesh.
/// Everything a face needs that does not vary within one shell.
struct FaceContext<'a> {
    brep: &'a BRep<NodeId>,
    graph: &'a axiolid_model::GeometryGraph,
    tolerance: axiolid_core::Tolerance,
    /// How far a chord may sit from a curved face or edge (#165).
    chord_error: Scalar,
    /// What the tessellated shells are: [`MeshClosure::Surface`] for a
    /// surface model, where a planar face that encloses no area covers
    /// nothing and is skipped (#171). In a solid it stays refused.
    closure: MeshClosure,
}

#[allow(clippy::too_many_arguments)]
fn append_face(
    mesh: &mut TriMesh,
    ctx: &FaceContext<'_>,
    face: &axiolid_topology::Face<NodeId>,
    flip: bool,
    welded: &mut std::collections::HashMap<axiolid_topology::VertexId, u32>,
    cache: &mut EdgeSamples,
    total_curved_records: &mut usize,
    deviation: &mut Deviation,
) -> GeomResult<()> {
    let (brep, graph) = (ctx.brep, ctx.graph);
    // A curved support cannot be tessellated by projecting the boundary onto
    // a plane: the interior curves away from it and the error is invisible in
    // the output. Sample the surface itself when the face states its boundary
    // in surface parameters, and refuse when it does not: a missing wall is
    // cheap, a wrong wall corrupts every downstream quantity.
    let support = face_surface(graph, face)?;
    if let Some(surface) = support {
        if !surface_is_planar(surface) {
            let mut bound = FaceBound::new(surface);
            with_curved_face_transaction(
                mesh,
                welded,
                cache,
                total_curved_records,
                ctx.chord_error,
                |state| append_curved_face(state, ctx, face, surface, flip, &mut bound),
            )?;
            deviation.add(
                DeviationPath::BRepFace,
                surface_family(surface),
                bound.finish(),
            );
            return Ok(());
        }
    }
    deviation.add(
        DeviationPath::BRepFace,
        "plane",
        planar_face_bound(ctx, face, support)?,
    );
    let mut rings: Vec<Vec<(axiolid_topology::VertexId, Vec3)>> = Vec::new();
    let mut outer_index = None;
    for bound in &face.bounds {
        let points = loop_points(brep, bound)?;
        if points.len() < 3 {
            return Err(GeomError::Degenerate(
                "planar face bound has fewer than three vertices".into(),
            ));
        }
        if bound.outer && outer_index.is_none() {
            outer_index = Some(rings.len());
        }
        rings.push(points);
    }
    if rings.is_empty() {
        return Err(GeomError::Degenerate(
            "planar face has no non-degenerate bounds".into(),
        ));
    }
    let outer_index = outer_index.ok_or_else(|| {
        GeomError::InvalidInput("planar face has no outer bound after topology audit".into())
    })?;
    rings.swap(0, outer_index);

    for (index, ring) in rings.iter().enumerate() {
        if plane_axes(newell_normal(ring.iter().map(|&(_, p)| p))).is_none() {
            // A surface model's outer bound that encloses no area adds no
            // surface: skipping it removes nothing the file drew. Real
            // exports (Nova pipe fittings) write each end cap as a strip of
            // bowtie quads through the pipe axis, whose signed area cancels
            // exactly, and every edge of such a strip is still used twice
            // by the fitting's real faces. A solid is not skipped: a face
            // it declares to bound a volume must enclose area, and a hole
            // with no area is refused in both, since skipping it would
            // leave real area of its face uncut.
            let finite = ring.iter().all(|&(_, p)| p.is_finite());
            if index == 0 && finite && ctx.closure == MeshClosure::Surface {
                // Its winding may still enclose lobes of area (a bowtie),
                // which the mesh then lacks.
                deviation.add(
                    DeviationPath::BRepFace,
                    "plane",
                    DeviationBound::Unbounded("planar face with no net area, skipped"),
                );
                return Ok(());
            }
            return Err(GeomError::Degenerate(
                "planar face bound has zero or non-finite area".into(),
            ));
        }
    }
    let normal = newell_normal(rings[0].iter().map(|&(_, p)| p));
    let (u, v) = plane_axes(normal).ok_or_else(|| {
        GeomError::Degenerate("planar face outer bound has no stable plane".into())
    })?;
    let origin = rings[0][0].1;

    let mut flat: Vec<[Scalar; 2]> = Vec::new();
    let mut positions: Vec<(axiolid_topology::VertexId, Vec3)> = Vec::new();
    let mut hole_starts: Vec<usize> = Vec::new();
    for (index, ring) in rings.iter().enumerate() {
        if index > 0 {
            hole_starts.push(flat.len());
        }
        for &(vertex, point) in ring {
            let d = point - origin;
            flat.push([d.dot(u), d.dot(v)]);
            positions.push((vertex, point));
        }
    }

    let indices = earcut_planar_face(&flat, &hole_starts, ctx.tolerance.linear());
    if indices.is_empty() || indices.len() % 3 != 0 {
        return Err(GeomError::Degenerate(format!(
            "face triangulation produced {} indices for {} vertices",
            indices.len(),
            flat.len()
        )));
    }
    checked_output_len(
        mesh.indices.len(),
        indices.len(),
        MAX_MESH_INDICES,
        "mesh indices",
    )?;

    // Weld by topological vertex. Adjacent facets already share interned
    // vertices upstream; emitting per-face copies would leave every edge
    // unshared, so the mesh would look correct yet fail a manifold check.
    let new_vertices: HashSet<_> = positions
        .iter()
        .map(|(vertex, _)| *vertex)
        .filter(|vertex| !welded.contains_key(vertex))
        .collect();
    checked_output_len(
        mesh.positions.len(),
        new_vertices.len(),
        MAX_MESH_POSITIONS,
        "mesh positions",
    )?;
    let mut local: Vec<u32> = Vec::with_capacity(positions.len());
    for &(vertex, point) in &positions {
        let index = *welded.entry(vertex).or_insert_with(|| {
            let next = mesh.positions.len() as u32;
            mesh.positions.push(point);
            next
        });
        local.push(index);
    }

    for triangle in indices.chunks_exact(3) {
        let (a, b, c) = (local[triangle[0]], local[triangle[1]], local[triangle[2]]);
        if flip {
            mesh.indices.extend([a, c, b]);
        } else {
            mesh.indices.extend([a, b, c]);
        }
    }
    Ok(())
}

/// Walk a bound loop and collect its vertex positions in order.
fn loop_points(
    brep: &BRep<NodeId>,
    bound: &axiolid_topology::FaceBound,
) -> GeomResult<Vec<(axiolid_topology::VertexId, Vec3)>> {
    let wire = brep
        .loops()
        .get(bound.loop_id.index())
        .ok_or_else(|| GeomError::InvalidInput("loop missing".to_string()))?;
    let mut points: Vec<(axiolid_topology::VertexId, Vec3)> = Vec::with_capacity(wire.edges.len());
    for use_ in &wire.edges {
        let edge = brep
            .edges()
            .get(use_.edge.index())
            .ok_or_else(|| GeomError::InvalidInput("edge missing".to_string()))?;
        // Each edge contributes its start under traversal orientation; the
        // loop is closed, so the final end repeats the first start.
        let vertex = if use_.orientation == Orientation::Forward {
            edge.start
        } else {
            edge.end
        };
        let position = brep
            .vertices()
            .get(vertex.index())
            .ok_or_else(|| GeomError::InvalidInput("vertex missing".to_string()))?
            .position;
        points.push((vertex, position));
    }
    if bound.orientation == Orientation::Reversed {
        points.reverse();
    }
    Ok(points)
}

/// Whether an exact support surface is planar.
///
/// A planar support adds nothing the boundary polygon does not already say,
/// so those faces stay on the fast path. Every other surface family
/// curves away from the boundary plane, and projecting it there is a
/// silent modelling error.
fn surface_is_planar(surface: &axiolid_surface::Surface) -> bool {
    matches!(surface, axiolid_surface::Surface::Plane(_))
}

/// Resolve a face's support surface from the graph, if it declares one.
///
/// A face may legitimately omit its surface: a planar polygon's boundary
/// already determines its plane. A declared handle that does not resolve to
/// a Surface node is a graph error, not an absent surface.
fn face_surface<'g>(
    graph: &'g axiolid_model::GeometryGraph,
    face: &axiolid_topology::Face<NodeId>,
) -> GeomResult<Option<&'g axiolid_surface::Surface>> {
    let Some(id) = face.surface else {
        return Ok(None);
    };
    match graph.get(id) {
        Some(axiolid_model::GeometryNode::Surface(s)) => Ok(Some(s)),
        Some(_) => Err(GeomError::InvalidInput(format!(
            "face support {id:?} is not a Surface node"
        ))),
        None => Err(GeomError::InvalidInput(format!(
            "face support {id:?} does not belong to this graph"
        ))),
    }
}

fn seed_grid_vertices(
    state: &mut CurvedMeshState<'_>,
    boundary: &CurvedBoundary,
    patch: &GridPatch,
    surface: &axiolid_surface::Surface,
) -> GeomResult<Option<(Vec<SurfaceVertex>, u32)>> {
    let columns = patch.nu.checked_add(1).ok_or(GeomError::BudgetExceeded {
        resource: "curved-face grid columns",
    })?;
    let rows = patch.nv.checked_add(1).ok_or(GeomError::BudgetExceeded {
        resource: "curved-face grid rows",
    })?;
    let total = columns
        .checked_mul(rows)
        .filter(|&count| count <= MAX_CURVED_FACE_VERTICES)
        .ok_or(GeomError::BudgetExceeded {
            resource: "curved-face grid vertices",
        })?;
    // Index boundary points by their grid cell so existing vertices win.
    let mut lookup: HashMap<(usize, usize), SurfaceVertex> =
        HashMap::with_capacity(boundary.uv.len());
    for (position, p) in boundary.uv.iter().enumerate() {
        let iu = ((p.x - patch.u_start) / patch.du).round();
        let iv = ((p.y - patch.v_start) / patch.dv).round();
        if iu < 0.0 || iv < 0.0 {
            continue;
        }
        let (iu, iv) = (iu as usize, iv as usize);
        if iu <= patch.nu && iv <= patch.nv {
            let current = SurfaceVertex {
                uv: *p,
                mesh: boundary.shared[position],
                local: u32::try_from(position).map_err(|_| GeomError::BudgetExceeded {
                    resource: "curved-face local vertices",
                })?,
            };
            if let Some(previous) = lookup.insert((iu, iv), current) {
                if previous.mesh != current.mesh {
                    return Ok(None);
                }
            }
        }
    }
    let interior = total.saturating_sub(lookup.len());
    state.reserve_face_vertices(interior)?;
    let mut next_local =
        u32::try_from(boundary.uv.len()).map_err(|_| GeomError::BudgetExceeded {
            resource: "curved-face local vertices",
        })?;
    let mut grid = Vec::with_capacity(total);
    for iv in 0..rows {
        for iu in 0..columns {
            if let Some(&existing) = lookup.get(&(iu, iv)) {
                grid.push(existing);
                continue;
            }
            let u = patch.u_start + patch.du * iu as Scalar;
            let v = patch.v_start + patch.dv * iv as Scalar;
            let uv = Point2::new(u, v);
            let point = axiolid_reference::surface::evaluate(surface, u, v)?;
            let mesh = state.push_face_position(point)?;
            let local = take_local_vertex(&mut next_local)?;
            grid.push(SurfaceVertex { uv, mesh, local });
        }
    }
    Ok(Some((grid, next_local)))
}

/// Tessellate a face whose support surface is curved.
///
/// The face's boundary must be stated in the surface's parameter domain:
/// every edge use needs a pcurve. Without them the trim is unknown, and
/// guessing it by inverting the surface is not generally solvable. A
/// missing pcurve is therefore refused, not approximated.
///
/// The parameter-space boundary gives a `(u, v)` polygon. Its bounding box
/// is the patch actually sampled, so a half-cylinder costs half a cylinder
/// rather than a full revolution clipped afterwards.
// Shared boundary samples, keyed canonically in the topological edge's
// start-to-end direction. Interior refinement is deliberately face-local.
type EdgeSamples = HashMap<axiolid_topology::EdgeId, (Vec<u32>, Orientation)>;
type WeldedVertices = HashMap<axiolid_topology::VertexId, u32>;
type LocalEdgeKey = (u32, u32);

/// Mutable state for one curved face plus the shell-wide boundary caches.
struct CurvedMeshState<'a> {
    mesh: &'a mut TriMesh,
    welded: &'a mut WeldedVertices,
    edge_samples: &'a mut EdgeSamples,
    chord_error: Scalar,
    /// Boundary occurrences plus refinement vertices for this face. This is
    /// deliberately independent of shell vertices created by earlier faces.
    face_local_vertices: usize,
    total_curved_records: &'a mut usize,
    new_welded: Vec<axiolid_topology::VertexId>,
    new_edge_samples: Vec<axiolid_topology::EdgeId>,
}

fn with_curved_face_transaction(
    mesh: &mut TriMesh,
    welded: &mut WeldedVertices,
    edge_samples: &mut EdgeSamples,
    total_curved_records: &mut usize,
    chord_error: Scalar,
    operation: impl FnOnce(&mut CurvedMeshState<'_>) -> GeomResult<()>,
) -> GeomResult<()> {
    let positions_start = mesh.positions.len();
    let indices_start = mesh.indices.len();
    let records_start = *total_curved_records;
    let mut state = CurvedMeshState {
        mesh,
        welded,
        edge_samples,
        chord_error,
        face_local_vertices: 0,
        total_curved_records,
        new_welded: Vec::new(),
        new_edge_samples: Vec::new(),
    };
    let result = operation(&mut state);
    if result.is_err() {
        state.mesh.positions.truncate(positions_start);
        state.mesh.indices.truncate(indices_start);
        for vertex in &state.new_welded {
            state.welded.remove(vertex);
        }
        for edge in &state.new_edge_samples {
            state.edge_samples.remove(edge);
        }
        *state.total_curved_records = records_start;
    }
    result
}

impl CurvedMeshState<'_> {
    fn reserve_face_vertices(&mut self, count: usize) -> GeomResult<()> {
        let next =
            self.face_local_vertices
                .checked_add(count)
                .ok_or(GeomError::BudgetExceeded {
                    resource: "curved-face local vertices",
                })?;
        if next > MAX_CURVED_FACE_VERTICES {
            return Err(GeomError::BudgetExceeded {
                resource: "curved-face local vertices",
            });
        }
        let total =
            self.total_curved_records
                .checked_add(count)
                .ok_or(GeomError::BudgetExceeded {
                    resource: "total curved-face records",
                })?;
        if total > MAX_TOTAL_CURVED_RECORDS {
            return Err(GeomError::BudgetExceeded {
                resource: "total curved-face records",
            });
        }
        self.face_local_vertices = next;
        *self.total_curved_records = total;
        Ok(())
    }

    fn push_face_position(&mut self, point: Vec3) -> GeomResult<u32> {
        if !point.is_finite() {
            return Err(GeomError::Degenerate(
                "curved-face position is non-finite".to_owned(),
            ));
        }
        let next =
            u32::try_from(self.mesh.positions.len()).map_err(|_| GeomError::BudgetExceeded {
                resource: "mesh vertex indices",
            })?;
        extend_output(
            &mut self.mesh.positions,
            &[point],
            MAX_MESH_POSITIONS,
            "mesh positions",
        )?;
        Ok(next)
    }
}

#[derive(Debug, Clone, Copy)]
struct SurfaceVertex {
    uv: Point2,
    mesh: u32,
    /// Face-local identity. Unlike `mesh`, periodic seam occurrences remain
    /// distinct while adjacent Earcut triangles share the same occurrence.
    local: u32,
}

#[derive(Debug, Clone, Copy)]
struct SurfaceTriangle {
    vertices: [SurfaceVertex; 3],
    /// Boundary provenance for edges AB, BC, and CA. Pcurve-sampled trim
    /// segments are immutable during support-surface interior refinement.
    boundary: [bool; 3],
    depth: u8,
}

const MAX_CURVED_EDGE_SEGMENTS: usize = 4_096;
const MAX_REFINEMENT_DEPTH: u8 = 20;
const MAX_CURVED_FACE_VERTICES: usize = 1 << 18;

/// Samples an edge needs so its 3D image meets the chord budget.
///
/// The trim is flattened in PARAMETER space, but tolerance is a claim about
/// 3D. A cylinder rim is the straight pcurve v = 0,
/// u in [0, TAU]: two points in parameter space, a full circle in space.
/// Subdivide until the mapped midpoint is within tolerance of the chord.
fn edge_sample_count(
    curve: &axiolid_curve::Curve2,
    surface: &axiolid_surface::Surface,
    tolerance: Scalar,
) -> GeomResult<usize> {
    if !(tolerance.is_finite() && tolerance > 0.0) {
        return Err(GeomError::InvalidInput(
            "curved-edge chord tolerance must be positive and finite".to_owned(),
        ));
    }
    let domain = axiolid_reference::curve::domain2(curve);
    let at = |t: Scalar| -> GeomResult<Vec3> {
        let p = axiolid_reference::curve::evaluate2(curve, t)?;
        axiolid_reference::surface::evaluate(surface, p.x, p.y)
    };
    // Powers of two keep the count stable: a face arriving later computes
    // the same value from the same inputs.
    let mut n = 1usize;
    loop {
        let mut worst: Scalar = 0.0;
        for i in 0..n {
            let t0 = domain.start + (domain.end - domain.start) * (i as Scalar) / (n as Scalar);
            let t1 =
                domain.start + (domain.end - domain.start) * ((i + 1) as Scalar) / (n as Scalar);
            let (a, b) = (at(t0)?, at(t1)?);
            let m = at(0.5 * (t0 + t1))?;
            worst = worst.max((m - 0.5 * (a + b)).length());
        }
        if worst <= tolerance {
            return Ok(certified_edge_count(curve, surface, tolerance, n));
        }
        if n >= MAX_CURVED_EDGE_SEGMENTS {
            return Err(GeomError::BudgetExceeded {
                resource: "curved edge samples",
            });
        }
        n *= 2;
    }
}

/// `n`, doubled until every span's certified 3D chord bound is within
/// `tolerance` too, where that is cheap (#232); `n` itself when the
/// bound cannot be had or would need more than the segment cap.
///
/// A span's image `S(g(t))` is within `G d + q(D) / 8` of its chord: `d`
/// the pcurve's certified chord bound in parameters and `G` the surface's
/// Lipschitz bound, plus the surface's interpolation error along the
/// parameter segment `D` between the span's ends (the one-dimensional case
/// of the triangle bound in the module notes). A boundary edge is never
/// split by interior refinement, so a coarse one would leave its
/// triangles' bounds wide.
fn certified_edge_count(
    curve: &axiolid_curve::Curve2,
    surface: &axiolid_surface::Surface,
    tolerance: Scalar,
    n: usize,
) -> usize {
    let Some(oracle) = SurfaceBoundOracle::new(surface) else {
        return n;
    };
    let domain = axiolid_reference::curve::domain2(curve);
    let fits = |count: usize| -> Option<bool> {
        let at = |i: usize| {
            domain.start + (domain.end - domain.start) * (i as Scalar) / (count as Scalar)
        };
        for i in 0..count {
            let (t0, t1) = (at(i), at(i + 1));
            let d = axiolid_reference::bound::chord_bound2(curve, t0, t1)?;
            let a = axiolid_reference::curve::evaluate2(curve, t0).ok()?;
            let b = axiolid_reference::curve::evaluate2(curve, t1).ok()?;
            let grow = Point2::splat(d);
            let (lo, hi) = (a.min(b) - grow, a.max(b) + grow);
            let bounds = oracle.bounds((lo.x, hi.x), (lo.y, hi.y))?;
            let span = b - a;
            let [x, y, z] = bounds.interpolation;
            let along = 0.125
                * (x * span.x * span.x + 2.0 * y * (span.x * span.y).abs() + z * span.y * span.y);
            let along = if bounds.smooth {
                along
            } else {
                bounds.du * span.x.abs() + bounds.dv * span.y.abs()
            };
            let lipschitz = (bounds.du * bounds.du + bounds.dv * bounds.dv).sqrt();
            if lipschitz * d + along > tolerance {
                return Some(false);
            }
        }
        Some(true)
    };
    let mut count = n;
    while count <= MAX_CURVED_EDGE_SEGMENTS {
        match fits(count) {
            Some(true) => return count,
            Some(false) => count *= 2,
            None => return n,
        }
    }
    n
}

/// Sample a trim as `n` segments, including both endpoints.
///
/// The complete sequence is cached so a reversed edge use can reverse it
/// without shifting the pcurve-to-vertex pairing. Callers omit the final
/// point when appending the edge to a loop; the next edge contributes that
/// junction.
fn trim_samples(curve: &axiolid_curve::Curve2, n: usize) -> GeomResult<Vec<axiolid_core::Point2>> {
    let domain = axiolid_reference::curve::domain2(curve);
    let mut out = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let t = domain.start + (domain.end - domain.start) * (i as Scalar) / (n as Scalar);
        out.push(axiolid_reference::curve::evaluate2(curve, t)?);
    }
    Ok(out)
}

/// Boundary of a curved face in parameter space, paired with the shared 3D
/// vertex index for each point.
///
/// The pairing is what makes a seam watertight: (u, v) drives triangulation,
/// and the index says which already-created vertex to emit.
struct CurvedBoundary {
    uv: Vec<axiolid_core::Point2>,
    shared: Vec<u32>,
    hole_starts: Vec<usize>,
    winding_reversed: bool,
}

type CurvedRing = (usize, usize, bool, bool);

fn validate_curved_rings(rings: &[CurvedRing]) -> GeomResult<usize> {
    if rings
        .iter()
        .any(|&(start, end, _, _)| end.checked_sub(start).is_none_or(|len| len < 3))
    {
        return Err(GeomError::Degenerate(
            "curved bound has fewer than three samples".to_owned(),
        ));
    }
    let mut outers = rings.iter().enumerate().filter(|(_, ring)| ring.2);
    let Some((index, _)) = outers.next() else {
        return Err(GeomError::Degenerate(
            "curved face must declare one outer bound".to_owned(),
        ));
    };
    if outers.next().is_some() {
        return Err(GeomError::Degenerate(
            "curved face declares multiple outer bounds".to_owned(),
        ));
    }
    Ok(index)
}

/// A face boundary recognised as a rectangle in parameter space.
///
/// Cylinders, cones, spheres and tori are trimmed by constant-u and
/// constant-v curves in the overwhelmingly common case, and such a patch
/// has a structure a polygon triangulator cannot see. Earcut works in UV,
/// where one unit of u and one unit of v are interchangeable; in 3D they
/// are not. On a cylinder of radius r, a step in u is r times longer than
/// the same step in v, so a triangulation that is perfectly reasonable in
/// parameter space can join a point on the bottom rim to a distant point
/// on the top rim, and that chord cuts THROUGH the tube rather than
/// running along it. The mesh stays closed and its area grows.
///
/// Recognising the rectangle lets the patch be meshed as a grid, where
/// every quad spans exactly one cell and no chord can cut across.
struct GridPatch {
    /// Cells along u and v.
    nu: usize,
    nv: usize,
    u_start: Scalar,
    v_start: Scalar,
    du: Scalar,
    dv: Scalar,
}

/// Recognise a boundary as a rectangular patch, or decline.
///
/// Declining is not a failure: a trimmed face with a hole, a slanted trim
/// or an irregular sample layout is genuinely not a grid, and earcut
/// remains the right tool for it. This only claims the cases it can prove.
fn recognise_grid(boundary: &CurvedBoundary) -> Option<GridPatch> {
    // A hole means the patch is not simply a rectangle.
    if !boundary.hole_starts.is_empty() || boundary.uv.len() < 4 {
        return None;
    }
    let (mut u_min, mut u_max) = (Scalar::INFINITY, Scalar::NEG_INFINITY);
    let (mut v_min, mut v_max) = (Scalar::INFINITY, Scalar::NEG_INFINITY);
    for p in &boundary.uv {
        if !p.x.is_finite() || !p.y.is_finite() {
            return None;
        }
        u_min = u_min.min(p.x);
        u_max = u_max.max(p.x);
        v_min = v_min.min(p.y);
        v_max = v_max.max(p.y);
    }
    let (u_span, v_span) = (u_max - u_min, v_max - v_min);
    if !(u_span > 0.0 && v_span > 0.0) {
        return None;
    }
    // Every boundary point must lie ON the rectangle's border, and the
    // spacing must be uniform. Counting distinct coordinates recovers the
    // cell counts without assuming which side the walk started on.
    let on_edge = |value: Scalar, lo: Scalar, hi: Scalar, span: Scalar| {
        (value - lo).abs() <= span * 1e-9 || (value - hi).abs() <= span * 1e-9
    };
    let mut us: Vec<Scalar> = Vec::new();
    let mut vs: Vec<Scalar> = Vec::new();
    for p in &boundary.uv {
        let u_on = on_edge(p.x, u_min, u_max, u_span);
        let v_on = on_edge(p.y, v_min, v_max, v_span);
        if !u_on && !v_on {
            return None;
        }
        if v_on {
            us.push(p.x);
        }
        if u_on {
            vs.push(p.y);
        }
    }
    let nu = distinct_steps(&mut us, u_min, u_span)?;
    let nv = distinct_steps(&mut vs, v_min, v_span)?;
    if nu == 0 || nv == 0 {
        return None;
    }
    let patch = GridPatch {
        nu,
        nv,
        u_start: u_min,
        v_start: v_min,
        du: u_span / nu as Scalar,
        dv: v_span / nv as Scalar,
    };
    if !ordered_grid_perimeter(boundary, &patch) {
        return None;
    }
    Some(patch)
}

fn grid_lattice_index(value: Scalar, start: Scalar, step: Scalar, cells: usize) -> Option<usize> {
    let raw = (value - start) / step;
    if !raw.is_finite() {
        return None;
    }
    let index = raw.round();
    if index < 0.0 || index > cells as Scalar {
        return None;
    }
    Some(index as usize)
}

fn ordered_grid_perimeter(boundary: &CurvedBoundary, patch: &GridPatch) -> bool {
    let Some(expected) = patch
        .nu
        .checked_mul(2)
        .and_then(|u| patch.nv.checked_mul(2).and_then(|v| u.checked_add(v)))
    else {
        return false;
    };
    if boundary.uv.len() != expected {
        return false;
    }
    let mut keys = Vec::with_capacity(expected);
    for point in &boundary.uv {
        let Some(i) = grid_lattice_index(point.x, patch.u_start, patch.du, patch.nu) else {
            return false;
        };
        let Some(j) = grid_lattice_index(point.y, patch.v_start, patch.dv, patch.nv) else {
            return false;
        };
        if i != 0 && i != patch.nu && j != 0 && j != patch.nv {
            return false;
        }
        keys.push((i, j));
    }
    if keys.iter().copied().collect::<HashSet<_>>().len() != expected {
        return false;
    }
    keys.iter()
        .zip(keys.iter().cycle().skip(1))
        .all(|(&(a, b), &(c, d))| a.abs_diff(c) + b.abs_diff(d) == 1)
}

/// Count uniform cells spanned by a set of coordinates, or decline.
///
/// Non-uniform spacing means the sides were sampled at different rates and
/// a grid would not line up with the boundary, so the caller must fall
/// back rather than weld mismatched vertices.
fn distinct_steps(values: &mut Vec<Scalar>, start: Scalar, span: Scalar) -> Option<usize> {
    values.sort_by(Scalar::total_cmp);
    values.dedup_by(|a, b| (*a - *b).abs() <= span * 1e-9);
    // At least two distinct values survive dedup whenever the span is
    // positive, which recognise_grid has already established, so cells
    // is at least one and the division below is safe.
    let cells = values.len().checked_sub(1)?;
    let step = span / cells as Scalar;
    for (index, value) in values.iter().enumerate() {
        let want = start + step * index as Scalar;
        if (value - want).abs() > span * 1e-6 {
            return None;
        }
    }
    Some(cells)
}

/// Collect a curved face's boundary, sampling each edge exactly once across
/// the whole shell.
fn curved_boundary(
    state: &mut CurvedMeshState<'_>,
    ctx: &FaceContext<'_>,
    face: &axiolid_topology::Face<NodeId>,
    surface: &axiolid_surface::Surface,
    bound: &mut FaceBound<'_>,
) -> GeomResult<CurvedBoundary> {
    let (brep, graph) = (ctx.brep, ctx.graph);
    let mut out = CurvedBoundary {
        uv: Vec::new(),
        shared: Vec::new(),
        hole_starts: Vec::new(),
        winding_reversed: false,
    };
    let mut rings: Vec<CurvedRing> = Vec::with_capacity(face.bounds.len());
    for face_bound in &face.bounds {
        let ring_start = out.uv.len();
        let wire = brep
            .loops()
            .get(face_bound.loop_id.index())
            .ok_or_else(|| GeomError::InvalidInput("loop missing".to_string()))?;
        // The previous pcurve's end and the ring's first start, for the
        // gaps between consecutive pcurves (#232).
        let mut previous_end: Option<Point2> = None;
        let mut ring_first: Option<Point2> = None;
        for use_ in &wire.edges {
            let Some(pcurve) = use_.pcurve else {
                return Err(GeomError::Unsupported {
                    backend: crate::BACKEND_ID,
                    operation: axiolid_contracts::Operation::Tessellation,
                });
            };
            let Some(axiolid_model::GeometryNode::Curve2(trim)) = graph.get(pcurve) else {
                return Err(GeomError::InvalidInput(
                    "edge pcurve must reference a Curve2 node".to_string(),
                ));
            };
            let n = match state.edge_samples.get(&use_.edge) {
                Some((vertices, _)) => vertices.len().checked_sub(1).ok_or_else(|| {
                    GeomError::Degenerate("cached edge has no segments".to_owned())
                })?,
                None => edge_sample_count(trim, surface, ctx.chord_error)?,
            };
            state.reserve_face_vertices(n)?;
            let params = trim_samples(trim, n)?;
            bound.trim(trim, n, &params);
            if let (Some(end), Some(&start)) = (previous_end, params.first()) {
                bound.gap(end, start, surface);
            }
            ring_first = ring_first.or_else(|| params.first().copied());
            previous_end = params.last().copied();
            // Evaluate the trim on the surface: these are the seam's
            // 3D points, and they are interned under the edge.
            let points: Vec<Vec3> = params
                .iter()
                .map(|p| axiolid_reference::surface::evaluate(surface, p.x, p.y))
                .collect::<GeomResult<_>>()?;
            let edge = brep
                .edges()
                .get(use_.edge.index())
                .ok_or_else(|| GeomError::InvalidInput("edge missing".to_string()))?;
            let (start_vertex, end_vertex) = if use_.orientation == Orientation::Forward {
                (edge.start, edge.end)
            } else {
                (edge.end, edge.start)
            };
            let shared = edge_samples(
                state,
                use_.edge,
                use_.orientation,
                start_vertex,
                end_vertex,
                &points,
            )?;
            // The next edge contributes this edge's final junction.
            out.uv.extend(params.into_iter().take(n));
            out.shared.extend(shared.into_iter().take(n));
        }
        if let (Some(end), Some(start)) = (previous_end, ring_first) {
            bound.gap(end, start, surface);
        }
        if face_bound.orientation == Orientation::Reversed {
            out.uv[ring_start..].reverse();
            out.shared[ring_start..].reverse();
        }
        rings.push((
            ring_start,
            out.uv.len(),
            face_bound.outer,
            face_bound.orientation == Orientation::Reversed,
        ));
    }
    let outer_index = validate_curved_rings(&rings)?;
    rings.swap(0, outer_index);
    out.winding_reversed = rings[0].3;
    let raw_uv = core::mem::take(&mut out.uv);
    let raw_shared = core::mem::take(&mut out.shared);
    for (index, &(start, end, _, _)) in rings.iter().enumerate() {
        let outer_anchor = if index == 0 {
            None
        } else {
            Some(parameter_chart_center(&out.uv).ok_or_else(|| {
                GeomError::Degenerate("curved outer trim has no parameter center".to_owned())
            })?)
        };
        if index > 0 {
            out.hole_starts.push(out.uv.len());
        }
        let ring_start = out.uv.len();
        out.uv.extend_from_slice(&raw_uv[start..end]);
        out.shared.extend_from_slice(&raw_shared[start..end]);
        unwrap_parameter_ring(surface, &mut out.uv[ring_start..], outer_anchor)?;
    }
    Ok(out)
}

/// Intern one edge's 3D samples, or return the existing ones.
///
/// Interning is keyed by EDGE, not by face. Whichever face reaches a seam
/// first creates the vertices; every later face reuses the identical
/// indices, so the seam is shared rather than merely coincident.
fn edge_samples(
    state: &mut CurvedMeshState<'_>,
    edge: axiolid_topology::EdgeId,
    sense: Orientation,
    start_vertex: axiolid_topology::VertexId,
    end_vertex: axiolid_topology::VertexId,
    points: &[Vec3],
) -> GeomResult<Vec<u32>> {
    if let Some((existing, stored)) = state.edge_samples.get(&edge) {
        let mut out = existing.clone();
        // Reverse the COMPLETE endpoint-inclusive sequence. Reversing an
        // endpoint-excluding sequence shifts every sample by one segment.
        if *stored != sense {
            out.reverse();
        }
        if out.len() != points.len() {
            return Err(GeomError::Degenerate(
                "shared edge pcurves produced different sample counts".to_owned(),
            ));
        }
        for (&index, point) in out.iter().zip(points) {
            let disagreement = (state.mesh.positions[index as usize] - *point).length();
            if !disagreement.is_finite() || disagreement > state.chord_error {
                return Err(GeomError::Degenerate(
                    "shared edge pcurves disagree in 3D".to_owned(),
                ));
            }
        }
        return Ok(out);
    }
    if points.len() < 2 {
        return Err(GeomError::Degenerate(
            "curved edge needs both endpoint samples".to_owned(),
        ));
    }

    let mut indices: Vec<u32> = Vec::with_capacity(points.len());
    for (i, point) in points.iter().enumerate() {
        let index = if i == 0 {
            welded_vertex(state, start_vertex, *point)?
        } else if i + 1 == points.len() {
            welded_vertex(state, end_vertex, *point)?
        } else {
            state.push_face_position(*point)?
        };
        indices.push(index);
    }

    state.edge_samples.insert(edge, (indices.clone(), sense));
    state.new_edge_samples.push(edge);
    Ok(indices)
}

fn welded_vertex(
    state: &mut CurvedMeshState<'_>,
    vertex: axiolid_topology::VertexId,
    point: Vec3,
) -> GeomResult<u32> {
    if let Some(&existing) = state.welded.get(&vertex) {
        let disagreement = (state.mesh.positions[existing as usize] - point).length();
        if !disagreement.is_finite() || disagreement > state.chord_error {
            return Err(GeomError::Degenerate(
                "pcurve endpoint disagrees with the welded topology vertex".to_owned(),
            ));
        }
        Ok(existing)
    } else {
        let next = state.push_face_position(point)?;
        state.welded.insert(vertex, next);
        state.new_welded.push(vertex);
        Ok(next)
    }
}

/// Tessellate a face whose support surface is curved.
///
/// The boundary comes from shared edge samples, so a seam between two curved
/// faces uses one set of vertices. Interior detail still comes from the
/// surface: the boundary alone would flatten the patch.
fn append_curved_face(
    state: &mut CurvedMeshState<'_>,
    ctx: &FaceContext<'_>,
    face: &axiolid_topology::Face<NodeId>,
    surface: &axiolid_surface::Surface,
    flip: bool,
    bound: &mut FaceBound<'_>,
) -> GeomResult<()> {
    let boundary = curved_boundary(state, ctx, face, surface, bound)?;
    if boundary.uv.len() < 3 {
        return Err(GeomError::Degenerate(
            "curved face boundary is underspecified".to_owned(),
        ));
    }
    if let Some(patch) = recognise_grid(&boundary) {
        let boundary_edges = local_boundary_edges(boundary.uv.len(), &boundary.hole_starts)?;
        if let Some((grid, next_local)) = seed_grid_vertices(state, &boundary, &patch, surface)? {
            let columns = patch.nu + 1;
            let capacity = patch.nu.saturating_mul(patch.nv).saturating_mul(2);
            let mut triangles = Vec::with_capacity(capacity);
            for iv in 0..patch.nv {
                for iu in 0..patch.nu {
                    let a = grid[iv * columns + iu];
                    let b = grid[iv * columns + iu + 1];
                    let c = grid[(iv + 1) * columns + iu + 1];
                    let d = grid[(iv + 1) * columns + iu];
                    triangles.push(surface_triangle([a, b, d], &boundary_edges));
                    triangles.push(surface_triangle([b, c, d], &boundary_edges));
                }
            }
            return refine_curved_face(
                state,
                surface,
                triangles,
                next_local,
                flip ^ boundary.winding_reversed,
                bound,
            );
        }
    }
    let flat: Vec<[Scalar; 2]> = boundary.uv.iter().map(|p| [p.x, p.y]).collect();
    let indices = earcut_projected(&flat, &boundary.hole_starts);
    if indices.is_empty() || indices.len() % 3 != 0 {
        return Err(GeomError::Degenerate(format!(
            "curved face trim triangulation produced {} indices for {} points",
            indices.len(),
            flat.len()
        )));
    }
    let ring_edges = local_boundary_edges(boundary.uv.len(), &boundary.hole_starts)?;
    let indices = restore_collinear_samples(&flat, indices);
    // Flipped in a metric where a unit of u and of v are about as long on
    // the surface: a cylinder's u is its radius times longer than its v.
    let metric = parameter_metric(surface, &boundary.uv);
    let scaled: Vec<[Scalar; 2]> = flat
        .iter()
        .map(|p| [p[0] * metric.0, p[1] * metric.1])
        .collect();
    let indices = delaunay_flips(&scaled, indices, &ring_edges);

    let boundary_edges = local_boundary_edges(boundary.uv.len(), &boundary.hole_starts)?;
    let mut triangles = Vec::with_capacity(indices.len() / 3);
    for t in indices.chunks_exact(3) {
        let vertex = |index: usize| -> GeomResult<SurfaceVertex> {
            Ok(SurfaceVertex {
                uv: boundary.uv[index],
                mesh: boundary.shared[index],
                local: u32::try_from(index).map_err(|_| GeomError::BudgetExceeded {
                    resource: "curved-face local vertices",
                })?,
            })
        };
        let vertices = [vertex(t[0])?, vertex(t[1])?, vertex(t[2])?];
        triangles.push(SurfaceTriangle {
            vertices,
            boundary: [
                boundary_edges.contains(&local_edge_key(vertices[0], vertices[1])),
                boundary_edges.contains(&local_edge_key(vertices[1], vertices[2])),
                boundary_edges.contains(&local_edge_key(vertices[2], vertices[0])),
            ],
            depth: 0,
        });
    }

    let next_local = u32::try_from(boundary.uv.len()).map_err(|_| GeomError::BudgetExceeded {
        resource: "curved-face local vertices",
    })?;
    // Earcut is trusted to cover the trim polygon; a triangulation whose
    // parameter area misses the polygon's would leave part of the face
    // unmeshed, so the bound is withdrawn rather than claimed.
    if !covers_polygon(&boundary, &triangles) {
        bound.unbounded("curved face triangulation does not cover its trim");
    }
    refine_curved_face(
        state,
        surface,
        triangles,
        next_local,
        flip ^ boundary.winding_reversed,
        bound,
    )
}

fn surface_triangle(
    vertices: [SurfaceVertex; 3],
    boundary_edges: &HashSet<LocalEdgeKey>,
) -> SurfaceTriangle {
    SurfaceTriangle {
        vertices,
        boundary: [
            boundary_edges.contains(&local_edge_key(vertices[0], vertices[1])),
            boundary_edges.contains(&local_edge_key(vertices[1], vertices[2])),
            boundary_edges.contains(&local_edge_key(vertices[2], vertices[0])),
        ],
        depth: 0,
    }
}

fn local_boundary_edges(
    vertex_count: usize,
    hole_starts: &[usize],
) -> GeomResult<HashSet<LocalEdgeKey>> {
    let mut ring_starts = Vec::with_capacity(hole_starts.len() + 2);
    ring_starts.push(0);
    ring_starts.extend_from_slice(hole_starts);
    ring_starts.push(vertex_count);
    let mut edges = HashSet::with_capacity(vertex_count);
    for range in ring_starts.windows(2) {
        let (start, end) = (range[0], range[1]);
        if end.saturating_sub(start) < 3 {
            return Err(GeomError::Degenerate(
                "curved trim ring has fewer than three samples".to_owned(),
            ));
        }
        for index in start..end {
            let next = if index + 1 == end { start } else { index + 1 };
            let a = u32::try_from(index).map_err(|_| GeomError::BudgetExceeded {
                resource: "curved-face local vertices",
            })?;
            let b = u32::try_from(next).map_err(|_| GeomError::BudgetExceeded {
                resource: "curved-face local vertices",
            })?;
            edges.insert(sorted_local_edge(a, b));
        }
    }
    Ok(edges)
}

/// Refine a curved face, also splitting where its certified bound misses
/// the budget (#232) while that stays within the refinement limits.
///
/// The certified splits are best effort, never a new refusal: if the pass
/// that makes them runs out of depth or vertices, everything it added is
/// discarded and the face is refined on the midpoint and centroid checks
/// alone, as before, its bound then reported as it is.
fn refine_curved_face(
    state: &mut CurvedMeshState<'_>,
    surface: &axiolid_surface::Surface,
    triangles: Vec<SurfaceTriangle>,
    next_local: u32,
    flip: bool,
    bound: &mut FaceBound<'_>,
) -> GeomResult<()> {
    let positions = state.mesh.positions.len();
    let indices = state.mesh.indices.len();
    let local = state.face_local_vertices;
    let records = *state.total_curved_records;
    match refine_pass(
        state,
        surface,
        triangles.clone(),
        next_local,
        flip,
        bound,
        true,
    ) {
        Err(GeomError::BudgetExceeded { resource }) => {
            let _ = resource;
            state.mesh.positions.truncate(positions);
            state.mesh.indices.truncate(indices);
            state.face_local_vertices = local;
            *state.total_curved_records = records;
            refine_pass(state, surface, triangles, next_local, flip, bound, false)
        }
        other => other,
    }
}

#[allow(clippy::too_many_arguments)]
fn refine_pass(
    state: &mut CurvedMeshState<'_>,
    surface: &axiolid_surface::Surface,
    mut triangles: Vec<SurfaceTriangle>,
    mut next_local: u32,
    flip: bool,
    bound: &mut FaceBound<'_>,
    certified: bool,
) -> GeomResult<()> {
    let mut midpoint_cache: HashMap<LocalEdgeKey, SurfaceVertex> = HashMap::new();

    loop {
        let mut requested_edges = HashSet::new();
        let mut self_masks = Vec::with_capacity(triangles.len());
        let mut centroid_failures = Vec::with_capacity(triangles.len());
        let mut pending = Vec::with_capacity(triangles.len());
        for triangle in triangles {
            let (edge_errors, centroid_error) =
                curved_triangle_errors(state.mesh, surface, triangle.vertices)?;
            let mut self_mask = 0_u8;
            for (edge, &edge_error) in edge_errors.iter().enumerate() {
                if !triangle.boundary[edge] && edge_error > state.chord_error {
                    let (a, b) = triangle_edge(triangle.vertices, edge);
                    requested_edges.insert(local_edge_key(a, b));
                    self_mask |= 1 << edge;
                }
            }
            // The midpoint and centroid checks above are measurements. Where
            // the certified bound still misses the budget, split the longest
            // edge that may be split, while it is cheap: shallow, within the
            // face's vertex budget, and with the bound not dominated by a
            // corner's own offset, which splitting cannot reduce (#232).
            if certified
                && self_mask == 0
                && triangle.depth < CERTIFIED_DEPTH
                && state.face_local_vertices < CERTIFIED_FACE_VERTICES
            {
                if let Some(certified) =
                    bound.triangle_bound(state.mesh, surface, triangle.vertices)?
                {
                    if certified > state.chord_error {
                        if let Some(edge) = bound.widest_free_edge(&triangle) {
                            let (a, b) = triangle_edge(triangle.vertices, edge);
                            requested_edges.insert(local_edge_key(a, b));
                            self_mask |= 1 << edge;
                        }
                    }
                }
            }
            self_masks.push(self_mask);
            centroid_failures.push(centroid_error > state.chord_error);
            pending.push(triangle);
        }

        let mut changed = false;
        let mut next = Vec::with_capacity(pending.len().saturating_mul(2));
        for ((triangle, centroid_fails), self_mask) in
            pending.into_iter().zip(centroid_failures).zip(self_masks)
        {
            let mut mask = 0_u8;
            for edge in 0..3 {
                let (a, b) = triangle_edge(triangle.vertices, edge);
                if requested_edges.contains(&local_edge_key(a, b)) {
                    mask |= 1 << edge;
                }
            }
            debug_assert_eq!(mask & self_mask, self_mask);
            if mask == 0 && !centroid_fails {
                next.push(triangle);
                continue;
            }
            let limit = if certified {
                MAX_CERTIFIED_REFINEMENT_DEPTH
            } else {
                MAX_REFINEMENT_DEPTH
            };
            if triangle.depth >= limit {
                return Err(GeomError::BudgetExceeded {
                    resource: "curved-face refinement depth",
                });
            }
            changed = true;
            if mask == 0 {
                split_at_centroid(state, surface, triangle, &mut next_local, &mut next)?;
            } else {
                split_requested_edges(
                    state,
                    surface,
                    triangle,
                    mask,
                    &mut midpoint_cache,
                    &mut next_local,
                    &mut next,
                )?;
            }
        }

        if !changed {
            let triangle_count = next
                .iter()
                .filter(|triangle| {
                    let [a, b, c] = triangle.vertices.map(|vertex| vertex.mesh);
                    a != b && b != c && c != a
                })
                .count();
            let index_count = triangle_count
                .checked_mul(3)
                .ok_or(GeomError::BudgetExceeded {
                    resource: "mesh indices",
                })?;
            checked_output_len(
                state.mesh.indices.len(),
                index_count,
                MAX_MESH_INDICES,
                "mesh indices",
            )?;
            bound.triangles(state.mesh, surface, &next)?;
            for triangle in next {
                let [a, b, c] = triangle.vertices.map(|vertex| vertex.mesh);
                if a == b || b == c || c == a {
                    continue;
                }
                if flip {
                    state.mesh.indices.extend([a, c, b]);
                } else {
                    state.mesh.indices.extend([a, b, c]);
                }
            }
            return Ok(());
        }
        triangles = next;
    }
}

fn split_at_centroid(
    state: &mut CurvedMeshState<'_>,
    surface: &axiolid_surface::Surface,
    triangle: SurfaceTriangle,
    next_local: &mut u32,
    out: &mut Vec<SurfaceTriangle>,
) -> GeomResult<()> {
    let [a, b, c] = triangle.vertices;
    let uv = (a.uv + b.uv + c.uv) / 3.0;
    state.reserve_face_vertices(1)?;
    let exact = axiolid_reference::surface::evaluate(surface, uv.x, uv.y)?;
    let center = SurfaceVertex {
        uv,
        mesh: state.push_face_position(exact)?,
        local: take_local_vertex(next_local)?,
    };
    let depth = triangle.depth + 1;
    out.extend([
        SurfaceTriangle {
            vertices: [a, b, center],
            boundary: [triangle.boundary[0], false, false],
            depth,
        },
        SurfaceTriangle {
            vertices: [b, c, center],
            boundary: [triangle.boundary[1], false, false],
            depth,
        },
        SurfaceTriangle {
            vertices: [c, a, center],
            boundary: [triangle.boundary[2], false, false],
            depth,
        },
    ]);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn split_requested_edges(
    state: &mut CurvedMeshState<'_>,
    surface: &axiolid_surface::Surface,
    triangle: SurfaceTriangle,
    mask: u8,
    cache: &mut HashMap<LocalEdgeKey, SurfaceVertex>,
    next_local: &mut u32,
    out: &mut Vec<SurfaceTriangle>,
) -> GeomResult<()> {
    let [a, b, c] = triangle.vertices;
    let [e0, e1, e2] = triangle.boundary;
    let m0 = if mask & 1 != 0 {
        Some(refinement_midpoint(
            state, surface, a, b, cache, next_local,
        )?)
    } else {
        None
    };
    let m1 = if mask & 2 != 0 {
        Some(refinement_midpoint(
            state, surface, b, c, cache, next_local,
        )?)
    } else {
        None
    };
    let m2 = if mask & 4 != 0 {
        Some(refinement_midpoint(
            state, surface, c, a, cache, next_local,
        )?)
    } else {
        None
    };
    let depth = triangle.depth + 1;
    let triangle = |vertices, boundary| SurfaceTriangle {
        vertices,
        boundary,
        depth,
    };
    match mask {
        1 => {
            let m0 = m0.unwrap();
            out.extend([
                triangle([a, m0, c], [e0, false, e2]),
                triangle([m0, b, c], [e0, e1, false]),
            ]);
        }
        2 => {
            let m1 = m1.unwrap();
            out.extend([
                triangle([b, m1, a], [e1, false, e0]),
                triangle([m1, c, a], [e1, e2, false]),
            ]);
        }
        4 => {
            let m2 = m2.unwrap();
            out.extend([
                triangle([c, m2, b], [e2, false, e1]),
                triangle([m2, a, b], [e2, e0, false]),
            ]);
        }
        3 => {
            let (m0, m1) = (m0.unwrap(), m1.unwrap());
            out.extend([
                triangle([m0, b, m1], [e0, e1, false]),
                triangle([a, m0, m1], [e0, false, false]),
                triangle([a, m1, c], [false, e1, e2]),
            ]);
        }
        6 => {
            let (m1, m2) = (m1.unwrap(), m2.unwrap());
            out.extend([
                triangle([m1, c, m2], [e1, e2, false]),
                triangle([b, m1, m2], [e1, false, false]),
                triangle([b, m2, a], [false, e2, e0]),
            ]);
        }
        5 => {
            let (m2, m0) = (m2.unwrap(), m0.unwrap());
            out.extend([
                triangle([m2, a, m0], [e2, e0, false]),
                triangle([c, m2, m0], [e2, false, false]),
                triangle([c, m0, b], [false, e0, e1]),
            ]);
        }
        7 => {
            let (m0, m1, m2) = (m0.unwrap(), m1.unwrap(), m2.unwrap());
            out.extend([
                triangle([a, m0, m2], [e0, false, e2]),
                triangle([m0, b, m1], [e0, e1, false]),
                triangle([m2, m1, c], [false, e1, e2]),
                triangle([m0, m1, m2], [false; 3]),
            ]);
        }
        _ => unreachable!("non-empty three-edge mask"),
    }
    Ok(())
}

fn triangle_edge(vertices: [SurfaceVertex; 3], edge: usize) -> (SurfaceVertex, SurfaceVertex) {
    match edge {
        0 => (vertices[0], vertices[1]),
        1 => (vertices[1], vertices[2]),
        2 => (vertices[2], vertices[0]),
        _ => unreachable!("triangle edge index"),
    }
}

fn local_edge_key(a: SurfaceVertex, b: SurfaceVertex) -> LocalEdgeKey {
    sorted_local_edge(a.local, b.local)
}

fn sorted_local_edge(a: u32, b: u32) -> LocalEdgeKey {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

fn take_local_vertex(next: &mut u32) -> GeomResult<u32> {
    let current = *next;
    *next = next.checked_add(1).ok_or(GeomError::BudgetExceeded {
        resource: "curved-face local vertices",
    })?;
    Ok(current)
}

fn surface_periods(surface: &axiolid_surface::Surface) -> (Option<Scalar>, Option<Scalar>) {
    let tau = core::f64::consts::TAU;
    match surface {
        axiolid_surface::Surface::Cylinder(_)
        | axiolid_surface::Surface::EllipticalCylinder(_)
        | axiolid_surface::Surface::Cone(_)
        | axiolid_surface::Surface::Sphere(_) => (Some(tau), None),
        axiolid_surface::Surface::Torus(_) => (Some(tau), Some(tau)),
        _ => (None, None),
    }
}

fn parameter_chart_center(points: &[Point2]) -> Option<Point2> {
    let first = *points.first()?;
    let mut min = first;
    let mut max = first;
    for &point in &points[1..] {
        min.x = min.x.min(point.x);
        min.y = min.y.min(point.y);
        max.x = max.x.max(point.x);
        max.y = max.y.max(point.y);
    }
    Some((min + max) * 0.5)
}

fn unwrap_parameter_ring(
    surface: &axiolid_surface::Surface,
    ring: &mut [Point2],
    outer_anchor: Option<Point2>,
) -> GeomResult<()> {
    let Some(first) = ring.first().copied() else {
        return Err(GeomError::Degenerate(
            "curved trim ring is empty".to_owned(),
        ));
    };
    let (u_period, v_period) = surface_periods(surface);
    if let Some(anchor) = outer_anchor {
        let aligned = Point2::new(
            localized_parameter(anchor.x, first.x, u_period),
            localized_parameter(anchor.y, first.y, v_period),
        );
        let shift = aligned - first;
        for point in ring.iter_mut() {
            *point += shift;
        }
    }
    for index in 1..ring.len() {
        ring[index] = Point2::new(
            localized_parameter(ring[index - 1].x, ring[index].x, u_period),
            localized_parameter(ring[index - 1].y, ring[index].y, v_period),
        );
        if !(ring[index].x.is_finite() && ring[index].y.is_finite()) {
            return Err(GeomError::InvalidInput(
                "periodic trim chart produced a non-finite parameter".to_owned(),
            ));
        }
    }
    Ok(())
}

fn localized_parameter(base: Scalar, value: Scalar, period: Option<Scalar>) -> Scalar {
    let Some(period) = period else {
        return value;
    };
    let mut offset = (value - base + period * 0.5).rem_euclid(period) - period * 0.5;
    // Exactly antipodal parameters have two equally short arcs. Choose the
    // same local chart after reversing an edge.
    if (offset + period * 0.5).abs() <= Scalar::EPSILON * period * 4.0 && value > base {
        offset = period * 0.5;
    }
    base + offset
}

fn curved_triangle_errors(
    mesh: &TriMesh,
    surface: &axiolid_surface::Surface,
    vertices: [SurfaceVertex; 3],
) -> GeomResult<([Scalar; 3], Scalar)> {
    let mut edge_errors = [0.0; 3];
    for (edge, edge_error) in edge_errors.iter_mut().enumerate() {
        let (a, b) = triangle_edge(vertices, edge);
        let uv = (a.uv + b.uv) * 0.5;
        let exact = axiolid_reference::surface::evaluate(surface, uv.x, uv.y)?;
        let chord = (mesh.positions[a.mesh as usize] + mesh.positions[b.mesh as usize]) * 0.5;
        *edge_error = (exact - chord).length();
        if !edge_error.is_finite() {
            return Err(GeomError::Degenerate(
                "curved-face edge error is non-finite".to_owned(),
            ));
        }
    }

    let uv = (vertices[0].uv + vertices[1].uv + vertices[2].uv) / 3.0;
    let exact = axiolid_reference::surface::evaluate(surface, uv.x, uv.y)?;
    let linear = (mesh.positions[vertices[0].mesh as usize]
        + mesh.positions[vertices[1].mesh as usize]
        + mesh.positions[vertices[2].mesh as usize])
        / 3.0;
    let centroid_error = (exact - linear).length();
    if !centroid_error.is_finite() {
        return Err(GeomError::Degenerate(
            "curved-face centroid error is non-finite".to_owned(),
        ));
    }
    Ok((edge_errors, centroid_error))
}

fn refinement_midpoint(
    state: &mut CurvedMeshState<'_>,
    surface: &axiolid_surface::Surface,
    a: SurfaceVertex,
    b: SurfaceVertex,
    cache: &mut HashMap<LocalEdgeKey, SurfaceVertex>,
    next_local: &mut u32,
) -> GeomResult<SurfaceVertex> {
    let uv = (a.uv + b.uv) * 0.5;
    let exact = axiolid_reference::surface::evaluate(surface, uv.x, uv.y)?;
    let key = local_edge_key(a, b);
    if let Some(existing) = cache.get(&key).copied() {
        let distance = (state.mesh.positions[existing.mesh as usize] - exact).length();
        if distance > state.chord_error {
            return Err(GeomError::Degenerate(format!(
                "face-local refined edge {key:?} disagrees by {distance}"
            )));
        }
        return Ok(SurfaceVertex {
            uv,
            mesh: existing.mesh,
            local: existing.local,
        });
    }

    state.reserve_face_vertices(1)?;
    let vertex = SurfaceVertex {
        uv,
        mesh: state.push_face_position(exact)?,
        local: take_local_vertex(next_local)?,
    };
    cache.insert(key, vertex);
    Ok(vertex)
}

/// The surface family's name, for a report.
fn surface_family(surface: &axiolid_surface::Surface) -> &'static str {
    match surface {
        axiolid_surface::Surface::Plane(_) => "plane",
        axiolid_surface::Surface::Cylinder(_) => "cylinder",
        axiolid_surface::Surface::EllipticalCylinder(_) => "elliptical cylinder",
        axiolid_surface::Surface::Cone(_) => "cone",
        axiolid_surface::Surface::Sphere(_) => "sphere",
        axiolid_surface::Surface::Torus(_) => "torus",
        axiolid_surface::Surface::BSpline(_) => "B-spline surface",
        _ => "surface",
    }
}

/// A planar face's bound: its vertices' distance off the declared plane,
/// when every edge is straight (#232).
fn planar_face_bound(
    ctx: &FaceContext<'_>,
    face: &axiolid_topology::Face<NodeId>,
    support: Option<&axiolid_surface::Surface>,
) -> GeomResult<DeviationBound> {
    let mut off_plane: Scalar = 0.0;
    let plane = match support {
        Some(axiolid_surface::Surface::Plane(plane)) => {
            let normal = plane.frame.x.cross(plane.frame.y).normalize_or_zero();
            (normal != Vec3::ZERO).then_some((plane.frame.origin, normal))
        }
        _ => None,
    };
    for bound in &face.bounds {
        let wire = ctx
            .brep
            .loops()
            .get(bound.loop_id.index())
            .ok_or_else(|| GeomError::InvalidInput("loop missing".to_owned()))?;
        for use_ in &wire.edges {
            let edge = ctx
                .brep
                .edges()
                .get(use_.edge.index())
                .ok_or_else(|| GeomError::InvalidInput("edge missing".to_owned()))?;
            let straight = match edge.curve {
                None => true,
                Some(id) => {
                    matches!(
                        ctx.graph.get(id),
                        Some(axiolid_model::GeometryNode::Curve3(
                            axiolid_curve::Curve3::Line(_)
                        ))
                    ) || matches!(
                        ctx.graph.get(id),
                        Some(axiolid_model::GeometryNode::Curve3(axiolid_curve::Curve3::Polyline(p)))
                            if p.points.len() == 2 && !p.closed
                    )
                }
            };
            if !straight {
                return Ok(DeviationBound::Unbounded(
                    "planar face with a curved edge, chorded vertex to vertex",
                ));
            }
            if let Some((origin, normal)) = plane {
                for vertex in [edge.start, edge.end] {
                    let position = ctx
                        .brep
                        .vertices()
                        .get(vertex.index())
                        .ok_or_else(|| GeomError::InvalidInput("vertex missing".to_owned()))?
                        .position;
                    off_plane = off_plane.max((position - origin).dot(normal).abs());
                }
            }
        }
    }
    Ok(if off_plane == 0.0 {
        DeviationBound::Proven(0.0)
    } else {
        DeviationBound::Certified(off_plane * (1.0 + 1e-12))
    })
}

/// The bound of one curved face, gathered while it is meshed (#232).
struct FaceBound<'s> {
    oracle: Option<SurfaceBoundOracle<'s>>,
    /// The worst triangle bound: interpolation plus vertex offsets.
    interior: Scalar,
    /// The worst trim lens or pcurve gap term.
    trim: Scalar,
    unbounded: Option<&'static str>,
}

impl<'s> FaceBound<'s> {
    fn new(surface: &'s axiolid_surface::Surface) -> Self {
        let oracle = SurfaceBoundOracle::new(surface);
        Self {
            unbounded: oracle
                .is_none()
                .then_some("curved face whose surface has no derivative bound"),
            oracle,
            interior: 0.0,
            trim: 0.0,
        }
    }

    fn unbounded(&mut self, reason: &'static str) {
        self.unbounded.get_or_insert(reason);
    }

    fn finish(&self) -> DeviationBound {
        match self.unbounded {
            Some(reason) => DeviationBound::Unbounded(reason),
            None => DeviationBound::Certified((self.interior + self.trim) * (1.0 + 1e-12)),
        }
    }

    /// `sqrt(|S_u|^2 + |S_v|^2)` over a box, the surface's Lipschitz bound
    /// in parameter distance.
    fn lipschitz(&self, lo: Point2, hi: Point2) -> Option<Scalar> {
        let b = self.oracle.as_ref()?.bounds((lo.x, hi.x), (lo.y, hi.y))?;
        Some((b.du * b.du + b.dv * b.dv).sqrt())
    }

    /// The lens terms of one pcurve sampled as `n` uniform spans.
    fn trim(&mut self, curve: &axiolid_curve::Curve2, n: usize, params: &[Point2]) {
        if self.unbounded.is_some() {
            return;
        }
        let domain = axiolid_reference::curve::domain2(curve);
        let at =
            |i: usize| domain.start + (domain.end - domain.start) * (i as Scalar) / (n as Scalar);
        for i in 0..n {
            let Some(d) = axiolid_reference::bound::chord_bound2(curve, at(i), at(i + 1)) else {
                self.unbounded("trim pcurve family with no chord bound");
                return;
            };
            if d == 0.0 {
                continue;
            }
            let (a, b) = (params[i], params[i + 1]);
            let grow = Point2::splat(d);
            let Some(g) = self.lipschitz(a.min(b) - grow, a.max(b) + grow) else {
                self.unbounded("surface has no derivative bound over a trim");
                return;
            };
            self.trim = self.trim.max(g * d);
        }
    }

    /// A gap between one pcurve's end and the next one's start.
    fn gap(&mut self, end: Point2, start: Point2, surface: &axiolid_surface::Surface) {
        if self.unbounded.is_some() || end == start {
            return;
        }
        let (u_period, v_period) = surface_periods(surface);
        let near = Point2::new(
            localized_parameter(end.x, start.x, u_period),
            localized_parameter(end.y, start.y, v_period),
        );
        let gap = (near - end).length();
        if gap == 0.0 {
            return;
        }
        match self.lipschitz(end.min(near), end.max(near)) {
            Some(g) => self.trim = self.trim.max(g * gap),
            None => self.unbounded("surface has no derivative bound over a trim gap"),
        }
    }

    /// The final triangles of one face.
    fn triangles(
        &mut self,
        mesh: &TriMesh,
        surface: &axiolid_surface::Surface,
        triangles: &[SurfaceTriangle],
    ) -> GeomResult<()> {
        if self.unbounded.is_some() {
            return Ok(());
        }
        // Edges of the emitted triangles, for the dropped ones.
        let mut emitted: HashSet<(u32, u32)> = HashSet::new();
        for triangle in triangles {
            let [a, b, c] = triangle.vertices.map(|vertex| vertex.mesh);
            if a != b && b != c && c != a {
                for (x, y) in [(a, b), (b, c), (c, a)] {
                    emitted.insert(sorted_local_edge(x, y));
                }
            }
        }
        // Keyed by vertex AND parameters: a seam vertex has two occurrences.
        let mut offsets: HashMap<(u32, u64, u64), Scalar> = HashMap::new();
        for triangle in triangles {
            let [a, b, c] = triangle.vertices.map(|vertex| vertex.mesh);
            if a == b || b == c || c == a {
                let pair = if a == b { (b, c) } else { (a, b) };
                let covered = if pair.0 == pair.1 {
                    emitted.iter().any(|&(x, y)| x == pair.0 || y == pair.0)
                } else {
                    emitted.contains(&sorted_local_edge(pair.0, pair.1))
                };
                if !covered {
                    self.unbounded("a dropped degenerate triangle is not covered by the mesh");
                    return Ok(());
                }
            }
            let mut offset: Scalar = 0.0;
            for vertex in triangle.vertices {
                let key = (vertex.mesh, vertex.uv.x.to_bits(), vertex.uv.y.to_bits());
                let value = match offsets.get(&key) {
                    Some(&value) => value,
                    None => {
                        let exact = axiolid_reference::surface::evaluate(
                            surface,
                            vertex.uv.x,
                            vertex.uv.y,
                        )?;
                        let value = (mesh.positions[vertex.mesh as usize] - exact).length();
                        offsets.insert(key, value);
                        value
                    }
                };
                offset = offset.max(value);
            }
            let Some(interp) = self.interpolation(triangle.vertices.map(|vertex| vertex.uv)) else {
                self.unbounded("surface has no derivative bound over a triangle");
                return Ok(());
            };
            self.interior = self.interior.max(interp + offset);
        }
        Ok(())
    }

    /// The bound on `|S - L|` over one parameter triangle.
    fn interpolation(&self, uv: [Point2; 3]) -> Option<Scalar> {
        let lo = uv[0].min(uv[1]).min(uv[2]);
        let hi = uv[0].max(uv[1]).max(uv[2]);
        let b = self.oracle.as_ref()?.bounds((lo.x, hi.x), (lo.y, hi.y))?;
        let first = [(0, 1), (1, 2), (2, 0)]
            .iter()
            .map(|&(i, j)| {
                let d = uv[j] - uv[i];
                b.du * d.x.abs() + b.dv * d.y.abs()
            })
            .fold(0.0, Scalar::max);
        if !b.smooth {
            return Some(first);
        }
        // `2 B |d_u d_v| <= B (alpha d_u^2 + d_v^2 / alpha)` for any
        // `alpha > 0`: each choice is a bound, so a few are tried, centred
        // on the triangle's own aspect (a ruled cone's triangles are long
        // along the rulings, where `C = 0`).
        let extent = hi - lo;
        let aspect = if extent.x > 0.0 && extent.y > 0.0 {
            extent.y / extent.x
        } else {
            1.0
        };
        let mut best = first;
        for alpha in [
            1.0,
            aspect / 4.0,
            aspect / 2.0,
            aspect,
            2.0 * aspect,
            4.0 * aspect,
        ] {
            if !(alpha > 0.0 && alpha.is_finite()) {
                continue;
            }
            // Either coefficient set bounds the error on its own (they
            // differ only for a rational spline); the better one is kept.
            for [a2, b2, c2] in [b.interpolation, [b.duu, b.duv, b.dvv]] {
                let (su, sv) = ((a2 + b2 * alpha).sqrt(), (c2 + b2 / alpha).sqrt());
                let z = uv.map(|p| Point2::new((p.x - lo.x) * su, (p.y - lo.y) * sv));
                let radius = enclosing_radius(z);
                best = best.min(0.5 * radius * radius);
            }
        }
        Some(best)
    }

    /// The triangle's widest edge in the metric its bound is measured in,
    /// when that edge is not a trim edge: splitting a narrower one leaves
    /// the bound as wide as it was.
    fn widest_free_edge(&self, triangle: &SurfaceTriangle) -> Option<usize> {
        let uv = triangle.vertices.map(|vertex| vertex.uv);
        let lo = uv[0].min(uv[1]).min(uv[2]);
        let hi = uv[0].max(uv[1]).max(uv[2]);
        let b = self.oracle.as_ref()?.bounds((lo.x, hi.x), (lo.y, hi.y))?;
        let [x, y, z] = b.interpolation;
        let width = |edge: usize| {
            let (p, q) = triangle_edge(triangle.vertices, edge);
            let d = q.uv - p.uv;
            x * d.x * d.x + 2.0 * y * (d.x * d.y).abs() + z * d.y * d.y
        };
        let widest = (0..3).max_by(|&i, &j| width(i).total_cmp(&width(j)))?;
        (!triangle.boundary[widest]).then_some(widest)
    }

    /// One triangle's bound while the face is refined: interpolation plus
    /// its corners' distance off the surface; `None` when unbounded.
    fn triangle_bound(
        &self,
        mesh: &TriMesh,
        surface: &axiolid_surface::Surface,
        vertices: [SurfaceVertex; 3],
    ) -> GeomResult<Option<Scalar>> {
        if self.unbounded.is_some() {
            return Ok(None);
        }
        let mut offset: Scalar = 0.0;
        for vertex in vertices {
            let exact = axiolid_reference::surface::evaluate(surface, vertex.uv.x, vertex.uv.y)?;
            offset = offset.max((mesh.positions[vertex.mesh as usize] - exact).length());
        }
        Ok(self
            .interpolation(vertices.map(|vertex| vertex.uv))
            .map(|interp| interp + offset))
    }
}

/// Refinement depth below which a triangle whose certified bound misses
/// the budget has its widest free edge split (#232); deeper ones are
/// accepted and their bound reported as it is.
const CERTIFIED_DEPTH: u8 = 24;
/// The depth guard of the pass that makes certified splits. Earcut fans
/// thin triangles out of a trim corner, and narrowing them by halving takes
/// more generations than the measurement checks ever asked for; the face's
/// vertex budget still bounds the work, and a pass that exhausts either is
/// discarded (see [`refine_curved_face`]).
const MAX_CERTIFIED_REFINEMENT_DEPTH: u8 = 40;
/// Vertices a face may reach before certified splits stop. Beyond it the
/// bound is reported as it stands rather than paid for: a spline's
/// derivative-net bounds can ask for many times the triangles the
/// measurement checks need.
const CERTIFIED_FACE_VERTICES: usize = 1 << 15;

/// Put back the trim samples earcut skipped as collinear.
///
/// Earcut drops a ring point that lies on the line through its neighbours,
/// so a straight trim sampled into many points (as a shared edge is, for
/// the curved face across it) comes back as one long triangle edge. The
/// neighbouring face still uses every sample, which leaves T-junctions in
/// the shell, and the long edge is not a trim edge any more, so refinement
/// splits it at points the neighbour does not have (#232). Each skipped
/// sample lying strictly inside a triangle edge splits that triangle in
/// two, which keeps the triangles' parameter area.
fn restore_collinear_samples(flat: &[[Scalar; 2]], mut indices: Vec<usize>) -> Vec<usize> {
    let mut used = vec![false; flat.len()];
    for &i in &indices {
        used[i] = true;
    }
    for (point, _) in used.iter().enumerate().filter(|(_, &u)| !u) {
        let p = flat[point];
        'triangles: for t in 0..indices.len() / 3 {
            for edge in 0..3 {
                let (a, b) = (indices[3 * t + edge], indices[3 * t + (edge + 1) % 3]);
                let (pa, pb) = (flat[a], flat[b]);
                let (dx, dy) = (pb[0] - pa[0], pb[1] - pa[1]);
                let (ex, ey) = (p[0] - pa[0], p[1] - pa[1]);
                let length2 = dx * dx + dy * dy;
                let along = ex * dx + ey * dy;
                let cross = dx * ey - dy * ex;
                if !(along > 0.0 && along < length2) || cross.abs() > 1e-12 * length2 {
                    continue;
                }
                let c = indices[3 * t + (edge + 2) % 3];
                indices[3 * t] = a;
                indices[3 * t + 1] = point;
                indices[3 * t + 2] = c;
                indices.extend([point, b, c]);
                break 'triangles;
            }
        }
    }
    indices
}

/// Scales for u and v that make them about equally long on the surface:
/// the bounds of `|S_u|` and `|S_v|` over the trim's box, or unit scales
/// when the surface has none.
fn parameter_metric(surface: &axiolid_surface::Surface, uv: &[Point2]) -> (Scalar, Scalar) {
    let Some(oracle) = SurfaceBoundOracle::new(surface) else {
        return (1.0, 1.0);
    };
    let lo = uv
        .iter()
        .fold(Point2::splat(Scalar::INFINITY), |m, p| m.min(*p));
    let hi = uv
        .iter()
        .fold(Point2::splat(Scalar::NEG_INFINITY), |m, p| m.max(*p));
    match oracle.bounds((lo.x, hi.x), (lo.y, hi.y)) {
        Some(b) if b.du > 0.0 && b.dv > 0.0 && b.du.is_finite() && b.dv.is_finite() => (b.du, b.dv),
        _ => (1.0, 1.0),
    }
}

/// Flip interior edges of a polygon triangulation until it is Delaunay in
/// the parameter plane (Lawson), never flipping a trim edge.
///
/// Putting skipped samples back fans thin triangles out of one corner;
/// slivers that long put their 3D triangles off the surface by more than
/// the surface itself bends, and leave the refinement nothing it can
/// split usefully (#232). A flip only replaces the diagonal of a convex
/// quad, so the triangles still tile the same polygon.
fn delaunay_flips(
    flat: &[[Scalar; 2]],
    mut indices: Vec<usize>,
    ring_edges: &HashSet<LocalEdgeKey>,
) -> Vec<usize> {
    let key = |a: usize, b: usize| sorted_local_edge(a as u32, b as u32);
    let orient = |a: usize, b: usize, c: usize| {
        let (pa, pb, pc) = (flat[a], flat[b], flat[c]);
        (pb[0] - pa[0]) * (pc[1] - pa[1]) - (pb[1] - pa[1]) * (pc[0] - pa[0])
    };
    let incircle = |a: usize, b: usize, c: usize, d: usize| {
        let row = |p: usize| {
            let (x, y) = (flat[p][0] - flat[d][0], flat[p][1] - flat[d][1]);
            [x, y, x * x + y * y]
        };
        let (ra, rb, rc) = (row(a), row(b), row(c));
        ra[0] * (rb[1] * rc[2] - rb[2] * rc[1]) - ra[1] * (rb[0] * rc[2] - rb[2] * rc[0])
            + ra[2] * (rb[0] * rc[1] - rb[1] * rc[0])
    };
    let triangles = indices.len() / 3;
    // A pass flips every edge that fails; Lawson's flips terminate, the
    // cap only guards against rounding cycling on cocircular points.
    for _ in 0..4 * triangles.max(1) {
        let mut owners: HashMap<LocalEdgeKey, Vec<(usize, usize)>> = HashMap::new();
        for t in 0..triangles {
            for e in 0..3 {
                let (a, b) = (indices[3 * t + e], indices[3 * t + (e + 1) % 3]);
                owners.entry(key(a, b)).or_default().push((t, e));
            }
        }
        let mut flipped = false;
        let mut touched = vec![false; triangles];
        for (edge, uses) in &owners {
            if uses.len() != 2 || ring_edges.contains(edge) {
                continue;
            }
            let ((t1, e1), (t2, e2)) = (uses[0], uses[1]);
            if touched[t1] || touched[t2] {
                continue;
            }
            let (a, b, c) = (
                indices[3 * t1 + e1],
                indices[3 * t1 + (e1 + 1) % 3],
                indices[3 * t1 + (e1 + 2) % 3],
            );
            let d = indices[3 * t2 + (e2 + 2) % 3];
            let side = orient(a, b, c);
            // A convex quad a, d, b, c: d across ab from c, and a, b on
            // opposite sides of cd.
            let convex = side * orient(a, b, d) < 0.0 && orient(c, d, a) * orient(c, d, b) < 0.0;
            if !convex || side * incircle(a, b, c, d) <= 0.0 {
                continue;
            }
            indices[3 * t1..3 * t1 + 3].copy_from_slice(&[a, d, c]);
            indices[3 * t2..3 * t2 + 3].copy_from_slice(&[d, b, c]);
            touched[t1] = true;
            touched[t2] = true;
            flipped = true;
        }
        if !flipped {
            break;
        }
    }
    indices
}

/// Radius of the smallest circle enclosing a triangle: half the longest
/// side when the triangle is right, obtuse or flat, else its circumradius.
fn enclosing_radius([a, b, c]: [Point2; 3]) -> Scalar {
    let (ab, bc, ca) = (
        (b - a).length_squared(),
        (c - b).length_squared(),
        (a - c).length_squared(),
    );
    let mut sides = [ab, bc, ca];
    sides.sort_by(Scalar::total_cmp);
    let [s0, s1, longest] = sides;
    if s0 + s1 <= longest {
        return 0.5 * longest.sqrt();
    }
    let twice_area = ((b - a).x * (c - a).y - (b - a).y * (c - a).x).abs();
    if twice_area <= 0.0 || twice_area.is_nan() {
        return 0.5 * longest.sqrt();
    }
    // R = |ab| |bc| |ca| / (4 area), never below half the longest side.
    ((ab * bc * ca).sqrt() / (2.0 * twice_area)).max(0.5 * longest.sqrt())
}

/// Whether the triangles' parameter area matches the trim polygon's, to a
/// relative `1e-9`.
fn covers_polygon(boundary: &CurvedBoundary, triangles: &[SurfaceTriangle]) -> bool {
    let signed = |ring: &[Point2]| -> Scalar {
        let n = ring.len();
        (0..n)
            .map(|i| {
                let (p, q) = (ring[i], ring[(i + 1) % n]);
                p.x * q.y - q.x * p.y
            })
            .sum::<Scalar>()
            * 0.5
    };
    let mut starts = vec![0];
    starts.extend_from_slice(&boundary.hole_starts);
    starts.push(boundary.uv.len());
    let mut polygon = 0.0;
    for (k, w) in starts.windows(2).enumerate() {
        let area = signed(&boundary.uv[w[0]..w[1]]).abs();
        polygon += if k == 0 { area } else { -area };
    }
    let covered: Scalar = triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.vertices.map(|vertex| vertex.uv);
            0.5 * ((b - a).x * (c - a).y - (b - a).y * (c - a).x).abs()
        })
        .sum();
    (covered - polygon).abs() <= 1e-9 * polygon.abs().max(covered)
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::{Frame3, Point3};
    use axiolid_curve::{Curve2, Polyline2};
    use axiolid_surface::{Cylinder, Plane, Surface};

    #[test]
    fn curved_edge_sampling_fails_closed_when_the_budget_cannot_meet_tolerance() {
        let trim = Curve2::Polyline(Polyline2 {
            points: vec![
                Point2::new(0.0, 0.0),
                Point2::new(core::f64::consts::TAU, 0.0),
            ],
            closed: false,
        });
        let surface = Surface::Cylinder(Cylinder {
            frame: Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius: 2.0,
        });

        let error = edge_sample_count(&trim, &surface, 1.0e-20)
            .expect_err("the boundary sampler must not return an unverified polyline");
        assert!(matches!(
            error,
            GeomError::BudgetExceeded {
                resource: "curved edge samples"
            }
        ));
    }

    #[test]
    fn adjacent_triangle_refinement_never_leaves_one_side_of_a_shared_edge() {
        let surface = Surface::Cylinder(Cylinder {
            frame: Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius: 2.0,
        });
        let uv = [
            Point2::new(0.0, 0.0),
            Point2::new(0.01, 1.0),
            Point2::new(core::f64::consts::FRAC_PI_2, 0.0),
            Point2::new(0.02, 0.0),
        ];
        let mut mesh = TriMesh {
            positions: uv
                .iter()
                .map(|p| axiolid_reference::surface::evaluate(&surface, p.x, p.y).unwrap())
                .collect(),
            indices: Vec::new(),
            normals: None,
            attributes: Vec::new(),
        };
        let mut welded = WeldedVertices::new();
        let mut edge_samples = EdgeSamples::new();
        let mut total_curved_records = 4;
        let mut state = CurvedMeshState {
            mesh: &mut mesh,
            welded: &mut welded,
            edge_samples: &mut edge_samples,
            chord_error: 1.0e-3,
            face_local_vertices: 4,
            total_curved_records: &mut total_curved_records,
            new_welded: Vec::new(),
            new_edge_samples: Vec::new(),
        };
        let vertex = |index: u32| SurfaceVertex {
            uv: uv[index as usize],
            mesh: index,
            local: index,
        };
        let triangles = vec![
            SurfaceTriangle {
                vertices: [vertex(0), vertex(1), vertex(2)],
                boundary: [false; 3],
                depth: 0,
            },
            SurfaceTriangle {
                vertices: [vertex(1), vertex(0), vertex(3)],
                boundary: [false; 3],
                depth: 0,
            },
        ];

        let mut bound = FaceBound::new(&surface);
        refine_curved_face(&mut state, &surface, triangles, 4, false, &mut bound).unwrap();

        let shared = (0_u32, 1_u32);
        let uses = state
            .mesh
            .indices
            .chunks_exact(3)
            .flat_map(|triangle| {
                [
                    (triangle[0], triangle[1]),
                    (triangle[1], triangle[2]),
                    (triangle[2], triangle[0]),
                ]
            })
            .filter(|&(a, b)| (if a < b { (a, b) } else { (b, a) }) == shared)
            .count();
        assert_ne!(
            uses, 1,
            "a shared edge must be retained by both triangles or split by both"
        );

        let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
        for triangle in state.mesh.indices.chunks_exact(3) {
            for (a, b) in [
                (triangle[0], triangle[1]),
                (triangle[1], triangle[2]),
                (triangle[2], triangle[0]),
            ] {
                *counts.entry(sorted_local_edge(a, b)).or_default() += 1;
            }
        }
        let boundary = [
            (uv[1], uv[2]),
            (uv[2], uv[0]),
            (uv[0], uv[3]),
            (uv[3], uv[1]),
        ];
        let parameter = |index: u32| {
            let point = state.mesh.positions[index as usize];
            Point2::new(point.y.atan2(point.x), point.z)
        };
        let on_segment = |point: Point2, start: Point2, end: Point2| {
            let edge = end - start;
            let offset = point - start;
            let cross = edge.x * offset.y - edge.y * offset.x;
            cross.abs() <= 1.0e-8
                && offset.dot(edge) >= -1.0e-8
                && offset.dot(edge) <= edge.length_squared() + 1.0e-8
        };
        for ((a, b), count) in counts {
            if count == 1 {
                let (pa, pb) = (parameter(a), parameter(b));
                assert!(
                    boundary.iter().any(|&(start, end)| {
                        on_segment(pa, start, end) && on_segment(pb, start, end)
                    }),
                    "unpaired interior edge ({a}, {b}) at {pa:?}->{pb:?}"
                );
            } else {
                assert_eq!(count, 2, "non-manifold edge ({a}, {b})");
            }
        }
    }

    #[test]
    fn refinement_midpoints_are_cached_only_within_one_face() {
        let frame = Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        };
        let cylinder = Surface::Cylinder(Cylinder { frame, radius: 2.0 });
        let plane = Surface::Plane(Plane { frame });
        let mut mesh = TriMesh {
            positions: vec![Point3::ZERO, Point3::new(1.0, 1.0, 1.0)],
            indices: Vec::new(),
            normals: None,
            attributes: Vec::new(),
        };
        let mut welded = WeldedVertices::new();
        let mut edge_samples = EdgeSamples::new();
        let mut total_curved_records = 2;
        let mut state = CurvedMeshState {
            mesh: &mut mesh,
            welded: &mut welded,
            edge_samples: &mut edge_samples,
            chord_error: 1.0e-3,
            face_local_vertices: 2,
            total_curved_records: &mut total_curved_records,
            new_welded: Vec::new(),
            new_edge_samples: Vec::new(),
        };
        let a = SurfaceVertex {
            uv: Point2::new(0.0, 0.0),
            mesh: 0,
            local: 0,
        };
        let b = SurfaceVertex {
            uv: Point2::new(core::f64::consts::FRAC_PI_2, 1.0),
            mesh: 1,
            local: 1,
        };
        let mut first_face_cache = HashMap::new();
        let mut first_next = 2;
        let first = refinement_midpoint(
            &mut state,
            &cylinder,
            a,
            b,
            &mut first_face_cache,
            &mut first_next,
        )
        .unwrap();
        let mut second_face_cache = HashMap::new();
        state.face_local_vertices = 2;
        let mut second_next = 2;
        let second = refinement_midpoint(
            &mut state,
            &plane,
            a,
            b,
            &mut second_face_cache,
            &mut second_next,
        )
        .unwrap();

        assert_ne!(first.mesh, second.mesh);
        assert!(
            (state.mesh.positions[first.mesh as usize]
                - state.mesh.positions[second.mesh as usize])
                .length()
                > 0.1
        );
    }

    #[test]
    fn the_enclosing_circle_of_an_acute_triangle_is_its_circumcircle() {
        let s = 3.0_f64.sqrt();
        let equilateral = [
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(1.0, s),
        ];
        assert!((enclosing_radius(equilateral) - 2.0 / s).abs() < 1e-12);
        let obtuse = [
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(1.0, 0.1),
        ];
        assert!((enclosing_radius(obtuse) - 1.0).abs() < 1e-12);
    }

    /// A seam vertex shared from the face across sits off this face's
    /// surface; that offset is part of the bound, since interpolation
    /// cannot reach it.
    #[test]
    fn a_vertex_off_its_surface_widens_the_triangle_bound() {
        let surface = Surface::Cylinder(Cylinder {
            frame: Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius: 2.0,
        });
        let uv = [
            Point2::new(0.0, 0.0),
            Point2::new(0.01, 0.0),
            Point2::new(0.0, 0.01),
        ];
        let mut positions: Vec<Vec3> = uv
            .iter()
            .map(|p| axiolid_reference::surface::evaluate(&surface, p.x, p.y).unwrap())
            .collect();
        positions[1] += Vec3::new(0.0, 0.0, 0.01);
        let mesh = TriMesh {
            positions,
            indices: vec![0, 1, 2],
            normals: None,
            attributes: Vec::new(),
        };
        let vertex = |i: u32| SurfaceVertex {
            uv: uv[i as usize],
            mesh: i,
            local: i,
        };
        let triangle = SurfaceTriangle {
            vertices: [vertex(0), vertex(1), vertex(2)],
            boundary: [false; 3],
            depth: 0,
        };
        let mut bound = FaceBound::new(&surface);
        bound.triangles(&mesh, &surface, &[triangle]).unwrap();
        assert!(bound.interior >= 0.01, "{}", bound.interior);
        assert!(bound.interior < 0.0101, "{}", bound.interior);
    }

    #[test]
    fn grid_rejects_bow_tie() {
        let uv = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        ];
        let b = CurvedBoundary {
            uv,
            shared: vec![0; 4],
            hole_starts: vec![],
            winding_reversed: false,
        };

        assert!(recognise_grid(&b).is_none());
    }

    #[test]
    fn curved_rings_require_one_complete_outer() {
        assert!(validate_curved_rings(&[(0, 2, true, false)]).is_err());
        assert!(validate_curved_rings(&[(0, 3, false, false)]).is_err());
        assert!(validate_curved_rings(&[(0, 3, true, false), (3, 6, true, false)]).is_err());
        assert_eq!(
            validate_curved_rings(&[(0, 3, false, false), (3, 6, true, true)]),
            Ok(1)
        );
    }

    #[test]
    fn curved_face_failure_rolls_back_partial_output() {
        let mut mesh = TriMesh {
            positions: vec![Point3::new(1.0, 2.0, 3.0)],
            indices: vec![0],
            normals: None,
            attributes: Vec::new(),
        };
        let before_positions = mesh.positions.clone();
        let before_indices = mesh.indices.clone();
        let mut welded = WeldedVertices::new();
        let mut edge_samples = EdgeSamples::new();
        let mut records = 7;
        let result = with_curved_face_transaction(
            &mut mesh,
            &mut welded,
            &mut edge_samples,
            &mut records,
            1.0e-3,
            |state| {
                state.reserve_face_vertices(1)?;
                state.push_face_position(Point3::new(4.0, 5.0, 6.0))?;
                state.mesh.indices.push(1);
                Err(GeomError::BudgetExceeded {
                    resource: "mesh positions",
                })
            },
        );
        assert!(result.is_err());
        assert_eq!(mesh.positions, before_positions);
        assert_eq!(mesh.indices, before_indices);
        assert_eq!(records, 7);
    }

    #[test]
    fn mesh_growth_is_checked_before_mutation() {
        let positions = checked_output_len(
            MAX_MESH_POSITIONS - 1,
            1,
            MAX_MESH_POSITIONS,
            "mesh positions",
        )
        .unwrap();
        assert_eq!(positions, MAX_MESH_POSITIONS);
        assert!(
            checked_output_len(MAX_MESH_POSITIONS, 1, MAX_MESH_POSITIONS, "mesh positions")
                .is_err()
        );
        assert!(checked_output_len(MAX_MESH_INDICES, 3, MAX_MESH_INDICES, "mesh indices").is_err());
        assert!(checked_output_len(usize::MAX, 1, usize::MAX, "mesh indices").is_err());
        let mut output = vec![7_u32];
        assert!(extend_output(&mut output, &[8], 1, "mesh indices").is_err());
        assert_eq!(output, vec![7]);
    }

    #[test]
    fn expanded_tessellation_work_is_bounded_transactionally() {
        let mut total = MAX_BREP_EDGE_USES - 1;
        consume_tessellation_work(&mut total, 1).expect("the exact budget is admitted");
        assert_eq!(total, MAX_BREP_EDGE_USES);
        let error = consume_tessellation_work(&mut total, 1)
            .expect_err("expanded shell references must not exceed the budget");
        assert!(matches!(
            error,
            GeomError::BudgetExceeded {
                resource: "expanded B-rep tessellation work"
            }
        ));
        assert_eq!(total, MAX_BREP_EDGE_USES);
    }

    #[test]
    fn curved_face_vertex_budget_ignores_prior_faces_and_bounds_local_work() {
        let mut mesh = TriMesh {
            positions: vec![Point3::ZERO; MAX_CURVED_FACE_VERTICES + 1],
            indices: Vec::new(),
            normals: None,
            attributes: Vec::new(),
        };
        let mut welded = WeldedVertices::new();
        let mut edge_samples = EdgeSamples::new();
        let mut total_curved_records = 0;
        let mut state = CurvedMeshState {
            mesh: &mut mesh,
            welded: &mut welded,
            edge_samples: &mut edge_samples,
            chord_error: 1.0e-3,
            face_local_vertices: 0,
            total_curved_records: &mut total_curved_records,
            new_welded: Vec::new(),
            new_edge_samples: Vec::new(),
        };
        state
            .reserve_face_vertices(1)
            .expect("prior shell geometry must not consume this face's budget");
        let error = state
            .reserve_face_vertices(MAX_CURVED_FACE_VERTICES)
            .expect_err("boundary plus refinement work must share one face-local budget");
        assert!(matches!(
            error,
            GeomError::BudgetExceeded {
                resource: "curved-face local vertices"
            }
        ));

        state.face_local_vertices = 0;
        *state.total_curved_records = MAX_TOTAL_CURVED_RECORDS;
        let error = state
            .reserve_face_vertices(1)
            .expect_err("curved-face records must also be bounded across the mesh");
        assert!(matches!(
            error,
            GeomError::BudgetExceeded {
                resource: "total curved-face records"
            }
        ));
        assert_eq!(*state.total_curved_records, MAX_TOTAL_CURVED_RECORDS);
    }

    #[test]
    fn an_elliptical_cylinder_is_periodic_in_u() {
        // The seam only closes if the u-period is known. A new surface that
        // falls into the catch-all arm silently reports (None, None), which
        // leaves the wall split open along its seam instead of wrapping.
        let frame = Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        };
        let surface = Surface::EllipticalCylinder(axiolid_surface::EllipticalCylinder {
            frame,
            semi_axis_x: 3.0,
            semi_axis_y: 1.0,
        });
        let (u, v) = surface_periods(&surface);
        assert_eq!(u, Some(core::f64::consts::TAU), "u must wrap");
        assert_eq!(v, None, "the axis direction does not wrap");
    }
}
