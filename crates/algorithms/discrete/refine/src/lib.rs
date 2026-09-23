//! Mesh refinement: more triangles, and optionally closer to the truth.
//!
//! Splitting a triangle is easy. The question this module exists to answer
//! is *where the new vertex goes*.
//!
//! A mesh-only kernel has one option: the edge midpoint. That subdivides
//! the approximation without improving it -- refining a faceted cylinder
//! forever leaves the same faceted cylinder, with more triangles.
//!
//! When the mesh came from a tessellated B-rep and the source surface is
//! still known, there is a better answer: invert the midpoint into the
//! surface's parameter domain and evaluate the surface there. The new
//! vertex lands on the ACTUAL cylinder. Refinement then converges on the
//! real geometry instead of preserving a facet forever.
//!
//! Keeping the analytic surface alongside the mesh is what makes that
//! possible, so it is the capability this module is really for.

pub mod smooth;

use ahash::AHashMap;

use axiolid_core::{Point3, Scalar, Tolerance};
use axiolid_mesh::{AttributeFate, TriMesh};
use axiolid_surface::Surface;

/// Why a refinement could not be performed.
#[derive(Debug, thiserror::Error, PartialEq)]
#[non_exhaustive]
pub enum RefineError {
    /// The index buffer is not a whole number of triangles.
    #[error("index buffer length {0} is not a multiple of 3")]
    RaggedIndices(usize),
    /// A triangle references a vertex that does not exist.
    #[error("triangle {0} references vertex {1}, which is out of range")]
    IndexOutOfRange(usize, u32),
    /// An edge-length target must be a positive, finite length.
    #[error("edge length target {0} is not a positive finite length")]
    InvalidTarget(Scalar),
    /// A surface refused to place a projected vertex.
    ///
    /// Propagated rather than absorbed. Falling back to the linear midpoint
    /// would return four times the triangles with none of the promised
    /// accuracy, and the caller could not tell the difference.
    #[error("surface-aware refinement refused: {0}")]
    SurfaceRefused(String),
    /// Refinement would exceed the triangle budget.
    ///
    /// Reported rather than silently truncated: a caller that asked for a
    /// 1mm edge on a building-sized model wants to know its request was
    /// impossible, not receive a partially refined mesh that looks fine.
    #[error("refinement would produce {produced} triangles, over the {limit} budget")]
    BudgetExceeded {
        /// Triangles the request would have produced.
        produced: usize,
        /// The cap that was not raised.
        limit: usize,
    },
}

/// How much to refine.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum RefineTarget {
    /// Split every triangle into four, `levels` times.
    Uniform {
        /// Number of subdivision passes.
        levels: u32,
    },
    /// Split edges until none is longer than this.
    EdgeLength {
        /// Maximum permitted edge length, in model units.
        max_edge: Scalar,
    },
}

/// What a refinement actually did.
///
/// `max_deviation` is the honest part. A planar refinement must report
/// exactly zero: a midpoint on a flat triangle lies in that triangle's
/// plane, so any movement means a defect. A surface-aware refinement
/// reports how far it MOVED the surface toward the analytic one, which is
/// the measure of what the caller gained.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct RefineReport {
    /// Triangles before.
    pub input_triangles: usize,
    /// Triangles after.
    pub output_triangles: usize,
    /// Vertices introduced.
    pub vertices_added: usize,
    /// Whether new vertices were placed on an analytic surface.
    ///
    /// `false` means linear midpoints: the result is a finer tessellation
    /// of the same approximation, not a better approximation.
    pub surface_aware: bool,
    /// Largest distance a new vertex sits from the linear midpoint it
    /// would otherwise have occupied, in model units.
    ///
    /// Zero for planar input even when surface-aware, because a plane's
    /// midpoint already lies on the plane.
    pub max_deviation: Scalar,
    /// What happened to each named attribute channel.
    pub attribute_fates: Vec<(String, AttributeFate)>,
}

impl RefineReport {
    /// Whether the mesh was left untouched.
    pub fn is_noop(&self) -> bool {
        self.vertices_added == 0
    }
}

/// Cap on output size, mirroring the budget discipline used elsewhere.
const MAX_TRIANGLES: usize = 20_000_000;

fn validate(mesh: &TriMesh) -> Result<(), RefineError> {
    if mesh.indices.len() % 3 != 0 {
        return Err(RefineError::RaggedIndices(mesh.indices.len()));
    }
    let vertex_count = mesh.positions.len();
    for (triangle, chunk) in mesh.indices.chunks_exact(3).enumerate() {
        for &index in chunk {
            if index as usize >= vertex_count {
                return Err(RefineError::IndexOutOfRange(triangle, index));
            }
        }
    }
    Ok(())
}

/// Refine a mesh, optionally snapping new vertices onto a known surface.
///
/// Passing `Some(surface)` is what turns subdivision into approximation
/// improvement. Passing `None` subdivides linearly and says so in the
/// report rather than implying an accuracy gain it did not deliver.
///
/// Deterministic: new vertices are numbered in the order edges are first
/// split during the triangle walk, which is fixed by the index buffer. The
/// midpoint cache is keyed by the ordered vertex pair and is only ever
/// queried by key, never iterated, so its internal ordering cannot reach
/// the output. The same input produces the same output vertex ordering on
/// every run and across processes.
///
/// # Errors
///
/// Refuses a ragged index buffer, out-of-range indices, a non-positive
/// edge-length target, and a request that would exceed the triangle budget.
pub fn refine(
    mesh: &TriMesh,
    target: RefineTarget,
    surface: Option<&Surface>,
    tolerance: Tolerance,
) -> Result<(TriMesh, RefineReport), RefineError> {
    validate(mesh)?;

    let levels = match target {
        RefineTarget::Uniform { levels } => levels,
        RefineTarget::EdgeLength { max_edge } => {
            if !max_edge.is_finite() || max_edge <= 0.0 {
                return Err(RefineError::InvalidTarget(max_edge));
            }
            passes_for_edge_length(mesh, max_edge)
        }
    };

    let input_triangles = mesh.triangle_count();
    // Each pass quadruples the triangle count. Checking the projection up
    // front turns an out-of-memory kill into a typed refusal.
    let projected = input_triangles
        .checked_mul(4usize.saturating_pow(levels))
        .unwrap_or(usize::MAX);
    if projected > MAX_TRIANGLES {
        return Err(RefineError::BudgetExceeded {
            produced: projected,
            limit: MAX_TRIANGLES,
        });
    }

    let mut positions = mesh.positions.clone();
    let mut indices = mesh.indices.clone();
    let mut max_deviation: Scalar = 0.0;

    for _ in 0..levels {
        let mut midpoints: AHashMap<(u32, u32), u32> = AHashMap::new();
        let mut next = Vec::with_capacity(indices.len() * 4);

        for chunk in indices.chunks_exact(3) {
            let [a, b, c] = [chunk[0], chunk[1], chunk[2]];
            let ab = split_edge(
                a,
                b,
                &mut positions,
                &mut midpoints,
                surface,
                tolerance,
                &mut max_deviation,
            )?;
            let bc = split_edge(
                b,
                c,
                &mut positions,
                &mut midpoints,
                surface,
                tolerance,
                &mut max_deviation,
            )?;
            let ca = split_edge(
                c,
                a,
                &mut positions,
                &mut midpoints,
                surface,
                tolerance,
                &mut max_deviation,
            )?;

            // Four children, each wound the same way as the parent so the
            // result keeps the input's orientation.
            next.extend_from_slice(&[a, ab, ca]);
            next.extend_from_slice(&[ab, b, bc]);
            next.extend_from_slice(&[ca, bc, c]);
            next.extend_from_slice(&[ab, bc, ca]);
        }
        indices = next;
    }

    let vertices_added = positions.len() - mesh.positions.len();
    let mut out = TriMesh::new(positions, indices);
    out.normals = None;

    // With no vertex created (zero levels, or no triangles to split) the
    // geometry is the input's, so its channels and normals are still exact
    // and are returned, not merely reported as surviving. Before this the
    // report said `Preserved` while the mesh came back without the channel.
    if vertices_added == 0 {
        out.normals = mesh.normals.clone();
        out.attributes = mesh.attributes.clone();
    }

    // A refinement creates vertices, so a channel survives only if its own
    // blend rule permits deriving a value. Unlike a boolean cut, the new
    // vertex HAS a preimage: it sits on a known edge between two vertices,
    // so a blendable channel is genuinely interpolatable here.
    let attribute_fates = mesh
        .attributes
        .iter()
        .map(|channel| {
            let fate = match channel.blend {
                // Checked first: nothing was derived, so even a channel that
                // forbids derivation came through untouched.
                _ if vertices_added == 0 => AttributeFate::Preserved,
                axiolid_mesh::Blend::None => {
                    AttributeFate::Dropped(axiolid_mesh::DropReason::NotBlendable)
                }
                _ => AttributeFate::Dropped(axiolid_mesh::DropReason::ProviderLimitation),
            };
            (channel.name.clone(), fate)
        })
        .collect();

    let report = RefineReport {
        input_triangles,
        output_triangles: out.triangle_count(),
        vertices_added,
        surface_aware: surface.is_some(),
        max_deviation,
        attribute_fates,
    };
    Ok((out, report))
}

/// Number of uniform passes needed to bring every edge under `max_edge`.
///
/// Each pass halves every edge, so the requirement is
/// `longest / 2^n <= max_edge`. Computed rather than iterated so the
/// budget check can happen before any memory is allocated.
fn passes_for_edge_length(mesh: &TriMesh, max_edge: Scalar) -> u32 {
    let mut longest: Scalar = 0.0;
    for chunk in mesh.indices.chunks_exact(3) {
        for (from, to) in [(0, 1), (1, 2), (2, 0)] {
            let a = mesh.positions[chunk[from] as usize];
            let b = mesh.positions[chunk[to] as usize];
            longest = longest.max((b - a).length());
        }
    }
    if longest <= max_edge || !longest.is_finite() {
        return 0;
    }
    (longest / max_edge).log2().ceil().max(0.0) as u32
}

/// Return the vertex splitting an edge, creating it on first encounter.
///
/// The edge key is ordered so both adjacent triangles find the same
/// midpoint. Without that the mesh would crack along every shared edge.
fn split_edge(
    a: u32,
    b: u32,
    positions: &mut Vec<Point3>,
    midpoints: &mut AHashMap<(u32, u32), u32>,
    surface: Option<&Surface>,
    tolerance: Tolerance,
    max_deviation: &mut Scalar,
) -> Result<u32, RefineError> {
    let key = if a < b { (a, b) } else { (b, a) };
    if let Some(&existing) = midpoints.get(&key) {
        return Ok(existing);
    }

    let linear = positions[a as usize].midpoint(positions[b as usize]);
    let placed = match surface {
        // PROJECT, not invert. The midpoint of a chord across a faceted
        // surface lies strictly off that surface, so inversion correctly
        // refuses it; projection is the question actually being asked.
        //
        // A refusal propagates instead of falling back to `linear`. A
        // silent fallback would return a mesh with four times the
        // triangles and none of the promised accuracy, which is worse
        // than an error because the caller cannot detect it.
        Some(surface) => {
            let (u, v) = axiolid_evaluate::surface::project(surface, linear, tolerance)
                .map_err(|error| RefineError::SurfaceRefused(error.to_string()))?;
            let on_surface = axiolid_evaluate::surface::evaluate(surface, u, v)
                .map_err(|error| RefineError::SurfaceRefused(error.to_string()))?;
            if !on_surface.is_finite() {
                return Err(RefineError::SurfaceRefused(
                    "projected midpoint is not finite".into(),
                ));
            }
            on_surface
        }
        None => linear,
    };

    *max_deviation = max_deviation.max((placed - linear).length());
    let index = positions.len() as u32;
    positions.push(placed);
    midpoints.insert(key, index);
    Ok(index)
}

/// What a smoothing pass actually did.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct SmoothReport {
    /// Vertices whose position changed.
    pub vertices_moved: usize,
    /// Vertices held fixed because they sit on an open border.
    pub boundary_vertices: usize,
    /// Largest distance any vertex moved, in model units.
    pub max_movement: Scalar,
    /// What happened to each named attribute channel.
    pub attribute_fates: Vec<(String, AttributeFate)>,
}

/// Report every channel as dropped.
///
/// Smoothing moves vertices without creating them, but the values it would
/// need to keep are only valid at the ORIGINAL positions. Rather than carry
/// stale data forward under its old name, every channel is dropped and said
/// to be dropped.
fn carry_attributes(mesh: &TriMesh) -> Vec<(String, AttributeFate)> {
    mesh.attributes
        .iter()
        .map(|channel| {
            (
                channel.name.clone(),
                AttributeFate::Dropped(axiolid_mesh::DropReason::ProviderLimitation),
            )
        })
        .collect()
}
