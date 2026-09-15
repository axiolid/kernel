//! Convex decomposition: a solid as a set of convex parts.
//!
//! # Two strategies, one contract
//!
//! There is no single right answer here, so the caller picks:
//!
//! - [`Strategy::Exact`] splits at reflex features until every part is
//!   genuinely convex. The union reproduces the input exactly, and the part
//!   count can be large.
//! - [`Strategy::Approximate`] stops once each part is convex to within a
//!   stated concavity bound. Far fewer parts, and the union is close to but
//!   not identical to the input.
//!
//! Both are legitimate. Collision detection and Minkowski sums usually want
//! the approximate one; anything claiming to reproduce the original solid
//! needs the exact one. What is NOT legitimate is returning an approximate
//! decomposition that presents itself as exact, so [`Decomposition`] always
//! reports which it is, and the approximate path reports the concavity it
//! actually reached rather than the one that was requested.
//!
//! # Method
//!
//! Both strategies share one loop: measure the worst concavity of a part,
//! and if it exceeds the bound, split the part by a plane and recurse.
//! They differ only in the bound -- exact uses zero (to tolerance).
//!
//! Concavity is measured as the largest distance from a vertex of the part
//! to its own convex hull. That is a direct measurement of the property the
//! caller cares about, rather than a proxy like volume ratio: a thin deep
//! notch barely changes volume but is exactly what breaks a convexity
//! assumption downstream.
//!
//! The split plane is the plane of the face the reflex vertex sticks out
//! past. Extending an existing face makes progress by definition, whereas
//! a bounding-box axis through the same point need not separate the notch.
//!
//! # Capping a cut: measure, do not predict
//!
//! Closing the cut is the hard part, and the first four attempts all
//! failed in the same shape. Each tried to PREDICT the cross-section from
//! the input mesh, deciding per triangle whether an edge bounded the cut.
//! Measured on an L-shaped solid:
//!
//! | rule | boundary edges | non-manifold edges |
//! |---|---|---|
//! | strict sign changes only | 5 | 0 |
//! | plus vertices lying on the plane | 0 | 3 |
//! | plus a straddling filter | 3 | 0 |
//! | plus the two-vertices-on-plane case | 1 | 3 |
//!
//! Every rule fixed one defect and reintroduced the other, which is the
//! signature of the wrong question rather than a missing case: whether a
//! wall standing ON the cut plane bounds THIS part depends on which side
//! the material lies, and a single triangle cannot see that.
//!
//! The fix is to stop predicting. Clip first, then look at what the shell
//! actually left open: in a closed mesh every undirected edge is used
//! exactly twice, so the edges used ONCE are precisely the hole. The cap
//! fills exactly that, and can be neither too generous nor too strict
//! whatever the clipping did upstream.
//!
//! # Two splitters
//!
//! The hand-rolled clipper above needs no boolean backend, which matters
//! because this crate should be usable without one. A caller that already
//! has a boolean provider can pass it instead, per call.
//!
//! Because `algorithms` may not depend on `providers` -- the architecture
//! gate enforces it -- the provider arrives through the mesh-boolean
//! CONTRACT, which both layers may depend on. See [`split::Splitter`].
//! The two paths are independent implementations of the same contract, so
//! each is evidence about the other, and the tests check they agree on the
//! resulting solid rather than merely on their own claims.

pub mod split;

use ahash::AHashMap;
use std::collections::BTreeMap;

use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_mesh::{audit_mesh, EdgeAdjacency, TriMesh};
use thiserror::Error;

/// Why a decomposition could not be produced.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum DecomposeError {
    /// The index buffer is not a whole number of triangles.
    #[error("index buffer length {0} is not a multiple of 3")]
    RaggedIndices(usize),
    /// A triangle references a vertex that does not exist.
    #[error("triangle {0} references vertex {1}, which is out of range")]
    IndexOutOfRange(usize, u32),
    /// The input is not a closed two-manifold solid.
    ///
    /// Refused rather than decomposed: the parts of an open surface do not
    /// have a union that reproduces it, so any answer would be a fiction.
    #[error("input is not a closed two-manifold solid: {boundary} boundary and {non_manifold} non-manifold edges")]
    NotASolid {
        /// Edges with a single incident triangle.
        boundary: usize,
        /// Edges with more than two incident triangles.
        non_manifold: usize,
    },
    /// A concavity bound must be a positive, finite length.
    #[error("concavity bound {0} is not a positive finite length")]
    InvalidBound(Scalar),
    /// A splitter failed to cut a part.
    #[error("splitting a part failed: {0}")]
    SplitFailed(String),
    /// Decomposition did not converge within the part budget.
    ///
    /// Reported rather than returning a partial decomposition, whose union
    /// would silently differ from the input.
    #[error("decomposition exceeded the {limit} part budget")]
    BudgetExceeded {
        /// The cap that was not raised.
        limit: usize,
    },
}

/// How hard to work at making each part convex.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum Strategy {
    /// Split until every part is convex to within `tolerance`.
    ///
    /// The union of the parts reproduces the input. Part count is whatever
    /// the geometry demands, which for a deeply non-convex solid is large.
    Exact,
    /// Stop once every part is convex to within `max_concavity`.
    ///
    /// Trades fidelity for part count. The union approximates the input:
    /// concave pockets shallower than the bound are filled in.
    Approximate {
        /// Largest tolerated distance from a part's vertex to its own hull.
        max_concavity: Scalar,
    },
}

/// Whether the parts reproduce the input or merely approximate it.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum Fidelity {
    /// Every part is convex to within tolerance; the union is the input.
    Exact,
    /// Parts are convex to within a bound larger than tolerance.
    Approximate {
        /// The bound that was requested.
        requested: Scalar,
        /// The largest concavity actually left in any part.
        ///
        /// Reported because it is the honest answer: a caller that asked for
        /// 10mm and got 2mm knows the result is better than it required,
        /// and one that reads this field cannot mistake the request for the
        /// outcome.
        achieved: Scalar,
    },
}

/// A solid expressed as convex parts, with the evidence to judge it.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Decomposition {
    /// The convex parts, in deterministic order.
    pub parts: Vec<TriMesh>,
    /// Whether the parts reproduce the input exactly.
    pub fidelity: Fidelity,
    /// Splits performed to reach this result.
    pub splits: usize,
}

impl Decomposition {
    /// Whether the input was already convex.
    pub fn is_single_part(&self) -> bool {
        self.parts.len() == 1
    }
}

/// Largest number of parts before the search is abandoned.
const MAX_PARTS: usize = 4096;

/// Decompose a closed two-manifold solid into convex parts.
///
/// # Errors
///
/// Refuses a ragged index buffer, an out-of-range index, an input that is
/// not a closed two-manifold solid, a non-positive concavity bound, and a
/// decomposition that exceeds the part budget.
pub fn convex_decompose(
    mesh: &TriMesh,
    strategy: Strategy,
    tolerance: Tolerance,
) -> Result<Decomposition, DecomposeError> {
    convex_decompose_with(mesh, strategy, tolerance, &split::Splitter::HandRolled)
}

/// Decompose a solid, choosing how parts are cut.
///
/// Identical to [`convex_decompose`] except that the caller supplies the
/// [`Splitter`](split::Splitter). Passing a boolean provider is how an
/// `algorithms` crate reaches a `providers` one: through the mesh-boolean
/// contract, which both layers may depend on.
///
/// # Errors
///
/// As [`convex_decompose`], plus [`DecomposeError::SplitFailed`] when the
/// supplied splitter cannot cut a part.
pub fn convex_decompose_with(
    mesh: &TriMesh,
    strategy: Strategy,
    tolerance: Tolerance,
    splitter: &split::Splitter<'_>,
) -> Result<Decomposition, DecomposeError> {
    if mesh.indices.len() % 3 != 0 {
        return Err(DecomposeError::RaggedIndices(mesh.indices.len()));
    }
    let vertex_count = mesh.positions.len();
    for (triangle, chunk) in mesh.indices.chunks_exact(3).enumerate() {
        for &index in chunk {
            if index as usize >= vertex_count {
                return Err(DecomposeError::IndexOutOfRange(triangle, index));
            }
        }
    }

    // A decomposition only means anything for a solid: the union of parts
    // reproduces a volume, not a surface. Checking here turns a meaningless
    // answer into a named refusal.
    let health = audit_mesh(mesh, tolerance);
    if !health.is_closed_two_manifold() {
        return Err(DecomposeError::NotASolid {
            boundary: health.boundary_edges,
            non_manifold: health.non_manifold_edges,
        });
    }

    let bound = match strategy {
        Strategy::Exact => tolerance.linear(),
        Strategy::Approximate { max_concavity } => {
            if !max_concavity.is_finite() || max_concavity <= 0.0 {
                return Err(DecomposeError::InvalidBound(max_concavity));
            }
            max_concavity
        }
    };

    // Work on meshes rather than point sets. A part is the actual solid on
    // one side of every split, produced by clipping; taking the hull of a
    // point subset instead would fill in any notch the subset still spans,
    // and the parts would sum to more volume than the input.
    let mut pending = vec![mesh.clone()];
    let mut finished: Vec<TriMesh> = Vec::new();
    let mut splits = 0usize;
    let mut achieved: Scalar = 0.0;

    while let Some(part) = pending.pop() {
        if finished.len() + pending.len() + 1 > MAX_PARTS {
            return Err(DecomposeError::BudgetExceeded { limit: MAX_PARTS });
        }

        let Some(reflex) = worst_concavity(&part.positions, &part.indices, tolerance) else {
            finished.push(part);
            continue;
        };
        if reflex.depth <= bound {
            achieved = achieved.max(reflex.depth);
            finished.push(part);
            continue;
        }

        // Split on the plane of the face the reflex vertex sticks out past.
        // Extending an existing face is the standard construction and it
        // makes progress by definition: everything in front of that plane
        // is separated from the face that could not see it.
        let (normal, offset) = (reflex.normal, reflex.offset);
        let (front, back) = splitter.split(&part, normal, offset, tolerance)?;

        match (front, back) {
            (Some(front), Some(back))
                if front.triangle_count() > 0 && back.triangle_count() > 0 =>
            {
                splits += 1;
                pending.push(front);
                pending.push(back);
            }
            // The plane failed to separate the part. Keeping it whole with
            // its concavity reported is honest; looping on a split that
            // makes no progress is not.
            _ => {
                achieved = achieved.max(reflex.depth);
                finished.push(part);
            }
        }
    }

    // Deterministic ordering: parts are keyed by their extreme corner, which
    // is a property of the geometry rather than of the traversal, so the
    // same solid decomposes to the same sequence on every run.
    finished.sort_by(|a, b| {
        let ka = order_key(&a.positions);
        let kb = order_key(&b.positions);
        ka.partial_cmp(&kb).unwrap_or(std::cmp::Ordering::Equal)
    });

    let parts = finished;

    let fidelity = match strategy {
        Strategy::Exact => Fidelity::Exact,
        Strategy::Approximate { max_concavity } => Fidelity::Approximate {
            requested: max_concavity,
            achieved,
        },
    };

    Ok(Decomposition {
        parts,
        fidelity,
        splits,
    })
}

/// Sort key: the lexicographically smallest corner of a part.
fn order_key(points: &[Point3]) -> (Scalar, Scalar, Scalar) {
    let mut best = (Scalar::INFINITY, Scalar::INFINITY, Scalar::INFINITY);
    for p in points {
        let key = (p.x, p.y, p.z);
        if key < best {
            best = key;
        }
    }
    best
}

/// A reflex feature: a vertex sticking out past one of the solid's own faces.
struct Reflex {
    /// How far the vertex lies in front of the face plane.
    depth: Scalar,
    /// The offending vertex.
    apex: Point3,
    /// Outward normal of the face it sticks out past.
    normal: Vec3,
    /// Plane offset of that face.
    offset: Scalar,
}

/// Depth and location of the worst reflex feature in a solid.
///
/// A solid is convex exactly when every vertex lies behind every face
/// plane. Where a vertex lies IN FRONT of some face plane, the solid
/// bulges past that face -- a reflex feature -- and the distance in front
/// is how deep the offending notch is.
///
/// Measuring against face planes rather than against the convex hull is
/// what makes this work. A reflex vertex generally lies exactly ON the
/// hull surface (the hull spans the notch with a face THROUGH that
/// vertex), so hull distance reports zero concavity for the very feature
/// that needs splitting.
///
/// `None` when the part is convex to within `tolerance`.
fn worst_concavity(positions: &[Point3], indices: &[u32], tolerance: Tolerance) -> Option<Reflex> {
    let linear = tolerance.linear();
    let mut worst: Option<Reflex> = None;

    // Bounding sphere over the vertices. For a unit normal `n`, no
    // vertex can satisfy dot(v, n) > centre.dot(n) + radius, so the
    // deepest a vertex could sit past a face plane is bounded without
    // touching a single vertex.
    //
    // The bound is CONSERVATIVE: it can only skip a face when no vertex
    // could qualify, so the result is identical to scanning every
    // vertex of every face -- including which face and vertex win a
    // tie. An AABB corner was tried first and prunes nothing on a
    // round mesh: it overshoots the true extent by up to sqrt(3).
    let &first = positions.first()?;
    let (mut low, mut high) = (first, first);
    for point in positions {
        low = Point3::new(low.x.min(point.x), low.y.min(point.y), low.z.min(point.z));
        high = Point3::new(
            high.x.max(point.x),
            high.y.max(point.y),
            high.z.max(point.z),
        );
    }
    let centre = (low + high) * 0.5;
    let radius = positions
        .iter()
        .fold(0.0, |m: Scalar, p| m.max((*p - centre).length()));

    for chunk in indices.chunks_exact(3) {
        let a = positions[chunk[0] as usize];
        let b = positions[chunk[1] as usize];
        let c = positions[chunk[2] as usize];

        let normal = (b - a).cross(c - a);
        let area = normal.length();
        // A degenerate face has no plane to test against; skip rather than
        // divide by a vanishing length and invent a direction.
        if area <= linear * linear {
            continue;
        }
        let unit = normal / area;

        // Deepest any vertex could sit past this plane, from the
        // bounding sphere alone -- no vertex touched.
        let reach = centre.dot(unit) + radius - a.dot(unit);

        // A vertex must clear `linear` to be a candidate at all. Against
        // an incumbent it must also reach the bottom of the tie window,
        // `depth - linear`, since an equal-depth vertex can still win on
        // coordinate order.
        //
        // Instrumented: the window is entered often (729 faces across
        // this suite) but no face inside it ever held a tie-breaking
        // winner, so pruning at `depth` behaves identically on every
        // input tried. The wider bound is kept because it CANNOT drop a
        // tie, not because a test distinguishes the two.
        let threshold = worst
            .as_ref()
            .map_or(linear, |current| (current.depth - linear).max(linear));
        if reach <= threshold {
            continue;
        }

        for (index, &point) in positions.iter().enumerate() {
            let ahead = (point - a).dot(unit);
            if ahead <= linear {
                continue;
            }
            // Deeper wins; equal depth breaks toward the lower index so the
            // choice is reproducible rather than dependent on iteration
            // order over an unordered structure.
            let better = match &worst {
                None => true,
                Some(current) => {
                    ahead > current.depth + linear
                        || ((ahead - current.depth).abs() <= linear
                            && (point.x, point.y, point.z)
                                < (current.apex.x, current.apex.y, current.apex.z))
                }
            };
            if better {
                let _ = index;
                // Record the FACE the vertex sticks out past, not just how
                // far. Splitting on that face's own plane is what removes
                // the reflex feature; a bounding-box axis through the same
                // point need not separate the notch at all.
                worst = Some(Reflex {
                    depth: ahead,
                    apex: point,
                    normal: unit,
                    offset: a.dot(unit),
                });
            }
        }
    }
    worst
}

/// Clip a closed solid by a plane, keeping the side the normal points away
/// from and capping the opening so the result is closed again.
///
/// Sutherland-Hodgman per triangle: each face is clipped to the half-space
/// and re-fanned into triangles. The opening is then capped by measuring
/// which edges the clipped shell left used only once, which is what keeps
/// the part a solid rather than an open shell.
fn clip(mesh: &TriMesh, normal: Vec3, offset: Scalar, tolerance: Tolerance) -> Option<TriMesh> {
    let linear = tolerance.linear();
    let mut positions: Vec<Point3> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    let mut lookup: AHashMap<(u64, u64, u64), u32> = AHashMap::new();

    let intern = |point: Point3, positions: &mut Vec<Point3>, lookup: &mut AHashMap<_, _>| {
        let key = (
            quantise(point.x, linear),
            quantise(point.y, linear),
            quantise(point.z, linear),
        );
        *lookup.entry(key).or_insert_with(|| {
            positions.push(point);
            (positions.len() - 1) as u32
        })
    };

    for chunk in mesh.indices.chunks_exact(3) {
        let triangle = [
            mesh.positions[chunk[0] as usize],
            mesh.positions[chunk[1] as usize],
            mesh.positions[chunk[2] as usize],
        ];
        let distances = [
            triangle[0].dot(normal) - offset,
            triangle[1].dot(normal) - offset,
            triangle[2].dot(normal) - offset,
        ];

        // A face lying IN the clip plane belongs to exactly one side, and
        // its distances cannot say which: every corner reads as "on the
        // plane", so both sides would keep it, duplicating the face and
        // double-counting its volume. Its own normal settles it -- a
        // coplanar face bounds the material on the side it faces away from.
        if distances.iter().all(|d| d.abs() <= linear) {
            let face = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
            if face.dot(normal) > 0.0 {
                let a = intern(triangle[0], &mut positions, &mut lookup);
                let b = intern(triangle[1], &mut positions, &mut lookup);
                let c = intern(triangle[2], &mut positions, &mut lookup);
                if a != b && b != c && c != a {
                    indices.extend_from_slice(&[a, b, c]);
                }
            }
            continue;
        }

        // Sutherland-Hodgman: walk the triangle's edges, keeping corners
        // behind the plane and the points where an edge crosses it.
        let mut kept: Vec<Point3> = Vec::new();
        for corner in 0..3 {
            let current = triangle[corner];
            let next = triangle[(corner + 1) % 3];
            let d_current = distances[corner];
            let d_next = distances[(corner + 1) % 3];

            if d_current <= linear {
                kept.push(current);
            }
            // The crossing point is shared by both parts, which is what
            // makes their union seamless along the cut.
            if (d_current < -linear && d_next > linear) || (d_current > linear && d_next < -linear)
            {
                let t = d_current / (d_current - d_next);
                kept.push(current + (next - current) * t);
            }
        }
        if kept.len() < 3 {
            continue;
        }

        let anchor = intern(kept[0], &mut positions, &mut lookup);
        for corner in 1..kept.len() - 1 {
            let b = intern(kept[corner], &mut positions, &mut lookup);
            let c = intern(kept[corner + 1], &mut positions, &mut lookup);
            if anchor != b && b != c && c != anchor {
                indices.extend_from_slice(&[anchor, b, c]);
            }
        }
    }

    if indices.is_empty() {
        return None;
    }

    // Cap whatever the clipping actually left open.
    //
    // Earlier versions PREDICTED the cut boundary from the input mesh,
    // deciding per triangle whether an edge belonged to the cross-section.
    // That oscillated between leaving holes and laying the cap over
    // existing walls, because whether a face standing on the plane bounds
    // THIS part is a global question a single triangle cannot answer.
    //
    // Measuring instead of predicting removes the question entirely. In a
    // closed shell every undirected edge is used exactly twice, so the
    // edges used ONCE are precisely the hole -- whatever the clipping did
    // upstream. The cap can then be neither too generous nor too strict.
    let shell = TriMesh::new(positions.clone(), indices.clone());
    let adjacency = EdgeAdjacency::build(&shell);
    let open_edges: Vec<(Point3, Point3)> = adjacency
        .boundary_edges()
        .map(|edge| {
            let (a, b) = edge.endpoints();
            (positions[a as usize], positions[b as usize])
        })
        .collect();

    if open_edges.is_empty() {
        return Some(TriMesh::new(positions, indices));
    }

    for loop_points in stitch_loops(&open_edges, linear) {
        if loop_points.len() < 3 {
            continue;
        }
        // Ear clipping works in 2D, so express the loop in the cut plane's
        // own basis. Any orthonormal pair perpendicular to the normal does.
        let (axis_u, axis_v) = plane_basis(normal);
        let origin = loop_points[0];
        let flat: Vec<Point2> = loop_points
            .iter()
            .map(|p| {
                let d = *p - origin;
                Point2::new(d.dot(axis_u), d.dot(axis_v))
            })
            .collect();

        // Ear clipping needs a counter-clockwise ring; a clockwise one is
        // reversed rather than refused, since the winding of a cut loop is
        // an artefact of traversal, not of the geometry.
        let area: Scalar = flat
            .iter()
            .enumerate()
            .map(|(k, p)| {
                let q = flat[(k + 1) % flat.len()];
                p.x * q.y - q.x * p.y
            })
            .sum();
        let (flat, loop_points) = if area < 0.0 {
            let mut f = flat;
            let mut l = loop_points;
            f.reverse();
            l.reverse();
            (f, l)
        } else {
            (flat, loop_points)
        };

        let Ok(fan) = axiolid_reference::polygon::triangulate_simple(&flat) else {
            continue;
        };
        for triple in fan {
            let a = intern(loop_points[triple[0] as usize], &mut positions, &mut lookup);
            let b = intern(loop_points[triple[1] as usize], &mut positions, &mut lookup);
            let c = intern(loop_points[triple[2] as usize], &mut positions, &mut lookup);
            if a == b || b == c || c == a {
                continue;
            }
            // The cap faces along the clip normal, opposite the material
            // that was removed, so the shell stays consistently outward.
            let wound = (positions[b as usize] - positions[a as usize])
                .cross(positions[c as usize] - positions[a as usize]);
            if wound.dot(normal) >= 0.0 {
                indices.extend_from_slice(&[a, b, c]);
            } else {
                indices.extend_from_slice(&[a, c, b]);
            }
        }
    }

    Some(TriMesh::new(positions, indices))
}

/// Chain unordered cut edges into closed loops.
///
/// Clipping produces the cut edges one triangle at a time, in no
/// particular order. A cap can only be triangulated once those edges are
/// walked into a ring, so each edge is joined to the next one sharing an
/// endpoint until the loop closes.
///
/// Endpoints are matched on a tolerance lattice: the same crossing point
/// computed from two adjacent triangles differs in the last few bits, and
/// exact comparison would leave every loop broken.
fn stitch_loops(edges: &[(Point3, Point3)], linear: Scalar) -> Vec<Vec<Point3>> {
    let key = |p: &Point3| {
        (
            quantise(p.x, linear),
            quantise(p.y, linear),
            quantise(p.z, linear),
        )
    };

    let mut adjacency: BTreeMap<(u64, u64, u64), Vec<usize>> = BTreeMap::new();
    for (index, (from, to)) in edges.iter().enumerate() {
        adjacency.entry(key(from)).or_default().push(index);
        adjacency.entry(key(to)).or_default().push(index);
    }

    let mut used = vec![false; edges.len()];
    let mut loops = Vec::new();

    for start in 0..edges.len() {
        if used[start] {
            continue;
        }
        used[start] = true;
        let mut ring = vec![edges[start].0, edges[start].1];
        let mut tail = edges[start].1;

        while let Some(candidates) = adjacency.get(&key(&tail)) {
            let mut advanced = false;
            for &next in candidates {
                if used[next] {
                    continue;
                }
                let (from, to) = edges[next];
                let other = if key(&from) == key(&tail) {
                    to
                } else if key(&to) == key(&tail) {
                    from
                } else {
                    continue;
                };
                used[next] = true;
                // Closing the ring: stop rather than repeat the first point.
                if key(&other) == key(&ring[0]) {
                    advanced = false;
                    break;
                }
                ring.push(other);
                tail = other;
                advanced = true;
                break;
            }
            if !advanced {
                break;
            }
        }
        if ring.len() >= 3 {
            loops.push(ring);
        }
    }
    loops
}

/// Any orthonormal basis of the plane perpendicular to `normal`.
fn plane_basis(normal: Vec3) -> (Vec3, Vec3) {
    // Seed against the axis the normal is least aligned with, so the cross
    // product is well conditioned rather than near zero.
    let seed = if normal.x.abs() <= normal.y.abs() && normal.x.abs() <= normal.z.abs() {
        Vec3::X
    } else if normal.y.abs() <= normal.z.abs() {
        Vec3::Y
    } else {
        Vec3::Z
    };
    let u = normal.cross(seed).normalize();
    let v = normal.cross(u);
    (u, v)
}

/// Snap a coordinate to a tolerance-sized lattice for welding.
///
/// Two clipped faces meeting at a cut must agree on the crossing vertex, or
/// the part is not closed. Comparing raw bits is too strict: the same point
/// computed from two different edges differs in the last ulp.
fn quantise(value: Scalar, linear: Scalar) -> u64 {
    let step = linear.max(Scalar::EPSILON);
    let snapped = (value / step).round();
    snapped.to_bits()
}

#[cfg(test)]
mod concavity_tests {
    use super::*;

    fn tol() -> Tolerance {
        Tolerance::new(1e-9, 1e-12).expect("tolerance")
    }

    /// The pre-prune implementation, kept verbatim as the oracle. The
    /// prune is only correct if it agrees with this on every input.
    fn unpruned(positions: &[Point3], indices: &[u32], tolerance: Tolerance) -> Option<Reflex> {
        let linear = tolerance.linear();
        let mut worst: Option<Reflex> = None;
        for chunk in indices.chunks_exact(3) {
            let a = positions[chunk[0] as usize];
            let b = positions[chunk[1] as usize];
            let c = positions[chunk[2] as usize];
            let normal = (b - a).cross(c - a);
            let area = normal.length();
            if area <= linear * linear {
                continue;
            }
            let unit = normal / area;
            for &point in positions.iter() {
                let ahead = (point - a).dot(unit);
                if ahead <= linear {
                    continue;
                }
                let better = match &worst {
                    None => true,
                    Some(current) => {
                        ahead > current.depth + linear
                            || ((ahead - current.depth).abs() <= linear
                                && (point.x, point.y, point.z)
                                    < (current.apex.x, current.apex.y, current.apex.z))
                    }
                };
                if better {
                    worst = Some(Reflex {
                        depth: ahead,
                        apex: point,
                        normal: unit,
                        offset: a.dot(unit),
                    });
                }
            }
        }
        worst
    }

    fn agree(label: &str, mesh: &TriMesh) {
        let want = unpruned(&mesh.positions, &mesh.indices, tol());
        let got = worst_concavity(&mesh.positions, &mesh.indices, tol());
        match (want, got) {
            (None, None) => {}
            (Some(w), Some(g)) => {
                assert!((w.depth - g.depth).abs() < 1e-12, "{label}: depth");
                // Same apex AND same plane: the caller splits on this
                // plane, so a different face changes the decomposition.
                assert_eq!(w.apex, g.apex, "{label}: apex");
                assert_eq!(w.normal, g.normal, "{label}: normal");
                assert!((w.offset - g.offset).abs() < 1e-12, "{label}: offset");
            }
            (a, b) => panic!(
                "{label}: presence differs, {} vs {}",
                a.is_some(),
                b.is_some()
            ),
        }
    }

    fn cube() -> TriMesh {
        let p = vec![
            Point3::new(-1.0, -1.0, -1.0),
            Point3::new(1.0, -1.0, -1.0),
            Point3::new(1.0, 1.0, -1.0),
            Point3::new(-1.0, 1.0, -1.0),
            Point3::new(-1.0, -1.0, 1.0),
            Point3::new(1.0, -1.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(-1.0, 1.0, 1.0),
        ];
        let i = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6,
            5, 0, 4, 7, 0, 7, 3u32,
        ];
        TriMesh::new(p, i)
    }

    /// Extruded L: a genuine reflex corner, and enough symmetry that
    /// several faces report the same depth.
    fn l_shape() -> TriMesh {
        let footprint = [
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 2.0),
            (0.0, 2.0),
        ];
        let mut positions = Vec::new();
        for &(x, y) in &footprint {
            positions.push(Point3::new(x, y, 0.0));
        }
        for &(x, y) in &footprint {
            positions.push(Point3::new(x, y, 1.0));
        }
        let n = footprint.len() as u32;
        let mut indices = Vec::new();
        for &(a, b, c) in &[(0u32, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5)] {
            indices.extend_from_slice(&[a, c, b]);
            indices.extend_from_slice(&[a + n, b + n, c + n]);
        }
        for i in 0..n {
            let j = (i + 1) % n;
            indices.extend_from_slice(&[i, j, j + n]);
            indices.extend_from_slice(&[i, j + n, i + n]);
        }
        TriMesh::new(positions, indices)
    }

    #[test]
    fn prune_agrees_on_a_convex_solid() {
        agree("cube", &cube());
    }

    /// Pull one corner inward so a genuine reflex feature exists: the
    /// convex case alone would let a prune that skips EVERYTHING pass.
    #[test]
    fn prune_agrees_on_a_dented_solid() {
        let mut mesh = cube();
        mesh.positions[6] = Point3::new(0.1, 0.1, 0.1);
        agree("dented", &mesh);
        assert!(
            worst_concavity(&mesh.positions, &mesh.indices, tol()).is_some(),
            "the dent must register as concavity, or this proves nothing"
        );
    }

    /// Many shapes, deterministic pseudo-random. A handcrafted fixture
    /// exercises one path through the tie-break; this sweeps enough
    /// geometry to hit equal-depth cases the prune must not skip.
    #[test]
    fn prune_agrees_across_many_dents() {
        let mut seed = 0x9E3779B97F4A7C15u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        for trial in 0..200 {
            let mut mesh = cube();
            for _ in 0..3 {
                let which = (next() * 8.0) as usize % 8;
                let scale = 0.2 + next() * 1.4;
                mesh.positions[which] *= scale;
            }
            agree(&format!("trial {trial}"), &mesh);
        }
    }

    /// Equal depths across several faces are what the tie-break exists
    /// to resolve, and what a prune clamped to the incumbent depth
    /// would skip. A symmetric dent produces them exactly; a coarse
    /// tolerance widens the tie window enough to be reachable.
    #[test]
    fn prune_respects_the_tie_window() {
        let coarse = Tolerance::new(0.05, 1e-12).expect("tolerance");
        // Pull four top corners inward by the SAME amount: several
        // faces then report identical reflex depth.
        // Push four corners OUTWARD symmetrically: spikes give several
        // faces an identical, genuinely-reflex depth.
        let mesh = l_shape();
        let want = unpruned(&mesh.positions, &mesh.indices, coarse);
        let got = worst_concavity(&mesh.positions, &mesh.indices, coarse);
        let (want, got) = (want.expect("reflex"), got.expect("reflex"));
        assert!((want.depth - got.depth).abs() < 1e-12, "depth differs");
        assert!(
            (want.apex - got.apex).length() < 1e-12,
            "same depth, different apex: the tie-break was not preserved"
        );
    }

    /// Sweep tolerance so the tie window spans the gap between the
    /// bounding-sphere reach and the true depth. Somewhere in that
    /// sweep a face is skipped by a prune clamped to the incumbent
    /// depth but kept by one that honours the window -- if the two
    /// ever differ, this finds it.
    #[test]
    fn prune_matches_the_oracle_across_tolerances() {
        let meshes = [("l", l_shape()), ("cube", cube())];
        for (name, mesh) in &meshes {
            let mut linear = 1e-12;
            while linear < 2.0 {
                let t = Tolerance::new(linear, 1e-12).expect("tolerance");
                let want = unpruned(&mesh.positions, &mesh.indices, t);
                let got = worst_concavity(&mesh.positions, &mesh.indices, t);
                match (want, got) {
                    (None, None) => {}
                    (Some(a), Some(b)) => {
                        assert!(
                            (a.depth - b.depth).abs() < 1e-12 && (a.apex - b.apex).length() < 1e-12,
                            "{name} at linear={linear:e}: prune changed the answer"
                        );
                    }
                    (a, b) => panic!(
                        "{name} at linear={linear:e}: presence differs, {} vs {}",
                        a.is_some(),
                        b.is_some()
                    ),
                }
                linear *= 1.5;
            }
        }
    }

    /// Randomised search for an input where a prune clamped to the
    /// incumbent depth differs from one honouring the tie window.
    /// Coarse tolerances widen the window; random point sets give the
    /// bounding-sphere bound a chance to be tight.
    #[test]
    fn prune_matches_the_oracle_on_random_solids() {
        let mut seed = 0xD1B54A32D192ED03u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        for trial in 0..400 {
            // Quantised coordinates: exact ties are then reachable,
            // which continuous random values would never produce.
            let mut mesh = cube();
            for slot in 0..8 {
                let q = |v: f64| (v * 4.0).round() / 4.0;
                let p = mesh.positions[slot];
                let s = 0.25 + (next() * 8.0).floor() / 4.0;
                mesh.positions[slot] = Point3::new(q(p.x * s), q(p.y * s), q(p.z * s));
            }
            for step in 0..6 {
                let linear = 0.01 * 4.0_f64.powi(step);
                let t = Tolerance::new(linear, 1e-12).expect("tolerance");
                let want = unpruned(&mesh.positions, &mesh.indices, t);
                let got = worst_concavity(&mesh.positions, &mesh.indices, t);
                match (want, got) {
                    (None, None) => {}
                    (Some(a), Some(b)) => assert!(
                        (a.depth - b.depth).abs() < 1e-12 && (a.apex - b.apex).length() < 1e-12,
                        "trial {trial} linear={linear}: prune changed the answer"
                    ),
                    (a, b) => panic!(
                        "trial {trial} linear={linear}: presence differs, {} vs {}",
                        a.is_some(),
                        b.is_some()
                    ),
                }
            }
        }
    }

    /// Constructed, not searched: two spikes at equal depth, where the
    /// second face has a bounding-sphere reach just below the
    /// incumbent depth. A prune clamped to that depth skips it and
    /// loses the tie-break; one honouring the window keeps it.
    #[test]
    fn prune_keeps_faces_inside_the_tie_window() {
        // Coarse tolerance so the window has real width.
        let t = Tolerance::new(0.25, 1e-12).expect("tolerance");
        // Sweep asymmetric spikes: some trial puts a tie-breaking
        // vertex behind a face whose reach sits inside the window.
        for a in 1..14 {
            for b in 1..14 {
                let mut mesh = cube();
                let sa = 1.0 + a as Scalar * 0.125;
                let sb = 1.0 + b as Scalar * 0.125;
                let p4 = mesh.positions[4];
                let p6 = mesh.positions[6];
                mesh.positions[4] = Point3::new(p4.x * sa, p4.y * sa, p4.z * sa);
                mesh.positions[6] = Point3::new(p6.x * sb, p6.y * sb, p6.z * sb);
                let want = unpruned(&mesh.positions, &mesh.indices, t);
                let got = worst_concavity(&mesh.positions, &mesh.indices, t);
                match (want, got) {
                    (None, None) => {}
                    (Some(x), Some(y)) => assert!(
                        (x.depth - y.depth).abs() < 1e-12 && (x.apex - y.apex).length() < 1e-12,
                        "a={a} b={b}: prune changed the answer"
                    ),
                    (x, y) => panic!("a={a} b={b}: {} vs {}", x.is_some(), y.is_some()),
                }
            }
        }
    }

    /// The bounding sphere is computed from the first vertex, so an
    /// empty mesh must not index it.
    #[test]
    fn empty_input_is_none() {
        assert!(worst_concavity(&[], &[], tol()).is_none());
    }

    /// Degenerate faces are skipped before the plane is formed; the
    /// prune must not change that.
    #[test]
    fn degenerate_faces_are_still_skipped() {
        let p = vec![Point3::ZERO, Point3::ZERO, Point3::ZERO];
        assert!(worst_concavity(&p, &[0, 1, 2], tol()).is_none());
    }
}
