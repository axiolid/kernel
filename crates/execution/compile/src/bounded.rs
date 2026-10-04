//! Curve-bounded planes (#192): a plane trimmed by closed boundary curves,
//! as IFC `IfcCurveBoundedPlane` space-boundary connection surfaces are.
//!
//! The boundaries live in the plane's parameter space: a point `(u, v)` of
//! a boundary is `origin + u x + v y` in the plane's frame. The first
//! boundary is the outer loop, the rest are holes. Straight boundaries
//! (lines, polylines) are kept exactly; curved ones are chorded within the
//! chord budget. The region is triangulated in its plane by the same path
//! authored polygon faces use, which refuses crossing rings, holes outside
//! the outer loop and loops with no area rather than filling them. The
//! result is a surface, never a solid.
//!
//! A boundary may also be a curve relation, a composite or a trim (#255):
//! it is resolved to points by the sweep directrix reader (see
//! [`relation`]) and must close, its ends no farther apart than the linear
//! tolerance, as composite joints may be. An open one is refused by name.

use axiolid_contracts::{BackendId, ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Frame3, Point2, Point3, Scalar, Tolerance};
use axiolid_curve::{Curve2, Curve3};
use axiolid_mesh::TriMesh;
use axiolid_model::{GeometryGraph, GeometryNode, NodeId};

pub(crate) mod relation;

/// A runaway guard on chord subdivision, as for profiles.
const MAX_SUBDIVISION_DEPTH: u32 = 24;

/// Triangulate a curve-bounded plane.
///
/// Also returns the largest distance a merged near-duplicate boundary point
/// lay from the point kept in its place, in the plane's parameters (#232):
/// merging moves a ring by no more than that.
pub(crate) fn curve_bounded(
    backend: BackendId,
    graph: &GeometryGraph,
    basis: NodeId,
    boundaries: &[NodeId],
    implicit_outer: bool,
    options: &ExecutionOptions,
) -> GeomResult<(TriMesh, Scalar)> {
    let tolerance = options.tolerance();
    if implicit_outer {
        // The outer loop would be the basis surface's own boundary, which a
        // plane does not have.
        return Err(GeomError::UnsupportedInput {
            backend,
            operation: Operation::Tessellation,
            input: "curve-bounded plane with an implicit outer boundary",
        });
    }
    let frame = match graph.get(basis) {
        Some(GeometryNode::Surface(axiolid_surface::Surface::Plane(plane))) => plane.frame,
        Some(GeometryNode::Surface(_)) => {
            return Err(GeomError::UnsupportedInput {
                backend,
                operation: Operation::SurfaceEvaluation,
                input: "curve-bounded surface over a non-planar basis",
            })
        }
        Some(_) => {
            return Err(GeomError::InvalidInput(format!(
                "curve-bounded basis {basis:?} is not a Surface node"
            )))
        }
        None => {
            return Err(GeomError::InvalidInput(format!(
                "curve-bounded basis {basis:?} does not belong to this graph"
            )))
        }
    };
    if boundaries.is_empty() {
        return Err(GeomError::InvalidInput(
            "curve-bounded plane has no boundary".to_owned(),
        ));
    }
    let local = boundary_options(&frame, options);
    let chord = crate::compiler::chord_error(&local);
    let mut rings: Vec<Vec<Point2>> = Vec::with_capacity(boundaries.len());
    let mut merged: Scalar = 0.0;
    for (index, &id) in boundaries.iter().enumerate() {
        let mut ring = boundary_ring(graph, id, &local, chord, tolerance, &mut merged)?;
        if ring.len() < 3 {
            return Err(GeomError::Degenerate(format!(
                "curve-bounded boundary {index} has {} distinct points, need at least 3",
                ring.len()
            )));
        }
        // The outer loop counter-clockwise in the plane's parameters, so the
        // triangles face along the plane's normal. Reversed about its first
        // point, so a loop read backwards (a composite segment against its
        // curve's sense, #255) triangulates as the loop read forwards.
        if index == 0 && signed_area(&ring) < 0.0 {
            ring[1..].reverse();
        }
        rings.push(ring);
    }
    let map = |p: Point2| frame.origin + frame.x * p.x + frame.y * p.y;
    let points: Vec<Vec<Point3>> = rings
        .iter()
        .map(|ring| ring.iter().map(|&p| map(p)).collect())
        .collect();
    let views: Vec<&[Point3]> = points.iter().map(Vec::as_slice).collect();
    let indices =
        crate::planar::triangulate_polygon(&views, tolerance.linear()).map_err(|refusal| {
            GeomError::InvalidInput(format!(
                "curve-bounded plane cannot be triangulated: {refusal}"
            ))
        })?;
    let positions: Vec<Point3> = points.into_iter().flatten().collect();
    let indices: Vec<u32> = indices
        .into_iter()
        .map(|i| {
            u32::try_from(i).map_err(|_| GeomError::BudgetExceeded {
                resource: "curve-bounded plane indices",
            })
        })
        .collect::<GeomResult<_>>()?;
    let mesh = TriMesh::new(positions, indices);
    mesh.validate_structure().map_err(|error| {
        GeomError::InvalidInput(format!("invalid curve-bounded plane mesh: {error}"))
    })?;
    Ok((mesh, merged))
}

/// The options boundaries are flattened with: the chord budget shrunk by
/// the plane frame's stretch (the boundaries live in its parameters), so
/// the deviation report's bound, that budget times the stretch, stays
/// within the budget asked for (#255). Unchanged when the stretch is not a
/// finite factor of at least one, or the shrunk budget is not positive.
pub(crate) fn boundary_options(frame: &Frame3, options: &ExecutionOptions) -> ExecutionOptions {
    let chord = crate::compiler::chord_error(options);
    let stretch = axiolid_reference::bound::frame_stretch3(frame);
    if !(stretch.is_finite() && stretch >= 1.0) {
        return options.clone();
    }
    let mut budget = chord / stretch;
    while budget > 0.0 && budget * stretch > chord {
        budget = budget.next_down();
    }
    options
        .clone()
        .with_chord_error(budget)
        .unwrap_or_else(|| options.clone())
}

/// A closed boundary as a ring in the plane's parameters.
fn boundary_ring(
    graph: &GeometryGraph,
    id: NodeId,
    options: &ExecutionOptions,
    chord: Scalar,
    tolerance: Tolerance,
    merged: &mut Scalar,
) -> GeomResult<Vec<Point2>> {
    let linear = tolerance.linear();
    let open = |closed: bool, first: Option<&Point2>, last: Option<&Point2>| {
        let meets = match (first, last) {
            (Some(a), Some(b)) => (a.x - b.x).abs() <= linear && (a.y - b.y).abs() <= linear,
            _ => false,
        };
        !closed && !meets
    };
    let raw: Vec<Point2> = match graph.get(id) {
        Some(GeometryNode::Curve2(Curve2::Polyline(p))) => {
            if open(p.closed, p.points.first(), p.points.last()) {
                return Err(GeomError::InvalidInput(format!(
                    "curve-bounded boundary {id:?} is an open polyline"
                )));
            }
            p.points.clone()
        }
        Some(GeometryNode::Curve2(curve)) => {
            let domain = axiolid_reference::curve::domain2(curve);
            axiolid_reference::curve::flatten2(curve, domain, chord, MAX_SUBDIVISION_DEPTH)?
        }
        Some(GeometryNode::CurveRelation(_)) => {
            let points = relation::points(graph, id, options)?;
            let gap = points[0].distance(points[points.len() - 1]);
            if gap.is_nan() || gap > linear {
                return Err(GeomError::InvalidInput(format!(
                    "curve-bounded boundary {id:?} is an open curve relation: its ends are \
                     {gap} apart, more than the linear tolerance {linear}"
                )));
            }
            in_parameter_plane(id, &points, linear)?
        }
        Some(GeometryNode::Curve3(curve)) => {
            let points = match curve {
                Curve3::Polyline(p) => {
                    let flat = |q: &Point3| Point2::new(q.x, q.y);
                    if open(
                        p.closed,
                        p.points.first().map(flat).as_ref(),
                        p.points.last().map(flat).as_ref(),
                    ) {
                        return Err(GeomError::InvalidInput(format!(
                            "curve-bounded boundary {id:?} is an open polyline"
                        )));
                    }
                    p.points.clone()
                }
                other => {
                    let domain = axiolid_reference::curve::domain3(other);
                    axiolid_reference::curve::flatten3(other, domain, chord, MAX_SUBDIVISION_DEPTH)?
                }
            };
            in_parameter_plane(id, &points, linear)?
        }
        Some(_) => {
            return Err(GeomError::InvalidInput(format!(
                "curve-bounded boundary {id:?} is not a curve node"
            )))
        }
        None => {
            return Err(GeomError::InvalidInput(format!(
                "curve-bounded boundary {id:?} does not belong to this graph"
            )))
        }
    };
    if raw.iter().any(|p| !p.is_finite()) {
        return Err(GeomError::InvalidInput(format!(
            "curve-bounded boundary {id:?} has a non-finite point"
        )));
    }
    // Drop repeated points and the closing repeat of the first. A repeat
    // within the rounding of the points' own evaluation (a full turn ends
    // `sin(2 pi) r` off its start) is the same point, welded exactly; only
    // a wider one moves the ring and is counted (#255).
    let weld = 2.0 * relation::JOINT_ROUNDING * relation::boundary_scale(graph, id, options);
    let mut record = |distance: Scalar| {
        if distance > weld {
            *merged = merged.max(distance);
        }
    };
    let near = |a: Point2, b: Point2| (a.x - b.x).abs() <= linear && (a.y - b.y).abs() <= linear;
    let mut ring: Vec<Point2> = Vec::with_capacity(raw.len());
    for p in raw {
        match ring.last() {
            Some(&q) if near(q, p) => record((q - p).length()),
            _ => ring.push(p),
        }
    }
    while ring.len() > 1 && near(ring[0], ring[ring.len() - 1]) {
        let dropped = ring.pop().expect("non-empty");
        record((ring[0] - dropped).length());
    }
    Ok(ring)
}

/// A 3D boundary is still in the plane's parameters: its points must lie
/// in the parameter plane, z = 0.
fn in_parameter_plane(id: NodeId, points: &[Point3], linear: Scalar) -> GeomResult<Vec<Point2>> {
    if let Some(off) = points.iter().find(|p| p.z.abs() > linear) {
        return Err(GeomError::InvalidInput(format!(
            "curve-bounded boundary {id:?} leaves the parameter plane (z = {})",
            off.z
        )));
    }
    Ok(points.iter().map(|p| Point2::new(p.x, p.y)).collect())
}

fn signed_area(ring: &[Point2]) -> Scalar {
    let n = ring.len();
    (0..n)
        .map(|i| ring[i].x * ring[(i + 1) % n].y - ring[(i + 1) % n].x * ring[i].y)
        .sum::<Scalar>()
        * 0.5
}
