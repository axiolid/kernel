//! Exact polyhedral lofts between planar polygon sections (#247).
//!
//! [`loft_polygon_stations_exact`] joins planar polygon sections placed in
//! their own frames, vertex `k` of every ring to vertex `k` of the same
//! ring at the next section, into a closed B-rep of planar faces: a cap at
//! the first and at the last section (holes included) and one wall per
//! ring edge between consecutive sections. It is the exact form of a
//! linear taper along a straight line, such as a station-placed spine whose
//! rectangle widens from one station to the next.
//!
//! Each wall is the ruled surface between two section edges. It is a
//! plane exactly when its four corners are coplanar; a wall whose fourth
//! corner leaves the plane of the other three by more than the linear
//! tolerance is a bilinear patch, which this builder does not make, and is
//! refused by name ([`WARPED_WALL`]), never flattened.
//!
//! The conventions are the exact extrusion's: every section's outer
//! ring runs counter-clockwise and each hole clockwise about its frame's
//! `z`, which points from one section to the next; the first cap faces
//! back along it, the last forward, and each wall's plane normal is the
//! right-hand normal of its loop `a -> b -> d -> c` (`a, b` an edge of one
//! section, `c, d` the same edge of the next), which points out of the
//! material.

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Frame3, Interval, Point2, Point3, Tolerance, Vec2, Vec3};
use axiolid_curve::{Curve2, Line2};
use axiolid_surface::{Plane, Surface};
use axiolid_topology::{EdgeId, EdgeUse, Face, FaceBound, Orientation, VertexId};

use crate::extrude_exact::{add_line_edge, add_loop, finish_closed, reserve};
use crate::BACKEND_ID;

/// The refusal of a wall whose corners are not coplanar.
pub const WARPED_WALL: &str =
    "an exact loft wall whose four corners are not coplanar: the ruled face between the two \
     section edges is a bilinear patch";

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BACKEND_ID,
        operation: Operation::Sweep,
        input,
    }
}

fn place(frame: &Frame3, p: Point2) -> Point3 {
    frame.origin + frame.x * p.x + frame.y * p.y
}

/// Loft planar polygon sections exactly; see the
/// [module documentation](self).
///
/// Each section is its plane's frame (orthonormal, `z` towards the next
/// section) and its rings in that frame: the outer first, counter-clockwise
/// about `z`, then each hole, clockwise. Every section has the same ring
/// structure.
///
/// # Errors
///
/// Fewer than two sections, sections whose ring structures differ, a ring
/// of fewer than three vertices, and a collapsed edge, by name; a wall
/// whose corners are not coplanar within `tolerance` as
/// [`GeomError::UnsupportedInput`] naming [`WARPED_WALL`]; and an assembly
/// that is not a closed manifold.
pub fn loft_polygon_stations_exact(
    sections: &[(Frame3, Vec<Vec<Point2>>)],
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    if sections.len() < 2 {
        return Err(GeomError::InvalidInput(format!(
            "an exact loft needs at least two sections, got {}",
            sections.len()
        )));
    }
    let shape: Vec<usize> = sections[0].1.iter().map(Vec::len).collect();
    if shape.is_empty() || shape.iter().any(|&n| n < 3) {
        return Err(GeomError::Degenerate(
            "an exact loft section needs rings of at least three vertices".into(),
        ));
    }
    if sections
        .iter()
        .any(|(_, rings)| rings.iter().map(Vec::len).ne(shape.iter().copied()))
    {
        return Err(GeomError::InvalidInput(
            "exact loft sections must share their ring structure".into(),
        ));
    }
    let count = sections.len();
    let per: usize = shape.iter().sum();
    let walls = (count - 1) * per;
    let mut builder = ExactBRepBuilder::default();
    reserve(
        &mut builder,
        count * per,
        count * per + walls,
        2 * shape.len() + walls,
        2 + walls,
        count * per + walls,
        2 * per + 4 * walls,
        2 + walls,
    )?;
    let points: Vec<Vec<Vec<Point3>>> = sections
        .iter()
        .map(|(frame, rings)| {
            rings
                .iter()
                .map(|ring| ring.iter().map(|p| place(frame, *p)).collect())
                .collect()
        })
        .collect();
    let vertices: Vec<Vec<Vec<VertexId>>> = points
        .iter()
        .map(|rings| {
            rings
                .iter()
                .map(|ring| {
                    ring.iter()
                        .map(|&position| {
                            builder
                                .topology_mut()
                                .add_vertex(axiolid_topology::Vertex { position })
                        })
                        .collect()
                })
                .collect()
        })
        .collect();
    let edge = |builder: &mut ExactBRepBuilder,
                (from, to): (VertexId, VertexId),
                (a, b): (Point3, Point3)|
     -> GeomResult<EdgeId> {
        if (b - a).length() == 0.0 {
            return Err(GeomError::Degenerate(
                "an exact loft has a collapsed edge".into(),
            ));
        }
        Ok(add_line_edge(builder, from, to, a, b - a))
    };
    // Each section's ring edges, and the edges joining consecutive sections.
    let mut rims: Vec<Vec<Vec<EdgeId>>> = Vec::with_capacity(count);
    for j in 0..count {
        let mut rings = Vec::with_capacity(shape.len());
        for (r, &n) in shape.iter().enumerate() {
            let mut ring = Vec::with_capacity(n);
            for k in 0..n {
                let m = (k + 1) % n;
                ring.push(edge(
                    &mut builder,
                    (vertices[j][r][k], vertices[j][r][m]),
                    (points[j][r][k], points[j][r][m]),
                )?);
            }
            rings.push(ring);
        }
        rims.push(rings);
    }
    let mut rails: Vec<Vec<Vec<EdgeId>>> = Vec::with_capacity(count - 1);
    for j in 0..count - 1 {
        let mut rings = Vec::with_capacity(shape.len());
        for (r, &n) in shape.iter().enumerate() {
            let mut ring = Vec::with_capacity(n);
            for k in 0..n {
                ring.push(edge(
                    &mut builder,
                    (vertices[j][r][k], vertices[j + 1][r][k]),
                    (points[j][r][k], points[j + 1][r][k]),
                )?);
            }
            rings.push(ring);
        }
        rails.push(rings);
    }
    let mut faces = Vec::with_capacity(2 + walls);
    // The caps, in their sections' planes.
    for (j, orientation) in [
        (0, Orientation::Reversed),
        (count - 1, Orientation::Forward),
    ] {
        let (frame, rings) = &sections[j];
        let surface = builder.add_surface(Surface::Plane(Plane { frame: *frame }));
        let mut bounds = Vec::with_capacity(shape.len());
        for (r, ring) in rings.iter().enumerate() {
            let n = ring.len();
            let mut uses = Vec::with_capacity(n);
            for k in 0..n {
                let (a, b) = (ring[k], ring[(k + 1) % n]);
                let pcurve = builder.add_curve2(Curve2::Line(Line2 {
                    origin: Vec2::new(a.x, a.y),
                    direction: b - a,
                }));
                uses.push(EdgeUse {
                    edge: rims[j][r][k],
                    orientation: Orientation::Forward,
                    pcurve: Some(pcurve),
                });
            }
            let loop_id = add_loop(&mut builder, uses, vec![Interval::UNIT; n]);
            bounds.push(FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: r == 0,
            });
        }
        faces.push(builder.topology_mut().add_face(Face {
            surface: Some(surface),
            bounds,
            orientation,
        }));
    }
    // The walls: one plane per ring edge and span.
    let linear = tolerance.linear();
    for j in 0..count - 1 {
        for (r, &n) in shape.iter().enumerate() {
            for k in 0..n {
                let m = (k + 1) % n;
                let (a, b) = (points[j][r][k], points[j][r][m]);
                let (c, d) = (points[j + 1][r][k], points[j + 1][r][m]);
                // Newell's normal of the loop a -> b -> d -> c.
                let corners = [a, b, d, c];
                let mut normal = Vec3::ZERO;
                for i in 0..4 {
                    let (p, q) = (corners[i], corners[(i + 1) % 4]);
                    normal += p.cross(q);
                }
                let length = normal.length();
                let x = (b - a).normalize_or_zero();
                if !(length.is_finite() && length > 0.0) || x == Vec3::ZERO {
                    return Err(GeomError::Degenerate(
                        "an exact loft wall has no plane".into(),
                    ));
                }
                let z = normal / length;
                let centre = 0.25 * (a + b + c + d);
                if corners.iter().any(|p| (*p - centre).dot(z).abs() > linear) {
                    return Err(unsupported(WARPED_WALL));
                }
                // The plane through the corners' centroid, with `x` along
                // the section edge made square to the normal.
                let x = (x - z * x.dot(z)).normalize();
                let y = z.cross(x);
                let frame = Frame3 {
                    origin: centre,
                    x,
                    y,
                    z,
                };
                let surface = builder.add_surface(Surface::Plane(Plane { frame }));
                let uv = |p: Point3| Vec2::new((p - centre).dot(x), (p - centre).dot(y));
                let line = |from: Point3, to: Point3| {
                    Curve2::Line(Line2 {
                        origin: uv(from),
                        direction: uv(to) - uv(from),
                    })
                };
                let sides = [
                    (rims[j][r][k], Orientation::Forward, line(a, b)),
                    (rails[j][r][m], Orientation::Forward, line(b, d)),
                    (rims[j + 1][r][k], Orientation::Reversed, line(c, d)),
                    (rails[j][r][k], Orientation::Reversed, line(a, c)),
                ];
                let mut uses = Vec::with_capacity(4);
                let mut intervals = Vec::with_capacity(4);
                for (edge, orientation, pcurve) in sides {
                    let pcurve = builder.add_curve2(pcurve);
                    uses.push(EdgeUse {
                        edge,
                        orientation,
                        pcurve: Some(pcurve),
                    });
                    intervals.push(match orientation {
                        Orientation::Forward => Interval::UNIT,
                        Orientation::Reversed => Interval::new(1.0, 0.0),
                    });
                }
                let loop_id = add_loop(&mut builder, uses, intervals);
                faces.push(builder.topology_mut().add_face(Face {
                    surface: Some(surface),
                    bounds: vec![FaceBound {
                        loop_id,
                        orientation: Orientation::Forward,
                        outer: true,
                    }],
                    orientation: Orientation::Forward,
                }));
            }
        }
    }
    finish_closed(builder, faces)
}
