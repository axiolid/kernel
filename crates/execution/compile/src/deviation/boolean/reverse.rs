//! The mesh's distance to the exact result (#329): the other direction.
//!
//! [`super::measure_faces`] certifies the distance from every point of the
//! exact result to the mesh, so it cannot see material the mesh has and
//! the exact result lacks: a clip whose mesh kept a strip the exact result
//! cut away reports the strip's neighbours' distance, not the strip's
//! height. The reverse quantity is the supremum, over the mesh's points, of
//! the distance to the exact result's surface.
//!
//! The exact surface has no distance query, so the mesh is certified
//! against a *witness*: a triangle mesh every point of which lies within
//! `error` of the exact surface. A point within `d` of a witness point is
//! within `d + error` of the exact surface. The witness does not have to
//! cover the surface, only to lie on it, so a sparse witness is sound and
//! only loosens the bound.
//!
//! # Scope: all-planar results
//!
//! Only an exact result whose faces are all planes is bounded
//! ([`super::of_boolean`]): a plane's witness is exact up to rounding and
//! the affine search runs to the measured distance. A curved result keeps
//! the one-sided bound (#348).
//!
//! # The witness
//!
//! Each face's domain polygon (the rings of [`super::Domain`], as its
//! segments run them) is triangulated by the certified ear clipper
//! ([`triangulate_with`], [`PinchPolicy::Accept`]), the triangles are
//! mapped through the plane's affine jet and turned outward by the face's
//! and its shell's orientation, and the vertices are welded across faces.
//! Welding and orientation are what let the search settle: the faces a
//! union splits out of one plane join into one flat region only across
//! shared, opposed edges (`certify::flat`), and without them the search
//! hits its cap.
//!
//! # Why `error` holds
//!
//! 1. A point `x` of the polygon `P` of a face lies within `2 * slack` of
//!    the face's domain `D`. If `d(x, bd P) > slack`, the parity argument in
//!    the "Its domain is the region its loops' pcurves enclose" bullet of
//!    the module's "Measured, not derived" puts `x` in `D`; otherwise
//!    some `b` in `bd P` lies within `slack` of `x`, and `b` lies within
//!    `slack` of `bd D`.
//! 2. The ear clipper's triangles lie in `P` (ADR 0083).
//! 3. A plane's parameterisation is affine, so
//!    `|T(x) - T(y)| <= (|du| + |dv|) * |x - y|`: the first term of the
//!    face's error is `2 * slack * (|du| + |dv|)`.
//! 4. Coordinates are relative to `origin`, the middle of the box the
//!    faces' jet centres span, so the search's rounding margins follow the
//!    result's extent, not its place in the world (a wall four million
//!    metres from the origin otherwise reads its planes' coplanarity
//!    against a margin of millimetres). A face's `q - origin` is rounded by
//!    at most half an ulp, and each mapped vertex by at most `9 * EPSILON`
//!    of `|q - origin| + reach`, `reach` the jet's span over the domain's
//!    box.
//! 5. Welding moves every point by at most `moved`; a triangle the weld
//!    drops only loses coverage.
//! 6. A face whose domain cannot be triangulated (no outer ring of three
//!    points, or an ear-clipper refusal) is skipped, not refused: the
//!    witness need not cover `S`, so a skipped face only loosens the bound.
//!    The mesh is subtracted from `origin` the same way, each vertex
//!    moving by less than `shift = EPSILON * |m - origin|`; the bound is
//!    widened by it, and by `error`.
//! 7. So `d(m, S) <= d(m, W) + error + shift` for a mesh point `m`, the
//!    surface `S` and the witness `W`. The search's sampled lower bound is only a
//!    stopping heuristic; every returned value is an upper bound, also when
//!    the search stops at its cell cap.

use axiolid_construct::profile::{triangulate_with, PinchPolicy, Rings};
use axiolid_core::{Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_mesh::TriMesh;

use super::{clip_half, widened, ExactFace};
use crate::certify::{certify_affine_capped, Cell, Coverage, Patch};

/// Cell evaluations the reverse search may spend: a quarter of a forward
/// search's.
const REVERSE_MAX_CELLS: usize = 1 << 17;

/// A mesh every point of which lies within `error` of the exact faces,
/// in coordinates relative to `origin`.
pub(super) struct Witness {
    pub(super) mesh: TriMesh,
    pub(super) error: Scalar,
    /// The point the mesh's coordinates are relative to.
    pub(super) origin: Point3,
}

/// The signed area of `ring`, positive when it runs counter-clockwise.
fn signed_area(ring: &[Point2]) -> Scalar {
    let n = ring.len();
    (0..n)
        .map(|k| ring[k].perp_dot(ring[(k + 1) % n]))
        .sum::<Scalar>()
        * 0.5
}

/// The triangulation of a face's domain polygon: its outer ring (the
/// bound the face marks outer) counter-clockwise, the others as holes. `None`
/// when the face marks no outer ring with at least three points or the
/// ear clipper refuses: the face then adds no witness triangle, which
/// only loses coverage (see "Why `error` holds", step 6).
fn triangulate_domain(
    rings: &[Vec<Point2>],
    outer: Option<usize>,
) -> Option<(Vec<Point2>, Vec<[u32; 3]>)> {
    let outer = outer.filter(|&k| rings.get(k).is_some_and(|r| r.len() >= 3))?;
    let oriented = |ring: &Vec<Point2>, ccw: bool| {
        let mut ring = ring.clone();
        if (signed_area(&ring) > 0.0) != ccw {
            ring.reverse();
        }
        ring
    };
    let polygon = Rings {
        outer: oriented(&rings[outer], true),
        holes: rings
            .iter()
            .enumerate()
            .filter(|&(k, r)| k != outer && r.len() >= 3)
            .map(|(_, r)| oriented(r, false))
            .collect(),
    };
    triangulate_with(&polygon, PinchPolicy::Accept).ok()
}

/// The witness of `faces`: their domain polygons triangulated, mapped onto
/// their planes relative to a common origin, turned outward and welded.
pub(super) fn witness(faces: &[ExactFace<'_>]) -> Result<Witness, &'static str> {
    // The jets are taken at the domain boxes' centres; the origin is the
    // middle of the box they span, so every coordinate below is relative to
    // the result's extent, not to its place in the world.
    let centres: Vec<Option<(Point2, Point3, Vec3, Vec3)>> = faces
        .iter()
        .map(|face| {
            let d = &face.domain;
            let c = Point2::new(0.5 * (d.lo.x + d.hi.x), 0.5 * (d.lo.y + d.hi.y));
            face.jet(c.x, c.y).map(|(q, du, dv)| (c, q, du, dv))
        })
        .collect();
    let (mut lo, mut hi) = (
        Vec3::splat(Scalar::INFINITY),
        Vec3::splat(Scalar::NEG_INFINITY),
    );
    for (_, q, _, _) in centres.iter().flatten() {
        lo = lo.min(*q);
        hi = hi.max(*q);
    }
    let origin = (lo + hi) * 0.5;
    if !origin.is_finite() {
        return Err("boolean: mesh to exact: an exact face cannot be evaluated");
    }
    let mut positions: Vec<Point3> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    let mut error: Scalar = 0.0;
    // The world scale of the coordinates the faces' planes were rounded at.
    let mut scale: Scalar = 0.0;
    for (face, jet) in faces.iter().zip(&centres) {
        let domain = &face.domain;
        let Some((c, q, du, dv)) = *jet else {
            return Err("boolean: mesh to exact: an exact face cannot be evaluated");
        };
        let Some((points, triangles)) = triangulate_domain(&domain.rings, domain.outer) else {
            continue;
        };
        let rel = q - origin;
        let (sx, sy) = (
            (c.x - domain.lo.x).max(domain.hi.x - c.x),
            (c.y - domain.lo.y).max(domain.hi.y - c.y),
        );
        let reach = du.length() * sx + dv.length() * sy;
        // The jet's origin moved by `origin` and each vertex rounded.
        let rounding = 9.0 * Scalar::EPSILON * (rel.abs().max_element() + reach);
        scale = scale.max(q.abs().max_element() + reach);
        let off = 2.0 * domain.slack * (du.length() + dv.length());
        if !(rounding.is_finite() && off.is_finite()) {
            return Err("boolean: mesh to exact: an exact face cannot be evaluated");
        }
        error = error.max(widened(off, rounding));
        let base = positions.len() as u32;
        positions.extend(
            points
                .iter()
                .map(|p| rel + du * (p.x - c.x) + dv * (p.y - c.y)),
        );
        for t in &triangles {
            let [a, b, c] = t.map(|i| base + i);
            indices.extend(if face.outward { [a, b, c] } else { [a, c, b] });
        }
    }
    if indices.is_empty() {
        return Err("boolean: mesh to exact: the exact result has no face to measure against");
    }
    let (mesh, moved) = weld(&positions, &indices, scale);
    if mesh.indices.is_empty() {
        return Err("boolean: mesh to exact: the exact result has no face to measure against");
    }
    Ok(Witness {
        mesh,
        error: widened(error, moved),
        origin,
    })
}

/// Ulps of the witness's scale within which two vertices are one: the
/// rounding two faces' independent evaluations of a shared corner differ by.
const WELD_ULPS: Scalar = 128.0;

/// `positions` and `indices` with vertices within [`WELD_ULPS`] of each
/// other merged into the first of them, and the farthest any vertex moved.
/// `scale` is the largest world coordinate the faces' planes were rounded
/// at, whatever origin the positions are relative to.
///
/// Two faces share a corner, but each reaches it from its own centre, so
/// the copies differ in the last places and the faces' triangles share no
/// edge. [`crate::certify`] covers a piece by a flat region only across
/// shared edges. A triangle's points move by at most the farthest vertex
/// moved, which the caller adds to the witness's error. A triangle that
/// loses a corner to a merge is dropped: its points lie on the edge it
/// keeps.
fn weld(positions: &[Point3], indices: &[u32], scale: Scalar) -> (TriMesh, Scalar) {
    use std::collections::HashMap;
    let pitch = WELD_ULPS * Scalar::EPSILON * scale;
    if !(pitch > 0.0 && pitch.is_finite()) {
        return (TriMesh::new(positions.to_vec(), indices.to_vec()), 0.0);
    }
    let key = |p: Point3| [p.x, p.y, p.z].map(|c| (c / pitch).floor() as i64);
    let mut kept: Vec<Point3> = Vec::new();
    let mut cells: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
    let mut moved: Scalar = 0.0;
    let mut map = Vec::with_capacity(positions.len());
    for &p in positions {
        let [kx, ky, kz] = key(p);
        let mut found = None;
        'search: for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(ids) = cells.get(&[kx + dx, ky + dy, kz + dz]) else {
                        continue;
                    };
                    for &id in ids {
                        if (kept[id as usize] - p).length() <= pitch {
                            found = Some(id);
                            break 'search;
                        }
                    }
                }
            }
        }
        let id = match found {
            Some(id) => {
                moved = moved.max((kept[id as usize] - p).length());
                id
            }
            None => {
                let id = kept.len() as u32;
                kept.push(p);
                cells.entry([kx, ky, kz]).or_default().push(id);
                id
            }
        };
        map.push(id);
    }
    let welded = indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| [map[t[0] as usize], map[t[1] as usize], map[t[2] as usize]])
        .filter(|t| t[0] != t[1] && t[1] != t[2] && t[2] != t[0])
        .flatten()
        .collect();
    (TriMesh::new(kept, welded), moved)
}

/// One triangle of the mesh as an affine patch over the unit parameter
/// triangle `s, t >= 0, s + t <= 1`.
struct MeshTriangle {
    a: Point3,
    ab: Vec3,
    ac: Vec3,
}

impl Patch for MeshTriangle {
    fn jet(&self, s: Scalar, t: Scalar) -> Option<(Point3, Vec3, Vec3)> {
        Some((self.a + self.ab * s + self.ac * t, self.ab, self.ac))
    }

    fn second(&self, _x: (Scalar, Scalar), _y: (Scalar, Scalar)) -> Option<[Scalar; 3]> {
        Some([0.0; 3])
    }

    fn coverage(&self, x: (Scalar, Scalar), y: (Scalar, Scalar)) -> Coverage {
        if x.1 + y.1 <= 1.0 {
            Coverage::Inside
        } else if x.0 + y.0 > 1.0 {
            Coverage::Outside
        } else {
            Coverage::Boundary
        }
    }

    fn clip(&self, x: (Scalar, Scalar), y: (Scalar, Scalar)) -> Option<Vec<Vec<Point2>>> {
        let (lo, hi) = (Point2::new(x.0, y.0), Point2::new(x.1, y.1));
        let cell = [lo, Point2::new(hi.x, lo.y), hi, Point2::new(lo.x, hi.y)];
        let piece = clip_half(
            &cell,
            Point2::new(1.0, 0.0),
            Vec2::new(-1.0, -1.0) / core::f64::consts::SQRT_2,
            4.0 * Scalar::EPSILON,
        );
        Some(vec![piece].into_iter().filter(|p| !p.is_empty()).collect())
    }
}

/// A certified bound on the distance from every point of `mesh` to the
/// surface of `faces` (all planar), or why there is none. `flattened` is
/// the mesh-side deviation the report already carries (it lets the budget
/// settle the search, as for the forward one) and `added` is what the
/// caller widens the result by whatever it finds, the flattening, the
/// floor and the readings; the search settles once the bound is within
/// what `budget` leaves after `added` and the errors below.
pub(super) fn mesh_to_exact(
    faces: &[ExactFace<'_>],
    mesh: &TriMesh,
    budget: Scalar,
    flattened: Scalar,
    added: Scalar,
) -> Result<Scalar, &'static str> {
    mesh_to_exact_capped(faces, mesh, budget, flattened, added, REVERSE_MAX_CELLS)
}

/// [`mesh_to_exact`] with a search of at most `max_cells` cell evaluations.
fn mesh_to_exact_capped(
    faces: &[ExactFace<'_>],
    mesh: &TriMesh,
    budget: Scalar,
    flattened: Scalar,
    added: Scalar,
    max_cells: usize,
) -> Result<Scalar, &'static str> {
    const NOT_FINITE: &str = "boolean: mesh to exact: the mesh has a triangle that is not finite";
    if mesh.indices.is_empty() {
        return Ok(0.0);
    }
    let wit = witness(faces)?;
    // The mesh relative to the witness's origin; each coordinate difference
    // is rounded by at most half an ulp of its result, so a vertex moves by
    // less than `EPSILON` of its largest relative coordinate.
    let mut shift: Scalar = 0.0;
    let mut relative = Vec::with_capacity(mesh.positions.len());
    for p in &mesh.positions {
        let r = *p - wit.origin;
        if !(p.is_finite() && r.is_finite()) {
            return Err(NOT_FINITE);
        }
        shift = shift.max(Scalar::EPSILON * r.abs().max_element());
        relative.push(r);
    }
    let mut triangles = Vec::with_capacity(mesh.indices.len() / 3);
    for corners in mesh.indices.as_chunks::<3>().0 {
        let at = |k: usize| relative.get(corners[k] as usize).copied();
        let (Some(a), Some(b), Some(c)) = (at(0), at(1), at(2)) else {
            return Err(NOT_FINITE);
        };
        let (ab, ac) = (b - a, c - a);
        if !(ab.is_finite() && ac.is_finite()) {
            return Err(NOT_FINITE);
        }
        triangles.push(MeshTriangle { a, ab, ac });
    }
    if triangles.is_empty() {
        return Ok(0.0);
    }
    // What the bound is widened by: the witness's error and the mesh's
    // recentring.
    let error = widened(wit.error, shift);
    let patches: Vec<&dyn Patch> = triangles.iter().map(|t| t as &dyn Patch).collect();
    let cells = (0..triangles.len())
        .map(|patch| Cell {
            patch,
            x: (0.0, 1.0),
            y: (0.0, 1.0),
        })
        .collect();
    let widening = widened(added, error);
    let settle = if widening < budget {
        budget - widening
    } else {
        budget
    };
    let distance = certify_affine_capped(&patches, cells, &wit.mesh, settle, flattened, max_cells)
        .ok_or("boolean: mesh to exact: the mesh could not be bounded against the exact result")?;
    Ok(widened(distance, error))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::deviation::boolean::{exact_faces, of_boolean};
    use crate::deviation::{Deviation, DeviationBound, DeviationPath, MESH_TO_EXACT};
    use crate::{ReferenceExactCompiler, ReferenceMeshCompiler};
    use axiolid_contracts::ExecutionOptions;
    use axiolid_core::{BooleanOperator, Plane3, Tolerance, Transform3};
    use axiolid_curve::{Curve2, Polyline2};
    use axiolid_exact_compile_contract::ExactCompiler;
    use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
    use axiolid_mesh_compile_contract::MeshCompiler;
    use axiolid_model::{
        GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
    };
    use axiolid_primitive::HalfSpace;
    use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

    /// A graph, the boolean's node, the uncut operand's node and the options.
    type Fixture = (GeometryGraph, NodeId, NodeId, ExecutionOptions);

    fn options() -> ExecutionOptions {
        ExecutionOptions::new(Tolerance::new(1e-3, 1e-9).unwrap())
            .with_chord_error(1e-3)
            .unwrap()
    }

    fn rectangle(g: &mut GeometryGraphBuilder, x: Scalar, y: Scalar) -> NodeId {
        g.push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x,
                y,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )))
        .unwrap()
    }

    fn extrusion(g: &mut GeometryGraphBuilder, profile: NodeId, depth: Scalar) -> NodeId {
        g.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
        .unwrap()
    }

    fn boolean(
        g: &mut GeometryGraphBuilder,
        left: NodeId,
        right: NodeId,
        operator: BooleanOperator,
    ) -> NodeId {
        g.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator,
        }))
        .unwrap()
    }

    /// A wall `6 x 0.35 x 7.76` clipped by the half-space above `z = 3`.
    fn clipped_wall() -> Fixture {
        let mut g = GeometryGraphBuilder::new();
        let profile = rectangle(&mut g, 6.0, 0.35);
        let wall = extrusion(&mut g, profile, 7.76);
        let half = g
            .push(GeometryNode::HalfSpace(HalfSpace {
                boundary: Plane3 {
                    origin: Point3::new(0.0, 0.0, 3.0),
                    normal: Vec3::Z,
                },
                agreement: true,
            }))
            .unwrap();
        let clip = boolean(&mut g, wall, half, BooleanOperator::Difference);
        (g.finish(vec![clip, wall]).unwrap(), clip, wall, options())
    }

    /// The wall of #325, in millimetres, less its bounded half-space strip.
    fn wall_325() -> Fixture {
        const WIDTH: Scalar = 238.549_352_729_28;
        let mut g = GeometryGraphBuilder::new();
        let profile = rectangle(&mut g, 10125.0, 350.0);
        let wall = extrusion(&mut g, profile, 16165.033064627303);
        let wall = g
            .push(GeometryNode::Instance(Instance {
                source: wall,
                transform: Transform3::from_cols(
                    -Vec3::X,
                    -Vec3::Y,
                    Vec3::Z,
                    Vec3::new(5062.5, 0.0, 0.0),
                ),
            }))
            .unwrap();
        let origin = Point3::new(9886.45064727072, -175.0, 11404.747431413967);
        let half_space = g
            .push(GeometryNode::HalfSpace(HalfSpace {
                boundary: Plane3 {
                    origin,
                    normal: Vec3::Z,
                },
                agreement: true,
            }))
            .unwrap();
        let boundary = g
            .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
                points: vec![
                    Point2::new(0.0, 0.0),
                    Point2::new(WIDTH, 0.0),
                    Point2::new(WIDTH, 350.0),
                    Point2::new(0.0, 350.0),
                    Point2::new(0.0, 0.0),
                ],
                closed: true,
            })))
            .unwrap();
        let clip = g
            .push(GeometryNode::SolidOperation(
                SolidOperation::BoundedHalfSpace {
                    half_space,
                    boundary,
                    placement: Transform3::from_translation(origin),
                },
            ))
            .unwrap();
        let cut = boolean(&mut g, wall, clip, BooleanOperator::Difference);
        let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
        (g.finish(vec![cut, wall]).unwrap(), cut, wall, options)
    }

    /// A wall `6 x 0.3 x 3` less a `1 x 1` window through it (turned across
    /// the wall) or, with `round`, less a round window of radius `0.4`.
    fn windowed_wall(round: bool) -> Fixture {
        let mut g = GeometryGraphBuilder::new();
        let profile = rectangle(&mut g, 6.0, 0.3);
        let wall = extrusion(&mut g, profile, 3.0);
        let profile = if round {
            g.push(GeometryNode::Profile(Profile::Circle(CircleProfile {
                radius: 0.4,
                thickness: None,
            })))
            .unwrap()
        } else {
            rectangle(&mut g, 1.0, 1.0)
        };
        let window = extrusion(&mut g, profile, 1.0);
        let window = g
            .push(GeometryNode::Instance(Instance {
                source: window,
                transform: Transform3::from_translation(Vec3::new(0.8, 0.5, 1.5))
                    * Transform3::from_rotation_x(core::f64::consts::FRAC_PI_2),
            }))
            .unwrap();
        let cut = boolean(&mut g, wall, window, BooleanOperator::Difference);
        (g.finish(vec![cut, wall]).unwrap(), cut, wall, options())
    }

    /// Two `4 x 0.3 x 3` walls under a building placement, the second turned
    /// a quarter and moved by `offset`.
    fn two_walls(offset: Vec3, operator: BooleanOperator) -> (GeometryGraph, NodeId) {
        let building = Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2))
            * Transform3::from_rotation_z(0.6);
        let mut g = GeometryGraphBuilder::new();
        let wall = |g: &mut GeometryGraphBuilder, transform: Transform3| {
            let profile = rectangle(g, 4.0, 0.3);
            let source = extrusion(g, profile, 3.0);
            g.push(GeometryNode::Instance(Instance { source, transform }))
                .unwrap()
        };
        let a = wall(&mut g, building);
        let b = wall(
            &mut g,
            building
                * Transform3::from_translation(offset)
                * Transform3::from_rotation_z(core::f64::consts::FRAC_PI_2),
        );
        let root = boolean(&mut g, a, b, operator);
        (g.finish(vec![root]).unwrap(), root)
    }

    fn mesh_of(graph: &GeometryGraph, node: NodeId, options: &ExecutionOptions) -> TriMesh {
        ReferenceMeshCompiler::new(BoolmeshBoolean::new())
            .compile_mesh_reported(graph, node, options)
            .unwrap()
            .mesh
    }

    /// The certified bound of the `MESH_TO_EXACT` contribution.
    fn reverse(deviation: &Deviation) -> Option<Scalar> {
        let found = deviation
            .0
            .iter()
            .find(|c| c.path == DeviationPath::Boolean && c.detail == MESH_TO_EXACT)?;
        match found.bound {
            DeviationBound::Certified(bound) => Some(bound),
            other => panic!("{other:?}"),
        }
    }

    /// The boolean's own contribution: the exact result's distance to the mesh.
    fn forward(deviation: &Deviation) -> Scalar {
        deviation.0[0].bound.value().expect("bounded")
    }

    #[test]
    fn extra_material_in_the_mesh_is_bounded() {
        let (graph, clip, wall, options) = clipped_wall();
        let deviation = of_boolean(&graph, clip, &options, &mesh_of(&graph, wall, &options));
        let height = reverse(&deviation).expect("the other direction");
        assert!((4.76..=1.2 * 4.76).contains(&height), "{height}");
        assert!(forward(&deviation) < 0.18, "{deviation:?}");
        let deviation = of_boolean(&graph, clip, &options, &mesh_of(&graph, clip, &options));
        let cut = reverse(&deviation).expect("the other direction");
        assert!(cut <= 1e-10, "{cut:e}");
    }

    #[test]
    fn the_325_wall_reports_its_uncut_strip() {
        const WIDTH: Scalar = 238.549_352_729_28;
        let (graph, cut, wall, options) = wall_325();
        let deviation = of_boolean(&graph, cut, &options, &mesh_of(&graph, wall, &options));
        let strip = reverse(&deviation).expect("the other direction");
        assert!((WIDTH..=1.2 * WIDTH).contains(&strip), "{strip}");
        assert!(forward(&deviation) < 175.001, "{deviation:?}");
    }

    #[test]
    fn an_uncut_window_shows_and_a_cut_one_is_at_rounding() {
        let (graph, cut, wall, options) = windowed_wall(false);
        let deviation = of_boolean(&graph, cut, &options, &mesh_of(&graph, wall, &options));
        let uncut = reverse(&deviation).expect("the other direction");
        assert!((0.5..=0.6).contains(&uncut), "{deviation:?}");
        let deviation = of_boolean(&graph, cut, &options, &mesh_of(&graph, cut, &options));
        let matching = reverse(&deviation).expect("the other direction");
        assert!(matching <= 1e-10, "{deviation:?}");
    }

    #[test]
    fn the_witness_lies_on_the_exact_faces_and_covers_them() {
        let (graph, clip, _, options) = clipped_wall();
        let exact = ReferenceExactCompiler::new()
            .compile_exact(&graph, clip, &options)
            .unwrap();
        let faces = exact_faces(&exact, 1e-3, 0.0).unwrap();
        let witness = witness(&faces).unwrap();
        let (lo, hi) = (Point3::new(-3.0, -0.175, 0.0), Point3::new(3.0, 0.175, 3.0));
        let centre = (lo + hi) * 0.5;
        let to_boundary = |p: Point3| {
            let outside = (lo - p).max(p - hi).max(Point3::ZERO).length();
            let inside = (p - lo).min(hi - p).min_element();
            if outside > 0.0 {
                outside
            } else {
                inside.max(0.0)
            }
        };
        let mesh = &witness.mesh;
        let mut area = 0.0;
        for t in mesh.indices.as_chunks::<3>().0 {
            let [a, b, c] = [0, 1, 2].map(|k| mesh.positions[t[k] as usize] + witness.origin);
            let normal = (b - a).cross(c - a);
            area += 0.5 * normal.length();
            let middle = (a + b + c) / 3.0;
            assert!(normal.dot(middle - centre) > 0.0, "inward at {middle:?}");
            for p in [a, b, c, middle] {
                assert!(
                    to_boundary(p) <= witness.error + 1e-15,
                    "{p:?}: {:e} over {:e}",
                    to_boundary(p),
                    witness.error
                );
            }
        }
        let box_area = 2.0 * (6.0 * 0.35 + 6.0 * 3.0 + 0.35 * 3.0);
        assert!((area - box_area).abs() <= 1e-9, "{area}");
    }

    #[test]
    fn the_witness_of_coplanar_split_faces_is_closed() {
        let options = options();
        let (graph, root) = two_walls(Vec3::new(1.5, 1.8, 0.0), BooleanOperator::Union);
        let exact = ReferenceExactCompiler::new()
            .compile_exact(&graph, root, &options)
            .unwrap();
        let faces = exact_faces(&exact, 1e-3, 0.0).unwrap();
        let witness = witness(&faces).unwrap();
        let mut edges = std::collections::HashMap::new();
        for t in witness.mesh.indices.as_chunks::<3>().0 {
            for k in 0..3 {
                *edges.entry((t[k], t[(k + 1) % 3])).or_insert(0) += 1;
            }
        }
        for (&(a, b), &n) in &edges {
            assert_eq!(n, 1, "edge {a}-{b} used {n} times");
            assert_eq!(edges.get(&(b, a)), Some(&1), "edge {a}-{b} has no opposite");
        }
        let deviation = of_boolean(&graph, root, &options, &mesh_of(&graph, root, &options));
        let union = reverse(&deviation).expect("the other direction");
        assert!(union <= 1e-10, "{deviation:?}");
    }

    #[test]
    fn the_reverse_settles_within_what_its_readings_leave_of_the_budget() {
        let options = options();
        let (graph, root) = two_walls(Vec3::new(1.5, 1.8, 3e-8), BooleanOperator::Union);
        let deviation = of_boolean(&graph, root, &options, &mesh_of(&graph, root, &options));
        let apart = reverse(&deviation).expect("the other direction");
        assert!((3e-8..=1e-3).contains(&apart), "{deviation:?}");
    }

    #[test]
    fn a_curved_result_adds_no_mesh_to_exact_contribution() {
        let (graph, cut, wall, options) = windowed_wall(true);
        let deviation = of_boolean(&graph, cut, &options, &mesh_of(&graph, wall, &options));
        assert_eq!(reverse(&deviation), None, "{deviation:?}");
        assert!(deviation.bound().is_some(), "{deviation:?}");
    }

    /// The wall of the far-placement regression: `5.524 x 0.2037 x 3.055`
    /// placed far from the origin and turned, less three openings each
    /// turned about its axis. With `clip`, instead one half-space cut above
    /// the wall's middle.
    fn far_wall(clip: bool) -> Fixture {
        let placement = Transform3::from_translation(Vec3::new(3910692.85, -2737484.99, 5.50))
            * Transform3::from_rotation_z(-1.0);
        let mut g = GeometryGraphBuilder::new();
        let profile = rectangle(&mut g, 5.524, 0.2037);
        let wall = extrusion(&mut g, profile, 3.055);
        let wall = g
            .push(GeometryNode::Instance(Instance {
                source: wall,
                transform: placement,
            }))
            .unwrap();
        let mut cut = wall;
        if clip {
            let half = g
                .push(GeometryNode::HalfSpace(HalfSpace {
                    boundary: Plane3 {
                        origin: Point3::new(3910692.85, -2737484.99, 5.50 + 1.5),
                        normal: Vec3::Z,
                    },
                    agreement: true,
                }))
                .unwrap();
            cut = boolean(&mut g, wall, half, BooleanOperator::Difference);
        } else {
            for (x, turn) in [(-1.8, 0.025), (0.0, -0.22), (1.8, 0.29)] {
                let profile = rectangle(&mut g, 0.9, 1.2);
                let window = extrusion(&mut g, profile, 1.0);
                let window = g
                    .push(GeometryNode::Instance(Instance {
                        source: window,
                        transform: placement
                            * Transform3::from_translation(Vec3::new(x, 0.5, 1.5))
                            * Transform3::from_rotation_x(core::f64::consts::FRAC_PI_2)
                            * Transform3::from_rotation_z(turn),
                    }))
                    .unwrap();
                cut = boolean(&mut g, cut, window, BooleanOperator::Difference);
            }
        }
        (g.finish(vec![cut, wall]).unwrap(), cut, wall, options())
    }

    #[test]
    fn a_far_placed_mesh_is_bounded_near_rounding() {
        let (graph, cut, _, options) = far_wall(false);
        let mesh = mesh_of(&graph, cut, &options);
        let (_, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
            .compile_mesh_with_deviation(&graph, cut, &options)
            .unwrap();
        let deviation = of_boolean(&graph, cut, &options, &mesh);
        let far = reverse(&deviation).expect("the other direction");
        // 1.6e-5 recentred; 3.2e-4 read at the world coordinates.
        assert!(far <= 1e-4, "{deviation:?}");
        assert!(report.meets_requested(), "{report:?}");
    }

    #[test]
    fn a_far_placed_uncut_strip_still_shows() {
        let (graph, cut, wall, options) = far_wall(true);
        let deviation = of_boolean(&graph, cut, &options, &mesh_of(&graph, wall, &options));
        let strip = reverse(&deviation).expect("the other direction");
        assert!((1.5..=1.2 * 1.55).contains(&strip), "{strip}");
    }

    /// The wall `6 x 0.3 x 3` united with a box at `(2.5, y, 0)` turned by
    /// `turn` about its axis.
    fn wall_and_box(
        size: (Scalar, Scalar, Scalar),
        y: Scalar,
        turn: Scalar,
    ) -> (GeometryGraph, NodeId) {
        let mut g = GeometryGraphBuilder::new();
        let profile = rectangle(&mut g, 6.0, 0.3);
        let wall = extrusion(&mut g, profile, 3.0);
        let profile = rectangle(&mut g, size.0, size.1);
        let block = extrusion(&mut g, profile, size.2);
        let block = g
            .push(GeometryNode::Instance(Instance {
                source: block,
                transform: Transform3::from_translation(Vec3::new(2.5, y, 0.0))
                    * Transform3::from_rotation_z(turn),
            }))
            .unwrap();
        let root = boolean(&mut g, wall, block, BooleanOperator::Union);
        (g.finish(vec![root]).unwrap(), root)
    }

    #[test]
    fn readings_do_not_inflate_the_reverse_bound() {
        let options = options();
        for (name, size, y, turn) in [
            ("turned 1e-5", (2.0, 0.4, 3.0), 0.2, 1e-5),
            ("turned 3e-4", (2.0, 0.4, 3.0), 0.2, 3e-4),
            ("1e-9 off the wall", (2.0, 1.0, 3.5), 0.65 - 1e-9, 0.0),
        ] {
            let (graph, root) = wall_and_box(size, y, turn);
            let deviation = of_boolean(&graph, root, &options, &mesh_of(&graph, root, &options));
            let back = reverse(&deviation).expect("the other direction");
            assert!(back < 1e-4, "{name}: {deviation:?}");
            assert!(
                back <= 2.0 * forward(&deviation) + 1e-9,
                "{name}: {deviation:?}"
            );
        }
    }

    #[test]
    fn a_face_that_cannot_be_triangulated_is_skipped_not_refused() {
        let (graph, clip, wall, options) = clipped_wall();
        let exact = ReferenceExactCompiler::new()
            .compile_exact(&graph, clip, &options)
            .unwrap();
        let mut faces = exact_faces(&exact, 1e-3, 0.0).unwrap();
        let whole = witness(&faces).unwrap().mesh.indices.len();
        // No outer ring, two points, three collinear points: all skipped.
        faces[0].domain.outer = None;
        faces[1].domain.rings = vec![vec![Point2::ZERO, Point2::X]];
        faces[2].domain.rings = vec![vec![Point2::ZERO, Point2::X, Point2::new(2.0, 0.0)]];
        faces[2].domain.outer = Some(0);
        let sparse = witness(&faces).unwrap().mesh.indices.len();
        assert!(sparse < whole, "{sparse} {whole}");
        // The sparse witness still bounds the uncut mesh from above.
        let mesh = mesh_of(&graph, wall, &options);
        let skipped = mesh_to_exact(&faces, &mesh, 1e-3, 0.0, 0.0).unwrap();
        assert!(skipped >= 4.76, "{skipped}");
        // Nothing left to triangulate is refused by name.
        for face in &mut faces {
            face.domain.rings.clear();
        }
        assert!(witness(&faces).is_err());
    }

    #[test]
    fn a_search_stopped_at_its_cap_is_still_an_upper_bound() {
        let (graph, clip, wall, options) = clipped_wall();
        let exact = ReferenceExactCompiler::new()
            .compile_exact(&graph, clip, &options)
            .unwrap();
        let faces = exact_faces(&exact, 1e-3, 0.0).unwrap();
        for mesh in [
            mesh_of(&graph, wall, &options),
            mesh_of(&graph, clip, &options),
        ] {
            let full = mesh_to_exact(&faces, &mesh, 1e-3, 0.0, 0.0).unwrap();
            for cap in [1, 8, 64] {
                let capped = mesh_to_exact_capped(&faces, &mesh, 1e-3, 0.0, 0.0, cap).unwrap();
                assert!(
                    capped >= full * (1.0 - 1e-9),
                    "cap {cap}: {capped} below {full}"
                );
            }
        }
    }

    #[test]
    fn the_witness_error_holds_the_slack_and_the_rounding_of_every_face() {
        let (graph, clip, _, options) = clipped_wall();
        let exact = ReferenceExactCompiler::new()
            .compile_exact(&graph, clip, &options)
            .unwrap();
        let faces = exact_faces(&exact, 1e-3, 0.0).unwrap();
        let error = witness(&faces).unwrap().error;
        for face in &faces {
            let d = &face.domain;
            let c = Point2::new(0.5 * (d.lo.x + d.hi.x), 0.5 * (d.lo.y + d.hi.y));
            let (_, du, dv) = face.jet(c.x, c.y).unwrap();
            let reach = du.length() * (d.hi.x - d.lo.x) + dv.length() * (d.hi.y - d.lo.y);
            let slack = 2.0 * d.slack * (du.length() + dv.length());
            assert!(
                error >= widened(slack, Scalar::EPSILON * reach),
                "{error:e} below slack {slack:e}, reach {reach}"
            );
        }
    }

    /// A mesh triangle as a patch: cells wholly inside, wholly outside and
    /// across its hypotenuse are told apart, and a boundary cell's pieces
    /// hold the cell's part in the triangle and reach no further than the
    /// clip's rounding.
    #[test]
    fn mesh_triangle_cells() {
        let triangle = MeshTriangle {
            a: Point3::ZERO,
            ab: Vec3::X,
            ac: Vec3::Y,
        };
        let (mut inside, mut outside, mut boundary) = (0, 0, 0);
        for n in [1_u32, 2, 4, 8] {
            let step = 1.0 / Scalar::from(n);
            for i in 0..n {
                for j in 0..n {
                    let x = (Scalar::from(i) * step, Scalar::from(i + 1) * step);
                    let y = (Scalar::from(j) * step, Scalar::from(j + 1) * step);
                    let on = |s: Scalar, t: Scalar| s + t <= 1.0;
                    match triangle.coverage(x, y) {
                        Coverage::Inside => {
                            inside += 1;
                            assert!(on(x.1, y.1), "{x:?} {y:?}");
                        }
                        Coverage::Outside => {
                            outside += 1;
                            assert!(!on(x.0, y.0), "{x:?} {y:?}");
                        }
                        Coverage::Boundary => {
                            boundary += 1;
                            let pieces = triangle.clip(x, y).unwrap();
                            for s in 0..=8 {
                                for t in 0..=8 {
                                    let p = Point2::new(
                                        x.0 + (x.1 - x.0) * Scalar::from(s) / 8.0,
                                        y.0 + (y.1 - y.0) * Scalar::from(t) / 8.0,
                                    );
                                    if on(p.x, p.y) {
                                        assert!(
                                            pieces.iter().any(|piece| inside_convex(p, piece)),
                                            "{p:?} is in no piece of {x:?} x {y:?}"
                                        );
                                    }
                                }
                            }
                            for v in pieces.iter().flatten() {
                                assert!(v.x + v.y <= 1.0 + 1e-14, "{v:?}");
                            }
                        }
                    }
                }
            }
        }
        assert!(inside > 0 && outside > 0 && boundary > 0);
    }

    /// Whether `p` is in the convex `piece`, up to rounding.
    fn inside_convex(p: Point2, piece: &[Point2]) -> bool {
        let n = piece.len();
        let sides: Vec<Scalar> = (0..n)
            .map(|k| (piece[(k + 1) % n] - piece[k]).perp_dot(p - piece[k]))
            .collect();
        sides.iter().all(|&s| s >= -1e-12) || sides.iter().all(|&s| s <= 1e-12)
    }
}
