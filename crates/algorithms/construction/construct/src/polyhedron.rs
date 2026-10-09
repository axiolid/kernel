//! Exact boolean over general planar-faced solids (#77).
//!
//! # Why this is not a BSP tree
//!
//! A BSP boolean constructs split points recursively, so each generation of
//! cuts is computed from coordinates that were themselves computed. Error
//! compounds with depth and the exactness claim decays silently.
//!
//! Here every fragment is carried as a polygon whose plane is one of the
//! ORIGINAL input planes, never a derived one. A face is split only against
//! input planes, so a vertex is at worst one intersection away from input
//! data. Classification then asks a certified predicate which side of the
//! other solid a fragment lies on.
//!
//! # Split points are correctly rounded
//!
//! Split points are stored in f64 (ADR 0045, `plane_crossing`), each
//! coordinate the double NEAREST to the exact crossing. The earlier f64
//! formula rounded twice, so after many operations a vertex two operands
//! should share landed a few ULPs apart, and a split through it emitted a
//! fragment enclosing no area (#199). Correct rounding is a function of the
//! exact point alone, so that drift cannot arise, and a crossing that is a
//! double -- every cut of an axis-aligned edge by an axis-aligned plane --
//! is returned exactly.
//!
//! # Thin fragments are classified exactly
//!
//! The input's own coordinates can put two planes one ULP apart (`2/3`
//! computed two ways), so a correct result may contain fragments one ULP
//! wide. No double lies strictly inside such a fragment, so it is
//! classified at an exact dyadic interior point (`Probe`); a fragment that
//! encloses no area has no interior and is refused by name.
//!
//! # Collapsed fragments are not dropped
//!
//! Deleting zero-area fragments was tried and reverted (`40b5069`): the
//! chain then completes, but the holes it leaves make the shell integrate to
//! a plausibly wrong volume. A refusal is actionable; a wrong volume is
//! silent. Do not reintroduce a drop-based fix (#199).

use crate::boolean_exact::unsupported;
use crate::polyhedron_exact::{self as exact, ExactPoint, RingPosition};
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Point3, Vec3};
use axiolid_guarantees::Sign;
use axiolid_mesh::TriMesh;
use axiolid_predicates::{orient2d, orient3d};
use std::collections::BTreeMap;

/// A closed solid bounded by planar polygonal faces.
///
/// Each face is a vertex ring wound counter-clockwise seen from outside, so
/// the outward normal follows the right-hand rule. That convention is what
/// makes containment decidable without a separate inside/outside oracle.
#[derive(Debug, Clone, PartialEq)]
pub struct Polyhedron {
    faces: Vec<Vec<Point3>>,
}

/// Which boolean to evaluate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOp {
    /// Everything in either solid.
    Union,
    /// Only what lies in both.
    Intersection,
    /// The subject with the tool removed.
    Difference,
}

impl Polyhedron {
    /// Build from outward-wound planar faces.
    ///
    /// Faces are validated as planar here rather than trusted, because every
    /// later decision assumes it. A non-planar ring has no single plane to
    /// classify against, so accepting one would make the exactness claim
    /// meaningless.
    pub fn new(faces: Vec<Vec<Point3>>) -> GeomResult<Self> {
        if faces.len() < 4 {
            return Err(GeomError::InvalidInput(
                "a closed solid needs at least 4 faces".to_owned(),
            ));
        }
        for face in &faces {
            if face.len() < 3 {
                return Err(GeomError::InvalidInput(
                    "a face needs at least 3 vertices".to_owned(),
                ));
            }
            if face.iter().any(|p| !p.is_finite()) {
                return Err(GeomError::InvalidInput(
                    "face vertices must be finite".to_owned(),
                ));
            }
            for &v in &face[3..] {
                if orient3d(face[0], face[1], face[2], v).sign() != Some(Sign::Zero) {
                    return Err(GeomError::InvalidInput(
                        "face is not planar; no single plane to classify against".to_owned(),
                    ));
                }
            }
        }
        Ok(Self { faces })
    }

    /// The bounding faces, each an outward-wound ring.
    #[must_use]
    pub fn faces(&self) -> &[Vec<Point3>] {
        &self.faces
    }
}

/// Which side of a face's plane a point lies on, decided exactly.
///
/// Returns `None` when the predicate cannot certify a sign, which is the
/// signal to refuse rather than guess.
fn side_of_face(face: &[Point3], point: Point3) -> Option<Sign> {
    orient3d(face[0], face[1], face[2], point).sign()
}

/// Where a point sits relative to a solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Containment {
    Inside,
    OnBoundary,
    Outside,
}

/// The point a fragment is classified at: one point of its relative
/// interior, since after splitting every such point has the same answer.
///
/// The f64 centroid serves whenever it is certified to lie in the
/// fragment's plane and strictly inside its ring. A fragment one ULP wide
/// has no double strictly inside it: its centroid rounds onto its own edge,
/// which is often the other solid's boundary as well, and the answer is
/// then about the wrong point (#199). Such a fragment is classified at an
/// exact dyadic interior point instead, with every predicate evaluated
/// exactly. Nothing is perturbed either way.
#[derive(Debug, Clone)]
enum Probe {
    Float(Point3),
    Exact(ExactPoint),
}

impl Probe {
    /// Choose the classification point of a fragment, or refuse by name when
    /// the fragment encloses no area and so has no interior to classify.
    fn of(fragment: &[Point3]) -> GeomResult<Self> {
        let centroid = centroid_of(fragment);
        if side_of_face(fragment, centroid) == Some(Sign::Zero)
            && ring_position(fragment, centroid) == Some(RingPosition::Inside)
        {
            return Ok(Self::Float(centroid));
        }
        // A collapsed ring corresponds to nothing in the modelled solid, so
        // neither keeping nor dropping it is known to be right: dropping was
        // tried and integrated to a wrong volume (see the module docs).
        exact::interior_point(fragment)
            .map(Self::Exact)
            .ok_or_else(|| {
                unsupported(
                    "a split fragment encloses no area (collapsed ring); cannot classify it",
                )
            })
    }

    /// A nearby double, used only to size the probe ray.
    fn approx(&self) -> Point3 {
        match self {
            Self::Float(p) => *p,
            Self::Exact(p) => exact::approx(p),
        }
    }

    /// Which side of `face`'s plane the probe point lies on.
    fn side(&self, face: &[Point3]) -> Option<Sign> {
        match self {
            Self::Float(p) => side_of_face(face, *p),
            Self::Exact(p) => exact::side(face, p),
        }
    }

    /// Which side of `face`'s plane `point` lies on, in the same arithmetic
    /// and sign convention as [`Probe::side`], so the two can be compared.
    fn side_of(&self, face: &[Point3], point: Point3) -> Option<Sign> {
        match self {
            Self::Float(_) => side_of_face(face, point),
            Self::Exact(_) => exact::side(face, &exact::exact(point)?),
        }
    }

    /// The turn of edge `a`-`b` seen along the ray to `far`.
    fn turn(&self, far: Point3, a: Point3, b: Point3) -> Option<Sign> {
        match self {
            Self::Float(p) => orient3d(*p, far, a, b).sign(),
            Self::Exact(p) => exact::turn(p, far, a, b),
        }
    }

    /// Where the probe point sits relative to a ring it is coplanar with.
    fn ring(&self, face: &[Point3]) -> Option<RingPosition> {
        match self {
            Self::Float(p) => ring_position(face, *p),
            Self::Exact(p) => exact::ring_position(face, p),
        }
    }
}

/// Whether `point` is inside `solid`, by exact ray crossing parity.
///
/// A convex all-faces test is wrong for non-convex solids: a point in the
/// notch of an L-shaped prism is on the inner side of every face plane and
/// would be called inside. Parity counting is correct for any closed
/// orientable solid, convex or not.
///
/// The ray direction is chosen so it misses every vertex and edge. Rather
/// than perturbing coordinates -- which would forfeit exactness -- a
/// degenerate hit makes the whole operation refuse.
fn contains(solid: &Polyhedron, probe: &Probe, direction: Vec3) -> Option<Containment> {
    // The ray is represented by a segment, so it must be long enough to
    // leave the solid: a unit-length direction would miss every crossing
    // beyond it and invert the parity. Scaling by the solid's own extent
    // keeps the far endpoint outside for any input size. The far endpoint
    // is any point beyond the solid on a ray from the probe, so building it
    // from the probe's f64 approximation loses nothing.
    let point = probe.approx();
    let reach = solid_reach(solid, point);
    let far = point + direction * reach;
    let mut crossings = 0usize;
    for face in solid.faces() {
        match ray_crosses_face(face, probe, far)? {
            RayHit::Miss => {}
            RayHit::Crosses => crossings += 1,
            RayHit::OnFace => return Some(Containment::OnBoundary),
        }
    }
    Some(if crossings % 2 == 1 {
        Containment::Inside
    } else {
        Containment::Outside
    })
}

/// A length that certainly carries a ray from `point` clear of `solid`.
fn solid_reach(solid: &Polyhedron, point: Point3) -> f64 {
    let mut furthest: f64 = 1.0;
    for face in solid.faces() {
        for &v in face {
            furthest = furthest.max((v - point).length());
        }
    }
    // Doubling leaves the far endpoint strictly outside even when the
    // furthest vertex lies exactly along the probe direction.
    furthest * 2.0
}

/// Outcome of testing one ray against one face.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RayHit {
    Miss,
    Crosses,
    OnFace,
}

/// Whether the ray from the probe point to `far` crosses `face`.
///
/// Decided with orientation signs alone. The ray is represented by two
/// points on it, the probe and `far`; a crossing requires the face to
/// separate them, and the hit point to fall inside the face ring. Both
/// questions are sign tests, so no intersection coordinate is constructed.
fn ray_crosses_face(face: &[Point3], probe: &Probe, far: Point3) -> Option<RayHit> {
    let near_side = probe.side(face)?;
    let far_side = probe.side_of(face, far)?;

    if near_side == Sign::Zero {
        // The origin lies in the face plane: it may be ON the face.
        return if probe.ring(face)? == RingPosition::Outside {
            Some(RayHit::Miss)
        } else {
            Some(RayHit::OnFace)
        };
    }
    if near_side == far_side || far_side == Sign::Zero {
        // Both endpoints on one side, or the segment ends exactly in the
        // plane: extend the segment rather than deciding on a tangency.
        return Some(RayHit::Miss);
    }
    ray_enters_ring(face, probe, far)
}

/// Whether the segment `origin`-`far` passes through the face's interior.
///
/// For each ring edge, the tetrahedron (origin, far, edge start, edge end)
/// has a sign. The segment passes inside the ring exactly when every such
/// sign agrees. A zero sign means the segment meets an edge or vertex --
/// the degenerate case this refuses on rather than resolving arbitrarily.
fn ray_enters_ring(face: &[Point3], probe: &Probe, far: Point3) -> Option<RayHit> {
    let mut sign: Option<Sign> = None;
    for i in 0..face.len() {
        let a = face[i];
        let b = face[(i + 1) % face.len()];
        match probe.turn(far, a, b)? {
            Sign::Zero => return None,
            s => match sign {
                None => sign = Some(s),
                Some(previous) if previous == s => {}
                Some(_) => return Some(RayHit::Miss),
            },
        }
    }
    Some(RayHit::Crosses)
}

/// Where a coplanar point lies relative to the face ring.
///
/// The face is dropped to 2D by discarding its largest-normal-component
/// axis, which keeps the projection non-degenerate, and containment is then
/// decided by exact crossing parity using `orient2d`.
///
/// Parity is required rather than an all-same-side test: a same-side test
/// is only valid for CONVEX rings, and silently reports "outside" for any
/// point in the concave region of an L-shaped face. That failure is
/// invisible -- it makes coplanar contact go undetected, and the boolean
/// then keeps duplicate faces from both operands.
///
/// On an edge is reported apart from inside: a contact test treats the two
/// alike, but a classification point on its own fragment's edge is not an
/// interior point and must not be used.
fn ring_position(face: &[Point3], point: Point3) -> Option<RingPosition> {
    let normal = face_normal(face);
    let (nx, ny, nz) = (normal.x.abs(), normal.y.abs(), normal.z.abs());
    let flatten = |p: Point3| {
        if nx >= ny && nx >= nz {
            Point2::new(p.y, p.z)
        } else if ny >= nz {
            Point2::new(p.z, p.x)
        } else {
            Point2::new(p.x, p.y)
        }
    };

    let ring: Vec<Point2> = face.iter().map(|&v| flatten(v)).collect();
    let q = flatten(point);

    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        if orient2d(a, b, q).sign()? == Sign::Zero
            && q.x >= a.x.min(b.x)
            && q.x <= a.x.max(b.x)
            && q.y >= a.y.min(b.y)
            && q.y <= a.y.max(b.y)
        {
            return Some(RingPosition::OnEdge);
        }
    }

    let mut inside = false;
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        if (a.y > q.y) != (b.y > q.y) {
            // The edge straddles the horizontal through `q`; the crossing is
            // to the right exactly when the triangle orientation says so, so
            // no intersection abscissa is constructed.
            let sign = orient2d(a, b, q).sign()?;
            let upward = b.y > a.y;
            let right = if upward {
                sign == Sign::Negative
            } else {
                sign == Sign::Positive
            };
            if right {
                inside = !inside;
            }
        }
    }
    Some(if inside {
        RingPosition::Inside
    } else {
        RingPosition::Outside
    })
}

/// Whether a coplanar fragment's outward normal agrees with the opposing
/// face it lies in.
///
/// Two solids touching along a shared plane either face the same way (one
/// surface, keep a single copy) or face each other (the surfaces cancel).
/// Distinguishing them is what stops a duplicate face entering the shell.
/// The normals of coplanar faces are parallel, so the sign of their dot
/// product, taken exactly, decides it.
fn coplanar_normals_agree(
    fragment: &[Point3],
    probe: &Probe,
    other: &Polyhedron,
) -> GeomResult<bool> {
    for face in other.faces() {
        let on_plane = probe
            .side(face)
            .ok_or_else(|| unsupported("coplanar classification undecidable"))?;
        if on_plane != Sign::Zero {
            continue;
        }
        let position = probe
            .ring(face)
            .ok_or_else(|| unsupported("coplanar containment undecidable"))?;
        if position != RingPosition::Outside {
            return exact::normals_agree(fragment, face)
                .ok_or_else(|| unsupported("coplanar normals undecidable"));
        }
    }
    // No opposing face carries this fragment, so there is nothing to
    // duplicate and the fragment stands on its own.
    Ok(true)
}

/// Unnormalised outward normal of a face.
fn face_normal(face: &[Point3]) -> Vec3 {
    (face[1] - face[0]).cross(face[2] - face[0])
}

/// The two sides a polygon falls into when cut by a plane; `None` on a
/// side means the polygon does not reach it.
type SplitParts = (Option<Vec<Point3>>, Option<Vec<Point3>>);

/// Split a polygon by a plane, returning the negative and positive parts.
///
/// The plane is given by three points of an input face, never a derived one,
/// so the crossing points computed here are one step from input data. A
/// polygon lying wholly on one side comes back whole, so a non-crossing
/// plane costs nothing and introduces no vertices.
fn split_polygon(polygon: &[Point3], plane: &[Point3]) -> Option<SplitParts> {
    let mut signs = Vec::with_capacity(polygon.len());
    for &v in polygon {
        signs.push(side_of_face(plane, v)?);
    }
    let has_negative = signs.contains(&Sign::Negative);
    let has_positive = signs.contains(&Sign::Positive);
    if !has_positive {
        return Some((Some(polygon.to_vec()), None));
    }
    if !has_negative {
        return Some((None, Some(polygon.to_vec())));
    }

    let mut negative = Vec::new();
    let mut positive = Vec::new();
    for i in 0..polygon.len() {
        let j = (i + 1) % polygon.len();
        let (vi, vj) = (polygon[i], polygon[j]);
        let (si, sj) = (signs[i], signs[j]);
        match si {
            Sign::Negative => negative.push(vi),
            Sign::Positive => positive.push(vi),
            Sign::Zero => {
                negative.push(vi);
                positive.push(vi);
            }
            _ => {}
        }
        let crosses = matches!(
            (si, sj),
            (Sign::Negative, Sign::Positive) | (Sign::Positive, Sign::Negative)
        );
        if crosses {
            let cut = plane_crossing(plane, vi, vj)?;
            negative.push(cut);
            positive.push(cut);
        }
    }
    Some((
        (negative.len() >= 3).then_some(negative),
        (positive.len() >= 3).then_some(positive),
    ))
}

/// Where segment `a`-`b` meets the plane through `plane`'s first 3 points.
///
/// This is the only place in the module that constructs a coordinate. It is
/// stored in f64 (ADR 0045), each coordinate the double nearest to the exact
/// crossing (#199): see `polyhedron_exact::plane_crossing`. The f64 formula
/// below only seeds that search.
fn plane_crossing(plane: &[Point3], a: Point3, b: Point3) -> Option<Point3> {
    let normal = face_normal(plane);
    let t = normal.dot(plane[0] - a) / normal.dot(b - a);
    exact::plane_crossing(plane, a, b, a + (b - a) * t)
}

/// Exact boolean over two planar-faced solids.
///
/// Each operand's faces are split against every plane of the other, so no
/// fragment straddles the other solid's boundary. Each fragment is then kept
/// or dropped by classifying one point of its relative interior, and
/// difference reverses the tool fragments so the result stays outward-wound.
///
/// Refuses rather than guessing whenever a certified predicate cannot decide
/// a classification. A refusal is a typed error, never an approximate mesh.
pub fn boolean_polyhedra_exact(
    subject: &Polyhedron,
    tool: &Polyhedron,
    op: BooleanOp,
) -> GeomResult<Polyhedron> {
    let subject_parts = split_all(subject.faces(), tool.faces())?;
    let tool_parts = split_all(tool.faces(), subject.faces())?;

    let mut faces = Vec::new();
    for fragment in subject_parts {
        let probe = Probe::of(&fragment)?;
        let keep = match classify_fragment(&probe, tool)? {
            Containment::Inside => matches!(op, BooleanOp::Intersection),
            Containment::Outside => matches!(op, BooleanOp::Union | BooleanOp::Difference),
            // Coplanar contact: this fragment lies IN the tool's surface, so
            // both operands carry a copy. Exactly one must survive or the
            // shell gains a duplicate face and stops being manifold.
            //
            // Keeping the subject's copy is only correct when the two faces
            // agree on which side is solid. When their outward normals
            // OPPOSE, the surfaces cancel: an intersection there has zero
            // thickness, and a union has interior contact, so neither keeps
            // a face. That distinction is what the tool-side loop cannot
            // make, which is why it is made here.
            // Coplanar contact. Both operands carry a copy of this surface,
            // so exactly one must survive or the shell gains a duplicate
            // face -- which reads as a self-intersection, not as a
            // manifold error, because the duplicate is geometrically
            // coincident rather than topologically loose.
            //
            // The tool-side loop drops all its boundary fragments, so the
            // subject's copy is the survivor whenever the two normals
            // agree. When they OPPOSE, the surfaces are interior contact:
            // union and intersection both drop them, and difference keeps
            // the subject's copy because that face becomes the cut wall.
            Containment::OnBoundary => {
                if coplanar_normals_agree(&fragment, &probe, tool)? {
                    !matches!(op, BooleanOp::Difference)
                } else {
                    matches!(op, BooleanOp::Difference)
                }
            }
        };
        if keep {
            faces.push(fragment);
        }
    }
    for fragment in tool_parts {
        let containment = classify_fragment(&Probe::of(&fragment)?, subject)?;
        // A tool fragment on the subject's boundary is the same surface the
        // subject loop already kept, so it is always dropped here.
        let keep = match op {
            BooleanOp::Union => containment == Containment::Outside,
            BooleanOp::Intersection | BooleanOp::Difference => containment == Containment::Inside,
        };
        if keep {
            // Difference turns the tool's surface into an inward-facing
            // cavity wall, so its winding must flip to stay outward.
            faces.push(if op == BooleanOp::Difference {
                fragment.into_iter().rev().collect()
            } else {
                fragment
            });
        }
    }

    if faces.len() < 4 {
        return Err(unsupported("boolean produced no closed solid"));
    }
    Polyhedron::new(faces)
}

/// Split every face against every plane of the other solid.
fn split_all(faces: &[Vec<Point3>], planes: &[Vec<Point3>]) -> GeomResult<Vec<Vec<Point3>>> {
    let mut current: Vec<Vec<Point3>> = faces.to_vec();
    for plane in planes {
        let mut next = Vec::with_capacity(current.len());
        for polygon in current {
            let (negative, positive) = split_polygon(&polygon, plane).ok_or_else(|| {
                unsupported("face not splittable exactly against an operand plane")
            })?;
            next.extend(negative);
            next.extend(positive);
        }
        current = next;
    }
    Ok(current)
}

/// Classify a fragment at its probe point.
///
/// After splitting, a fragment lies wholly inside or wholly outside the other
/// solid, so one point of its relative interior decides for the whole
/// fragment. That point landing exactly on the boundary means the fragment
/// is coplanar with an opposing face -- handled by its own arm rather than
/// resolved arbitrarily.
fn classify_fragment(probe: &Probe, other: &Polyhedron) -> GeomResult<Containment> {
    // A degenerate ray is an unlucky direction, not an unanswerable point:
    // containment is the same along every ray, so try the next direction
    // rather than refusing. Each attempt is exact; none perturbs coordinates.
    for direction in probe_directions() {
        if let Some(containment) = contains(other, probe, direction) {
            return Ok(containment);
        }
    }
    // Every direction in the family was degenerate. That is vanishingly
    // unlikely for real geometry, and refusing remains correct: guessing a
    // parity here would silently produce a wrong solid.
    Err(unsupported(
        "every probe direction met a vertex or edge exactly",
    ))
}

/// Average of a polygon's vertices.
///
/// A coordinate every vertex shares is returned as that coordinate rather
/// than as a rounded average of copies of it, so the centroid of an
/// axis-aligned face lies exactly in the face's plane and qualifies as a
/// float probe point.
fn centroid_of(polygon: &[Point3]) -> Point3 {
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    for &v in polygon {
        sum += v - Point3::new(0.0, 0.0, 0.0);
    }
    let mean = Point3::new(0.0, 0.0, 0.0) + sum / polygon.len() as f64;
    let first = polygon[0];
    let shared = |axis: fn(Point3) -> f64, averaged: f64| {
        if polygon.iter().all(|&v| axis(v) == axis(first)) {
            axis(first)
        } else {
            averaged
        }
    };
    Point3::new(
        shared(|p| p.x, mean.x),
        shared(|p| p.y, mean.y),
        shared(|p| p.z, mean.z),
    )
}

/// Ray directions tried in order when classifying a point.
///
/// Containment does not depend on the probe direction: a closed orientable
/// solid has the same inside/outside answer along every ray. So a ray that
/// meets a vertex or edge exactly is not an unanswerable input, only an
/// unlucky one, and trying another direction is exact rather than a fudge.
///
/// The family is fixed, not random, so the same input gives the same answer
/// on every run. The first entry is the long-standing direction, so inputs
/// that already worked keep taking the same path. The rest are chosen to be
/// mutually non-parallel with irrational-ish ratios, which is what keeps them
/// from lining up with the axis-aligned and diagonal features that made the
/// first one degenerate.
fn probe_directions() -> [Vec3; 4] {
    [
        Vec3::new(0.577_215_664_9, 0.313_724_518_3, 0.144_729_885_8),
        Vec3::new(0.211_324_865_4, 0.788_675_134_6, 0.366_025_403_8),
        Vec3::new(0.867_513_459_5, 0.132_486_540_5, 0.539_189_129_1),
        Vec3::new(0.404_508_497_2, 0.595_491_502_8, 0.951_056_516_3),
    ]
}

/// Triangulate a polyhedron for measurement and diagnosis.
///
/// Vertices are shared through exact-coordinate keying: emitting a fresh
/// vertex per face would leave every edge used once, so an audit would
/// report a cloud of boundary edges for a solid that is in fact closed.
/// Coordinates that meet do so bit-identically, because they come from the
/// same literal or the same split, so exact keying is correct and no welding
/// tolerance is invented.
///
/// Fanning assumes convex rings. A non-convex face fans into triangles that
/// leave the footprint, so callers measuring such a solid must supply a
/// closed-form oracle instead.
#[must_use]
pub fn triangulate(solid: &Polyhedron) -> TriMesh {
    let mut positions: Vec<Point3> = Vec::new();
    let mut indices = Vec::new();
    let mut lookup: BTreeMap<[u64; 3], u32> = BTreeMap::new();
    for face in solid.faces() {
        let ring: Vec<u32> = face
            .iter()
            .map(|&p| {
                // `+ 0.0` keys -0.0 as 0.0, the same corner (#269).
                let key = [
                    (p.x + 0.0).to_bits(),
                    (p.y + 0.0).to_bits(),
                    (p.z + 0.0).to_bits(),
                ];
                let next = u32::try_from(positions.len()).unwrap_or(u32::MAX);
                *lookup.entry(key).or_insert_with(|| {
                    positions.push(p);
                    next
                })
            })
            .collect();
        for i in 1..ring.len() - 1 {
            indices.extend([ring[0], ring[i], ring[i + 1]]);
        }
    }
    TriMesh::new(positions, indices)
}
