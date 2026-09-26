//! Exact intersection curves for elementary surface pairs.
//!
//! A traced-and-fitted spline is an approximation with an error bound. For
//! the surface pairs whose intersection has a closed-form conic or linear
//! answer, no fitting is needed: the curve is derived symbolically from the
//! operands and is exact in the same sense as the rest of the exact B-rep
//! path.
//!
//! This module covers only those pairs, and refuses everything else rather
//! than falling back to approximation. Each derivation states the identity
//! it relies on, so a reader can check the algebra rather than trust it.

use axiolid_core::{Frame3, Interval, Point3, Scalar, Vec3};
use axiolid_curve::{Circle3, Curve3, Ellipse3};
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;
use axiolid_surface::{Cylinder, Plane, Sphere, Surface};

// --- exact decisions ---------------------------------------------------------
//
// Which closed form applies (tangent or crossing, parallel or oblique,
// perpendicular or tilted) is a sign question about the operands' own
// doubles. It is decided here in exact dyadic arithmetic, so a tangency
// that holds exactly for the given numbers is never mistaken for a tiny
// crossing, and vice versa. The constructed curves are still `f64`: their
// sizes are rounded once, from exact numerators where that is cheap.

type D3 = [Dyadic; 3];

/// An exact copy of a vector; `None` when a component is not finite.
fn exact3(v: Vec3) -> Option<D3> {
    Some([
        Dyadic::try_from_f64(v.x)?,
        Dyadic::try_from_f64(v.y)?,
        Dyadic::try_from_f64(v.z)?,
    ])
}

fn edot(a: &D3, b: &D3) -> Dyadic {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}

fn esub(a: &D3, b: &D3) -> D3 {
    [a[0].sub(&b[0]), a[1].sub(&b[1]), a[2].sub(&b[2])]
}

fn ecross_is_zero(a: &D3, b: &D3) -> bool {
    let c = [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ];
    c.iter().all(|v| esign(v) == Sign::Zero)
}

fn esign(v: &Dyadic) -> Sign {
    v.sign().expect("dyadic signs are always decided")
}

/// `r^2 |n|^2 - (n . (c - o))^2`, exactly: positive when the point `c` lies
/// closer than `r` to the plane through `o` with (non-unit) normal `n`,
/// zero when exactly at distance `r`. Also returns `|n|^2`.
fn within_radius(
    radius: Scalar,
    normal: Vec3,
    point: Point3,
    plane_origin: Point3,
) -> Result<(Dyadic, Dyadic), ExactIntersectionRefusal> {
    let bad = ExactIntersectionRefusal::DegenerateFrame;
    let n = exact3(normal).ok_or(bad.clone())?;
    let c = exact3(point).ok_or(bad.clone())?;
    let o = exact3(plane_origin).ok_or(bad.clone())?;
    let r = Dyadic::try_from_f64(radius).ok_or(bad.clone())?;
    let nn = edot(&n, &n);
    if esign(&nn) == Sign::Zero {
        return Err(bad);
    }
    let nd = edot(&n, &esub(&c, &o));
    Ok((r.square().mul(&nn).sub(&nd.square()), nn))
}

/// `numerator / nn` rounded once, for a size whose square is known exactly.
fn rounded_square(numerator: &Dyadic, nn: &Dyadic) -> Result<Scalar, ExactIntersectionRefusal> {
    let value = numerator.to_f64() / nn.to_f64();
    if value > 0.0 && value.is_finite() {
        Ok(value)
    } else {
        Err(ExactIntersectionRefusal::DegenerateFrame)
    }
}

/// Why an elementary pair has no exact closed-form intersection curve here.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ExactIntersectionRefusal {
    /// The pair is not one of the supported elementary combinations.
    ///
    /// Not a statement about the geometry: the intersection may well be a
    /// nameable curve, just not one this module derives.
    UnsupportedPair,
    /// The surfaces are parallel or concentric and do not meet at all.
    Disjoint,
    /// The surfaces coincide or touch tangentially, so the intersection is
    /// not a regular curve.
    ///
    /// A single tangential point or a shared surface patch cannot be
    /// returned as a curve without inventing structure.
    NotRegularCurve,
    /// A required frame axis was degenerate, so no exact frame can be built.
    DegenerateFrame,
    /// The intersection is a parabola or hyperbola, which `Curve3` has no
    /// variant for.
    ///
    /// The curve is perfectly well defined and exactly derivable; it simply
    /// cannot be represented without adding a conic variant. Refusing names
    /// that representational gap instead of substituting a nearby ellipse or
    /// a fitted spline.
    UnrepresentableConic,
}

/// The exact intersection curve of two elementary surfaces, when one exists
/// in closed form.
///
/// Returns the curve together with the identity used to derive it, so a
/// caller can record provenance rather than re-deriving trust.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ExactIntersectionCurve {
    /// The branches of the intersection, in deterministic order.
    ///
    /// Usually one, but several elementary pairs genuinely meet in TWO
    /// disjoint components: equal-radius cylinders on intersecting axes
    /// cut two ellipses, and parallel cylinders cut two lines. Returning a
    /// single curve would have forced this code to pick one and discard the
    /// other, which is exactly the silent geometry loss the rest of this
    /// module refuses to do. Exact: every coordinate comes from the
    /// operands' own numbers through the stated identity, never from a fit.
    pub branches: Vec<Curve3>,
    /// Which closed-form identity produced `branches`.
    pub derivation: Derivation,
    /// The parameter span each branch exists on, aligned with `branches`:
    /// `None` for a curve defined on its whole natural domain (a line, a
    /// full circle or ellipse), `Some` for a piece of a ruled section
    /// (ADR 0076), which exists only where its discriminant is not
    /// negative. Two pieces over one span join at both ends into a loop.
    pub spans: Vec<Option<Interval>>,
}

impl ExactIntersectionCurve {
    /// Branches on their whole natural domains.
    pub(crate) fn whole(branches: Vec<Curve3>, derivation: Derivation) -> Self {
        let spans = vec![None; branches.len()];
        Self {
            branches,
            derivation,
            spans,
        }
    }

    /// Branches with explicit spans.
    pub(crate) fn with_spans(
        branches: Vec<Curve3>,
        spans: Vec<Option<Interval>>,
        derivation: Derivation,
    ) -> Self {
        Self {
            branches,
            derivation,
            spans,
        }
    }

    /// The sole branch, when the caller expects exactly one.
    ///
    /// Panics when there are several: a caller that assumes one branch and
    /// silently sees only the first would lose geometry, so this fails
    /// loudly instead.
    pub fn single(&self) -> &Curve3 {
        assert_eq!(self.branches.len(), 1, "expected one branch");
        &self.branches[0]
    }
}

/// The closed-form identity behind an exact intersection curve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Derivation {
    /// Two non-parallel planes meet in a line along `n1 x n2`.
    PlanePlaneLine,
    /// A plane perpendicular to a cylinder axis cuts a circle of the
    /// cylinder's own radius.
    CylinderPlanePerpendicularCircle,
    /// A plane oblique to a cylinder axis cuts an ellipse with semi-axes
    /// `r` and `r / cos(theta)`.
    CylinderPlaneObliqueEllipse,
    /// A plane at signed distance `d` from a sphere centre cuts a circle of
    /// radius `sqrt(r^2 - d^2)`.
    SpherePlaneCircle,
    /// Two spheres meet in a circle in their radical plane.
    SphereSphereCircle,
    /// Equal-radius cylinders on intersecting axes cut two ellipses.
    CylinderCylinderSteinmetzEllipses,
    /// Cylinders with parallel axes meet in one or two axis-parallel lines.
    ParallelCylinderLines,
    /// A plane parallel to a cylinder axis cuts one or two rulings.
    CylinderPlaneParallelRulings,
    /// Two coaxial surfaces of revolution meet in circles perpendicular to
    /// the shared axis, found by intersecting their meridian profiles.
    CoaxialRevolutionCircles,
    /// A quadric substituted into a ruled carrier (cylinder, elliptical
    /// cylinder, or a cone cut by a plane) is quadratic in the ruling
    /// parameter at every angle; the curve is a root branch of that
    /// quadratic (ADR 0076).
    RuledQuadricSection,
    /// A plane or sphere meets each circle of a torus about its axis where
    /// `A(v) cos u + B(v) sin u = C(v)`; the curve is `u` as a function of
    /// the tube angle `v` (ADR 0076).
    TorusAngleSection,
}

/// Derive the exact intersection curve of two elementary surfaces.
///
/// Returns `Err` with an explicit refusal for every pair this module does
/// not derive in closed form. Refusal is never a fallback to approximation:
/// a caller that needs those cases must use the certified numeric analysis
/// and decide for itself what to do with an unproven region.
pub fn exact_surface_intersection(
    first: &Surface,
    second: &Surface,
) -> Result<ExactIntersectionCurve, ExactIntersectionRefusal> {
    match closed_form(first, second) {
        Ok(curve) => Ok(curve),
        // Apart, or a frame that cannot be read: final either way.
        Err(
            refusal @ (ExactIntersectionRefusal::Disjoint
            | ExactIntersectionRefusal::DegenerateFrame),
        ) => Err(refusal),
        // No conic or line for this pair: a cylinder or cone cut by a
        // quadric is still exact as a ruled section (ADR 0076).
        Err(refusal) => {
            if let Some(curve) = crate::ruled_section::ruled_section(first, second)? {
                return Ok(curve);
            }
            match crate::torus_section::torus_section(first, second)? {
                Some(curve) => Ok(curve),
                None => Err(refusal),
            }
        }
    }
}

/// The pairs with a line, circle or ellipse in closed form.
fn closed_form(
    first: &Surface,
    second: &Surface,
) -> Result<ExactIntersectionCurve, ExactIntersectionRefusal> {
    match (first, second) {
        (Surface::Plane(a), Surface::Plane(b)) => plane_plane(a, b),
        (Surface::Cylinder(c), Surface::Plane(p)) => cylinder_plane(c, p),
        (Surface::Plane(p), Surface::Cylinder(c)) => cylinder_plane(c, p),
        (Surface::Sphere(s), Surface::Plane(p)) => sphere_plane(s, p),
        (Surface::Plane(p), Surface::Sphere(s)) => sphere_plane(s, p),
        // A plane cutting a cone gives a conic whose kind depends on the
        // tilt: circle and ellipse are representable, parabola and
        // hyperbola are not. Coaxial (tilt 0) is handled by the shared
        // revolution path; any other tilt is refused by kind.
        (Surface::Cone(c), Surface::Plane(p)) | (Surface::Plane(p), Surface::Cone(c)) => {
            cone_plane(c, p)
        }
        (Surface::Sphere(a), Surface::Sphere(b)) => sphere_sphere(a, b),
        (Surface::Cylinder(a), Surface::Cylinder(b)) => cylinder_cylinder(a, b),
        // Every remaining elementary pair is covered by the shared
        // surface-of-revolution identity when the axes coincide.
        _ => crate::revolution_profile::coaxial_revolution_intersection(first, second),
    }
}

/// Two planes meet in a line, unless their normals are parallel.
///
/// Identity: the direction is `n1 x n2`. With `di = ni . oi`, the point
/// `p = (d1 (n2 x dir) + d2 (dir x n1)) / |dir|^2` satisfies both plane
/// equations and lies nearest the origin, so it is a canonical choice that
/// depends only on the operands.
fn plane_plane(
    first: &Plane,
    second: &Plane,
) -> Result<ExactIntersectionCurve, ExactIntersectionRefusal> {
    let first_normal = first.frame.z;
    let second_normal = second.frame.z;
    let direction = first_normal.cross(second_normal);
    let direction_squared = direction.dot(direction);
    if direction_squared == 0.0 {
        // Parallel normals: the planes either coincide or never meet.
        // Neither outcome is a regular curve.
        return Err(ExactIntersectionRefusal::Disjoint);
    }
    let first_offset = first_normal.dot(first.frame.origin);
    let second_offset = second_normal.dot(second.frame.origin);
    let origin = (second_normal.cross(direction) * first_offset
        + direction.cross(first_normal) * second_offset)
        / direction_squared;
    let unit_direction = direction / direction_squared.sqrt();
    if !origin.is_finite() || !unit_direction.is_finite() {
        return Err(ExactIntersectionRefusal::DegenerateFrame);
    }
    Ok(ExactIntersectionCurve {
        spans: vec![None],
        branches: vec![Curve3::Line(axiolid_curve::Line3 {
            origin,
            direction: unit_direction,
        })],
        derivation: Derivation::PlanePlaneLine,
    })
}

/// A plane cuts a sphere in a circle.
///
/// Identity: with `d` the signed distance from the centre to the plane, the
/// section has radius `sqrt(r^2 - d^2)` and is centred at the centre's
/// projection onto the plane. `|d| >= r` is refused: `|d| > r` misses the
/// sphere entirely and `|d| == r` touches at one point, which is not a
/// regular curve.
fn sphere_plane(
    sphere: &Sphere,
    plane: &Plane,
) -> Result<ExactIntersectionCurve, ExactIntersectionRefusal> {
    let normal = plane.frame.z;
    let normal_squared = normal.dot(normal);
    if normal_squared == 0.0 {
        return Err(ExactIntersectionRefusal::DegenerateFrame);
    }
    let unit_normal = normal / normal_squared.sqrt();
    let centre = sphere.frame.origin;
    let signed_distance = unit_normal.dot(centre - plane.frame.origin);
    // Tangent or not is decided exactly; in `f64` a plane exactly tangent
    // along a normal like (3, 2, 6) came out as a circle of radius 1e-7.
    let (numerator, nn) = within_radius(sphere.radius, normal, centre, plane.frame.origin)?;
    match esign(&numerator) {
        Sign::Positive => {}
        // Exactly tangent: a touch is a point, not a curve.
        Sign::Zero => return Err(ExactIntersectionRefusal::NotRegularCurve),
        _ => return Err(ExactIntersectionRefusal::Disjoint),
    }
    let radius_squared = rounded_square(&numerator, &nn)?;
    let section_centre = centre - unit_normal * signed_distance;
    let frame = frame_from_normal(section_centre, unit_normal)?;
    Ok(ExactIntersectionCurve {
        spans: vec![None],
        branches: vec![Curve3::Circle(Circle3 {
            frame,
            radius: radius_squared.sqrt(),
        })],
        derivation: Derivation::SpherePlaneCircle,
    })
}

/// A plane cuts an infinite cylinder in a circle or an ellipse.
///
/// Identities, with `theta` the angle between the plane normal and the
/// cylinder axis:
/// - `theta == 0` (plane perpendicular to the axis): a circle of radius `r`.
/// - `0 < theta < pi/2`: an ellipse with minor semi-axis `r` across the
///   axis and major semi-axis `r / cos(theta)` along the tilt direction.
/// - `theta == pi/2` (plane parallel to the axis): refused. The section is
///   then two parallel lines, one line, or empty, none of which is a single
///   regular curve.
fn cylinder_plane(
    cylinder: &Cylinder,
    plane: &Plane,
) -> Result<ExactIntersectionCurve, ExactIntersectionRefusal> {
    let axis_squared = cylinder.frame.z.dot(cylinder.frame.z);
    let normal_squared = plane.frame.z.dot(plane.frame.z);
    if axis_squared == 0.0 || normal_squared == 0.0 {
        return Err(ExactIntersectionRefusal::DegenerateFrame);
    }
    let axis = cylinder.frame.z / axis_squared.sqrt();
    let normal = plane.frame.z / normal_squared.sqrt();
    // Parallel and perpendicular are decided exactly on the given axes: in
    // `f64` an axis (-3, -3, -3) against a normal (-3, 1, 2), whose dot
    // product is exactly zero, gave cos = 2.8e-17 and a 10^16-long ellipse.
    let bad = ExactIntersectionRefusal::DegenerateFrame;
    let exact_axis = exact3(cylinder.frame.z).ok_or(bad.clone())?;
    let exact_normal = exact3(plane.frame.z).ok_or(bad.clone())?;
    let along = edot(&exact_axis, &exact_normal);
    if esign(&along) == Sign::Zero {
        // Plane parallel to the axis: the section is a pair of rulings, one
        // ruling when the plane is tangent, or empty. Each is an exact line,
        // so this is derived rather than refused.
        return cylinder_plane_parallel(cylinder, plane, axis, normal);
    }
    let perpendicular = ecross_is_zero(&exact_axis, &exact_normal);
    // cos(theta) between axis and plane normal, from the exact dot product
    // rounded once. A tilted plane keeps a cosine below 1 even when its
    // tilt is below `f64` resolution, so it is still reported as the
    // ellipse it is.
    let cosine = if perpendicular {
        1.0
    } else {
        let scale = (edot(&exact_axis, &exact_axis).to_f64()
            * edot(&exact_normal, &exact_normal).to_f64())
        .sqrt();
        (along.to_f64().abs() / scale).min(1.0_f64.next_down())
    };
    // Non-zero exactly, but it can round to zero or NaN for extreme inputs.
    if cosine.is_nan() || cosine <= 0.0 {
        return Err(bad);
    }
    // The section centre is where the cylinder axis pierces the plane.
    let axis_origin = cylinder.frame.origin;
    let to_plane = normal.dot(plane.frame.origin - axis_origin);
    let centre = axis_origin + axis * (to_plane / axis.dot(normal));
    if !centre.is_finite() {
        return Err(ExactIntersectionRefusal::DegenerateFrame);
    }
    Ok(ExactIntersectionCurve {
        spans: vec![None],
        branches: vec![cylinder_section_curve(
            centre,
            axis,
            normal,
            cylinder.radius,
            cosine,
        )?],
        derivation: if perpendicular {
            Derivation::CylinderPlanePerpendicularCircle
        } else {
            Derivation::CylinderPlaneObliqueEllipse
        },
    })
}

/// Build the circle or ellipse a plane cuts from a cylinder.
///
/// The minor axis lies along `axis x normal`, which is perpendicular to the
/// tilt and so always spans the cylinder at its own radius. The major axis
/// completes the frame and is stretched by `1 / cos(theta)`.
fn cylinder_section_curve(
    centre: Point3,
    axis: Vec3,
    normal: Vec3,
    radius: Scalar,
    cosine: Scalar,
) -> Result<Curve3, ExactIntersectionRefusal> {
    if cosine == 1.0 {
        let frame = frame_from_normal(centre, normal)?;
        return Ok(Curve3::Circle(Circle3 { frame, radius }));
    }
    let across = axis.cross(normal);
    let across_squared = across.dot(across);
    if across_squared == 0.0 {
        return Err(ExactIntersectionRefusal::DegenerateFrame);
    }
    let minor = across / across_squared.sqrt();
    let major = normal.cross(minor);
    if !minor.is_finite() || !major.is_finite() {
        return Err(ExactIntersectionRefusal::DegenerateFrame);
    }
    let frame = Frame3 {
        origin: centre,
        x: minor,
        y: major,
        z: normal,
    };
    Ok(Curve3::Ellipse(Ellipse3 {
        frame,
        semi_axis_x: radius,
        semi_axis_y: radius / cosine,
    }))
}

/// Build an orthonormal frame whose `z` is the given unit normal.
///
/// The in-plane axes are otherwise arbitrary, so they are chosen
/// deterministically from the normal's own components: pick the coordinate
/// axis least aligned with the normal as a seed. A deterministic choice
/// matters because the frame ends up in the returned curve, and an
/// orientation that varied run to run would make results irreproducible.
pub(crate) fn frame_from_normal(
    origin: Point3,
    normal: Vec3,
) -> Result<Frame3, ExactIntersectionRefusal> {
    let seed = if normal.x.abs() <= normal.y.abs() && normal.x.abs() <= normal.z.abs() {
        Vec3::new(1.0, 0.0, 0.0)
    } else if normal.y.abs() <= normal.z.abs() {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let x_axis = normal.cross(seed);
    let x_squared = x_axis.dot(x_axis);
    if x_squared == 0.0 {
        return Err(ExactIntersectionRefusal::DegenerateFrame);
    }
    let x = x_axis / x_squared.sqrt();
    let y = normal.cross(x);
    if !x.is_finite() || !y.is_finite() {
        return Err(ExactIntersectionRefusal::DegenerateFrame);
    }
    Ok(Frame3 {
        origin,
        x,
        y,
        z: normal,
    })
}

/// Two spheres meet in a circle lying in their radical plane.
///
/// Identity: with `d` the centre distance and
/// `a = (d^2 + r1^2 - r2^2) / (2 d)` the distance from the first centre
/// along the centre line, the circle has radius `sqrt(r1^2 - a^2)` and is
/// centred at `c1 + a * u`. Both are built from the operands' own numbers.
fn sphere_sphere(
    first: &Sphere,
    second: &Sphere,
) -> Result<ExactIntersectionCurve, ExactIntersectionRefusal> {
    let separation = second.frame.origin - first.frame.origin;
    let distance = separation.length();
    if distance == 0.0 {
        // Concentric: identical spheres coincide, otherwise they never meet.
        return Err(if first.radius == second.radius {
            ExactIntersectionRefusal::NotRegularCurve
        } else {
            ExactIntersectionRefusal::Disjoint
        });
    }
    if distance > first.radius + second.radius {
        return Err(ExactIntersectionRefusal::Disjoint);
    }
    if distance < (first.radius - second.radius).abs() {
        // One sphere strictly encloses the other.
        return Err(ExactIntersectionRefusal::Disjoint);
    }
    let axis = separation / distance;
    let along = (distance * distance + first.radius * first.radius - second.radius * second.radius)
        / (2.0 * distance);
    let squared = first.radius * first.radius - along * along;
    if squared <= 0.0 {
        // Tangent spheres touch at a single point, which is not a curve.
        return Err(ExactIntersectionRefusal::NotRegularCurve);
    }
    let centre = first.frame.origin + axis * along;
    let frame = frame_from_normal(centre, axis)?;
    Ok(ExactIntersectionCurve {
        spans: vec![None],
        branches: vec![Curve3::Circle(Circle3 {
            frame,
            radius: squared.sqrt(),
        })],
        derivation: Derivation::SphereSphereCircle,
    })
}

/// Two cylinders, where the pair has a closed-form conic decomposition.
///
/// Only two configurations decompose into conics:
/// - parallel axes: the cross-section is two circles meeting in at most
///   two points, so the intersection is one or two lines along the axis;
/// - intersecting axes with EQUAL radii: the Steinmetz case, where the
///   solid identity `|p-P|^2 - (a.(p-P))^2 = |p-P|^2 - (b.(p-P))^2`
///   factors into the two planes with normals `a-b` and `a+b`, each
///   cutting the first cylinder in an ellipse.
///
/// Every other configuration -- notably unequal radii on intersecting
/// axes -- is a genuine space quartic that is NOT planar and therefore
/// has no exact conic representation. Verified numerically before this
/// was written: the third singular value of sampled points is 2.83, not
/// ~0, so no plane contains the curve. Those refuse rather than being
/// approximated by a fitted spline.
fn cylinder_cylinder(
    first: &Cylinder,
    second: &Cylinder,
) -> Result<ExactIntersectionCurve, ExactIntersectionRefusal> {
    let first_axis = first.frame.z;
    let second_axis = second.frame.z;
    let cross = first_axis.cross(second_axis);
    if cross.length() == 0.0 {
        return parallel_cylinders(first, second, first_axis);
    }
    if first.radius != second.radius {
        return Err(ExactIntersectionRefusal::NotRegularCurve);
    }
    // Steinmetz needs the axes to actually meet. Skew axes of equal
    // radius still give a quartic, so the common point is required, not
    // assumed: the shortest connecting segment must have zero length.
    let between = second.frame.origin - first.frame.origin;
    let unit_cross = cross / cross.length();
    if between.dot(unit_cross) != 0.0 {
        return Err(ExactIntersectionRefusal::NotRegularCurve);
    }
    // Solve for the crossing point on the first axis.
    let denominator = first_axis
        .dot(second_axis)
        .mul_add(-first_axis.dot(second_axis), 1.0);
    if denominator == 0.0 {
        return Err(ExactIntersectionRefusal::NotRegularCurve);
    }
    let along = (between.dot(first_axis) - first_axis.dot(second_axis) * between.dot(second_axis))
        / denominator;
    let meeting = first.frame.origin + first_axis * along;
    let mut branches = Vec::new();
    for normal in [first_axis - second_axis, first_axis + second_axis] {
        let length = normal.length();
        if length == 0.0 {
            continue;
        }
        let unit_normal = normal / length;
        let cosine = unit_normal.dot(first_axis).abs();
        if cosine == 0.0 {
            return Err(ExactIntersectionRefusal::NotRegularCurve);
        }
        branches.push(cylinder_section_curve(
            meeting,
            first_axis,
            unit_normal,
            first.radius,
            cosine,
        )?);
    }
    Ok(ExactIntersectionCurve::whole(
        branches,
        Derivation::CylinderCylinderSteinmetzEllipses,
    ))
}

/// Cylinders with parallel axes meet in lines parallel to those axes.
///
/// Identity: reduce to the plane perpendicular to the shared axis, where
/// the problem is two circles. With `d` the perpendicular centre offset,
/// the circles meet where `x = (d^2 + r1^2 - r2^2) / 2d` along the offset
/// direction and `y = +/- sqrt(r1^2 - x^2)` across it. Each solution
/// lifts to a full line along the axis.
fn parallel_cylinders(
    first: &Cylinder,
    second: &Cylinder,
    axis: Vec3,
) -> Result<ExactIntersectionCurve, ExactIntersectionRefusal> {
    let between = second.frame.origin - first.frame.origin;
    let offset = between - axis * between.dot(axis);
    let distance = offset.length();
    if distance == 0.0 {
        return Err(if first.radius == second.radius {
            ExactIntersectionRefusal::NotRegularCurve
        } else {
            ExactIntersectionRefusal::Disjoint
        });
    }
    if distance > first.radius + second.radius || distance < (first.radius - second.radius).abs() {
        return Err(ExactIntersectionRefusal::Disjoint);
    }
    let toward = offset / distance;
    let along = (distance * distance + first.radius * first.radius - second.radius * second.radius)
        / (2.0 * distance);
    let squared = first.radius * first.radius - along * along;
    let base = first.frame.origin + toward * along;
    if squared <= 0.0 {
        // Tangent cylinders share exactly one line.
        return Ok(ExactIntersectionCurve {
            spans: vec![None],
            branches: vec![Curve3::Line(axiolid_curve::Line3 {
                origin: base,
                direction: axis,
            })],
            derivation: Derivation::ParallelCylinderLines,
        });
    }
    let across = axis.cross(toward);
    let half = squared.sqrt();
    let mut branches = Vec::new();
    for sign in [1.0, -1.0] {
        branches.push(Curve3::Line(axiolid_curve::Line3 {
            origin: base + across * (sign * half),
            direction: axis,
        }));
    }
    Ok(ExactIntersectionCurve::whole(
        branches,
        Derivation::ParallelCylinderLines,
    ))
}

/// A plane cuts a cone in a conic whose kind follows the tilt.
///
/// With `phi` the angle between the plane and the cone axis and `alpha`
/// the semi-angle, the section is an ellipse while `phi > alpha`, a
/// parabola at `phi == alpha`, and a hyperbola below. Only the
/// perpendicular case reduces to a circle, and that is the coaxial case
/// handled by the shared revolution identity.
///
/// `Curve3` has no parabola or hyperbola variant, so those kinds are
/// refused by name. The ellipse case is genuinely derivable in closed
/// form but needs the apex-offset construction rather than the
/// cylinder's parallel-axis one, so it is refused as unsupported until
/// that derivation is written and tested. Naming the two gaps
/// differently keeps 'not representable' distinct from 'not implemented'.
fn cone_plane(
    cone: &axiolid_surface::Cone,
    plane: &axiolid_surface::Plane,
) -> Result<ExactIntersectionCurve, ExactIntersectionRefusal> {
    let axis = cone.frame.z;
    let normal = plane.frame.z;
    let alignment = axis.dot(normal).abs().clamp(0.0, 1.0);
    // Angle between the plane and the axis, from the axis/normal angle.
    let plane_axis_angle = alignment.asin();
    let semi = cone.semi_angle.abs();
    let perpendicular = match (exact3(axis), exact3(normal)) {
        (Some(a), Some(n)) => ecross_is_zero(&a, &n),
        _ => return Err(ExactIntersectionRefusal::DegenerateFrame),
    };
    if perpendicular {
        // Perpendicular plane (decided exactly): a circle, via the coaxial
        // profile path.
        return crate::revolution_profile::coaxial_revolution_intersection(
            &Surface::Cone(*cone),
            &Surface::Plane(*plane),
        );
    }
    if plane_axis_angle > semi {
        return Err(ExactIntersectionRefusal::UnsupportedPair);
    }
    Err(ExactIntersectionRefusal::UnrepresentableConic)
}

/// A plane parallel to a cylinder axis cuts rulings, not a conic.
///
/// Identity: with `d` the distance from the axis to the plane and the
/// half-chord `h = sqrt(r^2 - d^2)`, the plane meets the cylinder in the
/// two lines through `foot +/- t*h` along the axis, where `foot` is the
/// axis point projected onto the plane and `t = axis x normal` is the unit
/// in-plane direction perpendicular to the axis. A tangent plane gives one
/// line; a plane clear of the cylinder gives none.
fn cylinder_plane_parallel(
    cylinder: &Cylinder,
    plane: &Plane,
    axis: Vec3,
    normal: Vec3,
) -> Result<ExactIntersectionCurve, ExactIntersectionRefusal> {
    let distance = normal.dot(cylinder.frame.origin - plane.frame.origin);
    // Two rulings, one (tangent) or none, decided exactly.
    let (numerator, nn) = within_radius(
        cylinder.radius,
        plane.frame.z,
        cylinder.frame.origin,
        plane.frame.origin,
    )?;
    let tangent_plane = match esign(&numerator) {
        Sign::Positive => false,
        Sign::Zero => true,
        _ => return Err(ExactIntersectionRefusal::Disjoint),
    };
    // `axis` and `normal` are unit and perpendicular here, so their cross
    // product is already unit: no second normalisation is needed.
    let tangent = axis.cross(normal);
    let foot = cylinder.frame.origin - normal * distance;
    if !foot.is_finite() || !tangent.is_finite() {
        return Err(ExactIntersectionRefusal::DegenerateFrame);
    }
    let half_chord = if tangent_plane {
        0.0
    } else {
        rounded_square(&numerator, &nn)?.sqrt()
    };
    let mut branches = Vec::new();
    let offsets: &[Scalar] = if tangent_plane {
        &[0.0]
    } else {
        &[half_chord, -half_chord]
    };
    branches
        .try_reserve_exact(offsets.len())
        .map_err(|_| ExactIntersectionRefusal::DegenerateFrame)?;
    for offset in offsets {
        branches.push(Curve3::Line(axiolid_curve::Line3 {
            origin: foot + tangent * *offset,
            direction: axis,
        }));
    }
    Ok(ExactIntersectionCurve::whole(
        branches,
        Derivation::CylinderPlaneParallelRulings,
    ))
}
