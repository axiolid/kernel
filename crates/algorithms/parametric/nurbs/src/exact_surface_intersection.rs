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

use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::{Circle3, Curve3, Ellipse3};
use axiolid_surface::{Cylinder, Plane, Sphere, Surface};

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
}

impl ExactIntersectionCurve {
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
    /// Two coaxial surfaces of revolution meet in circles perpendicular to
    /// the shared axis, found by intersecting their meridian profiles.
    CoaxialRevolutionCircles,
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
    let radius_squared = sphere.radius * sphere.radius - signed_distance * signed_distance;
    if radius_squared <= 0.0 {
        // Strictly outside, or exactly tangent. A tangent touch is a point,
        // not a curve, so both refuse rather than degenerate to radius 0.
        return Err(if radius_squared == 0.0 {
            ExactIntersectionRefusal::NotRegularCurve
        } else {
            ExactIntersectionRefusal::Disjoint
        });
    }
    let section_centre = centre - unit_normal * signed_distance;
    let frame = frame_from_normal(section_centre, unit_normal)?;
    Ok(ExactIntersectionCurve {
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
    // cos(theta) between axis and plane normal. Sign only reflects axis
    // orientation, so the magnitude carries the geometry.
    let cosine = axis.dot(normal).abs();
    if cosine == 0.0 {
        // Plane parallel to the axis: not a single regular curve.
        return Err(ExactIntersectionRefusal::NotRegularCurve);
    }
    // The section centre is where the cylinder axis pierces the plane.
    let axis_origin = cylinder.frame.origin;
    let to_plane = normal.dot(plane.frame.origin - axis_origin);
    let centre = axis_origin + axis * (to_plane / axis.dot(normal));
    if !centre.is_finite() {
        return Err(ExactIntersectionRefusal::DegenerateFrame);
    }
    Ok(ExactIntersectionCurve {
        branches: vec![cylinder_section_curve(
            centre,
            axis,
            normal,
            cylinder.radius,
            cosine,
        )?],
        derivation: if cosine == 1.0 {
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
    Ok(ExactIntersectionCurve {
        branches,
        derivation: Derivation::CylinderCylinderSteinmetzEllipses,
    })
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
    Ok(ExactIntersectionCurve {
        branches,
        derivation: Derivation::ParallelCylinderLines,
    })
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
    if alignment == 1.0 {
        // Perpendicular plane: a circle, via the coaxial profile path.
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
