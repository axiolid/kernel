//! The exit check of a result: every vertex lies on its curves (#445).
//!
//! A boolean's report says what the tolerance decided, not where the
//! result's vertices ended up. A slack that grows with a parameter can leave
//! a vertex further from its edge's curve than any reading admits while the
//! report stays exact. This check measures the result itself: each vertex
//! against the 3D curve of every edge it bounds, and against each edge use's
//! pcurve lifted through its face's surface at the use's two ends. The
//! largest residue must lie within [`bound`]: the caller's linear
//! tolerance, and for an exact report the rounding floor too. A result that
//! exceeds it is refused ([`BooleanError::VertexOffCurve`]), never returned
//! with a report that misdescribes it. The check runs inside the cut, so
//! the finer, extending and rounding-end retries apply to its refusal.
//!
//! At [`Tolerance::ZERO`] only the rounding floor bounds an exact report.
//! Far from the origin that floor exceeds the caller's tolerance (ADR 0080);
//! the check caps what reaches a vertex but cannot tighten the floor.

use axiolid_brep::ExactBRep;
use axiolid_core::{Interval, Point3, Scalar, Tolerance};
use axiolid_curve::{Curve2, Curve3};
use axiolid_evaluate::{evaluate2, evaluate3, surface};
use axiolid_surface::Surface;
use axiolid_topology::Orientation;

use crate::report;
use crate::BooleanError;

/// The largest residue a result may have: the caller's tolerance `eps`
/// capped by the rounding `floor` for an exact report, the floor alone at a
/// zero tolerance, and `eps` for a report that recorded readings.
pub(crate) fn bound(caller: Tolerance, floor: Scalar, exact: bool) -> Scalar {
    let eps = caller.linear();
    match (exact, eps > 0.0) {
        (true, true) => floor.min(eps),
        (true, false) => floor,
        (false, _) => eps,
    }
}

/// Whether a residue is over its bound; a NaN residue is.
fn exceeds(residue: Scalar, bound: Scalar) -> bool {
    residue.is_nan() || residue > bound
}

/// How far `start` and `end` lie from `curve` at the ends of `span`.
fn edge_residue(
    curve: &Curve3,
    span: Interval,
    start: Point3,
    end: Point3,
) -> Result<Scalar, BooleanError> {
    let mut worst: Scalar = 0.0;
    for (parameter, vertex) in [(span.start, start), (span.end, end)] {
        let point = evaluate3(curve, parameter).map_err(|_| BooleanError::Evaluation)?;
        worst = worst.max(distance(point, vertex)?);
    }
    Ok(worst)
}

/// How far `start` and `end` lie from `pcurve` lifted through `surface` at
/// the ends of `span`.
fn use_residue(
    pcurve: &Curve2,
    span: Interval,
    surface: &Surface,
    start: Point3,
    end: Point3,
) -> Result<Scalar, BooleanError> {
    let mut worst: Scalar = 0.0;
    for (parameter, vertex) in [(span.start, start), (span.end, end)] {
        let uv = evaluate2(pcurve, parameter).map_err(|_| BooleanError::Evaluation)?;
        let point = surface::evaluate(surface, uv.x, uv.y).map_err(|_| BooleanError::Evaluation)?;
        worst = worst.max(distance(point, vertex)?);
    }
    Ok(worst)
}

fn distance(a: Point3, b: Point3) -> Result<Scalar, BooleanError> {
    let d = (a - b).length();
    if d.is_finite() {
        Ok(d)
    } else {
        Err(BooleanError::Evaluation)
    }
}

/// The largest residue of a B-rep's vertices against its edge curves and
/// its edge uses' lifted pcurves.
///
/// # Errors
///
/// [`BooleanError::DanglingReference`] for a missing handle or interval,
/// [`BooleanError::Evaluation`] for a failed evaluation or a non-finite
/// distance.
pub(crate) fn residue(brep: &ExactBRep) -> Result<Scalar, BooleanError> {
    let topology = brep.topology();
    let position = |vertex: axiolid_topology::VertexId| {
        topology
            .vertices()
            .get(vertex.index())
            .map(|v| v.position)
            .ok_or(BooleanError::DanglingReference)
    };
    let mut worst: Scalar = 0.0;
    for (index, edge) in topology.edges().iter().enumerate() {
        let (Some(curve), Some(id)) = (edge.curve, topology.edge_id_at(index)) else {
            return Err(BooleanError::DanglingReference);
        };
        let curve = brep
            .curves3()
            .get(curve.index())
            .ok_or(BooleanError::DanglingReference)?;
        let span = brep
            .edge_interval(id)
            .ok_or(BooleanError::DanglingReference)?;
        worst = worst.max(edge_residue(
            curve,
            span,
            position(edge.start)?,
            position(edge.end)?,
        )?);
    }
    for face in topology.faces() {
        let Some(surface) = face.surface else {
            continue;
        };
        let surface = brep
            .surfaces()
            .get(surface.index())
            .ok_or(BooleanError::DanglingReference)?;
        for bound in &face.bounds {
            let wire = topology
                .loops()
                .get(bound.loop_id.index())
                .ok_or(BooleanError::DanglingReference)?;
            for (index, edge_use) in wire.edges.iter().enumerate() {
                let Some(pcurve) = edge_use.pcurve else {
                    continue;
                };
                let pcurve = brep
                    .curves2()
                    .get(pcurve.index())
                    .ok_or(BooleanError::DanglingReference)?;
                let span = brep
                    .pcurve_interval(bound.loop_id, index)
                    .ok_or(BooleanError::DanglingReference)?;
                let edge = topology
                    .edges()
                    .get(edge_use.edge.index())
                    .ok_or(BooleanError::DanglingReference)?;
                let (start, end) = match edge_use.orientation {
                    Orientation::Forward => (edge.start, edge.end),
                    Orientation::Reversed => (edge.end, edge.start),
                };
                worst = worst.max(use_residue(
                    pcurve,
                    span,
                    surface,
                    position(start)?,
                    position(end)?,
                )?);
            }
        }
    }
    Ok(worst)
}

/// Refuse a result whose vertices lie further off their curves than the
/// open session's report allows.
///
/// # Errors
///
/// [`BooleanError::VertexOffCurve`] over the bound, else as [`residue`].
pub(crate) fn check(brep: &ExactBRep, tolerance: Tolerance) -> Result<(), BooleanError> {
    let caller = report::caller().unwrap_or(tolerance);
    let allowed = bound(caller, report::rounding(), report::exact());
    if exceeds(residue(brep)?, allowed) {
        return Err(BooleanError::VertexOffCurve);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_construct::extrude::extrude_profile_exact;
    use axiolid_core::Vec3;
    use axiolid_curve::{Line2, Line3};
    use axiolid_profile::{Profile, RectangleProfile};
    use axiolid_surface::Plane;

    #[test]
    fn the_bound_follows_the_report_and_the_tolerance() {
        let metre = Tolerance::METRE;
        assert_eq!(bound(metre, 5.09e-6, true), 1e-6);
        assert_eq!(bound(metre, 1.8e-12, true), 1.8e-12);
        assert_eq!(bound(Tolerance::ZERO, 5.09e-6, true), 5.09e-6);
        assert_eq!(bound(metre, 5.09e-6, false), 1e-6);
    }

    #[test]
    fn a_residue_over_its_bound_or_not_a_number_exceeds() {
        assert!(exceeds(2e-6, 1e-6));
        assert!(!exceeds(1e-6, 1e-6));
        assert!(exceeds(Scalar::NAN, 1.0));
    }

    #[test]
    fn an_edge_end_off_its_curve_is_measured() {
        let x_axis = Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::X,
        });
        let span = Interval::new(0.0, 1.0);
        let residue = edge_residue(&x_axis, span, Point3::ZERO, Point3::new(1.0, 0.0, 2e-6))
            .expect("evaluates");
        assert!((residue - 2e-6).abs() <= 1e-15, "{residue}");
    }

    #[test]
    fn a_use_end_off_its_lifted_pcurve_is_measured() {
        let plane = Surface::Plane(Plane {
            frame: axiolid_core::Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
        });
        let pcurve = Curve2::Line(Line2 {
            origin: axiolid_core::Point2::new(0.0, 0.0),
            direction: axiolid_core::Vec2::new(1.0, 0.0),
        });
        let residue = use_residue(
            &pcurve,
            Interval::new(0.0, 1.0),
            &plane,
            Point3::ZERO,
            Point3::new(1.0, 0.0, 2e-6),
        )
        .expect("evaluates");
        assert!((residue - 2e-6).abs() <= 1e-15, "{residue}");
    }

    #[test]
    fn an_exact_block_has_no_residue() {
        let profile = Profile::Rectangle(RectangleProfile {
            x: 0.4,
            y: 0.5,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        });
        let block =
            extrude_profile_exact(&profile, Vec3::Z, 6.0, Tolerance::METRE).expect("a block");
        assert!(residue(&block).expect("measured") <= 1e-15);
        let _guard = report::open(&[&block], Tolerance::METRE);
        assert_eq!(check(&block, Tolerance::METRE), Ok(()));
    }
}
