//! How far a node placed at a station lies from its exact placement (#405,
//! ADR 0082 amendment).
//!
//! [`placement`](super::placement) puts a source's local `x`, `y`, `z` on
//! the oriented tangent, lateral and up of the frame the kernel reads at
//! the station, its origin at the offsets read in the base frame. That
//! frame is exact on a line ([`ResolvedPlacement::exact`]); on any other
//! basis it is the exact frame law's at a distance read by an estimate, and
//! [`placement_accuracy`] states how far the placement lies from the exact
//! one: [`PlacementAccuracy::Exact`], [`PlacementAccuracy::Certified`] with
//! a [`PlacementBound`], or [`PlacementAccuracy::Unbounded`] by name.
//!
//! # The exact placement
//!
//! The source's local `x`, `y`, `z` on the oriented axes of the exact
//! frame law (`axiolid_reference::station::FrameLaw`, #316, #361) at the
//! station's distance, snapped by its window (#294) and read from its seam
//! side (#263), its origin at the offsets read in the exact base frame, a
//! section's or the plan's. The kernel's frame is compared with it, never
//! trusted: `FrameLaw::error_at` gives the error `e_p` of the read frame's
//! point and `e_a` of each of its axes.
//!
//! # The bound
//!
//! With `L1 = |lateral| + |vertical| + |longitudinal|` of the offsets, the
//! station's point `P = O + lateral L + vertical U + longitudinal T` is
//! within `e_p + e_a L1` of the exact one, plus the rounding of computing
//! it; the placement's axes are the base frame's, or its axes turned by
//! the orientation, each a combination of the three with unit components
//! `w`, `|w|_1 <= sqrt 3`: within `sqrt 3 e_a` of the exact ones, plus the
//! rounding of turning them. A local point `q` then lands within
//! `point_error + axis_error |q|_1` of its exact place
//! ([`PlacementBound::displacement`]); the bound is the same for every
//! point of a mesh, since `|q|_1` is convex. The rounding of applying the
//! transform in `f64` is not part of it: the mesh compiler adds it as it
//! adds it to every instance (#295).
//!
//! # Unbounded
//!
//! By name ([`UNCERTIFIED_STATION_FRAME`]), where `FrameLaw::error_at`
//! certifies nothing: frames read by quadrature (an ellipse, a B-spline, an
//! offset), chains, banked curves over intrinsic or chain plans, intrinsic
//! space curves with torsion or tilted, elevated curves with an intrinsic
//! profile or another plan, a skewed circle, a steep tangent, a basis
//! placed at a station whose own frame is not exact, and a station within
//! its distance's resolution of a seam it is not read on.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar};
use axiolid_model::{GeometryGraph, GeometryNode, NodeId, StationOffsets};
use axiolid_reference::station::POINT_ROUNDING;

use super::basis::Basis;
use super::bound::reading as frame_reading;
use super::resolve_parts;

/// Why a placement whose frame has no certified error is unbounded.
pub(crate) const UNCERTIFIED_STATION_FRAME: &str =
    "a node placed at a station whose frame has no certified error (a frame read by quadrature \
     on an ellipse, a B-spline or an offset, a chain, a banked curve over a chain \
     plan, an intrinsic space curve with torsion or tilted, an elevated curve with an intrinsic \
     profile or another plan, a skewed circle, a steep tangent, a basis placed at a station \
     whose frame is not exact, or a station within its distance's resolution of a seam it is \
     not read on)";

/// The rounding of turning the base frame by an orientation, per axis, for
/// an axis and reference direction at a right angle: turning is conditioned
/// by `1 / |axis x reference|` (Gram-Schmidt of two nearly parallel
/// vectors).
const TURN_ROUNDING: Scalar = 1e-14;

/// The margin on every composed term, over the rounding of composing it.
const MARGIN: Scalar = 1e-12;

/// How far a placement at a station lies from the exact one, see
/// [placing at a station](super#placing-at-a-station-264).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlacementBound {
    /// How far the placement's origin, the station's point, lies from the
    /// exact one.
    pub point_error: Scalar,
    /// How far each of the placement's axes lies from the exact one.
    pub axis_error: Scalar,
}

impl PlacementBound {
    /// How far the image of the source's local point `local` may lie from
    /// its exact image: `point_error + axis_error (|x| + |y| + |z|)`.
    #[must_use]
    pub fn displacement(&self, local: Point3) -> Scalar {
        self.point_error + self.axis_error * (local.x.abs() + local.y.abs() + local.z.abs())
    }
}

/// How far a node placed at a station lies from its exact placement.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlacementAccuracy {
    /// The frame is exact, rounding aside (`ResolvedPlacement::exact`): a
    /// line basis, itself placed, if at all, in exact frames.
    Exact,
    /// The placement lies within the bound of the exact one.
    Certified(PlacementBound),
    /// No bound is certified, for the reason named.
    Unbounded(&'static str),
}

/// How far the placement [`placement`](super::placement) resolves for the
/// [`InstanceAtStation`] node `id` lies from the exact one, its frame law's
/// quadratures certified to `tolerance` (see
/// [placing at a station](super#placing-at-a-station-264)).
///
/// [`PlacementAccuracy::Exact`] exactly when the placement is
/// ([`ResolvedPlacement::exact`](super::ResolvedPlacement::exact)); a
/// certified bound is finite and not negative. Both read the frame the
/// placement reads, bit for bit.
///
/// # Errors
///
/// A node that is not an [`InstanceAtStation`], a `tolerance` that is not
/// finite and positive, and every refusal of [`placement`](super::placement),
/// each by name.
///
/// [`InstanceAtStation`]: axiolid_model::InstanceAtStation
pub fn placement_accuracy(
    graph: &GeometryGraph,
    id: NodeId,
    tolerance: Scalar,
) -> GeomResult<PlacementAccuracy> {
    let Some(GeometryNode::InstanceAtStation(placed)) = graph.get(id) else {
        return Err(GeomError::InvalidInput(format!(
            "node {id:?} is not an instance at a station"
        )));
    };
    if !(tolerance.is_finite() && tolerance > 0.0) {
        return Err(GeomError::InvalidInput(format!(
            "station: a certified placement needs a finite, positive tolerance, not {tolerance}"
        )));
    }
    let station = &placed.station;
    let parts = resolve_parts(graph, station, 0)?;
    if parts.exact {
        return Ok(PlacementAccuracy::Exact);
    }
    let unbounded = PlacementAccuracy::Unbounded(UNCERTIFIED_STATION_FRAME);
    // A curve placed at a station whose frame is not exact has no exact
    // frame law to compare with.
    if matches!(
        parts.basis,
        Basis::Atomic {
            placed_exactly: false,
            ..
        }
    ) {
        return Ok(unbounded);
    }
    let Some(reading) = frame_reading(station.station.frame) else {
        return Ok(unbounded);
    };
    let law = parts.basis.frame_law(tolerance)?;
    let Some(error) = law.error_at(parts.distance, station.seam, reading) else {
        return Ok(unbounded);
    };
    // `None` for the base frame, else the sine of the angle between the
    // orientation's axis and reference direction.
    let turn = if station.orientation.is_base() {
        None
    } else {
        match station.orientation.unit_axes() {
            Ok((axis, reference)) => Some(axis.cross(reference).length()),
            Err(_) => return Ok(unbounded),
        }
    };
    if turn.is_some_and(|sin| !(sin.is_finite() && sin > 0.0)) {
        return Ok(unbounded);
    }
    let bound = compose(
        (error.point, error.axis),
        &station.station.station.offsets,
        turn,
        parts.resolved.point,
    );
    Ok(
        if bound.point_error.is_finite() && bound.axis_error.is_finite() {
            PlacementAccuracy::Certified(bound)
        } else {
            unbounded
        },
    )
}

/// The placement's bound from the errors `(read_point, read_axis)` of the
/// frame read: `turn` the sine of the angle between the orientation's axis
/// and reference direction, `None` when it leaves the base frame,
/// `offsets` read in it and `point` the station's point.
fn compose(
    (read_point, read_axis): (Scalar, Scalar),
    offsets: &StationOffsets,
    turn: Option<Scalar>,
    point: Point3,
) -> PlacementBound {
    let reach = offsets.lateral.abs() + offsets.vertical.abs() + offsets.longitudinal.abs();
    let scale = (point.abs().max_element() + reach).max(1.0);
    let point_error = (read_point + read_axis * reach + POINT_ROUNDING * scale) * (1.0 + MARGIN);
    let axis_error = match turn {
        None => read_axis,
        // Each turned axis combines the base axes with unit components.
        Some(sin) => 3.0_f64.sqrt() * read_axis + TURN_ROUNDING / sin,
    } * (1.0 + MARGIN);
    PlacementBound {
        point_error,
        axis_error,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::Vec3;

    /// A frame read: its point and its tangent, lateral and up.
    #[derive(Clone, Copy)]
    struct Read {
        point: Vec3,
        axes: [Vec3; 3],
    }

    fn base() -> Read {
        // Orthonormal, not axis-aligned.
        let t = Vec3::new(0.6, 0.8, 0.0);
        let l = Vec3::new(-0.8, 0.6, 0.0);
        Read {
            point: Vec3::new(1000.0, 2000.0, 30.0),
            axes: [t, l, Vec3::Z],
        }
    }

    /// The read frame off by `a` in its point and `b` in each axis, all in
    /// the direction `n`: the worst case the bound allows.
    fn perturbed(read: Read, a: Scalar, b: Scalar, n: Vec3) -> Read {
        Read {
            point: read.point + a * n,
            axes: read.axes.map(|axis| axis + b * n),
        }
    }

    fn offsets() -> StationOffsets {
        // Positive, so the axes' errors add.
        StationOffsets::new(2.0, 1.0, 0.5)
    }

    /// The station's point: the origin moved by the offsets in the frame.
    fn station(read: Read, offsets: &StationOffsets) -> Vec3 {
        read.point
            + offsets.longitudinal * read.axes[0]
            + offsets.lateral * read.axes[1]
            + offsets.vertical * read.axes[2]
    }

    /// A local point placed in the axes `axes`, at `origin`.
    fn placed(origin: Vec3, axes: [Vec3; 3], local: Point3) -> Vec3 {
        origin + local.x * axes[0] + local.y * axes[1] + local.z * axes[2]
    }

    /// Axes turned by the unit components `w` of each of the base axes.
    fn turned(axes: [Vec3; 3], w: Vec3) -> [Vec3; 3] {
        let mix = |u: Vec3| u.x * axes[0] + u.y * axes[1] + u.z * axes[2];
        [
            mix(w),
            mix(Vec3::new(w.y, w.z, w.x)),
            mix(Vec3::new(w.z, w.x, w.y)),
        ]
    }

    fn at() -> Point3 {
        let p = station(base(), &offsets());
        Point3::new(p.x, p.y, p.z)
    }

    const A: Scalar = 3e-7;
    const B: Scalar = 2e-9;

    fn n() -> Vec3 {
        Vec3::new(1.0, -2.0, 2.0) / 3.0
    }

    #[test]
    fn the_point_error_is_the_station_points_worst_displacement() {
        let error = (A, B);
        let offsets = offsets();
        let exact = base();
        let read = perturbed(exact, A, B, n());
        let moved = (station(read, &offsets) - station(exact, &offsets)).length();
        let reach = 3.5; // 2 + 1 + 0.5
                         // The worst case reaches `a + b L1` exactly.
        assert!((moved - (A + B * reach)).abs() <= 1e-12, "{moved}");
        let bound = compose(error, &offsets, None, at());
        assert!(bound.point_error >= moved, "{bound:?} < {moved}");
        assert!(
            bound.point_error <= moved * (1.0 + 1e-9) + 1e-9,
            "{bound:?}"
        );
        // No offsets, no lever on the axes.
        let flat = compose(error, &StationOffsets::ZERO, None, at());
        assert!(flat.point_error < A + 1e-9 && flat.point_error >= A);
    }

    #[test]
    fn a_turned_axis_moves_by_sqrt_three_of_the_axis_error() {
        let w = Vec3::ONE / 3.0_f64.sqrt();
        let exact = base();
        let read = perturbed(exact, A, B, n());
        let (e, r) = (turned(exact.axes, w), turned(read.axes, w));
        let moved = (0..3)
            .map(|i| (e[i] - r[i]).length())
            .fold(0.0, Scalar::max);
        assert!((moved - 3.0_f64.sqrt() * B).abs() <= 1e-15, "{moved}");
        let error = (A, B);
        let off = offsets();
        let turned_bound = compose(error, &off, Some(1.0), at());
        assert!(
            turned_bound.axis_error >= moved,
            "{turned_bound:?} < {moved}"
        );
        assert!(turned_bound.axis_error <= moved * (1.0 + 1e-6) + 1e-13);
        // The base frame's axes move by `b` and the bound is not looser.
        let base_bound = compose(error, &off, None, at());
        assert!(base_bound.axis_error >= B && base_bound.axis_error <= B * (1.0 + 1e-9));
        assert!(turned_bound.axis_error > base_bound.axis_error);
    }

    #[test]
    fn a_body_vertex_moves_by_the_one_norm_of_its_local_point() {
        let off = offsets();
        let exact = base();
        let read = perturbed(exact, A, B, n());
        let error = (A, B);
        let w = Vec3::ONE / 3.0_f64.sqrt();
        for k in [0.0, 1.0, -3.0, 6.0] {
            let local = Point3::new(k, k, k);
            // The base frame.
            let bound = compose(error, &off, None, at());
            let moved = (placed(station(read, &off), read.axes, local)
                - placed(station(exact, &off), exact.axes, local))
            .length();
            let shown = bound.displacement(local);
            assert!(shown >= moved, "k {k}: {shown} < {moved}");
            // `|q|_1` times the axes' error, and no looser than the margin
            // where the lever's signs agree with the offsets'.
            assert!(
                k < 0.0 || shown <= moved * (1.0 + 1e-9) + 1e-9,
                "k {k}: {shown} {moved}"
            );
            // A turned frame, components `(1, 1, 1) / sqrt 3`.
            let bound = compose(error, &off, Some(1.0), at());
            let moved = (placed(station(read, &off), turned(read.axes, w), local)
                - placed(station(exact, &off), turned(exact.axes, w), local))
            .length();
            let shown = bound.displacement(local);
            assert!(shown >= moved, "turned, k {k}: {shown} < {moved}");
            assert!(
                k < 0.0 || shown <= moved * (1.0 + 1e-6) + 1e-9,
                "turned, k {k}: {shown} {moved}"
            );
        }
    }
}
