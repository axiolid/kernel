//! The ADR 0014 near-degenerate halfspace case, as an executable regression.
//!
//! `test/fixtures/ifclite-geometry/issue_1155_halfspace_flyaway.ifc` clips a
//! millimetre-scale column with an `IfcHalfSpaceSolid` whose plane normal is
//! `(-1, 1.99999991118124e-9, 0)` — 2e-9 off axis. The dimensions below are
//! transcribed from ADR 0014, which records the expected result bounds:
//!
//! ```text
//! column bounds  min(-125, -125, 11940)  max(125, 125, 23880)
//! clipped bounds min( 124.9999995, -125, 11940)  max(125, 124.9999, 23880)
//! ```
//!
//! The kept region is therefore a SLIVER roughly 5e-7 mm wide in x, spanning
//! the full 250 mm of y — an aspect ratio near 5e9:1. That is the whole point
//! of the fixture: the near-degenerate plane is not incidental, it is what
//! produces a body whose width is at the edge of what f64 coordinates at
//! millimetre magnitude can express.
//!
//! The historical failure mode is "flyaway": the clip emits geometry far
//! outside the input bounds. That is the dangerous shape of this bug, because a
//! flyaway result is still structurally valid and still passes a manifold
//! check — it is only detectable by measuring where the geometry ended up.
//! Volume conservation does not catch it either, since escaped material can be
//! self-consistent.
//!
//! The mesh-level geometry is asserted here rather than through the IFC lowerer
//! so this gate stays meaningful before the parser path is wired up, matching
//! the `issue_2019` fixture next door.

mod support;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Point3, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
use support::{boxx, volume};

/// Column half-extent in millimetres, from the ADR bounds.
const HALF: f64 = 125.0;
/// Column base height in millimetres.
const Z0: f64 = 11_940.0;
/// Column top height in millimetres.
const Z1: f64 = 23_880.0;

/// The off-axis component of the clipping plane normal, verbatim from ADR 0014.
///
/// This exact value matters: it is small enough that a naive implementation
/// normalises it away or divides by it, and large enough that it is not zero.
const OFF_AXIS: f64 = 1.999_999_911_181_24e-9;

/// Width in millimetres of the kept sliver at its widest point.
///
/// From the ADR bounds: `max.x 125` minus `min.x 124.9999995`.
const SLIVER_WIDTH: f64 = 5e-7;

/// The column being clipped: 250 x 250 x 11940 mm, centred on the z axis.
fn column() -> TriMesh {
    boxx(0.0, 0.0, Z0, 2.0 * HALF, 2.0 * HALF, Z1 - Z0, 0.0)
}

/// Axis-aligned bounds of a mesh as `(min, max)`.
fn bounds(mesh: &TriMesh) -> (Point3, Point3) {
    let mut min = Point3::new(f64::MAX, f64::MAX, f64::MAX);
    let mut max = Point3::new(f64::MIN, f64::MIN, f64::MIN);
    for point in &mesh.positions {
        min = Point3::new(min.x.min(point.x), min.y.min(point.y), min.z.min(point.z));
        max = Point3::new(max.x.max(point.x), max.y.max(point.y), max.z.max(point.z));
    }
    (min, max)
}

/// A bounded stand-in for the half-space whose boundary plane is `x = c - t*y`.
///
/// A true `IfcHalfSpaceSolid` is unbounded; mesh booleans need it bounded.
/// Graph compilation bounds one in `axiolid_construct::half_space`; this
/// provider-level fixture has no graph, so it builds its own bound.
///
/// Sizing matters and is easy to get wrong. The bound must comfortably cover
/// the column so a flyaway cannot be masked by a tight cutter, but it must NOT
/// be enormous: the plane tilt is 2e-9, so over a cutter of extent `L` the
/// boundary sweeps `2e-9 * L` in x. At `L = 2e5` mm that sweep is 4e-4 mm,
/// which swamps the 5e-7 mm sliver this fixture exists to produce and simply
/// removes the whole column. The cutter is therefore sized to a few column
/// extents, where the tilt still resolves the sliver.
///
/// The kept sliver tapers from 5e-7 mm wide at `y = +125` to zero at
/// `y = -125`, reproducing the ADR's `min.x 124.9999995 / max.x 125`.
fn near_degenerate_cutter() -> TriMesh {
    let angle = (-OFF_AXIS).atan2(1.0);
    // Depth in x and height in z are generous; the y extent is what interacts
    // with the tilt, so it is kept at four column widths.
    let depth = 8.0 * HALF;
    let height = 4.0 * (Z1 - Z0);
    let width = 4.0 * (2.0 * HALF);
    // Plane offset placing the widest point of the kept sliver at x = HALF.
    let plane_offset = HALF + OFF_AXIS * HALF - SLIVER_WIDTH;
    // The rotated box's `+x` face passes through `plane_offset` on the axis.
    let cx = plane_offset - (depth / 2.0) / angle.cos();
    boxx(cx, 0.0, Z0 - height / 4.0, depth, width, height, angle)
}

/// A near-degenerate halfspace clip must stay inside the column, or refuse.
///
/// The provider may decline this input with a typed `Unsupported` error —
/// refusing is an honest answer for a case it cannot certify. What it must
/// never do is return a silently wrong mesh, so both accepted outcomes are
/// checked and any third outcome fails the test.
#[test]
fn near_degenerate_halfspace_clip_does_not_fly_away() {
    let column = column();
    let cutter = near_degenerate_cutter();
    let (column_min, column_max) = bounds(&column);

    let outcome = BoolmeshBoolean::new().boolean(
        &column,
        &cutter,
        BooleanOperator::Difference,
        &ExecutionOptions::new(Tolerance::METRE),
    );

    let result = match outcome {
        Ok(outcome) => outcome.mesh,
        Err(GeomError::Unsupported { .. }) => return,
        Err(error) => panic!("expected a result or a typed Unsupported refusal, got {error:?}"),
    };

    // An empty result would mean the clip removed the entire column, which for
    // this fixture is a FAILURE rather than a benign outcome: the sliver is the
    // subject under test. An earlier revision of this fixture returned empty
    // (the cutter was sized so large that the 2e-9 tilt swept past the sliver),
    // and an `is_empty` early return turned that into a silent vacuous pass.
    assert!(
        !result.positions.is_empty(),
        "the clip removed the whole column; the sliver under test never existed"
    );

    assert!(
        result.validate_structure().is_ok(),
        "result must be well formed"
    );

    // The load-bearing assertion. A flyaway escapes the input bounds, and no
    // structural or manifold check detects that.
    let (min, max) = bounds(&result);
    let slack = 1e-6;
    assert!(
        min.x >= column_min.x - slack
            && min.y >= column_min.y - slack
            && min.z >= column_min.z - slack
            && max.x <= column_max.x + slack
            && max.y <= column_max.y + slack
            && max.z <= column_max.z + slack,
        "clipped geometry escaped the column: got min({}, {}, {}) max({}, {}, {}), \
         column min({}, {}, {}) max({}, {}, {})",
        min.x,
        min.y,
        min.z,
        max.x,
        max.y,
        max.z,
        column_min.x,
        column_min.y,
        column_min.z,
        column_max.x,
        column_max.y,
        column_max.z
    );

    // A clip removes material or leaves it unchanged; it never adds any.
    let clipped_volume = volume(&result);
    let column_volume = volume(&column);
    assert!(
        clipped_volume <= column_volume * (1.0 + 1e-9),
        "clip grew the volume: {clipped_volume} > {column_volume}"
    );

    // The kept body is the near-degenerate sliver, not a fraction of the
    // column. Asserting thinness is what makes this the issue_1155 case rather
    // than an ordinary clip that happens to stay in bounds: a provider that
    // quietly widened the sliver to something numerically comfortable would
    // pass every check above and fail here.
    assert!(
        max.x - min.x <= 1e-3,
        "expected a sliver ~5e-7 mm wide in x, got {}",
        max.x - min.x
    );
}

/// The same clip with an exactly axis-aligned plane must also stay bounded.
///
/// This is the control. It isolates the 2e-9 perturbation: if this passed while
/// the perturbed case failed, the off-axis term would be proven to be the
/// trigger rather than an incidental detail of the fixture.
#[test]
fn the_axis_aligned_control_clip_also_stays_bounded() {
    let column = column();
    let span = 16.0 * (Z1 - Z0);
    let cutter = boxx(
        HALF - span / 2.0,
        0.0,
        Z0 - span,
        span,
        span,
        3.0 * span,
        0.0,
    );
    let (column_min, column_max) = bounds(&column);

    let result = BoolmeshBoolean::new()
        .boolean(
            &column,
            &cutter,
            BooleanOperator::Difference,
            &ExecutionOptions::new(Tolerance::METRE),
        )
        .expect("an axis-aligned clip is a supported case")
        .mesh;

    if result.positions.is_empty() {
        return;
    }

    let (min, max) = bounds(&result);
    let slack = 1e-6;
    assert!(
        min.x >= column_min.x - slack
            && max.x <= column_max.x + slack
            && min.y >= column_min.y - slack
            && max.y <= column_max.y + slack
            && min.z >= column_min.z - slack
            && max.z <= column_max.z + slack,
        "axis-aligned control escaped the column"
    );
}
