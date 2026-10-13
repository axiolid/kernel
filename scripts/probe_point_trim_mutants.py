"""Mutation probe for B-spline trims by points in the graph (#306).

A graph trim's point selector on a B-spline is resolved by the inversion a
curve path's point bound takes (#304, `StationCurve::parameter_at_bound`),
so that lowering and derivation cannot disagree. Each mutant breaks one
step of that sharing: in `axiolid-mesh-compile` the B-spline arm itself
(dropped, so the closed-form inversion refuses; or replaced by the
evaluator's other B-spline point search, `locate2`/`locate3`), the 2D
point's lift into the plane, and the keeping of the shared inversion's
refusal; and, run against the compiler's test alone, the shared inversion
in `axiolid-evaluate` (its tolerance either way, its Newton refinement, its
refusal of several places and its merging of one place). Each must turn
`compile/tests/point_trimmed_spline.rs` red: the last group shows that a
change to the shared inversion reaches the graph trims, so there is one
inversion, not two.

Not probed, because no test can tell them apart: which of several refused
point selectors' refusals is kept (each test trim names one point per end).
"""
import probe_lib

D = "crates/execution/compile/src/directrix.rs"
S = "crates/execution/compile/src/station/basis.rs"
B = "crates/algorithms/parametric/evaluate/src/station/path_bound.rs"
TRIM = ["-p", "axiolid-mesh-compile", "--test", "point_trimmed_spline"]

ARM = (
    "        StationCurve::Two(Curve2::BSpline(_)) | StationCurve::Three(Curve3::BSpline(_)) => {\n"
    "            curve.parameter_at_bound(PathBound::Point(point))\n"
    "        }\n"
)

MUTANTS = [
    # The compiler's B-spline arm.
    ('a B-spline point is inverted in closed form', D, ARM, '', TRIM),
    ('a B-spline point is located by a second search', D, ARM,
     "        StationCurve::Two(spline @ Curve2::BSpline(_)) => axiolid_reference::curve::locate2(\n"
     "            spline,\n"
     "            Point2::new(point.x, point.y),\n"
     "            tolerance,\n"
     "        ),\n"
     "        StationCurve::Three(spline @ Curve3::BSpline(_)) => {\n"
     "            axiolid_reference::curve::locate3(spline, point, tolerance)\n"
     "        }\n", TRIM),
    ('a 2D point is lifted off its plane', S,
     'Some(axiolid_core::Point3::new(p.x, p.y, 0.0))',
     'Some(axiolid_core::Point3::new(p.x, p.y, 1e-6))', TRIM),
    ('a B-spline refusal is not kept', D,
     '                if spline && refusal.is_none() {',
     '                if false && spline && refusal.is_none() {', TRIM),
    ('a 3D trim drops the kept refusal', D,
     '            Some(refusal) => refusal,',
     '            Some(_) => GeomError::InvalidInput(String::new()),', TRIM),
    ('a 2D trim drops the kept refusal', S,
     '            Some(refusal) => refusal,',
     '            Some(_) => GeomError::InvalidInput(String::new()),', TRIM),
    # The shared inversion, seen through the graph trims alone.
    ('the shared tolerance is a thousand times looser', B,
     '    JOINT_TOLERANCE * point.abs().max_element().max(1.0)',
     '    1e3 * JOINT_TOLERANCE * point.abs().max_element().max(1.0)', TRIM),
    ('the shared tolerance is a thousand times tighter', B,
     '    JOINT_TOLERANCE * point.abs().max_element().max(1.0)',
     '    1e-3 * JOINT_TOLERANCE * point.abs().max_element().max(1.0)', TRIM),
    ('the shared B-spline candidate is not refined', B,
     '            let refined = self.refine(at[i], a, b, point)?;',
     '            let refined = at[i];', TRIM),
    ('the shared refusal of several places resolves to the first', B,
     '        _ => Err(refused(PATH_POINT_AMBIGUOUS)),',
     '        [t, ..] => Ok(*t),', TRIM),
    ('the shared B-spline places always merge', B,
     '            Ok(between && self.distance(0.5 * (ta + tb), point)? <= tolerance)',
     '            Ok(true || between)', TRIM),
]

probe_lib.main(MUTANTS, timeout=1800)
