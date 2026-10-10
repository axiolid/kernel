"""Mutation probe for exact decisions in exact_surface_intersection (#119)."""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
PATH = ROOT / "crates/algorithms/parametric/nurbs/src/exact_surface_intersection.rs"

MUTANTS = [('sphere tangency read as a circle', '        Sign::Positive => {}\n        // Exactly tangent: a touch is a point, not a curve.\n        Sign::Zero => return Err(ExactIntersectionRefusal::NotRegularCurve),', '        Sign::Positive | Sign::Zero => {}\n        // Exactly tangent: a touch is a point, not a curve.\n        Sign::Negative => return Err(ExactIntersectionRefusal::NotRegularCurve),'), ('parallel cylinder test weakened', '    if esign(&along) == Sign::Zero {\n        // Plane parallel to the axis', '    if along.to_f64().abs() < 1e-12 {\n        // Plane parallel to the axis'), ('perpendicular cylinder test dropped', '    let perpendicular = ecross_is_zero(&exact_axis, &exact_normal);\n    // cos(theta)', '    let perpendicular = false;\n    // cos(theta)'), ('cone perpendicular test dropped', '        (Some(a), Some(n)) => ecross_is_zero(&a, &n),', '        (Some(a), Some(n)) => { let _ = (a, n); false }'), ('tangent ruling counted as two', '    let tangent_plane = match esign(&numerator) {\n        Sign::Positive => false,\n        Sign::Zero => true,', '    let tangent_plane = match esign(&numerator) {\n        Sign::Positive => false,\n        Sign::Zero => false,')]

probe_lib.main(MUTANTS, tests=[["-p", "axiolid-nurbs", "--test", "exact_intersection"]], path=PATH, timeout=900)
