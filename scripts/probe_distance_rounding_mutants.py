"""Mutation probe for the derived rounding margins of the certified
distance searches (#318).

The margins are ulps of the coordinates, not `1e-12` of them: a pad per
side in `exact_distance.rs`, the floor a search settles at, the touch
margin, and the Hausdorff matched bound, re-charted translate, cone pole
and displacement gate (`exact_hausdorff.rs`, `exact_hausdorff/translate.rs`).
Each mutant drops or weakens one, or certifies a foot witness that is not on
its face, and must turn a test red. The unit tests of `exact_distance.rs`
pin the constants and the far projection; the construct tests measure
bodies at (450000, 5600000), where the old margins left a width of `1e-5`.

Not probed, because their loss changes only how fast a search converges,
which no test can observe against a step budget: the witnesses counted when
their pair or node is queued, the first-in-first-out ties of `Nearest`, the
point queries' floor, and the Hausdorff outer loop's floor. Their tests
(`parallel_faces_close_on_a_foot_witness`, `far_boxes_apart_settle_at_the_floor`,
`a_point_under_an_edge_finds_the_edge_witness`) close either way. Nor the
displacement gate of `displacements` losing its rounding (it survived): a
displacement only seeds a support point to query, whose lower bound holds
for any point of the source, so a seed taken or missed changes how fast the
interval closes, and `square_translates_far_from_the_origin_close` closes
either way.
"""
import probe_lib

D = "crates/algorithms/query/measure/src/exact_distance.rs"
T = "crates/algorithms/query/measure/src/exact_distance/touch.rs"
W = "crates/algorithms/query/measure/src/exact_distance/witness.rs"
H = "crates/algorithms/query/measure/src/exact_hausdorff.rs"
X = "crates/algorithms/query/measure/src/exact_hausdorff/translate.rs"
UNIT = ["-p", "axiolid-measure", "--features", "exact", "--lib"]
HAUS = ["-p", "axiolid-construct", "--test", "boundary_hausdorff"]

MUTANTS = [
    ("the pad drops its rounding of the centre", D,
     "        + POINT_EPSILONS * Scalar::EPSILON * scale\n        + Scalar::MIN_POSITIVE",
     "        + 0.0 * scale\n        + Scalar::MIN_POSITIVE", UNIT),
    ("a projected range drops its pad", D,
     "        let pad = PROJECTION_EPSILONS * Scalar::EPSILON * self.scale;",
     "        let pad = 0.0 * self.scale;", UNIT),
    ("an interval closes only at the accuracy", D,
     "    upper - lower <= accuracy.max(floor)",
     "    upper - lower <= accuracy", UNIT),
    ("the floor drops its scale", D,
     "    FLOOR_EPSILONS * Scalar::EPSILON * scale + RATIO_EPSILONS * Scalar::EPSILON * upper",
     "    0.0 * scale + RATIO_EPSILONS * Scalar::EPSILON * upper", UNIT),
    ("touch margin back to 1e-12 of the magnitudes", T,
     "    d <= FLOOR_EPSILONS * Scalar::EPSILON * (p.length() + q.length())",
     "    d <= 1e-12 * (p.length() + q.length())", UNIT),
    ("a foot beyond the face certified", W,
     "    let certified = (inside && in_patch)\n"
     "        || match &side.domains[face] {\n"
     "            Some(domain) => domain.contains(at)? == Some(true),\n"
     "            None => false,\n"
     "        };",
     "    let certified = true;", UNIT),
    ("matched bound back to 1e-12 of the magnitudes", H,
     "    let bound = (centre.length() + spread) * (1.0 + RATIO_EPSILONS * Scalar::EPSILON)\n"
     "        + POINT_EPSILONS * Scalar::EPSILON * scale\n"
     "        + Scalar::MIN_POSITIVE;",
     "    let bound = centre.length() + spread + 1e-12 * scale + Scalar::MIN_POSITIVE;", HAUS),
    ("translate deviation back to 1e-12", X,
     "    let gap = gap * (1.0 + RATIO_EPSILONS * Scalar::EPSILON)\n"
     "        + POINT_EPSILONS * Scalar::EPSILON * (size + s.length());",
     "    let gap = gap + 1e-12 * (size + s.length());", HAUS),
    ("a plane's Lipschitz bound grown by its chart", X,
     "        Surface::Plane(_) | Surface::Cylinder(_) | Surface::EllipticalCylinder(_) => sum,",
     "        Surface::Plane(_) | Surface::Cylinder(_) | Surface::EllipticalCylinder(_) => {\n"
     "            sum * (1.0 + lo.y.abs().max(hi.y.abs()) + grow)\n"
     "        }", (UNIT, HAUS)),
]

probe_lib.main(MUTANTS, timeout=1200)
