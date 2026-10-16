"""Mutation probe for bounded half-spaces along their placement's z (#321,
ADR 0085 format 1.4) and for their whole frame (#322).

Each mutant breaks one decision a `BoundedHalfSpaceAlongZ` rests on: in
`axiolid-model` the wire format (the variant read in a 1.3 payload, a
payload never labelled 1.4, the feature test missing the variant) and the
variant's validation; in `axiolid-mesh-compile` the variant read along the
plane normal; in `axiolid-construct` the prism (the placement's z replaced by
the normal, the parallel-axis and degenerate-placement refusals, the
boundary carried along the normal instead of the axis, its origin left
off the plane, the far cap at the reach along the axis rather than across
the plane, the stations ordered by the side rather than the frame, a
clockwise ring left as wound); in `axiolid-mesh-compile` the exact tool
(the axis not dispatched, the plane's side not cut from the prism, the
axis not turned across the boundary's plane, the prism's margin meeting
the side's). Each must turn a test red.

Equivalent mutants, deliberately not listed: the prism's caps at the
envelope's own heights rather than a margin past them (every clip the
tests cut lies above the plane, which stands inside the envelope, so a
cap there removes nothing more); the tool's boolean report dropped from
the clip's (the tests' tools are cut exactly).

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import probe_lib

TIMEOUT = 1800
V = "crates/representations/modeling/graph/src/validation.rs"
B = "crates/execution/compile/src/half_space_boundary.rs"
G = "crates/representations/modeling/graph/src/wire.rs"
H = "crates/algorithms/construction/construct/src/half_space.rs"
C = "crates/execution/compile/src/exact/clip.rs"
P = "crates/algorithms/construction/construct/src/profile_lower.rs"
M = "crates/execution/compile/src/compiler.rs"

WIRE = ["-p", "axiolid-model", "--all-features", "--test", "wire"]
SLAB = ["-p", "axiolid-construct", "--test", "half_space"]
AXIS = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "bounded_half_space_prism_axis"]
WHOLE = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "bounded_half_space_whole_frame"]
DERIVED = ["-p", "axiolid-construct", "--test", "derived_composite"]

MUTANTS = [
    # The wire format (ADR 0085, format 1.4).
    ("the variant is read in a 1.3 payload", G,
     "    if version < ALONG_Z_SINCE && is_along_z(node) {",
     "    if false && is_along_z(node) {", [WIRE]),
    ("a payload is never labelled 1.4", G,
     "        version = version.max(ALONG_Z_SINCE);",
     "        version = version.max(BASIS_END_SINCE);", [WIRE]),
    ("the feature test misses the variant", G,
     "        GeometryNode::SolidOperation(SolidOperation::BoundedHalfSpaceAlongZ { .. })",
     "        GeometryNode::SolidOperation(SolidOperation::BoundedHalfSpace { .. })", [WIRE]),
    ("the variant's boundary is not validated", V,
     "            expect_reference(nodes, *boundary, ExpectedReference::HalfSpaceBoundary)",
     "            let _ = boundary;\n            Ok(())", [WIRE]),
    # Both compilers read the variant along the placement's z.
    ("the variant is read along the plane normal", B,
     "        } => Some((*half_space, *boundary, *placement, PrismAxis::PlacementZ)),",
     "        } => Some((*half_space, *boundary, *placement, PrismAxis::PlaneNormal)),",
     [AXIS]),
    # The prism.
    ("the placement's z is replaced by the normal", H,
     "            axis: placement.matrix3.z_axis,",
     "            axis: x.cross(y),", [SLAB, AXIS]),
    ("an axis in the plane is let through", H,
     "        if crossing.abs() <= unit_floor(tolerance) {",
     "        if crossing == 0.0 {", [SLAB, AXIS]),
    ("a degenerate placement is let through", H,
     "    if volume.abs() <= unit_floor(tolerance) * x.length() * y.length() * z.length() {",
     "    if false {", [SLAB, AXIS]),
    ("the boundary is carried along the normal", H,
     "    let carry = |v: Vec3| v - axis * (v.dot(normal) / crossing);",
     "    let carry = |v: Vec3| v - normal * v.dot(normal);", [SLAB, AXIS]),
    ("the boundary's origin is left off the plane", H,
     "    let origin = prism.origin - axis * ((prism.origin - plane.origin).dot(normal) / crossing);",
     "    let origin = prism.origin;", [SLAB, AXIS]),
    ("the far cap is at the reach along the axis", H,
     "    let step = toward * (reach / lean.abs());",
     "    let step = toward * reach;", [SLAB, AXIS]),
    ("the stations are ordered by the side alone", H,
     "    let ordered = if face.dot(toward) > 0.0 {",
     "    let ordered = if agreement {", [SLAB, AXIS]),
    ("a clockwise ring is left as wound", H,
     "    if twice_area >= 0.0 {",
     "    if true {", [SLAB]),
    # The exact tool.
    ("the exact compiler does not dispatch the axis", C,
     [("            Some((footprint, frame, PrismAxis::PlacementZ)) => {",
       "            Some((footprint, frame, PrismAxis::PlacementZ)) if false => {"),
      ("            Some((footprint, frame, PrismAxis::PlaneNormal)) => Some((footprint, frame)),",
       "            Some((footprint, frame, _)) => Some((footprint, frame)),")],
     [AXIS]),
    ("the plane's side is not cut from the prism", C,
     "            Ok((side, _)) => side,",
     "            Ok(_) => {\n"
     "                return long\n"
     "                    .transformed(placement)\n"
     "                    .map(|solid| Ok((solid, None)))\n"
     "                    .map_err(transform_error)\n"
     "            }", [AXIS]),
    ("the axis is not turned across the boundary's plane", C,
     "        if along.z < 0.0 {\n            along = -along;\n        }",
     "        if false {\n            along = -along;\n        }", [AXIS]),
    ("the prism's caps meet the side's", C,
     "        let margin = 0.5 * (local.hi - local.lo).length() + 8.0 * tolerance.linear();",
     "        let margin = 0.25 * (local.hi - local.lo).length() + 4.0 * tolerance.linear();",
     [AXIS]),
    # The whole frame (#322).
    ("a left-handed placement is read right-handed", H,
     "            handedness < 0.0,",
     "            false,",
     [SLAB, WHOLE]),
    ("a frame mirrored about the normal is read unmirrored", H,
     "            turn < 0.0,",
     "            false,", [SLAB]),
    ("a frame upright on the plane is let through", H,
     "        if turn.abs() <= unit_floor(tolerance) {",
     "        if false {", [SLAB]),
    ("the frame's scale is dropped", H,
     "            x: base.x * length,",
     "            x: base.x,", [SLAB, WHOLE]),
    ("the frame's shear is dropped", H,
     "            y: base.x * shear + across * rise,",
     "            y: across * rise,", [SLAB, WHOLE]),
    ("the frame's y is read at unit length", H,
     "            y: base.x * shear + across * rise,",
     "            y: base.x * shear + across,", [SLAB, WHOLE]),
    ("a derived arc is scaled twice", P,
     "                        x: transform.matrix2 * circle.frame.x / scale,",
     "                        x: transform.matrix2 * circle.frame.x,", [DERIVED, WHOLE]),
    ("the boundary is chorded at its own scale", M,
     "        let chord = chord_error(options) / stretch.max(1.0);",
     "        let chord = chord_error(options);", [WHOLE]),
    ("the exact compiler reads every frame as a unit one", C,
     "                let shape = if is_identity(&shape) {",
     "                let shape = if true {", [WHOLE]),
    ("a reflecting boundary map is not reflected", C,
     "        let reflects = into.matrix2.determinant() < 0.0;",
     "        let reflects = false;", [WHOLE]),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
