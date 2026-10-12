"""Mutation probe for certified joint distances and seam-snapping windows
(#294, ADR 0082 amendment).

Each mutant breaks one step the window or the certificate rests on: in
`axiolid-evaluate`'s window module the snap itself, the window's width,
both named refusals and what counts as inside; in its certify module the
bound on a chord's excess over its arc, the midpoint rule's error term,
reading a numerical length as exact, the grade a leaning frame adds to an
offset's plan speed, the curvature's sign, the tolerance check and the
certified interval a seam is snapped by; in the composite which seams
count; in the provider and the contract the query and its default
refusal; in the conformance suite the containment check; in the compiler
the window applied at all; and in `axiolid-model` the window's validation
and its wire format (omitted at zero, refused in a 1.0 payload, and
lifting the payload to 1.1 from either carrier, a `-0.0` window too,
#297). For run stations (#301): the window ignored by a section run, an
offset curve or an offset used as a basis, snapped stations out of order
read, an interior station's side not read, the first station's side read,
a station's side overriding the mitre at its knot, each carrier's window
validation, a default side written, and format 1.2's feature test (the
1.1 refusal, the label, the variant, a section's side, a `-0.0` window,
each carrier). Each must turn `evaluate/tests/station_certify.rs`,
`curve_evaluate_path.rs`, `compile/tests/station_window.rs`,
`compile/tests/station_run_window.rs` or the model's `station`/`wire`
tests red.

Not probed, because no test can tell them apart: halving a bound the
certificate holds with room to spare (the derivative suprema are loose by
design, so tightening one by a constant still contains the truth), the
compensated summation and the rounding of the chord chain (they move a
bound by a few units in the last place), and the panel budget.

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 1800
W = "crates/algorithms/parametric/evaluate/src/station/window.rs"
C = "crates/algorithms/parametric/evaluate/src/station/certify.rs"
K = "crates/algorithms/parametric/evaluate/src/station/composite.rs"
P = "crates/algorithms/parametric/evaluate/src/provider.rs"
T = "crates/contracts/operations/curve-evaluate/src/contract.rs"
F = "crates/contracts/operations/curve-evaluate/src/conformance.rs"
S = "crates/execution/compile/src/station.rs"
B = "crates/execution/compile/src/station/basis.rs"
V ="crates/representations/modeling/graph/src/validation.rs"
G = "crates/representations/modeling/graph/src/wire.rs"
M = "crates/representations/modeling/graph/src/station.rs"

CERT = ["-p", "axiolid-evaluate", "--all-features", "--test", "station_certify"]
PATH = ["-p", "axiolid-evaluate", "--all-features", "--test", "curve_evaluate_path"]
COMP = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_window"]
RUN = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_run_window"]
STATION = ["-p", "axiolid-model", "--test", "station"]
WIRE = ["-p", "axiolid-model", "--all-features", "--test", "wire"]

MUTANTS = [
    # The window.
    ('a seam within the window is not snapped onto', W,
     '    Ok(inside.unwrap_or(distance))', '    Ok(distance)', [CERT, COMP]),
    ('the window is half as wide', W,
     '    let (lo, hi) = (distance - window, distance + window);',
     '    let (lo, hi) = (distance - 0.5 * window, distance + 0.5 * window);', [CERT, COMP]),
    ('two seams in one window are not refused', W,
     '    if ambiguous {', '    if false {', [CERT, COMP]),
    ('a seam on the window\'s edge is not refused', W,
     '    if straddled {', '    if false {', [CERT, COMP]),
    ('a seam on the window\'s edge counts as inside', W,
     '        if interval.lower >= lo && interval.upper <= hi {',
     '        if interval.upper >= lo && interval.lower <= hi {', [CERT, COMP]),
    ('smooth seams are snapped onto', C,
     '        .filter(|seam| !seam.smooth)\n', '', [CERT]),
    ('a joint is snapped by the kernel\'s estimate, not its certificate', K,
     '                    distance: self.starts[k],\n                    interval: run.interval(),',
     '                    distance: self.starts[k],\n'
     '                    interval: DistanceInterval::point(self.starts[k]),', [CERT, COMP]),
    ('the compiler ignores the window', S,
     '    let distance = basis.snapped_distance(station.distance, *seam_window)?;',
     '    let distance = station.distance + 0.0 * *seam_window;', [COMP]),
    # The certificate.
    ('a chord is taken for its arc', C,
     '                chord * (1.0 + 4.0 * Scalar::EPSILON) * 2.0 * half * half / turn.cos()',
     '                0.0 * half', [CERT]),
    ('the midpoint rule claims no error', C,
     '                h * h * h / 24.0 * bend + VALUE_ROUNDING * h * speed_mid',
     '                0.0 * bend + VALUE_ROUNDING * h * speed_mid', [CERT]),
    ('a numerical offset\'s length is read as exact', C,
     '            if offset.measure_is_exact() {\n                return Ok(DistanceInterval::point(piece.length()));',
     '            if true {\n                return Ok(DistanceInterval::point(piece.length()));',
     [CERT]),
    ('a leaning frame adds no grade to the plan speed', C,
     '        let grade = elevation.filter(|_| leans && (vertical != zero || longitudinal != zero));',
     '        let grade = elevation.filter(|_| false && leans && (vertical != zero || longitudinal != zero));',
     [CERT]),
    ('the curvature\'s sign is flipped', C,
     '                .map(|k| sign * k)', '                .map(|k| -sign * k)', [CERT]),
    ('a zero tolerance is accepted', C,
     '    if tolerance.is_finite() && tolerance > 0.0 {', '    if tolerance.is_finite() {',
     [CERT]),
    # The contract.
    ('the provider does not certify', P,
     '        CompositeBasis::from_path(path)?.joint_intervals(tolerance)',
     '        CompositeBasis::from_path(path)?.joint_intervals(tolerance * 0.0 + 1.0)?;\n'
     '        Err(GeomError::InvalidInput("no".into()))', [CERT, PATH]),
    ('the default refuses as a path, not a joint query', T,
     '            input: PATH_JOINT_DISTANCES_UNSUPPORTED,',
     '            input: CURVE_PATH_UNSUPPORTED,', [PATH]),
    ('the conformance suite accepts an interval missing the truth', F,
     '    if interval.lower > truth + rounding || interval.upper < truth - rounding {',
     '    if false {', [PATH]),
    # The graph.
    ('a negative window is accepted', V,
     '    if window < 0.0 {', '    if false {', [STATION, WIRE]),
    ('a zero window is written', M,
     '    window.to_bits() == 0', '    false', [WIRE]),
    ('a window in a 1.0 payload is read', G,
     '    if version < SEAM_WINDOW_SINCE && uses_seam_window(node) {',
     '    if false && uses_seam_window(node) {', [WIRE]),
    ('the window of an oriented station node is ignored', G,
     '        GeometryNode::OrientedCurveStation(station) => core::slice::from_ref(station),',
     '        GeometryNode::OrientedCurveStation(_) => &[],', [WIRE]),
    ('the window of an instance at a station is ignored', G,
     '        GeometryNode::InstanceAtStation(placed) => core::slice::from_ref(&placed.station),',
     '        GeometryNode::InstanceAtStation(_) => &[],', [WIRE]),
    ('a -0.0 window is labelled 1.0', G,
     '        .any(|station| station.seam_window.to_bits() != 0)',
     '        .any(|station| station.seam_window != 0.0)', [WIRE]),
    # Run stations (#301): snapping, the order after it, the side, the
    # mitre, validation and format 1.2.
    ('a section run ignores its windows', S,
     '            .map(|section| (section.station.distance, section.seam_window)),',
     '            .map(|section| (section.station.distance, 0.0 * section.seam_window)),',
     [RUN]),
    ('an offset curve ignores its stations\' windows', S,
     '            .map(|station| (station.station.distance, station.seam_window)),',
     '            .map(|station| (station.station.distance, 0.0 * station.seam_window)),',
     [RUN]),
    ('an offset basis ignores its stations\' windows', B,
     '            .map(|station| (station.station.distance, station.seam_window)),',
     '            .map(|station| (station.station.distance, 0.0 * station.seam_window)),',
     [RUN]),
    ('snapped stations out of order are read', S,
     '    if distances.windows(2).any(|pair| pair[1] <= pair[0]) {\n'
     '        return Err(GeomError::InvalidInput(SNAPPED_NOT_INCREASING.into()));',
     '    if false {\n'
     '        return Err(GeomError::InvalidInput(SNAPPED_NOT_INCREASING.into()));', [RUN]),
    ('an interior station\'s side is not read', S,
     '                .map_or(SeamSide::Outgoing, |(_, side)| *side)',
     '                .map_or(SeamSide::Outgoing, |_| SeamSide::Outgoing)', [RUN]),
    ('the run\'s first station reads its own side', S,
     '        sides.retain(|(at, _)| *at > first);', '', [RUN]),
    ('a station\'s side overrides the mitre at its knot', S,
     '        if let Some((_, mitre)) = self.mitres.iter().find(|(knot, _)| *knot == s) {',
     '        if let Some((_, mitre)) = self\n'
     '            .mitres\n'
     '            .iter()\n'
     '            .find(|(knot, _)| *knot == s && !self.sides.iter().any(|(at, _)| *at == s))\n'
     '        {', [RUN]),
    ('an offset run\'s window is not validated', V,
     '                validate_seam_window(station.seam_window)?;\n', '', [STATION]),
    ('an open section\'s window is not validated', V,
     '                validate_seam_window(section.seam_window)?;\n'
     '            }\n'
     '            validate_station_run(sections.iter().map(|section| &section.station))?;\n'
     '            validate_section_tags(sections.iter().map(|section| section.tags.as_slice()), true)',
     '            }\n'
     '            validate_station_run(sections.iter().map(|section| &section.station))?;\n'
     '            validate_section_tags(sections.iter().map(|section| section.tags.as_slice()), true)',
     [STATION]),
    ('a closed section\'s window is not validated', V,
     '                validate_seam_window(section.seam_window)?;\n'
     '            }\n'
     '            validate_station_run(sections.iter().map(|section| &section.station))?;\n'
     '            validate_section_tags(\n',
     '            }\n'
     '            validate_station_run(sections.iter().map(|section| &section.station))?;\n'
     '            validate_section_tags(\n', [STATION]),
    ('a default run side is written', M,
     '    *side == SeamSide::Outgoing', '    false', [WIRE]),
    ('a 1.2 feature in a 1.1 payload is read', G,
     '    if version < RUN_SEAMS_SINCE {', '    if false {', [WIRE]),
    ('a payload is never labelled 1.2', G,
     '    if nodes.iter().any(|node| run_seams(node).is_some()) {', '    if false {', [WIRE]),
    ('an offset at run stations is not a 1.2 variant', G,
     '            CurveRelation::OffsetsAtStations { .. } => Some(RunSeams::Kind("OffsetsAtStations")),',
     '            CurveRelation::OffsetsAtStations { .. } => None,', [WIRE]),
    ('a section\'s side does not select 1.2', G,
     '        if section.seam != SeamSide::Outgoing {', '        if false {', [WIRE]),
    ('a -0.0 section window is labelled 1.1', G,
     '        } else if section.seam_window.to_bits() != 0 {',
     '        } else if section.seam_window != 0.0 {', [WIRE]),
    ('the sections of a solid carry nothing', G,
     '            SolidOperation::SectionsAtStations { sections, .. } => sections,',
     '            SolidOperation::SectionsAtStations { .. } => &[],', [WIRE]),
    ('the sections of a surface carry nothing', G,
     '            SurfaceRelation::OpenSectionsAtStations { sections, .. } => sections,',
     '            SurfaceRelation::OpenSectionsAtStations { .. } => &[],', [WIRE]),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
