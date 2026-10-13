"""Mutation probe for run stations at their basis's ends (#303, ADR 0082
amendment, ADR 0085 format 1.3).

Each mutant breaks one step a run station at its basis's start or end
rests on: in `axiolid-mesh-compile` the resolution of the end to the
kernel's own length (the end and the start moved, the length of an atomic
or a composite basis not the one its curve path states, in the swept
reading and in the flattened one), the order check after resolution, and
the window applied after it; in `axiolid-model` the validation (a start
not first, an end not last, an end stating a distance, the end ordered by
its unread distance, its offsets unchecked) and the wire format (`None`
written, a basis end read in a 1.2 payload, a payload never labelled 1.3,
the carrier's stations not looked at). Each must turn
`compile/tests/station_run_basis_end.rs` or the model's `station`/`wire`
tests red.

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import probe_lib

TIMEOUT = 1800
S = "crates/execution/compile/src/station.rs"
B = "crates/execution/compile/src/station/basis.rs"
V = "crates/representations/modeling/graph/src/validation.rs"
G = "crates/representations/modeling/graph/src/wire.rs"
M = "crates/representations/modeling/graph/src/station.rs"

END = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_run_basis_end"]
STATION = ["-p", "axiolid-model", "--test", "station"]
WIRE = ["-p", "axiolid-model", "--all-features", "--test", "wire"]

MUTANTS = [
    # Resolution against the kernel's own length.
    ('the end resolves short of the basis\'s length', S,
     '                Some(BasisEnd::End) => end,',
     '                Some(BasisEnd::End) => end - 1e-9,', [END]),
    ('the start resolves after the basis\'s start', S,
     '                Some(BasisEnd::Start) => 0.0,',
     '                Some(BasisEnd::Start) => 1e-9,', [END]),
    ('an atomic basis\'s end is not its whole curve\'s', B,
     '            Self::Atomic { curve, .. } => Ok(StationPiece::whole(*curve)?.length()),',
     '            Self::Atomic { curve, .. } => Ok(StationPiece::whole(*curve)?.length() * (1.0 - 1e-12)),',
     [END]),
    ('a composite basis\'s end is not its running sum', B,
     '            Self::Atomic { curve, .. } => Ok(StationPiece::whole(*curve)?.length()),\n'
     '            Self::Composite(composite) => Ok(composite.length()),',
     '            Self::Atomic { curve, .. } => Ok(StationPiece::whole(*curve)?.length()),\n'
     '            Self::Composite(composite) => Ok(composite.length() * (1.0 - 1e-12)),',
     [END]),
    ('a flattened run\'s atomic basis ends elsewhere', B,
     '            Self::Atomic(curve, ..) => Ok(StationPiece::whole(*curve)?.length()),',
     '            Self::Atomic(curve, ..) => Ok(StationPiece::whole(*curve)?.length() * (1.0 - 1e-12)),',
     [END]),
    ('a flattened run\'s composite basis ends elsewhere', B,
     '            Self::Atomic(curve, ..) => Ok(StationPiece::whole(*curve)?.length()),\n'
     '            Self::Composite(composite) => Ok(composite.length()),',
     '            Self::Atomic(curve, ..) => Ok(StationPiece::whole(*curve)?.length()),\n'
     '            Self::Composite(composite) => Ok(composite.length() * (1.0 - 1e-12)),',
     [END]),
    ('a flattened run reads its stated distances', B,
     '        super::resolved_run(|| base.length(), stations)?,',
     '        stations\n'
     '            .iter()\n'
     '            .map(|station| (station.station.distance, station.seam_window)),',
     [END]),
    # Order and windows after resolution.
    ('resolved stations out of order are not refused by name', S,
     '    if resolved.windows(2).any(|pair| pair[1].0 <= pair[0].0) {',
     '    if false {', [END]),
    ('the end\'s window is not applied after resolution', S,
     '            Ok((distance, station.seam_window))',
     '            Ok((distance, if station.basis_end.is_some() { 0.0 } else { station.seam_window }))',
     [END]),
    # Validation.
    ('a start station need not be first', V,
     '            BasisEnd::Start if index != 0 => {', '            BasisEnd::Start if false => {',
     [STATION]),
    ('an end station need not be last', V,
     '            BasisEnd::End if index != last => {', '            BasisEnd::End if false => {',
     [STATION]),
    ('an end station may state a distance', V,
     '        if station.station.distance != 0.0 {', '        if false {', [STATION, WIRE]),
    ('an end station is ordered by its unread distance', V,
     '                .filter(|station| station.basis_end != Some(BasisEnd::End));',
     '                .filter(|_| true);', [STATION, WIRE]),
    ('an end station\'s offsets are not checked', V,
     '        validate_station(&station.station)?;\n    }\n    Ok(())\n}',
     '    }\n    Ok(())\n}', [STATION]),
    # Format 1.3.
    ('an absent basis end is written', M,
     '        serde(default, skip_serializing_if = "Option::is_none")',
     '        serde(default)', [WIRE]),
    ('a basis end in a 1.2 payload is read', G,
     '    if version < BASIS_END_SINCE && uses_basis_end(node) {',
     '    if false && uses_basis_end(node) {', [WIRE]),
    ('a payload is never labelled 1.3', G,
     '    if nodes.iter().any(uses_basis_end) {', '    if false {', [WIRE]),
    ('an offset\'s run stations are not looked at', G,
     '            CurveRelation::OffsetsAtStations { stations, .. } => stations,',
     '            CurveRelation::OffsetsAtStations { .. } => &[],', [WIRE]),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
