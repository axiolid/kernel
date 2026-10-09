"""Mutation probe for nodes placed at stations (#264, ADR 0082 amendment).

Each mutant breaks one step a placement rests on: the axis mapping (local
x, y, z onto tangent, left lateral, up), the origin on the station's
offset point, the oriented (not base) frame, the seam side passed through,
the exactness claim and what the compilers do with it, the frame a
station along a placed curve carries and the refusal of a tilted
plan-measured one (#285 reads a tilted arc-length one in its own frame), the
directrix moved by the placement, and the graph's validation (family,
dimension, station checks, references). Each must turn a test red:
`evaluate/tests/station_placement.rs` checks the placement against closed
forms; `compile/tests/station_placement.rs` resolves placements on a line,
an arc, an elevated curve and a polyline seam against frames computed by
hand, meshes and exactly compiles placed solids and sweeps placed curves,
and checks the refusals; `graph/tests/station.rs` checks validation and
the round trip.

Not probed, because no test can tell them apart: the nesting depth limit
(no test nests 64 placements) and a 2D relation as a placed directrix's
source (its refusal only renames an error the reader already raises).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
E = "crates/algorithms/parametric/evaluate/src/station.rs"
C = "crates/execution/compile/src/station.rs"
B = "crates/execution/compile/src/station/basis.rs"
K = "crates/execution/compile/src/compiler.rs"
X = "crates/execution/compile/src/exact.rs"
D = "crates/execution/compile/src/directrix.rs"
V = "crates/representations/modeling/graph/src/validation.rs"
N = "crates/representations/modeling/graph/src/node.rs"
EVAL = ["-p", "axiolid-evaluate", "--test", "station_placement"]
PLACE = ["-p", "axiolid-mesh-compile", "--test", "station_placement"]
MODEL = ["-p", "axiolid-model", "--test", "station"]

MUTANTS = [
    ('local y onto up, z onto lateral', E,
     '        Transform3::from_cols(self.tangent, self.lateral, self.up, self.point)',
     '        Transform3::from_cols(self.tangent, self.up, self.lateral, self.point)',
     [EVAL, PLACE]),
    ('local y onto the right', E,
     '        Transform3::from_cols(self.tangent, self.lateral, self.up, self.point)',
     '        Transform3::from_cols(self.tangent, -self.lateral, self.up, self.point)',
     [EVAL, PLACE]),
    ('moved frame keeps its tangent', E,
     '            tangent: rigid.transform_vector3(self.tangent),',
     '            tangent: self.tangent,', [EVAL, PLACE]),
    ('every 2D frame claimed exact', E,
     '    matches!(curve, Curve2::Line(_))',
     '    true', [EVAL, PLACE]),
    ('every 3D frame claimed exact', E,
     '    matches!(curve, Curve3::Line(_))',
     '    true', [EVAL, PLACE]),
    ('placement ignores its seam side', C,
     '    let section = basis.section_on(station.distance, *frame, *seam)?;',
     '    let section = basis.section_on(station.distance, *frame, SeamSide::Outgoing)?;',
     [PLACE]),
    ('placement origin without the offsets', C,
     '    transform.translation = resolved.point;\n',
     '', [PLACE]),
    ('placement in the unturned frame', C,
     '    let mut transform = turned.placement();',
     '    let mut transform = resolved.section.placement();', [PLACE]),
    ('a placed basis is not carried', B,
     '                    Some(rigid) => section.carried(*rigid, curve.convention())?,',
     '                    Some(_) => section,', [PLACE]),
    ('a tilted plan-measured placed basis accepted', E,
     '        if measure == DistanceConvention::PlanDistance {',
     '        if false {', [PLACE]),
    ('a placed basis loses its seams', B,
     '                StationCurve::Two(curve) => exact_station_seams2(curve),',
     '                StationCurve::Two(_) => Ok(Vec::new()),', [PLACE]),
    ('an inexact placement reported bounded', K,
     '                if !placement.exact {',
     '                if false {', [PLACE]),
    ('an inexact placement compiled exactly', X,
     '                if !placement.exact {',
     '                if false {', [PLACE]),
    ('placed directrix not moved', D,
     '                    .map(|p| transform.transform_point3(p))',
     '                    .map(|p| p)', [PLACE]),
    ('a placed curve accepted in a 2D slot', V,
     '                if dimension != CurveDimension::Three\n                    || ',
     '                if ', [MODEL]),
    ('a placed non-curve accepted as a curve', V,
     '                    || !(curve_has_dimension(placed.source, nodes, CurveDimension::Two)\n'
     '                        || curve_has_dimension(placed.source, nodes, CurveDimension::Three))\n',
     '', [MODEL]),
    ('a placed profile still a profile', V,
     '                    if !matches!(self, Self::Profile | Self::OpenProfile) =>',
     '                    if true =>', [MODEL]),
    ('a placed solid not a solid', V,
     '                GeometryNode::InstanceAtStation(placed)\n'
     '                    if !matches!(self, Self::Profile | Self::OpenProfile) =>',
     '                GeometryNode::InstanceAtStation(placed)\n'
     '                    if false =>', [MODEL]),
    ('placement station not validated', V,
     '            validate_station(&station.station.station)?;\n',
     '', [MODEL]),
    ('placement orientation not validated', V,
     '            validate_orientation(&station.orientation)\n',
     '            Ok(())\n', [MODEL]),
    ('placement basis not checked', V,
     '            expect_reference(nodes, station.station.basis, ExpectedReference::Curve)?;\n',
     '', [MODEL]),
    ('placement does not reference its basis', N,
     '                references.extend([value.source, value.station.station.basis]);',
     '                references.push(value.source);', [MODEL]),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
    ).returncode


survivors = []
for name, rel, old, new, targets in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    code = 0
    try:
        for target in targets:
            try:
                code = run(target)
            except subprocess.TimeoutExpired:
                code = -1
            if code != 0:
                break
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
