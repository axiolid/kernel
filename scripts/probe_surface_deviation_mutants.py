"""Mutation probe for the surface-to-mesh bound of curved bodies (#231, #232,
#245).

Each mutant gives a direction of a doubly curved surface the whole chord
budget again, sizes a step at the wrong radius, drops a term of the bound
(the walls' twist, the taper, the section's reach, the directrix's own
sagitta), loosens the quad lemma, tilts a swept end, or clamps a budget it
must refuse, and must turn a test red. The pipe mutants (#232) do the same
to swept disks along chains of segments and arcs: they drop a term of a
piece's bound or the joint's measured shift, stop carrying the frame across
a joint, loosen the joint tolerance, cut fillets wrongly or let them
overrun their segments, and accept folded bends, horn tori or closed
polylines. The mitre mutants (#245) chord the ring without the mitre's
stretch, skip one side's projection, tilt the mitre plane, project
obliquely, drop the reach check or one end of it, and mitre reversals.
The #248 mutants turn the band beside a mitred arc the wrong way, let it
reach past its arc, share the wrong ring between two arcs, trim the
certified band or planar face the wrong way or drop the certificate, and
break loops: drop or flip the closing twist, forget the bends' share of
it, drop the closing joint or the closing fillet, leave the loft open,
drop the closing segment or accept a range on a closed polyline.
`surface_deviation.rs` samples the exact surfaces against the compiled
meshes; the `loft` unit tests check the quad lemma against a brute-force
distance.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/execution/compile/src/compiler.rs"
D = "crates/execution/compile/src/directrix.rs"
R = "crates/algorithms/construction/construct/src/revolve.rs"
S = "crates/algorithms/construction/construct/src/sweep.rs"
L = "crates/algorithms/construction/construct/src/loft.rs"
P = "crates/algorithms/reference/src/primitive.rs"
PIPE = "crates/algorithms/construction/construct/src/pipe.rs"
PIECES = "crates/execution/compile/src/directrix/pieces.rs"
BAND = "crates/algorithms/construction/construct/src/pipe/band.rs"
CERT = "crates/execution/compile/src/deviation/mitre.rs"
SURFACE = ["-p", "axiolid-mesh-compile", "--test", "surface_deviation"]
# The #248 tests alone, so a mutant there is not run against the whole file.
BESIDE = SURFACE + ["--", "beside", "after_a_bend"]
LOOPS = SURFACE + ["--", "closed", "twist"]
QUAD = ["-p", "axiolid-construct", "--lib", "loft"]

MUTANTS = [
    ('revolution profile chorded to the whole budget', C, '`axiolid_construct::revolve`).\n                let half = 0.5 * chord_error(options);', '`axiolid_construct::revolve`).\n                let half = chord_error(options);', [SURFACE]),
    ('revolution turned to the whole budget', C, '    Tolerance::new(0.5 * chord_error(options), options.tolerance().angular())', '    Tolerance::new(chord_error(options), options.tolerance().angular())', [SURFACE]),
    ('turn sized at the tube centre', R, '    let n = bounded_steps(rings, axis_origin, dir, angle, max_r, tolerance.linear())?;', '    let n = bounded_steps(rings, axis_origin, dir, angle, 0.8 * max_r, tolerance.linear())?;', [SURFACE]),
    ('revolution wall twist ignored', R, '        if sagitta + twist <= budget {', '        if sagitta <= budget {', [SURFACE]),
    # The revolution's step cap is guarded twice (the sagitta count and the
    # growth loop both refuse), so no single mutant there can pass the
    # refusal test; the primitive's single guard is probed below.
    ('tapered taper speed ignored', S, '    let curvature = angle * angle * max_r + 2.0 * angle.abs() * travel;', '    let curvature = angle * angle * max_r;', [SURFACE]),
    ('tapered wall twist ignored', S, '        if chord + twist <= budget {', '        if chord <= budget {', [SURFACE]),
    ('swept disk chorded to the whole budget', S, '    let rings = disk_rings(radius, inner_radius, 0.5 * chord)?;', '    let rings = disk_rings(radius, inner_radius, chord)?;', [SURFACE]),
    ('section reach ignored', S, '    directrix + reach * (1.0 - (0.5 * rotation).cos())', '    directrix', [SURFACE]),
    ('directrix sagitta ignored', S, '    let directrix = budget.min(0.5 * chord * (0.25 * turn).tan());', '    let directrix = 0.0;', [SURFACE]),
    ('end sections square to the end chords', S, '                ends.map_or(path[1] - path[0], |[start, _]| start)', '                path[1] - path[0]', [SURFACE]),
    ('smooth directrix never refined', S, '        if !smooth {\n            return Ok(stations);', '        if true {\n            return Ok(stations);', [SURFACE]),
    ('reversed trim keeps its tangents', D, '    (points, ends.map(|[start, end]| [-end, -start]))', '    (points, ends.map(|[start, end]| [end, start]))', [SURFACE]),
    ('conic end tangents not reported', D, '    Some([unit(start)?, unit(end)?])', '    Some([unit(start)?, unit(end)?]).filter(|_| false)', [SURFACE]),
    ('quad flattening bound quartered', L, '            best = best.min(g.abs());', '            best = best.min(0.25 * g.abs());', [QUAD, SURFACE]),
    ('quad mid-plane bound halved', L, '            best = best.min(2.0 * e.abs());', '            best = best.min(e.abs());', [QUAD, SURFACE]),
    ('quad convexity not checked', L, '        if convex_projection([a, b, d - n * g, c], n) {', '        if true {', [QUAD, SURFACE]),
    ('quad twist bound halved', L, '    let twist = (a - b - c + d).length() / 4.0;', '    let twist = (a - b - c + d).length() / 8.0;', [QUAD, SURFACE]),
    ('sphere stacks rounded down', P, '    let stacks = n.div_ceil(2).max(2);', '    let stacks = (n / 2).max(2);', [SURFACE]),
    ('sphere at the whole budget', P, '    let n = segments_within(r, 0.5 * tol)?;', '    let n = segments_within(r, tol)?;', [SURFACE]),
    ('torus equator at the whole budget', P, '    let n = segments_within(big + r, 0.5 * tol)?;', '    let n = segments_within(big + r, tol)?;', [SURFACE]),
    ('primitive segment cap clamps', P, '    if n == 4096 && radius', '    if false && radius', [SURFACE]),
    ('cylinder segments clamped', P, '"cylinder height")?;\n    let n = segments_within(r, tol)?;', '"cylinder height")?;\n    let n = segments(r, tol);', [SURFACE]),
    ('cone segments clamped', P, '"cone height")?;\n    let n = segments_within(r, tol)?;', '"cone height")?;\n    let n = segments(r, tol);', [SURFACE]),
    ('pipe disk chorded to the whole budget', PIPE, '            0.5 * self.chord * (0.5 * self.sharpest).cos(),', '            self.chord * (0.5 * self.sharpest).cos(),', [SURFACE]),
    ('pipe bends sampled to the whole budget', PIPE, '        let share = 0.5 * self.chord;\n        let band_tol', '        let share = self.chord;\n        let band_tol', [SURFACE]),
    ('pipe section reach ignored', PIPE, 'sweep::span_bound(&f[0], &f[1], b, reach)', 'sweep::span_bound(&f[0], &f[1], b, 0.0)', [SURFACE]),
    ('pipe bends never refined', PIPE, '        if worst <= share {', '        if true {', [SURFACE]),
    ('pipe joint shift ignored', PIPE, '                let moved = if k == 0 { shift } else { 0.0 };', '                let moved = 0.0;', [SURFACE]),
    ('pipe frame not carried across a joint', PIPE, '            |f| least_rotation(f.x, f.x.cross(f.y), tangent),', '            |_| seed_reference(tangent),', [SURFACE]),
    ('pipe corner tolerance quadrupled', PIPE, '    2.0 * (chord / (8.0 * radius)).min(1.0).asin()', '    2.0 * (chord / (2.0 * radius)).min(1.0).asin()', [SURFACE]),
    ('fillet cut by the sine', PIPE, '        let cut = fillet * (0.5 * theta).tan();', '        let cut = fillet * (0.5 * theta).sin();', [SURFACE]),
    ('fillet overrunning its segment accepted', PIPE, '        if needed > length + slack {', '        if false {', [SURFACE]),
    ('disk wider than the fillet accepted', PIPE, '        if radius > fillet {', '        if false {', [SURFACE]),
    ('disk as wide as the fillet not named', PIPE, '        if radius == fillet {', '        if false {', [SURFACE]),
    ('bend as tight as the disk accepted', PIPE, '                if bend == radius {', '                if false {', [SURFACE]),
    ('mitre ring without the stretch', PIPE, '            0.5 * self.chord * (0.5 * self.sharpest).cos(),', '            0.5 * self.chord,', [SURFACE]),
    ('mitre projected at the start only', PIPE, 'for (k, normal) in [(0, ends[0]), (last, ends[1])] {', 'for (k, normal) in [(0, ends[0]), (last, None)] {', [SURFACE]),
    ('mitre plane square to the incoming leg', PIPE, '        mitres.push(Some(bisector.normalize()));', '        mitres.push(Some(u));', [SURFACE]),
    ('mitre projected obliquely', PIPE, '.map(|q| *q - tangent * ((*q - frame.origin).dot(normal) / along))', '.map(|q| *q - tangent * (*q - frame.origin).dot(normal))', [SURFACE]),
    ('mitre reach not checked', PIPE, '        if reach >= length {', '        if false {', [SURFACE]),
    ('mitre reach of one end only', PIPE, '        let g = lean(ends[1]) - lean(ends[0]);', '        let g = lean(ends[1]) - u;', [SURFACE]),
    ('reversal mitred', PIPE, '        if theta >= core::f64::consts::PI - 1e-9 || bisector.length() <= 1e-9 {', '        if false {', [SURFACE]),
    ('reversed arc keeps its turn', PIPE, '                start: self.end_point(),\n                angle: -angle,', '                start: self.end_point(),\n                angle,', [SURFACE]),
    # #248: corners mitred beside an arc.
    ('band beside an outgoing arc turned away from the corner', PIPE, '(-axis * angle.signum(), path[i].start_point(), -n)', '(axis * angle.signum(), path[i].start_point(), -n)', BESIDE),
    ('band beside an incoming arc turned away from the corner', PIPE, '(axis * angle.signum(), path[i].end_point(), n)', '(-axis * angle.signum(), path[i].end_point(), n)', BESIDE),
    ('band continued past the corner without clearing the plane', BAND, '            >= margin\n    })?;', '            >= -radius\n    })?;', BESIDE),
    ('bands reaching past their arc accepted', PIPE, '            if need > angle.abs() {', '            if false {', BESIDE),
    ('two arcs share a projected ring', PIPE, '            (PathPiece::Arc { .. }, PathPiece::Arc { .. }) => {\n                self.bands[next][0].ok_or_else(not_an_arc)?.cut(&square)\n            }', '            (PathPiece::Arc { .. }, PathPiece::Arc { .. }) => Ok(mitred(&square, frame, normal)),', BESIDE),
    ('band vertices left on the far section', BAND, '.map(|(&p, &t)| self.turned(p, s * t))', '.map(|(&p, &t)| self.turned(p, 0.0 * s * t))', BESIDE),
    ('band trimmed to the wrong side of the plane', CERT, '        if height + moves <= 0.0 {\n            Coverage::Inside\n        } else if height - moves > 0.0 {', '        if height - moves > 0.0 {\n            Coverage::Inside\n        } else if height + moves <= 0.0 {', BESIDE),
    ('planar face taken where both tubes are', CERT, '            (Some(a), Some(b)) if a != b => Coverage::Inside,', '            (Some(a), Some(b)) if a == b => Coverage::Inside,', BESIDE),
    ('planar face reads a segment tube as a ball', CERT, '                within((q - *direction * q.dot(*direction)).length())', '                within(q.length())', BESIDE),
    ('mitre beside an arc not certified', CERT, '    if !plan.arc_mitres.is_empty() {', '    if false {', BESIDE),
    # #248: closed polylines.
    ('closing twist not turned back', PIPE, '                pipe.rate = pipe.twist / clear;', '                pipe.rate = 0.0;', LOOPS),
    ('closing twist turned the wrong way', PIPE, '        x.cross(x0).dot(t0).atan2(x.dot(x0))', '        x0.cross(x).dot(t0).atan2(x.dot(x0))', LOOPS),
    ('bends take no twist', PIPE, '            twist = self.rate * part.length();', '            twist = 0.0;', LOOPS),
    ('no closing joint', PIPE, '        let joints = if closed { m } else { m - 1 };', '        let joints = m - 1;', LOOPS),
    ('closing corner not filleted', PIPE, '    let first = if closed { 0 } else { 1 };', '    let first = 1;', LOOPS),
    ('loop lofted open', PIPE, '        let mut mesh = loft::loft(&rings, &stations, self.closed)?;', '        let mut mesh = loft::loft(&rings, &stations, false)?;', LOOPS),
    ('twist with no clear length accepted', PIPE, '                if clear <= 0.0 {', '                if false {', LOOPS),
    ('closing segment dropped', PIECES, '            let pieces: Vec<PathPiece> = (0..n)', '            let pieces: Vec<PathPiece> = (0..n - 1)', LOOPS),
    ('range on a closed polyline accepted', PIECES, '            if range.is_some() {', '            if false {', LOOPS),
]

probe_lib.main(MUTANTS, timeout=1200)
