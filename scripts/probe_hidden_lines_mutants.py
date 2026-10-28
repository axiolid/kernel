"""Mutation probe for the exact hidden-line provider (#138, ADR 0090).

Each mutant breaks one decision the provider owes the contract (what hides
what, which edge is which kind, how pieces are cut and rounded, what is
refused) and must turn a test of the provider red: the conformance suite, the
independent ray-cast oracle, the exactness test or a unit test.
"""
import probe_lib

OCC = "crates/algorithms/discrete/hidden-lines/src/occlusion.rs"
EDG = "crates/algorithms/discrete/hidden-lines/src/edges.rs"
PAR = "crates/algorithms/discrete/hidden-lines/src/param.rs"
DEC = "crates/algorithms/discrete/hidden-lines/src/decide.rs"
OUT = "crates/algorithms/discrete/hidden-lines/src/output.rs"
BRD = "crates/algorithms/discrete/hidden-lines/src/broad.rs"
PRO = "crates/algorithms/discrete/hidden-lines/src/provider.rs"
TESTS = ["-p", "axiolid-hidden-lines"]

MUTANTS = [
    # What hides what.
    ('touching occludes', OCC,
     '            lift: Lift::Own,\n        },\n        true,',
     '            lift: Lift::Own,\n        },\n        false,', TESTS),
    ('behind occludes', OCC,
     '            neg: !clockwise,\n            lift: Lift::Own,',
     '            neg: clockwise,\n            lift: Lift::Own,', TESTS),
    ('open image coverage', OCC,
     '                neg: clockwise,\n            },\n            false,\n        );\n        match solver.clip(span, &side)? {',
     '                neg: clockwise,\n            },\n            true,\n        );\n        match solver.clip(span, &side)? {', TESTS),
    ('image orientation ignored', OCC,
     '                neg: clockwise,\n            },\n            false,\n        );\n        match solver.clip(span, &side)? {',
     '                neg: false,\n            },\n            false,\n        );\n        match solver.clip(span, &side)? {', TESTS),
    ('edge-on faces never occlude', OCC,
     '            for span in edge_on(solver, scene, edge, tri, kept)? {',
     '            for span in edge_on(solver, scene, edge, tri, kept)?.into_iter().take(0) {', TESTS),
    ('edge-on fibre depth ignored', OCC,
     'let above = Constraint::new(Kind::SegmentAbove { from, to }, true);',
     'let above = Constraint::new(Kind::SegmentAtOrBelow { from, to, near: f64::MAX }, true);', TESTS),
    ('edge-on fibre reads one edge only', OCC,
     '        let Some(within) = within else { continue };',
     '        let Some(within) = within else { continue };\n        if !in_front.is_empty() {\n            break;\n        }', TESTS),
    ('clipped occluder still occludes', OCC,
     '    if let Some(near) = scene.near {\n        // The triangle\'s point must not be removed by the near plane.',
     '    if let Some(near) = scene.near.filter(|_| false) {\n        // The triangle\'s point must not be removed by the near plane.', TESTS),
    ('far ignored', OCC,
     '    if let Some(far) = scene.far {',
     '    if let Some(far) = scene.far.filter(|_| false) {', TESTS),
    ('region ignored', OCC,
     '    if let Some(region) = scene.region {',
     '    if let Some(region) = scene.region.filter(|_| false) {', TESTS),
    ('point edges drawn', OCC,
     '        if axis == 1 {\n            return Ok(Resolved::Point);',
     '        if axis == 7 {\n            return Ok(Resolved::Point);', TESTS),
    ('candidate limit unchecked', OCC,
     '        if *pairs >= max_pairs {',
     '        if false && *pairs >= max_pairs {', TESTS),
    # Which edge is which kind.
    ('silhouette needs strictly opposite facings', EDG,
     '                if tris[t1].sigma != tris[t2].sigma {',
     '                if tris[t1].sigma.reverse() == tris[t2].sigma\n                    && tris[t1].sigma != Ordering::Equal\n                {', TESTS),
    ('crease threshold flipped', EDG,
     '(decide(&crease, escalations).map_err(refuse)? == Ordering::Less)',
     '(decide(&crease, escalations).map_err(refuse)? == Ordering::Greater)', TESTS),
    ('silhouette drawn when deselected', EDG,
     '                    selection.silhouette.then_some(EdgeKind::Silhouette)',
     '                    Some(EdgeKind::Silhouette).filter(|_| true)', TESTS),
    # How the edge is cut.
    ('same-visibility pieces not merged', PAR,
     'while next < list.len() && self.cmp(list[next].lo, merged.hi)? != Ordering::Greater {',
     'while next < list.len() && self.cmp(list[next].lo, merged.hi)? == Ordering::Less {', TESTS),
    ('root compare ignores denominator signs', PAR,
     '            if diff == Ordering::Less {\n                order = order.reverse();',
     '            if false && diff == Ordering::Less {\n                order = order.reverse();', TESTS),
    ('identically zero constraint ignores strictness', PAR,
     '        (Ordering::Equal, Ordering::Equal) if constraint.strict => Holds::Nowhere,',
     '        (Ordering::Equal, Ordering::Equal) if false => Holds::Nowhere,', TESTS),
    ('undecided filter read as zero', DEC,
     '    *escalations += 1;\n    match expr.sign_in::<Dyadic>() {',
     '    *escalations += 1;\n    match Some(Sign::Zero).filter(|_| expr.sign_in::<Dyadic>().is_some()) {', TESTS),
    # What is written out.
    ('ties round up', OUT,
     'std::cmp::Ordering::Equal => sticky || kept.bit(0),',
     'std::cmp::Ordering::Equal => true,', TESTS),
    ('sticky remainder ignored', OUT,
     'std::cmp::Ordering::Equal => sticky || kept.bit(0),',
     'std::cmp::Ordering::Equal => kept.bit(0),', TESTS),
    ('output limit unchecked', PRO,
     '            if segments.len() >= request.limits.max_output_segments {',
     '            if false && segments.len() >= request.limits.max_output_segments {', TESTS),
    ('near clip through a cut solid accepted', PRO,
     '        if beyond && kept {',
     '        if false && beyond && kept {', TESTS),
    # No mutant makes `cut_by_near` report an uncut solid as cut: with caps
    # that only filters work. A solid wholly kept has no triangle beyond the
    # plane, and a solid wholly beyond it is crossed an even number of times
    # (the ray enters and leaves), so neither ever has an odd parity.
    ('collapsed neighbours not merged', PRO,
     '                if last.visibility == visibility {',
     '                if false && last.visibility == visibility {', TESTS),
    ('edge-on wall ignores near', OCC,
     '    let Some(near) = scene.near else {',
     '    let Some(near) = scene.near.filter(|_| false) else {', TESTS),
    # Near-plane caps.
    ('caps ignored', OCC,
     '        for span in hidden_by_caps(solver, scene, edge, near, kept, pairs, max_pairs)? {',
     '        for span in hidden_by_caps(solver, scene, edge, near, kept, pairs, max_pairs)?\n            .into_iter()\n            .filter(|_| false)\n        {', TESTS),
    ('perturbation replaced by closed coverage', PAR,
     '        Ok((side == Ordering::Greater).then_some(span))',
     '        Ok((side == Ordering::Greater || true).then_some(span))', TESTS),
    # No mutant for the tie "a triangle at the near plane counts as crossed":
    # over that triangle's image every point below the plane is hidden by the
    # triangle itself, and elsewhere it adds no toggle, so no drawing differs.
    ('perturbation ignores orientation', PAR,
     '        let side = if neg { side.reverse() } else { side };',
     '        let side = if neg && false { side.reverse() } else { side };', TESTS),
    ('caps applied to sheets', PRO,
     '            (Some(near), MeshRole::Solid) => {',
     '            (Some(near), _) => {', TESTS),
    # The broad phase.
    ('BVH probe misses touching boxes', BRD,
     '                    max: Point3::new(query.image.hi[0], query.image.hi[1], query.w_hi),',
     '                    max: Point3::new(query.image.hi[0].next_down(), query.image.hi[1], query.w_hi),', TESTS),
    ('BVH probe misses touching boxes at y max', BRD,
     '                    max: Point3::new(query.image.hi[0], query.image.hi[1], query.w_hi),',
     '                    max: Point3::new(query.image.hi[0], query.image.hi[1].next_down(), query.w_hi),', TESTS),
    ('BVH probe misses touching boxes at depth max', BRD,
     '                    max: Point3::new(query.image.hi[0], query.image.hi[1], query.w_hi),',
     '                    max: Point3::new(query.image.hi[0], query.image.hi[1], query.w_hi.next_down()),', TESTS),
]

probe_lib.main(MUTANTS, timeout=3600, reject_uncompiled=True)
