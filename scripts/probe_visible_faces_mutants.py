"""Mutation probe for the exact visible-face provider (#383, ADR 0094).

Each mutant breaks one decision the provider owes the contract (which face is
frontmost, how ties, caps, edge-on triangles, clips and the crop behave, how
regions are traced, snap-rounded, attributed and ordered, what is refused) and must
turn a test of the provider red: the conformance suite, the independent
ray-casting oracle or a unit test.
"""
import probe_lib

FAC = "crates/algorithms/discrete/hidden-lines/src/faces/"
MOD = "crates/algorithms/discrete/hidden-lines/src/faces.rs"
LAB = FAC + "label.rs"
SEG = FAC + "segments.rs"
INT = FAC + "intersect.rs"
REG = FAC + "regions.rs"
SID = FAC + "side.rs"
PTS = "crates/algorithms/planar/exact-arrangement/src/points.rs"
ARR = "crates/algorithms/planar/exact-arrangement/src/arrangement.rs"
PRE = "crates/algorithms/discrete/hidden-lines/src/prepare.rs"
PRO = FAC + "provider.rs"
GRP = "crates/algorithms/discrete/hidden-lines/src/groups.rs"
ARR_TESTS = ["-p", "axiolid-exact-arrangement"]
TESTS = ["-p", "axiolid-hidden-lines"]

MUTANTS = [
    # Which face is frontmost.
    ('equal depth occludes: one tied key kept', LAB,
     '            Ordering::Equal => best.push(surface),',
     '            Ordering::Equal => {}', TESTS),
    ('coverage sign ignores sigma', LAB,
     '        let sign = if clockwise { sign.reverse() } else { sign };',
     '        let sign = if clockwise && false { sign.reverse() } else { sign };', TESTS),
    ('edge-on triangles count as covering', MOD,
     '        self.planes[t].is_some() && !self.removed[t]',
     '        !self.removed[t]', TESTS),
    ('edge-on occludes in the prune', MOD,
     '            edge_on_occludes: false,',
     '            edge_on_occludes: true,', TESTS),
    ('near ignored in kept', LAB,
     '            kept &= sign_at(side, &plane.below(near), &mut work.escalations)? != Ordering::Less;',
     '            kept &= sign_at(side, &plane.below(near), &mut work.escalations)? != Ordering::Less || true;', TESTS),
    ('far ignored in kept', LAB,
     '            kept &= sign_at(side, &plane.above(far), &mut work.escalations)? != Ordering::Less;',
     '            kept &= sign_at(side, &plane.above(far), &mut work.escalations)? != Ordering::Less || true;', TESTS),
    ('crop ignored', LAB,
     '        if sign_at(side, form, &mut work.escalations)? == Ordering::Less {',
     '        if false && sign_at(side, form, &mut work.escalations)? == Ordering::Less {', TESTS),
    ('cap parity always even', LAB,
     '            odd = !odd;',
     '            odd = false;', TESTS),
    ('a triangle exactly at near counts as a cap crossing', LAB,
     '            && sign_at(side, &plane.below(near), &mut work.escalations)? == Ordering::Less',
     '            && sign_at(side, &plane.below(near), &mut work.escalations)? != Ordering::Greater', TESTS),
    ('lexicographic side point drops its third term', SID,
     '    for order in 0..3 {',
     '    for order in 0..2 {', TESTS),
    # The arrangement's input.
    ('intersection segments skipped', "crates/algorithms/discrete/hidden-lines/src/faces/provider.rs",
     '        &mut segments,\n    )?;\n    let arrangement = arrange(',
     '        &mut segments,\n    )?;\n    segments.truncate(outline_segments);\n    let arrangement = arrange(', TESTS),
    ('near-clip segments of edge-on solid triangles skipped', SEG,
     '        if !(faces.covers(t) || in_cap[t]) {',
     '        if !faces.covers(t) {', TESTS),
    ('pair limit unchecked for intersections', INT,
     '            work.budget.count().map_err(|e| faces.blame(a, e))?;',
     '            let _ = &work.budget;', TESTS),
    ('pair limit unchecked for labels', LAB,
     '        work.budget.count()?;\n        if !covers(faces, side, t, work)? {',
     '        if !covers(faces, side, t, work)? {', TESTS),
    # Regions.
    ('collinear-vertex rule applies at three or more pieces', REG,
     '        if let [a, b] = others[..] {',
     '        if let [a, b, ..] = others[..] {', TESTS),
    ('rounding through f64 division instead of round_ratio', PTS,
     '        match (round_ratio(x, w), round_ratio(y, w)) {',
     '        match (Some(x.to_f64() / w.to_f64()), Some(y.to_f64() / w.to_f64())) {', TESTS),
    ('the group is ignored', MOD,
     '            group_of.extend(resolve_groups(entry, options, &mut work.escalations)?);',
     '            group_of.extend(vec![FaceGroup(0); range.triangles.len()]);', TESTS),
    ('floating components not attributed', REG,
     '                    let placed = place(points, &outlines, at, work).map_err(blame)?;',
     '                    let _ = at;\n                    let placed = 0;', TESTS),
    ('region limit unchecked', REG,
     '            if built.regions.len() >= limits.max_regions {',
     '            if false && built.regions.len() >= limits.max_regions {', TESTS),
    ('region vertex limit unchecked', REG,
     '            if vertex_total > limits.max_region_vertices {',
     '            if false && vertex_total > limits.max_region_vertices {', TESTS),
    ('holes unsorted', REG,
     '    holes.sort_by(|a, b| sequence_order(a, b));',
     '    holes.sort_by(|_, _| std::cmp::Ordering::Equal);', TESTS),
    ('outer ring reversed', REG,
     '            let outer = positions(&rings[found.outer])?;',
     '            let mut outer = positions(&rings[found.outer])?;\n            outer.reverse();', TESTS),
    # Snap rounding of the regions (ADR 0100); the snap itself is probed by
    # `probe_snap_rounding_mutants.py`.
    ('a ring is netted with the rest of its class', REG,
     '            let own: Vec<(usize, bool)> = chain\n                .iter()\n                .zip(element_rings)\n                .filter(|&(_, &of)| of == ring)',
     '            let own: Vec<(usize, bool)> = chain\n                .iter()\n                .zip(element_rings)\n                .filter(|&(_, &of)| of == of || of == ring)', TESTS),
    ('snapped regions list no triangles', REG,
     '                        listed.extend(&source.triangles[region]);',
     '                        let _ = region;', TESTS),
    ('collapsed regions are not counted', REG,
     '                built.collapsed_regions += 1;',
     '                let _ = &built;', TESTS),
    ('collapsed holes are not counted', REG,
     '                    .filter(|&&hole| !ring_survives[hole])',
     '                    .filter(|&&hole| false && !ring_survives[hole])', TESTS),
    ('the displacement is not copied to the evidence', PRO,
     '    evidence.snap_displacement = built.snap_displacement;',
     '', TESTS),
    ('the rerouted count is not copied to the evidence', PRO,
     '    evidence.rerouted_edges = built.rerouted_edges;',
     '', TESTS),
    ('the collapsed counts are not copied to the evidence', PRO,
     '    evidence.collapsed_regions = built.collapsed_regions;',
     '', TESTS),
    # Planar face groups (#384, ADR 0106).
    ('planar groups ignore the facing test', GRP,
     '    Ok(decide(&dot, escalations)? == Ordering::Greater)',
     '    Ok(decide(&dot, escalations)? != Ordering::Less || true)', TESTS),
    ('planar groups join across every interior edge', GRP,
     '            || !same_plane_and_facing(&points(t1), &points(t2), escalations).map_err(refuse)?',
     '            || false && !same_plane_and_facing(&points(t1), &points(t2), escalations).map_err(refuse)?', TESTS),
    ('planar group id is the largest triangle', GRP,
     '        parent[r1.max(r2) as usize] = r1.min(r2);',
     '        parent[r1.min(r2) as usize] = r1.max(r2);', TESTS),
    ('planar groups never join', GRP,
     '        if t1 == t2\n',
     '        if true\n            || t1 == t2\n', TESTS),
    # What is refused.
    ('mesh limit unchecked', PRO,
     '            max_meshes: limits.max_meshes,',
     '            max_meshes: usize::MAX,', TESTS),
    ('vertex limit unchecked', PRO,
     '            max_source_vertices: limits.max_source_vertices,',
     '            max_source_vertices: usize::MAX,', TESTS),
    ('triangle limit unchecked', PRO,
     '            max_source_triangles: limits.max_source_triangles,',
     '            max_source_triangles: usize::MAX,', TESTS),
    ('mesh limit named as vertices', PRE,
     '        meshes: "visible-face meshes",',
     '        meshes: "visible-face source vertices",', TESTS),
    ('vertex limit named as triangles', PRE,
     '        vertices: "visible-face source vertices",',
     '        vertices: "visible-face source triangles",', TESTS),
    ('triangle limit named as meshes', PRE,
     '        triangles: "visible-face source triangles",',
     '        triangles: "visible-face meshes",', TESTS),
    ('intersection segments not charged as they come', INT,
     '                    if declared.saturating_add(charged) > budget {',
     '                    if false && declared.saturating_add(charged) > budget {', TESTS),
    # The arrangement.
    ('collinear segments of one end are not cut', ARR,
     '                if points.orient(shared, near, far)? == Sign::Zero {',
     '                if false && points.orient(shared, near, far)? == Sign::Zero {', ARR_TESTS),
    ('collinear segments of one end cut only the first way', ARR,
     '                    } else if points.strictly_between(shared, far, near)? {',
     '                    } else if false && points.strictly_between(shared, far, near)? {', ARR_TESTS),
    ('splits are free', ARR,
     'pub const SPLIT_SCRATCH_BYTES: usize = 64;',
     'pub const SPLIT_SCRATCH_BYTES: usize = 0;', ARR_TESTS),
]

probe_lib.main(MUTANTS, timeout=3600, reject_uncompiled=True)
