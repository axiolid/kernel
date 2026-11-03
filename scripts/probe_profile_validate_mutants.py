"""Mutation probe for the exact profile validator: each mutant must be caught."""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "crates/algorithms/construction/construct/src"

MUTANTS = [
    ("outer segment_of_edge replaced by the identity", "profile_validate.rs",
     "        ArcRingRef::Outer => lowered[0].segment_of_edge[edge],",
     "        ArcRingRef::Outer => edge,"),
    ("hole segment_of_edge read from the wrong hole", "profile_validate.rs",
     "        ArcRingRef::Hole(k) => lowered[k + 1].segment_of_edge[edge],",
     "        ArcRingRef::Hole(k) => lowered[k].segment_of_edge[edge],"),
    ("split arc maps to one edge only", "contour_lower.rs",
     "std::iter::repeat_n(index, lowered.vertices.len())",
     "std::iter::repeat_n(index, 1)"),
    ("broken contour does not short-circuit", "profile_validate.rs",
     "    if broken {\n        return report;\n    }",
     "    if false {\n        return report;\n    }"),
    ("max joint gap not tracked", "profile_validate.rs",
     "                report.max_joint_gap = report.max_joint_gap.max(done.max_gap);",
     "                let _ = done.max_gap;"),
    ("inner joint gaps not tracked", "contour_lower.rs",
     "            max_gap = max_gap.max(gap);\n        } else {",
     "        } else {"),
    ("closing joint gap not tracked", "contour_lower.rs",
     "        max_gap = max_gap.max(gap);\n    }\n\n    Ok(LoweredContour {",
     "    }\n\n    Ok(LoweredContour {"),
    ("not-lowered profile swallowed", "profile_validate.rs",
     "                findings: vec![ProfileFinding::NotLowered(error)],",
     "                findings: {\n                    let _ = error;\n                    vec![]\n                },"),
    ("holes not lowered", "profile_validate.rs",
     "                    .map(|(k, hole)| (ContourRef::Hole(k), hole)),",
     "                    .filter(|_| false)\n                    .map(|(k, hole)| (ContourRef::Hole(k), hole)),"),
    ("hole indices off by one", "profile_validate.rs",
     "(ContourRef::Hole(k), hole)",
     "(ContourRef::Hole(k + 1), hole)"),
    ("contact contours not mapped", "profile_validate.rs",
     "            contours: contact.rings.map(ContourRef::from),",
     "            contours: [ContourRef::Outer; 2],"),
    ("gap names the following segment", "contour_lower.rs",
     "                    segment: index - 1,",
     "                    segment: index,"),
    ("self-contact segments left unsorted", "profile_validate.rs",
     "                segments.sort_unstable();",
     "                segments.reverse();"),
]

probe_lib.main(MUTANTS, tests=[["-p", "axiolid-construct", "--test", "profile_validate"],
                           ["-p", "axiolid-construct", "--lib"]], root=SRC, timeout=900)
