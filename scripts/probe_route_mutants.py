"""Mutation probe for route visibility (#187).

Each mutant drops one condition a visible segment must meet -- the proper
crossing test -- and must turn a test red.

Equivalent mutants, deliberately not listed:

- testing a stretch along an obstacle edge by its midpoint instead of
  accepting it. The midpoint of a stretch on the boundary is on the
  boundary, which counts as inside.
- not cutting the segment at obstacle vertices lying on it. Since #189
  every such vertex is a graph vertex whose sectors are checked on each
  side of travel (scripts/probe_route_sector_mutants.py), and a gap
  outside the region beyond it is a sector that is not free; the cuts
  stay as a second, independent guard.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = "crates/algorithms/planar/route/src/lib.rs"
TESTS = ["-p", "axiolid-route"]

MUTANTS = [
    ('proper crossings ignored', R, '        if crosses(a, b, *p, *q)? {\n            return Ok(false);\n        }', '        if crosses(a, b, *p, *q)? && false {\n            return Ok(false);\n        }', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
