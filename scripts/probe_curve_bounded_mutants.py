"""Mutation probe for curve-bounded planes (#192).

Each mutant fills what must be refused, maps boundaries wrongly, or meshes
the plane facing the wrong way, and must turn a test red.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
B = "crates/execution/compile/src/bounded.rs"
TESTS = ["-p", "axiolid-mesh-compile", "--test", "curve_bounded"]

MUTANTS = [
    ('implicit outer boundary not refused', B, '    if implicit_outer {', '    if false && implicit_outer {', TESTS),
    ('non-planar basis taken as a plane', B, '        Some(GeometryNode::Surface(axiolid_surface::Surface::Plane(plane))) => plane.frame,', '        Some(GeometryNode::Surface(axiolid_surface::Surface::Plane(plane))) => plane.frame,\n        Some(GeometryNode::Surface(axiolid_surface::Surface::Cylinder(c))) => c.frame,', TESTS),
    ('clockwise outer loop kept', B, '        if index == 0 && signed_area(&ring) < 0.0 {', '        if false && index == 0 && signed_area(&ring) < 0.0 {', TESTS),
    ('parameters mapped with the axes swapped', B, '    let map = |p: Point2| frame.origin + frame.x * p.x + frame.y * p.y;', '    let map = |p: Point2| frame.origin + frame.y * p.x + frame.x * p.y;', TESTS),
    ('open polylines closed silently', B, '        !closed && !meets', '        false && !closed && !meets', TESTS),
    ('3D boundaries off the parameter plane flattened', B, '    if let Some(off) = points.iter().find(|p| p.z.abs() > linear) {', '    if let Some(off) = points.iter().find(|p| p.z.abs() > linear && false) {', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
