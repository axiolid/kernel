"""Mutation probe for B-rep faces on surface relations (#411, ADR 0103).

Each mutant weakens one decision of `crates/execution/compile/src/
surface_relation.rs` (the parameterisation, the closed-form derivative
bounds, the periods, the domain and the refusals), the view's arm for a
revolution, or the refusal of an edge use without a pcurve.
`tests/relation_faces.rs`, the module's unit tests, `tests/ruled_surface.rs`
or the tessellator's own tests must fail on it. A build error also counts as
killed.

Not listed, as not observable (no fixture, not equivalent):
- the `ROUNDING` inflation of a bound and the midpoint distance's rounding
  term: a relative 1e-12 widening of a bound that no fixture measures.
- the `DOMAIN_SLACK` clamp of a read inside the slack: the clamped point is
  the evaluated one to the last bit.
- the sort of the ruled sides' breaks: the smooth test reads them as a set.
- `Support::same` comparing two relations as always different: it only makes
  the edge pass locate shared samples again, with the same result.
- the `finite` check of a bound: the derivative bounds of the supported
  curve families are finite for finite data.
"""
import probe_lib

S = "crates/execution/compile/src/surface_relation.rs"
V = "crates/execution/compile/src/brep/view.rs"
I = "crates/execution/compile/src/brep/inversion.rs"
R = "crates/execution/compile/src/ruled.rs"
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--lib", "surface_relation"],
    ["-p", "axiolid-mesh-compile", "--test", "relation_faces"],
    ["-p", "axiolid-mesh-compile", "--test", "ruled_surface"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::"],
]

MUTANTS = [
    # Revolution bounds.
    ("the revolution's distance from the axis drops the width term", S,
     "                    + 0.5 * (v1 - v0) * perp\n",
     "                    + 0.0 * (v1 - v0) * perp\n"),
    ("the revolution's |S_uu| bound is 0", S,
     "                    duu: reach,\n",
     "                    duu: 0.0,\n"),
    ("the revolution's |S_uv| bound is 0", S,
     "                    duv: perp,\n",
     "                    duv: 0.0,\n"),
    ("the revolution's |S_vv| bound is 0", S,
     "                    dvv: d.second,\n",
     "                    dvv: 0.0,\n"),
    ("a line profile takes its whole speed as the speed normal to the axis", S,
     "                    Curve3::Line(line) => {\n"
     "                        (line.direction - *axis * line.direction.dot(*axis)).length()\n"
     "                    }\n",
     "                    Curve3::Line(_) => d.first,\n"),
    # Parameterisation.
    ("the rotation's sine term changes sign", S,
     "    along + (w - along) * u.cos() + a.cross(w) * u.sin()\n",
     "    along + (w - along) * u.cos() - a.cross(w) * u.sin()\n"),
    ("the axis is not normalised", S,
     "        let axis = axis_direction.normalize_or_zero();\n",
     "        let axis = axis_direction;\n"),
    ("an extrusion reads v in units of the unit direction", S,
     "            Self::Extrusion { curve, direction } => Ok(curve.point(u)? + *direction * v),\n",
     "            Self::Extrusion { curve, direction } => {\n"
     "                Ok(curve.point(u)? + direction.normalize() * v)\n"
     "            }\n"),
    # Extrusion bounds.
    ("the extrusion's |S_uu| bound is 0", S,
     "                    duu: d.second,\n",
     "                    duu: 0.0,\n"),
    # Periods.
    ("the revolution's period in u is None", S,
     "            Self::Revolution { curve, .. } => (Some(TAU), curve.is_conic().then_some(TAU)),\n",
     "            Self::Revolution { curve, .. } => (None, curve.is_conic().then_some(TAU)),\n"),
    ("a conic profile's period in v is None", S,
     "            Self::Revolution { curve, .. } => (Some(TAU), curve.is_conic().then_some(TAU)),\n",
     "            Self::Revolution { .. } => (Some(TAU), None),\n"),
    # Refusals and domain.
    ("the swept-family refusal admits a polyline", S,
     "            _ => return Err(refused(SWEPT_CURVE_FAMILY)),\n",
     "            Curve3::Polyline(_) => None,\n"
     "            _ => return Err(refused(SWEPT_CURVE_FAMILY)),\n"),
    ("the domain-slack check never refuses", S,
     "        if !(t >= domain.start - slack && t <= domain.end + slack) {\n",
     "        if false {\n"),
    ("the oracle returns None for a box reaching past the domain instead of clamping it", S,
     "                if lo > domain.end + slack || hi < domain.start - slack {\n",
     "                if lo > domain.end + slack || hi < domain.start - slack\n"
     "                    || lo < domain.start || hi > domain.end\n"
     "                {\n"),
    ("a box wholly within the slack past the end is refused", S,
     "                if lo > domain.end + slack || hi < domain.start - slack {\n",
     "                if lo > domain.end || hi < domain.start - slack {\n"),
    ("a box with a break inside is smooth", S,
     "            !self.breaks[axis].iter().any(|&b| b > lo && b < hi)\n",
     "            true\n"),
    # Ruled bounds.
    ("the ruled |S_uv| bound drops its second-derivative term", S,
     "                    + 0.5 * h * (a2 + b2)\n",
     "                    + 0.0 * h * (a2 + b2)\n"),
    ("the ruled |S_v| bound drops its Lipschitz term", S,
     "                let dv = (b - a).length() + 0.5 * h * duv + ROUNDING * (a.length() + b.length());\n",
     "                let dv = (b - a).length() + ROUNDING * (a.length() + b.length());\n"),
    ("the ruled |S_u| and |S_uu| bounds ignore the v blend", S,
     "                let (du, duu) = (blend(a1, b1), blend(a2, b2));\n",
     "                let (du, duu) = (a1, a2);\n"),
    # The view and the tessellator.
    ("the view's arm for a revolution is guarded off", V,
     "            node @ SurfaceRelation::Revolution {\n"
     "                swept_curve,\n"
     "                axis_origin,\n"
     "                axis_direction,\n"
     "            },\n"
     "        )) => Ok(Some(Support::Relation {\n",
     "            node @ SurfaceRelation::Revolution {\n"
     "                swept_curve,\n"
     "                axis_origin,\n"
     "                axis_direction,\n"
     "            },\n"
     "        )) if false => Ok(Some(Support::Relation {\n"),
    ("a ruled box across a break keeps the Lipschitz twist", S,
     "                let duv = if smooth_box {\n"
     "                    lipschitz.min(a1 + b1)\n"
     "                } else {\n"
     "                    a1 + b1\n"
     "                };\n",
     "                let duv = lipschitz;\n"),
    ("a ruled side's break is not mapped through its interval", R,
     "                    .map(|t| (t - start) / (end - start))\n",
     "                    .map(|t| t)\n"),
    ("the view's arm for a linear extrusion is guarded off", V,
     "                direction,\n"
     "            },\n"
     "        )) => Ok(Some(Support::Relation {\n",
     "                direction,\n"
     "            },\n"
     "        )) if false => Ok(Some(Support::Relation {\n"),
    ("the view's arm for a ruled surface is guarded off", V,
     "        Some(GeometryNode::SurfaceRelation(node @ SurfaceRelation::Ruled { first, second })) => {\n",
     "        Some(GeometryNode::SurfaceRelation(node @ SurfaceRelation::Ruled { first, second }))\n"
     "            if false =>\n"
     "        {\n"),
    ("the family of a revolution is named a surface", S,
     '            Self::Revolution { .. } => "surface of revolution",\n',
     '            Self::Revolution { .. } => "surface",\n'),
    ("an edge use without a pcurve on a relation is refused as degenerate", I,
     "            return Err(GeomError::UnsupportedInput {\n"
     "                backend: crate::BACKEND_ID,\n"
     "                operation: Operation::Tessellation,\n"
     "                input: RELATION_WITHOUT_PCURVE,\n"
     "            });\n",
     '            return Err(GeomError::Degenerate("a relation without a pcurve".into()));\n'),
]

probe_lib.main(MUTANTS, tests=TARGETS, command=probe_lib.CARGO_TEST + ("--all-features",), timeout=1800)
