"""Mutation probe for perspective hidden lines (#387, ADR 0107).

Each mutant breaks one part of how a perspective view is drawn: the exact
determinants of the sign questions (image side, plane in front, near-plane
lift, edge-on extent and segment), the depth and crop bounds, the homogeneous
vertex and its enclosure, the rounding of an endpoint, the cap query, the
view validation, or the flag that carries the projection through the
constraints. The hand scenes, the eye-ray oracle and the perspective
conformance fixtures must turn a test red for each.
"""
import probe_lib

HL = "crates/algorithms/discrete/hidden-lines/src/"
R = HL + "recipe.rs"
V = HL + "view.rs"
O = HL + "occlusion.rs"
CP = HL + "caps.rs"
OUT = HL + "output.rs"
PR = HL + "prepare.rs"
PV = HL + "provider.rs"
ED = HL + "edges.rs"
PA = HL + "param.rs"
BR = HL + "broad.rs"

TESTS = [
    "-p", "axiolid-hidden-lines",
    "--test", "perspective_lines",
    "--test", "perspective_oracle",
    "--test", "perspective",
]

BROAD = ["-p", "axiolid-hidden-lines", "--lib", "broad"]

ONE = "&Dyadic::from_parts(1.into(), 0)"

MUTANTS = [
    # The sign questions, as determinants of vectors from the eye.
    ("image side: the sign is flipped", R,
     "Kind::ImageSide { p, q, neg } => signed(neg, image_orient(verts, p, q, at)),",
     "Kind::ImageSide { p, q, neg } => signed(!neg, image_orient(verts, p, q, at)),"),
    ("image orientation: the weight column is second", R,
     "            T::from_dyadic(&exact[3]),\n            T::from_dyadic(&exact[0]),\n            T::from_dyadic(&exact[1]),\n        ]\n    };\n    dot3(&row(a)",
     "            T::from_dyadic(&exact[0]),\n            T::from_dyadic(&exact[3]),\n            T::from_dyadic(&exact[1]),\n        ]\n    };\n    dot3(&row(a)"),
    ("plane in front: P - t0 reversed", R,
     "Lift::Own => sub3(&v(at), &v(t0)),",
     "Lift::Own => sub3(&v(t0), &v(at)),"),
    ("near lift: t0 is not scaled by the depth", R,
     "sub3(&scale3(&T::from_f64(near), &p), &scale3(&p[2], &v(t0)))",
     "sub3(&scale3(&T::from_f64(near), &p), &v(t0))"),
    ("near lift: P is not scaled by the near depth", R,
     "sub3(&scale3(&T::from_f64(near), &p), &scale3(&p[2], &v(t0)))",
     "sub3(&p, &scale3(&p[2], &v(t0)))"),
    ("edge-on extent from: the cross product is reversed", R,
     "dot3(&cross3(&c, &v(at)), &cross3(&c, &d))",
     "dot3(&cross3(&v(at), &c), &cross3(&c, &d))"),
    ("edge-on extent to: the cross product is reversed", R,
     "dot3(&cross3(&v(at), &d), &cross3(&c, &d))",
     "dot3(&cross3(&d, &v(at)), &cross3(&c, &d))"),
    ("edge-on segment above: the norm is dropped", R,
     "s.sub(&dot3(&normal, &normal))",
     "s"),
    ("edge-on segment above: the norm is added", R,
     "s.sub(&dot3(&normal, &normal))",
     "s.add(&dot3(&normal, &normal))"),
    ("edge-on segment at the near plane: the depth term is added", R,
     ".sub(&dot3(&normal, &normal).mul(&v(at)[2]))",
     ".add(&dot3(&normal, &normal).mul(&v(at)[2]))"),
    ("edge-on segment at the near plane: the near factor is dropped", R,
     "                T::from_f64(near)\n                    .mul(&s)\n                    .sub(",
     "                s\n                    .sub("),
    # Depth and crop bounds.
    ("below near: the sign is flipped", R,
     "Kind::BelowNear(near) => T::from_f64(near).sub(&v(at)[2]),",
     "Kind::BelowNear(near) => v(at)[2].sub(&T::from_f64(near)),"),
    ("above far: the sign is flipped", R,
     "Kind::AboveFar(far) => v(at)[2].sub(&T::from_f64(far)),",
     "Kind::AboveFar(far) => T::from_f64(far).sub(&v(at)[2]),"),
    ("region: the sides are swapped", R,
     "                if upper {\n                    scaled.sub(&coordinate)",
     "                if !upper {\n                    scaled.sub(&coordinate)"),
    ("region: the bound is not scaled by the weight", R,
     "let scaled = T::from_f64(bound).mul(&T::from_dyadic(&exact[3]));",
     "let scaled = T::from_f64(bound);"),
    # Orientation and image difference of the vertices.
    ("orientation: the sign is flipped", R,
     "return image_orient::<T>(self.verts, self.a, self.b, self.c).sign();",
     "return image_orient::<T>(self.verts, self.a, self.b, self.c).neg().sign();"),
    ("image difference: the weights are added", R,
     ".sub(&T::from_dyadic(&t[3]).mul(&T::from_dyadic(&f[axis])))",
     ".add(&T::from_dyadic(&t[3]).mul(&T::from_dyadic(&f[axis])))"),
    ("a perspective constraint is read as orthographic", R,
     "        if self.perspective {\n            return self.eye_value(verts, at);\n        }",
     "        if false {\n            return self.eye_value(verts, at);\n        }"),
    # The homogeneous vertex and its enclosure.
    ("the weight is the depth, not its negative", V,
     "let h = w.neg();",
     "let h = w.clone();"),
    ("the weight is one", V,
     "let h = w.neg();",
     "let h = one();"),
    ("the image enclosure is not divided by the weight", V,
     "x.enclosure().quotient(h_enc)",
     "x.enclosure()"),
    ("the image enclosure divides by a negative weight", V,
     "            if h_enc.lo() > 0.0 {",
     "            if h_enc.lo() > 0.0 || h_enc.hi() < 0.0 {"),
    ("a perspective view builds orthographic vertices", V,
     "ViewProjection::Perspective => Self::perspective(p, &view.frame),",
     "ViewProjection::Perspective => Self::new(p, &view.frame),"),
    # Occlusion and the cap.
    ("the near plane's lift has the orthographic sign", O,
     "neg: clockwise != scene.perspective,",
     "neg: clockwise,"),
    ("a scene's constraints are orthographic", O,
     "Constraint::in_view(kind, strict, self.perspective)",
     "Constraint::new(kind, strict)"),
    ("the cap query keeps the orthographic depth range", CP,
     "        w_lo: if scene.perspective {\n            f64::NEG_INFINITY\n        } else {\n            near\n        },",
     "        w_lo: near,"),
    ("the cap crossing has the opposite sign", CP,
     "                    neg: !clockwise,\n                    lift: Lift::Near(near),",
     "                    neg: clockwise,\n                    lift: Lift::Near(near),"),
    # Endpoints: the exact rational, rounded once.
    ("an endpoint at the start ignores the weight", OUT,
     "Param::Zero => round_ratio(coord(start, axis), weight(start)),",
     "Param::Zero => round_ratio(coord(start, axis), " + ONE + "),"),
    ("an endpoint at the end ignores the weight", OUT,
     "Param::One => round_ratio(coord(end, axis), weight(end)),",
     "Param::One => round_ratio(coord(end, axis), " + ONE + "),"),
    ("a root's denominator ignores the weights", OUT,
     "let denominator = a.mul(weight(end)).sub(&b.mul(weight(start)));",
     "let denominator = a.sub(&b);"),
    ("a root's denominator swaps the weights", OUT,
     "let denominator = a.mul(weight(end)).sub(&b.mul(weight(start)));",
     "let denominator = a.mul(weight(start)).sub(&b.mul(weight(end)));"),
    # An edge-on triangle that contains the eye.
    ("an edge-on triangle never contains the eye", O,
     "let eye_inside = scene.perspective && eye_in_triangle(verts, tri, &mut solver.escalations)?;",
     "let eye_inside = false && eye_in_triangle(verts, tri, &mut solver.escalations)?;"),
    ("an orthographic edge-on triangle is tested for the eye", O,
     "let eye_inside = scene.perspective && eye_in_triangle(verts, tri, &mut solver.escalations)?;",
     "let eye_inside = eye_in_triangle(verts, tri, &mut solver.escalations)?;"),
    ("an edge-on triangle always contains the eye", O,
     "let eye_inside = scene.perspective && eye_in_triangle(verts, tri, &mut solver.escalations)?;",
     "let eye_inside = scene.perspective;"),
    ("the eye is in the triangle despite opposite signs", O,
     "    Ok(!(positive && negative))",
     "    Ok(!(positive && negative) || positive)"),
    ("a zero sign puts the eye outside the triangle", O,
     "            Ordering::Equal => {}\n        }\n    }\n    Ok(",
     "            Ordering::Equal => negative = true,\n        }\n    }\n    Ok("),
    ("the eye in the triangle keeps the nearer-than test", O,
     "        if eye_inside {\n            solver.insert(&mut in_front, within)?;",
     "        if false {\n            solver.insert(&mut in_front, within)?;"),
    # The broad phase: the hierarchy's depth is scaled in the build and the probe.
    ("the probe's depth is not scaled", BR,
     "                query.w_lo * self.depth_scale,",
     "                query.w_lo,", BROAD),
    ("the probe's upper depth is not scaled", BR,
     "                query.w_hi * self.depth_scale,",
     "                query.w_hi,", BROAD),
    ("the depth scale is negative", BR,
     "    if image.is_finite() && depth.is_finite() && scale.is_finite() && scale > 0.0 {\n        scale",
     "    if image.is_finite() && depth.is_finite() && scale.is_finite() && scale > 0.0 {\n        -scale", BROAD),
    ("an orthographic view scales the depth", BR,
     "        let depth_scale = if perspective {",
     "        let depth_scale = if true {", BROAD),
    # Validation and the flag through the pipeline.
    ("a near depth of zero is accepted", PR,
     "if !view.near.is_some_and(|near| near < 0.0) {",
     "if !view.near.is_some_and(|near| near <= 0.0) {"),
    ("a missing near depth is accepted", PR,
     "if !view.near.is_some_and(|near| near < 0.0) {",
     "if view.near.is_some_and(|near| !(near < 0.0)) {"),
    ("the scene is orthographic", PV,
     "perspective: request.view.projection == ViewProjection::Perspective,",
     "perspective: false,"),
    ("the triangle orientation is orthographic", ED,
     "                c: v[2],\n                perspective,",
     "                c: v[2],\n                perspective: false,"),
    ("the perturbed image difference is orthographic", PA,
     "perspective: constraint.perspective,",
     "perspective: false,"),
    ("the near cut test is orthographic", PR,
     "Constraint::in_view(Kind::BelowNear(near), false, perspective)",
     "Constraint::new(Kind::BelowNear(near), false)"),
]

probe_lib.main(MUTANTS, tests=TESTS, timeout=1800, baseline=True, reject_uncompiled=True)
