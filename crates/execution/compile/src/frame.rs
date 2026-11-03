//! The frame two placed operands share (#356).
//!
//! Each `Instance` is placed on its own: every vertex of a source is rounded
//! by the `f64` transform. Two faces that lie on one plane in exact
//! arithmetic (a column's and its corbel's under one turned placement) are
//! therefore a rounding residue apart once placed, and the snap, which lands
//! only on axis-aligned planes, cannot close them (ADR 0096). Operands
//! placed by the same linear part are axis-aligned alike in their own
//! frame, where the faces are exactly coplanar, or a residue apart that the
//! snap closes. These helpers decide when two transforms share a frame and
//! where the tool then lies in it; the boolean itself is cut by
//! `ReferenceMeshCompiler::boolean_in_shared_frame`.
//!
//! The linear parts must be equal bit for bit: parts that differ in the last
//! bits (a rotation composed another way) place the same face on slightly
//! different planes, and no frame holds both exactly (TODO(#391)).

use axiolid_core::{Transform3, Vec3};

/// Whether two transforms have the same linear part, bit for bit (`-0.0`
/// read as `0.0`), finite and invertible.
pub(crate) fn shared_linear(a: Transform3, b: Transform3) -> bool {
    let columns = a.matrix3.to_cols_array();
    let determinant = a.matrix3.determinant();
    columns.iter().all(|v| v.is_finite())
        && columns == b.matrix3.to_cols_array()
        && determinant.is_finite()
        && determinant != 0.0
}

/// The tool's offset in the subject's frame: `subject.matrix3^-1 *
/// (tool.translation - subject.translation)`; `Vec3::ZERO` exactly when the
/// translations are equal bit for bit; `None` when not finite.
pub(crate) fn relative_offset(subject: Transform3, tool: Transform3) -> Option<Vec3> {
    let same = subject
        .translation
        .to_array()
        .iter()
        .zip(tool.translation.to_array())
        .all(|(a, b)| *a == b);
    if same {
        return Some(Vec3::ZERO);
    }
    let offset = subject.matrix3.inverse() * (tool.translation - subject.translation);
    offset.is_finite().then_some(offset)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn turn(angle: f64, at: Vec3) -> Transform3 {
        Transform3::from_translation(at) * Transform3::from_rotation_z(angle)
    }

    #[test]
    fn bitwise_equal_linear_parts_are_shared() {
        let a = turn(0.3, Vec3::new(1.0, 2.0, 3.0));
        let b = turn(0.3, Vec3::new(-5.0, 0.0, 7.0));
        assert!(shared_linear(a, b));
        assert!(shared_linear(a, a));
    }

    #[test]
    fn a_negative_zero_is_shared() {
        let a = Transform3::from_scale(Vec3::new(1.0, 1.0, 1.0));
        let mut b = a;
        b.matrix3.x_axis.y = -0.0;
        assert!(shared_linear(a, b));
    }

    #[test]
    fn a_rotation_composed_another_way_is_not_shared() {
        let a = Transform3::from_rotation_z(0.3);
        let b = Transform3::from_rotation_z(0.15) * Transform3::from_rotation_z(0.15);
        assert_ne!(a.matrix3, b.matrix3, "the test needs different bits");
        assert!(!shared_linear(a, b));
    }

    #[test]
    fn a_singular_linear_part_is_not_shared() {
        let flat = Transform3::from_scale(Vec3::new(1.0, 1.0, 0.0));
        assert!(!shared_linear(flat, flat));
        let mut nan = Transform3::IDENTITY;
        nan.matrix3.x_axis.x = f64::NAN;
        assert!(!shared_linear(nan, nan));
    }

    #[test]
    fn equal_translations_give_a_zero_offset_exactly() {
        let t = turn(0.7, Vec3::new(6.0e5, 5.6e6, 204.3));
        assert_eq!(relative_offset(t, t), Some(Vec3::ZERO));
    }

    #[test]
    fn the_offset_is_in_the_subjects_frame() {
        let site = Vec3::new(6.0e5, 5.6e6, 204.3);
        let tl = turn(0.7, site);
        let tr = Transform3 {
            matrix3: tl.matrix3,
            translation: tl.translation + tl.matrix3 * Vec3::new(0.0, 0.0, 3.5),
        };
        let d = relative_offset(tl, tr).expect("finite");
        assert!((d - Vec3::new(0.0, 0.0, 3.5)).abs().max_element() < 1e-9);
        let p = Vec3::new(0.25, -0.5, 1.0);
        let by_subject = tl.transform_point3(p + d);
        let by_tool = tr.transform_point3(p);
        assert!(
            (by_subject - by_tool).abs().max_element() < 1e-15 * site.abs().max_element(),
            "{by_subject:?} {by_tool:?}"
        );
    }
}
