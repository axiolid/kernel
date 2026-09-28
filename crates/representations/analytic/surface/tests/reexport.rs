//! `BSplineSurface` still resolves from both of its published paths, with
//! every published field (docs/architecture/semver-exceptions.toml).

use axiolid_core::Point3;
use axiolid_curve::KnotSpec;

fn net() -> Vec<Vec<Point3>> {
    vec![
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
        vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)],
    ]
}

#[test]
fn both_published_paths_name_the_same_struct_with_every_field() {
    let a = axiolid_surface::BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        control_points: net(),
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![2, 2],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![2, 2],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    };
    let b: axiolid_surface::spline::BSplineSurface = a.clone();
    let c: axiolid_curve::BSplineSurface = b.clone();
    assert_eq!(a, c);
}
