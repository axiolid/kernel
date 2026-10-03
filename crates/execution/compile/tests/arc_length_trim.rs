//! `TrimSelector::ArcLength` on a sweep directrix (#239).
//!
//! The trim is resolved against the basis by arc length, in the trim's
//! sense, and compiles to the same solid as a parameter trim at the
//! parameter computed independently by hand.

use axiolid_contracts::{ExecutionOptions, GeomResult};
use axiolid_core::{Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::{BSplineCurve3, Curve3, KnotSpec, Line3};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    CurveRelation, GeometryGraphBuilder, GeometryNode, NodeId, SolidOperation, TrimSelector,
    TrimmingPreference,
};
use axiolid_profile::{Profile, RectangleProfile};

const AREA: Scalar = 0.4 * 0.6;

fn sweep(build: impl FnOnce(&mut GeometryGraphBuilder) -> NodeId) -> GeomResult<TriMesh> {
    let mut b = GeometryGraphBuilder::new();
    let profile = b
        .push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: 0.4,
                y: 0.6,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )))
        .unwrap();
    let directrix = build(&mut b);
    let swept = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::FixedReferenceSweep {
                profile,
                directrix,
                reference_direction: Vec3::Z,
                parameter_range: None,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![swept]).unwrap();
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(
        &graph,
        swept,
        &ExecutionOptions::new(Tolerance::MILLIMETRE),
    )
}

fn volume(mesh: &TriMesh) -> Scalar {
    volume_properties(mesh, Tolerance::MILLIMETRE)
        .expect("a swept solid is closed")
        .signed_volume
        .abs()
}

fn trim(
    b: &mut GeometryGraphBuilder,
    basis: NodeId,
    start: TrimSelector,
    end: TrimSelector,
    sense_agreement: bool,
) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
        basis,
        start: vec![start],
        end: vec![end],
        sense_agreement,
        preference: TrimmingPreference::Parameter,
    }))
    .unwrap()
}

/// A line whose parameter runs at speed 2: arc length 10 is parameter 5.
fn fast_line(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push(GeometryNode::Curve3(Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::new(2.0, 0.0, 0.0),
    })))
    .unwrap()
}

#[test]
fn an_arc_length_trim_of_a_line_sweeps_that_length_in_either_sense() {
    for sense in [true, false] {
        let mesh = sweep(|b| {
            let line = fast_line(b);
            trim(
                b,
                line,
                TrimSelector::ArcLength(0.0),
                TrimSelector::ArcLength(10.0),
                sense,
            )
        })
        .unwrap();
        // Read as a parameter, the same trim would sweep 20.
        let v = volume(&mesh);
        assert!((v - AREA * 10.0).abs() <= 1e-9, "sense {sense}: volume {v}");
        // Measured from parameter 0 in the trim's sense.
        let xs = mesh.positions.iter().map(|p| p.x);
        let (lo, hi) = xs.fold((Scalar::INFINITY, Scalar::NEG_INFINITY), |(lo, hi), x| {
            (lo.min(x), hi.max(x))
        });
        let expected = if sense { (0.0, 10.0) } else { (-10.0, 0.0) };
        assert!(
            (lo - expected.0).abs() <= 1e-9 && (hi - expected.1).abs() <= 1e-9,
            "sense {sense}: spans {lo}..{hi}"
        );
    }
}

/// `y = a x^3` over `x in [0, reach]` as an exact cubic Bezier in z = 0.
fn cubic(a: Scalar, reach: Scalar) -> Curve3 {
    Curve3::BSpline(BSplineCurve3 {
        degree: 3,
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(reach / 3.0, 0.0, 0.0),
            Point3::new(2.0 * reach / 3.0, 0.0, 0.0),
            Point3::new(reach, a * reach.powi(3), 0.0),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::PiecewiseBezier,
    })
}

/// Abscissa at arc length `length` along `y = a x^3`, by Newton on the
/// binomial series of `int sqrt(1 + 9 a^2 x^4) dx`.
fn cubic_abscissa(a: Scalar, length: Scalar) -> Scalar {
    let series = |x: Scalar| {
        let z = 9.0 * a * a * x.powi(4);
        let (mut binomial, mut power, mut sum) = (1.0, 1.0, 0.0);
        for n in 0..200 {
            sum += binomial * power / (4.0 * n as Scalar + 1.0);
            binomial *= (0.5 - n as Scalar) / (n as Scalar + 1.0);
            power *= z;
        }
        x * sum
    };
    let mut x = length;
    for _ in 0..50 {
        x -= (series(x) - length) / (1.0 + 9.0 * a * a * x.powi(4)).sqrt();
    }
    x
}

#[test]
fn an_arc_length_trim_of_a_cubic_matches_the_parameter_trim_at_its_series_abscissa() {
    let (radius, length) = (30.0, 20.0);
    let a = 1.0 / (6.0 * radius * length);
    let by_length = sweep(|b| {
        let basis = b.push(GeometryNode::Curve3(cubic(a, length))).unwrap();
        trim(
            b,
            basis,
            TrimSelector::ArcLength(0.0),
            TrimSelector::ArcLength(length),
            true,
        )
    })
    .unwrap();
    let t = cubic_abscissa(a, length) / length;
    let by_parameter = sweep(|b| {
        let basis = b.push(GeometryNode::Curve3(cubic(a, length))).unwrap();
        trim(
            b,
            basis,
            TrimSelector::Parameter(0.0),
            TrimSelector::Parameter(t),
            true,
        )
    })
    .unwrap();
    let (u, v) = (volume(&by_length), volume(&by_parameter));
    assert!((u - v).abs() <= 1e-9 * v, "{u} against {v}");
    // Pappus: a section square to a plane path sweeps area times length.
    assert!((u / (AREA * length) - 1.0).abs() <= 5e-3, "volume {u}");
}

#[test]
fn an_arc_length_trim_is_refused_by_name_where_it_cannot_be_measured() {
    // Past the end of a bounded basis.
    let error = sweep(|b| {
        let basis = b.push(GeometryNode::Curve3(cubic(1e-3, 10.0))).unwrap();
        trim(
            b,
            basis,
            TrimSelector::ArcLength(0.0),
            TrimSelector::ArcLength(100.0),
            true,
        )
    })
    .unwrap_err();
    assert!(error.to_string().contains("does not resolve"), "{error}");

    // On a relation basis there is no single curve to measure along.
    let error = sweep(|b| {
        let line = fast_line(b);
        let inner = trim(
            b,
            line,
            TrimSelector::Parameter(0.0),
            TrimSelector::Parameter(8.0),
            true,
        );
        trim(
            b,
            inner,
            TrimSelector::ArcLength(0.0),
            TrimSelector::ArcLength(4.0),
            true,
        )
    })
    .unwrap_err();
    assert!(error.to_string().contains("arc-length selector"), "{error}");
}
