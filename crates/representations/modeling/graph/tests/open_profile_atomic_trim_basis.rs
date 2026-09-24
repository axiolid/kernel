use axiolid_core::{Frame2, Vec2};
use axiolid_curve::{Circle2, Curve2, Ellipse2, Line2, Sinusoid2};
use axiolid_model::{
    CurveRelation, GeometryGraphBuilder, GraphError, OpenProfile, TrimSelector, TrimmingPreference,
};

fn wrap_trim(builder: &mut GeometryGraphBuilder, basis_curve: Curve2) -> Result<(), GraphError> {
    let basis = builder.push_value(basis_curve)?;
    let trimmed = builder.push_value(CurveRelation::Trimmed {
        basis,
        start: vec![TrimSelector::Parameter(0.0)],
        end: vec![TrimSelector::Parameter(1.0)],
        sense_agreement: true,
        preference: TrimmingPreference::Parameter,
    })?;
    builder.push_value(OpenProfile::new(trimmed)).map(|_| ())
}

fn frame() -> Frame2 {
    Frame2 {
        origin: Vec2::ZERO,
        x: Vec2::X,
        y: Vec2::Y,
    }
}

#[test]
fn open_profile_accepts_valid_trimmed_lines_and_conics() {
    for curve in [
        Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::X,
        }),
        Curve2::Circle(Circle2 {
            frame: frame(),
            radius: 2.0,
        }),
        Curve2::Ellipse(Ellipse2 {
            frame: frame(),
            semi_axis_x: 3.0,
            semi_axis_y: 2.0,
        }),
    ] {
        wrap_trim(&mut GeometryGraphBuilder::new(), curve).unwrap();
    }
}

#[test]
fn open_profile_rejects_malformed_atomic_trim_bases() {
    let malformed = [
        Curve2::Line(Line2 {
            origin: Vec2::splat(f64::NAN),
            direction: Vec2::X,
        }),
        Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::ZERO,
        }),
        Curve2::Circle(Circle2 {
            frame: frame(),
            radius: f64::NAN,
        }),
        Curve2::Circle(Circle2 {
            frame: frame(),
            radius: 0.0,
        }),
        Curve2::Ellipse(Ellipse2 {
            frame: frame(),
            semi_axis_x: -1.0,
            semi_axis_y: 2.0,
        }),
        Curve2::Ellipse(Ellipse2 {
            frame: Frame2 {
                origin: Vec2::ZERO,
                x: Vec2::X,
                y: Vec2::X,
            },
            semi_axis_x: 3.0,
            semi_axis_y: 2.0,
        }),
    ];

    for curve in malformed {
        let error = wrap_trim(&mut GeometryGraphBuilder::new(), curve).unwrap_err();
        assert!(matches!(
            error,
            GraphError::InvalidReferenceType {
                expected: "bounded open curve2",
                ..
            }
        ));
    }
}

/// A cylinder-rim wave (ADR 0071) is a valid trim basis when finite, and a
/// non-finite coefficient in any of its three terms is refused.
#[test]
fn open_profile_accepts_a_finite_wave_and_rejects_a_non_finite_one() {
    let wave = Sinusoid2 {
        mean: 1.0,
        cosine: 0.5,
        sine: -0.25,
    };
    wrap_trim(&mut GeometryGraphBuilder::new(), Curve2::Sinusoid(wave)).unwrap();
    for bad in [
        Sinusoid2 {
            mean: f64::NAN,
            ..wave
        },
        Sinusoid2 {
            cosine: f64::INFINITY,
            ..wave
        },
        Sinusoid2 {
            sine: f64::NAN,
            ..wave
        },
    ] {
        let error = wrap_trim(&mut GeometryGraphBuilder::new(), Curve2::Sinusoid(bad)).unwrap_err();
        assert!(matches!(
            error,
            GraphError::InvalidReferenceType {
                expected: "bounded open curve2",
                ..
            }
        ));
    }
}
