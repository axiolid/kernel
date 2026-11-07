//! Carriers with skewed frames, for the surface-bound and lifted-pcurve
//! tests (#414).

#![allow(dead_code)]

use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::{BSplineSurface, Carrier, KnotSpec, RuledCarrier, TorusCarrier};

/// A frame whose axes are neither unit nor orthogonal.
pub fn skewed() -> Frame3 {
    Frame3 {
        origin: Point3::new(0.3, -0.2, 0.5),
        x: Vec3::new(1.0, 0.2, -0.1),
        y: Vec3::new(0.15, 0.9, 0.3),
        z: Vec3::new(0.1, -0.2, 1.1),
    }
}

/// A rotated orthonormal frame.
pub fn rotated() -> Frame3 {
    let (s, c) = 0.7_f64.sin_cos();
    let (s2, c2) = 0.4_f64.sin_cos();
    let x = Vec3::new(c, s, 0.0);
    let y0 = Vec3::new(-s, c, 0.0);
    let y = y0 * c2 + Vec3::new(0.0, 0.0, s2);
    let z = x.cross(y);
    Frame3 {
        origin: Point3::new(0.4, -0.3, 0.2),
        x,
        y,
        z,
    }
}

/// A rational biquadratic patch over the unit square, bumped in the middle.
pub fn patch() -> BSplineSurface {
    let weights = vec![
        vec![1.0, 1.3, 1.0],
        vec![0.9, 1.6, 1.1],
        vec![1.0, 1.2, 1.0],
    ];
    BSplineSurface {
        u_degree: 2,
        v_degree: 2,
        control_points: (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| {
                        let z = if i == 1 && j == 1 {
                            0.6
                        } else {
                            0.1 * i as Scalar
                        };
                        Point3::new(i as Scalar * 0.5, j as Scalar * 0.5 + 0.1 * i as Scalar, z)
                    })
                    .collect()
            })
            .collect(),
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![3, 3],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![3, 3],
        weights: Some(weights),
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    }
}

/// A carrier with a parameter box inside its smooth part.
pub struct Case {
    pub name: &'static str,
    pub carrier: Carrier,
    pub u: (Scalar, Scalar),
    pub v: (Scalar, Scalar),
}

/// Every carrier family, in skewed or rotated frames.
pub fn cases() -> Vec<Case> {
    let ruled = |x_radius, y_radius, slope| {
        Carrier::Ruled(RuledCarrier {
            frame: skewed(),
            x_radius,
            y_radius,
            slope,
        })
    };
    vec![
        Case {
            name: "plane",
            carrier: Carrier::Plane(skewed()),
            u: (-1.0, 1.0),
            v: (-0.7, 0.9),
        },
        Case {
            name: "cylinder",
            carrier: ruled(0.6, 0.6, 0.0),
            u: (0.3, 2.0),
            v: (-0.5, 0.7),
        },
        Case {
            name: "elliptical cylinder",
            carrier: ruled(0.8, 0.4, 0.0),
            u: (-1.0, 0.8),
            v: (-0.2, 0.9),
        },
        Case {
            name: "cone",
            carrier: ruled(0.5, 0.5, 0.4),
            u: (0.1, 1.9),
            v: (0.0, 0.8),
        },
        Case {
            name: "sphere",
            carrier: Carrier::Sphere {
                frame: skewed(),
                radius: 0.7,
            },
            u: (0.2, 1.3),
            v: (-0.5, 0.6),
        },
        Case {
            name: "torus",
            carrier: Carrier::Torus(TorusCarrier {
                frame: rotated(),
                major_radius: 0.5,
                minor_radius: 0.2,
            }),
            u: (0.0, 1.2),
            v: (-1.0, 1.0),
        },
        // A negative minor radius: the ring is widest at `v = pi`.
        Case {
            name: "torus with a negative minor radius",
            carrier: Carrier::Torus(TorusCarrier {
                frame: rotated(),
                major_radius: 2.0,
                minor_radius: -1.0,
            }),
            u: (0.0, 0.6),
            v: (2.9, 3.4),
        },
        Case {
            name: "spline",
            carrier: Carrier::Spline(Box::new(patch())),
            u: (0.1, 0.9),
            v: (0.15, 0.85),
        },
    ]
}
